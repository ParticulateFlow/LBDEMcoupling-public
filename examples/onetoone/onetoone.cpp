
#include "palabos3D.h"
#include "palabos3D.hh"

#include "plb_ib.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>

// necessary LAMMPS/LIGGGHTS includes
#include "lammps.h"
#include "input.h"
#include "library.h"
#include "library_cfd_coupling.h"
#include "comm.h"
#include "atom.h"
#include "periodicPressureFunctionals3D.h"
#include "liggghtsCouplingWrapper.h"

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::ImmersedBoundaryD3Q19Descriptor
#define DYNAMICS IBdynamics<T, DESCRIPTOR>(parameters.getOmega())

void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
              IncomprFlowParam<T> const& parameters,
              PhysUnits3D<T> const& units, plint iter)
{
  
  T p_fact = units.getPhysForce(1)/pow(units.getPhysLength(1),2)/3.;
  
  std::string fname(createFileName("vtk", iter, 6));
  
  VtkImageOutput3D<T> vtkOut(fname, units.getPhysLength(1));
  vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", units.getPhysVel(1));
  vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", units.getPhysVel(1));  
  vtkOut.writeData<float>(*computeDensity(lattice), "density",units.getPhysDensity(1)); 
  
  MultiScalarField3D<T> p(*computeDensity(lattice));
  subtractInPlace(p,1.);
  vtkOut.writeData<float>(p,"pressure",p_fact ); 
  
 
  vtkOut.writeData<float>(*computeExternalScalar(lattice,DESCRIPTOR<T>::ExternalField::volumeFractionBeginsAt),
                          "SolidFraction",1);
  vtkOut.writeData<float>(*computeExternalScalar(lattice,DESCRIPTOR<T>::ExternalField::particleIdBeginsAt),
                          "PartId",1);
  vtkOut.writeData<float>(*computeExternalScalar(lattice,DESCRIPTOR<T>::ExternalField::hydrodynamicForceBeginsAt),
                          "fx",1);
  vtkOut.writeData<float>(*computeExternalScalar(lattice,DESCRIPTOR<T>::ExternalField::hydrodynamicForceBeginsAt+1),
                          "fy",1);
  vtkOut.writeData<float>(*computeExternalScalar(lattice,DESCRIPTOR<T>::ExternalField::hydrodynamicForceBeginsAt+2),
                          "fz",1);
  pcout << "wrote " << fname << std::endl;
}

void writeGif(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,plint iT)
{
    const plint imSize = 600;
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNy();
    const plint nz = lattice.getNz();
    Box3D slice(0, nx-1, (ny-1)/2, (ny-1)/2, 0, nz-1);
    //Box3D slice(0, nx-1, 0, ny-1, (nz-1)/2, (nz-1)/2);
    ImageWriter<T> imageWriter("leeloo.map");
    std::string fname(createFileName("u", iT, 6));
    imageWriter.writeScaledGif(fname,
                               *computeVelocityNorm(lattice, slice),
                               imSize, imSize);
    pcout << "wrote " << fname << std::endl;
}

void writePopulation(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,plint iPop, plint iT)
{
  std::stringstream fname_stream;
  fname_stream << global::directories().getOutputDir()
               << "f_" << setfill('0') << setw(2) << iPop << "_"
               << setfill('0') << setw(8) << iT << ".dat";

  plb_ofstream ofile(fname_stream.str().c_str());
  Box3D domain(lattice.getNx()/2,lattice.getNx()/2,
               0,lattice.getNy(),0,lattice.getNz());
  ofile << *computePopulation(lattice, domain, iPop);
}
void writeExternal(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint which, char const *prefix, plint iT)
{
  std::stringstream fname_stream;
  fname_stream << global::directories().getOutputDir()
               << prefix << setfill('0') << setw(8) << iT << ".dat";

  plb_ofstream ofile(fname_stream.str().c_str());
  Box3D domain(lattice.getNx()/2,lattice.getNx()/2,
               0,lattice.getNy(),0,lattice.getNz());
  ofile << *computeExternalScalar(lattice,which,domain);
}


class LatticeDecomposition {
public:
  LatticeDecomposition(plint nx_, plint ny_, plint nz_, LAMMPS_NS::LAMMPS *lmp_)
    : nx(nx_),ny(ny_),nz(nz_),lmp(*lmp_),
      blockStructure(0),threadAttribution(0)
  {
    npx = lmp.comm->procgrid[0];
    npy = lmp.comm->procgrid[1];
    npz = lmp.comm->procgrid[2];
    for(plint i=0;i<=npx;i++)
      xVal.push_back(round( lmp.comm->xsplit[i]*(T)nx ));
    for(plint i=0;i<=npy;i++)
      yVal.push_back(round( lmp.comm->ysplit[i]*(T)ny ));
    for(plint i=0;i<=npz;i++)
      zVal.push_back(round( lmp.comm->zsplit[i]*(T)nz ));

    blockStructure = new SparseBlockStructure3D(Box3D(xVal[0], xVal.back()-1, 
                                                yVal[0], yVal.back()-1, 
                                                zVal[0], zVal.back()-1) );
    
    threadAttribution = new ExplicitThreadAttribution();

    for (plint iX=0; iX<xVal.size()-1; ++iX) {
      for (plint iY=0; iY<yVal.size()-1; ++iY) {
        for (plint iZ=0; iZ<zVal.size()-1; ++iZ) {
          plint id = blockStructure->nextIncrementalId();
          blockStructure->addBlock (Box3D( xVal[iX], xVal[iX+1]-1, yVal[iY],
                                           yVal[iY+1]-1, zVal[iZ], zVal[iZ+1]-1 ),
                                id);
          threadAttribution->addBlock(id,(plint)lmp.comm->grid2proc);
        }
      }
    }
    
  }
  
  ~LatticeDecomposition()
  {
    if(blockStructure != 0) delete blockStructure;
    if(threadAttribution != 0) delete threadAttribution;
  }


  SparseBlockStructure3D getBlockDistribution()
  {
    return SparseBlockStructure3D(*blockStructure);
  }
  ExplicitThreadAttribution* getThreadAttribution()
  {
    return new ExplicitThreadAttribution(*threadAttribution);
  }
private:
  plint nx,ny,nz;
  LAMMPS_NS::LAMMPS &lmp;
  plint npx,npy,npz;
  std::vector<plint> xVal, yVal, zVal;
  SparseBlockStructure3D *blockStructure;
  ExplicitThreadAttribution *threadAttribution;
};

int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);

    T uMax;

    plint N;
    
    T nu_f,d_part,v_frac, v_inf;

    std::string outDir;
    
    try {
        global::argv(1).read(N);
    } catch(PlbIOException& exception) {
        pcout << exception.what() << endl;
        pcout << "Command line arguments:\n";
        pcout << "1 : N per particle diameter\n";
        exit(1);
    }

    std::string lbOutDir(outDir), demOutDir(outDir);
    lbOutDir.append("tmp/"); demOutDir.append("post/");
    global::directories().setOutputDir(lbOutDir);

    const T rho_f = 1000;

    char **argv_lmp = 0;
    argv_lmp = new char*[1];
    argv_lmp[0] = argv[0];

    LiggghtsCouplingWrapper wrapper(argv,global::mpi().getGlobalCommunicator());
    
    wrapper.setVariable("r_part",d_part/2);
    wrapper.setVariable("v_frac",v_frac);
    
    // domain.lbdem should set up everything until the create_box command
    wrapper.execFile("domain.lbdem");
    // pair.lbdem should contain pair style and integrate fix
    wrapper.execFile("pair.lbdem");
    // geom.lbdem contains all geometry definitions, stl imports, gravity, whatsoever
    wrapper.execFile("geom.lbdem");
    // ins.lbdem contains all insert commands
    wrapper.execFile("ins.lbdem");
    wrapper.run(1);

    // wrapper.allocateVariables();
    
    PhysUnits3D<T> units(1.,1.,1.,1.,1.,1.,N,0.02,1.);

    IncomprFlowParam<T> parameters(units.getLbParam());

    LatticeDecomposition lDec(parameters.getNx(),parameters.getNy(),parameters.getNz(),
                              wrapper.lmp);
    
    SparseBlockStructure3D blockStructure = lDec.getBlockDistribution();
    ExplicitThreadAttribution* threadAttribution = lDec.getThreadAttribution();
    plint demSubsteps = 10;
    
    const T maxT = 1.;
    const T vtkT = 100;
    const T logT = 0.0000001;

    const plint maxSteps = 10;//units.getLbSteps(maxT);
    const plint vtkSteps = max<plint>(units.getLbSteps(vtkT),1);
    const plint logSteps = max<plint>(units.getLbSteps(logT),1);

    writeLogFile(parameters, "rect channel showcase");

    plint nx = parameters.getNx(), ny = parameters.getNy(), nz = parameters.getNz()-1;

    T dt_phys = units.getPhysTime(1);
    pcout << "omega: " << parameters.getOmega() << "\n" 
          << "dt_phys: " << dt_phys << "\n"
          << "maxT: " << maxT << " | maxSteps: " << maxSteps << "\n"
          << "v_inf: " << v_inf << "\n"
          << "Re : " << parameters.getRe() << "\n"
          << "vtkSteps: " << vtkSteps << "\n"
          << "grid size: " << nx << " " << ny << " " << nz << std::endl;

    plint envelopeWidth = 1;
    MultiBlockLattice3D<T, DESCRIPTOR> 
      lattice (MultiBlockManagement3D (blockStructure, threadAttribution, envelopeWidth ),
               defaultMultiBlockPolicy3D().getBlockCommunicator(),
               defaultMultiBlockPolicy3D().getCombinedStatistics(),
               defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
               new DYNAMICS );
    
    lattice.initialize();

    T dt_dem = dt_phys/(T)demSubsteps;
    std::stringstream cmd;
    // cmd << "variable t_step equal " << dt_dem;
    // pcout << cmd.str() << std::endl;
    // wrapper.execCommand(cmd);
    // cmd.str("");

    wrapper.setVariable("t_step",dt_dem);
    wrapper.setVariable("dmp_stp",vtkSteps*demSubsteps);
    wrapper.setVariable("dmp_dir",demOutDir);

    wrapper.execFile("ts_thermo_dmp.lbdem");
    wrapper.runUpto(demSubsteps-1);

    clock_t start = clock();
    clock_t loop = clock();
    clock_t end = clock(); 
    // Loop over main time iteration.
    for (plint iT=0; iT<=maxSteps; ++iT) {
      
      plint nlocal = wrapper.lmp->atom->nlocal;
      plint nghost = wrapper.lmp->atom->nghost;
      std::cout << global::mpi().getRank() << " " << nlocal << " " << nghost << std::endl;
      for(plint i=0;i<nlocal+nghost;i++)
        std::cout << global::mpi().getRank() << " | "
                  << i << " | "
                  << wrapper.lmp->atom->x[i][0] << " "
                  << wrapper.lmp->atom->x[i][1] << " "
                  << wrapper.lmp->atom->x[i][2] << std::endl;


      //wrapper.dataFromLiggghts();

      bool initWithVel = false;
      //setSpheresOnLattice(lattice,wrapper,units,initWithVel);
      

      if(iT%vtkSteps == 0 && iT > 0) // LIGGGHTS does not write at timestep 0
        writeVTK(lattice,parameters,units,iT);

      lattice.collideAndStream();

      //getForcesFromLattice(lattice,wrapper,units);

      //wrapper.dataToLiggghts();
      wrapper.run(demSubsteps);


      if(iT%logSteps == 0){
        end = clock();
        T time = difftime(end,loop)/((T)CLOCKS_PER_SEC);
        T totaltime = difftime(end,start)/((T)CLOCKS_PER_SEC);
        T mlups = ((T) (lattice.getNx()*lattice.getNy()*lattice.getNz()*logSteps))/time/1e6;
        pcout << "time: " << time << " " ;
        pcout << "calculating at " << mlups << " MLU/s"
              << " | total time running: " << totaltime << std::endl;
        loop = clock();
      }
    }
    T totaltime = difftime(end,start)/((T)CLOCKS_PER_SEC);
    T totalmlups = ((T) (lattice.getNx()*lattice.getNy()*lattice.getNz()*(maxSteps+1)))/totaltime/1e6;
    pcout << " ********************** \n"
          << "total time: " << totaltime
          << " calculating at " << totalmlups << " MLU/s" << std::endl;

}
