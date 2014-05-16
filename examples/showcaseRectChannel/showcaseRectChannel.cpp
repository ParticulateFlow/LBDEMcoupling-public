
#include "palabos3D.h"
#include "palabos3D.hh"

#include "plb_ib.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>

#include "periodicPressureFunctionals3D.h"
#include "liggghtsCouplingWrapper.h"
#include "latticeDecomposition.h"

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
  // vtkOut.writeData<float>(*computeExternalScalar(lattice,DESCRIPTOR<T>::ExternalField::particleIdBeginsAt),
  //                         "PartId",1);
  // vtkOut.writeData<float>(*computeExternalScalar(lattice,DESCRIPTOR<T>::ExternalField::hydrodynamicForceBeginsAt),
  //                         "fx",1);
  // vtkOut.writeData<float>(*computeExternalScalar(lattice,DESCRIPTOR<T>::ExternalField::hydrodynamicForceBeginsAt+1),
  //                         "fy",1);
  // vtkOut.writeData<float>(*computeExternalScalar(lattice,DESCRIPTOR<T>::ExternalField::hydrodynamicForceBeginsAt+2),
  //                         "fz",1);
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

int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);

    const T uMax = 0.02;

    plint N;
    T deltaP;
    
    std::string outDir;
    
    try {
        global::argv(1).read(N);
        global::argv(2).read(deltaP);
        global::argv(3).read(outDir);
    } catch(PlbIOException& exception) {
        pcout << exception.what() << endl;
        pcout << "Command line arguments:\n";
        pcout << "1 : N\n";
        pcout << "2 : deltaP\n";
        pcout << "3 : outDir\n";
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
    wrapper.execFile("in.lbdem");

    const T nu_f = 1e-3;

    const T lx = 0.8, ly = 0.2, lz = 0.2;

    T gradP = deltaP/lz;
    T lx_eff = lx;//*(T)(N-1)/(T)N;
    T u_phys = gradP*lx_eff*lx_eff/(nu_f*8*rho_f);

    
    PhysUnits3D<T> units(lz,u_phys,nu_f,lx,ly,lz,N,uMax,rho_f);

    IncomprFlowParam<T> parameters(units.getLbParam());

    LatticeDecomposition lDec(parameters.getNx(),parameters.getNy(),parameters.getNz(),
                              wrapper.lmp);
    
    SparseBlockStructure3D blockStructure = lDec.getBlockDistribution();
    ExplicitThreadAttribution* threadAttribution = lDec.getThreadAttribution();
    plint demSubsteps = 10;
    
    plint envelopeWidth = 1;
    MultiBlockLattice3D<T, DESCRIPTOR> 
      lattice (MultiBlockManagement3D (blockStructure, threadAttribution, envelopeWidth ),
               defaultMultiBlockPolicy3D().getBlockCommunicator(),
               defaultMultiBlockPolicy3D().getCombinedStatistics(),
               defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
               new DYNAMICS );

    const T maxT = (T)1000.;
    const T vtkT = 0.1;
    const T gifT = 100;
    const T logT = 0.000000002;

    const plint maxSteps = units.getLbSteps(maxT);
    const plint vtkSteps = max<plint>(units.getLbSteps(vtkT),1);
    const plint gifSteps = units.getLbSteps(gifT);
    const plint logSteps = max<plint>(units.getLbSteps(logT),1);

    writeLogFile(parameters, "rect channel showcase");

    plint nx = parameters.getNx(), ny = parameters.getNy(), nz = parameters.getNz()-1;

    lattice.periodicity().toggle(0,true);

    Box3D inlet(0,0,1,ny-2,1,nz-2), outlet(nx-1,nx-1,1,ny-2,1,nz-2);
    
    T deltaRho = units.getLbRho(deltaP);
    T gradRho = units.getLbRho(gradRho);
    // T rhoHi = 1.+0.5*deltaRho, rhoLo = 1.-0.5*deltaRho;
    T rhoHi = 1., rhoLo = 1.-deltaRho;

    // initializeAtEquilibrium( lattice, lattice.getBoundingBox(), 
    //                          PressureGradient<T>(rhoHi,rhoLo,nz,0) );
    initializeAtEquilibrium( lattice, lattice.getBoundingBox(), 
                             PoiseuilleProfileAndPressureGradient<T>(rhoHi,rhoLo,uMax,nx,ny,nz,0) );

    lattice.initialize();
    lattice.toggleInternalStatistics(false);

    T dt_phys = units.getPhysTime(1);
    pcout << "omega: " << parameters.getOmega() << "\n" 
          << "dt_phys: " << dt_phys << "\n"
          << "u_phys: " << u_phys << "\n"
          << "Re : " << parameters.getRe() << "\n"
          << "deltaRho : " << deltaRho << "\n"
          << "vtkSteps: " << vtkSteps << "\n"
          << "grid size: " << nx << " " << ny << " " << nz << std::endl;

    T dt_dem = dt_phys/(T)demSubsteps;
    wrapper.setVariable("t_step",dt_dem);
    wrapper.setVariable("dmp_stp",vtkSteps*demSubsteps);
    wrapper.setVariable("dmp_dir",demOutDir);


    wrapper.execFile("in2.lbdem");
    wrapper.runUpto(demSubsteps-1);

    clock_t start = clock();
    clock_t loop = clock();
    clock_t end = clock(); 

    // Loop over main time iteration.
    for (plint iT=0; iT<=maxSteps; ++iT) {

      static bool initWithVel = true;
      setSpheresOnLatticeNew(lattice,wrapper,units,initWithVel);
      if(initWithVel) initWithVel = false;

      // if(iT%vtkSteps == 0 && iT > 3000)
      if(iT%vtkSteps == 0 && iT > 0) // LIGGGHTS does not write at timestep 0
      // if(iT%vtkSteps == 0)
        writeVTK(lattice,parameters,units,iT);

      PeriodicPressureManager<T,DESCRIPTOR> ppm(lattice,rhoHi,rhoLo,inlet,outlet,0,1,-1);
      ppm.preColl(lattice);
      lattice.collideAndStream();
      ppm.postColl(lattice);

      getForcesFromLatticeNew(lattice,wrapper,units);

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
