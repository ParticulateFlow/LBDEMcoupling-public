
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

int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);

    T uMax;

    plint N;
    
    T nu_f,d_part,v_frac, v_inf;

    std::string outDir;
    
    try {
        global::argv(1).read(d_part);
        global::argv(2).read(N);
        global::argv(3).read(v_frac);
        global::argv(4).read(nu_f);
        global::argv(5).read(v_inf);
        global::argv(6).read(uMax);
        global::argv(7).read(outDir);
    } catch(PlbIOException& exception) {
        pcout << exception.what() << endl;
        pcout << "Command line arguments:\n";
        pcout << "1 : d_part\n";
        pcout << "2 : N per particle diameter\n";
        pcout << "3 : particle volume fraction\n";
        pcout << "4 : nu_fluid\n";
        pcout << "5 : estimated v_inf\n";
        pcout << "6 : uMax\n";
        pcout << "7 : outDir\n";
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
    
    wrapper.execFile("in.lbdem");
    wrapper.allocateVariables();

    T g = 9.81;

    const T lx = 1., ly = 1., lz = 2.;

    wrapper.dataFromLiggghts();

    T r_ = wrapper.r[0][0];
    T rho_s = 1100.;
    T m = r_*r_*r_*4./3.*3.14*rho_s;

    // T v_inf_calc = 2.*(rho_s-rho_f)/rho_f*g*r_*r_/9./nu_f; // stokes
    T v_inf_calc = sqrt(4./3.*0.44*(rho_s-rho_f)/rho_f*g*2.*r_); // something else

    pcout << "v_inf: " << v_inf << " m: " << m << " r: " << r_ << std::endl;
    
    PhysUnits3D<T> units(2.*r_,v_inf,nu_f,lx,ly,lz,N,uMax,rho_f);

    IncomprFlowParam<T> parameters(units.getLbParam());

    plint demSubsteps = 10;
    
    const T maxT = ceil(3.*lz/v_inf);// (T)1.;
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

    MultiBlockLattice3D<T, DESCRIPTOR> lattice (nx,ny,nz, new DYNAMICS );


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

    // cmd << "variable dmp_dir string " << demOutDir;
    // pcout << cmd.str() << std::endl;
    // wrapper.execCommand(cmd);
    // cmd.str("");

    wrapper.execFile("in2.lbdem");
    wrapper.execCommand("run 9 upto");

    clock_t start = clock();    


    // Loop over main time iteration.
    for (plint iT=0; iT<=maxSteps; ++iT) {

      wrapper.dataFromLiggghts();

      bool initWithVel = false;
      setSpheresOnLattice(lattice,wrapper,units,initWithVel);
      

      if(iT%vtkSteps == 0 && iT > 0) // LIGGGHTS does not write at timestep 0
        writeVTK(lattice,parameters,units,iT);

      lattice.collideAndStream();

      getForcesFromLattice(lattice,wrapper,units);

      wrapper.dataToLiggghts();
      wrapper.execCommand("run 10");


      if(iT%logSteps == 0){
        clock_t end = clock();
        T time = difftime(end,start)/((T)CLOCKS_PER_SEC);
        T mlups = ((T) (lattice.getNx()*lattice.getNy()*lattice.getNz()*logSteps))/time/1e6;
        pcout << "time: " << time << " " ;
        pcout << "calculating at " << mlups << " MLU/s" << std::endl;
        start = clock();
      }
    }

}
