/*
  simulations of the experiments performed in

  Ten Cate, A., et al. "Particle imaging velocimetry experiments and 
  lattice-Boltzmann simulations on a single sphere settling under gravity." 
  Physics of Fluids (1994-present) 14.11 (2002): 4012-4025.
  
  parameters for the experiments in SI units : (rho_f, mu_f, u_inf)
  E1: 970 ; 0.373 ; 0.038
  E2: 965 ; 0.212 ; 0.6
  E3: 962 ; 0.113 ; 0.091
  E4: 960 ; 0.058 ; 0.128
  
 */


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
// #include "lammps.h"
// #include "input.h"
// #include "library.h"
// #include "library_cfd_coupling.h"

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
  
  MultiScalarField3D<T> tmp(lattice.getNx(),lattice.getNy(),lattice.getNz());
  
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
                          "SolidFraction",1.);
  vtkOut.writeData<float>(*computeExternalScalar(lattice,DESCRIPTOR<T>::ExternalField::particleIdBeginsAt),
                          "PartId",1.);
  vtkOut.writeData<float>(*computeExternalScalar(lattice,
                                                 DESCRIPTOR<T>::ExternalField::hydrodynamicForceBeginsAt),
                          "fx",1.);
  vtkOut.writeData<float>(*computeExternalScalar(lattice,
                                                 DESCRIPTOR<T>::ExternalField::hydrodynamicForceBeginsAt+1),
                          "fy",1.);
  vtkOut.writeData<float>(*computeExternalScalar(lattice,
                                                 DESCRIPTOR<T>::ExternalField::hydrodynamicForceBeginsAt+2),
                          "fz",1.);

  // for(plint i=0;i<DESCRIPTOR<T>::ExternalField::numScalars;i++){
  //   std::stringstream s;
  //   s << i;
  //   vtkOut.writeData<float>(*computeExternalScalar(lattice,i),s.str().c_str(),1.);
  // }

  pcout << "wrote " << fname << std::endl;
}

void writePopulation(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, Box3D slice, plint iPop, plint iter)
{
  std::stringstream ofname;
  ofname << global::directories().getOutputDir();
  ofname << "f_" << std::setfill('0') << std::setw(2) << iPop 
         << "_" << std::setfill('0') << std::setw(8) << iter << ".dat";
  
  plb_ofstream of(ofname.str().c_str());
  
  of << *computePopulation(lattice,slice,iPop);
 
  of.close();

}

void writeExternal(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, Box3D slice, plint which, 
                   char* const prefix, plint iter)
{
  std::stringstream ofname;
  ofname << global::directories().getOutputDir();
  ofname << prefix << "_" << std::setfill('0') << std::setw(8) << iter << ".dat";
  
  plb_ofstream of(ofname.str().c_str());
  
  of << *computeExternalScalar(lattice,which,slice);
 
  of.close();

}
int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);


    plint N(0);
    T uMax(0.),rho_f(0.),mu_f(0.),v_inf(0.);
    std::string outDir;

    try {
        global::argv(1).read(N);
        global::argv(2).read(uMax);
        global::argv(3).read(rho_f);
        global::argv(4).read(mu_f);
        global::argv(5).read(v_inf);
        global::argv(6).read(outDir);
    } catch(PlbIOException& exception) {
        pcout << "Error the parameters are wrong. The structure must be :\n";
        pcout << "1 : N\n";
        pcout << "2 : uMax\n";
        pcout << "3 : rho_f\n";
        pcout << "4 : mu_f\n";
        pcout << "5 : expected v_inf\n";
        pcout << "6 : outDir\n";
        exit(1);
    }


    std::string lbOutDir(outDir), demOutDir(outDir);
    lbOutDir.append("tmp/"); demOutDir.append("post/");
    global::directories().setOutputDir(lbOutDir);
    LiggghtsCouplingWrapper wrapper(argv,global::mpi().getGlobalCommunicator());
    
    wrapper.execFile("in.lbdem");

    const T g = 9.81;
    const T nu_f = mu_f/rho_f;
    const T lx = 0.1, ly = 0.1, lz = 0.16;

    T r_ = 0.015/2.;

    PhysUnits3D<T> units(2.*r_,v_inf,nu_f,lx,ly,lz,N,uMax,rho_f);

    IncomprFlowParam<T> parameters(units.getLbParam());

    
    
    const T maxT = (T)8.0;
    const T vtkT = 0.01;
    const T logT = 0.001;

    const plint maxSteps = units.getLbSteps(maxT);
    const plint vtkSteps = max<plint>(units.getLbSteps(vtkT),1);
    const plint logSteps = max<plint>(units.getLbSteps(logT),1);

    writeLogFile(parameters, "3D sedimenting sphere");

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


    lattice.periodicity().toggleAll(false);

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();


    T dt_phys = units.getPhysTime(1);
    pcout << "omega: " << parameters.getOmega() << "\n" 
          << "dt_phys: " << dt_phys << "\n"
          << "Re : " << parameters.getRe() << "\n"
          << "vtkSteps: " << vtkSteps << "\n"
          << "grid size: " 
          << parameters.getNx() << " "
          << parameters.getNy() << " "
          << parameters.getNz() << std::endl;
    
    lattice.initialize();

    T dt_dem = dt_phys/(T)demSubsteps;

    wrapper.setVariable("t_step",dt_dem);
    wrapper.setVariable("dmp_stp",vtkSteps*demSubsteps);
    wrapper.setVariable("dmp_dir",demOutDir);


    wrapper.execFile("in2.lbdem");
    wrapper.runUpto(demSubsteps-1);

    clock_t start = clock();    


    // Loop over main time iteration.
    for (plint iT=0; iT<maxSteps; ++iT) {

      setSpheresOnLatticeNew(lattice,wrapper,units,false);
    
      // if(iT%vtkSteps == 0)
      //   writeVTK(lattice,parameters,units,iT);

      lattice.collideAndStream();

      getForcesFromLatticeNew(lattice,wrapper,units);

      wrapper.run(demSubsteps);
      
      if(iT%logSteps == 0){
        clock_t end = clock();
        T time = difftime(end,start)/((T)CLOCKS_PER_SEC);
        T mlups = ((T) (lattice.getNx()*lattice.getNy()*lattice.getNz()*units.getLbSteps(logT)))/time/1e6;
        pcout << "time: " << time << " " ;
        pcout << "calculating at " << mlups << " MLU/s" << std::endl;
        start = clock();

      }

      
    }

    delete boundaryCondition;
}
