
#include "palabos3D.h"
#include "palabos3D.hh"

#include "periodicPressureFunctionals3D.h"

#include "plb_ib.h"
#include <ctime>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D3Q19Descriptor
#define DYNAMICS BGKdynamics<T, DESCRIPTOR>(parameters.getOmega())


void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
              IncomprFlowParam<T> const& parameters,
              PhysUnits3D<T> const& units, plint iter)
{
  
  MultiScalarField3D<T> tmp(lattice.getNx(),lattice.getNy(),lattice.getNz());
  
  T p_fact = units.getPhysForce(1)/pow(units.getPhysLength(1),2) / 3.;
  
  std::string fname(createFileName("vtk", iter, 6));
  
  VtkImageOutput3D<T> vtkOut(fname, units.getPhysLength(1));
  // vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", units.getPhysVel(1));
  // vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", units.getPhysVel(1));  
  // vtkOut.writeData<float>(*computeDensity(lattice), "density",units.getPhysDensity(1)); 

  vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", 1.);
  vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", 1.);  
  vtkOut.writeData<float>(*computeDensity(lattice), "density", 1.); 
  
  MultiScalarField3D<T> p(*computeDensity(lattice));
  subtractInPlace(p,1.);
  vtkOut.writeData<float>(p,"pressure",p_fact ); 

  for(plint i=0;i<DESCRIPTOR<T>::q;i++){
    std::stringstream s;
    s << "f " << i;
    vtkOut.writeData<float>(*computePopulation(lattice,i),s.str().c_str(),1.);
  }

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





int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    plint N;
    T uMax;
    T deltaP;
    T lFactor;
    try {
        global::argv(1).read(N);
        global::argv(2).read(deltaP);
        global::argv(3).read(uMax);
        global::argv(4).read(lFactor);
    } catch(PlbIOException& exception) {
        pcout << exception.what() << endl;
        pcout << "Command line arguments:\n";
        pcout << "1 : N\n";
        pcout << "2 : deltaP\n";
        pcout << "3 : uMax\n";
        pcout << "4 : lengthFactor\n";
        exit(1);
    }

    const T rho_f = 1000;

    const T nu_f = 1e-3;

    const T lx = 0.2, ly = 0.02, lz = 0.2*lFactor;
    
    T gradP = deltaP/lz;
    T lx_eff = lx;//*(T)(N-1)/(T)N;
    T u_phys = gradP*lx_eff*lx_eff/(nu_f*8*rho_f);
    
    PhysUnits3D<T> units(lx,u_phys,nu_f,lx,ly,lz,N,uMax,rho_f);

    IncomprFlowParam<T> parameters(units.getLbParam());

    const T maxT = (T)1;
    const T vtkT = 0.5;
    const T logT = 0.5;

    const plint maxSteps = units.getLbSteps(maxT);
    const plint vtkSteps = max<plint>(units.getLbSteps(vtkT),1);
    const plint logSteps = max<plint>(units.getLbSteps(logT),1);

    writeLogFile(parameters, "3D sedimenting sphere");

    plint nx = parameters.getNx(), ny = parameters.getNy(), nz = parameters.getNz()-1;

    // lattice setup
    MultiBlockLattice3D<T, DESCRIPTOR> lattice (nx,ny,nz,new DYNAMICS );

    lattice.periodicity().toggle(1,true);
    lattice.periodicity().toggle(2,true);

    Box3D inlet(0,nx-1,0,ny-1,0,0), outlet(0,nx-1,0,ny-1,nz-1,nz-1);
    
    T deltaRho = units.getLbRho(deltaP);//*(T)(nz)/(T)(nz+1);

    initializeAtEquilibrium( lattice, lattice.getBoundingBox(), 
                             PressureGradient<T>(1.+0.5*deltaRho,1.-0.5*deltaRho, nz, 2) );

    T dt_phys = units.getPhysTime(1);
    pcout << "omega: " << parameters.getOmega() << "\n" 
          << "dt_phys: " << dt_phys << "\n"
          << "u_phys: " << u_phys << "\n"
          << "Re : " << parameters.getRe() << "\n"
          << "deltaRho : " << deltaRho << "\n"
          << "vtkSteps: " << vtkSteps << "\n"
          << "grid size: " << nx << " " << ny << " " << nz << std::endl;
    
    clock_t start = clock();    
    for (plint iT=0; iT<maxSteps; ++iT) {

      if(iT%vtkSteps == 0)
        writeVTK(lattice,parameters,units,iT);


      pcerr << units.getPhysPress(computeAverageDensity(lattice,inlet)) << " " 
            << units.getPhysPress(computeAverageDensity(lattice,outlet)) << std::endl;

      PeriodicPressureManager<T,DESCRIPTOR> ppm(lattice,1.+0.5*deltaRho,1.-0.5*deltaRho,inlet,outlet,2,1,-1);
      ppm.preColl(lattice);
      lattice.collideAndStream();
      ppm.postColl(lattice);

      if(iT%logSteps == 0){
        clock_t end = clock();
        T time = difftime(end,start)/((T)CLOCKS_PER_SEC);
        T mlups = ((T) (lattice.getNx()*lattice.getNy()*lattice.getNz()*units.getLbSteps(logT)))/time/1e6;
        pcout << "time: " << time << " " ;
        pcout << "calculating at " << mlups << " MLU/s" << std::endl;
        start = clock();
      }
    }
}
