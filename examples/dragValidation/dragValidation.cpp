
#include "palabos3D.h"
#include "palabos3D.hh"

#include "plb_ib.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include <ctime>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::ImmersedBoundaryD3Q19Descriptor
#define DYNAMICS IBdynamics<T, DESCRIPTOR>(parameters.getOmega())

const T pi = 4.*atan(1.);
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
  
  
  applyProcessingFunctional
    (new GetExternalScalarFunctional3D<T,DESCRIPTOR,T>(DESCRIPTOR<T>::ExternalField::volumeFractionBeginsAt),
     lattice.getBoundingBox(), lattice, tmp);
  vtkOut.writeData<float>(tmp, "SolidFraction", 1.);
  
  pcout << "wrote " << fname << std::endl;
}


void writeToAscii(std::string fnamePrefix, MultiScalarField3D<T> &data, int iT)
{
  std::string fname(global::directories().getOutputDir());
  fname.append(createFileName(fnamePrefix, iT, 6));
  fname.append(".txt");
  plb_ofstream ofile(fname.c_str());
  ofile << setprecision(10) << data << endl;
  ofile.close();

}


void writeAscii(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,int iT)
{
  
  MultiScalarField3D<T> tmp(lattice.getNx(),lattice.getNy(),lattice.getNz());
                            
  applyProcessingFunctional
    (new GetExternalScalarFunctional3D<T,DESCRIPTOR,T>(DESCRIPTOR<T>::ExternalField::hydrodynamicForceBeginsAt),
     lattice.getBoundingBox(), lattice, tmp);
  
  writeToAscii("fx",tmp,iT);
  
  applyProcessingFunctional
    (new GetExternalScalarFunctional3D<T,DESCRIPTOR,T>(DESCRIPTOR<T>::ExternalField::volumeFractionBeginsAt),
     lattice.getBoundingBox(), lattice, tmp);
  
  writeToAscii("fs",tmp,iT);
  writeToAscii("ux",*computeVelocityComponent(lattice,lattice.getBoundingBox(),0),iT);
}

int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    if (argc != 6) {
        pcout << "Error the parameters are wrong. The structure must be :\n";
        pcout << "1 : N\n";
        pcout << "2 : u0\n";
        pcout << "3 : uMax\n";
        pcout << "4 : maxT\n";
        pcout << "5 : vtkT\n";
        exit(1);
    }

    const plint N = atoi(argv[1]);
    const T u0 = atof(argv[2]);
    const T uMax = atof(argv[3]);
    const T maxT = atof(argv[4]);
    const T vtkT = atof(argv[5]);

    const T rho_f = 1000;

    const T lx = 0.6, ly = 0.2, lz = 0.2;

    plint *id = new plint;
    *id = 0;
    T **x, **v, **f, **t, **omega;
    T r = 0.02;
    x = new T*[1]; x[0] = new T[3];
    v = new T*[1]; v[0] = new T[3];
    omega = new T*[1]; omega[0] = new T[3];
    f = new T*[1]; f[0] = new T[3];
    t = new T*[1]; t[0] = new T[3];


    x[0][0] = 0.2; x[0][1] = 0.1; x[0][2] = 0.1;
    v[0][0] = 0; v[0][1] = 0; v[0][2] = 0;
    omega[0][0] = 0; omega[0][1] = 0; omega[0][2] = 0;
    
    f[0][0] = 0; f[0][1] = 0; f[0][2] = 0;
    t[0][0] = 0; t[0][1] = 0; t[0][2] = 0;


    const T nu_f = 1e-4;

    PhysUnits3D<T> units(2.*r,u0,nu_f,lx,ly,lz,N,uMax,rho_f);

    IncomprFlowParam<T> parameters(units.getLbParam());

    
    const T logT = 0.05;

    const plint maxSteps = units.getLbSteps(maxT);
    const plint vtkSteps = max<plint>(units.getLbSteps(vtkT),1);
    const plint logSteps = max<plint>(units.getLbSteps(logT),1);

    writeLogFile(parameters, "3D sedimenting sphere");

    const plint nx = parameters.getNx(), ny = parameters.getNy(), nz = parameters.getNz();

    MultiBlockLattice3D<T, DESCRIPTOR> lattice (nx,ny,nz,new DYNAMICS);

    plint envelopeWidth = 1;

    lattice.periodicity().toggleAll(false);
    lattice.periodicity().toggle(1,true);
    lattice.periodicity().toggle(2,true);

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
      = createLocalBoundaryCondition3D<T,DESCRIPTOR>();


    Box3D interior(1,nx-2,1,ny-2,1,nz-2);
    Box3D inlet(0,0,0,ny-1,0,nz-1),
      outlet(nx-1,nx-1,0,ny-1,0,nz-1);
    
    Array<T,3> u0_vec(uMax,0.,0.);

    boundaryCondition->setVelocityConditionOnBlockBoundaries ( lattice, inlet );
    boundaryCondition->setPressureConditionOnBlockBoundaries (lattice, outlet );
    
    setBoundaryDensity(lattice,outlet,1.);
    setBoundaryVelocity(lattice, inlet, u0_vec );
    
    initializeAtEquilibrium(lattice,lattice.getBoundingBox(),1.,u0_vec);

    lattice.periodicity().toggleAll(false);
    lattice.periodicity().toggle(1,true);
    lattice.periodicity().toggle(2,true);

    T Re = parameters.getRe();

    pcout << "omega: " << parameters.getOmega() << "\n" 
          << "Re : " << Re << "\n"
          << "grid size: " << nx << " " << ny << " " << nz << "\n"
          << "physical time step: " << units.getPhysTime(1) << std::endl;

    T c_drag_sn = 24./Re*(1. + 0.15*pow(Re,0.687));

    T **x_lb, **v_lb, **omega_lb;
    x_lb = new T*[1]; x_lb[0] = new T[3];
    v_lb = new T*[1]; v_lb[0] = new T[3];
    omega_lb = new T*[1]; omega_lb[0] = new T[3];
    T r_lb = units.getLbLength(r);

    for(plint i=0;i<3;i++){
      x_lb[0][i] = units.getLbLength(x[0][i]);
      v_lb[0][i] = units.getLbVel(v[0][i]);
      omega_lb[0][i] = units.getLbFreq(omega[0][i]);
    }

    applyProcessingFunctional
      (new SetSpheres3D<T,DESCRIPTOR>(x_lb,v_lb,omega_lb,&r_lb,id,1),
       lattice.getBoundingBox(), lattice);
    pcout << "spheres set" << std::endl;
    clock_t start = clock();    

    // Loop over main time iteration.
    for (plint iT=0; iT<maxSteps; ++iT) {
      for(int i=0;i<3;i++){
        f[0][i] = 0; t[0][i] = 0;
      }
      
      lattice.collideAndStream();

      if(iT%vtkSteps == 0)
        writeVTK(lattice,parameters,units,iT);

      SumForceTorque3D<T,DESCRIPTOR> funct(1,x);
      if(iT > 1) 
        applyProcessingFunctional(funct, lattice.getBoundingBox(), lattice);
      
      
      std::vector<double> &v_ft = funct.getForceTorque();
      
      for(int i=0;i<3;i++){
        f[0][i] = units.getPhysForce(v_ft[i]);
        t[0][i] = units.getPhysTorque(v_ft[i+3]);
      }
      

      T c_drag = 2.*f[0][0]/(rho_f*u0*u0*r*r*3.1416);
      
      if(iT%logSteps == 0){
        clock_t end = clock();
        T time = difftime(end,start)/((T)CLOCKS_PER_SEC);
        T mlups = ((T) (lattice.getNx()*lattice.getNy()*lattice.getNz()*units.getLbSteps(logT)))/time/1e6;
        pcout << "time: " << time << " " ;
        pcout << "calculating at " << mlups << " MLU/s" << std::endl;
        start = clock();
        pcerr << units.getPhysTime((T)iT) << " | ";
        pcerr << f[0][0] << " " << f[0][1] << " " << f[0][2] << " | ";
        pcerr << t[0][0] << " " << t[0][1] << " " << t[0][2] << " | ";
        pcerr << c_drag << " " << c_drag_sn << " | ";
        pcerr << sqrt(2.*f[0][0]/(rho_f*u0*u0*3.1416*c_drag_sn)) << std::endl;
      }
      
    }

    delete boundaryCondition;
}
