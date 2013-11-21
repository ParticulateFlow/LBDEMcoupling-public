
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
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    std::string fname(createFileName("vtk", iter, 6));
    VtkImageOutput3D<T> vtkOut(fname, dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", dx/dt);   

    MultiScalarField3D<T> tmp(lattice.getNx(),lattice.getNy(),lattice.getNz());

    applyProcessingFunctional
      (new GetExternalScalarFunctional3D<T,DESCRIPTOR,T>(DESCRIPTOR<T>::ExternalField::volumeFractionBeginsAt),
       lattice.getBoundingBox(), lattice, tmp);
    vtkOut.writeData<float>(tmp, "SolidFraction", 1.);

    pcout << "wrote " << fname << std::endl;
}

void writeGif(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,int iT)
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

int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    if (argc != 3) {
        pcout << "Error the parameters are wrong. The structure must be :\n";
        pcout << "1 : N\n";
        pcout << "2 : rho_s/rho_f\n";
        exit(1);
    }

    const plint N = atoi(argv[1]);
    const T rho_fact = atof(argv[2]);
    const T uMax = 0.02;

    const T rho_f = 1000;

    plint *id = new plint;
    *id = 0;
    T **x, **v, **f, **t, **omega;
    T r = 0.01;
    T rho_s = rho_f*rho_fact;
    T m = 4./3.*r*r*r*3.14*rho_s;
    T I = 2./5.*m*r*r;
    x = new T*[1]; x[0] = new T[3];
    v = new T*[1]; v[0] = new T[3];
    omega = new T*[1]; omega[0] = new T[3];

    T g = 9.81;

    pcout << "m: " << m << " | I: " << I << std::endl;

    const T lx = 0.2, ly = 0.2, lz = 0.4;

    x[0][0] = lx/2; x[0][1] = ly/2; x[0][2] = lz-2*r;
    v[0][0] = 0; v[0][1] = 0; v[0][2] = 0;
    omega[0][0] = 0; omega[0][1] = 0; omega[0][2] = 0;

    T v_inf = sqrt(2.*m*g/(rho_f*r*r*3.14*0.47));

    pcout << "v_inf: " << v_inf << std::endl;


    const T nu_f = 1e-4;
    //    PhysUnits3D<T> units(lx,v_inf,nu_f,lx,ly,lz,N,uMax,rho_f);
    PhysUnits3D<T> units(2.*r,v_inf,nu_f,lx,ly,lz,N,uMax,rho_f);

    IncomprFlowParam<T> parameters(units.getLbParam());

    
    const T maxT = (T)10.0;
    const T vtkT = 0.05;
    const T gifT = 0.01;
    const T logT = 0.001;

    writeLogFile(parameters, "3D sedimenting sphere");

    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
        parameters.getNx(), parameters.getNy(), parameters.getNz(), 
        new DYNAMICS );

    // Use periodic boundary conditions.
    lattice.periodicity().toggle(2,true);

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();


    f = new T*[1]; f[0] = new T[3];
    t = new T*[1]; t[0] = new T[3];
    
    f[0][0] = 0; f[0][1] = 0; f[0][2] = 0;
    t[0][0] = 0; t[0][1] = 0; t[0][2] = 0;


    T dt_dem = units.getPhysTime(1);
    pcout << "omega: " << parameters.getOmega() << "\n" 
          << "dt_dem: " << dt_dem << "\n"
          << "Re : " << parameters.getRe() << std::endl;

    Array<T,3> x_(units.getLbLength(x[0][0]),
                  units.getLbLength(x[0][1]),
                  units.getLbLength(x[0][2]));
    Array<T,3> v_(units.getLbVel(v[0][0]),
                  units.getLbVel(v[0][1]),
                  units.getLbVel(v[0][2]));
    Array<T,3> omega_(units.getLbFreq(omega[0][0]),
                  units.getLbFreq(omega[0][1]),
                  units.getLbFreq(omega[0][2]));

    
    applyProcessingFunctional
      (new SetSphere3D<T,DESCRIPTOR>(x_,v_,omega_,units.getLbLength(r),*id),
       lattice.getBoundingBox(), lattice);
    

    clock_t start = clock();


    // Loop over main time iteration.
    //    for (plint iT=0; iT<parameters.nStep(maxT); ++iT) {
    pcout << units.getLbSteps(0.01) << std::endl;
    for (plint iT=0; iT<units.getLbSteps(maxT); ++iT) {

      v[0][0] += 0.5*dt_dem*f[0][0]/m;
      v[0][1] += 0.5*dt_dem*f[0][1]/m;
      v[0][2] += 0.5*dt_dem*(f[0][2]/m-g*(rho_fact-1));

      x[0][0] += v[0][0]*dt_dem + 0.5*f[0][0]*dt_dem*dt_dem/m;
      x[0][1] += v[0][1]*dt_dem + 0.5*f[0][1]*dt_dem*dt_dem/m;
      x[0][2] += v[0][2]*dt_dem + 0.5*dt_dem*dt_dem*(f[0][2]/m-g*(rho_fact-1));

      omega[0][0] += 0.5*dt_dem*t[0][0]/I;
      omega[0][1] += 0.5*dt_dem*t[0][1]/I;
      omega[0][2] += 0.5*dt_dem*t[0][2]/I;

      //      if(iT%parameters.nStep(vtkT) == 0)
      if(iT%units.getLbSteps(vtkT) == 0)
        writeVTK(lattice,parameters,iT);

      //      if(iT%parameters.nStep(gifT) == 0)
      if(iT%units.getLbSteps(gifT) == 0)
        writeGif(lattice,iT);
      
      for(int i=0;i<3;i++){
        f[0][i] = 0; t[0][i] = 0;
      }

      if(iT > 1) 
        applyProcessingFunctional(new SumForceTorque3D<T,DESCRIPTOR>(id,x,f,t),
                                  lattice.getBoundingBox(), lattice);
      for(int i=0;i<3;i++){
        f[0][i] = units.getPhysForce(f[0][i]);
        t[0][i] = units.getPhysTorque(t[0][i]);
      }

      if(iT%max<plint>(1,units.getLbSteps(logT)) == 0){
        clock_t end = clock();
        T time = difftime(end,start)/((T)CLOCKS_PER_SEC);
        T mlups = ((T) (lattice.getNx()*lattice.getNy()*lattice.getNz()*units.getLbSteps(logT)))/time/1e6;
        pcout << "time: " << time << " " ;
        pcout << "calculating at " << mlups << " MLU/s" << std::endl;
        start = clock();
        pcerr << x[0][0] << " " << x[0][1] << " " << x[0][2] << " | ";
        pcerr << v[0][0] << " " << v[0][1] << " " << v[0][2] << " | ";
        pcerr << omega[0][0] << " " << omega[0][1] << " " << omega[0][2] << " | ";
        pcerr << f[0][0] << " " << f[0][1] << " " << f[0][2] << " | ";
        pcerr << t[0][0] << " " << t[0][1] << " " << t[0][2] << std::endl;
      }
      
      
      Array<T,3> x_(units.getLbLength(x[0][0]),
                    units.getLbLength(x[0][1]),
                    units.getLbLength(x[0][2]));
      Array<T,3> v_(units.getLbVel(v[0][0]),
                    units.getLbVel(v[0][1]),
                    units.getLbVel(v[0][2]));
      Array<T,3> omega_(units.getLbFreq(omega[0][0]),
                        units.getLbFreq(omega[0][1]),
                        units.getLbFreq(omega[0][2]));
    
      applyProcessingFunctional
        (new SetSphere3D<T,DESCRIPTOR>(x_,v_,omega_,units.getLbLength(r),*id),
         lattice.getBoundingBox(), lattice);

        // Execute a time iteration.
      lattice.collideAndStream();


      v[0][0] += 0.5*dt_dem*f[0][0]/m;
      v[0][1] += 0.5*dt_dem*f[0][1]/m;
      v[0][2] += 0.5*dt_dem*(f[0][2]/m-g*(rho_fact-1));

      omega[0][0] += 0.5*dt_dem*t[0][0]/I;
      omega[0][1] += 0.5*dt_dem*t[0][1]/I;
      omega[0][2] += 0.5*dt_dem*t[0][2]/I;

    }

    delete boundaryCondition;
}
