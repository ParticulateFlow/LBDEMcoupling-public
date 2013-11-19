
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
    //Box3D slice(0, nx-1, (ny-1)/2, (ny-1)/2, 0, nz-1);
    Box3D slice(0, nx-1, 0, ny-1, (nz-1)/2, (nz-1)/2);
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

    if (argc != 2) {
        pcout << "Error the parameters are wrong. The structure must be :\n";
        pcout << "1 : N\n";
        exit(1);
    }

    const plint N = atoi(argv[1]);
    const T Re = 10.0;
    const T uMax = 0.01;

    IncomprFlowParam<T> parameters(
            uMax,
            Re,
            N,
            1.,        // lx
            1.,        // ly
            1.         // lz
    );

    const plint nx = parameters.getNx();
    const plint ny = parameters.getNx();
    const plint nz = parameters.getNx();

    const T maxT     = (T)10.0;
    const T vtkT = 0.05;
    const T gifT = 0.05;
    const T logT = 0.01;

    writeLogFile(parameters, "3D sedimenting sphere");

    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
        parameters.getNx(), parameters.getNy(), parameters.getNz(), 
        new DYNAMICS );

    // Use periodic boundary conditions.
    lattice.periodicity().toggle(2,true);

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    plint *id = new plint;
    *id = 0;
    T **x, **v, **f, **t, **omega;
    T r = nx/8;
    T rho = 1.01;
    T m = 4./3.*r*r*r*3.14*rho;
    T I = 2./5.*m*r*r;
    x = new T*[1]; x[0] = new T[3];
    v = new T*[1]; v[0] = new T[3];
    omega = new T*[1]; omega[0] = new T[3];

    f = new T*[1]; f[0] = new T[3];
    t = new T*[1]; t[0] = new T[3];

    x[0][0] = nx/2; x[0][1] = ny/2; x[0][2] = nz/2;
    v[0][0] = 0; v[0][1] = 0; v[0][2] = 0;
    omega[0][0] = 0; omega[0][1] = 0; omega[0][2] = 0;
    
    f[0][0] = 0; f[0][1] = 0; f[0][2] = 0;
    t[0][0] = 0; t[0][1] = 0; t[0][2] = 0;

    T g = 100;

    T dt_dem = parameters.getDeltaT();


    Array<T,3> x_(x[0][0],x[0][1],x[0][2]),v_(v[0][0],v[0][1],v[0][2]),
      omega_(omega[0][0],omega[0][1],omega[0][2]);
    
    applyProcessingFunctional
      (new SetSphere3D<T,DESCRIPTOR>(x_,v_,omega_,r,*id),
       lattice.getBoundingBox(), lattice);
    

    clock_t start = clock();


    // Loop over main time iteration.
    for (plint iT=0; iT<parameters.nStep(maxT); ++iT) {

      v[0][0] += 0.5*dt_dem*f[0][0]/m;
      v[0][1] += 0.5*dt_dem*f[0][1]/m;
      v[0][2] += 0.5*dt_dem*(f[0][2]-g)/m;

      x[0][0] += v[0][0]*dt_dem + 0.5*f[0][0]*dt_dem*dt_dem/m;
      x[0][1] += v[0][1]*dt_dem + 0.5*f[0][1]*dt_dem*dt_dem/m;
      x[0][2] += v[0][2]*dt_dem + 0.5*(f[0][2]-g)*dt_dem*dt_dem/m;

      omega[0][0] += 0.5*dt_dem*t[0][0]/I;
      omega[0][1] += 0.5*dt_dem*t[0][1]/I;
      omega[0][2] += 0.5*dt_dem*t[0][2]/I;

      if(iT%parameters.nStep(vtkT) == 0)
        writeVTK(lattice,parameters,iT);

      if(iT%parameters.nStep(gifT) == 0)
        writeGif(lattice,iT);
      
      if(iT%parameters.nStep(logT) == 0){
        clock_t end = clock();
        T time = difftime(end,start)/((T)CLOCKS_PER_SEC);
        T mlups = ((T) (lattice.getNx()*lattice.getNy()*lattice.getNz()*parameters.nStep(logT)))/time/1e6;
        pcout << "time: " << time << " " ;
        pcout << "calculating at " << mlups << " MLU/s" << std::endl;
        start = clock();
        pcout << "force: " << f[0][0] << " " << f[0][1] << " " << f[0][2] << " | ";
        pcout << "torque: " << t[0][0] << " " << t[0][1] << " " << t[0][2] << std::endl;
        pcout << x[0][0] << " " << x[0][1] << " " << x[0][2] << " | ";
        pcout << v[0][0]
              << " " << v[0][1] << " " << v[0][2] << " " << std::endl;
      }
      
      f[0][0] = 0;      f[0][1] = 0;      f[0][2] = 0;
      t[0][0] = 0;      t[0][1] = 0;      t[0][2] = 0;
      
      if(iT > 1) applyProcessingFunctional(new SumForceTorque3D<T,DESCRIPTOR>(id,x,f,t),
                                           lattice.getBoundingBox(), lattice);

      
      Array<T,3> x_(x[0][0],x[0][1],x[0][2]),v_(v[0][0],v[0][1],v[0][2]),
        omega_(omega[0][0],omega[0][1],omega[0][2]);
    
      applyProcessingFunctional
        (new SetSphere3D<T,DESCRIPTOR>(x_,v_,omega_,r,*id),
         lattice.getBoundingBox(), lattice);

        // Execute a time iteration.
      lattice.collideAndStream();


      v[0][0] += 0.5*dt_dem*f[0][0]/m;
      v[0][1] += 0.5*dt_dem*f[0][1]/m;
      v[0][2] += 0.5*dt_dem*(f[0][2]-g)/m;

      omega[0][0] += 0.5*dt_dem*t[0][0]/I;
      omega[0][1] += 0.5*dt_dem*t[0][1]/I;
      omega[0][2] += 0.5*dt_dem*t[0][2]/I;

    }

    delete boundaryCondition;
}
