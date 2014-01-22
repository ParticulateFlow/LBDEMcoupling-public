
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
  
  applyProcessingFunctional
    (new GetExternalScalarFunctional3D<T,DESCRIPTOR,T>(DESCRIPTOR<T>::ExternalField::particleIdBeginsAt),
     lattice.getBoundingBox(), lattice, tmp);
  vtkOut.writeData<float>(tmp, "PartId", 1.);

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
        pcout << "2 : uMax\n";
        exit(1);
    }

    const plint N = atoi(argv[1]);
    const T uMax = atof(argv[2]);

    const T rho_f = 1000;

    char **argv_lmp = 0;
    argv_lmp = new char*[1];
    argv_lmp[0] = argv[0];

    LAMMPS_NS::LAMMPS *lmp = new LAMMPS_NS::LAMMPS(1,argv_lmp,global::mpi().getGlobalCommunicator());

    lmp->input->file("in.lbdem");

    void *fix = locate_coupling_fix(lmp);

    int nAtoms = lammps_get_natoms(lmp);  

    double **x=0, **v=0, **f=0, **r=0;
    int **id=0;
    allocate_external_double(x,3,"nparticles",1.,lmp);
    allocate_external_double(v,3,"nparticles",0.,lmp);
    allocate_external_double(f,3,"nparticles",0.,lmp);
    allocate_external_double(r,1,"nparticles",0.,lmp);
    allocate_external_int(id,1,"nparticles",0,lmp);

    T g = 9.81;
    const T nu_f = 1e-4;


    const T lx = 0.4, ly = 0.2, lz = 0.2;



    data_liggghts_to_of("x","vector-atom",lmp,(void*&)x,"double");
    data_liggghts_to_of("v","vector-atom",lmp,(void*&)v,"double");
    data_liggghts_to_of("radius","scalar-atom",lmp,(void*&)r,"double");
    data_liggghts_to_of("id","scalar-atom",lmp,(void*&)id,"int");
     

    T r_ = r[0][0];
    T rho_s = 1100.;
    T m = r_*r_*r_*4./3.*3.14*rho_s;

    T v_inf = 2.*(rho_s-rho_f)/rho_f*g*r_*r_/9./nu_f; // stokes

    pcout << "v_inf: " << v_inf << " m: " << m << " r: " << r_ << std::endl;


    PhysUnits3D<T> units(2.*r_,v_inf,nu_f,lx,ly,lz,N,uMax,rho_f);

    IncomprFlowParam<T> parameters(units.getLbParam());

    
    
    const T maxT = (T)0.6;
    const T vtkT = 0.001;
    const T gifT = 100;
    const T logT = 0.05;

    const plint maxSteps = units.getLbSteps(maxT);
    const plint vtkSteps = max<plint>(units.getLbSteps(vtkT),1);
    const plint gifSteps = units.getLbSteps(gifT);
    const plint logSteps = max<plint>(units.getLbSteps(logT),1);

    writeLogFile(parameters, "3D sedimenting sphere");

    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
        parameters.getNx(), parameters.getNy(), parameters.getNz(), 
        new DYNAMICS );

    lattice.periodicity().toggleAll(true);

    lattice.initialize();
    lattice.collideAndStream();
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
    
    T **x_lb, **v_lb, **omega_lb;
    x_lb = new T*[nAtoms]; 
    v_lb = new T*[nAtoms]; 
    f = new T*[nAtoms];

    for(plint n=0;n<nAtoms;n++){
      x_lb[n] = new T[3];
      v_lb[n] = new T[3];
      f[n] = new T[3];
    }

    T *r_lb = new T[nAtoms];
    for(plint n=0;n<nAtoms;n++)
      r_lb[n] = units.getLbLength(r[n][0]);

    T dt_dem = dt_phys/10.;
    std::stringstream cmd;
    cmd << "variable t_step equal " << dt_dem;
    pcout << cmd.str() << std::endl;
    lmp->input->one(cmd.str().c_str()); cmd.str("");
    // cmd << "variable dmp_time equal " << vtkT;
    // pcout << cmd.str() << std::endl;
    // lmp->input->one(cmd.str().c_str()); cmd.str("");

    cmd << "variable dmp_stp equal " << vtkSteps*10;
    pcout << cmd.str() << std::endl;
    lmp->input->one(cmd.str().c_str()); cmd.str("");
    
    
    lmp->input->file("in2.lbdem");    
    
    lmp->input->one("run 9 upto");

    clock_t start = clock();    


    // Loop over main time iteration.
    for (plint iT=0; iT<maxSteps; ++iT) {

      nAtoms = lammps_get_natoms(lmp);  
      data_liggghts_to_of("x","vector-atom",lmp,(void*&)x,"double");
      data_liggghts_to_of("v","vector-atom",lmp,(void*&)v,"double");
      data_liggghts_to_of("radius","scalar-atom",lmp,(void*&)r,"double");
      data_liggghts_to_of("id","scalar-atom",lmp,(void*&)id,"int");

      for(plint n=0;n<nAtoms;n++){
        r_lb[n] = units.getLbLength(r[n][0]);
        for(plint i=0;i<3;i++){
          x_lb[n][i] = units.getLbLength(x[n][i]);
          v_lb[n][i] = units.getLbVel(v[n][i]);
        }
      }

      double **omega = 0;
      setSpheresOnLattice(lattice,x_lb,v_lb,omega,r_lb,id,nAtoms);
      // applyProcessingFunctional
      //   (new SetSpheres3D<T,DESCRIPTOR>(x_lb,v_lb,r_lb,id,nAtoms),
      //    lattice.getBoundingBox(), lattice);
    
      if(iT%vtkSteps == 0)
        writeVTK(lattice,parameters,units,iT);

      if(iT%gifSteps == 0)
        writeGif(lattice,iT);      

      lattice.collideAndStream();
      
      for(int n=0;n<nAtoms;n++)
        for(int i=0;i<3;i++){
          f[n][i] = 0;
        }
      SumForceTorque3D<T,DESCRIPTOR> sft(nAtoms,x_lb);
      //if(iT > 1) 
      applyProcessingFunctional(sft,lattice.getBoundingBox(), lattice);
      
      // for(plint i=0;i<sft.getForceTorque().size();i++)
      //   pcerr << sft.getForceTorque()[i] << " ";

      for(plint n=0;n<nAtoms;n++)
        for(int i=0;i<3;i++){
          f[n][i] = units.getPhysForce(sft.getForceTorque()[6*n+i]);
        }

      //      data_of_to_liggghts("dragforce","vector-atom",lmp,(void*)f,"double");
      lmp->input->one("run 10");

      for(plint n=0;n<nAtoms;n++)
        for(int i=0;i<3;i++){
          f[n][i] = units.getPhysForce(sft.getForceTorque()[6*n+i]);
        }
      if(iT%logSteps == 0){
        clock_t end = clock();
        T time = difftime(end,start)/((T)CLOCKS_PER_SEC);
        T mlups = ((T) (lattice.getNx()*lattice.getNy()*lattice.getNz()*units.getLbSteps(logT)))/time/1e6;
        pcout << "time: " << time << " " ;
        pcout << "calculating at " << mlups << " MLU/s" << std::endl;
        start = clock();
      }
      pcerr << x[0][0] << " " << x[0][1] << " " << x[0][2] << " "
            << v[0][0] << " " << v[0][1] << " " << v[0][2] << " "
            << f[0][0] << " " << f[0][1] << " " << f[0][2] << " "
            << x[1][0] << " " << x[1][1] << " " << x[1][2] << std::endl;
    }

    delete boundaryCondition;
}
