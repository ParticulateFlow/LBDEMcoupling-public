
#include "palabos3D.h"
#include "palabos3D.hh"


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
  vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", units.getPhysVel(1));
  vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", units.getPhysVel(1));  
  vtkOut.writeData<float>(*computeDensity(lattice), "density",units.getPhysDensity(1)); 
  
  MultiScalarField3D<T> p(*computeDensity(lattice));
  subtractInPlace(p,1.);
  vtkOut.writeData<float>(p,"pressure",p_fact ); 

  pcout << "wrote " << fname << std::endl;
}

template<typename T, template<typename U> class Descriptor>
class PeriodicPressureFunctional3D : public BoxProcessingFunctional3D_L<T,Descriptor> {
public:
  // dimension: 0,1,2 for x,y,z
  // direction: +1 for pos, -1 for neg
  PeriodicPressureFunctional3D(T const deltaRho, T const rhoAvg, 
                               plint const dimension, plint const direction)
    : rescaleFactor((1.+3.*deltaRho)/rhoAvg)
  {
    // to do: perform error checking here
    for(plint iPop=0;iPop<Descriptor<T>::q;iPop++){
      if(Descriptor<T>::c[iPop][dimension] == direction){
        rescalePop.push_back(iPop);      
      }
    }
  }
  
  void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice)
  {
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
          Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
          for(IndexVec::iterator it = rescalePop.begin();
              it != rescalePop.end(); it++){
            T fTmp = cell[*it] + Descriptor<T>::t[*it];
            fTmp *= rescaleFactor;
            cell[*it] = fTmp - Descriptor<T>::t[*it];
          }
        }
      }
    }
  }
  virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
  {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::dynamicVariables;
    modified[3] = modif::allVariables;
  }
  PeriodicPressureFunctional3D<T,Descriptor>* clone() const 
  { return new PeriodicPressureFunctional3D<T,Descriptor>(*this);}
  
private:
  typedef std::vector<plint> IndexVec;
  IndexVec rescalePop;
  T rescaleFactor;
};

int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    plint N;
    T uMax = 0.02;
    T deltaP;
    T lFactor;
    try {
        global::argv(1).read(N);
        global::argv(2).read(deltaP);
        global::argv(3).read(lFactor);
    } catch(PlbIOException& exception) {
        pcout << exception.what() << endl;
        pcout << "Command line arguments:\n";
        pcout << "1 : N\n";
        pcout << "2 : deltaP\n";
        pcout << "3 : lengthFactor\n";
        exit(1);
    }

    const T rho_f = 1000;

    const T nu_f = 1e-3;


    const T lx = 0.2, ly = 0.2, lz = 0.2*lFactor;
    
    T gradP = deltaP/lz/rho_f;
    T u_phys = gradP*lx*lx/(nu_f*8);

    
    PhysUnits3D<T> units(lx,u_phys,nu_f,lx,ly,lz,N,uMax,rho_f);

    IncomprFlowParam<T> parameters(units.getLbParam());

    const T maxT = (T)50.;
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

    Box3D inlet(1,nx-2,0,ny-1,0,0), outlet(1,nx-2,0,ny-1,nz-1,nz-1);
    Box3D backWall(0,0,0,ny-1,0,nz-1), frontWall(nx-1,nx-1,0,ny-1,0,nz-1);
    defineDynamics(lattice,backWall,new BounceBack<T,DESCRIPTOR>);
    defineDynamics(lattice,frontWall,new BounceBack<T,DESCRIPTOR>);
    //Box3D inlet(1,nx-2,1,ny-2,0,0), outlet(1,nx-2,1,ny-2,nz-1,nz-1);
    //Box3D inlet(1,nx-2,0,ny-1,0,0), outlet(1,nx-2,0,ny-1,nz-1,nz-1);

    
    T deltaRho = units.getLbRho(deltaP);//*((T)nz)/((T)(nz-2));//((T)(nx-1))/((T)nx);


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

      lattice.collideAndStream();

      T rhoAvgIn = computeAverageDensity(lattice,inlet);
      T rhoAvgOut = computeAverageDensity(lattice,outlet);

      applyProcessingFunctional(new PeriodicPressureFunctional3D<T,DESCRIPTOR>(0., rhoAvgOut,2, 1),
                                inlet,lattice);
      applyProcessingFunctional(new PeriodicPressureFunctional3D<T,DESCRIPTOR>(-deltaRho, rhoAvgIn,2, -1),
                                outlet,lattice);



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
