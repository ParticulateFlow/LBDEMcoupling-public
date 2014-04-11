#ifndef PERIODIC_PRESSURE_FUNCTIONALS_3D
#define PERIODIC_PRESSURE_FUNCTIONALS_3D

namespace plb {

  /**
   * this derives directly from BoxProcessingFunctional
   * it creates a coupling between three objects: a lattice,
   * a scalar field and a tensor field.
   */

  template<typename T1, template<typename U> class Descriptor, typename T2, typename T3, int nDim>
  class BoxProcessingFunctional3D_LST : public BoxProcessingFunctional3D {
  public:
    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor> &lattice,
                         ScalarField3D<T2> &scalar, TensorField3D<T3,nDim> &tensor) =0;
    void processGenericBlocks (Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);

  };

  /*
   * implements a periodic pressure gradient across a periodic boundary
   * following 
   * Kim and Pitsch, Phys. Fluids 19, 108101(2007)
   */

  template<typename T1, template<typename U> class Descriptor, typename T2, typename T3, int nDim>
  class KimPeriodicPressureFunctional3D : public BoxProcessingFunctional3D_LST<T1,Descriptor,T2,T3,nDim> {
  public:
    KimPeriodicPressureFunctional3D(T1 rhoTarget_, T1 rhoAvg_, plint dimension_, plint direction_);
    void process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice,
                 ScalarField3D<T2> &rho, TensorField3D<T3,nDim> &u);
    virtual KimPeriodicPressureFunctional3D<T1,Descriptor,T2,T3,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const;

  private:
    typedef std::vector<plint> IndexVec;
    IndexVec rescalePop;
    T1 rhoTarget,rhoAvg;
    plint dimension, direction;
    T1 computeEquilibrium(Cell<T1,Descriptor>& cell,
                          plint iEq,
                          T2 rho,
                          Array<T3,Descriptor<T1>::d> u);  

  };

  /*
   * convenience wrapper class that contains instantiation and application 
   * of KimPeriodicPressureFunctional3D
   * usage:
   * ppm.preColl(lattice);
   * lattice.collideAndStream();
   * ppm.postColl(lattice);
   * dimension: 0 = x, 1 = y, 2 = z
   * in/out direction: +1 if boundary
   */

  template<typename T, template<typename U> class Descriptor>
  class PeriodicPressureManager {
  public:
    PeriodicPressureManager(MultiBlockLattice3D<T,Descriptor> &lattice, 
                            T rhoIn_, T rhoOut_, Box3D inlet_, Box3D outlet_, 
                            plint dimension_, plint inDirection_, plint outDirection_);
    
    void preColl(MultiBlockLattice3D<T,Descriptor> &lattice);
    void postColl(MultiBlockLattice3D<T,Descriptor> &lattice);
  private:
    T rhoIn, rhoOut, rhoAvgIn, rhoAvgOut;
    Box3D inlet,outlet,tmp;
    plint dimension, inDirection, outDirection;
    MultiScalarField3D<T> rho;
    MultiTensorField3D<T,Descriptor<T>::d> u;
  };
  
  /*
   * class PressureGradient
   * helper class to initialize a domain with a constant pressure gradient
   * usage:
   * initializeAtEquilibrium( lattice, lattice.getBoundingBox(), 
   *                          PressureGradient<T>(rhoHi,rhoLo, nGrid, dimension) );
   * nGrid: grid points in direction of gradient
   * dimension: 0 = x, 1 = y, 2 = z
   */

  template<typename T>
  class PressureGradient {
  public:
    PressureGradient(T pHi_, T pLo_, plint n_, plint dimension_);
    void operator() (plint iX, plint iY, plint iZ, T& density, Array<T,3>& velocity) const;
  private:
    T pHi, pLo;
    plint n,dimension;
  };

}; /* namespace plb */

#include "periodicPressureFunctionals3D.hh"

#endif /* PERIODIC_PRESSURE_FUNCTIONALS_3D */