#ifndef PERIODIC_PRESSURE_FUNCTIONALS_3D
#define PERIODIC_PRESSURE_FUNCTIONALS_3D

namespace plb {
  /* 
     this implements the Zhang/Kwok periodic pressore boundary condition
     which only works for very low Mach numbers
  */
template<typename T, template<typename U> class Descriptor>
class PeriodicPressureFunctional3D : public BoxProcessingFunctional3D_L<T,Descriptor> {
public:
  // dimension: 0,1,2 for x,y,z
  // direction: +1 for pos, -1 for neg
  PeriodicPressureFunctional3D(T const deltaRho, T const rhoAvg, 
                               plint const dimension, plint const direction);
  
  void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
  virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
  PeriodicPressureFunctional3D<T,Descriptor>* clone() const;  
private:
  typedef std::vector<plint> IndexVec;
  IndexVec rescalePop;
  T rescaleFactor;

  void processCell(Cell<T,Descriptor>& cell);

};

// template<typename T, template<typename U> class Descriptor>
// class StripeOffEquilibriumFunctional3D : public BoxProcessingFunctional3D_L<T,Descriptor> {
// public:
//   // dimension: 0,1,2 for x,y,z
//   // direction: +1 for pos, -1 for neg
//   StripeOffEquilibriumFunctional3D(T const deltaRho, T const rhoAvg, 
//                                plint const dimension, plint const direction);
  
//   void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
//   virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
//   PeriodicPressureFunctional3D<T,Descriptor>* clone() const;  
// private:
//   typedef std::vector<plint> IndexVec;
//   IndexVec rescalePop;
//   plint xoff,yoff,zoff;
//   void processCell(Cell<T,Descriptor>& cell);

// };


}; /* namespace plb */

#include "periodicPressureFunctionals3D.hh"

#endif /* PERIODIC_PRESSURE_FUNCTIONALS_3D */
