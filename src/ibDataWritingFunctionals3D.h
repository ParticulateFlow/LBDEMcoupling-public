#ifndef IBDATAWRITINGFUNCTIONALS3D_H_LBDEM
#define IBDATAWRITINGFUNCTIONALS3D_H_LBDEM

namespace plb {

  enum IBscalarQuantity { SolidFraction, ParticleId };
  enum IBvectorQuantity { ParticleVelocity, HydrodynamicForce };

  template<typename T1, template<typename U> class Descriptor, typename T2>
  struct GetScalarQuantityFromDynamicsFunctional : public BoxProcessingFunctional3D_LS<T1,Descriptor,T2> {
  public:

    IBscalarQuantity which;
    GetScalarQuantityFromDynamicsFunctional(IBscalarQuantity const which_)
      : which(which_) {}
    
    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice,
                         ScalarField3D<T2>& data);

    GetScalarQuantityFromDynamicsFunctional<T1,Descriptor,T2>* clone() const;
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    
  };

  template<typename T1, template<typename U> class Descriptor, typename T2, int nDim>
  struct GetVectorQuantityFromDynamicsFunctional : public BoxProcessingFunctional3D_LT<T1,Descriptor,T2,nDim> {
  public:

    IBvectorQuantity which;
    GetVectorQuantityFromDynamicsFunctional(IBvectorQuantity const which_)
      : which(which_) {}
    
    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice,
                         TensorField3D<T2,nDim>& data);

    GetVectorQuantityFromDynamicsFunctional<T1,Descriptor,T2,nDim>* clone() const;
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
  };

}; /* namespace plb */

#include "ibDataWritingFunctionals3D.hh"

#endif /* IBDATAWRITINGFUNCTIONALS3D_H_LBDEM */
