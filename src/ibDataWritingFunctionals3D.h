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
                         ScalarField3D<T2>& data) 
    {
      for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
          for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
            
            Cell<T1,Descriptor>& cell = lattice.get(iX,iY,iZ);
            Dynamics<T1,Descriptor> *dyn = &(cell.getDynamics());
            T1 val = 0;
            if(dyn->isComposite()){
              IBcompositeDynamics<T1,Descriptor> *cDyn = 
                static_cast< IBcompositeDynamics<T1,Descriptor>* >( dyn );
              switch(which){
              case SolidFraction:
                val = (T2) cDyn->particleData.solidFraction;
                break;
              case ParticleId:
                val = (T2) cDyn->particleData.partId;
                break;

              }

            }
            data.get(iX,iY,iZ) = val;
          }
        }
      }
    }

    GetScalarQuantityFromDynamicsFunctional<T1,Descriptor,T2>* clone() const
    {
      return new GetScalarQuantityFromDynamicsFunctional<T1,Descriptor,T2>(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {
      modified[0] = modif::nothing;
      modified[1] = modif::staticVariables;
    }

    
  };

  template<typename T1, template<typename U> class Descriptor, typename T2, int nDim>
  struct GetVectorQuantityFromDynamicsFunctional : public BoxProcessingFunctional3D_LT<T1,Descriptor,T2,nDim> {
  public:

    IBvectorQuantity which;
    GetVectorQuantityFromDynamicsFunctional(IBvectorQuantity const which_)
      : which(which_) {}
    
    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice,
                         TensorField3D<T2,nDim>& solfrac) 
    {
      for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
          for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
            
            Cell<T1,Descriptor>& cell = lattice.get(iX,iY,iZ);
            Dynamics<T1,Descriptor> *dyn = &(cell.getDynamics());
            Array<T2,nDim> val;
            val.resetToZero();
            if(dyn->isComposite()){
              IBcompositeDynamics<T1,Descriptor> *cDyn = 
                static_cast< IBcompositeDynamics<T1,Descriptor>* >( dyn );
              switch(which){
              case ParticleVelocity:
                val = cDyn->particleData.uPart;
                break;
              case HydrodynamicForce:
                val = cDyn->particleData.hydrodynamicForce;
                break;

              }

            }
            solfrac.get(iX,iY,iZ) = val;
          }
        }
      }
    }

    GetVectorQuantityFromDynamicsFunctional<T1,Descriptor,T2,nDim>* clone() const
    {
      return new GetVectorQuantityFromDynamicsFunctional<T1,Descriptor,T2,nDim>(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {
      modified[0] = modif::nothing;
      modified[1] = modif::staticVariables;
    }
  };

}; /* namespace plb */

#endif /* IBDATAWRITINGFUNCTIONALS3D_H_LBDEM */
