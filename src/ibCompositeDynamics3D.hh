#ifndef IB_COMPOSITE_DYNAMICS_HH_LBDEM
#define IB_COMPOSITE_DYNAMICS_HH_LBDEM

#include "ibDef.h"

namespace plb {

  template<typename T, template<typename U> class Descriptor>
  void IBdynamicsParticleData<T,Descriptor>::serialize(HierarchicSerializer &serializer) const
  {

    for(plint i=0;i<Descriptor<T>::d;i++)
      serializer.addValue(uPart[i]);

    for(plint i=0;i<Descriptor<T>::d;i++)
      serializer.addValue(hydrodynamicForce[i]);
    
    serializer.addValue(solidFraction);
    serializer.addValue((T)partId);

  }

  template<typename T, template<typename U> class Descriptor>
  void IBdynamicsParticleData<T,Descriptor>::unserialize(HierarchicUnserializer &unserializer)
  {

    for(plint i=0;i<Descriptor<T>::d;i++)
      unserializer.readValue(uPart[i]);

    for(plint i=0;i<Descriptor<T>::d;i++)
      unserializer.readValue(hydrodynamicForce[i]);

    unserializer.readValue(solidFraction);


    // make sure no roundoff errors are introduced in
    // conversion back from float to int
    T partIdDummy;
    unserializer.readValue(partIdDummy);
    partId = (plint)floor(partIdDummy+0.01); 

  }


  template<typename T, template<typename U> class Descriptor>
  int IBcompositeDynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor, IBcompositeDynamics<T,Descriptor> >("IBcomposite");

  template<typename T, template<typename U> class Descriptor>
  IBcompositeDynamics<T,Descriptor>::IBcompositeDynamics(Dynamics<T,Descriptor>* baseDynamics_,
							 bool automaticPrepareCollision_ = true)
    : CompositeDynamics<T,Descriptor>(baseDynamics_,automaticPrepareCollision_), 
      particleData()
  { }
  
  template<typename T, template<typename U> class Descriptor>
  IBcompositeDynamics<T,Descriptor>::IBcompositeDynamics(HierarchicUnserializer &unserializer)
    : CompositeDynamics<T,Descriptor>(0,false)
  {
    // pcout << "entering serialize constructor" << std::endl;
    unserialize(unserializer);
  }

  template<typename T, template<typename U> class Descriptor>
  IBcompositeDynamics<T,Descriptor>::IBcompositeDynamics(const IBcompositeDynamics &orig)
    : CompositeDynamics<T,Descriptor>(orig),
      particleData(orig.particleData),
      fPre(orig.fPre)
  { }
  
  template<typename T, template<typename U> class Descriptor>
  IBcompositeDynamics<T,Descriptor>::~IBcompositeDynamics() {}
  
  template<typename T, template<typename U> class Descriptor>
  IBcompositeDynamics<T,Descriptor>* IBcompositeDynamics<T,Descriptor>::clone() const {
    return new IBcompositeDynamics<T,Descriptor>(*this);
  }
  
  template<typename T, template<typename U> class Descriptor>
  int IBcompositeDynamics<T,Descriptor>::getId() const
  {
    return id;
  }
  
  template<typename T, template<typename U> class Descriptor>
  void IBcompositeDynamics<T,Descriptor>::serialize(HierarchicSerializer &serializer) const
  {

    particleData.serialize(serializer);
    CompositeDynamics<T,Descriptor>::serialize(serializer);
  }

  template<typename T, template<typename U> class Descriptor>
  void IBcompositeDynamics<T,Descriptor>::unserialize(HierarchicUnserializer &unserializer)
  {
    PLB_PRECONDITION( unserializer.getId() == this->getId() );

    particleData.unserialize(unserializer);
    CompositeDynamics<T,Descriptor>::unserialize(unserializer);
  }
  
  template<typename T, template<typename U> class Descriptor>
  void IBcompositeDynamics<T,Descriptor>::prepareCollision(Cell<T,Descriptor>& cell)
  {

    if(particleData.solidFraction > SOLFRAC_MIN)
      fPre = cell.getRawPopulations();

  }

  template<typename T, template<typename U> class Descriptor>
  void IBcompositeDynamics<T,Descriptor>::defineVelocity(Cell<T,Descriptor>& cell, 
                                                         Array<T,Descriptor<T>::d> const& u)
  {
    Array<T,Descriptor<T>::q> fEq;
    T const rhoBar = 1.;
    T const uSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(u); 
    CompositeDynamics<T,Descriptor>::getBaseDynamics().computeEquilibria(fEq,0.,u,uSqr);
    for(plint i=0;i<Descriptor<T>::q;i++)
      cell[i] = fEq[i];
  }
  
  /// Implementation of the collision step
  template<typename T, template<typename U> class Descriptor>
  void IBcompositeDynamics<T,Descriptor>::collide(Cell<T,Descriptor>& cell,
                                                  BlockStatistics& statistics)
  {

    prepareCollision(cell);

    particleData.hydrodynamicForce.resetToZero();

    if(particleData.solidFraction <= SOLFRAC_MAX)
      CompositeDynamics<T,Descriptor>::getBaseDynamics().collide(cell,statistics);
    
    if(particleData.solidFraction < SOLFRAC_MIN)
      return;

    
    // compute equilibrium distributions for solid velocity
    Array<T,Descriptor<T>::q> fEqSolid;

    T const rhoBar = momentTemplates<T,Descriptor>::get_rhoBar(cell);
    Array<T,Descriptor<T>::d> uPart = particleData.uPart*(1.+rhoBar);
    T const uPartSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(uPart); 
    CompositeDynamics<T,Descriptor>::getBaseDynamics().computeEquilibria(fEqSolid,rhoBar,uPart,uPartSqr);

    // T const uPartSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(particleData.uPart); 
    // CompositeDynamics<T,Descriptor>::getBaseDynamics().computeEquilibria(fEq,0.,particleData.uPart,uPartSqr);
    
    if(particleData.solidFraction > SOLFRAC_MAX){
      cell[0] = fPre[0];
      for(plint iPop=1;iPop<Descriptor<T>::q;iPop++){
        plint const shift = Descriptor<T>::q/2;
        plint iOpposite = iPop <= shift ? iPop + shift : iPop - shift;
        
        T coll = 0.5*(fPre[iOpposite] - fEqSolid[iOpposite] + fEqSolid[iPop] - fPre[iPop]);
        
        cell[iPop] = fPre[iPop] + coll;

        for(plint iDim=0;iDim<Descriptor<T>::d;iDim++)
          particleData.hydrodynamicForce[iDim] -= Descriptor<T>::c[iPop][iDim]*coll;
      }
    } else {
      T const omega = CompositeDynamics<T,Descriptor>::getBaseDynamics().getOmega();
      T const ooo = 1./omega-0.5;
      T const B = particleData.solidFraction*ooo/((1.- particleData.solidFraction) + ooo);
      T const oneMinB = 1. - B;

      cell[0] = fPre[0] + oneMinB*(cell[0] - fPre[0]);
      
      for(plint iPop=1;iPop<Descriptor<T>::q;iPop++){
        plint const shift = Descriptor<T>::q/2;
        plint iOpposite = iPop <= shift ? iPop + shift : iPop - shift;

        T coll = -B*0.5*(fPre[iOpposite] - fEqSolid[iOpposite] + fEqSolid[iPop] - fPre[iPop]);

        cell[iPop] = fPre[iPop] + oneMinB*(cell[iPop] - fPre[iPop]) - coll;

        for(plint iDim=0;iDim<Descriptor<T>::d;iDim++)
          particleData.hydrodynamicForce[iDim] += Descriptor<T>::c[iPop][iDim]*coll;
      }
    }
    
  }

}; /* namespace plb */

#endif /* IB_COMPOSITE_DYNAMICS_HH_LBDEM */
