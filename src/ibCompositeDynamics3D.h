#ifndef IB_COMPOSITE_DYNAMICS_H_LBDEM
#define IB_COMPOSITE_DYNAMICS_H_LBDEM

namespace plb {

  template<typename T, template<typename U> class Descriptor>
  struct IBdynamicsParticleData {
  public:
    IBdynamicsParticleData() : partId(0),solidFraction(0.) 
    {
      uPart.resetToZero();
      hydrodynamicForce.resetToZero();
    }
    IBdynamicsParticleData(IBdynamicsParticleData const &orig)
      : partId(orig.partId), solidFraction(orig.solidFraction),
        uPart(orig.uPart), hydrodynamicForce(orig.hydrodynamicForce) {}
    plint partId;
    T solidFraction;
    Array<T,Descriptor<T>::d> uPart;
    Array<T,Descriptor<T>::d> hydrodynamicForce;
  };
  

  template<typename T, template<typename U> class Descriptor>
  class IBcompositeDynamics : public CompositeDynamics<T,Descriptor> {
  public:
    IBcompositeDynamics(Dynamics<T,Descriptor>* baseDynamics_, bool automaticPrepareCollision_);
    IBcompositeDynamics(const IBcompositeDynamics &orig);
    IBcompositeDynamics(HierarchicUnserializer &unserializer);
    ~IBcompositeDynamics();

    virtual IBcompositeDynamics<T,Descriptor>* clone() const;    
    
    virtual int getId() const;

    virtual void serialize(HierarchicSerializer& serializer) const;
    virtual void unserialize(HierarchicUnserializer& unserializer);

    virtual void prepareCollision(Cell<T,Descriptor>& cell);
    virtual void collide(Cell<T,Descriptor>& cell,  BlockStatistics& statistics);

    IBdynamicsParticleData<T,Descriptor> particleData;
    Array<T,Descriptor<T>::q> fPre; // pre-collision populations
    
    virtual void defineVelocity(Cell<T,Descriptor>& cell, 
                                Array<T,Descriptor<T>::d> const& u);

  private:
    static int id;    

  };

}; /* namespace plb */

#include "ibCompositeDynamics3D.hh"

#endif /* IB_COMPOSITE_DYNAMICS_H_LBDEM */
