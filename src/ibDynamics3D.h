/*
 * This file is part of the LBDEMcoupling software.
 *
 * LBDEMcoupling is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright 2014 Johannes Kepler University Linz
 *
 * Author: Philippe Seil (philippe.seil@jku.at)
 */

#ifndef IB_DYNAMICS_3D_H
#define IB_DYNAMICS_3D_H
#include "ibDef.h"

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

  // template<typename T, template<typename U> class Descriptor>
  // int IBcompositeDynamics<T,Descriptor>::id =
  //   meta::registerCompositeDynamics<T,Descriptor, IBcompositeDynamics<T,Descriptor> >("IBcomposite");
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

    for(plint i=0;i<Descriptor<T>::d;i++)
      serializer.addValue(particleData.uPart[i]);

    for(plint i=0;i<Descriptor<T>::d;i++)
      serializer.addValue(particleData.hydrodynamicForce[i]);
    
    serializer.addValue(particleData.solidFraction);

    T *partIdPtr = (T*) (&particleData.partId);
    serializer.addValue((T)particleData.partId);

    CompositeDynamics<T,Descriptor>::serialize(serializer);
  }

  template<typename T, template<typename U> class Descriptor>
  void IBcompositeDynamics<T,Descriptor>::unserialize(HierarchicUnserializer &unserializer)
  {
    PLB_PRECONDITION( unserializer.getId() == this->getId() );

    for(plint i=0;i<Descriptor<T>::d;i++)
      unserializer.readValue(particleData.uPart[i]);

    for(plint i=0;i<Descriptor<T>::d;i++)
      unserializer.readValue(particleData.hydrodynamicForce[i]);

    unserializer.readValue(particleData.solidFraction);

    T partIdDummy;
    unserializer.readValue(partIdDummy);
    plint *partIdPtr = (plint*) (&partIdDummy);
    particleData.partId = *partIdPtr;

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


  /* ************************************************************************* */
  /* IB dynamics standalone stuff                                              */
  /* ************************************************************************* */

  template<typename T, template<typename U> class Descriptor>
  class IBdynamics : public BGKdynamics<T,Descriptor> {
  public:
    /* *************** Construction / Destruction ************************ */
    IBdynamics(T omega_);
    
    /// Clone the object on its dynamic type.
    virtual IBdynamics<T,Descriptor>* clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

  /* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics& statistics_);
    // no need for equilibrium distribution, it gets inherited from BGKdynamics

  private:
    static int id;
  }; // class IBdynamics

  // implementation of immersed boundary collision
  template<typename T, template<typename U> class Descriptor>
  struct dynamicsTemplatesImplIB {
    static T ib_collision(Array<T,Descriptor<T>::q>& f, T rhoBar, T* externalScalars,
        Array<T,Descriptor<T>::d> const& j, T omega);
  };
  template<typename T, template<typename U> class Descriptor>
  struct dynamicsTemplatesIB {
    static T ib_collision(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, T omega)
    {
        return dynamicsTemplatesImplIB<T,Descriptor>
            ::ib_collision(cell.getRawPopulations(), rhoBar, cell.getExternal(0), j, omega);
    }
  };

  /* *************** Class IBdynamics *********************************************** */

  template<typename T, template<typename U> class Descriptor>
  int IBdynamics<T,Descriptor>::id =
      meta::registerOneParamDynamics<T,Descriptor,IBdynamics<T,Descriptor> >("IB");


  /** \param omega_ relaxation parameter, related to the dynamic viscosity
   */
  template<typename T, template<typename U> class Descriptor>
  IBdynamics<T,Descriptor>::IBdynamics(T omega_ )
      : BGKdynamics<T,Descriptor>(omega_)
  { }

  template<typename T, template<typename U> class Descriptor>
  IBdynamics<T,Descriptor>* IBdynamics<T,Descriptor>::clone() const {
    return new IBdynamics<T,Descriptor>(*this);
  }

  template<typename T, template<typename U> class Descriptor>
  int IBdynamics<T,Descriptor>::getId() const {
    return id;
  }

  template<typename T, template<typename U> class Descriptor>
  void IBdynamics<T,Descriptor>::collide (
          Cell<T,Descriptor>& cell,
          BlockStatistics& statistics )
  {
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T uSqr = dynamicsTemplatesIB<T,Descriptor>::ib_collision(cell, rhoBar, j, this->getOmega());
    if (cell.takesStatistics()) {
      gatherStatistics(statistics, rhoBar, uSqr);
    }

  }


  template<typename T, template<typename U> class Descriptor>
  T dynamicsTemplatesImplIB<T,Descriptor>::ib_collision(Array<T,Descriptor<T>::q>& f, T rhoBar,
    T* externalScalars, Array<T,Descriptor<T>::d> const& j, T omega)
  {

    T const solidFraction = externalScalars[Descriptor<T>::ExternalField::volumeFractionBeginsAt];

    T* force = &externalScalars[Descriptor<T>::ExternalField::hydrodynamicForceBeginsAt];
    for(int i=0;i<Descriptor<T>::d;i++) force[i] = 0.;


    if(solidFraction < SOLFRAC_MIN) // pure liquid --> bgk collision
      return dynamicsTemplatesImpl<T,Descriptor<T> >::bgk_inc_collision(f, rhoBar, j, omega);

    T const invRho = Descriptor<T>::invRho(rhoBar);
    T const jSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(j);

    Array<T,Descriptor<T>::d> uPart;
    uPart.from_cArray(&externalScalars[Descriptor<T>::ExternalField::boundaryVelocityBeginsAt]);

    for(int i=0;i<Descriptor<T>::d;i++) uPart[i] *= (1.+rhoBar);

    T const uPartSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(uPart); 
    T const ooo = 1./omega-0.5;
    T const B = solidFraction*ooo/((1.- solidFraction) + ooo);
    

    if(B > SOLFRAC_MAX){ // then we have pure solid and do not need any collision
      for(int iPop=1;iPop <= Descriptor<T>::q/2;iPop++){
        plint iOpposite = iPop+Descriptor<T>::q/2;
        T const bias = 
          dynamicsTemplatesImpl<T,Descriptor<T> >::bgk_ma2_equilibrium 
	  (iPop, rhoBar, invRho, uPart, uPartSqr ) -
	  dynamicsTemplatesImpl<T,Descriptor<T> >::bgk_ma2_equilibrium 
	  (iOpposite, rhoBar, invRho, uPart, uPartSqr );
	
        T const bounce = 0.5*(f[iOpposite] - f[iPop] + bias);

        
        f[iPop] += bounce;
        f[iOpposite] -= bounce;
        
        for(int i=0;i<Descriptor<T>::d;i++) 
          force[i] -= 2.*Descriptor<T>::c[iPop][i]*bounce; 

      }
    } else{
      T const omega_1minB = omega*(1.-B);

      f[0] -= 
	omega_1minB*( f[0] - 
		      dynamicsTemplatesImpl<T,Descriptor<T> >::bgk_ma2_equilibrium (0, rhoBar, invRho, j, jSqr ) );
      
      for(int iPop=1;iPop <= Descriptor<T>::q/2;iPop++){
        plint iOpposite = iPop+Descriptor<T>::q/2;
        
        T const bias =
	  dynamicsTemplatesImpl<T,Descriptor<T> >::bgk_ma2_equilibrium
          (iPop, rhoBar, invRho, uPart, uPartSqr ) -
	  dynamicsTemplatesImpl<T,Descriptor<T> >::bgk_ma2_equilibrium
          (iOpposite, rhoBar, invRho, uPart, uPartSqr );
        
        T const bounce = 0.5*B*(f[iOpposite] - f[iPop] + bias);
        
        T const deltaFColl = f[iPop] 
	  - dynamicsTemplatesImpl<T,Descriptor<T> >::bgk_ma2_equilibrium (iPop, rhoBar, invRho, j, jSqr );
        T const deltaFCollOp = f[iOpposite] 
	  - dynamicsTemplatesImpl<T,Descriptor<T> >::bgk_ma2_equilibrium (iOpposite, rhoBar, invRho, j, jSqr );
        
        
        f[iPop] -= omega_1minB*deltaFColl;           f[iPop] += bounce;
        f[iOpposite] -= omega_1minB*deltaFCollOp;    f[iOpposite] -= bounce;

        for(int i=0;i<Descriptor<T>::d;i++) 
          force[i] -= 2.*Descriptor<T>::c[iPop][i]*bounce;
      }
        
    }
    return jSqr*invRho*invRho;

  }

};

#endif // IB_DYNAMICS_3D_H
