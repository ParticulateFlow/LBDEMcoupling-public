/*
 * ibDynamics3D.h
 */  

#ifndef IB_DYNAMICS_3D_H
#define IB_DYNAMICS_3D_H

namespace plb {

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
      meta::registerOneParamDynamics<T,Descriptor,BGKdynamics<T,Descriptor> >("IB");


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
    T invRho = Descriptor<T>::invRho(rhoBar),
      jSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(j);

    T solidFraction = externalScalars[Descriptor<T>::ExternalField::volumeFractionBeginsAt];

    T* force = &externalScalars[Descriptor<T>::ExternalField::hydrodynamicForceBeginsAt];
    for(int i=0;i<Descriptor<T>::d;i++) force[i] = 0.;


    if(solidFraction < 0.001) // pure liquid --> bgk collision
      return dynamicsTemplatesImpl<T,Descriptor<T> >::bgk_inc_collision(f, rhoBar, j, omega);

    Array<T,Descriptor<T>::d> uPart;
    uPart.from_cArray(&externalScalars[Descriptor<T>::ExternalField::boundaryVelocityBeginsAt]);

    for(int i=0;i<Descriptor<T>::d;i++) uPart[i] /= invRho;

    T uPartSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(uPart); 
    T B = solidFraction;//*(1./omega-0.5)/(0.5 - solidFraction + 1./omega);

    if(B > 0.999){ // then we have pure solid and do not need any collision
      for(int iPop=1,iOpposite=Descriptor<T>::q/2+1;iPop <= Descriptor<T>::q/2;iPop++,iOpposite++){
        T bias = 
          dynamicsTemplatesImpl<T,Descriptor<T> >::bgk_ma2_equilibrium 
	  (iPop, rhoBar, invRho, uPart, uPartSqr ) -
	  dynamicsTemplatesImpl<T,Descriptor<T> >::bgk_ma2_equilibrium 
	  (iOpposite, rhoBar, invRho, uPart, uPartSqr );
	
        T bounce = 0.5*(f[iOpposite] - f[iPop] + bias);
        
        f[iPop] += bounce;
        f[iOpposite] -= bounce;
        
        for(int i=0;i<Descriptor<T>::d;i++) 
          force[i] -= 2.*Descriptor<T>::c[iPop][i]*bounce;	
      }
    } else{
      T omega_1minB = omega*(1.-B);
      
      f[0] -= 
	omega_1minB*( f[0] - 
		      dynamicsTemplatesImpl<T,Descriptor<T> >::bgk_ma2_equilibrium (0, rhoBar, invRho, j, jSqr ) );
      
      for(int iPop=1,iOpposite=Descriptor<T>::q/2+1;iPop <= Descriptor<T>::q/2;iPop++,iOpposite++){
        
        T bias =
	  dynamicsTemplatesImpl<T,Descriptor<T> >::bgk_ma2_equilibrium
          (iPop, rhoBar, invRho, uPart, uPartSqr ) -
	  dynamicsTemplatesImpl<T,Descriptor<T> >::bgk_ma2_equilibrium
          (iOpposite, rhoBar, invRho, uPart, uPartSqr );
        
        T bounce = 0.5*B*(f[iOpposite] - f[iPop] + bias);
        
        T deltaFColl = f[iPop] 
	  - dynamicsTemplatesImpl<T,Descriptor<T> >::bgk_ma2_equilibrium (iPop, rhoBar, invRho, j, jSqr );
        T deltaFCollOp = f[iOpposite] 
	  - dynamicsTemplatesImpl<T,Descriptor<T> >::bgk_ma2_equilibrium (iOpposite, rhoBar, invRho, j, jSqr );
        
        
        f[iPop] -= omega_1minB*deltaFColl;           f[iPop] += bounce;
        f[iOpposite] -= omega_1minB*deltaFCollOp;    f[iOpposite] -= bounce;

        for(int i=0;i<Descriptor<T>::d;i++) 
          force[i] -= 2.*Descriptor<T>::c[iPop][i]*bounce;
        //force[i] -= Descriptor<T>::c[iPop][i]*bounce;
      }
    }
   
    return jSqr*invRho*invRho;

  }

};

#endif // IB_DYNAMICS_3D_H
