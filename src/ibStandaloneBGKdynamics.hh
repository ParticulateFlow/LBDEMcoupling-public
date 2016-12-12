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

#ifndef IB_STANDALONE_BGK_DYNAMICS_HH_LBDEM
#define IB_STANDALONE_BGK_DYNAMICS_HH_LBDEM

#include "ibDef.h"

namespace plb {

  // leave this the same because some other parts rely on this information
  template<typename T, template<typename U> class Descriptor>
  int IBstandaloneBGKdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor, IBstandaloneBGKdynamics<T,Descriptor> >("IBstandaloneBGK");

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q> IBstandaloneBGKdynamics<T,Descriptor>::fEq =
    Array<T,Descriptor<T>::q>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q> IBstandaloneBGKdynamics<T,Descriptor>::fEqSolid =
    Array<T,Descriptor<T>::q>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q> IBstandaloneBGKdynamics<T,Descriptor>::fPostColl =
    Array<T,Descriptor<T>::q>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q/2+1> IBstandaloneBGKdynamics<T,Descriptor>::f_plus
    = Array<T,Descriptor<T>::q/2+1>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q/2+1> IBstandaloneBGKdynamics<T,Descriptor>::f_minus
    = Array<T,Descriptor<T>::q/2+1>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q/2+1> IBstandaloneBGKdynamics<T,Descriptor>::fEq_plus
    = Array<T,Descriptor<T>::q/2+1>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q/2+1> IBstandaloneBGKdynamics<T,Descriptor>::fEq_minus
    = Array<T,Descriptor<T>::q/2+1>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q/2+1> IBstandaloneBGKdynamics<T,Descriptor>::fEqSolid_plus
    = Array<T,Descriptor<T>::q/2+1>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q/2+1> IBstandaloneBGKdynamics<T,Descriptor>::fEqSolid_minus
    = Array<T,Descriptor<T>::q/2+1>();

  template<typename T, template<typename U> class Descriptor>
  IBstandaloneBGKdynamics<T,Descriptor>::IBstandaloneBGKdynamics(T const omega_, T const lambda_)
    : IsoThermalBulkDynamics<T,Descriptor>(omega_), 
    lambda(lambda_)
  { }
  
  template<typename T, template<typename U> class Descriptor>
  IBstandaloneBGKdynamics<T,Descriptor>::IBstandaloneBGKdynamics(HierarchicUnserializer &unserializer)
    : IsoThermalBulkDynamics<T,Descriptor>(1.)
  {
    unserialize(unserializer);
  }

  template<typename T, template<typename U> class Descriptor>
  IBstandaloneBGKdynamics<T,Descriptor>::IBstandaloneBGKdynamics(const IBstandaloneBGKdynamics &orig)
    : IsoThermalBulkDynamics<T,Descriptor>(orig),
    IBdynamicsParticleData<T,Descriptor>(orig),
    lambda(orig.lambda)
  { }
  
  template<typename T, template<typename U> class Descriptor>
  IBstandaloneBGKdynamics<T,Descriptor>::~IBstandaloneBGKdynamics() {}
  
  template<typename T, template<typename U> class Descriptor>
  IBstandaloneBGKdynamics<T,Descriptor>* IBstandaloneBGKdynamics<T,Descriptor>::clone() const {
    return new IBstandaloneBGKdynamics<T,Descriptor>(*this);
  }
  
  template<typename T, template<typename U> class Descriptor>
  int IBstandaloneBGKdynamics<T,Descriptor>::getId() const
  {
    return id;
  }
  
  template<typename T, template<typename U> class Descriptor>
  void IBstandaloneBGKdynamics<T,Descriptor>::serialize(HierarchicSerializer &serializer) const
  {

    this->particleData.serialize(serializer);
    serializer.addValue(lambda);
    IsoThermalBulkDynamics<T,Descriptor>::serialize(serializer);
  }

  template<typename T, template<typename U> class Descriptor>
  void IBstandaloneBGKdynamics<T,Descriptor>::unserialize(HierarchicUnserializer &unserializer)
  {
    PLB_PRECONDITION( unserializer.getId() == this->getId() );

    this->particleData.unserialize(unserializer);
    unserializer.readValue<T>(lambda);
    IsoThermalBulkDynamics<T,Descriptor>::unserialize(unserializer);
  }
  
  template<typename T, template<typename U> class Descriptor>
  void IBstandaloneBGKdynamics<T,Descriptor>::defineVelocity(Cell<T,Descriptor>& cell, 
                                                         Array<T,Descriptor<T>::d> const& u)
  {
    Array<T,Descriptor<T>::q> fEq;
    T const rhoBar = 1.;
    T const invRho = 1./rhoBar;
    T const uSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(u); 
    for(plint i=0;i<Descriptor<T>::q;i++)
      cell[i] = dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(i, rhoBar, invRho, u, uSqr);
  }
  
   /// Implementation of the collision step
  template<typename T, template<typename U> class Descriptor>
  void IBstandaloneBGKdynamics<T,Descriptor>::collide(Cell<T,Descriptor>& cell,
                                                  BlockStatistics& statistics)
  {
    // reset all forces to zero, regardless of what is going on
    this->particleData.hydrodynamicForce.resetToZero();

    Array<T,3> j;
    T rhoBar;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T jSqr = normSqr(j);
    T invRho = Descriptor<T>::invRho(rhoBar);

    dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria(rhoBar, invRho, j, jSqr, fEq);

    T const omega = this->getOmega();
    T const fs = this->particleData.solidFraction;
    
    // no solid --> regular BGK collision
    if(fs < SOLFRAC_MIN){
      for(plint i=0;i<Descriptor<T>::q;i++)
        cell[i] += -omega*(cell[i] - fEq[i]);

      return; // exiting because no solid treatment necessary
    }

    // for the moment, only implement "traditional" N-T collision
    Array<T,Descriptor<T>::d> jPart = this->particleData.uPart*(1.+rhoBar);
    T const jPartSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(jPart);
    dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria(rhoBar, invRho, jPart, jPartSqr, fEqSolid);
    
    T const b = calcB(omega,fs);
    T const one_min_b = 1. - b;
    {
      T const coll = b*( (cell[0] - fEq[0]) - (cell[0] - fEqSolid[0]) );
      cell[0] += -one_min_b*omega*(cell[0] - fEq[0]);// + coll;
    }

    for(plint i=1;i<=Descriptor<T>::q/2;i++){
      plint const iOpp = i+Descriptor<T>::q/2;

      T const coll = b*( (cell[iOpp] - fEq[iOpp]) - (cell[i] - fEqSolid[i]) );
      T const collOpp = b*( (cell[i] - fEq[i]) - (cell[iOpp] - fEqSolid[iOpp]) );

      cell[i] += -omega*one_min_b*(cell[i]-fEq[i]) + coll;
      cell[iOpp] += -omega*one_min_b*(cell[iOpp]-fEq[iOpp]) + collOpp;

      for(plint iDim=0;iDim<Descriptor<T>::d;iDim++)
        this->particleData.hydrodynamicForce[iDim]
          -= Descriptor<T>::c[i][iDim]*(coll-collOpp);
      
    }    
  }

  template<typename T, template<typename U> class Descriptor>
  T IBstandaloneBGKdynamics<T,Descriptor>::calcB(T const omega, T const fs)
  {
    return fs*(1./omega-0.5)/((1.-fs)+(1./omega-0.5));
  }
 
  template<typename T, template<typename U> class Descriptor>
  T IBstandaloneBGKdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar,
                                                              Array<T,Descriptor<T>::d> const& j,
                                                              T jSqr, T thetaBar) const
  {
      T invRho = Descriptor<T>::invRho(rhoBar);
      return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
  }
  
}; /* namespace plb */

#endif /* IB_STANDALONE_BGK_DYNAMICS_HH_LBDEM */
