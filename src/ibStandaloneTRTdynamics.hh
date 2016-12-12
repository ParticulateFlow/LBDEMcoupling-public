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

#ifndef IB_STANDALONE_TRT_DYNAMICS_HH_LBDEM
#define IB_STANDALONE_TRT_DYNAMICS_HH_LBDEM

#include "ibDef.h"

namespace plb {

  // leave this the same because some other parts rely on this information
  template<typename T, template<typename U> class Descriptor>
  int IBstandaloneTRTdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor, IBstandaloneTRTdynamics<T,Descriptor> >("IBstandaloneTRT");

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q> IBstandaloneTRTdynamics<T,Descriptor>::fEq =
    Array<T,Descriptor<T>::q>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q> IBstandaloneTRTdynamics<T,Descriptor>::fEqSolid =
    Array<T,Descriptor<T>::q>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q> IBstandaloneTRTdynamics<T,Descriptor>::fPostColl =
    Array<T,Descriptor<T>::q>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q/2+1> IBstandaloneTRTdynamics<T,Descriptor>::f_plus
    = Array<T,Descriptor<T>::q/2+1>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q/2+1> IBstandaloneTRTdynamics<T,Descriptor>::f_minus
    = Array<T,Descriptor<T>::q/2+1>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q/2+1> IBstandaloneTRTdynamics<T,Descriptor>::fEq_plus
    = Array<T,Descriptor<T>::q/2+1>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q/2+1> IBstandaloneTRTdynamics<T,Descriptor>::fEq_minus
    = Array<T,Descriptor<T>::q/2+1>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q/2+1> IBstandaloneTRTdynamics<T,Descriptor>::fEqSolid_plus
    = Array<T,Descriptor<T>::q/2+1>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q/2+1> IBstandaloneTRTdynamics<T,Descriptor>::fEqSolid_minus
    = Array<T,Descriptor<T>::q/2+1>();

  template<typename T, template<typename U> class Descriptor>
  IBstandaloneTRTdynamics<T,Descriptor>::IBstandaloneTRTdynamics(T const omega_, T const lambda_)
    : IsoThermalBulkDynamics<T,Descriptor>(omega_), 
    lambda(lambda_)
  { }
  
  template<typename T, template<typename U> class Descriptor>
  IBstandaloneTRTdynamics<T,Descriptor>::IBstandaloneTRTdynamics(HierarchicUnserializer &unserializer)
    : IsoThermalBulkDynamics<T,Descriptor>(1.)
  {
    unserialize(unserializer);
  }

  template<typename T, template<typename U> class Descriptor>
  IBstandaloneTRTdynamics<T,Descriptor>::IBstandaloneTRTdynamics(const IBstandaloneTRTdynamics &orig)
    : IsoThermalBulkDynamics<T,Descriptor>(orig),
    IBdynamicsParticleData<T,Descriptor>(orig),
    lambda(orig.lambda)
  { }
  
  template<typename T, template<typename U> class Descriptor>
  IBstandaloneTRTdynamics<T,Descriptor>::~IBstandaloneTRTdynamics() {}
  
  template<typename T, template<typename U> class Descriptor>
  IBstandaloneTRTdynamics<T,Descriptor>* IBstandaloneTRTdynamics<T,Descriptor>::clone() const {
    return new IBstandaloneTRTdynamics<T,Descriptor>(*this);
  }
  
  template<typename T, template<typename U> class Descriptor>
  int IBstandaloneTRTdynamics<T,Descriptor>::getId() const
  {
    return id;
  }
  
  template<typename T, template<typename U> class Descriptor>
  void IBstandaloneTRTdynamics<T,Descriptor>::serialize(HierarchicSerializer &serializer) const
  {

    this->particleData.serialize(serializer);
    serializer.addValue(lambda);
    IsoThermalBulkDynamics<T,Descriptor>::serialize(serializer);
  }

  template<typename T, template<typename U> class Descriptor>
  void IBstandaloneTRTdynamics<T,Descriptor>::unserialize(HierarchicUnserializer &unserializer)
  {
    PLB_PRECONDITION( unserializer.getId() == this->getId() );

    this->particleData.unserialize(unserializer);
    unserializer.readValue<T>(lambda);
    IsoThermalBulkDynamics<T,Descriptor>::unserialize(unserializer);
  }
  
  template<typename T, template<typename U> class Descriptor>
  void IBstandaloneTRTdynamics<T,Descriptor>::defineVelocity(Cell<T,Descriptor>& cell, 
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
  void IBstandaloneTRTdynamics<T,Descriptor>::collide(Cell<T,Descriptor>& cell,
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

    T const omega_plus = this->getOmega();
    T const omega_minus =  ( 4. - 2.*omega_plus )/( 4.*lambda*omega_plus + 2. - omega_plus);

    // pcout << "omega_plus " << omega_plus << " | omega_minus" << omega_minus << std::endl;
    
    T const fs = this->particleData.solidFraction;

    for(plint i=1;i<=Descriptor<T>::q/2;i++){
      plint const iOpp = i+Descriptor<T>::q/2;
      f_plus[i] = 0.5*(cell[i]+cell[iOpp]);
      f_minus[i] = 0.5*(cell[i]-cell[iOpp]);
      fEq_plus[i] = 0.5*(fEq[i]+fEq[iOpp]);
      fEq_minus[i] = 0.5*(fEq[i]-fEq[iOpp]);
    }
    
    // no solid --> regular TRT collision
    if(fs < SOLFRAC_MIN){
      cell[0] += -omega_plus*(cell[0]-fEq[0]);
      for(plint i=1;i<=Descriptor<T>::q/2;i++){
        plint const iOpp = i + Descriptor<T>::q/2;
        cell[i] += -omega_plus*(f_plus[i] - fEq_plus[i]) - omega_minus*(f_minus[i] - fEq_minus[i]);
        cell[iOpp] += -omega_plus*(f_plus[i] - fEq_plus[i]) + omega_minus*(f_minus[i] - fEq_minus[i]);
      }
      return; // exiting because no solid treatment necessary
    }

    // shear and bulk viscosity
    // from Kruger's book, p. 422
    T const nu = Descriptor<T>::cs2*(1./omega_plus - 0.5);
    T const eta = Descriptor<T>::cs2*(1./omega_minus - 0.5) - nu/3.;
    // pcout << "nu " << nu << " | eta " << eta << std::endl;
    // for the moment, only implement "traditional" N-T collision
    Array<T,Descriptor<T>::d> jPart = this->particleData.uPart*(1.+rhoBar);
    T const jPartSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(jPart);
    dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria(rhoBar, invRho, jPart, jPartSqr, fEqSolid);
    
    for(plint i=1;i<=Descriptor<T>::q/2;i++){
      plint const iOpp = i + Descriptor<T>::q/2;
      fEqSolid_plus[i] = 0.5*(fEqSolid[i]+fEqSolid[iOpp]);
      fEqSolid_minus[i] = 0.5*(fEqSolid[i]-fEqSolid[iOpp]);
    }

    T const b_plus = fs*nu*Descriptor<T>::invCs2 / ( (1.-fs) + nu*Descriptor<T>::invCs2 );
    T const b_minus = fs*eta*Descriptor<T>::invCs2 / ( (1.-fs) + eta*Descriptor<T>::invCs2 );
    // T const b_plus = fs*(1./omega_plus-0.5) / ( (1.-fs) + (1./omega_plus-0.5) );
    // T const b_minus = fs*(1./omega_minus-0.5) / ( (1.-fs) + (1./omega_minus-0.5) );
    T const one_min_b_plus = 1. - b_plus;
    T const one_min_b_minus = 1. - b_minus;

    // pcout << "fs " << fs << " | b_plus " << b_plus << " | b_minus " << b_minus << std::endl;
    
    {
      T const coll = b_plus*( (cell[0] - fEq[0]) - (cell[0] - fEqSolid[0]) );
      cell[0] += -one_min_b_plus*omega_plus*(cell[0] - fEq[0]);// + coll;
    }

    for(plint i=1;i<=Descriptor<T>::q/2;i++){
      plint const iOpp = i+Descriptor<T>::q/2;

      T const coll = b_plus*( (cell[iOpp] - fEq[iOpp]) - (cell[i] - fEqSolid[i]) );
      T const collOpp = b_plus*( (cell[i] - fEq[i]) - (cell[iOpp] - fEqSolid[iOpp]) );

      cell[i] += -omega_plus*one_min_b_plus*(f_plus[i] - fEq_plus[i])
        - omega_minus*one_min_b_plus*(f_minus[i] - fEq_minus[i]) + coll;
      cell[iOpp] += -omega_plus*one_min_b_plus*(f_plus[i] - fEq_plus[i])
        + omega_minus*one_min_b_plus*(f_minus[i] - fEq_minus[i]) + collOpp;

      for(plint iDim=0;iDim<Descriptor<T>::d;iDim++)
        this->particleData.hydrodynamicForce[iDim]
          -= Descriptor<T>::c[i][iDim]*(coll-collOpp);
      
    }    

    // {
    //   T const coll = b_plus*( fEqSolid[0] - fEq[0] );
    //   cell[0] += -one_min_b_plus*omega_plus*(cell[0] - fEq[0]) + coll;
    // }

    // for(plint i=1;i<=Descriptor<T>::q/2;i++){
    //   plint const iOpp = i + Descriptor<T>::q/2;

    //   // // collision part
    //   // cell[i] += -omega_plus*one_min_b_plus*(f_plus[i] - fEq_plus[i])
    //   //   - omega_minus*one_min_b_minus*(f_minus[i] - fEq_minus[i]);
    //   // cell[iOpp] += -omega_plus*one_min_b_plus*(f_plus[i] - fEq_plus[i])
    //   //   + omega_minus*one_min_b_minus*(f_minus[i] - fEq_minus[i]);


    //   // // bounceback part
    //   // T const coll_plus = fEqSolid_plus[i] - fEq_plus[i];
    //   // T const coll_minus = (fEqSolid_minus[i]  - f_minus[i]) + (fEq_minus[i]  - f_minus[i]);
    //   // cell[i] += b_plus*coll_plus + b_minus*coll_minus;
    //   // cell[iOpp] += b_plus*coll_plus - b_minus*coll_minus;
      
    //   // for(plint iDim=0;iDim<Descriptor<T>::d;iDim++)
    //   //   this->particleData.hydrodynamicForce[iDim]
    //   //     -= Descriptor<T>::c[i][iDim]*2.*coll_minus;
    //   //     // -= Descriptor<T>::c[i][iDim]*2.*coll_plus;
      
    // }    
  }

  template<typename T, template<typename U> class Descriptor>
  T IBstandaloneTRTdynamics<T,Descriptor>::calcB(T const omega, T const fs)
  {
    return fs*(1./omega-0.5)/((1.-fs)+(1./omega-0.5));
  }
 
  template<typename T, template<typename U> class Descriptor>
  T IBstandaloneTRTdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar,
                                                              Array<T,Descriptor<T>::d> const& j,
                                                              T jSqr, T thetaBar) const
  {
      T invRho = Descriptor<T>::invRho(rhoBar);
      return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
  }
  
}; /* namespace plb */

#endif /* IB_STANDALONE_TRT_DYNAMICS_HH_LBDEM */
