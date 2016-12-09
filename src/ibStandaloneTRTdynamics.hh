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

#ifndef IB_COMPOSITE_DYNAMICS_HH_LBDEM
#define IB_COMPOSITE_DYNAMICS_HH_LBDEM

#include "ibDef.h"

namespace plb {

  // leave this the same because some other parts rely on this information
  template<typename T, template<typename U> class Descriptor>
  int IBstandaloneTRTdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor, IBstandaloneTRTdynamics<T,Descriptor> >("IBcomposite");

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
    particleData(),
    lambda(lambda_)
  { }
  
  template<typename T, template<typename U> class Descriptor>
  IBstandaloneTRTdynamics<T,Descriptor>::IBstandaloneTRTdynamics(HierarchicUnserializer &unserializer)
    : IsoThermalBulkDynamics<T,Descriptor>(1.)
  {
    // pcout << "entering serialize constructor" << std::endl;
    unserialize(unserializer);
  }

  template<typename T, template<typename U> class Descriptor>
  IBstandaloneTRTdynamics<T,Descriptor>::IBstandaloneTRTdynamics(const IBstandaloneTRTdynamics &orig)
    : IsoThermalBulkDynamics<T,Descriptor>(orig),
    particleData(orig.particleData),
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

    particleData.serialize(serializer);
    serializer.addValue(lambda);
    IsoThermalBulkDynamics<T,Descriptor>::serialize(serializer);
  }

  template<typename T, template<typename U> class Descriptor>
  void IBstandaloneTRTdynamics<T,Descriptor>::unserialize(HierarchicUnserializer &unserializer)
  {
    PLB_PRECONDITION( unserializer.getId() == this->getId() );

    particleData.unserialize(unserializer);
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
    particleData.hydrodynamicForce.resetToZero();

    Array<T,3> j;
    T rhoBar;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T jSqr = normSqr(j);
    T invRho = Descriptor<T>::invRho(rhoBar);

    dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria(rhoBar, invRho, j, jSqr, fEq);

    for (plint i=1; i<=Descriptor<T>::q/2; ++i) {
      fEq_plus[i]  = 0.5*(fEq[i] + fEq[i+Descriptor<T>::q/2]);
      fEq_minus[i] = 0.5*(fEq[i] - fEq[i+Descriptor<T>::q/2]);
      f_plus[i]   = 0.5*(cell[i] + cell[i+Descriptor<T>::q/2]);
      f_minus[i]  = 0.5*(cell[i] - cell[i+Descriptor<T>::q/2]);
    }

    T const omega_plus = this->getOmega();
    T const omega_minus = ( 4. - 2.*omega_plus )/( 4.*lambda*omega_plus + 2. - omega_plus );
      
    T const fs = particleData.solidFraction;

    // no solid --> regular TRT collision
    if(particleData.solidFraction < SOLFRAC_MIN){
      cell[0] += -omega_plus*(cell[0] - fEq[0]);
      
      for (plint i=1; i<=Descriptor<T>::q/2; ++i) {
        cell[i] += -omega_plus*(f_plus[i]-fEq_plus[i]) - omega_minus*(f_minus[i]-fEq_minus[i]);
        cell[i+Descriptor<T>::q/2] += -omega_plus*(f_plus[i]-fEq_plus[i]) + omega_minus*(f_minus[i]-fEq_minus[i]);
      }
      return; // exiting because no solid treatment necessary
    }

    // for the moment, only implement "traditional" N-T collision
    Array<T,Descriptor<T>::d> uPart = particleData.uPart*(1.+rhoBar);
    T const uPartSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(uPart); 
    dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria(rhoBar, invRho, uPart, uPartSqr, fEqSolid);
    for (plint i=1; i<=Descriptor<T>::q/2; ++i) {
      fEqSolid_plus[i]  = 0.5*(fEqSolid[i] + fEqSolid[i+Descriptor<T>::q/2]);
      fEqSolid_minus[i] = 0.5*(fEqSolid[i] - fEqSolid[i+Descriptor<T>::q/2]);
    }

    T const b_plus = calcB(omega_plus,fs), b_minus = calcB(omega_minus,fs);
    T const one_min_b_plus = 1. - b_plus, one_min_b_minus = 1. - b_minus;

    // treat f[0] separately
    T coll = b_plus*(fEqSolid[0] - fEq[0]);
    cell[0] += -omega_plus*(cell[0] - fEq[0])*one_min_b_plus + coll;

    for (plint i=1; i<=Descriptor<T>::q/2; ++i) {
      plint const iOp = i+Descriptor<T>::q/2;
      T coll_plus = b_plus*(fEqSolid_plus[i] - fEqSolid[i]);
      T coll_minus = b_minus*(fEqSolid_minus[i] - f_minus[i]) + (fEq[i] - f_minus[i]);

      cell[i] += -omega_plus*(f_plus[i]-fEq_plus[i])*one_min_b_plus
        - omega_minus*(f_minus[i]-fEq_minus[i])*one_min_b_minus;
      cell[i] += coll_plus + coll_minus;
      
      cell[iOp] += -omega_plus*(f_plus[i]-fEq_plus[i])*one_min_b_plus
        + omega_minus*(f_minus[i]-fEq_minus[i])*one_min_b_minus;
      cell[iOp] += coll_plus - coll_minus;

      for(plint iDim=0;iDim<Descriptor<T>::d;iDim++)
        particleData.hydrodynamicForce[iDim] += Descriptor<T>::c[i][iDim]*coll_plus;

    }    
    
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

#endif /* IB_COMPOSITE_DYNAMICS_HH_LBDEM */
