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

/*
 * composite dynamics that implement the method of Noble and Torczinsky
 */

#ifndef IB_STANDALONE_TRT_DYNAMICS_H_LBDEM
#define IB_STANDALONE_TRT_DYNAMICS_H_LBDEM

#include "ibDynamicsParticleData.h"

namespace plb {

  template<typename T, template<typename U> class Descriptor>
  class IBstandaloneTRTdynamics : public IsoThermalBulkDynamics<T,Descriptor>,
                                  public IBdynamicsParticleData <T,Descriptor>
  {
  public:

    IBstandaloneTRTdynamics(T const omega_, T const lambda_= 3./16.);
    IBstandaloneTRTdynamics(const IBstandaloneTRTdynamics &orig);
    IBstandaloneTRTdynamics(HierarchicUnserializer &unserializer);
    ~IBstandaloneTRTdynamics();

    virtual IBstandaloneTRTdynamics<T,Descriptor>* clone() const;    
    
    virtual int getId() const;

    virtual void serialize(HierarchicSerializer& serializer) const;
    virtual void unserialize(HierarchicUnserializer& unserializer);

    virtual void collide(Cell<T,Descriptor>& cell,  BlockStatistics& statistics);

    virtual void defineVelocity(Cell<T,Descriptor>& cell, 
                                Array<T,Descriptor<T>::d> const& u);

    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar) const;
  private:
    static int id;
    T lambda;

    static Array<T,Descriptor<T>::q> fPostColl, fEq, fEqSolid;
    static Array<T,Descriptor<T>::q/2+1> fEq_plus, fEq_minus,
      fEqSolid_plus, fEqSolid_minus,
      f_plus, f_minus;

    T calcB(T const omega, T const fs) { return fs*(1./omega-0.5)/((1.-fs)-(1./omega-0.5)); }
  };

}; /* namespace plb */

#include "ibStandaloneTRTdynamics.hh"

#endif /* IB_STANDALONE_TRT_DYNAMICS_H_LBDEM */
