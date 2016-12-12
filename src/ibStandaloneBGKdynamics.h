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

#ifndef IB_STANDALONE_BGK_DYNAMICS_H_LBDEM
#define IB_STANDALONE_BGK_DYNAMICS_H_LBDEM

#include "ibDynamicsParticleData.h"

namespace plb {

  template<typename T, template<typename U> class Descriptor>
  class IBstandaloneBGKdynamics : public IsoThermalBulkDynamics<T,Descriptor>,
                                  public IBdynamicsParticleData <T,Descriptor>
  {
  public:

    IBstandaloneBGKdynamics(T const omega_);
    IBstandaloneBGKdynamics(const IBstandaloneBGKdynamics &orig);
    IBstandaloneBGKdynamics(HierarchicUnserializer &unserializer);
    ~IBstandaloneBGKdynamics();

    virtual IBstandaloneBGKdynamics<T,Descriptor>* clone() const;    
    
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

    static Array<T,Descriptor<T>::q> fEq, fEqSolid;

    T calcB(T const omega, T const fs);
  };

}; /* namespace plb */

#include "ibStandaloneBGKdynamics.hh"

#endif /* IB_STANDALONE_BGK_DYNAMICS_H_LBDEM */
