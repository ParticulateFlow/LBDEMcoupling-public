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


#ifndef UTILS_H_LBDEM
#define UTILS_H_LBDEM

#include "ibDynamicsParticleData.h"

namespace plb {
  template<typename T, template<typename U> class Descriptor>
  IBdynamicsParticleData<T,Descriptor>* getParticleDataFromCell(Cell<T,Descriptor> &cell);

  
}; /* namespace plb */

// --------------------------
// implementaton starts here

namespace plb {
  template<typename T, template<typename U> class Descriptor>
  IBdynamicsParticleData<T,Descriptor>* getParticleDataFromCell(Cell<T,Descriptor> &cell)
  {
    static IBcompositeDynamics<T,Descriptor> const cdyn(new NoDynamics<T,Descriptor>(),false);
    static plint const cdynId = cdyn.getId();

    Dynamics<T,Descriptor> *dyn = &( cell.getDynamics() );

    if(dyn->getId() == cdynId){
      return &( (static_cast< IBcompositeDynamics<T,Descriptor>* >(dyn))->particleData );
    } 

    while(dyn->isComposite()){
      dyn = &(static_cast<CompositeDynamics<T,Descriptor>* >(dyn))->getBaseDynamics();
      if(dyn->getId() == cdynId)
        return &( (static_cast< IBcompositeDynamics<T,Descriptor>* >(dyn))->particleData );
    }

    return 0;
  }

}; /* namespace plb */

#endif /* UTILS_H_LBDEM */
