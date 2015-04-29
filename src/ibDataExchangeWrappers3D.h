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
 * some wrapper functions for data exchange between LIGGGHTS and Palabos
 */

#ifndef IBDATAEXCHANGEWRAPPERS_H_LBDEM
#define IBDATAEXCHANGEWRAPPERS_H_LBDEM

#include "ibProcessors3D.h"

namespace plb {
  /*
   * if initVelFlag is set true, the populations at the grid nodes covered with 
   * particles are initialized with the particle velocity. This is useful to avoid
   * force spikes when a particle already has a certain velocity upon beginning of
   * the simulation. During normal simulation, it should be set to false
   */
  template<typename T, template<typename U> class Descriptor>
  void setSpheresOnLattice(MultiBlockLattice3D<T,Descriptor> &lattice,
                           LiggghtsCouplingWrapper &wrapper,
                           PhysUnits3D<T> const &units,
                           std::vector<plint> &excludeType,
                           bool initVelFlag);

  template<typename T, template<typename U> class Descriptor>
  void setSpheresOnLattice(MultiBlockLattice3D<T,Descriptor> &lattice,
                           LiggghtsCouplingWrapper &wrapper,
                           PhysUnits3D<T> const &units,
                           bool initVelFlag);

  template<typename T, template<typename U> class Descriptor>
  void getForcesFromLattice(MultiBlockLattice3D<T,Descriptor> &lattice,
                            LiggghtsCouplingWrapper &wrapper,
                            PhysUnits3D<T> const &units);
  
}; /* namespace plb */

#include "ibDataExchangeWrappers3D.hh"

#endif /* IBDATAEXCHANGEWRAPPERS_H_LBDEM */
