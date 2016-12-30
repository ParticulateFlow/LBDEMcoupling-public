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

#ifndef IB_DYNAMICS_PARTICLE_DATA_H_LBDEM
#define IB_DYNAMICS_PARTICLE_DATA_H_LBDEM

namespace plb {

  template<typename T, template<typename U> class Descriptor>
  struct IBdynamicsParticleData {
  public:
    IBdynamicsParticleData()
      : partId(0), solidFraction(0.),
        uPart(0.,0.,0.), hydrodynamicForce(0.,0.,0.)
    {}
    plint partId;
    T solidFraction;
    Array<T,Descriptor<T>::d> uPart;
    Array<T,Descriptor<T>::d> hydrodynamicForce;
  };

}; /* namespace plb */ 


#endif /* IB_DYNAMICS_PARTICLE_DATA_H_LBDEM */
