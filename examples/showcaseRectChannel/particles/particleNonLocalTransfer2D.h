/* This file is part of the Palabos library.
 *
 * The Palabos softare is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
 *
 * Contact:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
 *
 * The most recent release of Palabos can be downloaded at
 * <https://palabos.unige.ch/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PARTICLE_NON_LOCAL_TRANSFER_2D_H
#define PARTICLE_NON_LOCAL_TRANSFER_2D_H

#include "core/globalDefs.h"
#include "multiBlock/localMultiBlockInfo2D.h"
#include "particles/multiParticleField2D.h"

namespace plb {

template <class ParticleFieldT>
void copy(
    MultiParticleField2D<ParticleFieldT> const &from, Box2D const &fromDomain,
    MultiParticleField2D<ParticleFieldT> &to, Box2D const &toDomain);

/// Gather particles from a particle-field into a vector of particles,
/// on the main processor only.
template <typename T, template <typename U> class Descriptor, class ParticleFieldT>
void gatherParticles(
    MultiParticleField2D<ParticleFieldT> &particleField,
    std::vector<Particle2D<T, Descriptor> *> &particles, Box2D domain);

/// Inject particle which are only provided on the main processor into
/// a particle-field.
template <typename T, template <typename U> class Descriptor, class ParticleFieldT>
void injectParticlesAtMainProc(
    std::vector<Particle2D<T, Descriptor> *> &particles,
    MultiParticleField2D<ParticleFieldT> &particleField, Box2D domain);

}  // namespace plb

#endif  // PARTICLE_NON_LOCAL_TRANSFER_2D_H
