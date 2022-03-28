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

/** \file
 * BounceBack dynamics models in 2D -- header file.
 */

#ifndef BOUNCE_BACK_MODELS_2D_H
#define BOUNCE_BACK_MODELS_2D_H

#include "boundaryCondition/bounceBackModels.h"
#include "dataProcessors/dataInitializerFunctional2D.h"

namespace plb {

template <typename T, template <class U> class Descriptor>
void initializeMomentumExchange(BlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <class U> class Descriptor>
void initializeMomentumExchange(
    BlockLattice2D<T, Descriptor> &lattice, Box2D boundingBox, DomainFunctional2D *domain);

template <typename T, template <class U> class Descriptor>
void initializeMomentumExchange(BlockLattice2D<T, Descriptor> &lattice, DotList2D const &dotList);

template <typename T, template <class U> class Descriptor>
void initializeMomentumExchange(MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <class U> class Descriptor>
void initializeMomentumExchange(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D boundingBox, DomainFunctional2D *domain);

template <typename T, template <class U> class Descriptor>
void initializeMomentumExchange(
    MultiBlockLattice2D<T, Descriptor> &lattice, DotList2D const &dotList);

template <typename T, template <typename U> class Descriptor>
class SetAverageWallDensityOnVelocityBounceBack2D :
    public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual SetAverageWallDensityOnVelocityBounceBack2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

}  // namespace plb

#endif  // BOUNCE_BACK_MODELS_2D_H
