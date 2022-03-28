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
 * Helper classes for serial 2D multiblock lattice -- generic implementation.
 */
#ifndef SERIAL_MULTI_BLOCK_LATTICE_2D_HH
#define SERIAL_MULTI_BLOCK_LATTICE_2D_HH

#include "atomicBlock/blockLattice2D.h"
#include "multiBlock/multiBlockManagement2D.h"
#include "multiBlock/serialMultiBlockLattice2D.h"

namespace plb {

////////////////////// Class SerialCellAccess2D /////////////////////

template <typename T, template <typename U> class Descriptor>
SerialCellAccess2D<T, Descriptor>::SerialCellAccess2D() : locatedBlock(0)
{ }

template <typename T, template <typename U> class Descriptor>
void SerialCellAccess2D<T, Descriptor>::broadCastCell(
    Cell<T, Descriptor> &cell, plint fromBlock,
    MultiBlockManagement2D const &multiBlockManagement) const
{
    // Nothing to do in the serial case
}

template <typename T, template <typename U> class Descriptor>
Cell<T, Descriptor> &SerialCellAccess2D<T, Descriptor>::getDistributedCell(
    plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
    std::map<plint, BlockLattice2D<T, Descriptor> *> &lattices)
{
    plint localX, localY;
#ifdef PLB_DEBUG
    bool ok =
#endif
        multiBlockManagement.findInLocalBulk(iX, iY, locatedBlock, localX, localY);
    PLB_PRECONDITION(ok);
    return lattices[locatedBlock]->get(localX, localY);
}

template <typename T, template <typename U> class Descriptor>
Cell<T, Descriptor> const &SerialCellAccess2D<T, Descriptor>::getDistributedCell(
    plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
    std::map<plint, BlockLattice2D<T, Descriptor> *> const &lattices) const
{
    plint localX, localY;
#ifdef PLB_DEBUG
    bool ok =
#endif
        multiBlockManagement.findInLocalBulk(iX, iY, locatedBlock, localX, localY);
    PLB_PRECONDITION(ok);
    return lattices.find(locatedBlock)->second->get(localX, localY);
}

template <typename T, template <typename U> class Descriptor>
SerialCellAccess2D<T, Descriptor> *SerialCellAccess2D<T, Descriptor>::clone() const
{
    return new SerialCellAccess2D<T, Descriptor>;
}

}  // namespace plb

#endif  // SERIAL_MULTI_BLOCK_LATTICE_2D_HH
