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
 * Helper classes for serial 3D multiblock lattice -- generic implementation.
 */
#ifndef SERIAL_MULTI_BLOCK_LATTICE_3D_HH
#define SERIAL_MULTI_BLOCK_LATTICE_3D_HH

#include "atomicBlock/blockLattice3D.h"
#include "multiBlock/multiBlockManagement3D.h"
#include "multiBlock/serialMultiBlockLattice3D.h"

namespace plb {

////////////////////// Class SerialCellAccess3D /////////////////////

template <typename T, template <typename U> class Descriptor>
SerialCellAccess3D<T, Descriptor>::SerialCellAccess3D() : locatedBlock(0)
{ }

template <typename T, template <typename U> class Descriptor>
void SerialCellAccess3D<T, Descriptor>::broadCastCell(
    Cell<T, Descriptor> &cell, plint fromBlock,
    MultiBlockManagement3D const &multiBlockManagement) const
{
    // Nothing to do in the serial case
}

template <typename T, template <typename U> class Descriptor>
Cell<T, Descriptor> &SerialCellAccess3D<T, Descriptor>::getDistributedCell(
    plint iX, plint iY, plint iZ, MultiBlockManagement3D const &multiBlockManagement,
    std::map<plint, BlockLattice3D<T, Descriptor> *> &lattices)
{
    plint localX, localY, localZ;
#ifdef PLB_DEBUG
    bool ok =
#endif
        multiBlockManagement.findInLocalBulk(iX, iY, iZ, locatedBlock, localX, localY, localZ);
    PLB_PRECONDITION(ok);
    return lattices[locatedBlock]->get(localX, localY, localZ);
}

template <typename T, template <typename U> class Descriptor>
Cell<T, Descriptor> const &SerialCellAccess3D<T, Descriptor>::getDistributedCell(
    plint iX, plint iY, plint iZ, MultiBlockManagement3D const &multiBlockManagement,
    std::map<plint, BlockLattice3D<T, Descriptor> *> const &lattices) const
{
    plint localX, localY, localZ;
#ifdef PLB_DEBUG
    bool ok =
#endif
        multiBlockManagement.findInLocalBulk(iX, iY, iZ, locatedBlock, localX, localY, localZ);
    PLB_PRECONDITION(ok);
    return lattices.find(locatedBlock)->second->get(localX, localY, localZ);
}

template <typename T, template <typename U> class Descriptor>
SerialCellAccess3D<T, Descriptor> *SerialCellAccess3D<T, Descriptor>::clone() const
{
    return new SerialCellAccess3D<T, Descriptor>;
}

}  // namespace plb

#endif  // SERIAL_MULTI_BLOCK_LATTICE_3D_HH
