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
 * Helper classes for serial 2D multiblock lattice -- header file.
 */
#ifndef SERIAL_MULTI_BLOCK_LATTICE_2D_H
#define SERIAL_MULTI_BLOCK_LATTICE_2D_H

#include <vector>

#include "core/globalDefs.h"
#include "multiBlock/multiBlockLattice2D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class SerialCellAccess2D : public MultiCellAccess2D<T, Descriptor> {
public:
    SerialCellAccess2D();
    virtual void broadCastCell(
        Cell<T, Descriptor> &cell, plint fromBlock,
        MultiBlockManagement2D const &multiBlockManagement) const;
    Cell<T, Descriptor> &getDistributedCell(
        plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
        std::map<plint, BlockLattice2D<T, Descriptor> *> &lattices);
    Cell<T, Descriptor> const &getDistributedCell(
        plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
        std::map<plint, BlockLattice2D<T, Descriptor> *> const &lattices) const;
    virtual SerialCellAccess2D<T, Descriptor> *clone() const;

private:
    mutable plint locatedBlock;
};

}  // namespace plb

#endif  // SERIAL_MULTI_BLOCK_LATTICE_2D_H
