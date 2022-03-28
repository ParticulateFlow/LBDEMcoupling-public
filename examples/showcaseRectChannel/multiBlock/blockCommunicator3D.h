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
 * Block Communicator -- Abstract base class.
 */
#ifndef BLOCK_COMMUNICATOR_3D_H
#define BLOCK_COMMUNICATOR_3D_H

#include <vector>

#include "core/globalDefs.h"

namespace plb {

class MultiBlock3D;
class MultiBlockManagement3D;
class Overlap3D;

struct BlockCommunicator3D {
    virtual ~BlockCommunicator3D() { }
    virtual BlockCommunicator3D *clone() const = 0;
    /// Fill the overlaps (the "envelopes") with data from the corresponding bulks.
    /** The variable whichData specifies which type of content (static/dynamic/full dynamics object)
     *  is being transmitted.
     **/
    virtual void duplicateOverlaps(MultiBlock3D &multiBlock, modif::ModifT whichData) const = 0;
    /// Transmit data between two multi-blocks, according to a user-defined pattern.
    /** The variable whichData specifies which type of content (static/dynamic/full dynamics object)
     *  is being transmitted.
     **/
    virtual void communicate(
        std::vector<Overlap3D> const &overlaps, MultiBlock3D const &originMultiBlock,
        MultiBlock3D &destinationMultiBlock, modif::ModifT whichData) const = 0;
    virtual void signalPeriodicity() const = 0;
};

}  // namespace plb

#endif  // BLOCK_COMMUNICATOR_3D_H
