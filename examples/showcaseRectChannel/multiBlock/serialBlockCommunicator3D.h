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
 * Serial version of the 3D block communicator -- header file.
 */
#ifndef SERIAL_BLOCK_COMMUNICATOR_3D_H
#define SERIAL_BLOCK_COMMUNICATOR_3D_H

#include "core/globalDefs.h"
#include "multiBlock/blockCommunicator3D.h"
#include "multiBlock/multiBlockManagement3D.h"

namespace plb {

class MultiBlockManagement3D;

class SerialBlockCommunicator3D : public BlockCommunicator3D {
public:
    SerialBlockCommunicator3D();
    virtual SerialBlockCommunicator3D *clone() const;
    virtual void communicate(
        std::vector<Overlap3D> const &overlaps, MultiBlock3D const &originMultiBlock,
        MultiBlock3D &destinationMultiBlock, modif::ModifT whichData) const;
    virtual void duplicateOverlaps(MultiBlock3D &multiBlock, modif::ModifT whichData) const;
    virtual void signalPeriodicity() const;

private:
    void copyOverlap(
        Overlap3D const &overlap, MultiBlock3D const &fromMultiBlock, MultiBlock3D &toMultiBlock,
        modif::ModifT whichData) const;
};

}  // namespace plb

#endif  // SERIAL_BLOCK_COMMUNICATOR_3D_H
