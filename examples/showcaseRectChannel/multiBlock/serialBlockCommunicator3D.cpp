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
 * Serial version of the 3D block communicator -- implementation.
 */
#include "multiBlock/serialBlockCommunicator3D.h"

#include "atomicBlock/atomicBlock3D.h"
#include "core/plbDebug.h"
#include "multiBlock/multiBlock3D.h"

namespace plb {

////////////////////// Class SerialBlockCommunicator3D /////////////////////

SerialBlockCommunicator3D::SerialBlockCommunicator3D() { }

SerialBlockCommunicator3D *SerialBlockCommunicator3D::clone() const
{
    return new SerialBlockCommunicator3D;
}

void SerialBlockCommunicator3D::copyOverlap(
    Overlap3D const &overlap, MultiBlock3D const &fromMultiBlock, MultiBlock3D &toMultiBlock,
    modif::ModifT whichData) const
{
    MultiBlockManagement3D const &fromManagement = fromMultiBlock.getMultiBlockManagement();
    MultiBlockManagement3D const &toManagement = toMultiBlock.getMultiBlockManagement();
    plint fromEnvelopeWidth = fromManagement.getEnvelopeWidth();
    plint toEnvelopeWidth = toManagement.getEnvelopeWidth();
    SparseBlockStructure3D const &fromSparseBlock = fromManagement.getSparseBlockStructure();
    SparseBlockStructure3D const &toSparseBlock = toManagement.getSparseBlockStructure();
    plint originalId = overlap.getOriginalId();
    plint overlapId = overlap.getOverlapId();
    SmartBulk3D originalBulk(fromSparseBlock, fromEnvelopeWidth, originalId);
    SmartBulk3D overlapBulk(toSparseBlock, toEnvelopeWidth, overlapId);

    Box3D originalCoords(originalBulk.toLocal(overlap.getOriginalCoordinates()));
    Box3D overlapCoords(overlapBulk.toLocal(overlap.getOverlapCoordinates()));

    PLB_PRECONDITION(originalCoords.x1 - originalCoords.x0 == overlapCoords.x1 - overlapCoords.x0);
    PLB_PRECONDITION(originalCoords.y1 - originalCoords.y0 == overlapCoords.y1 - overlapCoords.y0);
    PLB_PRECONDITION(originalCoords.z1 - originalCoords.z0 == overlapCoords.z1 - overlapCoords.z0);

    AtomicBlock3D const *originalBlock = &fromMultiBlock.getComponent(originalId);
    AtomicBlock3D *overlapBlock = &toMultiBlock.getComponent(overlapId);
    plint deltaX = originalCoords.x0 - overlapCoords.x0;
    plint deltaY = originalCoords.y0 - overlapCoords.y0;
    plint deltaZ = originalCoords.z0 - overlapCoords.z0;

    overlapBlock->getDataTransfer().attribute(
        overlapCoords, deltaX, deltaY, deltaZ, *originalBlock, whichData);
}

void SerialBlockCommunicator3D::duplicateOverlaps(
    MultiBlock3D &multiBlock, modif::ModifT whichData) const
{
    MultiBlockManagement3D const &multiBlockManagement = multiBlock.getMultiBlockManagement();
    LocalMultiBlockInfo3D const &localInfo = multiBlockManagement.getLocalInfo();

    // Non-periodic communication
    for (pluint iOverlap = 0; iOverlap < localInfo.getNormalOverlaps().size(); ++iOverlap) {
        copyOverlap(localInfo.getNormalOverlaps()[iOverlap], multiBlock, multiBlock, whichData);
    }

    // Periodic communication
    PeriodicitySwitch3D const &periodicity = multiBlock.periodicity();
    for (pluint iOverlap = 0; iOverlap < localInfo.getPeriodicOverlaps().size(); ++iOverlap) {
        PeriodicOverlap3D const &pOverlap = localInfo.getPeriodicOverlaps()[iOverlap];
        if (periodicity.get(pOverlap.normalX, pOverlap.normalY, pOverlap.normalZ)) {
            copyOverlap(pOverlap.overlap, multiBlock, multiBlock, whichData);
        }
    }
}

void SerialBlockCommunicator3D::communicate(
    std::vector<Overlap3D> const &overlaps, MultiBlock3D const &originMultiBlock,
    MultiBlock3D &destinationMultiBlock, modif::ModifT whichData) const
{
    for (pluint iOverlap = 0; iOverlap < overlaps.size(); ++iOverlap) {
        copyOverlap(overlaps[iOverlap], originMultiBlock, destinationMultiBlock, whichData);
    }
}

void SerialBlockCommunicator3D::signalPeriodicity() const { }

}  // namespace plb
