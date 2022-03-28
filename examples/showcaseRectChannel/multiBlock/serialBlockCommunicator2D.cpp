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
 * Serial version of the 2D block communicator -- implementation.
 */
#include "multiBlock/serialBlockCommunicator2D.h"

#include "atomicBlock/atomicBlock2D.h"
#include "core/plbDebug.h"
#include "multiBlock/multiBlock2D.h"

namespace plb {

////////////////////// Class SerialBlockCommunicator2D /////////////////////

SerialBlockCommunicator2D::SerialBlockCommunicator2D() { }

SerialBlockCommunicator2D *SerialBlockCommunicator2D::clone() const
{
    return new SerialBlockCommunicator2D;
}

void SerialBlockCommunicator2D::copyOverlap(
    Overlap2D const &overlap, MultiBlock2D const &fromMultiBlock, MultiBlock2D &toMultiBlock,
    modif::ModifT whichData) const
{
    MultiBlockManagement2D const &fromManagement = fromMultiBlock.getMultiBlockManagement();
    MultiBlockManagement2D const &toManagement = toMultiBlock.getMultiBlockManagement();
    plint fromEnvelopeWidth = fromManagement.getEnvelopeWidth();
    plint toEnvelopeWidth = toManagement.getEnvelopeWidth();
    SparseBlockStructure2D const &fromSparseBlock = fromManagement.getSparseBlockStructure();
    SparseBlockStructure2D const &toSparseBlock = toManagement.getSparseBlockStructure();
    plint originalId = overlap.getOriginalId();
    plint overlapId = overlap.getOverlapId();
    SmartBulk2D originalBulk(fromSparseBlock, fromEnvelopeWidth, originalId);
    SmartBulk2D overlapBulk(toSparseBlock, toEnvelopeWidth, overlapId);

    Box2D originalCoords(originalBulk.toLocal(overlap.getOriginalCoordinates()));
    Box2D overlapCoords(overlapBulk.toLocal(overlap.getOverlapCoordinates()));

    PLB_PRECONDITION(originalCoords.x1 - originalCoords.x0 == overlapCoords.x1 - overlapCoords.x0);
    PLB_PRECONDITION(originalCoords.y1 - originalCoords.y0 == overlapCoords.y1 - overlapCoords.y0);

    AtomicBlock2D const *originalBlock = &fromMultiBlock.getComponent(originalId);
    AtomicBlock2D *overlapBlock = &toMultiBlock.getComponent(overlapId);
    plint deltaX = originalCoords.x0 - overlapCoords.x0;
    plint deltaY = originalCoords.y0 - overlapCoords.y0;

    overlapBlock->getDataTransfer().attribute(
        overlapCoords, deltaX, deltaY, *originalBlock, whichData);
}

void SerialBlockCommunicator2D::duplicateOverlaps(
    MultiBlock2D &multiBlock, modif::ModifT whichData) const
{
    MultiBlockManagement2D const &multiBlockManagement = multiBlock.getMultiBlockManagement();
    LocalMultiBlockInfo2D const &localInfo = multiBlockManagement.getLocalInfo();

    // Non-periodic communication
    for (pluint iOverlap = 0; iOverlap < localInfo.getNormalOverlaps().size(); ++iOverlap) {
        copyOverlap(localInfo.getNormalOverlaps()[iOverlap], multiBlock, multiBlock, whichData);
    }

    // Periodic communication
    PeriodicitySwitch2D const &periodicity = multiBlock.periodicity();
    for (pluint iOverlap = 0; iOverlap < localInfo.getPeriodicOverlaps().size(); ++iOverlap) {
        PeriodicOverlap2D const &pOverlap = localInfo.getPeriodicOverlaps()[iOverlap];
        if (periodicity.get(pOverlap.normalX, pOverlap.normalY)) {
            copyOverlap(pOverlap.overlap, multiBlock, multiBlock, whichData);
        }
    }
}

void SerialBlockCommunicator2D::communicate(
    std::vector<Overlap2D> const &overlaps, MultiBlock2D const &originMultiBlock,
    MultiBlock2D &destinationMultiBlock, modif::ModifT whichData) const
{
    for (pluint iOverlap = 0; iOverlap < overlaps.size(); ++iOverlap) {
        copyOverlap(overlaps[iOverlap], originMultiBlock, destinationMultiBlock, whichData);
    }
}

void SerialBlockCommunicator2D::signalPeriodicity() const { }

}  // namespace plb
