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
 * Geometry specifications for 3D multiblocks -- header file.
 */

#include "multiBlock/localMultiBlockInfo3D.h"

#include <algorithm>

#include "core/globalDefs.h"
#include "multiBlock/multiBlockManagement3D.h"

namespace plb {

PeriodicOverlap3D::PeriodicOverlap3D(
    Overlap3D const &overlap_, plint normalX_, plint normalY_, plint normalZ_) :
    overlap(overlap_), normalX(normalX_), normalY(normalY_), normalZ(normalZ_)
{ }

LocalMultiBlockInfo3D::LocalMultiBlockInfo3D(
    SparseBlockStructure3D const &sparseBlock, ThreadAttribution const &attribution,
    plint envelopeWidth_) :
    envelopeWidth(envelopeWidth_)
{
    computeMyBlocks(sparseBlock, attribution);
    computeAllNormalOverlaps(sparseBlock);
    computeAllPeriodicOverlaps(sparseBlock);
    // This is important: the overlaps must be sorted so they
    //   appear in the same order on different processors, to
    //   guarantee a match in the communication pattern.
    std::sort(normalOverlaps.begin(), normalOverlaps.end());
    std::sort(periodicOverlaps.begin(), periodicOverlaps.end());
}

std::vector<plint> const &LocalMultiBlockInfo3D::getBlocks() const
{
    return myBlocks;
}

std::vector<Overlap3D> const &LocalMultiBlockInfo3D::getNormalOverlaps() const
{
    return normalOverlaps;
}

std::vector<PeriodicOverlap3D> const &LocalMultiBlockInfo3D::getPeriodicOverlaps() const
{
    return periodicOverlaps;
}

std::vector<PeriodicOverlap3D> const &LocalMultiBlockInfo3D::getPeriodicOverlapWithRemoteData()
    const
{
    return periodicOverlapWithRemoteData;
}

void LocalMultiBlockInfo3D::swap(LocalMultiBlockInfo3D &rhs)
{
    std::swap(envelopeWidth, rhs.envelopeWidth);
    myBlocks.swap(rhs.myBlocks);
    normalOverlaps.swap(rhs.normalOverlaps);
    periodicOverlaps.swap(rhs.periodicOverlaps);
    periodicOverlapWithRemoteData.swap(rhs.periodicOverlapWithRemoteData);
}

void LocalMultiBlockInfo3D::computeMyBlocks(
    SparseBlockStructure3D const &sparseBlock, ThreadAttribution const &attribution)
{
    myBlocks = sparseBlock.getLocalBlocks(attribution);
}

void LocalMultiBlockInfo3D::computeAllNormalOverlaps(SparseBlockStructure3D const &sparseBlock)
{
    for (pluint iBlock = 0; iBlock < myBlocks.size(); ++iBlock) {
        plint blockId = myBlocks[iBlock];
        computeNormalOverlaps(sparseBlock, blockId);
    }
}

void LocalMultiBlockInfo3D::computeNormalOverlaps(
    SparseBlockStructure3D const &sparseBlock, plint blockId)
{
    Box3D intersection;
    SmartBulk3D bulk(sparseBlock, envelopeWidth, blockId);

    std::vector<plint> neighbors;
    sparseBlock.findNeighbors(blockId, envelopeWidth, neighbors);

    for (pluint iNeighbor = 0; iNeighbor < neighbors.size(); ++iNeighbor) {
        plint neighborId = neighbors[iNeighbor];
        SmartBulk3D neighborBulk(sparseBlock, envelopeWidth, neighborId);
        if (intersect(neighborBulk.getBulk(), bulk.computeNonPeriodicEnvelope(), intersection)) {
            normalOverlaps.push_back(Overlap3D(neighborId, blockId, intersection));
        }
        if (intersect(bulk.getBulk(), neighborBulk.computeNonPeriodicEnvelope(), intersection)) {
            normalOverlaps.push_back(Overlap3D(blockId, neighborId, intersection));
        }
    }
}

void LocalMultiBlockInfo3D::computeAllPeriodicOverlaps(SparseBlockStructure3D const &sparseBlock)
{
    for (pluint iBlock = 0; iBlock < myBlocks.size(); ++iBlock) {
        plint blockId = myBlocks[iBlock];
        Box3D bulk;
        sparseBlock.getBulk(blockId, bulk);
        // Speed optimization: execute the test for periodicity
        //   only for bulk-domains which touch the bounding box.
        if (!contained(bulk.enlarge(1), sparseBlock.getBoundingBox())) {
            computePeriodicOverlaps(sparseBlock, blockId);
        }
    }
}

void LocalMultiBlockInfo3D::computePeriodicOverlaps(
    SparseBlockStructure3D const &sparseBlock, plint blockId)
{
    Box3D intersection;            // Temporary variable.
    std::vector<plint> neighbors;  // Temporary variable.
    SmartBulk3D bulk(sparseBlock, envelopeWidth, blockId);

    for (plint dx = -1; dx <= +1; dx += 1) {
        for (plint dy = -1; dy <= +1; dy += 1) {
            for (plint dz = -1; dz <= +1; dz += 1) {
                if (dx != 0 || dy != 0 || dz != 0) {
                    // The new block is shifted by the length of the full multi block in each space
                    //   direction. Consequently, overlaps between the original multi block and the
                    //   shifted new block are identified as periodic overlaps.
                    plint shiftX = dx * sparseBlock.getBoundingBox().getNx();
                    plint shiftY = dy * sparseBlock.getBoundingBox().getNy();
                    plint shiftZ = dz * sparseBlock.getBoundingBox().getNz();
                    Box3D shiftedBulk(bulk.getBulk().shift(shiftX, shiftY, shiftZ));
                    Box3D shiftedEnvelope(bulk.computeEnvelope().shift(shiftX, shiftY, shiftZ));
                    // Speed optimization: perform following checks only if the shifted
                    //   domain touches the bounding box.
                    Box3D dummyIntersection;
                    if (intersect(shiftedEnvelope, sparseBlock.getBoundingBox(), dummyIntersection))
                    {
                        neighbors.clear();
                        sparseBlock.findNeighbors(shiftedBulk, envelopeWidth, neighbors);
                        // Check overlap with each existing block in the neighborhood, including
                        // with the newly added one.
                        for (pluint iNeighbor = 0; iNeighbor < neighbors.size(); ++iNeighbor) {
                            plint neighborId = neighbors[iNeighbor];
                            SmartBulk3D neighborBulk(sparseBlock, envelopeWidth, neighborId);
                            // Does the envelope of the shifted new block overlap with the bulk of a
                            // previous
                            //   block? If yes, add an overlap, in which the previous block has the
                            //   "original position", and the new block has the "overlap position".
                            if (intersect(neighborBulk.getBulk(), shiftedEnvelope, intersection)) {
                                PeriodicOverlap3D overlap(
                                    Overlap3D(
                                        neighborId, blockId, intersection, shiftX, shiftY, shiftZ),
                                    dx, dy, dz);
                                periodicOverlaps.push_back(overlap);
                                periodicOverlapWithRemoteData.push_back(overlap);
                            }
                            // Does the bulk of the shifted new block overlap with the envelope of a
                            // previous
                            //   block? If yes, add an overlap, in which the new block has the
                            //   "original position", and the previous block has the "overlap
                            //   position". If we are in the situation in which the newly added
                            //   block is periodic with itself, this step must be skipped, because
                            //   otherwise the overlap is counted twice.
                            if (!(neighborId == blockId)
                                && intersect(
                                    shiftedBulk, neighborBulk.computeEnvelope(), intersection))
                            {
                                intersection = intersection.shift(-shiftX, -shiftY, -shiftZ);
                                periodicOverlaps.push_back(PeriodicOverlap3D(
                                    Overlap3D(
                                        blockId, neighborId, intersection, -shiftX, -shiftY,
                                        -shiftZ),
                                    -dx, -dy, -dz));
                            }
                        }
                    }
                }
            }
        }
    }
}

}  // namespace plb
