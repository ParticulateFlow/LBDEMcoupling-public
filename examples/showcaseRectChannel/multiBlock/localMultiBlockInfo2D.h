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
 * Geometry specifications for 2D multiblocks -- header file.
 */

#ifndef LOCAL_MULTI_BLOCK_INFO_2D_H
#define LOCAL_MULTI_BLOCK_INFO_2D_H

#include <vector>

#include "core/geometry2D.h"
#include "core/globalDefs.h"
#include "multiBlock/sparseBlockStructure2D.h"
#include "multiBlock/threadAttribution.h"

namespace plb {

/// Hold extra information on the blocks which are local to the current MPI thread;
///  for example, overlaps with adjacent blocks.
class Overlap2D {
public:
    Overlap2D(plint originalId_, plint overlapId_, Box2D const &intersection_) :
        originalId(originalId_),
        overlapId(overlapId_),
        originalRegion(intersection_),
        overlapRegion(intersection_)
    { }
    Overlap2D(
        plint originalId_, plint overlapId_, Box2D const &originalRegion_, plint shiftX,
        plint shiftY) :
        originalId(originalId_),
        overlapId(overlapId_),
        originalRegion(originalRegion_),
        overlapRegion(originalRegion.shift(-shiftX, -shiftY))
    { }
    plint getOriginalId() const
    {
        return originalId;
    }
    plint getOverlapId() const
    {
        return overlapId;
    }
    /// Region (in absolute coordinates) on the original block.
    Box2D const &getOriginalCoordinates() const
    {
        return originalRegion;
    }
    /// Region (in absolute coordinates) on the overlapping block.
    /** This is usually identical with the region on the original block. An
     *  exception are periodic overlaps, in whick regions on opposite ends
     *  of the block are brought into relation.
     **/
    Box2D const &getOverlapCoordinates() const
    {
        return overlapRegion;
    }
    plint getShiftX() const
    {
        return originalRegion.x0 - overlapRegion.x0;
    }
    plint getShiftY() const
    {
        return originalRegion.y0 - overlapRegion.y0;
    }

private:
    plint originalId, overlapId;
    Box2D originalRegion, overlapRegion;
};

/// Define a global ordering for overlaps.
/** This can be used for example to guarantee that the MPI communication
 *  between a pair of processes is executed in the same order, and the
 *  communications don't cross.
 * */
inline bool operator<(Overlap2D const &overlap1, Overlap2D const &overlap2)
{
    return (overlap1.getOriginalId() < overlap2.getOriginalId())
           || ((overlap1.getOriginalId() == overlap2.getOriginalId())
               && ((overlap1.getOverlapId() < overlap2.getOverlapId())
                   || ((overlap1.getOverlapId() == overlap2.getOverlapId())
                       && (overlap1.getOriginalCoordinates()
                           < overlap2.getOriginalCoordinates()))));
}

/// This structure holds both overlap information and orientation of the boundary.
/** In case of periodic overlaps, it is important to know the orientation of the
 *  boundary, additionally to the coordinates of the overlap region. This is
 *  required when the communication step within a multi block is executed. Given
 *  that the user can selectively swith on/off periodicity, the multi block
 *  must be able to decide which periodic overlaps to communicate and which not.
 */
struct PeriodicOverlap2D {
    PeriodicOverlap2D(Overlap2D const &overlap_, plint normalX_, plint normalY_);
    Overlap2D overlap;
    plint normalX;
    plint normalY;
};

/// Define a global ordering for periodic overlaps.
inline bool operator<(PeriodicOverlap2D const &overlap1, PeriodicOverlap2D const &overlap2)
{
    return overlap1.overlap < overlap2.overlap;
}

/// Determine pairs of domains associated to a data transfer between two blocks.
std::vector<Overlap2D> copyAllDataTransfer(
    SparseBlockStructure2D const &block1, SparseBlockStructure2D const &block2);

/// Determine pairs of domains associated to a data transfer between domains on two blocks.
/** It is assumed that the two domains have the same extent.
 **/
std::vector<Overlap2D> copyDomainDataTransfer(
    SparseBlockStructure2D const &block1, Box2D block1Domain, SparseBlockStructure2D const &block2,
    Box2D block2Domain);

class LocalMultiBlockInfo2D {
public:
    LocalMultiBlockInfo2D(
        SparseBlockStructure2D const &sparseBlock, ThreadAttribution const &attribution,
        plint envelopeWidth_);
    /// Index of all blocks local to current processor
    std::vector<plint> const &getBlocks() const;
    /// Index of all overlaps for which original or overlap data are on current processor
    std::vector<Overlap2D> const &getNormalOverlaps() const;
    /// Index of all periodic overlaps for which original or overlap data are
    ///   on current processor.
    std::vector<PeriodicOverlap2D> const &getPeriodicOverlaps() const;
    /// Index of all periodic overlaps for which overlap data are on current processor
    std::vector<PeriodicOverlap2D> const &getPeriodicOverlapWithRemoteData() const;
    void swap(LocalMultiBlockInfo2D &rhs);

private:
    /// Determine all blocks which are associated to the current MPI thread.
    void computeMyBlocks(
        SparseBlockStructure2D const &sparseBlock, ThreadAttribution const &attribution);
    /// Compute normal overlaps for all local blocks.
    void computeAllNormalOverlaps(SparseBlockStructure2D const &sparseBlock);
    /// Compute normal overlaps for one local block.
    void computeNormalOverlaps(SparseBlockStructure2D const &sparseBlock, plint blockId);
    /// Compute periodic overlaps for all local blocks.
    void computeAllPeriodicOverlaps(SparseBlockStructure2D const &sparseBlock);
    /// Compute periodic overlaps for one local block.
    void computePeriodicOverlaps(SparseBlockStructure2D const &sparseBlock, plint blockId);

private:
    plint envelopeWidth;
    std::vector<plint> myBlocks;
    std::vector<Overlap2D> normalOverlaps;
    std::vector<PeriodicOverlap2D> periodicOverlaps;
    std::vector<PeriodicOverlap2D> periodicOverlapWithRemoteData;
};

}  // namespace plb

#endif  // LOCAL_MULTI_BLOCK_INFO_2D_H
