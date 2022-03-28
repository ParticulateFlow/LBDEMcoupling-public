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

#ifndef LOCAL_MULTI_BLOCK_INFO_3D_H
#define LOCAL_MULTI_BLOCK_INFO_3D_H

#include <vector>

#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "multiBlock/sparseBlockStructure3D.h"
#include "multiBlock/threadAttribution.h"

namespace plb {

/// Hold extra information on the blocks which are local to the current MPI thread;
///  for example, overlaps with adjacent blocks.
class Overlap3D {
public:
    Overlap3D(plint originalId_, plint overlapId_, Box3D const &intersection_) :
        originalId(originalId_),
        overlapId(overlapId_),
        originalRegion(intersection_),
        overlapRegion(intersection_)
    { }
    Overlap3D(
        plint originalId_, plint overlapId_, Box3D const &originalRegion_, plint shiftX,
        plint shiftY, plint shiftZ) :
        originalId(originalId_),
        overlapId(overlapId_),
        originalRegion(originalRegion_),
        overlapRegion(originalRegion.shift(-shiftX, -shiftY, -shiftZ))
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
    Box3D const &getOriginalCoordinates() const
    {
        return originalRegion;
    }
    /// Region (in absolute coordinates) on the overlapping block.
    /** This is usually identical with the region on the original block. An
     *  exception are periodic overlaps, in which regions on opposite ends
     *  of the block are brought into relation.
     **/
    Box3D const &getOverlapCoordinates() const
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
    plint getShiftZ() const
    {
        return originalRegion.z0 - overlapRegion.z0;
    }

private:
    plint originalId, overlapId;
    Box3D originalRegion, overlapRegion;
};

/// Define a global ordering for overlaps.
/** This can be used for example to guarantee that the MPI communication
 *  between a pair of processes is executed in the same order, and the
 *  communications don't cross.
 * */
inline bool operator<(Overlap3D const &overlap1, Overlap3D const &overlap2)
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
struct PeriodicOverlap3D {
    PeriodicOverlap3D(Overlap3D const &overlap_, plint normalX_, plint normalY_, plint normalZ_);
    Overlap3D overlap;
    plint normalX;
    plint normalY;
    plint normalZ;
};

/// Define a global ordering for periodic overlaps.
inline bool operator<(PeriodicOverlap3D const &overlap1, PeriodicOverlap3D const &overlap2)
{
    return overlap1.overlap < overlap2.overlap;
}

class LocalMultiBlockInfo3D {
public:
    LocalMultiBlockInfo3D(
        SparseBlockStructure3D const &sparseBlock, ThreadAttribution const &attribution,
        plint envelopeWidth_);
    /// Index of all blocks local to current processor
    std::vector<plint> const &getBlocks() const;
    /// Index of all overlaps for which original or overlap data are on current processor
    std::vector<Overlap3D> const &getNormalOverlaps() const;
    /// Index of all periodic overlaps for which original or overlap data are
    ///   on current processor.
    std::vector<PeriodicOverlap3D> const &getPeriodicOverlaps() const;
    /// Index of all periodic overlaps for which overlap data are on current processor
    std::vector<PeriodicOverlap3D> const &getPeriodicOverlapWithRemoteData() const;
    void swap(LocalMultiBlockInfo3D &rhs);

private:
    /// Determine all blocks which are associated to the current MPI thread.
    void computeMyBlocks(
        SparseBlockStructure3D const &sparseBlock, ThreadAttribution const &attribution);
    /// Compute normal overlaps for all local blocks.
    void computeAllNormalOverlaps(SparseBlockStructure3D const &sparseBlock);
    /// Compute normal overlaps for one local block.
    void computeNormalOverlaps(SparseBlockStructure3D const &sparseBlock, plint blockId);
    /// Compute periodic overlaps for all local blocks.
    void computeAllPeriodicOverlaps(SparseBlockStructure3D const &sparseBlock);
    /// Compute periodic overlaps for one local block.
    void computePeriodicOverlaps(SparseBlockStructure3D const &sparseBlock, plint blockId);

private:
    plint envelopeWidth;
    std::vector<plint> myBlocks;
    std::vector<Overlap3D> normalOverlaps;
    std::vector<PeriodicOverlap3D> periodicOverlaps;
    std::vector<PeriodicOverlap3D> periodicOverlapWithRemoteData;
};

}  // namespace plb

#endif  // LOCAL_MULTI_BLOCK_INFO_3D_H
