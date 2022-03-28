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
 * Fast algorithm for handling sparse multi-block structure -- header file.
 */

#ifndef SPARSE_BLOCK_STRUCTURE_2D_H
#define SPARSE_BLOCK_STRUCTURE_2D_H

#include <map>
#include <vector>

#include "core/geometry2D.h"
#include "core/globalDefs.h"
#include "multiBlock/threadAttribution.h"

namespace plb {

class SparseBlockStructure2D {
public:
    typedef std::map<Dot2D, std::vector<plint> > GridT;

public:
    /// Sparse grid structure with default internal implementation.
    SparseBlockStructure2D(plint nx, plint ny);
    /// Sparse grid structure with default internal implementation.
    SparseBlockStructure2D(Box2D boundingBox_);
    /// Sparse grid structure with explicit grid parameters.
    SparseBlockStructure2D(Box2D boundingBox_, plint gridNx_, plint gridNy_);
    /// Restrict an existing SparseBlockStructure2D to a sub-domain.
    SparseBlockStructure2D(SparseBlockStructure2D const &rhs, Box2D boundingBox_);
    /// Add a new block to the sparse block-structure.
    void addBlock(Box2D const &bulk, Box2D const &uniqueBulk, plint blockId);
    /// Add a new block to the sparse block-structure; uniqueBulk = bulk.
    void addBlock(Box2D const &bulk, plint blockId);
    /// Remove an existing block from the sparse block-structure.
    void removeBlock(plint blockId);
    /// Check if a block with the given ID already exists in the block-structure.
    bool exists(plint blockId);
    /// Return 1 + the maximum block ID currently found in the structure.
    plint nextIncrementalId() const;
    /// Return outer bounding box.
    Box2D getBoundingBox() const;
    /// Get the bulk of a given block. Returns false if the block does not exist.
    bool getBulk(plint blockId, Box2D &bulk) const;
    /// Get the unique bulk of a given block. Returns false if the block does not exist.
    bool getUniqueBulk(plint blockId, Box2D &uniqueBulk) const;
    /// Return the number of blocks contained in the structure.
    plint getNumBlocks() const;
    /// Get the number of bulk cells, cumulated over all blocks.
    plint getNumBulkCells() const;
    /// Get the id-to-bulk map for all blocks.
    std::map<plint, Box2D> const &getBulks() const;
    /// Enumerate all blocks attributed to the current MPI thread.
    std::vector<plint> getLocalBlocks(ThreadAttribution const &attribution) const;
    /// Find a block which contains given coordinates. Return -1 in case no block is found.
    plint locate(plint iX, plint iY) const;
    /// Intersect a given block with all blocks of the sparse block-structure.
    void intersect(
        Box2D const &bulk, std::vector<plint> &ids, std::vector<Box2D> &intersections) const;
    /// Intersect the envelope of block "blockId" with all neighbors
    /// in a given direction.
    void computeOverlaps(
        Box2D const &bulk, plint dx, plint dy, std::vector<plint> &ids,
        std::vector<Box2D> &overlapsOnBulk, std::vector<Box2D> &overlapsOnNeighbors) const;
    /// Intersect the envelope of block "blockId" with all neighbors.
    void computeOverlaps(
        plint blockId, plint envelopeWidth, std::vector<plint> &ids,
        std::vector<Box2D> &overlapsOnBulk, std::vector<Box2D> &overlapsOnNeighbors) const;
    /// Find all blocks adjacent to a given area (a specific block can be
    ///   excluded from the search.
    void findNeighbors(
        Box2D const &block, plint neighborhoodWidth, std::vector<plint> &neighbors,
        plint excludeId = -1) const;
    /// Find all blocks adjacent to a given block of the structure.
    void findNeighbors(plint blockId, plint neighborhoodWidth, std::vector<plint> &neighbors) const;
    void swap(SparseBlockStructure2D &rhs);
    bool equals(SparseBlockStructure2D const &rhs) const;

private:
    /// Default resolution of the sparse grid.
    void defaultGridN();
    /// Compute gridLx, gridLy.
    void iniGridParameters();
    /// Convert block x-coordinate into coordinate of the sparse-block grid.
    plint gridPosX(plint realX) const;
    /// Convert block y-coordinate into coordinate of the sparse-block grid.
    plint gridPosY(plint realY) const;
    /// Convert block coordinates into coordinates of the sparse-block grid.
    Box2D getGridBox(Box2D const &realBlock) const;
    /// Extend bulk by an envelope layer in a given direction, in view of
    ///   computing overlaps with neighbors.
    void computeEnvelopeTerm(
        plint block0, plint block1, plint &env0, plint &env1, plint delta) const;
    /// Integrate block into the sparse-block indirect addressing scheme.
    void integrateBlock(plint blockId, Box2D bulk);
    /// Remove block from the sparse-block indirect addressing scheme.
    void extractBlock(plint blockId);

private:
    Box2D boundingBox;
    plint gridLx, gridLy;
    plint gridNx, gridNy;
    GridT grid;
    // Attention: If replacing the map by a hashed_map, remember that
    // nextIncrementalId() uses the fact that elements are ordered inside
    // the map. Therefore, nextIncrementalId() must then be rewritten.
    std::map<plint, Box2D> bulks;
    /// UniqueBulks is never used for internal algorithms. It's just here to be
    ///   provided to the user when asked for.
    std::map<plint, Box2D> uniqueBulks;
};

/// The bounding box of the result is the intersection of bounding box and
///   domain if crop=true, and is identified with the original bounding box
///   otherwise.
SparseBlockStructure2D intersect(
    SparseBlockStructure2D const &sparseBlock, Box2D domain, bool crop);

/// Choose yourself the bounding box of the resulting block structure.
SparseBlockStructure2D intersect(
    SparseBlockStructure2D const &sparseBlock, Box2D domain, Box2D newBoundingBox);

/// The bounding box of the result is the intersection of the two original bounding
///   boxes if crop=true, and the bound of the two otherwise.
SparseBlockStructure2D intersect(
    SparseBlockStructure2D const &sparseBlock1, SparseBlockStructure2D const &sparseBlock2,
    bool crop);

/// Union of two sparse block structures.
SparseBlockStructure2D block_union(
    SparseBlockStructure2D const &sparseBlock1, SparseBlockStructure2D const &sparseBlock2,
    std::map<plint, std::vector<plint> > &remappedIds);

/// By default, the freshly added domain is associated to the main mpi process in
///   the areas which don't intersect with the original domain. The added block is
///   potentially split into non-overlapping components, and the ids of these new
///   components is returned in the vector newIds.
SparseBlockStructure2D extend(
    SparseBlockStructure2D const &sparseBlock, Box2D addedBulk, Box2D addedUniqueBulk,
    std::vector<plint> &newIds);

/// Make a hole in the original block structure.
SparseBlockStructure2D except(
    SparseBlockStructure2D const &sparseBlock, Box2D exceptedBlock,
    std::map<plint, std::vector<plint> > &remappedIds);

/// Align the distribution of originalStructure on partnerStructure, so that it
///   becomes possible to execute couplings between partnerStructure and the
///   newly created structure. The returned newIds refer to newly created blocks
///   (domains of originalStructure which don't overlap with partnerStructure) that
///   need to be parallelized manually. The returned remappedFromPartner refers to
///   newly created block which overlap with partnerStructure and should be parallelized
///   correspondingly.
SparseBlockStructure2D alignDistribution2D(
    SparseBlockStructure2D const &originalStructure, SparseBlockStructure2D const &partnerStructure,
    std::vector<plint> &newIds, std::map<plint, std::vector<plint> > &remappedFromPartner);

/// Iterate in a structured way over a sparse multi-block structure.
class EuclideanIterator2D {
public:
    EuclideanIterator2D(SparseBlockStructure2D const &sparseBlock_);
    /// Locate block and determine how much data can be
    ///     accessed in a x-contiguous manner. Return value indicates if
    ///     the memory corresponding to the chunk is allocated or not.
    bool getNextChunkX(plint iX, plint iY, plint &blockId, plint &chunkSize) const;
    /// Locate block and determine how much data can be
    ///     accessed in a y-contiguous manner. Return value indicates if
    ///     the memory corresponding to the chunk is allocated or not.
    bool getNextChunkY(plint iX, plint iY, plint &blockId, plint &chunkSize) const;

private:
    SparseBlockStructure2D const &sparseBlock;
};

}  // namespace plb

#endif  // SPARSE_BLOCK_STRUCTURE_2D_H
