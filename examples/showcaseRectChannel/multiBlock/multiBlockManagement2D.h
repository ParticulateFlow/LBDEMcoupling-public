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

#ifndef MULTI_BLOCK_MANAGEMENT_2D_H
#define MULTI_BLOCK_MANAGEMENT_2D_H

#include <vector>

#include "core/geometry2D.h"
#include "core/globalDefs.h"
#include "multiBlock/localMultiBlockInfo2D.h"
#include "multiBlock/sparseBlockStructure2D.h"
#include "multiBlock/threadAttribution.h"

namespace plb {

class MultiBlockManagement2D {
public:
    MultiBlockManagement2D(
        SparseBlockStructure2D const &sparseBlock_, ThreadAttribution *threadAttribution_,
        plint envelopeWidth_, plint refinementLevel_ = 0);
    MultiBlockManagement2D(MultiBlockManagement2D const &rhs);
    MultiBlockManagement2D &operator=(MultiBlockManagement2D const &rhs);
    void swap(MultiBlockManagement2D &rhs);
    ~MultiBlockManagement2D();
    plint getEnvelopeWidth() const;
    Box2D getBoundingBox() const;
    Box2D getBulk(plint blockId) const;
    Box2D getUniqueBulk(plint blockId) const;
    Box2D getEnvelope(plint blockId) const;
    SparseBlockStructure2D const &getSparseBlockStructure() const;
    LocalMultiBlockInfo2D const &getLocalInfo() const;
    ThreadAttribution const &getThreadAttribution() const;
    bool findInLocalBulk(plint iX, plint iY, plint &foundId, plint &localX, plint &localY) const;
    bool findAllLocalRepresentations(
        plint iX, plint iY, std::vector<plint> &foundId, std::vector<plint> &foundX,
        std::vector<plint> &foundY) const;
    plint getRefinementLevel() const;
    void setRefinementLevel(plint newLevel);
    void changeEnvelopeWidth(plint newEnvelopeWidth);
    // Same multi-block-management, except for envelope-width
    bool equivalentTo(MultiBlockManagement2D const &rhs) const;

private:
    plint envelopeWidth;
    SparseBlockStructure2D sparseBlock;
    ThreadAttribution *threadAttribution;
    LocalMultiBlockInfo2D localInfo;
    plint refinementLevel;
};

/// Create a new block-management, corresponding to a sub-domain of the old one.
/** If the parameter crop is true, the bounding-box of the new block-management
 * is equal to the specified sub-domain. If crop is false, the bounding-box is
 * the same as the bounding-box of the original block-management.
 */
MultiBlockManagement2D intersect(
    MultiBlockManagement2D const &originalManagement, Box2D subDomain, bool crop);

/// Choose yourself the bounding box of the resulting block-management.
MultiBlockManagement2D intersect(
    MultiBlockManagement2D const &originalManagement, Box2D subDomain, Box2D newBoundingBox);

/// Create a new block-management as an intersection of the old ones.
/** The bounding box of the result is the intersection of the two original bounding
 *  boxes if crop=true, and the bound of the two otherwise. The block-ids and
 *  the thread-attribution are inherited from management1.
 */
MultiBlockManagement2D intersect(
    MultiBlockManagement2D const &management1, MultiBlockManagement2D const &management2,
    bool crop);

/// Create a new block-management which extends the original one by a given block.
/** Added blocks are default-associated to the mpiProcess "bossId".
 **/
MultiBlockManagement2D extend(
    MultiBlockManagement2D const &management, Box2D addedBulk, Box2D addedUniqueBulk);

/// Create a new block-management from which a given block is excepted.
MultiBlockManagement2D except(MultiBlockManagement2D const &management, Box2D exceptedBlock);

/// Union of two block-management structures.
/** The block-ids and thread-attribution of management1 are replicated without
 *  modification in the result.
 **/
MultiBlockManagement2D block_union(
    MultiBlockManagement2D const &management1, MultiBlockManagement2D const &management2);

/// Align parallelization of originalManagement so it can be coupled through
///   data processors with partnerManagement.
MultiBlockManagement2D align(
    MultiBlockManagement2D const &originalManagement,
    MultiBlockManagement2D const &partnerManagement);

MultiBlockManagement2D align(
    std::vector<Box2D> const &originalDomain, MultiBlockManagement2D const &alignWith,
    plint envelopeWidth, plint refinementLevel, bool crop = true);

/// Re-create a block-management by covering the sparse structure with regular blocks.
/** The parameters blockLx and blockLy indicate the approximate size of the blocks.
 **/
MultiBlockManagement2D reparallelize(
    MultiBlockManagement2D const &management, plint blockLx, plint blockLy);

/// Compute envelope and things alike.
class SmartBulk2D {
public:
    /** The bulk has zero coordinates in case the block id does not exist. **/
    SmartBulk2D(MultiBlockManagement2D const &management, plint blockId);
    /** The bulk has zero coordinates in case the block id does not exist. **/
    SmartBulk2D(SparseBlockStructure2D const &sparseBlock_, plint envelopeWidth_, plint blockId);
    SmartBulk2D(
        SparseBlockStructure2D const &sparseBlock_, plint envelopeWidth_, Box2D const &bulk_);
    /// Access the bulk.
    Box2D getBulk() const;
    /// Compute envelope of a given block.
    Box2D computeEnvelope() const;
    /// Compute envelope of a given block, exluding margins of the outer domain.
    Box2D computeNonPeriodicEnvelope() const;
    /// Convert to local coordinates of a given block.
    Box2D toLocal(Box2D const &coord) const;
    /// Convert to local x-coordinate of a given block.
    plint toLocalX(plint iX) const;
    /// Convert to local y-coordinate of a given block.
    plint toLocalY(plint iY) const;

private:
    SparseBlockStructure2D const &sparseBlock;
    plint envelopeWidth;
    Box2D bulk;
};

}  // namespace plb

#endif  // MULTI_BLOCK_MANAGEMENT_2D_H
