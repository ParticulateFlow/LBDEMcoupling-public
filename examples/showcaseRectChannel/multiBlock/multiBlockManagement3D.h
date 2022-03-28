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

#ifndef MULTI_BLOCK_MANAGEMENT_3D_H
#define MULTI_BLOCK_MANAGEMENT_3D_H

#include <vector>

#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "multiBlock/localMultiBlockInfo3D.h"
#include "multiBlock/sparseBlockStructure3D.h"
#include "multiBlock/threadAttribution.h"
#include "multiGrid/multiScale.h"

namespace plb {

class MultiBlockManagement3D {
public:
    MultiBlockManagement3D(
        SparseBlockStructure3D const &sparseBlock_, ThreadAttribution *threadAttribution_,
        plint envelopeWidth_, plint refinementLevel_ = 0);
    MultiBlockManagement3D(MultiBlockManagement3D const &rhs);
    MultiBlockManagement3D &operator=(MultiBlockManagement3D const &rhs);
    void swap(MultiBlockManagement3D &rhs);
    ~MultiBlockManagement3D();
    plint getEnvelopeWidth() const;
    Box3D getBoundingBox() const;
    Box3D getBulk(plint blockId) const;
    Box3D getUniqueBulk(plint blockId) const;
    Box3D getEnvelope(plint blockId) const;
    SparseBlockStructure3D const &getSparseBlockStructure() const;
    LocalMultiBlockInfo3D const &getLocalInfo() const;
    ThreadAttribution const &getThreadAttribution() const;
    void setCoProcessors(std::map<plint, int> const &coProcessors);
    bool findInLocalBulk(
        plint iX, plint iY, plint iZ, plint &foundId, plint &localX, plint &localY,
        plint &localZ) const;
    bool findAllLocalRepresentations(
        plint iX, plint iY, plint iZ, std::vector<plint> &foundId, std::vector<plint> &foundX,
        std::vector<plint> &foundY, std::vector<plint> &foundZ) const;
    plint getRefinementLevel() const;
    void setRefinementLevel(plint newLevel);
    void changeEnvelopeWidth(plint newEnvelopeWidth);
    // Same multi-block-management, except for envelope-width
    bool equivalentTo(MultiBlockManagement3D const &rhs) const;

private:
    plint envelopeWidth;
    SparseBlockStructure3D sparseBlock;
    ThreadAttribution *threadAttribution;
    LocalMultiBlockInfo3D localInfo;
    plint refinementLevel;
};

MultiBlockManagement3D scale(MultiBlockManagement3D const &originalManagement, plint relativeLevel);

/// Create a new block-management, corresponding to a sub-domain of the old one.
/** If the parameter crop is true, the bounding-box of the new block-management
 * is equal to the specified sub-domain. If crop is false, the bounding-box is
 * the same as the bounding-box of the original block-management.
 */
MultiBlockManagement3D intersect(
    MultiBlockManagement3D const &originalManagement, Box3D subDomain, bool crop);

/// Choose yourself the bounding box of the resulting block-management.
MultiBlockManagement3D intersect(
    MultiBlockManagement3D const &originalManagement, Box3D subDomain, Box3D newBoundingBox);

/// Create a new block-management as an intersection of the old ones.
/** The bounding box of the result is the intersection of the two original bounding
 *  boxes if crop=true, and the bound of the two otherwise. The block-ids and
 *  the thread-attribution are inherited from management1.
 */
MultiBlockManagement3D intersect(
    MultiBlockManagement3D const &management1, MultiBlockManagement3D const &management2,
    bool crop);

/// Create a new block-management which extends the original one by a given block.
/** Added blocks are default-associated to the mpiProcess "bossId".
 **/
MultiBlockManagement3D extend(
    MultiBlockManagement3D const &management, Box3D addedBulk, Box3D addedUniqueBulk);

/// Create a new block-management from which a given block is excepted.
MultiBlockManagement3D except(MultiBlockManagement3D const &management, Box3D exceptedBlock);

/// Union of two block-management structures.
/** The block-ids and thread-attribution of management1 are replicated without
 *  modification in the result.
 **/
MultiBlockManagement3D block_union(
    MultiBlockManagement3D const &management1, MultiBlockManagement3D const &management2);

/// Align parallelization of originalManagement so it can be coupled through
///   data processors with partnerManagement.
MultiBlockManagement3D align(
    MultiBlockManagement3D const &originalManagement,
    MultiBlockManagement3D const &partnerManagement);

MultiBlockManagement3D align(
    std::vector<Box3D> const &originalDomain, MultiBlockManagement3D const &alignWith,
    plint envelopeWidth, plint refinementLevel, bool crop = true);

/// Re-create a block-management by covering the sparse structure with regular blocks.
/** The parameters blockLx, blockLy, and blockLz indicate the approximate size of the
 *  blocks.
 **/
MultiBlockManagement3D reparallelize(
    MultiBlockManagement3D const &management, plint blockLx, plint blockLy, plint blockLz);

/// Compute envelope and things alike.
class SmartBulk3D {
public:
    /** The bulk has zero coordinates in case the block id does not exist. **/
    SmartBulk3D(MultiBlockManagement3D const &management, plint blockId);
    /** The bulk has zero coordinates in case the block id does not exist. **/
    SmartBulk3D(SparseBlockStructure3D const &sparseBlock_, plint envelopeWidth_, plint blockId);
    SmartBulk3D(
        SparseBlockStructure3D const &sparseBlock_, plint envelopeWidth_, Box3D const &bulk_);
    /// Access the bulk.
    Box3D getBulk() const;
    /// Compute envelope of a given block.
    Box3D computeEnvelope() const;
    /// Compute envelope of a given block, exluding margins of the outer domain.
    Box3D computeNonPeriodicEnvelope() const;
    /// Convert to local coordinates of a given block.
    Box3D toLocal(Box3D const &coord) const;
    /// Convert to local x-coordinate of a given block.
    plint toLocalX(plint iX) const;
    /// Convert to local y-coordinate of a given block.
    plint toLocalY(plint iY) const;
    /// Convert to local z-coordinate of a given block.
    plint toLocalZ(plint iZ) const;

private:
    SparseBlockStructure3D const &sparseBlock;
    plint envelopeWidth;
    Box3D bulk;
};

}  // namespace plb

#endif  // MULTI_BLOCK_MANAGEMENT_3D_H
