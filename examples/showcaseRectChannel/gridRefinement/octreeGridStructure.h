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

#ifndef OCTREE_GRID_STRUCTURE_H
#define OCTREE_GRID_STRUCTURE_H

#include <map>
#include <vector>

#include "core/array.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "multiBlock/multiBlockManagement3D.h"

namespace plb {

class OctreeGridStructure {
public:
    OctreeGridStructure();
    OctreeGridStructure(std::string xmlFileName);

    void addBlock(plint blockId, Box3D const &bulk, plint level, plint processId, bool isOverlap);
    void addNeighborBlockIds(plint blockId, Array<plint, 26> const &ids);

    bool removeBlock(plint blockId, bool removeConnectivityInfo);
    bool removeNeighborBlockIds(plint blockId);

    plint getNumLevels() const;
    plint getNumProcesses() const;
    Box3D getBoundingBox(plint level) const;
    // Compute a closed cover (full bounding box of all levels) in the units of "level".
    Box3D getClosedCover(plint level) const;

    bool neighborIsAtSameLevel(plint blockId, plint direction) const;
    bool neighborIsAtCoarserLevel(plint blockId, plint direction) const;
    bool neighborIsAtFinerLevel(plint blockId, plint direction) const;
    bool neighborIsBoundary(plint blockId, plint direction) const;
    bool neighborIsAllocated(plint blockId, plint direction) const;

    void getBlock(plint blockId, Box3D &bulk, plint &level, plint &processId) const;
    Array<plint, 26> getNeighborBlockIds(plint blockId) const;
    bool isOverlap(plint blockId) const;

    std::vector<plint> getBlockIdsAtLevel(plint level, bool includeOverlaps) const;
    std::vector<plint> getOverlapBlockIdsAtLevel(plint level) const;

    MultiBlockManagement3D getMultiBlockManagement(
        plint level, Box3D const &boundingBox, plint envelopeWidth = 1) const;
    MultiBlockManagement3D getMultiBlockManagement(plint level, plint envelopeWidth = 1) const;

    MultiBlockManagement3D getMultiBlockManagementForOutput(
        plint level, Box3D const &boundingBox, bool crop, plint envelopeWidth = 1) const;
    MultiBlockManagement3D getMultiBlockManagementForOutput(
        plint level, bool crop, plint envelopeWidth = 1) const;

    void writeXML(std::string xmlFileName) const;

private:
    // This structure holds a bulk, its level,  its assigned process id and if it is an overlap
    // block or not. At each grid refinement level, the bulks are expressed in the coordinates of
    // the specific level.
    struct Block {
        Box3D bulk;
        plint level;
        plint processId;
        bool isOverlap;
    };

    std::vector<Box3D> computeBoxDifference(
        Box3D const &box, std::vector<Box3D> const &boxesToBeRemoved) const;

    std::vector<std::pair<plint, Box3D> > mergeBlocks(
        std::vector<std::pair<plint, Box3D> > const &blkToMerge) const;

    std::vector<std::pair<plint, Box3D> > getBlocksAndIdsAtLevelAndProcessorId(
        plint level, plint processId) const;

private:
    plint maxLevel;
    plint maxProcessId;
    std::map<plint, Box3D>
        boundingBoxes;  // Bounding box at each level (in the units of the same level).

    std::map<plint, Block> blocks;  // All blocks (normal and overlaps).
    std::map<plint, Array<plint, 26> >
        neighborBlockIds;  // Connectivity info for a block.
                           // Each block has 26 neighbors. If a neighbor is at the same
                           // or at a coarser grid level, then its id is stored. If there
                           // are many neighbors that belong to a finer grid level then the
                           // value OctreeTables::smaller() is stored. If there are no
                           // neighbors (domain boundary), then the value OctreeTables::border()
                           // is stored.
                           // Usually, we store the connectivity info for non-overlapping
                           // blocks, this is why we keep it in a separate map.
};

}  // namespace plb

#endif  // OCTREE_GRID_STRUCTURE_H
