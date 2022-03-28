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

#include "gridRefinement/octreeGridStructure.h"

#include "core/util.h"
#include "gridRefinement/octree.h"
#include "io/parallelIO.h"
#include "libraryInterfaces/TINYXML_xmlIO.h"
#include "libraryInterfaces/TINYXML_xmlIO.hh"
#include "multiBlock/sparseBlockStructure3D.h"
#include "multiBlock/threadAttribution.h"
// #include "interfaceGeneration/boxLogic3D.h"
// #include "interfaceGeneration/boxLogic3D.hh"

#include <algorithm>  // std::random_shuffle
#include <map>
#include <vector>

namespace plb {

OctreeGridStructure::OctreeGridStructure() : maxLevel(-1), maxProcessId(-1) { }

OctreeGridStructure::OctreeGridStructure(std::string xmlFileName) : maxLevel(-1), maxProcessId(-1)
{
    abortIfCannotOpenFileForReading(xmlFileName);
    XMLreader document(xmlFileName);

    XMLreaderProxy blocks = document["block"];
    for (; blocks.isValid(); blocks = blocks.iterId()) {
        plint blockId = blocks.getId();
        Array<plint, 6> domain;
        plint level;
        plint processId;
        blocks["bulk"].read<plint, 6>(domain);
        blocks["level"].read(level);
        blocks["processId"].read(processId);
        Box3D bulk(domain[0], domain[1], domain[2], domain[3], domain[4], domain[5]);
        bool isOverlap = false;  // By definition.
        addBlock(blockId, bulk, level, processId, isOverlap);
    }

    XMLreaderProxy overlaps = document["metadata"]["overlap"];
    for (; overlaps.isValid(); overlaps = overlaps.iterId()) {
        plint blockId = overlaps.getId();
        Array<plint, 6> domain;
        plint level;
        plint processId;
        overlaps["bulk"].read<plint, 6>(domain);
        overlaps["level"].read(level);
        overlaps["processId"].read(processId);
        Box3D bulk(domain[0], domain[1], domain[2], domain[3], domain[4], domain[5]);
        bool isOverlap = true;  // By definition.
        addBlock(blockId, bulk, level, processId, isOverlap);
    }

    XMLreaderProxy neighborBlockIds = document["metadata"]["neighborBlockIds"];
    for (; neighborBlockIds.isValid(); neighborBlockIds = neighborBlockIds.iterId()) {
        plint blockId = neighborBlockIds.getId();
        Array<plint, 26> ids;
        neighborBlockIds.read<plint, 26>(ids);
        addNeighborBlockIds(blockId, ids);
    }
}

void OctreeGridStructure::addBlock(
    plint blockId, Box3D const &bulk, plint level, plint processId, bool isOverlap)
{
    PLB_ASSERT(blockId >= 0 && level >= 0 && processId >= 0);

    maxLevel = std::max(maxLevel, level);
    maxProcessId = std::max(maxProcessId, processId);

    if (boundingBoxes.find(level) == boundingBoxes.end()) {
        boundingBoxes[level] = bulk;
    } else {
        boundingBoxes[level] = bound(boundingBoxes[level], bulk);
    }

    // PLB_ASSERT(blocks.find(blockId) == blocks.end());
    blocks[blockId].bulk = bulk;
    blocks[blockId].level = level;
    blocks[blockId].processId = processId;
    blocks[blockId].isOverlap = isOverlap;
}

void OctreeGridStructure::addNeighborBlockIds(plint blockId, Array<plint, 26> const &ids)
{
    PLB_ASSERT(blockId >= 0);
    // PLB_ASSERT(neighborBlockIds.find(blockId) == neighborBlockIds.end());
    neighborBlockIds[blockId] = ids;
}

bool OctreeGridStructure::removeBlock(plint blockId, bool removeConnectivityInfo)
{
    std::map<plint, Block>::const_iterator it = blocks.find(blockId);
    if (it == blocks.end()) {
        return (false);
    }

    plint levelOfRemovedBlock = it->second.level;

    plint numRemovedBlocks = blocks.erase(blockId);
    if (removeConnectivityInfo) {
        (void)removeNeighborBlockIds(blockId);
    }

    // If a block was indeed removed, then we need to update all the relevant information.
    if (numRemovedBlocks != 0) {
        maxLevel = -1;
        maxProcessId = -1;
        for (std::map<plint, Block>::const_iterator it = blocks.begin(); it != blocks.end(); ++it) {
            Block const &block = it->second;
            maxLevel = std::max(maxLevel, block.level);
            maxProcessId = std::max(maxProcessId, block.processId);
        }

#ifdef PLB_DEBUG
        plint numRemovedLevels = boundingBoxes.erase(levelOfRemovedBlock);
#endif
        PLB_ASSERT(numRemovedLevels == 1);
        for (std::map<plint, Block>::const_iterator it = blocks.begin(); it != blocks.end(); ++it) {
            Block const &block = it->second;
            if (block.level == levelOfRemovedBlock) {
                if (boundingBoxes.find(block.level) == boundingBoxes.end()) {
                    boundingBoxes[block.level] = block.bulk;
                } else {
                    boundingBoxes[block.level] = bound(boundingBoxes[block.level], block.bulk);
                }
            }
        }
    }

    return (numRemovedBlocks == 1);
}

bool OctreeGridStructure::removeNeighborBlockIds(plint blockId)
{
    plint numRemoved = neighborBlockIds.erase(blockId);
    return (numRemoved == 1);
}

plint OctreeGridStructure::getNumLevels() const
{
    return (maxLevel + 1);
}

plint OctreeGridStructure::getNumProcesses() const
{
    return (maxProcessId + 1);
}

Box3D OctreeGridStructure::getBoundingBox(plint level) const
{
    std::map<plint, Box3D>::const_iterator it = boundingBoxes.find(level);
    PLB_ASSERT(it != boundingBoxes.end());
    return (it->second);
}

Box3D OctreeGridStructure::getClosedCover(plint level) const
{
    PLB_ASSERT(level >= 0 && level <= maxLevel);

    plint levelCount = 0;
    Box3D openCoverAtMaxLevel;
    std::map<plint, Box3D>::const_iterator boundingBoxIt = boundingBoxes.begin();
    for (; boundingBoxIt != boundingBoxes.end(); ++boundingBoxIt) {
        plint boundingBoxLevel = boundingBoxIt->first;
        Box3D closedBoundingBoxAtLevel(boundingBoxIt->second);
        Box3D openBoundingBoxAtLevel(closedBoundingBoxAtLevel);
        openBoundingBoxAtLevel.x1++;
        openBoundingBoxAtLevel.y1++;
        openBoundingBoxAtLevel.z1++;
        Box3D openBoundingBoxAtMaxLevel(
            openBoundingBoxAtLevel.multiply(util::intTwoToThePower(maxLevel - boundingBoxLevel)));
        if (levelCount == 0) {
            openCoverAtMaxLevel = openBoundingBoxAtMaxLevel;
        } else {
            openCoverAtMaxLevel = bound(openCoverAtMaxLevel, openBoundingBoxAtMaxLevel);
        }
        levelCount++;
    }
    PLB_ASSERT(levelCount == (maxLevel + 1));

    Box3D openCoverAtLevel(openCoverAtMaxLevel.divide(util::intTwoToThePower(maxLevel - level)));
    Box3D closedCoverAtLevel(openCoverAtLevel);
    closedCoverAtLevel.x1--;
    closedCoverAtLevel.y1--;
    closedCoverAtLevel.z1--;

    return (closedCoverAtLevel);
}

bool OctreeGridStructure::neighborIsAtSameLevel(plint blockId, plint direction) const
{
    std::map<plint, Block>::const_iterator blockIt = blocks.find(blockId);
    std::map<plint, Array<plint, 26> >::const_iterator neighborBlockIdsIt =
        neighborBlockIds.find(blockId);

    PLB_ASSERT(blockIt != blocks.end());
    PLB_ASSERT(neighborBlockIdsIt != neighborBlockIds.end());
    PLB_ASSERT(direction >= 0 && direction < 26);

    plint neighborId = neighborBlockIdsIt->second[direction];
    if (neighborId >= 0) {
        std::map<plint, Block>::const_iterator neighborIt = blocks.find(neighborId);
        PLB_ASSERT(neighborIt != blocks.end());
        return (blockIt->second.level == neighborIt->second.level);
    }
    return (false);
}

bool OctreeGridStructure::neighborIsAtCoarserLevel(plint blockId, plint direction) const
{
    std::map<plint, Block>::const_iterator blockIt = blocks.find(blockId);
    std::map<plint, Array<plint, 26> >::const_iterator neighborBlockIdsIt =
        neighborBlockIds.find(blockId);

    PLB_ASSERT(blockIt != blocks.end());
    PLB_ASSERT(neighborBlockIdsIt != neighborBlockIds.end());
    PLB_ASSERT(direction >= 0 && direction < 26);

    plint neighborId = neighborBlockIdsIt->second[direction];
    if (neighborId >= 0) {
        std::map<plint, Block>::const_iterator neighborIt = blocks.find(neighborId);
        PLB_ASSERT(neighborIt != blocks.end());
        return (blockIt->second.level > neighborIt->second.level);
    }
    return (false);
}

bool OctreeGridStructure::neighborIsAtFinerLevel(plint blockId, plint direction) const
{
    typedef OctreeTables OT;

#ifdef PLB_DEBUG
    std::map<plint, Block>::const_iterator blockIt = blocks.find(blockId);
#endif
    std::map<plint, Array<plint, 26> >::const_iterator neighborBlockIdsIt =
        neighborBlockIds.find(blockId);

    PLB_ASSERT(blockIt != blocks.end());
    PLB_ASSERT(neighborBlockIdsIt != neighborBlockIds.end());
    PLB_ASSERT(direction >= 0 && direction < 26);

    plint neighborId = neighborBlockIdsIt->second[direction];
    return (neighborId == OT::smaller());
}

bool OctreeGridStructure::neighborIsBoundary(plint blockId, plint direction) const
{
    typedef OctreeTables OT;

#ifdef PLB_DEBUG
    std::map<plint, Block>::const_iterator blockIt = blocks.find(blockId);
#endif
    std::map<plint, Array<plint, 26> >::const_iterator neighborBlockIdsIt =
        neighborBlockIds.find(blockId);

    PLB_ASSERT(blockIt != blocks.end());
    PLB_ASSERT(neighborBlockIdsIt != neighborBlockIds.end());
    PLB_ASSERT(direction >= 0 && direction < 26);

    plint neighborId = neighborBlockIdsIt->second[direction];
    return (neighborId == OT::border());
}

bool OctreeGridStructure::neighborIsAllocated(plint blockId, plint direction) const
{
    typedef OctreeTables OT;

#ifdef PLB_DEBUG
    std::map<plint, Block>::const_iterator blockIt = blocks.find(blockId);
#endif
    std::map<plint, Array<plint, 26> >::const_iterator neighborBlockIdsIt =
        neighborBlockIds.find(blockId);

    PLB_ASSERT(blockIt != blocks.end());
    PLB_ASSERT(neighborBlockIdsIt != neighborBlockIds.end());
    PLB_ASSERT(direction >= 0 && direction < 26);

    plint neighborId = neighborBlockIdsIt->second[direction];
    return (!(neighborId == OT::undef()));
}

void OctreeGridStructure::getBlock(plint blockId, Box3D &bulk, plint &level, plint &processId) const
{
    std::map<plint, Block>::const_iterator it = blocks.find(blockId);
    PLB_ASSERT(it != blocks.end());
    Block const &block = it->second;
    bulk = block.bulk;
    level = block.level;
    processId = block.processId;
}

Array<plint, 26> OctreeGridStructure::getNeighborBlockIds(plint blockId) const
{
    std::map<plint, Array<plint, 26> >::const_iterator it = neighborBlockIds.find(blockId);
    PLB_ASSERT(it != neighborBlockIds.end());
    return (it->second);
}

bool OctreeGridStructure::isOverlap(plint blockId) const
{
    std::map<plint, Block>::const_iterator it = blocks.find(blockId);
    PLB_ASSERT(it != blocks.end());
    return (it->second.isOverlap);
}

std::vector<plint> OctreeGridStructure::getBlockIdsAtLevel(plint level, bool includeOverlaps) const
{
    PLB_ASSERT(level >= 0 && level <= maxLevel);

    std::vector<plint> blockIds;
    for (std::map<plint, Block>::const_iterator it = blocks.begin(); it != blocks.end(); ++it) {
        plint blockId = it->first;
        Block const &block = it->second;
        if (block.level == level && (includeOverlaps || !block.isOverlap)) {
            blockIds.push_back(blockId);
        }
    }

    return (blockIds);
}

std::vector<plint> OctreeGridStructure::getOverlapBlockIdsAtLevel(plint level) const
{
    PLB_ASSERT(level >= 0 && level <= maxLevel);

    std::vector<plint> blockIds;
    for (std::map<plint, Block>::const_iterator it = blocks.begin(); it != blocks.end(); ++it) {
        plint blockId = it->first;
        Block const &block = it->second;
        if (block.level == level && block.isOverlap) {
            blockIds.push_back(blockId);
        }
    }

    return (blockIds);
}

MultiBlockManagement3D OctreeGridStructure::getMultiBlockManagement(
    plint level, Box3D const &boundingBox, plint envelopeWidth) const
{
    PLB_ASSERT(level >= 0 && level <= maxLevel);
    PLB_ASSERT(envelopeWidth >= 0);

    SparseBlockStructure3D sparseBlock(boundingBox);
    ExplicitThreadAttribution threadAttribution;

    std::map<plint, Box3D>::const_iterator boundingBoxIt = boundingBoxes.find(level);
    PLB_ASSERT(boundingBoxIt != boundingBoxes.end());
    Box3D const &boundingBoxAtLevel = boundingBoxIt->second;

    if (doesIntersect(boundingBoxAtLevel, boundingBox)) {  // This is meant as an optimization.

        // for each process id ...
        for (plint iP = 0; iP < getNumProcesses(); ++iP) {
            // ... get all blocks at level "level"
            std::vector<std::pair<plint, Box3D> > minBlks =
                getBlocksAndIdsAtLevelAndProcessorId(level, iP);
            plint minNumBlks = (plint)minBlks.size();
            // pcout << "========Min Num Blocks = " << minNumBlks << std::endl;

            minBlks = mergeBlocks(minBlks);
            minNumBlks = (plint)minBlks.size();
            // pcout << "======= Old Num Blocks = " << minNumBlks << std::endl;

            do {
                minNumBlks = (plint)minBlks.size();
                minBlks = mergeBlocks(minBlks);
            } while (minNumBlks != (plint)minBlks.size());
            // minNumBlks = (plint)minBlks.size();
            // pcout << "======= Min Num Blocks = " << minNumBlks << std::endl;

            if (minNumBlks > 0) {
                for (std::vector<std::pair<plint, Box3D> >::iterator it = minBlks.begin();
                     it != minBlks.end(); ++it)
                {
                    // Although not all merged bulks may have intersections with the boundingBox,
                    // we retain their blockId (which does not correspond anymore to the blockId
                    // of the original grid structure file since blocks have been merged).
                    Box3D intersection;
                    if (intersect(it->second, boundingBox, intersection)) {
                        sparseBlock.addBlock(intersection, it->first);
                        threadAttribution.addBlock(it->first, iP);
                    }
                }
            }
        }
    }

    return (MultiBlockManagement3D(sparseBlock, threadAttribution.clone(), envelopeWidth, level));
}

MultiBlockManagement3D OctreeGridStructure::getMultiBlockManagement(
    plint level, plint envelopeWidth) const
{
    PLB_ASSERT(level >= 0 && level <= maxLevel);
    std::map<plint, Box3D>::const_iterator boundingBoxIt = boundingBoxes.find(level);
    PLB_ASSERT(boundingBoxIt != boundingBoxes.end());
    Box3D const &boundingBox = boundingBoxIt->second;
    return (getMultiBlockManagement(level, boundingBox, envelopeWidth));
}

std::vector<Box3D> OctreeGridStructure::computeBoxDifference(
    Box3D const &box, std::vector<Box3D> const &boxesToBeRemoved) const
{
    std::vector<Box3D> boxes;
    boxes.push_back(box);
    for (plint iBoxToBeRemoved = 0; iBoxToBeRemoved < (plint)boxesToBeRemoved.size();
         iBoxToBeRemoved++) {
        std::vector<Box3D> boxDifference;
        for (plint iBox = 0; iBox < (plint)boxes.size(); iBox++) {
            except(boxes[iBox], boxesToBeRemoved[iBoxToBeRemoved], boxDifference);
        }
        std::swap(boxes, boxDifference);
    }

    return (boxes);
}

// returns all the blocks that are located at a certain level located on a certain processor
std::vector<std::pair<plint, Box3D> > OctreeGridStructure::getBlocksAndIdsAtLevelAndProcessorId(
    plint level, plint processId) const
{
    std::vector<std::pair<plint, Box3D> > blkAndIdAtLvlAndProcId;
    for (std::map<plint, Block>::const_iterator it = blocks.begin(); it != blocks.end(); ++it) {
        plint blockId = it->first;
        Block const &block = it->second;
        if (block.level == level && block.processId == processId) {
            // but overlap blocks at the beginning. useful to merge them first
            if (!block.isOverlap) {
                blkAndIdAtLvlAndProcId.push_back(std::pair<plint, Box3D>(blockId, block.bulk));
            } else {
                blkAndIdAtLvlAndProcId.insert(
                    blkAndIdAtLvlAndProcId.begin(), std::pair<plint, Box3D>(blockId, block.bulk));
            }
        }
    }
    return blkAndIdAtLvlAndProcId;
}

// try to merge as many possible blocks in a vector of blocks.
// The blocks must be at the same level and at the same processorId.
// Comsumes the blkMerge.
std::vector<std::pair<plint, Box3D> > OctreeGridStructure::mergeBlocks(
    std::vector<std::pair<plint, Box3D> > const &blkToMerge) const
{
    std::vector<std::pair<plint, Box3D> > result = blkToMerge;
    std::sort(result.begin(), result.end());

    plint size = 0;
    while (size != (plint)result.size()) {
        size = (plint)result.size();
        for (plint iA = 0; iA < (plint)result.size(); ++iA) {
            for (plint iB = iA + 1; iB < (plint)result.size(); ++iB) {
                bool merged = merge(result[iA].second, result[iB].second);
                if (merged) {
                    result.erase(result.begin() + iB);
                    iB -= 1;
                }
            }
        }
    }

    return result;
}

MultiBlockManagement3D OctreeGridStructure::getMultiBlockManagementForOutput(
    plint level, Box3D const &boundingBox, bool crop, plint envelopeWidth) const
{
    typedef OctreeTables OT;

    PLB_ASSERT(level >= 0 && level <= maxLevel);
    PLB_ASSERT(envelopeWidth >= 0);

    SparseBlockStructure3D sparseBlock(boundingBox);
    ExplicitThreadAttribution threadAttribution;

    std::map<plint, Box3D>::const_iterator boundingBoxIt = boundingBoxes.find(level);
    PLB_ASSERT(boundingBoxIt != boundingBoxes.end());
    Box3D const &boundingBoxAtLevel = boundingBoxIt->second;

    if (doesIntersect(boundingBoxAtLevel, boundingBox)) {  // This is meant as an optimization.
        plint id = 0;
        for (std::map<plint, Block>::const_iterator it = blocks.begin(); it != blocks.end(); ++it) {
            plint blockId = it->first;
            Block const &block = it->second;
            if (block.level == level && !block.isOverlap) {
                plint xMargin = 0;
                plint yMargin = 0;
                plint zMargin = 0;
                plint h = 1;

                if (neighborIsAtCoarserLevel(blockId, OT::surface0P())
                    || (level != 0 && neighborIsBoundary(blockId, OT::surface0P()) && crop))
                {
                    xMargin = h;
                }
                if (neighborIsAtCoarserLevel(blockId, OT::surface1P())
                    || (level != 0 && neighborIsBoundary(blockId, OT::surface1P()) && crop))
                {
                    yMargin = h;
                }
                if (neighborIsAtCoarserLevel(blockId, OT::surface2P())
                    || (level != 0 && neighborIsBoundary(blockId, OT::surface2P()) && crop))
                {
                    zMargin = h;
                }

                Box3D bulk(block.bulk);
                std::vector<Box3D> domainsToBeRemoved;

                if (neighborIsAtCoarserLevel(blockId, OT::edge0PP())
                    || (level != 0 && neighborIsBoundary(blockId, OT::edge0PP()) && crop))
                {
                    domainsToBeRemoved.push_back(
                        Box3D(bulk.x0, bulk.x1, bulk.y1, bulk.y1, bulk.z1, bulk.z1));
                }
                if (neighborIsAtCoarserLevel(blockId, OT::edge1PP())
                    || (level != 0 && neighborIsBoundary(blockId, OT::edge1PP()) && crop))
                {
                    domainsToBeRemoved.push_back(
                        Box3D(bulk.x1, bulk.x1, bulk.y0, bulk.y1, bulk.z1, bulk.z1));
                }
                if (neighborIsAtCoarserLevel(blockId, OT::edge2PP())
                    || (level != 0 && neighborIsBoundary(blockId, OT::edge2PP()) && crop))
                {
                    domainsToBeRemoved.push_back(
                        Box3D(bulk.x1, bulk.x1, bulk.y1, bulk.y1, bulk.z0, bulk.z1));
                }

                if ((neighborIsAtCoarserLevel(blockId, OT::cornerPPP())
                     || (level != 0 && neighborIsBoundary(blockId, OT::cornerPPP())))
                    && crop)
                {
                    domainsToBeRemoved.push_back(
                        Box3D(bulk.x1, bulk.x1, bulk.y1, bulk.y1, bulk.z1, bulk.z1));
                }

                bulk.x1 -= xMargin;
                bulk.y1 -= yMargin;
                bulk.z1 -= zMargin;
                std::vector<Box3D> domains = computeBoxDifference(bulk, domainsToBeRemoved);

                for (plint iDomain = 0; iDomain < (plint)domains.size(); iDomain++) {
                    Box3D intersection;
                    if (intersect(domains[iDomain], boundingBox, intersection)) {
                        sparseBlock.addBlock(intersection, id);
                        threadAttribution.addBlock(id, block.processId);
                        id++;
                    }
                }
            }
        }
    }

    return (MultiBlockManagement3D(sparseBlock, threadAttribution.clone(), envelopeWidth, level));
}

MultiBlockManagement3D OctreeGridStructure::getMultiBlockManagementForOutput(
    plint level, bool crop, plint envelopeWidth) const
{
    PLB_ASSERT(level >= 0 && level <= maxLevel);
    std::map<plint, Box3D>::const_iterator boundingBoxIt = boundingBoxes.find(level);
    PLB_ASSERT(boundingBoxIt != boundingBoxes.end());
    Box3D const &boundingBox = boundingBoxIt->second;
    return (getMultiBlockManagementForOutput(level, boundingBox, crop, envelopeWidth));
}

void OctreeGridStructure::writeXML(std::string xmlFileName) const
{
    XMLwriter document;

    for (std::map<plint, Block>::const_iterator it = blocks.begin(); it != blocks.end(); ++it) {
        plint blockId = it->first;
        Block const &block = it->second;
        if (block.isOverlap) {
            continue;
        }
        Box3D bulk(block.bulk);
        plint level(block.level);
        plint processId(block.processId);
        Array<plint, 6> domain(bulk.x0, bulk.x1, bulk.y0, bulk.y1, bulk.z0, bulk.z1);
        document["block"][blockId]["bulk"].set<plint, 6>(domain);
        document["block"][blockId]["level"].set(level);
        document["block"][blockId]["processId"].set(processId);
    }

    for (std::map<plint, Block>::const_iterator it = blocks.begin(); it != blocks.end(); ++it) {
        plint blockId = it->first;
        Block const &block = it->second;
        if (!block.isOverlap) {
            continue;
        }
        Box3D bulk(block.bulk);
        plint level(block.level);
        plint processId(block.processId);
        Array<plint, 6> domain(bulk.x0, bulk.x1, bulk.y0, bulk.y1, bulk.z0, bulk.z1);
        document["metadata"]["overlap"][blockId]["bulk"].set<plint, 6>(domain);
        document["metadata"]["overlap"][blockId]["level"].set(level);
        document["metadata"]["overlap"][blockId]["processId"].set(processId);
    }

    for (std::map<plint, Array<plint, 26> >::const_iterator it = neighborBlockIds.begin();
         it != neighborBlockIds.end(); ++it)
    {
        plint blockId = it->first;
        Array<plint, 26> const &neighborIds = it->second;
        document["metadata"]["neighborBlockIds"][blockId].set<plint, 26>(neighborIds);
    }

    document.print(xmlFileName);
}

}  // namespace plb
