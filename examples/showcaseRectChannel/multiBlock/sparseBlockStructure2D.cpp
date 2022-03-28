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
 * Fast algorithm for handling sparse multi-block structure -- implementation.
 */

#include "multiBlock/sparseBlockStructure2D.h"

#include <set>

#include "core/plbDebug.h"
#include "multiBlock/defaultMultiBlockPolicy2D.h"

namespace plb {

SparseBlockStructure2D::SparseBlockStructure2D(plint nx, plint ny) :
    boundingBox(Box2D(0, nx - 1, 0, ny - 1))
{
    PLB_PRECONDITION(nx >= 1 && ny >= 1);
    defaultGridN();
    iniGridParameters();
}

SparseBlockStructure2D::SparseBlockStructure2D(Box2D boundingBox_) : boundingBox(boundingBox_)
{
    PLB_PRECONDITION(boundingBox.getNx() >= 1 && boundingBox.getNy() >= 1);
    defaultGridN();
    iniGridParameters();
}

SparseBlockStructure2D::SparseBlockStructure2D(Box2D boundingBox_, plint gridNx_, plint gridNy_) :
    boundingBox(boundingBox_), gridNx(gridNx_), gridNy(gridNy_)
{
    PLB_PRECONDITION(boundingBox.getNx() >= 1 && boundingBox.getNy() >= 1);
    iniGridParameters();
}

SparseBlockStructure2D::SparseBlockStructure2D(
    SparseBlockStructure2D const &rhs, Box2D boundingBox_) :
    boundingBox(boundingBox_)
{
    PLB_PRECONDITION(boundingBox.getNx() >= 1 && boundingBox.getNy() >= 1);
    gridNx =
        (plint)(0.5 + (double)rhs.gridNx * (double)boundingBox.getNx() / (double)rhs.boundingBox.getNx());
    if (gridNx < 1)
        gridNx = 1;
    gridNy =
        (plint)(0.5 + (double)rhs.gridNy * (double)boundingBox.getNy() / (double)rhs.boundingBox.getNy());
    if (gridNy < 1)
        gridNy = 1;
}

void SparseBlockStructure2D::addBlock(Box2D const &bulk, plint blockId)
{
    addBlock(bulk, bulk, blockId);
}

void SparseBlockStructure2D::addBlock(Box2D const &bulk, Box2D const &uniqueBulk, plint blockId)
{
    PLB_ASSERT(contained(uniqueBulk, bulk));
    // To avoid inconsistent state if the same blockId is added twice,
    //   remove a potentially existing instance of blockId.
#ifdef PLB_DEBUG
    PLB_ASSERT(!exists(blockId));
#else
    removeBlock(blockId);
#endif
    bulks[blockId] = bulk;
    uniqueBulks[blockId] = uniqueBulk;
    integrateBlock(blockId, bulk);
}

void SparseBlockStructure2D::removeBlock(plint blockId)
{
    std::map<plint, Box2D>::iterator it = bulks.find(blockId);
    if (it != bulks.end()) {
        extractBlock(blockId);
        bulks.erase(it);
        uniqueBulks.erase(blockId);
    }
}

bool SparseBlockStructure2D::exists(plint blockId)
{
    return bulks.find(blockId) != bulks.end();
}

plint SparseBlockStructure2D::nextIncrementalId() const
{
    if (bulks.empty()) {
        return 0;
    } else {
        return (--bulks.end())->first + 1;
    }
}

Box2D SparseBlockStructure2D::getBoundingBox() const
{
    return boundingBox;
}

bool SparseBlockStructure2D::getBulk(plint blockId, Box2D &bulk) const
{
    std::map<plint, Box2D>::const_iterator it = bulks.find(blockId);
    if (it != bulks.end()) {
        bulk = it->second;
        return true;
    }
    return false;
}

bool SparseBlockStructure2D::getUniqueBulk(plint blockId, Box2D &uniqueBulk) const
{
    std::map<plint, Box2D>::const_iterator it = uniqueBulks.find(blockId);
    if (it != uniqueBulks.end()) {
        uniqueBulk = it->second;
        return true;
    }
    return false;
}

plint SparseBlockStructure2D::getNumBlocks() const
{
    return (plint)bulks.size();
}

plint SparseBlockStructure2D::getNumBulkCells() const
{
    plint nCells = 0;
    std::map<plint, Box2D>::const_iterator it = bulks.begin();
    for (; it != bulks.end(); ++it) {
        nCells += it->second.nCells();
    }
    return nCells;
}

std::map<plint, Box2D> const &SparseBlockStructure2D::getBulks() const
{
    return bulks;
}

std::vector<plint> SparseBlockStructure2D::getLocalBlocks(
    ThreadAttribution const &attribution) const
{
    std::vector<plint> localBlocks;
    std::map<plint, Box2D>::const_iterator it = bulks.begin();
    for (; it != bulks.end(); ++it) {
        plint blockId = it->first;
        if (attribution.isLocal(blockId)) {
            localBlocks.push_back(blockId);
        }
    }
    return localBlocks;
}

plint SparseBlockStructure2D::locate(plint iX, plint iY) const
{
    GridT::const_iterator gridIter = grid.find(Dot2D(gridPosX(iX), gridPosY(iY)));
    if (gridIter != grid.end()) {
        std::vector<plint> const &blockList = gridIter->second;
        for (pluint iBlock = 0; iBlock < blockList.size(); ++iBlock) {
            Box2D const &bulk = bulks.find(blockList[iBlock])->second;
            if (contained(iX, iY, bulk)) {
                return blockList[iBlock];
            }
        }
    }
    return -1;
}

void SparseBlockStructure2D::defaultGridN()
{
    double uniformGridN = std::pow((double)defaultMultiBlockPolicy2D().getNumGridPoints(), 1. / 2.);
    double uniformNcell = std::pow((double)(boundingBox.nCells()), 1. / 2.);
    gridNx = (plint)(0.5 + (double)boundingBox.getNx() / uniformNcell * uniformGridN);
    if (gridNx < 1)
        gridNx = 1;
    gridNy = (plint)(0.5 + (double)boundingBox.getNy() / uniformNcell * uniformGridN);
    if (gridNy < 1)
        gridNy = 1;
}

void SparseBlockStructure2D::iniGridParameters()
{
    gridLx = boundingBox.getNx() / gridNx;
    if (boundingBox.getNx() % gridNx != 0) {
        ++gridLx;
    }
    gridLy = boundingBox.getNy() / gridNy;
    if (boundingBox.getNy() % gridNy != 0) {
        ++gridLy;
    }
}

plint SparseBlockStructure2D::gridPosX(plint realX) const
{
    return (realX - boundingBox.x0) / gridLx;
}

plint SparseBlockStructure2D::gridPosY(plint realY) const
{
    return (realY - boundingBox.y0) / gridLy;
}

Box2D SparseBlockStructure2D::getGridBox(Box2D const &realBlock) const
{
    return Box2D(
        gridPosX(realBlock.x0), gridPosX(realBlock.x1), gridPosY(realBlock.y0),
        gridPosY(realBlock.y1));
}

void SparseBlockStructure2D::intersect(
    Box2D const &bulk, std::vector<plint> &ids, std::vector<Box2D> &intersections) const
{
    Box2D gridBox = getGridBox(bulk);
    Box2D intersection;  // Temporary variable.

    std::set<plint> idsToTest;
    for (plint gridX = gridBox.x0; gridX <= gridBox.x1; ++gridX) {
        for (plint gridY = gridBox.y0; gridY <= gridBox.y1; ++gridY) {
            GridT::const_iterator gridIter = grid.find(Dot2D(gridX, gridY));
            if (gridIter != grid.end()) {
                std::vector<plint> const &blockList = gridIter->second;
                idsToTest.insert(blockList.begin(), blockList.end());
            }
        }
    }

    std::set<plint>::const_iterator it = idsToTest.begin();
    for (; it != idsToTest.end(); ++it) {
        Box2D testBlock = bulks.find(*it)->second;
        if (plb::intersect(bulk, testBlock, intersection)) {
            intersections.push_back(intersection);
            ids.push_back(*it);
        }
    }
}

void SparseBlockStructure2D::computeEnvelopeTerm(
    plint block0, plint block1, plint &env0, plint &env1, plint delta) const
{
    if (delta < 0) {
        env0 = block0 + delta;
        env1 = block0 - 1;
    } else if (delta > 0) {
        env0 = block1 + 1;
        env1 = block1 + delta;
    } else {
        env0 = block0;
        env1 = block1;
    }
}

void SparseBlockStructure2D::computeOverlaps(
    Box2D const &bulk, plint dx, plint dy, std::vector<plint> &ids,
    std::vector<Box2D> &overlapsOnBulk, std::vector<Box2D> &overlapsOnNeighbors) const
{
    Box2D envelope;
    std::vector<Box2D> intersections;
    computeEnvelopeTerm(bulk.x0, bulk.x1, envelope.x0, envelope.x1, dx);
    computeEnvelopeTerm(bulk.y0, bulk.y1, envelope.y0, envelope.y1, dy);
    intersect(envelope, ids, intersections);
    overlapsOnNeighbors.insert(
        overlapsOnNeighbors.end(), intersections.begin(), intersections.end());
    for (pluint iOverlap = 0; iOverlap < intersections.size(); ++iOverlap) {
        overlapsOnBulk.push_back(intersections[iOverlap].shift(-dx, -dy));
    }
}

void SparseBlockStructure2D::computeOverlaps(
    plint blockId, plint envelopeWidth, std::vector<plint> &ids, std::vector<Box2D> &overlapsOnBulk,
    std::vector<Box2D> &overlapsOnNeighbors) const
{
    Box2D bulk = bulks.find(blockId)->second;
    for (plint dx = -1; dx <= +1; ++dx) {
        for (plint dy = -1; dy <= +1; ++dy) {
            if (!(dx == 0 && dy == 0)) {
                computeOverlaps(
                    bulk, dx * envelopeWidth, dy * envelopeWidth, ids, overlapsOnBulk,
                    overlapsOnNeighbors);
            }
        }
    }
}

void SparseBlockStructure2D::findNeighbors(
    Box2D const &bulk, plint neighborhoodWidth, std::vector<plint> &neighbors,
    plint excludeId) const
{
    Box2D extendedBlock(bulk.enlarge(neighborhoodWidth));
    Box2D gridBox = getGridBox(extendedBlock);

    std::set<plint> idsToTest;
    for (plint gridX = gridBox.x0; gridX <= gridBox.x1; ++gridX) {
        for (plint gridY = gridBox.y0; gridY <= gridBox.y1; ++gridY) {
            // Test only on the surface of the grid-box.
            if (gridX == gridBox.x0 || gridX == gridBox.x1 || gridY == gridBox.y0
                || gridY == gridBox.y1) {
                GridT::const_iterator gridIter = grid.find(Dot2D(gridX, gridY));
                if (gridIter != grid.end()) {
                    std::vector<plint> const &blockList = gridIter->second;
                    idsToTest.insert(blockList.begin(), blockList.end());
                }
            }
        }
    }

    if (excludeId >= 0) {
        std::set<plint>::iterator excludeIter = idsToTest.find(excludeId);
        if (excludeIter != idsToTest.end()) {
            idsToTest.erase(excludeIter);
        }
    }

    std::set<plint>::const_iterator it = idsToTest.begin();
    for (; it != idsToTest.end(); ++it) {
        Box2D testBlock = bulks.find(*it)->second;
        if (plb::doesIntersect(extendedBlock, testBlock)) {
            neighbors.push_back(*it);
        }
    }
}

void SparseBlockStructure2D::findNeighbors(
    plint blockId, plint neighborhoodWidth, std::vector<plint> &neighbors) const
{
    Box2D bulk = bulks.find(blockId)->second;
    findNeighbors(bulk, neighborhoodWidth, neighbors, blockId);
}

void SparseBlockStructure2D::swap(SparseBlockStructure2D &rhs)
{
    std::swap(boundingBox, rhs.boundingBox);
    std::swap(gridLx, rhs.gridLx);
    std::swap(gridLy, rhs.gridLy);
    std::swap(gridNx, rhs.gridNx);
    std::swap(gridNy, rhs.gridNy);
    grid.swap(rhs.grid);
    bulks.swap(rhs.bulks);
    uniqueBulks.swap(rhs.uniqueBulks);
}

bool SparseBlockStructure2D::equals(SparseBlockStructure2D const &rhs) const
{
    return boundingBox == rhs.boundingBox && gridNx == rhs.gridNx && gridNy == rhs.gridNy
           && grid == rhs.grid && bulks == rhs.bulks;
}

void SparseBlockStructure2D::integrateBlock(plint blockId, Box2D bulk)
{
    Box2D gridBox = getGridBox(bulk);

    for (plint gridX = gridBox.x0; gridX <= gridBox.x1; ++gridX) {
        for (plint gridY = gridBox.y0; gridY <= gridBox.y1; ++gridY) {
            grid[Dot2D(gridX, gridY)].push_back(blockId);
        }
    }
}

void SparseBlockStructure2D::extractBlock(plint blockId)
{
    Box2D const &bulk = bulks[blockId];
    Box2D gridBox = getGridBox(bulk);

    for (plint gridX = gridBox.x0; gridX <= gridBox.x1; ++gridX) {
        for (plint gridY = gridBox.y0; gridY <= gridBox.y1; ++gridY) {
            std::vector<plint> &blockList = grid[Dot2D(gridX, gridY)];
            // Use remove-erase idiom (because blockList is a std::vector).
            blockList.erase(
                std::remove(blockList.begin(), blockList.end(), blockId), blockList.end());
        }
    }
}

SparseBlockStructure2D intersect(SparseBlockStructure2D const &sparseBlock, Box2D domain, bool crop)
{
    Box2D newBoundingBox;
    if (crop) {
        if (!intersect(domain, sparseBlock.getBoundingBox(), newBoundingBox)) {
            // Refuse to crop the domain if the domains don't intersect,
            //   to avoid creating a block in undefined state.
            newBoundingBox = sparseBlock.getBoundingBox();
        }
    } else {
        newBoundingBox = sparseBlock.getBoundingBox();
    }
    return intersect(sparseBlock, domain, newBoundingBox);
}

SparseBlockStructure2D intersect(
    SparseBlockStructure2D const &sparseBlock, Box2D domain, Box2D newBoundingBox)
{
    SparseBlockStructure2D newSparseBlock(newBoundingBox);
    std::vector<plint> ids;
    std::vector<Box2D> intersections;
    sparseBlock.intersect(domain, ids, intersections);
    for (pluint iBox = 0; iBox < ids.size(); ++iBox) {
        Box2D uniqueBulk, intersectedUniqueBulk;
#ifdef PLB_DEBUG
        bool ok =
#endif
            sparseBlock.getUniqueBulk(ids[iBox], uniqueBulk);
        PLB_ASSERT(ok);
#ifdef PLB_DEBUG
        bool doesIntersect =
#endif
            intersect(domain, uniqueBulk, intersectedUniqueBulk);
        PLB_ASSERT(doesIntersect);
        newSparseBlock.addBlock(intersections[iBox], intersectedUniqueBulk, ids[iBox]);
    }
    return newSparseBlock;
}

/** The block-ids are inherited from sparseBlocks1.
 **/
SparseBlockStructure2D intersect(
    SparseBlockStructure2D const &sparseBlock1, SparseBlockStructure2D const &sparseBlock2,
    bool crop)
{
    Box2D newBoundingBox;
    if (crop) {
        if (!intersect(
                sparseBlock1.getBoundingBox(), sparseBlock2.getBoundingBox(), newBoundingBox)) {
            // Avoid undefined state if there is no intersection.
            newBoundingBox = bound(sparseBlock1.getBoundingBox(), sparseBlock2.getBoundingBox());
        }
    } else {
        newBoundingBox = bound(sparseBlock1.getBoundingBox(), sparseBlock2.getBoundingBox());
    }
    SparseBlockStructure2D newSparseBlock(newBoundingBox);

    std::vector<plint> ids1;           // temporary variable.
    std::vector<Box2D> intersections;  // temporary variable.

    std::map<plint, Box2D> const &bulks2 = sparseBlock2.getBulks();
    std::map<plint, Box2D>::const_iterator it2 = bulks2.begin();
    for (; it2 != bulks2.end(); ++it2) {
        ids1.clear();
        intersections.clear();
        sparseBlock1.intersect(it2->second, ids1, intersections);
        for (pluint iBox = 0; iBox < ids1.size(); ++iBox) {
            Box2D uniqueBulk, intersectedUniqueBulk;
            sparseBlock1.getUniqueBulk(ids1[iBox], uniqueBulk);
#ifdef PLB_DEBUG
            bool doesIntersect =
#endif
                intersect(it2->second, uniqueBulk, intersectedUniqueBulk);
            PLB_ASSERT(doesIntersect);
            newSparseBlock.addBlock(intersections[iBox], intersectedUniqueBulk, ids1[iBox]);
        }
    }
    return newSparseBlock;
}

SparseBlockStructure2D extend(
    SparseBlockStructure2D const &sparseBlock, Box2D addedBulk, Box2D addedUniqueBulk,
    std::vector<plint> &newIds)
{
    Box2D newBoundingBox = bound(sparseBlock.getBoundingBox(), addedBulk);
    SparseBlockStructure2D newSparseBlock(newBoundingBox);

    // Simply take over all blocks from the original sparseBlock.
    std::map<plint, Box2D> const &bulks = sparseBlock.getBulks();
    std::map<plint, Box2D>::const_iterator it = bulks.begin();
    for (; it != bulks.end(); ++it) {
        Box2D uniqueBulk;
        sparseBlock.getUniqueBulk(it->first, uniqueBulk);
        newSparseBlock.addBlock(it->second, uniqueBulk, it->first);
    }

    // Remove from the added bulk all domains which are already in
    //   the sparseBlock.
    std::vector<Box2D> intersections;
    std::vector<plint> ids;
    sparseBlock.intersect(addedBulk, ids, intersections);
    std::vector<Box2D> newDomains;
    std::vector<Box2D> tmpNewDomains;
    newDomains.push_back(addedBulk);
    for (pluint iInters = 0; iInters < intersections.size(); ++iInters) {
        tmpNewDomains.clear();
        for (pluint iNew = 0; iNew < newDomains.size(); ++iNew) {
            except(newDomains[iNew], intersections[iInters], tmpNewDomains);
        }
        tmpNewDomains.swap(newDomains);
    }

    // Add the computed domains to the new sparseBlock.
    for (pluint iNew = 0; iNew < newDomains.size(); ++iNew) {
        plint nextId = newSparseBlock.nextIncrementalId();
        Box2D newBulk = newDomains[iNew];
        // The unique bulks are evaluated as the intersection between
        //   the computed domains and the user-provided unique bulk.
        Box2D newUniqueBulk;
#ifdef PLB_DEBUG
        bool doesIntersect =
#endif
            intersect(newBulk, addedUniqueBulk, newUniqueBulk);
        PLB_ASSERT(doesIntersect);

        newSparseBlock.addBlock(newBulk, newUniqueBulk, nextId);
        newIds.push_back(nextId);
    }
    return newSparseBlock;
}

/** The variable remappedIds maps blockIds from the original sparse block to blockIds in
 *  the newly created structure, for all blocks in the new structure which were created
 *  by exeption from the original sparseBlock.
 **/
SparseBlockStructure2D except(
    SparseBlockStructure2D const &sparseBlock, Box2D exceptedBlock,
    std::map<plint, std::vector<plint> > &remappedIds)
{
    // Simply take over all blocks and the bounding-box from the original sparseBlock.
    SparseBlockStructure2D newSparseBlock(sparseBlock);

    // Find all blocks that intersect with "exceptedBlock"...
    std::vector<Box2D> intersections;
    std::vector<plint> ids;
    sparseBlock.intersect(exceptedBlock, ids, intersections);
    std::vector<Box2D> fragments;
    for (pluint iBlock = 0; iBlock < intersections.size(); ++iBlock) {
        Box2D bulk, uniqueBulk;
        newSparseBlock.getBulk(ids[iBlock], bulk);
        newSparseBlock.getUniqueBulk(ids[iBlock], uniqueBulk);
        // ... remove them ...
        newSparseBlock.removeBlock(ids[iBlock]);

        fragments.clear();
        except(bulk, intersections[iBlock], fragments);
        std::vector<plint> newIds;
        // ... and put back evertything which is not in the intersection.
        for (pluint iFragment = 0; iFragment < fragments.size(); ++iFragment) {
            plint newId = newSparseBlock.nextIncrementalId();
            newIds.push_back(newId);
            Box2D fragmentUniqueBulk;
#ifdef PLB_DEBUG
            bool doesIntersect =
#endif
                intersect(uniqueBulk, fragments[iFragment], fragmentUniqueBulk);
            PLB_ASSERT(doesIntersect);
            newSparseBlock.addBlock(fragments[iFragment], fragmentUniqueBulk, newId);
        }
        remappedIds[ids[iBlock]] = newIds;
    }
    return newSparseBlock;
}

/** The bounding box of the created structure contains the two original bounding boxes.
 *  All blocks contained in sparseBlock1 are also contained in the created structure,
 *  with same blockIds. The blocks of sparseBlock2 are fractioned in order not to in-
 *  tersect with sparseBlock1's blocks, and the blockIds are remapped. The variable
 *  remappedIds maps blockIds from sparseBlock2 to blockIds in the newly created
 *  structure, for all blocks in the new structure which were created from sparseBlock2.
 **/
SparseBlockStructure2D block_union(
    SparseBlockStructure2D const &sparseBlock1, SparseBlockStructure2D const &sparseBlock2,
    std::map<plint, std::vector<plint> > &remappedIds)
{
    Box2D newBoundingBox = bound(sparseBlock1.getBoundingBox(), sparseBlock2.getBoundingBox());
    SparseBlockStructure2D newSparseBlock(newBoundingBox);

    // Simply take over all blocks from sparseBlock1.
    std::map<plint, Box2D> const &bulks1 = sparseBlock1.getBulks();
    std::map<plint, Box2D>::const_iterator it1 = bulks1.begin();
    for (; it1 != bulks1.end(); ++it1) {
        Box2D uniqueBulk1;
        sparseBlock1.getUniqueBulk(it1->first, uniqueBulk1);
        newSparseBlock.addBlock(it1->second, uniqueBulk1, it1->first);
    }

    // From sparseBlock2, take over everything with exception of what's already
    //   in newSparseBlock. This is not so easy, of course.
    std::map<plint, Box2D> const &bulks2 = sparseBlock2.getBulks();
    std::map<plint, Box2D>::const_iterator it2 = bulks2.begin();
    for (; it2 != bulks2.end(); ++it2) {
        Box2D uniqueBulk2;
        sparseBlock2.getUniqueBulk(it2->first, uniqueBulk2);
        std::vector<plint> ids;
        std::vector<Box2D> intersections;
        sparseBlock1.intersect(it2->second, ids, intersections);
        // The simplest case: if there's no intersection, the new block can
        //   simply be added.
        if (ids.empty()) {
            plint nextId = newSparseBlock.nextIncrementalId();
            newSparseBlock.addBlock(it2->second, uniqueBulk2, nextId);
            // Keep track of ID remapping as we go from sparseBlock2 to newSparseBlock.
            std::vector<plint> remappedId;
            remappedId.push_back(nextId);
            remappedIds[it2->first] = remappedId;
        } else {
            std::vector<Box2D> originalBlocks;
            std::vector<Box2D> tmp;
            originalBlocks.push_back(it2->second);
            for (pluint iInters = 0; iInters < intersections.size(); ++iInters) {
                tmp.clear();
                for (pluint iOrigin = 0; iOrigin < originalBlocks.size(); ++iOrigin) {
                    except(originalBlocks[iOrigin], intersections[iInters], tmp);
                }
                originalBlocks.swap(tmp);
            }
            std::vector<plint> nextIds;
            for (pluint iOrigin = 0; iOrigin < originalBlocks.size(); ++iOrigin) {
                Box2D intersectedUniqueBulk2;
#ifdef PLB_DEBUG
                bool doesIntersect =
#endif
                    intersect(uniqueBulk2, originalBlocks[iOrigin], intersectedUniqueBulk2);
                PLB_ASSERT(doesIntersect);
                plint nextId = newSparseBlock.nextIncrementalId();
                newSparseBlock.addBlock(originalBlocks[iOrigin], intersectedUniqueBulk2, nextId);
                nextIds.push_back(nextId);
            }
            // Keep track of ID remapping as we go from sparseBlock2 to newSparseBlock.
            remappedIds[it2->first] = nextIds;
        }
    }
    return newSparseBlock;
}

SparseBlockStructure2D alignDistribution2D(
    SparseBlockStructure2D const &originalStructure, SparseBlockStructure2D const &partnerStructure,
    std::vector<plint> &newIds, std::map<plint, std::vector<plint> > &remappedFromPartner)
{
    // Holds the return value.
    SparseBlockStructure2D newSparseBlock(originalStructure);

    std::vector<plint> ids;                                                    // Temporary.
    std::vector<Box2D> intersections, exceptedBlocks, intersectedUniqueBulks;  // Temporary.

    std::vector<plint> overlapRegionIds;
    std::vector<std::vector<Box2D> > overlapRegionBulks;
    std::vector<std::vector<Box2D> > overlapRegionUniqueBulks;

    // Iterate over blocks of partner, in order to remove corresponding areas
    //   in the original structure, and then add them with appropriate new IDs.
    std::map<plint, Box2D> const &partnerBulks = partnerStructure.getBulks();
    std::map<plint, Box2D>::const_iterator partnerIt = partnerBulks.begin();
    for (; partnerIt != partnerBulks.end(); ++partnerIt) {
        Box2D partnerBulk = partnerIt->second;
        ids.clear();
        intersections.clear();
        newSparseBlock.intersect(partnerBulk, ids, intersections);
        intersectedUniqueBulks.clear();
        // Remove areas of overlap (they will be added at the end of this function),
        //   but keep pieces of the original domain which don't overlap.
        for (pluint iInters = 0; iInters < intersections.size(); ++iInters) {
            Box2D intersection = intersections[iInters];
            plint originalId = ids[iInters];
            Box2D originalBulk, originalUniqueBulk;
            newSparseBlock.getBulk(originalId, originalBulk);
            newSparseBlock.getUniqueBulk(originalId, originalUniqueBulk);

            Box2D intersectedUniqueBulk;
#ifdef PLB_DEBUG
            bool doesIntersect =
#endif
                intersect(intersection, originalUniqueBulk, intersectedUniqueBulk);
            PLB_ASSERT(doesIntersect);
            intersectedUniqueBulks.push_back(intersectedUniqueBulk);

            // Remove blocks from the original structure which overlap with partner,
            //   fraction them, and add back components which don't overlap.
            newSparseBlock.removeBlock(originalId);
            exceptedBlocks.clear();
            except(originalBulk, intersection, exceptedBlocks);
            for (pluint iExc = 0; iExc < exceptedBlocks.size(); ++iExc) {
                Box2D exceptedBlock = exceptedBlocks[iExc];
                Box2D exceptedUniqueBulk;
#ifdef PLB_DEBUG
                bool doesIntersect =
#endif
                    intersect(exceptedBlock, originalUniqueBulk, exceptedUniqueBulk);
                PLB_ASSERT(doesIntersect);
                newSparseBlock.addBlock(
                    exceptedBlock, exceptedUniqueBulk, newSparseBlock.nextIncrementalId());
            }
        }
        overlapRegionIds.push_back(partnerIt->first);
        overlapRegionBulks.push_back(intersections);
        overlapRegionUniqueBulks.push_back(intersectedUniqueBulks);
    }

    // Keep track of IDs for left-overs from originalBlock which don't overlap
    //   with partnerBlock.
    std::map<plint, Box2D> const &originalBulks = newSparseBlock.getBulks();
    std::map<plint, Box2D>::const_iterator originalIt = originalBulks.begin();
    for (; originalIt != originalBulks.end(); ++originalIt) {
        newIds.push_back(originalIt->first);
    }

    std::vector<plint> remappedIds;
    for (pluint iOverlap = 0; iOverlap < overlapRegionIds.size(); ++iOverlap) {
        plint oldId = overlapRegionIds[iOverlap];
        remappedIds.clear();
        for (pluint iComp = 0; iComp < overlapRegionBulks[iOverlap].size(); ++iComp) {
            plint nextId = newSparseBlock.nextIncrementalId();
            newSparseBlock.addBlock(
                overlapRegionBulks[iOverlap][iComp], overlapRegionUniqueBulks[iOverlap][iComp],
                nextId);
            remappedIds.push_back(nextId);
        }
        remappedFromPartner[oldId] = remappedIds;
    }
    return newSparseBlock;
}

EuclideanIterator2D::EuclideanIterator2D(SparseBlockStructure2D const &sparseBlock_) :
    sparseBlock(sparseBlock_)
{ }

bool EuclideanIterator2D::getNextChunkX(plint iX, plint iY, plint &blockId, plint &chunkSize) const
{
    blockId = sparseBlock.locate(iX, iY);
    if (blockId == -1) {
        Box2D boundingBox = sparseBlock.getBoundingBox();
        plint exploreX = iX + 1;
        while (exploreX < boundingBox.getNx() && sparseBlock.locate(exploreX, iY) == -1) {
            ++exploreX;
        }
        chunkSize = exploreX - iX;
        return false;
    } else {
        Box2D bulk;
        sparseBlock.getBulk(blockId, bulk);
        chunkSize = bulk.x1 - iX + 1;
        return true;
    }
}

bool EuclideanIterator2D::getNextChunkY(plint iX, plint iY, plint &blockId, plint &chunkSize) const
{
    blockId = sparseBlock.locate(iX, iY);
    if (blockId == -1) {
        Box2D boundingBox = sparseBlock.getBoundingBox();
        plint exploreY = iY + 1;
        while (exploreY < boundingBox.getNy() && sparseBlock.locate(iX, exploreY) == -1) {
            ++exploreY;
        }
        chunkSize = exploreY - iY;
        return false;
    } else {
        Box2D bulk;
        sparseBlock.getBulk(blockId, bulk);
        chunkSize = bulk.y1 - iY + 1;
        return true;
    }
}

}  // namespace plb
