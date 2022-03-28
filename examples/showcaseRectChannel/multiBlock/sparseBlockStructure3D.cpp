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

#include "multiBlock/sparseBlockStructure3D.h"

#include <set>

#include "multiBlock/defaultMultiBlockPolicy3D.h"

namespace plb {

SparseBlockStructure3D::SparseBlockStructure3D(plint nx, plint ny, plint nz) :
    boundingBox(Box3D(0, nx - 1, 0, ny - 1, 0, nz - 1))
{
    PLB_PRECONDITION(nx >= 1 && ny >= 1 && nz >= 1);
    defaultGridN();
    iniGridParameters();
}

SparseBlockStructure3D::SparseBlockStructure3D(Box3D boundingBox_) : boundingBox(boundingBox_)
{
    PLB_PRECONDITION(
        boundingBox.getNx() >= 1 && boundingBox.getNy() >= 1 && boundingBox.getNz() >= 1);
    defaultGridN();
    iniGridParameters();
}

SparseBlockStructure3D::SparseBlockStructure3D(
    Box3D boundingBox_, plint gridNx_, plint gridNy_, plint gridNz_) :
    boundingBox(boundingBox_), gridNx(gridNx_), gridNy(gridNy_), gridNz(gridNz_)
{
    PLB_PRECONDITION(
        boundingBox.getNx() >= 1 && boundingBox.getNy() >= 1 && boundingBox.getNz() >= 1);
    iniGridParameters();
}

SparseBlockStructure3D::SparseBlockStructure3D(
    SparseBlockStructure3D const &rhs, Box3D boundingBox_) :
    boundingBox(boundingBox_)
{
    PLB_PRECONDITION(
        boundingBox.getNx() >= 1 && boundingBox.getNy() >= 1 && boundingBox.getNz() >= 1);
    gridNx =
        (plint)(0.5 + (double)rhs.gridNx * (double)boundingBox.getNx() / (double)rhs.boundingBox.getNx());
    if (gridNx < 1)
        gridNx = 1;
    gridNy =
        (plint)(0.5 + (double)rhs.gridNy * (double)boundingBox.getNy() / (double)rhs.boundingBox.getNy());
    if (gridNy < 1)
        gridNy = 1;
    gridNz =
        (plint)(0.5 + (double)rhs.gridNz * (double)boundingBox.getNz() / (double)rhs.boundingBox.getNz());
    if (gridNz < 1)
        gridNz = 1;
}

void SparseBlockStructure3D::addBlock(Box3D const &bulk, plint blockId)
{
    addBlock(bulk, bulk, blockId);
}

void SparseBlockStructure3D::addBlock(Box3D const &bulk, Box3D const &uniqueBulk, plint blockId)
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

void SparseBlockStructure3D::removeBlock(plint blockId)
{
    std::map<plint, Box3D>::iterator it = bulks.find(blockId);
    if (it != bulks.end()) {
        extractBlock(blockId);
        bulks.erase(it);
        uniqueBulks.erase(blockId);
    }
}

bool SparseBlockStructure3D::exists(plint blockId)
{
    return bulks.find(blockId) != bulks.end();
}

plint SparseBlockStructure3D::nextIncrementalId() const
{
    if (bulks.empty()) {
        return 0;
    } else {
        return (--bulks.end())->first + 1;
    }
}

Box3D SparseBlockStructure3D::getBoundingBox() const
{
    return boundingBox;
}

bool SparseBlockStructure3D::getBulk(plint blockId, Box3D &bulk) const
{
    std::map<plint, Box3D>::const_iterator it = bulks.find(blockId);
    if (it != bulks.end()) {
        bulk = it->second;
        return true;
    }
    return false;
}

bool SparseBlockStructure3D::getUniqueBulk(plint blockId, Box3D &uniqueBulk) const
{
    std::map<plint, Box3D>::const_iterator it = uniqueBulks.find(blockId);
    if (it != uniqueBulks.end()) {
        uniqueBulk = it->second;
        return true;
    }
    return false;
}

plint SparseBlockStructure3D::getNumBlocks() const
{
    return (plint)bulks.size();
}

plint SparseBlockStructure3D::getNumBulkCells() const
{
    plint nCells = 0;
    std::map<plint, Box3D>::const_iterator it = bulks.begin();
    for (; it != bulks.end(); ++it) {
        nCells += it->second.nCells();
    }
    return nCells;
}

std::map<plint, Box3D> const &SparseBlockStructure3D::getBulks() const
{
    return bulks;
}

std::vector<plint> SparseBlockStructure3D::getLocalBlocks(
    ThreadAttribution const &attribution) const
{
    std::vector<plint> localBlocks;
    std::map<plint, Box3D>::const_iterator it = bulks.begin();
    for (; it != bulks.end(); ++it) {
        plint blockId = it->first;
        if (attribution.isLocal(blockId)) {
            localBlocks.push_back(blockId);
        }
    }
    return localBlocks;
}

plint SparseBlockStructure3D::locate(plint iX, plint iY, plint iZ) const
{
    GridT::const_iterator gridIter = grid.find(Dot3D(gridPosX(iX), gridPosY(iY), gridPosZ(iZ)));
    if (gridIter != grid.end()) {
        std::vector<plint> const &blockList = gridIter->second;
        for (pluint iBlock = 0; iBlock < blockList.size(); ++iBlock) {
            Box3D const &bulk = bulks.find(blockList[iBlock])->second;
            if (contained(iX, iY, iZ, bulk)) {
                return blockList[iBlock];
            }
        }
    }
    return -1;
}

void SparseBlockStructure3D::defaultGridN()
{
    double uniformGridN = std::pow((double)defaultMultiBlockPolicy3D().getNumGridPoints(), 1. / 3.);
    double uniformNcell = std::pow((double)(boundingBox.nCells()), 1. / 3.);
    gridNx = (plint)(0.5 + (double)boundingBox.getNx() / uniformNcell * uniformGridN);
    if (gridNx < 1)
        gridNx = 1;
    gridNy = (plint)(0.5 + (double)boundingBox.getNy() / uniformNcell * uniformGridN);
    if (gridNy < 1)
        gridNy = 1;
    gridNz = (plint)(0.5 + (double)boundingBox.getNz() / uniformNcell * uniformGridN);
    if (gridNz < 1)
        gridNz = 1;
}

void SparseBlockStructure3D::iniGridParameters()
{
    gridLx = boundingBox.getNx() / gridNx;
    if (boundingBox.getNx() % gridNx != 0) {
        ++gridLx;
    }
    gridLy = boundingBox.getNy() / gridNy;
    if (boundingBox.getNy() % gridNy != 0) {
        ++gridLy;
    }
    gridLz = boundingBox.getNz() / gridNz;
    if (boundingBox.getNz() % gridNz != 0) {
        ++gridLz;
    }
}

plint SparseBlockStructure3D::gridPosX(plint realX) const
{
    return (realX - boundingBox.x0) / gridLx;
}

plint SparseBlockStructure3D::gridPosY(plint realY) const
{
    return (realY - boundingBox.y0) / gridLy;
}

plint SparseBlockStructure3D::gridPosZ(plint realZ) const
{
    return (realZ - boundingBox.z0) / gridLz;
}

Box3D SparseBlockStructure3D::getGridBox(Box3D const &realBlock) const
{
    return Box3D(
        gridPosX(realBlock.x0), gridPosX(realBlock.x1), gridPosY(realBlock.y0),
        gridPosY(realBlock.y1), gridPosZ(realBlock.z0), gridPosZ(realBlock.z1));
}

void SparseBlockStructure3D::intersect(
    Box3D const &bulk, std::vector<plint> &ids, std::vector<Box3D> &intersections) const
{
    Box3D gridBox = getGridBox(bulk);
    Box3D intersection;  // Temporary variable.

    std::set<plint> idsToTest;
    for (plint gridX = gridBox.x0; gridX <= gridBox.x1; ++gridX) {
        for (plint gridY = gridBox.y0; gridY <= gridBox.y1; ++gridY) {
            for (plint gridZ = gridBox.z0; gridZ <= gridBox.z1; ++gridZ) {
                GridT::const_iterator gridIter = grid.find(Dot3D(gridX, gridY, gridZ));
                if (gridIter != grid.end()) {
                    std::vector<plint> const &blockList = gridIter->second;
                    idsToTest.insert(blockList.begin(), blockList.end());
                }
            }
        }
    }

    std::set<plint>::const_iterator it = idsToTest.begin();
    for (; it != idsToTest.end(); ++it) {
        Box3D testBlock = bulks.find(*it)->second;
        if (plb::intersect(bulk, testBlock, intersection)) {
            intersections.push_back(intersection);
            ids.push_back(*it);
        }
    }
}

void SparseBlockStructure3D::computeEnvelopeTerm(
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

void SparseBlockStructure3D::computeOverlaps(
    Box3D const &bulk, plint dx, plint dy, plint dz, std::vector<plint> &ids,
    std::vector<Box3D> &overlapsOnBulk, std::vector<Box3D> &overlapsOnNeighbors) const
{
    Box3D envelope;
    std::vector<Box3D> intersections;
    computeEnvelopeTerm(bulk.x0, bulk.x1, envelope.x0, envelope.x1, dx);
    computeEnvelopeTerm(bulk.y0, bulk.y1, envelope.y0, envelope.y1, dy);
    computeEnvelopeTerm(bulk.z0, bulk.z1, envelope.z0, envelope.z1, dz);
    intersect(envelope, ids, intersections);
    overlapsOnNeighbors.insert(
        overlapsOnNeighbors.end(), intersections.begin(), intersections.end());
    for (pluint iOverlap = 0; iOverlap < intersections.size(); ++iOverlap) {
        overlapsOnBulk.push_back(intersections[iOverlap].shift(-dx, -dy, -dz));
    }
}

void SparseBlockStructure3D::computeOverlaps(
    plint blockId, plint envelopeWidth, std::vector<plint> &ids, std::vector<Box3D> &overlapsOnBulk,
    std::vector<Box3D> &overlapsOnNeighbors) const
{
    Box3D bulk = bulks.find(blockId)->second;
    for (plint dx = -1; dx <= +1; ++dx) {
        for (plint dy = -1; dy <= +1; ++dy) {
            for (plint dz = -1; dz <= +1; ++dz) {
                if (!(dx == 0 && dy == 0 && dz == 0)) {
                    computeOverlaps(
                        bulk, dx * envelopeWidth, dy * envelopeWidth, dz * envelopeWidth, ids,
                        overlapsOnBulk, overlapsOnNeighbors);
                }
            }
        }
    }
}

void SparseBlockStructure3D::findNeighbors(
    Box3D const &bulk, plint neighborhoodWidth, std::vector<plint> &neighbors,
    plint excludeId) const
{
    Box3D extendedBlock(bulk.enlarge(neighborhoodWidth));
    Box3D gridBox = getGridBox(extendedBlock);

    std::set<plint> idsToTest;
    for (plint gridX = gridBox.x0; gridX <= gridBox.x1; ++gridX) {
        for (plint gridY = gridBox.y0; gridY <= gridBox.y1; ++gridY) {
            for (plint gridZ = gridBox.z0; gridZ <= gridBox.z1; ++gridZ) {
                GridT::const_iterator gridIter = grid.find(Dot3D(gridX, gridY, gridZ));
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
        Box3D testBlock = bulks.find(*it)->second;
        if (plb::doesIntersect(extendedBlock, testBlock)) {
            neighbors.push_back(*it);
        }
    }
}

void SparseBlockStructure3D::findNeighbors(
    plint blockId, plint neighborhoodWidth, std::vector<plint> &neighbors) const
{
    Box3D bulk = bulks.find(blockId)->second;
    findNeighbors(bulk, neighborhoodWidth, neighbors, blockId);
}

void SparseBlockStructure3D::swap(SparseBlockStructure3D &rhs)
{
    std::swap(boundingBox, rhs.boundingBox);
    std::swap(gridLx, rhs.gridLx);
    std::swap(gridLy, rhs.gridLy);
    std::swap(gridLz, rhs.gridLz);
    std::swap(gridNx, rhs.gridNx);
    std::swap(gridNy, rhs.gridNy);
    std::swap(gridNz, rhs.gridNz);
    grid.swap(rhs.grid);
    bulks.swap(rhs.bulks);
    uniqueBulks.swap(rhs.uniqueBulks);
}

bool SparseBlockStructure3D::equals(SparseBlockStructure3D const &rhs) const
{
    return boundingBox == rhs.boundingBox && gridNx == rhs.gridNx && gridNy == rhs.gridNy
           && gridNz == rhs.gridNz && grid == rhs.grid && bulks == rhs.bulks;
}

void SparseBlockStructure3D::integrateBlock(plint blockId, Box3D bulk)
{
    Box3D gridBox = getGridBox(bulk);

    for (plint gridX = gridBox.x0; gridX <= gridBox.x1; ++gridX) {
        for (plint gridY = gridBox.y0; gridY <= gridBox.y1; ++gridY) {
            for (plint gridZ = gridBox.z0; gridZ <= gridBox.z1; ++gridZ) {
                grid[Dot3D(gridX, gridY, gridZ)].push_back(blockId);
            }
        }
    }
}

void SparseBlockStructure3D::extractBlock(plint blockId)
{
    Box3D const &bulk = bulks[blockId];
    Box3D gridBox = getGridBox(bulk);

    for (plint gridX = gridBox.x0; gridX <= gridBox.x1; ++gridX) {
        for (plint gridY = gridBox.y0; gridY <= gridBox.y1; ++gridY) {
            for (plint gridZ = gridBox.z0; gridZ <= gridBox.z1; ++gridZ) {
                std::vector<plint> &blockList = grid[Dot3D(gridX, gridY, gridZ)];
                // Use remove-erase idiom (because blockList is a std::vector).
                blockList.erase(
                    std::remove(blockList.begin(), blockList.end(), blockId), blockList.end());
            }
        }
    }
}

SparseBlockStructure3D scale(SparseBlockStructure3D const &sparseBlock, plint relativeLevel)
{
    SparseBlockStructure3D newSparseBlock(global::getDefaultMultiScaleManager().scaleBox(
        sparseBlock.getBoundingBox(), relativeLevel));
    std::map<plint, Box3D>::const_iterator bulks = sparseBlock.getBulks().begin();
    for (; bulks != sparseBlock.getBulks().end(); ++bulks) {
        plint id = bulks->first;
        Box3D bulk = bulks->second;
        Box3D uniqueBulk;
        sparseBlock.getUniqueBulk(id, uniqueBulk);
        newSparseBlock.addBlock(
            global::getDefaultMultiScaleManager().scaleBox(bulk, relativeLevel),
            global::getDefaultMultiScaleManager().scaleBox(uniqueBulk, relativeLevel), id);
    }
    return newSparseBlock;
}

SparseBlockStructure3D intersect(SparseBlockStructure3D const &sparseBlock, Box3D domain, bool crop)
{
    Box3D newBoundingBox;
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

SparseBlockStructure3D intersect(
    SparseBlockStructure3D const &sparseBlock, Box3D domain, Box3D newBoundingBox)
{
    SparseBlockStructure3D newSparseBlock(newBoundingBox);
    std::vector<plint> ids;
    std::vector<Box3D> intersections;
    sparseBlock.intersect(domain, ids, intersections);
    for (pluint iBox = 0; iBox < ids.size(); ++iBox) {
        Box3D uniqueBulk, intersectedUniqueBulk;
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
SparseBlockStructure3D intersect(
    SparseBlockStructure3D const &sparseBlock1, SparseBlockStructure3D const &sparseBlock2,
    bool crop)
{
    Box3D newBoundingBox;
    if (crop) {
        if (!intersect(
                sparseBlock1.getBoundingBox(), sparseBlock2.getBoundingBox(), newBoundingBox)) {
            // Avoid undefined state if there is no intersection.
            newBoundingBox = bound(sparseBlock1.getBoundingBox(), sparseBlock2.getBoundingBox());
        }
    } else {
        newBoundingBox = bound(sparseBlock1.getBoundingBox(), sparseBlock2.getBoundingBox());
    }
    SparseBlockStructure3D newSparseBlock(newBoundingBox);

    std::vector<plint> ids1;           // temporary variable.
    std::vector<Box3D> intersections;  // temporary variable.

    std::map<plint, Box3D> const &bulks2 = sparseBlock2.getBulks();
    std::map<plint, Box3D>::const_iterator it2 = bulks2.begin();
    for (; it2 != bulks2.end(); ++it2) {
        ids1.clear();
        intersections.clear();
        sparseBlock1.intersect(it2->second, ids1, intersections);
        for (pluint iBox = 0; iBox < ids1.size(); ++iBox) {
            Box3D uniqueBulk, intersectedUniqueBulk;
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

SparseBlockStructure3D extend(
    SparseBlockStructure3D const &sparseBlock, Box3D addedBulk, Box3D addedUniqueBulk,
    std::vector<plint> &newIds)
{
    Box3D newBoundingBox = bound(sparseBlock.getBoundingBox(), addedBulk);
    SparseBlockStructure3D newSparseBlock(newBoundingBox);

    // Simply take over all blocks from the original sparseBlock.
    std::map<plint, Box3D> const &bulks = sparseBlock.getBulks();
    std::map<plint, Box3D>::const_iterator it = bulks.begin();
    for (; it != bulks.end(); ++it) {
        Box3D uniqueBulk;
        sparseBlock.getUniqueBulk(it->first, uniqueBulk);
        newSparseBlock.addBlock(it->second, uniqueBulk, it->first);
    }

    // Remove from the added bulk all domains which are already in
    //   the sparseBlock.
    std::vector<Box3D> intersections;
    std::vector<plint> ids;
    sparseBlock.intersect(addedBulk, ids, intersections);
    std::vector<Box3D> newDomains;
    std::vector<Box3D> tmpNewDomains;
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
        Box3D newBulk = newDomains[iNew];
        // The unique bulks are evaluated as the intersection between
        //   the computed domains and the user-provided unique bulk.
        Box3D newUniqueBulk;
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
SparseBlockStructure3D except(
    SparseBlockStructure3D const &sparseBlock, Box3D exceptedBlock,
    std::map<plint, std::vector<plint> > &remappedIds)
{
    // Simply take over all blocks and the bounding-box from the original sparseBlock.
    SparseBlockStructure3D newSparseBlock(sparseBlock);

    // Find all blocks that intersect with "exceptedBlock"...
    std::vector<Box3D> intersections;
    std::vector<plint> ids;
    sparseBlock.intersect(exceptedBlock, ids, intersections);
    std::vector<Box3D> fragments;
    for (pluint iBlock = 0; iBlock < intersections.size(); ++iBlock) {
        Box3D bulk, uniqueBulk;
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
            Box3D fragmentUniqueBulk;
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
 *  All blocks contained in sparseBlock are also contained in the created structure,
 *  with same blockIds. The blocks of sparseBlock2 are fractioned in order not to in-
 *  tersect with sparseBlock1's blocks, and the blockIds are remapped. The variable
 *  remappedIds maps blockIds from sparseBlock2 to blockIds in the newly created
 *  structure, for all blocks in the new structure which were created from sparseBlock2.
 **/
SparseBlockStructure3D block_union(
    SparseBlockStructure3D const &sparseBlock1, SparseBlockStructure3D const &sparseBlock2,
    std::map<plint, std::vector<plint> > &remappedIds)
{
    Box3D newBoundingBox = bound(sparseBlock1.getBoundingBox(), sparseBlock2.getBoundingBox());
    SparseBlockStructure3D newSparseBlock(newBoundingBox);

    // Simply take over all blocks from sparseBlock1.
    std::map<plint, Box3D> const &bulks1 = sparseBlock1.getBulks();
    std::map<plint, Box3D>::const_iterator it1 = bulks1.begin();
    for (; it1 != bulks1.end(); ++it1) {
        Box3D uniqueBulk1;
        sparseBlock1.getUniqueBulk(it1->first, uniqueBulk1);
        newSparseBlock.addBlock(it1->second, uniqueBulk1, it1->first);
    }

    // From sparseBlock2, take over everything with exception of what's already
    //   in newSparseBlock. This is not so easy, of course.
    std::map<plint, Box3D> const &bulks2 = sparseBlock2.getBulks();
    std::map<plint, Box3D>::const_iterator it2 = bulks2.begin();
    for (; it2 != bulks2.end(); ++it2) {
        Box3D uniqueBulk2;
        sparseBlock2.getUniqueBulk(it2->first, uniqueBulk2);
        std::vector<plint> ids;
        std::vector<Box3D> intersections;
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
            std::vector<Box3D> originalBlocks;
            std::vector<Box3D> tmp;
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
                Box3D intersectedUniqueBulk2;
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

SparseBlockStructure3D alignDistribution3D(
    SparseBlockStructure3D const &originalStructure, SparseBlockStructure3D const &partnerStructure,
    std::vector<plint> &newIds, std::map<plint, std::vector<plint> > &remappedFromPartner)
{
    // Holds the return value.
    SparseBlockStructure3D newSparseBlock(originalStructure);

    std::vector<plint> ids;                                                    // Temporary.
    std::vector<Box3D> intersections, exceptedBlocks, intersectedUniqueBulks;  // Temporary.

    std::vector<plint> overlapRegionIds;
    std::vector<std::vector<Box3D> > overlapRegionBulks;
    std::vector<std::vector<Box3D> > overlapRegionUniqueBulks;

    // Iterate over blocks of partner, in order to remove corresponding areas
    //   in the original structure, and then add them with appropriate new IDs.
    std::map<plint, Box3D> const &partnerBulks = partnerStructure.getBulks();
    std::map<plint, Box3D>::const_iterator partnerIt = partnerBulks.begin();
    for (; partnerIt != partnerBulks.end(); ++partnerIt) {
        Box3D partnerBulk = partnerIt->second;
        ids.clear();
        intersections.clear();
        newSparseBlock.intersect(partnerBulk, ids, intersections);
        intersectedUniqueBulks.clear();
        // Remove areas of overlap (they will be added at the end of this function),
        //   but keep pieces of the original domain which don't overlap.
        for (pluint iInters = 0; iInters < intersections.size(); ++iInters) {
            Box3D intersection = intersections[iInters];
            plint originalId = ids[iInters];
            Box3D originalBulk, originalUniqueBulk;
            newSparseBlock.getBulk(originalId, originalBulk);
            newSparseBlock.getUniqueBulk(originalId, originalUniqueBulk);

            Box3D intersectedUniqueBulk;
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
                Box3D exceptedBlock = exceptedBlocks[iExc];
                Box3D exceptedUniqueBulk;
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
    std::map<plint, Box3D> const &originalBulks = newSparseBlock.getBulks();
    std::map<plint, Box3D>::const_iterator originalIt = originalBulks.begin();
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

EuclideanIterator3D::EuclideanIterator3D(SparseBlockStructure3D const &sparseBlock_) :
    sparseBlock(sparseBlock_)
{ }

bool EuclideanIterator3D::getNextChunkX(
    plint iX, plint iY, plint iZ, plint &blockId, plint &chunkSize) const
{
    blockId = sparseBlock.locate(iX, iY, iZ);
    if (blockId == -1) {
        Box3D boundingBox = sparseBlock.getBoundingBox();
        plint exploreX = iX + 1;
        while (exploreX < boundingBox.getNx() && sparseBlock.locate(exploreX, iY, iZ) == -1) {
            ++exploreX;
        }
        chunkSize = exploreX - iX;
        return false;
    } else {
        Box3D bulk;
        sparseBlock.getBulk(blockId, bulk);
        chunkSize = bulk.x1 - iX + 1;
        return true;
    }
}

bool EuclideanIterator3D::getNextChunkY(
    plint iX, plint iY, plint iZ, plint &blockId, plint &chunkSize) const
{
    blockId = sparseBlock.locate(iX, iY, iZ);
    if (blockId == -1) {
        Box3D boundingBox = sparseBlock.getBoundingBox();
        plint exploreY = iY + 1;
        while (exploreY < boundingBox.getNy() && sparseBlock.locate(iX, exploreY, iZ) == -1) {
            ++exploreY;
        }
        chunkSize = exploreY - iY;
        return false;
    } else {
        Box3D bulk;
        sparseBlock.getBulk(blockId, bulk);
        chunkSize = bulk.y1 - iY + 1;
        return true;
    }
}

bool EuclideanIterator3D::getNextChunkZ(
    plint iX, plint iY, plint iZ, plint &blockId, plint &chunkSize) const
{
    blockId = sparseBlock.locate(iX, iY, iZ);
    if (blockId == -1) {
        Box3D boundingBox = sparseBlock.getBoundingBox();
        plint exploreZ = iZ + 1;
        while (exploreZ < boundingBox.getNz() && sparseBlock.locate(iX, iY, exploreZ) == -1) {
            ++exploreZ;
        }
        chunkSize = exploreZ - iZ;
        return false;
    } else {
        Box3D bulk;
        sparseBlock.getBulk(blockId, bulk);
        chunkSize = bulk.z1 - iZ + 1;
        return true;
    }
}

}  // namespace plb
