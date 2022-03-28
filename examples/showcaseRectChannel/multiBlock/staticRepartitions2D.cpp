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
 * Utilities for 2D multi data distributions -- implementation.
 */

#include "multiBlock/staticRepartitions2D.h"

#include <algorithm>

#include "algorithm/basicAlgorithms.h"
#include "atomicBlock/dataField2D.hh"
#include "core/util.h"

namespace plb {

////////////////////// function createRegularDistribution2D /////////////////////

SparseBlockStructure2D createRegularDistribution2D(
    Box2D const &domain, plint numBlocksX, plint numBlocksY)
{
    SparseBlockStructure2D dataGeometry(domain);
    plint posX = domain.x0;
    for (plint iBlockX = 0; iBlockX < numBlocksX; ++iBlockX) {
        plint lx = domain.getNx() / numBlocksX;
        if (iBlockX < domain.getNx() % numBlocksX)
            ++lx;
        plint posY = domain.y0;
        for (plint iBlockY = 0; iBlockY < numBlocksY; ++iBlockY) {
            plint ly = domain.getNy() / numBlocksY;
            if (iBlockY < domain.getNy() % numBlocksY)
                ++ly;
            dataGeometry.addBlock(
                Box2D(posX, posX + lx - 1, posY, posY + ly - 1), dataGeometry.nextIncrementalId());
            posY += ly;
        }
        posX += lx;
    }
    return dataGeometry;
}

SparseBlockStructure2D createRegularDistribution2D(
    plint nx, plint ny, plint numBlocksX, plint numBlocksY)
{
    return createRegularDistribution2D(Box2D(0, nx - 1, 0, ny - 1), numBlocksX, numBlocksY);
}

SparseBlockStructure2D createRegularDistribution2D(Box2D const &domain, int numProc)
{
    std::vector<plint> repartition = algorithm::evenRepartition(numProc, 2);
    return createRegularDistribution2D(domain, repartition[0], repartition[1]);
}

SparseBlockStructure2D createRegularDistribution2D(plint nx, plint ny, int numProc)
{
    return createRegularDistribution2D(Box2D(0, nx - 1, 0, ny - 1), numProc);
}

void mergeIntersections(std::vector<Box2D> &intersections)
{
    std::vector<Box2D> merged;
    for (pluint iInters = 0; iInters < intersections.size(); ++iInters) {
        Box2D intersection = intersections[iInters];
        bool hasMerged = false;
        for (pluint iMerged = 0; iMerged < merged.size(); ++iMerged) {
            if (merge(merged[iMerged], intersection)) {
                hasMerged = true;
                break;
            }
        }
        if (!hasMerged) {
            merged.push_back(intersection);
        }
    }
    intersections.swap(merged);
}

plint cumNcells(std::vector<Box2D> const &domains)
{
    plint sum = 0;
    for (pluint iDomain = 0; iDomain < domains.size(); ++iDomain) {
        sum += domains[iDomain].nCells();
    }
    return sum;
}

SparseBlockStructure2D reparallelize(
    SparseBlockStructure2D const &originalStructure, plint blockLx, plint blockLy)
{
    std::vector<std::pair<plint, plint> > rangesX, rangesY;
    Box2D boundingBox = originalStructure.getBoundingBox();
    util::linearBlockRepartition(boundingBox.x0, boundingBox.x1, blockLx, rangesX);
    util::linearBlockRepartition(boundingBox.y0, boundingBox.y1, blockLy, rangesY);
    SparseBlockStructure2D newStructure(boundingBox);
    std::vector<plint> ids;
    std::vector<Box2D> intersections;
    for (pluint blockX = 0; blockX < rangesX.size(); ++blockX) {
        for (pluint blockY = 0; blockY < rangesY.size(); ++blockY) {
            Box2D currentBlock(
                rangesX[blockX].first, rangesX[blockX].second, rangesY[blockY].first,
                rangesY[blockY].second);
            ids.clear();
            intersections.clear();
            originalStructure.intersect(currentBlock, ids, intersections);
            // It is possible that the current block fully covers the domain of the old
            // distribution. In this case, simply add current block, in order to avoid
            // fragmentation. Note that this explicit test is really necessary, because
            // the function mergeIntersection, which is called below, is not always able
            // to reconstruct a full block from its fragments.
            if (currentBlock.nCells() == cumNcells(intersections)) {
                plint nextId = newStructure.nextIncrementalId();
                newStructure.addBlock(currentBlock, nextId);
            } else {
                // Construct bigger blocks if possible, in order to avoid fragmentation.
                mergeIntersections(intersections);
                for (pluint iInters = 0; iInters < intersections.size(); ++iInters) {
                    plint nextId = newStructure.nextIncrementalId();
                    newStructure.addBlock(intersections[iInters], nextId);
                }
            }
        }
    }
    return newStructure;
}

SparseBlockStructure2D reparallelize(SparseBlockStructure2D const &originalStructure)
{
    // Decide that atomic block-components have size 64.
    plint blockLx = 64;
    plint blockLy = 64;
    return reparallelize(originalStructure, blockLx, blockLy);
}

SparseBlockStructure2D createXSlicedMultiBlockDistribution2D(
    CellTypeField2D const &cellTypeField, plint numBlocks)
{
    plint nX = cellTypeField.getNx();
    plint nY = cellTypeField.getNy();

    std::vector<plint> numActivePerSlice;
    plint numActiveTotal = 0;
    for (plint iX = 0; iX < nX; iX++) {
        plint numActiveCurrentSlice = 0;
        for (plint iY = 0; iY < nY; iY++) {
            if (cellTypeField.get(iX, iY) > 0)
                numActiveCurrentSlice++;
        }
        numActivePerSlice.push_back(numActiveCurrentSlice);
        numActiveTotal += numActiveCurrentSlice;
    }
    plint numActivePerBlock = numActiveTotal / numBlocks;

    SparseBlockStructure2D dataGeometry(nX, nY);

    plint iX = 0;
    for (plint iBlock = 0; iBlock < numBlocks; ++iBlock) {
        plint posX = iX;
        plint numActiveCurrentBlock = 0;
        while (numActiveCurrentBlock < numActivePerBlock && iX < nX) {
            numActiveCurrentBlock += numActivePerSlice[iX];
            iX++;
        }
        dataGeometry.addBlock(Box2D(posX, iX - 1, 0, nY - 1), dataGeometry.nextIncrementalId());
    }
    return dataGeometry;
}

SparseBlockStructure2D createYSlicedMultiBlockDistribution2D(
    CellTypeField2D const &cellTypeField, plint numBlocks)
{
    plint nX = cellTypeField.getNx();
    plint nY = cellTypeField.getNy();

    std::vector<plint> numActivePerSlice;
    plint numActiveTotal = 0;
    for (plint iY = 0; iY < nY; iY++) {
        plint numActiveCurrentSlice = 0;
        for (plint iX = 0; iX < nX; iX++) {
            if (cellTypeField.get(iX, iY) > 0)
                numActiveCurrentSlice++;
        }
        numActivePerSlice.push_back(numActiveCurrentSlice);
        numActiveTotal += numActiveCurrentSlice;
    }
    plint numActivePerBlock = numActiveTotal / numBlocks;

    SparseBlockStructure2D dataGeometry(nX, nY);

    plint iY = 0;
    for (plint iBlock = 0; iBlock < numBlocks; ++iBlock) {
        plint posY = iY;
        plint numActiveCurrentBlock = 0;
        while (numActiveCurrentBlock < numActivePerBlock && iY < nY) {
            numActiveCurrentBlock += numActivePerSlice[iY];
            iY++;
        }
        dataGeometry.addBlock(Box2D(0, nX - 1, posY, iY - 1), dataGeometry.nextIncrementalId());
    }
    return dataGeometry;
}

SparseBlockStructure2D createXSlicedMultiBlockDistribution2D(CellTypeField2D const &cellTypeField)
{
    return createXSlicedMultiBlockDistribution2D(cellTypeField, global::mpi().getSize());
}

SparseBlockStructure2D createYSlicedMultiBlockDistribution2D(CellTypeField2D const &cellTypeField)
{
    return createYSlicedMultiBlockDistribution2D(cellTypeField, global::mpi().getSize());
}

}  // namespace plb
