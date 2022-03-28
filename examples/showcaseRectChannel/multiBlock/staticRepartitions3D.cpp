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
 * Utilities for 3D multi data distributions -- implementation.
 */

#include "multiBlock/staticRepartitions3D.h"

#include <algorithm>

#include "algorithm/basicAlgorithms.h"
#include "atomicBlock/dataField3D.hh"
#include "core/util.h"

namespace plb {

////////////////////// function createRegularDistribution3D /////////////////////

SparseBlockStructure3D createRegularDistribution3D(
    Box3D const &domain, plint numBlocksX, plint numBlocksY, plint numBlocksZ)
{
    SparseBlockStructure3D dataGeometry(domain);
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
            plint posZ = domain.z0;
            for (plint iBlockZ = 0; iBlockZ < numBlocksZ; ++iBlockZ) {
                plint lz = domain.getNz() / numBlocksZ;
                if (iBlockZ < domain.getNz() % numBlocksZ)
                    ++lz;
                dataGeometry.addBlock(
                    Box3D(posX, posX + lx - 1, posY, posY + ly - 1, posZ, posZ + lz - 1),
                    dataGeometry.nextIncrementalId());
                posZ += lz;
            }
            posY += ly;
        }
        posX += lx;
    }
    return dataGeometry;
}

SparseBlockStructure3D createRegularDistribution3D(
    plint nx, plint ny, plint nz, plint numBlocksX, plint numBlocksY, plint numBlocksZ)
{
    return createRegularDistribution3D(
        Box3D(0, nx - 1, 0, ny - 1, 0, nz - 1), numBlocksX, numBlocksY, numBlocksZ);
}

SparseBlockStructure3D createRegularDistribution3D(Box3D const &domain, int numProc)
{
    std::vector<plint> repartition = algorithm::evenRepartition(numProc, 3);
    std::vector<plint> newRepartition(3);
    if (domain.getNx() > domain.getNy()) {      // nx>ny
        if (domain.getNx() > domain.getNz()) {  // nx>nz
            newRepartition[0] = repartition[0];
            if (domain.getNy() > domain.getNz()) {  // ny>nz
                newRepartition[1] = repartition[1];
                newRepartition[2] = repartition[2];
            } else {  // nz>ny
                newRepartition[1] = repartition[2];
                newRepartition[2] = repartition[1];
            }
        } else {  // nz>nx
            newRepartition[2] = repartition[0];
            newRepartition[1] = repartition[2];
            newRepartition[0] = repartition[1];
        }
    } else {                                    // ny>nx
        if (domain.getNy() > domain.getNz()) {  // ny>nz
            newRepartition[1] = repartition[0];
            if (domain.getNx() > domain.getNz()) {  // nx>nz
                newRepartition[0] = repartition[1];
                newRepartition[2] = repartition[2];
            } else {  // nz>nx
                newRepartition[0] = repartition[2];
                newRepartition[2] = repartition[1];
            }
        } else {  // nz>ny
            newRepartition[2] = repartition[0];
            newRepartition[1] = repartition[1];
            newRepartition[0] = repartition[2];
        }
    }
    return createRegularDistribution3D(
        domain, newRepartition[0], newRepartition[1], newRepartition[2]);
}

SparseBlockStructure3D createRegularDistribution3D(plint nx, plint ny, plint nz, int numProc)
{
    return createRegularDistribution3D(Box3D(0, nx - 1, 0, ny - 1, 0, nz - 1), numProc);
}

SparseBlockStructure3D createRegularDistributionYZ3D(plint nx, plint ny, plint nz, int numProc)
{
    std::vector<plint> repartition = algorithm::evenRepartition(numProc, 2);
    plint repartitionMin = std::min(repartition[0], repartition[1]);
    plint repartitionMax = std::max(repartition[0], repartition[1]);
    plint numBlocksX = 1;
    plint numBlocksY = ny > nz ? repartitionMax : repartitionMin;
    plint numBlocksZ = nz >= ny ? repartitionMax : repartitionMin;

    return createRegularDistribution3D(
        Box3D(0, nx - 1, 0, ny - 1, 0, nz - 1), numBlocksX, numBlocksY, numBlocksZ);
}

SparseBlockStructure3D createRegularDistributionXZ3D(plint nx, plint ny, plint nz, int numProc)
{
    std::vector<plint> repartition = algorithm::evenRepartition(numProc, 2);
    plint repartitionMin = std::min(repartition[0], repartition[1]);
    plint repartitionMax = std::max(repartition[0], repartition[1]);
    plint numBlocksY = 1;
    plint numBlocksX = nx > nz ? repartitionMax : repartitionMin;
    plint numBlocksZ = nz >= nx ? repartitionMax : repartitionMin;

    return createRegularDistribution3D(
        Box3D(0, nx - 1, 0, ny - 1, 0, nz - 1), numBlocksX, numBlocksY, numBlocksZ);
}

SparseBlockStructure3D createRegularDistributionXY3D(plint nx, plint ny, plint nz, int numProc)
{
    std::vector<plint> repartition = algorithm::evenRepartition(numProc, 2);
    plint repartitionMin = std::min(repartition[0], repartition[1]);
    plint repartitionMax = std::max(repartition[0], repartition[1]);
    plint numBlocksZ = 1;
    plint numBlocksX = nx > ny ? repartitionMax : repartitionMin;
    plint numBlocksY = ny >= nx ? repartitionMax : repartitionMin;

    return createRegularDistribution3D(
        Box3D(0, nx - 1, 0, ny - 1, 0, nz - 1), numBlocksX, numBlocksY, numBlocksZ);
}

void mergeIntersections(std::vector<Box3D> &intersections)
{
    std::vector<Box3D> merged;
    for (pluint iInters = 0; iInters < intersections.size(); ++iInters) {
        Box3D intersection = intersections[iInters];
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

plint cumNcells(std::vector<Box3D> const &domains)
{
    plint sum = 0;
    for (pluint iDomain = 0; iDomain < domains.size(); ++iDomain) {
        sum += domains[iDomain].nCells();
    }
    return sum;
}

SparseBlockStructure3D reparallelize(
    SparseBlockStructure3D const &originalStructure, plint blockLx, plint blockLy, plint blockLz)
{
    std::vector<std::pair<plint, plint> > rangesX, rangesY, rangesZ;
    Box3D boundingBox = originalStructure.getBoundingBox();
    util::linearBlockRepartition(boundingBox.x0, boundingBox.x1, blockLx, rangesX);
    util::linearBlockRepartition(boundingBox.y0, boundingBox.y1, blockLy, rangesY);
    util::linearBlockRepartition(boundingBox.z0, boundingBox.z1, blockLz, rangesZ);
    SparseBlockStructure3D newStructure(boundingBox);
    std::vector<plint> ids;
    std::vector<Box3D> intersections;
    for (pluint blockX = 0; blockX < rangesX.size(); ++blockX) {
        for (pluint blockY = 0; blockY < rangesY.size(); ++blockY) {
            for (pluint blockZ = 0; blockZ < rangesZ.size(); ++blockZ) {
                Box3D currentBlock(
                    rangesX[blockX].first, rangesX[blockX].second, rangesY[blockY].first,
                    rangesY[blockY].second, rangesZ[blockZ].first, rangesZ[blockZ].second);
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
    }
    return newStructure;
}

SparseBlockStructure3D reparallelize(SparseBlockStructure3D const &originalStructure)
{
    // Decide that atomic block-components have size 16.
    plint blockLx = 16;
    plint blockLy = 16;
    plint blockLz = 16;
    return reparallelize(originalStructure, blockLx, blockLy, blockLz);
}

SparseBlockStructure3D createXSlicedDistribution3D(
    CellTypeField3D const &cellTypeField, plint numBlocks)
{
    plint nX = cellTypeField.getNx();
    plint nY = cellTypeField.getNy();
    plint nZ = cellTypeField.getNz();

    std::vector<plint> numActivePerSlice;
    plint numActiveTotal = 0;
    for (plint iX = 0; iX < nX; iX++) {
        plint numActiveCurrentSlice = 0;
        for (plint iY = 0; iY < nY; iY++) {
            for (plint iZ = 0; iZ < nZ; iZ++) {
                if (cellTypeField.get(iX, iY, iZ) > 0)
                    numActiveCurrentSlice++;
            }
        }
        numActivePerSlice.push_back(numActiveCurrentSlice);
        numActiveTotal += numActiveCurrentSlice;
    }
    plint numActivePerBlock = numActiveTotal / numBlocks;

    SparseBlockStructure3D dataGeometry(nX, nY, nZ);

    plint iX = 0;
    for (plint iBlock = 0; iBlock < numBlocks; ++iBlock) {
        plint posX = iX;
        plint numActiveCurrentBlock = 0;
        while (numActiveCurrentBlock < numActivePerBlock && iX < nX) {
            numActiveCurrentBlock += numActivePerSlice[iX];
            iX++;
        }
        dataGeometry.addBlock(
            Box3D(posX, iX - 1, 0, nY - 1, 0, nZ - 1), dataGeometry.nextIncrementalId());
    }
    return dataGeometry;
}

SparseBlockStructure3D createYSlicedDistribution3D(
    CellTypeField3D const &cellTypeField, plint numBlocks)
{
    plint nX = cellTypeField.getNx();
    plint nY = cellTypeField.getNy();
    plint nZ = cellTypeField.getNz();

    std::vector<plint> numActivePerSlice;
    plint numActiveTotal = 0;
    for (plint iY = 0; iY < nY; iY++) {
        plint numActiveCurrentSlice = 0;
        for (plint iX = 0; iX < nX; iX++) {
            for (plint iZ = 0; iZ < nZ; iZ++) {
                if (cellTypeField.get(iX, iY, iZ) > 0)
                    numActiveCurrentSlice++;
            }
        }
        numActivePerSlice.push_back(numActiveCurrentSlice);
        numActiveTotal += numActiveCurrentSlice;
    }
    plint numActivePerBlock = numActiveTotal / numBlocks;

    SparseBlockStructure3D dataGeometry(nX, nY, nZ);

    plint iY = 0;
    for (plint iBlock = 0; iBlock < numBlocks; ++iBlock) {
        plint posY = iY;
        plint numActiveCurrentBlock = 0;
        while (numActiveCurrentBlock < numActivePerBlock && iY < nY) {
            numActiveCurrentBlock += numActivePerSlice[iY];
            iY++;
        }
        dataGeometry.addBlock(
            Box3D(0, nX - 1, posY, iY - 1, 0, nZ - 1), dataGeometry.nextIncrementalId());
    }
    return dataGeometry;
}

SparseBlockStructure3D createZSlicedDistribution3D(
    CellTypeField3D const &cellTypeField, plint numBlocks)
{
    plint nX = cellTypeField.getNx();
    plint nY = cellTypeField.getNy();
    plint nZ = cellTypeField.getNz();

    std::vector<plint> numActivePerSlice;
    plint numActiveTotal = 0;
    for (plint iZ = 0; iZ < nZ; iZ++) {
        plint numActiveCurrentSlice = 0;
        for (plint iX = 0; iX < nX; iX++) {
            for (plint iY = 0; iY < nY; iY++) {
                if (cellTypeField.get(iX, iY, iZ) > 0)
                    numActiveCurrentSlice++;
            }
        }
        numActivePerSlice.push_back(numActiveCurrentSlice);
        numActiveTotal += numActiveCurrentSlice;
    }
    plint numActivePerBlock = numActiveTotal / numBlocks;

    SparseBlockStructure3D dataGeometry(nX, nY, nZ);

    plint iZ = 0;
    for (plint iBlock = 0; iBlock < numBlocks; ++iBlock) {
        plint posZ = iZ;
        plint numActiveCurrentBlock = 0;
        while (numActiveCurrentBlock < numActivePerBlock && iZ < nZ) {
            numActiveCurrentBlock += numActivePerSlice[iZ];
            iZ++;
        }
        dataGeometry.addBlock(
            Box3D(0, nX - 1, 0, nY - 1, posZ, iZ - 1), dataGeometry.nextIncrementalId());
    }
    return dataGeometry;
}

SparseBlockStructure3D createXSlicedDistribution3D(CellTypeField3D const &cellTypeField)
{
    return createXSlicedDistribution3D(cellTypeField, global::mpi().getSize());
}

SparseBlockStructure3D createYSlicedDistribution3D(CellTypeField3D const &cellTypeField)
{
    return createYSlicedDistribution3D(cellTypeField, global::mpi().getSize());
}

SparseBlockStructure3D createZSlicedDistribution3D(CellTypeField3D const &cellTypeField)
{
    return createZSlicedDistribution3D(cellTypeField, global::mpi().getSize());
}

}  // namespace plb
