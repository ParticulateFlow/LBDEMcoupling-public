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
 * Helper functions for domain initialization -- header file.
 */
#include "multiBlock/multiBlockInfo3D.h"

namespace plb {

bool getMultiBlockInfo(
    MultiBlock3D const &multiBlock, plint &nx, plint &ny, plint &nz, plint &numBlocks,
    Box3D &smallest, Box3D &largest, plint &numAllocatedCells)
{
    nx = multiBlock.getNx();
    ny = multiBlock.getNy();
    nz = multiBlock.getNz();

    MultiBlockManagement3D const &management = multiBlock.getMultiBlockManagement();
    SparseBlockStructure3D const &sparseBlock = management.getSparseBlockStructure();
    if (sparseBlock.getNumBlocks() == 0) {
        return false;
    }
    plint firstBulk = sparseBlock.getBulks().begin()->first;
    plint maxNumCells = management.getBulk(firstBulk).nCells();
    plint minNumCells = management.getBulk(firstBulk).nCells();
    numAllocatedCells = 0;
    numBlocks = sparseBlock.getNumBlocks();
    std::map<plint, Box3D>::const_iterator it = sparseBlock.getBulks().begin();
    plint smallestBlock = it->first;
    plint largestBlock = it->first;
    for (; it != sparseBlock.getBulks().end(); ++it) {
        plint numCells = it->second.nCells();
        numAllocatedCells += numCells;
        if (numCells > maxNumCells) {
            maxNumCells = numCells;
            largestBlock = it->first;
        }
        if (numCells < minNumCells) {
            minNumCells = numCells;
            smallestBlock = it->first;
        }
    }

    smallest = management.getBulk(smallestBlock);
    largest = management.getBulk(largestBlock);
    return true;
}

std::string getMultiBlockInfo(MultiBlock3D const &multiBlock)
{
    plint nx, ny, nz;
    plint numBlocks;
    plint numAllocatedCells;
    Box3D smallest, largest;
    if (!getMultiBlockInfo(multiBlock, nx, ny, nz, numBlocks, smallest, largest, numAllocatedCells))
    {
        return std::string("Empty multi-block\n");
    }
    std::stringstream blockInfo;
    blockInfo << "Size of the multi-block:     " << nx << "-by-" << ny << "-by-" << nz << "\n";
    blockInfo << "Number of atomic-blocks:     " << numBlocks << "\n";
    blockInfo << "Smallest atomic-block:       " << smallest.getNx() << "-by-" << smallest.getNy()
              << "-by-" << smallest.getNz() << "\n";
    blockInfo << "Largest atomic-block:        " << largest.getNx() << "-by-" << largest.getNy()
              << "-by-" << largest.getNz() << "\n";
    blockInfo << "Number of allocated cells:   " << (double)numAllocatedCells / 1.e6
              << " million\n";
    blockInfo << "Fraction of allocated domain: "
              << (double)numAllocatedCells / (double)(nx * ny * nz) * 100 << " percent\n";

    return blockInfo.str();
}

}  // namespace plb
