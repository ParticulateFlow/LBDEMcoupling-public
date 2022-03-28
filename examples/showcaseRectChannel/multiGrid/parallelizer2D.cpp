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

#include "multiGrid/parallelizer2D.h"

#include "io/parallelIO.h"
#include "multiGrid/multiScale.h"
#include "parallelism/mpiManager.h"

namespace plb {

void Parallelizer2D::parallelizeLevel(
    plint whichLevel, std::vector<std::vector<Box2D> > const &originalBlocks,
    std::vector<Box2D> const &parallelRegions, std::vector<plint> const &regionIDs)
{
    PLB_PRECONDITION(parallelRegions.size() == regionIDs.size());
    PLB_PRECONDITION(whichLevel < (plint)originalBlocks.size());
    std::vector<Box2D> newBlocks;
    // IDs are going to be reattributed at the level whichLevel.
    if (finalMpiDistribution.size() <= (pluint)whichLevel) {
        finalMpiDistribution.resize(whichLevel + 1);
    }
    for (pluint iRegion = 0; iRegion < parallelRegions.size(); ++iRegion) {
        plint currentId = regionIDs[iRegion];
        for (pluint iBlock = 0; iBlock < originalBlocks[whichLevel].size(); ++iBlock) {
            Box2D intersection;
            if (intersect(
                    originalBlocks[whichLevel][iBlock], parallelRegions[iRegion], intersection)) {
                newBlocks.push_back(intersection);
                finalMpiDistribution[whichLevel].push_back(currentId);
            }
        }
    }

    recomputedBlocks[whichLevel].insert(
        recomputedBlocks[whichLevel].end(), newBlocks.begin(), newBlocks.end());
}

plint Parallelizer2D::computeCost(std::vector<std::vector<Box2D> > const &originalBlocks, Box2D box)
{
    plint totalCost = 0;
    plint numLevels = originalBlocks.size();

    for (plint iLevel = (plint)originalBlocks.size() - 1; iLevel >= 0; --iLevel) {
        // convert the box to the current level
        Box2D levelBox =
            global::getDefaultMultiScaleManager().scaleBox(box, iLevel - (numLevels - 1));
        for (pluint iComp = 0; iComp < originalBlocks[iLevel].size(); ++iComp) {
            Box2D currentBox;
            if (intersect(originalBlocks[iLevel][iComp], levelBox, currentBox)) {
                plint volume = currentBox.getNx() * currentBox.getNy();
                totalCost += (plint)util::twoToThePower(iLevel) * volume;
            }
        }
    }

    return totalCost;
}

/* ************* ParallellizeBySquares2D **************** */
ParallellizeBySquares2D::ParallellizeBySquares2D(
    std::vector<std::vector<Box2D> > const &originalBlocks_, Box2D finestBoundingBox_,
    plint xTiles_, plint yTiles_) :
    originalBlocks(originalBlocks_),
    finestBoundingBox(finestBoundingBox_),
    processorNumber(global::mpi().getSize()),
    xTiles(xTiles_),
    yTiles(yTiles_)
{
    // divide the finest bounding box in xTiles by yTiles squares
    computeFinestDivision(xTiles, yTiles);
}

void ParallellizeBySquares2D::computeFinestDivision(plint xTiles, plint yTiles)
{
    std::vector<std::pair<plint, plint> > rangesX;
    std::vector<std::pair<plint, plint> > rangesY;

    plint nx = finestBoundingBox.x1;
    plint ny = finestBoundingBox.y1;

    util::linearRepartition(0, nx, xTiles, rangesX);
    util::linearRepartition(0, ny, yTiles, rangesY);

    finestDivision.resize(0);
    mpiDistribution.resize(xTiles * yTiles);

    for (plint iX = 0; iX < (plint)rangesX.size(); ++iX) {
        for (plint iY = 0; iY < (plint)rangesY.size(); ++iY) {
            // create a Box2D with coordinates for the new sector
            Box2D box(rangesX[iX].first, rangesX[iX].second, rangesY[iY].first, rangesY[iY].second);
            // put it in the finestDivision
            finestDivision.push_back(box);
        }
    }
}

void ParallellizeBySquares2D::parallelize()
{
    // we must have the same number of blocks as processors
    PLB_PRECONDITION(xTiles * yTiles == processorNumber);

    plint totalCost = computeCost(originalBlocks, finestBoundingBox);
    plint idealCostPerProcessor = totalCost / processorNumber;

    pcout << "Total cost of computations = " << totalCost << std::endl;
    pcout << "We are using " << processorNumber << " processors...\n";
    pcout << "Ideal cost per processor = " << idealCostPerProcessor << std::endl;

    std::vector<plint> totalCosts(processorNumber);

    plint total = 0;
    for (plint iProc = 0; iProc < processorNumber; ++iProc) {
        plint blockCost = computeCost(originalBlocks, finestDivision[iProc]);
        totalCosts[iProc] += blockCost;
        mpiDistribution[iProc] = iProc;
        total += blockCost;
    }

    pcout << "---- Costs Per Processor ----\n";

    for (pluint i = 0; i < totalCosts.size(); ++i) {
        pcout << i << " : " << totalCosts[i] << std::endl;
        // check if everyone is doing something
        if (totalCosts[i] == 0) {
            pcout << "\t >> processor " << i << " does not have work to do. Exiting.....\n";
            std::exit(1);
        }
    }

    pcout << "*******************************\n";
    pcout << "Sum of all costs = " << total << std::endl;
    pcout << "*******************************\n";

    // convert the original blocks to the new blocks
    recomputedBlocks.resize(originalBlocks.size());
    finalMpiDistribution.resize(originalBlocks.size());

    plint finestLevel = (plint)originalBlocks.size() - 1;
    for (plint iLevel = finestLevel; iLevel >= 0; --iLevel) {
        parallelizeLevel(iLevel, originalBlocks, finestDivision, mpiDistribution);
        // Adapt the regions to the next-coarser level.
        for (pluint iRegion = 0; iRegion < finestDivision.size(); ++iRegion) {
            finestDivision[iRegion] = finestDivision[iRegion].divideAndFitSmaller(2);
        }
    }
}

}  // namespace plb
