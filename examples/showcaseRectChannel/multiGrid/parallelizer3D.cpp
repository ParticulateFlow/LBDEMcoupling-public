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

#ifndef PARALLELIZER_3D_CPP
#define PARALLELIZER_3D_CPP

#include "multiGrid/parallelizer3D.h"

#include "io/parallelIO.h"
#include "multiGrid/multiScale.h"
#include "parallelism/mpiManager.h"

namespace plb {

void Parallelizer3D::parallelizeLevel(
    plint whichLevel, std::vector<std::vector<Box3D> > const &originalBlocks,
    std::vector<Box3D> const &parallelRegions, std::vector<plint> const &regionIDs)
{
    PLB_PRECONDITION(parallelRegions.size() == regionIDs.size());
    PLB_PRECONDITION(whichLevel < (plint)originalBlocks.size());
    std::vector<Box3D> newBlocks;
    // IDs are going to be reattributed at the level whichLevel.
    if (finalMpiDistribution.size() <= (pluint)whichLevel) {
        finalMpiDistribution.resize(whichLevel + 1);
    }
    for (pluint iRegion = 0; iRegion < parallelRegions.size(); ++iRegion) {
        plint currentId = regionIDs[iRegion];
        for (pluint iBlock = 0; iBlock < originalBlocks[whichLevel].size(); ++iBlock) {
            Box3D intersection;
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

plint Parallelizer3D::computeCost(std::vector<std::vector<Box3D> > const &originalBlocks, Box3D box)
{
    plint totalCost = 0;
    plint numLevels = originalBlocks.size();

    for (plint iLevel = (plint)originalBlocks.size() - 1; iLevel >= 0; --iLevel) {
        // convert the box to the current level
        Box3D levelBox =
            global::getDefaultMultiScaleManager().scaleBox(box, iLevel - (numLevels - 1));
        for (pluint iComp = 0; iComp < originalBlocks[iLevel].size(); ++iComp) {
            Box3D currentBox;
            if (intersect(originalBlocks[iLevel][iComp], levelBox, currentBox)) {
                plint volume = currentBox.getNx() * currentBox.getNy() * currentBox.getNz();
                totalCost += (plint)util::twoToThePower(iLevel) * volume;
            }
        }
    }

    return totalCost;
}

/* ************* ParallellizeByCubes3D ************************************** */
ParallellizeByCubes3D::ParallellizeByCubes3D(
    std::vector<std::vector<Box3D> > const &originalBlocks_, Box3D finestBoundingBox_,
    plint xTiles_, plint yTiles_, plint zTiles_, bool optimiseDivision) :
    originalBlocks(originalBlocks_),
    finestBoundingBox(finestBoundingBox_),
    processorNumber(global::mpi().getSize()),
    xTiles(xTiles_),
    yTiles(yTiles_),
    zTiles(zTiles_)
{
    // divide the finest bounding box in xTiles by yTiles by zTiles cubes
    if (!optimiseDivision) {
        computeFinestDivision();
    } else {
        computeDivision();
    }
}

void ParallellizeByCubes3D::computeFinestDivision()
{
    std::vector<std::pair<plint, plint> > rangesX;
    std::vector<std::pair<plint, plint> > rangesY;
    std::vector<std::pair<plint, plint> > rangesZ;

    plint nx = finestBoundingBox.x1;
    plint ny = finestBoundingBox.y1;
    plint nz = finestBoundingBox.z1;

    util::linearRepartition(0, nx, xTiles, rangesX);
    util::linearRepartition(0, ny, yTiles, rangesY);
    util::linearRepartition(0, nz, zTiles, rangesZ);

    finestDivision.resize(0);
    mpiDistribution.resize(xTiles * yTiles * zTiles);

    for (plint iX = 0; iX < (plint)rangesX.size(); ++iX) {
        for (plint iY = 0; iY < (plint)rangesY.size(); ++iY) {
            for (plint iZ = 0; iZ < (plint)rangesZ.size(); ++iZ) {
                // create a Box3D with coordinates for the new sector
                Box3D box(
                    rangesX[iX].first, rangesX[iX].second, rangesY[iY].first, rangesY[iY].second,
                    rangesZ[iZ].first, rangesZ[iZ].second);
                // put it in the finestDivision
                finestDivision.push_back(box);
            }
        }
    }
}

void ParallellizeByCubes3D::computeDivision()
{
    // compute a simple linear Repartition
    std::vector<std::pair<plint, plint> > rangesX;
    std::vector<std::pair<plint, plint> > rangesY;
    std::vector<std::pair<plint, plint> > rangesZ;

    plint nx = finestBoundingBox.x1;
    plint ny = finestBoundingBox.y1;
    plint nz = finestBoundingBox.z1;

    util::linearRepartition(0, nx, xTiles, rangesX);
    util::linearRepartition(0, ny, yTiles, rangesY);
    util::linearRepartition(0, nz, zTiles, rangesZ);

    finestDivision.resize(0);
    mpiDistribution.resize(xTiles * yTiles * zTiles);

    // move the ranges to fit the cost
    plint cost1, cost2;
    Box3D box1, box2;
    bool notDoneYet = true;

    // xrange
    while (notDoneYet) {
        notDoneYet = false;
        for (plint iX = 0; iX < (plint)rangesX.size() - 1; ++iX) {
            box1 = Box3D(
                rangesX[iX].first, rangesX[iX].second, 0, finestBoundingBox.y1, 0,
                finestBoundingBox.z1);
            box2 = Box3D(
                rangesX[iX + 1].first, rangesX[iX + 1].second, 0, finestBoundingBox.y1, 0,
                finestBoundingBox.z1);
            cost1 = computeCost(originalBlocks, box1);
            cost2 = computeCost(originalBlocks, box2);

            if (cost1 <= 1.05 * cost2 && cost1 >= 0.95 * cost2) {  // allow a tolerance of 10%
                // do nothing
            } else if (cost1 < cost2) {
                ++rangesX[iX].second;
                ++rangesX[iX + 1].first;
                notDoneYet = true;
            } else {
                --rangesX[iX].second;
                --rangesX[iX + 1].first;
                notDoneYet = true;
            }
        }
    }
    // yrange
    notDoneYet = true;
    while (notDoneYet) {
        notDoneYet = false;
        for (plint iY = 0; iY < (plint)rangesY.size() - 1; ++iY) {
            box1 = Box3D(
                0, finestBoundingBox.x1, rangesY[iY].first, rangesY[iY].second, 0,
                finestBoundingBox.z1);
            box2 = Box3D(
                0, finestBoundingBox.x1, rangesY[iY + 1].first, rangesY[iY + 1].second, 0,
                finestBoundingBox.z1);
            cost1 = computeCost(originalBlocks, box1);
            cost2 = computeCost(originalBlocks, box2);

            if (cost1 <= 1.05 * cost2 && cost1 >= 0.95 * cost2) {  // allow a tolerance of 10%
                // do nothing
            } else if (cost1 < cost2) {
                ++rangesY[iY].second;
                ++rangesY[iY + 1].first;
                notDoneYet = true;
            } else {
                --rangesY[iY].second;
                --rangesY[iY + 1].first;
                notDoneYet = true;
            }
        }
    }
    // zrange
    notDoneYet = true;
    while (notDoneYet) {
        notDoneYet = false;
        for (plint iZ = 0; iZ < (plint)rangesZ.size() - 1; ++iZ) {
            box1 = Box3D(
                0, finestBoundingBox.x1, 0, finestBoundingBox.y1, rangesZ[iZ].first,
                rangesZ[iZ].second);
            box2 = Box3D(
                0, finestBoundingBox.x1, 0, finestBoundingBox.y1, rangesZ[iZ + 1].first,
                rangesZ[iZ + 1].second);
            cost1 = computeCost(originalBlocks, box1);
            cost2 = computeCost(originalBlocks, box2);

            if (cost1 <= 1.05 * cost2 && cost1 >= 0.95 * cost2) {  // allow a tolerance of 10%
                // do nothing
            } else if (cost1 < cost2) {
                ++rangesZ[iZ].second;
                ++rangesZ[iZ + 1].first;
                notDoneYet = true;
            } else {
                --rangesZ[iZ].second;
                --rangesZ[iZ + 1].first;
                notDoneYet = true;
            }
        }
    }

    // check that we're not dividing the fineGridBoundaryDynamics
    // (because this causes problems)...
    for (plint iX = 0; iX < (plint)rangesX.size() - 1; ++iX) {
        for (pluint iLevel = 0; iLevel < originalBlocks.size(); ++iLevel) {
            for (pluint iBox = 0; iBox < originalBlocks[iLevel].size(); ++iBox) {
                plint x0 = originalBlocks[iLevel][iBox].x0;
                plint x1 = originalBlocks[iLevel][iBox].x1;
                if ((rangesX[iX].second > x0 - 5 && rangesX[iX].second < x0 + 5)
                    || (rangesX[iX].second > x1 - 5 && rangesX[iX].second < x1 + 5))
                {
                    rangesX[iX].second -= 10;
                    rangesX[iX + 1].first -= 10;
                }
            }
        }
    }
    for (plint iY = 0; iY < (plint)rangesY.size() - 1; ++iY) {
        for (pluint iLevel = 0; iLevel < originalBlocks.size(); ++iLevel) {
            for (pluint iBox = 0; iBox < originalBlocks[iLevel].size(); ++iBox) {
                plint y0 = originalBlocks[iLevel][iBox].y0;
                plint y1 = originalBlocks[iLevel][iBox].y1;
                if ((rangesY[iY].second > y0 - 5 && rangesY[iY].second < y0 + 5)
                    || (rangesY[iY].second > y1 - 5 && rangesY[iY].second < y1 + 5))
                {
                    rangesY[iY].second -= 10;
                    rangesY[iY + 1].first -= 10;
                }
            }
        }
    }
    for (plint iZ = 0; iZ < (plint)rangesZ.size() - 1; ++iZ) {
        for (pluint iLevel = 0; iLevel < originalBlocks.size(); ++iLevel) {
            for (pluint iBox = 0; iBox < originalBlocks[iLevel].size(); ++iBox) {
                plint z0 = originalBlocks[iLevel][iBox].z0;
                plint z1 = originalBlocks[iLevel][iBox].z1;
                if ((rangesZ[iZ].second > z0 - 5 && rangesZ[iZ].second < z0 + 5)
                    || (rangesZ[iZ].second > z1 - 5 && rangesZ[iZ].second < z1 + 5))
                {
                    rangesZ[iZ].second -= 10;
                    rangesZ[iZ + 1].first -= 10;
                }
            }
        }
    }

    // done shifting ranges - add boxes to finestDivision!
    finestDivision.resize(0);
    mpiDistribution.resize(xTiles * yTiles * zTiles);

    for (plint iX = 0; iX < (plint)rangesX.size(); ++iX) {
        for (plint iY = 0; iY < (plint)rangesY.size(); ++iY) {
            for (plint iZ = 0; iZ < (plint)rangesZ.size(); ++iZ) {
                // create a Box3D with coordinates for the new sector
                Box3D box(
                    rangesX[iX].first, rangesX[iX].second, rangesY[iY].first, rangesY[iY].second,
                    rangesZ[iZ].first, rangesZ[iZ].second);
                // put it in the finestDivision
                finestDivision.push_back(box);
            }
        }
    }
}

void ParallellizeByCubes3D::parallelize()
{
    // we must have the same number of blocks as processors
    PLB_PRECONDITION(xTiles * yTiles * zTiles == processorNumber);
    plint totalCost = computeCost(originalBlocks, finestBoundingBox);
    plint idealCostPerProcessor = totalCost / processorNumber;

    pcout << "Total cost of computations = " << totalCost << std::endl;
    pcout << "We are using " << processorNumber << " processors...\n";
    pcout << "Ideal cost per processor = " << idealCostPerProcessor << std::endl;

    std::vector<plint> totalCosts(processorNumber);

    // greedy load balancing part
    //     plint currentProcessor = 0;
    //     bool allAssigned = false;
    //     plint iBlock = 0;
    //     plint maxBlockCost = 0;
    //     //while ( (currentProcessor<processorNumber) && !allAssigned) {
    //     while (!allAssigned) {
    //         plint blockCost = computeCost(originalBlocks, finestDivision[iBlock]);
    //         if (blockCost > maxBlockCost) maxBlockCost = blockCost;
    //         if ( totalCosts[currentProcessor] < idealCostPerProcessor ){
    //             totalCosts[currentProcessor] += blockCost;
    //             mpiDistribution[iBlock] = currentProcessor;
    //             ++iBlock;
    //             currentProcessor++;
    //         }
    //         else {
    //             currentProcessor++;
    //         }
    //         if(currentProcessor>=processorNumber) currentProcessor = 0;
    //         if (iBlock>=(plint)finestDivision.size()){
    //             allAssigned = true;
    //         }
    //     }
    //
    //     if (maxBlockCost > idealCostPerProcessor){
    //         pcout << "There is a problem : maxBlockCost=" << maxBlockCost << " and ideal cost="
    //         << idealCostPerProcessor
    //               << std::endl;
    //     }

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

/**
  IMPLEMENTED BY HELEN MORRISON (2015) (doesn't seem to work properly...)
**/
///* ************* ParallellizeByLoadBalancing3D **************** */
// ParallellizeByLoadBalancing3D::ParallellizeByLoadBalancing3D(
//         std::vector<std::vector<Box3D> > const& originalBlocks_,
//         Box3D finestBoundingBox_)
//     : originalBlocks(originalBlocks_), finestBoundingBox(finestBoundingBox_),
//       processorNumber(global::mpi().getSize())
//{

//}

// void ParallellizeByLoadBalancing3D::parallelize(){
//     // we must have the same number of blocks as processors
//     plint totalCost = computeCost(originalBlocks, finestBoundingBox);
//     plint idealCostPerProcessor = totalCost/processorNumber;

//    pcout << "Total cost of computations = " << totalCost << std::endl;
//    pcout << "We are using " << processorNumber << " processors...\n";
//    pcout << "Ideal cost per processor = " << idealCostPerProcessor << std::endl;

//    std::vector<plint> totalCosts(processorNumber);

//    float deviationFromIdealCost = 1.0; // in percent

//    pcout << "costs of the original blocks: " << std::endl;
//    bool separationNeeded = true;
//    bool neededAtLeastOnce = false;
//    plint numBlocks = 0;
//    std::vector<plint> costPerBlock;
//    for(plint iLevel=0; iLevel<originalBlocks.size(); ++iLevel){
//        for(plint i=0; i<originalBlocks[iLevel].size(); ++i){
//            costPerBlock.push_back(
//                        originalBlocks[iLevel][i].getNx()*
//                        originalBlocks[iLevel][i].getNy()*
//                        originalBlocks[iLevel][i].getNz());
//            pcout << numBlocks << ": " << costPerBlock[numBlocks] << std::endl;
//            if(costPerBlock[numBlocks] > (1+deviationFromIdealCost)*idealCostPerProcessor){
//                neededAtLeastOnce = true;
//            }
//            numBlocks++;
//        }
//    }
//    if(!neededAtLeastOnce) separationNeeded = false;

//    std::vector<std::vector<Box3D> > newBlocks = originalBlocks;
//    std::vector<std::vector<Box3D> > newBlocksTmp;
//    Box3D box1, box2;

//    while(separationNeeded){
//        newBlocksTmp = newBlocks;
//        plint itBlock = 0;
//        for(plint iLevel=0; iLevel<newBlocks.size(); ++iLevel){
//           for(plint i=0; i<newBlocks[iLevel].size(); ++i){
//               if(costPerBlock[itBlock] > (1+deviationFromIdealCost)*idealCostPerProcessor){
//                   splitBlockInTwo(newBlocks[iLevel][i], box1, box2);
//                   newBlocksTmp[iLevel].push_back(box1);
//                   newBlocksTmp[iLevel].push_back(box2);
//               }
//               else{
//                   newBlocksTmp[iLevel].push_back(newBlocks[iLevel][i]);
//               }
//           }
//        }
//        newBlocks = newBlocksTmp;

//        costPerBlock.resize(newBlocks.size());
//        numBlocks = 0;
//        bool neededAtLeastOnce = false;
//        for(plint iLevel=0; iLevel<newBlocks.size(); ++iLevel){
//           for(plint i=0; i<newBlocks[iLevel].size(); ++i){
//               costPerBlock[numBlocks] =
//                       newBlocks[iLevel][i].getNx()*
//                       newBlocks[iLevel][i].getNy()*
//                       newBlocks[iLevel][i].getNz();
//               if(costPerBlock[numBlocks] > (1+deviationFromIdealCost)*idealCostPerProcessor){
//                   neededAtLeastOnce = true;
//               }
//               numBlocks++;
//           }
//        }
//        if(!neededAtLeastOnce) separationNeeded = false;
//    }

//    // check that there are at least as many blocks as processors
//    while(processorNumber > numBlocks){
//        plint maxBlock = 0;
//        plint maxCost = 0;
//        for(plint i=0; i<numBlocks;++i){
//            if(costPerBlock[i] > maxCost){
//                maxCost = costPerBlock[i];
//                maxBlock = i;
//            }
//        }
//        newBlocksTmp = newBlocks;
//        plint itBlock = 0;
//        for(plint iLevel=0; iLevel<newBlocks.size(); ++iLevel){
//           for(plint i=0; i<newBlocks[iLevel].size(); ++i){
//               if(itBlock == maxBlock){
//                   splitBlockInTwo(newBlocks[iLevel][i], box1, box2);
//                   newBlocksTmp[iLevel].push_back(box1);
//                   newBlocksTmp[iLevel].push_back(box2);
//               }
//               else{
//                   newBlocksTmp[iLevel].push_back(newBlocks[iLevel][i]);
//               }
//           }
//        }
//        newBlocks = newBlocksTmp;

//        costPerBlock.resize(newBlocks.size());
//        numBlocks = 0;
//        for(plint iLevel=0; iLevel<newBlocks.size(); ++iLevel){
//           for(plint i=0; i<newBlocks[iLevel].size(); ++i){
//               costPerBlock[numBlocks] =
//                       newBlocks[iLevel][i].getNx()*
//                       newBlocks[iLevel][i].getNy()*
//                       newBlocks[iLevel][i].getNz();
//               numBlocks++;
//           }
//        }
//    }

//    // THE BLOCKS ARE NOW THE RIGHT SIZE - TODO: ALLOCATE THEM APPROPRIATELY TO THE INDIVIDUAL
//    THREADS std::vector<bool> allocated; for(plint i=0; i<numBlocks; ++i){
//        allocated.push_back(false);
//    }
//    plint numBlocksAllocated = 0;

//    //first, allocate the largest blocks until each processor has at least one thing to do!
//    std::vector<plint> threads(costPerBlock.size());
//    for(plint iProc = 0; iProc < processorNumber; ++iProc){
//        plint maxBlock = 0;
//        plint maxCost = 0;
//        for(plint i=0; i<numBlocks;++i){
//            if(costPerBlock[i] > maxCost && !allocated[i]){
//                maxCost = costPerBlock[i];
//                maxBlock = i;
//            }
//        }
//        threads[maxBlock]=iProc;
//        allocated[maxBlock]=true;
//        ++numBlocksAllocated;
//        totalCosts[iProc] = maxCost;
//    }

//    // then, keep adding blocks to processors until idealCost is reached
//    // (starting with the proc that has the least amount to do)
//    while(numBlocksAllocated < numBlocks){
//        plint lowestCost = totalCosts[0];
//        plint whichProc = 0;
//        for(plint iProc=1; iProc<processorNumber; ++iProc){
//            if(totalCosts[iProc] < lowestCost){
//                lowestCost = totalCosts[iProc];
//                whichProc = iProc;
//            }
//        }

//        plint maxBlock = 0;
//        plint maxCost = 0;
//        for(plint i=0; i<numBlocks;++i){
//            if(costPerBlock[i] > maxCost && !allocated[i]){
//                maxCost = costPerBlock[i];
//                maxBlock = i;
//            }
//        }
//        threads[maxBlock]=whichProc;
//        allocated[maxBlock]=true;
//        ++numBlocksAllocated;
//        totalCosts[whichProc] += maxCost;
//    }

////    // update the final mpi distribution
////    this->finalMpiDistribution.resize(recomputedBlocks.size());
////    plint iBlock = 0;
////    for(plint iLevel=0; iLevel<recomputedBlocks.size(); ++iLevel){
////       for(plint i=0; i<recomputedBlocks[iLevel].size(); ++i){
////           this->finalMpiDistribution[iLevel].push_back(threads[iBlock]);
////           ++iBlock;
////       }
////    }

//    //print the total cost per processor:
//    pcout << "Total costs per proc: " << std::endl;
//    for(plint iProc=0; iProc < processorNumber; ++iProc){
//        pcout << iProc << ": " << totalCosts[iProc] << std::endl;
//    }

//    //the finestDivision can now be determined from newBlocks!
//    plint finestLevel= (plint)newBlocks.size()-1;
//    std::vector<Box3D> finestDivision;
//    std::vector<plint> mpiDistribution;
//    std::vector<Box3D> newBoxes;
//    plint iComp=0;
//    for(plint iLevel=0; iLevel<newBlocks.size(); ++iLevel){
//       for(plint i=0; i<newBlocks[iLevel].size(); ++i){
//           Box3D newBox = newBlocks[iLevel][i].multiply(pow(2,finestLevel-iLevel));
//           newBoxes.clear();
//           newBoxes.push_back(newBox);
//           if(finestDivision.size()!=0){
//               for(plint iBox=0; iBox<finestDivision.size(); ++iBox){
//                   std::vector<Box3D> newBoxesTmp;
//                   for(plint iNewBox=0; iNewBox<newBoxes.size(); ++iNewBox){
//                    except(newBoxes[iNewBox], finestDivision[iBox], newBoxesTmp);
//                   }
//                   newBoxes = newBoxesTmp;
//               }
//           }
//           for(plint iBox=0; iBox < newBoxes.size(); ++iBox){
//               finestDivision.push_back(newBoxes[iBox]);
//               mpiDistribution.push_back(threads[iComp]);
//           }
//           ++iComp;
//       }
//    }

//    // set the recomputedBlocks and the finalMpiDistribution appropriately
//    recomputedBlocks.resize(originalBlocks.size());
//    finalMpiDistribution.resize(originalBlocks.size());
//    for (plint iLevel=finestLevel; iLevel>=0; --iLevel) {
//        parallelizeLevel(iLevel, originalBlocks, finestDivision, mpiDistribution);
//        // Adapt the regions to the next-coarser level.
//        for (pluint iRegion=0; iRegion<finestDivision.size(); ++iRegion) {
//            finestDivision[iRegion] = finestDivision[iRegion].divideAndFitSmaller(2);
//        }
//    }

//}

// void ParallellizeByLoadBalancing3D::splitBlockInTwo(Box3D originalBox, Box3D& box1, Box3D& box2){
//     plint nx = originalBox.getNx();
//     plint ny = originalBox.getNy();
//     plint nz = originalBox.getNz();

//    if(nx > ny) {
//        if(nx > nz) {
//            // part in x
//            box1 = Box3D(originalBox.x0, originalBox.x0+nx/2.-1,
//                         originalBox.y0, originalBox.y1,
//                         originalBox.z0, originalBox.z1);
//            box2 = Box3D(originalBox.x0+nx/2., originalBox.x1,
//                         originalBox.y0, originalBox.y1,
//                         originalBox.z0, originalBox.z1);
//        }
//        else {
//            // part in z
//            box1 = Box3D(originalBox.x0, originalBox.x1,
//                         originalBox.y0, originalBox.y1,
//                         originalBox.z0, originalBox.z0+nz/2.-1);
//            box2 = Box3D(originalBox.x0, originalBox.x1,
//                         originalBox.y0, originalBox.y1,
//                         originalBox.z0+nz/2., originalBox.z1);
//        }
//    }
//    else {
//        if (ny > nz) {
//            // part in y
//            box1 = Box3D(originalBox.x0, originalBox.x1,
//                         originalBox.y0, originalBox.y0+ny/2.-1,
//                         originalBox.z0, originalBox.z1);
//            box2 = Box3D(originalBox.x0, originalBox.x1,
//                         originalBox.y0+ny/2., originalBox.y1,
//                         originalBox.z0, originalBox.z1);
//        }
//        else {
//            // part in z
//            box1 = Box3D(originalBox.x0, originalBox.x1,
//                         originalBox.y0, originalBox.y1,
//                         originalBox.z0, originalBox.z0+nz/2.-1);
//            box2 = Box3D(originalBox.x0, originalBox.x1,
//                         originalBox.y0, originalBox.y1,
//                         originalBox.z0+nz/2., originalBox.z1);
//        }
//    }
//}

}  // namespace plb

#endif  // PARALLELIZER_3D_CPP
