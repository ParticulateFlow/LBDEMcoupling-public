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
 * Operations on the 2D multigrid -- implementation
 */

#include "multiGrid/multiGridOperations2D.h"

#include <memory>
#include <vector>

#include "atomicBlock/dataProcessor2D.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "io/parallelIO.h"
#include "multiBlock/multiBlockInfo2D.h"
#include "multiBlock/multiBlockOperations2D.h"
#include "multiGrid/multiGridUtil.h"
#include "multiGrid/multiScale.h"

namespace plb {

/// Execute a data processor over several MultiGrid2D
void executeDataProcessor(
    DataProcessorGenerator2D const &generator, std::vector<MultiGrid2D *> multiGrids,
    plint referenceLevel)
{
    if (multiGrids.empty())
        return;

    for (plint iLevel = 0; iLevel < (plint)multiGrids[0]->getNumLevels(); ++iLevel) {
        std::unique_ptr<DataProcessorGenerator2D> localGenerator(generator.clone());
        plint dxScale = referenceLevel - iLevel;
        plint dtScale = dxScale;  // TODO: here, we assume convective scaling; general case should
                                  // be considered.

        int boxRescaleFactor =
            util::roundToInt(util::twoToThePower(std::abs(referenceLevel - iLevel)));
        if (dxScale > 0)  // if we go to a coarser grid
            localGenerator->divide(boxRescaleFactor);
        else  // otherwise we go to a finer grid
            localGenerator->multiply(boxRescaleFactor);

        localGenerator->setscale(dxScale, dtScale);

        std::vector<MultiBlock2D *> localBlocks(multiGrids.size());
        for (plint iBlock = 0; iBlock < (plint)localBlocks.size(); ++iBlock) {
            localBlocks[iBlock] = &multiGrids[iBlock]->getComponent(iLevel);
        }
        executeDataProcessor(*localGenerator, localBlocks);
    }
}

/// Execute a data processor over a single multiGrid2D
void executeDataProcessor(
    DataProcessorGenerator2D const &generator, MultiGrid2D &object, plint referenceLevel)
{
    std::vector<MultiGrid2D *> blocks;
    blocks.push_back(&object);
    executeDataProcessor(generator, blocks, referenceLevel);
}

/// Execute a data processor over two multiGrid2D
void executeDataProcessor(
    DataProcessorGenerator2D const &generator, MultiGrid2D &object1, MultiGrid2D &object2,
    plint referenceLevel)
{
    std::vector<MultiGrid2D *> blocks;
    blocks.push_back(&object1);
    blocks.push_back(&object2);
    executeDataProcessor(generator, blocks, referenceLevel);
}

void executeDataProcessor(
    ReductiveDataProcessorGenerator2D &generator, std::vector<MultiGrid2D *> multiGrids,
    plint referenceLevel)
{
    if (multiGrids.empty())
        return;
    plint numLevels = (plint)multiGrids[0]->getNumLevels();
    std::vector<int> dimensionsX, dimensionsT;
    generator.getDimensionsX(dimensionsX);
    generator.getDimensionsT(dimensionsT);
    std::vector<ReductiveDataProcessorGenerator2D *> localGenerators(numLevels);
    std::vector<BlockStatistics *> localStatistics(numLevels);
    for (plint iLevel = 0; iLevel < numLevels; ++iLevel) {
        int dxScale = referenceLevel - iLevel;
        int dtScale =
            dxScale;  // TODO: here, we assume convective scaling; general case could be considered.
        localGenerators[iLevel] = generator.clone();

        int boxRescaleFactor =
            util::roundToInt(util::twoToThePower(std::abs(referenceLevel - iLevel)));
        if (dxScale < 0)
            generator.divide(boxRescaleFactor);
        else
            generator.multiply(boxRescaleFactor);

        std::vector<MultiBlock2D *> localBlocks(multiGrids.size());
        for (plint iBlock = 0; iBlock < (plint)localBlocks.size(); ++iBlock) {
            localBlocks[iBlock] = &multiGrids[iBlock]->getComponent(iLevel);
        }
        executeDataProcessor(*localGenerators[iLevel], localBlocks);
        std::vector<double> scales(dimensionsX.size());
        for (pluint iScale = 0; iScale < scales.size(); ++iScale) {
            scales[iScale] =
                scaleToReference(dxScale, dimensionsX[iScale], dtScale, dimensionsT[iScale]);
        }
        localGenerators[iLevel]->getStatistics().rescale(scales);
        localStatistics[iLevel] = &(localGenerators[iLevel]->getStatistics());
    }
    combine(localStatistics, generator.getStatistics());
}

void executeDataProcessor(
    ReductiveDataProcessorGenerator2D &generator, MultiGrid2D &object, plint referenceLevel)
{
    std::vector<MultiGrid2D *> blocks;
    blocks.push_back(&object);
    executeDataProcessor(generator, blocks, referenceLevel);
}

void executeDataProcessor(
    ReductiveDataProcessorGenerator2D &generator, MultiGrid2D &object1, MultiGrid2D &object2,
    plint referenceLevel)
{
    std::vector<MultiGrid2D *> blocks;
    blocks.push_back(&object1);
    blocks.push_back(&object2);
    executeDataProcessor(generator, blocks, referenceLevel);
}

void addInternalProcessor(
    DataProcessorGenerator2D const &generator, std::vector<MultiGrid2D *> multiGrids,
    plint referenceLevel, plint processorLevel)
{
    if (multiGrids.empty())
        return;
    for (plint iLevel = 0; iLevel < (plint)multiGrids[0]->getNumLevels(); ++iLevel) {
        std::unique_ptr<DataProcessorGenerator2D> localGenerator(generator.clone());
        plint dxScale = referenceLevel - iLevel;
        plint dtScale = dxScale;  // TODO: here, we assume convective scaling; general case should
                                  // be considered.

        // TODO: correct for referenceLevel!=0
        localGenerator->setscale(dxScale, dtScale);
        int boxRescaleFactor =
            util::roundToInt(util::twoToThePower(std::abs(referenceLevel - iLevel)));
        if (referenceLevel - iLevel > 0)
            localGenerator->divide(boxRescaleFactor);
        else
            localGenerator->multiply(boxRescaleFactor);

        std::vector<MultiBlock2D *> localBlocks(multiGrids.size());
        for (plint iBlock = 0; iBlock < (plint)localBlocks.size(); ++iBlock) {
            localBlocks[iBlock] = &multiGrids[iBlock]->getComponent(iLevel);
        }
        addInternalProcessor(*localGenerator, localBlocks, processorLevel);
    }
}

void addInternalProcessor(
    DataProcessorGenerator2D const &generator, MultiGrid2D &object, plint referenceLevel,
    plint processorLevel)
{
    std::vector<MultiGrid2D *> blocks;
    blocks.push_back(&object);
    addInternalProcessor(generator, blocks, referenceLevel, processorLevel);
}

void addInternalProcessor(
    DataProcessorGenerator2D const &generator, MultiGrid2D &object1, MultiGrid2D &object2,
    plint referenceLevel, plint processorLevel)
{
    std::vector<MultiGrid2D *> blocks;
    blocks.push_back(&object1);
    blocks.push_back(&object2);
    addInternalProcessor(generator, blocks, referenceLevel, processorLevel);
}

}  // namespace plb
