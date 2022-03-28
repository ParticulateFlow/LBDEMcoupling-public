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

#include "dataProcessors/metaStuffWrapper3D.h"

#include <algorithm>
#include <random>

#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/dataProcessingFunctional3D.hh"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "atomicBlock/dataProcessorWrapper3D.hh"
#include "dataProcessors/dataInitializerFunctional3D.h"
#include "dataProcessors/dataInitializerFunctional3D.hh"
#include "dataProcessors/dataInitializerWrapper3D.h"
#include "dataProcessors/dataInitializerWrapper3D.hh"
#include "dataProcessors/metaStuffFunctional3D.h"
#include "multiBlock/multiContainerBlock3D.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"
#include "multiBlock/multiDataProcessorWrapper3D.hh"
#include "multiBlock/reductiveMultiDataProcessorWrapper3D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper3D.hh"

namespace plb {

bool allFlagsTrue(MultiBlock3D *multiBlock)
{
    std::vector<MultiBlock3D *> singleBlockVector;
    singleBlockVector.push_back(multiBlock);
    AllFlagsTrueFunctional3D functional;
    applyProcessingFunctional(functional, multiBlock->getBoundingBox(), singleBlockVector);
    return functional.allTrue();
}

void getThreadNum(MultiScalarField3D<int> &threadNum)
{
    applyProcessingFunctional(
        new GetThreadNumFunctional3D(), threadNum.getBoundingBox(), threadNum);
}

void getBlockNum(MultiScalarField3D<plint> &blockNum)
{
    std::vector<plint> const &blocks = blockNum.getLocalInfo().getBlocks();
    for (pluint iBlock = 0; iBlock < blocks.size(); ++iBlock) {
        plint blockId = blocks[iBlock];
        ScalarField3D<plint> &component =
            dynamic_cast<ScalarField3D<plint> &>(blockNum.getComponent(blockId));
        setToConstant(component, component.getBoundingBox(), blockId);
    }
}

void getRandomBlockNum(MultiScalarField3D<plint> &blockNum)
{
    std::vector<plint> const &blocks = blockNum.getLocalInfo().getBlocks();
    plint numIds = blockNum.getSparseBlockStructure().nextIncrementalId();
    std::vector<plint> ids(numIds);
    for (plint i = 0; i < numIds; ++i) {
        ids[i] = i;
    }
    std::random_device rd;
    std::mt19937 rndGenerator(rd());
    std::shuffle(ids.begin(), ids.end(), rndGenerator);
    for (pluint iBlock = 0; iBlock < blocks.size(); ++iBlock) {
        plint blockId = blocks[iBlock];
        ScalarField3D<plint> &component =
            dynamic_cast<ScalarField3D<plint> &>(blockNum.getComponent(blockId));
        setToConstant(component, component.getBoundingBox(), ids[blockId]);
    }
}

}  // namespace plb
