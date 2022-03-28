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

#include "dataProcessors/metaStuffWrapper2D.h"

#include "atomicBlock/dataProcessingFunctional2D.h"
#include "atomicBlock/dataProcessingFunctional2D.hh"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "atomicBlock/dataProcessorWrapper2D.hh"
#include "dataProcessors/dataInitializerFunctional2D.h"
#include "dataProcessors/dataInitializerFunctional2D.hh"
#include "dataProcessors/dataInitializerWrapper2D.h"
#include "dataProcessors/dataInitializerWrapper2D.hh"
#include "dataProcessors/metaStuffFunctional2D.h"
#include "multiBlock/multiContainerBlock2D.h"
#include "multiBlock/multiDataProcessorWrapper2D.h"
#include "multiBlock/multiDataProcessorWrapper2D.hh"
#include "multiBlock/reductiveMultiDataProcessorWrapper2D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper2D.hh"

namespace plb {

bool allFlagsTrue(MultiBlock2D *multiBlock)
{
    std::vector<MultiBlock2D *> singleBlockVector;
    singleBlockVector.push_back(multiBlock);
    AllFlagsTrueFunctional2D functional;
    applyProcessingFunctional(functional, multiBlock->getBoundingBox(), singleBlockVector);
    return functional.allTrue();
}

}  // namespace plb
