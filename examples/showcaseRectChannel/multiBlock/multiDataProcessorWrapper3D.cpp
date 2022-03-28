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

#include "multiBlock/multiDataProcessorWrapper3D.h"

#include "atomicBlock/dataProcessor3D.h"
#include "core/plbDebug.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiBlockOperations3D.h"
#include "multiBlock/multiDataField3D.h"

namespace plb {

/* *************** BoxProcessing3D, general case *************************** */

void applyProcessingFunctional(
    BoxProcessingFunctional3D *functional, Box3D domain, std::vector<MultiBlock3D *> multiBlocks)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), multiBlocks);
}

void integrateProcessingFunctional(
    BoxProcessingFunctional3D *functional, Box3D domain, std::vector<MultiBlock3D *> multiBlocks,
    plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), multiBlocks, level);
}

void integrateProcessingFunctional(
    BoxProcessingFunctional3D *functional, Box3D domain, MultiBlock3D &actor,
    std::vector<MultiBlock3D *> multiBlockArgs, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), actor, multiBlockArgs, level);
}

/* *************** DotProcessing, general case ***************************** */

void applyProcessingFunctional(
    DotProcessingFunctional3D *functional, DotList3D const &dotList,
    std::vector<MultiBlock3D *> multiBlocks)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), multiBlocks);
}

void integrateProcessingFunctional(
    DotProcessingFunctional3D *functional, DotList3D const &dotList,
    std::vector<MultiBlock3D *> multiBlocks, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), multiBlocks, level);
}

/* *************** BoundedBoxProcessing3D, general case *************************** */

void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D *functional, Box3D domain,
    std::vector<MultiBlock3D *> multiBlocks, plint boundaryWidth)
{
    std::vector<BoxProcessorGenerator3D *> generators;
    functional->getGenerators(domain, boundaryWidth, generators);
    delete functional;
    for (pluint iGen = 0; iGen < generators.size(); ++iGen) {
        executeDataProcessor(*generators[iGen], multiBlocks);
        delete generators[iGen];
    }
}

void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D *functional, Box3D domain,
    std::vector<MultiBlock3D *> multiBlocks, plint boundaryWidth, plint level)
{
    std::vector<BoxProcessorGenerator3D *> generators;
    functional->getGenerators(domain, boundaryWidth, generators);
    delete functional;
    for (pluint iGen = 0; iGen < generators.size(); ++iGen) {
        addInternalProcessor(*generators[iGen], multiBlocks, level);
        delete generators[iGen];
    }
}

}  // namespace plb
