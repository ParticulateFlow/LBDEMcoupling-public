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

#include "atomicBlock/dataProcessorWrapper3D.h"

#include "atomicBlock/atomicBlockOperations3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "core/plbDebug.h"

namespace plb {

/* *************** BoxProcessing3D, general case *************************** */

void applyProcessingFunctional(
    BoxProcessingFunctional3D *functional, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), atomicBlocks);
}

void integrateProcessingFunctional(
    BoxProcessingFunctional3D *functional, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks,
    plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), atomicBlocks, level);
}

/* *************** DotProcessing, general case ***************************** */

void applyProcessingFunctional(
    DotProcessingFunctional3D *functional, DotList3D const &dotList,
    std::vector<AtomicBlock3D *> atomicBlocks)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), atomicBlocks);
}

void integrateProcessingFunctional(
    DotProcessingFunctional3D *functional, DotList3D const &dotList,
    std::vector<AtomicBlock3D *> atomicBlocks, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), atomicBlocks, level);
}

/* *************** BoundedBoxProcessing3D, general case *************************** */

void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D *functional, Box3D domain,
    std::vector<AtomicBlock3D *> atomicBlocks, plint boundaryWidth)
{
    std::vector<BoxProcessorGenerator3D *> generators;
    functional->getGenerators(domain, boundaryWidth, generators);
    delete functional;
    for (pluint iGen = 0; iGen < generators.size(); ++iGen) {
        executeDataProcessor(*generators[iGen], atomicBlocks);
        delete generators[iGen];
    }
}

void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D *functional, Box3D domain,
    std::vector<AtomicBlock3D *> atomicBlocks, plint boundaryWidth, plint level)
{
    std::vector<BoxProcessorGenerator3D *> generators;
    functional->getGenerators(domain, boundaryWidth, generators);
    delete functional;
    for (pluint iGen = 0; iGen < generators.size(); ++iGen) {
        addInternalProcessor(*generators[iGen], atomicBlocks, level);
        delete generators[iGen];
    }
}

}  // namespace plb
