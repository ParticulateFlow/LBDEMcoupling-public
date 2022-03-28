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

#include "atomicBlock/reductiveDataProcessorWrapper3D.h"

#include "atomicBlock/atomicBlockOperations3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "core/plbDebug.h"
#include "multiBlock/combinedStatistics.h"

namespace plb {

/* *************** BoxProcessing3D, general case *************************** */

void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional3D &functional, Box3D domain,
    std::vector<AtomicBlock3D *> atomicBlocks)
{
    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveBoxProcessorGenerator3D generator(functional.clone(), domain);
    executeDataProcessor(generator, atomicBlocks);
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** DotProcessing, general case ***************************** */

void applyProcessingFunctional(
    ReductiveDotProcessingFunctional3D &functional, DotList3D const &dotList,
    std::vector<AtomicBlock3D *> atomicBlocks)
{
    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveDotProcessorGenerator3D generator(functional.clone(), dotList);
    executeDataProcessor(generator, atomicBlocks);
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** BoundedReductiveBoxProcessing3D, general case *********** */

void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional3D &functional, Box3D domain,
    std::vector<AtomicBlock3D *> atomicBlocks, plint boundaryWidth)
{
    std::vector<ReductiveBoxProcessorGenerator3D *> generators;
    functional.getGenerators(domain, boundaryWidth, generators);
    std::vector<BlockStatistics const *> individualStatistics(generators.size());
    for (pluint iGen = 0; iGen < generators.size(); ++iGen) {
        executeDataProcessor(*generators[iGen], atomicBlocks);
        individualStatistics[iGen] = &(generators[iGen]->getStatistics());
    }
    SerialCombinedStatistics().combine(individualStatistics, functional.getStatistics());
    for (pluint iGen = 0; iGen < generators.size(); ++iGen) {
        delete generators[iGen];
    }
}

}  // namespace plb
