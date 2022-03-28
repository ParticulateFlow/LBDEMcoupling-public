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
 * Operations on the 2D multiblock -- implementation.
 */

#include "multiBlock/multiBlockOperations2D.h"

#include "atomicBlock/atomicBlock2D.h"
#include "atomicBlock/atomicBlockOperations2D.h"
#include "atomicBlock/dataProcessor2D.h"
#include "core/plbDebug.h"
#include "multiBlock/domainManipulation2D.h"
#include "multiBlock/multiBlock2D.h"
#include "multiBlock/multiBlockOperations2D.hh"
#include "multiGrid/multiScale.h"

namespace plb {

void executeDataProcessor(
    DataProcessorGenerator2D const &generator, std::vector<MultiBlock2D *> multiBlocks)
{
    MultiProcessing2D<DataProcessorGenerator2D const, DataProcessorGenerator2D> multiProcessing(
        generator, multiBlocks);
    std::vector<DataProcessorGenerator2D *> const &retainedGenerators =
        multiProcessing.getRetainedGenerators();
    std::vector<std::vector<plint> > const &atomicBlockNumbers =
        multiProcessing.getAtomicBlockNumbers();

    for (pluint iGenerator = 0; iGenerator < retainedGenerators.size(); ++iGenerator) {
        std::vector<AtomicBlock2D *> extractedAtomicBlocks(multiBlocks.size());
        for (pluint iBlock = 0; iBlock < extractedAtomicBlocks.size(); ++iBlock) {
            extractedAtomicBlocks[iBlock] =
                &multiBlocks[iBlock]->getComponent(atomicBlockNumbers[iGenerator][iBlock]);
        }
        // Delegate to the "AtomicBlock version" of executeDataProcessor.
        plb::executeDataProcessor(*retainedGenerators[iGenerator], extractedAtomicBlocks);
    }
    // In the "executeProcessor" version, envelopes are updated right here, because the processor
    //   has already been executed. This behavior is unlike the behavior of the
    //   "addInternalProcessor" version, where envelopes are updated from within the multi-block,
    //   after the execution of internal processors.
    multiProcessing.updateEnvelopesWhereRequired();
}

void executeDataProcessor(DataProcessorGenerator2D const &generator, MultiBlock2D &object)
{
    std::vector<MultiBlock2D *> objects(1);
    objects[0] = &object;
    executeDataProcessor(generator, objects);
}

void executeDataProcessor(
    DataProcessorGenerator2D const &generator, MultiBlock2D &object1, MultiBlock2D &object2)
{
    std::vector<MultiBlock2D *> objects(2);
    objects[0] = &object1;
    objects[1] = &object2;
    executeDataProcessor(generator, objects);
}

void executeDataProcessor(
    ReductiveDataProcessorGenerator2D &generator, std::vector<MultiBlock2D *> multiBlocks)
{
    MultiProcessing2D<ReductiveDataProcessorGenerator2D, ReductiveDataProcessorGenerator2D>
        multiProcessing(generator, multiBlocks);
    std::vector<ReductiveDataProcessorGenerator2D *> const &retainedGenerators =
        multiProcessing.getRetainedGenerators();
    std::vector<std::vector<plint> > const &atomicBlockNumbers =
        multiProcessing.getAtomicBlockNumbers();

    std::vector<BlockStatistics const *> individualStatistics(retainedGenerators.size());
    for (pluint iGenerator = 0; iGenerator < retainedGenerators.size(); ++iGenerator) {
        std::vector<AtomicBlock2D *> extractedAtomicBlocks(multiBlocks.size());
        for (pluint iBlock = 0; iBlock < extractedAtomicBlocks.size(); ++iBlock) {
            extractedAtomicBlocks[iBlock] =
                &multiBlocks[iBlock]->getComponent(atomicBlockNumbers[iGenerator][iBlock]);
        }
        // Delegate to the "AtomicBlock Reductive version" of executeDataProcessor.
        plb::executeDataProcessor(*retainedGenerators[iGenerator], extractedAtomicBlocks);
        individualStatistics[iGenerator] = &(retainedGenerators[iGenerator]->getStatistics());
    }
    multiBlocks[0]->getCombinedStatistics().combine(
        individualStatistics, generator.getStatistics());
    // In the "executeProcessor" version, envelopes are updated right here, because the processor
    //   has already been executed. This behavior is unlike the behavior of the
    //   "addInternalProcessor" version, where envelopes are updated from within the multi-block,
    //   after the execution of internal processors.
    multiProcessing.updateEnvelopesWhereRequired();
}

void executeDataProcessor(ReductiveDataProcessorGenerator2D &generator, MultiBlock2D &object)
{
    std::vector<MultiBlock2D *> objects(1);
    objects[0] = &object;
    executeDataProcessor(generator, objects);
}

void executeDataProcessor(
    ReductiveDataProcessorGenerator2D &generator, MultiBlock2D &object1, MultiBlock2D &object2)
{
    std::vector<MultiBlock2D *> objects(2);
    objects[0] = &object1;
    objects[1] = &object2;
    executeDataProcessor(generator, objects);
}

void addInternalProcessor(
    DataProcessorGenerator2D const &generator, std::vector<MultiBlock2D *> multiBlocks, plint level)
{
    MultiProcessing2D<DataProcessorGenerator2D const, DataProcessorGenerator2D> multiProcessing(
        generator, multiBlocks);
    std::vector<DataProcessorGenerator2D *> const &retainedGenerators =
        multiProcessing.getRetainedGenerators();
    std::vector<std::vector<plint> > const &atomicBlockNumbers =
        multiProcessing.getAtomicBlockNumbers();

    for (pluint iGenerator = 0; iGenerator < retainedGenerators.size(); ++iGenerator) {
        std::vector<AtomicBlock2D *> extractedAtomicBlocks(multiBlocks.size());
        for (pluint iBlock = 0; iBlock < extractedAtomicBlocks.size(); ++iBlock) {
            extractedAtomicBlocks[iBlock] =
                &multiBlocks[iBlock]->getComponent(atomicBlockNumbers[iGenerator][iBlock]);
        }
        // Delegate to the "AtomicBlock version" of addInternal.
        plb::addInternalProcessor(*retainedGenerators[iGenerator], extractedAtomicBlocks, level);
    }
    // Subscribe the processor in the multi-block. This guarantees that the multi-block is aware
    //   of the maximal current processor level, and it instantiates the communication pattern
    //   for an update of envelopes after processor execution.
    std::vector<MultiBlock2D *> updatedMultiBlocks;
    std::vector<modif::ModifT> typeOfModification;
    multiProcessing.multiBlocksWhichRequireUpdate(updatedMultiBlocks, typeOfModification);
    multiBlocks[0]->subscribeProcessor(
        level, updatedMultiBlocks, typeOfModification,
        BlockDomain::usesEnvelope(generator.appliesTo()));
    multiBlocks[0]->storeProcessor(generator, multiBlocks, level);
}

void addInternalProcessor(
    DataProcessorGenerator2D const &generator, MultiBlock2D &object, plint level)
{
    std::vector<MultiBlock2D *> objects(1);
    objects[0] = &object;
    addInternalProcessor(generator, objects, level);
}

void addInternalProcessor(
    DataProcessorGenerator2D const &generator, MultiBlock2D &object1, MultiBlock2D &object2,
    plint level)
{
    std::vector<MultiBlock2D *> objects(2);
    objects[0] = &object1;
    objects[1] = &object2;
    addInternalProcessor(generator, objects, level);
}

}  // namespace plb
