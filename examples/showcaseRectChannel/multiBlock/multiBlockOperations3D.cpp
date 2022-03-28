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
 * Operations on the 3D multiblock -- implementation.
 */

#include "multiBlock/multiBlockOperations3D.h"

#include "core/plbDebug.h"
#include "multiBlock/multiBlockOperations3D.hh"

namespace plb {

void executeDataProcessor(
    DataProcessorGenerator3D const &generator, std::vector<MultiBlock3D *> multiBlocks)
{
    MultiProcessing3D<DataProcessorGenerator3D const, DataProcessorGenerator3D> multiProcessing(
        generator, multiBlocks);
    std::vector<DataProcessorGenerator3D *> const &retainedGenerators =
        multiProcessing.getRetainedGenerators();
    std::vector<std::vector<plint> > const &atomicBlockNumbers =
        multiProcessing.getAtomicBlockNumbers();

    for (pluint iGenerator = 0; iGenerator < retainedGenerators.size(); ++iGenerator) {
        std::vector<AtomicBlock3D *> extractedAtomicBlocks(multiBlocks.size());
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

void executeDataProcessor(DataProcessorGenerator3D const &generator, MultiBlock3D &object)
{
    std::vector<MultiBlock3D *> objects(1);
    objects[0] = &object;
    executeDataProcessor(generator, objects);
}

void executeDataProcessor(
    DataProcessorGenerator3D const &generator, MultiBlock3D &object1, MultiBlock3D &object2)
{
    std::vector<MultiBlock3D *> objects(2);
    objects[0] = &object1;
    objects[1] = &object2;
    executeDataProcessor(generator, objects);
}

void executeDataProcessor(
    ReductiveDataProcessorGenerator3D &generator, std::vector<MultiBlock3D *> multiBlocks)
{
    MultiProcessing3D<ReductiveDataProcessorGenerator3D, ReductiveDataProcessorGenerator3D>
        multiProcessing(generator, multiBlocks);
    std::vector<ReductiveDataProcessorGenerator3D *> const &retainedGenerators =
        multiProcessing.getRetainedGenerators();
    std::vector<std::vector<plint> > const &atomicBlockNumbers =
        multiProcessing.getAtomicBlockNumbers();

    std::vector<BlockStatistics const *> individualStatistics(retainedGenerators.size());
    for (pluint iGenerator = 0; iGenerator < retainedGenerators.size(); ++iGenerator) {
        std::vector<AtomicBlock3D *> extractedAtomicBlocks(multiBlocks.size());
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

void executeDataProcessor(ReductiveDataProcessorGenerator3D &generator, MultiBlock3D &object)
{
    std::vector<MultiBlock3D *> objects(1);
    objects[0] = &object;
    executeDataProcessor(generator, objects);
}

void executeDataProcessor(
    ReductiveDataProcessorGenerator3D &generator, MultiBlock3D &object1, MultiBlock3D &object2)
{
    std::vector<MultiBlock3D *> objects(2);
    objects[0] = &object1;
    objects[1] = &object2;
    executeDataProcessor(generator, objects);
}

void addInternalProcessor(
    DataProcessorGenerator3D const &generator, MultiBlock3D &actor,
    std::vector<MultiBlock3D *> multiBlockArgs, plint level)
{
    MultiProcessing3D<DataProcessorGenerator3D const, DataProcessorGenerator3D> multiProcessing(
        generator, multiBlockArgs);
    std::vector<DataProcessorGenerator3D *> const &retainedGenerators =
        multiProcessing.getRetainedGenerators();
    std::vector<std::vector<plint> > const &atomicBlockNumbers =
        multiProcessing.getAtomicBlockNumbers();

    for (pluint iGenerator = 0; iGenerator < retainedGenerators.size(); ++iGenerator) {
        std::vector<AtomicBlock3D *> extractedAtomicBlocks(multiBlockArgs.size());
        for (pluint iBlock = 0; iBlock < extractedAtomicBlocks.size(); ++iBlock) {
            extractedAtomicBlocks[iBlock] =
                &multiBlockArgs[iBlock]->getComponent(atomicBlockNumbers[iGenerator][iBlock]);
        }
        // It is assumed that the actor has the same distribution as block 0.
        PLB_ASSERT(!atomicBlockNumbers[iGenerator].empty());
        AtomicBlock3D &atomicActor = actor.getComponent(atomicBlockNumbers[iGenerator][0]);
        // Delegate to the "AtomicBlock version" of addInternal.
        plb::addInternalProcessor(
            *retainedGenerators[iGenerator], atomicActor, extractedAtomicBlocks, level);
    }
    // Subscribe the processor in the multi-block. This guarantees that the multi-block is aware
    //   of the maximal current processor level, and it instantiates the communication pattern
    //   for an update of envelopes after processor execution.
    std::vector<MultiBlock3D *> updatedMultiBlocks;
    std::vector<modif::ModifT> typeOfModification;
    multiProcessing.multiBlocksWhichRequireUpdate(updatedMultiBlocks, typeOfModification);
    actor.subscribeProcessor(
        level, updatedMultiBlocks, typeOfModification,
        BlockDomain::usesEnvelope(generator.appliesTo()));
    actor.storeProcessor(generator, multiBlockArgs, level);
}

void addInternalProcessor(
    DataProcessorGenerator3D const &generator, std::vector<MultiBlock3D *> multiBlocks, plint level)
{
    if (multiBlocks.size() > 0) {
        addInternalProcessor(generator, *multiBlocks[0], multiBlocks, level);
    }
}

void addInternalProcessor(
    DataProcessorGenerator3D const &generator, MultiBlock3D &object, plint level)
{
    std::vector<MultiBlock3D *> objects(1);
    objects[0] = &object;
    addInternalProcessor(generator, objects, level);
}

void addInternalProcessor(
    DataProcessorGenerator3D const &generator, MultiBlock3D &object1, MultiBlock3D &object2,
    plint level)
{
    std::vector<MultiBlock3D *> objects(2);
    objects[0] = &object1;
    objects[1] = &object2;
    addInternalProcessor(generator, objects, level);
}

}  // namespace plb
