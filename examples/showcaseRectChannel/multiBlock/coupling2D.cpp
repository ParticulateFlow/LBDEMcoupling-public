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
 * The 2D coupling -- implementation file.
 */
#include "multiBlock/coupling2D.h"

#include <algorithm>

#include "atomicBlock/atomicBlock2D.h"
#include "core/plbDebug.h"
#include "core/plbProfiler.h"
#include "multiBlock/multiBlockOperations2D.h"
#include "multiBlock/multiBlockOperations2D.hh"

namespace plb {

/* ********************  Coupling2D ************************** */

Coupling2D::Coupling2D(
    DataProcessorGenerator2D *generator_, std::vector<plint> const &multiBlocks_,
    std::vector<id_t> &allMultiBlocks) :
    generator(generator_), multiBlocks(multiBlocks_)
{
    generate(allMultiBlocks);
}

Coupling2D::~Coupling2D()
{
    delete generator;
    clearDataProcessors();
}

Coupling2D::Coupling2D(Coupling2D const &rhs)
{
    generator = rhs.generator->clone();
    multiBlocks = rhs.multiBlocks;
    for (pluint i = 0; i < dataProcessors.size(); ++i) {
        dataProcessors.push_back(rhs.dataProcessors[i]->clone());
    }
}

void Coupling2D::swap(Coupling2D &rhs)
{
    std::swap(generator, rhs.generator);
    multiBlocks.swap(rhs.multiBlocks);
    dataProcessors.swap(rhs.dataProcessors);
}

Coupling2D *Coupling2D::clone() const
{
    return new Coupling2D(*this);
}

void Coupling2D::execute()
{
    for (pluint i = 0; i < dataProcessors.size(); ++i) {
        dataProcessors[i]->process();
    }
}

void Coupling2D::generate(std::vector<id_t> &allMultiBlocks)
{
    std::vector<MultiBlock2D *> multiBlockArgs(multiBlocks.size());
    for (pluint i = 0; i < multiBlocks.size(); ++i) {
        id_t id = allMultiBlocks[multiBlocks[i]];
        MultiBlock2D *block = multiBlockRegistration2D().find(id);
        PLB_ASSERT(block);
        multiBlockArgs[i] = block;
    }
    MultiProcessing2D<DataProcessorGenerator2D const, DataProcessorGenerator2D> multiProcessing(
        *generator, multiBlockArgs);
    std::vector<DataProcessorGenerator2D *> const &retainedGenerators =
        multiProcessing.getRetainedGenerators();
    std::vector<std::vector<plint> > const &atomicBlockNumbers =
        multiProcessing.getAtomicBlockNumbers();

    clearDataProcessors();
    for (pluint iGenerator = 0; iGenerator < retainedGenerators.size(); ++iGenerator) {
        std::vector<AtomicBlock2D *> extractedAtomicBlocks(multiBlockArgs.size());
        for (pluint iBlock = 0; iBlock < extractedAtomicBlocks.size(); ++iBlock) {
            extractedAtomicBlocks[iBlock] =
                &multiBlockArgs[iBlock]->getComponent(atomicBlockNumbers[iGenerator][iBlock]);
        }
        dataProcessors.push_back(retainedGenerators[iGenerator]->generate(extractedAtomicBlocks));
    }
}

void Coupling2D::clearDataProcessors()
{
    for (pluint i = 0; i < dataProcessors.size(); ++i) {
        delete dataProcessors[i];
    }
    dataProcessors.clear();
}

/* ********************  CouplingAction2D ************************** */

CouplingAction2D::CouplingAction2D(Coupling2D *coupling_) : coupling(coupling_) { }

CouplingAction2D::CouplingAction2D(CouplingAction2D const &rhs) :
    coupling(rhs.coupling->clone()) { }

CouplingAction2D &CouplingAction2D::operator=(CouplingAction2D const &rhs)
{
    CouplingAction2D(rhs).swap(*this);
    return *this;
}

void CouplingAction2D::swap(CouplingAction2D &rhs)
{
    std::swap(coupling, rhs.coupling);
}

CouplingAction2D::~CouplingAction2D()
{
    delete coupling;
}

CouplingAction2D *CouplingAction2D::clone() const
{
    return new CouplingAction2D(*this);
}

void CouplingAction2D::execute(std::vector<id_t> &allMultiBlocks)
{
    coupling->execute();
}

void CouplingAction2D::regenerate(std::vector<id_t> &allMultiBlocks)
{
    coupling->generate(allMultiBlocks);
}

/* ********************  CommunicateAction2D ************************** */

CommunicateAction2D::CommunicateAction2D(plint blockId_, modif::ModifT whichData_) :
    blockId(blockId_), whichData(whichData_)
{ }

CommunicateAction2D *CommunicateAction2D::clone() const
{
    return new CommunicateAction2D(*this);
}

void CommunicateAction2D::execute(std::vector<id_t> &allMultiBlocks)
{
    PLB_ASSERT(blockId < (plint)allMultiBlocks.size());
    MultiBlock2D *block = multiBlockRegistration2D().find(allMultiBlocks[blockId]);
    PLB_ASSERT(block);
    block->duplicateOverlaps(whichData);
}

void CommunicateAction2D::regenerate(std::vector<id_t> &allMultiBlocks) { }

/* ********************  ExecuteInternalProcAction2D ************************** */

ExecuteInternalProcAction2D::ExecuteInternalProcAction2D(plint blockId_, plint level_) :
    blockId(blockId_), level(level_)
{ }

ExecuteInternalProcAction2D *ExecuteInternalProcAction2D::clone() const
{
    return new ExecuteInternalProcAction2D(*this);
}

void ExecuteInternalProcAction2D::execute(std::vector<id_t> &allMultiBlocks)
{
    PLB_ASSERT(blockId < (plint)allMultiBlocks.size());
    MultiBlock2D *block = multiBlockRegistration2D().find(allMultiBlocks[blockId]);
    PLB_ASSERT(block);
    bool communicate = false;
    block->executeInternalProcessors(level, communicate);
}

void ExecuteInternalProcAction2D::regenerate(std::vector<id_t> &allMultiBlocks) { }

/* ********************  EvaluateStatsAction2D ************************** */

EvaluateStatsAction2D::EvaluateStatsAction2D(plint blockId_) : blockId(blockId_) { }

EvaluateStatsAction2D *EvaluateStatsAction2D::clone() const
{
    return new EvaluateStatsAction2D(*this);
}

void EvaluateStatsAction2D::execute(std::vector<id_t> &allMultiBlocks)
{
    PLB_ASSERT(blockId < (plint)allMultiBlocks.size());
    MultiBlock2D *block = multiBlockRegistration2D().find(allMultiBlocks[blockId]);
    PLB_ASSERT(block);
    block->evaluateStatistics();
}

void EvaluateStatsAction2D::regenerate(std::vector<id_t> &allMultiBlocks) { }

/* ********************  Actions2D ************************** */

Actions2D::Actions2D() { }

Actions2D::~Actions2D()
{
    for (pluint i = 0; i < actions.size(); ++i) {
        delete actions[i];
    }
}

Actions2D::Actions2D(Actions2D const &rhs)
{
    actions.resize(rhs.actions.size());
    for (pluint i = 0; i < actions.size(); ++i) {
        actions[i] = rhs.actions[i]->clone();
    }
    allMultiBlocks = rhs.allMultiBlocks;
}

Actions2D &Actions2D::operator=(Actions2D const &rhs)
{
    Actions2D(rhs).swap(*this);
    return *this;
}

void Actions2D::swap(Actions2D &rhs)
{
    allMultiBlocks.swap(rhs.allMultiBlocks);
    actions.swap(rhs.actions);
}

Actions2D *Actions2D::clone() const
{
    return new Actions2D(*this);
}

plint Actions2D::addBlock(MultiBlock2D &block)
{
    allMultiBlocks.push_back(block.getId());
    return (plint)allMultiBlocks.size() - 1;
}

void Actions2D::replaceBlock(plint id, MultiBlock2D &block)
{
    PLB_ASSERT(id < (plint)allMultiBlocks.size());
    allMultiBlocks[id] = block.getId();
    for (pluint i = 0; i < actions.size(); ++i) {
        actions[i]->regenerate(allMultiBlocks);
    }
}

plint Actions2D::addProcessor(
    BoxProcessingFunctional2D *functional, std::vector<plint> blockNums, Box2D domain)
{
    for (pluint i = 0; i < blockNums.size(); ++i) {
        PLB_ASSERT(blockNums[i] < (plint)allMultiBlocks.size());
    }
    Coupling2D *coupling =
        new Coupling2D(new BoxProcessorGenerator2D(functional, domain), blockNums, allMultiBlocks);
    actions.push_back(new CouplingAction2D(coupling));
    return (plint)actions.size() - 1;
}

plint Actions2D::addProcessor(BoxProcessingFunctional2D *functional, plint blockNum1, Box2D domain)
{
    std::vector<plint> blockNums;
    blockNums.push_back(blockNum1);
    return addProcessor(functional, blockNums, domain);
}

plint Actions2D::addProcessor(
    BoxProcessingFunctional2D *functional, plint blockNum1, plint blockNum2, Box2D domain)
{
    std::vector<plint> blockNums;
    blockNums.push_back(blockNum1);
    blockNums.push_back(blockNum2);
    return addProcessor(functional, blockNums, domain);
}

plint Actions2D::addProcessor(
    BoxProcessingFunctional2D *functional, plint blockNum1, plint blockNum2, plint blockNum3,
    Box2D domain)
{
    std::vector<plint> blockNums;
    blockNums.push_back(blockNum1);
    blockNums.push_back(blockNum2);
    blockNums.push_back(blockNum3);
    return addProcessor(functional, blockNums, domain);
}

plint Actions2D::addProcessor(
    BoxProcessingFunctional2D *functional, plint blockNum1, plint blockNum2, plint blockNum3,
    plint blockNum4, Box2D domain)
{
    std::vector<plint> blockNums;
    blockNums.push_back(blockNum1);
    blockNums.push_back(blockNum2);
    blockNums.push_back(blockNum3);
    blockNums.push_back(blockNum4);
    return addProcessor(functional, blockNums, domain);
}

plint Actions2D::addProcessor(
    BoxProcessingFunctional2D *functional, plint blockNum1, plint blockNum2, plint blockNum3,
    plint blockNum4, plint blockNum5, Box2D domain)
{
    std::vector<plint> blockNums;
    blockNums.push_back(blockNum1);
    blockNums.push_back(blockNum2);
    blockNums.push_back(blockNum3);
    blockNums.push_back(blockNum4);
    blockNums.push_back(blockNum5);
    return addProcessor(functional, blockNums, domain);
}

plint Actions2D::addInternalProcessors(plint blockNum, plint level)
{
    PLB_ASSERT(blockNum < (plint)allMultiBlocks.size());
    actions.push_back(new ExecuteInternalProcAction2D(blockNum, level));
    return (plint)actions.size() - 1;
}

plint Actions2D::addCommunication(plint blockNum, modif::ModifT whichData)
{
    PLB_ASSERT(blockNum < (plint)allMultiBlocks.size());
    actions.push_back(new CommunicateAction2D(blockNum, whichData));
    return (plint)actions.size() - 1;
}

plint Actions2D::addEvaluateStats(plint blockNum)
{
    PLB_ASSERT(blockNum < (plint)allMultiBlocks.size());
    actions.push_back(new EvaluateStatsAction2D(blockNum));
    return (plint)actions.size() - 1;
}

void Actions2D::execute()
{
    for (pluint i = 0; i < actions.size(); ++i) {
        actions[i]->execute(allMultiBlocks);
    }
}

void Actions2D::execute(plint actionId)
{
    PLB_ASSERT(actionId < (plint)actions.size());
    actions[actionId]->execute(allMultiBlocks);
}

void Actions2D::execute(plint actionFrom, plint actionTo)
{
    PLB_ASSERT(actionFrom < (plint)actions.size());
    PLB_ASSERT(actionTo < (plint)actions.size());
    for (plint i = actionFrom; i <= actionTo; ++i) {
        actions[i]->execute(allMultiBlocks);
    }
}

}  // namespace plb
