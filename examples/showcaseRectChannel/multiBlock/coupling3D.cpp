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
 * The 3D coupling -- implementation file.
 */
#include "multiBlock/coupling3D.h"

#include <algorithm>

#include "atomicBlock/atomicBlock3D.h"
#include "core/plbDebug.h"
#include "core/plbProfiler.h"
#include "multiBlock/defaultMultiBlockPolicy3D.h"
#include "multiBlock/multiBlockOperations3D.h"
#include "multiBlock/multiBlockOperations3D.hh"

namespace plb {

/* ********************  Coupling3D ************************** */

Coupling3D::Coupling3D(
    DataProcessorGenerator3D *generator_, std::vector<plint> const &multiBlocks_,
    std::vector<id_t> &allMultiBlocks_) :
    generator(generator_)
{
    for (pluint i = 0; i < multiBlocks_.size(); ++i) {
        id_t id = allMultiBlocks_[multiBlocks_[i]];
        multiBlocks.push_back(id);
    }
    generate();
}

Coupling3D::~Coupling3D()
{
    delete generator;
    clearDataProcessors();
}

Coupling3D::Coupling3D(Coupling3D const &rhs)
{
    generator = rhs.generator->clone();
    multiBlocks = rhs.multiBlocks;
    generate();
}

void Coupling3D::swap(Coupling3D &rhs)
{
    std::swap(generator, rhs.generator);
    multiBlocks.swap(rhs.multiBlocks);
    dataProcessors.swap(rhs.dataProcessors);
}

Coupling3D *Coupling3D::clone() const
{
    return new Coupling3D(*this);
}

void Coupling3D::execute()
{
    for (pluint i = 0; i < dataProcessors.size(); ++i) {
        dataProcessors[i]->process();
    }
}

void Coupling3D::generate()
{
    std::vector<MultiBlock3D *> multiBlockArgs(multiBlocks.size());
    for (pluint i = 0; i < multiBlocks.size(); ++i) {
        id_t id = multiBlocks[i];
        MultiBlock3D *block = multiBlockRegistration3D().find(id);
        PLB_ASSERT(block);
        multiBlockArgs[i] = block;
    }
    MultiProcessing3D<DataProcessorGenerator3D const, DataProcessorGenerator3D> multiProcessing(
        *generator, multiBlockArgs);
    std::vector<DataProcessorGenerator3D *> const &retainedGenerators =
        multiProcessing.getRetainedGenerators();
    std::vector<std::vector<plint> > const &atomicBlockNumbers =
        multiProcessing.getAtomicBlockNumbers();

    clearDataProcessors();
    for (pluint iGenerator = 0; iGenerator < retainedGenerators.size(); ++iGenerator) {
        std::vector<AtomicBlock3D *> extractedAtomicBlocks(multiBlockArgs.size());
        for (pluint iBlock = 0; iBlock < extractedAtomicBlocks.size(); ++iBlock) {
            extractedAtomicBlocks[iBlock] =
                &multiBlockArgs[iBlock]->getComponent(atomicBlockNumbers[iGenerator][iBlock]);
        }
        dataProcessors.push_back(retainedGenerators[iGenerator]->generate(extractedAtomicBlocks));
    }
}

void Coupling3D::clearDataProcessors()
{
    for (pluint i = 0; i < dataProcessors.size(); ++i) {
        delete dataProcessors[i];
    }
    dataProcessors.clear();
}

/* ********************  ReductiveCoupling3D ************************** */

ReductiveCoupling3D::ReductiveCoupling3D(
    ReductiveDataProcessorGenerator3D *generator_, std::vector<plint> const &multiBlocks_,
    std::vector<id_t> &allMultiBlocks_) :
    generator(generator_), combinedStatistics(defaultMultiBlockPolicy3D().getCombinedStatistics())
{
    for (pluint i = 0; i < multiBlocks_.size(); ++i) {
        id_t id = allMultiBlocks_[multiBlocks_[i]];
        multiBlocks.push_back(id);
    }
    generate();
}

ReductiveCoupling3D::~ReductiveCoupling3D()
{
    delete generator;
    delete combinedStatistics;
    clearData();
}

ReductiveCoupling3D::ReductiveCoupling3D(ReductiveCoupling3D const &rhs) :
    generator(rhs.generator->clone()),
    multiBlocks(rhs.multiBlocks),
    combinedStatistics(rhs.combinedStatistics->clone())
{
    generate();
}

void ReductiveCoupling3D::swap(ReductiveCoupling3D &rhs)
{
    std::swap(generator, rhs.generator);
    multiBlocks.swap(rhs.multiBlocks);
    std::swap(combinedStatistics, rhs.combinedStatistics);
    retainedGenerators.swap(rhs.retainedGenerators);
    individualStatistics.swap(rhs.individualStatistics);
    dataProcessors.swap(rhs.dataProcessors);
}

ReductiveCoupling3D *ReductiveCoupling3D::clone() const
{
    return new ReductiveCoupling3D(*this);
}

void ReductiveCoupling3D::execute()
{
    for (pluint i = 0; i < dataProcessors.size(); ++i) {
        dataProcessors[i]->process();
    }
    combinedStatistics->combine(individualStatistics, generator->getStatistics());
}

void ReductiveCoupling3D::generate()
{
    clearData();

    std::vector<MultiBlock3D *> multiBlockArgs(multiBlocks.size());
    for (pluint i = 0; i < multiBlocks.size(); ++i) {
        id_t id = multiBlocks[i];
        MultiBlock3D *block = multiBlockRegistration3D().find(id);
        PLB_ASSERT(block);
        multiBlockArgs[i] = block;
    }
    MultiProcessing3D<ReductiveDataProcessorGenerator3D, ReductiveDataProcessorGenerator3D>
        multiProcessing(*generator, multiBlockArgs);
    std::vector<std::vector<plint> > const &atomicBlockNumbers =
        multiProcessing.getAtomicBlockNumbers();

    retainedGenerators = multiProcessing.releaseRetainedGenerators();
    individualStatistics.resize(retainedGenerators.size());

    for (pluint iGenerator = 0; iGenerator < retainedGenerators.size(); ++iGenerator) {
        std::vector<AtomicBlock3D *> extractedAtomicBlocks(multiBlockArgs.size());
        for (pluint iBlock = 0; iBlock < extractedAtomicBlocks.size(); ++iBlock) {
            extractedAtomicBlocks[iBlock] =
                &multiBlockArgs[iBlock]->getComponent(atomicBlockNumbers[iGenerator][iBlock]);
        }
        // Here, the retainedGenerators are of type ReductiveBoxProcessorGenerator3D, which
        // contain and own a ReductiveBoxProcessingFunctional3D, which owns the BlockStatistics.
        dataProcessors.push_back(retainedGenerators[iGenerator]->generate(extractedAtomicBlocks));
        individualStatistics[iGenerator] = &(retainedGenerators[iGenerator]->getStatistics());
    }
}

BlockStatistics const &ReductiveCoupling3D::statistics() const
{
    return generator->getStatistics();
}

void ReductiveCoupling3D::clearData()
{
    for (pluint i = 0; i < dataProcessors.size(); ++i) {
        delete dataProcessors[i];
    }
    dataProcessors.clear();
    for (pluint iGenerator = 0; iGenerator < retainedGenerators.size(); ++iGenerator) {
        delete retainedGenerators[iGenerator];
    }
    retainedGenerators.clear();
}

/* ********************  CouplingAction3D ************************** */

CouplingAction3D::CouplingAction3D(Coupling3D *coupling_) : coupling(coupling_) { }

CouplingAction3D::CouplingAction3D(CouplingAction3D const &rhs) :
    coupling(rhs.coupling->clone()) { }

CouplingAction3D &CouplingAction3D::operator=(CouplingAction3D const &rhs)
{
    CouplingAction3D(rhs).swap(*this);
    return *this;
}

void CouplingAction3D::swap(CouplingAction3D &rhs)
{
    std::swap(coupling, rhs.coupling);
}

CouplingAction3D::~CouplingAction3D()
{
    delete coupling;
}

CouplingAction3D *CouplingAction3D::clone() const
{
    return new CouplingAction3D(*this);
}

void CouplingAction3D::execute(std::vector<id_t> &allMultiBlocks)
{
    coupling->execute();
}

/* ********************  ReductiveCouplingAction3D ************************** */

ReductiveCouplingAction3D::ReductiveCouplingAction3D(ReductiveCoupling3D *coupling_) :
    coupling(coupling_)
{ }

ReductiveCouplingAction3D::ReductiveCouplingAction3D(ReductiveCouplingAction3D const &rhs) :
    coupling(rhs.coupling->clone())
{ }

ReductiveCouplingAction3D &ReductiveCouplingAction3D::operator=(
    ReductiveCouplingAction3D const &rhs)
{
    ReductiveCouplingAction3D(rhs).swap(*this);
    return *this;
}

void ReductiveCouplingAction3D::swap(ReductiveCouplingAction3D &rhs)
{
    std::swap(coupling, rhs.coupling);
}

ReductiveCouplingAction3D::~ReductiveCouplingAction3D()
{
    delete coupling;
}

ReductiveCouplingAction3D *ReductiveCouplingAction3D::clone() const
{
    return new ReductiveCouplingAction3D(*this);
}

void ReductiveCouplingAction3D::execute(std::vector<id_t> &allMultiBlocks)
{
    coupling->execute();
}

BlockStatistics const &ReductiveCouplingAction3D::statistics() const
{
    return coupling->statistics();
}

/* ********************  CommunicateAction3D ************************** */

CommunicateAction3D::CommunicateAction3D(plint blockId_, modif::ModifT whichData_) :
    blockId(blockId_), whichData(whichData_)
{ }

CommunicateAction3D *CommunicateAction3D::clone() const
{
    return new CommunicateAction3D(*this);
}

void CommunicateAction3D::execute(std::vector<id_t> &allMultiBlocks)
{
    PLB_ASSERT(blockId < (plint)allMultiBlocks.size());
    MultiBlock3D *block = multiBlockRegistration3D().find(allMultiBlocks[blockId]);
    PLB_ASSERT(block);
    block->duplicateOverlaps(whichData);
}

/* ********************  ExecuteInternalProcAction3D ************************** */

ExecuteInternalProcAction3D::ExecuteInternalProcAction3D(plint blockId_, plint level_) :
    blockId(blockId_), level(level_)
{ }

ExecuteInternalProcAction3D *ExecuteInternalProcAction3D::clone() const
{
    return new ExecuteInternalProcAction3D(*this);
}

void ExecuteInternalProcAction3D::execute(std::vector<id_t> &allMultiBlocks)
{
    PLB_ASSERT(blockId < (plint)allMultiBlocks.size());
    MultiBlock3D *block = multiBlockRegistration3D().find(allMultiBlocks[blockId]);
    PLB_ASSERT(block);
    bool communicate = false;
    block->executeInternalProcessors(level, communicate);
}

/* ********************  ExecuteAllInternalProcAction3D ************************** */

ExecuteAllInternalProcAction3D::ExecuteAllInternalProcAction3D(plint blockId_) : blockId(blockId_)
{ }

ExecuteAllInternalProcAction3D *ExecuteAllInternalProcAction3D::clone() const
{
    return new ExecuteAllInternalProcAction3D(*this);
}

void ExecuteAllInternalProcAction3D::execute(std::vector<id_t> &allMultiBlocks)
{
    PLB_ASSERT(blockId < (plint)allMultiBlocks.size());
    MultiBlock3D *block = multiBlockRegistration3D().find(allMultiBlocks[blockId]);
    PLB_ASSERT(block);
    block->executeInternalProcessors();
}

/* ********************  SubActionsAction3D ************************** */

SubActionsAction3D::SubActionsAction3D(Actions3D *actions_) : actions(actions_) { }

SubActionsAction3D::SubActionsAction3D(SubActionsAction3D const &rhs) :
    actions(rhs.actions->clone())
{ }

SubActionsAction3D &SubActionsAction3D::operator=(SubActionsAction3D const &rhs)
{
    SubActionsAction3D(rhs).swap(*this);
    return *this;
}

void SubActionsAction3D::swap(SubActionsAction3D &rhs)
{
    std::swap(actions, rhs.actions);
}

SubActionsAction3D::~SubActionsAction3D()
{
    delete actions;
}

SubActionsAction3D *SubActionsAction3D::clone() const
{
    return new SubActionsAction3D(*this);
}

void SubActionsAction3D::execute(std::vector<id_t> &allMultiBlocks)
{
    actions->execute();
}

Actions3D &SubActionsAction3D::get()
{
    return *actions;
}

/* ********************  EvaluateStatsAction3D ************************** */

EvaluateStatsAction3D::EvaluateStatsAction3D(plint blockId_) : blockId(blockId_) { }

EvaluateStatsAction3D *EvaluateStatsAction3D::clone() const
{
    return new EvaluateStatsAction3D(*this);
}

void EvaluateStatsAction3D::execute(std::vector<id_t> &allMultiBlocks)
{
    PLB_ASSERT(blockId < (plint)allMultiBlocks.size());
    MultiBlock3D *block = multiBlockRegistration3D().find(allMultiBlocks[blockId]);
    PLB_ASSERT(block);
    block->evaluateStatistics();
}

/* ********************  Actions3D ************************** */

Actions3D::Actions3D() { }

Actions3D::~Actions3D()
{
    for (pluint i = 0; i < actions.size(); ++i) {
        delete actions[i];
    }
}

Actions3D::Actions3D(Actions3D const &rhs)
{
    actions.resize(rhs.actions.size());
    for (pluint i = 0; i < actions.size(); ++i) {
        actions[i] = rhs.actions[i]->clone();
    }
    allMultiBlocks = rhs.allMultiBlocks;
}

Actions3D &Actions3D::operator=(Actions3D const &rhs)
{
    Actions3D(rhs).swap(*this);
    return *this;
}

void Actions3D::swap(Actions3D &rhs)
{
    allMultiBlocks.swap(rhs.allMultiBlocks);
    actions.swap(rhs.actions);
}

Actions3D *Actions3D::clone() const
{
    return new Actions3D(*this);
}

plint Actions3D::addBlock(MultiBlock3D &block)
{
    allMultiBlocks.push_back(block.getId());
    return (plint)allMultiBlocks.size() - 1;
}

plint Actions3D::addAction(Action3D *action)
{
    actions.push_back(action);
    return (plint)actions.size() - 1;
}

plint Actions3D::addProcessor(
    BoxProcessingFunctional3D *functional, std::vector<plint> blockNums, Box3D domain)
{
    for (pluint i = 0; i < blockNums.size(); ++i) {
        PLB_ASSERT(blockNums[i] < (plint)allMultiBlocks.size());
    }
    Coupling3D *coupling =
        new Coupling3D(new BoxProcessorGenerator3D(functional, domain), blockNums, allMultiBlocks);
    actions.push_back(new CouplingAction3D(coupling));
    return (plint)actions.size() - 1;
}

plint Actions3D::addProcessor(
    BoxProcessingFunctional3D *functional, std::vector<plint> blockNums,
    std::vector<Box3D> const &domains)
{
    for (pluint i = 0; i < blockNums.size(); ++i) {
        PLB_ASSERT(blockNums[i] < (plint)allMultiBlocks.size());
    }
    Coupling3D *coupling = new Coupling3D(
        new MultiBoxProcessorGenerator3D(functional, domains), blockNums, allMultiBlocks);
    actions.push_back(new CouplingAction3D(coupling));
    return (plint)actions.size() - 1;
}

plint Actions3D::addProcessor(BoxProcessingFunctional3D *functional, plint blockNum1, Box3D domain)
{
    std::vector<plint> blockNums;
    blockNums.push_back(blockNum1);
    return addProcessor(functional, blockNums, domain);
}

plint Actions3D::addProcessor(
    BoxProcessingFunctional3D *functional, plint blockNum1, std::vector<Box3D> const &domains)
{
    std::vector<plint> blockNums;
    blockNums.push_back(blockNum1);
    return addProcessor(functional, blockNums, domains);
}

plint Actions3D::addProcessor(
    BoxProcessingFunctional3D *functional, plint blockNum1, plint blockNum2, Box3D domain)
{
    std::vector<plint> blockNums;
    blockNums.push_back(blockNum1);
    blockNums.push_back(blockNum2);
    return addProcessor(functional, blockNums, domain);
}

plint Actions3D::addProcessor(
    BoxProcessingFunctional3D *functional, plint blockNum1, plint blockNum2,
    std::vector<Box3D> const &domains)
{
    std::vector<plint> blockNums;
    blockNums.push_back(blockNum1);
    blockNums.push_back(blockNum2);
    return addProcessor(functional, blockNums, domains);
}

plint Actions3D::addProcessor(
    BoxProcessingFunctional3D *functional, plint blockNum1, plint blockNum2, plint blockNum3,
    Box3D domain)
{
    std::vector<plint> blockNums;
    blockNums.push_back(blockNum1);
    blockNums.push_back(blockNum2);
    blockNums.push_back(blockNum3);
    return addProcessor(functional, blockNums, domain);
}

plint Actions3D::addProcessor(
    BoxProcessingFunctional3D *functional, plint blockNum1, plint blockNum2, plint blockNum3,
    std::vector<Box3D> const &domains)
{
    std::vector<plint> blockNums;
    blockNums.push_back(blockNum1);
    blockNums.push_back(blockNum2);
    blockNums.push_back(blockNum3);
    return addProcessor(functional, blockNums, domains);
}

plint Actions3D::addProcessor(
    BoxProcessingFunctional3D *functional, plint blockNum1, plint blockNum2, plint blockNum3,
    plint blockNum4, Box3D domain)
{
    std::vector<plint> blockNums;
    blockNums.push_back(blockNum1);
    blockNums.push_back(blockNum2);
    blockNums.push_back(blockNum3);
    blockNums.push_back(blockNum4);
    return addProcessor(functional, blockNums, domain);
}

plint Actions3D::addProcessor(
    BoxProcessingFunctional3D *functional, plint blockNum1, plint blockNum2, plint blockNum3,
    plint blockNum4, plint blockNum5, Box3D domain)
{
    std::vector<plint> blockNums;
    blockNums.push_back(blockNum1);
    blockNums.push_back(blockNum2);
    blockNums.push_back(blockNum3);
    blockNums.push_back(blockNum4);
    blockNums.push_back(blockNum5);
    return addProcessor(functional, blockNums, domain);
}

plint Actions3D::addProcessor(
    ReductiveBoxProcessingFunctional3D *functional, std::vector<plint> blockNums, Box3D domain)
{
    for (pluint i = 0; i < blockNums.size(); ++i) {
        PLB_ASSERT(blockNums[i] < (plint)allMultiBlocks.size());
    }
    ReductiveCoupling3D *coupling = new ReductiveCoupling3D(
        new ReductiveBoxProcessorGenerator3D(functional, domain), blockNums, allMultiBlocks);
    actions.push_back(new ReductiveCouplingAction3D(coupling));
    return (plint)actions.size() - 1;
}

plint Actions3D::addProcessor(
    ReductiveBoxProcessingFunctional3D *functional, plint blockNum1, Box3D domain)
{
    std::vector<plint> blockNums;
    blockNums.push_back(blockNum1);
    return addProcessor(functional, blockNums, domain);
}

plint Actions3D::addProcessor(
    ReductiveBoxProcessingFunctional3D *functional, plint blockNum1, plint blockNum2, Box3D domain)
{
    std::vector<plint> blockNums;
    blockNums.push_back(blockNum1);
    blockNums.push_back(blockNum2);
    return addProcessor(functional, blockNums, domain);
}

plint Actions3D::addProcessor(
    ReductiveBoxProcessingFunctional3D *functional, plint blockNum1, plint blockNum2,
    plint blockNum3, Box3D domain)
{
    std::vector<plint> blockNums;
    blockNums.push_back(blockNum1);
    blockNums.push_back(blockNum2);
    blockNums.push_back(blockNum3);
    return addProcessor(functional, blockNums, domain);
}

plint Actions3D::addActions(Actions3D *subActions)
{
    actions.push_back(new SubActionsAction3D(subActions));
    return (plint)actions.size() - 1;
}

Actions3D &Actions3D::getActions(plint id)
{
    PLB_ASSERT(id < (plint)actions.size());
    SubActionsAction3D *action = dynamic_cast<SubActionsAction3D *>(actions[id]);
    PLB_ASSERT(action);
    return action->get();
}

plint Actions3D::addInternalProcessors(plint blockNum, plint level)
{
    PLB_ASSERT(blockNum < (plint)allMultiBlocks.size());
    actions.push_back(new ExecuteInternalProcAction3D(blockNum, level));
    return (plint)actions.size() - 1;
}

plint Actions3D::addAllInternalProcessorsAndComm(plint blockNum)
{
    PLB_ASSERT(blockNum < (plint)allMultiBlocks.size());
    actions.push_back(new ExecuteAllInternalProcAction3D(blockNum));
    return (plint)actions.size() - 1;
}

plint Actions3D::addCommunication(plint blockNum, modif::ModifT whichData)
{
    PLB_ASSERT(blockNum < (plint)allMultiBlocks.size());
    actions.push_back(new CommunicateAction3D(blockNum, whichData));
    return (plint)actions.size() - 1;
}

plint Actions3D::addEvaluateStats(plint blockNum)
{
    PLB_ASSERT(blockNum < (plint)allMultiBlocks.size());
    actions.push_back(new EvaluateStatsAction3D(blockNum));
    return (plint)actions.size() - 1;
}

plint Actions3D::addSumReductionToVariable(plint actionID, double &result)
{
    return addAction(new CopySumReductionToVariable3D(*this, actionID, result));
}

plint Actions3D::addNOP()
{
    return addAction(new NOPaction());
}

void Actions3D::execute()
{
    for (pluint i = 0; i < actions.size(); ++i) {
        actions[i]->execute(allMultiBlocks);
    }
}

void Actions3D::execute(plint actionId)
{
    PLB_ASSERT(actionId < (plint)actions.size());
    actions[actionId]->execute(allMultiBlocks);
}

void Actions3D::execute(plint actionFrom, plint actionTo)
{
    PLB_ASSERT(actionFrom < (plint)actions.size());
    PLB_ASSERT(actionTo < (plint)actions.size());
    for (plint i = actionFrom; i <= actionTo; ++i) {
        actions[i]->execute(allMultiBlocks);
    }
}

BlockStatistics const &Actions3D::statistics(plint actionId) const
{
    PLB_ASSERT(actionId < (plint)actions.size());
    ReductiveCouplingAction3D *action =
        dynamic_cast<ReductiveCouplingAction3D *>(actions[actionId]);
    PLB_ASSERT(action);
    return action->statistics();
}

CopySumReductionToVariable3D::CopySumReductionToVariable3D(
    Actions3D const &actions_, plint actionID_, double &result_) :
    actions(actions_), actionID(actionID_), result(result_)
{ }

CopySumReductionToVariable3D *CopySumReductionToVariable3D::clone() const
{
    return new CopySumReductionToVariable3D(*this);
}

void CopySumReductionToVariable3D::execute(std::vector<id_t> &allMultiBlocks)
{
    result = actions.statistics(actionID).getSum(0);
}

NOPaction *NOPaction::clone() const
{
    return new NOPaction(*this);
}

void NOPaction::execute(std::vector<id_t> &allMultiBlocks) { }

}  // namespace plb
