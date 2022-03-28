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
 * The 2D multiblock -- implementation file.
 */
#include "multiBlock/multiBlock2D.h"

#include <algorithm>
#include <cmath>

#include "atomicBlock/atomicBlock2D.h"
#include "core/plbDebug.h"
#include "core/plbProfiler.h"
#include "multiBlock/defaultMultiBlockPolicy2D.h"
#include "multiBlock/multiBlockOperations2D.h"
#include "multiBlock/multiBlockSerializer2D.h"

namespace plb {

/* *************** Class MultiStatSubscriber2D ****************************** */

MultiStatSubscriber2D::MultiStatSubscriber2D(MultiBlock2D &multiBlock_) :
    multiBlock(multiBlock_) { }

plint MultiStatSubscriber2D::subscribeAverage()
{
    std::vector<plint> const &blocks = multiBlock.getLocalInfo().getBlocks();
    for (pluint iBlock = 0; iBlock < blocks.size(); ++iBlock) {
        plint blockId = blocks[iBlock];
        multiBlock.getComponent(blockId).internalStatSubscription().subscribeAverage();
    }
    return multiBlock.getInternalStatistics().subscribeAverage();
}

plint MultiStatSubscriber2D::subscribeSum()
{
    std::vector<plint> const &blocks = multiBlock.getLocalInfo().getBlocks();
    for (pluint iBlock = 0; iBlock < blocks.size(); ++iBlock) {
        plint blockId = blocks[iBlock];
        multiBlock.getComponent(blockId).internalStatSubscription().subscribeSum();
    }
    return multiBlock.getInternalStatistics().subscribeSum();
}

plint MultiStatSubscriber2D::subscribeMax()
{
    std::vector<plint> const &blocks = multiBlock.getLocalInfo().getBlocks();
    for (pluint iBlock = 0; iBlock < blocks.size(); ++iBlock) {
        plint blockId = blocks[iBlock];
        multiBlock.getComponent(blockId).internalStatSubscription().subscribeMax();
    }
    return multiBlock.getInternalStatistics().subscribeMax();
}

plint MultiStatSubscriber2D::subscribeIntSum()
{
    std::vector<plint> const &blocks = multiBlock.getLocalInfo().getBlocks();
    for (pluint iBlock = 0; iBlock < blocks.size(); ++iBlock) {
        plint blockId = blocks[iBlock];
        multiBlock.getComponent(blockId).internalStatSubscription().subscribeIntSum();
    }
    return multiBlock.getInternalStatistics().subscribeIntSum();
}

/* *************** Class PeriodicitySwitch2D ******************************** */

PeriodicitySwitch2D::PeriodicitySwitch2D(MultiBlock2D &block_) :
    // Default all boundaries to non-periodic. Note: it is not possible to send a
    // signal "signalPeriodicity" to "block" at this point, because "block" might
    // not be fully constructed, and virtual function calls not available.
    periodicity(false, false),
    block(block_)
{ }

PeriodicitySwitch2D::PeriodicitySwitch2D(MultiBlock2D &block_, PeriodicitySwitch2D const &rhs) :
    periodicity(rhs.periodicity), block(block_)
{ }

PeriodicitySwitch2D &PeriodicitySwitch2D::operator=(PeriodicitySwitch2D const &rhs)
{
    // Don't modify the MultiBlock2D reference, as it is immutable.
    periodicity = rhs.periodicity;
    return *this;
}

void PeriodicitySwitch2D::toggle(plint direction, bool periodic)
{
    PLB_PRECONDITION(direction == 0 || direction == 1);
    periodicity[direction] = periodic;
    // Make sure periodic envelope is filled with data, else all
    //   subsequent operations are wrong.
    block.signalPeriodicity();
    block.duplicateOverlaps(modif::dataStructure);
}

void PeriodicitySwitch2D::toggleAll(bool periodic)
{
    periodicity[0] = periodic;
    periodicity[1] = periodic;
    // Make sure periodic envelope is filled with data, else all
    //   subsequent operations are wrong.
    block.signalPeriodicity();
    block.duplicateOverlaps(modif::dataStructure);
}

bool PeriodicitySwitch2D::get(plint direction) const
{
    PLB_PRECONDITION(direction == 0 || direction == 1);

    return periodicity[direction];
}

bool PeriodicitySwitch2D::get(plint normalX, plint normalY) const
{
    bool testX = normalX != 0;
    bool testY = normalY != 0;
    return ((!testX || (testX && periodicity[0])) && (!testY || (testY && periodicity[1])));
}

Box2D PeriodicitySwitch2D::getPeriodicEnvelope(Box2D const &bulk, plint envelopeWidth) const
{
    Box2D envelope(bulk);
    if (get(0)) {  // If periodic in x, extend bulk in x-direction
        envelope.x0 -= envelopeWidth;
        envelope.x1 += envelopeWidth;
    }
    if (get(1)) {  // If periodic in y, extend bulk in y-direction
        envelope.y0 -= envelopeWidth;
        envelope.y1 += envelopeWidth;
    }
    return envelope;
}

/* *************** Class MultiBlock2D::ProcessorStorage2D ******************* */

MultiBlock2D::ProcessorStorage2D::ProcessorStorage2D(
    DataProcessorGenerator2D const &generator_, std::vector<MultiBlock2D *> const &multiBlocks_,
    plint level_) :
    generator(generator_.clone()), multiBlockIds(multiBlocks_.size()), level(level_)
{
    for (pluint iBlock = 0; iBlock < multiBlockIds.size(); ++iBlock) {
        multiBlockIds[iBlock] = multiBlocks_[iBlock]->getId();
    }
}

MultiBlock2D::ProcessorStorage2D::~ProcessorStorage2D()
{
    delete generator;
}

MultiBlock2D::ProcessorStorage2D::ProcessorStorage2D(MultiBlock2D::ProcessorStorage2D const &rhs) :
    generator(rhs.generator->clone()), multiBlockIds(rhs.multiBlockIds), level(rhs.level)
{ }

MultiBlock2D::ProcessorStorage2D &MultiBlock2D::ProcessorStorage2D::operator=(
    MultiBlock2D::ProcessorStorage2D const &rhs)
{
    MultiBlock2D::ProcessorStorage2D(rhs).swap(*this);
    return *this;
}

void MultiBlock2D::ProcessorStorage2D::swap(MultiBlock2D::ProcessorStorage2D &rhs)
{
    std::swap(generator, rhs.generator);
    multiBlockIds.swap(rhs.multiBlockIds);
    std::swap(level, rhs.level);
}

MultiBlock2D::ProcessorStorage2D *MultiBlock2D::ProcessorStorage2D::clone() const
{
    return new MultiBlock2D::ProcessorStorage2D(*this);
}

DataProcessorGenerator2D const &MultiBlock2D::ProcessorStorage2D::getGenerator() const
{
    return *generator;
}

std::vector<id_t> const &MultiBlock2D::ProcessorStorage2D::getMultiBlockIds() const
{
    return multiBlockIds;
}

std::vector<MultiBlock2D *> MultiBlock2D::ProcessorStorage2D::getMultiBlocks() const
{
    std::vector<MultiBlock2D *> multiBlocks(multiBlockIds.size());
    for (pluint iBlock = 0; iBlock < multiBlockIds.size(); ++iBlock) {
        multiBlocks[iBlock] = multiBlockRegistration2D().find(multiBlockIds[iBlock]);
        PLB_ASSERT(multiBlocks[iBlock]);
    }
    return multiBlocks;
}

plint MultiBlock2D::ProcessorStorage2D::getLevel() const
{
    return level;
}

void MultiBlock2D::ProcessorStorage2D::replace(id_t oldBlock, id_t newBlock)
{
    for (pluint iBlock = 0; iBlock < multiBlockIds.size(); ++iBlock) {
        if (multiBlockIds[iBlock] == oldBlock) {
            multiBlockIds[iBlock] = newBlock;
        }
    }
}

/* *************** Class MultiBlock2D *************************************** */

MultiBlock2D::MultiBlock2D(
    MultiBlockManagement2D const &multiBlockManagement_, BlockCommunicator2D *blockCommunicator_,
    CombinedStatistics *combinedStatistics_) :
    multiBlockManagement(multiBlockManagement_),
    maxProcessorLevel(-1),
    blockCommunicator(blockCommunicator_),
    internalStatistics(),
    combinedStatistics(combinedStatistics_),
    statSubscriber(*this),
    statisticsOn(true),
    periodicitySwitch(*this),
    internalModifT(modif::staticVariables)
{
    id = multiBlockRegistration2D().announce(*this);
}

MultiBlock2D::MultiBlock2D(plint nx, plint ny, plint envelopeWidth) :
    multiBlockManagement(defaultMultiBlockPolicy2D().getMultiBlockManagement(
        Box2D(0, nx - 1, 0, ny - 1), envelopeWidth)),
    maxProcessorLevel(-1),
    blockCommunicator(defaultMultiBlockPolicy2D().getBlockCommunicator()),
    internalStatistics(),
    combinedStatistics(defaultMultiBlockPolicy2D().getCombinedStatistics()),
    statSubscriber(*this),
    statisticsOn(true),
    periodicitySwitch(*this),
    internalModifT(modif::staticVariables)
{
    id = multiBlockRegistration2D().announce(*this);
}

MultiBlock2D::MultiBlock2D(MultiBlock2D const &rhs) :
    multiBlockManagement(rhs.multiBlockManagement),
    multiBlocksChangedByManualProcessors(rhs.multiBlocksChangedByManualProcessors),
    multiBlocksChangedByAutomaticProcessors(rhs.multiBlocksChangedByAutomaticProcessors),
    maxProcessorLevel(rhs.maxProcessorLevel),
    storedProcessors(rhs.storedProcessors),
    blockCommunicator(rhs.blockCommunicator->clone()),
    internalStatistics(rhs.internalStatistics),
    combinedStatistics(rhs.combinedStatistics->clone()),
    statSubscriber(*this),
    statisticsOn(rhs.statisticsOn),
    periodicitySwitch(*this, rhs.periodicitySwitch),
    internalModifT(rhs.internalModifT)
{
    id = multiBlockRegistration2D().announce(*this);
}

MultiBlock2D::MultiBlock2D(MultiBlock2D const &rhs, Box2D subDomain, bool crop) :
    multiBlockManagement(intersect(rhs.getMultiBlockManagement(), subDomain, crop)),
    maxProcessorLevel(-1),
    storedProcessors(rhs.storedProcessors),
    blockCommunicator(rhs.blockCommunicator->clone()),
    internalStatistics(),
    combinedStatistics(rhs.combinedStatistics->clone()),
    statSubscriber(*this),
    statisticsOn(true),
    periodicitySwitch(*this),
    internalModifT(rhs.internalModifT)
{
    id = multiBlockRegistration2D().announce(*this);
}

void MultiBlock2D::swap(MultiBlock2D &rhs)
{
    multiBlockManagement.swap(rhs.multiBlockManagement);
    multiBlocksChangedByManualProcessors.swap(rhs.multiBlocksChangedByManualProcessors);
    multiBlocksChangedByAutomaticProcessors.swap(rhs.multiBlocksChangedByAutomaticProcessors);
    std::swap(maxProcessorLevel, rhs.maxProcessorLevel);
    storedProcessors.swap(rhs.storedProcessors);
    std::swap(blockCommunicator, rhs.blockCommunicator);
    std::swap(internalStatistics, rhs.internalStatistics);
    std::swap(combinedStatistics, rhs.combinedStatistics);
    std::swap(statisticsOn, rhs.statisticsOn);
    std::swap(periodicitySwitch, rhs.periodicitySwitch);
    std::swap(internalModifT, rhs.internalModifT);
}

MultiBlock2D::~MultiBlock2D()
{
    delete blockCommunicator;
    delete combinedStatistics;
    multiBlockRegistration2D().release(*this);
}

id_t MultiBlock2D::getId() const
{
    return id;
}

void MultiBlock2D::initialize()
{
    executeInternalProcessors();
}

Box2D MultiBlock2D::getBoundingBox() const
{
    return multiBlockManagement.getBoundingBox();
}

plint MultiBlock2D::getNx() const
{
    return multiBlockManagement.getBoundingBox().getNx();
}

plint MultiBlock2D::getNy() const
{
    return multiBlockManagement.getBoundingBox().getNy();
}

MultiBlockManagement2D const &MultiBlock2D::getMultiBlockManagement() const
{
    return multiBlockManagement;
}

LocalMultiBlockInfo2D const &MultiBlock2D::getLocalInfo() const
{
    return multiBlockManagement.getLocalInfo();
}

SparseBlockStructure2D const &MultiBlock2D::getSparseBlockStructure() const
{
    return multiBlockManagement.getSparseBlockStructure();
}

BlockStatistics &MultiBlock2D::getInternalStatistics()
{
    return internalStatistics;
}

BlockStatistics const &MultiBlock2D::getInternalStatistics() const
{
    return internalStatistics;
}

CombinedStatistics const &MultiBlock2D::getCombinedStatistics() const
{
    return *combinedStatistics;
}

StatSubscriber &MultiBlock2D::internalStatSubscription()
{
    return statSubscriber;
}

BlockCommunicator2D const &MultiBlock2D::getBlockCommunicator() const
{
    return *blockCommunicator;
}

void MultiBlock2D::duplicateOverlaps(modif::ModifT whichData)
{
    this->getBlockCommunicator().duplicateOverlaps(*this, whichData);
}

void MultiBlock2D::signalPeriodicity()
{
    getBlockCommunicator().signalPeriodicity();
}

DataSerializer *MultiBlock2D::getBlockSerializer(
    Box2D const &domain, IndexOrdering::OrderingT ordering) const
{
    return new MultiBlockFastSerializer2D(*this, domain, ordering);
}

DataUnSerializer *MultiBlock2D::getBlockUnSerializer(
    Box2D const &domain, IndexOrdering::OrderingT ordering)
{
    return new MultiBlockFastUnSerializer2D(*this, domain, ordering);
}

void MultiBlock2D::resetFlags()
{
    std::vector<plint> const &blocks = getLocalInfo().getBlocks();
    for (pluint iBlock = 0; iBlock < blocks.size(); ++iBlock) {
        plint blockId = blocks[iBlock];
        getComponent(blockId).setFlag(false);
    }
}

void MultiBlock2D::getDynamicsDict(Box2D domain, std::map<std::string, int> &dict)
{
    return dict.clear();
}

void MultiBlock2D::evaluateStatistics()
{
    std::vector<plint> const &blocks = getLocalInfo().getBlocks();
    for (pluint iBlock = 0; iBlock < blocks.size(); ++iBlock) {
        plint blockId = blocks[iBlock];
        getComponent(blockId).evaluateStatistics();
    }
    if (isInternalStatisticsOn())
        reduceStatistics();
}

void MultiBlock2D::reduceStatistics()
{
    std::vector<plint> const &blocks = getLocalInfo().getBlocks();
    std::vector<BlockStatistics const *> individualStatistics;
    // Prepare a vector containing the BlockStatistics of all components
    for (pluint iBlock = 0; iBlock < blocks.size(); ++iBlock) {
        plint blockId = blocks[iBlock];
        individualStatistics.push_back(&getComponent(blockId).getInternalStatistics());
    }

    // Execute reduction operation on all individual statistics and store result into
    //   statistics of current MultiBlock.
    combinedStatistics->combine(individualStatistics, this->getInternalStatistics());
    // Copy result to each individual statistics
    for (pluint iBlock = 0; iBlock < blocks.size(); ++iBlock) {
        plint blockId = blocks[iBlock];
        (getComponent(blockId).getInternalStatistics()) = (this->getInternalStatistics());
    }
}

void MultiBlock2D::toggleInternalStatistics(bool statisticsOn_)
{
    statisticsOn = statisticsOn_;
}

bool MultiBlock2D::isInternalStatisticsOn() const
{
    return statisticsOn;
}

PeriodicitySwitch2D const &MultiBlock2D::periodicity() const
{
    return periodicitySwitch;
}

PeriodicitySwitch2D &MultiBlock2D::periodicity()
{
    return periodicitySwitch;
}

modif::ModifT MultiBlock2D::getInternalTypeOfModification() const
{
    return internalModifT;
}

void MultiBlock2D::setInternalTypeOfModification(modif::ModifT internalModifT_)
{
    internalModifT = internalModifT_;
}

void MultiBlock2D::setRefinementLevel(plint newLevel)
{
    multiBlockManagement.setRefinementLevel(newLevel);
}

void MultiBlock2D::executeInternalProcessors()
{
    global::profiler().start("dataProcessor");
    // Execute all automatic internal processors.
    for (plint iLevel = 0; iLevel <= maxProcessorLevel; ++iLevel) {
        executeInternalProcessors(iLevel);
    }
    // Duplicate boundaries at least once in case there is no automatic processor.
    if (maxProcessorLevel == -1) {
        global::profiler().start("envelope-update");
        this->duplicateOverlaps(internalModifT);
        global::profiler().stop("envelope-update");
    }
    global::profiler().stop("dataProcessor");
}

void MultiBlock2D::executeInternalProcessors(plint level, bool communicate)
{
    std::vector<plint> const &blocks = getLocalInfo().getBlocks();
    for (pluint iBlock = 0; iBlock < blocks.size(); ++iBlock) {
        plint blockId = blocks[iBlock];
        getComponent(blockId).executeInternalProcessors(level);
    }
    if (communicate) {
        duplicateOverlapsInModifiedMultiBlocks(level);
    }
}

void MultiBlock2D::subscribeProcessor(
    plint level, std::vector<MultiBlock2D *> modifiedBlocks,
    std::vector<modif::ModifT> typeOfModification, bool includesEnvelope)
{
    maxProcessorLevel = std::max(level, maxProcessorLevel);

    if (level >= 0) {
        addModifiedBlocks(
            level, modifiedBlocks, typeOfModification, multiBlocksChangedByAutomaticProcessors,
            includesEnvelope);
    } else if (level < 0) {
        addModifiedBlocks(
            -level, modifiedBlocks, typeOfModification, multiBlocksChangedByManualProcessors,
            includesEnvelope);
    }
}

void MultiBlock2D::storeProcessor(
    DataProcessorGenerator2D const &generator, std::vector<MultiBlock2D *> multiBlocks, plint level)
{
    storedProcessors.push_back(ProcessorStorage2D(generator, multiBlocks, level));
}

std::vector<MultiBlock2D::ProcessorStorage2D> const &MultiBlock2D::getStoredProcessors() const
{
    return storedProcessors;
}

void MultiBlock2D::addModifiedBlocks(
    plint level, std::vector<MultiBlock2D *> modifiedBlocks,
    std::vector<modif::ModifT> typeOfModification,
    std::vector<std::vector<BlockAndModif> > &multiBlockCollection, bool includesEnvelope)
{
    PLB_PRECONDITION(modifiedBlocks.size() == typeOfModification.size());
    // Resize vector which collects modified blocks (resize needs to be
    //   done even when no block is added, to avoid memory violations
    //   during read access to the vector).
    if ((pluint)level >= multiBlockCollection.size()) {
        multiBlockCollection.resize(level + 1);
    }
    // Unless envelope is already included in the domain of application of the data
    //   processor, subscribe modified blocks for an update of the envelope.
    if (!includesEnvelope) {
        for (pluint iNewBlock = 0; iNewBlock < modifiedBlocks.size(); ++iNewBlock) {
            bool alreadyAdded = false;
            // Check if the block is already in the collection.
            // Note: It may seem stupid to use a linear-complexity algorithm to
            // find existing blocks. However, it would be a mistake to store the data
            // processors in a sorted structure. Indeed, sorting pointers leads to
            // an unpredictable order (because the pointers have an unpredictable value).
            // This is a problem in parallel programs, because different threads may find a
            // different order, and the communication pattern between processes gets messed up.
            for (pluint iOriginal = 0; iOriginal < multiBlockCollection[level].size(); ++iOriginal)
            {
                MultiBlock2D *existingBlock = multiBlockCollection[level][iOriginal].first;
                if (existingBlock == modifiedBlocks[iNewBlock]) {
                    // If the block already exists, simply update the type of modification
                    // that applies to it.
                    modif::ModifT &existingModif = multiBlockCollection[level][iOriginal].second;
                    existingModif = combine(existingModif, typeOfModification[iNewBlock]);
                    alreadyAdded = true;
                }
            }
            // If the block is not already in the structure, add it.
            if (!alreadyAdded) {
                multiBlockCollection[level].push_back(
                    BlockAndModif(modifiedBlocks[iNewBlock], typeOfModification[iNewBlock]));
            }
        }
    }
}

void MultiBlock2D::duplicateOverlapsInModifiedMultiBlocks(plint level)
{
    if (level == 0) {
        duplicateOverlapsAtLevelZero(multiBlocksChangedByAutomaticProcessors[level]);
    } else if (level > 0) {
        if (level < (plint)multiBlocksChangedByAutomaticProcessors.size()) {
            duplicateOverlapsInModifiedMultiBlocks(multiBlocksChangedByAutomaticProcessors[level]);
        }
    } else {  // level < 0
        if (-level < (plint)multiBlocksChangedByManualProcessors.size()) {
            duplicateOverlapsInModifiedMultiBlocks(multiBlocksChangedByManualProcessors[-level]);
        }
    }
}

void MultiBlock2D::duplicateOverlapsInModifiedMultiBlocks(std::vector<BlockAndModif> &multiBlocks)
{
    for (pluint iBlock = 0; iBlock < multiBlocks.size(); ++iBlock) {
        multiBlocks[iBlock].first->duplicateOverlaps(multiBlocks[iBlock].second);
    }
}

void MultiBlock2D::duplicateOverlapsAtLevelZero(std::vector<BlockAndModif> &multiBlocks)
{
    bool treatedThis = false;
    for (pluint iBlock = 0; iBlock < multiBlocks.size(); ++iBlock) {
        MultiBlock2D *modifiedBlock = multiBlocks[iBlock].first;
        modif::ModifT modificationType = multiBlocks[iBlock].second;
        if (modifiedBlock == this) {
            treatedThis = true;
            // If it's the current multi-block we are treating, make sure
            //   type of modification is equal to internalModifT or stronger.
            this->duplicateOverlaps(combine(modificationType, internalModifT));
        } else {
            modifiedBlock->duplicateOverlaps(modificationType);
        }
    }
    // If current multi-block has not already been treated, duplicate
    //   overlaps explicitly (because overlaps are expected to be duplicated
    //   in any case at level 0).
    if (!treatedThis) {
        this->duplicateOverlaps(internalModifT);
    }
}

/* *************** Class MultiBlockRegistration2D ******************************** */

MultiBlockRegistration2D::MultiBlockRegistration2D() { }

id_t MultiBlockRegistration2D::announce(MultiBlock2D &block)
{
    id_t id = uniqueId.getId();
    multiBlocks.insert(std::pair<id_t, MultiBlock2D *>(id, &block));
    return id;
}

void MultiBlockRegistration2D::release(MultiBlock2D &block)
{
    std::map<id_t, MultiBlock2D *>::iterator it = multiBlocks.find(block.getId());
    if (it == multiBlocks.end()) {
        throw PlbLogicException("Releasing a block which is not registered.");
    } else {
        uniqueId.releaseId(it->first);
        multiBlocks.erase(it);
    }
}

MultiBlock2D *MultiBlockRegistration2D::find(id_t id)
{
    std::map<id_t, MultiBlock2D *>::iterator it = multiBlocks.find(id);
    if (it == multiBlocks.end()) {
        return 0;
    } else {
        return it->second;
    }
}

// TODO: Why are copy const and equal doing nothing?
MultiBlockRegistration2D::MultiBlockRegistration2D(MultiBlockRegistration2D const &rhs) { }

MultiBlockRegistration2D &MultiBlockRegistration2D::operator=(MultiBlockRegistration2D const &rhs)
{
    return *this;
}

MultiBlockRegistration2D &multiBlockRegistration2D()
{
    static MultiBlockRegistration2D instance;
    return instance;
}

}  // namespace plb
