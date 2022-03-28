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
 * The 3D multiblock -- implementation file.
 */
#include "multiBlock/multiBlock3D.h"

#include <algorithm>
#include <cmath>

#include "atomicBlock/atomicBlock3D.h"
#include "core/plbDebug.h"
#include "core/plbProfiler.h"
#include "multiBlock/defaultMultiBlockPolicy3D.h"
#include "multiBlock/multiBlockOperations3D.h"
#include "multiBlock/multiBlockSerializer3D.h"

namespace plb {

/* *************** Class MultiStatSubscriber3D ****************************** */

MultiStatSubscriber3D::MultiStatSubscriber3D(MultiBlock3D &multiBlock_) :
    multiBlock(multiBlock_) { }

plint MultiStatSubscriber3D::subscribeAverage()
{
    std::vector<plint> const &blocks = multiBlock.getLocalInfo().getBlocks();
    for (pluint iBlock = 0; iBlock < blocks.size(); ++iBlock) {
        plint blockId = blocks[iBlock];
        multiBlock.getComponent(blockId).internalStatSubscription().subscribeAverage();
    }
    return multiBlock.getInternalStatistics().subscribeAverage();
}

plint MultiStatSubscriber3D::subscribeSum()
{
    std::vector<plint> const &blocks = multiBlock.getLocalInfo().getBlocks();
    for (pluint iBlock = 0; iBlock < blocks.size(); ++iBlock) {
        plint blockId = blocks[iBlock];
        multiBlock.getComponent(blockId).internalStatSubscription().subscribeSum();
    }
    return multiBlock.getInternalStatistics().subscribeSum();
}

plint MultiStatSubscriber3D::subscribeMax()
{
    std::vector<plint> const &blocks = multiBlock.getLocalInfo().getBlocks();
    for (pluint iBlock = 0; iBlock < blocks.size(); ++iBlock) {
        plint blockId = blocks[iBlock];
        multiBlock.getComponent(blockId).internalStatSubscription().subscribeMax();
    }
    return multiBlock.getInternalStatistics().subscribeMax();
}

plint MultiStatSubscriber3D::subscribeIntSum()
{
    std::vector<plint> const &blocks = multiBlock.getLocalInfo().getBlocks();
    for (pluint iBlock = 0; iBlock < blocks.size(); ++iBlock) {
        plint blockId = blocks[iBlock];
        multiBlock.getComponent(blockId).internalStatSubscription().subscribeIntSum();
    }
    return multiBlock.getInternalStatistics().subscribeIntSum();
}

/* *************** Class PeriodicitySwitch3D ******************************** */

PeriodicitySwitch3D::PeriodicitySwitch3D(MultiBlock3D &block_) :
    // Default all boundaries to non-periodic. Note: it is not possible to send a
    // signal "signalPeriodicity" to "block" at this point, because "block" might
    // not be fully constructed, and virtual function calls not available.
    periodicity(false, false, false),
    block(block_)
{ }

PeriodicitySwitch3D::PeriodicitySwitch3D(MultiBlock3D &block_, PeriodicitySwitch3D const &rhs) :
    periodicity(rhs.periodicity), block(block_)
{ }

PeriodicitySwitch3D &PeriodicitySwitch3D::operator=(PeriodicitySwitch3D const &rhs)
{
    // Don't modify the MultiBlock3D reference, as it is immutable.
    periodicity = rhs.periodicity;
    return *this;
}

void PeriodicitySwitch3D::toggle(plint direction, bool periodic)
{
    PLB_PRECONDITION(direction == 0 || direction == 1 || direction == 2);
    periodicity[direction] = periodic;
    // Make sure periodic envelope is filled with data, else all
    //   subsequent operations are wrong.
    block.signalPeriodicity();
    block.duplicateOverlaps(modif::dataStructure);
}

void PeriodicitySwitch3D::toggleAll(bool periodic)
{
    periodicity[0] = periodic;
    periodicity[1] = periodic;
    periodicity[2] = periodic;
    // Make sure periodic envelope is filled with data, else all
    //   subsequent operations are wrong.
    block.signalPeriodicity();
    block.duplicateOverlaps(modif::dataStructure);
}

bool PeriodicitySwitch3D::get(plint direction) const
{
    PLB_PRECONDITION(direction == 0 || direction == 1 || direction == 2);

    return periodicity[direction];
}

bool PeriodicitySwitch3D::get(plint normalX, plint normalY, plint normalZ) const
{
    bool testX = normalX != 0;
    bool testY = normalY != 0;
    bool testZ = normalZ != 0;
    return (
        (!testX || (testX && periodicity[0])) && (!testY || (testY && periodicity[1]))
        && (!testZ || (testZ && periodicity[2])));
}

Box3D PeriodicitySwitch3D::getPeriodicEnvelope(Box3D const &bulk, plint envelopeWidth) const
{
    Box3D envelope(bulk);
    if (get(0)) {  // If periodic in x, extend bulk in x-direction
        envelope.x0 -= envelopeWidth;
        envelope.x1 += envelopeWidth;
    }
    if (get(1)) {  // If periodic in y, extend bulk in y-direction
        envelope.y0 -= envelopeWidth;
        envelope.y1 += envelopeWidth;
    }
    if (get(2)) {  // If periodic in z, extend bulk in z-direction
        envelope.z0 -= envelopeWidth;
        envelope.z1 += envelopeWidth;
    }
    return envelope;
}

/* *************** Class MultiBlock3D::ProcessorStorage3D ******************* */

MultiBlock3D::ProcessorStorage3D::ProcessorStorage3D(
    DataProcessorGenerator3D const &generator_, std::vector<MultiBlock3D *> const &multiBlocks_,
    plint level_) :
    generator(generator_.clone()), multiBlockIds(multiBlocks_.size()), level(level_)
{
    for (pluint iBlock = 0; iBlock < multiBlockIds.size(); ++iBlock) {
        multiBlockIds[iBlock] = multiBlocks_[iBlock]->getId();
    }
}

MultiBlock3D::ProcessorStorage3D::~ProcessorStorage3D()
{
    delete generator;
}

MultiBlock3D::ProcessorStorage3D::ProcessorStorage3D(MultiBlock3D::ProcessorStorage3D const &rhs) :
    generator(rhs.generator->clone()), multiBlockIds(rhs.multiBlockIds), level(rhs.level)
{ }

MultiBlock3D::ProcessorStorage3D &MultiBlock3D::ProcessorStorage3D::operator=(
    MultiBlock3D::ProcessorStorage3D const &rhs)
{
    MultiBlock3D::ProcessorStorage3D(rhs).swap(*this);
    return *this;
}

void MultiBlock3D::ProcessorStorage3D::swap(MultiBlock3D::ProcessorStorage3D &rhs)
{
    std::swap(generator, rhs.generator);
    multiBlockIds.swap(rhs.multiBlockIds);
    std::swap(level, rhs.level);
}

MultiBlock3D::ProcessorStorage3D *MultiBlock3D::ProcessorStorage3D::clone() const
{
    return new MultiBlock3D::ProcessorStorage3D(*this);
}

DataProcessorGenerator3D const &MultiBlock3D::ProcessorStorage3D::getGenerator() const
{
    return *generator;
}

std::vector<id_t> const &MultiBlock3D::ProcessorStorage3D::getMultiBlockIds() const
{
    return multiBlockIds;
}

std::vector<MultiBlock3D *> MultiBlock3D::ProcessorStorage3D::getMultiBlocks() const
{
    std::vector<MultiBlock3D *> multiBlocks(multiBlockIds.size());
    for (pluint iBlock = 0; iBlock < multiBlockIds.size(); ++iBlock) {
        multiBlocks[iBlock] = multiBlockRegistration3D().find(multiBlockIds[iBlock]);
        PLB_ASSERT(multiBlocks[iBlock]);
    }
    return multiBlocks;
}

plint MultiBlock3D::ProcessorStorage3D::getLevel() const
{
    return level;
}

void MultiBlock3D::ProcessorStorage3D::replace(id_t oldBlock, id_t newBlock)
{
    for (pluint iBlock = 0; iBlock < multiBlockIds.size(); ++iBlock) {
        if (multiBlockIds[iBlock] == oldBlock) {
            multiBlockIds[iBlock] = newBlock;
        }
    }
}

/* *************** Class MultiBlock3D *************************************** */

MultiBlock3D::MultiBlock3D(
    MultiBlockManagement3D const &multiBlockManagement_, BlockCommunicator3D *blockCommunicator_,
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
    id = multiBlockRegistration3D().announce(*this);
}

MultiBlock3D::MultiBlock3D(plint nx, plint ny, plint nz, plint envelopeWidth) :
    multiBlockManagement(defaultMultiBlockPolicy3D().getMultiBlockManagement(
        Box3D(0, nx - 1, 0, ny - 1, 0, nz - 1), envelopeWidth)),
    maxProcessorLevel(-1),
    blockCommunicator(defaultMultiBlockPolicy3D().getBlockCommunicator()),
    internalStatistics(),
    combinedStatistics(defaultMultiBlockPolicy3D().getCombinedStatistics()),
    statSubscriber(*this),
    statisticsOn(true),
    periodicitySwitch(*this),
    internalModifT(modif::staticVariables)
{
    id = multiBlockRegistration3D().announce(*this);
}

MultiBlock3D::MultiBlock3D(MultiBlock3D const &rhs) :
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
    id = multiBlockRegistration3D().announce(*this);
}

MultiBlock3D::MultiBlock3D(MultiBlock3D const &rhs, Box3D subDomain, bool crop) :
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
    id = multiBlockRegistration3D().announce(*this);
}

void MultiBlock3D::swap(MultiBlock3D &rhs)
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

MultiBlock3D::~MultiBlock3D()
{
    delete blockCommunicator;
    delete combinedStatistics;
    multiBlockRegistration3D().release(*this);
}

id_t MultiBlock3D::getId() const
{
    return id;
}

void MultiBlock3D::initialize()
{
    executeInternalProcessors();
}

Box3D MultiBlock3D::getBoundingBox() const
{
    return multiBlockManagement.getBoundingBox();
}

plint MultiBlock3D::getNx() const
{
    return multiBlockManagement.getBoundingBox().getNx();
}

plint MultiBlock3D::getNy() const
{
    return multiBlockManagement.getBoundingBox().getNy();
}

plint MultiBlock3D::getNz() const
{
    return multiBlockManagement.getBoundingBox().getNz();
}

MultiBlockManagement3D const &MultiBlock3D::getMultiBlockManagement() const
{
    return multiBlockManagement;
}

void MultiBlock3D::setCoProcessors(std::map<plint, int> const &coProcessors)
{
    multiBlockManagement.setCoProcessors(coProcessors);
}

LocalMultiBlockInfo3D const &MultiBlock3D::getLocalInfo() const
{
    return multiBlockManagement.getLocalInfo();
}

SparseBlockStructure3D const &MultiBlock3D::getSparseBlockStructure() const
{
    return multiBlockManagement.getSparseBlockStructure();
}

BlockStatistics &MultiBlock3D::getInternalStatistics()
{
    return internalStatistics;
}

BlockStatistics const &MultiBlock3D::getInternalStatistics() const
{
    return internalStatistics;
}

CombinedStatistics const &MultiBlock3D::getCombinedStatistics() const
{
    return *combinedStatistics;
}

StatSubscriber &MultiBlock3D::internalStatSubscription()
{
    return statSubscriber;
}

BlockCommunicator3D const &MultiBlock3D::getBlockCommunicator() const
{
    return *blockCommunicator;
}

void MultiBlock3D::duplicateOverlaps(modif::ModifT whichData)
{
    this->getBlockCommunicator().duplicateOverlaps(*this, whichData);
}

void MultiBlock3D::signalPeriodicity()
{
    getBlockCommunicator().signalPeriodicity();
}

DataSerializer *MultiBlock3D::getBlockSerializer(
    Box3D const &domain, IndexOrdering::OrderingT ordering) const
{
    return new MultiBlockFastSerializer3D(*this, domain, ordering);
}

DataUnSerializer *MultiBlock3D::getBlockUnSerializer(
    Box3D const &domain, IndexOrdering::OrderingT ordering)
{
    return new MultiBlockFastUnSerializer3D(*this, domain, ordering);
}

void MultiBlock3D::resetFlags()
{
    std::vector<plint> const &blocks = getLocalInfo().getBlocks();
    for (pluint iBlock = 0; iBlock < blocks.size(); ++iBlock) {
        plint blockId = blocks[iBlock];
        getComponent(blockId).setFlag(false);
    }
}

void MultiBlock3D::getDynamicsDict(Box3D domain, std::map<std::string, int> &dict)
{
    return dict.clear();
}

void MultiBlock3D::evaluateStatistics()
{
    std::vector<plint> const &blocks = getLocalInfo().getBlocks();
    for (pluint iBlock = 0; iBlock < blocks.size(); ++iBlock) {
        plint blockId = blocks[iBlock];
        getComponent(blockId).evaluateStatistics();
    }
    if (isInternalStatisticsOn())
        reduceStatistics();
}

void MultiBlock3D::reduceStatistics()
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

void MultiBlock3D::toggleInternalStatistics(bool statisticsOn_)
{
    statisticsOn = statisticsOn_;
}

bool MultiBlock3D::isInternalStatisticsOn() const
{
    return statisticsOn;
}

PeriodicitySwitch3D const &MultiBlock3D::periodicity() const
{
    return periodicitySwitch;
}

PeriodicitySwitch3D &MultiBlock3D::periodicity()
{
    return periodicitySwitch;
}

modif::ModifT MultiBlock3D::getInternalTypeOfModification() const
{
    return internalModifT;
}

void MultiBlock3D::setInternalTypeOfModification(modif::ModifT internalModifT_)
{
    internalModifT = internalModifT_;
}

void MultiBlock3D::setRefinementLevel(plint newLevel)
{
    multiBlockManagement.setRefinementLevel(newLevel);
}

void MultiBlock3D::setDataTransfer(BlockDataTransfer3D *newDataTransfer)
{
    std::vector<plint> const &blocks = getLocalInfo().getBlocks();
    for (pluint iBlock = 0; iBlock < blocks.size(); ++iBlock) {
        plint blockId = blocks[iBlock];
        getComponent(blockId).setDataTransfer(newDataTransfer->clone());
    }
    delete newDataTransfer;
}

void MultiBlock3D::executeInternalProcessors()
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

void MultiBlock3D::executeInternalProcessors(plint level, bool communicate)
{
    if (level < 0) {
        global::timer("execute_dp").start();
    }
    std::vector<plint> const &blocks = getLocalInfo().getBlocks();
    for (pluint iBlock = 0; iBlock < blocks.size(); ++iBlock) {
        plint blockId = blocks[iBlock];
        getComponent(blockId).executeInternalProcessors(level);
    }
    if (level < 0) {
        global::timer("execute_dp").stop();
        global::timer("communicate_dp").start();
    }
    if (communicate) {
        duplicateOverlapsInModifiedMultiBlocks(level);
    }
    if (level < 0) {
        global::timer("communicate_dp").stop();
    }
}

void MultiBlock3D::subscribeProcessor(
    plint level, std::vector<MultiBlock3D *> modifiedBlocks,
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

void MultiBlock3D::storeProcessor(
    DataProcessorGenerator3D const &generator, std::vector<MultiBlock3D *> multiBlocks, plint level)
{
    storedProcessors.push_back(ProcessorStorage3D(generator, multiBlocks, level));
}

std::vector<MultiBlock3D::ProcessorStorage3D> const &MultiBlock3D::getStoredProcessors() const
{
    return storedProcessors;
}

void MultiBlock3D::addModifiedBlocks(
    plint level, std::vector<MultiBlock3D *> modifiedBlocks,
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
                MultiBlock3D *existingBlock = multiBlockCollection[level][iOriginal].first;
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

void MultiBlock3D::duplicateOverlapsInModifiedMultiBlocks(plint level)
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

void MultiBlock3D::duplicateOverlapsInModifiedMultiBlocks(std::vector<BlockAndModif> &multiBlocks)
{
    for (pluint iBlock = 0; iBlock < multiBlocks.size(); ++iBlock) {
        multiBlocks[iBlock].first->duplicateOverlaps(multiBlocks[iBlock].second);
    }
}

void MultiBlock3D::duplicateOverlapsAtLevelZero(std::vector<BlockAndModif> &multiBlocks)
{
    bool treatedThis = false;
    for (pluint iBlock = 0; iBlock < multiBlocks.size(); ++iBlock) {
        MultiBlock3D *modifiedBlock = multiBlocks[iBlock].first;
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

/* *************** Class MultiBlockRegistration3D ******************************** */

MultiBlockRegistration3D::MultiBlockRegistration3D() { }

id_t MultiBlockRegistration3D::announce(MultiBlock3D &block)
{
    id_t id = uniqueId.getId();
    multiBlocks.insert(std::pair<id_t, MultiBlock3D *>(id, &block));
    return id;
}

void MultiBlockRegistration3D::release(MultiBlock3D &block)
{
    std::map<id_t, MultiBlock3D *>::iterator it = multiBlocks.find(block.getId());
    if (it == multiBlocks.end()) {
        throw PlbLogicException("Releasing a block which is not registered.");
    } else {
        uniqueId.releaseId(it->first);
        multiBlocks.erase(it);
    }
}

MultiBlock3D *MultiBlockRegistration3D::find(id_t id)
{
    std::map<id_t, MultiBlock3D *>::iterator it = multiBlocks.find(id);
    if (it == multiBlocks.end()) {
        return 0;
    } else {
        return it->second;
    }
}

// TODO: WHy copy const and equal do nothing?
MultiBlockRegistration3D::MultiBlockRegistration3D(MultiBlockRegistration3D const &rhs) { }

MultiBlockRegistration3D &MultiBlockRegistration3D::operator=(MultiBlockRegistration3D const &rhs)
{
    return *this;
}

MultiBlockRegistration3D &multiBlockRegistration3D()
{
    static MultiBlockRegistration3D instance;
    return instance;
}

}  // namespace plb
