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
 * The 2D multiblock -- header file.
 */
#ifndef MULTI_BLOCK_2D_H
#define MULTI_BLOCK_2D_H

#include <string>
#include <utility>
#include <vector>

#include "atomicBlock/dataProcessor2D.h"
#include "core/array.h"
#include "core/block2D.h"
#include "core/blockStatistics.h"
#include "core/geometry2D.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "multiBlock/blockCommunicator2D.h"
#include "multiBlock/combinedStatistics.h"
#include "multiBlock/localMultiBlockInfo2D.h"
#include "multiBlock/multiBlockManagement2D.h"

namespace plb {

class AtomicBlock2D;
class MultiBlock2D;
class MultiBlockRegistration2D;
template <typename T>
class TypedAtomicBlock2D;
template <typename T>
class EulerianAtomicBlock2D;

/// Handles statistics subscriptions for the MultiBlockLattice2D
class MultiStatSubscriber2D : public StatSubscriber {
public:
    MultiStatSubscriber2D(MultiBlock2D &multiBlock_);
    /// Subscribe a new observable for which the average value is computed.
    virtual plint subscribeAverage();
    /// Subscribe a new observable for which the sum is computed.
    virtual plint subscribeSum();
    /// Subscribe a new observable for which the maximum is computed.
    virtual plint subscribeMax();
    /// Subscribe a new integer observable for which the sum is computed.
    virtual plint subscribeIntSum();

private:
    MultiBlock2D &multiBlock;
};

class PeriodicitySwitch2D {
public:
    /// The constructor defaults all directions to false (i.e., non-periodic).
    PeriodicitySwitch2D(MultiBlock2D &block_);
    PeriodicitySwitch2D(MultiBlock2D &block_, PeriodicitySwitch2D const &rhs);
    PeriodicitySwitch2D &operator=(PeriodicitySwitch2D const &rhs);

    /// Set periodicity status of a direction (direction=0 means x-direction etc.)
    void toggle(plint direction, bool periodic);
    /// Set periodicity status of all directions synchronously
    void toggleAll(bool periodic);
    /// Get periodicity status of a direction;
    bool get(plint direction) const;
    /// Get periodicity along a general direction;
    bool get(plint normalX, plint normalY) const;
    /// Extend the bulk in each periodic direction, and return the result.
    Box2D getPeriodicEnvelope(Box2D const &bulk, plint envelopeWidth) const;

private:
    Array<bool, 2> periodicity;
    MultiBlock2D &block;
};

class MultiBlock2D : public Block2D {
public:
    typedef std::pair<MultiBlock2D *, modif::ModifT> BlockAndModif;

public:
    class ProcessorStorage2D {
    public:
        ProcessorStorage2D(
            DataProcessorGenerator2D const &generator_,
            std::vector<MultiBlock2D *> const &multiBlocks_, plint level_);
        ~ProcessorStorage2D();
        ProcessorStorage2D(ProcessorStorage2D const &rhs);
        ProcessorStorage2D &operator=(ProcessorStorage2D const &rhs);
        void swap(ProcessorStorage2D &rhs);
        ProcessorStorage2D *clone() const;
        DataProcessorGenerator2D const &getGenerator() const;
        std::vector<id_t> const &getMultiBlockIds() const;
        std::vector<MultiBlock2D *> getMultiBlocks() const;
        plint getLevel() const;
        void replace(id_t oldBlock, id_t newBlock);

    private:
        DataProcessorGenerator2D *generator;
        std::vector<id_t> multiBlockIds;
        plint level;
    };

public:
    MultiBlock2D(
        MultiBlockManagement2D const &multiBlockManagement_,
        BlockCommunicator2D *blockCommunicator_, CombinedStatistics *combinedStatistics_);
    MultiBlock2D(plint nx, plint ny, plint envelopeWidth);
    MultiBlock2D(MultiBlock2D const &rhs);
    MultiBlock2D(MultiBlock2D const &rhs, Box2D subDomain, bool crop);
    void swap(MultiBlock2D &rhs);
    virtual ~MultiBlock2D();
    id_t getId() const;
    virtual int getStaticId() const = 0;
    virtual MultiBlock2D *clone() const = 0;
    virtual MultiBlock2D *clone(MultiBlockManagement2D const &newManagement) const = 0;
    /// Initialize block content by executing internal processors once.
    void initialize();

public:
    virtual Box2D getBoundingBox() const;
    /// Get number of cells in x-direction.
    plint getNx() const;
    /// Get number of cells in y-direction.
    plint getNy() const;
    /// Execute all internal dataProcessors at positive or zero level.
    void executeInternalProcessors();
    /// Execute all internal dataProcessors at a given level.
    void executeInternalProcessors(plint level, bool communicate = true);
    /// After adding an internal processor to the atomic-blocks, subscribe it
    /// in the multi-block to guarantee it will be executed.
    void subscribeProcessor(
        plint level, std::vector<MultiBlock2D *> modifiedBlocks,
        std::vector<modif::ModifT> typeOfModification, bool includesEnvelope);
    void storeProcessor(
        DataProcessorGenerator2D const &generator, std::vector<MultiBlock2D *> multiBlocks,
        plint level);
    std::vector<ProcessorStorage2D> const &getStoredProcessors() const;

public:
    MultiBlockManagement2D const &getMultiBlockManagement() const;
    LocalMultiBlockInfo2D const &getLocalInfo() const;
    SparseBlockStructure2D const &getSparseBlockStructure() const;
    /// Get a handle to internal statistics. Don't use this to subscribe new
    /// statistics. Use the method internalStatSubscription() instead.
    BlockStatistics &getInternalStatistics();
    /// Get a constant handle to internal statistics.
    BlockStatistics const &getInternalStatistics() const;
    /// Get object to subscribe new internal statistics.
    StatSubscriber &internalStatSubscription();
    /// Copy running statistics to public statistics, and reset running stats.
    void evaluateStatistics();
    CombinedStatistics const &getCombinedStatistics() const;
    void toggleInternalStatistics(bool statisticsOn_);
    bool isInternalStatisticsOn() const;
    PeriodicitySwitch2D const &periodicity() const;
    PeriodicitySwitch2D &periodicity();
    /// Returns: which kind of data is modified by level-0 processors and by
    ///   dynamics objects.
    modif::ModifT getInternalTypeOfModification() const;
    /// Specify which kind of data is modifified by level-0 processors and by
    ///   dynamics objects.
    void setInternalTypeOfModification(modif::ModifT internalModifT_);
    void setRefinementLevel(plint newLevel);

public:
    virtual AtomicBlock2D &getComponent(plint blockId) = 0;
    virtual AtomicBlock2D const &getComponent(plint blockId) const = 0;
    virtual plint sizeOfCell() const = 0;
    virtual plint getCellDim() const = 0;

private:
    void addModifiedBlocks(
        plint level, std::vector<MultiBlock2D *> modifiedBlocks,
        std::vector<modif::ModifT> typeOfModification,
        std::vector<std::vector<BlockAndModif> > &multiBlockCollection, bool includesEnvelope);
    void duplicateOverlapsInModifiedMultiBlocks(plint level);
    void duplicateOverlapsInModifiedMultiBlocks(std::vector<BlockAndModif> &multiBlocks);
    void duplicateOverlapsAtLevelZero(std::vector<BlockAndModif> &multiBlocks);
    void reduceStatistics();

public:
    BlockCommunicator2D const &getBlockCommunicator() const;
    virtual void copyReceive(
        MultiBlock2D const &fromBlock, Box2D const &fromDomain, Box2D const &toDomain,
        modif::ModifT whichData = modif::dataStructure) = 0;
    void duplicateOverlaps(modif::ModifT whichData);
    void signalPeriodicity();
    virtual DataSerializer *getBlockSerializer(
        Box2D const &domain, IndexOrdering::OrderingT ordering) const;
    virtual DataUnSerializer *getBlockUnSerializer(
        Box2D const &domain, IndexOrdering::OrderingT ordering);
    /// Set all flags of the atomic-blocks to false (to be used by data processors).
    void resetFlags();
    /// Return a unique map for the IDs of dynamics objects present in the multi-block, if any.
    virtual void getDynamicsDict(Box2D domain, std::map<std::string, int> &dict);
    /// Get a string identifier for the type of block. E.g. "lattice2d"
    virtual std::string getBlockName() const = 0;
    /// Get one or two string identifiers for the template parameters of the block.
    ///   E.g. "double" and "d2q9"
    virtual std::vector<std::string> getTypeInfo() const = 0;

private:
    MultiBlockManagement2D multiBlockManagement;
    /// List of MultiBlocks which are modified by the manual processors and require
    /// an update of their envelope.
    std::vector<std::vector<BlockAndModif> > multiBlocksChangedByManualProcessors;
    /// List of MultiBlocks which are modified by the automatic processors and require
    /// an update of their envelope.
    std::vector<std::vector<BlockAndModif> > multiBlocksChangedByAutomaticProcessors;
    plint maxProcessorLevel;
    std::vector<ProcessorStorage2D> storedProcessors;
    BlockCommunicator2D *blockCommunicator;
    BlockStatistics internalStatistics;
    CombinedStatistics *combinedStatistics;
    MultiStatSubscriber2D statSubscriber;
    bool statisticsOn;
    PeriodicitySwitch2D periodicitySwitch;
    modif::ModifT internalModifT;
    id_t id;
};

class MultiBlockRegistration2D {
public:
    id_t announce(MultiBlock2D &block);
    void release(MultiBlock2D &block);
    MultiBlock2D *find(id_t id);

private:
    MultiBlockRegistration2D();
    MultiBlockRegistration2D(MultiBlockRegistration2D const &rhs);
    MultiBlockRegistration2D &operator=(MultiBlockRegistration2D const &rhs);

private:
    util::UniqueId uniqueId;
    std::map<id_t, MultiBlock2D *> multiBlocks;
    friend MultiBlockRegistration2D &multiBlockRegistration2D();
};

MultiBlockRegistration2D &multiBlockRegistration2D();

}  // namespace plb

#endif  // MULTI_BLOCK_2D_H
