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
 * The 3D multiblock -- header file.
 */
#ifndef MULTI_BLOCK_3D_H
#define MULTI_BLOCK_3D_H

#include <string>
#include <utility>
#include <vector>

#include "atomicBlock/dataProcessor3D.h"
#include "core/array.h"
#include "core/block3D.h"
#include "core/blockStatistics.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "multiBlock/blockCommunicator3D.h"
#include "multiBlock/combinedStatistics.h"
#include "multiBlock/multiBlockManagement3D.h"

namespace plb {

class AtomicBlock3D;
class MultiBlock3D;
class MultiBlockRegistration3D;
template <typename T>
class TypedAtomicBlock3D;
template <typename T>
class EulerianAtomicBlock3D;

/// Handles statistics subscriptions for the MultiBlockLattice3D
class MultiStatSubscriber3D : public StatSubscriber {
public:
    MultiStatSubscriber3D(MultiBlock3D &multiBlock_);
    /// Subscribe a new observable for which the average value is computed.
    virtual plint subscribeAverage();
    /// Subscribe a new observable for which the sum is computed.
    virtual plint subscribeSum();
    /// Subscribe a new observable for which the maximum is computed.
    virtual plint subscribeMax();
    /// Subscribe a new integer observable for which the sum is computed.
    virtual plint subscribeIntSum();

private:
    MultiBlock3D &multiBlock;
};

class PeriodicitySwitch3D {
public:
    /// The constructor defaults all directions to false (i.e., non-periodic).
    PeriodicitySwitch3D(MultiBlock3D &block_);
    PeriodicitySwitch3D(MultiBlock3D &block_, PeriodicitySwitch3D const &rhs);
    PeriodicitySwitch3D &operator=(PeriodicitySwitch3D const &rhs);

    /// Set periodicity status of a direction (direction=0 means x-direction etc.)
    void toggle(plint direction, bool periodic);
    /// Set periodicity status of all directions synchronously
    void toggleAll(bool periodic);
    /// Get periodicity status of a direction;
    bool get(plint direction) const;
    /// Get periodicity along a general direction;
    bool get(plint normalX, plint normalY, plint normalZ) const;
    /// Extend the bulk in each periodic direction, and return the result.
    Box3D getPeriodicEnvelope(Box3D const &bulk, plint envelopeWidth) const;

private:
    Array<bool, 3> periodicity;
    MultiBlock3D &block;
};

struct BlockDataTransfer3D;

class MultiBlock3D : public Block3D {
public:
    typedef std::pair<MultiBlock3D *, modif::ModifT> BlockAndModif;

public:
    class ProcessorStorage3D {
    public:
        ProcessorStorage3D(
            DataProcessorGenerator3D const &generator_,
            std::vector<MultiBlock3D *> const &multiBlocks_, plint level_);
        ~ProcessorStorage3D();
        ProcessorStorage3D(ProcessorStorage3D const &rhs);
        ProcessorStorage3D &operator=(ProcessorStorage3D const &rhs);
        void swap(ProcessorStorage3D &rhs);
        ProcessorStorage3D *clone() const;
        DataProcessorGenerator3D const &getGenerator() const;
        std::vector<id_t> const &getMultiBlockIds() const;
        std::vector<MultiBlock3D *> getMultiBlocks() const;
        plint getLevel() const;
        void replace(id_t oldBlock, id_t newBlock);

    private:
        DataProcessorGenerator3D *generator;
        std::vector<id_t> multiBlockIds;
        plint level;
    };

public:
    MultiBlock3D(
        MultiBlockManagement3D const &multiBlockManagement_,
        BlockCommunicator3D *blockCommunicator_, CombinedStatistics *combinedStatistics_);
    MultiBlock3D(plint nx, plint ny, plint nz, plint envelopeWidth);
    MultiBlock3D(MultiBlock3D const &rhs);
    MultiBlock3D(MultiBlock3D const &rhs, Box3D subDomain, bool crop);
    void swap(MultiBlock3D &rhs);
    virtual ~MultiBlock3D();
    id_t getId() const;
    virtual int getStaticId() const = 0;
    virtual MultiBlock3D *clone() const = 0;
    virtual MultiBlock3D *clone(MultiBlockManagement3D const &newManagement) const = 0;
    /// Initialize block content by executing internal processors once.
    void initialize();

public:
    virtual Box3D getBoundingBox() const;
    /// Get number of cells in x-direction.
    plint getNx() const;
    /// Get number of cells in y-direction.
    plint getNy() const;
    /// Get number of cells in z-direction.
    plint getNz() const;
    /// Execute all internal dataProcessors at positive or zero level.
    void executeInternalProcessors();
    /// Execute all internal dataProcessors at a given level.
    void executeInternalProcessors(plint level, bool communicate = true);
    /// After adding an internal processor to the atomic-blocks, subscribe it
    /// in the multi-block to guarantee it will be executed.
    void subscribeProcessor(
        plint level, std::vector<MultiBlock3D *> modifiedBlocks,
        std::vector<modif::ModifT> typeOfModification, bool includesEnvelope);
    void storeProcessor(
        DataProcessorGenerator3D const &generator, std::vector<MultiBlock3D *> multiBlocks,
        plint level);
    std::vector<ProcessorStorage3D> const &getStoredProcessors() const;

public:
    MultiBlockManagement3D const &getMultiBlockManagement() const;
    void setCoProcessors(std::map<plint, int> const &coProcessors);
    LocalMultiBlockInfo3D const &getLocalInfo() const;
    SparseBlockStructure3D const &getSparseBlockStructure() const;
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
    PeriodicitySwitch3D const &periodicity() const;
    PeriodicitySwitch3D &periodicity();
    /// Returns: which kind of data is modified by level-0 processors and by
    ///   dynamics objects.
    modif::ModifT getInternalTypeOfModification() const;
    /// Specify which kind of data is modifified by level-0 processors and by
    ///   dynamics objects.
    void setInternalTypeOfModification(modif::ModifT internalModifT_);
    void setRefinementLevel(plint newLevel);
    /// Assign a data transfer policy to all atomic-blocks.
    void setDataTransfer(BlockDataTransfer3D *newDataTransfer);

public:
    virtual AtomicBlock3D &getComponent(plint blockId) = 0;
    virtual AtomicBlock3D const &getComponent(plint blockId) const = 0;
    virtual plint sizeOfCell() const = 0;
    virtual plint getCellDim() const = 0;

private:
    void addModifiedBlocks(
        plint level, std::vector<MultiBlock3D *> modifiedBlocks,
        std::vector<modif::ModifT> typeOfModification,
        std::vector<std::vector<BlockAndModif> > &multiBlockCollection, bool includesEnvelope);
    void duplicateOverlapsInModifiedMultiBlocks(plint level);
    void duplicateOverlapsInModifiedMultiBlocks(std::vector<BlockAndModif> &multiBlocks);
    void duplicateOverlapsAtLevelZero(std::vector<BlockAndModif> &multiBlocks);
    void reduceStatistics();

public:
    BlockCommunicator3D const &getBlockCommunicator() const;
    virtual void copyReceive(
        MultiBlock3D const &fromBlock, Box3D const &fromDomain, Box3D const &toDomain,
        modif::ModifT whichData = modif::dataStructure) = 0;
    void duplicateOverlaps(modif::ModifT whichData);
    void signalPeriodicity();
    virtual DataSerializer *getBlockSerializer(
        Box3D const &domain, IndexOrdering::OrderingT ordering) const;
    virtual DataUnSerializer *getBlockUnSerializer(
        Box3D const &domain, IndexOrdering::OrderingT ordering);
    /// Set all flags of the atomic-blocks to false (to be used by data processors).
    void resetFlags();
    /// Return a unique map for the IDs of dynamics objects present in the multi-block, if any.
    virtual void getDynamicsDict(Box3D domain, std::map<std::string, int> &dict);
    /// Get a string identifier for the type of block. E.g. "lattice3d"
    virtual std::string getBlockName() const = 0;
    /// Get one or two string identifiers for the template parameters of the block.
    ///   E.g. "double" and "d3q19"
    virtual std::vector<std::string> getTypeInfo() const = 0;

private:
    MultiBlockManagement3D multiBlockManagement;
    /// List of MultiBlocks which are modified by the manual processors and require
    /// an update of their envelope.
    std::vector<std::vector<BlockAndModif> > multiBlocksChangedByManualProcessors;
    /// List of MultiBlocks which are modified by the automatic processors and require
    /// an update of their envelope.
    std::vector<std::vector<BlockAndModif> > multiBlocksChangedByAutomaticProcessors;
    plint maxProcessorLevel;
    std::vector<ProcessorStorage3D> storedProcessors;
    BlockCommunicator3D *blockCommunicator;
    BlockStatistics internalStatistics;
    CombinedStatistics *combinedStatistics;
    MultiStatSubscriber3D statSubscriber;
    bool statisticsOn;
    PeriodicitySwitch3D periodicitySwitch;
    modif::ModifT internalModifT;
    id_t id;
};

class MultiBlockRegistration3D {
public:
    id_t announce(MultiBlock3D &block);
    void release(MultiBlock3D &block);
    MultiBlock3D *find(id_t id);

private:
    MultiBlockRegistration3D();
    MultiBlockRegistration3D(MultiBlockRegistration3D const &rhs);
    MultiBlockRegistration3D &operator=(MultiBlockRegistration3D const &rhs);

private:
    util::UniqueId uniqueId;
    std::map<id_t, MultiBlock3D *> multiBlocks;
    friend MultiBlockRegistration3D &multiBlockRegistration3D();
};

MultiBlockRegistration3D &multiBlockRegistration3D();

}  // namespace plb

#endif  // MULTI_BLOCK_3D_H
