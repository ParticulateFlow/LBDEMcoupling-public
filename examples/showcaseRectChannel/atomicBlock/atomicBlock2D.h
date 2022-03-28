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

#ifndef ATOMIC_BLOCK_2D_H
#define ATOMIC_BLOCK_2D_H

#include <algorithm>
#include <map>

#include "atomicBlock/dataProcessorWrapper2D.h"
#include "atomicBlock/reductiveDataProcessorWrapper2D.h"
#include "core/block2D.h"
#include "core/blockIdentifiers.h"
#include "core/blockStatistics.h"
#include "core/geometry2D.h"
#include "core/globalDefs.h"

namespace plb {

// Forward declarations
class AtomicBlock2D;

/// Handles statistics subscriptions for the BlockLattice2D
class StatSubscriber2D : public StatSubscriber {
public:
    StatSubscriber2D(AtomicBlock2D &block_);
    /// Subscribe a new observable for which the average value is computed.
    virtual plint subscribeAverage();
    /// Subscribe a new observable for which the sum is computed.
    virtual plint subscribeSum();
    /// Subscribe a new observable for which the maximum is computed.
    virtual plint subscribeMax();
    /// Subscribe a new integer observable for which the sum is computed.
    virtual plint subscribeIntSum();

private:
    AtomicBlock2D &block;
};

struct BlockDataTransfer2D {
    virtual ~BlockDataTransfer2D() { }
    virtual plint staticCellSize() const = 0;
    virtual void send(Box2D domain, std::vector<char> &buffer, modif::ModifT kind) const = 0;
    virtual void receive(Box2D domain, std::vector<char> const &buffer, modif::ModifT kind) = 0;
    /// Receive data from a byte-stream into the block, and adjust coordinates if the block contains
    /// abolute coordinate data.
    /** By default, offset information is ignored. **/
    virtual void receive(
        Box2D domain, std::vector<char> const &buffer, modif::ModifT kind, Dot2D offset)
    {
        receive(domain, buffer, kind);
    }
    /// Receive data from a byte-stream into the block, and re-map IDs for dynamics if exist.
    virtual void receive(
        Box2D domain, std::vector<char> const &buffer, modif::ModifT kind,
        std::map<int, std::string> const &foreignIds) = 0;
    virtual void attribute(
        Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D const &from,
        modif::ModifT kind) = 0;
    /// Attribute data between two blocks, and adjust coordinates if the block contains absolute
    /// coordinate data.
    virtual void attribute(
        Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D const &from, modif::ModifT kind,
        Dot2D absoluteOffset)
    {
        attribute(toDomain, deltaX, deltaY, from, kind);
    }
};

class AtomicBlock2D : public Block2D {
public:
    AtomicBlock2D(plint nx_, plint ny_);
    AtomicBlock2D(AtomicBlock2D const &rhs);
    virtual ~AtomicBlock2D();
    void swap(AtomicBlock2D &rhs);
    /// Initialize block content by executing internal processors once.
    void initialize();

public:
    /// Get bounding box of the lattice.
    virtual Box2D getBoundingBox() const;
    /// Get number of cells in x-direction.
    plint getNx() const
    {
        return nx;
    }
    /// Get number of cells in y-direction.
    plint getNy() const
    {
        return ny;
    }
    /// Execute all internal dataProcessors at positive or zero level.
    void executeInternalProcessors();
    /// Execute all internal dataProcessors at a given level.
    void executeInternalProcessors(plint level);
    /// Add a dataProcessor, which is executed after each iteration.
    void integrateDataProcessor(DataProcessor2D *processor, plint level);
    /// Remove all data processors.
    void clearDataProcessors();
    /// Remove all data processors with a given ID.
    void removeDataProcessors(int staticId);
    /// Get access to data transfer between blocks.
    virtual BlockDataTransfer2D &getDataTransfer() = 0;
    /// Get access to data transfer between blocks (const version)
    virtual BlockDataTransfer2D const &getDataTransfer() const = 0;
    /// Get an object through which the atomic block can be serialized
    virtual DataSerializer *getBlockSerializer(
        Box2D const &domain, IndexOrdering::OrderingT ordering) const;
    /// Get an object through which the atomic block can be un-serialized
    virtual DataUnSerializer *getBlockUnSerializer(
        Box2D const &domain, IndexOrdering::OrderingT ordering);

public:
    /// Get object to subscribe new internal statistics.
    StatSubscriber &internalStatSubscription();
    /// Copy running statistics to public statistics, and reset running stats.
    void evaluateStatistics();
    /// Get a handle to internal statistics. Don't use this to subscribe new
    /// statistics. Use the method internalStatSubscription() instead.
    BlockStatistics &getInternalStatistics();
    /// Get a constant handle to internal statistics.
    BlockStatistics const &getInternalStatistics() const;

public:
    /// Specify the relative location of the atomic-block in a
    /// global coordinate system.
    void setLocation(Dot2D const &location_);
    /// Get the relative location of the atomic-block in a
    /// global coordinate system.
    Dot2D getLocation() const;
    /// Set the flag value (to be exploited internally by data processors).
    void setFlag(bool value);
    /// Get the flag value (to be exploited internally by data processors).
    bool getFlag() const;

private:
    typedef std::vector<std::vector<DataProcessor2D *> > DataProcessorVector;

private:
    /// Common implementation for explicit/automatic processors.
    void integrateDataProcessor(
        DataProcessor2D *processor, plint level, DataProcessorVector &processors);
    /// Common implementation for explicit/automatic processors.
    void executeInternalProcessors(plint level, DataProcessorVector &processors);
    /// Copy processors from one vector to another.
    void copyDataProcessors(DataProcessorVector const &from, DataProcessorVector &to);
    /// Release memory for a given species of lattice processors.
    void clearDataProcessors(DataProcessorVector &processors);

private:
    plint nx, ny;    /// Dimensions of the lattice.
    Dot2D location;  /// Absolute, real-space coordinate of the local (0,0) position.
    bool flag;       /// Flag for internal use by data processors.
    BlockStatistics internalStatistics;
    StatSubscriber2D statisticsSubscriber;
    DataProcessorVector explicitInternalProcessors;
    DataProcessorVector automaticInternalProcessors;
};

Dot2D computeRelativeDisplacement(AtomicBlock2D const &block1, AtomicBlock2D const &block2);

}  // namespace plb

#endif  // ATOMIC_BLOCK_2D
