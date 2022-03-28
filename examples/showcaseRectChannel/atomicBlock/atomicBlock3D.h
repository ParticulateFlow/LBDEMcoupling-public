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

#ifndef ATOMIC_BLOCK_3D_H
#define ATOMIC_BLOCK_3D_H

#include <algorithm>
#include <map>
#include <string>

#include "atomicBlock/dataProcessorWrapper3D.h"
#include "atomicBlock/reductiveDataProcessorWrapper3D.h"
#include "core/block3D.h"
#include "core/blockIdentifiers.h"
#include "core/blockStatistics.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"

namespace plb {

// Forward declarations
class AtomicBlock3D;

/// Handles statistics subscriptions for the BlockLattice3D
class StatSubscriber3D : public StatSubscriber {
public:
    StatSubscriber3D(AtomicBlock3D &block_);
    /// Subscribe a new observable for which the average value is computed.
    virtual plint subscribeAverage();
    /// Subscribe a new observable for which the sum is computed.
    virtual plint subscribeSum();
    /// Subscribe a new observable for which the maximum is computed.
    virtual plint subscribeMax();
    /// Subscribe a new integer observable for which the sum is computed.
    virtual plint subscribeIntSum();

private:
    AtomicBlock3D &block;
};

struct BlockDataTransfer3D {
    virtual ~BlockDataTransfer3D() { }
    virtual void setBlock(AtomicBlock3D &block_) = 0;
    virtual void setConstBlock(AtomicBlock3D const &block_) = 0;
    virtual BlockDataTransfer3D *clone() const = 0;
    virtual plint staticCellSize() const = 0;
    /// Send data from the block into a byte-stream.
    /** This method automatically resizes the buffer so it holds the proper amount
     *  of data. This is particularly important if the transmitted data is dynamic,
     *  in which case the size of the buffer cannot be predicted exactly. You can
     *  always preallocate data for the buffer though (through resize or reserve),
     *  to avoid reallocation and improve performance.
     **/
    virtual void send(Box3D domain, std::vector<char> &buffer, modif::ModifT kind) const = 0;
    /// Receive data from a byte-stream into the block.
    virtual void receive(Box3D domain, std::vector<char> const &buffer, modif::ModifT kind) = 0;
    /// Receive data from a byte-stream into the block, and adjust coordinates if the block contains
    /// abolute coordinate data.
    /** By default, offset information is ignored. **/
    virtual void receive(
        Box3D domain, std::vector<char> const &buffer, modif::ModifT kind, Dot3D absoluteOffset)
    {
        receive(domain, buffer, kind);
    }
    /// Receive data from a byte-stream into the block, and re-map IDs for dynamics if exist.
    virtual void receive(
        Box3D domain, std::vector<char> const &buffer, modif::ModifT kind,
        std::map<int, std::string> const &foreignIds) = 0;
    /// Attribute data between two blocks.
    virtual void attribute(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
        modif::ModifT kind) = 0;
    /// Attribute data between two blocks, and adjust coordinates if the block contains absolute
    /// coordinate data.
    virtual void attribute(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
        modif::ModifT kind, Dot3D absoluteOffset)
    {
        attribute(toDomain, deltaX, deltaY, deltaZ, from, kind);
    }
};

class AtomicBlock3D : public Block3D {
public:
    AtomicBlock3D(plint nx_, plint ny_, plint nz_, BlockDataTransfer3D *defaultDataTransfer);
    AtomicBlock3D(AtomicBlock3D const &rhs);
    AtomicBlock3D(AtomicBlock3D const &rhs, BlockDataTransfer3D *defaultDataTransfer);
    virtual ~AtomicBlock3D();
    void swap(AtomicBlock3D &rhs);
    /// Initialize block content by executing internal processors once.
    void initialize();

public:
    /// Get bounding box of the lattice.
    virtual Box3D getBoundingBox() const;
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
    /// Get number of cells in z-direction.
    plint getNz() const
    {
        return nz;
    }
    /// Execute all internal dataProcessors at positive or zero level.
    void executeInternalProcessors();
    /// Execute all internal dataProcessors at a given level.
    void executeInternalProcessors(plint level);
    /// Add a dataProcessor, which is executed after each iteration.
    void integrateDataProcessor(DataProcessor3D *processor, plint level);
    /// Remove all data processors.
    void clearDataProcessors();
    /// Remove all data processors with a given ID.
    void removeDataProcessors(int staticId);
    /// Reset data transfer policy to a new one.
    void setDataTransfer(BlockDataTransfer3D *newDataTransfer);
    /// Get access to data transfer between blocks.
    BlockDataTransfer3D &getDataTransfer();
    /// Get access to data transfer between blocks (const version)
    BlockDataTransfer3D const &getDataTransfer() const;
    /// Get an object through which the atomic block can be serialized
    virtual DataSerializer *getBlockSerializer(
        Box3D const &domain, IndexOrdering::OrderingT ordering) const;
    /// Get an object through which the atomic block can be un-serialized
    virtual DataUnSerializer *getBlockUnSerializer(
        Box3D const &domain, IndexOrdering::OrderingT ordering);

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
    void setLocation(Dot3D const &location_);
    /// Get the relative location of the atomic-block in a
    /// global coordinate system.
    Dot3D getLocation() const;
    /// Set the flag value (to be exploited internally by data processors).
    void setFlag(bool value);
    /// Get the flag value (to be exploited internally by data processors).
    bool getFlag() const;

private:
    typedef std::vector<std::vector<DataProcessor3D *> > DataProcessorVector;

private:
    /// Common implementation for explicit/automatic processors.
    void integrateDataProcessor(
        DataProcessor3D *processor, plint level, DataProcessorVector &processors);
    /// Common implementation for explicit/automatic processors.
    void executeInternalProcessors(plint level, DataProcessorVector &processors);
    /// Copy processors from one vector to another.
    void copyDataProcessors(DataProcessorVector const &from, DataProcessorVector &to);
    /// Release memory for a given species of lattice processors.
    void clearDataProcessors(DataProcessorVector &processors);

private:
    plint nx, ny, nz;  /// Dimensions of the lattice.
    Dot3D location;    /// Absolute, real-space coordinate of the local (0,0,0) position.
    bool flag;         /// Flag for internal use by data processors.
    BlockStatistics internalStatistics;
    StatSubscriber3D statisticsSubscriber;
    DataProcessorVector explicitInternalProcessors;
    DataProcessorVector automaticInternalProcessors;
    mutable BlockDataTransfer3D *dataTransfer;
};

Dot3D computeRelativeDisplacement(AtomicBlock3D const &block1, AtomicBlock3D const &block2);

}  // namespace plb

#endif  // ATOMIC_BLOCK_3D
