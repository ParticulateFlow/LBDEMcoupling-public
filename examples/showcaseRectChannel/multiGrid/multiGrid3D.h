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

#ifndef MULTI_GRID_3D_H
#define MULTI_GRID_3D_H

#include "core/array.h"
#include "core/block3D.h"
#include "core/blockStatistics.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "multiBlock/blockCommunicator3D.h"
#include "multiBlock/combinedStatistics.h"
#include "multiBlock/multiBlock3D.h"
#include "multiBlock/multiBlockManagement3D.h"
#include "multiGrid/multiGridManagement3D.h"
#include "multiGrid/multiScale.h"

namespace plb {

class MultiGrid3D;

/// Class that allows periodicity to be toggled on/off in the multigrid
class MultiGridPeriodicitySwitch3D {
public:
    MultiGridPeriodicitySwitch3D(MultiGrid3D *block_);
    MultiGridPeriodicitySwitch3D(MultiGridPeriodicitySwitch3D const &rhs);
    MultiGridPeriodicitySwitch3D &operator=(MultiGridPeriodicitySwitch3D const &rhs);

    void swap(MultiGridPeriodicitySwitch3D &rhs);

    void toggle(int direction, bool periodicity);
    void toggleAll(bool periodicity);

    Array<bool, 3> const &getPeriodicityArray() const;

private:
    MultiGrid3D *block;
    Array<bool, 3> periodicityArray;
};

/// Handles statistics subscriptions for the MultiGridLattice3D
/** This class will provide the same functionalities as the MultiStatSubscriber3D, but
 *  each subscription comes with two parameters which will allow the rescaling of quantities
 *  among grids of different resolution.
 */
class MultiGridStatSubscriber3D {
public:
    MultiGridStatSubscriber3D(MultiGrid3D *multiGrid_);
    MultiGridStatSubscriber3D(MultiGridStatSubscriber3D const &rhs);

    void swap(MultiGridStatSubscriber3D &rhs);
    MultiGridStatSubscriber3D &operator=(MultiGridStatSubscriber3D const &rhs);

    /// Subscribe a new observable for which the average value is computed.
    plint subscribeAverage(plint dimDx, plint dimDt);
    /// Subscribe a new observable for which the sum is computed.
    plint subscribeSum(plint dimDx, plint dimDt);
    /// Subscribe a new observable for which the maximum is computed.
    plint subscribeMax(plint dimDx, plint dimDt);
    /// Subscribe a new integer observable for which the sum is computed.
    plint subscribeIntSum(plint dimDx, plint dimDt);
    /// Initialize the default subscriptions (rho, uSqr, et max uSqr)
    void initialize();

    /// Retrieve the dimensions of dx and dt for all the components
    std::vector<int> const &getDimensionsX() const;
    std::vector<int> const &getDimensionsT() const;

private:
    MultiGrid3D *multiGrid;
    std::vector<int> dimensionsX;
    std::vector<int> dimensionsT;
};

/// Base non-typed object that represents a multigrid
class MultiGrid3D : public Block3D {
public:
    MultiGrid3D(MultiGridManagement3D management, plint behaviorLevel_);
    // Copy constructor
    MultiGrid3D(const MultiGrid3D &rhs);
    MultiGrid3D(MultiGrid3D const &rhs, Box3D subDomain, bool crop);
    virtual ~MultiGrid3D();

    void swap(MultiGrid3D &rhs);

    /// Retrieving the components of the multigrid
    virtual MultiBlock3D const &getComponent(plint level) const = 0;
    virtual MultiBlock3D &getComponent(plint level) = 0;

    MultiGridManagement3D const &getMultiGridManagement() const;
    MultiGridManagement3D &getMultiGridManagement();

    /// "Sizes" of the multigrid, according to a definition made by the end-user class.
    plint getNx() const;
    plint getNy() const;
    plint getNz() const;

    /// Retrieve the bounding box of the domain according to a definition made by user
    virtual Box3D getBoundingBox() const;

    /// Retrieve the multigrid informations (for the implementing classes)
    plint getReferenceLevel() const;
    plint getNumLevels() const;
    plint getBehaviorLevel() const;
    void setBehaviorLevel(plint behaviorLevel_);

    /// Execute all processors one
    void initialize();

    /// Execute data processors
    void executeInternalProcessors();
    void executeInternalProcessors(plint level);

    /// Subscription of a Data Processor
    void subscribeProcessor(plint level);

    /// Retrieve the multigrid statistics (this contains all the statistics rescaled)
    BlockStatistics &getInternalStatistics();
    BlockStatistics const &getInternalStatistics() const;

    /// statistics related manipulations like evaluation
    void reduceStatistics();
    void evaluateStatistics();
    void toggleInternalStatistics(bool statisticsOn_);
    bool isInternalStatisticsOn() const;

    /// Periodicity control
    MultiGridPeriodicitySwitch3D const &periodicity() const;
    /// Periodicity control
    MultiGridPeriodicitySwitch3D &periodicity();
    void signalPeriodicity();

    /// Retrieve the stats subscriber
    MultiGridStatSubscriber3D &internalStatSubscription();

    /// Retrieve the scale manager
    MultiScaleManager const &getScaleManager() const;

    virtual DataSerializer *getBlockSerializer(
        Box3D const &domain, IndexOrdering::OrderingT ordering) const;
    virtual DataUnSerializer *getBlockUnSerializer(
        Box3D const &domain, IndexOrdering::OrderingT ordering);

    virtual int getBlockId() const = 0;

private:
    // reference level management
    MultiGridManagement3D management;

    plint behaviorLevel;
    MultiGridPeriodicitySwitch3D periodicitySwitch;

    plint maxProcessorLevel;
    bool statisticsOn;

    MultiGridStatSubscriber3D statsSubscriber;
    BlockStatistics internalStatistics;

    // Used to easen the scale conversion in this class
    MultiScaleManager *scaleManager;
};

// checkpointing for the multiGrid
void saveBinaryGrid(MultiGrid3D const &block, std::string fName, bool enforceUint = false);
void loadBinaryGrid(MultiGrid3D &block, std::string fName, bool enforceUint = false);

}  // namespace plb

#endif  // MULTI_GRID_3D_H
