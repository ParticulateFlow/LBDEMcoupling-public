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

#ifndef MULTI_GRID_2D_H
#define MULTI_GRID_2D_H

#include "core/array.h"
#include "core/block2D.h"
#include "core/blockStatistics.h"
#include "core/geometry2D.h"
#include "core/globalDefs.h"
#include "multiBlock/blockCommunicator2D.h"
#include "multiBlock/combinedStatistics.h"
#include "multiBlock/multiBlock2D.h"
#include "multiBlock/multiBlockManagement2D.h"
#include "multiGrid/multiGridManagement2D.h"
#include "multiGrid/multiScale.h"

namespace plb {

class MultiGrid2D;

/// Class that allows periodicity to be toggled on/off in the multigrid
class MultiGridPeriodicitySwitch2D {
public:
    MultiGridPeriodicitySwitch2D(MultiGrid2D *block_);
    MultiGridPeriodicitySwitch2D(MultiGridPeriodicitySwitch2D const &rhs);
    MultiGridPeriodicitySwitch2D &operator=(MultiGridPeriodicitySwitch2D const &rhs);
    void swap(MultiGridPeriodicitySwitch2D &rhs);

    void toggle(int direction, bool periodicity);
    void toggleAll(bool periodicity);

    Array<bool, 2> const &getPeriodicityArray() const;

private:
    MultiGrid2D *block;
    Array<bool, 2> periodicityArray;
};

/// Handles statistics subscriptions for the MultiGridLattice2D
/** This class will provide the same functionalities as the MultiStatSubscriber2D, but
 *  each subscription comes with two parameters which will allow the rescaling of quantities
 *  among grids of different resolution.
 */
class MultiGridStatSubscriber2D {
public:
    MultiGridStatSubscriber2D(MultiGrid2D *multiGrid_);
    MultiGridStatSubscriber2D(MultiGridStatSubscriber2D const &rhs);

    void swap(MultiGridStatSubscriber2D &rhs);
    MultiGridStatSubscriber2D &operator=(MultiGridStatSubscriber2D const &rhs);

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
    MultiGrid2D *multiGrid;
    std::vector<int> dimensionsX;
    std::vector<int> dimensionsT;
};

/// Base non-typed object that represents a multigrid
class MultiGrid2D : public Block2D {
public:
    MultiGrid2D(MultiGridManagement2D management, plint behaviorLevel_);
    // Copy constructor
    MultiGrid2D(const MultiGrid2D &rhs);
    MultiGrid2D(MultiGrid2D const &rhs, Box2D subDomain, bool crop);
    virtual ~MultiGrid2D();

    void swap(MultiGrid2D &rhs);

    /// Retrieving the components of the multigrid
    virtual MultiBlock2D const &getComponent(plint level) const = 0;
    virtual MultiBlock2D &getComponent(plint level) = 0;

    MultiGridManagement2D const &getMultiGridManagement() const;
    MultiGridManagement2D &getMultiGridManagement();

    /// "Sizes" of the multigrid, according to a definition made by the end-user class.
    plint getNx() const;
    plint getNy() const;

    /// Retrieve the bounding box of the domain according to a definition made by user
    virtual Box2D getBoundingBox() const;

    /// Retrieve the multigrid informations (for the implementing classes)
    plint getReferenceLevel() const;
    plint getNumLevels() const;
    plint getBehaviorLevel() const;
    void setBehaviorLevel(plint behaviorLevel_);

    /// Execute all processors one
    virtual void initialize();

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
    MultiGridPeriodicitySwitch2D const &periodicity() const;
    /// Periodicity control
    MultiGridPeriodicitySwitch2D &periodicity();
    void signalPeriodicity();

    /// Retrieve the stats subscriber
    MultiGridStatSubscriber2D &internalStatSubscription();

    /// Retrieve the scale manager
    MultiScaleManager const &getScaleManager() const;

    virtual DataSerializer *getBlockSerializer(
        Box2D const &domain, IndexOrdering::OrderingT ordering) const;
    virtual DataUnSerializer *getBlockUnSerializer(
        Box2D const &domain, IndexOrdering::OrderingT ordering);

    virtual int getBlockId() const = 0;

private:
    // reference level management
    MultiGridManagement2D management;

    plint behaviorLevel;
    MultiGridPeriodicitySwitch2D periodicitySwitch;

    plint maxProcessorLevel;
    bool statisticsOn;

    MultiGridStatSubscriber2D statsSubscriber;
    BlockStatistics internalStatistics;

    // Used to easen the scale conversion in this class
    MultiScaleManager *scaleManager;
};

// checkpointing for the multiGrid
void saveBinaryGrid(MultiGrid2D const &block, std::string fName, bool enforceUint = false);
void loadBinaryGrid(MultiGrid2D &block, std::string fName, bool enforceUint = false);

}  // namespace plb

#endif  // MULTI_GRID_2D_H
