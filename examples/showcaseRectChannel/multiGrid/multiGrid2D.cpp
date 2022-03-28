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

#include "multiGrid/multiGrid2D.h"

#include <fstream>
#include <istream>
#include <ostream>

#include "core/plbDebug.h"
#include "core/runTimeDiagnostics.h"
#include "io/serializerIO_2D.h"
#include "multiBlock/multiBlockSerializer2D.h"
#include "multiGrid/multiGridUtil.h"
#include "parallelism/mpiManager.h"

namespace plb {

/* ******************* MultiGridPeriodicity ******************** */

/// Constructors
MultiGridPeriodicitySwitch2D::MultiGridPeriodicitySwitch2D(MultiGrid2D *block_) : block(block_) { }

MultiGridPeriodicitySwitch2D::MultiGridPeriodicitySwitch2D(
    MultiGridPeriodicitySwitch2D const &rhs) :
    block(rhs.block), periodicityArray(rhs.periodicityArray)
{ }

MultiGridPeriodicitySwitch2D &MultiGridPeriodicitySwitch2D::operator=(
    MultiGridPeriodicitySwitch2D const &rhs)
{
    block = rhs.block;
    periodicityArray = rhs.periodicityArray;
    return *this;
}

void MultiGridPeriodicitySwitch2D::swap(MultiGridPeriodicitySwitch2D &rhs)
{
    std::swap(block, rhs.block);
    std::swap(periodicityArray[0], rhs.periodicityArray[0]);
    std::swap(periodicityArray[1], rhs.periodicityArray[1]);
}

Array<bool, 2> const &MultiGridPeriodicitySwitch2D::getPeriodicityArray() const
{
    return periodicityArray;
}

/// toggle on/off the periodicity in one direction
void MultiGridPeriodicitySwitch2D::toggle(int direction, bool periodicity)
{
    PLB_PRECONDITION(direction == 0 || direction == 1);
    periodicityArray[direction] = periodicity;
    block->signalPeriodicity();
}

/// toggle on/off the periodicity in all directions
void MultiGridPeriodicitySwitch2D::toggleAll(bool periodicity)
{
    periodicityArray[0] = periodicity;
    periodicityArray[1] = periodicity;
    block->signalPeriodicity();
}

/* ****************** MultiGridStatsSubscriber2D ******************* */

MultiGridStatSubscriber2D::MultiGridStatSubscriber2D(MultiGrid2D *multiGrid_) :
    multiGrid(multiGrid_)
{ }

MultiGridStatSubscriber2D::MultiGridStatSubscriber2D(MultiGridStatSubscriber2D const &rhs) :
    multiGrid(rhs.multiGrid), dimensionsX(rhs.dimensionsX), dimensionsT(rhs.dimensionsT)

{ }

/** In order to swap two MultiGridStatsSubscriber2D, we need to swap all
 *  the scale information (dx and dt) and the internal statistics of each
 *  block as well as the multigrid itself.
 */
void MultiGridStatSubscriber2D::swap(MultiGridStatSubscriber2D &rhs)
{
    dimensionsX.swap(rhs.dimensionsX);
    dimensionsT.swap(rhs.dimensionsT);
    std::swap(multiGrid, rhs.multiGrid);
}

/// Copy
MultiGridStatSubscriber2D &MultiGridStatSubscriber2D::operator=(
    MultiGridStatSubscriber2D const &rhs)
{
    multiGrid = rhs.multiGrid;  // assign the pointer to the multigrid
    dimensionsX.assign(
        rhs.dimensionsX.begin(), rhs.dimensionsX.end());  // copy all spatial dimensions
    dimensionsT.assign(
        rhs.dimensionsT.begin(), rhs.dimensionsT.end());  // copy all temporal dimensions
    return *this;
}

/// initialize everything ONCE all the multi grid matrices have been created
/** This is necessary because the multiBlockLattice2D creates automatically the
 *  3 statistics. It is not included in the constructor cause the multiBlockLattice2D
 *  might not be allocated.
 */
void MultiGridStatSubscriber2D::initialize()
{
    // Subscribe the automatic quantities from the MultiBlockLattice2D
    multiGrid->getInternalStatistics().subscribeAverage();  // rho
    dimensionsX.push_back(0);
    dimensionsT.push_back(0);

    multiGrid->getInternalStatistics().subscribeAverage();  // uSqr
    dimensionsX.push_back(2);                               // dx^2
    dimensionsT.push_back(-2);                              // dt^-2

    multiGrid->getInternalStatistics().subscribeMax();  // max uSqr
    dimensionsX.push_back(2);
    dimensionsT.push_back(-2);
}

/// Subscribe a new observable for which the average value is computed.
plint MultiGridStatSubscriber2D::subscribeAverage(plint dimDx, plint dimDt)
{
    // we subscribe an average for each level
    for (plint iLevel = 0; iLevel < multiGrid->getNumLevels(); ++iLevel) {
        multiGrid->getComponent(iLevel).internalStatSubscription().subscribeAverage();
    }
    // adding the corresponding rescaling parameters
    dimensionsX.push_back(dimDx);
    dimensionsT.push_back(dimDt);

    // finally subscribe the statistic in the multiGrid internal BlockStatistics
    return multiGrid->getInternalStatistics().subscribeAverage();
}

/// Subscribe a new observable for which the sum is computed.
plint MultiGridStatSubscriber2D::subscribeSum(plint dimDx, plint dimDt)
{
    // we subscribe a sum for each level
    for (plint iLevel = 0; iLevel < multiGrid->getNumLevels(); ++iLevel) {
        multiGrid->getComponent(iLevel).internalStatSubscription().subscribeSum();
    }
    // adding the corresponding rescaling parameters
    dimensionsX.push_back(dimDx);
    dimensionsT.push_back(dimDt);

    // finally subscribe the statistic in the multiGrid internal BlockStatistics
    return multiGrid->getInternalStatistics().subscribeSum();
}

/// Subscribe a new observable for which the maximum is computed.
plint MultiGridStatSubscriber2D::subscribeMax(plint dimDx, plint dimDt)
{
    // we subscribe a max for each level
    for (plint iLevel = 0; iLevel < multiGrid->getNumLevels(); ++iLevel) {
        multiGrid->getComponent(iLevel).internalStatSubscription().subscribeMax();
    }
    // adding the corresponding rescaling parameters
    dimensionsX.push_back(dimDx);
    dimensionsT.push_back(dimDt);

    // finally subscribe the statistic in the multiGrid internal BlockStatistics
    return multiGrid->getInternalStatistics().subscribeMax();
}

/// Subscribe a new integer observable for which the sum is computed.
plint MultiGridStatSubscriber2D::subscribeIntSum(plint dimDx, plint dimDt)
{
    // we subscribe a sum for each level
    for (plint iLevel = 0; iLevel < multiGrid->getNumLevels(); ++iLevel) {
        multiGrid->getComponent(iLevel).internalStatSubscription().subscribeIntSum();
    }
    // adding the corresponding rescaling parameters
    dimensionsX.push_back(dimDx);
    dimensionsT.push_back(dimDt);

    // finally subscribe the statistic in the multiGrid internal BlockStatistics
    return multiGrid->getInternalStatistics().subscribeIntSum();
}

std::vector<int> const &MultiGridStatSubscriber2D::getDimensionsX() const
{
    return dimensionsX;
}
std::vector<int> const &MultiGridStatSubscriber2D::getDimensionsT() const
{
    return dimensionsT;
}

/* *********************** MultiGrid2D ************************* */

MultiGrid2D::MultiGrid2D(MultiGridManagement2D management_, plint behaviorLevel_)

    :
    management(management_),
    behaviorLevel(behaviorLevel_),
    periodicitySwitch(this),
    maxProcessorLevel(-1),
    statisticsOn(true),
    statsSubscriber(this),
    scaleManager(global::getDefaultMultiScaleManager().clone())
{ }

/// Copy constructor
MultiGrid2D::MultiGrid2D(const MultiGrid2D &rhs) :
    Block2D(rhs),
    management(rhs.management),
    behaviorLevel(rhs.behaviorLevel),
    periodicitySwitch(rhs.periodicitySwitch),
    maxProcessorLevel(rhs.maxProcessorLevel),
    statisticsOn(rhs.statisticsOn),
    statsSubscriber(this),
    scaleManager(rhs.scaleManager->clone())
{ }

MultiGrid2D::MultiGrid2D(MultiGrid2D const &rhs, Box2D subDomain, bool crop) :
    Block2D(rhs),
    management(extractManagement(rhs.getMultiGridManagement(), subDomain, crop)),
    behaviorLevel(rhs.behaviorLevel),
    periodicitySwitch(rhs.periodicitySwitch),
    maxProcessorLevel(rhs.maxProcessorLevel),
    statisticsOn(rhs.statisticsOn),
    statsSubscriber(this),
    scaleManager(rhs.scaleManager->clone())
{ }

MultiGrid2D::~MultiGrid2D()
{
    // we only delete the scaleManager
    delete scaleManager;
    // the internal combined statistics are deleted by each MultiBlock2D
}

// Swapping two MultiGrid2D
void MultiGrid2D::swap(MultiGrid2D &rhs)
{
    management.swap(rhs.management);
    std::swap(behaviorLevel, rhs.behaviorLevel);

    periodicitySwitch.swap(rhs.periodicitySwitch);

    std::swap(maxProcessorLevel, rhs.maxProcessorLevel);
    std::swap(statisticsOn, rhs.statisticsOn);
    statsSubscriber.swap(rhs.statsSubscriber);
}

/// Retrieving the MultiGridManagement2D
MultiGridManagement2D const &MultiGrid2D::getMultiGridManagement() const
{
    return management;
}
MultiGridManagement2D &MultiGrid2D::getMultiGridManagement()
{
    return management;
}

/// Execute all processors one
void MultiGrid2D::initialize()
{
    executeInternalProcessors();
}

void MultiGrid2D::subscribeProcessor(plint level)
{
    maxProcessorLevel = std::max(level, maxProcessorLevel);
}

/// Execute data processors
void MultiGrid2D::executeInternalProcessors()
{
    // for each MultiBlock in the multigrid
    for (int iLevel = 0; iLevel < this->getNumLevels(); ++iLevel) {
        getComponent(iLevel).executeInternalProcessors();
    }
}

/// Execute data processors at a certain level
void MultiGrid2D::executeInternalProcessors(plint pLevel)
{
    // for each MultiBlock in the multigrid execute the processors of a given level
    for (int iLevel = 0; iLevel < this->getNumLevels(); ++iLevel) {
        getComponent(iLevel).executeInternalProcessors(pLevel);
    }
}

/// Retrieve the reference level
plint MultiGrid2D::getReferenceLevel() const
{
    return management.getReferenceLevel();
}

plint MultiGrid2D::getNumLevels() const
{
    return management.getNumLevels();
}

plint MultiGrid2D::getBehaviorLevel() const
{
    return behaviorLevel;
}

/** Get the object charged for periodicity (non-modifiable).
 *  Note that periodicity is affected at the reference level only. If you
 *  want your simulation to be overall periodic, it is important that
 *  the reference block reaches all borders.
 **/
MultiGridPeriodicitySwitch2D const &MultiGrid2D::periodicity() const
{
    return periodicitySwitch;
}

/** Get the object charged for periodicity (modifiable).
 *  Note that periodicity is affected at the reference level only. If you
 *  want your simulation to be overall periodic, it is important that
 *  the reference block reaches all borders.
 **/
MultiGridPeriodicitySwitch2D &MultiGrid2D::periodicity()
{
    return periodicitySwitch;
}

/// Retrieve the multigrid statistics (this contains stats of all the levels )
BlockStatistics &MultiGrid2D::getInternalStatistics()
{
    return internalStatistics;
}

BlockStatistics const &MultiGrid2D::getInternalStatistics() const
{
    return internalStatistics;
}

MultiGridStatSubscriber2D &MultiGrid2D::internalStatSubscription()
{
    return statsSubscriber;
}

plint MultiGrid2D::getNx() const
{
    return management.getBoundingBox(behaviorLevel).getNx();
}

plint MultiGrid2D::getNy() const
{
    return management.getBoundingBox(behaviorLevel).getNy();
}

Box2D MultiGrid2D::getBoundingBox() const
{
    return management.getBoundingBox(behaviorLevel);
}

void MultiGrid2D::setBehaviorLevel(plint behaviorLevel_)
{
    behaviorLevel = behaviorLevel_;
}

/// This method rescales each statistic and combines all the results
void MultiGrid2D::reduceStatistics()
{
    std::vector<BlockStatistics *> levelStatistics(this->getNumLevels());
    std::vector<int> dimensionsX, dimensionsT;

    dimensionsX = statsSubscriber.getDimensionsX();
    dimensionsT = statsSubscriber.getDimensionsT();

    for (plint iLevel = 0; iLevel < this->getNumLevels(); ++iLevel) {
        int dxScale = this->getReferenceLevel() - iLevel;
        int dtScale = dxScale;  // TODO: here, we assume convective scaling; general case should be
                                // considered.
        std::vector<double> scales(dimensionsX.size());
        for (pluint iScale = 0; iScale < scales.size(); ++iScale) {
            scales[iScale] =
                scaleToReference(dxScale, dimensionsX[iScale], dtScale, dimensionsT[iScale]);
        }
        // copy the statistics
        levelStatistics[iLevel] = new BlockStatistics(getComponent(iLevel).getInternalStatistics());
        // rescale the statistics to the reference level
        levelStatistics[iLevel]->rescale(scales);
    }
    combine(levelStatistics, getInternalStatistics());

    internalStatistics.incrementStats();

    for (plint iLevel = 0; iLevel < this->getNumLevels(); ++iLevel) {
        delete levelStatistics[iLevel];
    }
}

/// Compute statistics. For each level, first compute their statistics, then rescale and join
/// everything
void MultiGrid2D::evaluateStatistics()
{
    // join the stats
    if (statisticsOn) {
        reduceStatistics();
    }
}

/** Simply use the periodicity of the reference level lattice which should at least
 *  include all the borders of the simulation
 */
void MultiGrid2D::signalPeriodicity()
{
    Array<bool, 2> periodicity = periodicitySwitch.getPeriodicityArray();
    getComponent(this->getReferenceLevel()).periodicity().toggle(0, periodicity[0]);
    getComponent(this->getReferenceLevel()).periodicity().toggle(1, periodicity[1]);
    getComponent(this->getReferenceLevel()).signalPeriodicity();
}

void MultiGrid2D::toggleInternalStatistics(bool statisticsOn_)
{
    statisticsOn = statisticsOn_;
}

bool MultiGrid2D::isInternalStatisticsOn() const
{
    return statisticsOn;
}

MultiScaleManager const &MultiGrid2D::getScaleManager() const
{
    return *scaleManager;
}

/// I/O
DataSerializer *MultiGrid2D::getBlockSerializer(
    Box2D const &domain, IndexOrdering::OrderingT ordering) const
{
    PLB_ASSERT(false);
    return 0;
}

DataUnSerializer *MultiGrid2D::getBlockUnSerializer(
    Box2D const &domain, IndexOrdering::OrderingT ordering)
{
    PLB_ASSERT(false);
    return 0;
}

/// Checkpointing
void saveBinaryGrid(MultiGrid2D const &block, std::string fName, bool enforceUint)
{
    std::ofstream *ostr = 0;
    bool isOK = true;
    if (global::mpi().isMainProcessor()) {
        ostr = new std::ofstream(fName.c_str());
        isOK = (bool)(*ostr);
    }
    plbMainProcIOError(
        !isOK, std::string("Could not open binary file ") + fName + std::string(" for saving"));
    for (plint iLevel = 0; iLevel < block.getNumLevels(); ++iLevel) {
        serializerToBase64Stream(
            block.getComponent(iLevel).getBlockSerializer(
                block.getComponent(iLevel).getBoundingBox(),
                global::IOpolicy().getIndexOrderingForStreams()),
            ostr, enforceUint);
    }

    delete ostr;
}

void loadBinaryGrid(MultiGrid2D &block, std::string fName, bool enforceUint)
{
    std::ifstream *istr = 0;
    bool isOK = true;
    if (global::mpi().isMainProcessor()) {
        istr = new std::ifstream(fName.c_str());
        isOK = (bool)(*istr);
    }
    plbMainProcIOError(
        !isOK, std::string("Could not open binary file ") + fName + std::string(" for reading"));
    for (plint iLevel = 0; iLevel < block.getNumLevels(); ++iLevel) {
        base64StreamToUnSerializer(
            istr,
            block.getComponent(iLevel).getBlockUnSerializer(
                block.getComponent(iLevel).getBoundingBox(),
                global::IOpolicy().getIndexOrderingForStreams()),
            enforceUint);
    }

    delete istr;
}

}  // namespace plb
