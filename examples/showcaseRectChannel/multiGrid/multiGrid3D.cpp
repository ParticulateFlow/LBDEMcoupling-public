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

#include "multiGrid/multiGrid3D.h"

#include <fstream>
#include <istream>
#include <ostream>

#include "core/plbDebug.h"
#include "core/runTimeDiagnostics.h"
#include "io/serializerIO_3D.h"
#include "multiBlock/multiBlockSerializer3D.h"
#include "multiGrid/multiGridUtil.h"
#include "parallelism/mpiManager.h"

// Who deletes the communicators?
// How can you swap statsSubscriber if it should point to *this?

namespace plb {

/* ******************* MultiGridPeriodicity ******************** */

/// Constructors
MultiGridPeriodicitySwitch3D::MultiGridPeriodicitySwitch3D(MultiGrid3D *block_) :
    block(block_), periodicityArray(Array<bool, 3>(false, false, false))
{ }

MultiGridPeriodicitySwitch3D::MultiGridPeriodicitySwitch3D(
    MultiGridPeriodicitySwitch3D const &rhs) :
    block(rhs.block), periodicityArray(rhs.periodicityArray)
{ }

MultiGridPeriodicitySwitch3D &MultiGridPeriodicitySwitch3D::operator=(
    MultiGridPeriodicitySwitch3D const &rhs)
{
    block = rhs.block;
    periodicityArray = rhs.periodicityArray;
    return *this;
}

void MultiGridPeriodicitySwitch3D::swap(MultiGridPeriodicitySwitch3D &rhs)
{
    std::swap(block, rhs.block);
    std::swap(periodicityArray[0], rhs.periodicityArray[0]);
    std::swap(periodicityArray[1], rhs.periodicityArray[1]);
    std::swap(periodicityArray[2], rhs.periodicityArray[2]);
}

Array<bool, 3> const &MultiGridPeriodicitySwitch3D::getPeriodicityArray() const
{
    return periodicityArray;
}

/// toggle on/off the periodicity in one direction
void MultiGridPeriodicitySwitch3D::toggle(int direction, bool periodicity)
{
    PLB_PRECONDITION(direction == 0 || direction == 1 || direction == 2);
    periodicityArray[direction] = periodicity;
    block->signalPeriodicity();
}

/// toggle on/off the periodicity in all directions
void MultiGridPeriodicitySwitch3D::toggleAll(bool periodicity)
{
    periodicityArray[0] = periodicity;
    periodicityArray[1] = periodicity;
    periodicityArray[2] = periodicity;
    block->signalPeriodicity();
}

/* ****************** MultiGridStatsSubscriber ******************* */

MultiGridStatSubscriber3D::MultiGridStatSubscriber3D(MultiGrid3D *multiGrid_) :
    multiGrid(multiGrid_)
{ }

MultiGridStatSubscriber3D::MultiGridStatSubscriber3D(MultiGridStatSubscriber3D const &rhs) :
    multiGrid(rhs.multiGrid), dimensionsX(rhs.dimensionsX), dimensionsT(rhs.dimensionsT)
{ }

/** In order to swap two MultiGridStatsSubscriber3D, we need to swap all
 *  the scale information (dx and dt) and the internal statistics of each
 *  block. The reference to each block does not change.
 */
void MultiGridStatSubscriber3D::swap(MultiGridStatSubscriber3D &rhs)
{
    dimensionsX.swap(rhs.dimensionsX);
    dimensionsT.swap(rhs.dimensionsT);
    std::swap(multiGrid, rhs.multiGrid);
}

/// Copy
MultiGridStatSubscriber3D &MultiGridStatSubscriber3D::operator=(
    MultiGridStatSubscriber3D const &rhs)
{
    multiGrid = rhs.multiGrid;  // assign the pointer to the multigrid
    dimensionsX.assign(
        rhs.dimensionsX.begin(), rhs.dimensionsX.end());  // copy all spatial dimensions
    dimensionsT.assign(
        rhs.dimensionsT.begin(), rhs.dimensionsT.end());  // copy all temporal dimensions
    return *this;
}

/// initialize everything ONCE all the multi grid matrices have been created
/** This is necessary because the multiBlockLattice3D creates automatically the
 *  3 statistics. It is not included in the constructor cause the multiBlockLattice3D
 *  might not be allocated.
 */
void MultiGridStatSubscriber3D::initialize()
{
    // Subscribe the automatic quantities from the MultiBlockLattice3D
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
plint MultiGridStatSubscriber3D::subscribeAverage(plint dimDx, plint dimDt)
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
plint MultiGridStatSubscriber3D::subscribeSum(plint dimDx, plint dimDt)
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
plint MultiGridStatSubscriber3D::subscribeMax(plint dimDx, plint dimDt)
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
plint MultiGridStatSubscriber3D::subscribeIntSum(plint dimDx, plint dimDt)
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

std::vector<int> const &MultiGridStatSubscriber3D::getDimensionsX() const
{
    return dimensionsX;
}
std::vector<int> const &MultiGridStatSubscriber3D::getDimensionsT() const
{
    return dimensionsT;
}

/* *********************** MultiGrid3D ************************* */

// Constructor
MultiGrid3D::MultiGrid3D(MultiGridManagement3D management_, plint behaviorLevel_)

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
MultiGrid3D::MultiGrid3D(const MultiGrid3D &rhs) :
    Block3D(rhs),
    management(rhs.management),
    behaviorLevel(rhs.behaviorLevel),
    periodicitySwitch(rhs.periodicitySwitch),
    maxProcessorLevel(rhs.maxProcessorLevel),
    statisticsOn(rhs.statisticsOn),
    statsSubscriber(this),
    scaleManager(rhs.scaleManager->clone())
{ }

MultiGrid3D::MultiGrid3D(MultiGrid3D const &rhs, Box3D subDomain, bool crop) :
    Block3D(rhs),
    management(extractManagement(rhs.getMultiGridManagement(), subDomain, crop)),
    // management(rhs.getMultiGridManagement()),
    behaviorLevel(rhs.behaviorLevel),
    periodicitySwitch(rhs.periodicitySwitch),
    maxProcessorLevel(rhs.maxProcessorLevel),
    statisticsOn(rhs.statisticsOn),
    statsSubscriber(this),  // TODO
    scaleManager(rhs.scaleManager->clone())
{ }

MultiGrid3D::~MultiGrid3D()
{
    // we only delete the scaleManager
    delete scaleManager;
}

// Swapping two MultiGrid3D TODO
void MultiGrid3D::swap(MultiGrid3D &rhs)
{
    management.swap(rhs.management);
    std::swap(behaviorLevel, rhs.behaviorLevel);

    // Swap cannot be used because Periodicity switch
    // holds a reference to the block it is modifying.
    MultiGridPeriodicitySwitch3D tmp = periodicitySwitch;
    periodicitySwitch = rhs.periodicitySwitch;
    rhs.periodicitySwitch = tmp;

    std::swap(maxProcessorLevel, rhs.maxProcessorLevel);
    std::swap(statisticsOn, rhs.statisticsOn);
    statsSubscriber.swap(rhs.statsSubscriber);
}

/// Retrieving the MultiGridManagement3D
MultiGridManagement3D const &MultiGrid3D::getMultiGridManagement() const
{
    return management;
}
MultiGridManagement3D &MultiGrid3D::getMultiGridManagement()
{
    return management;
}

/// Execute all processors one
void MultiGrid3D::initialize()
{
    executeInternalProcessors();
}

void MultiGrid3D::subscribeProcessor(plint level)
{
    maxProcessorLevel = std::max(level, maxProcessorLevel);
}

/// Execute data processors
void MultiGrid3D::executeInternalProcessors()
{
    // for each MultiBlock in the multigrid
    for (int iLevel = 0; iLevel < this->getNumLevels(); ++iLevel) {
        getComponent(iLevel).executeInternalProcessors();
    }
}

/// Execute data processors at a certain level
void MultiGrid3D::executeInternalProcessors(plint pLevel)
{
    // for each MultiBlock in the multigrid execute the processors of a given level
    for (int iLevel = 0; iLevel < this->getNumLevels(); ++iLevel) {
        getComponent(iLevel).executeInternalProcessors(pLevel);
    }
}

/// Retrieve the reference level
plint MultiGrid3D::getReferenceLevel() const
{
    return management.getReferenceLevel();
}

plint MultiGrid3D::getNumLevels() const
{
    return management.getNumLevels();
}

plint MultiGrid3D::getBehaviorLevel() const
{
    return behaviorLevel;
}

/** Get the object charged for periodicity (non-modifiable).
 *  Note that periodicity is affected at the reference level only. If you
 *  want your simulation to be overall periodic, it is important that
 *  the reference block reaches all borders.
 **/
MultiGridPeriodicitySwitch3D const &MultiGrid3D::periodicity() const
{
    return periodicitySwitch;
}

/** Get the object charged for periodicity (modifiable).
 *  Note that periodicity is affected at the reference level only. If you
 *  want your simulation to be overall periodic, it is important that
 *  the reference block reaches all borders.
 **/
MultiGridPeriodicitySwitch3D &MultiGrid3D::periodicity()
{
    return periodicitySwitch;
}

/// Retrieve the multigrid statistics (this contains all the level statistics)
BlockStatistics &MultiGrid3D::getInternalStatistics()
{
    return internalStatistics;
}

BlockStatistics const &MultiGrid3D::getInternalStatistics() const
{
    return internalStatistics;
}

MultiGridStatSubscriber3D &MultiGrid3D::internalStatSubscription()
{
    return statsSubscriber;
}

plint MultiGrid3D::getNx() const
{
    return management.getBoundingBox(behaviorLevel).getNx();
}

plint MultiGrid3D::getNy() const
{
    return management.getBoundingBox(behaviorLevel).getNy();
}

plint MultiGrid3D::getNz() const
{
    return management.getBoundingBox(behaviorLevel).getNz();
}

Box3D MultiGrid3D::getBoundingBox() const
{
    return management.getBoundingBox(behaviorLevel);
}

void MultiGrid3D::setBehaviorLevel(plint behaviorLevel_)
{
    behaviorLevel = behaviorLevel_;
}

/// This method rescales each statistic and combines all the results
void MultiGrid3D::reduceStatistics()
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
        levelStatistics[iLevel] =
            new BlockStatistics(getComponent(iLevel).getInternalStatistics());  // copy
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
void MultiGrid3D::evaluateStatistics()
{
    // join the stats
    if (statisticsOn) {
        reduceStatistics();
    }
}

/** Simply use the periodicity of every multiBlock representing each level
 */
void MultiGrid3D::signalPeriodicity()
{
    Array<bool, 3> periodicity = periodicitySwitch.getPeriodicityArray();
    for (plint iLevel = 0; iLevel < this->getNumLevels(); ++iLevel) {
        getComponent(iLevel).periodicity().toggle(0, periodicity[0]);
        getComponent(iLevel).periodicity().toggle(1, periodicity[1]);
        getComponent(iLevel).periodicity().toggle(2, periodicity[2]);
        getComponent(iLevel).signalPeriodicity();
    }
}

void MultiGrid3D::toggleInternalStatistics(bool statisticsOn_)
{
    statisticsOn = statisticsOn_;
}

bool MultiGrid3D::isInternalStatisticsOn() const
{
    return statisticsOn;
}

MultiScaleManager const &MultiGrid3D::getScaleManager() const
{
    return *scaleManager;
}

// TODO think how to do a good serialization
/// I/O
DataSerializer *MultiGrid3D::getBlockSerializer(
    Box3D const &domain, IndexOrdering::OrderingT ordering) const
{
    PLB_ASSERT(false);
    return 0;
}

DataUnSerializer *MultiGrid3D::getBlockUnSerializer(
    Box3D const &domain, IndexOrdering::OrderingT ordering)
{
    PLB_ASSERT(false);
    return 0;
}

/// Checkpointing
void saveBinaryGrid(MultiGrid3D const &block, std::string fName, bool enforceUint)
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

void loadBinaryGrid(MultiGrid3D &block, std::string fName, bool enforceUint)
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
