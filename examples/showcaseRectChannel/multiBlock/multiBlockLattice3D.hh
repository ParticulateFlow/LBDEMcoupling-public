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
 * A 3D multiblock lattice -- generic implementation.
 */
#ifndef MULTI_BLOCK_LATTICE_3D_HH
#define MULTI_BLOCK_LATTICE_3D_HH

#include <algorithm>
#include <cmath>
#include <limits>

#include "atomicBlock/blockLattice3D.h"
#include "coProcessors/coProcessor3D.h"
#include "core/dynamicsIdentifiers.h"
#include "core/latticeStatistics.h"
#include "core/multiBlockIdentifiers3D.h"
#include "core/plbProfiler.h"
#include "core/plbTypenames.h"
#include "dataProcessors/metaStuffWrapper3D.h"
#include "multiBlock/defaultMultiBlockPolicy3D.h"
#include "multiBlock/multiBlockGenerator3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/nonLocalTransfer3D.h"

namespace plb {

////////////////////// Class MultiBlockLattice3D /////////////////////////

template <typename T, template <typename U> class Descriptor>
const int MultiBlockLattice3D<T, Descriptor>::staticId = meta::registerMultiBlock3D(
    MultiBlockLattice3D<T, Descriptor>::basicType(),
    MultiBlockLattice3D<T, Descriptor>::descriptorType(),
    MultiBlockLattice3D<T, Descriptor>::blockName(),
    defaultGenerateMultiBlockLattice3D<T, Descriptor>);

template <typename T, template <typename U> class Descriptor>
MultiBlockLattice3D<T, Descriptor>::MultiBlockLattice3D(
    MultiBlockManagement3D const &multiBlockManagement_, BlockCommunicator3D *blockCommunicator_,
    CombinedStatistics *combinedStatistics_, MultiCellAccess3D<T, Descriptor> *multiCellAccess_,
    Dynamics<T, Descriptor> *backgroundDynamics_) :
    MultiBlock3D(multiBlockManagement_, blockCommunicator_, combinedStatistics_),
    backgroundDynamics(backgroundDynamics_),
    multiCellAccess(multiCellAccess_)
{
    allocateAndInitialize();
    eliminateStatisticsInEnvelope();
    this->evaluateStatistics();  // Reset statistics to default.
}

template <typename T, template <typename U> class Descriptor>
MultiBlockLattice3D<T, Descriptor>::MultiBlockLattice3D(
    plint nx, plint ny, plint nz, Dynamics<T, Descriptor> *backgroundDynamics_) :
    MultiBlock3D(nx, ny, nz, Descriptor<T>::vicinity),
    backgroundDynamics(backgroundDynamics_),
    multiCellAccess(defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>())
{
    allocateAndInitialize();
    eliminateStatisticsInEnvelope();
    this->evaluateStatistics();  // Reset statistics to default.
}

template <typename T, template <typename U> class Descriptor>
MultiBlockLattice3D<T, Descriptor>::~MultiBlockLattice3D()
{
    for (typename BlockMap::iterator it = blockLattices.begin(); it != blockLattices.end(); ++it) {
        delete it->second;
    }
    delete backgroundDynamics;
    delete multiCellAccess;
}

template <typename T, template <typename U> class Descriptor>
MultiBlockLattice3D<T, Descriptor>::MultiBlockLattice3D(
    MultiBlockLattice3D<T, Descriptor> const &rhs) :
    BlockLatticeBase3D<T, Descriptor>(rhs),
    MultiBlock3D(rhs),
    backgroundDynamics(rhs.backgroundDynamics->clone()),
    multiCellAccess(rhs.multiCellAccess->clone())
{
    for (typename BlockMap::const_iterator it = rhs.blockLattices.begin();
         it != rhs.blockLattices.end(); ++it)
    {
        blockLattices[it->first] = new BlockLattice3D<T, Descriptor>(*it->second);
    }
}

template <typename T, template <typename U> class Descriptor>
MultiBlockLattice3D<T, Descriptor>::MultiBlockLattice3D(MultiBlock3D const &rhs)
    // Use MultiBlock's sub-domain constructor to avoid that the data-processors are copied
    :
    MultiBlock3D(rhs, rhs.getBoundingBox(), false),
    backgroundDynamics(new NoDynamics<T, Descriptor>),
    multiCellAccess(defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>())
{
    allocateAndInitialize();
    eliminateStatisticsInEnvelope();
}

template <typename T, template <typename U> class Descriptor>
MultiBlockLattice3D<T, Descriptor>::MultiBlockLattice3D(
    MultiBlock3D const &rhs, Box3D subDomain, bool crop) :
    MultiBlock3D(rhs, subDomain, crop),
    backgroundDynamics(new NoDynamics<T, Descriptor>),
    multiCellAccess(defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>())
{
    allocateAndInitialize();
    eliminateStatisticsInEnvelope();
}

template <typename T, template <typename U> class Descriptor>
void MultiBlockLattice3D<T, Descriptor>::swap(MultiBlockLattice3D<T, Descriptor> &rhs)
{
    BlockLatticeBase3D<T, Descriptor>::swap(rhs);
    MultiBlock3D::swap(rhs);
    std::swap(backgroundDynamics, rhs.backgroundDynamics);
    std::swap(multiCellAccess, rhs.multiCellAccess);
    blockLattices.swap(rhs.blockLattices);
}

template <typename T, template <typename U> class Descriptor>
MultiBlockLattice3D<T, Descriptor> &MultiBlockLattice3D<T, Descriptor>::operator=(
    MultiBlockLattice3D<T, Descriptor> const &rhs)
{
    MultiBlockLattice3D<T, Descriptor> tmp(rhs);
    swap(tmp);
    return *this;
}

template <typename T, template <typename U> class Descriptor>
MultiBlockLattice3D<T, Descriptor> *MultiBlockLattice3D<T, Descriptor>::clone() const
{
    return new MultiBlockLattice3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
MultiBlockLattice3D<T, Descriptor> *MultiBlockLattice3D<T, Descriptor>::clone(
    MultiBlockManagement3D const &newManagement) const
{
    MultiBlockLattice3D<T, Descriptor> *newLattice = new MultiBlockLattice3D<T, Descriptor>(
        newManagement, this->getBlockCommunicator().clone(), this->getCombinedStatistics().clone(),
        multiCellAccess->clone(), getBackgroundDynamics().clone());
    // Make sure background dynamics is not used in the cloned lattice, and instead every cell
    // has an individual dynamics instance. This is the behavior we expect most of the time
    // when cloning a lattice. As a matter of fact, the background dynamics is now mostly
    // considered to be dangereous, and is kept as a default behavior of the multi-block lattice
    // for legacy reasons only.
    newLattice->resetDynamics(getBackgroundDynamics());
    // Use the same domain in the "from" and "to" argument, so that the data is not shifted
    // in space during the creation of the new block.
    copy(
        *this, newLattice->getBoundingBox(), *newLattice, newLattice->getBoundingBox(),
        modif::dataStructure);
    return newLattice;
}

template <typename T, template <typename U> class Descriptor>
void MultiBlockLattice3D<T, Descriptor>::resetDynamics(Dynamics<T, Descriptor> const &dynamics)
{
    for (typename BlockMap::const_iterator it = blockLattices.begin(); it != blockLattices.end();
         ++it) {
        blockLattices[it->first]->resetDynamics(dynamics);
    }
}

template <typename T, template <typename U> class Descriptor>
Dynamics<T, Descriptor> const &MultiBlockLattice3D<T, Descriptor>::getBackgroundDynamics() const
{
    return *backgroundDynamics;
}

template <typename T, template <typename U> class Descriptor>
Cell<T, Descriptor> &MultiBlockLattice3D<T, Descriptor>::get(plint iX, plint iY, plint iZ)
{
    return multiCellAccess->getDistributedCell(
        iX, iY, iZ, this->getMultiBlockManagement(), blockLattices);
}

template <typename T, template <typename U> class Descriptor>
Cell<T, Descriptor> const &MultiBlockLattice3D<T, Descriptor>::get(
    plint iX, plint iY, plint iZ) const
{
    return multiCellAccess->getDistributedCell(
        iX, iY, iZ, this->getMultiBlockManagement(), blockLattices);
}

template <typename T, template <typename U> class Descriptor>
void MultiBlockLattice3D<T, Descriptor>::specifyStatisticsStatus(Box3D domain, bool status)
{
    Box3D inters;
    for (typename BlockMap::iterator it = blockLattices.begin(); it != blockLattices.end(); ++it) {
        SmartBulk3D bulk(this->getMultiBlockManagement(), it->first);
        if (intersect(domain, bulk.getBulk(), inters)) {
            inters = bulk.toLocal(inters);
            it->second->specifyStatisticsStatus(inters, status);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void MultiBlockLattice3D<T, Descriptor>::collide(Box3D domain)
{
    Box3D inters;
    for (typename BlockMap::iterator it = blockLattices.begin(); it != blockLattices.end(); ++it) {
        SmartBulk3D bulk(this->getMultiBlockManagement(), it->first);
        if (intersect(domain, bulk.computeEnvelope(), inters)) {
            inters = bulk.toLocal(inters);
            it->second->collide(inters);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void MultiBlockLattice3D<T, Descriptor>::collide()
{
    for (typename BlockMap::iterator it = blockLattices.begin(); it != blockLattices.end(); ++it) {
        it->second->collide();
    }
}

template <typename T, template <typename U> class Descriptor>
void MultiBlockLattice3D<T, Descriptor>::stream(Box3D domain)
{
    Box3D inters;
    for (typename BlockMap::iterator it = blockLattices.begin(); it != blockLattices.end(); ++it) {
        SmartBulk3D bulk(this->getMultiBlockManagement(), it->first);
        if (intersect(domain, bulk.computeNonPeriodicEnvelope(), inters)) {
            inters = bulk.toLocal(inters);
            it->second->stream(inters);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
Box3D MultiBlockLattice3D<T, Descriptor>::extendPeriodic(
    Box3D const &box, plint envelopeWidth) const
{
    Box3D boundingBox(this->getBoundingBox());
    Box3D periodicBox(box);
    bool periodicX = this->periodicity().get(0);
    bool periodicY = this->periodicity().get(1);
    bool periodicZ = this->periodicity().get(2);
    if (periodicX) {
        if (periodicBox.x0 == boundingBox.x0) {
            periodicBox.x0 -= envelopeWidth;
        }
        if (periodicBox.x1 == boundingBox.x1) {
            periodicBox.x1 += envelopeWidth;
        }
    }
    if (periodicY) {
        if (periodicBox.y0 == boundingBox.y0) {
            periodicBox.y0 -= envelopeWidth;
        }
        if (periodicBox.y1 == boundingBox.y1) {
            periodicBox.y1 += envelopeWidth;
        }
    }
    if (periodicZ) {
        if (periodicBox.z0 == boundingBox.z0) {
            periodicBox.z0 -= envelopeWidth;
        }
        if (periodicBox.z1 == boundingBox.z1) {
            periodicBox.z1 += envelopeWidth;
        }
    }
    return periodicBox;
}

template <typename T, template <typename U> class Descriptor>
void MultiBlockLattice3D<T, Descriptor>::stream()
{
    streamImplementation();
    this->executeInternalProcessors();
    this->evaluateStatistics();
    this->incrementTime();
}

template <typename T, template <typename U> class Descriptor>
void MultiBlockLattice3D<T, Descriptor>::externalStream()
{
    streamImplementation();
}

template <typename T, template <typename U> class Descriptor>
void MultiBlockLattice3D<T, Descriptor>::streamImplementation()
{
    for (typename BlockMap::iterator it = blockLattices.begin(); it != blockLattices.end(); ++it) {
        SmartBulk3D bulk(this->getMultiBlockManagement(), it->first);
        // Stream must be applied to full domain, including currently active envelopes.
        Box3D domain = extendPeriodic(
            bulk.computeNonPeriodicEnvelope(), this->getMultiBlockManagement().getEnvelopeWidth());
        it->second->stream(bulk.toLocal(domain));
    }
}

template <typename T, template <typename U> class Descriptor>
void MultiBlockLattice3D<T, Descriptor>::collideAndStream(Box3D domain)
{
    Box3D inters;
    for (typename BlockMap::iterator it = blockLattices.begin(); it != blockLattices.end(); ++it) {
        SmartBulk3D bulk(this->getMultiBlockManagement(), it->first);
        if (intersect(domain, bulk.computeNonPeriodicEnvelope(), inters)) {
            inters = bulk.toLocal(inters);
            it->second->collideAndStream(inters);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void MultiBlockLattice3D<T, Descriptor>::collideAndStream()
{
    global::profiler().start("cycle");
    collideAndStreamImplementation();
    this->executeInternalProcessors();
    this->evaluateStatistics();
    this->incrementTime();
    global::profiler().stop("cycle");
    if (global::profiler().cyclingIsAutomatic()) {
        global::profiler().cycle();
    }
}

template <typename T, template <typename U> class Descriptor>
void MultiBlockLattice3D<T, Descriptor>::externalCollideAndStream()
{
    global::profiler().start("cycle");
    collideAndStreamImplementation();
    if (global::profiler().cyclingIsAutomatic()) {
        global::profiler().cycle();
    }
    global::profiler().stop("cycle");
}

template <typename T, template <typename U> class Descriptor>
void MultiBlockLattice3D<T, Descriptor>::collideAndStreamImplementation()
{
    ThreadAttribution const &threadAttribution =
        this->getMultiBlockManagement().getThreadAttribution();
    if (threadAttribution.hasCoProcessors()) {
        for (typename BlockMap::iterator it = blockLattices.begin(); it != blockLattices.end();
             ++it) {
            plint blockId = it->first;
            int handle = threadAttribution.getCoProcessorHandle(blockId);
            if (handle >= 0) {
                global::defaultCoProcessor3D<T>().collideAndStream(handle);
            } else {
                SmartBulk3D bulk(this->getMultiBlockManagement(), blockId);
                // CollideAndStream must be applied to full domain,
                //   including currently active envelopes.
                Box3D domain = extendPeriodic(
                    bulk.computeNonPeriodicEnvelope(),
                    this->getMultiBlockManagement().getEnvelopeWidth());
                it->second->collideAndStream(bulk.toLocal(domain));
            }
        }
    } else {
        for (typename BlockMap::iterator it = blockLattices.begin(); it != blockLattices.end();
             ++it) {
            SmartBulk3D bulk(this->getMultiBlockManagement(), it->first);
            // CollideAndStream must be applied to full domain,
            //   including currently active envelopes.
            Box3D domain = extendPeriodic(
                bulk.computeNonPeriodicEnvelope(),
                this->getMultiBlockManagement().getEnvelopeWidth());
            it->second->collideAndStream(bulk.toLocal(domain));
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void MultiBlockLattice3D<T, Descriptor>::incrementTime()
{
    for (typename BlockMap::iterator it = blockLattices.begin(); it != blockLattices.end(); ++it) {
        it->second->incrementTime();
    }
    this->getTimeCounter().incrementTime();
}

template <typename T, template <typename U> class Descriptor>
void MultiBlockLattice3D<T, Descriptor>::resetTime(pluint value)
{
    for (typename BlockMap::iterator it = blockLattices.begin(); it != blockLattices.end(); ++it) {
        it->second->getTimeCounter().resetTime(value);
    }
    this->getTimeCounter().resetTime(value);
}

template <typename T, template <typename U> class Descriptor>
void MultiBlockLattice3D<T, Descriptor>::allocateAndInitialize()
{
    this->getInternalStatistics().subscribeAverage();  // Subscribe average rho-bar
    this->getInternalStatistics().subscribeAverage();  // Subscribe average uSqr
    this->getInternalStatistics().subscribeMax();      // Subscribe max uSqr

    for (pluint iBlock = 0; iBlock < this->getLocalInfo().getBlocks().size(); ++iBlock) {
        plint blockId = this->getLocalInfo().getBlocks()[iBlock];
        SmartBulk3D bulk(this->getMultiBlockManagement(), blockId);
        Box3D envelope = bulk.computeEnvelope();
        BlockLattice3D<T, Descriptor> *newLattice = new BlockLattice3D<T, Descriptor>(
            envelope.getNx(), envelope.getNy(), envelope.getNz(), backgroundDynamics->clone());
        newLattice->setLocation(Dot3D(envelope.x0, envelope.y0, envelope.z0));
        blockLattices[blockId] = newLattice;
    }
}

template <typename T, template <typename U> class Descriptor>
void MultiBlockLattice3D<T, Descriptor>::eliminateStatisticsInEnvelope()
{
    for (typename BlockMap::iterator it = blockLattices.begin(); it != blockLattices.end(); ++it) {
        plint envelopeWidth = this->getMultiBlockManagement().getEnvelopeWidth();
        BlockLattice3D<T, Descriptor> &block = *it->second;
        plint maxX = block.getNx() - 1;
        plint maxY = block.getNy() - 1;
        plint maxZ = block.getNz() - 1;

        block.specifyStatisticsStatus(Box3D(0, maxX, 0, maxY, 0, envelopeWidth - 1), false);
        block.specifyStatisticsStatus(
            Box3D(0, maxX, 0, maxY, maxZ - envelopeWidth + 1, maxZ), false);
        block.specifyStatisticsStatus(Box3D(0, maxX, 0, envelopeWidth - 1, 0, maxZ), false);
        block.specifyStatisticsStatus(
            Box3D(0, maxX, maxY - envelopeWidth + 1, maxY, 0, maxZ), false);
        block.specifyStatisticsStatus(Box3D(0, envelopeWidth - 1, 0, maxY, 0, maxZ), false);
        block.specifyStatisticsStatus(
            Box3D(maxX - envelopeWidth + 1, maxX, 0, maxY, 0, maxZ), false);
    }
}

template <typename T, template <typename U> class Descriptor>
std::map<plint, BlockLattice3D<T, Descriptor> *>
    &MultiBlockLattice3D<T, Descriptor>::getBlockLattices()
{
    return blockLattices;
}

template <typename T, template <typename U> class Descriptor>
std::map<plint, BlockLattice3D<T, Descriptor> *> const &
    MultiBlockLattice3D<T, Descriptor>::getBlockLattices() const
{
    return blockLattices;
}

template <typename T, template <typename U> class Descriptor>
void MultiBlockLattice3D<T, Descriptor>::getDynamicsDict(
    Box3D domain, std::map<std::string, int> &dict)
{
    std::vector<int> ids;
    uniqueDynamicsIds(*this, domain, ids);
    dict.clear();
    for (pluint i = 0; i < ids.size(); ++i) {
        int id = ids[i];
        std::string name = meta::dynamicsRegistration<T, Descriptor>().getName(id);
        dict.insert(std::pair<std::string, int>(name, id));
    }
}

template <typename T, template <typename U> class Descriptor>
std::string MultiBlockLattice3D<T, Descriptor>::getBlockName() const
{
    return std::string("BlockLattice3D");
}

template <typename T, template <typename U> class Descriptor>
std::vector<std::string> MultiBlockLattice3D<T, Descriptor>::getTypeInfo() const
{
    std::vector<std::string> info;
    info.push_back(basicType());
    info.push_back(descriptorType());
    return info;
}

template <typename T, template <typename U> class Descriptor>
std::string MultiBlockLattice3D<T, Descriptor>::blockName()
{
    return std::string("BlockLattice3D");
}

template <typename T, template <typename U> class Descriptor>
std::string MultiBlockLattice3D<T, Descriptor>::basicType()
{
    return std::string(NativeType<T>::getName());
}

template <typename T, template <typename U> class Descriptor>
std::string MultiBlockLattice3D<T, Descriptor>::descriptorType()
{
    return std::string(Descriptor<T>::name);
}

template <typename T, template <typename U> class Descriptor>
BlockLattice3D<T, Descriptor> &MultiBlockLattice3D<T, Descriptor>::getComponent(plint blockId)
{
    typename BlockMap::iterator it = blockLattices.find(blockId);
    PLB_ASSERT(it != blockLattices.end());
    return *it->second;
}

template <typename T, template <typename U> class Descriptor>
BlockLattice3D<T, Descriptor> const &MultiBlockLattice3D<T, Descriptor>::getComponent(
    plint blockId) const
{
    typename BlockMap::const_iterator it = blockLattices.find(blockId);
    PLB_ASSERT(it != blockLattices.end());
    return *it->second;
}

template <typename T, template <typename U> class Descriptor>
plint MultiBlockLattice3D<T, Descriptor>::sizeOfCell() const
{
    return sizeof(T) * (Descriptor<T>::numPop + Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor>
plint MultiBlockLattice3D<T, Descriptor>::getCellDim() const
{
    return Descriptor<T>::numPop + Descriptor<T>::ExternalField::numScalars;
}

template <typename T, template <typename U> class Descriptor>
int MultiBlockLattice3D<T, Descriptor>::getStaticId() const
{
    return staticId;
}

template <typename T, template <typename U> class Descriptor>
void MultiBlockLattice3D<T, Descriptor>::copyReceive(
    MultiBlock3D const &fromBlock, Box3D const &fromDomain, Box3D const &toDomain,
    modif::ModifT whichData)
{
    MultiBlockLattice3D<T, Descriptor> const *fromLattice =
        dynamic_cast<MultiBlockLattice3D<T, Descriptor> const *>(&fromBlock);
    PLB_ASSERT(fromLattice);
    copy(*fromLattice, fromDomain, *this, toDomain, whichData);
}

/////////// Free Functions //////////////////////////////

template <typename T, template <typename U> class Descriptor>
MultiBlockLattice3D<T, Descriptor> &findMultiBlockLattice3D(id_t id)
{
    MultiBlock3D *multiBlock = multiBlockRegistration3D().find(id);
    if (!multiBlock || multiBlock->getStaticId() != MultiBlockLattice3D<T, Descriptor>::staticId) {
        throw PlbLogicException("Trying to access a multi block lattice that is not registered.");
    }
    return (MultiBlockLattice3D<T, Descriptor> &)(*multiBlock);
}

template <typename T, template <typename U> class Descriptor>
double getStoredAverageDensity(MultiBlockLattice3D<T, Descriptor> const &blockLattice)
{
    return Descriptor<T>::fullRho(
        blockLattice.getInternalStatistics().getAverage(LatticeStatistics::avRhoBar));
}

template <typename T, template <typename U> class Descriptor>
double getStoredAverageEnergy(MultiBlockLattice3D<T, Descriptor> const &blockLattice)
{
    return 0.5 * blockLattice.getInternalStatistics().getAverage(LatticeStatistics::avUSqr);
}

template <typename T, template <typename U> class Descriptor>
double getStoredMaxVelocity(MultiBlockLattice3D<T, Descriptor> const &blockLattice)
{
    return std::sqrt(blockLattice.getInternalStatistics().getMax(LatticeStatistics::maxUSqr));
}

}  // namespace plb

#endif  // MULTI_BLOCK_LATTICE_3D_HH
