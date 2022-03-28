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
 * Atomic block -- implementation file.
 */
#include "atomicBlock/atomicBlock2D.h"

#include "atomicBlock/atomicBlockSerializer2D.h"

namespace plb {

/* *************** Class StatSubscriber2D *********************************** */

StatSubscriber2D::StatSubscriber2D(AtomicBlock2D &block_) : block(block_) { }

plint StatSubscriber2D::subscribeAverage()
{
    return block.getInternalStatistics().subscribeAverage();
}

plint StatSubscriber2D::subscribeSum()
{
    return block.getInternalStatistics().subscribeSum();
}

plint StatSubscriber2D::subscribeMax()
{
    return block.getInternalStatistics().subscribeMax();
}

plint StatSubscriber2D::subscribeIntSum()
{
    return block.getInternalStatistics().subscribeIntSum();
}

/* *************** Class AtomicBlock2D ************************************** */

AtomicBlock2D::AtomicBlock2D(plint nx_, plint ny_) :
    nx(nx_), ny(ny_), location(0, 0), flag(false), internalStatistics(), statisticsSubscriber(*this)
{ }

AtomicBlock2D::AtomicBlock2D(AtomicBlock2D const &rhs) :
    nx(rhs.nx),
    ny(rhs.ny),
    location(rhs.location),
    flag(rhs.flag),
    internalStatistics(rhs.internalStatistics),
    statisticsSubscriber(*this)
{
    copyDataProcessors(rhs.explicitInternalProcessors, explicitInternalProcessors);
    copyDataProcessors(rhs.automaticInternalProcessors, automaticInternalProcessors);
}

AtomicBlock2D::~AtomicBlock2D()
{
    clearDataProcessors();
}

void AtomicBlock2D::swap(AtomicBlock2D &rhs)
{
    std::swap(nx, rhs.nx);
    std::swap(ny, rhs.ny);
    std::swap(location, rhs.location);
    std::swap(flag, rhs.flag);
    std::swap(internalStatistics, rhs.internalStatistics);
    explicitInternalProcessors.swap(rhs.explicitInternalProcessors);
    automaticInternalProcessors.swap(rhs.automaticInternalProcessors);
}

void AtomicBlock2D::initialize()
{
    executeInternalProcessors();
}

BlockStatistics &AtomicBlock2D::getInternalStatistics()
{
    return internalStatistics;
}

BlockStatistics const &AtomicBlock2D::getInternalStatistics() const
{
    return internalStatistics;
}

void AtomicBlock2D::setLocation(Dot2D const &location_)
{
    location = location_;
}

Dot2D AtomicBlock2D::getLocation() const
{
    return location;
}

void AtomicBlock2D::setFlag(bool value)
{
    flag = value;
}

bool AtomicBlock2D::getFlag() const
{
    return flag;
}

Box2D AtomicBlock2D::getBoundingBox() const
{
    return Box2D(0, getNx() - 1, 0, getNy() - 1);
}

void AtomicBlock2D::integrateDataProcessor(DataProcessor2D *processor, plint level)
{
    // Negative level numbers account for explicit internal BlockProcessors
    if (level < 0) {
        integrateDataProcessor(processor, -level - 1, explicitInternalProcessors);
    }
    // Positive-or-zero level numbers account for automatic internal BlockProcessors
    else
    {
        integrateDataProcessor(processor, level, automaticInternalProcessors);
    }
}

void AtomicBlock2D::integrateDataProcessor(
    DataProcessor2D *processor, plint level, DataProcessorVector &processors)
{
    if (level >= (plint)processors.size()) {
        processors.resize(level + 1);
    }
    processors[level].push_back(processor);
}

void AtomicBlock2D::copyDataProcessors(DataProcessorVector const &from, DataProcessorVector &to)
{
    clearDataProcessors(to);
    to.resize(from.size());
    for (pluint iLevel = 0; iLevel < from.size(); ++iLevel) {
        to[iLevel].resize(from[iLevel].size());
        for (pluint iProc = 0; iProc < from[iLevel].size(); ++iProc) {
            to[iLevel][iProc] = from[iLevel][iProc]->clone();
        }
    }
}

void AtomicBlock2D::clearDataProcessors()
{
    clearDataProcessors(explicitInternalProcessors);
    clearDataProcessors(automaticInternalProcessors);
}

void AtomicBlock2D::removeDataProcessors(int staticId)
{
    for (pluint iLevel = 0; iLevel < explicitInternalProcessors.size(); ++iLevel) {
        std::vector<DataProcessor2D *>::iterator it = explicitInternalProcessors[iLevel].begin();
        for (; it != explicitInternalProcessors[iLevel].end(); ++it) {
            if ((*it)->getStaticId() == staticId) {
                delete *it;
                it = explicitInternalProcessors[iLevel].erase(it);
            }
        }
    }
    for (pluint iLevel = 0; iLevel < automaticInternalProcessors.size(); ++iLevel) {
        std::vector<DataProcessor2D *>::iterator it = automaticInternalProcessors[iLevel].begin();
        for (; it != automaticInternalProcessors[iLevel].end(); ++it) {
            if ((*it)->getStaticId() == staticId) {
                delete *it;
                it = automaticInternalProcessors[iLevel].erase(it);
            }
        }
    }
}

void AtomicBlock2D::clearDataProcessors(DataProcessorVector &processors)
{
    for (pluint iLevel = 0; iLevel < processors.size(); ++iLevel) {
        for (pluint iProc = 0; iProc < processors[iLevel].size(); ++iProc) {
            delete processors[iLevel][iProc];
        }
    }
    processors.clear();
}

void AtomicBlock2D::executeInternalProcessors()
{
    for (pluint iLevel = 0; iLevel < automaticInternalProcessors.size(); ++iLevel) {
        executeInternalProcessors(iLevel, automaticInternalProcessors);
    }
}

void AtomicBlock2D::executeInternalProcessors(plint level)
{
    // Negative level numbers account for explicit internal BlockProcessors
    if (level < 0) {
        executeInternalProcessors(-level - 1, explicitInternalProcessors);
    }
    // Positive-or-zero level numbers account for automatic internal BlockProcessors
    else
    {
        executeInternalProcessors(level, automaticInternalProcessors);
    }
}

void AtomicBlock2D::executeInternalProcessors(plint level, DataProcessorVector &processors)
{
    if (level < (plint)processors.size()) {
        for (pluint iProc = 0; iProc < processors[level].size(); ++iProc) {
            processors[level][iProc]->process();
        }
    }
}

DataSerializer *AtomicBlock2D::getBlockSerializer(
    Box2D const &domain, IndexOrdering::OrderingT ordering) const
{
    return new AtomicBlockSerializer2D(*this, domain, ordering);
}

DataUnSerializer *AtomicBlock2D::getBlockUnSerializer(
    Box2D const &domain, IndexOrdering::OrderingT ordering)
{
    return new AtomicBlockUnSerializer2D(*this, domain, ordering);
}

StatSubscriber &AtomicBlock2D::internalStatSubscription()
{
    return statisticsSubscriber;
}

void AtomicBlock2D::evaluateStatistics()
{
    getInternalStatistics().evaluate();
}

Dot2D computeRelativeDisplacement(AtomicBlock2D const &block1, AtomicBlock2D const &block2)
{
    return Dot2D(
        block1.getLocation().x - block2.getLocation().x,
        block1.getLocation().y - block2.getLocation().y);
}

}  // namespace plb
