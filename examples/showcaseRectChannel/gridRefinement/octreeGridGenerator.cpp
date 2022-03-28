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

#include "gridRefinement/octreeGridGenerator.h"

#include <cmath>
#include <cstdio>
#include <limits>

#include "gridRefinement/octree.h"
#include "gridRefinement/octree.hh"
#include "gridRefinement/octreeGridGenerator.hh"
#include "io/parallelIO.h"

namespace plb {
OctreeProcessLoads::OctreeProcessLoads(
    plint numProcesses_, plint numCellsPerBlock_, int minLeafLevel_) :
    numProcesses(numProcesses_), numCellsPerBlock(numCellsPerBlock_), minLeafLevel(minLeafLevel_)
{
    PLB_ASSERT(numProcesses > 0);
    PLB_ASSERT(numCellsPerBlock > 0);
    PLB_ASSERT(minLeafLevel >= 0);

    loads.resize(numProcesses, 0);
}

void OctreeProcessLoads::clearLoads()
{
    std::fill(loads.begin(), loads.end(), 0);
}

void OctreeProcessLoads::clearLoad(plint processId)
{
    PLB_ASSERT(processId >= 0 && processId < numProcesses);
    loads[processId] = 0;
}

plint OctreeProcessLoads::getNumProcesses() const
{
    return (numProcesses);
}

plint OctreeProcessLoads::getNumCellsPerBlock() const
{
    return (numCellsPerBlock);
}

int OctreeProcessLoads::getMinLeafLevel() const
{
    return (minLeafLevel);
}

plint OctreeProcessLoads::getLoad(plint processId) const
{
    PLB_ASSERT(processId >= 0 && processId < numProcesses);
    return (loads[processId]);
}

void OctreeProcessLoads::addLoad(plint processId, int level)
{
    PLB_ASSERT(processId >= 0 && processId < numProcesses);
    PLB_ASSERT(level >= minLeafLevel);
    loads[processId] += numCellsPerBlock * util::intTwoToThePower(level - minLeafLevel);
}

plint OctreeProcessLoads::getProcessWithMinLoad() const
{
    /*
    plint minLoad = loads[0];
    plint processId = 0;
    for (plint iProcess = 1; iProcess < numProcesses; iProcess++) {
        plint load = loads[iProcess];
        if (load < minLoad) {
            minLoad = load;
            processId = iProcess;
        }
    }
    return(processId);
    */

    // return((plint) std::distance(loads.begin(), std::min_element(loads.begin(), loads.end())));
    return ((plint)(std::min_element(loads.begin(), loads.end()) - loads.begin()));
}

}  // namespace plb
