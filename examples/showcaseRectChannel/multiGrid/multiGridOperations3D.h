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
 * Operations on the 3D multiGrid -- header file.
 */
#ifndef MULTI_GRID_OPERATIONS_3D_H
#define MULTI_GRID_OPERATIONS_3D_H

#include <vector>

#include "core/blockStatistics.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "multiGrid/multiGrid3D.h"

namespace plb {

class MultiBlock3D;
struct DataProcessorGenerator3D;
class ReductiveDataProcessorGenerator3D;

void executeDataProcessor(
    DataProcessorGenerator3D const &generator, std::vector<MultiGrid3D *> multiGrids,
    plint referenceLevel);

void executeDataProcessor(
    DataProcessorGenerator3D const &generator, MultiGrid3D &object, plint referenceLevel);

void executeDataProcessor(
    DataProcessorGenerator3D const &generator, MultiGrid3D &object1, MultiGrid3D &object2,
    plint referenceLevel);

void executeDataProcessor(
    ReductiveDataProcessorGenerator3D &generator, std::vector<MultiGrid3D *> multiGrids,
    plint referenceLevel);

void executeDataProcessor(
    ReductiveDataProcessorGenerator3D &generator, MultiGrid3D &object, plint referenceLevel);

void executeDataProcessor(
    ReductiveDataProcessorGenerator3D &generator, MultiGrid3D &object1, MultiGrid3D &object2,
    plint referenceLevel);

void addInternalProcessor(
    DataProcessorGenerator3D const &generator, std::vector<MultiGrid3D *> multiGrids,
    plint referenceLevel, plint processorLevel = 0);

void addInternalProcessor(
    DataProcessorGenerator3D const &generator, MultiGrid3D &object, plint referenceLevel,
    plint processorLevel = 0);

void addInternalProcessor(
    DataProcessorGenerator3D const &generator, MultiGrid3D &object1, MultiGrid3D &object2,
    plint referenceLevel, plint processorLevel = 0);

}  // namespace plb

#endif  // MULTI_GRID_OPERATIONS_3D_H
