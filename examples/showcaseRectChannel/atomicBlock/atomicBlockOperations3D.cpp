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
 * Operations on the 3D multiblock -- implementation.
 */

#include "atomicBlock/atomicBlockOperations3D.h"

#include "core/plbDebug.h"

namespace plb {

void executeDataProcessor(
    DataProcessorGenerator3D const &generator, std::vector<AtomicBlock3D *> objects)
{
    DataProcessor3D *processor = generator.generate(objects);
    processor->process();
    delete processor;
}

void executeDataProcessor(DataProcessorGenerator3D const &generator, AtomicBlock3D &object)
{
    std::vector<AtomicBlock3D *> objects(1);
    objects[0] = &object;
    executeDataProcessor(generator, objects);
}

void executeDataProcessor(
    DataProcessorGenerator3D const &generator, AtomicBlock3D &object1, AtomicBlock3D &object2)
{
    std::vector<AtomicBlock3D *> objects(2);
    objects[0] = &object1;
    objects[1] = &object2;
    executeDataProcessor(generator, objects);
}

void executeDataProcessor(
    ReductiveDataProcessorGenerator3D &generator, std::vector<AtomicBlock3D *> objects)
{
    DataProcessor3D *processor = generator.generate(objects);
    processor->process();
    delete processor;
}

void executeDataProcessor(ReductiveDataProcessorGenerator3D &generator, AtomicBlock3D &object)
{
    std::vector<AtomicBlock3D *> objects(1);
    objects[0] = &object;
    executeDataProcessor(generator, objects);
}

void executeDataProcessor(
    ReductiveDataProcessorGenerator3D &generator, AtomicBlock3D &object1, AtomicBlock3D &object2)
{
    std::vector<AtomicBlock3D *> objects(2);
    objects[0] = &object1;
    objects[1] = &object2;
    executeDataProcessor(generator, objects);
}

void addInternalProcessor(
    DataProcessorGenerator3D const &generator, AtomicBlock3D &actor,
    std::vector<AtomicBlock3D *> objects, plint level)
{
    PLB_PRECONDITION(!objects.empty());
    actor.integrateDataProcessor(generator.generate(objects), level);
}

void addInternalProcessor(
    DataProcessorGenerator3D const &generator, std::vector<AtomicBlock3D *> objects, plint level)
{
    PLB_PRECONDITION(!objects.empty());
    addInternalProcessor(generator, *objects[0], objects, level);
}

void addInternalProcessor(
    DataProcessorGenerator3D const &generator, AtomicBlock3D &object, plint level)
{
    std::vector<AtomicBlock3D *> objects(1);
    objects[0] = &object;
    addInternalProcessor(generator, objects, level);
}

void addInternalProcessor(
    DataProcessorGenerator3D const &generator, AtomicBlock3D &object1, AtomicBlock3D &object2,
    plint level)
{
    std::vector<AtomicBlock3D *> objects(2);
    objects[0] = &object1;
    objects[1] = &object2;
    addInternalProcessor(generator, objects, level);
}

}  // namespace plb
