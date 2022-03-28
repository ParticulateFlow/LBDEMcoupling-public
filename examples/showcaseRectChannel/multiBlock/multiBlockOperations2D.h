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
 * Operations on the 2D multiblock -- header file.
 */
#ifndef MULTI_BLOCK_OPERATIONS_2D_H
#define MULTI_BLOCK_OPERATIONS_2D_H

#include <vector>

#include "core/blockStatistics.h"
#include "core/geometry2D.h"
#include "core/globalDefs.h"

namespace plb {

class MultiBlock2D;
struct DataProcessorGenerator2D;
class ReductiveDataProcessorGenerator2D;

void executeDataProcessor(
    DataProcessorGenerator2D const &generator, std::vector<MultiBlock2D *> multiBlocks);

void executeDataProcessor(DataProcessorGenerator2D const &generator, MultiBlock2D &object);

void executeDataProcessor(
    DataProcessorGenerator2D const &generator, MultiBlock2D &object1, MultiBlock2D &object2);

void executeDataProcessor(
    ReductiveDataProcessorGenerator2D &generator, std::vector<MultiBlock2D *> multiBlocks);

void executeDataProcessor(ReductiveDataProcessorGenerator2D &generator, MultiBlock2D &object);

void executeDataProcessor(
    ReductiveDataProcessorGenerator2D &generator, MultiBlock2D &object1, MultiBlock2D &object2);

void addInternalProcessor(
    DataProcessorGenerator2D const &generator, std::vector<MultiBlock2D *> multiBlocks,
    plint level = 0);

void addInternalProcessor(
    DataProcessorGenerator2D const &generator, MultiBlock2D &object, plint level = 0);

void addInternalProcessor(
    DataProcessorGenerator2D const &generator, MultiBlock2D &object1, MultiBlock2D &object2,
    plint level = 0);

template <class OriginalGenerator, class MutableGenerator>
class MultiProcessing2D {
public:
    MultiProcessing2D(OriginalGenerator &generator_, std::vector<MultiBlock2D *> multiBlocks_);
    ~MultiProcessing2D();
    void extractProcessorsOnFirstBlock(BlockDomain::DomainT appliesTo);
    void intersectWithRemainingBlocks(BlockDomain::DomainT appliesTo);
    void subdivideGenerator();
    void adjustCoordinates();
    std::vector<MutableGenerator *> const &getRetainedGenerators() const;
    std::vector<std::vector<plint> > const &getAtomicBlockNumbers() const;
    void multiBlocksWhichRequireUpdate(
        std::vector<MultiBlock2D *> &multiBlocksModifiedByProcessor,
        std::vector<modif::ModifT> &typesOfModification) const;
    void updateEnvelopesWhereRequired();
    std::vector<MutableGenerator *> releaseRetainedGenerators();

private:
    void extractGeneratorOnBlocks(
        std::vector<Box2D> const &finalDomains, std::vector<std::vector<plint> > const &finalIds,
        plint shiftX = 0, plint shiftY = 0);

private:
    OriginalGenerator &generator;
    std::vector<MultiBlock2D *> multiBlocks;
    MultiBlock2D *firstMultiBlock;
    std::vector<MutableGenerator *> retainedGenerators;
    std::vector<std::vector<plint> > atomicBlockNumbers;
};

}  // namespace plb

#endif  // MULTI_BLOCK_OPERATIONS_2D_H
