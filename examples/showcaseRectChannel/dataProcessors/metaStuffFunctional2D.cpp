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

#include "dataProcessors/metaStuffFunctional2D.h"

#include <cmath>

#include "atomicBlock/atomicBlock2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include "core/blockStatistics.h"
#include "core/dynamicsIdentifiers.h"
#include "core/plbDebug.h"
#include "core/util.h"

namespace plb {

/* ******** IterateDynamicsFunctional2D ************************************ */

IterateDynamicsFunctional2D::IterateDynamicsFunctional2D(std::vector<int> previousMaximum_) :
    previousMaximum(previousMaximum_), maxIds(previousMaximum.size())
{
    for (pluint i = 0; i < maxIds.size(); ++i) {
        maxIds[i] = this->getStatistics().subscribeMax();
    }
}

void IterateDynamicsFunctional2D::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    AtomicContainerBlock2D &container = *dynamic_cast<AtomicContainerBlock2D *>(blocks[0]);
    StoreDynamicsID *storeId = dynamic_cast<StoreDynamicsID *>(container.getData());
    std::vector<int> nextIDs(previousMaximum.size());
    for (pluint i = 0; i < nextIDs.size(); ++i) {
        nextIDs[i] = -1;
    }
    if (!storeId->empty()) {
        nextIDs = storeId->getCurrent();
        util::extendVectorSize(nextIDs, previousMaximum.size());
        if (vectorEquals(nextIDs, previousMaximum)) {
            nextIDs = storeId->iterate();
            util::extendVectorSize(nextIDs, previousMaximum.size());
        }
    }
    for (pluint i = 0; i < nextIDs.size(); ++i) {
        if (i > 0 && this->getStatistics().getMax(maxIds[i - 1]) >= (double)nextIDs[i - 1]) {
            this->getStatistics().gatherMax(maxIds[i], (double)nextIDs[i]);
        }
    }
}

IterateDynamicsFunctional2D *IterateDynamicsFunctional2D::clone() const
{
    return new IterateDynamicsFunctional2D(*this);
}

void IterateDynamicsFunctional2D::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

BlockDomain::DomainT IterateDynamicsFunctional2D::appliesTo() const
{
    return BlockDomain::bulk;
}

std::vector<int> IterateDynamicsFunctional2D::getNextMaximum() const
{
    std::vector<int> maximum(maxIds.size());
    for (pluint i = 0; i < maxIds.size(); ++i) {
        double d_maximum = this->getStatistics().getMax(maxIds[i]);
        if (d_maximum < 0.5) {
            maximum[i] = -1;
        } else {
            maximum[i] = (int)(.5 + d_maximum);
        }
    }
    return maximum;
}

/* ******** AllFlagsTrueFunctional2D ************************************ */

AllFlagsTrueFunctional2D::AllFlagsTrueFunctional2D() :
    numFalseId(this->getStatistics().subscribeIntSum())
{ }

void AllFlagsTrueFunctional2D::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    if (!blocks[0]->getFlag()) {
        this->getStatistics().gatherIntSum(numFalseId, 1);
    }
}

AllFlagsTrueFunctional2D *AllFlagsTrueFunctional2D::clone() const
{
    return new AllFlagsTrueFunctional2D(*this);
}

void AllFlagsTrueFunctional2D::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
}

BlockDomain::DomainT AllFlagsTrueFunctional2D::appliesTo() const
{
    return BlockDomain::bulk;
}

bool AllFlagsTrueFunctional2D::allTrue() const
{
    return this->getStatistics().getIntSum(numFalseId) == 0;
}

}  // namespace plb
