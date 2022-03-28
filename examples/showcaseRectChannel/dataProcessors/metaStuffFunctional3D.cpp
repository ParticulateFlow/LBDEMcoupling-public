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

#include "dataProcessors/metaStuffFunctional3D.h"

#include <cmath>

#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/dataProcessingFunctional3D.hh"
#include "core/blockStatistics.h"
#include "core/plbDebug.h"
#include "core/util.h"
#include "dataProcessors/metaStuffHelper.h"
#include "parallelism/mpiManager.h"

namespace plb {

/* ******** IterateDynamicsFunctional3D ************************************ */

IterateDynamicsFunctional3D::IterateDynamicsFunctional3D(std::vector<int> previousMaximum_) :
    previousMaximum(previousMaximum_), maxIds(previousMaximum.size())
{
    for (pluint i = 0; i < maxIds.size(); ++i) {
        maxIds[i] = this->getStatistics().subscribeMax();
    }
}

void IterateDynamicsFunctional3D::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    AtomicContainerBlock3D &container = *dynamic_cast<AtomicContainerBlock3D *>(blocks[0]);
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

IterateDynamicsFunctional3D *IterateDynamicsFunctional3D::clone() const
{
    return new IterateDynamicsFunctional3D(*this);
}

void IterateDynamicsFunctional3D::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

BlockDomain::DomainT IterateDynamicsFunctional3D::appliesTo() const
{
    return BlockDomain::bulk;
}

std::vector<int> IterateDynamicsFunctional3D::getNextMaximum() const
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

/* ******** AllFlagsTrueFunctional3D ************************************ */

AllFlagsTrueFunctional3D::AllFlagsTrueFunctional3D() :
    numFalseId(this->getStatistics().subscribeIntSum())
{ }

void AllFlagsTrueFunctional3D::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    if (!blocks[0]->getFlag()) {
        this->getStatistics().gatherIntSum(numFalseId, 1);
    }
}

AllFlagsTrueFunctional3D *AllFlagsTrueFunctional3D::clone() const
{
    return new AllFlagsTrueFunctional3D(*this);
}

void AllFlagsTrueFunctional3D::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
}

BlockDomain::DomainT AllFlagsTrueFunctional3D::appliesTo() const
{
    return BlockDomain::bulk;
}

bool AllFlagsTrueFunctional3D::allTrue() const
{
    return this->getStatistics().getIntSum(numFalseId) == 0;
}

/* ******** GetThreadNumFunctional3D ************************************ */

void GetThreadNumFunctional3D::process(Box3D domain, ScalarField3D<int> &threadNum)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                threadNum.get(iX, iY, iZ) = global::mpi().getRank() + 1;
            }
        }
    }
}

GetThreadNumFunctional3D *GetThreadNumFunctional3D::clone() const
{
    return new GetThreadNumFunctional3D(*this);
}

void GetThreadNumFunctional3D::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

}  // namespace plb
