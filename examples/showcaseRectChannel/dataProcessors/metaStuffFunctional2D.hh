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
 * Helper functions for domain initialization -- header file.
 */
#ifndef META_STUFF_FUNCTIONAL_2D_HH
#define META_STUFF_FUNCTIONAL_2D_HH

#include <cmath>

#include "atomicBlock/atomicBlock2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include "core/blockStatistics.h"
#include "core/dynamicsIdentifiers.h"
#include "core/plbDebug.h"
#include "core/util.h"
#include "dataProcessors/metaStuffFunctional2D.h"

namespace plb {

/* ******** StoreDynamicsFunctional2D ************************************ */

template <typename T, template <typename U> class Descriptor>
StoreDynamicsFunctional2D<T, Descriptor>::StoreDynamicsFunctional2D() :
    maxChainLengthId(this->getStatistics().subscribeMax())
{ }

template <typename T, template <typename U> class Descriptor>
void StoreDynamicsFunctional2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    BlockLattice2D<T, Descriptor> &lattice =
        *dynamic_cast<BlockLattice2D<T, Descriptor> *>(blocks[0]);
    AtomicContainerBlock2D &container = *dynamic_cast<AtomicContainerBlock2D *>(blocks[1]);
    StoreDynamicsID *storeID = new StoreDynamicsID;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            std::vector<int> chain;
            constructIdChain(lattice.get(iX, iY).getDynamics(), chain);
            storeID->addIdChain(chain);
            this->getStatistics().gatherMax(maxChainLengthId, (double)chain.size());
        }
    }
    storeID->startIterations();
    container.setData(storeID);
}

template <typename T, template <typename U> class Descriptor>
StoreDynamicsFunctional2D<T, Descriptor> *StoreDynamicsFunctional2D<T, Descriptor>::clone() const
{
    return new StoreDynamicsFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void StoreDynamicsFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT StoreDynamicsFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
pluint StoreDynamicsFunctional2D<T, Descriptor>::getMaxChainLength() const
{
    double maximum = this->getStatistics().getMax(maxChainLengthId);
    return (pluint)(.5 + maximum);
}

/* ******** ExtractDynamicsChainFunctional2D ************************************ */

template <typename T, template <typename U> class Descriptor>
ExtractDynamicsChainFunctional2D<T, Descriptor>::ExtractDynamicsChainFunctional2D(
    typename ExtractDynamicsChainFunctional2D<T, Descriptor>::DMap const &dynamicsMap_,
    pluint maxChainSize_) :
    dynamicsMap(dynamicsMap_), maxChainSize(maxChainSize_)
{ }

template <typename T, template <typename U> class Descriptor>
void ExtractDynamicsChainFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<int> &mask)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            std::vector<int> chain;
            constructIdChain(lattice.get(iX, iY).getDynamics(), chain);
            util::extendVectorSize(chain, maxChainSize);
            mask.get(iX, iY) = dynamicsMap[chain];
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ExtractDynamicsChainFunctional2D<T, Descriptor>
    *ExtractDynamicsChainFunctional2D<T, Descriptor>::clone() const
{
    return new ExtractDynamicsChainFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void ExtractDynamicsChainFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT ExtractDynamicsChainFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** ExtractTopMostDynamicsFunctional2D ************************************ */

template <typename T, template <typename U> class Descriptor>
void ExtractTopMostDynamicsFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<int> &mask)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            mask.get(iX, iY) = lattice.get(iX, iY).getDynamics().getId();
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ExtractTopMostDynamicsFunctional2D<T, Descriptor>
    *ExtractTopMostDynamicsFunctional2D<T, Descriptor>::clone() const
{
    return new ExtractTopMostDynamicsFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void ExtractTopMostDynamicsFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT ExtractTopMostDynamicsFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** ExtractBottomMostDynamicsFunctional2D ************************************ */

template <typename T, template <typename U> class Descriptor>
void ExtractBottomMostDynamicsFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<int> &mask)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            mask.get(iX, iY) = getBottomMostDynamics(lattice.get(iX, iY).getDynamics()).getId();
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ExtractBottomMostDynamicsFunctional2D<T, Descriptor>
    *ExtractBottomMostDynamicsFunctional2D<T, Descriptor>::clone() const
{
    return new ExtractBottomMostDynamicsFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void ExtractBottomMostDynamicsFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT ExtractBottomMostDynamicsFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** AssignEntireCellFunctional2D ************************************ */

template <typename T, template <typename U> class Descriptor>
void AssignEntireCellFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &sourceLattice,
    BlockLattice2D<T, Descriptor> &destinationLattice)
{
    std::vector<char> data;
    Dot2D offset = computeRelativeDisplacement(sourceLattice, destinationLattice);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Cell<T, Descriptor> const &cell = sourceLattice.get(iX, iY);

            pluint pos = data.size();
            data.resize(pos + CellInfo<T, Descriptor>::n * sizeof(T));
            cell.serialize(&data[pos]);

            HierarchicSerializer serializer(data, cell.getDynamics().getId());
            cell.getDynamics().serialize(serializer);
        }
    }

    pluint pos = 0;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint iX_ = iX + offset.x;
            plint iY_ = iY + offset.y;
            Cell<T, Descriptor> &cell = destinationLattice.get(iX_, iY_);

            cell.unSerialize(&data[pos]);
            pos += CellInfo<T, Descriptor>::n * sizeof(T);

            HierarchicUnserializer unserializer(data, pos);
            destinationLattice.attributeDynamics(
                iX_, iY_, meta::dynamicsRegistration<T, Descriptor>().generate(unserializer));
            pos = unserializer.getCurrentPos();
        }
    }
}

template <typename T, template <typename U> class Descriptor>
AssignEntireCellFunctional2D<T, Descriptor> *AssignEntireCellFunctional2D<T, Descriptor>::clone()
    const
{
    return new AssignEntireCellFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void AssignEntireCellFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT AssignEntireCellFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

}  // namespace plb

#endif  // META_STUFF_FUNCTIONAL_2D_HH
