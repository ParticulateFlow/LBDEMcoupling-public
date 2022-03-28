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
#ifndef META_STUFF_FUNCTIONAL_3D_HH
#define META_STUFF_FUNCTIONAL_3D_HH

#include <cmath>

#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "core/blockStatistics.h"
#include "core/dynamicsIdentifiers.h"
#include "core/plbDebug.h"
#include "core/util.h"
#include "dataProcessors/metaStuffFunctional3D.h"
#include "dataProcessors/metaStuffHelper.h"

namespace plb {

/* ******** StoreDynamicsFunctional3D ************************************ */

template <typename T, template <typename U> class Descriptor>
StoreDynamicsFunctional3D<T, Descriptor>::StoreDynamicsFunctional3D() :
    maxChainLengthId(this->getStatistics().subscribeMax())
{ }

template <typename T, template <typename U> class Descriptor>
void StoreDynamicsFunctional3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    BlockLattice3D<T, Descriptor> &lattice =
        *dynamic_cast<BlockLattice3D<T, Descriptor> *>(blocks[0]);
    AtomicContainerBlock3D &container = *dynamic_cast<AtomicContainerBlock3D *>(blocks[1]);
    StoreDynamicsID *storeID = new StoreDynamicsID;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                std::vector<int> chain;
                constructIdChain(lattice.get(iX, iY, iZ).getDynamics(), chain);
                storeID->addIdChain(chain);
                this->getStatistics().gatherMax(maxChainLengthId, (double)chain.size());
            }
        }
    }
    storeID->startIterations();
    container.setData(storeID);
}

template <typename T, template <typename U> class Descriptor>
StoreDynamicsFunctional3D<T, Descriptor> *StoreDynamicsFunctional3D<T, Descriptor>::clone() const
{
    return new StoreDynamicsFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void StoreDynamicsFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT StoreDynamicsFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
pluint StoreDynamicsFunctional3D<T, Descriptor>::getMaxChainLength() const
{
    double maximum = this->getStatistics().getMax(maxChainLengthId);
    return (pluint)(.5 + maximum);
}

/* ******** ExtractDynamicsChainFunctional3D ************************************ */

template <typename T, template <typename U> class Descriptor>
ExtractDynamicsChainFunctional3D<T, Descriptor>::ExtractDynamicsChainFunctional3D(
    typename ExtractDynamicsChainFunctional3D<T, Descriptor>::DMap const &dynamicsMap_,
    pluint maxChainSize_) :
    dynamicsMap(dynamicsMap_), maxChainSize(maxChainSize_)
{ }

template <typename T, template <typename U> class Descriptor>
void ExtractDynamicsChainFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &mask)
{
    Dot3D offset = computeRelativeDisplacement(lattice, mask);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                std::vector<int> chain;
                constructIdChain(lattice.get(iX, iY, iZ).getDynamics(), chain);
                util::extendVectorSize(chain, maxChainSize);
                mask.get(iX + offset.x, iY + offset.y, iZ + offset.z) = dynamicsMap[chain];
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ExtractDynamicsChainFunctional3D<T, Descriptor>
    *ExtractDynamicsChainFunctional3D<T, Descriptor>::clone() const
{
    return new ExtractDynamicsChainFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void ExtractDynamicsChainFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT ExtractDynamicsChainFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** ExtractTopMostDynamicsFunctional3D ************************************ */

template <typename T, template <typename U> class Descriptor>
void ExtractTopMostDynamicsFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &mask)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                mask.get(iX, iY, iZ) = lattice.get(iX, iY, iZ).getDynamics().getId();
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ExtractTopMostDynamicsFunctional3D<T, Descriptor>
    *ExtractTopMostDynamicsFunctional3D<T, Descriptor>::clone() const
{
    return new ExtractTopMostDynamicsFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void ExtractTopMostDynamicsFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT ExtractTopMostDynamicsFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** ExtractBottomMostDynamicsFunctional3D ************************************ */

template <typename T, template <typename U> class Descriptor>
void ExtractBottomMostDynamicsFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &mask)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                mask.get(iX, iY, iZ) =
                    getBottomMostDynamics(lattice.get(iX, iY, iZ).getDynamics()).getId();
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ExtractBottomMostDynamicsFunctional3D<T, Descriptor>
    *ExtractBottomMostDynamicsFunctional3D<T, Descriptor>::clone() const
{
    return new ExtractBottomMostDynamicsFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void ExtractBottomMostDynamicsFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT ExtractBottomMostDynamicsFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** AssignEntireCellFunctional3D ************************************ */

template <typename T, template <typename U> class Descriptor>
void AssignEntireCellFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &sourceLattice,
    BlockLattice3D<T, Descriptor> &destinationLattice)
{
    std::vector<char> data;
    Dot3D offset = computeRelativeDisplacement(sourceLattice, destinationLattice);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> const &cell = sourceLattice.get(iX, iY, iZ);

                pluint pos = data.size();
                data.resize(pos + CellInfo<T, Descriptor>::n * sizeof(T));
                cell.serialize(&data[pos]);

                HierarchicSerializer serializer(data, cell.getDynamics().getId());
                cell.getDynamics().serialize(serializer);
            }
        }
    }

    pluint pos = 0;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iX_ = iX + offset.x;
                plint iY_ = iY + offset.y;
                plint iZ_ = iZ + offset.z;
                Cell<T, Descriptor> &cell = destinationLattice.get(iX_, iY_, iZ_);

                cell.unSerialize(&data[pos]);
                pos += CellInfo<T, Descriptor>::n * sizeof(T);

                HierarchicUnserializer unserializer(data, pos);
                destinationLattice.attributeDynamics(
                    iX_, iY_, iZ_,
                    meta::dynamicsRegistration<T, Descriptor>().generate(unserializer));
                pos = unserializer.getCurrentPos();
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
AssignEntireCellFunctional3D<T, Descriptor> *AssignEntireCellFunctional3D<T, Descriptor>::clone()
    const
{
    return new AssignEntireCellFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void AssignEntireCellFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT AssignEntireCellFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

}  // namespace plb

#endif  // META_STUFF_FUNCTIONAL_3D_HH
