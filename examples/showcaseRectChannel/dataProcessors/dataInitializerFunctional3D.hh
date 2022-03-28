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
 * Functionals for domain initialization -- generic implementation.
 */
#ifndef DATA_INITIALIZER_FUNCTIONAL_3D_HH
#define DATA_INITIALIZER_FUNCTIONAL_3D_HH

#include <algorithm>
#include <limits>

#include "atomicBlock/blockLattice3D.h"
#include "core/cell.h"
#include "dataProcessors/dataInitializerFunctional3D.h"
#include "latticeBoltzmann/externalFieldAccess.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "multiGrid/multiGridUtil.h"
#include "sitmo/prng_engine.hpp"

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Initialization of the block-lattice *************** */
/* ******************************************************************* */

template <typename T, template <class U> class Descriptor>
OneCellFunctional3D<T, Descriptor>::~OneCellFunctional3D()
{ }

template <typename T, template <class U> class Descriptor>
BlockDomain::DomainT OneCellFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <class U> class Descriptor>
void OneCellFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, template <class U> class Descriptor>
void OneCellFunctional3D<T, Descriptor>::setscale(int dxScale, int dtScale)
{ }

template <typename T, template <class U> class Descriptor>
OneCellIndexedFunctional3D<T, Descriptor>::~OneCellIndexedFunctional3D()
{ }

template <typename T, template <class U> class Descriptor>
BlockDomain::DomainT OneCellIndexedFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <class U> class Descriptor>
void OneCellIndexedFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, template <class U> class Descriptor>
void OneCellIndexedFunctional3D<T, Descriptor>::setscale(int dxScale, int dtScale)
{ }

template <typename T, template <class U> class Descriptor>
OneCellIndexedWithRandFunctional3D<T, Descriptor>::~OneCellIndexedWithRandFunctional3D()
{ }

template <typename T, template <class U> class Descriptor>
BlockDomain::DomainT OneCellIndexedWithRandFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <class U> class Descriptor>
void OneCellIndexedWithRandFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, template <class U> class Descriptor>
void OneCellIndexedWithRandFunctional3D<T, Descriptor>::setscale(int dxScale, int dtScale)
{ }

template <typename T, template <class U> class Descriptor>
GenericLatticeFunctional3D<T, Descriptor>::GenericLatticeFunctional3D(
    OneCellFunctional3D<T, Descriptor> *f_) :
    f(f_)
{ }

template <typename T, template <class U> class Descriptor>
GenericLatticeFunctional3D<T, Descriptor>::GenericLatticeFunctional3D(
    GenericLatticeFunctional3D<T, Descriptor> const &rhs) :
    f(rhs.f->clone())
{ }

template <typename T, template <class U> class Descriptor>
GenericLatticeFunctional3D<T, Descriptor>::~GenericLatticeFunctional3D()
{
    delete f;
}

template <typename T, template <class U> class Descriptor>
GenericLatticeFunctional3D<T, Descriptor> &GenericLatticeFunctional3D<T, Descriptor>::operator=(
    GenericLatticeFunctional3D<T, Descriptor> const &rhs)
{
    delete f;
    f = rhs.f->clone();
    return *this;
}

template <typename T, template <class U> class Descriptor>
void GenericLatticeFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                f->execute(lattice.get(iX, iY, iZ));
            }
        }
    }
}

template <typename T, template <class U> class Descriptor>
GenericLatticeFunctional3D<T, Descriptor> *GenericLatticeFunctional3D<T, Descriptor>::clone() const
{
    return new GenericLatticeFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <class U> class Descriptor>
void GenericLatticeFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    f->getTypeOfModification(modified);
}

template <typename T, template <class U> class Descriptor>
BlockDomain::DomainT GenericLatticeFunctional3D<T, Descriptor>::appliesTo() const
{
    return f->appliesTo();
}

template <typename T, template <class U> class Descriptor>
void GenericLatticeFunctional3D<T, Descriptor>::setscale(int dxScale, int dtScale)
{
    f->setscale(dxScale, dtScale);
}

template <typename T, template <class U> class Descriptor>
GenericIndexedLatticeFunctional3D<T, Descriptor>::GenericIndexedLatticeFunctional3D(
    OneCellIndexedFunctional3D<T, Descriptor> *f_) :
    f(f_)
{ }

template <typename T, template <class U> class Descriptor>
GenericIndexedLatticeFunctional3D<T, Descriptor>::GenericIndexedLatticeFunctional3D(
    GenericIndexedLatticeFunctional3D<T, Descriptor> const &rhs) :
    f(rhs.f->clone())
{ }

template <typename T, template <class U> class Descriptor>
GenericIndexedLatticeFunctional3D<T, Descriptor>::~GenericIndexedLatticeFunctional3D()
{
    delete f;
}

template <typename T, template <class U> class Descriptor>
GenericIndexedLatticeFunctional3D<T, Descriptor>
    &GenericIndexedLatticeFunctional3D<T, Descriptor>::operator=(
        GenericIndexedLatticeFunctional3D<T, Descriptor> const &rhs)
{
    delete f;
    f = rhs.f->clone();
    return *this;
}

template <typename T, template <class U> class Descriptor>
void GenericIndexedLatticeFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    Dot3D relativeOffset = lattice.getLocation();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                f->execute(
                    iX + relativeOffset.x, iY + relativeOffset.y, iZ + relativeOffset.z,
                    lattice.get(iX, iY, iZ));
            }
        }
    }
}

template <typename T, template <class U> class Descriptor>
GenericIndexedLatticeFunctional3D<T, Descriptor>
    *GenericIndexedLatticeFunctional3D<T, Descriptor>::clone() const
{
    return new GenericIndexedLatticeFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <class U> class Descriptor>
void GenericIndexedLatticeFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    f->getTypeOfModification(modified);
}

template <typename T, template <class U> class Descriptor>
BlockDomain::DomainT GenericIndexedLatticeFunctional3D<T, Descriptor>::appliesTo() const
{
    return f->appliesTo();
}

template <typename T, template <class U> class Descriptor>
void GenericIndexedLatticeFunctional3D<T, Descriptor>::setscale(int dxScale, int dtScale)
{
    f->setscale(dxScale, dtScale);
}

/* *************** Class GenericIndexedWithRandLatticeFunctional3D ************* */

template <typename T, template <class U> class Descriptor>
GenericIndexedWithRandLatticeFunctional3D<T, Descriptor>::GenericIndexedWithRandLatticeFunctional3D(
    OneCellIndexedWithRandFunctional3D<T, Descriptor> *f_, Box3D boundingBox_,
    uint32_t const *seed_) :
    f(f_), nY(boundingBox_.getNy()), nZ(boundingBox_.getNz()), seed(seed_)
{ }

template <typename T, template <class U> class Descriptor>
GenericIndexedWithRandLatticeFunctional3D<T, Descriptor>::GenericIndexedWithRandLatticeFunctional3D(
    GenericIndexedWithRandLatticeFunctional3D<T, Descriptor> const &rhs) :
    f(rhs.f->clone()), nY(rhs.nY), nZ(rhs.nZ), seed(rhs.seed)
{ }

template <typename T, template <class U> class Descriptor>
GenericIndexedWithRandLatticeFunctional3D<
    T, Descriptor>::~GenericIndexedWithRandLatticeFunctional3D()
{
    delete f;
}

template <typename T, template <class U> class Descriptor>
GenericIndexedWithRandLatticeFunctional3D<T, Descriptor>
    &GenericIndexedWithRandLatticeFunctional3D<T, Descriptor>::operator=(
        GenericIndexedWithRandLatticeFunctional3D<T, Descriptor> const &rhs)
{
    delete f;
    f = rhs.f->clone();
    return *this;
}

template <typename T, template <class U> class Descriptor>
void GenericIndexedWithRandLatticeFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    Dot3D relativeOffset = lattice.getLocation();
    sitmo::prng_engine eng(*seed);
    plint rng_index = 0;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint globalX = nY * (iX + relativeOffset.x);
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint globalY = nZ * (iY + relativeOffset.y + globalX);
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint globalZ = iZ + relativeOffset.z + globalY;
                PLB_ASSERT(globalZ >= rng_index);
                if (globalZ > rng_index) {
                    eng.discard(globalZ - rng_index);
                    rng_index = globalZ;
                }
                T rand_val = (T)eng() / ((T)sitmo::prng_engine::max() + 1.0);
                ++rng_index;
                f->execute(
                    iX + relativeOffset.x, iY + relativeOffset.y, iZ + relativeOffset.z, rand_val,
                    lattice.get(iX, iY, iZ));
            }
        }
    }
}

template <typename T, template <class U> class Descriptor>
GenericIndexedWithRandLatticeFunctional3D<T, Descriptor>
    *GenericIndexedWithRandLatticeFunctional3D<T, Descriptor>::clone() const
{
    return new GenericIndexedWithRandLatticeFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <class U> class Descriptor>
void GenericIndexedWithRandLatticeFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    f->getTypeOfModification(modified);
}

template <typename T, template <class U> class Descriptor>
BlockDomain::DomainT GenericIndexedWithRandLatticeFunctional3D<T, Descriptor>::appliesTo() const
{
    return f->appliesTo();
}

template <typename T, template <class U> class Descriptor>
void GenericIndexedWithRandLatticeFunctional3D<T, Descriptor>::setscale(int dxScale, int dtScale)
{
    f->setscale(dxScale, dtScale);
}

/* *************** Class InstantiateDynamicsFunctional3D ************* */

template <typename T, template <typename U> class Descriptor>
InstantiateDynamicsFunctional3D<T, Descriptor>::InstantiateDynamicsFunctional3D(
    Dynamics<T, Descriptor> *dynamics_) :
    dynamics(dynamics_)
{ }

template <typename T, template <typename U> class Descriptor>
InstantiateDynamicsFunctional3D<T, Descriptor>::InstantiateDynamicsFunctional3D(
    InstantiateDynamicsFunctional3D<T, Descriptor> const &rhs) :
    dynamics(rhs.dynamics->clone())
{ }

template <typename T, template <typename U> class Descriptor>
InstantiateDynamicsFunctional3D<T, Descriptor>
    &InstantiateDynamicsFunctional3D<T, Descriptor>::operator=(
        InstantiateDynamicsFunctional3D<T, Descriptor> const &rhs)
{
    delete dynamics;
    dynamics = rhs.dynamics->clone();
    return *this;
}

template <typename T, template <typename U> class Descriptor>
InstantiateDynamicsFunctional3D<T, Descriptor>::~InstantiateDynamicsFunctional3D()
{
    delete dynamics;
}

template <typename T, template <typename U> class Descriptor>
void InstantiateDynamicsFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                lattice.attributeDynamics(iX, iY, iZ, dynamics->clone());
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT InstantiateDynamicsFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void InstantiateDynamicsFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dataStructure;
}

template <typename T, template <typename U> class Descriptor>
InstantiateDynamicsFunctional3D<T, Descriptor>
    *InstantiateDynamicsFunctional3D<T, Descriptor>::clone() const
{
    return new InstantiateDynamicsFunctional3D<T, Descriptor>(*this);
}

/* *************** Class InstantiateDynamicsInBulkAndEnvelopeFunctional3D ************* */

template <typename T, template <typename U> class Descriptor>
InstantiateDynamicsInBulkAndEnvelopeFunctional3D<T, Descriptor>::
    InstantiateDynamicsInBulkAndEnvelopeFunctional3D(Dynamics<T, Descriptor> *dynamics_) :
    dynamics(dynamics_)
{ }

template <typename T, template <typename U> class Descriptor>
InstantiateDynamicsInBulkAndEnvelopeFunctional3D<T, Descriptor>::
    InstantiateDynamicsInBulkAndEnvelopeFunctional3D(
        InstantiateDynamicsInBulkAndEnvelopeFunctional3D<T, Descriptor> const &rhs) :
    dynamics(rhs.dynamics->clone())
{ }

template <typename T, template <typename U> class Descriptor>
InstantiateDynamicsInBulkAndEnvelopeFunctional3D<T, Descriptor>
    &InstantiateDynamicsInBulkAndEnvelopeFunctional3D<T, Descriptor>::operator=(
        InstantiateDynamicsInBulkAndEnvelopeFunctional3D<T, Descriptor> const &rhs)
{
    delete dynamics;
    dynamics = rhs.dynamics->clone();
    return *this;
}

template <typename T, template <typename U> class Descriptor>
InstantiateDynamicsInBulkAndEnvelopeFunctional3D<
    T, Descriptor>::~InstantiateDynamicsInBulkAndEnvelopeFunctional3D()
{
    delete dynamics;
}

template <typename T, template <typename U> class Descriptor>
void InstantiateDynamicsInBulkAndEnvelopeFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                lattice.attributeDynamics(iX, iY, iZ, dynamics->clone());
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT InstantiateDynamicsInBulkAndEnvelopeFunctional3D<T, Descriptor>::appliesTo()
    const
{
    return BlockDomain::bulkAndEnvelope;  // This is the only difference with the
                                          // InstantiateDynamicsFunctional3D.
}

template <typename T, template <typename U> class Descriptor>
void InstantiateDynamicsInBulkAndEnvelopeFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dataStructure;
}

template <typename T, template <typename U> class Descriptor>
InstantiateDynamicsInBulkAndEnvelopeFunctional3D<T, Descriptor>
    *InstantiateDynamicsInBulkAndEnvelopeFunctional3D<T, Descriptor>::clone() const
{
    return new InstantiateDynamicsInBulkAndEnvelopeFunctional3D<T, Descriptor>(*this);
}

/* ************* Class InstantiateComplexDomainDynamicsFunctional3D ** */

template <typename T, template <typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional3D<T, Descriptor>::
    InstantiateComplexDomainDynamicsFunctional3D(
        Dynamics<T, Descriptor> *dynamics_, DomainFunctional3D *domain_) :
    dynamics(dynamics_), domain(domain_)
{ }

template <typename T, template <typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional3D<T, Descriptor>::
    InstantiateComplexDomainDynamicsFunctional3D(
        InstantiateComplexDomainDynamicsFunctional3D<T, Descriptor> const &rhs) :
    dynamics(rhs.dynamics->clone()), domain(rhs.domain->clone())
{ }

template <typename T, template <typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional3D<T, Descriptor>
    &InstantiateComplexDomainDynamicsFunctional3D<T, Descriptor>::operator=(
        InstantiateComplexDomainDynamicsFunctional3D<T, Descriptor> const &rhs)
{
    delete dynamics;
    dynamics = rhs.dynamics->clone();
    delete domain;
    domain = rhs.domain->clone();
    return *this;
}

template <typename T, template <typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional3D<
    T, Descriptor>::~InstantiateComplexDomainDynamicsFunctional3D()
{
    delete dynamics;
    delete domain;
}

template <typename T, template <typename U> class Descriptor>
void InstantiateComplexDomainDynamicsFunctional3D<T, Descriptor>::process(
    Box3D boundingBox, BlockLattice3D<T, Descriptor> &lattice)
{
    Dot3D relativeOffset = lattice.getLocation();
    for (plint iX = boundingBox.x0; iX <= boundingBox.x1; ++iX) {
        for (plint iY = boundingBox.y0; iY <= boundingBox.y1; ++iY) {
            for (plint iZ = boundingBox.z0; iZ <= boundingBox.z1; ++iZ) {
                if ((*domain)(iX + relativeOffset.x, iY + relativeOffset.y, iZ + relativeOffset.z))
                {
                    lattice.attributeDynamics(iX, iY, iZ, dynamics->clone());
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT InstantiateComplexDomainDynamicsFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void InstantiateComplexDomainDynamicsFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dataStructure;
}

template <typename T, template <typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional3D<T, Descriptor>
    *InstantiateComplexDomainDynamicsFunctional3D<T, Descriptor>::clone() const
{
    return new InstantiateComplexDomainDynamicsFunctional3D<T, Descriptor>(*this);
}

/* ************* Class InstantiateDotDynamicsFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
InstantiateDotDynamicsFunctional3D<T, Descriptor>::InstantiateDotDynamicsFunctional3D(
    Dynamics<T, Descriptor> *dynamics_) :
    dynamics(dynamics_)
{ }

template <typename T, template <typename U> class Descriptor>
InstantiateDotDynamicsFunctional3D<T, Descriptor>::InstantiateDotDynamicsFunctional3D(
    InstantiateDotDynamicsFunctional3D<T, Descriptor> const &rhs) :
    dynamics(rhs.dynamics->clone())
{ }

template <typename T, template <typename U> class Descriptor>
InstantiateDotDynamicsFunctional3D<T, Descriptor>
    &InstantiateDotDynamicsFunctional3D<T, Descriptor>::operator=(
        InstantiateDotDynamicsFunctional3D<T, Descriptor> const &rhs)
{
    delete dynamics;
    dynamics = rhs.dynamics->clone();
    return *this;
}

template <typename T, template <typename U> class Descriptor>
InstantiateDotDynamicsFunctional3D<T, Descriptor>::~InstantiateDotDynamicsFunctional3D()
{
    delete dynamics;
}

template <typename T, template <typename U> class Descriptor>
void InstantiateDotDynamicsFunctional3D<T, Descriptor>::process(
    DotList3D const &dotList, BlockLattice3D<T, Descriptor> &lattice)
{
    for (plint iDot = 0; iDot < dotList.getN(); ++iDot) {
        Dot3D const &dot = dotList.getDot(iDot);
        lattice.attributeDynamics(dot.x, dot.y, dot.z, dynamics->clone());
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT InstantiateDotDynamicsFunctional3D<T, Descriptor>::appliesTo() const
{
    // Dynamics needs to be instantiated everywhere, including envelope.
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void InstantiateDotDynamicsFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dataStructure;
}

template <typename T, template <typename U> class Descriptor>
InstantiateDotDynamicsFunctional3D<T, Descriptor>
    *InstantiateDotDynamicsFunctional3D<T, Descriptor>::clone() const
{
    return new InstantiateDotDynamicsFunctional3D<T, Descriptor>(*this);
}

/* ************* Class DynamicsFromMaskFunctional3D ************************ */

template <typename T, template <typename U> class Descriptor>
DynamicsFromMaskFunctional3D<T, Descriptor>::DynamicsFromMaskFunctional3D(
    Dynamics<T, Descriptor> *dynamics_, bool whichFlag_) :
    dynamics(dynamics_), whichFlag(whichFlag_)
{ }

template <typename T, template <typename U> class Descriptor>
DynamicsFromMaskFunctional3D<T, Descriptor>::DynamicsFromMaskFunctional3D(
    DynamicsFromMaskFunctional3D<T, Descriptor> const &rhs) :
    dynamics(rhs.dynamics->clone()), whichFlag(rhs.whichFlag)
{ }

template <typename T, template <typename U> class Descriptor>
DynamicsFromMaskFunctional3D<T, Descriptor> &DynamicsFromMaskFunctional3D<T, Descriptor>::operator=(
    DynamicsFromMaskFunctional3D<T, Descriptor> const &rhs)
{
    delete dynamics;
    dynamics = rhs.dynamics->clone();
    whichFlag = rhs.whichFlag;
    return *this;
}

template <typename T, template <typename U> class Descriptor>
DynamicsFromMaskFunctional3D<T, Descriptor>::~DynamicsFromMaskFunctional3D()
{
    delete dynamics;
}

template <typename T, template <typename U> class Descriptor>
void DynamicsFromMaskFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<bool> &mask)
{
    Dot3D offset = computeRelativeDisplacement(lattice, mask);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                bool flag = mask.get(iX + offset.x, iY + offset.y, iZ + offset.z);
                if (util::boolIsEqual(flag, whichFlag)) {
                    lattice.attributeDynamics(iX, iY, iZ, dynamics->clone());
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT DynamicsFromMaskFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void DynamicsFromMaskFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dataStructure;
    modified[1] = modif::nothing;
}

template <typename T, template <typename U> class Descriptor>
DynamicsFromMaskFunctional3D<T, Descriptor> *DynamicsFromMaskFunctional3D<T, Descriptor>::clone()
    const
{
    return new DynamicsFromMaskFunctional3D<T, Descriptor>(*this);
}

/* ************* Class DynamicsFromIntMaskFunctional3D ************************ */

template <typename T, template <typename U> class Descriptor>
DynamicsFromIntMaskFunctional3D<T, Descriptor>::DynamicsFromIntMaskFunctional3D(
    Dynamics<T, Descriptor> *dynamics_, int whichFlag_) :
    dynamics(dynamics_), whichFlag(whichFlag_)
{ }

template <typename T, template <typename U> class Descriptor>
DynamicsFromIntMaskFunctional3D<T, Descriptor>::DynamicsFromIntMaskFunctional3D(
    DynamicsFromIntMaskFunctional3D<T, Descriptor> const &rhs) :
    dynamics(rhs.dynamics->clone()), whichFlag(rhs.whichFlag)
{ }

template <typename T, template <typename U> class Descriptor>
DynamicsFromIntMaskFunctional3D<T, Descriptor>
    &DynamicsFromIntMaskFunctional3D<T, Descriptor>::operator=(
        DynamicsFromIntMaskFunctional3D<T, Descriptor> const &rhs)
{
    delete dynamics;
    dynamics = rhs.dynamics->clone();
    whichFlag = rhs.whichFlag;
    return *this;
}

template <typename T, template <typename U> class Descriptor>
DynamicsFromIntMaskFunctional3D<T, Descriptor>::~DynamicsFromIntMaskFunctional3D()
{
    delete dynamics;
}

template <typename T, template <typename U> class Descriptor>
void DynamicsFromIntMaskFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &mask)
{
    Dot3D offset = computeRelativeDisplacement(lattice, mask);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                int flag = mask.get(iX + offset.x, iY + offset.y, iZ + offset.z);
                if (flag == whichFlag) {
                    lattice.attributeDynamics(iX, iY, iZ, dynamics->clone());
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT DynamicsFromIntMaskFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void DynamicsFromIntMaskFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dataStructure;
    modified[1] = modif::nothing;
}

template <typename T, template <typename U> class Descriptor>
DynamicsFromIntMaskFunctional3D<T, Descriptor>
    *DynamicsFromIntMaskFunctional3D<T, Descriptor>::clone() const
{
    return new DynamicsFromIntMaskFunctional3D<T, Descriptor>(*this);
}

/* ************* Class RecomposeFromOrderZeroVariablesFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
void RecomposeFromOrderZeroVariablesFunctional3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    BlockLattice3D<T, Descriptor> &lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> &>(*atomicBlocks[0]);
    ScalarField3D<T> const &rhoField = dynamic_cast<ScalarField3D<T> const &>(*atomicBlocks[1]);
    TensorField3D<T, 3> const &uField = dynamic_cast<TensorField3D<T, 3> const &>(*atomicBlocks[2]);
    TensorField3D<T, Descriptor<T>::q> const &fNeqField =
        dynamic_cast<TensorField3D<T, Descriptor<T>::q> const &>(*atomicBlocks[3]);

    Dot3D offset1 = computeRelativeDisplacement(lattice, rhoField);
    Dot3D offset2 = computeRelativeDisplacement(lattice, uField);
    Dot3D offset3 = computeRelativeDisplacement(lattice, fNeqField);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
                T rho = rhoField.get(iX + offset1.x, iY + offset1.y, iZ + offset1.z);
                Array<T, 3> const &u = uField.get(iX + offset2.x, iY + offset2.y, iZ + offset2.z);
                Array<T, Descriptor<T>::q> const &fNeq =
                    fNeqField.get(iX + offset3.x, iY + offset3.y, iZ + offset3.z);

                plint recomposeOrder = 0;
                std::vector<T> rawData(cell.getDynamics().numDecomposedVariables(recomposeOrder));

                // Convert rho --> rhoBar.
                rawData[0] = Descriptor<T>::rhoBar(rho);

                // Convert u --> j
                rawData[1] = rho * u[0];
                rawData[2] = rho * u[1];
                rawData[3] = rho * u[2];

                for (pluint iPop = 1 + Descriptor<T>::d;
                     iPop < 1 + Descriptor<T>::d + Descriptor<T>::q; ++iPop) {
                    rawData[iPop] = fNeq[iPop - (1 + Descriptor<T>::d)];
                }

                for (pluint iPop = 1 + Descriptor<T>::d + Descriptor<T>::q; iPop < rawData.size();
                     ++iPop) {
                    rawData[iPop] =
                        *cell.getExternal(iPop - (1 + Descriptor<T>::d + Descriptor<T>::q));
                }

                // Recompose the cell.
                cell.getDynamics().recompose(cell, rawData, recomposeOrder);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
RecomposeFromOrderZeroVariablesFunctional3D<T, Descriptor>
    *RecomposeFromOrderZeroVariablesFunctional3D<T, Descriptor>::clone() const
{
    return new RecomposeFromOrderZeroVariablesFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT RecomposeFromOrderZeroVariablesFunctional3D<T, Descriptor>::appliesTo() const
{
    // We could directly apply to the envelope too, but let's keep it
    //   bulk-only for future compatibility.
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void RecomposeFromOrderZeroVariablesFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
    modified[3] = modif::nothing;
}

/* ************* Class RecomposeFromFlowVariablesFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
void RecomposeFromFlowVariablesFunctional3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    BlockLattice3D<T, Descriptor> &lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> &>(*atomicBlocks[0]);
    ScalarField3D<T> const &rhoField = dynamic_cast<ScalarField3D<T> const &>(*atomicBlocks[1]);
    TensorField3D<T, 3> const &uField = dynamic_cast<TensorField3D<T, 3> const &>(*atomicBlocks[2]);
    TensorField3D<T, 6> const &SField = dynamic_cast<TensorField3D<T, 6> const &>(*atomicBlocks[3]);

    Dot3D offset1 = computeRelativeDisplacement(lattice, rhoField);
    Dot3D offset2 = computeRelativeDisplacement(lattice, uField);
    Dot3D offset3 = computeRelativeDisplacement(lattice, SField);

    std::vector<T> rawData(10 + Descriptor<T>::ExternalField::numScalars);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
                T rho = rhoField.get(iX + offset1.x, iY + offset1.y, iZ + offset1.z);
                Array<T, 3> const &u = uField.get(iX + offset2.x, iY + offset2.y, iZ + offset2.z);
                Array<T, 6> const &S = SField.get(iX + offset3.x, iY + offset3.y, iZ + offset3.z);

                // Convert rho --> rhoBar.
                rawData[0] = Descriptor<T>::rhoBar(rho);

                // Convert u --> j
                rawData[1] = rho * u[0];
                rawData[2] = rho * u[1];
                rawData[3] = rho * u[2];

                // Convert S --> PiNeq.
                T omega = cell.getDynamics().getOmega();
                T prefactor = -Descriptor<T>::cs2 * rho * (T)2 / omega;
                rawData[4] = S[0] * prefactor;
                rawData[5] = S[1] * prefactor;
                rawData[6] = S[2] * prefactor;
                rawData[7] = S[3] * prefactor;
                rawData[8] = S[4] * prefactor;
                rawData[9] = S[5] * prefactor;

                // Recompose the cell.
                plint recomposeOrder = 1;
                cell.getDynamics().recompose(cell, rawData, recomposeOrder);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
RecomposeFromFlowVariablesFunctional3D<T, Descriptor>
    *RecomposeFromFlowVariablesFunctional3D<T, Descriptor>::clone() const
{
    return new RecomposeFromFlowVariablesFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT RecomposeFromFlowVariablesFunctional3D<T, Descriptor>::appliesTo() const
{
    // We could directly apply to the envelope too, but let's keep it
    //   bulk-only for future compatibility.
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void RecomposeFromFlowVariablesFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
    modified[3] = modif::nothing;
}

/* ************* Class AssignOmegaFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
AssignOmegaFunctional3D<T, Descriptor>::AssignOmegaFunctional3D(T omega_) : omega(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
void AssignOmegaFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    // Define dimensions of a viscosity.
    int dimDx = 2;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx, this->getDtScale(), dimDt);
    T nu_cs2 = (T)1 / omega - (T)1 / (T)2;
    T scaledOmega = (T)1 / (scaleFactor * nu_cs2 + (T)1 / (T)2);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                lattice.get(iX, iY, iZ).getDynamics().setOmega(scaledOmega);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
AssignOmegaFunctional3D<T, Descriptor> *AssignOmegaFunctional3D<T, Descriptor>::clone() const
{
    return new AssignOmegaFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT AssignOmegaFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void AssignOmegaFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;
}

/* ************* Class AssignScalarFieldOmegaFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
AssignScalarFieldOmegaFunctional3D<T, Descriptor>::AssignScalarFieldOmegaFunctional3D()
{ }

template <typename T, template <typename U> class Descriptor>
void AssignScalarFieldOmegaFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &omega)
{
    // Define dimensions of a viscosity.
    int dimDx = 2;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx, this->getDtScale(), dimDt);

    Dot3D off = computeRelativeDisplacement(lattice, omega);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = iX + off.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY = iY + off.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint oZ = iZ + off.z;

                T nu_cs2 = (T)1 / omega.get(oX, oY, oZ) - (T)1 / (T)2;
                T scaledOmega = (T)1 / (scaleFactor * nu_cs2 + (T)1 / (T)2);

                lattice.get(iX, iY, iZ).getDynamics().setOmega(scaledOmega);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
AssignScalarFieldOmegaFunctional3D<T, Descriptor>
    *AssignScalarFieldOmegaFunctional3D<T, Descriptor>::clone() const
{
    return new AssignScalarFieldOmegaFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT AssignScalarFieldOmegaFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void AssignScalarFieldOmegaFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;
    modified[1] = modif::nothing;
}

/* ************* Class SetConstBoundaryVelocityFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
SetConstBoundaryVelocityFunctional3D<T, Descriptor>::SetConstBoundaryVelocityFunctional3D(
    Array<T, Descriptor<T>::d> velocity) :
    u(velocity)
{ }

template <typename T, template <typename U> class Descriptor>
void SetConstBoundaryVelocityFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    int dimDx = 1;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx, this->getDtScale(), dimDt);
    Array<T, 3> scaledU = u * scaleFactor;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                lattice.get(iX, iY, iZ).defineVelocity(scaledU);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
SetConstBoundaryVelocityFunctional3D<T, Descriptor>
    *SetConstBoundaryVelocityFunctional3D<T, Descriptor>::clone() const
{
    return new SetConstBoundaryVelocityFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT SetConstBoundaryVelocityFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void SetConstBoundaryVelocityFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::allVariables;
}

/* ************* Class SetConstBoundaryVelocityWithTensorForceFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
SetConstBoundaryVelocityWithTensorForceFunctional3D<T, Descriptor>::
    SetConstBoundaryVelocityWithTensorForceFunctional3D(Array<T, Descriptor<T>::d> velocity) :
    u(velocity)
{ }

template <typename T, template <typename U> class Descriptor>
void SetConstBoundaryVelocityWithTensorForceFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, TensorField3D<T, Descriptor<T>::d> &force)
{
    Dot3D ofs = computeRelativeDisplacement(lattice, force);
    int dimDx = 1;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx, this->getDtScale(), dimDt);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, 3> scaledU =
                    scaleFactor * (u - (T)0.5 * force.get(iX + ofs.x, iY + ofs.y, iZ + ofs.z));
                lattice.get(iX, iY, iZ).defineVelocity(scaledU);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
SetConstBoundaryVelocityWithTensorForceFunctional3D<T, Descriptor>
    *SetConstBoundaryVelocityWithTensorForceFunctional3D<T, Descriptor>::clone() const
{
    return new SetConstBoundaryVelocityWithTensorForceFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT SetConstBoundaryVelocityWithTensorForceFunctional3D<T, Descriptor>::appliesTo()
    const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void SetConstBoundaryVelocityWithTensorForceFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::allVariables;  // Lattice.
    modified[1] = modif::nothing;       // Force.
}

/* ************* Class SetConstBoundaryVelocityWithForceFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
SetConstBoundaryVelocityWithForceFunctional3D<T, Descriptor>::
    SetConstBoundaryVelocityWithForceFunctional3D(
        Array<T, Descriptor<T>::d> force, Array<T, Descriptor<T>::d> velocity) :
    f(force), u(velocity)
{ }

template <typename T, template <typename U> class Descriptor>
void SetConstBoundaryVelocityWithForceFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    int dimDx = 1;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx, this->getDtScale(), dimDt);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, Descriptor<T>::d> scaledU = scaleFactor * (u - (T)0.5 * f);
                lattice.get(iX, iY, iZ).defineVelocity(scaledU);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
SetConstBoundaryVelocityWithForceFunctional3D<T, Descriptor>
    *SetConstBoundaryVelocityWithForceFunctional3D<T, Descriptor>::clone() const
{
    return new SetConstBoundaryVelocityWithForceFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT SetConstBoundaryVelocityWithForceFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void SetConstBoundaryVelocityWithForceFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::allVariables;  // Lattice.
}

/* ************* Class SetCustomBoundaryVelocityWithTensorForceFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class VelocityFunction>
SetCustomBoundaryVelocityWithTensorForceFunctional3D<T, Descriptor, VelocityFunction>::
    SetCustomBoundaryVelocityWithTensorForceFunctional3D(VelocityFunction f_) :
    f(f_)
{ }

template <typename T, template <typename U> class Descriptor, class VelocityFunction>
void SetCustomBoundaryVelocityWithTensorForceFunctional3D<T, Descriptor, VelocityFunction>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, TensorField3D<T, Descriptor<T>::d> &force)
{
    Dot3D offset = computeRelativeDisplacement(lattice, force);
    Dot3D location = lattice.getLocation();
    int dimDx = 1;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx, this->getDtScale(), dimDt);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint x = iX + location.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint y = iY + location.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint z = iZ + location.z;

                Array<T, 3> u;
                f(x, y, z, u);
                Array<T, 3> scaledU =
                    scaleFactor
                    * (u - (T)0.5 * force.get(iX + offset.x, iY + offset.y, iZ + offset.z));
                lattice.get(iX, iY, iZ).defineVelocity(scaledU);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, class VelocityFunction>
SetCustomBoundaryVelocityWithTensorForceFunctional3D<T, Descriptor, VelocityFunction>
    *SetCustomBoundaryVelocityWithTensorForceFunctional3D<T, Descriptor, VelocityFunction>::clone()
        const
{
    return new SetCustomBoundaryVelocityWithTensorForceFunctional3D<
        T, Descriptor, VelocityFunction>(*this);
}

template <typename T, template <typename U> class Descriptor, class VelocityFunction>
BlockDomain::DomainT SetCustomBoundaryVelocityWithTensorForceFunctional3D<
    T, Descriptor, VelocityFunction>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor, class VelocityFunction>
void SetCustomBoundaryVelocityWithTensorForceFunctional3D<T, Descriptor, VelocityFunction>::
    getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::allVariables;  // Lattice.
    modified[1] = modif::nothing;       // Force.
}

/* ************* Class SetCustomBoundaryVelocityWithForceFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class VelocityFunction>
SetCustomBoundaryVelocityWithForceFunctional3D<T, Descriptor, VelocityFunction>::
    SetCustomBoundaryVelocityWithForceFunctional3D(
        Array<T, Descriptor<T>::d> force_, VelocityFunction f_) :
    force(force_), f(f_)
{ }

template <typename T, template <typename U> class Descriptor, class VelocityFunction>
void SetCustomBoundaryVelocityWithForceFunctional3D<T, Descriptor, VelocityFunction>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    Dot3D location = lattice.getLocation();
    int dimDx = 1;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx, this->getDtScale(), dimDt);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint x = iX + location.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint y = iY + location.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint z = iZ + location.z;

                Array<T, Descriptor<T>::d> u;
                f(x, y, z, u);
                Array<T, Descriptor<T>::d> scaledU = scaleFactor * (u - (T)0.5 * force);
                lattice.get(iX, iY, iZ).defineVelocity(scaledU);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, class VelocityFunction>
SetCustomBoundaryVelocityWithForceFunctional3D<T, Descriptor, VelocityFunction>
    *SetCustomBoundaryVelocityWithForceFunctional3D<T, Descriptor, VelocityFunction>::clone() const
{
    return new SetCustomBoundaryVelocityWithForceFunctional3D<T, Descriptor, VelocityFunction>(
        *this);
}

template <typename T, template <typename U> class Descriptor, class VelocityFunction>
BlockDomain::DomainT SetCustomBoundaryVelocityWithForceFunctional3D<
    T, Descriptor, VelocityFunction>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor, class VelocityFunction>
void SetCustomBoundaryVelocityWithForceFunctional3D<T, Descriptor, VelocityFunction>::
    getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::allVariables;  // Lattice.
}

/* ************* Class SetCustomBoundaryVelocityWithCustomForceFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class ForceVelocityFunction>
SetCustomBoundaryVelocityWithCustomForceFunctional3D<T, Descriptor, ForceVelocityFunction>::
    SetCustomBoundaryVelocityWithCustomForceFunctional3D(ForceVelocityFunction f_) :
    f(f_)
{ }

template <typename T, template <typename U> class Descriptor, class ForceVelocityFunction>
void SetCustomBoundaryVelocityWithCustomForceFunctional3D<T, Descriptor, ForceVelocityFunction>::
    process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    Dot3D location = lattice.getLocation();
    int dimDx = 1;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx, this->getDtScale(), dimDt);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint x = iX + location.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint y = iY + location.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint z = iZ + location.z;

                Array<T, Descriptor<T>::d> force;
                Array<T, Descriptor<T>::d> u;
                f(x, y, z, force, u);
                Array<T, Descriptor<T>::d> scaledU = scaleFactor * (u - (T)0.5 * force);
                lattice.get(iX, iY, iZ).defineVelocity(scaledU);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, class ForceVelocityFunction>
SetCustomBoundaryVelocityWithCustomForceFunctional3D<T, Descriptor, ForceVelocityFunction>
    *SetCustomBoundaryVelocityWithCustomForceFunctional3D<
        T, Descriptor, ForceVelocityFunction>::clone() const
{
    return new SetCustomBoundaryVelocityWithCustomForceFunctional3D<
        T, Descriptor, ForceVelocityFunction>(*this);
}

template <typename T, template <typename U> class Descriptor, class ForceVelocityFunction>
BlockDomain::DomainT SetCustomBoundaryVelocityWithCustomForceFunctional3D<
    T, Descriptor, ForceVelocityFunction>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor, class ForceVelocityFunction>
void SetCustomBoundaryVelocityWithCustomForceFunctional3D<T, Descriptor, ForceVelocityFunction>::
    getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::allVariables;  // Lattice.
}

/* ************* Class SetConstBoundaryDensityFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
SetConstBoundaryDensityFunctional3D<T, Descriptor>::SetConstBoundaryDensityFunctional3D(T rho_) :
    rho(rho_)
{ }

template <typename T, template <typename U> class Descriptor>
void SetConstBoundaryDensityFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                lattice.get(iX, iY, iZ).defineDensity(rho);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
SetConstBoundaryDensityFunctional3D<T, Descriptor>
    *SetConstBoundaryDensityFunctional3D<T, Descriptor>::clone() const
{
    return new SetConstBoundaryDensityFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT SetConstBoundaryDensityFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void SetConstBoundaryDensityFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::allVariables;
}

/* ************* Class IniConstEquilibriumFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
IniConstEquilibriumFunctional3D<T, Descriptor>::IniConstEquilibriumFunctional3D(
    T density, Array<T, Descriptor<T>::d> velocity, T temperature) :
    rho(density), rhoBar(Descriptor<T>::rhoBar(density)), u(velocity), thetaBar(temperature - (T)1)
{ }

template <typename T, template <typename U> class Descriptor>
void IniConstEquilibriumFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    int dimDx = 1;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx, this->getDtScale(), dimDt);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, Descriptor<T>::d> f;
                f[0] = getExternalForceComponent(lattice.get(iX, iY, iZ), 0);
                f[1] = getExternalForceComponent(lattice.get(iX, iY, iZ), 1);
                f[2] = getExternalForceComponent(lattice.get(iX, iY, iZ), 2);
                Array<T, Descriptor<T>::d> scaledJ;
                scaledJ[0] = scaleFactor * rho * (u[0] - (T)0.5 * f[0]);
                scaledJ[1] = scaleFactor * rho * (u[1] - (T)0.5 * f[1]);
                scaledJ[2] = scaleFactor * rho * (u[2] - (T)0.5 * f[2]);
                T scaledJsqr = normSqr(scaledJ);
                lattice.get(iX, iY, iZ)
                    .getDynamics()
                    .computeEquilibria(
                        lattice.get(iX, iY, iZ).getRawPopulations(), rhoBar, scaledJ, scaledJsqr,
                        thetaBar);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
IniConstEquilibriumFunctional3D<T, Descriptor>
    *IniConstEquilibriumFunctional3D<T, Descriptor>::clone() const
{
    return new IniConstEquilibriumFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT IniConstEquilibriumFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void IniConstEquilibriumFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************* Class MaskedIniConstEquilibriumFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
MaskedIniConstEquilibriumFunctional3D<T, Descriptor>::MaskedIniConstEquilibriumFunctional3D(
    T density, Array<T, Descriptor<T>::d> velocity, T temperature, int whichFlag_) :
    rho(density),
    rhoBar(Descriptor<T>::rhoBar(density)),
    u(velocity),
    thetaBar(temperature - (T)1),
    whichFlag(whichFlag_)
{ }

template <typename T, template <typename U> class Descriptor>
void MaskedIniConstEquilibriumFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &mask)
{
    int dimDx = 1;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx, this->getDtScale(), dimDt);
    Dot3D offset = computeRelativeDisplacement(lattice, mask);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (mask.get(iX + offset.x, iY + offset.y, iZ + offset.z) == whichFlag) {
                    Array<T, Descriptor<T>::d> f;
                    f[0] = getExternalForceComponent(lattice.get(iX, iY, iZ), 0);
                    f[1] = getExternalForceComponent(lattice.get(iX, iY, iZ), 1);
                    f[2] = getExternalForceComponent(lattice.get(iX, iY, iZ), 2);
                    Array<T, Descriptor<T>::d> scaledJ;
                    scaledJ[0] = scaleFactor * rho * (u[0] - (T)0.5 * f[0]);
                    scaledJ[1] = scaleFactor * rho * (u[1] - (T)0.5 * f[1]);
                    scaledJ[2] = scaleFactor * rho * (u[2] - (T)0.5 * f[2]);
                    T scaledJsqr = normSqr(scaledJ);
                    lattice.get(iX, iY, iZ)
                        .getDynamics()
                        .computeEquilibria(
                            lattice.get(iX, iY, iZ).getRawPopulations(), rhoBar, scaledJ,
                            scaledJsqr, thetaBar);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
MaskedIniConstEquilibriumFunctional3D<T, Descriptor>
    *MaskedIniConstEquilibriumFunctional3D<T, Descriptor>::clone() const
{
    return new MaskedIniConstEquilibriumFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT MaskedIniConstEquilibriumFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void MaskedIniConstEquilibriumFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Lattice.
    modified[1] = modif::nothing;          // Mask.
}

/* ************* Class IniConstTensorForceEquilibriumFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
IniConstTensorForceEquilibriumFunctional3D<T, Descriptor>::
    IniConstTensorForceEquilibriumFunctional3D(
        T density, Array<T, Descriptor<T>::d> velocity, T temperature) :
    rho(density), rhoBar(Descriptor<T>::rhoBar(density)), u(velocity), thetaBar(temperature - (T)1)
{ }

template <typename T, template <typename U> class Descriptor>
void IniConstTensorForceEquilibriumFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, TensorField3D<T, Descriptor<T>::d> &force)
{
    Dot3D offset = computeRelativeDisplacement(lattice, force);
    int dimDx = 1;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx, this->getDtScale(), dimDt);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, Descriptor<T>::d> const &f =
                    force.get(iX + offset.x, iY + offset.y, iZ + offset.z);
                Array<T, Descriptor<T>::d> scaledJ;
                scaledJ[0] = scaleFactor * rho * (u[0] - (T)0.5 * f[0]);
                scaledJ[1] = scaleFactor * rho * (u[1] - (T)0.5 * f[1]);
                scaledJ[2] = scaleFactor * rho * (u[2] - (T)0.5 * f[2]);
                T scaledJsqr = normSqr(scaledJ);
                lattice.get(iX, iY, iZ)
                    .getDynamics()
                    .computeEquilibria(
                        lattice.get(iX, iY, iZ).getRawPopulations(), rhoBar, scaledJ, scaledJsqr,
                        thetaBar);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
IniConstTensorForceEquilibriumFunctional3D<T, Descriptor>
    *IniConstTensorForceEquilibriumFunctional3D<T, Descriptor>::clone() const
{
    return new IniConstTensorForceEquilibriumFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT IniConstTensorForceEquilibriumFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void IniConstTensorForceEquilibriumFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Lattice.
    modified[1] = modif::nothing;          // Force.
}

/* ************* Class IniCustomTensorForceEquilibriumFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
IniCustomTensorForceEquilibriumFunctional3D<T, Descriptor, RhoUFunction>::
    IniCustomTensorForceEquilibriumFunctional3D(RhoUFunction f_, T temperature) :
    f(f_), thetaBar(temperature - (T)1)
{ }

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
void IniCustomTensorForceEquilibriumFunctional3D<T, Descriptor, RhoUFunction>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, TensorField3D<T, Descriptor<T>::d> &force)
{
    Dot3D offset = computeRelativeDisplacement(lattice, force);
    Dot3D location = lattice.getLocation();
    int dimDx = 1;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx, this->getDtScale(), dimDt);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint x = iX + location.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint y = iY + location.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint z = iZ + location.z;

                T rho;
                Array<T, Descriptor<T>::d> u;
                f(x, y, z, rho, u);
                T rhoBar = Descriptor<T>::rhoBar(rho);

                Array<T, Descriptor<T>::d> const &localForce =
                    force.get(iX + offset.x, iY + offset.y, iZ + offset.z);
                Array<T, Descriptor<T>::d> scaledJ;
                scaledJ[0] = scaleFactor * rho * (u[0] - (T)0.5 * localForce[0]);
                scaledJ[1] = scaleFactor * rho * (u[1] - (T)0.5 * localForce[1]);
                scaledJ[2] = scaleFactor * rho * (u[2] - (T)0.5 * localForce[2]);

                T scaledJsqr = normSqr(scaledJ);
                lattice.get(iX, iY, iZ)
                    .getDynamics()
                    .computeEquilibria(
                        lattice.get(iX, iY, iZ).getRawPopulations(), rhoBar, scaledJ, scaledJsqr,
                        thetaBar);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
IniCustomTensorForceEquilibriumFunctional3D<T, Descriptor, RhoUFunction>
    *IniCustomTensorForceEquilibriumFunctional3D<T, Descriptor, RhoUFunction>::clone() const
{
    return new IniCustomTensorForceEquilibriumFunctional3D<T, Descriptor, RhoUFunction>(*this);
}

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
BlockDomain::DomainT
    IniCustomTensorForceEquilibriumFunctional3D<T, Descriptor, RhoUFunction>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
void IniCustomTensorForceEquilibriumFunctional3D<
    T, Descriptor, RhoUFunction>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Lattice.
    modified[1] = modif::nothing;          // Force.
}

/* ************* Class IniCustomTensorForceRandomEquilibriumFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
IniCustomTensorForceRandomEquilibriumFunctional3D<T, Descriptor, RhoUFunction>::
    IniCustomTensorForceRandomEquilibriumFunctional3D(RhoUFunction f_, T temperature) :
    f(f_), thetaBar(temperature - (T)1)
{ }

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
void IniCustomTensorForceRandomEquilibriumFunctional3D<T, Descriptor, RhoUFunction>::
    processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() >= 3);

    plint iBlock = 0;

    BlockLattice3D<T, Descriptor> *lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> *>(blocks[iBlock]);
    PLB_ASSERT(lattice);
    iBlock++;

    TensorField3D<T, Descriptor<T>::d> *force =
        dynamic_cast<TensorField3D<T, Descriptor<T>::d> *>(blocks[iBlock]);
    PLB_ASSERT(force);
    Dot3D ofsF = computeRelativeDisplacement(*lattice, *force);
    iBlock++;

    plint numRandomFields = (plint)blocks.size() - 2;
    std::vector<ScalarField3D<T> *> randomFields(numRandomFields);
    std::vector<Dot3D> ofsR(numRandomFields);
    for (plint i = 0; i < numRandomFields; i++) {
        randomFields[i] = dynamic_cast<ScalarField3D<T> *>(blocks[iBlock]);
        PLB_ASSERT(randomFields[i]);
        ofsR[i] = computeRelativeDisplacement(*lattice, *randomFields[i]);
        iBlock++;
    }

    std::vector<T> randomValues(numRandomFields);

    Dot3D location = lattice->getLocation();

    int dimDx = 1;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx, this->getDtScale(), dimDt);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint x = iX + location.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint y = iY + location.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint z = iZ + location.z;

                T rho;
                Array<T, Descriptor<T>::d> u;

                for (plint i = 0; i < numRandomFields; i++) {
                    randomValues[i] =
                        randomFields[i]->get(iX + ofsR[i].x, iY + ofsR[i].y, iZ + ofsR[i].z);
                }

                f(x, y, z, randomValues, rho, u);

                T rhoBar = Descriptor<T>::rhoBar(rho);

                Array<T, Descriptor<T>::d> const &localForce =
                    force->get(iX + ofsF.x, iY + ofsF.y, iZ + ofsF.z);
                Array<T, Descriptor<T>::d> scaledJ;
                scaledJ[0] = scaleFactor * rho * (u[0] - (T)0.5 * localForce[0]);
                scaledJ[1] = scaleFactor * rho * (u[1] - (T)0.5 * localForce[1]);
                scaledJ[2] = scaleFactor * rho * (u[2] - (T)0.5 * localForce[2]);

                T scaledJsqr = normSqr(scaledJ);
                lattice->get(iX, iY, iZ)
                    .getDynamics()
                    .computeEquilibria(
                        lattice->get(iX, iY, iZ).getRawPopulations(), rhoBar, scaledJ, scaledJsqr,
                        thetaBar);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
IniCustomTensorForceRandomEquilibriumFunctional3D<T, Descriptor, RhoUFunction>
    *IniCustomTensorForceRandomEquilibriumFunctional3D<T, Descriptor, RhoUFunction>::clone() const
{
    return new IniCustomTensorForceRandomEquilibriumFunctional3D<T, Descriptor, RhoUFunction>(
        *this);
}

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
BlockDomain::DomainT IniCustomTensorForceRandomEquilibriumFunctional3D<
    T, Descriptor, RhoUFunction>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
void IniCustomTensorForceRandomEquilibriumFunctional3D<
    T, Descriptor, RhoUFunction>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::staticVariables;  // Lattice.
}

/* ************* Class IniConstForceEquilibriumFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
IniConstForceEquilibriumFunctional3D<T, Descriptor>::IniConstForceEquilibriumFunctional3D(
    Array<T, Descriptor<T>::d> force, T density, Array<T, Descriptor<T>::d> velocity,
    T temperature) :
    f(force),
    rho(density),
    rhoBar(Descriptor<T>::rhoBar(density)),
    u(velocity),
    thetaBar(temperature - (T)1)
{ }

template <typename T, template <typename U> class Descriptor>
void IniConstForceEquilibriumFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    int dimDx = 1;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx, this->getDtScale(), dimDt);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, Descriptor<T>::d> scaledJ;
                scaledJ[0] = scaleFactor * rho * (u[0] - (T)0.5 * f[0]);
                scaledJ[1] = scaleFactor * rho * (u[1] - (T)0.5 * f[1]);
                scaledJ[2] = scaleFactor * rho * (u[2] - (T)0.5 * f[2]);
                T scaledJsqr = normSqr(scaledJ);
                lattice.get(iX, iY, iZ)
                    .getDynamics()
                    .computeEquilibria(
                        lattice.get(iX, iY, iZ).getRawPopulations(), rhoBar, scaledJ, scaledJsqr,
                        thetaBar);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
IniConstForceEquilibriumFunctional3D<T, Descriptor>
    *IniConstForceEquilibriumFunctional3D<T, Descriptor>::clone() const
{
    return new IniConstForceEquilibriumFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT IniConstForceEquilibriumFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void IniConstForceEquilibriumFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Lattice.
}

/* ************* Class IniCustomForceEquilibriumFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
IniCustomForceEquilibriumFunctional3D<T, Descriptor, RhoUFunction>::
    IniCustomForceEquilibriumFunctional3D(
        Array<T, Descriptor<T>::d> force_, RhoUFunction f_, T temperature) :
    force(force_), f(f_), thetaBar(temperature - (T)1)
{ }

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
void IniCustomForceEquilibriumFunctional3D<T, Descriptor, RhoUFunction>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    Dot3D location = lattice.getLocation();
    int dimDx = 1;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx, this->getDtScale(), dimDt);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint x = iX + location.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint y = iY + location.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint z = iZ + location.z;

                T rho;
                Array<T, Descriptor<T>::d> u;
                f(x, y, z, rho, u);
                T rhoBar = Descriptor<T>::rhoBar(rho);

                Array<T, Descriptor<T>::d> scaledJ;
                scaledJ[0] = scaleFactor * rho * (u[0] - (T)0.5 * force[0]);
                scaledJ[1] = scaleFactor * rho * (u[1] - (T)0.5 * force[1]);
                scaledJ[2] = scaleFactor * rho * (u[2] - (T)0.5 * force[2]);

                T scaledJsqr = normSqr(scaledJ);
                lattice.get(iX, iY, iZ)
                    .getDynamics()
                    .computeEquilibria(
                        lattice.get(iX, iY, iZ).getRawPopulations(), rhoBar, scaledJ, scaledJsqr,
                        thetaBar);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
IniCustomForceEquilibriumFunctional3D<T, Descriptor, RhoUFunction>
    *IniCustomForceEquilibriumFunctional3D<T, Descriptor, RhoUFunction>::clone() const
{
    return new IniCustomForceEquilibriumFunctional3D<T, Descriptor, RhoUFunction>(*this);
}

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
BlockDomain::DomainT IniCustomForceEquilibriumFunctional3D<T, Descriptor, RhoUFunction>::appliesTo()
    const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
void IniCustomForceEquilibriumFunctional3D<T, Descriptor, RhoUFunction>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Lattice.
}

/* ************* Class IniCustomForceRandomEquilibriumFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
IniCustomForceRandomEquilibriumFunctional3D<T, Descriptor, RhoUFunction>::
    IniCustomForceRandomEquilibriumFunctional3D(
        Array<T, Descriptor<T>::d> force_, RhoUFunction f_, T temperature) :
    force(force_), f(f_), thetaBar(temperature - (T)1)
{ }

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
void IniCustomForceRandomEquilibriumFunctional3D<T, Descriptor, RhoUFunction>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() >= 2);

    plint iBlock = 0;

    BlockLattice3D<T, Descriptor> *lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> *>(blocks[iBlock]);
    PLB_ASSERT(lattice);
    iBlock++;

    plint numRandomFields = (plint)blocks.size() - 2;
    std::vector<ScalarField3D<T> *> randomFields(numRandomFields);
    std::vector<Dot3D> ofsR(numRandomFields);
    for (plint i = 0; i < numRandomFields; i++) {
        randomFields[i] = dynamic_cast<ScalarField3D<T> *>(blocks[iBlock]);
        PLB_ASSERT(randomFields[i]);
        ofsR[i] = computeRelativeDisplacement(*lattice, *randomFields[i]);
        iBlock++;
    }

    std::vector<T> randomValues(numRandomFields);

    Dot3D location = lattice->getLocation();

    int dimDx = 1;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx, this->getDtScale(), dimDt);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint x = iX + location.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint y = iY + location.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint z = iZ + location.z;

                T rho;
                Array<T, Descriptor<T>::d> u;

                for (plint i = 0; i < numRandomFields; i++) {
                    randomValues[i] =
                        randomFields[i]->get(iX + ofsR[i].x, iY + ofsR[i].y, iZ + ofsR[i].z);
                }

                f(x, y, z, randomValues, rho, u);

                T rhoBar = Descriptor<T>::rhoBar(rho);

                Array<T, Descriptor<T>::d> scaledJ;
                scaledJ[0] = scaleFactor * rho * (u[0] - (T)0.5 * force[0]);
                scaledJ[1] = scaleFactor * rho * (u[1] - (T)0.5 * force[1]);
                scaledJ[2] = scaleFactor * rho * (u[2] - (T)0.5 * force[2]);

                T scaledJsqr = normSqr(scaledJ);
                lattice->get(iX, iY, iZ)
                    .getDynamics()
                    .computeEquilibria(
                        lattice->get(iX, iY, iZ).getRawPopulations(), rhoBar, scaledJ, scaledJsqr,
                        thetaBar);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
IniCustomForceRandomEquilibriumFunctional3D<T, Descriptor, RhoUFunction>
    *IniCustomForceRandomEquilibriumFunctional3D<T, Descriptor, RhoUFunction>::clone() const
{
    return new IniCustomForceRandomEquilibriumFunctional3D<T, Descriptor, RhoUFunction>(*this);
}

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
BlockDomain::DomainT
    IniCustomForceRandomEquilibriumFunctional3D<T, Descriptor, RhoUFunction>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
void IniCustomForceRandomEquilibriumFunctional3D<
    T, Descriptor, RhoUFunction>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::staticVariables;  // Lattice.
}

/* ************* Class IniConstEquilibriumOnDomainFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class DomainFunctional>
IniConstEquilibriumOnDomainFunctional3D<T, Descriptor, DomainFunctional>::
    IniConstEquilibriumOnDomainFunctional3D(
        T density, Array<T, Descriptor<T>::d> velocity, T temperature,
        DomainFunctional const &domain_) :
    rho(density),
    rhoBar(Descriptor<T>::rhoBar(density)),
    u(velocity),
    thetaBar(temperature - (T)1),
    domain(domain_)
{ }

template <typename T, template <typename U> class Descriptor, class DomainFunctional>
void IniConstEquilibriumOnDomainFunctional3D<T, Descriptor, DomainFunctional>::process(
    Box3D bbox, BlockLattice3D<T, Descriptor> &lattice)
{
    int dimDx = 1;
    int dimDt = -1;
    Dot3D absPos = lattice.getLocation();
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx, this->getDtScale(), dimDt);
    for (plint iX = bbox.x0; iX <= bbox.x1; ++iX) {
        for (plint iY = bbox.y0; iY <= bbox.y1; ++iY) {
            for (plint iZ = bbox.z0; iZ <= bbox.z1; ++iZ) {
                if (domain(iX + absPos.x, iY + absPos.y, iZ + absPos.z)) {
                    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                        Array<T, Descriptor<T>::d> f;
                        f[0] = getExternalForceComponent(lattice.get(iX, iY, iZ), 0);
                        f[1] = getExternalForceComponent(lattice.get(iX, iY, iZ), 1);
                        f[2] = getExternalForceComponent(lattice.get(iX, iY, iZ), 2);
                        Array<T, Descriptor<T>::d> scaledJ;
                        scaledJ[0] = scaleFactor * rho * (u[0] - (T)0.5 * f[0]);
                        scaledJ[1] = scaleFactor * rho * (u[1] - (T)0.5 * f[1]);
                        scaledJ[2] = scaleFactor * rho * (u[2] - (T)0.5 * f[2]);
                        T scaledJsqr = normSqr(scaledJ);
                        lattice.get(iX, iY, iZ)[iPop] =
                            lattice.get(iX, iY, iZ)
                                .computeEquilibrium(iPop, rhoBar, scaledJ, scaledJsqr, thetaBar);
                    }
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, class DomainFunctional>
IniConstEquilibriumOnDomainFunctional3D<T, Descriptor, DomainFunctional>
    *IniConstEquilibriumOnDomainFunctional3D<T, Descriptor, DomainFunctional>::clone() const
{
    return new IniConstEquilibriumOnDomainFunctional3D<T, Descriptor, DomainFunctional>(*this);
}

template <typename T, template <typename U> class Descriptor, class DomainFunctional>
void IniConstEquilibriumOnDomainFunctional3D<T, Descriptor, DomainFunctional>::
    getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************* Class StripeOffDensityOffsetFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
StripeOffDensityOffsetFunctional3D<T, Descriptor>::StripeOffDensityOffsetFunctional3D(T deltaRho_) :
    deltaRho(deltaRho_)
{ }

template <typename T, template <typename U> class Descriptor>
void StripeOffDensityOffsetFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
                Dynamics<T, Descriptor> &dynamics = cell.getDynamics();
                plint orderOfDecomposition = 0;
                std::vector<T> rawData;
                dynamics.decompose(cell, rawData, orderOfDecomposition);
                T &rhoBar = rawData[0];
                rhoBar -= deltaRho;
                dynamics.recompose(cell, rawData, orderOfDecomposition);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
StripeOffDensityOffsetFunctional3D<T, Descriptor>
    *StripeOffDensityOffsetFunctional3D<T, Descriptor>::clone() const
{
    return new StripeOffDensityOffsetFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT StripeOffDensityOffsetFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void StripeOffDensityOffsetFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************* Class InstantiateCompositeDynamicsFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional3D<T, Descriptor>::InstantiateCompositeDynamicsFunctional3D(
    CompositeDynamics<T, Descriptor> *compositeDynamics_) :
    compositeDynamics(compositeDynamics_)
{ }

template <typename T, template <typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional3D<T, Descriptor>::InstantiateCompositeDynamicsFunctional3D(
    InstantiateCompositeDynamicsFunctional3D<T, Descriptor> const &rhs) :
    compositeDynamics(rhs.compositeDynamics->clone())
{ }

template <typename T, template <typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional3D<T, Descriptor>
    &InstantiateCompositeDynamicsFunctional3D<T, Descriptor>::operator=(
        InstantiateCompositeDynamicsFunctional3D<T, Descriptor> const &rhs)
{
    delete compositeDynamics;
    compositeDynamics = rhs.compositeDynamics->clone();
    return *this;
}

template <typename T, template <typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional3D<T, Descriptor>::~InstantiateCompositeDynamicsFunctional3D()
{
    delete compositeDynamics;
}

template <typename T, template <typename U> class Descriptor>
void InstantiateCompositeDynamicsFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                lattice.attributeDynamics(
                    iX, iY, iZ,
                    cloneAndInsertAtTopDynamics(
                        lattice.get(iX, iY, iZ).getDynamics(), compositeDynamics->clone()));
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void InstantiateCompositeDynamicsFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dataStructure;
}

template <typename T, template <typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional3D<T, Descriptor>
    *InstantiateCompositeDynamicsFunctional3D<T, Descriptor>::clone() const
{
    return new InstantiateCompositeDynamicsFunctional3D<T, Descriptor>(*this);
}

/* ************* Class SetExternalScalarFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
SetExternalScalarFunctional3D<T, Descriptor>::SetExternalScalarFunctional3D(
    int whichScalar_, T externalScalar_) :
    whichScalar(whichScalar_), externalScalar(externalScalar_)
{
    PLB_ASSERT(whichScalar < Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor>
void SetExternalScalarFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                *lattice.get(iX, iY, iZ).getExternal(whichScalar) = externalScalar;
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT SetExternalScalarFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void SetExternalScalarFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
SetExternalScalarFunctional3D<T, Descriptor> *SetExternalScalarFunctional3D<T, Descriptor>::clone()
    const
{
    return new SetExternalScalarFunctional3D<T, Descriptor>(*this);
}

/* ************* Class SetExternalScalarFromScalarFieldFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
SetExternalScalarFromScalarFieldFunctional3D<
    T, Descriptor>::SetExternalScalarFromScalarFieldFunctional3D(int whichScalar_) :
    whichScalar(whichScalar_)
{
    PLB_ASSERT(whichScalar < Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor>
void SetExternalScalarFromScalarFieldFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalar)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalar);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint oZ = iZ + offset.z;
                *lattice.get(iX, iY, iZ).getExternal(whichScalar) = scalar.get(oX, oY, oZ);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT SetExternalScalarFromScalarFieldFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void SetExternalScalarFromScalarFieldFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
SetExternalScalarFromScalarFieldFunctional3D<T, Descriptor>
    *SetExternalScalarFromScalarFieldFunctional3D<T, Descriptor>::clone() const
{
    return new SetExternalScalarFromScalarFieldFunctional3D<T, Descriptor>(*this);
}

/* ************* Class MaskedSetExternalScalarFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
MaskedSetExternalScalarFunctional3D<T, Descriptor>::MaskedSetExternalScalarFunctional3D(
    int flag_, int whichScalar_, T externalScalar_) :
    flag(flag_), whichScalar(whichScalar_), externalScalar(externalScalar_)
{
    PLB_ASSERT(whichScalar < Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor>
void MaskedSetExternalScalarFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &mask)
{
    Dot3D offset = computeRelativeDisplacement(lattice, mask);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint oZ = iZ + offset.z;
                if (mask.get(oX, oY, oZ) == flag) {
                    *lattice.get(iX, iY, iZ).getExternal(whichScalar) = externalScalar;
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT MaskedSetExternalScalarFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void MaskedSetExternalScalarFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Lattice
    modified[1] = modif::nothing;          // Mask
}

template <typename T, template <typename U> class Descriptor>
MaskedSetExternalScalarFunctional3D<T, Descriptor>
    *MaskedSetExternalScalarFunctional3D<T, Descriptor>::clone() const
{
    return new MaskedSetExternalScalarFunctional3D<T, Descriptor>(*this);
}

/* ************* Class SetGenericExternalScalarFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class Functional>
SetGenericExternalScalarFunctional3D<T, Descriptor, Functional>::
    SetGenericExternalScalarFunctional3D(int whichScalar_, Functional const &functional_) :
    whichScalar(whichScalar_), functional(functional_)
{
    PLB_ASSERT(whichScalar < Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor, class Functional>
void SetGenericExternalScalarFunctional3D<T, Descriptor, Functional>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    Dot3D absOffset = lattice.getLocation();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                *lattice.get(iX, iY, iZ).getExternal(whichScalar) =
                    functional(iX + absOffset.x, iY + absOffset.y, iZ + absOffset.z);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, class Functional>
void SetGenericExternalScalarFunctional3D<T, Descriptor, Functional>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor, class Functional>
SetGenericExternalScalarFunctional3D<T, Descriptor, Functional>
    *SetGenericExternalScalarFunctional3D<T, Descriptor, Functional>::clone() const
{
    return new SetGenericExternalScalarFunctional3D<T, Descriptor, Functional>(*this);
}

/* ************* Class MaskedSetGenericExternalScalarFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class Functional>
MaskedSetGenericExternalScalarFunctional3D<T, Descriptor, Functional>::
    MaskedSetGenericExternalScalarFunctional3D(
        int flag_, int whichScalar_, Functional const &functional_) :
    flag(flag_), whichScalar(whichScalar_), functional(functional_)
{
    PLB_ASSERT(whichScalar < Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor, class Functional>
void MaskedSetGenericExternalScalarFunctional3D<T, Descriptor, Functional>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &mask)
{
    Dot3D offset = computeRelativeDisplacement(lattice, mask);
    Dot3D absOffset = lattice.getLocation();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint oZ = iZ + offset.z;

                if (mask.get(oX, oY, oZ) == flag) {
                    *lattice.get(iX, iY, iZ).getExternal(whichScalar) =
                        functional(iX + absOffset.x, iY + absOffset.y, iZ + absOffset.z);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, class Functional>
void MaskedSetGenericExternalScalarFunctional3D<T, Descriptor, Functional>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Lattice.
    modified[1] = modif::nothing;          // Mask.
}

template <typename T, template <typename U> class Descriptor, class Functional>
MaskedSetGenericExternalScalarFunctional3D<T, Descriptor, Functional>
    *MaskedSetGenericExternalScalarFunctional3D<T, Descriptor, Functional>::clone() const
{
    return new MaskedSetGenericExternalScalarFunctional3D<T, Descriptor, Functional>(*this);
}

/* ************* Class AddToExternalScalarFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
AddToExternalScalarFunctional3D<T, Descriptor>::AddToExternalScalarFunctional3D(
    int whichScalar_, T externalScalar_) :
    whichScalar(whichScalar_), externalScalar(externalScalar_)
{
    PLB_ASSERT(whichScalar < Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor>
void AddToExternalScalarFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                *lattice.get(iX, iY, iZ).getExternal(whichScalar) += externalScalar;
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT AddToExternalScalarFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void AddToExternalScalarFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
AddToExternalScalarFunctional3D<T, Descriptor>
    *AddToExternalScalarFunctional3D<T, Descriptor>::clone() const
{
    return new AddToExternalScalarFunctional3D<T, Descriptor>(*this);
}

/* ************* Class SetExternalVectorFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
SetExternalVectorFunctional3D<T, Descriptor>::SetExternalVectorFunctional3D(
    int vectorStartsAt_, Array<T, Descriptor<T>::d> &externalVector_) :
    vectorStartsAt(vectorStartsAt_), externalVector(externalVector_)
{
    PLB_ASSERT(vectorStartsAt + Descriptor<T>::d <= Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor>
void SetExternalVectorFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                externalVector.to_cArray(lattice.get(iX, iY, iZ).getExternal(vectorStartsAt));
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT SetExternalVectorFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void SetExternalVectorFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
SetExternalVectorFunctional3D<T, Descriptor> *SetExternalVectorFunctional3D<T, Descriptor>::clone()
    const
{
    return new SetExternalVectorFunctional3D<T, Descriptor>(*this);
}

/* ************* Class MaskedSetExternalVectorFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
MaskedSetExternalVectorFunctional3D<T, Descriptor>::MaskedSetExternalVectorFunctional3D(
    int flag_, int vectorStartsAt_, Array<T, Descriptor<T>::d> &externalVector_) :
    flag(flag_), vectorStartsAt(vectorStartsAt_), externalVector(externalVector_)
{
    PLB_ASSERT(vectorStartsAt + Descriptor<T>::d <= Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor>
void MaskedSetExternalVectorFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &mask)
{
    Dot3D offset = computeRelativeDisplacement(lattice, mask);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint oZ = iZ + offset.z;
                if (mask.get(oX, oY, oZ) == flag) {
                    externalVector.to_cArray(lattice.get(iX, iY, iZ).getExternal(vectorStartsAt));
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void MaskedSetExternalVectorFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Lattice.
    modified[1] = modif::nothing;          // Mask.
}

template <typename T, template <typename U> class Descriptor>
MaskedSetExternalVectorFunctional3D<T, Descriptor>
    *MaskedSetExternalVectorFunctional3D<T, Descriptor>::clone() const
{
    return new MaskedSetExternalVectorFunctional3D<T, Descriptor>(*this);
}

/* ************* Class SetGenericExternalVectorFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class Functional>
SetGenericExternalVectorFunctional3D<T, Descriptor, Functional>::
    SetGenericExternalVectorFunctional3D(int vectorBeginsAt_, Functional const &functional_) :
    vectorBeginsAt(vectorBeginsAt_), functional(functional_)
{
    PLB_ASSERT(vectorBeginsAt < Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor, class Functional>
void SetGenericExternalVectorFunctional3D<T, Descriptor, Functional>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    Dot3D absOffset = lattice.getLocation();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, Descriptor<T>::d> u;
                functional(iX + absOffset.x, iY + absOffset.y, iZ + absOffset.z, u);
                u.to_cArray(lattice.get(iX, iY, iZ).getExternal(vectorBeginsAt));
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, class Functional>
void SetGenericExternalVectorFunctional3D<T, Descriptor, Functional>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor, class Functional>
SetGenericExternalVectorFunctional3D<T, Descriptor, Functional>
    *SetGenericExternalVectorFunctional3D<T, Descriptor, Functional>::clone() const
{
    return new SetGenericExternalVectorFunctional3D<T, Descriptor, Functional>(*this);
}

/* ************* Class MaskedSetGenericExternalVectorFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class Functional>
MaskedSetGenericExternalVectorFunctional3D<T, Descriptor, Functional>::
    MaskedSetGenericExternalVectorFunctional3D(
        int flag_, int vectorBeginsAt_, Functional const &functional_) :
    flag(flag_), vectorBeginsAt(vectorBeginsAt_), functional(functional_)
{
    PLB_ASSERT(vectorBeginsAt < Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor, class Functional>
void MaskedSetGenericExternalVectorFunctional3D<T, Descriptor, Functional>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &mask)
{
    Dot3D offset = computeRelativeDisplacement(lattice, mask);
    Dot3D absOffset = lattice.getLocation();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint oZ = iZ + offset.z;

                if (mask.get(oX, oY, oZ) == flag) {
                    Array<T, Descriptor<T>::d> u;
                    functional(iX + absOffset.x, iY + absOffset.y, iZ + absOffset.z, u);
                    u.to_cArray(lattice.get(iX, iY, iZ).getExternal(vectorBeginsAt));
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, class Functional>
void MaskedSetGenericExternalVectorFunctional3D<T, Descriptor, Functional>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Lattice.
    modified[1] = modif::nothing;          // Mask.
}

template <typename T, template <typename U> class Descriptor, class Functional>
MaskedSetGenericExternalVectorFunctional3D<T, Descriptor, Functional>
    *MaskedSetGenericExternalVectorFunctional3D<T, Descriptor, Functional>::clone() const
{
    return new MaskedSetGenericExternalVectorFunctional3D<T, Descriptor, Functional>(*this);
}

/* ************* Class SetExternalVectorFromTensorFieldFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, int nDim>
SetExternalVectorFromTensorFieldFunctional3D<
    T, Descriptor, nDim>::SetExternalVectorFromTensorFieldFunctional3D(int vectorStartsAt_) :
    vectorStartsAt(vectorStartsAt_)
{
    PLB_ASSERT(vectorStartsAt + nDim <= Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor, int nDim>
void SetExternalVectorFromTensorFieldFunctional3D<T, Descriptor, nDim>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, TensorField3D<T, nDim> &tensor)
{
    Dot3D offset = computeRelativeDisplacement(lattice, tensor);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint oZ = iZ + offset.z;
                Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
                Array<T, nDim> externalVector = tensor.get(oX, oY, oZ);

                for (plint iD = 0; iD < nDim; ++iD) {
                    *cell.getExternal(vectorStartsAt + iD) = externalVector[iD];
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, int nDim>
BlockDomain::DomainT SetExternalVectorFromTensorFieldFunctional3D<T, Descriptor, nDim>::appliesTo()
    const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor, int nDim>
void SetExternalVectorFromTensorFieldFunctional3D<T, Descriptor, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, template <typename U> class Descriptor, int nDim>
SetExternalVectorFromTensorFieldFunctional3D<T, Descriptor, nDim>
    *SetExternalVectorFromTensorFieldFunctional3D<T, Descriptor, nDim>::clone() const
{
    return new SetExternalVectorFromTensorFieldFunctional3D<T, Descriptor, nDim>(*this);
}

/* ************* Class InterpolatePopulationsFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
InterpolatePopulationsFunctional3D<T, Descriptor>::InterpolatePopulationsFunctional3D(
    plint minIter_, plint maxIter_) :
    minIter(minIter_), maxIter(maxIter_)
{ }

template <typename T, template <typename U> class Descriptor>
void InterpolatePopulationsFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice1, BlockLattice3D<T, Descriptor> &lattice2)
{
    plint currentIter = (plint)lattice1.getTimeCounter().getTime();
    T lambda = (T)0.5 + 0.5 * (T)(currentIter - minIter) / (T)(maxIter - minIter);
    if (lambda < (T)0.5)
        lambda = (T)0.5;
    if (lambda > T(1.))
        lambda = T(1.);
    Dot3D offset = computeRelativeDisplacement(lattice1, lattice2);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint oZ = iZ + offset.z;
                Cell<T, Descriptor> &cell1 = lattice1.get(iX, iY, iZ);
                Cell<T, Descriptor> &cell2 = lattice2.get(oX, oY, oZ);
                for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                    cell1[iPop] = lambda * cell1[iPop] + (1. - lambda) * cell2[iPop];
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT InterpolatePopulationsFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void InterpolatePopulationsFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, template <typename U> class Descriptor>
InterpolatePopulationsFunctional3D<T, Descriptor>
    *InterpolatePopulationsFunctional3D<T, Descriptor>::clone() const
{
    return new InterpolatePopulationsFunctional3D<T, Descriptor>(*this);
}

/* *************** PART II ******************************************* */
/* *************** Initialization of scalar- and tensor-fields ******* */
/* ******************************************************************* */

/* ************** Class IniConstScalarFunctional3D ******************* */

template <typename T>
IniConstScalarFunctional3D<T>::IniConstScalarFunctional3D(T value_) : value(value_)
{ }

template <typename T>
void IniConstScalarFunctional3D<T>::process(Box3D domain, ScalarField3D<T> &field)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                field.get(iX, iY, iZ) = value;
            }
        }
    }
}

template <typename T>
IniConstScalarFunctional3D<T> *IniConstScalarFunctional3D<T>::clone() const
{
    return new IniConstScalarFunctional3D<T>(*this);
}

template <typename T>
BlockDomain::DomainT IniConstScalarFunctional3D<T>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T>
void IniConstScalarFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************** Class MaskedIniConstScalarFunctional3D ******************* */

template <typename T>
MaskedIniConstScalarFunctional3D<T>::MaskedIniConstScalarFunctional3D(int flag_, T value_) :
    flag(flag_), value(value_)
{ }

template <typename T>
void MaskedIniConstScalarFunctional3D<T>::process(
    Box3D domain, ScalarField3D<T> &field, ScalarField3D<int> &mask)
{
    Dot3D offset = computeRelativeDisplacement(field, mask);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (mask.get(iX + offset.x, iY + offset.y, iZ + offset.z) == flag) {
                    field.get(iX, iY, iZ) = value;
                }
            }
        }
    }
}

template <typename T>
MaskedIniConstScalarFunctional3D<T> *MaskedIniConstScalarFunctional3D<T>::clone() const
{
    return new MaskedIniConstScalarFunctional3D<T>(*this);
}

template <typename T>
BlockDomain::DomainT MaskedIniConstScalarFunctional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T>
void MaskedIniConstScalarFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Scalar-Field.
    modified[1] = modif::nothing;          // Mask.
}

/* ************** Class MaskedIniConstScalarFunctional3D_N ******************* */

template <typename T>
MaskedIniConstScalarFunctional3D_N<T>::MaskedIniConstScalarFunctional3D_N(int flag_, T value_) :
    flag(flag_), value(value_)
{ }

template <typename T>
void MaskedIniConstScalarFunctional3D_N<T>::process(
    Box3D domain, ScalarField3D<T> &field, NTensorField3D<int> &mask)
{
    Dot3D offset = computeRelativeDisplacement(field, mask);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (*mask.get(iX + offset.x, iY + offset.y, iZ + offset.z) == flag) {
                    field.get(iX, iY, iZ) = value;
                }
            }
        }
    }
}

template <typename T>
MaskedIniConstScalarFunctional3D_N<T> *MaskedIniConstScalarFunctional3D_N<T>::clone() const
{
    return new MaskedIniConstScalarFunctional3D_N<T>(*this);
}

template <typename T>
BlockDomain::DomainT MaskedIniConstScalarFunctional3D_N<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T>
void MaskedIniConstScalarFunctional3D_N<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Scalar-Field.
    modified[1] = modif::nothing;          // Mask.
}

/* ************** Class IniConstTensorFunctional3D ******************* */

template <typename T, int nDim>
IniConstTensorFunctional3D<T, nDim>::IniConstTensorFunctional3D(Array<T, nDim> const &value_) :
    value(value_)
{ }

template <typename T, int nDim>
void IniConstTensorFunctional3D<T, nDim>::process(Box3D domain, TensorField3D<T, nDim> &field)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                field.get(iX, iY, iZ) = value;
            }
        }
    }
}

template <typename T, int nDim>
IniConstTensorFunctional3D<T, nDim> *IniConstTensorFunctional3D<T, nDim>::clone() const
{
    return new IniConstTensorFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
BlockDomain::DomainT IniConstTensorFunctional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, int nDim>
void IniConstTensorFunctional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************** Class MaskedIniConstTensorFunctional3D ******************* */

template <typename T, int nDim>
MaskedIniConstTensorFunctional3D<T, nDim>::MaskedIniConstTensorFunctional3D(
    int flag_, Array<T, nDim> const &value_) :
    flag(flag_), value(value_)
{ }

template <typename T, int nDim>
void MaskedIniConstTensorFunctional3D<T, nDim>::process(
    Box3D domain, ScalarField3D<int> &mask, TensorField3D<T, nDim> &field)
{
    Dot3D offset = computeRelativeDisplacement(mask, field);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (mask.get(iX, iY, iZ) == flag) {
                    field.get(iX + offset.x, iY + offset.y, iZ + offset.z) = value;
                }
            }
        }
    }
}

template <typename T, int nDim>
MaskedIniConstTensorFunctional3D<T, nDim> *MaskedIniConstTensorFunctional3D<T, nDim>::clone() const
{
    return new MaskedIniConstTensorFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
BlockDomain::DomainT MaskedIniConstTensorFunctional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void MaskedIniConstTensorFunctional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Mask.
    modified[1] = modif::staticVariables;  // Tensor-Field.
}

/* ************** Class MaskedIniConstTensorFunctional3D_N ******************* */

template <typename T, int nDim>
MaskedIniConstTensorFunctional3D_N<T, nDim>::MaskedIniConstTensorFunctional3D_N(
    int flag_, Array<T, nDim> const &value_) :
    flag(flag_), value(value_)
{ }

template <typename T, int nDim>
void MaskedIniConstTensorFunctional3D_N<T, nDim>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    NTensorField3D<int> *mask = dynamic_cast<NTensorField3D<int> *>(atomicBlocks[0]);
    PLB_ASSERT(mask);
    TensorField3D<T, nDim> *field = dynamic_cast<TensorField3D<T, nDim> *>(atomicBlocks[1]);
    PLB_ASSERT(field);

    Dot3D offset = computeRelativeDisplacement(*mask, *field);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (*mask->get(iX, iY, iZ) == flag) {
                    field->get(iX + offset.x, iY + offset.y, iZ + offset.z) = value;
                }
            }
        }
    }
}

template <typename T, int nDim>
MaskedIniConstTensorFunctional3D_N<T, nDim> *MaskedIniConstTensorFunctional3D_N<T, nDim>::clone()
    const
{
    return new MaskedIniConstTensorFunctional3D_N<T, nDim>(*this);
}

template <typename T, int nDim>
BlockDomain::DomainT MaskedIniConstTensorFunctional3D_N<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void MaskedIniConstTensorFunctional3D_N<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Mask.
    modified[1] = modif::staticVariables;  // Tensor-Field.
}

/* ************** Class SwapValuesBulkAndEnvelope3D_N ******************* */

template <typename T>
void SwapValuesBulkAndEnvelope3D_N<T>::process(
    Box3D domain, NTensorField3D<T> &A, NTensorField3D<T> &B)
{
    Dot3D offset = computeRelativeDisplacement(A, B);
    plint ndim = A.getNdim();
    PLB_ASSERT(ndim == B.getNdim());
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                for (plint i = 0; i < ndim; ++i) {
                    std::swap(
                        A.get(iX, iY, iZ)[i],
                        B.get(iX + offset.x, iY + offset.y, iZ + offset.z)[i]);
                }
            }
        }
    }
}

template <typename T>
BlockDomain::DomainT SwapValuesBulkAndEnvelope3D_N<T>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T>
SwapValuesBulkAndEnvelope3D_N<T> *SwapValuesBulkAndEnvelope3D_N<T>::clone() const
{
    return new SwapValuesBulkAndEnvelope3D_N<T>(*this);
}

template <typename T>
void SwapValuesBulkAndEnvelope3D_N<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************** Class SetToCoordinateFunctional3D ****************** */

template <typename T>
SetToCoordinateFunctional3D<T>::SetToCoordinateFunctional3D(plint index_) : index(index_)
{
    PLB_ASSERT(index >= 0 && index <= 2);
}

template <typename T>
void SetToCoordinateFunctional3D<T>::process(Box3D domain, ScalarField3D<T> &field)
{
    Dot3D relativeOffset = field.getLocation();
    Array<plint, 3> ofs(relativeOffset.x, relativeOffset.y, relativeOffset.z);
    Array<plint, 3> pos;
    for (pos[0] = domain.x0; pos[0] <= domain.x1; ++pos[0]) {
        for (pos[1] = domain.y0; pos[1] <= domain.y1; ++pos[1]) {
            for (pos[2] = domain.z0; pos[2] <= domain.z1; ++pos[2]) {
                field.get(pos[0], pos[1], pos[2]) = (T)(pos[index] + ofs[index]);
            }
        }
    }
}

template <typename T>
SetToCoordinateFunctional3D<T> *SetToCoordinateFunctional3D<T>::clone() const
{
    return new SetToCoordinateFunctional3D<T>(*this);
}

template <typename T>
BlockDomain::DomainT SetToCoordinateFunctional3D<T>::appliesTo() const
{
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template <typename T>
void SetToCoordinateFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************** Class SetToCoordinatesFunctional3D ***************** */

template <typename T>
SetToCoordinatesFunctional3D<T>::SetToCoordinatesFunctional3D()
{ }

template <typename T>
void SetToCoordinatesFunctional3D<T>::process(Box3D domain, TensorField3D<T, 3> &field)
{
    Dot3D relativeOffset = field.getLocation();
    Array<plint, 3> ofs(relativeOffset.x, relativeOffset.y, relativeOffset.z);
    Array<plint, 3> pos;
    for (pos[0] = domain.x0; pos[0] <= domain.x1; ++pos[0]) {
        for (pos[1] = domain.y0; pos[1] <= domain.y1; ++pos[1]) {
            for (pos[2] = domain.z0; pos[2] <= domain.z1; ++pos[2]) {
                Array<T, 3> &cell = field.get(pos[0], pos[1], pos[2]);
                cell[0] = (T)(pos[0] + ofs[0]);
                cell[1] = (T)(pos[1] + ofs[1]);
                cell[2] = (T)(pos[2] + ofs[2]);
            }
        }
    }
}

template <typename T>
SetToCoordinatesFunctional3D<T> *SetToCoordinatesFunctional3D<T>::clone() const
{
    return new SetToCoordinatesFunctional3D<T>(*this);
}

template <typename T>
BlockDomain::DomainT SetToCoordinatesFunctional3D<T>::appliesTo() const
{
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template <typename T>
void SetToCoordinatesFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************** Class SetToRandomFunctional3D ***************** */

template <typename T>
SetToRandomFunctional3D<T>::SetToRandomFunctional3D(Box3D boundingBox, uint32_t const *seed_) :
    nY(boundingBox.getNy()), nZ(boundingBox.getNz()), seed(seed_)
{ }

template <typename T>
void SetToRandomFunctional3D<T>::process(Box3D domain, ScalarField3D<T> &field)
{
    Dot3D relativeOffset = field.getLocation();
    sitmo::prng_engine eng(*seed);
    plint rng_index = 0;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint globalX = nY * (iX + relativeOffset.x);
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint globalY = nZ * (iY + relativeOffset.y + globalX);
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint globalZ = iZ + relativeOffset.z + globalY;
                PLB_ASSERT(globalZ >= rng_index);
                if (globalZ > rng_index) {
                    eng.discard(globalZ - rng_index);
                    rng_index = globalZ;
                }
                field.get(iX, iY, iZ) = (T)eng() / ((T)sitmo::prng_engine::max() + 1.0);
                ++rng_index;
            }
        }
    }
}

template <typename T>
SetToRandomFunctional3D<T> *SetToRandomFunctional3D<T>::clone() const
{
    return new SetToRandomFunctional3D<T>(*this);
}

template <typename T>
void SetToRandomFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************** Class SetTensorComponentFunctional3D ***************** */

template <typename T, int nDim>
SetTensorComponentFunctional3D<T, nDim>::SetTensorComponentFunctional3D(int whichDim_) :
    whichDim(whichDim_)
{ }

template <typename T, int nDim>
void SetTensorComponentFunctional3D<T, nDim>::process(
    Box3D domain, ScalarField3D<T> &scalarField, TensorField3D<T, nDim> &tensorField)
{
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                tensorField.get(iX + offset.x, iY + offset.y, iZ + offset.z)[whichDim] =
                    scalarField.get(iX, iY, iZ);
            }
        }
    }
}

template <typename T, int nDim>
SetTensorComponentFunctional3D<T, nDim> *SetTensorComponentFunctional3D<T, nDim>::clone() const
{
    return new SetTensorComponentFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
BlockDomain::DomainT SetTensorComponentFunctional3D<T, nDim>::appliesTo() const
{
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void SetTensorComponentFunctional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

/* ************** Class PropagateInZdirection3D ***************** */

template <typename T>
void PropagateInZdirection3D<T>::process(Box3D domain, ScalarField3D<T> &field)
{
    if (field.getFlag())
        return;
    Dot3D relativeOffset = field.getLocation();
    plint iniZ = domain.z0;
    if (iniZ + relativeOffset.z == 0) {
        ++iniZ;
    }
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = iniZ; iZ <= domain.z1; ++iZ) {
                T nextData = field.get(iX, iY, iZ - 1);
                if (nextData < (T)(-0.5))
                    return;
                field.get(iX, iY, iZ) = nextData;
            }
        }
    }
    // Indicate that all data has been propagated.
    field.setFlag(true);
}

template <typename T>
PropagateInZdirection3D<T> *PropagateInZdirection3D<T>::clone() const
{
    return new PropagateInZdirection3D<T>(*this);
}

template <typename T>
void PropagateInZdirection3D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ******** GrowDomainFunctional3D ************************************* */

template <typename T>
GrowDomainFunctional3D<T>::GrowDomainFunctional3D(T flag_) : flag(flag_)
{ }

template <typename T>
void GrowDomainFunctional3D<T>::process(Box3D domain, ScalarField3D<T> &voxels)
{
    ScalarField3D<T> tmpVoxels(voxels);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                for (plint dx = -1; dx <= 1; ++dx)
                    for (plint dy = -1; dy <= 1; ++dy)
                        for (plint dz = -1; dz <= 1; ++dz)
                            if (!(dx == 0 && dy == 0 && dz == 0)) {
                                plint nextX = iX + dx;
                                plint nextY = iY + dy;
                                plint nextZ = iZ + dz;
                                if (tmpVoxels.get(nextX, nextY, nextZ) == flag) {
                                    voxels.get(iX, iY, iZ) = flag;
                                }
                            }
            }
        }
    }
}

template <typename T>
GrowDomainFunctional3D<T> *GrowDomainFunctional3D<T>::clone() const
{
    return new GrowDomainFunctional3D<T>(*this);
}

template <typename T>
void GrowDomainFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT GrowDomainFunctional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ************** Class ScalarNeumannToUnusedEnvelopes3D ******************* */

template <typename T>
void ScalarNeumannToUnusedEnvelopes3D<T>::process(Box3D domain, ScalarField3D<T> &field)
{
    Box3D edomain(domain.enlarge(1));
    Box3D bulk(field.getBoundingBox().enlarge(-1));
    for (plint iX = edomain.x0; iX <= edomain.x1; ++iX) {
        for (plint iY = edomain.y0; iY <= edomain.y1; ++iY) {
            for (plint iZ = edomain.z0; iZ <= edomain.z1; ++iZ) {
                if (!contained(Dot3D(iX, iY, iZ), bulk)) {
                    T average = T();
                    plint numcells = 0;
                    for (plint x = -1; x <= +1; ++x) {
                        for (plint y = -1; y <= +1; ++y) {
                            for (plint z = -1; z <= +1; ++z) {
                                if (contained(Dot3D(iX + x, iY + y, iZ + z), bulk)) {
                                    average += field.get(iX + x, iY + y, iZ + z);
                                    ++numcells;
                                }
                            }
                        }
                    }
                    PLB_ASSERT(numcells > 0);
                    average /= (T)numcells;

                    field.get(iX, iY, iZ) = average;
                }
            }
        }
    }
}

template <typename T>
ScalarNeumannToUnusedEnvelopes3D<T> *ScalarNeumannToUnusedEnvelopes3D<T>::clone() const
{
    return new ScalarNeumannToUnusedEnvelopes3D<T>(*this);
}

template <typename T>
BlockDomain::DomainT ScalarNeumannToUnusedEnvelopes3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T>
void ScalarNeumannToUnusedEnvelopes3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************** Class TensorNeumannToUnusedEnvelopes3D ******************* */

template <typename T, int nDim>
void TensorNeumannToUnusedEnvelopes3D<T, nDim>::process(Box3D domain, TensorField3D<T, nDim> &field)
{
    Box3D edomain(domain.enlarge(1));
    Box3D bulk(field.getBoundingBox().enlarge(-1));
    for (plint iX = edomain.x0; iX <= edomain.x1; ++iX) {
        for (plint iY = edomain.y0; iY <= edomain.y1; ++iY) {
            for (plint iZ = edomain.z0; iZ <= edomain.z1; ++iZ) {
                if (!contained(Dot3D(iX, iY, iZ), bulk)) {
                    Array<T, nDim> average(Array<T, nDim>::zero());
                    plint numcells = 0;
                    for (plint x = -1; x <= +1; ++x) {
                        for (plint y = -1; y <= +1; ++y) {
                            for (plint z = -1; z <= +1; ++z) {
                                if (contained(Dot3D(iX + x, iY + y, iZ + z), bulk)) {
                                    average += field.get(iX + x, iY + y, iZ + z);
                                    ++numcells;
                                }
                            }
                        }
                    }
                    PLB_ASSERT(numcells > 0);
                    average /= (T)numcells;

                    field.get(iX, iY, iZ) = average;
                }
            }
        }
    }
}

template <typename T, int nDim>
TensorNeumannToUnusedEnvelopes3D<T, nDim> *TensorNeumannToUnusedEnvelopes3D<T, nDim>::clone() const
{
    return new TensorNeumannToUnusedEnvelopes3D<T, nDim>(*this);
}

template <typename T, int nDim>
BlockDomain::DomainT TensorNeumannToUnusedEnvelopes3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void TensorNeumannToUnusedEnvelopes3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

}  // namespace plb

#endif  // DATA_INITIALIZER_FUNCTIONAL_3D_HH
