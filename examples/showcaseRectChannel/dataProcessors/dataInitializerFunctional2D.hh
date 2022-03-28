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
#ifndef DATA_INITIALIZER_FUNCTIONAL_2D_HH
#define DATA_INITIALIZER_FUNCTIONAL_2D_HH

#include <limits>

#include "atomicBlock/blockLattice2D.h"
#include "core/cell.h"
#include "dataProcessors/dataInitializerFunctional2D.h"
#include "latticeBoltzmann/externalFieldAccess.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "multiGrid/multiGridUtil.h"
#include "sitmo/prng_engine.hpp"

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Initialization of the block-lattice *************** */
/* ******************************************************************* */

template <typename T, template <class U> class Descriptor>
OneCellFunctional2D<T, Descriptor>::~OneCellFunctional2D()
{ }

template <typename T, template <class U> class Descriptor>
BlockDomain::DomainT OneCellFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <class U> class Descriptor>
void OneCellFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, template <class U> class Descriptor>
void OneCellFunctional2D<T, Descriptor>::setscale(int dxScale, int dtScale)
{ }

template <typename T, template <class U> class Descriptor>
OneCellIndexedFunctional2D<T, Descriptor>::~OneCellIndexedFunctional2D()
{ }

template <typename T, template <class U> class Descriptor>
BlockDomain::DomainT OneCellIndexedFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <class U> class Descriptor>
void OneCellIndexedFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, template <class U> class Descriptor>
void OneCellIndexedFunctional2D<T, Descriptor>::setscale(int dxScale, int dtScale)
{ }

template <typename T, template <class U> class Descriptor>
OneCellIndexedWithRandFunctional2D<T, Descriptor>::~OneCellIndexedWithRandFunctional2D()
{ }

template <typename T, template <class U> class Descriptor>
BlockDomain::DomainT OneCellIndexedWithRandFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <class U> class Descriptor>
void OneCellIndexedWithRandFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, template <class U> class Descriptor>
void OneCellIndexedWithRandFunctional2D<T, Descriptor>::setscale(int dxScale, int dtScale)
{ }

template <typename T, template <class U> class Descriptor>
GenericLatticeFunctional2D<T, Descriptor>::GenericLatticeFunctional2D(
    OneCellFunctional2D<T, Descriptor> *f_) :
    f(f_)
{ }

template <typename T, template <class U> class Descriptor>
GenericLatticeFunctional2D<T, Descriptor>::GenericLatticeFunctional2D(
    GenericLatticeFunctional2D<T, Descriptor> const &rhs) :
    f(rhs.f->clone())
{ }

template <typename T, template <class U> class Descriptor>
GenericLatticeFunctional2D<T, Descriptor>::~GenericLatticeFunctional2D()
{
    delete f;
}

template <typename T, template <class U> class Descriptor>
GenericLatticeFunctional2D<T, Descriptor> &GenericLatticeFunctional2D<T, Descriptor>::operator=(
    GenericLatticeFunctional2D<T, Descriptor> const &rhs)
{
    delete f;
    f = rhs.f->clone();
    return *this;
}

template <typename T, template <class U> class Descriptor>
void GenericLatticeFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            f->execute(lattice.get(iX, iY));
        }
    }
}

template <typename T, template <class U> class Descriptor>
GenericLatticeFunctional2D<T, Descriptor> *GenericLatticeFunctional2D<T, Descriptor>::clone() const
{
    return new GenericLatticeFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <class U> class Descriptor>
void GenericLatticeFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    f->getModificationPattern(modified);
}

template <typename T, template <class U> class Descriptor>
BlockDomain::DomainT GenericLatticeFunctional2D<T, Descriptor>::appliesTo() const
{
    return f->appliesTo();
}

template <typename T, template <class U> class Descriptor>
void GenericLatticeFunctional2D<T, Descriptor>::setscale(int dxScale, int dtScale)
{
    f->setscale(dxScale, dtScale);
}

template <typename T, template <class U> class Descriptor>
GenericIndexedLatticeFunctional2D<T, Descriptor>::GenericIndexedLatticeFunctional2D(
    OneCellIndexedFunctional2D<T, Descriptor> *f_) :
    f(f_)
{ }

template <typename T, template <class U> class Descriptor>
GenericIndexedLatticeFunctional2D<T, Descriptor>::GenericIndexedLatticeFunctional2D(
    GenericIndexedLatticeFunctional2D<T, Descriptor> const &rhs) :
    f(rhs.f->clone())
{ }

template <typename T, template <class U> class Descriptor>
GenericIndexedLatticeFunctional2D<T, Descriptor>::~GenericIndexedLatticeFunctional2D()
{
    delete f;
}

template <typename T, template <class U> class Descriptor>
GenericIndexedLatticeFunctional2D<T, Descriptor>
    &GenericIndexedLatticeFunctional2D<T, Descriptor>::operator=(
        GenericIndexedLatticeFunctional2D<T, Descriptor> const &rhs)
{
    delete f;
    f = rhs.f->clone();
    return *this;
}

template <typename T, template <class U> class Descriptor>
void GenericIndexedLatticeFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    Dot2D relativeOffset = lattice.getLocation();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            f->execute(iX + relativeOffset.x, iY + relativeOffset.y, lattice.get(iX, iY));
        }
    }
}

template <typename T, template <class U> class Descriptor>
GenericIndexedLatticeFunctional2D<T, Descriptor>
    *GenericIndexedLatticeFunctional2D<T, Descriptor>::clone() const
{
    return new GenericIndexedLatticeFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <class U> class Descriptor>
void GenericIndexedLatticeFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    f->getTypeOfModification(modified);
}

template <typename T, template <class U> class Descriptor>
BlockDomain::DomainT GenericIndexedLatticeFunctional2D<T, Descriptor>::appliesTo() const
{
    return f->appliesTo();
}

template <typename T, template <class U> class Descriptor>
void GenericIndexedLatticeFunctional2D<T, Descriptor>::setscale(int dxScale, int dtScale)
{
    f->setscale(dxScale, dtScale);
}

/* *************** Class GenericIndexedWithRandLatticeFunctional2D ************* */

template <typename T, template <class U> class Descriptor>
GenericIndexedWithRandLatticeFunctional2D<T, Descriptor>::GenericIndexedWithRandLatticeFunctional2D(
    OneCellIndexedWithRandFunctional2D<T, Descriptor> *f_, Box2D boundingBox_,
    uint32_t const *seed_) :
    f(f_), nY(boundingBox_.getNy()), seed(seed_)
{ }

template <typename T, template <class U> class Descriptor>
GenericIndexedWithRandLatticeFunctional2D<T, Descriptor>::GenericIndexedWithRandLatticeFunctional2D(
    GenericIndexedWithRandLatticeFunctional2D<T, Descriptor> const &rhs) :
    f(rhs.f->clone()), nY(rhs.nY), seed(rhs.seed)
{ }

template <typename T, template <class U> class Descriptor>
GenericIndexedWithRandLatticeFunctional2D<
    T, Descriptor>::~GenericIndexedWithRandLatticeFunctional2D()
{
    delete f;
}

template <typename T, template <class U> class Descriptor>
GenericIndexedWithRandLatticeFunctional2D<T, Descriptor>
    &GenericIndexedWithRandLatticeFunctional2D<T, Descriptor>::operator=(
        GenericIndexedWithRandLatticeFunctional2D<T, Descriptor> const &rhs)
{
    delete f;
    f = rhs.f->clone();
    return *this;
}

template <typename T, template <class U> class Descriptor>
void GenericIndexedWithRandLatticeFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    Dot2D relativeOffset = lattice.getLocation();
    sitmo::prng_engine eng(*seed);
    plint rng_index = 0;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint globalX = nY * (iX + relativeOffset.x);
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint globalY = iY + relativeOffset.y + globalX;
            PLB_ASSERT(globalY >= rng_index);
            if (globalY > rng_index) {
                eng.discard(globalY - rng_index);
                rng_index = globalY;
            }
            T rand_val = (T)eng() / ((T)sitmo::prng_engine::max() + 1.0);
            ++rng_index;
            f->execute(iX + relativeOffset.x, iY + relativeOffset.y, rand_val, lattice.get(iX, iY));
        }
    }
}

template <typename T, template <class U> class Descriptor>
GenericIndexedWithRandLatticeFunctional2D<T, Descriptor>
    *GenericIndexedWithRandLatticeFunctional2D<T, Descriptor>::clone() const
{
    return new GenericIndexedWithRandLatticeFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <class U> class Descriptor>
void GenericIndexedWithRandLatticeFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    f->getTypeOfModification(modified);
}

template <typename T, template <class U> class Descriptor>
BlockDomain::DomainT GenericIndexedWithRandLatticeFunctional2D<T, Descriptor>::appliesTo() const
{
    return f->appliesTo();
}

template <typename T, template <class U> class Descriptor>
void GenericIndexedWithRandLatticeFunctional2D<T, Descriptor>::setscale(int dxScale, int dtScale)
{
    f->setscale(dxScale, dtScale);
}

/* *************** Class InstantiateDynamicsFunctional2D ************* */

template <typename T, template <typename U> class Descriptor>
InstantiateDynamicsFunctional2D<T, Descriptor>::InstantiateDynamicsFunctional2D(
    Dynamics<T, Descriptor> *dynamics_) :
    dynamics(dynamics_)
{ }

template <typename T, template <typename U> class Descriptor>
InstantiateDynamicsFunctional2D<T, Descriptor>::InstantiateDynamicsFunctional2D(
    InstantiateDynamicsFunctional2D<T, Descriptor> const &rhs) :
    dynamics(rhs.dynamics->clone())
{ }

template <typename T, template <typename U> class Descriptor>
InstantiateDynamicsFunctional2D<T, Descriptor>
    &InstantiateDynamicsFunctional2D<T, Descriptor>::operator=(
        InstantiateDynamicsFunctional2D<T, Descriptor> const &rhs)
{
    delete dynamics;
    dynamics = rhs.dynamics->clone();
    return *this;
}

template <typename T, template <typename U> class Descriptor>
InstantiateDynamicsFunctional2D<T, Descriptor>::~InstantiateDynamicsFunctional2D()
{
    delete dynamics;
}

template <typename T, template <typename U> class Descriptor>
void InstantiateDynamicsFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            lattice.attributeDynamics(iX, iY, dynamics->clone());
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT InstantiateDynamicsFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void InstantiateDynamicsFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dataStructure;
}

template <typename T, template <typename U> class Descriptor>
InstantiateDynamicsFunctional2D<T, Descriptor>
    *InstantiateDynamicsFunctional2D<T, Descriptor>::clone() const
{
    return new InstantiateDynamicsFunctional2D<T, Descriptor>(*this);
}

/* ************* Class InstantiateComplexDomainDynamicsFunctional2D ** */

template <typename T, template <typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional2D<T, Descriptor>::
    InstantiateComplexDomainDynamicsFunctional2D(
        Dynamics<T, Descriptor> *dynamics_, DomainFunctional2D *domain_) :
    dynamics(dynamics_), domain(domain_)
{ }

template <typename T, template <typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional2D<T, Descriptor>::
    InstantiateComplexDomainDynamicsFunctional2D(
        InstantiateComplexDomainDynamicsFunctional2D<T, Descriptor> const &rhs) :
    dynamics(rhs.dynamics->clone()), domain(rhs.domain->clone())
{ }

template <typename T, template <typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional2D<T, Descriptor>
    &InstantiateComplexDomainDynamicsFunctional2D<T, Descriptor>::operator=(
        InstantiateComplexDomainDynamicsFunctional2D<T, Descriptor> const &rhs)
{
    delete dynamics;
    dynamics = rhs.dynamics->clone();
    delete domain;
    domain = rhs.domain->clone();
    return *this;
}

template <typename T, template <typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional2D<
    T, Descriptor>::~InstantiateComplexDomainDynamicsFunctional2D()
{
    delete dynamics;
    delete domain;
}

template <typename T, template <typename U> class Descriptor>
void InstantiateComplexDomainDynamicsFunctional2D<T, Descriptor>::process(
    Box2D boundingBox, BlockLattice2D<T, Descriptor> &lattice)
{
    Dot2D relativeOffset = lattice.getLocation();
    for (plint iX = boundingBox.x0; iX <= boundingBox.x1; ++iX) {
        for (plint iY = boundingBox.y0; iY <= boundingBox.y1; ++iY) {
            if ((*domain)(iX + relativeOffset.x, iY + relativeOffset.y)) {
                lattice.attributeDynamics(iX, iY, dynamics->clone());
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT InstantiateComplexDomainDynamicsFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void InstantiateComplexDomainDynamicsFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dataStructure;
}

template <typename T, template <typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional2D<T, Descriptor>
    *InstantiateComplexDomainDynamicsFunctional2D<T, Descriptor>::clone() const
{
    return new InstantiateComplexDomainDynamicsFunctional2D<T, Descriptor>(*this);
}

/* ************* Class InstantiateDotDynamicsFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
InstantiateDotDynamicsFunctional2D<T, Descriptor>::InstantiateDotDynamicsFunctional2D(
    Dynamics<T, Descriptor> *dynamics_) :
    dynamics(dynamics_)
{ }

template <typename T, template <typename U> class Descriptor>
InstantiateDotDynamicsFunctional2D<T, Descriptor>::InstantiateDotDynamicsFunctional2D(
    InstantiateDotDynamicsFunctional2D<T, Descriptor> const &rhs) :
    dynamics(rhs.dynamics->clone())
{ }

template <typename T, template <typename U> class Descriptor>
InstantiateDotDynamicsFunctional2D<T, Descriptor>
    &InstantiateDotDynamicsFunctional2D<T, Descriptor>::operator=(
        InstantiateDotDynamicsFunctional2D<T, Descriptor> const &rhs)
{
    delete dynamics;
    dynamics = rhs.dynamics->clone();
    return *this;
}

template <typename T, template <typename U> class Descriptor>
InstantiateDotDynamicsFunctional2D<T, Descriptor>::~InstantiateDotDynamicsFunctional2D()
{
    delete dynamics;
}

template <typename T, template <typename U> class Descriptor>
void InstantiateDotDynamicsFunctional2D<T, Descriptor>::process(
    DotList2D const &dotList, BlockLattice2D<T, Descriptor> &lattice)
{
    for (plint iDot = 0; iDot < dotList.getN(); ++iDot) {
        Dot2D const &dot = dotList.getDot(iDot);
        lattice.attributeDynamics(dot.x, dot.y, dynamics->clone());
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT InstantiateDotDynamicsFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void InstantiateDotDynamicsFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dataStructure;
}

template <typename T, template <typename U> class Descriptor>
InstantiateDotDynamicsFunctional2D<T, Descriptor>
    *InstantiateDotDynamicsFunctional2D<T, Descriptor>::clone() const
{
    return new InstantiateDotDynamicsFunctional2D<T, Descriptor>(*this);
}

/* ************* Class DynamicsFromMaskFunctional2D ************************ */

template <typename T, template <typename U> class Descriptor>
DynamicsFromMaskFunctional2D<T, Descriptor>::DynamicsFromMaskFunctional2D(
    Dynamics<T, Descriptor> *dynamics_, bool whichFlag_) :
    dynamics(dynamics_), whichFlag(whichFlag_)
{ }

template <typename T, template <typename U> class Descriptor>
DynamicsFromMaskFunctional2D<T, Descriptor>::DynamicsFromMaskFunctional2D(
    DynamicsFromMaskFunctional2D<T, Descriptor> const &rhs) :
    dynamics(rhs.dynamics->clone()), whichFlag(rhs.whichFlag)
{ }

template <typename T, template <typename U> class Descriptor>
DynamicsFromMaskFunctional2D<T, Descriptor> &DynamicsFromMaskFunctional2D<T, Descriptor>::operator=(
    DynamicsFromMaskFunctional2D<T, Descriptor> const &rhs)
{
    delete dynamics;
    dynamics = rhs.dynamics->clone();
    whichFlag = rhs.whichFlag;
    return *this;
}

template <typename T, template <typename U> class Descriptor>
DynamicsFromMaskFunctional2D<T, Descriptor>::~DynamicsFromMaskFunctional2D()
{
    delete dynamics;
}

template <typename T, template <typename U> class Descriptor>
void DynamicsFromMaskFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<bool> &mask)
{
    Dot2D offset = computeRelativeDisplacement(lattice, mask);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            bool flag = mask.get(iX + offset.x, iY + offset.y);
            if (util::boolIsEqual(flag, whichFlag)) {
                lattice.attributeDynamics(iX, iY, dynamics->clone());
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT DynamicsFromMaskFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void DynamicsFromMaskFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dataStructure;
    modified[1] = modif::nothing;
}

template <typename T, template <typename U> class Descriptor>
DynamicsFromMaskFunctional2D<T, Descriptor> *DynamicsFromMaskFunctional2D<T, Descriptor>::clone()
    const
{
    return new DynamicsFromMaskFunctional2D<T, Descriptor>(*this);
}

/* ************* Class DynamicsFromIntMaskFunctional2D ************************ */

template <typename T, template <typename U> class Descriptor>
DynamicsFromIntMaskFunctional2D<T, Descriptor>::DynamicsFromIntMaskFunctional2D(
    Dynamics<T, Descriptor> *dynamics_, int whichFlag_) :
    dynamics(dynamics_), whichFlag(whichFlag_)
{ }

template <typename T, template <typename U> class Descriptor>
DynamicsFromIntMaskFunctional2D<T, Descriptor>::DynamicsFromIntMaskFunctional2D(
    DynamicsFromIntMaskFunctional2D<T, Descriptor> const &rhs) :
    dynamics(rhs.dynamics->clone()), whichFlag(rhs.whichFlag)
{ }

template <typename T, template <typename U> class Descriptor>
DynamicsFromIntMaskFunctional2D<T, Descriptor>
    &DynamicsFromIntMaskFunctional2D<T, Descriptor>::operator=(
        DynamicsFromIntMaskFunctional2D<T, Descriptor> const &rhs)
{
    delete dynamics;
    dynamics = rhs.dynamics->clone();
    whichFlag = rhs.whichFlag;
    return *this;
}

template <typename T, template <typename U> class Descriptor>
DynamicsFromIntMaskFunctional2D<T, Descriptor>::~DynamicsFromIntMaskFunctional2D()
{
    delete dynamics;
}

template <typename T, template <typename U> class Descriptor>
void DynamicsFromIntMaskFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<int> &mask)
{
    Dot2D offset = computeRelativeDisplacement(lattice, mask);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            int flag = mask.get(iX + offset.x, iY + offset.y);
            if (flag == whichFlag) {
                lattice.attributeDynamics(iX, iY, dynamics->clone());
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT DynamicsFromIntMaskFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void DynamicsFromIntMaskFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dataStructure;
    modified[1] = modif::nothing;
}

template <typename T, template <typename U> class Descriptor>
DynamicsFromIntMaskFunctional2D<T, Descriptor>
    *DynamicsFromIntMaskFunctional2D<T, Descriptor>::clone() const
{
    return new DynamicsFromIntMaskFunctional2D<T, Descriptor>(*this);
}

/* ************* Class RecomposeFromFlowVariablesFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
void RecomposeFromFlowVariablesFunctional2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    BlockLattice2D<T, Descriptor> &lattice =
        dynamic_cast<BlockLattice2D<T, Descriptor> &>(*atomicBlocks[0]);
    ScalarField2D<T> const &rhoField = dynamic_cast<ScalarField2D<T> const &>(*atomicBlocks[1]);
    TensorField2D<T, 2> const &uField = dynamic_cast<TensorField2D<T, 2> const &>(*atomicBlocks[2]);
    TensorField2D<T, 3> const &SField = dynamic_cast<TensorField2D<T, 3> const &>(*atomicBlocks[3]);

    Dot2D offset1 = computeRelativeDisplacement(lattice, rhoField);
    Dot2D offset2 = computeRelativeDisplacement(lattice, uField);
    Dot2D offset3 = computeRelativeDisplacement(lattice, SField);

    std::vector<T> rawData(
        1 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n
        + Descriptor<T>::ExternalField::numScalars);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Cell<T, Descriptor> &cell = lattice.get(iX, iY);
            T rho = rhoField.get(iX + offset1.x, iY + offset1.y);
            Array<T, 2> const &u = uField.get(iX + offset2.x, iY + offset2.y);
            Array<T, 3> const &S = SField.get(iX + offset3.x, iY + offset3.y);

            // Convert rho --> rhoBar.
            rawData[0] = Descriptor<T>::rhoBar(rho);

            // Convert u --> j
            rawData[1] = rho * u[0];
            rawData[2] = rho * u[1];

            // Convert S --> PiNeq.
            T omega = cell.getDynamics().getOmega();
            T prefactor = -Descriptor<T>::cs2 * rho * (T)2 / omega;
            rawData[3] = S[0] * prefactor;
            rawData[4] = S[1] * prefactor;
            rawData[5] = S[2] * prefactor;

            // Recompose the cell.
            plint recomposeOrder = 1;
            cell.getDynamics().recompose(cell, rawData, recomposeOrder);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
RecomposeFromFlowVariablesFunctional2D<T, Descriptor>
    *RecomposeFromFlowVariablesFunctional2D<T, Descriptor>::clone() const
{
    return new RecomposeFromFlowVariablesFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT RecomposeFromFlowVariablesFunctional2D<T, Descriptor>::appliesTo() const
{
    // We could directly apply to the envelope too, but let's keep it
    //   bulk-only for future compatibility.
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void RecomposeFromFlowVariablesFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
    modified[3] = modif::nothing;
}

/* ************* Class AssignOmegaFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
AssignOmegaFunctional2D<T, Descriptor>::AssignOmegaFunctional2D(T omega_) : omega(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
void AssignOmegaFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    // Define dimensions of a viscosity.
    int dimDx = 2;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx, this->getDtScale(), dimDt);
    T nu_cs2 = (T)1 / omega - (T)1 / (T)2;
    T scaledOmega = (T)1 / (scaleFactor * nu_cs2 + (T)1 / (T)2);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            lattice.get(iX, iY).getDynamics().setOmega(scaledOmega);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
AssignOmegaFunctional2D<T, Descriptor> *AssignOmegaFunctional2D<T, Descriptor>::clone() const
{
    return new AssignOmegaFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT AssignOmegaFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void AssignOmegaFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;
}

/* ************* Class SetConstBoundaryVelocityFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
SetConstBoundaryVelocityFunctional2D<T, Descriptor>::SetConstBoundaryVelocityFunctional2D(
    Array<T, Descriptor<T>::d> velocity) :
    u(velocity)
{ }

template <typename T, template <typename U> class Descriptor>
void SetConstBoundaryVelocityFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    int dimDx = 1;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx, this->getDtScale(), dimDt);
    Array<T, 2> scaledU = u * scaleFactor;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            lattice.get(iX, iY).defineVelocity(scaledU);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
SetConstBoundaryVelocityFunctional2D<T, Descriptor>
    *SetConstBoundaryVelocityFunctional2D<T, Descriptor>::clone() const
{
    return new SetConstBoundaryVelocityFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT SetConstBoundaryVelocityFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void SetConstBoundaryVelocityFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::allVariables;
}

/* ************* Class SetConstBoundaryDensityFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
SetConstBoundaryDensityFunctional2D<T, Descriptor>::SetConstBoundaryDensityFunctional2D(T rho_) :
    rho(rho_)
{ }

template <typename T, template <typename U> class Descriptor>
void SetConstBoundaryDensityFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            lattice.get(iX, iY).defineDensity(rho);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
SetConstBoundaryDensityFunctional2D<T, Descriptor>
    *SetConstBoundaryDensityFunctional2D<T, Descriptor>::clone() const
{
    return new SetConstBoundaryDensityFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT SetConstBoundaryDensityFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void SetConstBoundaryDensityFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::allVariables;
}

/* ************* Class SetConstBoundaryTemperatureFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
SetConstBoundaryTemperatureFunctional2D<T, Descriptor>::SetConstBoundaryTemperatureFunctional2D(
    T temperature_) :
    temperature(temperature_)
{ }

template <typename T, template <typename U> class Descriptor>
void SetConstBoundaryTemperatureFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            lattice.get(iX, iY).defineTemperature(temperature);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
SetConstBoundaryTemperatureFunctional2D<T, Descriptor>
    *SetConstBoundaryTemperatureFunctional2D<T, Descriptor>::clone() const
{
    return new SetConstBoundaryTemperatureFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT SetConstBoundaryTemperatureFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void SetConstBoundaryTemperatureFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::allVariables;
}

/* ************* Class IniConstEquilibriumFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
IniConstEquilibriumFunctional2D<T, Descriptor>::IniConstEquilibriumFunctional2D(
    T density, Array<T, Descriptor<T>::d> velocity, T temperature) :
    rho(density), rhoBar(Descriptor<T>::rhoBar(density)), u(velocity), thetaBar(temperature - (T)1)
{ }

template <typename T, template <typename U> class Descriptor>
void IniConstEquilibriumFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    int dimDx = 1;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx, this->getDtScale(), dimDt);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Array<T, Descriptor<T>::d> f;
            f[0] = getExternalForceComponent(lattice.get(iX, iY), 0);
            f[1] = getExternalForceComponent(lattice.get(iX, iY), 1);
            Array<T, Descriptor<T>::d> scaledJ;
            scaledJ[0] = scaleFactor * rho * (u[0] - (T)0.5 * f[0]);
            scaledJ[1] = scaleFactor * rho * (u[1] - (T)0.5 * f[1]);
            T scaledJsqr = normSqr(scaledJ);
            for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                lattice.get(iX, iY)[iPop] = lattice.get(iX, iY).computeEquilibrium(
                    iPop, rhoBar, scaledJ, scaledJsqr, thetaBar);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
IniConstEquilibriumFunctional2D<T, Descriptor>
    *IniConstEquilibriumFunctional2D<T, Descriptor>::clone() const
{
    return new IniConstEquilibriumFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT IniConstEquilibriumFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void IniConstEquilibriumFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************* Class IniConstEquilibriumComplexDomainFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
IniConstEquilibriumComplexDomainFunctional2D<T, Descriptor>::
    IniConstEquilibriumComplexDomainFunctional2D(
        DomainFunctional2D *domain_, T density, Array<T, Descriptor<T>::d> velocity,
        T temperature) :
    domain(domain_),
    rho(density),
    rhoBar(Descriptor<T>::rhoBar(density)),
    u(velocity),
    thetaBar(temperature - (T)1)
{ }

template <typename T, template <typename U> class Descriptor>
void IniConstEquilibriumComplexDomainFunctional2D<T, Descriptor>::process(
    Box2D box, BlockLattice2D<T, Descriptor> &lattice)
{
    int dimDx = 1;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx, this->getDtScale(), dimDt);
    Dot2D relativeOffset = lattice.getLocation();
    for (plint iX = box.x0; iX <= box.x1; ++iX) {
        for (plint iY = box.y0; iY <= box.y1; ++iY) {
            if ((*domain)(iX + relativeOffset.x, iY + relativeOffset.y)) {
                Array<T, Descriptor<T>::d> f;
                f[0] = getExternalForceComponent(lattice.get(iX, iY), 0);
                f[1] = getExternalForceComponent(lattice.get(iX, iY), 1);
                Array<T, Descriptor<T>::d> scaledJ;
                scaledJ[0] = scaleFactor * rho * (u[0] - (T)0.5 * f[0]);
                scaledJ[1] = scaleFactor * rho * (u[1] - (T)0.5 * f[1]);
                T scaledJsqr = normSqr(scaledJ);
                for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                    lattice.get(iX, iY)[iPop] = lattice.get(iX, iY).computeEquilibrium(
                        iPop, rhoBar, scaledJ, scaledJsqr, thetaBar);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
IniConstEquilibriumComplexDomainFunctional2D<T, Descriptor>
    *IniConstEquilibriumComplexDomainFunctional2D<T, Descriptor>::clone() const
{
    return new IniConstEquilibriumComplexDomainFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT IniConstEquilibriumComplexDomainFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void IniConstEquilibriumComplexDomainFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************* Class StripeOffDensityOffsetFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
StripeOffDensityOffsetFunctional2D<T, Descriptor>::StripeOffDensityOffsetFunctional2D(T deltaRho_) :
    deltaRho(deltaRho_)
{ }

template <typename T, template <typename U> class Descriptor>
void StripeOffDensityOffsetFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Cell<T, Descriptor> &cell = lattice.get(iX, iY);
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

template <typename T, template <typename U> class Descriptor>
StripeOffDensityOffsetFunctional2D<T, Descriptor>
    *StripeOffDensityOffsetFunctional2D<T, Descriptor>::clone() const
{
    return new StripeOffDensityOffsetFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT StripeOffDensityOffsetFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void StripeOffDensityOffsetFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************* Class InstantiateCompositeDynamicsFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional2D<T, Descriptor>::InstantiateCompositeDynamicsFunctional2D(
    CompositeDynamics<T, Descriptor> *compositeDynamics_) :
    compositeDynamics(compositeDynamics_)
{ }

template <typename T, template <typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional2D<T, Descriptor>::InstantiateCompositeDynamicsFunctional2D(
    InstantiateCompositeDynamicsFunctional2D<T, Descriptor> const &rhs) :
    compositeDynamics(rhs.compositeDynamics->clone())
{ }

template <typename T, template <typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional2D<T, Descriptor>
    &InstantiateCompositeDynamicsFunctional2D<T, Descriptor>::operator=(
        InstantiateCompositeDynamicsFunctional2D<T, Descriptor> const &rhs)
{
    delete compositeDynamics;
    compositeDynamics = rhs.compositeDynamics->clone();
    return *this;
}

template <typename T, template <typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional2D<T, Descriptor>::~InstantiateCompositeDynamicsFunctional2D()
{
    delete compositeDynamics;
}

template <typename T, template <typename U> class Descriptor>
void InstantiateCompositeDynamicsFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            lattice.attributeDynamics(
                iX, iY,
                cloneAndInsertAtTopDynamics(
                    lattice.get(iX, iY).getDynamics(), compositeDynamics->clone()));
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT InstantiateCompositeDynamicsFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void InstantiateCompositeDynamicsFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dataStructure;
}

template <typename T, template <typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional2D<T, Descriptor>
    *InstantiateCompositeDynamicsFunctional2D<T, Descriptor>::clone() const
{
    return new InstantiateCompositeDynamicsFunctional2D<T, Descriptor>(*this);
}

/* ************* Class SetExternalScalarFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
SetExternalScalarFunctional2D<T, Descriptor>::SetExternalScalarFunctional2D(
    int whichScalar_, T externalScalar_) :
    whichScalar(whichScalar_), externalScalar(externalScalar_)
{
    PLB_ASSERT(whichScalar < Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor>
void SetExternalScalarFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            *lattice.get(iX, iY).getExternal(whichScalar) = externalScalar;
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT SetExternalScalarFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void SetExternalScalarFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
SetExternalScalarFunctional2D<T, Descriptor> *SetExternalScalarFunctional2D<T, Descriptor>::clone()
    const
{
    return new SetExternalScalarFunctional2D<T, Descriptor>(*this);
}

/* ************* Class SetGenericExternalScalarFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor, class Functional>
SetGenericExternalScalarFunctional2D<T, Descriptor, Functional>::
    SetGenericExternalScalarFunctional2D(int whichScalar_, Functional const &functional_) :
    whichScalar(whichScalar_), functional(functional_)
{
    PLB_ASSERT(whichScalar < Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor, class Functional>
void SetGenericExternalScalarFunctional2D<T, Descriptor, Functional>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    Dot2D absOffset = lattice.getLocation();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            *lattice.get(iX, iY).getExternal(whichScalar) =
                functional(iX + absOffset.x, iY + absOffset.y);
        }
    }
}

template <typename T, template <typename U> class Descriptor, class Functional>
void SetGenericExternalScalarFunctional2D<T, Descriptor, Functional>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor, class Functional>
SetGenericExternalScalarFunctional2D<T, Descriptor, Functional>
    *SetGenericExternalScalarFunctional2D<T, Descriptor, Functional>::clone() const
{
    return new SetGenericExternalScalarFunctional2D<T, Descriptor, Functional>(*this);
}

/* ************* Class AddToExternalScalarFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
AddToExternalScalarFunctional2D<T, Descriptor>::AddToExternalScalarFunctional2D(
    int whichScalar_, T externalScalar_) :
    whichScalar(whichScalar_), externalScalar(externalScalar_)
{
    PLB_ASSERT(whichScalar < Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor>
void AddToExternalScalarFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            *lattice.get(iX, iY).getExternal(whichScalar) += externalScalar;
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT AddToExternalScalarFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void AddToExternalScalarFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
AddToExternalScalarFunctional2D<T, Descriptor>
    *AddToExternalScalarFunctional2D<T, Descriptor>::clone() const
{
    return new AddToExternalScalarFunctional2D<T, Descriptor>(*this);
}

/* ************* Class SetExternalScalarFromScalarFieldFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
SetExternalScalarFromScalarFieldFunctional2D<
    T, Descriptor>::SetExternalScalarFromScalarFieldFunctional2D(int whichScalar_) :
    whichScalar(whichScalar_)
{
    PLB_ASSERT(whichScalar < Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor>
void SetExternalScalarFromScalarFieldFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &field)
{
    Dot2D offset = computeRelativeDisplacement(lattice, field);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = offset.x + iX;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY = offset.y + iY;
            *lattice.get(iX, iY).getExternal(whichScalar) = field.get(oX, oY);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT SetExternalScalarFromScalarFieldFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void SetExternalScalarFromScalarFieldFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, template <typename U> class Descriptor>
SetExternalScalarFromScalarFieldFunctional2D<T, Descriptor>
    *SetExternalScalarFromScalarFieldFunctional2D<T, Descriptor>::clone() const
{
    return new SetExternalScalarFromScalarFieldFunctional2D<T, Descriptor>(*this);
}

/* ************* Class SetExternalVectorFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
SetExternalVectorFunctional2D<T, Descriptor>::SetExternalVectorFunctional2D(
    int vectorStartsAt_, Array<T, Descriptor<T>::d> const &externalVector_) :
    vectorStartsAt(vectorStartsAt_), externalVector(externalVector_)
{
    PLB_ASSERT(vectorStartsAt + Descriptor<T>::d <= Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor>
void SetExternalVectorFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Cell<T, Descriptor> &cell = lattice.get(iX, iY);
            for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
                *cell.getExternal(vectorStartsAt + iD) = externalVector[iD];
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT SetExternalVectorFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void SetExternalVectorFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
SetExternalVectorFunctional2D<T, Descriptor> *SetExternalVectorFunctional2D<T, Descriptor>::clone()
    const
{
    return new SetExternalVectorFunctional2D<T, Descriptor>(*this);
}

/* ************* Class SetCustomExternalVectorFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor, class VectorFunction>
SetCustomExternalVectorFunctional2D<T, Descriptor, VectorFunction>::
    SetCustomExternalVectorFunctional2D(int vectorStartsAt_, VectorFunction f_) :
    vectorStartsAt(vectorStartsAt_), f(f_)
{
    PLB_ASSERT(vectorStartsAt + Descriptor<T>::d <= Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor, class VectorFunction>
void SetCustomExternalVectorFunctional2D<T, Descriptor, VectorFunction>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    Dot2D offset = lattice.getLocation();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY = iY + offset.y;
            Cell<T, Descriptor> &cell = lattice.get(iX, iY);
            for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
                Array<T, Descriptor<T>::d> externalVector;
                f(oX, oY, externalVector);
                *cell.getExternal(vectorStartsAt + iD) = externalVector[iD];
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, class VectorFunction>
BlockDomain::DomainT SetCustomExternalVectorFunctional2D<T, Descriptor, VectorFunction>::appliesTo()
    const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor, class VectorFunction>
void SetCustomExternalVectorFunctional2D<T, Descriptor, VectorFunction>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor, class VectorFunction>
SetCustomExternalVectorFunctional2D<T, Descriptor, VectorFunction>
    *SetCustomExternalVectorFunctional2D<T, Descriptor, VectorFunction>::clone() const
{
    return new SetCustomExternalVectorFunctional2D<T, Descriptor, VectorFunction>(*this);
}

/* ************* Class SetExternalVectorFromTensorFieldFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor, int nDim>
SetExternalVectorFromTensorFieldFunctional2D<
    T, Descriptor, nDim>::SetExternalVectorFromTensorFieldFunctional2D(int vectorStartsAt_) :
    vectorStartsAt(vectorStartsAt_)
{
    PLB_ASSERT(vectorStartsAt + nDim <= Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor, int nDim>
void SetExternalVectorFromTensorFieldFunctional2D<T, Descriptor, nDim>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, TensorField2D<T, nDim> &tensor)
{
    Dot2D offset = computeRelativeDisplacement(lattice, tensor);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY = iY + offset.y;
            Cell<T, Descriptor> &cell = lattice.get(iX, iY);
            Array<T, nDim> externalVector = tensor.get(oX, oY);

            for (plint iD = 0; iD < nDim; ++iD) {
                *cell.getExternal(vectorStartsAt + iD) = externalVector[iD];
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, int nDim>
BlockDomain::DomainT SetExternalVectorFromTensorFieldFunctional2D<T, Descriptor, nDim>::appliesTo()
    const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor, int nDim>
void SetExternalVectorFromTensorFieldFunctional2D<T, Descriptor, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[0] = modif::nothing;
}

template <typename T, template <typename U> class Descriptor, int nDim>
SetExternalVectorFromTensorFieldFunctional2D<T, Descriptor, nDim>
    *SetExternalVectorFromTensorFieldFunctional2D<T, Descriptor, nDim>::clone() const
{
    return new SetExternalVectorFromTensorFieldFunctional2D<T, Descriptor, nDim>(*this);
}

/* *************** PART II ******************************************* */
/* *************** Initialization of scalar- and tensor-fields ******* */
/* ******************************************************************* */

/* ************** Class IniConstScalarFunctional2D ******************* */

template <typename T>
IniConstScalarFunctional2D<T>::IniConstScalarFunctional2D(T value_) : value(value_)
{ }

template <typename T>
void IniConstScalarFunctional2D<T>::process(Box2D domain, ScalarField2D<T> &field)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            field.get(iX, iY) = value;
        }
    }
}

template <typename T>
IniConstScalarFunctional2D<T> *IniConstScalarFunctional2D<T>::clone() const
{
    return new IniConstScalarFunctional2D<T>(*this);
}

template <typename T>
BlockDomain::DomainT IniConstScalarFunctional2D<T>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T>
void IniConstScalarFunctional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************** Class MaskedIniConstScalarFunctional2D ******************* */

template <typename T>
MaskedIniConstScalarFunctional2D<T>::MaskedIniConstScalarFunctional2D(int flag_, T value_) :
    flag(flag_), value(value_)
{ }

template <typename T>
void MaskedIniConstScalarFunctional2D<T>::process(
    Box2D domain, ScalarField2D<T> &field, ScalarField2D<int> &mask)
{
    Dot2D offset = computeRelativeDisplacement(field, mask);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            if (mask.get(iX + offset.x, iY + offset.y) == flag) {
                field.get(iX, iY) = value;
            }
        }
    }
}

template <typename T>
MaskedIniConstScalarFunctional2D<T> *MaskedIniConstScalarFunctional2D<T>::clone() const
{
    return new MaskedIniConstScalarFunctional2D<T>(*this);
}

template <typename T>
BlockDomain::DomainT MaskedIniConstScalarFunctional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T>
void MaskedIniConstScalarFunctional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Scalar-Field.
    modified[1] = modif::nothing;          // Mask.
}

/* ************** Class IniConstTensorFunctional2D ******************* */

template <typename T, int nDim>
IniConstTensorFunctional2D<T, nDim>::IniConstTensorFunctional2D(Array<T, nDim> const &value_) :
    value(value_)
{ }

template <typename T, int nDim>
void IniConstTensorFunctional2D<T, nDim>::process(Box2D domain, TensorField2D<T, nDim> &field)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            field.get(iX, iY) = value;
        }
    }
}

template <typename T, int nDim>
IniConstTensorFunctional2D<T, nDim> *IniConstTensorFunctional2D<T, nDim>::clone() const
{
    return new IniConstTensorFunctional2D<T, nDim>(*this);
}

template <typename T, int nDim>
BlockDomain::DomainT IniConstTensorFunctional2D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, int nDim>
void IniConstTensorFunctional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************** Class MaskedIniConstTensorFunctional2D ******************* */

template <typename T, int nDim>
MaskedIniConstTensorFunctional2D<T, nDim>::MaskedIniConstTensorFunctional2D(
    int flag_, Array<T, nDim> const &value_) :
    flag(flag_), value(value_)
{ }

template <typename T, int nDim>
void MaskedIniConstTensorFunctional2D<T, nDim>::process(
    Box2D domain, ScalarField2D<int> &mask, TensorField2D<T, nDim> &field)
{
    Dot2D offset = computeRelativeDisplacement(mask, field);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            if (mask.get(iX, iY) == flag) {
                field.get(iX + offset.x, iY + offset.y) = value;
            }
        }
    }
}

template <typename T, int nDim>
MaskedIniConstTensorFunctional2D<T, nDim> *MaskedIniConstTensorFunctional2D<T, nDim>::clone() const
{
    return new MaskedIniConstTensorFunctional2D<T, nDim>(*this);
}

template <typename T, int nDim>
BlockDomain::DomainT MaskedIniConstTensorFunctional2D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void MaskedIniConstTensorFunctional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Mask.
    modified[1] = modif::staticVariables;  // Tensor-Field.
}

/* ************** Class SetToCoordinateFunctional2D ****************** */

template <typename T>
SetToCoordinateFunctional2D<T>::SetToCoordinateFunctional2D(plint index_) : index(index_)
{
    PLB_ASSERT(index >= 0 && index <= 1);
}

template <typename T>
void SetToCoordinateFunctional2D<T>::process(Box2D domain, ScalarField2D<T> &field)
{
    Dot2D relativeOffset = field.getLocation();
    Array<plint, 2> ofs(relativeOffset.x, relativeOffset.y);
    Array<plint, 2> pos;
    for (pos[0] = domain.x0; pos[0] <= domain.x1; ++pos[0]) {
        for (pos[1] = domain.y0; pos[1] <= domain.y1; ++pos[1]) {
            field.get(pos[0], pos[1]) = (T)(pos[index] + ofs[index]);
        }
    }
}

template <typename T>
SetToCoordinateFunctional2D<T> *SetToCoordinateFunctional2D<T>::clone() const
{
    return new SetToCoordinateFunctional2D<T>(*this);
}

template <typename T>
BlockDomain::DomainT SetToCoordinateFunctional2D<T>::appliesTo() const
{
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template <typename T>
void SetToCoordinateFunctional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************** Class SetToCoordinatesFunctional2D ***************** */

template <typename T>
SetToCoordinatesFunctional2D<T>::SetToCoordinatesFunctional2D()
{ }

template <typename T>
void SetToCoordinatesFunctional2D<T>::process(Box2D domain, TensorField2D<T, 2> &field)
{
    Dot2D relativeOffset = field.getLocation();
    Array<plint, 2> ofs(relativeOffset.x, relativeOffset.y);
    Array<plint, 2> pos;
    for (pos[0] = domain.x0; pos[0] <= domain.x1; ++pos[0]) {
        for (pos[1] = domain.y0; pos[1] <= domain.y1; ++pos[1]) {
            Array<T, 2> &cell = field.get(pos[0], pos[1]);
            cell[0] = (T)(pos[0] + ofs[0]);
            cell[1] = (T)(pos[1] + ofs[1]);
        }
    }
}

template <typename T>
SetToCoordinatesFunctional2D<T> *SetToCoordinatesFunctional2D<T>::clone() const
{
    return new SetToCoordinatesFunctional2D<T>(*this);
}

template <typename T>
BlockDomain::DomainT SetToCoordinatesFunctional2D<T>::appliesTo() const
{
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template <typename T>
void SetToCoordinatesFunctional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************** Class SetToRandomFunctional2D ***************** */

template <typename T>
SetToRandomFunctional2D<T>::SetToRandomFunctional2D(
    Box2D const &boundingBox_, uint32_t const *seed_) :
    nY(boundingBox_.getNy()), seed(seed_)
{ }

template <typename T>
void SetToRandomFunctional2D<T>::process(Box2D domain, ScalarField2D<T> &field)
{
    Dot2D relativeOffset = field.getLocation();
    sitmo::prng_engine eng(*seed);
    plint rng_index = 0;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint globalX = nY * (iX + relativeOffset.x);
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint globalY = iY + relativeOffset.y + globalX;
            PLB_ASSERT(globalY >= rng_index);
            if (globalY > rng_index) {
                eng.discard(globalY - rng_index);
                rng_index = globalY;
            }
            field.get(iX, iY) = (T)eng() / ((T)sitmo::prng_engine::max() + 1.0);
            ++rng_index;
        }
    }
}

template <typename T>
SetToRandomFunctional2D<T> *SetToRandomFunctional2D<T>::clone() const
{
    return new SetToRandomFunctional2D<T>(*this);
}

template <typename T>
void SetToRandomFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************** Class SetTensorComponentFunctional2D ***************** */

template <typename T, int nDim>
SetTensorComponentFunctional2D<T, nDim>::SetTensorComponentFunctional2D(int whichDim_) :
    whichDim(whichDim_)
{ }

template <typename T, int nDim>
void SetTensorComponentFunctional2D<T, nDim>::process(
    Box2D domain, ScalarField2D<T> &scalarField, TensorField2D<T, nDim> &tensorField)
{
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            tensorField.get(iX + offset.x, iY + offset.y)[whichDim] = scalarField.get(iX, iY);
        }
    }
}

template <typename T, int nDim>
SetTensorComponentFunctional2D<T, nDim> *SetTensorComponentFunctional2D<T, nDim>::clone() const
{
    return new SetTensorComponentFunctional2D<T, nDim>(*this);
}

template <typename T, int nDim>
BlockDomain::DomainT SetTensorComponentFunctional2D<T, nDim>::appliesTo() const
{
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void SetTensorComponentFunctional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

/* ******** GrowDomainFunctional2D ************************************* */

template <typename T>
GrowDomainFunctional2D<T>::GrowDomainFunctional2D(T flag_) : flag(flag_)
{ }

template <typename T>
void GrowDomainFunctional2D<T>::process(Box2D domain, ScalarField2D<T> &voxels)
{
    ScalarField2D<T> tmpVoxels(voxels);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint dx = -1; dx <= 1; ++dx)
                for (plint dy = -1; dy <= 1; ++dy)
                    if (!(dx == 0 && dy == 0)) {
                        plint nextX = iX + dx;
                        plint nextY = iY + dy;
                        if (tmpVoxels.get(nextX, nextY) == flag) {
                            voxels.get(iX, iY) = flag;
                        }
                    }
        }
    }
}

template <typename T>
GrowDomainFunctional2D<T> *GrowDomainFunctional2D<T>::clone() const
{
    return new GrowDomainFunctional2D<T>(*this);
}

template <typename T>
void GrowDomainFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT GrowDomainFunctional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

}  // namespace plb

#endif  // DATA_INITIALIZER_FUNCTIONAL_2D_HH
