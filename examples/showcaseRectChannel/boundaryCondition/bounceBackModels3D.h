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
 * BounceBack dynamics models in 3D -- header file.
 */

#ifndef BOUNCE_BACK_MODELS_3D_H
#define BOUNCE_BACK_MODELS_3D_H

#include "atomicBlock/reductiveDataProcessorWrapper3D.h"
#include "boundaryCondition/bounceBackModels.h"
#include "core/vectorFunction3D.h"
#include "dataProcessors/dataInitializerFunctional3D.h"

namespace plb {

template <typename T, template <class U> class Descriptor>
void initializeMomentumExchange(BlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <class U> class Descriptor>
void initializeMomentumExchange(
    BlockLattice3D<T, Descriptor> &lattice, Box3D boundingBox, DomainFunctional3D *domain);

template <typename T, template <class U> class Descriptor>
void initializeMomentumExchange(BlockLattice3D<T, Descriptor> &lattice, DotList3D const &dotList);

template <typename T, template <class U> class Descriptor>
void initializeMomentumExchange(MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <class U> class Descriptor>
void initializeMomentumExchange(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D boundingBox, DomainFunctional3D *domain);

template <typename T, template <class U> class Descriptor>
void initializeMomentumExchange(
    MultiBlockLattice3D<T, Descriptor> &lattice, DotList3D const &dotList);

template <typename T, template <typename U> class Descriptor>
class CountBBNeighborsFunctional3D : public BoxProcessingFunctional3D_LS<T, Descriptor, int> {
public:
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &neighbors);
    virtual CountBBNeighborsFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
class CountBBNeighbors_NTensor3D : public BoxProcessingFunctional3D_LN<T, Descriptor, int> {
public:
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, NTensorField3D<int> &neighbors);
    virtual CountBBNeighbors_NTensor3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
class MaskedCountBBNeighbors_NTensor3D :
    public MaskedBoxProcessingFunctional3D_LN<T, Descriptor, int> {
public:
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, NTensorField3D<int> &neighbors,
        NTensorField3D<int> &mask);
    virtual MaskedCountBBNeighbors_NTensor3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
void countBBNeighbors(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiNTensorField3D<int> &neighbors, Box3D domain);

template <typename T, template <typename U> class Descriptor>
void maskedCountBBNeighbors(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiNTensorField3D<int> &neighbors,
    MultiNTensorField3D<int> &mask, Box3D domain);

template <typename T, template <typename U> class Descriptor>
MultiNTensorField3D<int> *countBBNeighbors(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
MultiNTensorField3D<int> *maskedCountBBNeighbors(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiNTensorField3D<int> &mask, Box3D domain);

template <typename T, template <typename U> class Descriptor>
class SetAverageWallDensityOnVelocityBounceBack3D :
    public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual SetAverageWallDensityOnVelocityBounceBack3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
class SetWallVelocityOnVelocityBounceBack3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    SetWallVelocityOnVelocityBounceBack3D(VectorFunction3D<T> const &velocityFunction_);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual SetWallVelocityOnVelocityBounceBack3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    VectorFunction3D<T> const &velocityFunction;
};

/* *************** Class ComputeMomentumExchangeFunctional3D ************* */

template <typename T, template <typename U> class Descriptor>
class ComputeMomentumExchangeFunctional3D :
    public ReductiveBoxProcessingFunctional3D_L<T, Descriptor> {
public:
    ComputeMomentumExchangeFunctional3D();
    ComputeMomentumExchangeFunctional3D(Array<T, 3> const &center_);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual ComputeMomentumExchangeFunctional3D<T, Descriptor> *clone() const;
    Array<T, 3> getForce() const;
    Array<T, 3> getTorque() const;

private:
    bool computeTorque;
    Array<T, 3> center;
    Array<plint, 3> forceIds;
    Array<plint, 3> torqueIds;
};

template <typename T, template <typename U> class Descriptor>
class MaskedComputeMomentumExchangeFunctional3D :
    public ReductiveBoxProcessingFunctional3D_LS<T, Descriptor, int> {
public:
    MaskedComputeMomentumExchangeFunctional3D();
    MaskedComputeMomentumExchangeFunctional3D(Array<T, 3> const &center_);
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &mask);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual MaskedComputeMomentumExchangeFunctional3D<T, Descriptor> *clone() const;
    Array<T, 3> getForce() const;
    Array<T, 3> getTorque() const;

private:
    bool computeTorque;
    Array<T, 3> center;
    Array<plint, 3> forceIds;
    Array<plint, 3> torqueIds;
};

}  // namespace plb

#endif  // BOUNCE_BACK_MODELS_3D_H
