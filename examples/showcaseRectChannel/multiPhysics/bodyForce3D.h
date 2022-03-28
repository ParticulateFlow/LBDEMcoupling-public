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

#ifndef BODY_FORCE_3D_H
#define BODY_FORCE_3D_H

#include <vector>

#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "multiBlock/defaultMultiBlockPolicy3D.h"
#include "multiBlock/multiBlock3D.h"
#include "multiBlock/multiDataField3D.h"

namespace plb {

// Implementation of the momentum correction algorithm for applying a constant body force.
template <typename T, template <typename U> class Descriptor>
class AddConstForceToMomentum3D : public BoxProcessingFunctional3D {
public:
    AddConstForceToMomentum3D(Array<T, 3> const &force_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual AddConstForceToMomentum3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    Array<T, 3> force;
};

template <typename T, template <typename U> class Descriptor>
void addConstForceToMomentum(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &rhoBar,
    MultiTensorField3D<T, 3> &j, Array<T, 3> const &force, Box3D const &domain);

// Implementation of the momentum correction algorithm for applying a body force given by a
// user-provided function.
template <typename T, template <typename U> class Descriptor, class ForceFunction>
class AddCustomForceToMomentum3D : public BoxProcessingFunctional3D {
public:
    AddCustomForceToMomentum3D(ForceFunction f_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual AddCustomForceToMomentum3D<T, Descriptor, ForceFunction> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    ForceFunction f;
};

template <typename T, template <class U> class Descriptor, class ForceFunction>
void addCustomForceToMomentum(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &rhoBar,
    MultiTensorField3D<T, 3> &j, ForceFunction f, Box3D const &domain);

// Implementation of the momentum correction algorithm for applying a body force.
template <typename T, template <typename U> class Descriptor>
class AddForceToMomentum3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual AddForceToMomentum3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
void addForceToMomentum(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &rhoBar,
    MultiTensorField3D<T, 3> &j, MultiTensorField3D<T, 3> &force, Box3D const &domain);

// Implementation of the momentum correction algorithm for applying a constant body force only to
// the wet nodes.
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceAddConstForceToMomentum3D : public BoxProcessingFunctional3D {
public:
    FreeSurfaceAddConstForceToMomentum3D(Array<T, 3> const &force_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual FreeSurfaceAddConstForceToMomentum3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    Array<T, 3> force;
};

template <typename T, template <typename U> class Descriptor>
void freeSurfaceAddConstForceToMomentum(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &rhoBar,
    MultiTensorField3D<T, 3> &j, MultiScalarField3D<int> &flag, Array<T, 3> const &force,
    Box3D const &domain);

// Implementation of the momentum correction algorithm for applying a body force only to the wet
// nodes.
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceAddForceToMomentum3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual FreeSurfaceAddForceToMomentum3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
void freeSurfaceAddForceToMomentum(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &rhoBar,
    MultiTensorField3D<T, 3> &j, MultiScalarField3D<int> &flag, MultiTensorField3D<T, 3> &force,
    Box3D const &domain);

// Transform a constant force in a rotating reference frame, and add the Coriolis and centripetal
// forces. The Coriolis force depends on the local velocity. To compute this velocity we need to
// know the body force. For this, we use the body force at the previous iteration, which is expected
// to be contained in the "force" block passed to processGenericBlocks (fourth block in the vector).
// This data processor overwrites the "force" block with the new body force computed.
template <typename T, template <typename U> class Descriptor>
class ComputeRotatingFrameForce3D : public BoxProcessingFunctional3D {
public:
    ComputeRotatingFrameForce3D(
        Array<T, 3> const &constantForce_, Array<T, 3> const &angularVelocity_,
        Array<T, 3> const &origin_, bool incompressibleModel_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual ComputeRotatingFrameForce3D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Array<T, 3> constantForce;
    Array<T, 3> angularVelocity;
    Array<T, 3> origin;
    bool incompressibleModel;
};

template <typename T, template <typename U> class Descriptor>
void computeRotatingFrameForce(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &rhoBar,
    MultiTensorField3D<T, 3> &j, MultiTensorField3D<T, 3> &force, Array<T, 3> const &constantForce,
    Array<T, 3> const &angularVelocity, Array<T, 3> const &origin, bool incompressibleModel,
    Box3D domain);

// Transform a constant force in a rotating reference frame, and add the Coriolis and centripetal
// forces only on wet nodes (on the rest of the cells, the force is set to zero). The Coriolis force
// depends on the local velocity. To compute this velocity we need to know the body force. For this,
// we use the body force at the previous iteration, which is expected to be contained in the "force"
// block passed to processGenericBlocks (fifth block in the vector). This data processor overwrites
// the "force" block with the new body force computed.
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceComputeRotatingFrameForce3D : public BoxProcessingFunctional3D {
public:
    FreeSurfaceComputeRotatingFrameForce3D(
        Array<T, 3> const &constantForce_, Array<T, 3> const &angularVelocity_,
        Array<T, 3> const &origin_, bool incompressibleModel_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual FreeSurfaceComputeRotatingFrameForce3D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Array<T, 3> constantForce;
    Array<T, 3> angularVelocity;
    Array<T, 3> origin;
    bool incompressibleModel;
};

template <typename T, template <typename U> class Descriptor>
void freeSurfaceComputeRotatingFrameForce(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &rhoBar,
    MultiTensorField3D<T, 3> &j, MultiScalarField3D<int> &flag, MultiTensorField3D<T, 3> &force,
    Array<T, 3> const &constantForce, Array<T, 3> const &angularVelocity, Array<T, 3> const &origin,
    bool incompressibleModel, Box3D domain);

template <typename T>
class ComputeTorqueFromBodyForce3D : public ReductiveBoxProcessingFunctional3D_T<T, 3> {
public:
    ComputeTorqueFromBodyForce3D(Array<T, 3> const &center_);
    virtual void process(Box3D domain, TensorField3D<T, 3> &force);
    virtual ComputeTorqueFromBodyForce3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    Array<T, 3> getTorque() const;

private:
    Array<T, 3> center;
    Array<plint, 3> torqueIds;
};

template <typename T>
class MaskedComputeTorqueFromBodyForce3D : public ReductiveBoxProcessingFunctional3D_ST<int, T, 3> {
public:
    MaskedComputeTorqueFromBodyForce3D(Array<T, 3> const &center_, int flag_);
    virtual void process(Box3D domain, ScalarField3D<int> &mask, TensorField3D<T, 3> &force);
    virtual MaskedComputeTorqueFromBodyForce3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    Array<T, 3> getTorque() const;

private:
    Array<T, 3> center;
    int flag;
    Array<plint, 3> torqueIds;
};

template <typename T>
Array<T, 3> computeTorqueFromBodyForce(
    MultiTensorField3D<T, 3> &force, Array<T, 3> const &center, Box3D domain);

template <typename T>
Array<T, 3> computeTorqueFromBodyForce(
    MultiTensorField3D<T, 3> &force, MultiScalarField3D<int> &mask, Array<T, 3> const &center,
    int flag, Box3D domain);

}  // namespace plb

#endif  // BODY_FORCE_3D_H
