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
 * Functionals for domain initialization -- header file.
 */
#ifndef DATA_INITIALIZER_FUNCTIONAL_3D_H
#define DATA_INITIALIZER_FUNCTIONAL_3D_H

#include "atomicBlock/dataProcessingFunctional3D.h"
#include "core/dynamics.h"
#include "core/globalDefs.h"

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Initialization of the block-lattice *************** */
/* ******************************************************************* */

template <typename T, template <class U> class Descriptor>
struct OneCellFunctional3D {
    virtual ~OneCellFunctional3D();
    virtual OneCellFunctional3D<T, Descriptor> *clone() const = 0;
    virtual void execute(Cell<T, Descriptor> &cell) const = 0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual void setscale(int dxScale, int dtScale);
};

template <typename T, template <class U> class Descriptor>
struct OneCellIndexedFunctional3D {
    virtual ~OneCellIndexedFunctional3D();
    virtual OneCellIndexedFunctional3D<T, Descriptor> *clone() const = 0;
    virtual void execute(plint iX, plint iY, plint iZ, Cell<T, Descriptor> &cell) const = 0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual void setscale(int dxScale, int dtScale);
};

template <typename T, template <class U> class Descriptor>
struct OneCellIndexedWithRandFunctional3D {
    virtual ~OneCellIndexedWithRandFunctional3D();
    virtual OneCellIndexedWithRandFunctional3D<T, Descriptor> *clone() const = 0;
    virtual void execute(
        plint iX, plint iY, plint iZ, T randVal, Cell<T, Descriptor> &cell) const = 0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual void setscale(int dxScale, int dtScale);
};

struct DomainFunctional3D {
    virtual ~DomainFunctional3D() { }
    virtual bool operator()(plint iX, plint iY, plint iZ) const = 0;
    virtual DomainFunctional3D *clone() const = 0;
};

template <typename T, template <class U> class Descriptor>
class GenericLatticeFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    GenericLatticeFunctional3D(OneCellFunctional3D<T, Descriptor> *f_);
    GenericLatticeFunctional3D(GenericLatticeFunctional3D<T, Descriptor> const &rhs);
    virtual ~GenericLatticeFunctional3D();
    GenericLatticeFunctional3D<T, Descriptor> &operator=(
        GenericLatticeFunctional3D<T, Descriptor> const &rhs);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual GenericLatticeFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void setscale(int dxScale, int dtScale);

private:
    OneCellFunctional3D<T, Descriptor> *f;
};

template <typename T, template <class U> class Descriptor>
class GenericIndexedLatticeFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    GenericIndexedLatticeFunctional3D(OneCellIndexedFunctional3D<T, Descriptor> *f_);
    GenericIndexedLatticeFunctional3D(GenericIndexedLatticeFunctional3D<T, Descriptor> const &rhs);
    virtual ~GenericIndexedLatticeFunctional3D();
    GenericIndexedLatticeFunctional3D<T, Descriptor> &operator=(
        GenericIndexedLatticeFunctional3D<T, Descriptor> const &rhs);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual GenericIndexedLatticeFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void setscale(int dxScale, int dtScale);

private:
    OneCellIndexedFunctional3D<T, Descriptor> *f;
};

template <typename T, template <class U> class Descriptor>
class GenericIndexedWithRandLatticeFunctional3D :
    public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    GenericIndexedWithRandLatticeFunctional3D(
        OneCellIndexedWithRandFunctional3D<T, Descriptor> *f_, Box3D boundingBox_,
        uint32_t const *seed_);
    GenericIndexedWithRandLatticeFunctional3D(
        GenericIndexedWithRandLatticeFunctional3D<T, Descriptor> const &rhs);
    virtual ~GenericIndexedWithRandLatticeFunctional3D();
    GenericIndexedWithRandLatticeFunctional3D<T, Descriptor> &operator=(
        GenericIndexedWithRandLatticeFunctional3D<T, Descriptor> const &rhs);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual GenericIndexedWithRandLatticeFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void setscale(int dxScale, int dtScale);

private:
    OneCellIndexedWithRandFunctional3D<T, Descriptor> *f;
    plint nY, nZ;
    uint32_t const *seed;
};

/* *************** Class InstantiateDynamicsFunctional3D ************* */

template <typename T, template <typename U> class Descriptor>
class InstantiateDynamicsFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    InstantiateDynamicsFunctional3D(Dynamics<T, Descriptor> *dynamics_);
    InstantiateDynamicsFunctional3D(InstantiateDynamicsFunctional3D<T, Descriptor> const &rhs);
    InstantiateDynamicsFunctional3D<T, Descriptor> &operator=(
        InstantiateDynamicsFunctional3D<T, Descriptor> const &rhs);
    virtual ~InstantiateDynamicsFunctional3D();
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual InstantiateDynamicsFunctional3D<T, Descriptor> *clone() const;

private:
    Dynamics<T, Descriptor> *dynamics;
};

/* *************** Class InstantiateDynamicsInBulkAndEnvelopeFunctional3D ************* */

// The following data processor is a "low-level" one to be used only in very specific cases, where
// one wants to assign the dynamics inside the envelope, and the standard
// "InstantiateDynamicsFunctional3D" functional (which is applied only in the bulk) will not be
// adequate because of an "irregular" sparse block decomposition. What we mean by that is, that
// sometimes the sparse block structure is such, that after communication, the envelopes will not
// have the desired values. In such cases the standard "InstantiateDynamicsFunctional3D" functional
// will not be sufficient, and the use of the following
// "InstantiateDynamicsInBulkAndEnvelopeFunctional3D" will be preferred.
// This data processor is the same as InstantiateDynamicsFunctional3D, except from the
// implementation of the "appliesTo" function.
template <typename T, template <typename U> class Descriptor>
class InstantiateDynamicsInBulkAndEnvelopeFunctional3D :
    public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    InstantiateDynamicsInBulkAndEnvelopeFunctional3D(Dynamics<T, Descriptor> *dynamics_);
    InstantiateDynamicsInBulkAndEnvelopeFunctional3D(
        InstantiateDynamicsInBulkAndEnvelopeFunctional3D<T, Descriptor> const &rhs);
    InstantiateDynamicsInBulkAndEnvelopeFunctional3D<T, Descriptor> &operator=(
        InstantiateDynamicsInBulkAndEnvelopeFunctional3D<T, Descriptor> const &rhs);
    virtual ~InstantiateDynamicsInBulkAndEnvelopeFunctional3D();
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual InstantiateDynamicsInBulkAndEnvelopeFunctional3D<T, Descriptor> *clone() const;

private:
    Dynamics<T, Descriptor> *dynamics;
};

/* ************* Class InstantiateComplexDomainDynamicsFunctional3D ** */

template <typename T, template <typename U> class Descriptor>
class InstantiateComplexDomainDynamicsFunctional3D :
    public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    InstantiateComplexDomainDynamicsFunctional3D(
        Dynamics<T, Descriptor> *dynamics_, DomainFunctional3D *domain_);
    InstantiateComplexDomainDynamicsFunctional3D(
        InstantiateComplexDomainDynamicsFunctional3D<T, Descriptor> const &rhs);
    InstantiateComplexDomainDynamicsFunctional3D<T, Descriptor> &operator=(
        InstantiateComplexDomainDynamicsFunctional3D<T, Descriptor> const &rhs);
    virtual ~InstantiateComplexDomainDynamicsFunctional3D();
    virtual void process(Box3D boundingBox, BlockLattice3D<T, Descriptor> &lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual InstantiateComplexDomainDynamicsFunctional3D<T, Descriptor> *clone() const;

private:
    Dynamics<T, Descriptor> *dynamics;
    DomainFunctional3D *domain;
};

/* ************* Class InstantiateDotDynamicsFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
class InstantiateDotDynamicsFunctional3D : public DotProcessingFunctional3D_L<T, Descriptor> {
public:
    InstantiateDotDynamicsFunctional3D(Dynamics<T, Descriptor> *dynamics_);
    InstantiateDotDynamicsFunctional3D(
        InstantiateDotDynamicsFunctional3D<T, Descriptor> const &rhs);
    InstantiateDotDynamicsFunctional3D<T, Descriptor> &operator=(
        InstantiateDotDynamicsFunctional3D<T, Descriptor> const &rhs);
    virtual ~InstantiateDotDynamicsFunctional3D();
    virtual void process(DotList3D const &dotList, BlockLattice3D<T, Descriptor> &lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual InstantiateDotDynamicsFunctional3D<T, Descriptor> *clone() const;

private:
    Dynamics<T, Descriptor> *dynamics;
};

/* ************* Class DynamicsFromMaskFunctional3D ************************ */

/// Assign dynamics to nodes specified by a boolean mask.
template <typename T, template <typename U> class Descriptor>
class DynamicsFromMaskFunctional3D : public BoxProcessingFunctional3D_LS<T, Descriptor, bool> {
public:
    DynamicsFromMaskFunctional3D(Dynamics<T, Descriptor> *dynamics_, bool whichFlag_);
    DynamicsFromMaskFunctional3D(DynamicsFromMaskFunctional3D<T, Descriptor> const &rhs);
    DynamicsFromMaskFunctional3D<T, Descriptor> &operator=(
        DynamicsFromMaskFunctional3D<T, Descriptor> const &rhs);
    virtual ~DynamicsFromMaskFunctional3D();
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<bool> &mask);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual DynamicsFromMaskFunctional3D<T, Descriptor> *clone() const;

private:
    Dynamics<T, Descriptor> *dynamics;
    bool whichFlag;
};

/* ************* Class DynamicsFromIntMaskFunctional3D ************************ */

/// Assign dynamics to nodes specified by a boolean mask.
template <typename T, template <typename U> class Descriptor>
class DynamicsFromIntMaskFunctional3D : public BoxProcessingFunctional3D_LS<T, Descriptor, int> {
public:
    DynamicsFromIntMaskFunctional3D(Dynamics<T, Descriptor> *dynamics_, int whichFlag_);
    DynamicsFromIntMaskFunctional3D(DynamicsFromIntMaskFunctional3D<T, Descriptor> const &rhs);
    DynamicsFromIntMaskFunctional3D<T, Descriptor> &operator=(
        DynamicsFromIntMaskFunctional3D<T, Descriptor> const &rhs);
    virtual ~DynamicsFromIntMaskFunctional3D();
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &mask);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual DynamicsFromIntMaskFunctional3D<T, Descriptor> *clone() const;

private:
    Dynamics<T, Descriptor> *dynamics;
    int whichFlag;
};

/* *************** Class RecomposeFromOrderZeroVariablesFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
class RecomposeFromOrderZeroVariablesFunctional3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual RecomposeFromOrderZeroVariablesFunctional3D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

/* *************** Class RecomposeFromFlowVariablesFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
class RecomposeFromFlowVariablesFunctional3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual RecomposeFromFlowVariablesFunctional3D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

/* ************* Class AssignOmegaFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
class AssignOmegaFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    AssignOmegaFunctional3D(T omega);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual AssignOmegaFunctional3D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T omega;
};

/* ************* Class AssignScalarFieldOmegaFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
class AssignScalarFieldOmegaFunctional3D : public BoxProcessingFunctional3D_LS<T, Descriptor, T> {
public:
    AssignScalarFieldOmegaFunctional3D();
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &omega);
    virtual AssignScalarFieldOmegaFunctional3D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

/* ************* Class SetConstBoundaryVelocityFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
class SetConstBoundaryVelocityFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    SetConstBoundaryVelocityFunctional3D(Array<T, Descriptor<T>::d> velocity);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual SetConstBoundaryVelocityFunctional3D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Array<T, Descriptor<T>::d> u;
};

/* ************* Class SetConstBoundaryVelocityWithTensorForceFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
class SetConstBoundaryVelocityWithTensorForceFunctional3D :
    public BoxProcessingFunctional3D_LT<T, Descriptor, T, Descriptor<T>::d> {
public:
    SetConstBoundaryVelocityWithTensorForceFunctional3D(Array<T, Descriptor<T>::d> velocity);
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        TensorField3D<T, Descriptor<T>::d> &force);
    virtual SetConstBoundaryVelocityWithTensorForceFunctional3D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Array<T, Descriptor<T>::d> u;
};

/* ************* Class SetConstBoundaryVelocityWithForceFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
class SetConstBoundaryVelocityWithForceFunctional3D :
    public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    SetConstBoundaryVelocityWithForceFunctional3D(
        Array<T, Descriptor<T>::d> force, Array<T, Descriptor<T>::d> velocity);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual SetConstBoundaryVelocityWithForceFunctional3D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Array<T, Descriptor<T>::d> f;
    Array<T, Descriptor<T>::d> u;
};

/* ************* Class SetCustomBoundaryVelocityFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class VelocityFunction>
class SetCustomBoundaryVelocityFunctional3D : public OneCellIndexedFunctional3D<T, Descriptor> {
public:
    SetCustomBoundaryVelocityFunctional3D(VelocityFunction f_);
    virtual SetCustomBoundaryVelocityFunctional3D<T, Descriptor, VelocityFunction> *clone() const;
    virtual void execute(plint iX, plint iY, plint iZ, Cell<T, Descriptor> &cell) const;
    virtual void setscale(int dxScale, int dtScale);

private:
    VelocityFunction f;
    T velocityScale;
};

/* ************* Class SetCustomBoundaryVelocityWithTensorForceFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class VelocityFunction>
class SetCustomBoundaryVelocityWithTensorForceFunctional3D :
    public BoxProcessingFunctional3D_LT<T, Descriptor, T, Descriptor<T>::d> {
public:
    SetCustomBoundaryVelocityWithTensorForceFunctional3D(VelocityFunction f_);
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        TensorField3D<T, Descriptor<T>::d> &force);
    virtual SetCustomBoundaryVelocityWithTensorForceFunctional3D<T, Descriptor, VelocityFunction>
        *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    VelocityFunction f;
};

/* ************* Class SetCustomBoundaryVelocityWithForceFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class VelocityFunction>
class SetCustomBoundaryVelocityWithForceFunctional3D :
    public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    SetCustomBoundaryVelocityWithForceFunctional3D(
        Array<T, Descriptor<T>::d> force_, VelocityFunction f_);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual SetCustomBoundaryVelocityWithForceFunctional3D<T, Descriptor, VelocityFunction> *clone()
        const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Array<T, Descriptor<T>::d> force;
    VelocityFunction f;
};

/* ************* Class SetCustomBoundaryVelocityWithCustomForceFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class ForceVelocityFunction>
class SetCustomBoundaryVelocityWithCustomForceFunctional3D :
    public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    SetCustomBoundaryVelocityWithCustomForceFunctional3D(ForceVelocityFunction f_);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual SetCustomBoundaryVelocityWithCustomForceFunctional3D<
        T, Descriptor, ForceVelocityFunction>
        *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    ForceVelocityFunction f;
};

/* ************* Class SetConstBoundaryDensityFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
class SetConstBoundaryDensityFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    SetConstBoundaryDensityFunctional3D(T rho_);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual SetConstBoundaryDensityFunctional3D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T rho;
};

/* ************* Class SetCustomBoundaryDensityFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class DensityFunction>
class SetCustomBoundaryDensityFunctional3D : public OneCellIndexedFunctional3D<T, Descriptor> {
public:
    SetCustomBoundaryDensityFunctional3D(DensityFunction f_);
    virtual void execute(plint iX, plint iY, plint iZ, Cell<T, Descriptor> &cell) const;
    virtual SetCustomBoundaryDensityFunctional3D<T, Descriptor, DensityFunction> *clone() const;

private:
    DensityFunction f;
};

/* ************* Class IniConstEquilibriumFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
class IniConstEquilibriumFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    IniConstEquilibriumFunctional3D(T density, Array<T, Descriptor<T>::d> velocity, T temperature);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual IniConstEquilibriumFunctional3D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T rho;
    T rhoBar;
    Array<T, Descriptor<T>::d> u;
    T thetaBar;
};

/* ************* Class MaskedIniConstEquilibriumFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
class MaskedIniConstEquilibriumFunctional3D :
    public BoxProcessingFunctional3D_LS<T, Descriptor, int> {
public:
    MaskedIniConstEquilibriumFunctional3D(
        T density, Array<T, Descriptor<T>::d> velocity, T temperature, int whichFlag_);
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &mask);
    virtual MaskedIniConstEquilibriumFunctional3D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T rho;
    T rhoBar;
    Array<T, Descriptor<T>::d> u;
    T thetaBar;
    int whichFlag;
};

/* ************* Class IniConstTensorForceEquilibriumFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
class IniConstTensorForceEquilibriumFunctional3D :
    public BoxProcessingFunctional3D_LT<T, Descriptor, T, Descriptor<T>::d> {
public:
    IniConstTensorForceEquilibriumFunctional3D(
        T density, Array<T, Descriptor<T>::d> velocity, T temperature);
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        TensorField3D<T, Descriptor<T>::d> &force);
    virtual IniConstTensorForceEquilibriumFunctional3D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T rho;
    T rhoBar;
    Array<T, Descriptor<T>::d> u;
    T thetaBar;
};

/* ************* Class IniCustomTensorForceEquilibriumFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
class IniCustomTensorForceEquilibriumFunctional3D :
    public BoxProcessingFunctional3D_LT<T, Descriptor, T, Descriptor<T>::d> {
public:
    IniCustomTensorForceEquilibriumFunctional3D(RhoUFunction f_, T temperature);
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        TensorField3D<T, Descriptor<T>::d> &force);
    virtual IniCustomTensorForceEquilibriumFunctional3D<T, Descriptor, RhoUFunction> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    RhoUFunction f;
    T thetaBar;
};

/* ************* Class IniCustomTensorForceRandomEquilibriumFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
class IniCustomTensorForceRandomEquilibriumFunctional3D : public BoxProcessingFunctional3D {
public:
    IniCustomTensorForceRandomEquilibriumFunctional3D(RhoUFunction f_, T temperature);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual IniCustomTensorForceRandomEquilibriumFunctional3D<T, Descriptor, RhoUFunction> *clone()
        const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    RhoUFunction f;
    T thetaBar;
};

/* ************* Class IniConstForceEquilibriumFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
class IniConstForceEquilibriumFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    IniConstForceEquilibriumFunctional3D(
        Array<T, Descriptor<T>::d> force, T density, Array<T, Descriptor<T>::d> velocity,
        T temperature);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual IniConstForceEquilibriumFunctional3D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Array<T, Descriptor<T>::d> f;
    T rho;
    T rhoBar;
    Array<T, Descriptor<T>::d> u;
    T thetaBar;
};

/* ************* Class IniCustomForceEquilibriumFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
class IniCustomForceEquilibriumFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    IniCustomForceEquilibriumFunctional3D(
        Array<T, Descriptor<T>::d> force_, RhoUFunction f_, T temperature);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual IniCustomForceEquilibriumFunctional3D<T, Descriptor, RhoUFunction> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Array<T, Descriptor<T>::d> force;
    RhoUFunction f;
    T thetaBar;
};

/* ************* Class IniCustomForceRandomEquilibriumFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
class IniCustomForceRandomEquilibriumFunctional3D : public BoxProcessingFunctional3D {
public:
    IniCustomForceRandomEquilibriumFunctional3D(
        Array<T, Descriptor<T>::d> force_, RhoUFunction f_, T temperature);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual IniCustomForceRandomEquilibriumFunctional3D<T, Descriptor, RhoUFunction> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Array<T, Descriptor<T>::d> force;
    RhoUFunction f;
    T thetaBar;
};

/* ************* Class IniConstEquilibriumOnDomainFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class DomainFunctional>
class IniConstEquilibriumOnDomainFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    IniConstEquilibriumOnDomainFunctional3D(
        T density, Array<T, Descriptor<T>::d> velocity, T temperature,
        DomainFunctional const &domain_);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual IniConstEquilibriumOnDomainFunctional3D<T, Descriptor, DomainFunctional> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T rho;
    T rhoBar;
    Array<T, Descriptor<T>::d> u;
    T thetaBar;
    DomainFunctional domain;
};

/* ************* Class IniCustomEquilibriumFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
class IniCustomEquilibriumFunctional3D : public OneCellIndexedFunctional3D<T, Descriptor> {
public:
    IniCustomEquilibriumFunctional3D(RhoUFunction f_);
    virtual void execute(plint iX, plint iY, plint iZ, Cell<T, Descriptor> &cell) const;
    virtual IniCustomEquilibriumFunctional3D<T, Descriptor, RhoUFunction> *clone() const;
    virtual void setscale(int dxScale, int dtScale);

private:
    RhoUFunction f;
    T velocityScale;
};

/* ************* Class IniCustomRandomEquilibriumFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
class IniCustomRandomEquilibriumFunctional3D :
    public OneCellIndexedWithRandFunctional3D<T, Descriptor> {
public:
    IniCustomRandomEquilibriumFunctional3D(RhoUFunction f_);
    virtual void execute(plint iX, plint iY, plint iZ, T randVal, Cell<T, Descriptor> &cell) const;
    virtual IniCustomRandomEquilibriumFunctional3D<T, Descriptor, RhoUFunction> *clone() const;
    virtual void setscale(int dxScale, int dtScale);

private:
    RhoUFunction f;
    T velocityScale;
};

/* ************* Class SetCustomOmegaFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class OmegaFunction>
class SetCustomOmegaFunctional3D : public OneCellIndexedFunctional3D<T, Descriptor> {
public:
    SetCustomOmegaFunctional3D(OmegaFunction f_);
    virtual void execute(plint iX, plint iY, plint iZ, Cell<T, Descriptor> &cell) const;
    virtual SetCustomOmegaFunctional3D<T, Descriptor, OmegaFunction> *clone() const;

private:
    OmegaFunction f;
};

/* ************* Class IniCustomThermalEquilibriumFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class RhoVelTempFunction>
class IniCustomThermalEquilibriumFunctional3D : public OneCellIndexedFunctional3D<T, Descriptor> {
public:
    IniCustomThermalEquilibriumFunctional3D(RhoVelTempFunction f_);
    virtual void execute(plint iX, plint iY, plint iZ, Cell<T, Descriptor> &cell) const;
    virtual IniCustomThermalEquilibriumFunctional3D<T, Descriptor, RhoVelTempFunction> *clone()
        const;
    virtual void setscale(int dxScale, int dtScale);

private:
    RhoVelTempFunction f;
    T velocityScale;
};

/* ************* Class StripeOffDensityOffsetFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
class StripeOffDensityOffsetFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    StripeOffDensityOffsetFunctional3D(T deltaRho_);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual StripeOffDensityOffsetFunctional3D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T deltaRho;
};

/* ************* Class InstantiateCompositeDynamicsFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
class InstantiateCompositeDynamicsFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    InstantiateCompositeDynamicsFunctional3D(CompositeDynamics<T, Descriptor> *compositeDynamics_);
    InstantiateCompositeDynamicsFunctional3D(
        InstantiateCompositeDynamicsFunctional3D<T, Descriptor> const &rhs);
    InstantiateCompositeDynamicsFunctional3D<T, Descriptor> &operator=(
        InstantiateCompositeDynamicsFunctional3D<T, Descriptor> const &rhs);
    virtual ~InstantiateCompositeDynamicsFunctional3D();
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual InstantiateCompositeDynamicsFunctional3D<T, Descriptor> *clone() const;

private:
    CompositeDynamics<T, Descriptor> *compositeDynamics;
};

/* ************* Class SetExternalScalarFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
class SetExternalScalarFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    SetExternalScalarFunctional3D(int whichScalar_, T externalScalar_);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual SetExternalScalarFunctional3D<T, Descriptor> *clone() const;

private:
    int whichScalar;
    T externalScalar;
};

/* ************* Class SetExternalScalarFromScalarFieldFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
class SetExternalScalarFromScalarFieldFunctional3D :
    public BoxProcessingFunctional3D_LS<T, Descriptor, T> {
public:
    SetExternalScalarFromScalarFieldFunctional3D(int whichScalar_);
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalar);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual SetExternalScalarFromScalarFieldFunctional3D<T, Descriptor> *clone() const;

private:
    int whichScalar;
};

/* ************* Class MaskedSetExternalScalarFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
class MaskedSetExternalScalarFunctional3D :
    public BoxProcessingFunctional3D_LS<T, Descriptor, int> {
public:
    MaskedSetExternalScalarFunctional3D(int flag_, int whichScalar_, T externalScalar_);
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &mask);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual MaskedSetExternalScalarFunctional3D<T, Descriptor> *clone() const;

private:
    int flag;
    int whichScalar;
    T externalScalar;
};

/* ************* Class SetGenericExternalScalarFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class Functional>
class SetGenericExternalScalarFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    SetGenericExternalScalarFunctional3D(int whichScalar_, Functional const &functional_);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual SetGenericExternalScalarFunctional3D<T, Descriptor, Functional> *clone() const;

private:
    int whichScalar;
    Functional functional;
};

/* ************* Class MaskedSetGenericExternalScalarFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class Functional>
class MaskedSetGenericExternalScalarFunctional3D :
    public BoxProcessingFunctional3D_LS<T, Descriptor, int> {
public:
    MaskedSetGenericExternalScalarFunctional3D(
        int flag_, int whichScalar_, Functional const &functional_);
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &mask);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual MaskedSetGenericExternalScalarFunctional3D<T, Descriptor, Functional> *clone() const;

private:
    int flag;
    int whichScalar;
    Functional functional;
};

/* ************* Class AddToExternalScalarFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
class AddToExternalScalarFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    AddToExternalScalarFunctional3D(int whichScalar_, T externalScalar_);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual AddToExternalScalarFunctional3D<T, Descriptor> *clone() const;

private:
    int whichScalar;
    T externalScalar;
};

/* ************* Class SetExternalVectorFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
class SetExternalVectorFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    SetExternalVectorFunctional3D(int vectorStartsAt_, Array<T, Descriptor<T>::d> &externalVector_);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual SetExternalVectorFunctional3D<T, Descriptor> *clone() const;

private:
    int vectorStartsAt;
    Array<T, Descriptor<T>::d> externalVector;
};

/* ************* Class MaskedSetExternalVectorFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
class MaskedSetExternalVectorFunctional3D :
    public BoxProcessingFunctional3D_LS<T, Descriptor, int> {
public:
    MaskedSetExternalVectorFunctional3D(
        int flag_, int vectorStartsAt_, Array<T, Descriptor<T>::d> &externalVector_);
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &mask);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual MaskedSetExternalVectorFunctional3D<T, Descriptor> *clone() const;

private:
    int flag;
    int vectorStartsAt;
    Array<T, Descriptor<T>::d> externalVector;
};

/* ************* Class SetGenericExternalVectorFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class Functional>
class SetGenericExternalVectorFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    SetGenericExternalVectorFunctional3D(int vectorBeginsAt_, Functional const &functional_);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual SetGenericExternalVectorFunctional3D<T, Descriptor, Functional> *clone() const;

private:
    int vectorBeginsAt;
    Functional functional;
};

/* ************* Class MaskedSetGenericExternalVectorFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, class Functional>
class MaskedSetGenericExternalVectorFunctional3D :
    public BoxProcessingFunctional3D_LS<T, Descriptor, int> {
public:
    MaskedSetGenericExternalVectorFunctional3D(
        int flag_, int vectorBeginsAt_, Functional const &functional_);
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &mask);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual MaskedSetGenericExternalVectorFunctional3D<T, Descriptor, Functional> *clone() const;

private:
    int flag;
    int vectorBeginsAt;
    Functional functional;
};

/* ************* Class SetExternalVectorFromTensorFieldFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor, int nDim>
class SetExternalVectorFromTensorFieldFunctional3D :
    public BoxProcessingFunctional3D_LT<T, Descriptor, T, nDim> {
public:
    SetExternalVectorFromTensorFieldFunctional3D(int vectorStartsAt_);
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, TensorField3D<T, nDim> &tensor);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual SetExternalVectorFromTensorFieldFunctional3D<T, Descriptor, nDim> *clone() const;

private:
    int vectorStartsAt;
};

/* ************* Class InterpolatePopulationsFunctional3D ******************* */

template <typename T, template <typename U> class Descriptor>
class InterpolatePopulationsFunctional3D :
    public BoxProcessingFunctional3D_LL<T, Descriptor, T, Descriptor> {
public:
    InterpolatePopulationsFunctional3D(plint minIter_, plint maxIter_);
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice1,
        BlockLattice3D<T, Descriptor> &lattice2);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual InterpolatePopulationsFunctional3D<T, Descriptor> *clone() const;

private:
    plint minIter, maxIter;
};

/* *************** PART II ******************************************* */
/* *************** Initialization of scalar- and tensor-fields ******* */
/* ******************************************************************* */

template <typename T>
class IniConstScalarFunctional3D : public BoxProcessingFunctional3D_S<T> {
public:
    IniConstScalarFunctional3D(T value_);
    virtual void process(Box3D domain, ScalarField3D<T> &field);
    virtual IniConstScalarFunctional3D<T> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T value;
};

template <typename T>
class MaskedIniConstScalarFunctional3D : public BoxProcessingFunctional3D_SS<T, int> {
public:
    MaskedIniConstScalarFunctional3D(int flag_, T value_);
    virtual void process(Box3D domain, ScalarField3D<T> &field, ScalarField3D<int> &mask);
    virtual MaskedIniConstScalarFunctional3D<T> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    int flag;
    T value;
};

template <typename T>
class MaskedIniConstScalarFunctional3D_N : public BoxProcessingFunctional3D_SN<T, int> {
public:
    MaskedIniConstScalarFunctional3D_N(int flag_, T value_);
    virtual void process(Box3D domain, ScalarField3D<T> &field, NTensorField3D<int> &mask);
    virtual MaskedIniConstScalarFunctional3D_N<T> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    int flag;
    T value;
};

template <typename T, class Function>
class SetToScalarFunctionFunctional3D : public BoxProcessingFunctional3D_S<T> {
public:
    SetToScalarFunctionFunctional3D(Function f_);
    virtual void process(Box3D domain, ScalarField3D<T> &field);
    virtual SetToScalarFunctionFunctional3D<T, Function> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Function f;
};

template <typename T, class Function>
class SetToNTensorFunctionFunctional3D : public BoxProcessingFunctional3D_N<T> {
public:
    SetToNTensorFunctionFunctional3D(Function f_);
    virtual void process(Box3D domain, NTensorField3D<T> &field);
    virtual SetToNTensorFunctionFunctional3D<T, Function> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Function f;
};

template <typename T, class Function>
class AnalyticalSetRhoBarJFunctional3D : public BoxProcessingFunctional3D_N<T> {
public:
    AnalyticalSetRhoBarJFunctional3D(Function const &function_);
    virtual void process(Box3D domain, NTensorField3D<T> &rhoBarJ);
    virtual AnalyticalSetRhoBarJFunctional3D<T, Function> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Function function;
};

template <typename T, int nDim>
class IniConstTensorFunctional3D : public BoxProcessingFunctional3D_T<T, nDim> {
public:
    IniConstTensorFunctional3D(Array<T, nDim> const &value_);
    virtual void process(Box3D domain, TensorField3D<T, nDim> &field);
    virtual IniConstTensorFunctional3D<T, nDim> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Array<T, nDim> value;
};

template <typename T, int nDim>
class MaskedIniConstTensorFunctional3D : public BoxProcessingFunctional3D_ST<int, T, nDim> {
public:
    MaskedIniConstTensorFunctional3D(int flag_, Array<T, nDim> const &value_);
    virtual void process(Box3D domain, ScalarField3D<int> &mask, TensorField3D<T, nDim> &field);
    virtual MaskedIniConstTensorFunctional3D<T, nDim> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    int flag;
    Array<T, nDim> value;
};

template <typename T, int nDim>
class MaskedIniConstTensorFunctional3D_N : public BoxProcessingFunctional3D {
public:
    MaskedIniConstTensorFunctional3D_N(int flag_, Array<T, nDim> const &value_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual MaskedIniConstTensorFunctional3D_N<T, nDim> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    int flag;
    Array<T, nDim> value;
};

template <typename T>
class SwapValuesBulkAndEnvelope3D_N : public BoxProcessingFunctional3D_NN<T, T> {
public:
    virtual void process(Box3D domain, NTensorField3D<T> &A, NTensorField3D<T> &B);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual SwapValuesBulkAndEnvelope3D_N<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

template <typename T, int nDim, class Function>
class SetToTensorFunctionFunctional3D : public BoxProcessingFunctional3D_T<T, nDim> {
public:
    SetToTensorFunctionFunctional3D(Function f_);
    virtual void process(Box3D domain, TensorField3D<T, nDim> &field);
    virtual SetToTensorFunctionFunctional3D<T, nDim, Function> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Function f;
};

template <typename T>
class SetToCoordinateFunctional3D : public BoxProcessingFunctional3D_S<T> {
public:
    SetToCoordinateFunctional3D(plint index_);
    virtual void process(Box3D domain, ScalarField3D<T> &field);
    virtual SetToCoordinateFunctional3D<T> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    plint index;
};

template <typename T>
class SetToCoordinatesFunctional3D : public BoxProcessingFunctional3D_T<T, 3> {
public:
    SetToCoordinatesFunctional3D();
    virtual void process(Box3D domain, TensorField3D<T, 3> &field);
    virtual SetToCoordinatesFunctional3D<T> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

template <typename T>
class SetToRandomFunctional3D : public BoxProcessingFunctional3D_S<T> {
public:
    SetToRandomFunctional3D(Box3D boundingBox, uint32_t const *seed_);
    virtual void process(Box3D domain, ScalarField3D<T> &field);
    virtual SetToRandomFunctional3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    plint nY, nZ;
    uint32_t const *seed;
};

template <typename T, int nDim>
class SetTensorComponentFunctional3D : public BoxProcessingFunctional3D_ST<T, T, nDim> {
public:
    SetTensorComponentFunctional3D(int whichDim_);
    virtual void process(
        Box3D domain, ScalarField3D<T> &scalarField, TensorField3D<T, nDim> &tensorField);
    virtual SetTensorComponentFunctional3D<T, nDim> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    int whichDim;
};

template <typename T>
class PropagateInZdirection3D : public BoxProcessingFunctional3D_S<T> {
public:
    virtual void process(Box3D domain, ScalarField3D<T> &field);
    virtual PropagateInZdirection3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

template <typename T>
class GrowDomainFunctional3D : public BoxProcessingFunctional3D_S<T> {
public:
    GrowDomainFunctional3D(T flag_);
    virtual void process(Box3D domain, ScalarField3D<T> &voxels);
    virtual GrowDomainFunctional3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    T flag;
};

// This functional fills the envelopes (of a scalar field) which cannot be updated by
// communication (due to an irregular block structure), with values extrapolated from
// the bulk. This functional must be applied on the full bounding box of the
// corresponding MultiScalarField3D. The envelopeWidth is assumed to be 1.
template <typename T>
class ScalarNeumannToUnusedEnvelopes3D : public BoxProcessingFunctional3D_S<T> {
public:
    virtual void process(Box3D domain, ScalarField3D<T> &field);
    virtual ScalarNeumannToUnusedEnvelopes3D<T> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

// This functional fills the envelopes (of a tensor field) which cannot be updated by
// communication (due to an irregular block structure), with values extrapolated from
// the bulk. This functional must be applied on the full bounding box of the
// corresponding MultiTensorField3D. The envelopeWidth is assumed to be 1.
template <typename T, int nDim>
class TensorNeumannToUnusedEnvelopes3D : public BoxProcessingFunctional3D_T<T, nDim> {
public:
    virtual void process(Box3D domain, TensorField3D<T, nDim> &field);
    virtual TensorNeumannToUnusedEnvelopes3D<T, nDim> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

}  // namespace plb

#endif  // DATA_INITIALIZER_FUNCTIONAL_3D_H

// Explicitly include generic algorithms which are never precompiled (not even in precompiled
// version)
#include "dataProcessors/dataInitializerGenerics3D.h"
