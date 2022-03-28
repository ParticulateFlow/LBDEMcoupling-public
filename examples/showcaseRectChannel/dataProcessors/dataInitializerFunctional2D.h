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
#ifndef DATA_INITIALIZER_FUNCTIONAL_2D_H
#define DATA_INITIALIZER_FUNCTIONAL_2D_H

#include "atomicBlock/dataProcessingFunctional2D.h"
#include "core/dynamics.h"
#include "core/globalDefs.h"

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Initialization of the block-lattice *************** */
/* ******************************************************************* */

template <typename T, template <class U> class Descriptor>
struct OneCellFunctional2D {
    virtual ~OneCellFunctional2D();
    virtual OneCellFunctional2D<T, Descriptor> *clone() const = 0;
    virtual void execute(Cell<T, Descriptor> &cell) const = 0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual void setscale(int dxScale, int dtScale);
};

template <typename T, template <class U> class Descriptor>
struct OneCellIndexedFunctional2D {
    virtual ~OneCellIndexedFunctional2D();
    virtual OneCellIndexedFunctional2D<T, Descriptor> *clone() const = 0;
    virtual void execute(plint iX, plint iY, Cell<T, Descriptor> &cell) const = 0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual void setscale(int dxScale, int dtScale);
};

template <typename T, template <class U> class Descriptor>
struct OneCellIndexedWithRandFunctional2D {
    virtual ~OneCellIndexedWithRandFunctional2D();
    virtual OneCellIndexedWithRandFunctional2D<T, Descriptor> *clone() const = 0;
    virtual void execute(plint iX, plint iY, T randVal, Cell<T, Descriptor> &cell) const = 0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual void setscale(int dxScale, int dtScale);
};

struct DomainFunctional2D {
    virtual ~DomainFunctional2D() { }
    virtual bool operator()(plint iX, plint iY) const = 0;
    virtual DomainFunctional2D *clone() const = 0;
};

template <typename T, template <class U> class Descriptor>
class GenericLatticeFunctional2D : public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    GenericLatticeFunctional2D(OneCellFunctional2D<T, Descriptor> *f_);
    GenericLatticeFunctional2D(GenericLatticeFunctional2D<T, Descriptor> const &rhs);
    virtual ~GenericLatticeFunctional2D();
    GenericLatticeFunctional2D<T, Descriptor> &operator=(
        GenericLatticeFunctional2D<T, Descriptor> const &rhs);
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual GenericLatticeFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void setscale(int dxScale, int dtScale);

private:
    OneCellFunctional2D<T, Descriptor> *f;
};

template <typename T, template <class U> class Descriptor>
class GenericIndexedLatticeFunctional2D : public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    GenericIndexedLatticeFunctional2D(OneCellIndexedFunctional2D<T, Descriptor> *f_);
    GenericIndexedLatticeFunctional2D(GenericIndexedLatticeFunctional2D<T, Descriptor> const &rhs);
    virtual ~GenericIndexedLatticeFunctional2D();
    GenericIndexedLatticeFunctional2D<T, Descriptor> &operator=(
        GenericIndexedLatticeFunctional2D<T, Descriptor> const &rhs);
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual GenericIndexedLatticeFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void setscale(int dxScale, int dtScale);

private:
    OneCellIndexedFunctional2D<T, Descriptor> *f;
};

template <typename T, template <class U> class Descriptor>
class GenericIndexedWithRandLatticeFunctional2D :
    public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    GenericIndexedWithRandLatticeFunctional2D(
        OneCellIndexedWithRandFunctional2D<T, Descriptor> *f_, Box2D boundingBox_,
        uint32_t const *seed_);
    GenericIndexedWithRandLatticeFunctional2D(
        GenericIndexedWithRandLatticeFunctional2D<T, Descriptor> const &rhs);
    virtual ~GenericIndexedWithRandLatticeFunctional2D();
    GenericIndexedWithRandLatticeFunctional2D<T, Descriptor> &operator=(
        GenericIndexedWithRandLatticeFunctional2D<T, Descriptor> const &rhs);
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual GenericIndexedWithRandLatticeFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void setscale(int dxScale, int dtScale);

private:
    OneCellIndexedWithRandFunctional2D<T, Descriptor> *f;
    plint nY;
    uint32_t const *seed;
};

/* *************** Class InstantiateDynamicsFunctional2D ************* */

template <typename T, template <typename U> class Descriptor>
class InstantiateDynamicsFunctional2D : public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    InstantiateDynamicsFunctional2D(Dynamics<T, Descriptor> *dynamics_);
    InstantiateDynamicsFunctional2D(InstantiateDynamicsFunctional2D<T, Descriptor> const &rhs);
    InstantiateDynamicsFunctional2D<T, Descriptor> &operator=(
        InstantiateDynamicsFunctional2D<T, Descriptor> const &rhs);
    virtual ~InstantiateDynamicsFunctional2D();
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual InstantiateDynamicsFunctional2D<T, Descriptor> *clone() const;

private:
    Dynamics<T, Descriptor> *dynamics;
};

/* ************* Class InstantiateComplexDomainDynamicsFunctional2D ** */

template <typename T, template <typename U> class Descriptor>
class InstantiateComplexDomainDynamicsFunctional2D :
    public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    InstantiateComplexDomainDynamicsFunctional2D(
        Dynamics<T, Descriptor> *dynamics_, DomainFunctional2D *domain_);
    InstantiateComplexDomainDynamicsFunctional2D(
        InstantiateComplexDomainDynamicsFunctional2D<T, Descriptor> const &rhs);
    InstantiateComplexDomainDynamicsFunctional2D<T, Descriptor> &operator=(
        InstantiateComplexDomainDynamicsFunctional2D<T, Descriptor> const &rhs);
    virtual ~InstantiateComplexDomainDynamicsFunctional2D();
    virtual void process(Box2D boundingBox, BlockLattice2D<T, Descriptor> &lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual InstantiateComplexDomainDynamicsFunctional2D<T, Descriptor> *clone() const;

private:
    Dynamics<T, Descriptor> *dynamics;
    DomainFunctional2D *domain;
};

/* ************* Class InstantiateDotDynamicsFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
class InstantiateDotDynamicsFunctional2D : public DotProcessingFunctional2D_L<T, Descriptor> {
public:
    InstantiateDotDynamicsFunctional2D(Dynamics<T, Descriptor> *dynamics_);
    InstantiateDotDynamicsFunctional2D(
        InstantiateDotDynamicsFunctional2D<T, Descriptor> const &rhs);
    InstantiateDotDynamicsFunctional2D<T, Descriptor> &operator=(
        InstantiateDotDynamicsFunctional2D<T, Descriptor> const &rhs);
    virtual ~InstantiateDotDynamicsFunctional2D();
    virtual void process(DotList2D const &dotList, BlockLattice2D<T, Descriptor> &lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual InstantiateDotDynamicsFunctional2D<T, Descriptor> *clone() const;

private:
    Dynamics<T, Descriptor> *dynamics;
};

/* ************* Class DynamicsFromMaskFunctional2D ************************ */

/// Assign dynamics to nodes specified by a boolean mask.
template <typename T, template <typename U> class Descriptor>
class DynamicsFromMaskFunctional2D : public BoxProcessingFunctional2D_LS<T, Descriptor, bool> {
public:
    DynamicsFromMaskFunctional2D(Dynamics<T, Descriptor> *dynamics_, bool whichFlag_);
    DynamicsFromMaskFunctional2D(DynamicsFromMaskFunctional2D<T, Descriptor> const &rhs);
    DynamicsFromMaskFunctional2D<T, Descriptor> &operator=(
        DynamicsFromMaskFunctional2D<T, Descriptor> const &rhs);
    virtual ~DynamicsFromMaskFunctional2D();
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<bool> &mask);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual DynamicsFromMaskFunctional2D<T, Descriptor> *clone() const;

private:
    Dynamics<T, Descriptor> *dynamics;
    bool whichFlag;
};

/* ************* Class DynamicsFromIntMaskFunctional2D ************************ */

/// Assign dynamics to nodes specified by an integer mask.
template <typename T, template <typename U> class Descriptor>
class DynamicsFromIntMaskFunctional2D : public BoxProcessingFunctional2D_LS<T, Descriptor, int> {
public:
    DynamicsFromIntMaskFunctional2D(Dynamics<T, Descriptor> *dynamics_, int whichFlag_);
    DynamicsFromIntMaskFunctional2D(DynamicsFromIntMaskFunctional2D<T, Descriptor> const &rhs);
    DynamicsFromIntMaskFunctional2D<T, Descriptor> &operator=(
        DynamicsFromIntMaskFunctional2D<T, Descriptor> const &rhs);
    virtual ~DynamicsFromIntMaskFunctional2D();
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<int> &mask);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual DynamicsFromIntMaskFunctional2D<T, Descriptor> *clone() const;

private:
    Dynamics<T, Descriptor> *dynamics;
    int whichFlag;
};

/* *************** Class RecomposeFromFlowVariablesFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
class RecomposeFromFlowVariablesFunctional2D : public BoxProcessingFunctional2D {
public:
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual RecomposeFromFlowVariablesFunctional2D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

/* ************* Class AssignOmegaFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
class AssignOmegaFunctional2D : public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    AssignOmegaFunctional2D(T omega_);
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual AssignOmegaFunctional2D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T omega;
};

/* ************* Class SetConstBoundaryVelocityFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
class SetConstBoundaryVelocityFunctional2D : public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    SetConstBoundaryVelocityFunctional2D(Array<T, Descriptor<T>::d> velocity);
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual SetConstBoundaryVelocityFunctional2D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Array<T, Descriptor<T>::d> u;
};

/* ************* Class SetCustomBoundaryVelocityFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor, class VelocityFunction>
class SetCustomBoundaryVelocityFunctional2D : public OneCellIndexedFunctional2D<T, Descriptor> {
public:
    SetCustomBoundaryVelocityFunctional2D(VelocityFunction f_);
    virtual SetCustomBoundaryVelocityFunctional2D<T, Descriptor, VelocityFunction> *clone() const;
    virtual void execute(plint iX, plint iY, Cell<T, Descriptor> &cell) const;
    virtual void setscale(int dxScale, int dtScale);

private:
    VelocityFunction f;
    T velocityScale;
};

/* ************* Class SetConstBoundaryDensityFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
class SetConstBoundaryDensityFunctional2D : public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    SetConstBoundaryDensityFunctional2D(T rho_);
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual SetConstBoundaryDensityFunctional2D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T rho;
};

/* ************* Class SetConstBoundaryTemperatureFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
class SetConstBoundaryTemperatureFunctional2D : public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    SetConstBoundaryTemperatureFunctional2D(T temperature_);
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual SetConstBoundaryTemperatureFunctional2D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T temperature;
};

/* ************* Class SetCustomBoundaryDensityFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor, class DensityFunction>
class SetCustomBoundaryDensityFunctional2D : public OneCellIndexedFunctional2D<T, Descriptor> {
public:
    SetCustomBoundaryDensityFunctional2D(DensityFunction f_);
    virtual void execute(plint iX, plint iY, Cell<T, Descriptor> &cell) const;
    virtual SetCustomBoundaryDensityFunctional2D<T, Descriptor, DensityFunction> *clone() const;

private:
    DensityFunction f;
};

/* ************* Class SetCustomBoundaryTemperatureFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor, class TemperatureFunction>
class SetCustomBoundaryTemperatureFunctional2D : public OneCellIndexedFunctional2D<T, Descriptor> {
public:
    SetCustomBoundaryTemperatureFunctional2D(TemperatureFunction f_);
    virtual void execute(plint iX, plint iY, Cell<T, Descriptor> &cell) const;
    virtual SetCustomBoundaryTemperatureFunctional2D<T, Descriptor, TemperatureFunction> *clone()
        const;

private:
    TemperatureFunction f;
};

/* ************* Class IniConstEquilibriumFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
class IniConstEquilibriumFunctional2D : public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    IniConstEquilibriumFunctional2D(T density, Array<T, Descriptor<T>::d> velocity, T temperature);
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual IniConstEquilibriumFunctional2D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T rho;
    T rhoBar;
    Array<T, Descriptor<T>::d> u;
    T thetaBar;
};

/* ************* Class IniCustomEquilibriumFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
class IniCustomEquilibriumFunctional2D : public OneCellIndexedFunctional2D<T, Descriptor> {
public:
    IniCustomEquilibriumFunctional2D(RhoUFunction f_);
    virtual void execute(plint iX, plint iY, Cell<T, Descriptor> &cell) const;
    virtual IniCustomEquilibriumFunctional2D<T, Descriptor, RhoUFunction> *clone() const;
    virtual void setscale(int dxScale, int dtScale);

private:
    RhoUFunction f;
    T velocityScale;
};

/* ************* Class IniConstEquilibriumComplexDomainFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
class IniConstEquilibriumComplexDomainFunctional2D :
    public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    IniConstEquilibriumComplexDomainFunctional2D(
        DomainFunctional2D *domain_, T density, Array<T, Descriptor<T>::d> velocity, T temperature);
    virtual void process(Box2D box, BlockLattice2D<T, Descriptor> &lattice);
    virtual IniConstEquilibriumComplexDomainFunctional2D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    DomainFunctional2D *domain;
    T rho;
    T rhoBar;
    Array<T, Descriptor<T>::d> u;
    T thetaBar;
};

/* ************* Class IniCustomThermalEquilibriumFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor, class RhoVelThetaFunction>
class IniCustomThermalEquilibriumFunctional2D : public OneCellIndexedFunctional2D<T, Descriptor> {
public:
    IniCustomThermalEquilibriumFunctional2D(RhoVelThetaFunction f_);
    virtual void execute(plint iX, plint iY, Cell<T, Descriptor> &cell) const;
    virtual IniCustomThermalEquilibriumFunctional2D<T, Descriptor, RhoVelThetaFunction> *clone()
        const;
    virtual void setscale(int dxScale, int dtScale);

private:
    RhoVelThetaFunction f;
    T velocityScale;
};

/* ************* Class SetCustomOmegaFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor, class OmegaFunction>
class SetCustomOmegaFunctional2D : public OneCellIndexedFunctional2D<T, Descriptor> {
public:
    SetCustomOmegaFunctional2D(OmegaFunction f_);
    virtual void execute(plint iX, plint iY, Cell<T, Descriptor> &cell) const;
    virtual SetCustomOmegaFunctional2D<T, Descriptor, OmegaFunction> *clone() const;

private:
    OmegaFunction f;
};

/* ************* Class StripeOffDensityOffsetFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
class StripeOffDensityOffsetFunctional2D : public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    StripeOffDensityOffsetFunctional2D(T deltaRho_);
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual StripeOffDensityOffsetFunctional2D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T deltaRho;
};

/* ************* Class InstantiateCompositeDynamicsFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
class InstantiateCompositeDynamicsFunctional2D : public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    InstantiateCompositeDynamicsFunctional2D(CompositeDynamics<T, Descriptor> *compositeDynamics_);
    InstantiateCompositeDynamicsFunctional2D(
        InstantiateCompositeDynamicsFunctional2D<T, Descriptor> const &rhs);
    InstantiateCompositeDynamicsFunctional2D<T, Descriptor> &operator=(
        InstantiateCompositeDynamicsFunctional2D<T, Descriptor> const &rhs);
    virtual ~InstantiateCompositeDynamicsFunctional2D();
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual InstantiateCompositeDynamicsFunctional2D<T, Descriptor> *clone() const;

private:
    CompositeDynamics<T, Descriptor> *compositeDynamics;
};

/* ************* Class SetExternalScalarFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
class SetExternalScalarFunctional2D : public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    SetExternalScalarFunctional2D(int whichScalar_, T externalScalar_);
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual SetExternalScalarFunctional2D<T, Descriptor> *clone() const;

private:
    int whichScalar;
    T externalScalar;
};

/* ************* Class SetGenericExternalScalarFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor, class Functional>
class SetGenericExternalScalarFunctional2D : public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    SetGenericExternalScalarFunctional2D(int whichScalar_, Functional const &functional_);
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual SetGenericExternalScalarFunctional2D<T, Descriptor, Functional> *clone() const;

private:
    int whichScalar;
    Functional functional;
};

/* ************* Class AddToExternalScalarFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
class AddToExternalScalarFunctional2D : public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    AddToExternalScalarFunctional2D(int whichScalar_, T externalScalar_);
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual AddToExternalScalarFunctional2D<T, Descriptor> *clone() const;

private:
    int whichScalar;
    T externalScalar;
};

/* ************* Class SetExternalScalarFromScalarFieldFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
class SetExternalScalarFromScalarFieldFunctional2D :
    public BoxProcessingFunctional2D_LS<T, Descriptor, T> {
public:
    SetExternalScalarFromScalarFieldFunctional2D(int whichScalar_);
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &field);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual SetExternalScalarFromScalarFieldFunctional2D<T, Descriptor> *clone() const;

private:
    int whichScalar;
};

/* ************* Class SetExternalVectorFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor>
class SetExternalVectorFunctional2D : public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    SetExternalVectorFunctional2D(
        int vectorStartsAt_, Array<T, Descriptor<T>::d> const &externalVector_);
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual SetExternalVectorFunctional2D<T, Descriptor> *clone() const;

private:
    int vectorStartsAt;
    Array<T, Descriptor<T>::d> externalVector;
};

/* ************* Class SetCustomExternalVectorFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor, class VectorFunction>
class SetCustomExternalVectorFunctional2D : public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    SetCustomExternalVectorFunctional2D(int vectorStartsAt_, VectorFunction f_);
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual SetCustomExternalVectorFunctional2D<T, Descriptor, VectorFunction> *clone() const;

private:
    int vectorStartsAt;
    VectorFunction f;
};

/* ************* Class SetExternalVectorFromTensorFieldFunctional2D ******************* */

template <typename T, template <typename U> class Descriptor, int nDim>
class SetExternalVectorFromTensorFieldFunctional2D :
    public BoxProcessingFunctional2D_LT<T, Descriptor, T, nDim> {
public:
    SetExternalVectorFromTensorFieldFunctional2D(int vectorStartsAt_);
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, TensorField2D<T, nDim> &tensor);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual SetExternalVectorFromTensorFieldFunctional2D<T, Descriptor, nDim> *clone() const;

private:
    int vectorStartsAt;
};

/* *************** PART II ******************************************* */
/* *************** Initialization of scalar- and tensor-fields ******* */
/* ******************************************************************* */

template <typename T>
class IniConstScalarFunctional2D : public BoxProcessingFunctional2D_S<T> {
public:
    IniConstScalarFunctional2D(T value_);
    virtual void process(Box2D domain, ScalarField2D<T> &field);
    virtual IniConstScalarFunctional2D<T> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T value;
};

template <typename T>
class MaskedIniConstScalarFunctional2D : public BoxProcessingFunctional2D_SS<T, int> {
public:
    MaskedIniConstScalarFunctional2D(int flag_, T value_);
    virtual void process(Box2D domain, ScalarField2D<T> &field, ScalarField2D<int> &mask);
    virtual MaskedIniConstScalarFunctional2D<T> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    int flag;
    T value;
};

template <typename T, class Function>
class SetToScalarFunctionFunctional2D : public BoxProcessingFunctional2D_S<T> {
public:
    SetToScalarFunctionFunctional2D(Function f_);
    virtual void process(Box2D domain, ScalarField2D<T> &field);
    virtual SetToScalarFunctionFunctional2D<T, Function> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Function f;
};

template <typename T, int nDim>
class IniConstTensorFunctional2D : public BoxProcessingFunctional2D_T<T, nDim> {
public:
    IniConstTensorFunctional2D(Array<T, nDim> const &value_);
    virtual void process(Box2D domain, TensorField2D<T, nDim> &field);
    virtual IniConstTensorFunctional2D<T, nDim> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Array<T, nDim> value;
};

template <typename T, int nDim>
class MaskedIniConstTensorFunctional2D : public BoxProcessingFunctional2D_ST<int, T, nDim> {
public:
    MaskedIniConstTensorFunctional2D(int flag_, Array<T, nDim> const &value_);
    virtual void process(Box2D domain, ScalarField2D<int> &mask, TensorField2D<T, nDim> &field);
    virtual MaskedIniConstTensorFunctional2D<T, nDim> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    int flag;
    Array<T, nDim> value;
};

template <typename T, int nDim, class Function>
class SetToTensorFunctionFunctional2D : public BoxProcessingFunctional2D_T<T, nDim> {
public:
    SetToTensorFunctionFunctional2D(Function f_);
    virtual void process(Box2D domain, TensorField2D<T, nDim> &field);
    virtual SetToTensorFunctionFunctional2D<T, nDim, Function> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Function f;
};

template <typename T>
class SetToCoordinateFunctional2D : public BoxProcessingFunctional2D_S<T> {
public:
    SetToCoordinateFunctional2D(plint index_);
    virtual void process(Box2D domain, ScalarField2D<T> &field);
    virtual SetToCoordinateFunctional2D<T> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    plint index;
};

template <typename T>
class SetToRandomFunctional2D : public BoxProcessingFunctional2D_S<T> {
public:
    SetToRandomFunctional2D(Box2D const &boundingBox_, uint32_t const *seed_);
    virtual void process(Box2D domain, ScalarField2D<T> &field);
    virtual SetToRandomFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    plint nY;
    uint32_t const *seed;
};

template <typename T>
class SetToCoordinatesFunctional2D : public BoxProcessingFunctional2D_T<T, 2> {
public:
    SetToCoordinatesFunctional2D();
    virtual void process(Box2D domain, TensorField2D<T, 2> &field);
    virtual SetToCoordinatesFunctional2D<T> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

template <typename T, int nDim>
class SetTensorComponentFunctional2D : public BoxProcessingFunctional2D_ST<T, T, nDim> {
public:
    SetTensorComponentFunctional2D(int whichDim_);
    virtual void process(
        Box2D domain, ScalarField2D<T> &scalarField, TensorField2D<T, nDim> &tensorField);
    virtual SetTensorComponentFunctional2D<T, nDim> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    int whichDim;
};

template <typename T>
class GrowDomainFunctional2D : public BoxProcessingFunctional2D_S<T> {
public:
    GrowDomainFunctional2D(T flag_);
    virtual void process(Box2D domain, ScalarField2D<T> &voxels);
    virtual GrowDomainFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    T flag;
};

}  // namespace plb

#endif  // DATA_INITIALIZER_FUNCTIONAL_2D_H

// Explicitly include generic algorithms which are never precompiled (not even in precompiled
// version)
#include "dataProcessors/dataInitializerGenerics2D.h"
