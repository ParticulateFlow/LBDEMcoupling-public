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

#ifndef ADVECTION_DIFFUSION_3D_H
#define ADVECTION_DIFFUSION_3D_H

#include <memory>

#include "atomicBlock/dataProcessor3D.h"
#include "core/block3D.h"
#include "core/globalDefs.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"

namespace plb {

/**
 * Multiphysics class for one-way coupling between Navier-Stokes and
 * advection-diffusion equations: the fluid velocity is copied
 * to the advection-diffusion field, which is advected passively.
 */
template <typename T, template <typename U> class TemperatureDescriptor>
class VelocityToPassiveAdvDiff3D :
    public BoxProcessingFunctional3D_LT<T, TemperatureDescriptor, T, 3> {
public:
    virtual void process(
        Box3D domain, BlockLattice3D<T, TemperatureDescriptor> &temperature,
        TensorField3D<T, 3> &velocity);
    virtual VelocityToPassiveAdvDiff3D<T, TemperatureDescriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

template <typename T, template <typename U> class TemperatureDescriptor>
void velocityToPassiveAdvDiff(
    MultiBlockLattice3D<T, TemperatureDescriptor> &temperature, MultiTensorField3D<T, 3> &velocity,
    Box3D domain);

template <typename T, template <typename U> class TemperatureDescriptor>
class N_VelocityToPassiveAdvDiff3D :
    public BoxProcessingFunctional3D_LN<T, TemperatureDescriptor, T> {
public:
    virtual void process(
        Box3D domain, BlockLattice3D<T, TemperatureDescriptor> &temperature,
        NTensorField3D<T> &velocity);
    virtual N_VelocityToPassiveAdvDiff3D<T, TemperatureDescriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

template <typename T, template <typename U> class TemperatureDescriptor>
void NVelocityToPassiveAdvDiff(
    MultiBlockLattice3D<T, TemperatureDescriptor> &temperature, MultiNTensorField3D<T> &velocity,
    Box3D domain);

/**
 * Multiphysics class for one-way coupling between Navier-Stokes and
 * advection-diffusion equations: the fluid velocity is copied
 * to the advection-diffusion field, which is advected passively.
 */
template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
class LatticeToPassiveAdvDiff3D :
    public BoxProcessingFunctional3D_LL<T, FluidDescriptor, T, ScalarDescriptor> {
public:
    LatticeToPassiveAdvDiff3D(T scaling_ = 1.);
    virtual void process(
        Box3D domain, BlockLattice3D<T, FluidDescriptor> &fluid,
        BlockLattice3D<T, ScalarDescriptor> &scalar);
    virtual LatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T scaling;
};

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
void latticeToPassiveAdvDiff(
    MultiBlockLattice3D<T, FluidDescriptor> &fluid,
    MultiBlockLattice3D<T, ScalarDescriptor> &scalar, Box3D domain);

/**
 * Multiphysics class for one-way coupling between Navier-Stokes and
 * advection-diffusion equations in the complete regularized advection diffusion case:
 * the fluid density, velocity and stress is copied
 * to the advection-diffusion field, which is advected passively.
 */
template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
class LatticeToPassiveComplRegAdvDiff3D :
    public BoxProcessingFunctional3D_LL<T, FluidDescriptor, T, ScalarDescriptor> {
public:
    LatticeToPassiveComplRegAdvDiff3D(T scaling_ = 1.);
    virtual void process(
        Box3D domain, BlockLattice3D<T, FluidDescriptor> &fluid,
        BlockLattice3D<T, ScalarDescriptor> &scalar);
    virtual LatticeToPassiveComplRegAdvDiff3D<T, FluidDescriptor, ScalarDescriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T scaling;
};

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
void latticeToPassiveComplRegAdvDiff(
    MultiBlockLattice3D<T, FluidDescriptor> &fluid,
    MultiBlockLattice3D<T, ScalarDescriptor> &scalar, Box3D domain);

/**
 * Multiphysics class for one-way coupling between Navier-Stokes and
 * advection-diffusion equations: the fluid velocity is copied
 * to the advection-diffusion field, which is advected passively.
 * The relaxation parameter of the advection-diffusion equation is
 * also changed according to an eddy-diffusivity-type model by using
 * the provided value of the turbulent Prandtl number.
 */
template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
class TurbulentLatticeToPassiveAdvDiff3D :
    public BoxProcessingFunctional3D_LL<T, FluidDescriptor, T, ScalarDescriptor> {
public:
    // Pr_t is the turbulent Prandtl number (the ratio between the eddy kinematic viscosity and the
    // eddy diffusivity).
    TurbulentLatticeToPassiveAdvDiff3D(T Pr_t);
    virtual void process(
        Box3D domain, BlockLattice3D<T, FluidDescriptor> &fluid,
        BlockLattice3D<T, ScalarDescriptor> &scalar);
    virtual TurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T C;
};

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
void turbulentLatticeToPassiveAdvDiff(
    MultiBlockLattice3D<T, FluidDescriptor> &fluid,
    MultiBlockLattice3D<T, ScalarDescriptor> &scalar, T Pr_t, Box3D domain);

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
class MaskedTurbulentLatticeToPassiveAdvDiff3D : public BoxProcessingFunctional3D {
public:
    // Pr_t is the turbulent Prandtl number (the ratio between the eddy kinematic viscosity and the
    // eddy diffusivity).
    MaskedTurbulentLatticeToPassiveAdvDiff3D(T Pr_t, int whichFlag_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual MaskedTurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor> *clone()
        const;
    virtual BlockDomain::DomainT appliesTo() const;
    void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T C;
    int whichFlag;
};

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
void maskedTurbulentLatticeToPassiveAdvDiff(
    MultiBlockLattice3D<T, FluidDescriptor> &fluid,
    MultiBlockLattice3D<T, ScalarDescriptor> &scalar, MultiScalarField3D<int> &mask, T Pr_t,
    int whichFlag, Box3D domain);

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
class N_MaskedTurbulentLatticeToPassiveAdvDiff3D : public BoxProcessingFunctional3D {
public:
    // Pr_t is the turbulent Prandtl number (the ratio between the eddy kinematic viscosity and the
    // eddy diffusivity).
    N_MaskedTurbulentLatticeToPassiveAdvDiff3D(T Pr_t, int whichFlag_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual N_MaskedTurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>
        *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T C;
    int whichFlag;
};

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
void NMaskedTurbulentLatticeToPassiveAdvDiff(
    MultiBlockLattice3D<T, FluidDescriptor> &fluid,
    MultiBlockLattice3D<T, ScalarDescriptor> &scalar, MultiNTensorField3D<int> &mask, T Pr_t,
    int whichFlag, Box3D domain);

template <typename T, template <typename U> class Descriptor>
class CrystallizeAndAggregate : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    CrystallizeAndAggregate(T Ncr_, T Nag_);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual CrystallizeAndAggregate<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T Ncr, Nag;
};

template <typename T, template <typename U> class Descriptor>
void crystallizeAndAggregate(
    MultiBlockLattice3D<T, Descriptor> &lattice, T Ncr, T Nag, Box3D domain);

/* Coupling between the free-surface model and an advection-diffusion equation.
 * The "empty space" of the free-surface simulation is approximated by an adiabatic
 * medium in the advection-diffusion simulation.
 */
template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
class AdvectionDiffusionFreeSurfaceCoupling3D : public BoxProcessingFunctional3D {
public:
    // Pr_t is the turbulent Prandtl number (the ratio between the eddy kinematic viscosity and the
    // eddy diffusivity).
    AdvectionDiffusionFreeSurfaceCoupling3D(
        Dynamics<T, AD_Descriptor> *conductiveDynamics_,
        Dynamics<T, AD_Descriptor> *adiabaticDynamics_, T iniVal_, T Pr_t);
    ~AdvectionDiffusionFreeSurfaceCoupling3D();
    AdvectionDiffusionFreeSurfaceCoupling3D(
        AdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> const &rhs);
    AdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> &operator=(
        AdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> const &rhs);
    void swap(AdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> &rhs);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual AdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Dynamics<T, AD_Descriptor> *conductiveDynamics;
    Dynamics<T, AD_Descriptor> *adiabaticDynamics;
    T iniVal;
    T C;
};

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
class MaskedAdvectionDiffusionFreeSurfaceCoupling3D : public BoxProcessingFunctional3D {
public:
    // Pr_t is the turbulent Prandtl number (the ratio between the eddy kinematic viscosity and the
    // eddy diffusivity).
    MaskedAdvectionDiffusionFreeSurfaceCoupling3D(
        Dynamics<T, AD_Descriptor> *conductiveDynamics_,
        Dynamics<T, AD_Descriptor> *adiabaticDynamics_, T iniVal_, T Pr_t, int maskValue_);
    ~MaskedAdvectionDiffusionFreeSurfaceCoupling3D();
    MaskedAdvectionDiffusionFreeSurfaceCoupling3D(
        MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> const &rhs);
    MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> &operator=(
        MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> const &rhs);
    void swap(MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> &rhs);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> *clone()
        const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Dynamics<T, AD_Descriptor> *conductiveDynamics;
    Dynamics<T, AD_Descriptor> *adiabaticDynamics;
    T iniVal;
    T C;
    int maskValue;
};

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
class N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D : public BoxProcessingFunctional3D {
public:
    // Pr_t is the turbulent Prandtl number (the ratio between the eddy kinematic viscosity and the
    // eddy diffusivity).
    N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D(
        Dynamics<T, AD_Descriptor> *conductiveDynamics_,
        Dynamics<T, AD_Descriptor> *adiabaticDynamics_, T iniVal_, T Pr_t, int maskValue_);
    ~N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D();
    N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D(
        N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> const
            &rhs);
    N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> &operator=(
        N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> const
            &rhs);
    void swap(
        N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> &rhs);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>
        *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Dynamics<T, AD_Descriptor> *conductiveDynamics;
    Dynamics<T, AD_Descriptor> *adiabaticDynamics;
    T iniVal;
    T C;
    int maskValue;
};

template <typename T>
class AdvectionDiffusionFd3D : public BoxProcessingFunctional3D {
public:
    AdvectionDiffusionFd3D(T d_, bool upwind_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual AdvectionDiffusionFd3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T d;
    bool upwind;
};

template <typename T>
class VariableDiffusivityAdvectionDiffusionFd3D : public BoxProcessingFunctional3D {
public:
    VariableDiffusivityAdvectionDiffusionFd3D(bool upwind_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual VariableDiffusivityAdvectionDiffusionFd3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    bool upwind;
};

// The class ADCellT describes an "advection-diffusion cell type". It needs to implement
// four functions:
//
// bool isSingular(int flag);
// bool isFixed(int flag, T& phi);
// bool isZeroFlux(int flag);
// bool isFixedFlux(int flag, T phi, T& dphi);
//
// where "flag" is an integer number which characterizes a cell, "phi" is the dependent
// variable on the specific cell, and "dphi" is the normal derivative of the dependent
// variable on the specific cell (the normal is the outward unit normal). If "isSingular"
// returns true on a cell, then on this cell the value "defaultValue" is simply imposed
// on "phi" for plotting purposes and it is not used for any meaningful computation
// inside the data processor. If a singular cell is next to a normal cell (where the AD
// equation needs to be solved), then a "zero flux" condition is implemented. The "phi"
// and "dphi" values can be used by the caller only if the corresponding functions return
// true, otherwise they are undefined. The "isFixedFlux" function takes as an argument
// the value "phi" and returns the value of its normal derivative (this is used to
// implement boundary conditions of Robin type among others).
// TODO: This implementation uses the ADCellT to execute a lot of tests, so it is a bit
//       inefficient. Optimization is needed.
template <typename T, class ADCellT>
class MaskedVariableDiffusivityAdvectionDiffusionFd3D : public BoxProcessingFunctional3D {
public:
    MaskedVariableDiffusivityAdvectionDiffusionFd3D(
        ADCellT adCellT_, T defaultValue_, Box3D const &fullClosedDomain_, bool xPeriodic_,
        bool yPeriodic_, bool zPeriodic_, bool upwind_ = true);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual MaskedVariableDiffusivityAdvectionDiffusionFd3D<T, ADCellT> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    ADCellT adCellT;
    T defaultValue;
    Box3D fullExtendedDomain;
    bool upwind;
};

template <typename T, template <typename U> class FluidDescriptor>
class TurbulentDiffusivityFunctional3D :
    public BoxProcessingFunctional3D_LS<T, FluidDescriptor, T> {
public:
    // D0 is the base diffusivity.
    // Pr_t is the turbulent Prandtl number (the ratio between the eddy kinematic viscosity and the
    // eddy diffusivity).
    TurbulentDiffusivityFunctional3D(T D0_, T Pr_t);
    virtual void process(
        Box3D domain, BlockLattice3D<T, FluidDescriptor> &fluid, ScalarField3D<T> &diffusivity);
    virtual TurbulentDiffusivityFunctional3D<T, FluidDescriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T D0, C;
};

template <typename T, template <typename U> class FluidDescriptor>
void computeTurbulentDiffusivity(
    T D0, T Pr_t, MultiBlockLattice3D<T, FluidDescriptor> &fluid,
    MultiScalarField3D<T> &diffusivity, Box3D domain);

template <typename T, template <typename U> class FluidDescriptor>
std::unique_ptr<MultiScalarField3D<T> > computeTurbulentDiffusivity(
    T D0, T Pr_t, MultiBlockLattice3D<T, FluidDescriptor> &fluid, Box3D domain);

template <typename T, template <typename U> class FluidDescriptor>
std::unique_ptr<MultiScalarField3D<T> > computeTurbulentDiffusivity(
    T D0, T Pr_t, MultiBlockLattice3D<T, FluidDescriptor> &fluid);

template <typename T, template <typename U> class FluidDescriptor>
class TurbulentDiffusivityFromFlagFunctional3D : public BoxProcessingFunctional3D {
public:
    // D0 is the vector of base diffusivities. Its values are in one-to-one correspondence with the
    // flags. Pr_t is the turbulent Prandtl number (the ratio between the eddy kinematic viscosity
    // and the eddy diffusivity).
    TurbulentDiffusivityFromFlagFunctional3D(std::vector<T> const &D0_, T Pr_t);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual TurbulentDiffusivityFromFlagFunctional3D<T, FluidDescriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    std::vector<T> D0;
    T C;
};

template <typename T, template <typename U> class FluidDescriptor>
void computeTurbulentDiffusivityFromFlag(
    std::vector<T> const &D0, T Pr_t, MultiBlockLattice3D<T, FluidDescriptor> &fluid,
    MultiScalarField3D<int> &flags, MultiScalarField3D<T> &diffusivity, Box3D domain);

template <typename T, template <typename U> class FluidDescriptor>
std::unique_ptr<MultiScalarField3D<T> > computeTurbulentDiffusivityFromFlag(
    std::vector<T> const &D0, T Pr_t, MultiBlockLattice3D<T, FluidDescriptor> &fluid,
    MultiScalarField3D<int> &flags, Box3D domain);

template <typename T, template <typename U> class FluidDescriptor>
std::unique_ptr<MultiScalarField3D<T> > computeTurbulentDiffusivityFromFlag(
    std::vector<T> const &D0, T Pr_t, MultiBlockLattice3D<T, FluidDescriptor> &fluid,
    MultiScalarField3D<int> &flags);

template <typename T, template <typename U> class FluidDescriptor>
class TurbulentDiffusivityFromScalarFunctional3D : public BoxProcessingFunctional3D {
public:
    // Pr_t is the turbulent Prandtl number (the ratio between the eddy kinematic viscosity and the
    // eddy diffusivity).
    TurbulentDiffusivityFromScalarFunctional3D(T Pr_t);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual TurbulentDiffusivityFromScalarFunctional3D<T, FluidDescriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T C;
};

template <typename T, template <typename U> class FluidDescriptor>
void computeTurbulentDiffusivityFromScalar(
    T Pr_t, MultiBlockLattice3D<T, FluidDescriptor> &fluid, MultiScalarField3D<T> &D0,
    MultiScalarField3D<T> &diffusivity, Box3D domain);

template <typename T, template <typename U> class FluidDescriptor>
std::unique_ptr<MultiScalarField3D<T> > computeTurbulentDiffusivityFromScalar(
    T Pr_t, MultiBlockLattice3D<T, FluidDescriptor> &fluid, MultiScalarField3D<T> &D0,
    Box3D domain);

template <typename T, template <typename U> class FluidDescriptor>
std::unique_ptr<MultiScalarField3D<T> > computeTurbulentDiffusivityFromScalar(
    T Pr_t, MultiBlockLattice3D<T, FluidDescriptor> &fluid, MultiScalarField3D<T> &D0);

// This data processor implements a chemical reaction coupling based on the Arrhenius equation.
// It couples a number of species (nSpecies) according to a number of chemical reactions
// (nReactions). It is very important to note that all data provided to this data processor are in
// physical units (as opposed to lattice units). This is because the units of the reaction rate
// contstant depend on the corresponding global order of reaction, so it is more convenient to do
// all calculations in physical units, compute the source terms of the involved advection-diffusion
// equations and then scale back to lattice units only these source terms. The units of the
// different quantities are as follows:
//
// preExponentialFactors: unspecified (depend on the corresponding global order of reaction)
// activationEnergies   : J / mol
// reactionEnthalpies   : J / mol
//
// The time step (dt) is provided in seconds, in order to transform the source terms from physical
// to lattice units.
//
// The stoichiometry matrix defines the stoichiometries for each chemical reaction.
// It has a size: nSpecies x nReactions.
//
// The reaction rate matrix defines the partial orders of reaction rates for each chemical reaction.
// It has a size: nReactions x nSpecies.
//
// As usual, the user needs to provide a set of blocks to this data processor. In the following, the
// units of the quantities provided as blocks are listed:
//
// mass concentrations (for each species): mol / m^3
// temperature                           : K
// volumetric heat capacity              : J / (m^3 * K)
//
// The volumetric heat capacity is the product of the density (kg / m^3) times the specific heat
// capacity (J / (kg * K)). The source terms of the advection-diffusion equations (which are
// provided as blocks and are modified by this data processor) are in lattice units.
template <typename T>
class ArrheniusChemicalReactionCoupling3D : public BoxProcessingFunctional3D {
public:
    ArrheniusChemicalReactionCoupling3D(
        plint nSpecies_, plint nReactions_, T dt_, std::vector<T> const &preExponentialFactors_,
        std::vector<T> const &activationEnergies_, std::vector<T> const &reactionEnthalpies_,
        std::vector<std::vector<int> > const &stoichiometryMatrix_,
        std::vector<std::vector<T> > const &reactionRateMatrix_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual ArrheniusChemicalReactionCoupling3D<T> *clone() const;
    void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    plint nSpecies;
    plint nReactions;
    T dt;
    std::vector<T> preExponentialFactors;
    std::vector<T> activationEnergies;
    std::vector<T> reactionEnthalpies;
    std::vector<std::vector<int> > stoichiometryMatrix;
    std::vector<std::vector<T> > reactionRateMatrix;
    static T R;  // Universal gas constant [J / (mol * K)].
};

}  // namespace plb

#endif  // ADVECTION_DIFFUSION_3D_H
