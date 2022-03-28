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
 * Helper functions for domain initialization -- generic implementation.
 */
#ifndef DATA_INITIALIZER_WRAPPER_3D_HH
#define DATA_INITIALIZER_WRAPPER_3D_HH

#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "core/cell.h"
#include "dataProcessors/dataInitializerWrapper3D.h"
#include "dataProcessors/metaStuffWrapper3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Initialization of the block-lattice: atomic-block * */
/* ******************************************************************* */

template <typename T, template <class U> class Descriptor>
void apply(
    BlockLattice3D<T, Descriptor> &lattice, Box3D domain, OneCellFunctional3D<T, Descriptor> *f)
{
    applyProcessingFunctional(new GenericLatticeFunctional3D<T, Descriptor>(f), domain, lattice);
}

template <typename T, template <class U> class Descriptor>
void applyIndexed(
    BlockLattice3D<T, Descriptor> &lattice, Box3D domain,
    OneCellIndexedFunctional3D<T, Descriptor> *f)
{
    applyProcessingFunctional(
        new GenericIndexedLatticeFunctional3D<T, Descriptor>(f), domain, lattice);
}

template <typename T, template <class U> class Descriptor>
void defineDynamics(
    BlockLattice3D<T, Descriptor> &lattice, Box3D domain, Dynamics<T, Descriptor> *dynamics)
{
    applyProcessingFunctional(
        new InstantiateDynamicsFunctional3D<T, Descriptor>(dynamics), domain, lattice);
}

template <typename T, template <class U> class Descriptor>
void defineDynamics(
    BlockLattice3D<T, Descriptor> &lattice, Box3D boundingBox, DomainFunctional3D *domain,
    Dynamics<T, Descriptor> *dynamics)
{
    applyProcessingFunctional(
        new InstantiateComplexDomainDynamicsFunctional3D<T, Descriptor>(dynamics, domain),
        boundingBox, lattice);
}

template <typename T, template <class U> class Descriptor>
void defineDynamics(
    BlockLattice3D<T, Descriptor> &lattice, DotList3D const &dotList,
    Dynamics<T, Descriptor> *dynamics)
{
    applyProcessingFunctional(
        new InstantiateDotDynamicsFunctional3D<T, Descriptor>(dynamics), dotList, lattice);
}

template <typename T, template <class U> class Descriptor>
void defineDynamics(
    BlockLattice3D<T, Descriptor> &lattice, plint iX, plint iY, plint iZ,
    Dynamics<T, Descriptor> *dynamics)
{
    DotList3D pos;
    pos.addDot(Dot3D(iX, iY, iZ));
    defineDynamics(lattice, pos, dynamics);
}

template <typename T, template <typename U> class Descriptor>
void defineDynamics(
    BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<bool> &boolMask, Box3D domain,
    Dynamics<T, Descriptor> *dynamics, bool whichFlag)
{
    applyProcessingFunctional(
        new DynamicsFromMaskFunctional3D<T, Descriptor>(dynamics, whichFlag), domain, lattice,
        boolMask);
}

template <typename T, template <typename U> class Descriptor>
void defineDynamics(
    BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<bool> &boolMask,
    Dynamics<T, Descriptor> *dynamics, bool whichFlag)
{
    defineDynamics(lattice, boolMask, lattice.getBoundingBox(), dynamics, whichFlag);
}

template <typename T, template <typename U> class Descriptor>
void defineDynamics(
    BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &intMask, Box3D domain,
    Dynamics<T, Descriptor> *dynamics, int whichFlag)
{
    applyProcessingFunctional(
        new DynamicsFromIntMaskFunctional3D<T, Descriptor>(dynamics, whichFlag), domain, lattice,
        intMask);
}

template <typename T, template <typename U> class Descriptor>
void defineDynamics(
    BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &intMask,
    Dynamics<T, Descriptor> *dynamics, int whichFlag)
{
    defineDynamics(lattice, intMask, lattice.getBoundingBox(), dynamics, whichFlag);
}

template <typename T, template <typename U> class Descriptor>
void recomposeFromFlowVariables(
    BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &density,
    TensorField3D<T, 3> &velocity, TensorField3D<T, 6> &strainRate, Box3D domain)
{
    std::vector<AtomicBlock3D *> atomicBlocks(4);
    atomicBlocks[0] = &lattice;
    atomicBlocks[1] = &density;
    atomicBlocks[2] = &velocity;
    atomicBlocks[3] = &strainRate;
    applyProcessingFunctional(
        new RecomposeFromFlowVariablesFunctional3D<T, Descriptor>, domain, atomicBlocks);
}

template <typename T, template <typename U> class Descriptor>
void recomposeFromFlowVariables(
    BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &density,
    TensorField3D<T, 3> &velocity, TensorField3D<T, 6> &strainRate)
{
    recomposeFromFlowVariables(lattice, density, velocity, strainRate, lattice.getBoundingBox());
}

template <typename T, template <class U> class Descriptor>
void setOmega(BlockLattice3D<T, Descriptor> &lattice, Box3D domain, T omega)
{
    applyProcessingFunctional(new AssignOmegaFunctional3D<T, Descriptor>(omega), domain, lattice);
}

template <typename T, template <class U> class Descriptor>
void setOmega(BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &omega, Box3D domain)
{
    applyProcessingFunctional(
        new AssignScalarFieldOmegaFunctional3D<T, Descriptor>(), domain, lattice, omega);
}

template <typename T, template <class U> class Descriptor>
void setBoundaryVelocity(BlockLattice3D<T, Descriptor> &lattice, Box3D domain, Array<T, 3> velocity)
{
    applyProcessingFunctional(
        new SetConstBoundaryVelocityFunctional3D<T, Descriptor>(velocity), domain, lattice);
}

template <typename T, template <class U> class Descriptor>
void setBoundaryVelocity(
    BlockLattice3D<T, Descriptor> &lattice, TensorField3D<T, Descriptor<T>::d> &force, Box3D domain,
    Array<T, Descriptor<T>::d> velocity)
{
    applyProcessingFunctional(
        new SetConstBoundaryVelocityWithTensorForceFunctional3D<T, Descriptor>(velocity), domain,
        lattice, force);
}

template <typename T, template <class U> class Descriptor>
void setBoundaryVelocity(
    BlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force, Box3D domain,
    Array<T, Descriptor<T>::d> velocity)
{
    applyProcessingFunctional(
        new SetConstBoundaryVelocityWithForceFunctional3D<T, Descriptor>(force, velocity), domain,
        lattice);
}

template <typename T, template <class U> class Descriptor, class VelocityFunction>
void setBoundaryVelocity(
    BlockLattice3D<T, Descriptor> &lattice, TensorField3D<T, Descriptor<T>::d> &force, Box3D domain,
    VelocityFunction f)
{
    applyProcessingFunctional(
        new SetCustomBoundaryVelocityWithTensorForceFunctional3D<T, Descriptor, VelocityFunction>(
            f),
        domain, lattice, force);
}

template <typename T, template <class U> class Descriptor, class VelocityFunction>
void setBoundaryVelocity(
    BlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force, Box3D domain,
    VelocityFunction f)
{
    applyProcessingFunctional(
        new SetCustomBoundaryVelocityWithForceFunctional3D<T, Descriptor, VelocityFunction>(
            force, f),
        domain, lattice);
}

template <typename T, template <class U> class Descriptor>
void setBoundaryDensity(BlockLattice3D<T, Descriptor> &lattice, Box3D domain, T rho)
{
    applyProcessingFunctional(
        new SetConstBoundaryDensityFunctional3D<T, Descriptor>(rho), domain, lattice);
}

template <typename T, template <class U> class Descriptor>
void initializeAtEquilibrium(
    BlockLattice3D<T, Descriptor> &lattice, Box3D domain, T rho, Array<T, 3> velocity,
    T temperature)
{
    applyProcessingFunctional(
        new IniConstEquilibriumFunctional3D<T, Descriptor>(rho, velocity, temperature), domain,
        lattice);
}

template <typename T, template <class U> class Descriptor>
void stripeOffDensityOffset(BlockLattice3D<T, Descriptor> &lattice, Box3D domain, T deltaRho)
{
    applyProcessingFunctional(
        new StripeOffDensityOffsetFunctional3D<T, Descriptor>(deltaRho), domain, lattice);
}

template <typename T, template <class U> class Descriptor>
void setCompositeDynamics(
    BlockLattice3D<T, Descriptor> &lattice, Box3D domain,
    CompositeDynamics<T, Descriptor> *compositeDynamics)
{
    applyProcessingFunctional(
        new InstantiateCompositeDynamicsFunctional3D<T, Descriptor>(compositeDynamics), domain,
        lattice);
}

template <typename T, template <class U> class Descriptor>
void setExternalScalar(
    BlockLattice3D<T, Descriptor> &lattice, Box3D domain, int whichScalar, T externalScalar)
{
    applyProcessingFunctional(
        new SetExternalScalarFunctional3D<T, Descriptor>(whichScalar, externalScalar), domain,
        lattice);
}

template <typename T, template <class U> class Descriptor>
void setExternalScalar(
    BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &mask, int flag, Box3D domain,
    int whichScalar, T externalScalar)
{
    applyProcessingFunctional(
        new MaskedSetExternalScalarFunctional3D<T, Descriptor>(flag, whichScalar, externalScalar),
        domain, lattice, mask);
}

template <typename T, template <class U> class Descriptor>
void setExternalVector(
    BlockLattice3D<T, Descriptor> &lattice, Box3D domain, int vectorStartsAt,
    Array<T, Descriptor<T>::d> externalVector)
{
    applyProcessingFunctional(
        new SetExternalVectorFunctional3D<T, Descriptor>(vectorStartsAt, externalVector), domain,
        lattice);
}

/* *************** PART II ******************************************* */
/* *************** Initialization of the block-lattice: multi-block * */
/* ******************************************************************* */

template <typename T, template <class U> class Descriptor>
void apply(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain,
    OneCellFunctional3D<T, Descriptor> *f)
{
    applyProcessingFunctional(new GenericLatticeFunctional3D<T, Descriptor>(f), domain, lattice);
}

template <typename T, template <class U> class Descriptor>
void applyIndexed(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain,
    OneCellIndexedFunctional3D<T, Descriptor> *f)
{
    applyProcessingFunctional(
        new GenericIndexedLatticeFunctional3D<T, Descriptor>(f), domain, lattice);
}

template <typename T, template <class U> class Descriptor>
void applyIndexed(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain,
    OneCellIndexedWithRandFunctional3D<T, Descriptor> *f, uint32_t seed)
{
    applyProcessingFunctional(
        new GenericIndexedWithRandLatticeFunctional3D<T, Descriptor>(
            f, lattice.getBoundingBox(), &seed),
        domain, lattice);
}

template <typename T, template <class U> class Descriptor>
void defineDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, Dynamics<T, Descriptor> *dynamics)
{
    applyProcessingFunctional(
        new InstantiateDynamicsFunctional3D<T, Descriptor>(dynamics), domain, lattice);
}

template <typename T, template <class U> class Descriptor>
void defineDynamicsInBulkAndEnvelope(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, Dynamics<T, Descriptor> *dynamics)
{
    applyProcessingFunctional(
        new InstantiateDynamicsInBulkAndEnvelopeFunctional3D<T, Descriptor>(dynamics), domain,
        lattice);
}

template <typename T, template <class U> class Descriptor>
void defineDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D boundingBox, DomainFunctional3D *domain,
    Dynamics<T, Descriptor> *dynamics)
{
    applyProcessingFunctional(
        new InstantiateComplexDomainDynamicsFunctional3D<T, Descriptor>(dynamics, domain),
        boundingBox, lattice);
}

template <typename T, template <class U> class Descriptor>
void defineDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, DotList3D const &dotList,
    Dynamics<T, Descriptor> *dynamics)
{
    applyProcessingFunctional(
        new InstantiateDotDynamicsFunctional3D<T, Descriptor>(dynamics), dotList, lattice);
}

template <typename T, template <class U> class Descriptor>
void defineDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, plint iX, plint iY, plint iZ,
    Dynamics<T, Descriptor> *dynamics)
{
    DotList3D pos;
    pos.addDot(Dot3D(iX, iY, iZ));
    defineDynamics(lattice, pos, dynamics);
}

template <typename T, template <typename U> class Descriptor>
void defineDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<bool> &boolMask, Box3D domain,
    Dynamics<T, Descriptor> *dynamics, bool whichFlag)
{
    applyProcessingFunctional(
        new DynamicsFromMaskFunctional3D<T, Descriptor>(dynamics, whichFlag), domain, lattice,
        boolMask);
}

template <typename T, template <typename U> class Descriptor>
void defineDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<bool> &boolMask,
    Dynamics<T, Descriptor> *dynamics, bool whichFlag)
{
    defineDynamics(lattice, boolMask, lattice.getBoundingBox(), dynamics, whichFlag);
}

template <typename T, template <typename U> class Descriptor>
void defineDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<int> &intMask, Box3D domain,
    Dynamics<T, Descriptor> *dynamics, int whichFlag)
{
    applyProcessingFunctional(
        new DynamicsFromIntMaskFunctional3D<T, Descriptor>(dynamics, whichFlag), domain, lattice,
        intMask);
}

template <typename T, template <typename U> class Descriptor>
void defineDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<int> &intMask,
    Dynamics<T, Descriptor> *dynamics, int whichFlag)
{
    defineDynamics(lattice, intMask, lattice.getBoundingBox(), dynamics, whichFlag);
}

template <typename T, template <typename U> class Descriptor>
void recomposeFromFlowVariables(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &density,
    MultiTensorField3D<T, 3> &velocity, MultiTensorField3D<T, 6> &strainRate, Box3D domain)
{
    std::vector<MultiBlock3D *> multiBlocks(4);
    multiBlocks[0] = &lattice;
    multiBlocks[1] = &density;
    multiBlocks[2] = &velocity;
    multiBlocks[3] = &strainRate;
    applyProcessingFunctional(
        new RecomposeFromFlowVariablesFunctional3D<T, Descriptor>, domain, multiBlocks);
}

template <typename T, template <typename U> class Descriptor>
void recomposeFromOrderZeroVariables(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &density,
    MultiTensorField3D<T, 3> &velocity, MultiTensorField3D<T, Descriptor<T>::q> &fNeq)
{
    recomposeFromOrderZeroVariables(lattice, density, velocity, fNeq, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void recomposeFromOrderZeroVariables(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &density,
    MultiTensorField3D<T, 3> &velocity, MultiTensorField3D<T, Descriptor<T>::q> &fNeq, Box3D domain)
{
    std::vector<MultiBlock3D *> multiBlocks(4);
    multiBlocks[0] = &lattice;
    multiBlocks[1] = &density;
    multiBlocks[2] = &velocity;
    multiBlocks[3] = &fNeq;
    applyProcessingFunctional(
        new RecomposeFromOrderZeroVariablesFunctional3D<T, Descriptor>, domain, multiBlocks);
}

template <typename T, template <typename U> class Descriptor>
void recomposeFromFlowVariables(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &density,
    MultiTensorField3D<T, 3> &velocity, MultiTensorField3D<T, 6> &strainRate)
{
    recomposeFromFlowVariables(lattice, density, velocity, strainRate, lattice.getBoundingBox());
}

template <typename T, template <class U> class Descriptor>
void setOmega(MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, T omega)
{
    applyProcessingFunctional(new AssignOmegaFunctional3D<T, Descriptor>(omega), domain, lattice);
}

template <typename T, template <class U> class Descriptor, class Function>
void setOmega(MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, Function f)
{
    applyIndexed(lattice, domain, new SetCustomOmegaFunctional3D<T, Descriptor, Function>(f));
}

template <typename T, template <class U> class Descriptor>
void setOmega(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &omega, Box3D domain)
{
    applyProcessingFunctional(
        new AssignScalarFieldOmegaFunctional3D<T, Descriptor>(), domain, lattice, omega);
}

template <typename T, template <class U> class Descriptor>
void setBoundaryVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, Array<T, 3> velocity)
{
    applyProcessingFunctional(
        new SetConstBoundaryVelocityFunctional3D<T, Descriptor>(velocity), domain, lattice);
}

template <typename T, template <class U> class Descriptor>
void setBoundaryVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    Box3D domain, Array<T, Descriptor<T>::d> velocity)
{
    applyProcessingFunctional(
        new SetConstBoundaryVelocityWithTensorForceFunctional3D<T, Descriptor>(velocity), domain,
        lattice, force);
}

template <typename T, template <class U> class Descriptor>
void setBoundaryVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force, Box3D domain,
    Array<T, Descriptor<T>::d> velocity)
{
    applyProcessingFunctional(
        new SetConstBoundaryVelocityWithForceFunctional3D<T, Descriptor>(force, velocity), domain,
        lattice);
}

template <typename T, template <class U> class Descriptor, class VelocityFunction>
void setBoundaryVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    Box3D domain, VelocityFunction f)
{
    applyProcessingFunctional(
        new SetCustomBoundaryVelocityWithTensorForceFunctional3D<T, Descriptor, VelocityFunction>(
            f),
        domain, lattice, force);
}

template <typename T, template <class U> class Descriptor, class VelocityFunction>
void setBoundaryVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force, Box3D domain,
    VelocityFunction f)
{
    applyProcessingFunctional(
        new SetCustomBoundaryVelocityWithForceFunctional3D<T, Descriptor, VelocityFunction>(
            force, f),
        domain, lattice);
}

template <typename T, template <class U> class Descriptor>
void setBoundaryDensity(MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, T rho)
{
    applyProcessingFunctional(
        new SetConstBoundaryDensityFunctional3D<T, Descriptor>(rho), domain, lattice);
}

template <typename T, template <class U> class Descriptor>
void initializeAtEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, T rho, Array<T, 3> velocity,
    T temperature)
{
    applyProcessingFunctional(
        new IniConstEquilibriumFunctional3D<T, Descriptor>(rho, velocity, temperature), domain,
        lattice);
}

template <typename T, template <class U> class Descriptor, class DomainFunctional>
void maskedInitializeAtEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D boundingBox, DomainFunctional const &domain,
    T density, Array<T, 3> velocity, T temperature)
{
    applyProcessingFunctional(
        new IniConstEquilibriumOnDomainFunctional3D<T, Descriptor, DomainFunctional>(
            density, velocity, temperature, domain),
        boundingBox, lattice);
}

template <typename T, template <class U> class Descriptor>
void maskedInitializeAtEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<int> &mask, Box3D boundingBox,
    T density, Array<T, 3> velocity, int whichFlag, T temperature)
{
    applyProcessingFunctional(
        new MaskedIniConstEquilibriumFunctional3D<T, Descriptor>(
            density, velocity, temperature, whichFlag),
        boundingBox, lattice, mask);
}

template <typename T, template <class U> class Descriptor>
void initializeAtEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    Box3D domain, T density, Array<T, Descriptor<T>::d> velocity, T temperature)
{
    applyProcessingFunctional(
        new IniConstTensorForceEquilibriumFunctional3D<T, Descriptor>(
            density, velocity, temperature),
        domain, lattice, force);
}

template <typename T, template <class U> class Descriptor, class RhoUFunction>
void initializeAtEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    Box3D domain, RhoUFunction f, T temperature)
{
    applyProcessingFunctional(
        new IniCustomTensorForceEquilibriumFunctional3D<T, Descriptor, RhoUFunction>(
            f, temperature),
        domain, lattice, force);
}

template <typename T, template <class U> class Descriptor, class RhoUFunction>
void initializeAtEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    std::vector<MultiScalarField3D<T> *> randomFields, Box3D domain, RhoUFunction f, T temperature)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&lattice);
    args.push_back(&force);
    for (size_t i = 0; i < randomFields.size(); i++) {
        args.push_back(randomFields[i]);
    }

    applyProcessingFunctional(
        new IniCustomTensorForceRandomEquilibriumFunctional3D<T, Descriptor, RhoUFunction>(
            f, temperature),
        domain, args);
}

template <typename T, template <class U> class Descriptor>
void initializeAtEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force, Box3D domain,
    T density, Array<T, Descriptor<T>::d> velocity, T temperature)
{
    applyProcessingFunctional(
        new IniConstForceEquilibriumFunctional3D<T, Descriptor>(
            force, density, velocity, temperature),
        domain, lattice);
}

template <typename T, template <class U> class Descriptor, class RhoUFunction>
void initializeAtEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force, Box3D domain,
    RhoUFunction f, T temperature)
{
    applyProcessingFunctional(
        new IniCustomForceEquilibriumFunctional3D<T, Descriptor, RhoUFunction>(
            force, f, temperature),
        domain, lattice);
}

template <typename T, template <class U> class Descriptor, class RhoUFunction>
void initializeAtEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force,
    std::vector<MultiScalarField3D<T> *> randomFields, Box3D domain, RhoUFunction f, T temperature)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&lattice);
    for (size_t i = 0; i < randomFields.size(); i++) {
        args.push_back(randomFields[i]);
    }

    applyProcessingFunctional(
        new IniCustomForceRandomEquilibriumFunctional3D<T, Descriptor, RhoUFunction>(
            force, f, temperature),
        domain, args);
}

template <typename T, template <class U> class Descriptor>
void stripeOffDensityOffset(MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, T deltaRho)
{
    applyProcessingFunctional(
        new StripeOffDensityOffsetFunctional3D<T, Descriptor>(deltaRho), domain, lattice);
}

template <typename T, template <class U> class Descriptor>
void setCompositeDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain,
    CompositeDynamics<T, Descriptor> *compositeDynamics)
{
    applyProcessingFunctional(
        new InstantiateCompositeDynamicsFunctional3D<T, Descriptor>(compositeDynamics), domain,
        lattice);
}

template <typename T, template <class U> class Descriptor>
void setExternalScalar(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, int whichScalar, T externalScalar)
{
    applyProcessingFunctional(
        new SetExternalScalarFunctional3D<T, Descriptor>(whichScalar, externalScalar), domain,
        lattice);
}

template <typename T, template <class U> class Descriptor>
void setExternalScalar(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, int whichScalar,
    MultiScalarField3D<T> &scalar)
{
    applyProcessingFunctional(
        new SetExternalScalarFromScalarFieldFunctional3D<T, Descriptor>(whichScalar), domain,
        lattice, scalar);
}

template <typename T, template <class U> class Descriptor>
void setExternalScalar(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<int> &mask, int flag,
    Box3D domain, int whichScalar, T externalScalar)
{
    applyProcessingFunctional(
        new MaskedSetExternalScalarFunctional3D<T, Descriptor>(flag, whichScalar, externalScalar),
        domain, lattice, mask);
}

template <typename T, template <class U> class Descriptor, class Functional>
void setGenericExternalScalar(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, int whichScalar,
    Functional const &functional)
{
    applyProcessingFunctional(
        new SetGenericExternalScalarFunctional3D<T, Descriptor, Functional>(
            whichScalar, functional),
        domain, lattice);
}

template <typename T, template <class U> class Descriptor, class Functional>
void setGenericExternalScalar(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<int> &mask, int flag,
    Box3D domain, int whichScalar, Functional const &functional)
{
    applyProcessingFunctional(
        new MaskedSetGenericExternalScalarFunctional3D<T, Descriptor, Functional>(
            flag, whichScalar, functional),
        domain, lattice, mask);
}

template <typename T, template <class U> class Descriptor>
void setExternalVector(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, int vectorStartsAt,
    Array<T, Descriptor<T>::d> externalVector)
{
    applyProcessingFunctional(
        new SetExternalVectorFunctional3D<T, Descriptor>(vectorStartsAt, externalVector), domain,
        lattice);
}

template <typename T, template <class U> class Descriptor>
void setExternalVector(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<int> &mask, int flag,
    Box3D domain, int vectorStartsAt, Array<T, Descriptor<T>::d> externalVector)
{
    applyProcessingFunctional(
        new MaskedSetExternalVectorFunctional3D<T, Descriptor>(
            flag, vectorStartsAt, externalVector),
        domain, lattice, mask);
}

template <typename T, template <class U> class Descriptor, int nDim>
void setExternalVector(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, int vectorStartsAt,
    MultiTensorField3D<T, nDim> &tensor)
{
    applyProcessingFunctional(
        new SetExternalVectorFromTensorFieldFunctional3D<T, Descriptor, nDim>(vectorStartsAt),
        domain, lattice, tensor);
}

template <typename T, template <class U> class Descriptor, class Functional>
void setExternalVector(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, int vectorStartsAt,
    Functional const &functional)
{
    applyProcessingFunctional(
        new SetGenericExternalVectorFunctional3D<T, Descriptor, Functional>(
            vectorStartsAt, functional),
        domain, lattice);
}

template <typename T, template <class U> class Descriptor, class Functional>
void setExternalVector(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<int> &mask, int flag,
    Box3D domain, int vectorStartsAt, Functional const &functional)
{
    applyProcessingFunctional(
        new MaskedSetGenericExternalVectorFunctional3D<T, Descriptor, Functional>(
            flag, vectorStartsAt, functional),
        domain, lattice, mask);
}

template <typename T, template <class U> class Descriptor>
void interpolatePopulations(
    MultiBlockLattice3D<T, Descriptor> &modifiedLattice,
    MultiBlockLattice3D<T, Descriptor> &constLattice, plint minIter, plint maxIter,
    Box3D const &domain)
{
    applyProcessingFunctional(
        new InterpolatePopulationsFunctional3D<T, Descriptor>(minIter, maxIter), domain,
        modifiedLattice, constLattice);
}

template <typename T, template <class U> class Descriptor>
void interpolatePopulations(
    MultiBlockLattice3D<T, Descriptor> &modifiedLattice,
    MultiBlockLattice3D<T, Descriptor> &constLattice, plint minIter, plint maxIter)
{
    interpolatePopulations(
        modifiedLattice, constLattice, minIter, maxIter, modifiedLattice.getBoundingBox());
}

/* *************** PART III ******************************************* */
/* *************** Initialization of scalar- and tensor-fields: ****** */
/* *************** Atomic-Block  ************************************* */
/* ******************************************************************* */

template <typename T>
void setToConstant(ScalarField3D<T> &field, Box3D domain, T value)
{
    applyProcessingFunctional(new IniConstScalarFunctional3D<T>(value), domain, field);
}

template <typename T>
void setToConstant(
    ScalarField3D<T> &field, ScalarField3D<int> &mask, int flag, Box3D domain, T value)
{
    applyProcessingFunctional(
        new MaskedIniConstScalarFunctional3D<T>(flag, value), domain, field, mask);
}

template <typename T, int nDim>
void setToConstant(TensorField3D<T, nDim> &field, Box3D domain, Array<T, nDim> const &value)
{
    applyProcessingFunctional(new IniConstTensorFunctional3D<T, nDim>(value), domain, field);
}

template <typename T, int nDim>
void setToConstant(
    TensorField3D<T, nDim> &field, ScalarField3D<int> &mask, int flag, Box3D domain,
    Array<T, nDim> const &value)
{
    applyProcessingFunctional(
        new MaskedIniConstTensorFunctional3D<T, nDim>(flag, value), domain, mask, field);
}

template <typename T>
void setToCoordinate(ScalarField3D<T> &field, Box3D domain, plint index)
{
    applyProcessingFunctional(new SetToCoordinateFunctional3D<T>(index), domain, field);
}

template <typename T>
void setToCoordinates(TensorField3D<T, 3> &field, Box3D domain)
{
    applyProcessingFunctional(new SetToCoordinatesFunctional3D<T>, domain, field);
}

template <typename T, int nDim>
void assignComponent(
    TensorField3D<T, nDim> &tensorField, int whichComponent, ScalarField3D<T> &scalarField,
    Box3D domain)
{
    applyProcessingFunctional(
        new SetTensorComponentFunctional3D<T, nDim>(whichComponent), domain, scalarField,
        tensorField);
}

/* *************** PART IV ******************************************* */
/* *************** Initialization of scalar- and tensor-fields: ****** */
/* *************** Multi-Block  ************************************** */
/* ******************************************************************* */

template <typename T>
void setToConstant(MultiScalarField3D<T> &field, Box3D domain, T value)
{
    applyProcessingFunctional(new IniConstScalarFunctional3D<T>(value), domain, field);
}

template <typename T>
void setToConstant(
    MultiScalarField3D<T> &field, MultiScalarField3D<int> &mask, int flag, Box3D domain, T value)
{
    applyProcessingFunctional(
        new MaskedIniConstScalarFunctional3D<T>(flag, value), domain, field, mask);
}

template <typename T>
void setToConstant(
    MultiScalarField3D<T> &field, MultiNTensorField3D<int> &mask, int flag, Box3D domain, T value)
{
    applyProcessingFunctional(
        new MaskedIniConstScalarFunctional3D_N<T>(flag, value), domain, field, mask);
}

template <typename T, int nDim>
void setToConstant(MultiTensorField3D<T, nDim> &field, Box3D domain, Array<T, nDim> const &value)
{
    applyProcessingFunctional(new IniConstTensorFunctional3D<T, nDim>(value), domain, field);
}

template <typename T, int nDim>
void setToConstant(
    MultiTensorField3D<T, nDim> &field, MultiScalarField3D<int> &mask, int flag, Box3D domain,
    Array<T, nDim> const &value)
{
    applyProcessingFunctional(
        new MaskedIniConstTensorFunctional3D<T, nDim>(flag, value), domain, mask, field);
}

template <typename T, int nDim>
void setToConstant(
    MultiTensorField3D<T, nDim> &field, MultiNTensorField3D<int> &mask, int flag, Box3D domain,
    Array<T, nDim> const &value)
{
    std::vector<MultiBlock3D *> blocks;
    blocks.push_back(&mask);
    blocks.push_back(&field);
    applyProcessingFunctional(
        new MaskedIniConstTensorFunctional3D_N<T, nDim>(flag, value), domain, blocks);
}

template <typename T>
void setToCoordinate(MultiScalarField3D<T> &field, Box3D domain, plint index)
{
    applyProcessingFunctional(new SetToCoordinateFunctional3D<T>(index), domain, field);
}

template <typename T>
void setToCoordinates(MultiTensorField3D<T, 3> &field, Box3D domain)
{
    applyProcessingFunctional(new SetToCoordinatesFunctional3D<T>, domain, field);
}

template <typename T>
void setToRandom(MultiScalarField3D<T> &field, Box3D domain, uint32_t seed)
{
    applyProcessingFunctional(
        new SetToRandomFunctional3D<T>(field.getBoundingBox(), &seed), domain, field);
}

template <typename T, int nDim>
void assignComponent(
    MultiTensorField3D<T, nDim> &tensorField, int whichComponent,
    MultiScalarField3D<T> &scalarField, Box3D domain)
{
    applyProcessingFunctional(
        new SetTensorComponentFunctional3D<T, nDim>(whichComponent), domain, scalarField,
        tensorField);
}

template <typename T>
void propagateInZdirection(MultiScalarField3D<T> &field)
{
    field.resetFlags();
    plint maxIterations = field.getNz();
    plint i = 0;
    while (!allFlagsTrue(&field) && i < maxIterations) {
        applyProcessingFunctional(new PropagateInZdirection3D<T>(), field.getBoundingBox(), field);
        ++i;
    }
    if (i == maxIterations) {
        pcout << "Propagation in z-direction failed." << std::endl;
    }
}

template <typename T>
void growDomain(MultiScalarField3D<T> &field, T flag, int nCells, Box3D domain)
{
    for (int i = 0; i < nCells; ++i) {
        applyProcessingFunctional(new GrowDomainFunctional3D<T>(flag), domain, field);
    }
}

template <typename T>
void growDomain(MultiScalarField3D<T> &field, T flag, int nCells)
{
    growDomain(field, flag, nCells, field.getBoundingBox());
}

}  // namespace plb

#endif  // DATA_INITIALIZER_WRAPPER_3D_HH
