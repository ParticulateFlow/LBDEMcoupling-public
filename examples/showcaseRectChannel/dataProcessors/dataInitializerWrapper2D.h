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
#ifndef DATA_INITIALIZER_WRAPPER_2D_H
#define DATA_INITIALIZER_WRAPPER_2D_H

#include "atomicBlock/dataProcessorWrapper2D.h"
#include "core/dynamics.h"
#include "core/globalDefs.h"
#include "dataProcessors/dataInitializerFunctional2D.h"
#include "multiBlock/multiDataProcessorWrapper2D.h"

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Initialization of the block-lattice: atomic-blocks  */
/* ******************************************************************* */

template <typename T, template <class U> class Descriptor>
void apply(
    BlockLattice2D<T, Descriptor> &lattice, Box2D domain, OneCellFunctional2D<T, Descriptor> *f);

template <typename T, template <class U> class Descriptor>
void applyIndexed(
    BlockLattice2D<T, Descriptor> &lattice, Box2D domain,
    OneCellIndexedFunctional2D<T, Descriptor> *f);

template <typename T, template <class U> class Descriptor>
void defineDynamics(
    BlockLattice2D<T, Descriptor> &lattice, Box2D domain, Dynamics<T, Descriptor> *dynamics);

template <typename T, template <class U> class Descriptor>
void defineDynamics(
    BlockLattice2D<T, Descriptor> &lattice, Box2D boundingBox, DomainFunctional2D *domain,
    Dynamics<T, Descriptor> *dynamics);

template <typename T, template <class U> class Descriptor>
void defineDynamics(
    BlockLattice2D<T, Descriptor> &lattice, plint iX, plint iY, Dynamics<T, Descriptor> *dynamics);

template <typename T, template <class U> class Descriptor>
void defineDynamics(
    BlockLattice2D<T, Descriptor> &lattice, DotList2D const &dotList,
    Dynamics<T, Descriptor> *dynamics);

template <typename T, template <typename U> class Descriptor>
void defineDynamics(
    BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<bool> &boolMask, Box2D domain,
    Dynamics<T, Descriptor> *dynamics, bool whichFlag);

template <typename T, template <typename U> class Descriptor>
void defineDynamics(
    BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<bool> &boolMask,
    Dynamics<T, Descriptor> *dynamics, bool whichFlag);

template <typename T, template <typename U> class Descriptor>
void defineDynamics(
    BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<int> &intMask, Box2D domain,
    Dynamics<T, Descriptor> *dynamics, int whichFlag);

template <typename T, template <typename U> class Descriptor>
void defineDynamics(
    BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<int> &intMask,
    Dynamics<T, Descriptor> *dynamics, int whichFlag);

template <typename T, template <typename U> class Descriptor>
void recomposeFromFlowVariables(
    BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &density,
    TensorField2D<T, 2> &velocity, TensorField2D<T, 3> &strainRate, Box2D domain);

template <typename T, template <typename U> class Descriptor>
void recomposeFromFlowVariables(
    BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &density,
    TensorField2D<T, 2> &velocity, TensorField2D<T, 3> &strainRate);

template <typename T, template <class U> class Descriptor>
void setOmega(BlockLattice2D<T, Descriptor> &lattice, Box2D domain, T omega);

template <typename T, template <class U> class Descriptor>
void setBoundaryVelocity(
    BlockLattice2D<T, Descriptor> &lattice, Box2D domain, Array<T, 2> velocity);

// This function is implemented in-place, because it cannot be precompiled due to its generic
// nature.
template <typename T, template <class U> class Descriptor, class VelocityFunction>
void setBoundaryVelocity(BlockLattice2D<T, Descriptor> &lattice, Box2D domain, VelocityFunction f)
{
    applyIndexed(
        lattice, domain,
        new SetCustomBoundaryVelocityFunctional2D<T, Descriptor, VelocityFunction>(f));
}

template <typename T, template <class U> class Descriptor>
void setBoundaryDensity(BlockLattice2D<T, Descriptor> &lattice, Box2D domain, T rho);

// This function is implemented in-place, because it cannot be precompiled due to its generic
// nature.
template <typename T, template <class U> class Descriptor, class DensityFunction>
void setBoundaryDensity(BlockLattice2D<T, Descriptor> &lattice, Box2D domain, DensityFunction f)
{
    applyIndexed(
        lattice, domain,
        new SetCustomBoundaryDensityFunctional2D<T, Descriptor, DensityFunction>(f));
}

template <typename T, template <class U> class Descriptor>
void setBoundaryTemperature(BlockLattice2D<T, Descriptor> &lattice, Box2D domain, T temperature);

// This function is implemented in-place, because it cannot be precompiled due to its generic
// nature.
template <typename T, template <class U> class Descriptor, class TemperatureFunction>
void setBoundaryTemperature(
    BlockLattice2D<T, Descriptor> &lattice, Box2D domain, TemperatureFunction f)
{
    applyIndexed(
        lattice, domain,
        new SetCustomBoundaryTemperatureFunctional2D<T, Descriptor, TemperatureFunction>(f));
}

template <typename T, template <class U> class Descriptor>
void initializeAtEquilibrium(
    BlockLattice2D<T, Descriptor> &lattice, Box2D domain, T density, Array<T, 2> velocity,
    T temperature = (T)1);

// This function is implemented in-place, because it cannot be precompiled due to its generic
// nature.
template <typename T, template <class U> class Descriptor, class RhoUFunction>
void initializeAtEquilibrium(BlockLattice2D<T, Descriptor> &lattice, Box2D domain, RhoUFunction f)
{
    applyIndexed(
        lattice, domain, new IniCustomEquilibriumFunctional2D<T, Descriptor, RhoUFunction>(f));
}

// This function is implemented in-place, because it cannot be precompiled due to its generic
// nature.
template <typename T, template <class U> class Descriptor, class RhoVelTempFunction>
void initializeAtThermalEquilibrium(
    BlockLattice2D<T, Descriptor> &lattice, Box2D domain, RhoVelTempFunction f)
{
    applyIndexed(
        lattice, domain,
        new IniCustomThermalEquilibriumFunctional2D<T, Descriptor, RhoVelTempFunction>(f));
}

template <typename T, template <class U> class Descriptor>
void stripeOffDensityOffset(BlockLattice2D<T, Descriptor> &lattice, Box2D domain, T deltaRho);

template <typename T, template <class U> class Descriptor>
void setCompositeDynamics(
    BlockLattice2D<T, Descriptor> &lattice, Box2D domain,
    CompositeDynamics<T, Descriptor> *compositeDynamics);

template <typename T, template <class U> class Descriptor>
void setExternalScalar(
    BlockLattice2D<T, Descriptor> &lattice, Box2D domain, int whichScalar, T externalScalar);

template <typename T, template <class U> class Descriptor>
void setExternalVector(
    BlockLattice2D<T, Descriptor> &lattice, Box2D domain, int vectorStartsAt,
    Array<T, Descriptor<T>::d> externalVector);

/* *************** PART II ******************************************* */
/* *************** Initialization of the block-lattice: multi-block ** */
/* ******************************************************************* */

template <typename T, template <class U> class Descriptor>
void apply(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain,
    OneCellFunctional2D<T, Descriptor> *f);

template <typename T, template <class U> class Descriptor>
void applyIndexed(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain,
    OneCellIndexedFunctional2D<T, Descriptor> *f);

template <typename T, template <class U> class Descriptor>
void applyIndexed(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain,
    OneCellIndexedWithRandFunctional2D<T, Descriptor> *f, uint32_t seed = 0);

template <typename T, template <class U> class Descriptor>
void defineDynamics(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, Dynamics<T, Descriptor> *dynamics);

template <typename T, template <class U> class Descriptor>
void defineDynamics(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D boundingBox, DomainFunctional2D *domain,
    Dynamics<T, Descriptor> *dynamics);

template <typename T, template <class U> class Descriptor>
void defineDynamics(
    MultiBlockLattice2D<T, Descriptor> &lattice, plint iX, plint iY,
    Dynamics<T, Descriptor> *dynamics);

template <typename T, template <class U> class Descriptor>
void defineDynamics(
    MultiBlockLattice2D<T, Descriptor> &lattice, DotList2D const &dotList,
    Dynamics<T, Descriptor> *dynamics);

template <typename T, template <typename U> class Descriptor>
void defineDynamics(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<bool> &boolMask, Box2D domain,
    Dynamics<T, Descriptor> *dynamics, bool whichFlag);

template <typename T, template <typename U> class Descriptor>
void defineDynamics(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<bool> &boolMask,
    Dynamics<T, Descriptor> *dynamics, bool whichFlag);

template <typename T, template <typename U> class Descriptor>
void defineDynamics(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<int> &intMask, Box2D domain,
    Dynamics<T, Descriptor> *dynamics, int whichFlag);

template <typename T, template <typename U> class Descriptor>
void defineDynamics(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<int> &intMask,
    Dynamics<T, Descriptor> *dynamics, int whichFlag);

template <typename T, template <typename U> class Descriptor>
void recomposeFromFlowVariables(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &density,
    MultiTensorField2D<T, 2> &velocity, MultiTensorField2D<T, 3> &strainRate, Box2D domain);

template <typename T, template <typename U> class Descriptor>
void recomposeFromFlowVariables(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &density,
    MultiTensorField2D<T, 2> &velocity, MultiTensorField2D<T, 3> &strainRate);

template <typename T, template <class U> class Descriptor>
void setOmega(MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, T omega);

template <typename T, template <class U> class Descriptor, class Function>
void setOmega(MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, Function f);

template <typename T, template <class U> class Descriptor>
void setBoundaryVelocity(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, Array<T, 2> velocity);

// This function is implemented in-place, because it cannot be precompiled due to its generic
// nature.
template <typename T, template <class U> class Descriptor, class VelocityFunction>
void setBoundaryVelocity(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, VelocityFunction f)
{
    applyIndexed(
        lattice, domain,
        new SetCustomBoundaryVelocityFunctional2D<T, Descriptor, VelocityFunction>(f));
}

template <typename T, template <class U> class Descriptor>
void setBoundaryDensity(MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, T rho);

// This function is implemented in-place, because it cannot be precompiled due to its generic
// nature.
template <typename T, template <class U> class Descriptor, class DensityFunction>
void setBoundaryDensity(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, DensityFunction f)
{
    applyIndexed(
        lattice, domain,
        new SetCustomBoundaryDensityFunctional2D<T, Descriptor, DensityFunction>(f));
}

template <typename T, template <class U> class Descriptor>
void setBoundaryTemperature(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, T temperature);

// This function is implemented in-place, because it cannot be precompiled due to its generic
// nature.
template <typename T, template <class U> class Descriptor, class TemperatureFunction>
void setBoundaryTemperature(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, TemperatureFunction f)
{
    applyIndexed(
        lattice, domain,
        new SetCustomBoundaryTemperatureFunctional2D<T, Descriptor, TemperatureFunction>(f));
}

template <typename T, template <class U> class Descriptor>
void initializeAtEquilibrium(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, T density, Array<T, 2> velocity,
    T temperature = (T)1);

// This function is implemented in-place, because it cannot be precompiled due to its generic
// nature.
template <typename T, template <class U> class Descriptor, class RhoUFunction>
void initializeAtEquilibrium(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, RhoUFunction f)
{
    applyIndexed(
        lattice, domain, new IniCustomEquilibriumFunctional2D<T, Descriptor, RhoUFunction>(f));
}

// This function is implemented in-place, because it cannot be precompiled due to its generic
// nature.
template <typename T, template <class U> class Descriptor, class RhoVelTempFunction>
void initializeAtThermalEquilibrium(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, RhoVelTempFunction f)
{
    applyIndexed(
        lattice, domain,
        new IniCustomThermalEquilibriumFunctional2D<T, Descriptor, RhoVelTempFunction>(f));
}

template <typename T, template <class U> class Descriptor>
void stripeOffDensityOffset(MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, T deltaRho);

template <typename T, template <class U> class Descriptor>
void setCompositeDynamics(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain,
    CompositeDynamics<T, Descriptor> *compositeDynamics);

template <typename T, template <class U> class Descriptor>
void setExternalScalar(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, int whichScalar, T externalScalar);

template <typename T, template <class U> class Descriptor>
void setExternalScalar(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, int whichScalar,
    MultiScalarField2D<T> &field);

template <typename T, template <class U> class Descriptor, class Functional>
void setGenericExternalScalar(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, int whichScalar,
    Functional const &functional);

template <typename T, template <class U> class Descriptor>
void setExternalVector(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, int vectorStartsAt,
    Array<T, Descriptor<T>::d> externalVector);

template <typename T, template <class U> class Descriptor, int nDim>
void setExternalVector(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, int vectorStartsAt,
    MultiTensorField2D<T, nDim> &tensor);

template <typename T, template <class U> class Descriptor, class VectorFunction>
void setExternalVector(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, int vectorStartsAt,
    VectorFunction f);

/* *************** PART III ******************************************* */
/* *************** Initialization of scalar- and tensor-fields: ****** */
/* *************** Atomic-Block  ************************************* */
/* ******************************************************************* */

/// Initialize scalar-field with the same constant value on each cell.
template <typename T>
void setToConstant(ScalarField2D<T> &field, Box2D domain, T value);

/// Initialize scalar-field with the same constant value on each cell on which
///   the mask is equal to the flag.
template <typename T>
void setToConstant(
    ScalarField2D<T> &field, ScalarField2D<int> &mask, int flag, Box2D domain, T value);

/// Initialize tensor-field with the same constant tensor/vector on each cell.
template <typename T, int nDim>
void setToConstant(TensorField2D<T, nDim> &field, Box2D domain, Array<T, nDim> const &value);

/// Initialize tensor-field with the same constant tensor/vector on each cell on which
///   the mask is equal to the flag.
template <typename T, int nDim>
void setToConstant(
    TensorField2D<T, nDim> &field, ScalarField2D<int> &mask, int flag, Box2D domain,
    Array<T, nDim> const &value);

/// Initialize scalar-field with the a value from a function.
// This function is implemented in-place, because it cannot be precompiled due to its generic
// nature.
template <typename T, class Function>
void setToFunction(ScalarField2D<T> &field, Box2D domain, Function f)
{
    applyProcessingFunctional(new SetToScalarFunctionFunctional2D<T, Function>(f), domain, field);
}

/// Initialize tensor-field with a vector/tensor value from a function.
// This function is implemented in-place, because it cannot be precompiled due to its generic
// nature.
template <typename T, int nDim, class Function>
void setToFunction(TensorField2D<T, nDim> &field, Box2D domain, Function f)
{
    applyProcessingFunctional(
        new SetToTensorFunctionFunctional2D<T, nDim, Function>(f), domain, field);
}

/// Assign the component "index" of its space coordinate to each cell.
template <typename T>
void setToCoordinate(ScalarField2D<T> &field, Box2D domain, plint index);

/// Assign its space coordinate to each cell.
template <typename T>
void setToCoordinates(TensorField2D<T, 2> &field, Box2D domain);

/// Assign scalar-field to one component of a tensor-field.
template <typename T, int nDim>
void assignComponent(
    TensorField2D<T, nDim> &tensorField, int whichComponent, ScalarField2D<T> &scalarField,
    Box2D domain);

/* *************** PART IV ******************************************* */
/* *************** Initialization of scalar- and tensor-fields: ****** */
/* *************** Multi-Block *************************************** */
/* ******************************************************************* */

/// Initialize scalar-field with the same constant value on each cell.
template <typename T>
void setToConstant(MultiScalarField2D<T> &field, Box2D domain, T value);

/// Initialize scalar-field with the same constant value on each cell on which
///   the mask is equal to the flag.
template <typename T>
void setToConstant(
    MultiScalarField2D<T> &field, MultiScalarField2D<int> &mask, int flag, Box2D domain, T value);

/// Initialize tensor-field with the same constant tensor/vector on each cell.
template <typename T, int nDim>
void setToConstant(MultiTensorField2D<T, nDim> &field, Box2D domain, Array<T, nDim> const &value);

/// Initialize tensor-field with the same constant tensor/vector on each cell on which
///   the mask is equal to the flag.
template <typename T, int nDim>
void setToConstant(
    MultiTensorField2D<T, nDim> &field, MultiScalarField2D<int> &mask, int flag, Box2D domain,
    Array<T, nDim> const &value);

/// Initialize scalar-field with the a value from a function.
// This function is implemented in-place, because it cannot be precompiled due to its generic
// nature.
template <typename T, class Function>
void setToFunction(MultiScalarField2D<T> &field, Box2D domain, Function f)
{
    applyProcessingFunctional(new SetToScalarFunctionFunctional2D<T, Function>(f), domain, field);
}

/// Initialize tensor-field with a vector/tensor value from a function.
// This function is implemented in-place, because it cannot be precompiled due to its generic
// nature.
template <typename T, int nDim, class Function>
void setToFunction(MultiTensorField2D<T, nDim> &field, Box2D domain, Function f)
{
    applyProcessingFunctional(
        new SetToTensorFunctionFunctional2D<T, nDim, Function>(f), domain, field);
}

/// Assign the component "index" of its space coordinate to each cell.
template <typename T>
void setToCoordinate(MultiScalarField2D<T> &field, Box2D domain, plint index);

/// Assign its space coordinate to each cell.
template <typename T>
void setToCoordinates(MultiTensorField2D<T, 2> &field, Box2D domain);

template <typename T>
void setToRandom(MultiScalarField2D<T> &field, Box2D domain, uint32_t seed = 0);

/// Assign scalar-field to one component of a tensor-field.
template <typename T, int nDim>
void assignComponent(
    MultiTensorField2D<T, nDim> &tensorField, int whichComponent,
    MultiScalarField2D<T> &scalarField, Box2D domain);

/// Grow a domain defined by the flag "flag" by n cells.
template <typename T>
void growDomain(MultiScalarField2D<T> &field, T flag, int nCells, Box2D domain);

template <typename T>
void growDomain(MultiScalarField2D<T> &field, T flag, int nCells);

}  // namespace plb

#endif  // DATA_INITIALIZER_WRAPPER_2D_H
