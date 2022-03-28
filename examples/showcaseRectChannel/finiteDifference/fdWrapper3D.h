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
 * Helper functions for finite differences -- header file.
 */

#ifndef FINITE_DIFFERENCE_WRAPPER_3D_H
#define FINITE_DIFFERENCE_WRAPPER_3D_H

#include <memory>

#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "core/globalDefs.h"
#include "dataProcessors/dataAnalysisFunctional3D.h"
#include "finiteDifference/fdWeights.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"

namespace plb {

template <typename T>
void computeXderivative(
    MultiScalarField3D<T> &value, MultiScalarField3D<T> &derivative, Box3D const &domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeXderivative(
    MultiScalarField3D<T> &value, Box3D const &domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeXderivative(MultiScalarField3D<T> &value);

template <typename T>
void computeYderivative(
    MultiScalarField3D<T> &value, MultiScalarField3D<T> &derivative, Box3D const &domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeYderivative(
    MultiScalarField3D<T> &value, Box3D const &domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeYderivative(MultiScalarField3D<T> &value);

template <typename T>
void computeZderivative(
    MultiScalarField3D<T> &value, MultiScalarField3D<T> &derivative, Box3D const &domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeZderivative(
    MultiScalarField3D<T> &value, Box3D const &domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeZderivative(MultiScalarField3D<T> &value);

template <typename T>
void computeGradientNorm(
    MultiScalarField3D<T> &value, MultiScalarField3D<T> &derivative, Box3D const &domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeGradientNorm(
    MultiScalarField3D<T> &value, Box3D const &domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeGradientNorm(MultiScalarField3D<T> &value);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computePoissonRHS(
    MultiTensorField3D<T, 3> &velocity, Box3D const &domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computePoissonRHS(MultiTensorField3D<T, 3> &velocity);

template <typename T>
void poissonIterate(
    MultiScalarField3D<T> &oldPressure, MultiScalarField3D<T> &newPressure,
    MultiScalarField3D<T> &rhs, T beta, Box3D const &domain, plint boundaryWidth = 1);

template <typename T>
T computePoissonResidue(
    MultiScalarField3D<T> &pressure, MultiScalarField3D<T> &rhs, Box3D const &domain);

// ========================================================================= //
// PERIODIC VERSIONS OF THE DERIVATIVES AND POISSON SCHEMES //
// ========================================================================= //

template <typename T>
void computeXperiodicDerivative(
    MultiScalarField3D<T> &value, MultiScalarField3D<T> &derivative, Box3D const &domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeXperiodicDerivative(
    MultiScalarField3D<T> &value, Box3D const &domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeXperiodicDerivative(MultiScalarField3D<T> &value);

template <typename T>
void computeYperiodicDerivative(
    MultiScalarField3D<T> &value, MultiScalarField3D<T> &derivative, Box3D const &domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeYperiodicDerivative(
    MultiScalarField3D<T> &value, Box3D const &domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeYperiodicDerivative(MultiScalarField3D<T> &value);

template <typename T>
void computeZperiodicDerivative(
    MultiScalarField3D<T> &value, MultiScalarField3D<T> &derivative, Box3D const &domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeZperiodicDerivative(
    MultiScalarField3D<T> &value, Box3D const &domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeZperiodicDerivative(MultiScalarField3D<T> &value);

template <typename T>
void computePeriodicGradientNorm(
    MultiScalarField3D<T> &value, MultiScalarField3D<T> &derivative, Box3D const &domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computePeriodicGradientNorm(
    MultiScalarField3D<T> &value, Box3D const &domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computePeriodicGradientNorm(MultiScalarField3D<T> &value);

template <typename T>
void computePeriodicGradient(
    MultiScalarField3D<T> &value, MultiTensorField3D<T, 3> &derivative, Box3D const &domain);

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 2> > computePeriodicGradient(
    MultiScalarField3D<T> &value, Box3D const &domain);

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 2> > computePeriodicGradient(MultiScalarField3D<T> &value);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computePeriodicPoissonRHS(
    MultiTensorField3D<T, 3> &velocity, Box3D const &domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computePeriodicPoissonRHS(
    MultiTensorField3D<T, 3> &velocity);

template <typename T>
void periodicPoissonIterate(
    MultiScalarField3D<T> &oldPressure, MultiScalarField3D<T> &newPressure,
    MultiScalarField3D<T> &rhs, T beta, Box3D const &domain);

// General Stencils.

// Compute the X/Y/Z derivative of a data field. The order of the derivative to be approximated
// is provided by the template parameter "order". The specific stencil is defined by "width"
// and "position". The value of the argument "width" cannot be greater than that of the template
// parameter "maxWidth". For best efficiency, one should use these functions with a fixed
// "maxWidth". This means that if one wants to use many widths, he should fix "maxWidth", and vary
// "width". The returned derivative value is computed on the local grid point (iX, iY, iZ) of the
// atomic-block.

template <typename T, int order, int maxWidth>
T computeScalarXderivative(
    ScalarField3D<T> const &scalar, int width, int position, plint iX, plint iY, plint iZ);

template <typename T, int order, int maxWidth>
T computeScalarYderivative(
    ScalarField3D<T> const &scalar, int width, int position, plint iX, plint iY, plint iZ);

template <typename T, int order, int maxWidth>
T computeScalarZderivative(
    ScalarField3D<T> const &scalar, int width, int position, plint iX, plint iY, plint iZ);

template <typename T, int nDim, int order, int maxWidth>
Array<T, nDim> computeTensorXderivative(
    TensorField3D<T, nDim> const &tensor, int width, int position, plint iX, plint iY, plint iZ);

template <typename T, int nDim, int order, int maxWidth>
Array<T, nDim> computeTensorYderivative(
    TensorField3D<T, nDim> const &tensor, int width, int position, plint iX, plint iY, plint iZ);

template <typename T, int nDim, int order, int maxWidth>
Array<T, nDim> computeTensorZderivative(
    TensorField3D<T, nDim> const &tensor, int width, int position, plint iX, plint iY, plint iZ);

}  // namespace plb

#endif  // FINITE_DIFFERENCE_WRAPPER_3D_H
