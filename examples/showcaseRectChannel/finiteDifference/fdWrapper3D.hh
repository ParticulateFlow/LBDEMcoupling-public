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

#ifndef FINITE_DIFFERENCE_WRAPPER_3D_HH
#define FINITE_DIFFERENCE_WRAPPER_3D_HH

#include "atomicBlock/dataProcessorWrapper3D.h"
#include "atomicBlock/reductiveDataProcessorWrapper3D.h"
#include "finiteDifference/fdFunctional3D.h"
#include "finiteDifference/fdWrapper3D.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper3D.h"

namespace plb {

template <typename T>
void computeXderivative(
    MultiScalarField3D<T> &value, MultiScalarField3D<T> &derivative, Box3D const &domain)
{
    plint boundaryWidth = 1;
    applyProcessingFunctional(
        new BoxXderivativeFunctional3D<T>, domain, value, derivative, boundaryWidth);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeXderivative(
    MultiScalarField3D<T> &value, Box3D const &domain)
{
    MultiScalarField3D<T> *derivative = new MultiScalarField3D<T>(value, domain);
    computeXderivative(value, *derivative, domain);
    return std::unique_ptr<MultiScalarField3D<T> >(derivative);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeXderivative(MultiScalarField3D<T> &value)
{
    return computeXderivative(value, value.getBoundingBox());
}

template <typename T>
void computeYderivative(
    MultiScalarField3D<T> &value, MultiScalarField3D<T> &derivative, Box3D const &domain)
{
    plint boundaryWidth = 1;
    applyProcessingFunctional(
        new BoxYderivativeFunctional3D<T>, domain, value, derivative, boundaryWidth);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeYderivative(
    MultiScalarField3D<T> &value, Box3D const &domain)
{
    MultiScalarField3D<T> *derivative = new MultiScalarField3D<T>(value, domain);
    computeYderivative(value, *derivative, domain);
    return std::unique_ptr<MultiScalarField3D<T> >(derivative);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeYderivative(MultiScalarField3D<T> &value)
{
    return computeYderivative(value, value.getBoundingBox());
}

template <typename T>
void computeZderivative(
    MultiScalarField3D<T> &value, MultiScalarField3D<T> &derivative, Box3D const &domain)
{
    plint boundaryWidth = 1;
    applyProcessingFunctional(
        new BoxZderivativeFunctional3D<T>, domain, value, derivative, boundaryWidth);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeZderivative(
    MultiScalarField3D<T> &value, Box3D const &domain)
{
    MultiScalarField3D<T> *derivative = new MultiScalarField3D<T>(value, domain);
    computeZderivative(value, *derivative, domain);
    return std::unique_ptr<MultiScalarField3D<T> >(derivative);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeZderivative(MultiScalarField3D<T> &value)
{
    return computeZderivative(value, value.getBoundingBox());
}

template <typename T>
void computeGradientNorm(
    MultiScalarField3D<T> &value, MultiScalarField3D<T> &derivative, Box3D const &domain)
{
    plint boundaryWidth = 1;
    applyProcessingFunctional(
        new BoxGradientNormFunctional3D<T>, domain, value, derivative, boundaryWidth);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeGradientNorm(
    MultiScalarField3D<T> &value, Box3D const &domain)
{
    MultiScalarField3D<T> *derivative = new MultiScalarField3D<T>(value, domain);
    computeGradientNorm(value, *derivative, domain);
    return std::unique_ptr<MultiScalarField3D<T> >(derivative);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeGradientNorm(MultiScalarField3D<T> &value)
{
    return computeGradientNorm(value, value.getBoundingBox());
}

template <typename T>
void computePeriodicGradient(
    MultiScalarField3D<T> &value, MultiTensorField3D<T, 3> &derivative, Box3D const &domain)
{
    applyProcessingFunctional(new BoxPeriodicGradientFunctional3D<T>, domain, value, derivative);
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computePeriodicGradient(
    MultiScalarField3D<T> &value, Box3D const &domain)
{
    MultiTensorField3D<T, 3> *derivative = new MultiTensorField3D<T, 2>(value, domain);
    computePeriodicGradient(value, *derivative, domain);
    return std::unique_ptr<MultiTensorField3D<T, 3> >(derivative);
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computePeriodicGradient(MultiScalarField3D<T> &value)
{
    return computePeriodicGradient(value, value.getBoundingBox());
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computePoissonRHS(
    MultiTensorField3D<T, 3> &velocity, Box3D const &domain)
{
    std::unique_ptr<MultiScalarField3D<T> > ux = extractComponent(velocity, domain, 0);
    std::unique_ptr<MultiScalarField3D<T> > uy = extractComponent(velocity, domain, 1);
    std::unique_ptr<MultiScalarField3D<T> > uz = extractComponent(velocity, domain, 2);

    std::unique_ptr<MultiScalarField3D<T> > dx_ux = computeXderivative(*ux, domain);
    std::unique_ptr<MultiScalarField3D<T> > dy_ux = computeYderivative(*ux, domain);
    std::unique_ptr<MultiScalarField3D<T> > dz_ux = computeZderivative(*ux, domain);
    std::unique_ptr<MultiScalarField3D<T> > dx_uy = computeXderivative(*uy, domain);
    std::unique_ptr<MultiScalarField3D<T> > dy_uy = computeYderivative(*uy, domain);
    std::unique_ptr<MultiScalarField3D<T> > dz_uy = computeZderivative(*uy, domain);
    std::unique_ptr<MultiScalarField3D<T> > dx_uz = computeXderivative(*uz, domain);
    std::unique_ptr<MultiScalarField3D<T> > dy_uz = computeYderivative(*uz, domain);
    std::unique_ptr<MultiScalarField3D<T> > dz_uz = computeZderivative(*uz, domain);

    std::unique_ptr<MultiScalarField3D<T> > term1 = multiply(*dx_ux, *dx_ux, domain);
    std::unique_ptr<MultiScalarField3D<T> > term2 = multiply(*dy_uy, *dy_uy, domain);
    std::unique_ptr<MultiScalarField3D<T> > term3 = multiply(*dz_uz, *dz_uz, domain);

    std::unique_ptr<MultiScalarField3D<T> > term4 =
        multiply((T)2, *multiply(*dx_uy, *dy_ux, domain), domain);
    std::unique_ptr<MultiScalarField3D<T> > term5 =
        multiply((T)2, *multiply(*dx_uz, *dz_ux, domain), domain);
    std::unique_ptr<MultiScalarField3D<T> > term6 =
        multiply((T)2, *multiply(*dy_uz, *dz_uy, domain), domain);

    std::unique_ptr<MultiScalarField3D<T> > rhs =
        add(*term6, *add(*term5, *add(*term4, *add(*term1, *add(*term2, *term3)))));
    return rhs;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computePoissonRHS(MultiTensorField3D<T, 3> &velocity)
{
    return computePoissonRHS(velocity, velocity.getBoundingBox());
}

template <typename T>
void poissonIterate(
    MultiScalarField3D<T> &oldPressure, MultiScalarField3D<T> &newPressure,
    MultiScalarField3D<T> &rhs, T beta, Box3D const &domain, plint boundaryWidth)
{
    std::vector<MultiScalarField3D<T> *> fields;
    fields.push_back(&oldPressure);
    fields.push_back(&newPressure);
    fields.push_back(&rhs);
    applyProcessingFunctional(new BoxPoissonIteration3D<T>(beta), domain, fields, boundaryWidth);
}

template <typename T>
T computePoissonResidue(
    MultiScalarField3D<T> &pressure, MultiScalarField3D<T> &rhs, Box3D const &domain)
{
    BoxPoissonResidueFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, pressure, rhs);
    return functional.getMaxResidue();
}

// ========================================================================= //
// PERIODIC VERSIONS OF THE DERIVATIVES AND POISSON SCHEMES //
// ========================================================================= //

template <typename T>
void computeXperiodicDerivative(
    MultiScalarField3D<T> &value, MultiScalarField3D<T> &derivative, Box3D const &domain)
{
    plint boundaryWidth = 1;
    applyProcessingFunctional(
        new BoxXperiodicDerivativeFunctional3D<T>, domain, value, derivative, boundaryWidth);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeXperiodicDerivative(
    MultiScalarField3D<T> &value, Box3D const &domain)
{
    MultiScalarField3D<T> *derivative = new MultiScalarField3D<T>(value, domain);
    computeXperiodicDerivative(value, *derivative, domain);
    return std::unique_ptr<MultiScalarField3D<T> >(derivative);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeXperiodicDerivative(MultiScalarField3D<T> &value)
{
    return computeXperiodicDerivative(value, value.getBoundingBox());
}

template <typename T>
void computeYperiodicDerivative(
    MultiScalarField3D<T> &value, MultiScalarField3D<T> &derivative, Box3D const &domain)
{
    plint boundaryWidth = 1;
    applyProcessingFunctional(
        new BoxYperiodicDerivativeFunctional3D<T>, domain, value, derivative, boundaryWidth);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeYperiodicDerivative(
    MultiScalarField3D<T> &value, Box3D const &domain)
{
    MultiScalarField3D<T> *derivative = new MultiScalarField3D<T>(value, domain);
    computeYperiodicDerivative(value, *derivative, domain);
    return std::unique_ptr<MultiScalarField3D<T> >(derivative);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeYperiodicDerivative(MultiScalarField3D<T> &value)
{
    return computeYperiodicDerivative(value, value.getBoundingBox());
}

template <typename T>
void computeZperiodicDerivative(
    MultiScalarField3D<T> &value, MultiScalarField3D<T> &derivative, Box3D const &domain)
{
    plint boundaryWidth = 1;
    applyProcessingFunctional(
        new BoxZperiodicDerivativeFunctional3D<T>, domain, value, derivative, boundaryWidth);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeZperiodicDerivative(
    MultiScalarField3D<T> &value, Box3D const &domain)
{
    MultiScalarField3D<T> *derivative = new MultiScalarField3D<T>(value, domain);
    computeZperiodicDerivative(value, *derivative, domain);
    return std::unique_ptr<MultiScalarField3D<T> >(derivative);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeZperiodicDerivative(MultiScalarField3D<T> &value)
{
    return computeZperiodicDerivative(value, value.getBoundingBox());
}

template <typename T>
void computePeriodicGradientNorm(
    MultiScalarField3D<T> &value, MultiScalarField3D<T> &derivative, Box3D const &domain)
{
    plint boundaryWidth = 1;
    applyProcessingFunctional(
        new BoxGradientNormFunctional3D<T>, domain, value, derivative, boundaryWidth);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computePeriodicGradientNorm(
    MultiScalarField3D<T> &value, Box3D const &domain)
{
    MultiScalarField3D<T> *derivative = new MultiScalarField3D<T>(value, domain);
    computePeriodicGradientNorm(value, *derivative, domain);
    return std::unique_ptr<MultiScalarField3D<T> >(derivative);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computePeriodicGradientNorm(MultiScalarField3D<T> &value)
{
    return computePeriodicGradientNorm(value, value.getBoundingBox());
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computePeriodicPoissonRHS(
    MultiTensorField3D<T, 3> &velocity, Box3D const &domain)
{
    std::unique_ptr<MultiScalarField3D<T> > ux = extractComponent(velocity, domain, 0);
    ux->periodicity().toggleAll(true);
    std::unique_ptr<MultiScalarField3D<T> > uy = extractComponent(velocity, domain, 1);
    uy->periodicity().toggleAll(true);
    std::unique_ptr<MultiScalarField3D<T> > uz = extractComponent(velocity, domain, 2);
    uz->periodicity().toggleAll(true);

    std::unique_ptr<MultiScalarField3D<T> > dx_ux = computeXperiodicDerivative(*ux, domain);
    std::unique_ptr<MultiScalarField3D<T> > dy_ux = computeYperiodicDerivative(*ux, domain);
    std::unique_ptr<MultiScalarField3D<T> > dz_ux = computeZperiodicDerivative(*ux, domain);
    std::unique_ptr<MultiScalarField3D<T> > dx_uy = computeXperiodicDerivative(*uy, domain);
    std::unique_ptr<MultiScalarField3D<T> > dy_uy = computeYperiodicDerivative(*uy, domain);
    std::unique_ptr<MultiScalarField3D<T> > dz_uy = computeZperiodicDerivative(*uy, domain);
    std::unique_ptr<MultiScalarField3D<T> > dx_uz = computeXperiodicDerivative(*uz, domain);
    std::unique_ptr<MultiScalarField3D<T> > dy_uz = computeYperiodicDerivative(*uz, domain);
    std::unique_ptr<MultiScalarField3D<T> > dz_uz = computeZperiodicDerivative(*uz, domain);

    std::unique_ptr<MultiScalarField3D<T> > term1 = multiply(*dx_ux, *dx_ux, domain);
    std::unique_ptr<MultiScalarField3D<T> > term2 = multiply(*dy_uy, *dy_uy, domain);
    std::unique_ptr<MultiScalarField3D<T> > term3 = multiply(*dz_uz, *dz_uz, domain);

    std::unique_ptr<MultiScalarField3D<T> > term4 =
        multiply((T)2, *multiply(*dx_uy, *dy_ux, domain), domain);
    std::unique_ptr<MultiScalarField3D<T> > term5 =
        multiply((T)2, *multiply(*dx_uz, *dz_ux, domain), domain);
    std::unique_ptr<MultiScalarField3D<T> > term6 =
        multiply((T)2, *multiply(*dy_uz, *dz_uy, domain), domain);

    std::unique_ptr<MultiScalarField3D<T> > rhs =
        add(*term6, *add(*term5, *add(*term4, *add(*term1, *add(*term2, *term3)))));
    return rhs;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computePeriodicPoissonRHS(
    MultiTensorField3D<T, 3> &velocity)
{
    return computePeriodicPoissonRHS(velocity, velocity.getBoundingBox());
}

template <typename T>
void periodicPoissonIterate(
    MultiScalarField3D<T> &oldPressure, MultiScalarField3D<T> &newPressure,
    MultiScalarField3D<T> &rhs, T beta, Box3D const &domain)
{
    std::vector<MultiScalarField3D<T> *> fields;
    fields.push_back(&oldPressure);
    fields.push_back(&newPressure);
    fields.push_back(&rhs);
    applyProcessingFunctional(new BoxPeriodicPoissonIteration3D<T>(beta), domain, fields);
}

// General Stencils.

template <typename T, int order, int maxWidth>
T computeScalarXderivative(
    ScalarField3D<T> const &scalar, int width, int position, plint iX, plint iY, plint iZ)
{
    T const *w = plb::fdWeights<T, order, maxWidth>().getWeights(width, position);
    T d = 0.0;
    for (plint i = 0, x = iX - (plint)position; i < (plint)width; i++, x++) {
        d += w[i] * scalar.get(x, iY, iZ);
    }
    return d;
}

template <typename T, int order, int maxWidth>
T computeScalarYderivative(
    ScalarField3D<T> const &scalar, int width, int position, plint iX, plint iY, plint iZ)
{
    T const *w = plb::fdWeights<T, order, maxWidth>().getWeights(width, position);
    T d = 0.0;
    for (plint i = 0, y = iY - (plint)position; i < (plint)width; i++, y++) {
        d += w[i] * scalar.get(iX, y, iZ);
    }
    return d;
}

template <typename T, int order, int maxWidth>
T computeScalarZderivative(
    ScalarField3D<T> const &scalar, int width, int position, plint iX, plint iY, plint iZ)
{
    T const *w = plb::fdWeights<T, order, maxWidth>().getWeights(width, position);
    T d = 0.0;
    for (plint i = 0, z = iZ - (plint)position; i < (plint)width; i++, z++) {
        d += w[i] * scalar.get(iX, iY, z);
    }
    return d;
}

template <typename T, int nDim, int order, int maxWidth>
Array<T, nDim> computeTensorXderivative(
    TensorField3D<T, nDim> const &tensor, int width, int position, plint iX, plint iY, plint iZ)
{
    T const *w = plb::fdWeights<T, order, maxWidth>().getWeights(width, position);
    Array<T, nDim> d;
    d.resetToZero();
    for (plint i = 0, x = iX - (plint)position; i < (plint)width; i++, x++) {
        d += w[i] * tensor.get(x, iY, iZ);
    }
    return d;
}

template <typename T, int nDim, int order, int maxWidth>
Array<T, nDim> computeTensorYderivative(
    TensorField3D<T, nDim> const &tensor, int width, int position, plint iX, plint iY, plint iZ)
{
    T const *w = plb::fdWeights<T, order, maxWidth>().getWeights(width, position);
    Array<T, nDim> d;
    d.resetToZero();
    for (plint i = 0, y = iY - (plint)position; i < (plint)width; i++, y++) {
        d += w[i] * tensor.get(iX, y, iZ);
    }
    return d;
}

template <typename T, int nDim, int order, int maxWidth>
Array<T, nDim> computeTensorZderivative(
    TensorField3D<T, nDim> const &tensor, int width, int position, plint iX, plint iY, plint iZ)
{
    T const *w = plb::fdWeights<T, order, maxWidth>().getWeights(width, position);
    Array<T, nDim> d;
    d.resetToZero();
    for (plint i = 0, z = iZ - (plint)position; i < (plint)width; i++, z++) {
        d += w[i] * tensor.get(iX, iY, z);
    }
    return d;
}

}  // namespace plb

#endif  // FINITE_DIFFERENCE_WRAPPER_3D_HH
