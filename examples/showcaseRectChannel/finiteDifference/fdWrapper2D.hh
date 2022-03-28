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

#ifndef FINITE_DIFFERENCE_WRAPPER_2D_HH
#define FINITE_DIFFERENCE_WRAPPER_2D_HH

#include "atomicBlock/dataProcessorWrapper2D.h"
#include "atomicBlock/reductiveDataProcessorWrapper2D.h"
#include "finiteDifference/fdFunctional2D.h"
#include "finiteDifference/fdWrapper2D.h"
#include "multiBlock/multiDataProcessorWrapper2D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper2D.h"
#include "multiGrid/gridConversion2D.h"

namespace plb {

template <typename T>
void computeLaplacian(
    MultiScalarField2D<T> &value, MultiScalarField2D<T> &laplacian, Box2D const &domain)
{
    applyProcessingFunctional(new BoxLaplacianFunctional2D<T>, domain, value, laplacian);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeLaplacian(
    MultiScalarField2D<T> &value, Box2D const &domain)
{
    MultiScalarField2D<T> *laplacian = new MultiScalarField2D<T>(value, domain);
    computeLaplacian(value, *laplacian, domain);
    return std::unique_ptr<MultiScalarField2D<T> >(laplacian);
}

template <typename T>
void computeXderivative(
    MultiScalarField2D<T> &value, MultiScalarField2D<T> &derivative, Box2D const &domain)
{
    plint boundaryWidth = 1;
    applyProcessingFunctional(
        new BoxXderivativeFunctional2D<T>, domain, value, derivative, boundaryWidth);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeXderivative(
    MultiScalarField2D<T> &value, Box2D const &domain)
{
    MultiScalarField2D<T> *derivative = new MultiScalarField2D<T>(value, domain);
    computeXderivative(value, *derivative, domain);
    return std::unique_ptr<MultiScalarField2D<T> >(derivative);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeXderivative(MultiScalarField2D<T> &value)
{
    return computeXderivative(value, value.getBoundingBox());
}

template <typename T>
void computeYderivative(
    MultiScalarField2D<T> &value, MultiScalarField2D<T> &derivative, Box2D const &domain)
{
    plint boundaryWidth = 1;
    applyProcessingFunctional(
        new BoxYderivativeFunctional2D<T>, domain, value, derivative, boundaryWidth);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeYderivative(
    MultiScalarField2D<T> &value, Box2D const &domain)
{
    MultiScalarField2D<T> *derivative = new MultiScalarField2D<T>(value, domain);
    computeYderivative(value, *derivative, domain);
    return std::unique_ptr<MultiScalarField2D<T> >(derivative);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeYderivative(MultiScalarField2D<T> &value)
{
    return computeYderivative(value, value.getBoundingBox());
}

template <typename T>
void computeGradientNorm(
    MultiScalarField2D<T> &value, MultiScalarField2D<T> &derivative, Box2D const &domain)
{
    plint boundaryWidth = 1;
    applyProcessingFunctional(
        new BoxGradientNormFunctional2D<T>, domain, value, derivative, boundaryWidth);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeGradientNorm(
    MultiScalarField2D<T> &value, Box2D const &domain)
{
    MultiScalarField2D<T> *derivative = new MultiScalarField2D<T>(value, domain);
    computeGradientNorm(value, *derivative, domain);
    return std::unique_ptr<MultiScalarField2D<T> >(derivative);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeGradientNorm(MultiScalarField2D<T> &value)
{
    return computeGradientNorm(value, value.getBoundingBox());
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computePoissonRHS(
    MultiTensorField2D<T, 2> &velocity, Box2D const &domain)
{
    std::unique_ptr<MultiScalarField2D<T> > ux = extractComponent(velocity, domain, 0);
    std::unique_ptr<MultiScalarField2D<T> > uy = extractComponent(velocity, domain, 1);

    std::unique_ptr<MultiScalarField2D<T> > dx_ux = computeXderivative(*ux, domain);
    std::unique_ptr<MultiScalarField2D<T> > dy_ux = computeYderivative(*ux, domain);
    std::unique_ptr<MultiScalarField2D<T> > dx_uy = computeXderivative(*uy, domain);
    std::unique_ptr<MultiScalarField2D<T> > dy_uy = computeYderivative(*uy, domain);

    std::unique_ptr<MultiScalarField2D<T> > term1 = multiply(*dx_ux, *dx_ux, domain);
    std::unique_ptr<MultiScalarField2D<T> > term2 =
        multiply((T)2, *multiply(*dx_uy, *dy_ux, domain), domain);
    std::unique_ptr<MultiScalarField2D<T> > term3 = multiply(*dy_uy, *dy_uy, domain);

    std::unique_ptr<MultiScalarField2D<T> > rhs = add(*term1, *add(*term2, *term3));
    return rhs;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computePoissonRHS(MultiTensorField2D<T, 2> &velocity)
{
    return computePoissonRHS(velocity, velocity.getBoundingBox());
}

template <typename T>
void poissonIterate(
    MultiScalarField2D<T> &oldPressure, MultiScalarField2D<T> &newPressure,
    MultiScalarField2D<T> &rhs, T beta, Box2D const &domain)
{
    std::vector<MultiScalarField2D<T> *> fields;
    fields.push_back(&oldPressure);
    fields.push_back(&newPressure);
    fields.push_back(&rhs);
    plint boundaryWidth = 1;
    applyProcessingFunctional(new BoxPoissonIteration2D<T>(beta), domain, fields, boundaryWidth);
}

template <typename T>
T computePoissonResidue(
    MultiScalarField2D<T> &pressure, MultiScalarField2D<T> &rhs, Box2D const &domain)
{
    BoxPoissonResidueFunctional2D<T> functional;
    applyProcessingFunctional(functional, domain, pressure, rhs);
    return functional.getMaxResidue();
}

/* ************ Wrapper for one Jacobi iteration *************** */
template <typename T>
void JacobiIteration(
    MultiScalarField2D<T> &u_h, MultiScalarField2D<T> &new_u_h, MultiScalarField2D<T> &rhs,
    Box2D const &domain)
{
    std::vector<MultiScalarField2D<T> *> fields;
    fields.push_back(&u_h);
    fields.push_back(&new_u_h);
    fields.push_back(&rhs);
    plint boundaryWidth = 1;
    applyProcessingFunctional(new JacobiIteration2D<T>(), domain, fields, boundaryWidth);
}

/* ************ Wrapper for one Gauss-Seidel iteration *************** */
template <typename T>
void GaussSeidelIteration(
    MultiScalarField2D<T> &u_h, MultiScalarField2D<T> &jacobi_u_h, MultiScalarField2D<T> &new_u_h,
    MultiScalarField2D<T> &rhs, Box2D const &domain)
{
    std::vector<MultiScalarField2D<T> *> fields;
    fields.push_back(&u_h);
    fields.push_back(&jacobi_u_h);
    fields.push_back(&new_u_h);
    fields.push_back(&rhs);
    plint boundaryWidth = 1;
    applyProcessingFunctional(new GaussSeidelIteration2D<T>(), domain, fields, boundaryWidth);
}

/* ************ Wrapper for Gauss-Seidel defect computation *************** */
template <typename T>
MultiScalarField2D<T> *computeGaussSeidelDefect(
    MultiScalarField2D<T> &u_h, MultiScalarField2D<T> &rhs, Box2D const &domain)
{
    MultiScalarField2D<T> *residual = new MultiScalarField2D<T>(u_h);
    std::vector<MultiScalarField2D<T> *> fields;
    fields.push_back(&u_h);
    fields.push_back(residual);
    fields.push_back(&rhs);
    plint boundaryWidth = 1;
    applyProcessingFunctional(new GaussSeidelDefect2D<T>(), domain, fields, boundaryWidth);

    return residual;
}

template <typename T>
T computeEuclidianNorm(MultiScalarField2D<T> &matrix, Box2D const &domain)
{
    T av = computeAverage(*multiply(matrix, matrix), domain);
    return std::sqrt(av);
}

/* ************ Gauss-Seidel Solver *************** */
template <typename T>
void GaussSeidelSolver(
    MultiScalarField2D<T> &initialValue, MultiScalarField2D<T> &result, MultiScalarField2D<T> &rhs,
    Box2D const &domain, T tolerance, plint maxIter)
{
    T originalNorm;
    T newNorm;
    plint infoIt = 100;

    // to contain the Jacobi result
    MultiScalarField2D<T> jacobiValue(initialValue);
    jacobiValue.reset();

    // computation of the initial residual. We will stop when
    // currentResidual=tolerance*initialResidual or when we have made maxIter iterations
    MultiScalarField2D<T> *defect = computeGaussSeidelDefect<T>(initialValue, rhs, domain);
    originalNorm = computeEuclidianNorm<T>(*defect, domain);
    delete defect;

    // MAIN LOOP
    for (plint iT = 0; iT < maxIter; ++iT) {
        // use one Jacobi iteration
        JacobiIteration<T>(initialValue, jacobiValue, rhs, domain);

        // with this jacobiValue, make one gauss-sidel iteration (eliminate the asynchronie for
        // parallelism)
        GaussSeidelIteration<T>(initialValue, jacobiValue, result, rhs, domain);

        // compute the new residual norm to know if we continue
        defect = computeGaussSeidelDefect<T>(result, rhs, domain);
        newNorm = computeEuclidianNorm<T>(*defect, domain);
        delete defect;

        if (newNorm < tolerance * originalNorm) {
            pcout << "Gauss-Seidel iterations: " << iT << std::endl;
            break;
        }

        // result becomes the new initial value
        result.swap(initialValue);

        if (iT % infoIt == 0 && iT > 0) {
            pcout << "Gauss-Seidel iteration " << iT << " : error=" << newNorm / originalNorm
                  << " : NORM=" << newNorm << std::endl;
        }
    }
    pcout << "Gauss-Seidel iterations: " << maxIter << std::endl;
}

/* ************ MultiGrid method *************** */
/// Iterate Gauss-Sidel for a given number of iterations
template <typename T>
MultiScalarField2D<T> *smooth(
    MultiScalarField2D<T> &initialValue, MultiScalarField2D<T> &rhs, Box2D const &domain,
    plint smoothIters)
{
    // create the solution for Jacobi
    MultiScalarField2D<T> jacobiValue(initialValue);
    jacobiValue.reset();

    // copy initialValue
    MultiScalarField2D<T> initialCopy(initialValue);
    // the result holder
    MultiScalarField2D<T> *newValue = new MultiScalarField2D<T>(initialValue);

    // iteration over the original system: one Jacobi + one Gauss-Seidel
    plint v1 = smoothIters;
    for (plint iV1 = 0; iV1 < v1; ++iV1) {
        newValue->swap(initialCopy);
        JacobiIteration<T>(initialCopy, jacobiValue, rhs, domain);
        GaussSeidelIteration<T>(initialCopy, jacobiValue, *newValue, rhs, domain);
    }

    return newValue;
}

/// Iterate Gauss-Sidel and then interpolate the result (go up in the V)
template <typename T>
MultiScalarField2D<T> *smoothAndInterpolate(
    MultiScalarField2D<T> &initialValue, MultiScalarField2D<T> &rhs, Box2D const &domain,
    plint smoothIters)
{
    // smooth
    MultiScalarField2D<T> *coarseNewValue = smooth(initialValue, rhs, domain, smoothIters);
    // interpolate
    MultiScalarField2D<T> *fineNewValue =
        new MultiScalarField2D<T>(*refine<T>(*coarseNewValue, 1, -1, -1, -1));
    delete coarseNewValue;

    return fineNewValue;
}

/// Iterate Gauss-Sidel and then compute the error (the finest level in the V)
template <typename T>
MultiScalarField2D<T> *smoothAndComputeError(
    MultiScalarField2D<T> &initialValue, MultiScalarField2D<T> &rhs, Box2D const &domain,
    T &newError, plint smoothIters)
{
    MultiScalarField2D<T> *newValue = smooth(initialValue, rhs, domain, smoothIters);
    // compute the new defect norm
    MultiScalarField2D<T> *newDefect = computeGaussSeidelDefect<T>(*newValue, rhs, domain);
    newError = computeEuclidianNorm<T>(*newDefect, domain);
    delete newDefect;

    return newValue;
}

/// Iterate Gauss-Sidel and then compute the defect (the last level before the coarsest grid)
template <typename T>
MultiScalarField2D<T> *smoothAndComputeCoarseDefect(
    MultiScalarField2D<T> &initialValue, MultiScalarField2D<T> &rhs, Box2D const &domain,
    plint smoothIters)
{
    MultiScalarField2D<T> *newValue = smooth(initialValue, rhs, domain, smoothIters);
    // compute the new defect norm
    MultiScalarField2D<T> *newDefect = computeGaussSeidelDefect<T>(*newValue, rhs, domain);
    MultiScalarField2D<T> *coarseDefect =
        new MultiScalarField2D<T>(*coarsen<T>(*newDefect, 1, -1, 1, 1));
    delete newDefect;

    return coarseDefect;
}

template <typename T>
void generateRHS(MultiScalarField2D<T> &originalRHS, std::vector<MultiScalarField2D<T> *> &rhs)
{
    plint vectorSize = (plint)rhs.size();
    rhs[vectorSize - 1] = new MultiScalarField2D<T>(originalRHS);
    for (plint iLevel = vectorSize - 2; iLevel >= 0; iLevel--) {
        rhs[iLevel] = new MultiScalarField2D<T>(*coarsen<T>(*rhs[iLevel + 1], 1, -1, 1, 1));
    }
}
// TODO: This "domain" is misleading. Unused.
template <typename T>
T multiGridVCycle(
    MultiScalarField2D<T> &initialValue, MultiScalarField2D<T> &newValue,
    MultiScalarField2D<T> &rhs, Box2D const &domain, plint depth)
{
    PLB_PRECONDITION(depth >= 1);

    // containers of the values to not lose anything
    std::vector<MultiScalarField2D<T> *> newValues(depth + 1);
    std::vector<MultiScalarField2D<T> *> defects(depth + 1);

    plint smoothIters1 = 5;
    plint smoothIters2 = 5;
    // we go down until the level before the coarsest
    newValues[depth] = smooth<T>(initialValue, rhs, initialValue.getBoundingBox(), smoothIters1);
    defects[depth] =
        computeGaussSeidelDefect<T>(*newValues[depth], rhs, initialValue.getBoundingBox());
    defects[depth - 1] = new MultiScalarField2D<T>(*coarsen<T>(*defects[depth], 1, -1, 1, 1));
    *defects[depth - 1] = *plb::multiply(-1.0, *defects[depth - 1]);

    /// GOING DOWN THE V
    for (plint iLevel = depth - 1; iLevel >= 1; iLevel--) {
        // create initial solution = 0
        MultiScalarField2D<T> initialSolution(*defects[iLevel]);
        initialSolution.reset();
        // smooth the new system
        newValues[iLevel] = smooth<T>(
            initialSolution, *defects[iLevel], initialSolution.getBoundingBox(), smoothIters1);
        MultiScalarField2D<T> *tempDefect = computeGaussSeidelDefect<T>(
            *newValues[iLevel], *defects[iLevel], newValues[iLevel]->getBoundingBox());
        defects[iLevel - 1] = new MultiScalarField2D<T>(*coarsen<T>(*tempDefect, 1, -1, 1, 1));
        *defects[iLevel - 1] = *plb::multiply(-1.0, *defects[iLevel - 1]);  // change the sign
        delete tempDefect;
    }

    /// THE EDGE OF THE V
    // resolve the system exactly for the coarse system
    MultiScalarField2D<T> initialValueCoarse(*defects[0]);
    MultiScalarField2D<T> resultCoarse(
        *defects[0]);            // container for the coarse solution to the correction scheme
    initialValueCoarse.reset();  // for the correction 0 is a good first approximation

    T tolerance = 1e-5;
    plint maxIter = 100;
    GaussSeidelSolver<T>(
        initialValueCoarse, resultCoarse, *defects[0], resultCoarse.getBoundingBox(), tolerance,
        maxIter);

    newValues[0] = new MultiScalarField2D<T>(resultCoarse);
    // we interpolate the error computed in the coarse grid and add it to the original result from
    // level 1
    MultiScalarField2D<T> *v_h =
        new MultiScalarField2D<T>(*refine<T>(*newValues[0], 1, -1, -1, -1));
    // compute new approximation as newValue = u_h + v_h
    *newValues[1] = *plb::add(*newValues[1], *v_h);
    delete v_h;

    /// GOING UP THE V
    // go up interpolating and smoothing up to the level before the finest
    for (plint iLevel = 1; iLevel < depth; ++iLevel) {
        MultiScalarField2D<T> *tempNewVal = smoothAndInterpolate<T>(
            *newValues[iLevel], *defects[iLevel], newValues[iLevel]->getBoundingBox(),
            smoothIters2);
        *newValues[iLevel + 1] = *plb::add(*newValues[iLevel + 1], *tempNewVal);
        delete tempNewVal;
    }

    // we smooth and compute the new error in the finest level
    T newError = 0.0;
    MultiScalarField2D<T> *result = smoothAndComputeError(
        *newValues[depth], rhs, newValues[depth]->getBoundingBox(), newError, smoothIters2);
    // copy new initial value
    newValue = *result;
    delete result;

    // CLEAN-UP
    for (plint iLevel = 0; iLevel <= depth; ++iLevel) {
        delete newValues[iLevel];
        delete defects[iLevel];
    }

    return newError;
}
// TODO: Unused initialValue. We should remove it.
template <typename T>
std::vector<MultiScalarField2D<T> *> fullMultiGrid(
    MultiScalarField2D<T> &initialValue, MultiScalarField2D<T> &originalRhs, Box2D const &domain,
    plint gridLevels, plint ncycles)
{
    PLB_PRECONDITION(gridLevels >= 2 && ncycles >= 1);

    std::vector<MultiScalarField2D<T> *> rhs(gridLevels);
    std::vector<MultiScalarField2D<T> *> solutions(gridLevels);

    // create all the right-hand sides from restriction over rhs
    generateRHS(originalRhs, rhs);

    // Initial solution computed exactly on the coarsest grid
    MultiScalarField2D<T> initialSolution(*rhs[0]);
    MultiScalarField2D<T> initialValueCoarse(*rhs[0]);
    initialValueCoarse.reset();
    GaussSeidelSolver<T>(
        initialValueCoarse, initialSolution, *rhs[0], initialSolution.getBoundingBox());
    solutions[0] = new MultiScalarField2D<T>(initialSolution);

    // MAIN LOOP (at each time we add one more level)
    for (plint iLevel = 1; iLevel < gridLevels; ++iLevel) {
        pcout << "Level: " << iLevel << std::endl;
        // interpolate the iLevel-1 solution
        MultiScalarField2D<T> *initialSolutionLevel =
            new MultiScalarField2D<T>(*refine<T>(*solutions[iLevel - 1], 1, -1, -1, -1));
        MultiScalarField2D<T> newValue(*initialSolutionLevel);
        // for each level we itere ncycles
        for (plint iNcycle = 0; iNcycle < ncycles * iLevel; ++iNcycle) {
            // we use a V cycle with the current number of grids
            T error;
            error =
                multiGridVCycle<T>(*initialSolutionLevel, newValue, *rhs[iLevel], domain, iLevel);
            pcout << "\tError=" << error << std::endl;
            newValue.swap(*initialSolutionLevel);
        }

        // save the current level solution
        solutions[iLevel] = new MultiScalarField2D<T>(*initialSolutionLevel);
        delete initialSolutionLevel;
    }

    // CLEAN-UP
    for (plint iLevel = 0; iLevel < gridLevels; ++iLevel) {
        delete rhs[iLevel];
    }

    return solutions;
}
// TODO: Unused initalValue. Misleading. We should remove it. Same for domain.
template <typename T>
std::vector<MultiScalarField2D<T> *> simpleMultiGrid(
    MultiScalarField2D<T> &initialValue, MultiScalarField2D<T> &originalRhs, Box2D const &domain,
    plint gridLevels)
{
    PLB_PRECONDITION(gridLevels >= 2);

    std::vector<MultiScalarField2D<T> *> rhs(gridLevels);
    std::vector<MultiScalarField2D<T> *> solutions(gridLevels);

    // create all the right-hand sides from restriction over rhs
    generateRHS(originalRhs, rhs);

    // Initial solution computed exactly on the coarsest grid
    pcout << "Level: 0" << std::endl;
    pcout << "Size of the field : " << rhs[0]->getNx() << " x " << rhs[0]->getNy() << std::endl;
    MultiScalarField2D<T> initialSolution(*rhs[0]);
    MultiScalarField2D<T> initialValueCoarse(*rhs[0]);
    initialValueCoarse.reset();
    GaussSeidelSolver<T>(
        initialValueCoarse, initialSolution, *rhs[0], initialSolution.getBoundingBox());
    solutions[0] = new MultiScalarField2D<T>(initialSolution);

    // MAIN LOOP (at each time we add one more level)
    for (plint iLevel = 1; iLevel < gridLevels; ++iLevel) {
        pcout << "Level: " << iLevel << std::endl;

        // interpolate the iLevel-1 solution
        MultiScalarField2D<T> *initialSolutionLevel =
            new MultiScalarField2D<T>(*refine<T>(*solutions[iLevel - 1], 1, -1, -1, -1));
        MultiScalarField2D<T> newValue(*initialSolutionLevel);
        pcout << "Size of the field : " << newValue.getNx() << " x " << newValue.getNy()
              << std::endl;
        GaussSeidelSolver<T>(
            *initialSolutionLevel, newValue, *rhs[iLevel], initialSolutionLevel->getBoundingBox());

        // save the current level solution
        solutions[iLevel] = new MultiScalarField2D<T>(*initialSolutionLevel);
        delete initialSolutionLevel;
    }

    // CLEAN-UP
    for (plint iLevel = 0; iLevel < gridLevels; ++iLevel) {
        delete rhs[iLevel];
    }

    return solutions;
}

// General Stencils.

template <typename T, int order, int maxWidth>
T computeScalarXderivative(
    ScalarField2D<T> const &scalar, int width, int position, plint iX, plint iY)
{
    T const *w = plb::fdWeights<T, order, maxWidth>().getWeights(width, position);
    T d = 0.0;
    for (plint i = 0, x = iX - (plint)position; i < (plint)width; i++, x++) {
        d += w[i] * scalar.get(x, iY);
    }
    return d;
}

template <typename T, int order, int maxWidth>
T computeScalarYderivative(
    ScalarField2D<T> const &scalar, int width, int position, plint iX, plint iY)
{
    T const *w = plb::fdWeights<T, order, maxWidth>().getWeights(width, position);
    T d = 0.0;
    for (plint i = 0, y = iY - (plint)position; i < (plint)width; i++, y++) {
        d += w[i] * scalar.get(iX, y);
    }
    return d;
}

template <typename T, int nDim, int order, int maxWidth>
Array<T, nDim> computeTensorXderivative(
    TensorField2D<T, nDim> const &tensor, int width, int position, plint iX, plint iY)
{
    T const *w = plb::fdWeights<T, order, maxWidth>().getWeights(width, position);
    Array<T, nDim> d;
    d.resetToZero();
    for (plint i = 0, x = iX - (plint)position; i < (plint)width; i++, x++) {
        d += w[i] * tensor.get(x, iY);
    }
    return d;
}

template <typename T, int nDim, int order, int maxWidth>
Array<T, nDim> computeTensorYderivative(
    TensorField2D<T, nDim> const &tensor, int width, int position, plint iX, plint iY)
{
    T const *w = plb::fdWeights<T, order, maxWidth>().getWeights(width, position);
    Array<T, nDim> d;
    d.resetToZero();
    for (plint i = 0, y = iY - (plint)position; i < (plint)width; i++, y++) {
        d += w[i] * tensor.get(iX, y);
    }
    return d;
}

}  // namespace plb

#endif  // FINITE_DIFFERENCE_WRAPPER_2D_HH
