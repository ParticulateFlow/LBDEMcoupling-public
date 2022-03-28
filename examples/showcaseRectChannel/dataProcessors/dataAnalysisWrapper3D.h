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
 * Helper functions for data analysis -- header file.
 */

#ifndef DATA_ANALYSIS_WRAPPER_3D_H
#define DATA_ANALYSIS_WRAPPER_3D_H

#include <memory>

#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "core/globalDefs.h"
#include "dataProcessors/dataAnalysisFunctional3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"

namespace plb {

/* ******************************************************************* */
/* *************** PART I. Atomic-block wrappers: Block-Lattice ****** */
/* ******************************************************************* */

/* *************** Reductive operations ****************************** */

template <typename T, template <typename U> class Descriptor>
T computeAverageDensity(BlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
T computeAverageDensity(BlockLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(BlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(BlockLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
T computeAverageEnergy(BlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
T computeAverageEnergy(BlockLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor, class BoolMask>
plint count(BlockLattice3D<T, Descriptor> &lattice, Box3D domain, BoolMask boolMask);

template <typename T, template <typename U> class Descriptor, class BoolMask>
plint count(BlockLattice3D<T, Descriptor> &lattice, BoolMask boolMask);

/* *************** Density ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeDensity(BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &density);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField3D<T> > computeDensity(BlockLattice3D<T, Descriptor> &lattice);

/* *************** RhoBar ******************************************** */

template <typename T, template <typename U> class Descriptor>
void computeRhoBar(BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &rhoBar);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField3D<T> > computeRhoBar(BlockLattice3D<T, Descriptor> &lattice);

/* *************** Kinetic Energy ************************************ */

template <typename T, template <typename U> class Descriptor>
void computeKineticEnergy(BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &energy);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField3D<T> > computeKineticEnergy(BlockLattice3D<T, Descriptor> &lattice);

/* *************** Velocity Norm ************************************* */

template <typename T, template <typename U> class Descriptor>
void computeVelocityNorm(BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &velocityNorm);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField3D<T> > computeVelocityNorm(BlockLattice3D<T, Descriptor> &lattice);

/* *************** Velocity Component ************************************ */

template <typename T, template <typename U> class Descriptor>
void computeVelocityComponent(
    BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &velocityComponent, plint iComponent);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField3D<T> > computeVelocityComponent(BlockLattice3D<T, Descriptor> &lattice);

/* *************** Velocity ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocity(
    BlockLattice3D<T, Descriptor> &lattice, TensorField3D<T, Descriptor<T>::d> &velocity);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<TensorField3D<T, Descriptor<T>::d> > computeVelocity(
    BlockLattice3D<T, Descriptor> &lattice);

/* *************** Pi Neq ********************************* */

template <typename T, template <typename U> class Descriptor>
void computePiNeq(
    BlockLattice3D<T, Descriptor> &lattice,
    TensorField3D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<TensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computePiNeq(
    BlockLattice3D<T, Descriptor> &lattice);

/* *************** Shear Stress ********************************* */

template <typename T, template <typename U> class Descriptor>
void computeShearStress(
    BlockLattice3D<T, Descriptor> &lattice,
    TensorField3D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<TensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computeShearStress(
    BlockLattice3D<T, Descriptor> &lattice);

/* *************** Strain Rate from Stress *************************** */

template <typename T, template <typename U> class Descriptor>
void computeStrainRateFromStress(
    BlockLattice3D<T, Descriptor> &lattice, TensorField3D<T, SymmetricTensor<T, Descriptor>::n> &S);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<TensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computeStrainRateFromStress(
    BlockLattice3D<T, Descriptor> &lattice);

/* *************** Population **************************************** */

template <typename T, template <typename U> class Descriptor>
void computePopulation(
    BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &population, plint iPop);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField3D<T> > computePopulation(
    BlockLattice3D<T, Descriptor> &lattice, plint iPop);

template <typename T, template <typename U> class Descriptor>
void computeAllPopulations(
    BlockLattice3D<T, Descriptor> &lattice, TensorField3D<T, Descriptor<T>::q> &populations);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<TensorField3D<T, Descriptor<T>::q> > computeAllPopulations(
    BlockLattice3D<T, Descriptor> &lattice);

/* ******************************************************************* */
/* *************** PART II. Atomic-block wrappers: Scalar-Field ****** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template <typename T>
T computeSum(ScalarField3D<T> &scalarField, Box3D domain);

template <typename T>
T computeSum(ScalarField3D<T> &scalarField);

template <typename T>
T computeSum(ScalarField3D<T> &scalarField, ScalarField3D<int> &mask, int flag, Box3D domain);

template <typename T>
T computeSum(ScalarField3D<T> &scalarField, ScalarField3D<int> &mask, int flag);

template <typename T>
T computeAverage(ScalarField3D<T> &scalarField, Box3D domain);

template <typename T>
T computeAverage(ScalarField3D<T> &scalarField);

template <typename T>
T computeAverage(ScalarField3D<T> &scalarField, ScalarField3D<int> &mask, int flag, Box3D domain);

template <typename T>
T computeAverage(ScalarField3D<T> &scalarField, ScalarField3D<int> &mask, int flag);

template <typename T>
T computeMin(ScalarField3D<T> &scalarField, Box3D domain);

template <typename T>
T computeMin(ScalarField3D<T> &scalarField);

template <typename T>
T computeMax(ScalarField3D<T> &scalarField, Box3D domain);

template <typename T>
T computeMax(ScalarField3D<T> &scalarField);

template <typename T>
T computeBoundedAverage(ScalarField3D<T> &scalarField, Box3D domain);

template <typename T>
T computeBoundedAverage(ScalarField3D<T> &scalarField);

template <typename T, class BoolMask>
plint count(ScalarField3D<T> &field, Box3D domain, BoolMask boolMask);

template <typename T, class BoolMask>
plint count(ScalarField3D<T> &field, BoolMask boolMask);

/* *************** Generic operations *************** */

template <typename T, class Function>
void apply(Function f, ScalarField3D<T> &field)
{
    applyProcessingFunctional(
        new ApplyScalarFunctional3D<T, Function>(f), field.getBoundingBox(), field);
}

template <typename T, class Function>
void evaluate(Function f, ScalarField3D<T> &field, ScalarField3D<T> &result)
{
    applyProcessingFunctional(
        new EvaluateScalarFunctional3D<T, Function>(f), field.getBoundingBox(), field, result);
}

template <typename T, class Function>
std::unique_ptr<ScalarField3D<T> > evaluate(Function f, ScalarField3D<T> &field)
{
    ScalarField3D<T> *result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    evaluate(f, field, result);
    return std::unique_ptr<ScalarField3D<T> >(result);
}

/* *************** ScalarField - Scalar operations *************** */

template <typename T>
void add(T scalar, ScalarField3D<T> &field, ScalarField3D<T> &result);

template <typename T>
std::unique_ptr<ScalarField3D<T> > add(T scalar, ScalarField3D<T> &field);

template <typename T>
void add(ScalarField3D<T> &field, T scalar, ScalarField3D<T> &result);

template <typename T>
std::unique_ptr<ScalarField3D<T> > add(ScalarField3D<T> &field, T scalar);

template <typename T>
void subtract(T scalar, ScalarField3D<T> &field, ScalarField3D<T> &result);

template <typename T>
std::unique_ptr<ScalarField3D<T> > subtract(T scalar, ScalarField3D<T> &field);

template <typename T>
void subtract(ScalarField3D<T> &field, T scalar, ScalarField3D<T> &result);

template <typename T>
std::unique_ptr<ScalarField3D<T> > subtract(ScalarField3D<T> &field, T scalar);

template <typename T>
void multiply(T scalar, ScalarField3D<T> &field, ScalarField3D<T> &result);

template <typename T>
std::unique_ptr<ScalarField3D<T> > multiply(T scalar, ScalarField3D<T> &field);

template <typename T>
void multiply(ScalarField3D<T> &field, T scalar, ScalarField3D<T> &result);

template <typename T>
std::unique_ptr<ScalarField3D<T> > multiply(ScalarField3D<T> &field, T scalar);

template <typename T>
void divide(T scalar, ScalarField3D<T> &field, ScalarField3D<T> &result);

template <typename T>
std::unique_ptr<ScalarField3D<T> > divide(T scalar, ScalarField3D<T> &field);

template <typename T>
void divide(ScalarField3D<T> &field, T scalar, ScalarField3D<T> &result);

template <typename T>
std::unique_ptr<ScalarField3D<T> > divide(ScalarField3D<T> &field, T scalar);

/* *************** ScalarField - Scalar inplace operations *************** */

template <typename T>
void addInPlace(ScalarField3D<T> &field, T scalar);

template <typename T>
void subtractInPlace(ScalarField3D<T> &field, T scalar);

template <typename T>
void multiplyInPlace(ScalarField3D<T> &field, T scalar);

template <typename T>
void divideInPlace(ScalarField3D<T> &field, T scalar);

/* *************** ScalarField - ScalarField operations *************** */

template <typename T1, typename T2>
void copy(ScalarField3D<T1> &field, ScalarField3D<T2> &convertedField);

template <typename T1, typename T2>
std::unique_ptr<ScalarField3D<T2> > copyConvert(ScalarField3D<T1> &field);

template <typename T>
void add(ScalarField3D<T> &A, ScalarField3D<T> &B, ScalarField3D<T> &result);

template <typename T>
std::unique_ptr<ScalarField3D<T> > add(ScalarField3D<T> &A, ScalarField3D<T> &B);

template <typename T>
void subtract(ScalarField3D<T> &A, ScalarField3D<T> &B, ScalarField3D<T> &result);

template <typename T>
std::unique_ptr<ScalarField3D<T> > subtract(ScalarField3D<T> &A, ScalarField3D<T> &B);

template <typename T>
void multiply(ScalarField3D<T> &A, ScalarField3D<T> &B, ScalarField3D<T> &result);

template <typename T>
std::unique_ptr<ScalarField3D<T> > multiply(ScalarField3D<T> &A, ScalarField3D<T> &B);

template <typename T>
void divide(ScalarField3D<T> &A, ScalarField3D<T> &B, ScalarField3D<T> &result);

template <typename T>
std::unique_ptr<ScalarField3D<T> > divide(ScalarField3D<T> &A, ScalarField3D<T> &B);

/* *************** ScalarField operations *************** */

template <typename T>
void computeSqrt(ScalarField3D<T> &field, ScalarField3D<T> &result, Box3D domain);

template <typename T>
std::unique_ptr<ScalarField3D<T> > computeSqrt(ScalarField3D<T> &field, Box3D domain);

template <typename T>
std::unique_ptr<ScalarField3D<T> > computeSqrt(ScalarField3D<T> &field);

template <typename T>
void computeAbsoluteValue(ScalarField3D<T> &field, ScalarField3D<T> &result, Box3D domain);

template <typename T>
std::unique_ptr<ScalarField3D<T> > computeAbsoluteValue(ScalarField3D<T> &field, Box3D domain);

template <typename T>
std::unique_ptr<ScalarField3D<T> > computeAbsoluteValue(ScalarField3D<T> &field);

/* *************** ScalarField - ScalarField inplace operations *************** */

template <typename T>
void addInPlace(ScalarField3D<T> &A, ScalarField3D<T> &B);

template <typename T>
void subtractInPlace(ScalarField3D<T> &A, ScalarField3D<T> &B);

template <typename T>
void multiplyInPlace(ScalarField3D<T> &A, ScalarField3D<T> &B);

template <typename T>
void divideInPlace(ScalarField3D<T> &A, ScalarField3D<T> &B);

/* ******************************************************************* */
/* *************** PART III. Atomic-block wrappers: Tensor-Field ***** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template <typename T, int nDim, class BoolMask>
plint count(TensorField3D<T, nDim> &field, Box3D domain, BoolMask boolMask);

template <typename T, int nDim, class BoolMask>
plint count(TensorField3D<T, nDim> &field, BoolMask boolMask);

template <typename T, int nDim>
Array<T, nDim> computeSum(TensorField3D<T, nDim> &tensorField, Box3D domain);

template <typename T, int nDim>
Array<T, nDim> computeSum(TensorField3D<T, nDim> &tensorField);

template <typename T, int nDim>
Array<T, nDim> computeSum(
    TensorField3D<T, nDim> &tensorField, ScalarField3D<int> &mask, int flag, Box3D domain);

template <typename T, int nDim>
Array<T, nDim> computeSum(TensorField3D<T, nDim> &tensorField, ScalarField3D<int> &mask, int flag);

template <typename T, int nDim>
Array<T, nDim> computeAverage(TensorField3D<T, nDim> &tensorField, Box3D domain);

template <typename T, int nDim>
Array<T, nDim> computeAverage(TensorField3D<T, nDim> &tensorField);

template <typename T, int nDim>
Array<T, nDim> computeAverage(
    TensorField3D<T, nDim> &tensorField, ScalarField3D<int> &mask, int flag, Box3D domain);

template <typename T, int nDim>
Array<T, nDim> computeAverage(
    TensorField3D<T, nDim> &tensorField, ScalarField3D<int> &mask, int flag);

/* *************** Component (scalar-field) out of a tensor-field ****** */

template <typename T, int nDim>
void extractComponent(
    TensorField3D<T, nDim> &tensorField, ScalarField3D<T> &component, int iComponent);

template <typename T, int nDim>
std::unique_ptr<ScalarField3D<T> > extractComponent(
    TensorField3D<T, nDim> &tensorField, int iComponent);

/* *************** Copy-convert a tensor-field *************** */

template <typename T1, typename T2, int nDim>
void copy(TensorField3D<T1, nDim> &field, TensorField3D<T2, nDim> &convertedField);

template <typename T1, typename T2, int nDim>
std::unique_ptr<TensorField3D<T2, nDim> > copyConvert(TensorField3D<T1, nDim> &field);

/* *************** Vector-norm of each cell in the field *************** */

template <typename T, int nDim>
void computeNorm(TensorField3D<T, nDim> &tensorField, ScalarField3D<T> &norm);

template <typename T, int nDim>
std::unique_ptr<ScalarField3D<T> > computeNorm(TensorField3D<T, nDim> &norm);

/* *************** Squared vector-norm of each cell in the field ******** */

template <typename T, int nDim>
void computeNormSqr(TensorField3D<T, nDim> &tensorField, ScalarField3D<T> &norm);

template <typename T, int nDim>
std::unique_ptr<ScalarField3D<T> > computeNormSqr(TensorField3D<T, nDim> &norm);

/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template <typename T>
void computeSymmetricTensorNorm(TensorField3D<T, 6> &tensorField, ScalarField3D<T> &norm);

template <typename T>
std::unique_ptr<ScalarField3D<T> > computeSymmetricTensorNorm(TensorField3D<T, 6> &norm);

/* *************** Squared Tensor-norm of each symmetric tensor of a field ************* */

template <typename T>
void computeSymmetricTensorNormSqr(TensorField3D<T, 6> &tensorField, ScalarField3D<T> &normSqr);

template <typename T>
std::unique_ptr<ScalarField3D<T> > computeSymmetricTensorNormSqr(TensorField3D<T, 6> &normSqr);

/* *************** Trace of each symmetric tensor of a field ************ */

template <typename T>
void computeSymmetricTensorTrace(TensorField3D<T, 6> &tensorField, ScalarField3D<T> &trace);

template <typename T>
std::unique_ptr<ScalarField3D<T> > computeSymmetricTensorTrace(TensorField3D<T, 6> &trace);

/* *************** Trace of each symmetric tensor of a field ************ */

template <typename T>
void computeSymmetricTensorTrace(TensorField3D<T, 6> &tensorField, ScalarField3D<T> &trace);

template <typename T>
std::unique_ptr<ScalarField3D<T> > computeSymmetricTensorTrace(TensorField3D<T, 6> &trace);

/* *************** Vorticity from Velocity field *********************** */

template <typename T>
void computeVorticity(TensorField3D<T, 3> &velocity, TensorField3D<T, 3> &vorticity);

template <typename T>
std::unique_ptr<TensorField3D<T, 3> > computeVorticity(TensorField3D<T, 3> &velocity);

/* *************** Vorticity, witout boundary treatment, from Velocity field ************** */

template <typename T>
void computeBulkVorticity(TensorField3D<T, 3> &velocity, TensorField3D<T, 3> &vorticity);

template <typename T>
std::unique_ptr<TensorField3D<T, 3> > computeBulkVorticity(TensorField3D<T, 3> &velocity);

/* *************** Helicity from Velocity field *********************** */

template <typename T>
void computeHelicity(TensorField3D<T, 3> &velocity, ScalarField3D<T> &helicity);

template <typename T>
std::unique_ptr<ScalarField3D<T> > computeHelicity(TensorField3D<T, 3> &velocity);

/* *************** Helicity, witout boundary treatment, from Velocity field ************** */

template <typename T>
void computeBulkHelicity(TensorField3D<T, 3> &velocity, ScalarField3D<T> &helicity);

template <typename T>
std::unique_ptr<ScalarField3D<T> > computeBulkHelicity(TensorField3D<T, 3> &velocity);

/* *************** Divergence, witout boundary treatment, from Velocity field ************* */

template <typename T>
void computeBulkDivergence(TensorField3D<T, 3> &velocity, ScalarField3D<T> &divergence);

template <typename T>
std::unique_ptr<ScalarField3D<T> > computeBulkDivergence(TensorField3D<T, 3> &velocity);

/* *************** Strain rate from Velocity field ********************* */

template <typename T>
void computeStrainRate(TensorField3D<T, 3> &velocity, TensorField3D<T, 6> &S);

template <typename T>
std::unique_ptr<TensorField3D<T, 6> > computeStrainRate(TensorField3D<T, 3> &velocity);

/* *************** Strain rate, witout boundary treatment, from Velocity field ************ */

template <typename T>
void computeBulkStrainRate(TensorField3D<T, 3> &velocity, TensorField3D<T, 6> &S);

template <typename T>
std::unique_ptr<TensorField3D<T, 6> > computeBulkStrainRate(TensorField3D<T, 3> &velocity);

/* *************** TensorField - TensorField operations *************** */

template <typename T, int nDim>
void add(TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B, TensorField3D<T, nDim> &result);

template <typename T, int nDim>
std::unique_ptr<TensorField3D<T, nDim> > add(TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B);

template <typename T, int nDim>
void subtract(TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B, TensorField3D<T, nDim> &result);

template <typename T, int nDim>
std::unique_ptr<TensorField3D<T, nDim> > subtract(
    TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B);

template <typename T, int nDim>
void multiply(TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B, TensorField3D<T, nDim> &result);

template <typename T, int nDim>
std::unique_ptr<TensorField3D<T, nDim> > multiply(
    TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B);

template <typename T, int nDim>
void multiply(
    T scalar, TensorField3D<T, nDim> &field, TensorField3D<T, nDim> &result, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<TensorField3D<T, nDim> > multiply(T scalar, TensorField3D<T, nDim> &field);

template <typename T, int nDim>
void multiply(
    TensorField3D<T, nDim> &field, T scalar, TensorField3D<T, nDim> &result, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<TensorField3D<T, nDim> > multiply(TensorField3D<T, nDim> &field, T scalar);

template <typename T, int nDim>
void divide(TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B, TensorField3D<T, nDim> &result);

template <typename T, int nDim>
std::unique_ptr<TensorField3D<T, nDim> > divide(
    TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B);

/* *************** TensorField operations *************** */

template <typename T, int nDim>
void computeSqrt(TensorField3D<T, nDim> &field, TensorField3D<T, nDim> &result, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<TensorField3D<T, nDim> > computeSqrt(TensorField3D<T, nDim> &field, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<TensorField3D<T, nDim> > computeSqrt(TensorField3D<T, nDim> &field);

/* *************** TensorField - TensorField inplace operations *************** */

template <typename T, int nDim>
void addInPlace(TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B);

template <typename T, int nDim>
void subtractInPlace(TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B);

template <typename T, int nDim>
void multiplyInPlace(TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B);

template <typename T, int nDim>
void multiplyInPlace(TensorField3D<T, nDim> &A, T alpha);

template <typename T, int nDim>
void divideInPlace(TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B);

/* ******************************************************************* */
/* *************** PART IV. Multi-block wrappers: Block-Lattice ****** */
/* ******************************************************************* */

/* *************** Reductive Functions ******************************* */

template <typename T, template <typename U> class Descriptor>
T computeAverageDensity(MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
T computeAverageDensity(MultiBlockLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(MultiBlockLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
T computeAverageEnergy(MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
T computeAverageEnergy(MultiBlockLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
T computeAverageForcedEnergy(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    Box3D domain);

template <typename T, template <typename U> class Descriptor>
T computeAverageForcedEnergy(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force);

template <typename T, template <typename U> class Descriptor>
T computeAverageForcedEnergy(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force, Box3D domain);

template <typename T, template <typename U> class Descriptor>
T computeAverageForcedEnergy(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force);

template <typename T, template <typename U> class Descriptor, class ForceFunction>
T computeAverageForcedEnergy(
    MultiBlockLattice3D<T, Descriptor> &lattice, ForceFunction f, Box3D domain);

template <typename T, template <typename U> class Descriptor, class ForceFunction>
T computeAverageForcedEnergy(MultiBlockLattice3D<T, Descriptor> &lattice, ForceFunction f);

template <typename T, template <typename U> class Descriptor, class BoolMask>
plint count(MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, BoolMask boolMask);

template <typename T, template <typename U> class Descriptor, class BoolMask>
plint count(MultiBlockLattice3D<T, Descriptor> &lattice, BoolMask boolMask);

template <typename T>
std::vector<T> scalarSingleProbes(
    MultiScalarField3D<T> &scalarField, Box3D domain, std::vector<Array<T, 3> > const &positions);

template <typename T, template <typename U> class Descriptor>
std::vector<T> densitySingleProbes(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain,
    std::vector<Array<T, 3> > const &positions);

template <typename T, template <typename U> class Descriptor>
std::vector<T> densitySingleProbes(
    MultiBlockLattice3D<T, Descriptor> &lattice, std::vector<Array<T, 3> > const &positions);

template <typename T, template <typename U> class Descriptor>
std::vector<Array<T, 3> > velocitySingleProbes(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain,
    std::vector<Array<T, 3> > const &positions);

template <typename T, template <typename U> class Descriptor>
std::vector<Array<T, 3> > velocitySingleProbes(
    MultiBlockLattice3D<T, Descriptor> &lattice, std::vector<Array<T, 3> > const &positions);

template <typename T, template <typename U> class Descriptor>
std::vector<Array<T, 3> > vorticitySingleProbes(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain,
    std::vector<Array<T, 3> > const &positions);

template <typename T, template <typename U> class Descriptor>
std::vector<Array<T, 3> > vorticitySingleProbes(
    MultiBlockLattice3D<T, Descriptor> &lattice, std::vector<Array<T, 3> > const &positions);

/* *************** Extract Sub-Lattice ******************************* */

template <typename T, template <typename U> class Descriptor>
void extractSubDomain(
    MultiBlockLattice3D<T, Descriptor> &lattice,
    MultiBlockLattice3D<T, Descriptor> &extractedLattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > extractSubDomain(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

/* *************** Density ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeDensity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &density, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeDensity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeDensity(MultiBlockLattice3D<T, Descriptor> &lattice);

/* *************** RhoBar ******************************************** */

template <typename T, template <typename U> class Descriptor>
void computeRhoBar(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &rhoBar, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeRhoBar(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeRhoBar(MultiBlockLattice3D<T, Descriptor> &lattice);

/* *************** RhoBarJ ******************************************** */

template <typename T, template <typename U> class Descriptor>
void computeRhoBarJ(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &rhoBar,
    MultiTensorField3D<T, 3> &j, Box3D domain);

template <typename T, template <typename U> class Descriptor>
void maskedComputeRhoBarJ(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &rhoBar,
    MultiTensorField3D<T, 3> &j, MultiScalarField3D<int> &mask, int whichFlag, Box3D domain);

/* *************** RhoBarJPiNeq ******************************************** */

template <typename T, template <typename U> class Descriptor>
void computeRhoBarJPiNeq(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &rhoBar,
    MultiTensorField3D<T, 3> &j, TensorField3D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq,
    Box3D domain);

/* *************** Packed RhoBar J *********************************** */

template <typename T, template <typename U> class Descriptor>
void computePackedRhoBarJ(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiNTensorField3D<T> &rhoBarJ, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiNTensorField3D<T> > computePackedRhoBarJ(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiNTensorField3D<T> > computePackedRhoBarJ(
    MultiBlockLattice3D<T, Descriptor> &lattice);

template <typename T>
void computeDensityFromRhoBarJ(
    MultiNTensorField3D<T> &rhoBarJ, MultiScalarField3D<T> &density, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeDensityFromRhoBarJ(
    MultiNTensorField3D<T> &rhoBarJ, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeDensityFromRhoBarJ(MultiNTensorField3D<T> &rhoBarJ);

template <typename T>
void computeVelocityFromRhoBarJ(
    MultiNTensorField3D<T> &rhoBarJ, MultiTensorField3D<T, 3> &velocity, Box3D domain,
    bool velIsJ = false);

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeVelocityFromRhoBarJ(
    MultiNTensorField3D<T> &rhoBarJ, Box3D domain, bool velIsJ = false);

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeVelocityFromRhoBarJ(
    MultiNTensorField3D<T> &rhoBarJ, bool velIsJ = false);

template <typename T>
void computeVelocityFromRhoBarAndJ(
    MultiScalarField3D<T> &rhoBar, MultiTensorField3D<T, 3> &j, MultiTensorField3D<T, 3> &velocity,
    Box3D domain, bool velIsJ = false);

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeVelocityFromRhoBarAndJ(
    MultiScalarField3D<T> &rhoBar, MultiTensorField3D<T, 3> &j, Box3D domain, bool velIsJ = false);

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeVelocityFromRhoBarAndJ(
    MultiScalarField3D<T> &rhoBar, MultiTensorField3D<T, 3> &j, bool velIsJ = false);

/* *************** Kinetic Energy ************************************ */

template <typename T, template <typename U> class Descriptor>
void computeKineticEnergy(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &energy, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeKineticEnergy(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeKineticEnergy(
    MultiBlockLattice3D<T, Descriptor> &lattice);

/* *************** Velocity Norm ************************************* */

template <typename T, template <typename U> class Descriptor>
void computeVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &velocityNorm, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
void computeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<T> &velocityNorm, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force);

template <typename T, template <typename U> class Descriptor>
void computeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force,
    MultiScalarField3D<T> &velocityNorm, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force);

template <typename T, template <typename U> class Descriptor, class ForceFunction>
void computeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, ForceFunction f,
    MultiScalarField3D<T> &velocityNorm, Box3D domain);

template <typename T, template <typename U> class Descriptor, class ForceFunction>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, ForceFunction f, Box3D domain);

template <typename T, template <typename U> class Descriptor, class ForceFunction>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, ForceFunction f);

/* *************** Velocity Component ******************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &velocityComponent,
    Box3D domain, plint iComponent);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, plint iComponent);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, plint iComponent);

template <typename T, template <typename U> class Descriptor>
void computeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<T> &velocityComponent, Box3D domain, plint iComponent);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    Box3D domain, plint iComponent);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    plint iComponent);

template <typename T, template <typename U> class Descriptor>
void computeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force,
    MultiScalarField3D<T> &velocityComponent, Box3D domain, plint iComponent);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force, Box3D domain,
    plint iComponent);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force,
    plint iComponent);

template <typename T, template <typename U> class Descriptor, class ForceFunction>
void computeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, ForceFunction f,
    MultiScalarField3D<T> &velocityComponent, Box3D domain, plint iComponent);

template <typename T, template <typename U> class Descriptor, class ForceFunction>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, ForceFunction f, Box3D domain, plint iComponent);

template <typename T, template <typename U> class Descriptor, class ForceFunction>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, ForceFunction f, plint iComponent);

/* *************** Velocity ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &velocity,
    Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
void computeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiTensorField3D<T, Descriptor<T>::d> &velocity, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force);

template <typename T, template <typename U> class Descriptor>
void computeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force,
    MultiTensorField3D<T, Descriptor<T>::d> &velocity, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force);

template <typename T, template <typename U> class Descriptor, class ForceFunction>
void computeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, ForceFunction f,
    MultiTensorField3D<T, Descriptor<T>::d> &velocity, Box3D domain);

template <typename T, template <typename U> class Descriptor, class ForceFunction>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, ForceFunction f, Box3D domain);

template <typename T, template <typename U> class Descriptor, class ForceFunction>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, ForceFunction f);

/* *************** Temperature ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeTemperature(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &temperature, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeTemperature(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeTemperature(
    MultiBlockLattice3D<T, Descriptor> &lattice);

/* *************** PiNeq ********************************* */

template <typename T, template <typename U> class Descriptor>
void computePiNeq(
    MultiBlockLattice3D<T, Descriptor> &lattice,
    MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computePiNeq(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computePiNeq(
    MultiBlockLattice3D<T, Descriptor> &lattice);

/* *************** ShearStress ********************************* */

template <typename T, template <typename U> class Descriptor>
void computeShearStress(
    MultiBlockLattice3D<T, Descriptor> &lattice,
    MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computeShearStress(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computeShearStress(
    MultiBlockLattice3D<T, Descriptor> &lattice);

/* *************** Stress ********************************* */

template <typename T, template <typename U> class Descriptor>
void computeStress(
    MultiBlockLattice3D<T, Descriptor> &lattice,
    MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> &stress, T rho0, bool isCompressible,
    Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computeStress(
    MultiBlockLattice3D<T, Descriptor> &lattice, T rho0, bool isCompressible, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computeStress(
    MultiBlockLattice3D<T, Descriptor> &lattice, T rho0, bool isCompressible);

/* *************** Strain Rate from Stress *************************** */

template <typename T, template <typename U> class Descriptor>
void computeStrainRateFromStress(
    MultiBlockLattice3D<T, Descriptor> &lattice,
    MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> &S, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> >
    computeStrainRateFromStress(MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> >
    computeStrainRateFromStress(MultiBlockLattice3D<T, Descriptor> &lattice);

/* *************** Shear Rate ******************************************** */

template <typename T, template <typename U> class Descriptor>
void computeShearRate(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &shearRate, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeShearRate(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeShearRate(
    MultiBlockLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
void computeShearRate_N(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiNTensorField3D<T> &shearRate, Box3D domain);

template <typename T, template <typename U> class Descriptor>
MultiNTensorField3D<T> *computeShearRate_N(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

/* *************** Population **************************************** */

template <typename T, template <typename U> class Descriptor>
void computePopulation(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &population, Box3D domain,
    plint iPop);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computePopulation(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, plint iPop);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computePopulation(
    MultiBlockLattice3D<T, Descriptor> &lattice, plint iPop);

template <typename T, template <typename U> class Descriptor>
void computeEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &equilibrium, Box3D domain,
    plint iPop);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, plint iPop);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice, plint iPop);

template <typename T, template <typename U> class Descriptor>
void computeEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice,
    MultiTensorField3D<T, Descriptor<T>::q> &equilibrium, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::q> > computeEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::q> > computeEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
void computeNonEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice,
    MultiTensorField3D<T, Descriptor<T>::q> &nonEquilibrium, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::q> > computeNonEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::q> > computeNonEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
void computeAllPopulations(
    MultiBlockLattice3D<T, Descriptor> &lattice,
    MultiTensorField3D<T, Descriptor<T>::q> &populations);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::q> > computeAllPopulations(
    MultiBlockLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::q> > computeAllPopulations(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
void computeAllPopulationsFromTensorField(
    MultiBlockLattice3D<T, Descriptor> &lattice,
    MultiTensorField3D<T, Descriptor<T>::q> &populations);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::q> > computeAllPopulationsFromTensorField(
    MultiBlockLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::q> > computeAllPopulationsFromTensorField(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
void copyPopulations(
    MultiBlockLattice3D<T, Descriptor> &latticeFrom, MultiBlockLattice3D<T, Descriptor> &latticeTo,
    Box3D domain);

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void copyConvertPopulations(
    MultiBlockLattice3D<T1, Descriptor1> &latticeFrom,
    MultiBlockLattice3D<T2, Descriptor2> &latticeTo, Box3D domain);

template <typename T, template <typename U> class Descriptor>
void copyAll(
    MultiBlockLattice3D<T, Descriptor> &latticeFrom, MultiBlockLattice3D<T, Descriptor> &latticeTo,
    Box3D domain);

template <typename T, template <typename U> class Descriptor>
void copyRegenerate(
    MultiBlockLattice3D<T, Descriptor> &latticeFrom, MultiBlockLattice3D<T, Descriptor> &latticeTo,
    Box3D domain);

/* *************** Omega ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeOmega(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &omega, Box3D domain);

template <typename T, template <typename U> class Descriptor>
void computeOmega(MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &omega);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeOmega(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeOmega(MultiBlockLattice3D<T, Descriptor> &lattice);

/* *************** KinematicViscosity ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeKinematicViscosity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &omega, Box3D domain);

template <typename T, template <typename U> class Descriptor>
void computeKinematicViscosity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &omega);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeKinematicViscosity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeKinematicViscosity(
    MultiBlockLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
void computeKinematicViscosity_N(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiNTensorField3D<T> &nu, Box3D domain);

template <typename T, template <typename U> class Descriptor>
MultiNTensorField3D<T> *computeKinematicViscosity_N(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
MultiNTensorField3D<T> *computeKinematicViscosity_N(MultiBlockLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
void computeKinematicEddyViscosity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &omega, Box3D domain);

template <typename T, template <typename U> class Descriptor>
void computeKinematicEddyViscosity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &omega);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeKinematicEddyViscosity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeKinematicEddyViscosity(
    MultiBlockLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
void computeKinematicEddyViscosity_N(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiNTensorField3D<T> &nu, Box3D domain);

template <typename T, template <typename U> class Descriptor>
MultiNTensorField3D<T> *computeKinematicEddyViscosity_N(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
MultiNTensorField3D<T> *computeKinematicEddyViscosity_N(
    MultiBlockLattice3D<T, Descriptor> &lattice);

/* *************** ExternalForce ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeExternalForce(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeExternalForce(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeExternalForce(
    MultiBlockLattice3D<T, Descriptor> &lattice);

/* *************** ExternalScalar ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeExternalScalar(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &scalar, int whichScalar,
    Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeExternalScalar(
    MultiBlockLattice3D<T, Descriptor> &lattice, int whichScalar, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeExternalScalar(
    MultiBlockLattice3D<T, Descriptor> &lattice, int whichScalar);

/* *************** ExternalVector ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeExternalVector(
    MultiBlockLattice3D<T, Descriptor> &lattice,
    MultiTensorField3D<T, Descriptor<T>::d> &tensorField, int vectorBeginsAt, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeExternalVector(
    MultiBlockLattice3D<T, Descriptor> &lattice, int vectorBeginsAt, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeExternalVector(
    MultiBlockLattice3D<T, Descriptor> &lattice, int vectorBeginsAt);

/* *************** DynamicParameter ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeDynamicParameter(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &scalar,
    plint whichParameter, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeDynamicParameter(
    MultiBlockLattice3D<T, Descriptor> &lattice, plint whichParameter, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeDynamicParameter(
    MultiBlockLattice3D<T, Descriptor> &lattice, plint whichParameter);

/* *************** DynamicViscosity ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeDynamicViscosity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &scalar, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeDynamicViscosity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeDynamicViscosity(
    MultiBlockLattice3D<T, Descriptor> &lattice);

/* ******************************************************************* */
/* *************** PART V. Multi-block wrappers: Scalar-Field ******** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template <typename T>
T computeSum(MultiScalarField3D<T> &scalarField, Box3D domain);

template <typename T>
T computeSum(MultiScalarField3D<T> &scalarField);

template <typename T>
T computeSum(
    MultiScalarField3D<T> &scalarField, MultiScalarField3D<int> &mask, int flag, Box3D domain);

template <typename T>
T computeSum(MultiScalarField3D<T> &scalarField, MultiScalarField3D<int> &mask, int flag);

template <typename T>
plint computeIntSum(MultiScalarField3D<T> &scalarField, Box3D domain);

template <typename T>
T computeAverage(MultiScalarField3D<T> &scalarField, Box3D domain);

template <typename T>
T computeAverage(MultiScalarField3D<T> &scalarField);

template <typename T>
T computeAverage(
    MultiScalarField3D<T> &scalarField, MultiScalarField3D<int> &mask, int flag, Box3D domain);

template <typename T>
T computeAverage(MultiScalarField3D<T> &scalarField, MultiScalarField3D<int> &mask, int flag);

template <typename T>
T computeMin(MultiScalarField3D<T> &scalarField, Box3D domain);

template <typename T>
T computeMin(MultiScalarField3D<T> &scalarField);

template <typename T>
T computeMin(
    MultiScalarField3D<T> &scalarField, MultiScalarField3D<int> &mask, int flag, Box3D domain);

template <typename T>
T computeMin(MultiScalarField3D<T> &scalarField, MultiScalarField3D<int> &mask, int flag);

template <typename T>
T computeMax(MultiScalarField3D<T> &scalarField, Box3D domain);

template <typename T>
T computeMax(MultiScalarField3D<T> &scalarField);

template <typename T>
T computeMax(
    MultiScalarField3D<T> &scalarField, MultiScalarField3D<int> &mask, int flag, Box3D domain);

template <typename T>
T computeMax(MultiScalarField3D<T> &scalarField, MultiScalarField3D<int> &mask, int flag);

template <typename T>
T computeBoundedAverage(MultiScalarField3D<T> &scalarField, Box3D domain);

template <typename T>
T computeBoundedAverage(MultiScalarField3D<T> &scalarField);

template <typename T, class BoolMask>
plint count(MultiScalarField3D<T> &field, Box3D domain, BoolMask boolMask);

template <typename T, class BoolMask>
plint count(MultiScalarField3D<T> &field, BoolMask boolMask);

/* *************** Generic operations *************** */

template <typename T, class Function>
void apply(Function f, MultiScalarField3D<T> &field, Box3D domain)
{
    applyProcessingFunctional(new ApplyScalarFunctional3D<T, Function>(f), domain, field);
}

template <typename T, class Function>
void apply(Function f, MultiScalarField3D<T> &field)
{
    apply(f, field, field.getBoundingBox());
}

template <typename T, class Function>
void evaluate(Function f, MultiScalarField3D<T> &field, MultiScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(
        new EvaluateScalarFunctional3D<T, Function>(f), domain, field, result);
}

template <typename T, class Function>
std::unique_ptr<MultiScalarField3D<T> > evaluate(
    Function f, MultiScalarField3D<T> &field, Box3D domain)
{
    MultiScalarField3D<T> *result = new MultiScalarField3D<T>(field, domain);
    evaluate(f, field, *result, domain);
    return std::unique_ptr<MultiScalarField3D<T> >(result);
}

template <typename T, class Function>
std::unique_ptr<MultiScalarField3D<T> > evaluate(Function f, MultiScalarField3D<T> &field)
{
    return evaluate(f, field, field.getBoundingBox());
}

/* *************** Extract Sub-ScalarField *************************** */

template <typename T>
void extractSubDomain(
    MultiScalarField3D<T> &field, MultiScalarField3D<T> &extractedField, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > extractSubDomain(
    MultiScalarField3D<T> &field, Box3D domain);

/* *************** Copy and Convert ScalarField *************************** */

template <typename T1, typename T2>
void copy(MultiScalarField3D<T1> &field, MultiScalarField3D<T2> &convertedField, Box3D domain);

template <typename T1, typename T2>
std::unique_ptr<MultiScalarField3D<T2> > copyConvert(MultiScalarField3D<T1> &field, Box3D domain);

template <typename T1, typename T2>
std::unique_ptr<MultiScalarField3D<T2> > copyConvert(MultiScalarField3D<T1> &field);

/* *************** MultiScalarField - Scalar operations *************** */

template <typename T>
void lessThan(
    MultiScalarField3D<T> &field, T scalar, MultiScalarField3D<int> &result, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<int> > lessThan(
    MultiScalarField3D<T> &field, T scalar, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<int> > lessThan(MultiScalarField3D<T> &field, T scalar);

template <typename T>
void greaterThan(
    MultiScalarField3D<T> &field, T scalar, MultiScalarField3D<int> &result, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<int> > greaterThan(
    MultiScalarField3D<T> &field, T scalar, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<int> > greaterThan(MultiScalarField3D<T> &field, T scalar);

template <typename T>
void add(T scalar, MultiScalarField3D<T> &field, MultiScalarField3D<T> &result, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > add(T scalar, MultiScalarField3D<T> &field, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > add(T scalar, MultiScalarField3D<T> &field);

template <typename T>
void add(MultiScalarField3D<T> &field, T scalar, MultiScalarField3D<T> &result, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > add(MultiScalarField3D<T> &field, T scalar, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > add(MultiScalarField3D<T> &field, T scalar);

template <typename T>
void subtract(T scalar, MultiScalarField3D<T> &field, MultiScalarField3D<T> &result, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > subtract(
    T scalar, MultiScalarField3D<T> &field, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > subtract(T scalar, MultiScalarField3D<T> &field);

template <typename T>
void subtract(MultiScalarField3D<T> &field, T scalar, MultiScalarField3D<T> &result, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > subtract(
    MultiScalarField3D<T> &field, T scalar, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > subtract(MultiScalarField3D<T> &field, T scalar);

template <typename T>
void multiply(T scalar, MultiScalarField3D<T> &field, MultiScalarField3D<T> &result, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > multiply(
    T scalar, MultiScalarField3D<T> &field, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > multiply(T scalar, MultiScalarField3D<T> &field);

template <typename T>
void multiply(MultiScalarField3D<T> &field, T scalar, MultiScalarField3D<T> &result, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > multiply(
    MultiScalarField3D<T> &field, T scalar, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > multiply(MultiScalarField3D<T> &field, T scalar);

template <typename T>
void divide(T scalar, MultiScalarField3D<T> &field, MultiScalarField3D<T> &result, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > divide(
    T scalar, MultiScalarField3D<T> &field, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > divide(T scalar, MultiScalarField3D<T> &field);

template <typename T>
void divide(MultiScalarField3D<T> &field, T scalar, MultiScalarField3D<T> &result, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > divide(
    MultiScalarField3D<T> &field, T scalar, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > divide(MultiScalarField3D<T> &field, T scalar);

/* *************** MultiScalarField operations *************** */

template <typename T>
void computeSqrt(MultiScalarField3D<T> &field, MultiScalarField3D<T> &result, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeSqrt(MultiScalarField3D<T> &field, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeSqrt(MultiScalarField3D<T> &field);

template <typename T>
void computeAbsoluteValue(
    MultiScalarField3D<T> &field, MultiScalarField3D<T> &result, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeAbsoluteValue(
    MultiScalarField3D<T> &field, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeAbsoluteValue(MultiScalarField3D<T> &field);

/* *************** MultiScalarField - Scalar inplace operations *************** */

template <typename T>
void addInPlace(MultiScalarField3D<T> &field, T scalar, Box3D domain);

template <typename T>
void addInPlace(MultiScalarField3D<T> &field, T scalar);

template <typename T>
void subtractInPlace(MultiScalarField3D<T> &field, T scalar, Box3D domain);

template <typename T>
void subtractInPlace(MultiScalarField3D<T> &field, T scalar);

template <typename T>
void multiplyInPlace(MultiScalarField3D<T> &field, T scalar, Box3D domain);

template <typename T>
void multiplyInPlace(MultiScalarField3D<T> &field, T scalar);

template <typename T>
void divideInPlace(MultiScalarField3D<T> &field, T scalar, Box3D domain);

template <typename T>
void divideInPlace(MultiScalarField3D<T> &field, T scalar);

/* *************** MultiScalarField - MultiScalarField operations *************** */

template <typename T>
void lessThan(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, MultiScalarField3D<int> &result,
    Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<int> > lessThan(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<int> > lessThan(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B);

template <typename T>
void greaterThan(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, MultiScalarField3D<int> &result,
    Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<int> > greaterThan(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<int> > greaterThan(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B);

template <typename T>
void add(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, MultiScalarField3D<T> &result,
    Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > add(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > add(MultiScalarField3D<T> &A, MultiScalarField3D<T> &B);

template <typename T>
void subtract(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, MultiScalarField3D<T> &result,
    Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > subtract(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > subtract(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B);

template <typename T>
void multiply(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, MultiScalarField3D<T> &result,
    Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > multiply(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > multiply(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B);

template <typename T>
void divide(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, MultiScalarField3D<T> &result,
    Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > divide(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > divide(MultiScalarField3D<T> &A, MultiScalarField3D<T> &B);

/* *************** MultiScalarField - MultiScalarField inplace operations *************** */

template <typename T>
void addInPlace(MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, Box3D domain);

template <typename T>
void addInPlace(MultiScalarField3D<T> &A, MultiScalarField3D<T> &B);

template <typename T>
void subtractInPlace(MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, Box3D domain);

template <typename T>
void subtractInPlace(MultiScalarField3D<T> &A, MultiScalarField3D<T> &B);

template <typename T>
void multiplyInPlace(MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, Box3D domain);

template <typename T>
void multiplyInPlace(MultiScalarField3D<T> &A, MultiScalarField3D<T> &B);

template <typename T>
void divideInPlace(MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, Box3D domain);

template <typename T>
void divideInPlace(MultiScalarField3D<T> &A, MultiScalarField3D<T> &B);

template <typename T>
void uniformlyBoundScalarField(MultiScalarField3D<T> &data, Box3D domain, T bound);

template <typename T>
void uniformlyBoundScalarField(MultiScalarField3D<T> &data, T bound);

template <typename T>
void boundScalarField(MultiScalarField3D<T> &data, Box3D domain, T lowerBound, T upperBound);

template <typename T>
void boundScalarField(MultiScalarField3D<T> &data, T lowerBound, T upperBound);

/* *************** LBMsmoothen3D ******************************************* */
// This is a "bulk" version. It has no special boundary treatment.

template <typename T, template <typename U> class Descriptor>
void lbmSmoothen(MultiScalarField3D<T> &data, MultiScalarField3D<T> &result, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > lbmSmoothen(MultiScalarField3D<T> &data, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > lbmSmoothen(MultiScalarField3D<T> &data);

/* *************** LBMsmoothenInPlace3D ******************************************* */
// This is a "bulk" version. It has no special boundary treatment.

template <typename T, template <typename U> class Descriptor>
void lbmSmoothenInPlace(MultiScalarField3D<T> &data, Box3D domain);

template <typename T, template <typename U> class Descriptor>
void lbmSmoothenInPlace(MultiScalarField3D<T> &data);

/* *************** Smoothen3D ******************************************* */
// This is a "bulk" version. It has no special boundary treatment.

template <typename T>
void smoothen(MultiScalarField3D<T> &data, MultiScalarField3D<T> &result, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > smoothen(MultiScalarField3D<T> &data, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > smoothen(MultiScalarField3D<T> &data);

/* *************** MollifyScalar3D ******************************************* */
// This is a "bulk" version. It has no special boundary treatment.

template <typename T>
void mollifyScalar(
    T l, plint d, Box3D globalDomain, int exclusionFlag, MultiScalarField3D<T> &data,
    MultiScalarField3D<int> &flag, MultiScalarField3D<T> &result, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > mollifyScalar(
    T l, plint d, Box3D globalDomain, int exclusionFlag, MultiScalarField3D<T> &data,
    MultiScalarField3D<int> &flag, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > mollifyScalar(
    T l, plint d, Box3D globalDomain, int exclusionFlag, MultiScalarField3D<T> &data,
    MultiScalarField3D<int> &flag);

/* *************** LBMcomputeGradient3D ******************************************* */
// This is a "bulk" version. It has no special boundary treatment.

template <typename T, template <typename U> class Descriptor>
void lbmComputeGradient(
    MultiScalarField3D<T> &scalarField, MultiTensorField3D<T, 3> &gradient, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, 3> > lbmComputeGradient(
    MultiScalarField3D<T> &scalarField, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, 3> > lbmComputeGradient(MultiScalarField3D<T> &scalarField);

/* ******************************************************************* */
/* *************** PART VI. Multi-block wrappers: Tensor-field ******* */
/* ******************************************************************* */

/* *************** Reductive functions ************** */

template <typename T, int nDim, class BoolMask>
plint count(MultiTensorField3D<T, nDim> &field, Box3D domain, BoolMask boolMask);

template <typename T, int nDim, class BoolMask>
plint count(MultiTensorField3D<T, nDim> &field, BoolMask boolMask);

template <typename T, int nDim>
Array<T, nDim> computeSum(MultiTensorField3D<T, nDim> &tensorField, Box3D domain);

template <typename T, int nDim>
Array<T, nDim> computeSum(MultiTensorField3D<T, nDim> &tensorField);

template <typename T, int nDim>
Array<T, nDim> computeSum(
    MultiTensorField3D<T, nDim> &tensorField, MultiScalarField3D<int> &mask, int flag,
    Box3D domain);

template <typename T, int nDim>
Array<T, nDim> computeSum(
    MultiTensorField3D<T, nDim> &tensorField, MultiScalarField3D<int> &mask, int flag);

template <typename T, int nDim>
Array<T, nDim> computeAverage(MultiTensorField3D<T, nDim> &tensorField, Box3D domain);

template <typename T, int nDim>
Array<T, nDim> computeAverage(MultiTensorField3D<T, nDim> &tensorField);

template <typename T, int nDim>
Array<T, nDim> computeAverage(
    MultiTensorField3D<T, nDim> &tensorField, MultiScalarField3D<int> &mask, int flag,
    Box3D domain);

template <typename T, int nDim>
Array<T, nDim> computeAverage(
    MultiTensorField3D<T, nDim> &tensorField, MultiScalarField3D<int> &mask, int flag);

/* *************** Copy-convert a tensor-field *************** */

template <typename T1, typename T2, int nDim>
void copy(
    MultiTensorField3D<T1, nDim> &field, MultiTensorField3D<T2, nDim> &convertedField,
    Box3D domain);

template <typename T1, typename T2, int nDim>
std::unique_ptr<MultiTensorField3D<T2, nDim> > copyConvert(
    MultiTensorField3D<T1, nDim> &field, Box3D domain);

template <typename T1, typename T2, int nDim>
std::unique_ptr<MultiTensorField3D<T2, nDim> > copyConvert(MultiTensorField3D<T1, nDim> &field);

/* *************** Extract Sub-TensorField *************************** */

template <typename T, int nDim>
void extractSubDomain(
    MultiTensorField3D<T, nDim> &field, MultiTensorField3D<T, nDim> &extractedField, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > extractSubDomain(
    MultiTensorField3D<T, nDim> &field, Box3D domain);

/* *************** Component (scalar-field) out of a tensor-field ****** */

template <typename T, int nDim>
void extractComponent(
    MultiTensorField3D<T, nDim> &tensorField, MultiScalarField3D<T> &component, Box3D domain,
    int iComponent);

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T> > extractComponent(
    MultiTensorField3D<T, nDim> &tensorField, Box3D domain, int iComponent);

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T> > extractComponent(
    MultiTensorField3D<T, nDim> &tensorField, int iComponent);

/* *************** Vector-norm of each cell in the field *************** */

template <typename T, int nDim>
void computeNorm(
    MultiTensorField3D<T, nDim> &tensorField, MultiScalarField3D<T> &norm, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T> > computeNorm(
    MultiTensorField3D<T, nDim> &tensorField, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T> > computeNorm(MultiTensorField3D<T, nDim> &tensorField);

/* *************** Squared vector-norm of each cell in the field ******** */

template <typename T, int nDim>
void computeNormSqr(
    MultiTensorField3D<T, nDim> &tensorField, MultiScalarField3D<T> &normSqr, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T> > computeNormSqr(
    MultiTensorField3D<T, nDim> &tensorField, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T> > computeNormSqr(MultiTensorField3D<T, nDim> &tensorField);

/* *************** Max element of each array of each cell *************** */

template <typename T, int nDim>
void computeMaximumElement(
    MultiTensorField3D<T, nDim> &A, MultiScalarField3D<T> &result, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T> > computeMaximumElement(
    MultiTensorField3D<T, nDim> &A, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T> > computeMaximumElement(MultiTensorField3D<T, nDim> &A);

/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template <typename T>
void computeSymmetricTensorNorm(
    MultiTensorField3D<T, 6> &tensorField, MultiScalarField3D<T> &norm, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeSymmetricTensorNorm(
    MultiTensorField3D<T, 6> &tensorField, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeSymmetricTensorNorm(
    MultiTensorField3D<T, 6> &tensorField);

/* *************** Squared Tensor-norm of each symmetric tensor of a field **************** */

template <typename T>
void computeSymmetricTensorNormSqr(
    MultiTensorField3D<T, 6> &tensorField, MultiScalarField3D<T> &normSqr, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeSymmetricTensorNormSqr(
    MultiTensorField3D<T, 6> &tensorField, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeSymmetricTensorNormSqr(
    MultiTensorField3D<T, 6> &tensorField);

/* *************** Trace of each symmetric tensor of a field ************ */

template <typename T>
void computeSymmetricTensorTrace(
    MultiTensorField3D<T, 6> &tensorField, MultiScalarField3D<T> &trace, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeSymmetricTensorTrace(
    MultiTensorField3D<T, 6> &tensorField, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeSymmetricTensorTrace(
    MultiTensorField3D<T, 6> &tensorField);

/* *************** Gradient from scalar field *********************** */

template <typename T>
void computeGradient(MultiScalarField3D<T> &phi, MultiTensorField3D<T, 3> &gradient, Box3D domain);

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeGradient(
    MultiScalarField3D<T> &phi, Box3D domain);

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeGradient(MultiScalarField3D<T> &phi);

/* *************** Gradient, witout boundary treatment, from scalar field **************** */

template <typename T>
void computeBulkGradient(
    MultiScalarField3D<T> &phi, MultiTensorField3D<T, 3> &gradient, Box3D domain);

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeBulkGradient(
    MultiScalarField3D<T> &phi, Box3D domain);

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeBulkGradient(MultiScalarField3D<T> &phi);

/* *************** Vorticity from Velocity field *********************** */

template <typename T>
void computeVorticity(
    MultiTensorField3D<T, 3> &velocity, MultiTensorField3D<T, 3> &vorticity, Box3D domain);

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeVorticity(
    MultiTensorField3D<T, 3> &velocity, Box3D domain);

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeVorticity(MultiTensorField3D<T, 3> &velocity);

/* *************** Vorticity, witout boundary treatment, from Velocity field ************** */

template <typename T>
void computeBulkVorticity(
    MultiTensorField3D<T, 3> &velocity, MultiTensorField3D<T, 3> &vorticity, Box3D domain);

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeBulkVorticity(
    MultiTensorField3D<T, 3> &velocity, Box3D domain);

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeBulkVorticity(MultiTensorField3D<T, 3> &velocity);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeBulkVorticity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeBulkVorticity(
    MultiBlockLattice3D<T, Descriptor> &lattice);

/* *************** 4th order Vorticity, witout boundary treatment, from Velocity field
 * ************** */

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > computeBulkVorticityOrderFour(
    MultiTensorField3D<T, nDim> &vel);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeBulkVorticityOrderFour(
    MultiBlockLattice3D<T, Descriptor> &lattice);

/* *************** 6th order Vorticity, witout boundary treatment, from Velocity field
 * ************** */

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > computeBulkVorticityOrderSix(
    MultiTensorField3D<T, nDim> &vel);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeBulkVorticityOrderSix(
    MultiBlockLattice3D<T, Descriptor> &lattice);

/* *************** Helicity from Velocity field *********************** */

template <typename T>
void computeHelicity(
    MultiTensorField3D<T, 3> &velocity, MultiScalarField3D<T> &helicity, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeHelicity(
    MultiTensorField3D<T, 3> &velocity, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeHelicity(MultiTensorField3D<T, 3> &velocity);

/* *************** Helicity, witout boundary treatment, from Velocity field ************** */

template <typename T>
void computeBulkHelicity(
    MultiTensorField3D<T, 3> &velocity, MultiScalarField3D<T> &helicity, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeBulkHelicity(
    MultiTensorField3D<T, 3> &velocity, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeBulkHelicity(MultiTensorField3D<T, 3> &velocity);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeBulkHelicity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeBulkHelicity(
    MultiBlockLattice3D<T, Descriptor> &lattice);

/* *************** Divergence, witout boundary treatment, from Velocity field *************** */

template <typename T>
void computeBulkDivergence(
    MultiTensorField3D<T, 3> &velocity, MultiScalarField3D<T> &divergence, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeBulkDivergence(
    MultiTensorField3D<T, 3> &velocity, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeBulkDivergence(MultiTensorField3D<T, 3> &velocity);

/* *************** Strain Rate from Velocity field ********************* */

template <typename T>
void computeStrainRate(
    MultiTensorField3D<T, 3> &velocity, MultiTensorField3D<T, 6> &S, Box3D domain);

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 6> > computeStrainRate(
    MultiTensorField3D<T, 3> &velocity, Box3D domain);

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 6> > computeStrainRate(MultiTensorField3D<T, 3> &velocity);

/* *************** Strain rate, witout boundary treatment, from Velocity field *************** */

template <typename T>
void computeBulkStrainRate(
    MultiTensorField3D<T, 3> &velocity, MultiTensorField3D<T, 6> &S, Box3D domain);

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 6> > computeBulkStrainRate(
    MultiTensorField3D<T, 3> &velocity, Box3D domain);

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 6> > computeBulkStrainRate(
    MultiTensorField3D<T, 3> &velocity);

/* *************** Q-criterion from vorticity and strain rate fields ******************** */
// This is useful for vortex identification. Vortex regions are where Q > 0.

template <typename T>
void computeQcriterion(
    MultiTensorField3D<T, 3> &vorticity, MultiTensorField3D<T, 6> &S,
    MultiScalarField3D<T> &qCriterion, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeQcriterion(
    MultiTensorField3D<T, 3> &vorticity, MultiTensorField3D<T, 6> &S, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeQcriterion(
    MultiTensorField3D<T, 3> &vorticity, MultiTensorField3D<T, 6> &S);

/* *************** Instantaneous reynolds stress computation from avg vel and vel
 * ******************** */
template <typename T>
void computeInstantaneousReynoldsStress(
    MultiTensorField3D<T, 3> &vel, MultiTensorField3D<T, 3> &avgVel, MultiTensorField3D<T, 6> &tau,
    Box3D domain);

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 6> > computeInstantaneousReynoldsStress(
    MultiTensorField3D<T, 3> &vel, MultiTensorField3D<T, 3> &avgVel, Box3D domain);

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 6> > computeInstantaneousReynoldsStress(
    MultiTensorField3D<T, 3> &vel, MultiTensorField3D<T, 3> &avgVel);

/* *************** lambda2-criterion from vorticity and strain rate fields ******************** */
// This is useful for vortex identification. Vortex regions are where lambda2 < 0.
#ifndef PLB_BGP
#ifdef PLB_USE_EIGEN
template <typename T>
void computeLambda2(
    MultiTensorField3D<T, 3> &vorticity, MultiTensorField3D<T, 6> &S,
    MultiScalarField3D<T> &lambda2, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeLambda2(
    MultiTensorField3D<T, 3> &vorticity, MultiTensorField3D<T, 6> &S, Box3D domain);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeLambda2(
    MultiTensorField3D<T, 3> &vorticity, MultiTensorField3D<T, 6> &S);
#endif
#endif

/* *************** MultiTensorField - MultiTensorField operations *************** */

template <typename T, int nDim>
void add(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B,
    MultiTensorField3D<T, nDim> &result, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > add(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > add(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B);

template <typename T, int nDim>
void subtract(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B,
    MultiTensorField3D<T, nDim> &result, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > subtract(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > subtract(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B);

template <typename T, int nDim>
void multiply(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B,
    MultiTensorField3D<T, nDim> &result, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > multiply(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > multiply(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B);

template <typename T, int nDim>
void multiply(
    T scalar, MultiTensorField3D<T, nDim> &field, MultiTensorField3D<T, nDim> &result,
    Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > multiply(
    T scalar, MultiTensorField3D<T, nDim> &field, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > multiply(
    T scalar, MultiTensorField3D<T, nDim> &field);

template <typename T, int nDim>
void multiply(
    MultiTensorField3D<T, nDim> &field, T scalar, MultiTensorField3D<T, nDim> &result,
    Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > multiply(
    MultiTensorField3D<T, nDim> &field, T scalar, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > multiply(
    MultiTensorField3D<T, nDim> &field, T scalar);

template <typename T, int nDim>
void divide(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B,
    MultiTensorField3D<T, nDim> &result, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > divide(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > divide(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B);

template <typename T, int nDim>
void normalize(
    MultiTensorField3D<T, nDim> &data, MultiTensorField3D<T, nDim> &result, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > normalize(
    MultiTensorField3D<T, nDim> &data, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > normalize(MultiTensorField3D<T, nDim> &data);

template <typename T, int nDim>
void symmetricTensorProduct(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, SymmetricTensorImpl<T, nDim>::n> &result,
    Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, SymmetricTensorImpl<T, nDim>::n> > symmetricTensorProduct(
    MultiTensorField3D<T, nDim> &A, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, SymmetricTensorImpl<T, nDim>::n> > symmetricTensorProduct(
    MultiTensorField3D<T, nDim> &A);

/* *************** MultiScalarField - MultiTensorField operations *************** */

template <typename T, int nDim>
void multiply(
    MultiScalarField3D<T> &A, MultiTensorField3D<T, nDim> &B, MultiTensorField3D<T, nDim> &result,
    Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > multiply(
    MultiScalarField3D<T> &A, MultiTensorField3D<T, nDim> &B, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > multiply(
    MultiScalarField3D<T> &A, MultiTensorField3D<T, nDim> &B);

template <typename T, int nDim>
void fullIndexContractionOfSymmetricTensors(
    MultiTensorField3D<T, SymmetricTensorImpl<T, nDim>::n> &A,
    MultiTensorField3D<T, SymmetricTensorImpl<T, nDim>::n> &B, MultiScalarField3D<T> &result,
    Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T> > fullIndexContractionOfSymmetricTensors(
    MultiTensorField3D<T, SymmetricTensorImpl<T, nDim>::n> &A,
    MultiTensorField3D<T, SymmetricTensorImpl<T, nDim>::n> &B, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T> > fullIndexContractionOfSymmetricTensors(
    MultiTensorField3D<T, SymmetricTensorImpl<T, nDim>::n> &A,
    MultiTensorField3D<T, SymmetricTensorImpl<T, nDim>::n> &B);

/* *************** MultiTensorField operations *************** */

template <typename T, int nDim>
void computeSqrt(
    MultiTensorField3D<T, nDim> &field, MultiTensorField3D<T, nDim> &result, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > computeSqrt(
    MultiTensorField3D<T, nDim> &field, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > computeSqrt(MultiTensorField3D<T, nDim> &field);

/* *************** MultiTensorField - MultiTensorField inplace operations *************** */

template <typename T, int nDim>
void addInPlace(MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B, Box3D domain);

template <typename T, int nDim>
void addInPlace(MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B);

template <typename T, int nDim>
void subtractInPlace(MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B, Box3D domain);

template <typename T, int nDim>
void subtractInPlace(MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B);

template <typename T, int nDim>
void multiplyInPlace(MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B, Box3D domain);

template <typename T, int nDim>
void multiplyInPlace(MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B);

template <typename T, int nDim>
void multiplyInPlace(MultiTensorField3D<T, nDim> &A, T alpha, Box3D domain);

template <typename T, int nDim>
void multiplyInPlace(MultiTensorField3D<T, nDim> &A, T alpha);

template <typename T, int nDim>
void divideInPlace(MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B, Box3D domain);

template <typename T, int nDim>
void divideInPlace(MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B);

template <typename T, int nDim>
void divideInPlace(MultiTensorField3D<T, nDim> &A, MultiScalarField3D<T> &B, Box3D domain);

template <typename T, int nDim>
void divideInPlace(MultiTensorField3D<T, nDim> &A, MultiScalarField3D<T> &B);

template <typename T, int nDim>
void divideInPlace(
    MultiTensorField3D<T, nDim> &A, MultiScalarField3D<T> &B, MultiScalarField3D<int> &mask,
    int flag, Box3D domain);

template <typename T, int nDim>
void divideInPlace(
    MultiTensorField3D<T, nDim> &A, MultiScalarField3D<T> &B, MultiScalarField3D<int> &mask,
    int flag);

template <typename T, int nDim>
void normalizeInPlace(MultiTensorField3D<T, nDim> &data, Box3D domain);

template <typename T, int nDim>
void normalizeInPlace(MultiTensorField3D<T, nDim> &data);

/* *************** LBMsmoothenTensor3D ******************************************* */
// This is a "bulk" version. It has no special boundary treatment.

template <typename T, int nDim, template <typename U> class Descriptor>
void lbmSmoothenTensor(
    MultiTensorField3D<T, nDim> &data, MultiTensorField3D<T, nDim> &result, Box3D domain);

template <typename T, int nDim, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, nDim> > lbmSmoothenTensor(
    MultiTensorField3D<T, nDim> &data, Box3D domain);

template <typename T, int nDim, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, nDim> > lbmSmoothenTensor(MultiTensorField3D<T, nDim> &data);

/* *************** LBMsmoothenTensorInPlace3D ******************************************* */
// This is a "bulk" version. It has no special boundary treatment.

template <typename T, int nDim, template <typename U> class Descriptor>
void lbmSmoothenTensorInPlace(MultiTensorField3D<T, nDim> &data, Box3D domain);

template <typename T, int nDim, template <typename U> class Descriptor>
void lbmSmoothenTensorInPlace(MultiTensorField3D<T, nDim> &data);

/* *************** SmoothenTensor3D ******************************************* */
// This is a "bulk" version. It has no special boundary treatment.

template <typename T, int nDim>
void smoothenTensor(
    MultiTensorField3D<T, nDim> &data, MultiTensorField3D<T, nDim> &result, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > smoothenTensor(
    MultiTensorField3D<T, nDim> &data, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > smoothenTensor(MultiTensorField3D<T, nDim> &data);

/* *************** MollifyTensor3D ******************************************* */
// This is a "bulk" version. It has no special boundary treatment.

template <typename T, int nDim>
void mollifyTensor(
    T l, plint d, Box3D globalDomain, int exclusionFlag, MultiTensorField3D<T, nDim> &data,
    MultiScalarField3D<int> &flag, MultiTensorField3D<T, nDim> &result, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > mollifyTensor(
    T l, plint d, Box3D globalDomain, int exclusionFlag, MultiTensorField3D<T, nDim> &data,
    MultiScalarField3D<int> &flag, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > mollifyTensor(
    T l, plint d, Box3D globalDomain, int exclusionFlag, MultiTensorField3D<T, nDim> &data,
    MultiScalarField3D<int> &flag);

/* *************** LBMcomputeDivergence3D ******************************************* */
// This is a "bulk" version. It has no special boundary treatment.

template <typename T, template <typename U> class Descriptor>
void lbmComputeDivergence(
    MultiScalarField3D<T> &divergence, MultiTensorField3D<T, 3> &vectorField, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > lbmComputeDivergence(
    MultiTensorField3D<T, 3> &vectorField, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > lbmComputeDivergence(MultiTensorField3D<T, 3> &vectorField);

/* ********************************************************************* */
/* *************** PART VII. Multi-block wrappers: NTensor-field ******* */
/* ********************************************************************* */

/* *************** Extract Sub-NTensorField **************************** */

template <typename T>
void extractSubDomain(
    MultiNTensorField3D<T> &field, MultiNTensorField3D<T> &extractedField, Box3D domain);

template <typename T>
std::unique_ptr<MultiNTensorField3D<T> > extractSubDomain(
    MultiNTensorField3D<T> &field, Box3D domain);

}  // namespace plb

#endif  // DATA_ANALYSIS_WRAPPER_3D_H
