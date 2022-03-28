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

#ifndef DATA_ANALYSIS_WRAPPER_2D_H
#define DATA_ANALYSIS_WRAPPER_2D_H

#include <memory>

#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include "core/globalDefs.h"
#include "dataProcessors/dataAnalysisFunctional2D.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiDataField2D.h"

namespace plb {

/* ******************************************************************* */
/* *************** PART I. Atomic-block wrappers: Block-Lattice ****** */
/* ******************************************************************* */

/* *************** Reductive operations ****************************** */

template <typename T, template <typename U> class Descriptor>
T computeAverageDensity(BlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
T computeAverageDensity(BlockLattice2D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(BlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(BlockLattice2D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(BlockLattice2D<T, Descriptor> &lattice, DotList2D dotList);

template <typename T, template <typename U> class Descriptor>
T computeAverageEnergy(BlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
T computeAverageEnergy(BlockLattice2D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
T computeAverageVelocityComponent(
    BlockLattice2D<T, Descriptor> &lattice, DotList2D dotList, plint iComponent);

template <typename T, template <typename U> class Descriptor, class BoolMask>
plint count(BlockLattice2D<T, Descriptor> &lattice, Box2D domain, BoolMask boolMask);

template <typename T, template <typename U> class Descriptor, class BoolMask>
plint count(BlockLattice2D<T, Descriptor> &lattice, BoolMask boolMask);

/* *************** Density ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeDensity(BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &density);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField2D<T> > computeDensity(BlockLattice2D<T, Descriptor> &lattice);

/* *************** RhoBar ******************************************** */

template <typename T, template <typename U> class Descriptor>
void computeRhoBar(BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &rhoBar);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField2D<T> > computeRhoBar(BlockLattice2D<T, Descriptor> &lattice);

/* *************** Kinetic Energy ************************************ */

template <typename T, template <typename U> class Descriptor>
void computeKineticEnergy(BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &energy);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField2D<T> > computeKineticEnergy(BlockLattice2D<T, Descriptor> &lattice);

/* *************** Velocity Norm ************************************* */

template <typename T, template <typename U> class Descriptor>
void computeVelocityNorm(BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &velocityNorm);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField2D<T> > computeVelocityNorm(BlockLattice2D<T, Descriptor> &lattice);

/* *************** Velocity Component ************************************ */

template <typename T, template <typename U> class Descriptor>
void computeVelocityComponent(
    BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &velocityComponent, plint iComponent);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField2D<T> > computeVelocityComponent(
    BlockLattice2D<T, Descriptor> &lattice, plint iComponent);

/* *************** Velocity ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocity(
    BlockLattice2D<T, Descriptor> &lattice, TensorField2D<T, Descriptor<T>::d> &velocity);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<TensorField2D<T, Descriptor<T>::d> > computeVelocity(
    BlockLattice2D<T, Descriptor> &lattice);

/* *************** Pi Neq ********************************* */

template <typename T, template <typename U> class Descriptor>
void computePiNeq(
    BlockLattice2D<T, Descriptor> &lattice,
    TensorField2D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<TensorField2D<T, SymmetricTensor<T, Descriptor>::n> > computePiNeq(
    BlockLattice2D<T, Descriptor> &lattice);

/* *************** Shear Stress ********************************* */

template <typename T, template <typename U> class Descriptor>
void computeShearStress(
    BlockLattice2D<T, Descriptor> &lattice,
    TensorField2D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<TensorField2D<T, SymmetricTensor<T, Descriptor>::n> > computeShearStress(
    BlockLattice2D<T, Descriptor> &lattice);

/* *************** Strain Rate from Stress *************************** */

template <typename T, template <typename U> class Descriptor>
void computeStrainRateFromStress(
    BlockLattice2D<T, Descriptor> &lattice, TensorField2D<T, SymmetricTensor<T, Descriptor>::n> &S);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<TensorField2D<T, SymmetricTensor<T, Descriptor>::n> > computeStrainRateFromStress(
    BlockLattice2D<T, Descriptor> &lattice);

/* *************** Temperature ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeTemperature(BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &temperature);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField2D<T> > computeTemperature(BlockLattice2D<T, Descriptor> &lattice);

/* *************** SoundSpeed ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeSoundSpeed(BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &soundSpeed);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField2D<T> > computeSoundSpeed(BlockLattice2D<T, Descriptor> &lattice);

/* *************** Population **************************************** */

template <typename T, template <typename U> class Descriptor>
void computePopulation(
    BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &population, plint iPop);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField2D<T> > computePopulation(
    BlockLattice2D<T, Descriptor> &lattice, plint iPop);

template <typename T, template <typename U> class Descriptor>
void computeAllPopulations(
    BlockLattice2D<T, Descriptor> &lattice, TensorField2D<T, Descriptor<T>::q> &populations);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<TensorField2D<T, Descriptor<T>::q> > computeAllPopulations(
    BlockLattice2D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
void copyPopulations(
    BlockLattice2D<T, Descriptor> &latticeFrom, BlockLattice2D<T, Descriptor> &latticeTo,
    Box2D domain);

/* ******************************************************************* */
/* *************** PART II. Atomic-block wrappers: Scalar-Field ****** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template <typename T>
T computeSum(ScalarField2D<T> &scalarField, Box2D domain);

template <typename T>
T computeSum(ScalarField2D<T> &scalarField);

template <typename T>
T computeAverage(ScalarField2D<T> &scalarField, Box2D domain);

template <typename T>
T computeAverage(ScalarField2D<T> &scalarField);

template <typename T>
T computeAverage(ScalarField2D<T> &scalarField, ScalarField2D<int> &mask, int flag, Box2D domain);

template <typename T>
T computeAverage(ScalarField2D<T> &scalarField, ScalarField2D<int> &mask, int flag);

template <typename T>
T computeMin(ScalarField2D<T> &scalarField, Box2D domain);

template <typename T>
T computeMin(ScalarField2D<T> &scalarField);

template <typename T>
T computeMax(ScalarField2D<T> &scalarField, Box2D domain);

template <typename T>
T computeMax(ScalarField2D<T> &scalarField);

template <typename T>
T computeBoundedAverage(ScalarField2D<T> &scalarField, Box2D domain);

template <typename T>
T computeBoundedAverage(ScalarField2D<T> &scalarField);

template <typename T, class BoolMask>
plint count(ScalarField2D<T> &field, Box2D domain, BoolMask boolMask);

template <typename T, class BoolMask>
plint count(ScalarField2D<T> &field, BoolMask boolMask);

/* *************** Generic operations *************** */

template <typename T, class Function>
void apply(Function f, ScalarField2D<T> &field)
{
    applyProcessingFunctional(
        new ApplyScalarFunctional2D<T, Function>(f), field.getBoundingBox(), field);
}

template <typename T, class Function>
void evaluate(Function f, ScalarField2D<T> &field, ScalarField2D<T> &result)
{
    applyProcessingFunctional(
        new EvaluateScalarFunctional2D<T, Function>(f), field.getBoundingBox(), field, result);
}

template <typename T, class Function>
std::unique_ptr<ScalarField2D<T> > evaluate(Function f, ScalarField2D<T> &field)
{
    ScalarField2D<T> *result = new ScalarField2D<T>(field.getNx(), field.getNy());
    evaluate(f, field, result);
    return std::unique_ptr<ScalarField2D<T> >(result);
}

/* *************** ScalarField - Scalar operations *************** */

template <typename T>
void add(T scalar, ScalarField2D<T> &field, ScalarField2D<T> &result);

template <typename T>
std::unique_ptr<ScalarField2D<T> > add(T scalar, ScalarField2D<T> &field);

template <typename T>
void add(ScalarField2D<T> &field, T scalar, ScalarField2D<T> &result);

template <typename T>
std::unique_ptr<ScalarField2D<T> > add(ScalarField2D<T> &field, T scalar);

template <typename T>
void subtract(T scalar, ScalarField2D<T> &field, ScalarField2D<T> &result);

template <typename T>
std::unique_ptr<ScalarField2D<T> > subtract(T scalar, ScalarField2D<T> &field);

template <typename T>
void subtract(ScalarField2D<T> &field, T scalar, ScalarField2D<T> &result);

template <typename T>
std::unique_ptr<ScalarField2D<T> > subtract(ScalarField2D<T> &field, T scalar);

template <typename T>
void multiply(T scalar, ScalarField2D<T> &field, ScalarField2D<T> &result);

template <typename T>
std::unique_ptr<ScalarField2D<T> > multiply(T scalar, ScalarField2D<T> &field);

template <typename T>
void multiply(ScalarField2D<T> &field, T scalar, ScalarField2D<T> &result);

template <typename T>
std::unique_ptr<ScalarField2D<T> > multiply(ScalarField2D<T> &field, T scalar);

template <typename T>
void divide(T scalar, ScalarField2D<T> &field, ScalarField2D<T> &result);

template <typename T>
std::unique_ptr<ScalarField2D<T> > divide(T scalar, ScalarField2D<T> &field);

template <typename T>
void divide(ScalarField2D<T> &field, T scalar, ScalarField2D<T> &result);

template <typename T>
std::unique_ptr<ScalarField2D<T> > divide(ScalarField2D<T> &field, T scalar);

/* *************** ScalarField operations *************** */

template <typename T>
void computeSqrt(ScalarField2D<T> &field, ScalarField2D<T> &result, Box2D domain);

template <typename T>
std::unique_ptr<ScalarField2D<T> > computeSqrt(ScalarField2D<T> &field, Box2D domain);

template <typename T>
std::unique_ptr<ScalarField2D<T> > computeSqrt(ScalarField2D<T> &field);

template <typename T>
void computeAbsoluteValue(ScalarField2D<T> &field, ScalarField2D<T> &result, Box2D domain);

template <typename T>
std::unique_ptr<ScalarField2D<T> > computeAbsoluteValue(ScalarField2D<T> &field, Box2D domain);

template <typename T>
std::unique_ptr<ScalarField2D<T> > computeAbsoluteValue(ScalarField2D<T> &field);

/* *************** ScalarField - Scalar inplace operations *************** */

template <typename T>
void addInPlace(ScalarField2D<T> &field, T scalar);

template <typename T>
void subtractInPlace(ScalarField2D<T> &field, T scalar);

template <typename T>
void multiplyInPlace(ScalarField2D<T> &field, T scalar);

template <typename T>
void divideInPlace(ScalarField2D<T> &field, T scalar);

/* *************** ScalarField - ScalarField operations *************** */

template <typename T1, typename T2>
void copy(ScalarField2D<T1> &field, ScalarField2D<T2> &convertedField);

template <typename T1, typename T2>
std::unique_ptr<ScalarField2D<T2> > copyConvert(ScalarField2D<T1> &field);

template <typename T>
void add(ScalarField2D<T> &A, ScalarField2D<T> &B, ScalarField2D<T> &result);

template <typename T>
std::unique_ptr<ScalarField2D<T> > add(ScalarField2D<T> &A, ScalarField2D<T> &B);

template <typename T>
void subtract(ScalarField2D<T> &A, ScalarField2D<T> &B, ScalarField2D<T> &result);

template <typename T>
std::unique_ptr<ScalarField2D<T> > subtract(ScalarField2D<T> &A, ScalarField2D<T> &B);

template <typename T>
void multiply(ScalarField2D<T> &A, ScalarField2D<T> &B, ScalarField2D<T> &result);

template <typename T>
std::unique_ptr<ScalarField2D<T> > multiply(ScalarField2D<T> &A, ScalarField2D<T> &B);

template <typename T>
void divide(ScalarField2D<T> &A, ScalarField2D<T> &B, ScalarField2D<T> &result);

template <typename T>
std::unique_ptr<ScalarField2D<T> > divide(ScalarField2D<T> &A, ScalarField2D<T> &B);

/* *************** ScalarField - ScalarField inplace operations *************** */

template <typename T>
void addInPlace(ScalarField2D<T> &A, ScalarField2D<T> &B);

template <typename T>
void subtractInPlace(ScalarField2D<T> &A, ScalarField2D<T> &B);

template <typename T>
void multiplyInPlace(ScalarField2D<T> &A, ScalarField2D<T> &B);

template <typename T>
void divideInPlace(ScalarField2D<T> &A, ScalarField2D<T> &B);

/* ******************************************************************* */
/* *************** PART III. Atomic-block wrappers: Tensor-field ***** */
/* ******************************************************************* */

/* *************** Reductive functions ********************************* */

template <typename T, int nDim, class BoolMask>
plint count(TensorField2D<T, nDim> &field, Box2D domain, BoolMask boolMask);

template <typename T, int nDim, class BoolMask>
plint count(TensorField2D<T, nDim> &field, BoolMask boolMask);

/* *************** Component (scalar-field) out of a tensor-field ****** */

template <typename T, int nDim>
void extractComponent(
    TensorField2D<T, nDim> &tensorField, ScalarField2D<T> &component, int iComponent);

template <typename T, int nDim>
std::unique_ptr<ScalarField2D<T> > extractComponent(
    TensorField2D<T, nDim> &tensorField, int iComponent);

/* *************** Copy-convert a tensor-field *************** */

template <typename T1, typename T2, int nDim>
void copy(TensorField2D<T1, nDim> &field, TensorField2D<T2, nDim> &convertedField);

template <typename T1, typename T2, int nDim>
std::unique_ptr<TensorField2D<T2, nDim> > copyConvert(TensorField2D<T1, nDim> &field);

/* *************** Vector-norm of each cell in the field *************** */

template <typename T, int nDim>
void computeNorm(TensorField2D<T, nDim> &tensorField, ScalarField2D<T> &norm);

template <typename T, int nDim>
std::unique_ptr<ScalarField2D<T> > computeNorm(TensorField2D<T, nDim> &norm);

/* *************** Squared vector-norm of each cell in the field ******** */

template <typename T, int nDim>
void computeNormSqr(TensorField2D<T, nDim> &tensorField, ScalarField2D<T> &normSqr);

template <typename T, int nDim>
std::unique_ptr<ScalarField2D<T> > computeNormSqr(TensorField2D<T, nDim> &normSqr);

/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template <typename T>
void computeSymmetricTensorNorm(TensorField2D<T, 3> &tensorField, ScalarField2D<T> &norm);

template <typename T>
std::unique_ptr<ScalarField2D<T> > computeSymmetricTensorNorm(TensorField2D<T, 3> &norm);

/* *************** Squared Tensor-norm of each symmetric tensor of a field*/

template <typename T>
void computeSymmetricTensorNormSqr(TensorField2D<T, 3> &tensorField, ScalarField2D<T> &normSqr);

template <typename T>
std::unique_ptr<ScalarField2D<T> > computeSymmetricTensorNormSqr(TensorField2D<T, 3> &normSqr);

/* *************** Trace of each symmetric tensor of a field ************ */

template <typename T>
void computeSymmetricTensorTrace(TensorField2D<T, 3> &tensorField, ScalarField2D<T> &trace);

template <typename T>
std::unique_ptr<ScalarField2D<T> > computeSymmetricTensorTrace(TensorField2D<T, 3> &trace);

/* *************** Vorticity from Velocity field *********************** */

template <typename T>
void computeVorticity(TensorField2D<T, 2> &velocity, ScalarField2D<T> &vorticity);

template <typename T>
std::unique_ptr<ScalarField2D<T> > computeVorticity(TensorField2D<T, 2> &velocity);

/* *************** Vorticity, witout boundary treatment, from Velocity field  */

template <typename T>
void computeBulkVorticity(TensorField2D<T, 2> &velocity, ScalarField2D<T> &vorticity);

template <typename T>
std::unique_ptr<ScalarField2D<T> > computeBulkVorticity(TensorField2D<T, 2> &velocity);

/* *************** Strain Rate from Velocity field ********************* */

template <typename T>
void computeStrainRate(TensorField2D<T, 2> &velocity, TensorField2D<T, 3> &S);

template <typename T>
std::unique_ptr<TensorField2D<T, 3> > computeStrainRate(TensorField2D<T, 2> &velocity);

/* *************** Str. rate, witout boundary treatment, from Velocity field  */

template <typename T>
void computeBulkStrainRate(TensorField2D<T, 2> &velocity, TensorField2D<T, 3> &S);

template <typename T>
std::unique_ptr<TensorField2D<T, 3> > computeBulkStrainRate(TensorField2D<T, 2> &velocity);

/* *************** TensorField - TensorField operations *************** */

template <typename T, int nDim>
void add(TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B, TensorField2D<T, nDim> &result);

template <typename T, int nDim>
std::unique_ptr<TensorField2D<T, nDim> > add(TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B);

template <typename T, int nDim>
void subtract(TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B, TensorField2D<T, nDim> &result);

template <typename T, int nDim>
std::unique_ptr<TensorField2D<T, nDim> > subtract(
    TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B);

template <typename T, int nDim>
void multiply(TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B, TensorField2D<T, nDim> &result);

template <typename T, int nDim>
std::unique_ptr<TensorField2D<T, nDim> > multiply(
    TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B);

template <typename T, int nDim>
void multiply(
    T scalar, TensorField2D<T, nDim> &field, TensorField2D<T, nDim> &result, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<TensorField2D<T, nDim> > multiply(T scalar, TensorField2D<T, nDim> &field);

template <typename T, int nDim>
void multiply(
    TensorField2D<T, nDim> &field, T scalar, TensorField2D<T, nDim> &result, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<TensorField2D<T, nDim> > multiply(TensorField2D<T, nDim> &field, T scalar);

template <typename T, int nDim>
void divide(TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B, TensorField2D<T, nDim> &result);

template <typename T, int nDim>
std::unique_ptr<TensorField2D<T, nDim> > divide(
    TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B);

/* *************** TensorField - TensorField inplace operations *************** */

template <typename T, int nDim>
void addInPlace(TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B);

template <typename T, int nDim>
void subtractInPlace(TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B);

template <typename T, int nDim>
void multiplyInPlace(TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B);

template <typename T, int nDim>
void multiplyInPlace(TensorField2D<T, nDim> &A, T alpha);

template <typename T, int nDim>
void divideInPlace(TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B);

/* ******************************************************************* */
/* *************** PART IV : Multi-block wrappers: Block-Lattice ***** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template <typename T, template <typename U> class Descriptor>
T computeAverageDensity(MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
T computeAverageDensity(MultiBlockLattice2D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(MultiBlockLattice2D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(MultiBlockLattice2D<T, Descriptor> &lattice, DotList2D dotList);

template <typename T, template <typename U> class Descriptor>
T computeAverageEnergy(MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
T computeAverageEnergy(MultiBlockLattice2D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
T computeAverageVelocityComponent(
    MultiBlockLattice2D<T, Descriptor> &lattice, DotList2D dotList, plint iComponent);

template <typename T, template <typename U> class Descriptor, class BoolMask>
plint count(MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, BoolMask boolMask);

template <typename T, template <typename U> class Descriptor, class BoolMask>
plint count(MultiBlockLattice2D<T, Descriptor> &lattice, BoolMask boolMask);

/* *************** Extract Sub-Lattice ******************************* */

template <typename T, template <typename U> class Descriptor>
void extractSubDomain(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiBlockLattice2D<T, Descriptor> &extractedLattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > extractSubDomain(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain);

/* *************** Density ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeDensity(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &density, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeDensity(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeDensity(MultiBlockLattice2D<T, Descriptor> &lattice);

/* *************** RhoBar ******************************************** */

template <typename T, template <typename U> class Descriptor>
void computeRhoBar(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &rhoBar, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeRhoBar(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeRhoBar(MultiBlockLattice2D<T, Descriptor> &lattice);

/* *************** RhoBarJ ******************************************** */

template <typename T, template <typename U> class Descriptor>
void computeRhoBarJ(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &rhoBar,
    MultiTensorField2D<T, 2> &j, Box2D domain);

/* *************** Packed RhoBar J *********************************** */

template <typename T, template <typename U> class Descriptor>
void computePackedRhoBarJ(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiNTensorField2D<T> &rhoBarJ, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiNTensorField2D<T> > computePackedRhoBarJ(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiNTensorField2D<T> > computePackedRhoBarJ(
    MultiBlockLattice2D<T, Descriptor> &lattice);

/* *************** Kinetic Energy ************************************ */

template <typename T, template <typename U> class Descriptor>
void computeKineticEnergy(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &energy, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeKineticEnergy(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeKineticEnergy(
    MultiBlockLattice2D<T, Descriptor> &lattice);

/* *************** Velocity Norm ************************************* */

template <typename T, template <typename U> class Descriptor>
void computeVelocityNorm(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &velocityNorm, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeVelocityNorm(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeVelocityNorm(
    MultiBlockLattice2D<T, Descriptor> &lattice);

/* *************** Velocity Component ******************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocityComponent(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &velocityComponent,
    Box2D domain, plint iComponent);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeVelocityComponent(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, plint iComponent);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeVelocityComponent(
    MultiBlockLattice2D<T, Descriptor> &lattice, plint iComponent);

/* *************** Velocity ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocity(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiTensorField2D<T, Descriptor<T>::d> &velocity,
    Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, Descriptor<T>::d> > computeVelocity(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, Descriptor<T>::d> > computeVelocity(
    MultiBlockLattice2D<T, Descriptor> &lattice);

/* *************** Pi Neq ********************************* */

template <typename T, template <typename U> class Descriptor>
void computePiNeq(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, SymmetricTensor<T, Descriptor>::n> > computePiNeq(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, SymmetricTensor<T, Descriptor>::n> > computePiNeq(
    MultiBlockLattice2D<T, Descriptor> &lattice);

/* *************** ShearStress ********************************* */

template <typename T, template <typename U> class Descriptor>
void computeShearStress(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, SymmetricTensor<T, Descriptor>::n> > computeShearStress(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, SymmetricTensor<T, Descriptor>::n> > computeShearStress(
    MultiBlockLattice2D<T, Descriptor> &lattice);

/* *************** Strain Rate from Stress *************************** */

template <typename T, template <typename U> class Descriptor>
void computeStrainRateFromStress(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, SymmetricTensor<T, Descriptor>::n> &S, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, SymmetricTensor<T, Descriptor>::n> >
    computeStrainRateFromStress(MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, SymmetricTensor<T, Descriptor>::n> >
    computeStrainRateFromStress(MultiBlockLattice2D<T, Descriptor> &lattice);

/* *************** Temperature ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeTemperature(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &temperature, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeTemperature(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeTemperature(
    MultiBlockLattice2D<T, Descriptor> &lattice);

/* *************** SoundSpeed ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeSoundSpeed(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &soundSpeed, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeSoundSpeed(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeSoundSpeed(
    MultiBlockLattice2D<T, Descriptor> &lattice);

/* *************** Population **************************************** */

template <typename T, template <typename U> class Descriptor>
void computePopulation(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &population, Box2D domain,
    plint iPop);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computePopulation(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, plint iPop);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computePopulation(
    MultiBlockLattice2D<T, Descriptor> &lattice, plint iPop);

template <typename T, template <typename U> class Descriptor>
void computeEquilibrium(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &equilibrium, Box2D domain,
    plint iPop);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeEquilibrium(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, plint iPop);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeEquilibrium(
    MultiBlockLattice2D<T, Descriptor> &lattice, plint iPop);

template <typename T, template <typename U> class Descriptor>
void computeAllPopulations(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::q> &populations, Box2D domain);

template <typename T, template <typename U> class Descriptor>
void computeAllPopulations(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::q> &populations);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, Descriptor<T>::q> > computeAllPopulations(
    MultiBlockLattice2D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, Descriptor<T>::q> > computeAllPopulations(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
void copyPopulations(
    MultiBlockLattice2D<T, Descriptor> &latticeFrom, MultiBlockLattice2D<T, Descriptor> &latticeTo,
    Box2D domain);

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void copyConvertPopulations(
    BlockLattice2D<T1, Descriptor1> &latticeFrom, BlockLattice2D<T2, Descriptor2> &latticeTo,
    Box2D domain);

template <typename T, template <typename U> class Descriptor>
void copyAll(
    MultiBlockLattice2D<T, Descriptor> &latticeFrom, MultiBlockLattice2D<T, Descriptor> &latticeTo,
    Box2D domain);

template <typename T, template <typename U> class Descriptor>
void copyRegenerate(
    MultiBlockLattice2D<T, Descriptor> &latticeFrom, MultiBlockLattice2D<T, Descriptor> &latticeTo,
    Box2D domain);

/* *************** Omega ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeOmega(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &omega, Box2D domain);

template <typename T, template <typename U> class Descriptor>
void computeOmega(MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &omega);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeOmega(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeOmega(MultiBlockLattice2D<T, Descriptor> &lattice);

/* *************** ExternalForce ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeExternalForce(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiTensorField2D<T, Descriptor<T>::d> &force,
    Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, Descriptor<T>::d> > computeExternalForce(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, Descriptor<T>::d> > computeExternalForce(
    MultiBlockLattice2D<T, Descriptor> &lattice);

/* *************** ExternalScalar ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeExternalScalar(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &scalar, int whichScalar,
    Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeExternalScalar(
    MultiBlockLattice2D<T, Descriptor> &lattice, int whichScalar, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeExternalScalar(
    MultiBlockLattice2D<T, Descriptor> &lattice, int whichScalar);

/* *************** ExternalVector ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeExternalVector(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::d> &tensorField, int vectorBeginsAt, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, Descriptor<T>::d> > computeExternalVector(
    MultiBlockLattice2D<T, Descriptor> &lattice, int vectorBeginsAt, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, Descriptor<T>::d> > computeExternalVector(
    MultiBlockLattice2D<T, Descriptor> &lattice, int vectorBeginsAt);

/* ******************************************************************* */
/* *************** PART V : Multi-block wrappers: Scalar-Field ******* */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template <typename T>
T computeSum(MultiScalarField2D<T> &scalarField, Box2D domain);

template <typename T>
T computeSum(MultiScalarField2D<T> &scalarField);

template <typename T>
T computeAverage(MultiScalarField2D<T> &scalarField, Box2D domain);

template <typename T>
T computeAverage(MultiScalarField2D<T> &scalarField);

template <typename T>
T computeAverage(
    MultiScalarField2D<T> &scalarField, MultiScalarField2D<int> &mask, int flag, Box2D domain);

template <typename T>
T computeAverage(MultiScalarField2D<T> &scalarField, MultiScalarField2D<int> &mask, int flag);

template <typename T>
T computeMin(MultiScalarField2D<T> &scalarField, Box2D domain);

template <typename T>
T computeMin(MultiScalarField2D<T> &scalarField);

template <typename T>
T computeMin(
    MultiScalarField2D<T> &scalarField, MultiScalarField2D<int> &mask, int flag, Box2D domain);

template <typename T>
T computeMin(MultiScalarField2D<T> &scalarField, MultiScalarField2D<int> &mask, int flag);

template <typename T>
T computeMax(MultiScalarField2D<T> &scalarField, Box2D domain);

template <typename T>
T computeMax(MultiScalarField2D<T> &scalarField);

template <typename T>
T computeMax(
    MultiScalarField2D<T> &scalarField, MultiScalarField2D<int> &mask, int flag, Box2D domain);

template <typename T>
T computeMax(MultiScalarField2D<T> &scalarField, MultiScalarField2D<int> &mask, int flag);

template <typename T>
T computeBoundedAverage(MultiScalarField2D<T> &scalarField, Box2D domain);

template <typename T>
T computeBoundedAverage(MultiScalarField2D<T> &scalarField);

template <typename T, class BoolMask>
plint count(MultiScalarField2D<T> &field, Box2D domain, BoolMask boolMask);

template <typename T, class BoolMask>
plint count(MultiScalarField2D<T> &field, BoolMask boolMask);

/* *************** Generic operations *************** */

template <typename T, class Function>
void apply(Function f, MultiScalarField2D<T> &field, Box2D domain)
{
    applyProcessingFunctional(new ApplyScalarFunctional2D<T, Function>(f), domain, field);
}

template <typename T, class Function>
void apply(Function f, MultiScalarField2D<T> &field)
{
    apply(f, field, field.getBoundingBox());
}

template <typename T, class Function>
void evaluate(Function f, MultiScalarField2D<T> &field, MultiScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(
        new EvaluateScalarFunctional2D<T, Function>(f), domain, field, result);
}

template <typename T, class Function>
std::unique_ptr<MultiScalarField2D<T> > evaluate(
    Function f, MultiScalarField2D<T> &field, Box2D domain)
{
    MultiScalarField2D<T> *result = new MultiScalarField2D<T>(field, domain);
    evaluate(f, field, *result, domain);
    return std::unique_ptr<MultiScalarField2D<T> >(result);
}

template <typename T, class Function>
std::unique_ptr<MultiScalarField2D<T> > evaluate(Function f, MultiScalarField2D<T> &field)
{
    return evaluate(f, field, field.getBoundingBox());
}

/* *************** Extract Sub-ScalarField *************************** */

template <typename T>
void extractSubDomain(
    MultiScalarField2D<T> &field, MultiScalarField2D<T> &extractedField, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > extractSubDomain(
    MultiScalarField2D<T> &field, Box2D domain);

/* *************** Copy and Convert ScalarField *************************** */

template <typename T1, typename T2>
void copy(MultiScalarField2D<T1> &field, MultiScalarField2D<T2> &convertedField, Box2D domain);

template <typename T1, typename T2>
std::unique_ptr<MultiScalarField2D<T2> > copyConvert(MultiScalarField2D<T1> &field, Box2D domain);

template <typename T1, typename T2>
std::unique_ptr<MultiScalarField2D<T2> > copyConvert(MultiScalarField2D<T1> &field);

/* *************** MultiScalarField - Scalar operations *************** */

template <typename T>
void lessThan(
    MultiScalarField2D<T> &field, T scalar, MultiScalarField2D<int> &result, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<int> > lessThan(
    MultiScalarField2D<T> &field, T scalar, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<int> > lessThan(MultiScalarField2D<T> &field, T scalar);

template <typename T>
void greaterThan(
    MultiScalarField2D<T> &field, T scalar, MultiScalarField2D<int> &result, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<int> > greaterThan(
    MultiScalarField2D<T> &field, T scalar, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<int> > greaterThan(MultiScalarField2D<T> &field, T scalar);

template <typename T>
void add(T scalar, MultiScalarField2D<T> &field, MultiScalarField2D<T> &result, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > add(T scalar, MultiScalarField2D<T> &field, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > add(T scalar, MultiScalarField2D<T> &field);

template <typename T>
void add(MultiScalarField2D<T> &field, T scalar, MultiScalarField2D<T> &result, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > add(MultiScalarField2D<T> &field, T scalar, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > add(MultiScalarField2D<T> &field, T scalar);

template <typename T>
void subtract(T scalar, MultiScalarField2D<T> &field, MultiScalarField2D<T> &result, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > subtract(
    T scalar, MultiScalarField2D<T> &field, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > subtract(T scalar, MultiScalarField2D<T> &field);

template <typename T>
void subtract(MultiScalarField2D<T> &field, T scalar, MultiScalarField2D<T> &result, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > subtract(
    MultiScalarField2D<T> &field, T scalar, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > subtract(MultiScalarField2D<T> &field, T scalar);

template <typename T>
void multiply(T scalar, MultiScalarField2D<T> &field, MultiScalarField2D<T> &result, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > multiply(
    T scalar, MultiScalarField2D<T> &field, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > multiply(T scalar, MultiScalarField2D<T> &field);

template <typename T>
void multiply(MultiScalarField2D<T> &field, T scalar, MultiScalarField2D<T> &result, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > multiply(
    MultiScalarField2D<T> &field, T scalar, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > multiply(MultiScalarField2D<T> &field, T scalar);

template <typename T>
void divide(T scalar, MultiScalarField2D<T> &field, MultiScalarField2D<T> &result, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > divide(
    T scalar, MultiScalarField2D<T> &field, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > divide(T scalar, MultiScalarField2D<T> &field);

template <typename T>
void divide(MultiScalarField2D<T> &field, T scalar, MultiScalarField2D<T> &result, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > divide(
    MultiScalarField2D<T> &field, T scalar, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > divide(MultiScalarField2D<T> &field, T scalar);

/* *************** MultiScalarField operations *************** */

template <typename T>
void computeSqrt(MultiScalarField2D<T> &field, MultiScalarField2D<T> &result, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeSqrt(MultiScalarField2D<T> &field, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeSqrt(MultiScalarField2D<T> &field);

template <typename T>
void computeLog(MultiScalarField2D<T> &field, MultiScalarField2D<T> &result, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeLog(MultiScalarField2D<T> &field, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeLog(MultiScalarField2D<T> &field);

template <typename T>
void computeAbsoluteValue(
    MultiScalarField2D<T> &field, MultiScalarField2D<T> &result, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeAbsoluteValue(
    MultiScalarField2D<T> &field, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeAbsoluteValue(MultiScalarField2D<T> &field);

template <typename T>
void computePower(
    MultiScalarField2D<T> &field, MultiScalarField2D<T> &result, T power, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeSqrt(
    MultiScalarField2D<T> &field, T power, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeSqrt(MultiScalarField2D<T> &field, T power);

/* *************** MultiScalarField - Scalar inplace operations *************** */

template <typename T>
void addInPlace(MultiScalarField2D<T> &field, T scalar, Box2D domain);

template <typename T>
void addInPlace(MultiScalarField2D<T> &field, T scalar);

template <typename T>
void subtractInPlace(MultiScalarField2D<T> &field, T scalar, Box2D domain);

template <typename T>
void subtractInPlace(MultiScalarField2D<T> &field, T scalar);

template <typename T>
void multiplyInPlace(MultiScalarField2D<T> &field, T scalar, Box2D domain);

template <typename T>
void multiplyInPlace(MultiScalarField2D<T> &field, T scalar);

template <typename T>
void divideInPlace(MultiScalarField2D<T> &field, T scalar, Box2D domain);

template <typename T>
void divideInPlace(MultiScalarField2D<T> &field, T scalar);

/* *************** MultiScalarField - MultiScalarField operations *************** */

template <typename T>
void lessThan(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, MultiScalarField2D<int> &result,
    Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<int> > lessThan(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<int> > lessThan(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B);

template <typename T>
void greaterThan(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, MultiScalarField2D<int> &result,
    Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<int> > greaterThan(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<int> > greaterThan(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B);

template <typename T>
void add(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, MultiScalarField2D<T> &result,
    Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > add(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > add(MultiScalarField2D<T> &A, MultiScalarField2D<T> &B);

template <typename T>
void subtract(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, MultiScalarField2D<T> &result,
    Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > subtract(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > subtract(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B);

template <typename T>
void multiply(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, MultiScalarField2D<T> &result,
    Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > multiply(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > multiply(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B);

template <typename T>
void divide(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, MultiScalarField2D<T> &result,
    Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > divide(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > divide(MultiScalarField2D<T> &A, MultiScalarField2D<T> &B);

/* *************** MultiScalarField - MultiScalarField inplace operations *************** */

template <typename T>
void addInPlace(MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, Box2D domain);

template <typename T>
void addInPlace(MultiScalarField2D<T> &A, MultiScalarField2D<T> &B);

template <typename T>
void subtractInPlace(MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, Box2D domain);

template <typename T>
void subtractInPlace(MultiScalarField2D<T> &A, MultiScalarField2D<T> &B);

template <typename T>
void multiplyInPlace(MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, Box2D domain);

template <typename T>
void multiplyInPlace(MultiScalarField2D<T> &A, MultiScalarField2D<T> &B);

template <typename T>
void divideInPlace(MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, Box2D domain);

template <typename T>
void divideInPlace(MultiScalarField2D<T> &A, MultiScalarField2D<T> &B);

/* ******************************************************************* */
/* *************** PART VI : Multi-block wrappers: Tensor-field ****** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template <typename T, int nDim, class BoolMask>
plint count(MultiTensorField2D<T, nDim> &field, Box2D domain, BoolMask boolMask);

template <typename T, int nDim, class BoolMask>
plint count(MultiTensorField2D<T, nDim> &field, BoolMask boolMask);

template <typename T, int nDim>
Array<T, nDim> computeAverage(MultiTensorField2D<T, nDim> &field, Box2D domain);

template <typename T, int nDim>
Array<T, nDim> computeAverage(MultiTensorField2D<T, nDim> &field);

/* *************** Copy-convert a tensor-field *************** */

template <typename T1, typename T2, int nDim>
void copy(
    MultiTensorField2D<T1, nDim> &field, MultiTensorField2D<T2, nDim> &convertedField,
    Box2D domain);

template <typename T1, typename T2, int nDim>
std::unique_ptr<MultiTensorField2D<T2, nDim> > copyConvert(
    MultiTensorField2D<T1, nDim> &field, Box2D domain);

template <typename T1, typename T2, int nDim>
std::unique_ptr<MultiTensorField2D<T2, nDim> > copyConvert(MultiTensorField2D<T1, nDim> &field);

/* *************** Extract Sub-TensorField *************************** */

template <typename T, int nDim>
void extractSubDomain(
    MultiTensorField2D<T, nDim> &field, MultiTensorField2D<T, nDim> &extractedField, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > extractSubDomain(
    MultiTensorField2D<T, nDim> &field, Box2D domain);

/* *************** Component (scalar-field) out of a tensor-field ****** */

template <typename T, int nDim>
void extractComponent(
    MultiTensorField2D<T, nDim> &tensorField, MultiScalarField2D<T> &component, Box2D domain,
    int iComponent);

template <typename T, int nDim>
std::unique_ptr<MultiScalarField2D<T> > extractComponent(
    MultiTensorField2D<T, nDim> &tensorField, Box2D domain, int iComponent);

template <typename T, int nDim>
std::unique_ptr<MultiScalarField2D<T> > extractComponent(
    MultiTensorField2D<T, nDim> &tensorField, int iComponent);

/* *************** Vector-norm of each cell in the field *************** */

template <typename T, int nDim>
void computeNorm(
    MultiTensorField2D<T, nDim> &tensorField, MultiScalarField2D<T> &norm, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiScalarField2D<T> > computeNorm(
    MultiTensorField2D<T, nDim> &tensorField, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiScalarField2D<T> > computeNorm(MultiTensorField2D<T, nDim> &tensorField);

/* *************** Squared vector-norm of each cell in the field ******** */

template <typename T, int nDim>
void computeNormSqr(
    MultiTensorField2D<T, nDim> &tensorField, MultiScalarField2D<T> &normSqr, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiScalarField2D<T> > computeNormSqr(
    MultiTensorField2D<T, nDim> &tensorField, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiScalarField2D<T> > computeNormSqr(MultiTensorField2D<T, nDim> &tensorField);

/* *************** Max element of each array of each cell *************** */

template <typename T, int nDim>
void computeMaximumElement(
    MultiTensorField2D<T, nDim> &A, MultiScalarField2D<T> &result, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiScalarField2D<T> > computeMaximumElement(
    MultiTensorField2D<T, nDim> &A, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiScalarField2D<T> > computeMaximumElement(MultiTensorField2D<T, nDim> &A);

/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template <typename T>
void computeSymmetricTensorNorm(
    MultiTensorField2D<T, 3> &tensorField, MultiScalarField2D<T> &norm, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeSymmetricTensorNorm(
    MultiTensorField2D<T, 3> &tensorField, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeSymmetricTensorNorm(
    MultiTensorField2D<T, 3> &tensorField);

/* *************** Squared Tensor-norm of each symmetric tensor of a field*/

template <typename T>
void computeSymmetricTensorNormSqr(
    MultiTensorField2D<T, 3> &tensorField, MultiScalarField2D<T> &normSqr, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeSymmetricTensorNormSqr(
    MultiTensorField2D<T, 3> &tensorField, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeSymmetricTensorNormSqr(
    MultiTensorField2D<T, 3> &tensorField);

/* *************** Trace of each symmetric tensor of a field ************ */

template <typename T>
void computeSymmetricTensorTrace(
    MultiTensorField2D<T, 3> &tensorField, MultiScalarField2D<T> &trace, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeSymmetricTensorTrace(
    MultiTensorField2D<T, 3> &tensorField, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeSymmetricTensorTrace(
    MultiTensorField2D<T, 3> &tensorField);

/* *************** Gradient from scalar field *********************** */

template <typename T>
void computeGradient(MultiScalarField2D<T> &phi, MultiTensorField2D<T, 2> &gradient, Box2D domain);

template <typename T>
std::unique_ptr<MultiTensorField2D<T, 2> > computeGradient(
    MultiScalarField2D<T> &phi, Box2D domain);

template <typename T>
std::unique_ptr<MultiTensorField2D<T, 2> > computeGradient(MultiScalarField2D<T> &phi);

/* *************** Vorticity from Velocity field *********************** */

template <typename T>
void computeVorticity(
    MultiTensorField2D<T, 2> &velocity, MultiScalarField2D<T> &vorticity, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeVorticity(
    MultiTensorField2D<T, 2> &velocity, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeVorticity(MultiTensorField2D<T, 2> &velocity);

/* *************** Vorticity, witout boundary treatment, from Velocity field  */

template <typename T>
void computeBulkVorticity(
    MultiTensorField2D<T, 2> &velocity, MultiScalarField2D<T> &vorticity, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeBulkVorticity(
    MultiTensorField2D<T, 2> &velocity, Box2D domain);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeBulkVorticity(MultiTensorField2D<T, 2> &velocity);

/* *************** Strain rate from Velocity field ********************* */

template <typename T>
void computeStrainRate(
    MultiTensorField2D<T, 2> &velocity, MultiTensorField2D<T, 3> &S, Box2D domain);

template <typename T>
std::unique_ptr<MultiTensorField2D<T, 3> > computeStrainRate(
    MultiTensorField2D<T, 2> &velocity, Box2D domain);

template <typename T>
std::unique_ptr<MultiTensorField2D<T, 3> > computeStrainRate(MultiTensorField2D<T, 2> &velocity);

/* *************** Str. rate, witout boundary treatment, from Velocity field  */

template <typename T>
void computeBulkStrainRate(
    MultiTensorField2D<T, 2> &velocity, MultiTensorField2D<T, 3> &S, Box2D domain);

template <typename T>
std::unique_ptr<MultiTensorField2D<T, 3> > computeBulkStrainRate(
    MultiTensorField2D<T, 2> &velocity, Box2D domain);

template <typename T>
std::unique_ptr<MultiTensorField2D<T, 3> > computeBulkStrainRate(
    MultiTensorField2D<T, 2> &velocity);

/* *************** MultiTensorField - MultiTensorField operations *************** */

template <typename T, int nDim>
void add(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B,
    MultiTensorField2D<T, nDim> &result, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > add(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > add(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B);

template <typename T, int nDim>
void subtract(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B,
    MultiTensorField2D<T, nDim> &result, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > subtract(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > subtract(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B);

template <typename T, int nDim>
void multiply(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B,
    MultiTensorField2D<T, nDim> &result, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > multiply(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > multiply(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B);

template <typename T, int nDim>
void multiply(
    T scalar, MultiTensorField2D<T, nDim> &field, MultiTensorField2D<T, nDim> &result,
    Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > multiply(
    T scalar, MultiTensorField2D<T, nDim> &field, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > multiply(
    T scalar, MultiTensorField2D<T, nDim> &field);

template <typename T, int nDim>
void multiply(
    MultiTensorField2D<T, nDim> &field, T scalar, MultiTensorField2D<T, nDim> &result,
    Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > multiply(
    MultiTensorField2D<T, nDim> &field, T scalar, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > multiply(
    MultiTensorField2D<T, nDim> &field, T scalar);

template <typename T, int nDim>
void divide(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B,
    MultiTensorField2D<T, nDim> &result, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > divide(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > divide(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B);

/* *************** MultiTensorField - MultiTensorField inplace operations *************** */

template <typename T, int nDim>
void addInPlace(MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B, Box2D domain);

template <typename T, int nDim>
void addInPlace(MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B);

template <typename T, int nDim>
void subtractInPlace(MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B, Box2D domain);

template <typename T, int nDim>
void subtractInPlace(MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B);

template <typename T, int nDim>
void multiplyInPlace(MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B, Box2D domain);

template <typename T, int nDim>
void multiplyInPlace(MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B);

template <typename T, int nDim>
void multiplyInPlace(MultiTensorField2D<T, nDim> &A, T alpha, Box2D domain);

template <typename T, int nDim>
void multiplyInPlace(MultiTensorField2D<T, nDim> &A, T alpha);

template <typename T, int nDim>
void divideInPlace(MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B, Box2D domain);

template <typename T, int nDim>
void divideInPlace(MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B);

}  // namespace plb

#endif  // DATA_ANALYSIS_WRAPPER_2D_H
