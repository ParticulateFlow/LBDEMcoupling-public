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

#ifndef MULTI_GRID_DATA_ANALYSIS_WRAPPER_2D_H
#define MULTI_GRID_DATA_ANALYSIS_WRAPPER_2D_H

#include "dataProcessors/dataAnalysisWrapper2D.h"
#include "multiGrid/multiGridDataField2D.h"
#include "multiGrid/multiGridLattice2D.h"

namespace plb {

/* *************** Reductive functions ******************************* */

template <typename T, template <typename U> class Descriptor>
T computeAverageDensity(MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
T computeAverageDensity(MultiGridLattice2D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(MultiGridLattice2D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
T computeAverageEnergy(MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
T computeAverageEnergy(MultiGridLattice2D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor, class BoolMask>
plint count(MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain, BoolMask boolMask);

template <typename T, template <typename U> class Descriptor, class BoolMask>
plint count(MultiGridLattice2D<T, Descriptor> &lattice, BoolMask boolMask);

/* *************** Extract Sub-Lattice ******************************* */

template <typename T, template <typename U> class Descriptor>
void extractSubDomain(
    MultiGridLattice2D<T, Descriptor> &lattice, MultiGridLattice2D<T, Descriptor> &extractedLattice,
    Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridLattice2D<T, Descriptor> > extractSubDomain(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain);

/* *************** Density ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeDensity(
    MultiGridLattice2D<T, Descriptor> &lattice, MultiGridScalarField2D<T> &density, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeDensity(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeDensity(
    MultiGridLattice2D<T, Descriptor> &lattice);

/* *************** RhoBar ******************************************** */

template <typename T, template <typename U> class Descriptor>
void computeDensity(
    MultiGridLattice2D<T, Descriptor> &lattice, MultiGridScalarField2D<T> &rhoBar, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeDensity(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeDensity(
    MultiGridLattice2D<T, Descriptor> &lattice);

/* *************** Kinetic Energy ************************************ */

template <typename T, template <typename U> class Descriptor>
void computeKineticEnergy(
    MultiGridLattice2D<T, Descriptor> &lattice, MultiGridScalarField2D<T> &energy, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeKineticEnergy(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeKineticEnergy(
    MultiGridLattice2D<T, Descriptor> &lattice);

/* *************** Velocity Norm ************************************* */

template <typename T, template <typename U> class Descriptor>
void computeVelocityNorm(
    MultiGridLattice2D<T, Descriptor> &lattice, MultiGridScalarField2D<T> &velocityNorm,
    Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeVelocityNorm(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeVelocityNorm(
    MultiGridLattice2D<T, Descriptor> &lattice);

/* *************** Velocity Component ******************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocityComponent(
    MultiGridLattice2D<T, Descriptor> &lattice, MultiGridScalarField2D<T> &velocityComponent,
    Box2D domain, plint iComponent);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeVelocityComponent(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain, plint iComponent);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeVelocityComponent(
    MultiGridLattice2D<T, Descriptor> &lattice, plint iComponent);

/* *************** Velocity ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocity(
    MultiGridLattice2D<T, Descriptor> &lattice,
    MultiGridTensorField2D<T, Descriptor<T>::d> &velocity, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField2D<T, Descriptor<T>::d> > computeVelocity(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField2D<T, Descriptor<T>::d> > computeVelocity(
    MultiGridLattice2D<T, Descriptor> &lattice);

/* *************** Pi Neq ********************************* */

template <typename T, template <typename U> class Descriptor>
void computePiNeq(
    MultiGridLattice2D<T, Descriptor> &lattice,
    MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> > computePiNeq(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> > computePiNeq(
    MultiGridLattice2D<T, Descriptor> &lattice);

/* *************** Shear Stress ********************************* */

template <typename T, template <typename U> class Descriptor>
void computeShearStress(
    MultiGridLattice2D<T, Descriptor> &lattice,
    MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> > computeShearStress(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> > computeShearStress(
    MultiGridLattice2D<T, Descriptor> &lattice);

/* *************** Strain Rate from Stress *************************** */

template <typename T, template <typename U> class Descriptor>
void computeStrainRateFromStress(
    MultiGridLattice2D<T, Descriptor> &lattice,
    MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> &S, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> >
    computeStrainRateFromStress(MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> >
    computeStrainRateFromStress(MultiGridLattice2D<T, Descriptor> &lattice);

/* *************** Temperature ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeTemperature(
    MultiGridLattice2D<T, Descriptor> &lattice, MultiGridScalarField2D<T> &temperature,
    Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeTemperature(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeTemperature(
    MultiGridLattice2D<T, Descriptor> &lattice);

/* *************** SoundSpeed ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeSoundSpeed(
    MultiGridLattice2D<T, Descriptor> &lattice, MultiGridScalarField2D<T> &soundSpeed,
    Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeSoundSpeed(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeSoundSpeed(
    MultiGridLattice2D<T, Descriptor> &lattice);

/* *************** Population **************************************** */

template <typename T, template <typename U> class Descriptor>
void computePopulation(
    MultiGridLattice2D<T, Descriptor> &lattice, MultiGridScalarField2D<T> &population, Box2D domain,
    plint iPop);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computePopulation(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain, plint iPop);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computePopulation(
    MultiGridLattice2D<T, Descriptor> &lattice, plint iPop);

template <typename T, template <typename U> class Descriptor>
void computeAllPopulations(
    MultiGridLattice2D<T, Descriptor> &lattice,
    MultiGridTensorField2D<T, Descriptor<T>::q> &populations);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField2D<T, Descriptor<T>::q> > computeAllPopulations(
    MultiGridLattice2D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField2D<T, Descriptor<T>::q> > computeAllPopulations(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
void copyPopulations(
    MultiGridLattice2D<T, Descriptor> &latticeFrom, MultiGridLattice2D<T, Descriptor> &latticeTo,
    Box2D domain);

/* *************** Reductive functions ******************************* */

template <typename T>
T computeAverage(MultiGridScalarField2D<T> &scalarField, Box2D domain);

template <typename T>
T computeAverage(MultiGridScalarField2D<T> &scalarField);

template <typename T>
T computeMin(MultiGridScalarField2D<T> &scalarField, Box2D domain);

template <typename T>
T computeMin(MultiGridScalarField2D<T> &scalarField);

template <typename T>
T computeMax(MultiGridScalarField2D<T> &scalarField, Box2D domain);

template <typename T>
T computeMax(MultiGridScalarField2D<T> &scalarField);

template <typename T>
T computeBoundedAverage(MultiGridScalarField2D<T> &scalarField, Box2D domain);

template <typename T>
T computeBoundedAverage(MultiGridScalarField2D<T> &scalarField);

template <typename T, class BoolMask>
plint count(MultiGridScalarField2D<T> &field, Box2D domain, BoolMask boolMask);

template <typename T, class BoolMask>
plint count(MultiGridScalarField2D<T> &field, BoolMask boolMask);

/* *************** Generic operations *************** */

template <typename T, class Function>
void apply(Function f, MultiGridScalarField2D<T> &field, Box2D domain)
{
    applyProcessingFunctional(new ApplyScalarFunctional2D<T, Function>(f), domain, field);
}

template <typename T, class Function>
void apply(Function f, MultiGridScalarField2D<T> &field)
{
    apply(f, field, field.getBoundingBox());
}

template <typename T, class Function>
void evaluate(
    Function f, MultiGridScalarField2D<T> &field, MultiGridScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(
        new EvaluateScalarFunctional2D<T, Function>(f), domain, field, result);
}

template <typename T, class Function>
std::unique_ptr<MultiGridScalarField2D<T> > evaluate(
    Function f, MultiGridScalarField2D<T> &field, Box2D domain)
{
    MultiGridScalarField2D<T> *result = new MultiGridScalarField2D<T>(field, domain);
    evaluate(f, field, *result, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(result);
}

template <typename T, class Function>
std::unique_ptr<MultiGridScalarField2D<T> > evaluate(Function f, MultiGridScalarField2D<T> &field)
{
    return evaluate(f, field, field.getBoundingBox());
}

/* *************** Extract Sub-ScalarField *************************** */

template <typename T>
void extractSubDomain(
    MultiGridScalarField2D<T> &field, MultiGridScalarField2D<T> &extractedField, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > extractSubDomain(
    MultiGridScalarField2D<T> &field, Box2D domain);

/* *************** Copy and Convert ScalarField *************************** */

template <typename T1, typename T2>
void copyConvert(
    MultiGridScalarField2D<T1> &field, MultiGridScalarField2D<T2> &convertedField, Box2D domain);

template <typename T1, typename T2>
std::unique_ptr<MultiGridScalarField2D<T2> > copyConvert(
    MultiGridScalarField2D<T1> &field, Box2D domain);

template <typename T1, typename T2>
std::unique_ptr<MultiGridScalarField2D<T2> > copyConvert(MultiGridScalarField2D<T1> &field);

/* *************** MultiScalarField - Scalar operations *************** */

template <typename T>
void add(
    T scalar, MultiGridScalarField2D<T> &field, MultiGridScalarField2D<T> &result, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > add(
    T scalar, MultiGridScalarField2D<T> &field, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > add(T scalar, MultiGridScalarField2D<T> &field);

template <typename T>
void add(
    MultiGridScalarField2D<T> &field, T scalar, MultiGridScalarField2D<T> &result, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > add(
    MultiGridScalarField2D<T> &field, T scalar, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > add(MultiGridScalarField2D<T> &field, T scalar);

template <typename T>
void subtract(
    T scalar, MultiGridScalarField2D<T> &field, MultiGridScalarField2D<T> &result, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > subtract(
    T scalar, MultiGridScalarField2D<T> &field, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > subtract(T scalar, MultiGridScalarField2D<T> &field);

template <typename T>
void subtract(
    MultiGridScalarField2D<T> &field, T scalar, MultiGridScalarField2D<T> &result, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > subtract(
    MultiGridScalarField2D<T> &field, T scalar, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > subtract(MultiGridScalarField2D<T> &field, T scalar);

template <typename T>
void multiply(
    T scalar, MultiGridScalarField2D<T> &field, MultiGridScalarField2D<T> &result, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > multiply(
    T scalar, MultiGridScalarField2D<T> &field, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > multiply(T scalar, MultiGridScalarField2D<T> &field);

template <typename T>
void multiply(
    MultiGridScalarField2D<T> &field, T scalar, MultiGridScalarField2D<T> &result, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > multiply(
    MultiGridScalarField2D<T> &field, T scalar, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > multiply(MultiGridScalarField2D<T> &field, T scalar);

template <typename T>
void divide(
    T scalar, MultiGridScalarField2D<T> &field, MultiGridScalarField2D<T> &result, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > divide(
    T scalar, MultiGridScalarField2D<T> &field, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > divide(T scalar, MultiGridScalarField2D<T> &field);

template <typename T>
void divide(
    MultiGridScalarField2D<T> &field, T scalar, MultiGridScalarField2D<T> &result, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > divide(
    MultiGridScalarField2D<T> &field, T scalar, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > divide(MultiGridScalarField2D<T> &field, T scalar);

/* *************** MultiScalarField - Scalar inplace operations *************** */

template <typename T>
void addInPlace(MultiGridScalarField2D<T> &field, T scalar, Box2D domain);

template <typename T>
void addInPlace(MultiGridScalarField2D<T> &field, T scalar);

template <typename T>
void subtractInPlace(MultiGridScalarField2D<T> &field, T scalar, Box2D domain);

template <typename T>
void subtractInPlace(MultiGridScalarField2D<T> &field, T scalar);

template <typename T>
void multiplyInPlace(MultiGridScalarField2D<T> &field, T scalar, Box2D domain);

template <typename T>
void multiplyInPlace(MultiGridScalarField2D<T> &field, T scalar);

template <typename T>
void divideInPlace(MultiGridScalarField2D<T> &field, T scalar, Box2D domain);

template <typename T>
void divideInPlace(MultiGridScalarField2D<T> &field, T scalar);

/* *************** MultiScalarField - MultiScalarField operations *************** */

template <typename T>
void add(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, MultiGridScalarField2D<T> &result,
    Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > add(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > add(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B);

template <typename T>
void subtract(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, MultiGridScalarField2D<T> &result,
    Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > subtract(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > subtract(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B);

template <typename T>
void multiply(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, MultiGridScalarField2D<T> &result,
    Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > multiply(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > multiply(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B);

template <typename T>
void divide(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, MultiGridScalarField2D<T> &result,
    Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > divide(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > divide(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B);

/* *************** MultiScalarField - MultiScalarField inplace operations *************** */

template <typename T>
void addInPlace(MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, Box2D domain);

template <typename T>
void addInPlace(MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B);

template <typename T>
void subtractInPlace(MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, Box2D domain);

template <typename T>
void subtractInPlace(MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B);

template <typename T>
void multiplyInPlace(MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, Box2D domain);

template <typename T>
void multiplyInPlace(MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B);

template <typename T>
void divideInPlace(MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, Box2D domain);

template <typename T>
void divideInPlace(MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B);

/* ******************************************************************* */
/* ***************           Tensor-field                       ****** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template <typename T, int nDim, class BoolMask>
plint count(MultiGridTensorField2D<T, nDim> &field, Box2D domain, BoolMask boolMask);

template <typename T, int nDim, class BoolMask>
plint count(MultiGridTensorField2D<T, nDim> &field, BoolMask boolMask);

/* *************** Copy-convert a tensor-field *************** */

template <typename T1, typename T2, int nDim>
void copyConvert(
    MultiGridTensorField2D<T1, nDim> &field, MultiGridTensorField2D<T2, nDim> &convertedField,
    Box2D domain);

template <typename T1, typename T2, int nDim>
std::unique_ptr<MultiGridTensorField2D<T2, nDim> > copyConvert(
    MultiGridTensorField2D<T1, nDim> &field, Box2D domain);

template <typename T1, typename T2, int nDim>
std::unique_ptr<MultiGridTensorField2D<T2, nDim> > copyConvert(
    MultiGridTensorField2D<T1, nDim> &field);

/* *************** Extract Sub-TensorField *************************** */

template <typename T, int nDim>
void extractSubDomain(
    MultiGridTensorField2D<T, nDim> &field, MultiGridTensorField2D<T, nDim> &extractedField,
    Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField2D<T, nDim> > extractSubDomain(
    MultiGridTensorField2D<T, nDim> &field, Box2D domain);

/* *************** Component (scalar-field) out of a tensor-field ****** */

template <typename T, int nDim>
void extractComponent(
    MultiGridTensorField2D<T, nDim> &tensorField, MultiGridScalarField2D<T> &component,
    Box2D domain, int iComponent);

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField2D<T> > extractComponent(
    MultiGridTensorField2D<T, nDim> &tensorField, Box2D domain, int iComponent);

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField2D<T> > extractComponent(
    MultiGridTensorField2D<T, nDim> &tensorField, int iComponent);

/* *************** Vector-norm of each cell in the field *************** */

template <typename T, int nDim>
void computeNorm(
    MultiGridTensorField2D<T, nDim> &tensorField, MultiGridScalarField2D<T> &norm, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField2D<T> > computeNorm(
    MultiGridTensorField2D<T, nDim> &tensorField, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField2D<T> > computeNorm(
    MultiGridTensorField2D<T, nDim> &tensorField);

/* *************** Squared vector-norm of each cell in the field ******** */

template <typename T, int nDim>
void computeNormSqr(
    MultiGridTensorField2D<T, nDim> &tensorField, MultiGridScalarField2D<T> &normSqr, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField2D<T> > computeNormSqr(
    MultiGridTensorField2D<T, nDim> &tensorField, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField2D<T> > computeNormSqr(
    MultiGridTensorField2D<T, nDim> &tensorField);

/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template <typename T>
void computeSymmetricTensorNorm(
    MultiGridTensorField2D<T, 3> &tensorField, MultiGridScalarField2D<T> &norm, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > computeSymmetricTensorNorm(
    MultiGridTensorField2D<T, 3> &tensorField, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > computeSymmetricTensorNorm(
    MultiGridTensorField2D<T, 3> &tensorField);

/* *************** Squared Tensor-norm of each symmetric tensor of a field*/

template <typename T>
void computeSymmetricTensorNormSqr(
    MultiGridTensorField2D<T, 3> &tensorField, MultiGridScalarField2D<T> &normSqr, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > computeSymmetricTensorNormSqr(
    MultiGridTensorField2D<T, 3> &tensorField, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > computeSymmetricTensorNormSqr(
    MultiGridTensorField2D<T, 3> &tensorField);

/* *************** Trace of each symmetric tensor of a field ************ */

template <typename T>
void computeSymmetricTensorTrace(
    MultiGridTensorField2D<T, 3> &tensorField, MultiGridScalarField2D<T> &trace, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > computeSymmetricTensorTrace(
    MultiGridTensorField2D<T, 3> &tensorField, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > computeSymmetricTensorTrace(
    MultiGridTensorField2D<T, 3> &tensorField);

/* *************** Vorticity from Velocity field *********************** */

template <typename T>
void computeVorticity(
    MultiGridTensorField2D<T, 2> &velocity, MultiGridScalarField2D<T> &vorticity, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > computeVorticity(
    MultiGridTensorField2D<T, 2> &velocity, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > computeVorticity(
    MultiGridTensorField2D<T, 2> &velocity);

/* *************** Vorticity, witout boundary treatment, from Velocity field  */

template <typename T>
void computeBulkVorticity(
    MultiGridTensorField2D<T, 2> &velocity, MultiGridScalarField2D<T> &vorticity, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > computeBulkVorticity(
    MultiGridTensorField2D<T, 2> &velocity, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > computeBulkVorticity(
    MultiGridTensorField2D<T, 2> &velocity);

/* *************** Strain rate from Velocity field ********************* */

template <typename T>
void computeStrainRate(
    MultiGridTensorField2D<T, 2> &velocity, MultiGridTensorField2D<T, 3> &S, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridTensorField2D<T, 3> > computeStrainRate(
    MultiGridTensorField2D<T, 2> &velocity, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridTensorField2D<T, 3> > computeStrainRate(
    MultiGridTensorField2D<T, 2> &velocity);

/* *************** Str. rate, witout boundary treatment, from Velocity field  */

template <typename T>
void computeBulkStrainRate(
    MultiGridTensorField2D<T, 2> &velocity, MultiGridTensorField2D<T, 3> &S, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridTensorField2D<T, 3> > computeBulkStrainRate(
    MultiGridTensorField2D<T, 2> &velocity, Box2D domain);

template <typename T>
std::unique_ptr<MultiGridTensorField2D<T, 3> > computeBulkStrainRate(
    MultiGridTensorField2D<T, 2> &velocity);

/* *************** MultiTensorField - MultiTensorField operations *************** */

template <typename T, int nDim>
void add(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B,
    MultiGridTensorField2D<T, nDim> &result, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField2D<T, nDim> > add(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField2D<T, nDim> > add(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B);

template <typename T, int nDim>
void subtract(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B,
    MultiGridTensorField2D<T, nDim> &result, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField2D<T, nDim> > subtract(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField2D<T, nDim> > subtract(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B);

template <typename T, int nDim>
void multiply(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B,
    MultiGridTensorField2D<T, nDim> &result, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField2D<T, nDim> > multiply(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField2D<T, nDim> > multiply(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B);

template <typename T, int nDim>
void divide(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B,
    MultiGridTensorField2D<T, nDim> &result, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField2D<T, nDim> > divide(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B, Box2D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField2D<T, nDim> > divide(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B);

/* *************** MultiTensorField - MultiTensorField inplace operations *************** */

template <typename T, int nDim>
void addInPlace(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B, Box2D domain);

template <typename T, int nDim>
void addInPlace(MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B);

template <typename T, int nDim>
void subtractInPlace(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B, Box2D domain);

template <typename T, int nDim>
void subtractInPlace(MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B);

template <typename T, int nDim>
void multiplyInPlace(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B, Box2D domain);

template <typename T, int nDim>
void multiplyInPlace(MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B);

template <typename T, int nDim>
void multiplyInPlace(MultiGridTensorField2D<T, nDim> &A, T alpha, Box2D domain);

template <typename T, int nDim>
void multiplyInPlace(MultiGridTensorField2D<T, nDim> &A, T alpha);

template <typename T, int nDim>
void divideInPlace(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B, Box2D domain);

template <typename T, int nDim>
void divideInPlace(MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B);

}  // namespace plb

#endif  // MULTI_GRID_DATA_ANALYSIS_WRAPPER_2D_H
