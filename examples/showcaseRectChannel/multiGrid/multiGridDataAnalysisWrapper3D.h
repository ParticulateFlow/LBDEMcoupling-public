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

#ifndef MULTI_GRID_DATA_ANALYSIS_WRAPPER_3D_H
#define MULTI_GRID_DATA_ANALYSIS_WRAPPER_3D_H

#include "dataProcessors/dataAnalysisWrapper3D.h"
#include "multiGrid/multiGridDataField3D.h"
#include "multiGrid/multiGridLattice3D.h"

namespace plb {

/* *************** Reductive functions ******************************* */

template <typename T, template <typename U> class Descriptor>
T computeAverageDensity(MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
T computeAverageDensity(MultiGridLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(MultiGridLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
T computeAverageEnergy(MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
T computeAverageEnergy(MultiGridLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor, class BoolMask>
plint count(MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain, BoolMask boolMask);

template <typename T, template <typename U> class Descriptor, class BoolMask>
plint count(MultiGridLattice3D<T, Descriptor> &lattice, BoolMask boolMask);

/* *************** Extract Sub-Lattice ******************************* */

template <typename T, template <typename U> class Descriptor>
void extractSubDomain(
    MultiGridLattice3D<T, Descriptor> &lattice, MultiGridLattice3D<T, Descriptor> &extractedLattice,
    Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridLattice3D<T, Descriptor> > extractSubDomain(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain);

/* *************** Density ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeDensity(
    MultiGridLattice3D<T, Descriptor> &lattice, MultiGridScalarField3D<T> &density, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeDensity(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeDensity(
    MultiGridLattice3D<T, Descriptor> &lattice);

/* *************** RhoBar ******************************************** */

template <typename T, template <typename U> class Descriptor>
void computeRhoBar(
    MultiGridLattice3D<T, Descriptor> &lattice, MultiGridScalarField3D<T> &rhoBar, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeRhoBar(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeRhoBar(
    MultiGridLattice3D<T, Descriptor> &lattice);

/* *************** Kinetic Energy ************************************ */

template <typename T, template <typename U> class Descriptor>
void computeKineticEnergy(
    MultiGridLattice3D<T, Descriptor> &lattice, MultiGridScalarField3D<T> &energy, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeKineticEnergy(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeKineticEnergy(
    MultiGridLattice3D<T, Descriptor> &lattice);

/* *************** Velocity Norm ************************************* */

template <typename T, template <typename U> class Descriptor>
void computeVelocityNorm(
    MultiGridLattice3D<T, Descriptor> &lattice, MultiGridScalarField3D<T> &velocityNorm,
    Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeVelocityNorm(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeVelocityNorm(
    MultiGridLattice3D<T, Descriptor> &lattice);

/* *************** Velocity Component ******************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocityComponent(
    MultiGridLattice3D<T, Descriptor> &lattice, MultiGridScalarField3D<T> &velocityComponent,
    Box3D domain, plint iComponent);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeVelocityComponent(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain, plint iComponent);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeVelocityComponent(
    MultiGridLattice3D<T, Descriptor> &lattice, plint iComponent);

/* *************** Velocity ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocity(
    MultiGridLattice3D<T, Descriptor> &lattice,
    MultiGridTensorField3D<T, Descriptor<T>::d> &velocity, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField3D<T, Descriptor<T>::d> > computeVelocity(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField3D<T, Descriptor<T>::d> > computeVelocity(
    MultiGridLattice3D<T, Descriptor> &lattice);

/* *************** Pi Neq ********************************* */

template <typename T, template <typename U> class Descriptor>
void computePiNeq(
    MultiGridLattice3D<T, Descriptor> &lattice,
    MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computePiNeq(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computePiNeq(
    MultiGridLattice3D<T, Descriptor> &lattice);

/* *************** Shear Stress ********************************* */

template <typename T, template <typename U> class Descriptor>
void computeShearStress(
    MultiGridLattice3D<T, Descriptor> &lattice,
    MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computeShearStress(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computeShearStress(
    MultiGridLattice3D<T, Descriptor> &lattice);

/* *************** Strain Rate from Stress *************************** */

template <typename T, template <typename U> class Descriptor>
void computeStrainRateFromStress(
    MultiGridLattice3D<T, Descriptor> &lattice,
    MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> &S, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> >
    computeStrainRateFromStress(MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> >
    computeStrainRateFromStress(MultiGridLattice3D<T, Descriptor> &lattice);

/* *************** Temperature ******************************************* */
// TODO

/* *************** SoundSpeed ******************************************* */
// TODO

/* *************** Population **************************************** */

template <typename T, template <typename U> class Descriptor>
void computePopulation(
    MultiGridLattice3D<T, Descriptor> &lattice, MultiGridScalarField3D<T> &population, Box3D domain,
    plint iPop);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computePopulation(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain, plint iPop);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computePopulation(
    MultiGridLattice3D<T, Descriptor> &lattice, plint iPop);

template <typename T, template <typename U> class Descriptor>
void computeAllPopulations(
    MultiGridLattice3D<T, Descriptor> &lattice,
    MultiGridTensorField3D<T, Descriptor<T>::q> &populations);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField3D<T, Descriptor<T>::q> > computeAllPopulations(
    MultiGridLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField3D<T, Descriptor<T>::q> > computeAllPopulations(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
void copyPopulations(
    MultiGridLattice3D<T, Descriptor> &latticeFrom, MultiGridLattice3D<T, Descriptor> &latticeTo,
    Box3D domain);

/* *************** Omega ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeOmega(
    MultiGridLattice3D<T, Descriptor> &lattice, MultiGridScalarField3D<T> &omega, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeOmega(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeOmega(
    MultiGridLattice3D<T, Descriptor> &lattice);

/* *************** Reductive functions ******************************* */

template <typename T>
T computeAverage(MultiGridScalarField3D<T> &scalarField, Box3D domain);

template <typename T>
T computeAverage(MultiGridScalarField3D<T> &scalarField);

template <typename T>
T computeMin(MultiGridScalarField3D<T> &scalarField, Box3D domain);

template <typename T>
T computeMin(MultiGridScalarField3D<T> &scalarField);

template <typename T>
T computeMax(MultiGridScalarField3D<T> &scalarField, Box3D domain);

template <typename T>
T computeMax(MultiGridScalarField3D<T> &scalarField);

template <typename T>
T computeBoundedAverage(MultiGridScalarField3D<T> &scalarField, Box3D domain);

template <typename T>
T computeBoundedAverage(MultiGridScalarField3D<T> &scalarField);

template <typename T, class BoolMask>
plint count(MultiGridScalarField3D<T> &field, Box3D domain, BoolMask boolMask);

template <typename T, class BoolMask>
plint count(MultiGridScalarField3D<T> &field, BoolMask boolMask);

/* *************** Generic operations *************** */

template <typename T, class Function>
void apply(Function f, MultiGridScalarField3D<T> &field, Box3D domain)
{
    applyProcessingFunctional(new ApplyScalarFunctional3D<T, Function>(f), domain, field);
}

template <typename T, class Function>
void apply(Function f, MultiGridScalarField3D<T> &field)
{
    apply(f, field, field.getBoundingBox());
}

template <typename T, class Function>
void evaluate(
    Function f, MultiGridScalarField3D<T> &field, MultiGridScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(
        new EvaluateScalarFunctional3D<T, Function>(f), domain, field, result);
}

template <typename T, class Function>
std::unique_ptr<MultiGridScalarField3D<T> > evaluate(
    Function f, MultiGridScalarField3D<T> &field, Box3D domain)
{
    MultiGridScalarField3D<T> *result = new MultiGridScalarField3D<T>(field, domain);
    evaluate(f, field, *result, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(result);
}

template <typename T, class Function>
std::unique_ptr<MultiGridScalarField3D<T> > evaluate(Function f, MultiGridScalarField3D<T> &field)
{
    return evaluate(f, field, field.getBoundingBox());
}

/* *************** Extract Sub-ScalarField *************************** */

template <typename T>
void extractSubDomain(
    MultiGridScalarField3D<T> &field, MultiGridScalarField3D<T> &extractedField, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > extractSubDomain(
    MultiGridScalarField3D<T> &field, Box3D domain);

/* *************** Copy and Convert ScalarField *************************** */

template <typename T1, typename T2>
void copyConvert(
    MultiGridScalarField3D<T1> &field, MultiGridScalarField3D<T2> &convertedField, Box3D domain);

template <typename T1, typename T2>
std::unique_ptr<MultiGridScalarField3D<T2> > copyConvert(
    MultiGridScalarField3D<T1> &field, Box3D domain);

template <typename T1, typename T2>
std::unique_ptr<MultiGridScalarField3D<T2> > copyConvert(MultiGridScalarField3D<T1> &field);

/* *************** MultiScalarField - Scalar operations *************** */

template <typename T>
void add(
    T scalar, MultiGridScalarField3D<T> &field, MultiGridScalarField3D<T> &result, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > add(
    T scalar, MultiGridScalarField3D<T> &field, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > add(T scalar, MultiGridScalarField3D<T> &field);

template <typename T>
void add(
    MultiGridScalarField3D<T> &field, T scalar, MultiGridScalarField3D<T> &result, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > add(
    MultiGridScalarField3D<T> &field, T scalar, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > add(MultiGridScalarField3D<T> &field, T scalar);

template <typename T>
void subtract(
    T scalar, MultiGridScalarField3D<T> &field, MultiGridScalarField3D<T> &result, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > subtract(
    T scalar, MultiGridScalarField3D<T> &field, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > subtract(T scalar, MultiGridScalarField3D<T> &field);

template <typename T>
void subtract(
    MultiGridScalarField3D<T> &field, T scalar, MultiGridScalarField3D<T> &result, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > subtract(
    MultiGridScalarField3D<T> &field, T scalar, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > subtract(MultiGridScalarField3D<T> &field, T scalar);

template <typename T>
void multiply(
    T scalar, MultiGridScalarField3D<T> &field, MultiGridScalarField3D<T> &result, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > multiply(
    T scalar, MultiGridScalarField3D<T> &field, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > multiply(T scalar, MultiGridScalarField3D<T> &field);

template <typename T>
void multiply(
    MultiGridScalarField3D<T> &field, T scalar, MultiGridScalarField3D<T> &result, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > multiply(
    MultiGridScalarField3D<T> &field, T scalar, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > multiply(MultiGridScalarField3D<T> &field, T scalar);

template <typename T>
void divide(
    T scalar, MultiGridScalarField3D<T> &field, MultiGridScalarField3D<T> &result, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > divide(
    T scalar, MultiGridScalarField3D<T> &field, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > divide(T scalar, MultiGridScalarField3D<T> &field);

template <typename T>
void divide(
    MultiGridScalarField3D<T> &field, T scalar, MultiGridScalarField3D<T> &result, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > divide(
    MultiGridScalarField3D<T> &field, T scalar, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > divide(MultiGridScalarField3D<T> &field, T scalar);

/* *************** MultiScalarField - Scalar inplace operations *************** */

template <typename T>
void addInPlace(MultiGridScalarField3D<T> &field, T scalar, Box3D domain);

template <typename T>
void addInPlace(MultiGridScalarField3D<T> &field, T scalar);

template <typename T>
void subtractInPlace(MultiGridScalarField3D<T> &field, T scalar, Box3D domain);

template <typename T>
void subtractInPlace(MultiGridScalarField3D<T> &field, T scalar);

template <typename T>
void multiplyInPlace(MultiGridScalarField3D<T> &field, T scalar, Box3D domain);

template <typename T>
void multiplyInPlace(MultiGridScalarField3D<T> &field, T scalar);

template <typename T>
void divideInPlace(MultiGridScalarField3D<T> &field, T scalar, Box3D domain);

template <typename T>
void divideInPlace(MultiGridScalarField3D<T> &field, T scalar);

/* *************** MultiScalarField - MultiScalarField operations *************** */

template <typename T>
void add(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, MultiGridScalarField3D<T> &result,
    Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > add(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > add(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B);

template <typename T>
void subtract(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, MultiGridScalarField3D<T> &result,
    Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > subtract(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > subtract(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B);

template <typename T>
void multiply(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, MultiGridScalarField3D<T> &result,
    Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > multiply(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > multiply(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B);

template <typename T>
void divide(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, MultiGridScalarField3D<T> &result,
    Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > divide(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > divide(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B);

/* *************** MultiScalarField - MultiScalarField inplace operations *************** */

template <typename T>
void addInPlace(MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, Box3D domain);

template <typename T>
void addInPlace(MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B);

template <typename T>
void subtractInPlace(MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, Box3D domain);

template <typename T>
void subtractInPlace(MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B);

template <typename T>
void multiplyInPlace(MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, Box3D domain);

template <typename T>
void multiplyInPlace(MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B);

template <typename T>
void divideInPlace(MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, Box3D domain);

template <typename T>
void divideInPlace(MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B);

/* ******************************************************************* */
/* ***************           Tensor-field                       ****** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template <typename T, int nDim, class BoolMask>
plint count(MultiGridTensorField3D<T, nDim> &field, Box3D domain, BoolMask boolMask);

template <typename T, int nDim, class BoolMask>
plint count(MultiGridTensorField3D<T, nDim> &field, BoolMask boolMask);

/* *************** Copy-convert a tensor-field *************** */

template <typename T1, typename T2, int nDim>
void copyConvert(
    MultiGridTensorField3D<T1, nDim> &field, MultiGridTensorField3D<T2, nDim> &convertedField,
    Box3D domain);

template <typename T1, typename T2, int nDim>
std::unique_ptr<MultiGridTensorField3D<T2, nDim> > copyConvert(
    MultiGridTensorField3D<T1, nDim> &field, Box3D domain);

template <typename T1, typename T2, int nDim>
std::unique_ptr<MultiGridTensorField3D<T2, nDim> > copyConvert(
    MultiGridTensorField3D<T1, nDim> &field);

/* *************** Extract Sub-TensorField *************************** */

template <typename T, int nDim>
void extractSubDomain(
    MultiGridTensorField3D<T, nDim> &field, MultiGridTensorField3D<T, nDim> &extractedField,
    Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField3D<T, nDim> > extractSubDomain(
    MultiGridTensorField3D<T, nDim> &field, Box3D domain);

/* *************** Component (scalar-field) out of a tensor-field ****** */

template <typename T, int nDim>
void extractComponent(
    MultiGridTensorField3D<T, nDim> &tensorField, MultiGridScalarField3D<T> &component,
    Box3D domain, int iComponent);

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField3D<T> > extractComponent(
    MultiGridTensorField3D<T, nDim> &tensorField, Box3D domain, int iComponent);

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField3D<T> > extractComponent(
    MultiGridTensorField3D<T, nDim> &tensorField, int iComponent);

/* *************** Vector-norm of each cell in the field *************** */

template <typename T, int nDim>
void computeNorm(
    MultiGridTensorField3D<T, nDim> &tensorField, MultiGridScalarField3D<T> &norm, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField3D<T> > computeNorm(
    MultiGridTensorField3D<T, nDim> &tensorField, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField3D<T> > computeNorm(
    MultiGridTensorField3D<T, nDim> &tensorField);

/* *************** Squared vector-norm of each cell in the field ******** */

template <typename T, int nDim>
void computeNormSqr(
    MultiGridTensorField3D<T, nDim> &tensorField, MultiGridScalarField3D<T> &normSqr, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField3D<T> > computeNormSqr(
    MultiGridTensorField3D<T, nDim> &tensorField, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField3D<T> > computeNormSqr(
    MultiGridTensorField3D<T, nDim> &tensorField);

/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template <typename T>
void computeSymmetricTensorNorm(
    MultiGridTensorField3D<T, 3> &tensorField, MultiGridScalarField3D<T> &norm, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > computeSymmetricTensorNorm(
    MultiGridTensorField3D<T, 3> &tensorField, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > computeSymmetricTensorNorm(
    MultiGridTensorField3D<T, 3> &tensorField);

/* *************** Squared Tensor-norm of each symmetric tensor of a field*/

template <typename T>
void computeSymmetricTensorNormSqr(
    MultiGridTensorField3D<T, 3> &tensorField, MultiGridScalarField3D<T> &normSqr, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > computeSymmetricTensorNormSqr(
    MultiGridTensorField3D<T, 3> &tensorField, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > computeSymmetricTensorNormSqr(
    MultiGridTensorField3D<T, 3> &tensorField);

/* *************** Trace of each symmetric tensor of a field ************ */

template <typename T>
void computeSymmetricTensorTrace(
    MultiGridTensorField3D<T, 3> &tensorField, MultiGridScalarField3D<T> &trace, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > computeSymmetricTensorTrace(
    MultiGridTensorField3D<T, 3> &tensorField, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > computeSymmetricTensorTrace(
    MultiGridTensorField3D<T, 3> &tensorField);

/* *************** Vorticity from Velocity field *********************** */

template <typename T>
void computeVorticity(
    MultiGridTensorField3D<T, 3> &velocity, MultiGridTensorField3D<T, 3> &vorticity, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridTensorField3D<T, 3> > computeVorticity(
    MultiGridTensorField3D<T, 3> &velocity, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridTensorField3D<T, 3> > computeVorticity(
    MultiGridTensorField3D<T, 3> &velocity);

/* *************** Vorticity, witout boundary treatment, from Velocity field  */

template <typename T>
void computeBulkVorticity(
    MultiGridTensorField3D<T, 3> &velocity, MultiGridTensorField3D<T, 3> &vorticity, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridTensorField3D<T, 3> > computeBulkVorticity(
    MultiGridTensorField3D<T, 3> &velocity, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridTensorField3D<T, 3> > computeBulkVorticity(
    MultiGridTensorField3D<T, 3> &velocity);

/* *************** Strain rate from Velocity field ********************* */

template <typename T>
void computeStrainRate(
    MultiGridTensorField3D<T, 3> &velocity, MultiGridTensorField3D<T, 6> &S, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridTensorField3D<T, 6> > computeStrainRate(
    MultiGridTensorField3D<T, 3> &velocity, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridTensorField3D<T, 6> > computeStrainRate(
    MultiGridTensorField3D<T, 3> &velocity);

/* *************** Str. rate, witout boundary treatment, from Velocity field  */

template <typename T>
void computeBulkStrainRate(
    MultiGridTensorField3D<T, 3> &velocity, MultiGridTensorField3D<T, 6> &S, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridTensorField3D<T, 6> > computeBulkStrainRate(
    MultiGridTensorField3D<T, 3> &velocity, Box3D domain);

template <typename T>
std::unique_ptr<MultiGridTensorField3D<T, 6> > computeBulkStrainRate(
    MultiGridTensorField3D<T, 3> &velocity);

/* *************** MultiTensorField - MultiTensorField operations *************** */

template <typename T, int nDim>
void add(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B,
    MultiGridTensorField3D<T, nDim> &result, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField3D<T, nDim> > add(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField3D<T, nDim> > add(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B);

template <typename T, int nDim>
void subtract(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B,
    MultiGridTensorField3D<T, nDim> &result, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField3D<T, nDim> > subtract(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField3D<T, nDim> > subtract(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B);

template <typename T, int nDim>
void multiply(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B,
    MultiGridTensorField3D<T, nDim> &result, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField3D<T, nDim> > multiply(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField3D<T, nDim> > multiply(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B);

template <typename T, int nDim>
void divide(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B,
    MultiGridTensorField3D<T, nDim> &result, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField3D<T, nDim> > divide(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField3D<T, nDim> > divide(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B);

/* *************** MultiTensorField - MultiTensorField inplace operations *************** */

template <typename T, int nDim>
void addInPlace(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B, Box3D domain);

template <typename T, int nDim>
void addInPlace(MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B);

template <typename T, int nDim>
void subtractInPlace(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B, Box3D domain);

template <typename T, int nDim>
void subtractInPlace(MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B);

template <typename T, int nDim>
void multiplyInPlace(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B, Box3D domain);

template <typename T, int nDim>
void multiplyInPlace(MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B);

template <typename T, int nDim>
void multiplyInPlace(MultiGridTensorField3D<T, nDim> &A, T alpha, Box3D domain);

template <typename T, int nDim>
void multiplyInPlace(MultiGridTensorField3D<T, nDim> &A, T alpha);

template <typename T, int nDim>
void divideInPlace(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B, Box3D domain);

template <typename T, int nDim>
void divideInPlace(MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B);

}  // namespace plb

#endif  // MULTI_GRID_DATA_ANALYSIS_WRAPPER_3D_H
