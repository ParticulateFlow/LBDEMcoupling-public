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

#ifndef MULTI_GRID_DATA_ANALYSIS_WRAPPER_3D_HH
#define MULTI_GRID_DATA_ANALYSIS_WRAPPER_3D_HH

#include "dataProcessors/dataAnalysisWrapper3D.h"
#include "multiGrid/multiGridDataAnalysisWrapper3D.h"
#include "multiGrid/multiGridDataProcessorWrapper3D.h"

namespace plb {

/* *************** Reductive functions ******************************* */

template <typename T>
T computeAverage(MultiGridScalarField3D<T> &scalarField, Box3D domain)
{
    BoxScalarSumFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField, scalarField.getReferenceLevel());
    return functional.getSumScalar() / (T)domain.nCells();
}

template <typename T>
T computeAverage(MultiGridScalarField3D<T> &scalarField)
{
    return computeAverage(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeMin(MultiGridScalarField3D<T> &scalarField, Box3D domain)
{
    BoxScalarMinFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField, scalarField.getReferenceLevel());
    return functional.getMinScalar();
}

template <typename T>
T computeMin(MultiGridScalarField3D<T> &scalarField)
{
    return computeMin(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeMax(MultiGridScalarField3D<T> &scalarField, Box3D domain)
{
    BoxScalarMaxFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField, scalarField.getReferenceLevel());
    return functional.getMaxScalar();
}

template <typename T>
T computeMax(MultiGridScalarField3D<T> &scalarField)
{
    return computeMax(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeBoundedAverage(MultiGridScalarField3D<T> &scalarField, Box3D domain)
{
    BoundedBoxScalarSumFunctional3D<T> functional;
    plint envelopeWidth = 1;
    applyProcessingFunctional(
        functional, domain, scalarField, scalarField.getReferenceLevel(), envelopeWidth);
    return functional.getSumScalar()
           / (T)((domain.getNx() - 1) * (domain.getNy() - 1) * (domain.getNz() - 1));
}

template <typename T>
T computeBoundedAverage(MultiGridScalarField3D<T> &scalarField)
{
    return computeBoundedAverage(scalarField, scalarField.getBoundingBox());
}

template <typename T, class BoolMask>
plint count(MultiGridScalarField3D<T> &field, Box3D domain, BoolMask boolMask)
{
    CountScalarElementsFunctional3D<T, BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, field, field.getReferenceLevel());
    return functional.getCount();
}

template <typename T, class BoolMask>
plint count(MultiGridScalarField3D<T> &field, BoolMask boolMask)
{
    return count(field, field.getBoundingBox(), boolMask);
}

/* *************** Extract Sub-Lattice ******************************* */

template <typename T, template <typename U> class Descriptor>
void extractSubDomain(
    MultiGridLattice3D<T, Descriptor> &lattice, MultiGridLattice3D<T, Descriptor> &extractedLattice,
    Box3D domain)
{
    applyProcessingFunctional(
        new LatticeRegenerateFunctional3D<T, Descriptor>, domain, lattice, extractedLattice,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridLattice3D<T, Descriptor> > extractSubDomain(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    MultiGridLattice3D<T, Descriptor> *extractedLattice =
        new MultiGridLattice3D<T, Descriptor>(lattice, domain);
    extractSubDomain(lattice, *extractedLattice, domain);
    return std::unique_ptr<MultiGridLattice3D<T, Descriptor> >(extractedLattice);
}

/* *************** Density ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeDensity(
    MultiGridLattice3D<T, Descriptor> &lattice, MultiGridScalarField3D<T> &density, Box3D domain)
{
    applyProcessingFunctional(
        new BoxDensityFunctional3D<T, Descriptor>, domain, lattice, density,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeDensity(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    MultiGridScalarField3D<T> *density = new MultiGridScalarField3D<T>(lattice, domain);
    computeDensity(lattice, *density, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(density);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeDensity(
    MultiGridLattice3D<T, Descriptor> &lattice)
{
    return computeDensity(lattice, lattice.getBoundingBox());
}

/* *************** RhoBar ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeRhoBar(
    MultiGridLattice3D<T, Descriptor> &lattice, MultiGridScalarField3D<T> &rhoBar, Box3D domain)
{
    applyProcessingFunctional(
        new BoxRhoBarFunctional3D<T, Descriptor>, domain, lattice, rhoBar,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeRhoBar(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    MultiGridScalarField3D<T> *rhoBar = new MultiGridScalarField3D<T>(lattice, domain);
    computeRhoBar(lattice, *rhoBar, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(rhoBar);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeRhoBar(
    MultiGridLattice3D<T, Descriptor> &lattice)
{
    return computeRhoBar(lattice, lattice.getBoundingBox());
}

/* *************** Kinetic Energy ************************************ */

template <typename T, template <typename U> class Descriptor>
void computeKineticEnergy(
    MultiGridLattice3D<T, Descriptor> &lattice, MultiGridScalarField3D<T> &energy, Box3D domain)
{
    applyProcessingFunctional(
        new BoxKineticEnergyFunctional3D<T, Descriptor>, domain, lattice, energy,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeKineticEnergy(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    MultiGridScalarField3D<T> *energy = new MultiGridScalarField3D<T>(lattice, domain);
    computeKineticEnergy(lattice, *energy, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(energy);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeKineticEnergy(
    MultiGridLattice3D<T, Descriptor> &lattice)
{
    return computeKineticEnergy(lattice, lattice.getBoundingBox());
}

/* *************** Velocity Norm ************************************* */

template <typename T, template <typename U> class Descriptor>
void computeVelocityNorm(
    MultiGridLattice3D<T, Descriptor> &lattice, MultiGridScalarField3D<T> &velocityNorm,
    Box3D domain)
{
    applyProcessingFunctional(
        new BoxVelocityNormFunctional3D<T, Descriptor>, domain, lattice, velocityNorm,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeVelocityNorm(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    MultiGridScalarField3D<T> *velocityNorm = new MultiGridScalarField3D<T>(lattice, domain);
    computeVelocityNorm(lattice, *velocityNorm, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(velocityNorm);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeVelocityNorm(
    MultiGridLattice3D<T, Descriptor> &lattice)
{
    return computeVelocityNorm(lattice, lattice.getBoundingBox());
}

/* *************** Velocity Component ******************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocityComponent(
    MultiGridLattice3D<T, Descriptor> &lattice, MultiGridScalarField3D<T> &velocityComponent,
    Box3D domain, plint iComponent)
{
    applyProcessingFunctional(
        new BoxVelocityComponentFunctional3D<T, Descriptor>(iComponent), domain, lattice,
        velocityComponent, lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeVelocityComponent(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain, plint iComponent)
{
    MultiGridScalarField3D<T> *velocityComponent = new MultiGridScalarField3D<T>(lattice, domain);
    computeVelocityComponent(lattice, *velocityComponent, domain, iComponent);
    return std::unique_ptr<MultiGridScalarField3D<T> >(velocityComponent);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeVelocityComponent(
    MultiGridLattice3D<T, Descriptor> &lattice, plint iComponent)
{
    return computeVelocityComponent(lattice, lattice.getBoundingBox(), iComponent);
}

/* *************** Velocity ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocity(
    MultiGridLattice3D<T, Descriptor> &lattice,
    MultiGridTensorField3D<T, Descriptor<T>::d> &velocity, Box3D domain)
{
    applyProcessingFunctional(
        new BoxVelocityFunctional3D<T, Descriptor>, domain, lattice, velocity,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField3D<T, Descriptor<T>::d> > computeVelocity(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    MultiGridTensorField3D<T, Descriptor<T>::d> *velocity =
        new MultiGridTensorField3D<T, Descriptor<T>::d>(lattice, domain);
    computeVelocity(lattice, *velocity, domain);
    return std::unique_ptr<MultiGridTensorField3D<T, Descriptor<T>::d> >(velocity);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField3D<T, Descriptor<T>::d> > computeVelocity(
    MultiGridLattice3D<T, Descriptor> &lattice)
{
    return computeVelocity(lattice, lattice.getBoundingBox());
}

/* *************** PiNeq ********************************* */

template <typename T, template <typename U> class Descriptor>
void computePiNeq(
    MultiGridLattice3D<T, Descriptor> &lattice,
    MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq, Box3D domain)
{
    applyProcessingFunctional(
        new BoxPiNeqFunctional3D<T, Descriptor>, domain, lattice, PiNeq,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computePiNeq(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> *PiNeq =
        new MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n>(lattice, domain);
    computePiNeq(lattice, *PiNeq, domain);
    return std::unique_ptr<MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> >(PiNeq);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computePiNeq(
    MultiGridLattice3D<T, Descriptor> &lattice)
{
    return computePiNeq(lattice, lattice.getBoundingBox());
}

/* *************** ShearStress ********************************* */

template <typename T, template <typename U> class Descriptor>
void computeShearStress(
    MultiGridLattice3D<T, Descriptor> &lattice,
    MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> &stress, Box3D domain)
{
    applyProcessingFunctional(
        new BoxShearStressFunctional3D<T, Descriptor>, domain, lattice, stress,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computeShearStress(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> *stress =
        new MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n>(lattice, domain);
    computeShearStress(lattice, *stress, domain);
    return std::unique_ptr<MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> >(stress);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computeShearStress(
    MultiGridLattice3D<T, Descriptor> &lattice)
{
    return computeShearStress(lattice, lattice.getBoundingBox());
}

/* *************** Strain Rate from Stress *************************** */

template <typename T, template <typename U> class Descriptor>
void computeStrainRateFromStress(
    MultiGridLattice3D<T, Descriptor> &lattice,
    MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> &S, Box3D domain)
{
    applyProcessingFunctional(
        new BoxStrainRateFromStressFunctional3D<T, Descriptor>, domain, lattice, S,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> >
    computeStrainRateFromStress(MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> *S =
        new MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n>(lattice, domain);
    computeStrainRateFromStress(lattice, *S, domain);
    return std::unique_ptr<MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> >(S);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField3D<T, SymmetricTensor<T, Descriptor>::n> >
    computeStrainRateFromStress(MultiGridLattice3D<T, Descriptor> &lattice)
{
    return computeStrainRateFromStress(lattice, lattice.getBoundingBox());
}

/* *************** Temperature ******************************************* */
// TODO

/* *************** SoundSpeed ******************************************* */
// TODO

/* *************** Population **************************************** */

template <typename T, template <typename U> class Descriptor>
void computePopulation(
    MultiGridLattice3D<T, Descriptor> &lattice, MultiGridScalarField3D<T> &population, Box3D domain,
    plint iPop)
{
    applyProcessingFunctional(
        new BoxPopulationFunctional3D<T, Descriptor>(iPop), domain, lattice, population,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computePopulation(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain, plint iPop)
{
    MultiGridScalarField3D<T> *population = new MultiGridScalarField3D<T>(lattice, domain);
    computePopulation(lattice, *population, domain, iPop);
    return std::unique_ptr<MultiGridScalarField3D<T> >(population);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computePopulation(
    MultiGridLattice3D<T, Descriptor> &lattice, plint iPop)
{
    return computePopulation(lattice, lattice.getBoundingBox(), iPop);
}

template <typename T, template <typename U> class Descriptor>
void computeAllPopulations(
    MultiGridLattice3D<T, Descriptor> &lattice,
    MultiGridTensorField3D<T, Descriptor<T>::q> &populations, Box3D domain)
{
    applyProcessingFunctional(
        new BoxAllPopulationsFunctional3D<T, Descriptor>(), domain, lattice, populations,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeAllPopulations(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    MultiGridTensorField3D<T, Descriptor<T>::q> *populations =
        new MultiGridTensorField3D<T, Descriptor<T>::q>(lattice, domain);
    computeAllPopulations(lattice, *populations, domain);
    return std::unique_ptr<MultiGridTensorField3D<T, Descriptor<T>::q> >(populations);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField3D<T, Descriptor<T>::q> > computeAllPopulations(
    MultiGridLattice3D<T, Descriptor> &lattice)
{
    return computeAllPopulations(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void copyPopulations(
    MultiGridLattice3D<T, Descriptor> &latticeFrom, MultiGridLattice3D<T, Descriptor> &latticeTo,
    Box3D domain)
{
    applyProcessingFunctional(
        new CopyPopulationsFunctional3D<T, Descriptor>(), domain, latticeFrom, latticeTo);
}

/* *************** Omega ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeOmega(
    MultiGridLattice3D<T, Descriptor> &lattice, MultiGridScalarField3D<T> &omega, Box3D domain)
{
    applyProcessingFunctional(
        new BoxOmegaFunctional3D<T, Descriptor>, domain, lattice, omega,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeOmega(
    MultiGridLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    MultiGridScalarField3D<T> *omega = new MultiGridScalarField3D<T>(lattice, domain);
    computeOmega(lattice, *omega, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(omega);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField3D<T> > computeOmega(MultiGridLattice3D<T, Descriptor> &lattice)
{
    return computeOmega(lattice, lattice.getBoundingBox());
}

/* *************** Extract Sub-ScalarField *************************** */

template <typename T>
void extractSubDomain(
    MultiGridScalarField3D<T> &field, MultiGridScalarField3D<T> &extractedField, Box3D domain)
{
    applyProcessingFunctional(
        new ExtractScalarSubDomainFunctional3D<T>, domain, field, extractedField,
        field.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > extractSubDomain(
    MultiGridScalarField3D<T> &field, Box3D domain)
{
    MultiGridScalarField3D<T> *extractedField = new MultiGridScalarField3D<T>(field, domain);
    extractSubDomain(field, *extractedField, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(extractedField);
}

/* *************** MultiScalarField - Scalar operations *************** */

template <typename T>
void add(
    MultiGridScalarField3D<T> &field, T scalar, MultiGridScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(
        new A_plus_alpha_functional3D<T>(scalar), domain, field, result, field.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > add(
    MultiGridScalarField3D<T> &field, T scalar, Box3D domain)
{
    MultiGridScalarField3D<T> *result = new MultiGridScalarField3D<T>(field, domain);
    add(field, scalar, *result, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > add(MultiGridScalarField3D<T> &field, T scalar)
{
    return add(field, scalar, field.getBoundingBox());
}

template <typename T>
void add(
    T scalar, MultiGridScalarField3D<T> &field, MultiGridScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(
        new A_plus_alpha_functional3D<T>(scalar), domain, field, result, field.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > add(
    T scalar, MultiGridScalarField3D<T> &field, Box3D domain)
{
    MultiGridScalarField3D<T> *result = new MultiGridScalarField3D<T>(field, domain);
    add(scalar, field, *result, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > add(T scalar, MultiGridScalarField3D<T> &field)
{
    return add(scalar, field, field.getBoundingBox());
}

template <typename T>
void subtract(
    MultiGridScalarField3D<T> &field, T scalar, MultiGridScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(
        new A_minus_alpha_functional3D<T>(scalar), domain, field, result,
        field.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > subtract(
    MultiGridScalarField3D<T> &field, T scalar, Box3D domain)
{
    MultiGridScalarField3D<T> *result = new MultiGridScalarField3D<T>(field, domain);
    subtract(field, scalar, *result, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > subtract(MultiGridScalarField3D<T> &field, T scalar)
{
    return subtract(field, scalar, field.getBoundingBox());
}

template <typename T>
void subtract(
    T scalar, MultiGridScalarField3D<T> &field, MultiGridScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(
        new Alpha_minus_A_functional3D<T>(scalar), domain, field, result,
        field.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > subtract(
    T scalar, MultiGridScalarField3D<T> &field, Box3D domain)
{
    MultiGridScalarField3D<T> *result = new MultiGridScalarField3D<T>(field, domain);
    subtract(scalar, field, *result, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > subtract(T scalar, MultiGridScalarField3D<T> &field)
{
    return subtract(scalar, field, field.getBoundingBox());
}

template <typename T>
void multiply(
    MultiGridScalarField3D<T> &field, T scalar, MultiGridScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(
        new A_times_alpha_functional3D<T>(scalar), domain, field, result,
        field.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > multiply(
    MultiGridScalarField3D<T> &field, T scalar, Box3D domain)
{
    MultiGridScalarField3D<T> *result = new MultiGridScalarField3D<T>(field, domain);
    multiply(field, scalar, *result, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > multiply(MultiGridScalarField3D<T> &field, T scalar)
{
    return multiply(field, scalar, field.getBoundingBox());
}

template <typename T>
void multiply(
    T scalar, MultiGridScalarField3D<T> &field, MultiGridScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(
        new A_times_alpha_functional3D<T>(scalar), domain, field, result,
        field.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > multiply(
    T scalar, MultiGridScalarField3D<T> &field, Box3D domain)
{
    MultiGridScalarField3D<T> *result = new MultiGridScalarField3D<T>(field, domain);
    multiply(scalar, field, *result, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > multiply(T scalar, MultiGridScalarField3D<T> &field)
{
    return multiply(scalar, field, field.getBoundingBox());
}

template <typename T>
void divide(
    MultiGridScalarField3D<T> &field, T scalar, MultiGridScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(
        new A_dividedBy_alpha_functional3D<T>(scalar), domain, field, result,
        field.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > divide(
    MultiGridScalarField3D<T> &field, T scalar, Box3D domain)
{
    MultiGridScalarField3D<T> *result = new MultiGridScalarField3D<T>(field, domain);
    divide(field, scalar, *result, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > divide(MultiGridScalarField3D<T> &field, T scalar)
{
    return divide(field, scalar, field.getBoundingBox());
}

template <typename T>
void divide(
    T scalar, MultiGridScalarField3D<T> &field, MultiGridScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(
        new Alpha_dividedBy_A_functional3D<T>(scalar), domain, field, result,
        field.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > divide(
    T scalar, MultiGridScalarField3D<T> &field, Box3D domain)
{
    MultiGridScalarField3D<T> *result = new MultiGridScalarField3D<T>(field, domain);
    divide(scalar, field, *result, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > divide(T scalar, MultiGridScalarField3D<T> &field)
{
    return divide(scalar, field, field.getBoundingBox());
}

/* *************** MultiScalarField - Scalar inplace operations *************** */

template <typename T>
void addInPlace(MultiGridScalarField3D<T> &field, T scalar, Box3D domain)
{
    applyProcessingFunctional(
        new A_plus_alpha_inplace_functional3D<T>(scalar), domain, field, field.getReferenceLevel());
}

template <typename T>
void addInPlace(MultiGridScalarField3D<T> &field, T scalar)
{
    addInPlace(field, scalar, field.getBoundingBox());
}

template <typename T>
void subtractInPlace(MultiGridScalarField3D<T> &field, T scalar, Box3D domain)
{
    applyProcessingFunctional(
        new A_minus_alpha_inplace_functional3D<T>(scalar), domain, field,
        field.getReferenceLevel());
}

template <typename T>
void subtractInPlace(MultiGridScalarField3D<T> &field, T scalar)
{
    subtractInPlace(field, scalar, field.getBoundingBox());
}

template <typename T>
void multiplyInPlace(MultiGridScalarField3D<T> &field, T scalar, Box3D domain)
{
    applyProcessingFunctional(
        new A_times_alpha_inplace_functional3D<T>(scalar), domain, field,
        field.getReferenceLevel());
}

template <typename T>
void multiplyInPlace(MultiGridScalarField3D<T> &field, T scalar)
{
    multiplyInPlace(field, scalar, field.getBoundingBox());
}

template <typename T>
void divideInPlace(MultiGridScalarField3D<T> &field, T scalar, Box3D domain)
{
    applyProcessingFunctional(
        new A_dividedBy_alpha_inplace_functional3D<T>(scalar), domain, field,
        field.getReferenceLevel());
}

template <typename T>
void divideInPlace(MultiGridScalarField3D<T> &field, T scalar)
{
    divideInPlace(field, scalar, field.getBoundingBox());
}

/* *************** MultiScalarField - MultiScalarField operations *************** */

template <typename T1, typename T2>
void copyConvert(
    MultiGridScalarField3D<T1> &field, MultiGridScalarField3D<T2> &convertedField, Box3D domain)
{
    applyProcessingFunctional(
        new CopyConvertScalarFunctional3D<T1, T2>, domain, field, convertedField,
        field.getReferenceLevel());
}

template <typename T1, typename T2>
std::unique_ptr<MultiGridScalarField3D<T2> > copyConvert(
    MultiGridScalarField3D<T1> &field, Box3D domain)
{
    MultiGridScalarField3D<T2> *convertedField = new MultiGridScalarField3D<T2>(field, domain);
    copyConvert(field, *convertedField, domain);
    return std::unique_ptr<MultiGridScalarField3D<T2> >(convertedField);
}

template <typename T1, typename T2>
std::unique_ptr<MultiGridScalarField3D<T2> > copyConvert(MultiGridScalarField3D<T1> &field)
{
    return copyConvert<T1, T2>(field, field.getBoundingBox());
}

template <typename T>
void add(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, MultiGridScalarField3D<T> &result,
    Box3D domain)
{
    std::vector<MultiGridScalarField3D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_plus_B_functional3D<T>, domain, fields, A.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > add(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, Box3D domain)
{
    MultiGridScalarField3D<T> *result = new MultiGridScalarField3D<T>(A, domain);
    add(A, B, *result, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > add(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B)
{
    return add(A, B, A.getBoundingBox());
}

template <typename T>
void subtract(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, MultiGridScalarField3D<T> &result,
    Box3D domain)
{
    std::vector<MultiGridScalarField3D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_minus_B_functional3D<T>, domain, fields, A.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > subtract(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, Box3D domain)
{
    MultiGridScalarField3D<T> *result = new MultiGridScalarField3D<T>(A, domain);
    subtract(A, B, *result, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > subtract(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B)
{
    return subtract(A, B, A.getBoundingBox());
}

template <typename T>
void multiply(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, MultiGridScalarField3D<T> &result,
    Box3D domain)
{
    std::vector<MultiGridScalarField3D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_times_B_functional3D<T>, domain, fields, A.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > multiply(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, Box3D domain)
{
    MultiGridScalarField3D<T> *result = new MultiGridScalarField3D<T>(A, domain);
    multiply(A, B, *result, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > multiply(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B)
{
    return multiply(A, B, A.getBoundingBox());
}

template <typename T>
void divide(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, MultiGridScalarField3D<T> &result,
    Box3D domain)
{
    std::vector<MultiGridScalarField3D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(
        new A_dividedBy_B_functional3D<T>, domain, fields, A.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > divide(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, Box3D domain)
{
    MultiGridScalarField3D<T> *result = new MultiGridScalarField3D<T>(A, domain);
    divide(A, B, *result, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > divide(
    MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B)
{
    return divide(A, B, A.getBoundingBox());
}

/* *************** ScalarField - ScalarField inplace operations *************** */

template <typename T>
void addInPlace(MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, Box3D domain)
{
    applyProcessingFunctional(
        new A_plus_B_inplace_functional3D<T>, domain, A, B, A.getReferenceLevel());
}

template <typename T>
void addInPlace(MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B)
{
    addInPlace(A, B, A.getBoundingBox());
}

template <typename T>
void subtractInPlace(MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, Box3D domain)
{
    applyProcessingFunctional(
        new A_minus_B_inplace_functional3D<T>, domain, A, B, A.getReferenceLevel());
}

template <typename T>
void subtractInPlace(MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B)
{
    subtractInPlace(A, B, A.getBoundingBox());
}

template <typename T>
void multiplyInPlace(MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, Box3D domain)
{
    applyProcessingFunctional(
        new A_times_B_inplace_functional3D<T>, domain, A, B, A.getReferenceLevel());
}

template <typename T>
void multiplyInPlace(MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B)
{
    multiplyInPlace(A, B, A.getBoundingBox());
}

template <typename T>
void divideInPlace(MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B, Box3D domain)
{
    applyProcessingFunctional(
        new A_dividedBy_B_inplace_functional3D<T>, domain, A, B, A.getReferenceLevel());
}

template <typename T>
void divideInPlace(MultiGridScalarField3D<T> &A, MultiGridScalarField3D<T> &B)
{
    divideInPlace(A, B, A.getBoundingBox());
}

/* ******************************************************************* */
/* *************** PART VI : Multi-block wrappers: Tensor-field ****** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template <typename T, plint nDim, class BoolMask>
plint count(MultiGridTensorField3D<T, nDim> &field, Box3D domain, BoolMask boolMask)
{
    CountTensorElementsFunctional3D<T, nDim, BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, field, field.getReferenceLevel());
    return functional.getCount();
}

template <typename T, plint nDim, class BoolMask>
plint count(MultiGridTensorField3D<T, nDim> &field, BoolMask boolMask)
{
    return count(field, field.getBoundingBox(), boolMask);
}

/* *************** Extract Sub-TensorField *************************** */

template <typename T, int nDim>
void extractSubDomain(
    MultiGridTensorField3D<T, nDim> &field, MultiGridTensorField3D<T, nDim> &extractedField,
    Box3D domain)
{
    applyProcessingFunctional(
        new ExtractTensorSubDomainFunctional3D<T, nDim>, domain, field, extractedField,
        field.getReferenceLevel());
}

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField3D<T, nDim> > extractSubDomain(
    MultiGridTensorField3D<T, nDim> &field, Box3D domain)
{
    MultiGridTensorField3D<T, nDim> *extractedField =
        new MultiGridTensorField3D<T, nDim>(field, domain);
    extractSubDomain(field, *extractedField, domain);
    return std::unique_ptr<MultiGridTensorField3D<T, nDim> >(extractedField);
}

/* *************** Component (scalar-field) out of a tensor-field ****** */

template <typename T, int nDim>
void extractComponent(
    MultiGridTensorField3D<T, nDim> &tensorField, MultiGridScalarField3D<T> &component,
    Box3D domain, int iComponent)
{
    applyProcessingFunctional(
        new ExtractTensorComponentFunctional3D<T, nDim>(iComponent), domain, component, tensorField,
        tensorField.getReferenceLevel());
}

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField3D<T> > extractComponent(
    MultiGridTensorField3D<T, nDim> &tensorField, Box3D domain, int iComponent)
{
    MultiGridScalarField3D<T> *component = new MultiGridScalarField3D<T>(tensorField, domain);
    extractComponent(tensorField, *component, domain, iComponent);
    return std::unique_ptr<MultiGridScalarField3D<T> >(component);
}

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField3D<T> > extractComponent(
    MultiGridTensorField3D<T, nDim> &tensorField, int iComponent)
{
    return extractComponent(tensorField, tensorField.getBoundingBox(), iComponent);
}

/* *************** Vector-norm of each cell in the field *************** */

template <typename T, int nDim>
void computeNorm(
    MultiGridTensorField3D<T, nDim> &tensorField, MultiGridScalarField3D<T> &component,
    Box3D domain)
{
    applyProcessingFunctional(
        new ComputeNormFunctional3D<T, nDim>, domain, component, tensorField,
        tensorField.getReferenceLevel());
}

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField3D<T> > computeNorm(
    MultiGridTensorField3D<T, nDim> &tensorField, Box3D domain)
{
    MultiGridScalarField3D<T> *component = new MultiGridScalarField3D<T>(tensorField, domain);
    computeNorm(tensorField, *component, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(component);
}

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField3D<T> > computeNorm(
    MultiGridTensorField3D<T, nDim> &tensorField)
{
    return computeNorm(tensorField, tensorField.getBoundingBox());
}

/* *************** Squared vector-norm of each cell in the field ******** */

template <typename T, int nDim>
void computeNormSqr(
    MultiGridTensorField3D<T, nDim> &tensorField, MultiGridScalarField3D<T> &component,
    Box3D domain)
{
    applyProcessingFunctional(
        new ComputeNormSqrFunctional3D<T, nDim>, domain, component, tensorField,
        tensorField.getReferenceLevel());
}

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField3D<T> > computeNormSqr(
    MultiGridTensorField3D<T, nDim> &tensorField, Box3D domain)
{
    MultiGridScalarField3D<T> *component = new MultiGridScalarField3D<T>(tensorField, domain);
    computeNormSqr(tensorField, *component, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(component);
}

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField3D<T> > computeNormSqr(
    MultiGridTensorField3D<T, nDim> &tensorField)
{
    return computeNormSqr(tensorField, tensorField.getBoundingBox());
}

/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template <typename T>
void computeSymmetricTensorNorm(
    MultiGridTensorField3D<T, 3> &tensorField, MultiGridScalarField3D<T> &norm, Box3D domain)
{
    applyProcessingFunctional(
        new ComputeSymmetricTensorNormFunctional3D<T>, domain, norm, tensorField,
        tensorField.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > computeSymmetricTensorNorm(
    MultiGridTensorField3D<T, 3> &tensorField, Box3D domain)
{
    MultiGridScalarField3D<T> *norm = new MultiGridScalarField3D<T>(tensorField, domain);
    computeSymmetricTensorNorm(tensorField, *norm, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(norm);
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > computeSymmetricTensorNorm(
    MultiGridTensorField3D<T, 3> &tensorField)
{
    return computeSymmetricTensorNorm(tensorField, tensorField.getBoundingBox());
}

/* *************** Squared Tensor-norm of each symmetric tensor of a field*/

template <typename T>
void computeSymmetricTensorNormSqr(
    MultiGridTensorField3D<T, 3> &tensorField, MultiGridScalarField3D<T> &normSqr, Box3D domain)
{
    applyProcessingFunctional(
        new ComputeSymmetricTensorNormSqrFunctional3D<T>, domain, normSqr, tensorField,
        tensorField.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > computeSymmetricTensorNormSqr(
    MultiGridTensorField3D<T, 3> &tensorField, Box3D domain)
{
    MultiGridScalarField3D<T> *normSqr = new MultiGridScalarField3D<T>(tensorField, domain);
    computeSymmetricTensorNormSqr(tensorField, *normSqr, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(normSqr);
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > computeSymmetricTensorNormSqr(
    MultiGridTensorField3D<T, 3> &tensorField)
{
    return computeSymmetricTensorNormSqr(tensorField, tensorField.getBoundingBox());
}

/* *************** Trace of each symmetric tensor of a field ************ */

template <typename T>
void computeSymmetricTensorTrace(
    MultiGridTensorField3D<T, 3> &tensorField, MultiGridScalarField3D<T> &trace, Box3D domain)
{
    applyProcessingFunctional(
        new ComputeSymmetricTensorTraceFunctional3D<T>, domain, trace, tensorField,
        tensorField.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > computeSymmetricTensorTrace(
    MultiGridTensorField3D<T, 3> &tensorField, Box3D domain)
{
    MultiGridScalarField3D<T> *trace = new MultiGridScalarField3D<T>(tensorField, domain);
    computeSymmetricTensorTrace(tensorField, *trace, domain);
    return std::unique_ptr<MultiGridScalarField3D<T> >(trace);
}

template <typename T>
std::unique_ptr<MultiGridScalarField3D<T> > computeSymmetricTensorTrace(
    MultiGridTensorField3D<T, 3> &tensorField)
{
    return computeSymmetricTensorTrace(tensorField, tensorField.getBoundingBox());
}

/* *************** Vorticity from Velocity field *********************** */

template <typename T>
void computeVorticity(
    MultiGridTensorField3D<T, 3> &velocity, MultiGridTensorField3D<T, 3> &vorticity, Box3D domain)
{
    plint envelopeWidth = 1;
    applyProcessingFunctional(
        new BoxVorticityFunctional3D<T, 3>, domain, vorticity, velocity,
        velocity.getReferenceLevel(), envelopeWidth);
}

template <typename T>
std::unique_ptr<MultiGridTensorField3D<T, 3> > computeVorticity(
    MultiGridTensorField3D<T, 3> &velocity, Box3D domain)
{
    MultiGridTensorField3D<T, 3> *vorticity = new MultiGridTensorField3D<T, 3>(velocity, domain);
    computeVorticity(velocity, *vorticity, domain);
    return std::unique_ptr<MultiGridTensorField3D<T, 3> >(vorticity);
}

template <typename T>
std::unique_ptr<MultiGridTensorField3D<T, 3> > computeVorticity(
    MultiGridTensorField3D<T, 3> &velocity)
{
    return computeVorticity(velocity, velocity.getBoundingBox());
}

/* *************** Vorticity, witout boundary treatment, from Velocity field  */

template <typename T>
void computeBulkVorticity(
    MultiGridTensorField3D<T, 3> &velocity, MultiGridTensorField3D<T, 3> &vorticity, Box3D domain)
{
    applyProcessingFunctional(
        new BoxBulkVorticityFunctional3D<T, 3>, domain, vorticity, velocity,
        velocity.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridTensorField3D<T, 3> > computeBulkVorticity(
    MultiGridTensorField3D<T, 3> &velocity, Box3D domain)
{
    MultiGridTensorField3D<T, 3> *vorticity = new MultiGridTensorField3D<T, 3>(velocity, domain);
    computeBulkVorticity(velocity, *vorticity, domain);
    return std::unique_ptr<MultiGridTensorField3D<T, 3> >(vorticity);
}

template <typename T>
std::unique_ptr<MultiGridTensorField3D<T, 3> > computeBulkVorticity(
    MultiGridTensorField3D<T, 3> &velocity)
{
    return computeBulkVorticity(velocity, velocity.getBoundingBox());
}

/* *************** Strain Rate from Velocity field ********************* */

template <typename T>
void computeStrainRate(
    MultiGridTensorField3D<T, 3> &velocity, MultiGridTensorField3D<T, 6> &S, Box3D domain)
{
    plint envelopeWidth = 1;
    applyProcessingFunctional(
        new BoxStrainRateFunctional3D<T, 3>, domain, velocity, S, velocity.getReferenceLevel(),
        envelopeWidth);
}

template <typename T>
std::unique_ptr<MultiGridTensorField3D<T, 6> > computeStrainRate(
    MultiGridTensorField3D<T, 3> &velocity, Box3D domain)
{
    MultiGridTensorField3D<T, 3> *S = new MultiGridTensorField3D<T, 3>(velocity, domain);
    computeStrainRate(velocity, *S, domain);
    return std::unique_ptr<MultiGridTensorField3D<T, 3> >(S);
}

template <typename T>
std::unique_ptr<MultiGridTensorField3D<T, 6> > computeStrainRate(
    MultiGridTensorField3D<T, 3> &velocity)
{
    return computeStrainRate(velocity, velocity.getBoundingBox());
}

/* *************** Str. rate, witout boundary treatment, from Velocity field  */

template <typename T>
void computeBulkStrainRate(
    MultiGridTensorField3D<T, 3> &velocity, MultiGridTensorField3D<T, 6> &S, Box3D domain)
{
    applyProcessingFunctional(
        new BoxBulkStrainRateFunctional3D<T, 3>, domain, velocity, S, velocity.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridTensorField3D<T, 6> > computeBulkStrainRate(
    MultiGridTensorField3D<T, 3> &velocity, Box3D domain)
{
    MultiGridTensorField3D<T, 6> *S = new MultiGridTensorField3D<T, 3>(velocity, domain);
    computeBulkStrainRate(velocity, *S, domain);
    return std::unique_ptr<MultiGridTensorField3D<T, 6> >(S);
}

template <typename T>
std::unique_ptr<MultiGridTensorField3D<T, 6> > computeBulkStrainRate(
    MultiGridTensorField3D<T, 3> &velocity)
{
    return computeBulkStrainRate(velocity, velocity.getBoundingBox());
}

/* *************** MultiTensorField - MultiTensorField operations *************** */

template <typename T1, typename T2, int nDim>
void copyConvert(
    MultiGridTensorField3D<T1, nDim> &field, MultiGridTensorField3D<T2, nDim> &convertedField,
    Box3D domain)
{
    applyProcessingFunctional(
        new CopyConvertTensorFunctional3D<T1, T2, nDim>, domain, field, convertedField,
        field.getReferenceLevel());
}

template <typename T1, typename T2, int nDim>
std::unique_ptr<MultiGridTensorField3D<T2, nDim> > copyConvert(
    MultiGridTensorField3D<T1, nDim> &field, Box3D domain)
{
    MultiGridTensorField3D<T2, nDim> *convertedField =
        new MultiGridTensorField3D<T2, nDim>(field, domain);
    copyConvert<T1, T2, nDim>(field, *convertedField, domain);
    return std::unique_ptr<MultiGridTensorField3D<T2, nDim> >(convertedField);
}

template <typename T1, typename T2, int nDim>
std::unique_ptr<MultiGridTensorField3D<T2, nDim> > copyConvert(
    MultiGridTensorField3D<T1, nDim> &field)
{
    return copyConvert<T1, T2, nDim>(field, field.getBoundingBox());
}

template <typename T, int nDim>
void add(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B,
    MultiGridTensorField3D<T, nDim> &result, Box3D domain)
{
    std::vector<MultiGridTensorField3D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(
        new Tensor_A_plus_B_functional3D<T, nDim>, domain, fields, A.getReferenceLevel());
}

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField3D<T, nDim> > add(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B, Box3D domain)
{
    MultiGridTensorField3D<T, nDim> *result = new MultiGridTensorField3D<T, nDim>(A, domain);
    add(A, B, *result, domain);
    return std::unique_ptr<MultiGridTensorField3D<T, nDim> >(result);
}

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField3D<T, nDim> > add(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B)
{
    return add(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void subtract(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B,
    MultiGridTensorField3D<T, nDim> &result, Box3D domain)
{
    std::vector<MultiGridTensorField3D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(
        new Tensor_A_minus_B_functional3D<T, nDim>, domain, fields, A.getReferenceLevel());
}

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField3D<T, nDim> > subtract(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B, Box3D domain)
{
    MultiGridTensorField3D<T, nDim> *result = new MultiGridTensorField3D<T, nDim>(A, domain);
    subtract(A, B, *result, domain);
    return std::unique_ptr<MultiGridTensorField3D<T, nDim> >(result);
}

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField3D<T, nDim> > subtract(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B)
{
    return subtract(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void multiply(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B,
    MultiGridTensorField3D<T, nDim> &result, Box3D domain)
{
    std::vector<MultiGridTensorField3D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(
        new Tensor_A_times_B_functional3D<T, nDim>, domain, fields, A.getReferenceLevel());
}

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField3D<T, nDim> > multiply(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B, Box3D domain)
{
    MultiGridTensorField3D<T, nDim> *result = new MultiGridTensorField3D<T, nDim>(A, domain);
    multiply(A, B, *result, domain);
    return std::unique_ptr<MultiGridTensorField3D<T, nDim> >(result);
}

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField3D<T, nDim> > multiply(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B)
{
    return multiply(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void divide(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B,
    MultiGridTensorField3D<T, nDim> &result, Box3D domain)
{
    std::vector<MultiGridTensorField3D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(
        new Tensor_A_dividedBy_B_functional3D<T, nDim>, domain, fields, A.getReferenceLevel());
}

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField3D<T, nDim> > divide(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B, Box3D domain)
{
    MultiGridTensorField3D<T, nDim> *result = new MultiGridTensorField3D<T, nDim>(A, domain);
    divide(A, B, *result, domain);
    return std::unique_ptr<MultiGridTensorField3D<T, nDim> >(result);
}

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField3D<T, nDim> > divide(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B)
{
    return divide(A, B, A.getBoundingBox());
}

/* *************** TensorField - TensorField inplace operations *************** */

template <typename T, int nDim>
void addInPlace(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B, Box3D domain)
{
    applyProcessingFunctional(
        new Tensor_A_plus_B_inplace_functional3D<T, nDim>, domain, A, B, A.getReferenceLevel());
}

template <typename T, int nDim>
void addInPlace(MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B)
{
    addInPlace(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void subtractInPlace(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B, Box3D domain)
{
    applyProcessingFunctional(new Tensor_A_minus_B_inplace_functional3D<T, nDim>, domain, A, B);
}

template <typename T, int nDim>
void subtractInPlace(MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B)
{
    subtractInPlace(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void multiplyInPlace(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B, Box3D domain)
{
    applyProcessingFunctional(
        new Tensor_A_times_B_inplace_functional3D<T, nDim>, domain, A, B, A.getReferenceLevel());
}

template <typename T, int nDim>
void multiplyInPlace(MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B)
{
    multiplyInPlace(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void multiplyInPlace(MultiGridTensorField3D<T, nDim> &A, T alpha, Box3D domain)
{
    applyProcessingFunctional(
        new Tensor_A_times_alpha_inplace_functional3D<T, nDim>(alpha), domain, A,
        A.getReferenceLevel());
}

template <typename T, int nDim>
void multiplyInPlace(MultiGridTensorField3D<T, nDim> &A, T alpha)
{
    multiplyInPlace(A, alpha, A.getBoundingBox());
}

template <typename T, int nDim>
void divideInPlace(
    MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B, Box3D domain)
{
    applyProcessingFunctional(
        new Tensor_A_dividedBy_B_inplace_functional3D<T, nDim>, domain, A, B,
        A.getReferenceLevel());
}

template <typename T, int nDim>
void divideInPlace(MultiGridTensorField3D<T, nDim> &A, MultiGridTensorField3D<T, nDim> &B)
{
    divideInPlace(A, B, A.getBoundingBox());
}

}  // namespace plb

#endif  // MULTI_GRID_DATA_ANALYSIS_WRAPPER_3D_HH
