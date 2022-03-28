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

#ifndef MULTI_GRID_DATA_ANALYSIS_WRAPPER_2D_HH
#define MULTI_GRID_DATA_ANALYSIS_WRAPPER_2D_HH

#include "dataProcessors/dataAnalysisWrapper2D.h"
#include "multiGrid/multiGridDataAnalysisWrapper2D.h"
#include "multiGrid/multiGridDataProcessorWrapper2D.h"

namespace plb {

/* *************** Reductive functions ******************************* */

template <typename T>
T computeAverage(MultiGridScalarField2D<T> &scalarField, Box2D domain)
{
    BoxScalarSumFunctional2D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField, scalarField.getReferenceLevel());
    return functional.getSumScalar() / (T)domain.nCells();
}

template <typename T>
T computeAverage(MultiGridScalarField2D<T> &scalarField)
{
    return computeAverage(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeMin(MultiGridScalarField2D<T> &scalarField, Box2D domain)
{
    BoxScalarMinFunctional2D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField, scalarField.getReferenceLevel());
    return functional.getMinScalar();
}

template <typename T>
T computeMin(MultiGridScalarField2D<T> &scalarField)
{
    return computeMin(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeMax(MultiGridScalarField2D<T> &scalarField, Box2D domain)
{
    BoxScalarMaxFunctional2D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField, scalarField.getReferenceLevel());
    return functional.getMaxScalar();
}

template <typename T>
T computeMax(MultiGridScalarField2D<T> &scalarField)
{
    return computeMax(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeBoundedAverage(MultiGridScalarField2D<T> &scalarField, Box2D domain)
{
    BoundedBoxScalarSumFunctional2D<T> functional;
    plint envelopeWidth = 1;
    applyProcessingFunctional(
        functional, domain, scalarField, scalarField.getReferenceLevel(), envelopeWidth);
    return functional.getSumScalar() / (T)((domain.getNx() - 1) * (domain.getNy() - 1));
}

template <typename T>
T computeBoundedAverage(MultiGridScalarField2D<T> &scalarField)
{
    return computeBoundedAverage(scalarField, scalarField.getBoundingBox());
}

template <typename T, class BoolMask>
plint count(MultiGridScalarField2D<T> &field, Box2D domain, BoolMask boolMask)
{
    CountScalarElementsFunctional2D<T, BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, field, field.getReferenceLevel());
    return functional.getCount();
}

template <typename T, class BoolMask>
plint count(MultiGridScalarField2D<T> &field, BoolMask boolMask)
{
    return count(field, field.getBoundingBox(), boolMask);
}

/* *************** Extract Sub-Lattice ******************************* */

template <typename T, template <typename U> class Descriptor>
void extractSubDomain(
    MultiGridLattice2D<T, Descriptor> &lattice, MultiGridLattice2D<T, Descriptor> &extractedLattice,
    Box2D domain)
{
    applyProcessingFunctional(
        new LatticeRegenerateFunctional2D<T, Descriptor>, domain, lattice, extractedLattice,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridLattice2D<T, Descriptor> > extractSubDomain(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    MultiGridLattice2D<T, Descriptor> *extractedLattice =
        new MultiGridLattice2D<T, Descriptor>(lattice, domain);
    extractSubDomain(lattice, *extractedLattice, domain);
    return std::unique_ptr<MultiGridLattice2D<T, Descriptor> >(extractedLattice);
}

/* *************** Density ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeDensity(
    MultiGridLattice2D<T, Descriptor> &lattice, MultiGridScalarField2D<T> &density, Box2D domain)
{
    applyProcessingFunctional(
        new BoxDensityFunctional2D<T, Descriptor>, domain, lattice, density,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeDensity(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    MultiGridScalarField2D<T> *density = new MultiGridScalarField2D<T>(lattice, domain);
    computeDensity(lattice, *density, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(density);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeDensity(
    MultiGridLattice2D<T, Descriptor> &lattice)
{
    return computeDensity(lattice, lattice.getBoundingBox());
}

/* *************** RhoBar ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeRhoBar(
    MultiGridLattice2D<T, Descriptor> &lattice, MultiGridScalarField2D<T> &rhoBar, Box2D domain)
{
    applyProcessingFunctional(
        new BoxRhoBarFunctional2D<T, Descriptor>, domain, lattice, rhoBar,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeRhoBar(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    MultiGridScalarField2D<T> *rhoBar = new MultiGridScalarField2D<T>(lattice, domain);
    computeRhoBar(lattice, *rhoBar, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(rhoBar);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeRhoBar(
    MultiGridLattice2D<T, Descriptor> &lattice)
{
    return computeRhoBar(lattice, lattice.getBoundingBox());
}

/* *************** Kinetic Energy ************************************ */

template <typename T, template <typename U> class Descriptor>
void computeKineticEnergy(
    MultiGridLattice2D<T, Descriptor> &lattice, MultiGridScalarField2D<T> &energy, Box2D domain)
{
    applyProcessingFunctional(
        new BoxKineticEnergyFunctional2D<T, Descriptor>, domain, lattice, energy,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeKineticEnergy(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    MultiGridScalarField2D<T> *energy = new MultiGridScalarField2D<T>(lattice, domain);
    computeKineticEnergy(lattice, *energy, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(energy);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeKineticEnergy(
    MultiGridLattice2D<T, Descriptor> &lattice)
{
    return computeKineticEnergy(lattice, lattice.getBoundingBox());
}

/* *************** Velocity Norm ************************************* */

template <typename T, template <typename U> class Descriptor>
void computeVelocityNorm(
    MultiGridLattice2D<T, Descriptor> &lattice, MultiGridScalarField2D<T> &velocityNorm,
    Box2D domain)
{
    applyProcessingFunctional(
        new BoxVelocityNormFunctional2D<T, Descriptor>, domain, lattice, velocityNorm,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeVelocityNorm(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    MultiGridScalarField2D<T> *velocityNorm = new MultiGridScalarField2D<T>(lattice, domain);
    computeVelocityNorm(lattice, *velocityNorm, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(velocityNorm);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeVelocityNorm(
    MultiGridLattice2D<T, Descriptor> &lattice)
{
    return computeVelocityNorm(lattice, lattice.getBoundingBox());
}

/* *************** Velocity Component ******************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocityComponent(
    MultiGridLattice2D<T, Descriptor> &lattice, MultiGridScalarField2D<T> &velocityComponent,
    Box2D domain, plint iComponent)
{
    applyProcessingFunctional(
        new BoxVelocityComponentFunctional2D<T, Descriptor>(iComponent), domain, lattice,
        velocityComponent, lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeVelocityComponent(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain, plint iComponent)
{
    MultiGridScalarField2D<T> *velocityComponent = new MultiGridScalarField2D<T>(lattice, domain);
    computeVelocityComponent(lattice, *velocityComponent, domain, iComponent);
    return std::unique_ptr<MultiGridScalarField2D<T> >(velocityComponent);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeVelocityComponent(
    MultiGridLattice2D<T, Descriptor> &lattice, plint iComponent)
{
    return computeVelocityComponent(lattice, lattice.getBoundingBox(), iComponent);
}

/* *************** Velocity ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocity(
    MultiGridLattice2D<T, Descriptor> &lattice,
    MultiGridTensorField2D<T, Descriptor<T>::d> &velocity, Box2D domain)
{
    applyProcessingFunctional(
        new BoxVelocityFunctional2D<T, Descriptor>, domain, lattice, velocity,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField2D<T, Descriptor<T>::d> > computeVelocity(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    MultiGridTensorField2D<T, Descriptor<T>::d> *velocity =
        new MultiGridTensorField2D<T, Descriptor<T>::d>(lattice, domain);
    computeVelocity(lattice, *velocity, domain);
    return std::unique_ptr<MultiGridTensorField2D<T, Descriptor<T>::d> >(velocity);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField2D<T, Descriptor<T>::d> > computeVelocity(
    MultiGridLattice2D<T, Descriptor> &lattice)
{
    return computeVelocity(lattice, lattice.getBoundingBox());
}

/* *************** PiNeq ********************************* */

template <typename T, template <typename U> class Descriptor>
void computePiNeq(
    MultiGridLattice2D<T, Descriptor> &lattice,
    MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq, Box2D domain)
{
    applyProcessingFunctional(
        new BoxPiNeqFunctional2D<T, Descriptor>, domain, lattice, PiNeq,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> > computePiNeq(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> *PiNeq =
        new MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n>(lattice, domain);
    computePiNeq(lattice, *PiNeq, domain);
    return std::unique_ptr<MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> >(PiNeq);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> > computePiNeq(
    MultiGridLattice2D<T, Descriptor> &lattice)
{
    return computePiNeq(lattice, lattice.getBoundingBox());
}

/* *************** Shear Stress ********************************* */

template <typename T, template <typename U> class Descriptor>
void computeShearStress(
    MultiGridLattice2D<T, Descriptor> &lattice,
    MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> &stress, Box2D domain)
{
    applyProcessingFunctional(
        new BoxShearStressFunctional2D<T, Descriptor>, domain, lattice, stress,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> > computeShearStress(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> *stress =
        new MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n>(lattice, domain);
    computeShearStress(lattice, *stress, domain);
    return std::unique_ptr<MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> >(stress);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> > computeShearStress(
    MultiGridLattice2D<T, Descriptor> &lattice)
{
    return computeShearStress(lattice, lattice.getBoundingBox());
}

/* *************** Strain Rate from Stress *************************** */

template <typename T, template <typename U> class Descriptor>
void computeStrainRateFromStress(
    MultiGridLattice2D<T, Descriptor> &lattice,
    MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> &S, Box2D domain)
{
    applyProcessingFunctional(
        new BoxStrainRateFromStressFunctional2D<T, Descriptor>, domain, lattice, S,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> >
    computeStrainRateFromStress(MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> *S =
        new MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n>(lattice, domain);
    computeStrainRateFromStress(lattice, *S, domain);
    return std::unique_ptr<MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> >(S);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField2D<T, SymmetricTensor<T, Descriptor>::n> >
    computeStrainRateFromStress(MultiGridLattice2D<T, Descriptor> &lattice)
{
    return computeStrainRateFromStress(lattice, lattice.getBoundingBox());
}

/* *************** Temperature ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeTemperature(
    MultiGridLattice2D<T, Descriptor> &lattice, MultiGridScalarField2D<T> &temperature,
    Box2D domain)
{
    applyProcessingFunctional(
        new BoxTemperatureFunctional2D<T, Descriptor>, domain, lattice, temperature,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeTemperature(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    MultiGridScalarField2D<T> *temperature = new MultiGridScalarField2D<T>(lattice, domain);
    computeTemperature(lattice, *temperature, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(temperature);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeTemperature(
    MultiGridLattice2D<T, Descriptor> &lattice)
{
    return computeTemperature(lattice, lattice.getBoundingBox());
}

/* *************** SoundSpeed ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeSoundSpeed(
    MultiGridLattice2D<T, Descriptor> &lattice, MultiGridScalarField2D<T> &soundSpeed, Box2D domain)
{
    applyProcessingFunctional(
        new BoxSoundSpeedFunctional2D<T, Descriptor>, domain, lattice, soundSpeed,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeSoundSpeed(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    MultiGridScalarField2D<T> *soundSpeed = new MultiGridScalarField2D<T>(lattice, domain);
    computeSoundSpeed(lattice, *soundSpeed, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(soundSpeed);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeSoundSpeed(
    MultiGridLattice2D<T, Descriptor> &lattice)
{
    return computeSoundSpeed(lattice, lattice.getBoundingBox());
}

/* *************** Population **************************************** */

template <typename T, template <typename U> class Descriptor>
void computePopulation(
    MultiGridLattice2D<T, Descriptor> &lattice, MultiGridScalarField2D<T> &population, Box2D domain,
    plint iPop)
{
    applyProcessingFunctional(
        new BoxPopulationFunctional2D<T, Descriptor>(iPop), domain, lattice, population,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computePopulation(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain, plint iPop)
{
    MultiGridScalarField2D<T> *population = new MultiGridScalarField2D<T>(lattice, domain);
    computePopulation(lattice, *population, domain, iPop);
    return std::unique_ptr<MultiGridScalarField2D<T> >(population);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computePopulation(
    MultiGridLattice2D<T, Descriptor> &lattice, plint iPop)
{
    return computePopulation(lattice, lattice.getBoundingBox(), iPop);
}

template <typename T, template <typename U> class Descriptor>
void computeAllPopulations(
    MultiGridLattice2D<T, Descriptor> &lattice,
    MultiGridTensorField2D<T, Descriptor<T>::q> &populations, Box2D domain)
{
    applyProcessingFunctional(
        new BoxAllPopulationsFunctional2D<T, Descriptor>(), domain, lattice, populations,
        lattice.getReferenceLevel());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridScalarField2D<T> > computeAllPopulations(
    MultiGridLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    MultiGridTensorField2D<T, Descriptor<T>::q> *populations =
        new MultiGridTensorField2D<T, Descriptor<T>::q>(lattice, domain);
    computeAllPopulations(lattice, *populations, domain);
    return std::unique_ptr<MultiGridTensorField2D<T, Descriptor<T>::q> >(populations);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridTensorField2D<T, Descriptor<T>::q> > computeAllPopulations(
    MultiGridLattice2D<T, Descriptor> &lattice)
{
    return computeAllPopulations(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void copyPopulations(
    MultiGridLattice2D<T, Descriptor> &latticeFrom, MultiGridLattice2D<T, Descriptor> &latticeTo,
    Box2D domain)
{
    applyProcessingFunctional(
        new CopyPopulationsFunctional2D<T, Descriptor>(), domain, latticeFrom, latticeTo);
}

/* *************** Extract Sub-ScalarField *************************** */

template <typename T>
void extractSubDomain(
    MultiGridScalarField2D<T> &field, MultiGridScalarField2D<T> &extractedField, Box2D domain)
{
    applyProcessingFunctional(
        new ExtractScalarSubDomainFunctional2D<T>, domain, field, extractedField,
        field.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > extractSubDomain(
    MultiGridScalarField2D<T> &field, Box2D domain)
{
    MultiGridScalarField2D<T> *extractedField = new MultiGridScalarField2D<T>(field, domain);
    extractSubDomain(field, *extractedField, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(extractedField);
}

/* *************** MultiScalarField - Scalar operations *************** */

template <typename T>
void add(
    MultiGridScalarField2D<T> &field, T scalar, MultiGridScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(
        new A_plus_alpha_functional2D<T>(scalar), domain, field, result, field.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > add(
    MultiGridScalarField2D<T> &field, T scalar, Box2D domain)
{
    MultiGridScalarField2D<T> *result = new MultiGridScalarField2D<T>(field, domain);
    add(field, scalar, *result, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > add(MultiGridScalarField2D<T> &field, T scalar)
{
    return add(field, scalar, field.getBoundingBox());
}

template <typename T>
void add(
    T scalar, MultiGridScalarField2D<T> &field, MultiGridScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(
        new A_plus_alpha_functional2D<T>(scalar), domain, field, result, field.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > add(
    T scalar, MultiGridScalarField2D<T> &field, Box2D domain)
{
    MultiGridScalarField2D<T> *result = new MultiGridScalarField2D<T>(field, domain);
    add(scalar, field, *result, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > add(T scalar, MultiGridScalarField2D<T> &field)
{
    return add(scalar, field, field.getBoundingBox());
}

template <typename T>
void subtract(
    MultiGridScalarField2D<T> &field, T scalar, MultiGridScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(
        new A_minus_alpha_functional2D<T>(scalar), domain, field, result,
        field.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > subtract(
    MultiGridScalarField2D<T> &field, T scalar, Box2D domain)
{
    MultiGridScalarField2D<T> *result = new MultiGridScalarField2D<T>(field, domain);
    subtract(field, scalar, *result, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > subtract(MultiGridScalarField2D<T> &field, T scalar)
{
    return subtract(field, scalar, field.getBoundingBox());
}

template <typename T>
void subtract(
    T scalar, MultiGridScalarField2D<T> &field, MultiGridScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(
        new Alpha_minus_A_functional2D<T>(scalar), domain, field, result,
        field.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > subtract(
    T scalar, MultiGridScalarField2D<T> &field, Box2D domain)
{
    MultiGridScalarField2D<T> *result = new MultiGridScalarField2D<T>(field, domain);
    subtract(scalar, field, *result, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > subtract(T scalar, MultiGridScalarField2D<T> &field)
{
    return subtract(scalar, field, field.getBoundingBox());
}

template <typename T>
void multiply(
    MultiGridScalarField2D<T> &field, T scalar, MultiGridScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(
        new A_times_alpha_functional2D<T>(scalar), domain, field, result,
        field.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > multiply(
    MultiGridScalarField2D<T> &field, T scalar, Box2D domain)
{
    MultiGridScalarField2D<T> *result = new MultiGridScalarField2D<T>(field, domain);
    multiply(field, scalar, *result, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > multiply(MultiGridScalarField2D<T> &field, T scalar)
{
    return multiply(field, scalar, field.getBoundingBox());
}

template <typename T>
void multiply(
    T scalar, MultiGridScalarField2D<T> &field, MultiGridScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(
        new A_times_alpha_functional2D<T>(scalar), domain, field, result,
        field.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > multiply(
    T scalar, MultiGridScalarField2D<T> &field, Box2D domain)
{
    MultiGridScalarField2D<T> *result = new MultiGridScalarField2D<T>(field, domain);
    multiply(scalar, field, *result, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > multiply(T scalar, MultiGridScalarField2D<T> &field)
{
    return multiply(scalar, field, field.getBoundingBox());
}

template <typename T>
void divide(
    MultiGridScalarField2D<T> &field, T scalar, MultiGridScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(
        new A_dividedBy_alpha_functional2D<T>(scalar), domain, field, result,
        field.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > divide(
    MultiGridScalarField2D<T> &field, T scalar, Box2D domain)
{
    MultiGridScalarField2D<T> *result = new MultiGridScalarField2D<T>(field, domain);
    divide(field, scalar, *result, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > divide(MultiGridScalarField2D<T> &field, T scalar)
{
    return divide(field, scalar, field.getBoundingBox());
}

template <typename T>
void divide(
    T scalar, MultiGridScalarField2D<T> &field, MultiGridScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(
        new Alpha_dividedBy_A_functional2D<T>(scalar), domain, field, result,
        field.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > divide(
    T scalar, MultiGridScalarField2D<T> &field, Box2D domain)
{
    MultiGridScalarField2D<T> *result = new MultiGridScalarField2D<T>(field, domain);
    divide(scalar, field, *result, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > divide(T scalar, MultiGridScalarField2D<T> &field)
{
    return divide(scalar, field, field.getBoundingBox());
}

/* *************** MultiScalarField - Scalar inplace operations *************** */

template <typename T>
void addInPlace(MultiGridScalarField2D<T> &field, T scalar, Box2D domain)
{
    applyProcessingFunctional(
        new A_plus_alpha_inplace_functional2D<T>(scalar), domain, field, field.getReferenceLevel());
}

template <typename T>
void addInPlace(MultiGridScalarField2D<T> &field, T scalar)
{
    addInPlace(field, scalar, field.getBoundingBox());
}

template <typename T>
void subtractInPlace(MultiGridScalarField2D<T> &field, T scalar, Box2D domain)
{
    applyProcessingFunctional(
        new A_minus_alpha_inplace_functional2D<T>(scalar), domain, field,
        field.getReferenceLevel());
}

template <typename T>
void subtractInPlace(MultiGridScalarField2D<T> &field, T scalar)
{
    subtractInPlace(field, scalar, field.getBoundingBox());
}

template <typename T>
void multiplyInPlace(MultiGridScalarField2D<T> &field, T scalar, Box2D domain)
{
    applyProcessingFunctional(
        new A_times_alpha_inplace_functional2D<T>(scalar), domain, field,
        field.getReferenceLevel());
}

template <typename T>
void multiplyInPlace(MultiGridScalarField2D<T> &field, T scalar)
{
    multiplyInPlace(field, scalar, field.getBoundingBox());
}

template <typename T>
void divideInPlace(MultiGridScalarField2D<T> &field, T scalar, Box2D domain)
{
    applyProcessingFunctional(
        new A_dividedBy_alpha_inplace_functional2D<T>(scalar), domain, field,
        field.getReferenceLevel());
}

template <typename T>
void divideInPlace(MultiGridScalarField2D<T> &field, T scalar)
{
    divideInPlace(field, scalar, field.getBoundingBox());
}

/* *************** MultiScalarField - MultiScalarField operations *************** */

template <typename T1, typename T2>
void copyConvert(
    MultiGridScalarField2D<T1> &field, MultiGridScalarField2D<T2> &convertedField, Box2D domain)
{
    applyProcessingFunctional(
        new CopyConvertScalarFunctional2D<T1, T2>, domain, field, convertedField,
        field.getReferenceLevel());
}

template <typename T1, typename T2>
std::unique_ptr<MultiGridScalarField2D<T2> > copyConvert(
    MultiGridScalarField2D<T1> &field, Box2D domain)
{
    MultiGridScalarField2D<T2> *convertedField = new MultiGridScalarField2D<T2>(field, domain);
    copyConvert(field, *convertedField, domain);
    return std::unique_ptr<MultiGridScalarField2D<T2> >(convertedField);
}

template <typename T1, typename T2>
std::unique_ptr<MultiGridScalarField2D<T2> > copyConvert(MultiGridScalarField2D<T1> &field)
{
    return copyConvert<T1, T2>(field, field.getBoundingBox());
}

template <typename T>
void add(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, MultiGridScalarField2D<T> &result,
    Box2D domain)
{
    std::vector<MultiGridScalarField2D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_plus_B_functional2D<T>, domain, fields, A.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > add(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, Box2D domain)
{
    MultiGridScalarField2D<T> *result = new MultiGridScalarField2D<T>(A, domain);
    add(A, B, *result, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > add(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B)
{
    return add(A, B, A.getBoundingBox());
}

template <typename T>
void subtract(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, MultiGridScalarField2D<T> &result,
    Box2D domain)
{
    std::vector<MultiGridScalarField2D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_minus_B_functional2D<T>, domain, fields, A.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > subtract(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, Box2D domain)
{
    MultiGridScalarField2D<T> *result = new MultiGridScalarField2D<T>(A, domain);
    subtract(A, B, *result, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > subtract(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B)
{
    return subtract(A, B, A.getBoundingBox());
}

template <typename T>
void multiply(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, MultiGridScalarField2D<T> &result,
    Box2D domain)
{
    std::vector<MultiGridScalarField2D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_times_B_functional2D<T>, domain, fields, A.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > multiply(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, Box2D domain)
{
    MultiGridScalarField2D<T> *result = new MultiGridScalarField2D<T>(A, domain);
    multiply(A, B, *result, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > multiply(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B)
{
    return multiply(A, B, A.getBoundingBox());
}

template <typename T>
void divide(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, MultiGridScalarField2D<T> &result,
    Box2D domain)
{
    std::vector<MultiGridScalarField2D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(
        new A_dividedBy_B_functional2D<T>, domain, fields, A.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > divide(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, Box2D domain)
{
    MultiGridScalarField2D<T> *result = new MultiGridScalarField2D<T>(A, domain);
    divide(A, B, *result, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(result);
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > divide(
    MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B)
{
    return divide(A, B, A.getBoundingBox());
}

/* *************** ScalarField - ScalarField inplace operations *************** */

template <typename T>
void addInPlace(MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, Box2D domain)
{
    applyProcessingFunctional(
        new A_plus_B_inplace_functional2D<T>, domain, A, B, A.getReferenceLevel());
}

template <typename T>
void addInPlace(MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B)
{
    addInPlace(A, B, A.getBoundingBox());
}

template <typename T>
void subtractInPlace(MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, Box2D domain)
{
    applyProcessingFunctional(
        new A_minus_B_inplace_functional2D<T>, domain, A, B, A.getReferenceLevel());
}

template <typename T>
void subtractInPlace(MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B)
{
    subtractInPlace(A, B, A.getBoundingBox());
}

template <typename T>
void multiplyInPlace(MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, Box2D domain)
{
    applyProcessingFunctional(
        new A_times_B_inplace_functional2D<T>, domain, A, B, A.getReferenceLevel());
}

template <typename T>
void multiplyInPlace(MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B)
{
    multiplyInPlace(A, B, A.getBoundingBox());
}

template <typename T>
void divideInPlace(MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B, Box2D domain)
{
    applyProcessingFunctional(
        new A_dividedBy_B_inplace_functional2D<T>, domain, A, B, A.getReferenceLevel());
}

template <typename T>
void divideInPlace(MultiGridScalarField2D<T> &A, MultiGridScalarField2D<T> &B)
{
    divideInPlace(A, B, A.getBoundingBox());
}

/* ******************************************************************* */
/* *************** PART VI : Multi-block wrappers: Tensor-field ****** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template <typename T, plint nDim, class BoolMask>
plint count(MultiGridTensorField2D<T, nDim> &field, Box2D domain, BoolMask boolMask)
{
    CountTensorElementsFunctional2D<T, nDim, BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, field, field.getReferenceLevel());
    return functional.getCount();
}

template <typename T, plint nDim, class BoolMask>
plint count(MultiGridTensorField2D<T, nDim> &field, BoolMask boolMask)
{
    return count(field, field.getBoundingBox(), boolMask);
}

/* *************** Extract Sub-TensorField *************************** */

template <typename T, int nDim>
void extractSubDomain(
    MultiGridTensorField2D<T, nDim> &field, MultiGridTensorField2D<T, nDim> &extractedField,
    Box2D domain)
{
    applyProcessingFunctional(
        new ExtractTensorSubDomainFunctional2D<T, nDim>, domain, field, extractedField,
        field.getReferenceLevel());
}

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField2D<T, nDim> > extractSubDomain(
    MultiGridTensorField2D<T, nDim> &field, Box2D domain)
{
    MultiGridTensorField2D<T, nDim> *extractedField =
        new MultiGridTensorField2D<T, nDim>(field, domain);
    extractSubDomain(field, *extractedField, domain);
    return std::unique_ptr<MultiGridTensorField2D<T, nDim> >(extractedField);
}

/* *************** Component (scalar-field) out of a tensor-field ****** */

template <typename T, int nDim>
void extractComponent(
    MultiGridTensorField2D<T, nDim> &tensorField, MultiGridScalarField2D<T> &component,
    Box2D domain, int iComponent)
{
    applyProcessingFunctional(
        new ExtractTensorComponentFunctional2D<T, nDim>(iComponent), domain, component, tensorField,
        tensorField.getReferenceLevel());
}

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField2D<T> > extractComponent(
    MultiGridTensorField2D<T, nDim> &tensorField, Box2D domain, int iComponent)
{
    MultiGridScalarField2D<T> *component = new MultiGridScalarField2D<T>(tensorField, domain);
    extractComponent(tensorField, *component, domain, iComponent);
    return std::unique_ptr<MultiGridScalarField2D<T> >(component);
}

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField2D<T> > extractComponent(
    MultiGridTensorField2D<T, nDim> &tensorField, int iComponent)
{
    return extractComponent(tensorField, tensorField.getBoundingBox(), iComponent);
}

/* *************** Vector-norm of each cell in the field *************** */

template <typename T, int nDim>
void computeNorm(
    MultiGridTensorField2D<T, nDim> &tensorField, MultiGridScalarField2D<T> &component,
    Box2D domain)
{
    applyProcessingFunctional(
        new ComputeNormFunctional2D<T, nDim>, domain, component, tensorField,
        tensorField.getReferenceLevel());
}

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField2D<T> > computeNorm(
    MultiGridTensorField2D<T, nDim> &tensorField, Box2D domain)
{
    MultiGridScalarField2D<T> *component = new MultiGridScalarField2D<T>(tensorField, domain);
    computeNorm(tensorField, *component, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(component);
}

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField2D<T> > computeNorm(
    MultiGridTensorField2D<T, nDim> &tensorField)
{
    return computeNorm(tensorField, tensorField.getBoundingBox());
}

/* *************** Squared vector-norm of each cell in the field ******** */

template <typename T, int nDim>
void computeNormSqr(
    MultiGridTensorField2D<T, nDim> &tensorField, MultiGridScalarField2D<T> &component,
    Box2D domain)
{
    applyProcessingFunctional(
        new ComputeNormSqrFunctional2D<T, nDim>, domain, component, tensorField,
        tensorField.getReferenceLevel());
}

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField2D<T> > computeNormSqr(
    MultiGridTensorField2D<T, nDim> &tensorField, Box2D domain)
{
    MultiGridScalarField2D<T> *component = new MultiGridScalarField2D<T>(tensorField, domain);
    computeNormSqr(tensorField, *component, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(component);
}

template <typename T, int nDim>
std::unique_ptr<MultiGridScalarField2D<T> > computeNormSqr(
    MultiGridTensorField2D<T, nDim> &tensorField)
{
    return computeNormSqr(tensorField, tensorField.getBoundingBox());
}

/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template <typename T>
void computeSymmetricTensorNorm(
    MultiGridTensorField2D<T, 3> &tensorField, MultiGridScalarField2D<T> &norm, Box2D domain)
{
    applyProcessingFunctional(
        new ComputeSymmetricTensorNormFunctional2D<T>, domain, norm, tensorField,
        tensorField.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > computeSymmetricTensorNorm(
    MultiGridTensorField2D<T, 3> &tensorField, Box2D domain)
{
    MultiGridScalarField2D<T> *norm = new MultiGridScalarField2D<T>(tensorField, domain);
    computeSymmetricTensorNorm(tensorField, *norm, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(norm);
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > computeSymmetricTensorNorm(
    MultiGridTensorField2D<T, 3> &tensorField)
{
    return computeSymmetricTensorNorm(tensorField, tensorField.getBoundingBox());
}

/* *************** Squared Tensor-norm of each symmetric tensor of a field*/

template <typename T>
void computeSymmetricTensorNormSqr(
    MultiGridTensorField2D<T, 3> &tensorField, MultiGridScalarField2D<T> &normSqr, Box2D domain)
{
    applyProcessingFunctional(
        new ComputeSymmetricTensorNormSqrFunctional2D<T>, domain, normSqr, tensorField,
        tensorField.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > computeSymmetricTensorNormSqr(
    MultiGridTensorField2D<T, 3> &tensorField, Box2D domain)
{
    MultiGridScalarField2D<T> *normSqr = new MultiGridScalarField2D<T>(tensorField, domain);
    computeSymmetricTensorNormSqr(tensorField, *normSqr, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(normSqr);
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > computeSymmetricTensorNormSqr(
    MultiGridTensorField2D<T, 3> &tensorField)
{
    return computeSymmetricTensorNormSqr(tensorField, tensorField.getBoundingBox());
}

/* *************** Trace of each symmetric tensor of a field ************ */

template <typename T>
void computeSymmetricTensorTrace(
    MultiGridTensorField2D<T, 3> &tensorField, MultiGridScalarField2D<T> &trace, Box2D domain)
{
    applyProcessingFunctional(
        new ComputeSymmetricTensorTraceFunctional2D<T>, domain, trace, tensorField,
        tensorField.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > computeSymmetricTensorTrace(
    MultiGridTensorField2D<T, 3> &tensorField, Box2D domain)
{
    MultiGridScalarField2D<T> *trace = new MultiGridScalarField2D<T>(tensorField, domain);
    computeSymmetricTensorTrace(tensorField, *trace, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(trace);
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > computeSymmetricTensorTrace(
    MultiGridTensorField2D<T, 3> &tensorField)
{
    return computeSymmetricTensorTrace(tensorField, tensorField.getBoundingBox());
}

/* *************** Vorticity from Velocity field *********************** */

template <typename T>
void computeVorticity(
    MultiGridTensorField2D<T, 2> &velocity, MultiGridScalarField2D<T> &vorticity, Box2D domain)
{
    plint envelopeWidth = 1;
    applyProcessingFunctional(
        new BoxVorticityFunctional2D<T, 2>, domain, vorticity, velocity,
        velocity.getReferenceLevel(), envelopeWidth);
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > computeVorticity(
    MultiGridTensorField2D<T, 2> &velocity, Box2D domain)
{
    MultiGridScalarField2D<T> *vorticity = new MultiGridScalarField2D<T>(velocity, domain);
    computeVorticity(velocity, *vorticity, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(vorticity);
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > computeVorticity(MultiGridTensorField2D<T, 2> &velocity)
{
    return computeVorticity(velocity, velocity.getBoundingBox());
}

/* *************** Vorticity, witout boundary treatment, from Velocity field  */

template <typename T>
void computeBulkVorticity(
    MultiGridTensorField2D<T, 2> &velocity, MultiGridScalarField2D<T> &vorticity, Box2D domain)
{
    applyProcessingFunctional(
        new BoxBulkVorticityFunctional2D<T, 2>, domain, vorticity, velocity,
        velocity.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > computeBulkVorticity(
    MultiGridTensorField2D<T, 2> &velocity, Box2D domain)
{
    MultiGridScalarField2D<T> *vorticity = new MultiGridScalarField2D<T>(velocity, domain);
    computeBulkVorticity(velocity, *vorticity, domain);
    return std::unique_ptr<MultiGridScalarField2D<T> >(vorticity);
}

template <typename T>
std::unique_ptr<MultiGridScalarField2D<T> > computeBulkVorticity(
    MultiGridTensorField2D<T, 2> &velocity)
{
    return computeBulkVorticity(velocity, velocity.getBoundingBox());
}

/* *************** Strain Rate from Velocity field ********************* */

template <typename T>
void computeStrainRate(
    MultiGridTensorField2D<T, 2> &velocity, MultiGridTensorField2D<T, 3> &S, Box2D domain)
{
    plint envelopeWidth = 1;
    applyProcessingFunctional(
        new BoxStrainRateFunctional2D<T, 2>, domain, velocity, S, velocity.getReferenceLevel(),
        envelopeWidth);
}

template <typename T>
std::unique_ptr<MultiGridTensorField2D<T, 3> > computeStrainRate(
    MultiGridTensorField2D<T, 2> &velocity, Box2D domain)
{
    MultiGridTensorField2D<T, 3> *S = new MultiGridTensorField2D<T, 3>(velocity, domain);
    computeStrainRate(velocity, *S, domain);
    return std::unique_ptr<MultiGridTensorField2D<T, 3> >(S);
}

template <typename T>
std::unique_ptr<MultiGridTensorField2D<T, 3> > computeStrainRate(
    MultiGridTensorField2D<T, 2> &velocity)
{
    return computeStrainRate(velocity, velocity.getBoundingBox());
}

/* *************** Str. rate, witout boundary treatment, from Velocity field  */

template <typename T>
void computeBulkStrainRate(
    MultiGridTensorField2D<T, 2> &velocity, MultiGridTensorField2D<T, 3> &S, Box2D domain)
{
    applyProcessingFunctional(
        new BoxBulkStrainRateFunctional2D<T, 2>, domain, velocity, S, velocity.getReferenceLevel());
}

template <typename T>
std::unique_ptr<MultiGridTensorField2D<T, 3> > computeBulkStrainRate(
    MultiGridTensorField2D<T, 2> &velocity, Box2D domain)
{
    MultiGridTensorField2D<T, 3> *S = new MultiGridTensorField2D<T, 3>(velocity, domain);
    computeBulkStrainRate(velocity, *S, domain);
    return std::unique_ptr<MultiGridTensorField2D<T, 3> >(S);
}

template <typename T>
std::unique_ptr<MultiGridTensorField2D<T, 3> > computeBulkStrainRate(
    MultiGridTensorField2D<T, 2> &velocity)
{
    return computeBulkStrainRate(velocity, velocity.getBoundingBox());
}

/* *************** MultiTensorField - MultiTensorField operations *************** */

template <typename T1, typename T2, int nDim>
void copyConvert(
    MultiGridTensorField2D<T1, nDim> &field, MultiGridTensorField2D<T2, nDim> &convertedField,
    Box2D domain)
{
    applyProcessingFunctional(
        new CopyConvertTensorFunctional2D<T1, T2, nDim>, domain, field, convertedField,
        field.getReferenceLevel());
}

template <typename T1, typename T2, int nDim>
std::unique_ptr<MultiGridTensorField2D<T2, nDim> > copyConvert(
    MultiGridTensorField2D<T1, nDim> &field, Box2D domain)
{
    MultiGridTensorField2D<T2, nDim> *convertedField =
        new MultiGridTensorField2D<T2, nDim>(field, domain);
    copyConvert<T1, T2, nDim>(field, *convertedField, domain);
    return std::unique_ptr<MultiGridTensorField2D<T2, nDim> >(convertedField);
}

template <typename T1, typename T2, int nDim>
std::unique_ptr<MultiGridTensorField2D<T2, nDim> > copyConvert(
    MultiGridTensorField2D<T1, nDim> &field)
{
    return copyConvert<T1, T2, nDim>(field, field.getBoundingBox());
}

template <typename T, int nDim>
void add(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B,
    MultiGridTensorField2D<T, nDim> &result, Box2D domain)
{
    std::vector<MultiGridTensorField2D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(
        new Tensor_A_plus_B_functional2D<T, nDim>, domain, fields, A.getReferenceLevel());
}

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField2D<T, nDim> > add(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B, Box2D domain)
{
    MultiGridTensorField2D<T, nDim> *result = new MultiGridTensorField2D<T, nDim>(A, domain);
    add(A, B, *result, domain);
    return std::unique_ptr<MultiGridTensorField2D<T, nDim> >(result);
}

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField2D<T, nDim> > add(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B)
{
    return add(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void subtract(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B,
    MultiGridTensorField2D<T, nDim> &result, Box2D domain)
{
    std::vector<MultiGridTensorField2D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(
        new Tensor_A_minus_B_functional2D<T, nDim>, domain, fields, A.getReferenceLevel());
}

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField2D<T, nDim> > subtract(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B, Box2D domain)
{
    MultiGridTensorField2D<T, nDim> *result = new MultiGridTensorField2D<T, nDim>(A, domain);
    subtract(A, B, *result, domain);
    return std::unique_ptr<MultiGridTensorField2D<T, nDim> >(result);
}

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField2D<T, nDim> > subtract(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B)
{
    return subtract(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void multiply(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B,
    MultiGridTensorField2D<T, nDim> &result, Box2D domain)
{
    std::vector<MultiGridTensorField2D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(
        new Tensor_A_times_B_functional2D<T, nDim>, domain, fields, A.getReferenceLevel());
}

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField2D<T, nDim> > multiply(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B, Box2D domain)
{
    MultiGridTensorField2D<T, nDim> *result = new MultiGridTensorField2D<T, nDim>(A, domain);
    multiply(A, B, *result, domain);
    return std::unique_ptr<MultiGridTensorField2D<T, nDim> >(result);
}

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField2D<T, nDim> > multiply(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B)
{
    return multiply(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void divide(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B,
    MultiGridTensorField2D<T, nDim> &result, Box2D domain)
{
    std::vector<MultiGridTensorField2D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(
        new Tensor_A_dividedBy_B_functional2D<T, nDim>, domain, fields, A.getReferenceLevel());
}

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField2D<T, nDim> > divide(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B, Box2D domain)
{
    MultiGridTensorField2D<T, nDim> *result = new MultiGridTensorField2D<T, nDim>(A, domain);
    divide(A, B, *result, domain);
    return std::unique_ptr<MultiGridTensorField2D<T, nDim> >(result);
}

template <typename T, int nDim>
std::unique_ptr<MultiGridTensorField2D<T, nDim> > divide(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B)
{
    return divide(A, B, A.getBoundingBox());
}

/* *************** TensorField - TensorField inplace operations *************** */

template <typename T, int nDim>
void addInPlace(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B, Box2D domain)
{
    applyProcessingFunctional(
        new Tensor_A_plus_B_inplace_functional2D<T, nDim>, domain, A, B, A.getReferenceLevel());
}

template <typename T, int nDim>
void addInPlace(MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B)
{
    addInPlace(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void subtractInPlace(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B, Box2D domain)
{
    applyProcessingFunctional(new Tensor_A_minus_B_inplace_functional2D<T, nDim>, domain, A, B);
}

template <typename T, int nDim>
void subtractInPlace(MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B)
{
    subtractInPlace(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void multiplyInPlace(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B, Box2D domain)
{
    applyProcessingFunctional(
        new Tensor_A_times_B_inplace_functional2D<T, nDim>, domain, A, B, A.getReferenceLevel());
}

template <typename T, int nDim>
void multiplyInPlace(MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B)
{
    multiplyInPlace(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void multiplyInPlace(MultiGridTensorField2D<T, nDim> &A, T alpha, Box2D domain)
{
    applyProcessingFunctional(
        new Tensor_A_times_alpha_inplace_functional2D<T, nDim>(alpha), domain, A,
        A.getReferenceLevel());
}

template <typename T, int nDim>
void multiplyInPlace(MultiGridTensorField2D<T, nDim> &A, T alpha)
{
    multiplyInPlace(A, alpha, A.getBoundingBox());
}

template <typename T, int nDim>
void divideInPlace(
    MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B, Box2D domain)
{
    applyProcessingFunctional(
        new Tensor_A_dividedBy_B_inplace_functional2D<T, nDim>, domain, A, B,
        A.getReferenceLevel());
}

template <typename T, int nDim>
void divideInPlace(MultiGridTensorField2D<T, nDim> &A, MultiGridTensorField2D<T, nDim> &B)
{
    divideInPlace(A, B, A.getBoundingBox());
}

}  // namespace plb

#endif  // MULTI_GRID_DATA_ANALYSIS_WRAPPER_2D_HH
