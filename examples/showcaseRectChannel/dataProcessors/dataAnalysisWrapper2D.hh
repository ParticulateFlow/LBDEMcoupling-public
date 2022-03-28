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

#ifndef DATA_ANALYSIS_WRAPPER_2D_HH
#define DATA_ANALYSIS_WRAPPER_2D_HH

#include "atomicBlock/dataProcessorWrapper2D.h"
#include "atomicBlock/reductiveDataProcessorWrapper2D.h"
#include "dataProcessors/dataAnalysisFunctional2D.h"
#include "dataProcessors/dataAnalysisWrapper2D.h"
#include "multiBlock/multiBlockGenerator2D.h"
#include "multiBlock/multiDataProcessorWrapper2D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper2D.h"

namespace plb {

/* ******************************************************************* */
/* *************** PART I. Atomic-block wrappers: Block-Lattice ****** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template <typename T, template <typename U> class Descriptor>
T computeAverageDensity(BlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    BoxSumRhoBarFunctional2D<T, Descriptor> functional;
    applyProcessingFunctional(functional, domain, lattice);
    return Descriptor<T>::fullRho(functional.getSumRhoBar() / (T)domain.nCells());
}

template <typename T, template <typename U> class Descriptor>
T computeAverageDensity(BlockLattice2D<T, Descriptor> &lattice)
{
    return computeAverageDensity(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(BlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    BoxSumRhoBarFunctional2D<T, Descriptor> functional;
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getSumRhoBar() / (T)domain.nCells();
}

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(BlockLattice2D<T, Descriptor> &lattice)
{
    return computeAverageRhoBar(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(BlockLattice2D<T, Descriptor> &lattice, DotList2D dotList)
{
    DotSumRhoBarFunctional2D<T, Descriptor> functional;
    applyProcessingFunctional(functional, dotList, lattice);
    return functional.getSumRhoBar() / (T)dotList.getN();
}

template <typename T, template <typename U> class Descriptor>
T computeAverageEnergy(BlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    BoxSumEnergyFunctional2D<T, Descriptor> functional;
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getSumEnergy() / (T)domain.nCells();
    ;
}

template <typename T, template <typename U> class Descriptor>
T computeAverageEnergy(BlockLattice2D<T, Descriptor> &lattice)
{
    return computeAverageEnergy(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
T computeAverageVelocityComponent(
    BlockLattice2D<T, Descriptor> &lattice, DotList2D dotList, plint iComponent)
{
    DotSumVelocityComponentFunctional2D<T, Descriptor> functional(iComponent);
    applyProcessingFunctional(functional, dotList, lattice);
    return functional.getSumVelocityComponent() / (T)dotList.getN();
}

template <typename T, template <typename U> class Descriptor, class BoolMask>
plint count(BlockLattice2D<T, Descriptor> &lattice, Box2D domain, BoolMask boolMask)
{
    CountLatticeElementsFunctional2D<T, Descriptor, BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getCount();
}

template <typename T, template <typename U> class Descriptor, class BoolMask>
plint count(BlockLattice2D<T, Descriptor> &lattice, BoolMask boolMask)
{
    return count(lattice, lattice.getBoundingBox(), boolMask);
}

/* *************** Density ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeDensity(BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &density)
{
    applyProcessingFunctional(
        new BoxDensityFunctional2D<T, Descriptor>, lattice.getBoundingBox(), lattice, density);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField2D<T> > computeDensity(BlockLattice2D<T, Descriptor> &lattice)
{
    ScalarField2D<T> *density = new ScalarField2D<T>(lattice.getNx(), lattice.getNy());
    computeDensity(lattice, *density);
    return std::unique_ptr<ScalarField2D<T> >(density);
}

/* *************** RhoBar ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeRhoBar(BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &rhoBar)
{
    applyProcessingFunctional(
        new BoxRhoBarFunctional2D<T, Descriptor>, lattice.getBoundingBox(), lattice, rhoBar);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField2D<T> > computeRhoBar(BlockLattice2D<T, Descriptor> &lattice)
{
    ScalarField2D<T> *rhoBar = new ScalarField2D<T>(lattice.getNx(), lattice.getNy());
    computeRhoBar(lattice, *rhoBar);
    return std::unique_ptr<ScalarField2D<T> >(rhoBar);
}

/* *************** Kinetic Energy ************************************ */

template <typename T, template <typename U> class Descriptor>
void computeKineticEnergy(BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &energy)
{
    applyProcessingFunctional(
        new BoxKineticEnergyFunctional2D<T, Descriptor>, lattice.getBoundingBox(), lattice, energy);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField2D<T> > computeKineticEnergy(BlockLattice2D<T, Descriptor> &lattice)
{
    ScalarField2D<T> *energy = new ScalarField2D<T>(lattice.getNx(), lattice.getNy());
    computeKineticEnergy(lattice, *energy);
    return std::unique_ptr<ScalarField2D<T> >(energy);
}

/* *************** Velocity Norm ************************************* */

template <typename T, template <typename U> class Descriptor>
void computeVelocityNorm(BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &velocityNorm)
{
    applyProcessingFunctional(
        new BoxVelocityNormFunctional2D<T, Descriptor>, lattice.getBoundingBox(), lattice,
        velocityNorm);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField2D<T> > computeVelocityNorm(BlockLattice2D<T, Descriptor> &lattice)
{
    ScalarField2D<T> *velocityNorm = new ScalarField2D<T>(lattice.getNx(), lattice.getNy());
    computeVelocityNorm(lattice, *velocityNorm);
    return std::unique_ptr<ScalarField2D<T> >(velocityNorm);
}

/* *************** Velocity Component ******************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocityComponent(
    BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &velocityComponent, plint iComponent)
{
    applyProcessingFunctional(
        new BoxVelocityComponentFunctional2D<T, Descriptor>(iComponent), lattice.getBoundingBox(),
        lattice, velocityComponent);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField2D<T> > computeVelocityComponent(
    BlockLattice2D<T, Descriptor> &lattice, plint iComponent)
{
    ScalarField2D<T> *velocityComponent = new ScalarField2D<T>(lattice.getNx(), lattice.getNy());
    computeVelocityComponent(lattice, *velocityComponent, iComponent);
    return std::unique_ptr<ScalarField2D<T> >(velocityComponent);
}

/* *************** Velocity ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocity(
    BlockLattice2D<T, Descriptor> &lattice, TensorField2D<T, Descriptor<T>::d> &velocity)
{
    applyProcessingFunctional(
        new BoxVelocityFunctional2D<T, Descriptor>, lattice.getBoundingBox(), lattice, velocity);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<TensorField2D<T, Descriptor<T>::d> > computeVelocity(
    BlockLattice2D<T, Descriptor> &lattice)
{
    TensorField2D<T, Descriptor<T>::d> *velocity =
        new TensorField2D<T, Descriptor<T>::d>(lattice.getNx(), lattice.getNy());
    computeVelocity(lattice, *velocity);
    return std::unique_ptr<TensorField2D<T, Descriptor<T>::d> >(velocity);
}

/* *************** Pi Neq ********************************* */

template <typename T, template <typename U> class Descriptor>
void computePiNeq(
    BlockLattice2D<T, Descriptor> &lattice,
    TensorField2D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq)
{
    applyProcessingFunctional(
        new BoxPiNeqFunctional2D<T, Descriptor>, lattice.getBoundingBox(), lattice, PiNeq);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<TensorField2D<T, SymmetricTensor<T, Descriptor>::n> > computePiNeq(
    BlockLattice2D<T, Descriptor> &lattice)
{
    TensorField2D<T, SymmetricTensor<T, Descriptor>::n> *PiNeq =
        new TensorField2D<T, SymmetricTensor<T, Descriptor>::n>(lattice.getNx(), lattice.getNy());
    computePiNeq(lattice, *PiNeq);
    return std::unique_ptr<TensorField2D<T, SymmetricTensor<T, Descriptor>::n> >(PiNeq);
}

/* *************** Shear Stress ********************************* */

template <typename T, template <typename U> class Descriptor>
void computeShearStress(
    BlockLattice2D<T, Descriptor> &lattice,
    TensorField2D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq)
{
    applyProcessingFunctional(
        new BoxShearStressFunctional2D<T, Descriptor>, lattice.getBoundingBox(), lattice, PiNeq);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<TensorField2D<T, SymmetricTensor<T, Descriptor>::n> > computeShearStress(
    BlockLattice2D<T, Descriptor> &lattice)
{
    TensorField2D<T, SymmetricTensor<T, Descriptor>::n> *PiNeq =
        new TensorField2D<T, SymmetricTensor<T, Descriptor>::n>(lattice.getNx(), lattice.getNy());
    computeShearStress(lattice, *PiNeq);
    return std::unique_ptr<TensorField2D<T, SymmetricTensor<T, Descriptor>::n> >(PiNeq);
}

/* *************** Strain Rate from Stress *************************** */

template <typename T, template <typename U> class Descriptor>
void computeStrainRateFromStress(
    BlockLattice2D<T, Descriptor> &lattice, TensorField2D<T, SymmetricTensor<T, Descriptor>::n> &S)
{
    applyProcessingFunctional(
        new BoxStrainRateFromStressFunctional2D<T, Descriptor>, lattice.getBoundingBox(), lattice,
        S);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<TensorField2D<T, SymmetricTensor<T, Descriptor>::n> > computeStrainRateFromStress(
    BlockLattice2D<T, Descriptor> &lattice)
{
    TensorField2D<T, SymmetricTensor<T, Descriptor>::n> *S =
        new TensorField2D<T, SymmetricTensor<T, Descriptor>::n>(lattice.getNx(), lattice.getNy());
    computeStrainRateFromStress(lattice, *S);
    return std::unique_ptr<TensorField2D<T, SymmetricTensor<T, Descriptor>::n> >(S);
}

/* *************** Temperature ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeTemperature(BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &temperature)
{
    applyProcessingFunctional(
        new BoxTemperatureFunctional2D<T, Descriptor>, lattice.getBoundingBox(), lattice,
        temperature);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField2D<T> > computeTemperature(BlockLattice2D<T, Descriptor> &lattice)
{
    ScalarField2D<T> *temperature = new ScalarField2D<T>(lattice.getNx(), lattice.getNy());
    computeTemperature(lattice, *temperature);
    return std::unique_ptr<ScalarField2D<T> >(temperature);
}

/* *************** SoundSpeed ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeSoundSpeed(BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &soundSpeed)
{
    applyProcessingFunctional(
        new BoxSoundSpeedFunctional2D<T, Descriptor>, lattice.getBoundingBox(), lattice,
        soundSpeed);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField2D<T> > computeSoundSpeed(BlockLattice2D<T, Descriptor> &lattice)
{
    ScalarField2D<T> *soundSpeed = new ScalarField2D<T>(lattice.getNx(), lattice.getNy());
    computeSoundSpeed(lattice, *soundSpeed);
    return std::unique_ptr<ScalarField2D<T> >(soundSpeed);
}

/* *************** Population *************************************** */

template <typename T, template <typename U> class Descriptor>
void computePopulation(
    BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &population, plint iPop)
{
    applyProcessingFunctional(
        new BoxPopulationFunctional2D<T, Descriptor>(iPop), lattice.getBoundingBox(), lattice,
        population);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField2D<T> > computePopulation(
    BlockLattice2D<T, Descriptor> &lattice, plint iPop)
{
    ScalarField2D<T> *population = new ScalarField2D<T>(lattice.getNx(), lattice.getNy());
    computePopulation(lattice, *population, iPop);
    return std::unique_ptr<ScalarField2D<T> >(population);
}

template <typename T, template <typename U> class Descriptor>
void computeAllPopulations(
    BlockLattice2D<T, Descriptor> &lattice, TensorField2D<T, Descriptor<T>::q> &populations,
    Box2D domain)
{
    applyProcessingFunctional(
        new BoxAllPopulationsFunctional2D<T, Descriptor>(), domain, lattice, populations);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField2D<T> > computeAllPopulations(BlockLattice2D<T, Descriptor> &lattice)
{
    TensorField2D<T, Descriptor<T>::q> *populations =
        new TensorField2D<T, Descriptor<T>::q>(lattice);
    computeAllPopulations(lattice, *populations);
    return std::unique_ptr<TensorField2D<T, Descriptor<T>::q> >(populations);
}

template <typename T, template <typename U> class Descriptor>
void copyPopulations(
    BlockLattice2D<T, Descriptor> &latticeFrom, BlockLattice2D<T, Descriptor> &latticeTo,
    Box2D domain)
{
    applyProcessingFunctional(
        new CopyPopulationsFunctional2D<T, Descriptor>(), domain, latticeFrom, latticeTo);
}

template <typename T, template <typename U> class Descriptor>
void copyAll(
    BlockLattice2D<T, Descriptor> &latticeFrom, BlockLattice2D<T, Descriptor> &latticeTo,
    Box2D domain)
{
    applyProcessingFunctional(
        new LatticeCopyAllFunctional2D<T, Descriptor>(), domain, latticeFrom, latticeTo);
}

template <typename T, template <typename U> class Descriptor>
void copyRegenerate(
    BlockLattice2D<T, Descriptor> &latticeFrom, BlockLattice2D<T, Descriptor> &latticeTo,
    Box2D domain)
{
    applyProcessingFunctional(
        new LatticeRegenerateFunctional2D<T, Descriptor>(), domain, latticeFrom, latticeTo);
}

/* ******************************************************************* */
/* *************** PART II. Atomic-block wrappers: Scalar-Field ****** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template <typename T>
T computeSum(ScalarField2D<T> &scalarField, Box2D domain)
{
    BoxScalarSumFunctional2D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getSumScalar();
}

template <typename T>
T computeSum(ScalarField2D<T> &scalarField)
{
    return computeSum(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeAverage(ScalarField2D<T> &scalarField, Box2D domain)
{
    BoxScalarSumFunctional2D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getSumScalar() / (T)domain.nCells();
}

template <typename T>
T computeAverage(ScalarField2D<T> &scalarField)
{
    return computeAverage(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeAverage(ScalarField2D<T> &scalarField, ScalarField2D<int> &mask, int flag, Box2D domain)
{
    MaskedBoxScalarAverageFunctional2D<T> functional(flag);
    applyProcessingFunctional(functional, domain, scalarField, mask);
    return functional.getAverageScalar();
}

template <typename T>
T computeAverage(ScalarField2D<T> &scalarField, ScalarField2D<int> &mask, int flag)
{
    return computeAverage(scalarField, mask, flag, scalarField.getBoundingBox());
}

template <typename T>
T computeMin(ScalarField2D<T> &scalarField, Box2D domain)
{
    BoxScalarMinFunctional2D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getMinScalar();
}

template <typename T>
T computeMin(ScalarField2D<T> &scalarField)
{
    return computeMin(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeMax(ScalarField2D<T> &scalarField, Box2D domain)
{
    BoxScalarMaxFunctional2D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getMaxScalar();
}

template <typename T>
T computeMax(ScalarField2D<T> &scalarField)
{
    return computeMax(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeBoundedSum(ScalarField2D<T> &scalarField, Box2D domain)
{
    BoundedBoxScalarSumFunctional2D<T> functional;
    plint envelopeWidth = 1;
    applyProcessingFunctional(functional, domain, scalarField, envelopeWidth);
    return functional.getSumScalar();
}

template <typename T>
T computeBoundedSum(ScalarField2D<T> &scalarField)
{
    return computeBoundedSum(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeBoundedAverage(ScalarField2D<T> &scalarField, Box2D domain)
{
    BoundedBoxScalarSumFunctional2D<T> functional;
    plint envelopeWidth = 1;
    applyProcessingFunctional(functional, domain, scalarField, envelopeWidth);
    return functional.getSumScalar() / (T)((domain.getNx() - 1) * (domain.getNy() - 1));
}

template <typename T>
T computeBoundedAverage(ScalarField2D<T> &scalarField)
{
    return computeBoundedAverage(scalarField, scalarField.getBoundingBox());
}

template <typename T, class BoolMask>
plint count(ScalarField2D<T> &field, Box2D domain, BoolMask boolMask)
{
    CountScalarElementsFunctional2D<T, BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, field);
    return functional.getCount();
}

template <typename T, class BoolMask>
plint count(ScalarField2D<T> &field, BoolMask boolMask)
{
    return count(field, field.getBoundingBox(), boolMask);
}

/* *************** ScalarField - Scalar operations *************** */

template <typename T>
void add(ScalarField2D<T> &field, T scalar, ScalarField2D<T> &result)
{
    applyProcessingFunctional(
        new A_plus_alpha_functional2D<T>(scalar), field.getBoundingBox(), field, result);
}

template <typename T>
std::unique_ptr<ScalarField2D<T> > add(ScalarField2D<T> &field, T scalar)
{
    ScalarField2D<T> *result = new ScalarField2D<T>(field.getNx(), field.getNy());
    add(field, scalar, *result);
    return std::unique_ptr<ScalarField2D<T> >(result);
}

template <typename T>
void add(T scalar, ScalarField2D<T> &field, ScalarField2D<T> &result)
{
    applyProcessingFunctional(
        new A_plus_alpha_functional2D<T>(scalar), field.getBoundingBox(), field, result);
}

template <typename T>
std::unique_ptr<ScalarField2D<T> > add(T scalar, ScalarField2D<T> &field)
{
    ScalarField2D<T> *result = new ScalarField2D<T>(field.getNx(), field.getNy());
    add(scalar, field, *result);
    return std::unique_ptr<ScalarField2D<T> >(result);
}

template <typename T>
void subtract(ScalarField2D<T> &field, T scalar, ScalarField2D<T> &result)
{
    applyProcessingFunctional(
        new A_minus_alpha_functional2D<T>(scalar), field.getBoundingBox(), field, result);
}

template <typename T>
std::unique_ptr<ScalarField2D<T> > subtract(ScalarField2D<T> &field, T scalar)
{
    ScalarField2D<T> *result = new ScalarField2D<T>(field.getNx(), field.getNy());
    subtract(field, scalar, *result);
    return std::unique_ptr<ScalarField2D<T> >(result);
}

template <typename T>
void subtract(T scalar, ScalarField2D<T> &field, ScalarField2D<T> &result)
{
    applyProcessingFunctional(
        new Alpha_minus_A_functional2D<T>(scalar), field.getBoundingBox(), field, result);
}

template <typename T>
std::unique_ptr<ScalarField2D<T> > subtract(T scalar, ScalarField2D<T> &field)
{
    ScalarField2D<T> *result = new ScalarField2D<T>(field.getNx(), field.getNy());
    subtract(scalar, field, *result);
    return std::unique_ptr<ScalarField2D<T> >(result);
}

template <typename T>
void multiply(ScalarField2D<T> &field, T scalar, ScalarField2D<T> &result)
{
    applyProcessingFunctional(
        new A_times_alpha_functional2D<T>(scalar), field.getBoundingBox(), field, result);
}

template <typename T>
std::unique_ptr<ScalarField2D<T> > multiply(ScalarField2D<T> &field, T scalar)
{
    ScalarField2D<T> *result = new ScalarField2D<T>(field.getNx(), field.getNy());
    multiply(field, scalar, *result);
    return std::unique_ptr<ScalarField2D<T> >(result);
}

template <typename T>
void multiply(T scalar, ScalarField2D<T> &field, ScalarField2D<T> &result)
{
    applyProcessingFunctional(
        new A_times_alpha_functional2D<T>(scalar), field.getBoundingBox(), field, result);
}

template <typename T>
std::unique_ptr<ScalarField2D<T> > multiply(T scalar, ScalarField2D<T> &field)
{
    ScalarField2D<T> *result = new ScalarField2D<T>(field.getNx(), field.getNy());
    multiply(scalar, field, *result);
    return std::unique_ptr<ScalarField2D<T> >(result);
}

template <typename T>
void divide(ScalarField2D<T> &field, T scalar, ScalarField2D<T> &result)
{
    applyProcessingFunctional(
        new A_dividedBy_alpha_functional2D<T>(scalar), field.getBoundingBox(), field, result);
}

template <typename T>
std::unique_ptr<ScalarField2D<T> > divide(ScalarField2D<T> &field, T scalar)
{
    ScalarField2D<T> *result = new ScalarField2D<T>(field.getNx(), field.getNy());
    divide(field, scalar, *result);
    return std::unique_ptr<ScalarField2D<T> >(result);
}

template <typename T>
void divide(T scalar, ScalarField2D<T> &field, ScalarField2D<T> &result)
{
    applyProcessingFunctional(
        new Alpha_dividedBy_A_functional2D<T>(scalar), field.getBoundingBox(), field, result);
}

template <typename T>
std::unique_ptr<ScalarField2D<T> > divide(T scalar, ScalarField2D<T> &field)
{
    ScalarField2D<T> *result = new ScalarField2D<T>(field.getNx(), field.getNy());
    divide(scalar, field, *result);
    return std::unique_ptr<ScalarField2D<T> >(result);
}

/* *************** ScalarField operations *************** */

template <typename T>
void computeSqrt(ScalarField2D<T> &A, ScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(new ComputeScalarSqrtFunctional2D<T>, domain, A, result);
}

template <typename T>
std::unique_ptr<ScalarField2D<T> > computeSqrt(ScalarField2D<T> &A, Box2D domain)
{
    ScalarField2D<T> *result = new ScalarField2D<T>(A.getNx(), A.getNy());
    computeSqrt(A, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<ScalarField2D<T> > computeSqrt(ScalarField2D<T> &A)
{
    return computeSqrt(A, A.getBoundingBox());
}

template <typename T>
void computeAbsoluteValue(ScalarField2D<T> &A, ScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(new ComputeAbsoluteValueFunctional2D<T>, domain, A, result);
}

template <typename T>
std::unique_ptr<ScalarField2D<T> > computeAbsoluteValue(ScalarField2D<T> &A, Box2D domain)
{
    ScalarField2D<T> *result = new ScalarField2D<T>(A.getNx(), A.getNy());
    computeAbsoluteValue(A, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<ScalarField2D<T> > computeAbsoluteValue(ScalarField2D<T> &A)
{
    return computeAbsoluteValue(A, A.getBoundingBox());
}

/* *************** ScalarField - Scalar inplace operations *************** */

template <typename T>
void addInPlace(ScalarField2D<T> &field, T scalar)
{
    applyProcessingFunctional(
        new A_plus_alpha_inplace_functional2D<T>(scalar), field.getBoundingBox(), field);
}

template <typename T>
void subtractInPlace(ScalarField2D<T> &field, T scalar)
{
    applyProcessingFunctional(
        new A_minus_alpha_inplace_functional2D<T>(scalar), field.getBoundingBox(), field);
}

template <typename T>
void multiplyInPlace(ScalarField2D<T> &field, T scalar)
{
    applyProcessingFunctional(
        new A_times_alpha_inplace_functional2D<T>(scalar), field.getBoundingBox(), field);
}

template <typename T>
void divideInPlace(ScalarField2D<T> &field, T scalar)
{
    applyProcessingFunctional(
        new A_dividedBy_alpha_inplace_functional2D<T>(scalar), field.getBoundingBox(), field);
}

/* *************** ScalarField - ScalarField operations *************** */

template <typename T1, typename T2>
void copy(ScalarField2D<T1> &field, ScalarField2D<T2> &convertedField)
{
    applyProcessingFunctional(
        new CopyConvertScalarFunctional2D<T1, T2>, field.getBoundingBox(), field, convertedField);
}

template <typename T1, typename T2>
std::unique_ptr<ScalarField2D<T2> > copyConvert(ScalarField2D<T1> &field)
{
    ScalarField2D<T2> *convertedField = new ScalarField2D<T2>(field.getNx(), field.getNy());
    plb::copy(field, *convertedField);
    return std::unique_ptr<ScalarField2D<T2> >(convertedField);
}

template <typename T>
void add(ScalarField2D<T> &A, ScalarField2D<T> &B, ScalarField2D<T> &result)
{
    std::vector<ScalarField2D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_plus_B_functional2D<T>, A.getBoundingBox(), fields);
}

template <typename T>
std::unique_ptr<ScalarField2D<T> > add(ScalarField2D<T> &A, ScalarField2D<T> &B)
{
    ScalarField2D<T> *result = new ScalarField2D<T>(A.getNx(), A.getNy());
    add(A, B, *result);
    return std::unique_ptr<ScalarField2D<T> >(result);
}

template <typename T>
void subtract(ScalarField2D<T> &A, ScalarField2D<T> &B, ScalarField2D<T> &result)
{
    std::vector<ScalarField2D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_minus_B_functional2D<T>, A.getBoundingBox(), fields);
}

template <typename T>
std::unique_ptr<ScalarField2D<T> > subtract(ScalarField2D<T> &A, ScalarField2D<T> &B)
{
    ScalarField2D<T> *result = new ScalarField2D<T>(A.getNx(), A.getNy());
    subtract(A, B, *result);
    return std::unique_ptr<ScalarField2D<T> >(result);
}

template <typename T>
void multiply(ScalarField2D<T> &A, ScalarField2D<T> &B, ScalarField2D<T> &result)
{
    std::vector<ScalarField2D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_times_B_functional2D<T>, A.getBoundingBox(), fields);
}

template <typename T>
std::unique_ptr<ScalarField2D<T> > multiply(ScalarField2D<T> &A, ScalarField2D<T> &B)
{
    ScalarField2D<T> *result = new ScalarField2D<T>(A.getNx(), A.getNy());
    multiply(A, B, *result);
    return std::unique_ptr<ScalarField2D<T> >(result);
}

template <typename T>
void divide(ScalarField2D<T> &A, ScalarField2D<T> &B, ScalarField2D<T> &result)
{
    std::vector<ScalarField2D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_dividedBy_B_functional2D<T>, A.getBoundingBox(), fields);
}

template <typename T>
std::unique_ptr<ScalarField2D<T> > divide(ScalarField2D<T> &A, ScalarField2D<T> &B)
{
    ScalarField2D<T> *result = new ScalarField2D<T>(A.getNx(), A.getNy());
    divide(A, B, *result);
    return std::unique_ptr<ScalarField2D<T> >(result);
}

/* *************** ScalarField - ScalarField inplace operations *************** */

template <typename T>
void addInPlace(ScalarField2D<T> &A, ScalarField2D<T> &B)
{
    applyProcessingFunctional(new A_plus_B_inplace_functional2D<T>, A.getBoundingBox(), A, B);
}

template <typename T>
void subtractInPlace(ScalarField2D<T> &A, ScalarField2D<T> &B)
{
    applyProcessingFunctional(new A_minus_B_inplace_functional2D<T>, A.getBoundingBox(), A, B);
}

template <typename T>
void multiplyInPlace(ScalarField2D<T> &A, ScalarField2D<T> &B)
{
    applyProcessingFunctional(new A_times_B_inplace_functional2D<T>, A.getBoundingBox(), A, B);
}

template <typename T>
void divideInPlace(ScalarField2D<T> &A, ScalarField2D<T> &B)
{
    applyProcessingFunctional(new A_dividedBy_B_inplace_functional2D<T>, A.getBoundingBox(), A, B);
}

/* ******************************************************************* */
/* *************** PART III. Atomic-block wrappers: Tensor-field ***** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template <typename T, plint nDim, class BoolMask>
plint count(TensorField2D<T, nDim> &field, Box2D domain, BoolMask boolMask)
{
    CountTensorElementsFunctional2D<T, nDim, BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, field);
    return functional.getCount();
}

template <typename T, plint nDim, class BoolMask>
plint count(TensorField2D<T, nDim> &field, BoolMask boolMask)
{
    return count(field, field.getBoundingBox(), boolMask);
}

/* *************** Component (scalar-field) out of a tensor-field ****** */

template <typename T, int nDim>
void extractComponent(
    TensorField2D<T, nDim> &tensorField, ScalarField2D<T> &component, int iComponent)
{
    applyProcessingFunctional(
        new ExtractTensorComponentFunctional2D<T, nDim>(iComponent), tensorField.getBoundingBox(),
        component, tensorField);
}

template <typename T, int nDim>
std::unique_ptr<ScalarField2D<T> > extractComponent(
    TensorField2D<T, nDim> &tensorField, int iComponent)
{
    ScalarField2D<T> *component = new ScalarField2D<T>(tensorField.getNx(), tensorField.getNy());
    extractComponent(tensorField, *component, iComponent);
    return std::unique_ptr<ScalarField2D<T> >(component);
}

/* *************** Vector-norm of each cell in the field *************** */

template <typename T, int nDim>
void computeNorm(TensorField2D<T, nDim> &tensorField, ScalarField2D<T> &norm)
{
    applyProcessingFunctional(
        new ComputeNormFunctional2D<T, nDim>, tensorField.getBoundingBox(), norm, tensorField);
}

template <typename T, int nDim>
std::unique_ptr<ScalarField2D<T> > computeNorm(TensorField2D<T, nDim> &tensorField)
{
    ScalarField2D<T> *norm = new ScalarField2D<T>(tensorField.getNx(), tensorField.getNy());
    computeNorm(tensorField, *norm);
    return std::unique_ptr<ScalarField2D<T> >(norm);
}

/* *************** Sqrt operation on each component of each cell *************** */

template <typename T, int nDim>
void computeSqrt(TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &result, Box2D domain)
{
    applyProcessingFunctional(new ComputeTensorSqrtFunctional2D<T, nDim>, domain, A, result);
}

template <typename T, int nDim>
std::unique_ptr<TensorField2D<T, nDim> > computeSqrt(TensorField2D<T, nDim> &A, Box2D domain)
{
    TensorField2D<T, nDim> *result = new TensorField2D<T, nDim>(A.getNx(), A.getNy());
    computeSqrt(A, *result, domain);
    return result;
}

template <typename T, int nDim>
std::unique_ptr<TensorField2D<T, nDim> > computeSqrt(TensorField2D<T, nDim> &A)
{
    return computeSqrt(A, A.getBoundingBox());
}

/* *************** Squared vector-norm of each cell in the field ******** */

template <typename T, int nDim>
void computeNormSqr(TensorField2D<T, nDim> &tensorField, ScalarField2D<T> &normSqr)
{
    applyProcessingFunctional(
        new ComputeNormSqrFunctional2D<T, nDim>, tensorField.getBoundingBox(), normSqr,
        tensorField);
}

template <typename T, int nDim>
std::unique_ptr<ScalarField2D<T> > computeNormSqr(TensorField2D<T, nDim> &tensorField)
{
    ScalarField2D<T> *normSqr = new ScalarField2D<T>(tensorField.getNx(), tensorField.getNy());
    computeNormSqr(tensorField, *normSqr);
    return std::unique_ptr<ScalarField2D<T> >(normSqr);
}

/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template <typename T>
void computeSymmetricTensorNorm(TensorField2D<T, 3> &tensorField, ScalarField2D<T> &norm)
{
    applyProcessingFunctional(
        new ComputeSymmetricTensorNormFunctional2D<T>, tensorField.getBoundingBox(), norm,
        tensorField);
}

template <typename T>
std::unique_ptr<ScalarField2D<T> > computeSymmetricTensorNorm(TensorField2D<T, 3> &tensorField)
{
    ScalarField2D<T> *norm = new ScalarField2D<T>(tensorField.getNx(), tensorField.getNy());
    computeSymmetricTensorNorm(tensorField, *norm);
    return std::unique_ptr<ScalarField2D<T> >(norm);
}

/* *************** Squared Tensor-norm of each symmetric tensor of a field*/

template <typename T>
void computeSymmetricTensorNormSqr(TensorField2D<T, 3> &tensorField, ScalarField2D<T> &normSqr)
{
    applyProcessingFunctional(
        new ComputeSymmetricTensorNormSqrFunctional2D<T>, tensorField.getBoundingBox(), normSqr,
        tensorField);
}

template <typename T>
std::unique_ptr<ScalarField2D<T> > computeSymmetricTensorNormSqr(TensorField2D<T, 3> &tensorField)
{
    ScalarField2D<T> *normSqr = new ScalarField2D<T>(tensorField.getNx(), tensorField.getNy());
    computeSymmetricTensorNormSqr(tensorField, *normSqr);
    return std::unique_ptr<ScalarField2D<T> >(normSqr);
}

/* *************** Trace of each symmetric tensor of a field ************* */

template <typename T>
void computeSymmetricTensorTrace(TensorField2D<T, 3> &tensorField, ScalarField2D<T> &trace)
{
    applyProcessingFunctional(
        new ComputeSymmetricTensorTraceFunctional2D<T>, tensorField.getBoundingBox(), trace,
        tensorField);
}

template <typename T>
std::unique_ptr<ScalarField2D<T> > computeSymmetricTensorTrace(TensorField2D<T, 3> &tensorField)
{
    ScalarField2D<T> *trace = new ScalarField2D<T>(tensorField.getNx(), tensorField.getNy());
    computeSymmetricTensorTrace(tensorField, *trace);
    return std::unique_ptr<ScalarField2D<T> >(trace);
}

/* *************** Vorticity from Velocity field *********************** */

template <typename T>
void computeVorticity(TensorField2D<T, 2> &velocity, ScalarField2D<T> &vorticity)
{
    plint envelopeWidth = 1;
    applyProcessingFunctional(
        new BoxVorticityFunctional2D<T, 2>, velocity.getBoundingBox(), vorticity, velocity,
        envelopeWidth);
}

template <typename T>
std::unique_ptr<ScalarField2D<T> > computeVorticity(TensorField2D<T, 2> &velocity)
{
    ScalarField2D<T> *vorticity = new ScalarField2D<T>(velocity.getNx(), velocity.getNy());
    computeVorticity(velocity, *vorticity);
    return std::unique_ptr<ScalarField2D<T> >(vorticity);
}

/* *************** Vorticity, witout boundary treatment, from Velocity field  */

template <typename T>
void computeBulkVorticity(TensorField2D<T, 2> &velocity, ScalarField2D<T> &vorticity)
{
    applyProcessingFunctional(
        new BoxBulkVorticityFunctional2D<T, 2>, velocity.getBoundingBox(), vorticity, velocity);
}

template <typename T>
std::unique_ptr<ScalarField2D<T> > computeBulkVorticity(TensorField2D<T, 2> &velocity)
{
    ScalarField2D<T> *vorticity = new ScalarField2D<T>(velocity.getNx(), velocity.getNy());
    computeBulkVorticity(velocity, *vorticity);
    return std::unique_ptr<ScalarField2D<T> >(vorticity);
}

/* *************** Strain rate from Velocity field ********************* */

template <typename T>
void computeStrainRate(TensorField2D<T, 2> &velocity, TensorField2D<T, 3> &S)
{
    plint envelopeWidth = 1;
    applyProcessingFunctional(
        new BoxStrainRateFunctional2D<T, 2>, velocity.getBoundingBox(), velocity, S, envelopeWidth);
}

template <typename T>
std::unique_ptr<TensorField2D<T, 3> > computeStrainRate(TensorField2D<T, 2> &velocity)
{
    TensorField2D<T, 3> *S = new TensorField2D<T, 3>(velocity.getNx(), velocity.getNy());
    computeStrainRate(velocity, *S);
    return std::unique_ptr<TensorField2D<T, 3> >(S);
}

/* *************** Vorticity, witout boundary treatment, from Velocity field  */

template <typename T>
void computeBulkStrainRate(TensorField2D<T, 2> &velocity, TensorField2D<T, 3> &S)
{
    applyProcessingFunctional(
        new BoxBulkStrainRateFunctional2D<T, 2>, velocity.getBoundingBox(), velocity, S);
}

template <typename T>
std::unique_ptr<TensorField2D<T, 3> > computeBulkVorticity(TensorField2D<T, 2> &velocity)
{
    TensorField2D<T, 3> *S = new TensorField2D<T, 3>(velocity.getNx(), velocity.getNy());
    computeBulkStrainRate(velocity, *S);
    return std::unique_ptr<TensorField2D<T, 3> >(S);
}

/* *************** TensorField - TensorField operations *************** */

template <typename T1, typename T2, int nDim>
void copy(TensorField2D<T1, nDim> &field, TensorField2D<T2, nDim> &convertedField)
{
    applyProcessingFunctional(
        new CopyConvertTensorFunctional2D<T1, T2, nDim>, field.getBoundingBox(), field,
        convertedField);
}

template <typename T1, typename T2, int nDim>
std::unique_ptr<TensorField2D<T2, nDim> > copyConvert(TensorField2D<T1, nDim> &field)
{
    TensorField2D<T2, nDim> *convertedField =
        new TensorField2D<T2, nDim>(field.getNx(), field.getNy());
    plb::copy(field, *convertedField);
    return std::unique_ptr<TensorField2D<T2, nDim> >(convertedField);
}

template <typename T, int nDim>
void add(TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B, TensorField2D<T, nDim> &result)
{
    std::vector<TensorField2D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(
        new Tensor_A_plus_B_functional2D<T, nDim>, A.getBoundingBox(), fields);
}

template <typename T, int nDim>
std::unique_ptr<TensorField2D<T, nDim> > add(TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B)
{
    TensorField2D<T, nDim> *result = new TensorField2D<T, nDim>(A.getNx(), A.getNy());
    add(A, B, *result);
    return std::unique_ptr<TensorField2D<T, nDim> >(result);
}

template <typename T, int nDim>
void subtract(TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B, TensorField2D<T, nDim> &result)
{
    std::vector<TensorField2D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(
        new Tensor_A_minus_B_functional2D<T, nDim>, A.getBoundingBox(), fields);
}

template <typename T, int nDim>
std::unique_ptr<TensorField2D<T, nDim> > subtract(
    TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B)
{
    TensorField2D<T, nDim> *result = new TensorField2D<T, nDim>(A.getNx(), A.getNy());
    subtract(A, B, *result);
    return std::unique_ptr<TensorField2D<T, nDim> >(result);
}

template <typename T, int nDim>
void multiply(TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B, TensorField2D<T, nDim> &result)
{
    std::vector<TensorField2D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(
        new Tensor_A_times_B_functional2D<T, nDim>, A.getBoundingBox(), fields);
}

template <typename T, int nDim>
std::unique_ptr<TensorField2D<T, nDim> > multiply(
    TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B)
{
    TensorField2D<T, nDim> *result = new TensorField2D<T, nDim>(A.getNx(), A.getNy());
    multiply(A, B, *result);
    return std::unique_ptr<TensorField2D<T, nDim> >(result);
}

template <typename T, int nDim>
void multiply(TensorField2D<T, nDim> &field, T scalar, TensorField2D<T, nDim> &result, Box2D domain)
{
    applyProcessingFunctional(new A_times_alpha_functional2D<T>(scalar), domain, field, result);
}

template <typename T, int nDim>
std::unique_ptr<TensorField2D<T, nDim> > multiply(TensorField2D<T, nDim> &field, T scalar)
{
    TensorField2D<T, nDim> *result = new TensorField2D<T, nDim>(field.getNx(), field.getNy());
    multiply(field, scalar, *result, field.getBoundingBox());
    return std::unique_ptr<TensorField2D<T, nDim> >(result);
}

template <typename T, int nDim>
void multiply(T scalar, TensorField2D<T, nDim> &field, TensorField2D<T, nDim> &result)
{
    applyProcessingFunctional(
        new Tensor_A_times_alpha_functional2D<T, nDim>(scalar), field.getBoundingBox(), field,
        result);
}

template <typename T, int nDim>
std::unique_ptr<TensorField2D<T, nDim> > multiply(T scalar, TensorField2D<T, nDim> &field)
{
    TensorField2D<T, nDim> *result = new TensorField2D<T, nDim>(field.getNx(), field.getNy());
    multiply(scalar, field, *result);
    return std::unique_ptr<TensorField2D<T, nDim> >(result);
}

template <typename T, int nDim>
void divide(TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B, TensorField2D<T, nDim> &result)
{
    std::vector<TensorField2D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(
        new Tensor_A_dividedBy_B_functional2D<T, nDim>, A.getBoundingBox(), fields);
}

template <typename T, int nDim>
std::unique_ptr<TensorField2D<T, nDim> > divide(
    TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B)
{
    TensorField2D<T, nDim> *result = new TensorField2D<T, nDim>(A.getNx(), A.getNy());
    divide(A, B, *result);
    return std::unique_ptr<TensorField2D<T, nDim> >(result);
}

/* *************** TensorField - TensorField inplace operations *************** */

template <typename T, int nDim>
void addInPlace(TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B)
{
    applyProcessingFunctional(
        new Tensor_A_plus_B_inplace_functional2D<T, nDim>, A.getBoundingBox(), A, B);
}

template <typename T, int nDim>
void subtractInPlace(TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B)
{
    applyProcessingFunctional(
        new Tensor_A_minus_B_inplace_functional2D<T, nDim>, A.getBoundingBox(), A, B);
}

template <typename T, int nDim>
void multiplyInPlace(TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B)
{
    applyProcessingFunctional(
        new Tensor_A_times_B_inplace_functional2D<T, nDim>, A.getBoundingBox(), A, B);
}

template <typename T, int nDim>
void multiplyInPlace(TensorField2D<T, nDim> &A, T alpha)
{
    applyProcessingFunctional(
        new Tensor_A_times_alpha_inplace_functional2D<T, nDim>(alpha), A.getBoundingBox(), A);
}

template <typename T, int nDim>
void divideInPlace(TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B)
{
    applyProcessingFunctional(
        new Tensor_A_dividedBy_B_inplace_functional2D<T, nDim>, A.getBoundingBox(), A, B);
}

/* ******************************************************************* */
/* *************** PART IV : Multi-block wrappers: Block-Lattice ***** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template <typename T, template <typename U> class Descriptor>
T computeAverageDensity(MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    BoxSumRhoBarFunctional2D<T, Descriptor> functional;
    applyProcessingFunctional(functional, domain, lattice);
    return Descriptor<T>::fullRho(functional.getSumRhoBar() / (T)domain.nCells());
}

template <typename T, template <typename U> class Descriptor>
T computeAverageDensity(MultiBlockLattice2D<T, Descriptor> &lattice)
{
    return computeAverageDensity(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    BoxSumRhoBarFunctional2D<T, Descriptor> functional;
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getSumRhoBar() / (T)domain.nCells();
}

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(MultiBlockLattice2D<T, Descriptor> &lattice)
{
    return computeAverageRhoBar(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(MultiBlockLattice2D<T, Descriptor> &lattice, DotList2D dotList)
{
    DotSumRhoBarFunctional2D<T, Descriptor> functional;
    applyProcessingFunctional(functional, dotList, lattice);
    return functional.getSumRhoBar() / (T)dotList.getN();
}

template <typename T, template <typename U> class Descriptor>
T computeAverageEnergy(MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    BoxSumEnergyFunctional2D<T, Descriptor> functional;
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getSumEnergy() / (T)domain.nCells();
    ;
}

template <typename T, template <typename U> class Descriptor>
T computeAverageEnergy(MultiBlockLattice2D<T, Descriptor> &lattice)
{
    return computeAverageEnergy(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
T computeAverageVelocityComponent(
    MultiBlockLattice2D<T, Descriptor> &lattice, DotList2D dotList, plint iComponent)
{
    DotSumVelocityComponentFunctional2D<T, Descriptor> functional(iComponent);
    applyProcessingFunctional(functional, dotList, lattice);
    return functional.getSumVelocityComponent() / (T)dotList.getN();
}

template <typename T, template <typename U> class Descriptor, class BoolMask>
plint count(MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, BoolMask boolMask)
{
    CountLatticeElementsFunctional2D<T, Descriptor, BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getCount();
}

template <typename T, template <typename U> class Descriptor, class BoolMask>
plint count(MultiBlockLattice2D<T, Descriptor> &lattice, BoolMask boolMask)
{
    return count(lattice, lattice.getBoundingBox(), boolMask);
}

/* ******************************************************************* */
/* *************** PART V : Multi-block wrappers: Scalar-Field ******* */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template <typename T>
T computeSum(MultiScalarField2D<T> &scalarField, Box2D domain)
{
    BoxScalarSumFunctional2D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getSumScalar();
}

template <typename T>
T computeSum(MultiScalarField2D<T> &scalarField)
{
    return computeSum(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeAverage(MultiScalarField2D<T> &scalarField, Box2D domain)
{
    BoxScalarSumFunctional2D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getSumScalar() / (T)domain.nCells();
}

template <typename T>
T computeAverage(MultiScalarField2D<T> &scalarField)
{
    return computeAverage(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeAverage(
    MultiScalarField2D<T> &scalarField, MultiScalarField2D<int> &mask, int flag, Box2D domain)
{
    MaskedBoxScalarAverageFunctional2D<T> functional(flag);
    applyProcessingFunctional(functional, domain, scalarField, mask);
    return functional.getAverageScalar();
}

template <typename T>
T computeAverage(MultiScalarField2D<T> &scalarField, MultiScalarField2D<int> &mask, int flag)
{
    return computeAverage(scalarField, mask, flag, scalarField.getBoundingBox());
}

template <typename T>
T computeMin(MultiScalarField2D<T> &scalarField, Box2D domain)
{
    BoxScalarMinFunctional2D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getMinScalar();
}

template <typename T>
T computeMin(MultiScalarField2D<T> &scalarField)
{
    return computeMin(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeMin(
    MultiScalarField2D<T> &scalarField, MultiScalarField2D<int> &mask, int flag, Box2D domain)
{
    MaskedBoxScalarMinFunctional2D<T> functional(flag);
    applyProcessingFunctional(functional, domain, scalarField, mask);
    return functional.getMinScalar();
}

template <typename T>
T computeMin(MultiScalarField2D<T> &scalarField, MultiScalarField2D<int> &mask, int flag)
{
    return computeMin(scalarField, mask, flag, scalarField.getBoundingBox());
}

template <typename T>
T computeMax(MultiScalarField2D<T> &scalarField, Box2D domain)
{
    BoxScalarMaxFunctional2D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getMaxScalar();
}

template <typename T>
T computeMax(MultiScalarField2D<T> &scalarField)
{
    return computeMax(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeMax(
    MultiScalarField2D<T> &scalarField, MultiScalarField2D<int> &mask, int flag, Box2D domain)
{
    MaskedBoxScalarMaxFunctional2D<T> functional(flag);
    applyProcessingFunctional(functional, domain, scalarField, mask);
    return functional.getMaxScalar();
}

template <typename T>
T computeMax(MultiScalarField2D<T> &scalarField, MultiScalarField2D<int> &mask, int flag)
{
    return computeMax(scalarField, mask, flag, scalarField.getBoundingBox());
}

template <typename T>
T computeBoundedSum(MultiScalarField2D<T> &scalarField, Box2D domain)
{
    BoundedBoxScalarSumFunctional2D<T> functional;
    plint envelopeWidth = 1;
    applyProcessingFunctional(functional, domain, scalarField, envelopeWidth);
    return functional.getSumScalar();
}

template <typename T>
T computeBoundedSum(MultiScalarField2D<T> &scalarField)
{
    return computeBoundedSum(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeBoundedAverage(MultiScalarField2D<T> &scalarField, Box2D domain)
{
    BoundedBoxScalarSumFunctional2D<T> functional;
    plint envelopeWidth = 1;
    applyProcessingFunctional(functional, domain, scalarField, envelopeWidth);
    return functional.getSumScalar() / (T)((domain.getNx() - 1) * (domain.getNy() - 1));
}

template <typename T>
T computeBoundedAverage(MultiScalarField2D<T> &scalarField)
{
    return computeBoundedAverage(scalarField, scalarField.getBoundingBox());
}

template <typename T, class BoolMask>
plint count(MultiScalarField2D<T> &field, Box2D domain, BoolMask boolMask)
{
    CountScalarElementsFunctional2D<T, BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, field);
    return functional.getCount();
}

template <typename T, class BoolMask>
plint count(MultiScalarField2D<T> &field, BoolMask boolMask)
{
    return count(field, field.getBoundingBox(), boolMask);
}

/* *************** Extract Sub-Lattice ******************************* */

template <typename T, template <typename U> class Descriptor>
void extractSubDomain(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiBlockLattice2D<T, Descriptor> &extractedLattice, Box2D domain)
{
    applyProcessingFunctional(
        new LatticeRegenerateFunctional2D<T, Descriptor>, domain, lattice, extractedLattice);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > extractSubDomain(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > extractedLattice =
        generateMultiBlockLattice<T, Descriptor>(lattice, domain);
    extractSubDomain(lattice, *extractedLattice, domain);
    return extractedLattice;
}

/* *************** Density ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeDensity(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &density, Box2D domain)
{
    applyProcessingFunctional(new BoxDensityFunctional2D<T, Descriptor>, domain, lattice, density);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeDensity(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > density = generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the BoxDensityFunctional2D()
    // acts on both bulk and envelope, you would expect the envelope layer around the domain, on the
    // density multi-block, to be assigned some proper values too. By default, this is however not
    // what happens, because the physical space occupied by these envelopes does not intersect with
    // the domain "domain". We work around this issue by extending the domain. There's no problem if
    // the enlarged domain gets beyond the actual extent of the lattice, because Palabos handles
    // these situations properly.

    computeDensity(
        lattice, *density, domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return density;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeDensity(MultiBlockLattice2D<T, Descriptor> &lattice)
{
    return computeDensity(lattice, lattice.getBoundingBox());
}

/* *************** RhoBar ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeRhoBar(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &rhoBar, Box2D domain)
{
    applyProcessingFunctional(new BoxRhoBarFunctional2D<T, Descriptor>, domain, lattice, rhoBar);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeRhoBar(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > rhoBar = generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the BoxRhoBarFunctional2D() acts
    // on both bulk and envelope, you would expect the envelope layer around the domain, on the
    // rhoBar multi-block, to be assigned some proper values too. By default, this is however not
    // what happens, because the physical space occupied by these envelopes does not intersect with
    // the domain "domain". We work around this issue by extending the domain. There's no problem if
    // the enlarged domain gets beyond the actual extent of the lattice, because Palabos handles
    // these situations properly.

    computeRhoBar(
        lattice, *rhoBar, domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return rhoBar;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeRhoBar(MultiBlockLattice2D<T, Descriptor> &lattice)
{
    return computeRhoBar(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void computeRhoBarJ(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &rhoBar,
    MultiTensorField2D<T, 2> &j, Box2D domain)
{
    std::vector<MultiBlock2D *> fields;
    fields.push_back(&lattice);
    fields.push_back(&rhoBar);
    fields.push_back(&j);
    applyProcessingFunctional(new BoxRhoBarJfunctional2D<T, Descriptor>, domain, fields);
}

/* *************** Packed RhoBar J *********************************** */

template <typename T, template <typename U> class Descriptor>
void computePackedRhoBarJ(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiNTensorField2D<T> &rhoBarJ, Box2D domain)
{
    applyProcessingFunctional(
        new PackedRhoBarJfunctional2D<T, Descriptor>, domain, lattice, rhoBarJ);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiNTensorField2D<T> > computePackedRhoBarJ(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    plint nDim = 1 + Descriptor<T>::d;
    std::unique_ptr<MultiNTensorField2D<T> > rhoBarJ(
        generateMultiNTensorField<T>(lattice, domain, nDim));

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the PackedRhoBarJFunctional2D()
    // acts on both bulk and envelope, you would expect the envelope layer around the domain, on the
    // rhoBarJ multi-block, to be assigned some proper values too. By default, this is however not
    // what happens, because the physical space occupied by these envelopes does not intersect with
    // the domain "domain". We work around this issue by extending the domain. There's no problem if
    // the enlarged domain gets beyond the actual extent of the lattice, because Palabos handles
    // these situations properly.

    computePackedRhoBarJ(
        lattice, *rhoBarJ, domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return rhoBarJ;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computePackedRhoBarJ(
    MultiBlockLattice2D<T, Descriptor> &lattice)
{
    return computePackedRhoBarJ(lattice, lattice.getBoundingBox());
}

/* *************** Kinetic Energy ************************************ */

template <typename T, template <typename U> class Descriptor>
void computeKineticEnergy(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &energy, Box2D domain)
{
    applyProcessingFunctional(
        new BoxKineticEnergyFunctional2D<T, Descriptor>, domain, lattice, energy);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeKineticEnergy(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > energy = generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxKineticEnergyFunctional2D() acts on both bulk and envelope, you would expect the envelope
    // layer around the domain, on the energy multi-block, to be assigned some proper values too. By
    // default, this is however not what happens, because the physical space occupied by these
    // envelopes does not intersect with the domain "domain". We work around this issue by extending
    // the domain. There's no problem if the enlarged domain gets beyond the actual extent of the
    // lattice, because Palabos handles these situations properly.

    computeKineticEnergy(
        lattice, *energy, domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return energy;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeKineticEnergy(
    MultiBlockLattice2D<T, Descriptor> &lattice)
{
    return computeKineticEnergy(lattice, lattice.getBoundingBox());
}

/* *************** Velocity Norm ************************************* */

template <typename T, template <typename U> class Descriptor>
void computeVelocityNorm(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &velocityNorm, Box2D domain)
{
    applyProcessingFunctional(
        new BoxVelocityNormFunctional2D<T, Descriptor>, domain, lattice, velocityNorm);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeVelocityNorm(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > velocityNorm =
        generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxVelocityNormFunctional2D() acts on both bulk and envelope, you would expect the envelope
    // layer around the domain, on the velocityNorm multi-block, to be assigned some proper values
    // too. By default, this is however not what happens, because the physical space occupied by
    // these envelopes does not intersect with the domain "domain". We work around this issue by
    // extending the domain. There's no problem if the enlarged domain gets beyond the actual extent
    // of the lattice, because Palabos handles these situations properly.

    computeVelocityNorm(
        lattice, *velocityNorm,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return velocityNorm;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeVelocityNorm(
    MultiBlockLattice2D<T, Descriptor> &lattice)
{
    return computeVelocityNorm(lattice, lattice.getBoundingBox());
}

/* *************** Velocity Component ******************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocityComponent(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &velocityComponent,
    Box2D domain, plint iComponent)
{
    applyProcessingFunctional(
        new BoxVelocityComponentFunctional2D<T, Descriptor>(iComponent), domain, lattice,
        velocityComponent);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeVelocityComponent(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, plint iComponent)
{
    std::unique_ptr<MultiScalarField2D<T> > velocityComponent =
        generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxVelocityComponentFunctional2D() acts on both bulk and envelope, you would expect the
    // envelope layer around the domain, on the velocityComponent multi-block, to be assigned some
    // proper values too. By default, this is however not what happens, because the physical space
    // occupied by these envelopes does not intersect with the domain "domain". We work around this
    // issue by extending the domain. There's no problem if the enlarged domain gets beyond the
    // actual extent of the lattice, because Palabos handles these situations properly.

    computeVelocityComponent(
        lattice, *velocityComponent,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()), iComponent);
    return velocityComponent;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeVelocityComponent(
    MultiBlockLattice2D<T, Descriptor> &lattice, plint iComponent)
{
    return computeVelocityComponent(lattice, lattice.getBoundingBox(), iComponent);
}

/* *************** Velocity ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocity(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiTensorField2D<T, Descriptor<T>::d> &velocity,
    Box2D domain)
{
    applyProcessingFunctional(
        new BoxVelocityFunctional2D<T, Descriptor>, domain, lattice, velocity);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, Descriptor<T>::d> > computeVelocity(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    std::unique_ptr<MultiTensorField2D<T, Descriptor<T>::d> > velocity =
        generateMultiTensorField<T, Descriptor<T>::d>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the BoxVelocityFunctional2D()
    // acts on both bulk and envelope, you would expect the envelope layer around the domain, on the
    // velocity multi-block, to be assigned some proper values too. By default, this is however not
    // what happens, because the physical space occupied by these envelopes does not intersect with
    // the domain "domain". We work around this issue by extending the domain. There's no problem if
    // the enlarged domain gets beyond the actual extent of the lattice, because Palabos handles
    // these situations properly.

    computeVelocity(
        lattice, *velocity, domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return velocity;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, Descriptor<T>::d> > computeVelocity(
    MultiBlockLattice2D<T, Descriptor> &lattice)
{
    return computeVelocity(lattice, lattice.getBoundingBox());
}

/* *************** Pi Neq ********************************* */

template <typename T, template <typename U> class Descriptor>
void computePiNeq(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq, Box2D domain)
{
    applyProcessingFunctional(new BoxPiNeqFunctional2D<T, Descriptor>, domain, lattice, PiNeq);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, SymmetricTensor<T, Descriptor>::n> > computePiNeq(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    std::unique_ptr<MultiTensorField2D<T, SymmetricTensor<T, Descriptor>::n> > PiNeq =
        generateMultiTensorField<T, SymmetricTensor<T, Descriptor>::n>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the BoxPiNeqFunctional2D() acts
    // on both bulk and envelope, you would expect the envelope layer around the domain, on the
    // PiNeq multi-block, to be assigned some proper values too. By default, this is however not
    // what happens, because the physical space occupied by these envelopes does not intersect with
    // the domain "domain". We work around this issue by extending the domain. There's no problem if
    // the enlarged domain gets beyond the actual extent of the lattice, because Palabos handles
    // these situations properly.

    computePiNeq(
        lattice, *PiNeq, domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return PiNeq;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, SymmetricTensor<T, Descriptor>::n> > computePiNeq(
    MultiBlockLattice2D<T, Descriptor> &lattice)
{
    return computePiNeq(lattice, lattice.getBoundingBox());
}

/* *************** Shear Stress ********************************* */

template <typename T, template <typename U> class Descriptor>
void computeShearStress(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq, Box2D domain)
{
    applyProcessingFunctional(
        new BoxShearStressFunctional2D<T, Descriptor>, domain, lattice, PiNeq);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, SymmetricTensor<T, Descriptor>::n> > computeShearStress(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    std::unique_ptr<MultiTensorField2D<T, SymmetricTensor<T, Descriptor>::n> > PiNeq =
        generateMultiTensorField<T, SymmetricTensor<T, Descriptor>::n>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the BoxShearStressFunctional2D()
    // acts on both bulk and envelope, you would expect the envelope layer around the domain, on the
    // PiNeq multi-block, to be assigned some proper values too. By default, this is however not
    // what happens, because the physical space occupied by these envelopes does not intersect with
    // the domain "domain". We work around this issue by extending the domain. There's no problem if
    // the enlarged domain gets beyond the actual extent of the lattice, because Palabos handles
    // these situations properly.

    computeShearStress(
        lattice, *PiNeq, domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return PiNeq;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, SymmetricTensor<T, Descriptor>::n> > computeShearStress(
    MultiBlockLattice2D<T, Descriptor> &lattice)
{
    return computeShearStress(lattice, lattice.getBoundingBox());
}

/* *************** Strain Rate from Stress *************************** */

template <typename T, template <typename U> class Descriptor>
void computeStrainRateFromStress(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, SymmetricTensor<T, Descriptor>::n> &S, Box2D domain)
{
    applyProcessingFunctional(
        new BoxStrainRateFromStressFunctional2D<T, Descriptor>, domain, lattice, S);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, SymmetricTensor<T, Descriptor>::n> >
    computeStrainRateFromStress(MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    std::unique_ptr<MultiTensorField2D<T, SymmetricTensor<T, Descriptor>::n> > S =
        generateMultiTensorField<T, SymmetricTensor<T, Descriptor>::n>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxStrainRateFromStressFunctional2D() acts on both bulk and envelope, you would expect the
    // envelope layer around the domain, on the S multi-block, to be assigned some proper values
    // too. By default, this is however not what happens, because the physical space occupied by
    // these envelopes does not intersect with the domain "domain". We work around this issue by
    // extending the domain. There's no problem if the enlarged domain gets beyond the actual extent
    // of the lattice, because Palabos handles these situations properly.

    computeStrainRateFromStress(
        lattice, *S, domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return S;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, SymmetricTensor<T, Descriptor>::n> >
    computeStrainRateFromStress(MultiBlockLattice2D<T, Descriptor> &lattice)
{
    return computeStrainRateFromStress(lattice, lattice.getBoundingBox());
}

/* *************** Temperature ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeTemperature(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &temperature, Box2D domain)
{
    applyProcessingFunctional(
        new BoxTemperatureFunctional2D<T, Descriptor>, domain, lattice, temperature);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeTemperature(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > temperature =
        generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the BoxTemperatureFunctional2D()
    // acts on both bulk and envelope, you would expect the envelope layer around the domain, on the
    // temperature multi-block, to be assigned some proper values too. By default, this is however
    // not what happens, because the physical space occupied by these envelopes does not intersect
    // with the domain "domain". We work around this issue by extending the domain. There's no
    // problem if the enlarged domain gets beyond the actual extent of the lattice, because Palabos
    // handles these situations properly.

    computeTemperature(
        lattice, *temperature,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return temperature;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeTemperature(
    MultiBlockLattice2D<T, Descriptor> &lattice)
{
    return computeTemperature(lattice, lattice.getBoundingBox());
}

/* *************** SoundSpeed ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeSoundSpeed(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &soundSpeed, Box2D domain)
{
    applyProcessingFunctional(
        new BoxSoundSpeedFunctional2D<T, Descriptor>, domain, lattice, soundSpeed);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeSoundSpeed(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > soundSpeed =
        generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the BoxSoundSpeedFunctional2D()
    // acts on both bulk and envelope, you would expect the envelope layer around the domain, on the
    // soundSpeed multi-block, to be assigned some proper values too. By default, this is however
    // not what happens, because the physical space occupied by these envelopes does not intersect
    // with the domain "domain". We work around this issue by extending the domain. There's no
    // problem if the enlarged domain gets beyond the actual extent of the lattice, because Palabos
    // handles these situations properly.

    computeSoundSpeed(
        lattice, *soundSpeed, domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return soundSpeed;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeSoundSpeed(
    MultiBlockLattice2D<T, Descriptor> &lattice)
{
    return computeSoundSpeed(lattice, lattice.getBoundingBox());
}

/* *************** Population **************************************** */

template <typename T, template <typename U> class Descriptor>
void computePopulation(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &population, Box2D domain,
    plint iPop)
{
    applyProcessingFunctional(
        new BoxPopulationFunctional2D<T, Descriptor>(iPop), domain, lattice, population);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computePopulation(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, plint iPop)
{
    std::unique_ptr<MultiScalarField2D<T> > population =
        generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the BoxPopulationFunctional2D()
    // acts on both bulk and envelope, you would expect the envelope layer around the domain, on the
    // population multi-block, to be assigned some proper values too. By default, this is however
    // not what happens, because the physical space occupied by these envelopes does not intersect
    // with the domain "domain". We work around this issue by extending the domain. There's no
    // problem if the enlarged domain gets beyond the actual extent of the lattice, because Palabos
    // handles these situations properly.

    computePopulation(
        lattice, *population, domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()),
        iPop);
    return population;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computePopulation(
    MultiBlockLattice2D<T, Descriptor> &lattice, plint iPop)
{
    return computePopulation(lattice, lattice.getBoundingBox(), iPop);
}

template <typename T, template <typename U> class Descriptor>
void computeEquilibrium(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &equilibrium, Box2D domain,
    plint iPop)
{
    applyProcessingFunctional(
        new BoxEquilibriumFunctional2D<T, Descriptor>(iPop), domain, lattice, equilibrium);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeEquilibrium(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, plint iPop)
{
    std::unique_ptr<MultiScalarField2D<T> > equilibrium =
        generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the BoxEquilibriumFunctional2D()
    // acts on both bulk and envelope, you would expect the envelope layer around the domain, on the
    // equilibrium multi-block, to be assigned some proper values too. By default, this is however
    // not what happens, because the physical space occupied by these envelopes does not intersect
    // with the domain "domain". We work around this issue by extending the domain. There's no
    // problem if the enlarged domain gets beyond the actual extent of the lattice, because Palabos
    // handles these situations properly.

    computeEquilibrium(
        lattice, *equilibrium, domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()),
        iPop);
    return equilibrium;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeEquilibrium(
    MultiBlockLattice2D<T, Descriptor> &lattice, plint iPop)
{
    return computeEquilibrium(lattice, lattice.getBoundingBox(), iPop);
}

template <typename T, template <typename U> class Descriptor>
void computeAllPopulations(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::q> &populations, Box2D domain)
{
    applyProcessingFunctional(
        new BoxAllPopulationsFunctional2D<T, Descriptor>(), domain, lattice, populations);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, Descriptor<T>::q> > computeAllPopulations(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    std::unique_ptr<MultiTensorField2D<T, Descriptor<T>::q> > populations =
        generateMultiTensorField<T, Descriptor<T>::q>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxAllPopulationsFunctional2D() acts on both bulk and envelope, you would expect the envelope
    // layer around the domain, on the populations multi-block, to be assigned some proper values
    // too. By default, this is however not what happens, because the physical space occupied by
    // these envelopes does not intersect with the domain "domain". We work around this issue by
    // extending the domain. There's no problem if the enlarged domain gets beyond the actual extent
    // of the lattice, because Palabos handles these situations properly.

    computeAllPopulations(
        lattice, *populations,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return populations;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, Descriptor<T>::q> > computeAllPopulations(
    MultiBlockLattice2D<T, Descriptor> &lattice)
{
    return computeAllPopulations(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void copyPopulations(
    MultiBlockLattice2D<T, Descriptor> &latticeFrom, MultiBlockLattice2D<T, Descriptor> &latticeTo,
    Box2D domain)
{
    applyProcessingFunctional(
        new CopyPopulationsFunctional2D<T, Descriptor>(), domain, latticeFrom, latticeTo);
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void copyConvertPopulations(
    MultiBlockLattice2D<T1, Descriptor1> &latticeFrom,
    MultiBlockLattice2D<T2, Descriptor2> &latticeTo, Box2D domain)
{
    applyProcessingFunctional(
        new CopyConvertPopulationsFunctional2D<T1, Descriptor1, T2, Descriptor2>(), domain,
        latticeFrom, latticeTo);
}

template <typename T, template <typename U> class Descriptor>
void copyAll(
    MultiBlockLattice2D<T, Descriptor> &latticeFrom, MultiBlockLattice2D<T, Descriptor> &latticeTo,
    Box2D domain)
{
    applyProcessingFunctional(
        new LatticeCopyAllFunctional2D<T, Descriptor>(), domain, latticeFrom, latticeTo);
}

template <typename T, template <typename U> class Descriptor>
void copyRegenerate(
    MultiBlockLattice2D<T, Descriptor> &latticeFrom, MultiBlockLattice2D<T, Descriptor> &latticeTo,
    Box2D domain)
{
    applyProcessingFunctional(
        new LatticeRegenerateFunctional2D<T, Descriptor>(), domain, latticeFrom, latticeTo);
}

/* *************** Omega ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeOmega(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &omega, Box2D domain)
{
    applyProcessingFunctional(new BoxOmegaFunctional2D<T, Descriptor>, domain, lattice, omega);
}

template <typename T, template <typename U> class Descriptor>
void computeOmega(MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &omega)
{
    computeOmega(lattice, omega, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeOmega(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > omega = generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the BoxOmegaFunctional2D() acts
    // on both bulk and envelope, you would expect the envelope layer around the domain, on the
    // omega multi-block, to be assigned some proper values too. By default, this is however not
    // what happens, because the physical space occupied by these envelopes does not intersect with
    // the domain "domain". We work around this issue by extending the domain. There's no problem if
    // the enlarged domain gets beyond the actual extent of the lattice, because Palabos handles
    // these situations properly.

    computeOmega(
        lattice, *omega, domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return omega;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeOmega(MultiBlockLattice2D<T, Descriptor> &lattice)
{
    return computeOmega(lattice, lattice.getBoundingBox());
}

/* *************** ExternalForce ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeExternalForce(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiTensorField2D<T, Descriptor<T>::d> &force,
    Box2D domain)
{
    applyProcessingFunctional(
        new BoxExternalForceFunctional2D<T, Descriptor>, domain, lattice, force);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, Descriptor<T>::d> > computeExternalForce(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    std::unique_ptr<MultiTensorField2D<T, Descriptor<T>::d> > force =
        generateMultiTensorField<T, Descriptor<T>::d>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxExternalForceFunctional2D() acts on both bulk and envelope, you would expect the envelope
    // layer around the domain, on the force multi-block, to be assigned some proper values too. By
    // default, this is however not what happens, because the physical space occupied by these
    // envelopes does not intersect with the domain "domain". We work around this issue by extending
    // the domain. There's no problem if the enlarged domain gets beyond the actual extent of the
    // lattice, because Palabos handles these situations properly.

    computeExternalForce(
        lattice, *force, domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return force;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, Descriptor<T>::d> > computeExternalForce(
    MultiBlockLattice2D<T, Descriptor> &lattice)
{
    return computeExternalForce(lattice, lattice.getBoundingBox());
}

/* *************** ExternalScalar ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeExternalScalar(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &scalar, int whichScalar,
    Box2D domain)
{
    applyProcessingFunctional(
        new BoxExternalScalarFunctional2D<T, Descriptor>(whichScalar), domain, lattice, scalar);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeExternalScalar(
    MultiBlockLattice2D<T, Descriptor> &lattice, int whichScalar, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > scalar = generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxExternalScalarFunctional2D() acts on both bulk and envelope, you would expect the envelope
    // layer around the domain, on the scalar multi-block, to be assigned some proper values too. By
    // default, this is however not what happens, because the physical space occupied by these
    // envelopes does not intersect with the domain "domain". We work around this issue by extending
    // the domain. There's no problem if the enlarged domain gets beyond the actual extent of the
    // lattice, because Palabos handles these situations properly.

    computeExternalScalar(
        lattice, *scalar, whichScalar,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return scalar;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<T> > computeExternalScalar(
    MultiBlockLattice2D<T, Descriptor> &lattice, int whichScalar)
{
    return computeExternalScalar(lattice, whichScalar, lattice.getBoundingBox());
}

/* *************** ExternalVector ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeExternalVector(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::d> &tensorField, int vectorBeginsAt, Box2D domain)
{
    applyProcessingFunctional(
        new BoxExternalVectorFunctional2D<T, Descriptor>(vectorBeginsAt), domain, lattice,
        tensorField);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, Descriptor<T>::d> > computeExternalVector(
    MultiBlockLattice2D<T, Descriptor> &lattice, int vectorBeginsAt, Box2D domain)
{
    std::unique_ptr<MultiTensorField2D<T, Descriptor<T>::d> > tensorField =
        generateMultiTensorField<T, Descriptor<T>::d>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxExternalVectorFunctional2D() acts on both bulk and envelope, you would expect the envelope
    // layer around the domain, on the tensor multi-block, to be assigned some proper values too. By
    // default, this is however not what happens, because the physical space occupied by these
    // envelopes does not intersect with the domain "domain". We work around this issue by extending
    // the domain. There's no problem if the enlarged domain gets beyond the actual extent of the
    // lattice, because Palabos handles these situations properly.

    computeExternalVector(
        lattice, *tensorField, vectorBeginsAt,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return tensorField;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField2D<T, Descriptor<T>::d> > computeExternalVector(
    MultiBlockLattice2D<T, Descriptor> &lattice, int vectorBeginsAt)
{
    return computeExternalVector(lattice, vectorBeginsAt, lattice.getBoundingBox());
}

/* *************** Extract Sub-ScalarField *************************** */

template <typename T>
void extractSubDomain(
    MultiScalarField2D<T> &field, MultiScalarField2D<T> &extractedField, Box2D domain)
{
    applyProcessingFunctional(
        new ExtractScalarSubDomainFunctional2D<T>, domain, field, extractedField);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > extractSubDomain(MultiScalarField2D<T> &field, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > extractedField =
        generateMultiScalarField<T>(field, domain);
    extractSubDomain(field, *extractedField, domain);
    return extractedField;
}

/* *************** MultiScalarField - Scalar operations *************** */

template <typename T>
void lessThan(MultiScalarField2D<T> &field, T scalar, MultiScalarField2D<int> &result, Box2D domain)
{
    applyProcessingFunctional(new A_lt_alpha_functional2D<T>(scalar), domain, field, result);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<int> > lessThan(
    MultiScalarField2D<T> &field, T scalar, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<int> > result = generateMultiScalarField<int>(field, domain);
    lessThan(field, scalar, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<int> > lessThan(MultiScalarField2D<T> &field, T scalar)
{
    return lessThan(field, scalar, field.getBoundingBox());
}

template <typename T>
void greaterThan(
    MultiScalarField2D<T> &field, T scalar, MultiScalarField2D<int> &result, Box2D domain)
{
    applyProcessingFunctional(new A_gt_alpha_functional2D<T>(scalar), domain, field, result);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<int> > greaterThan(
    MultiScalarField2D<T> &field, T scalar, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<int> > result = generateMultiScalarField<int>(field, domain);
    greaterThan(field, scalar, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<int> > greaterThan(MultiScalarField2D<T> &field, T scalar)
{
    return greaterThan(field, scalar, field.getBoundingBox());
}

template <typename T>
void add(MultiScalarField2D<T> &field, T scalar, MultiScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(new A_plus_alpha_functional2D<T>(scalar), domain, field, result);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > add(MultiScalarField2D<T> &field, T scalar, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > result = generateMultiScalarField<T>(field, domain);
    add(field, scalar, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > add(MultiScalarField2D<T> &field, T scalar)
{
    return add(field, scalar, field.getBoundingBox());
}

template <typename T>
void add(T scalar, MultiScalarField2D<T> &field, MultiScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(new A_plus_alpha_functional2D<T>(scalar), domain, field, result);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > add(T scalar, MultiScalarField2D<T> &field, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > result = generateMultiScalarField<T>(field, domain);
    add(scalar, field, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > add(T scalar, MultiScalarField2D<T> &field)
{
    return add(scalar, field, field.getBoundingBox());
}

template <typename T>
void subtract(MultiScalarField2D<T> &field, T scalar, MultiScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(new A_minus_alpha_functional2D<T>(scalar), domain, field, result);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > subtract(
    MultiScalarField2D<T> &field, T scalar, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > result = generateMultiScalarField<T>(field, domain);
    subtract(field, scalar, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > subtract(MultiScalarField2D<T> &field, T scalar)
{
    return subtract(field, scalar, field.getBoundingBox());
}

template <typename T>
void subtract(T scalar, MultiScalarField2D<T> &field, MultiScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(new Alpha_minus_A_functional2D<T>(scalar), domain, field, result);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > subtract(
    T scalar, MultiScalarField2D<T> &field, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > result = generateMultiScalarField<T>(field, domain);
    subtract(scalar, field, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > subtract(T scalar, MultiScalarField2D<T> &field)
{
    return subtract(scalar, field, field.getBoundingBox());
}

template <typename T>
void multiply(MultiScalarField2D<T> &field, T scalar, MultiScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(new A_times_alpha_functional2D<T>(scalar), domain, field, result);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > multiply(
    MultiScalarField2D<T> &field, T scalar, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > result = generateMultiScalarField<T>(field, domain);
    multiply(field, scalar, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > multiply(MultiScalarField2D<T> &field, T scalar)
{
    return multiply(field, scalar, field.getBoundingBox());
}

template <typename T>
void multiply(T scalar, MultiScalarField2D<T> &field, MultiScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(new A_times_alpha_functional2D<T>(scalar), domain, field, result);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > multiply(
    T scalar, MultiScalarField2D<T> &field, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > result = generateMultiScalarField<T>(field, domain);
    multiply(scalar, field, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > multiply(T scalar, MultiScalarField2D<T> &field)
{
    return multiply(scalar, field, field.getBoundingBox());
}

template <typename T>
void divide(MultiScalarField2D<T> &field, T scalar, MultiScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(new A_dividedBy_alpha_functional2D<T>(scalar), domain, field, result);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > divide(MultiScalarField2D<T> &field, T scalar, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > result = generateMultiScalarField<T>(field, domain);
    divide(field, scalar, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > divide(MultiScalarField2D<T> &field, T scalar)
{
    return divide(field, scalar, field.getBoundingBox());
}

template <typename T>
void divide(T scalar, MultiScalarField2D<T> &field, MultiScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(new Alpha_dividedBy_A_functional2D<T>(scalar), domain, field, result);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > divide(T scalar, MultiScalarField2D<T> &field, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > result = generateMultiScalarField<T>(field, domain);
    divide(scalar, field, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > divide(T scalar, MultiScalarField2D<T> &field)
{
    return divide(scalar, field, field.getBoundingBox());
}

/* *************** MultiScalarField - Scalar inplace operations *************** */

template <typename T>
void addInPlace(MultiScalarField2D<T> &field, T scalar, Box2D domain)
{
    applyProcessingFunctional(new A_plus_alpha_inplace_functional2D<T>(scalar), domain, field);
}

template <typename T>
void addInPlace(MultiScalarField2D<T> &field, T scalar)
{
    addInPlace(field, scalar, field.getBoundingBox());
}

template <typename T>
void subtractInPlace(MultiScalarField2D<T> &field, T scalar, Box2D domain)
{
    applyProcessingFunctional(new A_minus_alpha_inplace_functional2D<T>(scalar), domain, field);
}

template <typename T>
void subtractInPlace(MultiScalarField2D<T> &field, T scalar)
{
    subtractInPlace(field, scalar, field.getBoundingBox());
}

template <typename T>
void multiplyInPlace(MultiScalarField2D<T> &field, T scalar, Box2D domain)
{
    applyProcessingFunctional(new A_times_alpha_inplace_functional2D<T>(scalar), domain, field);
}

template <typename T>
void multiplyInPlace(MultiScalarField2D<T> &field, T scalar)
{
    multiplyInPlace(field, scalar, field.getBoundingBox());
}

template <typename T>
void divideInPlace(MultiScalarField2D<T> &field, T scalar, Box2D domain)
{
    applyProcessingFunctional(new A_dividedBy_alpha_inplace_functional2D<T>(scalar), domain, field);
}

template <typename T>
void divideInPlace(MultiScalarField2D<T> &field, T scalar)
{
    divideInPlace(field, scalar, field.getBoundingBox());
}

/* *************** MultiScalarField operations *************** */

template <typename T>
void computeSqrt(MultiScalarField2D<T> &A, MultiScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(new ComputeScalarSqrtFunctional2D<T>, domain, A, result);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeSqrt(MultiScalarField2D<T> &A, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > result = generateMultiScalarField<T>(A, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // ComputeScalarSqrtFunctional2D() acts on both bulk and envelope, you would expect the envelope
    // layer around the domain, on the result multi-block, to be assigned some proper values too. By
    // default, this is however not what happens, because the physical space occupied by these
    // envelopes does not intersect with the domain "domain". We work around this issue by extending
    // the domain. There's no problem if the enlarged domain gets beyond the actual extent of the
    // lattice, because Palabos handles these situations properly.

    computeSqrt(A, *result, domain.enlarge(A.getMultiBlockManagement().getEnvelopeWidth()));
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeSqrt(MultiScalarField2D<T> &A)
{
    return computeSqrt(A, A.getBoundingBox());
}

// =============  compute pow of each cell ========================= //

template <typename T>
void computePower(MultiScalarField2D<T> &A, MultiScalarField2D<T> &result, T power, Box2D domain)
{
    applyProcessingFunctional(new ComputeScalarPowFunctional2D<T>(power), domain, A, result);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computePower(
    MultiScalarField2D<T> &A, T power, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > result = generateMultiScalarField<T>(A, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // ComputeScalarSqrtFunctional2D() acts on both bulk and envelope, you would expect the envelope
    // layer around the domain, on the result multi-block, to be assigned some proper values too. By
    // default, this is however not what happens, because the physical space occupied by these
    // envelopes does not intersect with the domain "domain". We work around this issue by extending
    // the domain. There's no problem if the enlarged domain gets beyond the actual extent of the
    // lattice, because Palabos handles these situations properly.

    computePower(A, *result, power, domain.enlarge(A.getMultiBlockManagement().getEnvelopeWidth()));
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computePower(MultiScalarField2D<T> &A, T power)
{
    return computePower(A, power, A.getBoundingBox());
}

template <typename T>
void computeLog(MultiScalarField2D<T> &A, MultiScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(new ComputeScalarLogFunctional2D<T>, domain, A, result);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeLog(MultiScalarField2D<T> &A, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > result = generateMultiScalarField<T>(A, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // ComputeScalarLogFunctional2D() acts on both bulk and envelope, you would expect the envelope
    // layer around the domain, on the result multi-block, to be assigned some proper values too. By
    // default, this is however not what happens, because the physical space occupied by these
    // envelopes does not intersect with the domain "domain". We work around this issue by extending
    // the domain. There's no problem if the enlarged domain gets beyond the actual extent of the
    // lattice, because Palabos handles these situations properly.

    computeLog(A, *result, domain.enlarge(A.getMultiBlockManagement().getEnvelopeWidth()));
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeLog(MultiScalarField2D<T> &A)
{
    return computeLog(A, A.getBoundingBox());
}

template <typename T>
void computeAbsoluteValue(MultiScalarField2D<T> &A, MultiScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(new ComputeAbsoluteValueFunctional2D<T>, domain, A, result);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeAbsoluteValue(MultiScalarField2D<T> &A, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > result = generateMultiScalarField<T>(A, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // ComputeAbsoluteValueFunctional2D() acts on both bulk and envelope, you would expect the
    // envelope layer around the domain, on the result multi-block, to be assigned some proper
    // values too. By default, this is however not what happens, because the physical space occupied
    // by these envelopes does not intersect with the domain "domain". We work around this issue by
    // extending the domain. There's no problem if the enlarged domain gets beyond the actual extent
    // of the lattice, because Palabos handles these situations properly.

    computeAbsoluteValue(
        A, *result, domain.enlarge(A.getMultiBlockManagement().getEnvelopeWidth()));
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeAbsoluteValue(MultiScalarField2D<T> &A)
{
    return computeAbsoluteValue(A, A.getBoundingBox());
}

/* *************** MultiScalarField - MultiScalarField operations *************** */

template <typename T1, typename T2>
void copy(MultiScalarField2D<T1> &field, MultiScalarField2D<T2> &convertedField, Box2D domain)
{
    applyProcessingFunctional(
        new CopyConvertScalarFunctional2D<T1, T2>, domain, field, convertedField);
}

template <typename T1, typename T2>
std::unique_ptr<MultiScalarField2D<T2> > copyConvert(MultiScalarField2D<T1> &field, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T2> > convertedField =
        generateMultiScalarField<T2>(field, domain);
    plb::copy(field, *convertedField, domain);
    return convertedField;
}

template <typename T1, typename T2>
std::unique_ptr<MultiScalarField2D<T2> > copyConvert(MultiScalarField2D<T1> &field)
{
    return copyConvert<T1, T2>(field, field.getBoundingBox());
}

template <typename T>
void lessThan(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, MultiScalarField2D<int> &result,
    Box2D domain)
{
    std::vector<MultiBlock2D *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_lt_B_functional2D<T>, domain, fields);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<int> > lessThan(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<int> > result =
        generateIntersectMultiScalarField<int>(A, B, domain);
    lessThan(A, B, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<int> > lessThan(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B)
{
    return lessThan(A, B, A.getBoundingBox());
}

template <typename T>
void greaterThan(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, MultiScalarField2D<int> &result,
    Box2D domain)
{
    std::vector<MultiBlock2D *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_gt_B_functional2D<T>, domain, fields);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<int> > greaterThan(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<int> > result =
        generateIntersectMultiScalarField<int>(A, B, domain);
    greaterThan(A, B, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<int> > greaterThan(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B)
{
    return greaterThan(A, B, A.getBoundingBox());
}

template <typename T>
void add(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, MultiScalarField2D<T> &result, Box2D domain)
{
    std::vector<MultiScalarField2D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_plus_B_functional2D<T>, domain, fields);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > add(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > result =
        generateIntersectMultiScalarField<T>(A, B, domain);
    add(A, B, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > add(MultiScalarField2D<T> &A, MultiScalarField2D<T> &B)
{
    return add(A, B, A.getBoundingBox());
}

template <typename T>
void subtract(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, MultiScalarField2D<T> &result, Box2D domain)
{
    std::vector<MultiScalarField2D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_minus_B_functional2D<T>, domain, fields);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > subtract(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > result =
        generateIntersectMultiScalarField<T>(A, B, domain);
    subtract(A, B, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > subtract(MultiScalarField2D<T> &A, MultiScalarField2D<T> &B)
{
    return subtract(A, B, A.getBoundingBox());
}

template <typename T>
void multiply(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, MultiScalarField2D<T> &result, Box2D domain)
{
    std::vector<MultiScalarField2D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_times_B_functional2D<T>, domain, fields);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > multiply(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > result =
        generateIntersectMultiScalarField<T>(A, B, domain);
    multiply(A, B, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > multiply(MultiScalarField2D<T> &A, MultiScalarField2D<T> &B)
{
    return multiply(A, B, A.getBoundingBox());
}

template <typename T>
void divide(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, MultiScalarField2D<T> &result, Box2D domain)
{
    std::vector<MultiScalarField2D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_dividedBy_B_functional2D<T>, domain, fields);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > divide(
    MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > result =
        generateIntersectMultiScalarField<T>(A, B, domain);
    divide(A, B, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > divide(MultiScalarField2D<T> &A, MultiScalarField2D<T> &B)
{
    return divide(A, B, A.getBoundingBox());
}

/* *************** ScalarField - ScalarField inplace operations *************** */

template <typename T>
void addInPlace(MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, Box2D domain)
{
    applyProcessingFunctional(new A_plus_B_inplace_functional2D<T>, domain, A, B);
}

template <typename T>
void addInPlace(MultiScalarField2D<T> &A, MultiScalarField2D<T> &B)
{
    addInPlace(A, B, A.getBoundingBox());
}

template <typename T>
void subtractInPlace(MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, Box2D domain)
{
    applyProcessingFunctional(new A_minus_B_inplace_functional2D<T>, domain, A, B);
}

template <typename T>
void subtractInPlace(MultiScalarField2D<T> &A, MultiScalarField2D<T> &B)
{
    subtractInPlace(A, B, A.getBoundingBox());
}

template <typename T>
void multiplyInPlace(MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, Box2D domain)
{
    applyProcessingFunctional(new A_times_B_inplace_functional2D<T>, domain, A, B);
}

template <typename T>
void multiplyInPlace(MultiScalarField2D<T> &A, MultiScalarField2D<T> &B)
{
    multiplyInPlace(A, B, A.getBoundingBox());
}

template <typename T>
void divideInPlace(MultiScalarField2D<T> &A, MultiScalarField2D<T> &B, Box2D domain)
{
    applyProcessingFunctional(new A_dividedBy_B_inplace_functional2D<T>, domain, A, B);
}

template <typename T>
void divideInPlace(MultiScalarField2D<T> &A, MultiScalarField2D<T> &B)
{
    divideInPlace(A, B, A.getBoundingBox());
}

/* ******************************************************************* */
/* *************** PART VI : Multi-block wrappers: Tensor-field ****** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template <typename T, plint nDim, class BoolMask>
plint count(MultiTensorField2D<T, nDim> &field, Box2D domain, BoolMask boolMask)
{
    CountTensorElementsFunctional2D<T, nDim, BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, field);
    return functional.getCount();
}

template <typename T, plint nDim, class BoolMask>
plint count(MultiTensorField2D<T, nDim> &field, BoolMask boolMask)
{
    return count(field, field.getBoundingBox(), boolMask);
}

template <typename T, int nDim>
Array<T, nDim> computeAverage(MultiTensorField2D<T, nDim> &field, Box2D domain)
{
    BoxTensorSumFunctional2D<T, nDim> functional;
    applyProcessingFunctional(functional, domain, field);
    return functional.getSumTensor() / (T)domain.nCells();
}

template <typename T, int nDim>
Array<T, nDim> computeAverage(MultiTensorField2D<T, nDim> &field)
{
    return computeAverage(field, field.getBoundingBox());
}

/* *************** Extract Sub-TensorField *************************** */

template <typename T, int nDim>
void extractSubDomain(
    MultiTensorField2D<T, nDim> &field, MultiTensorField2D<T, nDim> &extractedField, Box2D domain)
{
    applyProcessingFunctional(
        new ExtractTensorSubDomainFunctional2D<T, nDim>, domain, field, extractedField);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > extractSubDomain(
    MultiTensorField2D<T, nDim> &field, Box2D domain)
{
    std::unique_ptr<MultiTensorField2D<T, nDim> > extractedField =
        generateMultiTensorField<T, nDim>(field, domain);
    extractSubDomain(field, *extractedField, domain);
    return extractedField;
}

/* *************** Component (scalar-field) out of a tensor-field ****** */

template <typename T, int nDim>
void extractComponent(
    MultiTensorField2D<T, nDim> &tensorField, MultiScalarField2D<T> &component, Box2D domain,
    int iComponent)
{
    applyProcessingFunctional(
        new ExtractTensorComponentFunctional2D<T, nDim>(iComponent), domain, component,
        tensorField);
}

template <typename T, int nDim>
std::unique_ptr<MultiScalarField2D<T> > extractComponent(
    MultiTensorField2D<T, nDim> &tensorField, Box2D domain, int iComponent)
{
    std::unique_ptr<MultiScalarField2D<T> > component =
        generateMultiScalarField<T>(tensorField, domain);
    extractComponent(tensorField, *component, domain, iComponent);
    return component;
}

template <typename T, int nDim>
std::unique_ptr<MultiScalarField2D<T> > extractComponent(
    MultiTensorField2D<T, nDim> &tensorField, int iComponent)
{
    return extractComponent(tensorField, tensorField.getBoundingBox(), iComponent);
}

/* *************** Vector-norm of each cell in the field *************** */

template <typename T, int nDim>
void computeNorm(
    MultiTensorField2D<T, nDim> &tensorField, MultiScalarField2D<T> &component, Box2D domain)
{
    applyProcessingFunctional(new ComputeNormFunctional2D<T, nDim>, domain, component, tensorField);
}

template <typename T, int nDim>
std::unique_ptr<MultiScalarField2D<T> > computeNorm(
    MultiTensorField2D<T, nDim> &tensorField, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > component =
        generateMultiScalarField<T>(tensorField, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the ComputeNormFunctional2D()
    // acts on both bulk and envelope, you would expect the envelope layer around the domain, on the
    // component multi-block, to be assigned some proper values too. By default, this is however not
    // what happens, because the physical space occupied by these envelopes does not intersect with
    // the domain "domain". We work around this issue by extending the domain. There's no problem if
    // the enlarged domain gets beyond the actual extent of the lattice, because Palabos handles
    // these situations properly.

    computeNorm(
        tensorField, *component,
        domain.enlarge(tensorField.getMultiBlockManagement().getEnvelopeWidth()));
    return component;
}

template <typename T, int nDim>
std::unique_ptr<MultiScalarField2D<T> > computeNorm(MultiTensorField2D<T, nDim> &tensorField)
{
    return computeNorm(tensorField, tensorField.getBoundingBox());
}

/* *************** Sqrt operation on each component of each cell *************** */

template <typename T, int nDim>
void computeSqrt(MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &result, Box2D domain)
{
    applyProcessingFunctional(new ComputeTensorSqrtFunctional2D<T, nDim>, domain, A, result);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > computeSqrt(
    MultiTensorField2D<T, nDim> &A, Box2D domain)
{
    std::unique_ptr<MultiTensorField2D<T, nDim> > result =
        generateMultiTensorField<T, nDim>(A, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // ComputeTensorSqrtFunctional2D() acts on both bulk and envelope, you would expect the envelope
    // layer around the domain, on the result multi-block, to be assigned some proper values too. By
    // default, this is however not what happens, because the physical space occupied by these
    // envelopes does not intersect with the domain "domain". We work around this issue by extending
    // the domain. There's no problem if the enlarged domain gets beyond the actual extent of the
    // lattice, because Palabos handles these situations properly.

    computeSqrt(A, *result, domain.enlarge(A.getMultiBlockManagement().getEnvelopeWidth()));
    return result;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > computeSqrt(MultiTensorField2D<T, nDim> &A)
{
    return computeSqrt(A, A.getBoundingBox());
}

/* *************** Squared vector-norm of each cell in the field ******** */

template <typename T, int nDim>
void computeNormSqr(
    MultiTensorField2D<T, nDim> &tensorField, MultiScalarField2D<T> &component, Box2D domain)
{
    applyProcessingFunctional(
        new ComputeNormSqrFunctional2D<T, nDim>, domain, component, tensorField);
}

template <typename T, int nDim>
std::unique_ptr<MultiScalarField2D<T> > computeNormSqr(
    MultiTensorField2D<T, nDim> &tensorField, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > component =
        generateMultiScalarField<T>(tensorField, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the ComputeNormSqrFunctional2D()
    // acts on both bulk and envelope, you would expect the envelope layer around the domain, on the
    // component multi-block, to be assigned some proper values too. By default, this is however not
    // what happens, because the physical space occupied by these envelopes does not intersect with
    // the domain "domain". We work around this issue by extending the domain. There's no problem if
    // the enlarged domain gets beyond the actual extent of the lattice, because Palabos handles
    // these situations properly.

    computeNormSqr(
        tensorField, *component,
        domain.enlarge(tensorField.getMultiBlockManagement().getEnvelopeWidth()));
    return component;
}

template <typename T, int nDim>
std::unique_ptr<MultiScalarField2D<T> > computeNormSqr(MultiTensorField2D<T, nDim> &tensorField)
{
    return computeNormSqr(tensorField, tensorField.getBoundingBox());
}

/* *************** Max element of each array of each cell *************** */

template <typename T, int nDim>
void computeMaximumElement(
    MultiTensorField2D<T, nDim> &A, MultiScalarField2D<T> &result, Box2D domain)
{
    applyProcessingFunctional(
        new BoxLocalMaximumPerComponentFunctional2D<T, nDim>(), domain, result, A);
}

template <typename T, int nDim>
std::unique_ptr<MultiScalarField2D<T> > computeMaximumElement(
    MultiTensorField2D<T, nDim> &A, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > result = generateMultiScalarField<T>(A, domain);

    computeMaximumElement(A, *result, domain);
    return result;
}

template <typename T, int nDim>
std::unique_ptr<MultiScalarField2D<T> > computeMaximumElement(MultiTensorField2D<T, nDim> &A)
{
    return computeMaximumElement(A, A.getBoundingBox());
}

/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template <typename T>
void computeSymmetricTensorNorm(
    MultiTensorField2D<T, 3> &tensorField, MultiScalarField2D<T> &norm, Box2D domain)
{
    applyProcessingFunctional(
        new ComputeSymmetricTensorNormFunctional2D<T>, domain, norm, tensorField);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeSymmetricTensorNorm(
    MultiTensorField2D<T, 3> &tensorField, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > norm = generateMultiScalarField<T>(tensorField, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // ComputesymmetricTensorNormFunctional2D() acts on both bulk and envelope, you would expect the
    // envelope layer around the domain, on the norm multi-block, to be assigned some proper values
    // too. By default, this is however not what happens, because the physical space occupied by
    // these envelopes does not intersect with the domain "domain". We work around this issue by
    // extending the domain. There's no problem if the enlarged domain gets beyond the actual extent
    // of the lattice, because Palabos handles these situations properly.

    computeSymmetricTensorNorm(
        tensorField, *norm,
        domain.enlarge(tensorField.getMultiBlockManagement().getEnvelopeWidth()));
    return norm;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeSymmetricTensorNorm(
    MultiTensorField2D<T, 3> &tensorField)
{
    return computeSymmetricTensorNorm(tensorField, tensorField.getBoundingBox());
}

/* *************** Squared Tensor-norm of each symmetric tensor of a field*/

template <typename T>
void computeSymmetricTensorNormSqr(
    MultiTensorField2D<T, 3> &tensorField, MultiScalarField2D<T> &normSqr, Box2D domain)
{
    applyProcessingFunctional(
        new ComputeSymmetricTensorNormSqrFunctional2D<T>, domain, normSqr, tensorField);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeSymmetricTensorNormSqr(
    MultiTensorField2D<T, 3> &tensorField, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > normSqr =
        generateMultiScalarField<T>(tensorField, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // ComputeSymmetricTensorNormSqrFunctional2D() acts on both bulk and envelope, you would expect
    // the envelope layer around the domain, on the normSqr multi-block, to be assigned some proper
    // values too. By default, this is however not what happens, because the physical space occupied
    // by these envelopes does not intersect with the domain "domain". We work around this issue by
    // extending the domain. There's no problem if the enlarged domain gets beyond the actual extent
    // of the lattice, because Palabos handles these situations properly.

    computeSymmetricTensorNormSqr(
        tensorField, *normSqr,
        domain.enlarge(tensorField.getMultiBlockManagement().getEnvelopeWidth()));
    return normSqr;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeSymmetricTensorNormSqr(
    MultiTensorField2D<T, 3> &tensorField)
{
    return computeSymmetricTensorNormSqr(tensorField, tensorField.getBoundingBox());
}

/* *************** Trace of each symmetric tensor of a field ************ */

template <typename T>
void computeSymmetricTensorTrace(
    MultiTensorField2D<T, 3> &tensorField, MultiScalarField2D<T> &trace, Box2D domain)
{
    applyProcessingFunctional(
        new ComputeSymmetricTensorTraceFunctional2D<T>, domain, trace, tensorField);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeSymmetricTensorTrace(
    MultiTensorField2D<T, 3> &tensorField, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > trace =
        generateMultiScalarField<T>(tensorField, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // ComputeSymmetricTensorTraceFunctional2D() acts on both bulk and envelope, you would expect
    // the envelope layer around the domain, on the trace multi-block, to be assigned some proper
    // values too. By default, this is however not what happens, because the physical space occupied
    // by these envelopes does not intersect with the domain "domain". We work around this issue by
    // extending the domain. There's no problem if the enlarged domain gets beyond the actual extent
    // of the lattice, because Palabos handles these situations properly.

    computeSymmetricTensorTrace(
        tensorField, *trace,
        domain.enlarge(tensorField.getMultiBlockManagement().getEnvelopeWidth()));
    return trace;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeSymmetricTensorTrace(
    MultiTensorField2D<T, 3> &tensorField)
{
    return computeSymmetricTensorTrace(tensorField, tensorField.getBoundingBox());
}

/* *************** Gradient from scalar field *********************** */

template <typename T>
void computeGradient(MultiScalarField2D<T> &phi, MultiTensorField2D<T, 2> &gradient, Box2D domain)
{
    plint envelopeWidth = 1;
    applyProcessingFunctional(new BoxGradientFunctional2D<T>, domain, phi, gradient, envelopeWidth);
}

template <typename T>
std::unique_ptr<MultiTensorField2D<T, 2> > computeGradient(MultiScalarField2D<T> &phi, Box2D domain)
{
    std::unique_ptr<MultiTensorField2D<T, 2> > gradient =
        generateMultiTensorField<T, 2>(phi, domain);
    computeGradient(phi, *gradient, domain);
    return gradient;
}

template <typename T>
std::unique_ptr<MultiTensorField2D<T, 2> > computeGradient(MultiScalarField2D<T> &phi)
{
    return computeGradient(phi, phi.getBoundingBox());
}

/* *************** Vorticity from Velocity field *********************** */

template <typename T>
void computeVorticity(
    MultiTensorField2D<T, 2> &velocity, MultiScalarField2D<T> &vorticity, Box2D domain)
{
    plint envelopeWidth = 1;
    applyProcessingFunctional(
        new BoxVorticityFunctional2D<T, 2>, domain, vorticity, velocity, envelopeWidth);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeVorticity(
    MultiTensorField2D<T, 2> &velocity, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > vorticity =
        generateMultiScalarField<T>(velocity, domain);
    computeVorticity(velocity, *vorticity, domain);
    return vorticity;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeVorticity(MultiTensorField2D<T, 2> &velocity)
{
    return computeVorticity(velocity, velocity.getBoundingBox());
}

/* *************** Vorticity, witout boundary treatment, from Velocity field  */

template <typename T>
void computeBulkVorticity(
    MultiTensorField2D<T, 2> &velocity, MultiScalarField2D<T> &vorticity, Box2D domain)
{
    applyProcessingFunctional(new BoxBulkVorticityFunctional2D<T, 2>, domain, vorticity, velocity);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeBulkVorticity(
    MultiTensorField2D<T, 2> &velocity, Box2D domain)
{
    std::unique_ptr<MultiScalarField2D<T> > vorticity =
        generateMultiScalarField<T>(velocity, domain);
    computeBulkVorticity(velocity, *vorticity, domain);
    return vorticity;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > computeBulkVorticity(MultiTensorField2D<T, 2> &velocity)
{
    return computeBulkVorticity(velocity, velocity.getBoundingBox());
}

/* *************** Strain Rate from Velocity field ********************* */

template <typename T>
void computeStrainRate(
    MultiTensorField2D<T, 2> &velocity, MultiTensorField2D<T, 3> &S, Box2D domain)
{
    plint envelopeWidth = 1;
    applyProcessingFunctional(
        new BoxStrainRateFunctional2D<T, 2>, domain, velocity, S, envelopeWidth);
}

template <typename T>
std::unique_ptr<MultiTensorField2D<T, 3> > computeStrainRate(
    MultiTensorField2D<T, 2> &velocity, Box2D domain)
{
    std::unique_ptr<MultiTensorField2D<T, 3> > S = generateMultiTensorField<T, 3>(velocity, domain);
    computeStrainRate(velocity, *S, domain);
    return S;
}

template <typename T>
std::unique_ptr<MultiTensorField2D<T, 3> > computeStrainRate(MultiTensorField2D<T, 2> &velocity)
{
    return computeStrainRate(velocity, velocity.getBoundingBox());
}

/* *************** Str. rate, witout boundary treatment, from Velocity field  */

template <typename T>
void computeBulkStrainRate(
    MultiTensorField2D<T, 2> &velocity, MultiTensorField2D<T, 3> &S, Box2D domain)
{
    applyProcessingFunctional(new BoxBulkStrainRateFunctional2D<T, 2>, domain, velocity, S);
}

template <typename T>
std::unique_ptr<MultiTensorField2D<T, 3> > computeBulkStrainRate(
    MultiTensorField2D<T, 2> &velocity, Box2D domain)
{
    std::unique_ptr<MultiTensorField2D<T, 3> > S = generateMultiTensorField<T, 3>(velocity, domain);
    computeBulkStrainRate(velocity, *S, domain);
    return S;
}

template <typename T>
std::unique_ptr<MultiTensorField2D<T, 3> > computeBulkStrainRate(MultiTensorField2D<T, 2> &velocity)
{
    return computeBulkStrainRate(velocity, velocity.getBoundingBox());
}

/* *************** MultiTensorField - MultiTensorField operations *************** */

template <typename T1, typename T2, int nDim>
void copy(
    MultiTensorField2D<T1, nDim> &field, MultiTensorField2D<T2, nDim> &convertedField, Box2D domain)
{
    applyProcessingFunctional(
        new CopyConvertTensorFunctional2D<T1, T2, nDim>, domain, field, convertedField);
}

template <typename T1, typename T2, int nDim>
std::unique_ptr<MultiTensorField2D<T2, nDim> > copyConvert(
    MultiTensorField2D<T1, nDim> &field, Box2D domain)
{
    std::unique_ptr<MultiTensorField2D<T2, nDim> > convertedField =
        generateMultiTensorField<T2, nDim>(field, domain);
    plb::copy<T1, T2, nDim>(field, *convertedField, domain);
    return convertedField;
}

template <typename T1, typename T2, int nDim>
std::unique_ptr<MultiTensorField2D<T2, nDim> > copyConvert(MultiTensorField2D<T1, nDim> &field)
{
    return copyConvert<T1, T2, nDim>(field, field.getBoundingBox());
}

template <typename T, int nDim>
void add(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B,
    MultiTensorField2D<T, nDim> &result, Box2D domain)
{
    std::vector<MultiTensorField2D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new Tensor_A_plus_B_functional2D<T, nDim>, domain, fields);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > add(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B, Box2D domain)
{
    std::unique_ptr<MultiTensorField2D<T, nDim> > result =
        generateIntersectMultiTensorField<T, nDim>(A, B, domain);
    add(A, B, *result, domain);
    return result;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > add(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B)
{
    return add(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void subtract(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B,
    MultiTensorField2D<T, nDim> &result, Box2D domain)
{
    std::vector<MultiTensorField2D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new Tensor_A_minus_B_functional2D<T, nDim>, domain, fields);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > subtract(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B, Box2D domain)
{
    std::unique_ptr<MultiTensorField2D<T, nDim> > result =
        generateIntersectMultiTensorField<T, nDim>(A, B, domain);
    subtract(A, B, *result, domain);
    return result;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > subtract(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B)
{
    return subtract(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void multiply(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B,
    MultiTensorField2D<T, nDim> &result, Box2D domain)
{
    std::vector<MultiTensorField2D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new Tensor_A_times_B_functional2D<T, nDim>, domain, fields);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > multiply(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B, Box2D domain)
{
    std::unique_ptr<MultiTensorField2D<T, nDim> > result =
        generateIntersectMultiTensorField<T, nDim>(A, B, domain);
    multiply(A, B, *result, domain);
    return result;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > multiply(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B)
{
    return multiply(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void divide(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B,
    MultiTensorField2D<T, nDim> &result, Box2D domain)
{
    std::vector<MultiTensorField2D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new Tensor_A_dividedBy_B_functional2D<T, nDim>, domain, fields);
}

template <typename T, int nDim>
void multiply(
    MultiTensorField2D<T, nDim> &field, T scalar, MultiTensorField2D<T, nDim> &result, Box2D domain)
{
    applyProcessingFunctional(
        new Tensor_A_times_alpha_functional2D<T, nDim>(scalar), domain, field, result);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > multiply(
    MultiTensorField2D<T, nDim> &field, T scalar, Box2D domain)
{
    std::unique_ptr<MultiTensorField2D<T, nDim> > result =
        generateMultiTensorField<T, nDim>(field, domain);
    multiply(field, scalar, *result, domain);
    return result;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > multiply(MultiTensorField2D<T, nDim> &field, T scalar)
{
    return multiply(field, scalar, field.getBoundingBox());
}

template <typename T, int nDim>
void multiply(
    T scalar, MultiTensorField2D<T, nDim> &field, MultiTensorField2D<T, nDim> &result, Box2D domain)
{
    applyProcessingFunctional(
        new Tensor_A_times_alpha_functional2D<T, nDim>(scalar), domain, field, result);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > multiply(
    T scalar, MultiTensorField2D<T, nDim> &field, Box2D domain)
{
    std::unique_ptr<MultiTensorField2D<T, nDim> > result =
        generateMultiTensorField<T, nDim>(field, domain);
    multiply(scalar, field, *result, domain);
    return result;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > multiply(T scalar, MultiTensorField2D<T, nDim> &field)
{
    return multiply(scalar, field, field.getBoundingBox());
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > divide(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B, Box2D domain)
{
    std::unique_ptr<MultiTensorField2D<T, nDim> > result =
        generateIntersectMultiTensorField<T, nDim>(A, B, domain);
    divide(A, B, *result, domain);
    return result;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > divide(
    MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B)
{
    return divide(A, B, A.getBoundingBox());
}

/* *************** TensorField - TensorField inplace operations *************** */

template <typename T, int nDim>
void addInPlace(MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B, Box2D domain)
{
    applyProcessingFunctional(new Tensor_A_plus_B_inplace_functional2D<T, nDim>, domain, A, B);
}

template <typename T, int nDim>
void addInPlace(MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B)
{
    addInPlace(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void subtractInPlace(MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B, Box2D domain)
{
    applyProcessingFunctional(new Tensor_A_minus_B_inplace_functional2D<T, nDim>, domain, A, B);
}

template <typename T, int nDim>
void subtractInPlace(MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B)
{
    subtractInPlace(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void multiplyInPlace(MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B, Box2D domain)
{
    applyProcessingFunctional(new Tensor_A_times_B_inplace_functional2D<T, nDim>, domain, A, B);
}

template <typename T, int nDim>
void multiplyInPlace(MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B)
{
    multiplyInPlace(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void multiplyInPlace(MultiTensorField2D<T, nDim> &A, T alpha, Box2D domain)
{
    applyProcessingFunctional(
        new Tensor_A_times_alpha_inplace_functional2D<T, nDim>(alpha), domain, A);
}

template <typename T, int nDim>
void multiplyInPlace(MultiTensorField2D<T, nDim> &A, T alpha)
{
    multiplyInPlace(A, alpha, A.getBoundingBox());
}

template <typename T, int nDim>
void divideInPlace(MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B, Box2D domain)
{
    applyProcessingFunctional(new Tensor_A_dividedBy_B_inplace_functional2D<T, nDim>, domain, A, B);
}

template <typename T, int nDim>
void divideInPlace(MultiTensorField2D<T, nDim> &A, MultiTensorField2D<T, nDim> &B)
{
    divideInPlace(A, B, A.getBoundingBox());
}

}  // namespace plb

#endif  // DATA_ANALYSIS_WRAPPER_2D_HH
