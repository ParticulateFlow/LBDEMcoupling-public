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

#ifndef DATA_ANALYSIS_WRAPPER_3D_HH
#define DATA_ANALYSIS_WRAPPER_3D_HH

#include "atomicBlock/dataProcessorWrapper3D.h"
#include "atomicBlock/reductiveDataProcessorWrapper3D.h"
#include "dataProcessors/dataAnalysisFunctional3D.h"
#include "dataProcessors/dataAnalysisWrapper3D.h"
#include "multiBlock/multiBlockGenerator3D.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper3D.h"

namespace plb {

/* ******************************************************************* */
/* *************** PART I. Atomic-block wrappers: Block-Lattice ****** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template <typename T, template <typename U> class Descriptor>
T computeAverageDensity(BlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    BoxSumRhoBarFunctional3D<T, Descriptor> functional;
    applyProcessingFunctional(functional, domain, lattice);
    return Descriptor<T>::fullRho(functional.getSumRhoBar() / (T)domain.nCells());
}

template <typename T, template <typename U> class Descriptor>
T computeAverageDensity(BlockLattice3D<T, Descriptor> &lattice)
{
    return computeAverageDensity(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(BlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    BoxSumRhoBarFunctional3D<T, Descriptor> functional;
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getSumRhoBar() / (T)domain.nCells();
}

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(BlockLattice3D<T, Descriptor> &lattice)
{
    return computeAverageRhoBar(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
T computeAverageEnergy(BlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    BoxSumEnergyFunctional3D<T, Descriptor> functional;
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getSumEnergy() / (T)domain.nCells();
    ;
}

template <typename T, template <typename U> class Descriptor>
T computeAverageEnergy(BlockLattice3D<T, Descriptor> &lattice)
{
    return computeAverageEnergy(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor, class BoolMask>
plint count(BlockLattice3D<T, Descriptor> &lattice, Box3D domain, BoolMask boolMask)
{
    CountLatticeElementsFunctional3D<T, Descriptor, BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getCount();
}

template <typename T, template <typename U> class Descriptor, class BoolMask>
plint count(BlockLattice3D<T, Descriptor> &lattice, BoolMask boolMask)
{
    return count(lattice, lattice.getBoundingBox(), boolMask);
}

/* *************** Density ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeDensity(BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &density)
{
    applyProcessingFunctional(
        new BoxDensityFunctional3D<T, Descriptor>, lattice.getBoundingBox(), lattice, density);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField3D<T> > computeDensity(BlockLattice3D<T, Descriptor> &lattice)
{
    ScalarField3D<T> *density =
        new ScalarField3D<T>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeDensity(lattice, *density);
    return std::unique_ptr<ScalarField3D<T> >(density);
}

/* *************** RhoBar ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeRhoBar(BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &rhoBar)
{
    applyProcessingFunctional(
        new BoxRhoBarFunctional3D<T, Descriptor>, lattice.getBoundingBox(), lattice, rhoBar);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField3D<T> > computeRhoBar(BlockLattice3D<T, Descriptor> &lattice)
{
    ScalarField3D<T> *rhoBar =
        new ScalarField3D<T>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeRhoBar(lattice, *rhoBar);
    return std::unique_ptr<ScalarField3D<T> >(rhoBar);
}

/* *************** Kinetic Energy ************************************ */

template <typename T, template <typename U> class Descriptor>
void computeKineticEnergy(BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &energy)
{
    applyProcessingFunctional(
        new BoxKineticEnergyFunctional3D<T, Descriptor>, lattice.getBoundingBox(), lattice, energy);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField3D<T> > computeKineticEnergy(BlockLattice3D<T, Descriptor> &lattice)
{
    ScalarField3D<T> *energy =
        new ScalarField3D<T>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeKineticEnergy(lattice, *energy);
    return std::unique_ptr<ScalarField3D<T> >(energy);
}

/* *************** Velocity Norm ************************************* */

template <typename T, template <typename U> class Descriptor>
void computeVelocityNorm(BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &velocityNorm)
{
    applyProcessingFunctional(
        new BoxVelocityNormFunctional3D<T, Descriptor>, lattice.getBoundingBox(), lattice,
        velocityNorm);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField3D<T> > computeVelocityNorm(BlockLattice3D<T, Descriptor> &lattice)
{
    ScalarField3D<T> *velocityNorm =
        new ScalarField3D<T>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeVelocityNorm(lattice, *velocityNorm);
    return std::unique_ptr<ScalarField3D<T> >(velocityNorm);
}

/* *************** Velocity Component ******************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocityComponent(
    BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &velocityComponent, plint iComponent)
{
    applyProcessingFunctional(
        new BoxVelocityComponentFunctional3D<T, Descriptor>(iComponent), lattice.getBoundingBox(),
        lattice, velocityComponent);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField3D<T> > computeVelocityComponent(
    BlockLattice3D<T, Descriptor> &lattice, plint iComponent)
{
    ScalarField3D<T> *velocityComponent =
        new ScalarField3D<T>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeVelocityComponent(lattice, *velocityComponent, iComponent);
    return std::unique_ptr<ScalarField3D<T> >(velocityComponent);
}

/* *************** Velocity ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocity(
    BlockLattice3D<T, Descriptor> &lattice, TensorField3D<T, Descriptor<T>::d> &velocity)
{
    applyProcessingFunctional(
        new BoxVelocityFunctional3D<T, Descriptor>, lattice.getBoundingBox(), lattice, velocity);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<TensorField3D<T, Descriptor<T>::d> > computeVelocity(
    BlockLattice3D<T, Descriptor> &lattice)
{
    TensorField3D<T, Descriptor<T>::d> *velocity =
        new TensorField3D<T, Descriptor<T>::d>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeVelocity(lattice, *velocity);
    return std::unique_ptr<TensorField3D<T, Descriptor<T>::d> >(velocity);
}

/* *************** Pi Neq ********************************* */

template <typename T, template <typename U> class Descriptor>
void computePiNeq(
    BlockLattice3D<T, Descriptor> &lattice,
    TensorField3D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq)
{
    applyProcessingFunctional(
        new BoxPiNeqFunctional3D<T, Descriptor>, lattice.getBoundingBox(), lattice, PiNeq);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<TensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computePiNeq(
    BlockLattice3D<T, Descriptor> &lattice)
{
    TensorField3D<T, SymmetricTensor<T, Descriptor>::n> *PiNeq =
        new TensorField3D<T, SymmetricTensor<T, Descriptor>::n>(
            lattice.getNx(), lattice.getNy(), lattice.getNz());
    computePiNeq(lattice, *PiNeq);
    return std::unique_ptr<TensorField3D<T, SymmetricTensor<T, Descriptor>::n> >(PiNeq);
}

/* *************** Shear Stress ********************************* */

template <typename T, template <typename U> class Descriptor>
void computeShearStress(
    BlockLattice3D<T, Descriptor> &lattice,
    TensorField3D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq)
{
    applyProcessingFunctional(
        new BoxShearStressFunctional3D<T, Descriptor>, lattice.getBoundingBox(), lattice, PiNeq);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<TensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computeShearStress(
    BlockLattice3D<T, Descriptor> &lattice)
{
    TensorField3D<T, SymmetricTensor<T, Descriptor>::n> *PiNeq =
        new TensorField3D<T, SymmetricTensor<T, Descriptor>::n>(
            lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeShearStress(lattice, *PiNeq);
    return std::unique_ptr<TensorField3D<T, SymmetricTensor<T, Descriptor>::n> >(PiNeq);
}

/* *************** Strain Rate from Stress *************************** */

template <typename T, template <typename U> class Descriptor>
void computeStrainRateFromStress(
    BlockLattice3D<T, Descriptor> &lattice, TensorField3D<T, SymmetricTensor<T, Descriptor>::n> &S)
{
    applyProcessingFunctional(
        new BoxStrainRateFromStressFunctional3D<T, Descriptor>, lattice.getBoundingBox(), lattice,
        S);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<TensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computeStrainRateFromStress(
    BlockLattice3D<T, Descriptor> &lattice)
{
    TensorField3D<T, SymmetricTensor<T, Descriptor>::n> *S =
        new TensorField3D<T, SymmetricTensor<T, Descriptor>::n>(
            lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeStrainRateFromStress(lattice, *S);
    return std::unique_ptr<TensorField3D<T, SymmetricTensor<T, Descriptor>::n> >(S);
}

/* *************** Population *************************************** */

template <typename T, template <typename U> class Descriptor>
void computePopulation(
    BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &population, plint iPop)
{
    applyProcessingFunctional(
        new BoxPopulationFunctional3D<T, Descriptor>(iPop), lattice.getBoundingBox(), lattice,
        population);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField3D<T> > computePopulation(
    BlockLattice3D<T, Descriptor> &lattice, plint iPop)
{
    ScalarField3D<T> *population =
        new ScalarField3D<T>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computePopulation(lattice, *population, iPop);
    return std::unique_ptr<ScalarField3D<T> >(population);
}

template <typename T, template <typename U> class Descriptor>
void computeAllPopulations(
    BlockLattice3D<T, Descriptor> &lattice, TensorField3D<T, Descriptor<T>::q> &populations,
    Box3D domain)
{
    applyProcessingFunctional(
        new BoxAllPopulationsFunctional3D<T, Descriptor>(), domain, lattice, populations);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<ScalarField3D<T> > computeAllPopulations(BlockLattice3D<T, Descriptor> &lattice)
{
    TensorField3D<T, Descriptor<T>::q> *populations =
        new TensorField3D<T, Descriptor<T>::q>(lattice);
    computeAllPopulations(lattice, *populations);
    return std::unique_ptr<TensorField3D<T, Descriptor<T>::q> >(populations);
}

template <typename T, template <typename U> class Descriptor>
void copyPopulations(
    BlockLattice3D<T, Descriptor> &latticeFrom, BlockLattice3D<T, Descriptor> &latticeTo,
    Box3D domain)
{
    applyProcessingFunctional(
        new CopyPopulationsFunctional3D<T, Descriptor>(), domain, latticeFrom, latticeTo);
}

template <typename T, template <typename U> class Descriptor>
void copyAll(
    BlockLattice3D<T, Descriptor> &latticeFrom, BlockLattice3D<T, Descriptor> &latticeTo,
    Box3D domain)
{
    applyProcessingFunctional(
        new LatticeCopyAllFunctional3D<T, Descriptor>(), domain, latticeFrom, latticeTo);
}

template <typename T, template <typename U> class Descriptor>
void copyRegenerate(
    BlockLattice3D<T, Descriptor> &latticeFrom, BlockLattice3D<T, Descriptor> &latticeTo,
    Box3D domain)
{
    applyProcessingFunctional(
        new LatticeRegenerateFunctional3D<T, Descriptor>(), domain, latticeFrom, latticeTo);
}

/* ******************************************************************* */
/* *************** PART II. Atomic-block wrappers: Scalar-Field ****** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template <typename T>
T computeSum(ScalarField3D<T> &scalarField, Box3D domain)
{
    BoxScalarSumFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getSumScalar();
}

template <typename T>
T computeSum(ScalarField3D<T> &scalarField)
{
    return computeSum(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeSum(ScalarField3D<T> &scalarField, ScalarField3D<int> &mask, int flag, Box3D domain)
{
    MaskedBoxScalarSumFunctional3D<T> functional(flag);
    applyProcessingFunctional(functional, domain, scalarField, mask);
    return functional.getSumScalar();
}

template <typename T>
T computeSum(ScalarField3D<T> &scalarField, ScalarField3D<int> &mask, int flag)
{
    return computeSum(scalarField, mask, flag, scalarField.getBoundingBox());
}

template <typename T>
T computeBoundedSum(ScalarField3D<T> &scalarField, Box3D domain)
{
    BoundedBoxScalarSumFunctional3D<T> functional;
    plint envelopeWidth = 1;
    applyProcessingFunctional(functional, domain, scalarField, envelopeWidth);
    return functional.getSumScalar();
}

template <typename T>
T computeBoundedSum(ScalarField3D<T> &scalarField)
{
    return computeBoundedSum(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeAverage(ScalarField3D<T> &scalarField, Box3D domain)
{
    BoxScalarSumFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getSumScalar() / (T)domain.nCells();
}

template <typename T>
T computeAverage(ScalarField3D<T> &scalarField)
{
    return computeAverage(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeAverage(ScalarField3D<T> &scalarField, ScalarField3D<int> &mask, int flag, Box3D domain)
{
    MaskedBoxScalarAverageFunctional3D<T> functional(flag);
    applyProcessingFunctional(functional, domain, scalarField, mask);
    return functional.getAverageScalar();
}

template <typename T>
T computeAverage(ScalarField3D<T> &scalarField, ScalarField3D<int> &mask, int flag)
{
    return computeAverage(scalarField, mask, flag, scalarField.getBoundingBox());
}

template <typename T>
T computeMin(ScalarField3D<T> &scalarField, Box3D domain)
{
    BoxScalarMinFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getMinScalar();
}

template <typename T>
T computeMin(ScalarField3D<T> &scalarField)
{
    return computeMin(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeMax(ScalarField3D<T> &scalarField, Box3D domain)
{
    BoxScalarMaxFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getMaxScalar();
}

template <typename T>
T computeMax(ScalarField3D<T> &scalarField)
{
    return computeMax(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeBoundedAverage(ScalarField3D<T> &scalarField, Box3D domain)
{
    BoundedBoxScalarSumFunctional3D<T> functional;
    plint envelopeWidth = 1;
    applyProcessingFunctional(functional, domain, scalarField, envelopeWidth);
    return functional.getSumScalar()
           / (T)((domain.getNx() - 1) * (domain.getNy() - 1) * (domain.getNz() - 1));
}

template <typename T>
T computeBoundedAverage(ScalarField3D<T> &scalarField)
{
    return computeBoundedAverage(scalarField, scalarField.getBoundingBox());
}

template <typename T, class BoolMask>
plint count(ScalarField3D<T> &field, Box3D domain, BoolMask boolMask)
{
    CountScalarElementsFunctional3D<T, BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, field);
    return functional.getCount();
}

template <typename T, class BoolMask>
plint count(ScalarField3D<T> &field, BoolMask boolMask)
{
    return count(field, field.getBoundingBox(), boolMask);
}

/* *************** ScalarField - Scalar operations *************** */

template <typename T>
void add(ScalarField3D<T> &field, T scalar, ScalarField3D<T> &result)
{
    applyProcessingFunctional(
        new A_plus_alpha_functional3D<T>(scalar), field.getBoundingBox(), field, result);
}

template <typename T>
std::unique_ptr<ScalarField3D<T> > add(ScalarField3D<T> &field, T scalar)
{
    ScalarField3D<T> *result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    add(field, scalar, *result);
    return std::unique_ptr<ScalarField3D<T> >(result);
}

template <typename T>
void add(T scalar, ScalarField3D<T> &field, ScalarField3D<T> &result)
{
    applyProcessingFunctional(
        new A_plus_alpha_functional3D<T>(scalar), field.getBoundingBox(), field, result);
}

template <typename T>
std::unique_ptr<ScalarField3D<T> > add(T scalar, ScalarField3D<T> &field)
{
    ScalarField3D<T> *result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    add(scalar, field, *result);
    return std::unique_ptr<ScalarField3D<T> >(result);
}

template <typename T>
void subtract(ScalarField3D<T> &field, T scalar, ScalarField3D<T> &result)
{
    applyProcessingFunctional(
        new A_minus_alpha_functional3D<T>(scalar), field.getBoundingBox(), field, result);
}

template <typename T>
std::unique_ptr<ScalarField3D<T> > subtract(ScalarField3D<T> &field, T scalar)
{
    ScalarField3D<T> *result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    subtract(field, scalar, *result);
    return std::unique_ptr<ScalarField3D<T> >(result);
}

template <typename T>
void subtract(T scalar, ScalarField3D<T> &field, ScalarField3D<T> &result)
{
    applyProcessingFunctional(
        new Alpha_minus_A_functional3D<T>(scalar), field.getBoundingBox(), field, result);
}

template <typename T>
std::unique_ptr<ScalarField3D<T> > subtract(T scalar, ScalarField3D<T> &field)
{
    ScalarField3D<T> *result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    subtract(scalar, field, *result);
    return std::unique_ptr<ScalarField3D<T> >(result);
}

template <typename T>
void multiply(ScalarField3D<T> &field, T scalar, ScalarField3D<T> &result)
{
    applyProcessingFunctional(
        new A_times_alpha_functional3D<T>(scalar), field.getBoundingBox(), field, result);
}

template <typename T>
std::unique_ptr<ScalarField3D<T> > multiply(ScalarField3D<T> &field, T scalar)
{
    ScalarField3D<T> *result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    multiply(field, scalar, *result);
    return std::unique_ptr<ScalarField3D<T> >(result);
}

template <typename T>
void multiply(T scalar, ScalarField3D<T> &field, ScalarField3D<T> &result)
{
    applyProcessingFunctional(
        new A_times_alpha_functional3D<T>(scalar), field.getBoundingBox(), field, result);
}

template <typename T>
std::unique_ptr<ScalarField3D<T> > multiply(T scalar, ScalarField3D<T> &field)
{
    ScalarField3D<T> *result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    multiply(scalar, field, *result);
    return std::unique_ptr<ScalarField3D<T> >(result);
}

template <typename T>
void divide(ScalarField3D<T> &field, T scalar, ScalarField3D<T> &result)
{
    applyProcessingFunctional(
        new A_dividedBy_alpha_functional3D<T>(scalar), field.getBoundingBox(), field, result);
}

template <typename T>
std::unique_ptr<ScalarField3D<T> > divide(ScalarField3D<T> &field, T scalar)
{
    ScalarField3D<T> *result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    divide(field, scalar, *result);
    return std::unique_ptr<ScalarField3D<T> >(result);
}

template <typename T>
void divide(T scalar, ScalarField3D<T> &field, ScalarField3D<T> &result)
{
    applyProcessingFunctional(
        new Alpha_dividedBy_A_functional3D<T>(scalar), field.getBoundingBox(), field, result);
}

template <typename T>
std::unique_ptr<ScalarField3D<T> > divide(T scalar, ScalarField3D<T> &field)
{
    ScalarField3D<T> *result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    divide(scalar, field, *result);
    return std::unique_ptr<ScalarField3D<T> >(result);
}

/* *************** ScalarField - Scalar inplace operations *************** */

template <typename T>
void addInPlace(ScalarField3D<T> &field, T scalar)
{
    applyProcessingFunctional(
        new A_plus_alpha_inplace_functional3D<T>(scalar), field.getBoundingBox(), field);
}

template <typename T>
void subtractInPlace(ScalarField3D<T> &field, T scalar)
{
    applyProcessingFunctional(
        new A_minus_alpha_inplace_functional3D<T>(scalar), field.getBoundingBox(), field);
}

template <typename T>
void multiplyInPlace(ScalarField3D<T> &field, T scalar)
{
    applyProcessingFunctional(
        new A_times_alpha_inplace_functional3D<T>(scalar), field.getBoundingBox(), field);
}

template <typename T>
void divideInPlace(ScalarField3D<T> &field, T scalar)
{
    applyProcessingFunctional(
        new A_dividedBy_alpha_inplace_functional3D<T>(scalar), field.getBoundingBox(), field);
}

/* *************** ScalarField - ScalarField operations *************** */

template <typename T1, typename T2>
void copy(ScalarField3D<T1> &field, ScalarField3D<T2> &convertedField)
{
    applyProcessingFunctional(
        new CopyConvertScalarFunctional3D<T1, T2>, field.getBoundingBox(), field, convertedField);
}

template <typename T1, typename T2>
std::unique_ptr<ScalarField3D<T2> > copyConvert(ScalarField3D<T1> &field)
{
    ScalarField3D<T2> *convertedField =
        new ScalarField3D<T2>(field.getNx(), field.getNy(), field.getNz());
    plb::copy(field, *convertedField);
    return std::unique_ptr<ScalarField3D<T2> >(convertedField);
}

template <typename T>
void add(ScalarField3D<T> &A, ScalarField3D<T> &B, ScalarField3D<T> &result)
{
    std::vector<ScalarField3D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_plus_B_functional3D<T>, A.getBoundingBox(), fields);
}

template <typename T>
std::unique_ptr<ScalarField3D<T> > add(ScalarField3D<T> &A, ScalarField3D<T> &B)
{
    ScalarField3D<T> *result = new ScalarField3D<T>(A.getNx(), A.getNy(), A.getNz());
    add(A, B, *result);
    return std::unique_ptr<ScalarField3D<T> >(result);
}

template <typename T>
void subtract(ScalarField3D<T> &A, ScalarField3D<T> &B, ScalarField3D<T> &result)
{
    std::vector<ScalarField3D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_minus_B_functional3D<T>, A.getBoundingBox(), fields);
}

template <typename T>
std::unique_ptr<ScalarField3D<T> > subtract(ScalarField3D<T> &A, ScalarField3D<T> &B)
{
    ScalarField3D<T> *result = new ScalarField3D<T>(A.getNx(), A.getNy(), A.getNz());
    subtract(A, B, *result);
    return std::unique_ptr<ScalarField3D<T> >(result);
}

template <typename T>
void multiply(ScalarField3D<T> &A, ScalarField3D<T> &B, ScalarField3D<T> &result)
{
    std::vector<ScalarField3D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_times_B_functional3D<T>, A.getBoundingBox(), fields);
}

template <typename T>
std::unique_ptr<ScalarField3D<T> > multiply(ScalarField3D<T> &A, ScalarField3D<T> &B)
{
    ScalarField3D<T> *result = new ScalarField3D<T>(A.getNx(), A.getNy(), A.getNz());
    multiply(A, B, *result);
    return std::unique_ptr<ScalarField3D<T> >(result);
}

template <typename T>
void divide(ScalarField3D<T> &A, ScalarField3D<T> &B, ScalarField3D<T> &result)
{
    std::vector<ScalarField3D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_dividedBy_B_functional3D<T>, A.getBoundingBox(), fields);
}

template <typename T>
std::unique_ptr<ScalarField3D<T> > divide(ScalarField3D<T> &A, ScalarField3D<T> &B)
{
    ScalarField3D<T> *result = new ScalarField3D<T>(A.getNx(), A.getNy(), A.getNz());
    divide(A, B, *result);
    return std::unique_ptr<ScalarField3D<T> >(result);
}

/* *************** ScalarField operations *************** */

template <typename T>
void computeSqrt(ScalarField3D<T> &A, ScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(new ComputeScalarSqrtFunctional3D<T>, domain, A, result);
}

template <typename T>
std::unique_ptr<ScalarField3D<T> > computeSqrt(ScalarField3D<T> &A, Box3D domain)
{
    ScalarField3D<T> *result = new ScalarField3D<T>(A.getNx(), A.getNy(), A.getNz());
    computeSqrt(A, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<ScalarField3D<T> > computeSqrt(ScalarField3D<T> &A)
{
    return computeSqrt(A, A.getBoundingBox());
}

template <typename T>
void computeAbsoluteValue(ScalarField3D<T> &A, ScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(new ComputeAbsoluteValueFunctional3D<T>, domain, A, result);
}

template <typename T>
std::unique_ptr<ScalarField3D<T> > computeAbsoluteValue(ScalarField3D<T> &A, Box3D domain)
{
    ScalarField3D<T> *result = new ScalarField3D<T>(A.getNx(), A.getNy(), A.getNz());
    computeAbsoluteValue(A, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<ScalarField3D<T> > computeAbsoluteValue(ScalarField3D<T> &A)
{
    return computeAbsoluteValue(A, A.getBoundingBox());
}

/* *************** ScalarField - ScalarField inplace operations *************** */

template <typename T>
void addInPlace(ScalarField3D<T> &A, ScalarField3D<T> &B)
{
    applyProcessingFunctional(new A_plus_B_inplace_functional3D<T>, A.getBoundingBox(), A, B);
}

template <typename T>
void subtractInPlace(ScalarField3D<T> &A, ScalarField3D<T> &B)
{
    applyProcessingFunctional(new A_minus_B_inplace_functional3D<T>, A.getBoundingBox(), A, B);
}

template <typename T>
void multiplyInPlace(ScalarField3D<T> &A, ScalarField3D<T> &B)
{
    applyProcessingFunctional(new A_times_B_inplace_functional3D<T>, A.getBoundingBox(), A, B);
}

template <typename T>
void divideInPlace(ScalarField3D<T> &A, ScalarField3D<T> &B)
{
    applyProcessingFunctional(new A_dividedBy_B_inplace_functional3D<T>, A.getBoundingBox(), A, B);
}

/* ******************************************************************* */
/* *************** PART III. Atomic-block wrappers: Tensor-Field ***** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template <typename T, plint nDim, class BoolMask>
plint count(TensorField3D<T, nDim> &field, Box3D domain, BoolMask boolMask)
{
    CountTensorElementsFunctional3D<T, nDim, BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, field);
    return functional.getCount();
}

template <typename T, plint nDim, class BoolMask>
plint count(TensorField3D<T, nDim> &field, BoolMask boolMask)
{
    return count(field, field.getBoundingBox(), boolMask);
}

template <typename T, int nDim>
Array<T, nDim> computeSum(TensorField3D<T, nDim> &tensorField, Box3D domain)
{
    BoxTensorSumFunctional3D<T, nDim> functional;
    applyProcessingFunctional(functional, domain, tensorField);
    return functional.getSumTensor();
}

template <typename T, int nDim>
Array<T, nDim> computeSum(TensorField3D<T, nDim> &tensorField)
{
    return computeSum(tensorField, tensorField.getBoundingBox());
}

template <typename T, int nDim>
Array<T, nDim> computeSum(
    TensorField3D<T, nDim> &tensorField, ScalarField3D<int> &mask, int flag, Box3D domain)
{
    MaskedBoxTensorSumFunctional3D<T, nDim> functional(flag);
    applyProcessingFunctional(functional, domain, mask, tensorField);
    return functional.getSumTensor();
}

template <typename T, int nDim>
Array<T, nDim> computeSum(TensorField3D<T, nDim> &tensorField, ScalarField3D<int> &mask, int flag)
{
    return computeSum(tensorField, mask, flag, tensorField.getBoundingBox());
}

template <typename T, int nDim>
Array<T, nDim> computeAverage(TensorField3D<T, nDim> &tensorField, Box3D domain)
{
    BoxTensorSumFunctional3D<T, nDim> functional;
    applyProcessingFunctional(functional, domain, tensorField);
    return functional.getSumTensor() / (T)domain.nCells();
}

template <typename T, int nDim>
Array<T, nDim> computeAverage(TensorField3D<T, nDim> &tensorField)
{
    return computeAverage(tensorField, tensorField.getBoundingBox());
}

template <typename T, int nDim>
Array<T, nDim> computeAverage(
    TensorField3D<T, nDim> &tensorField, ScalarField3D<int> &mask, int flag, Box3D domain)
{
    MaskedBoxTensorAverageFunctional3D<T, nDim> functional(flag);
    applyProcessingFunctional(functional, domain, mask, tensorField);
    return functional.getAverageTensor();
}

template <typename T, int nDim>
Array<T, nDim> computeAverage(
    TensorField3D<T, nDim> &tensorField, ScalarField3D<int> &mask, int flag)
{
    return computeAverage(tensorField, mask, flag, tensorField.getBoundingBox());
}

/* *************** Component (scalar-field) out of a tensor-field ****** */

template <typename T, int nDim>
void extractComponent(
    TensorField3D<T, nDim> &tensorField, ScalarField3D<T> &component, int iComponent)
{
    applyProcessingFunctional(
        new ExtractTensorComponentFunctional3D<T, nDim>(iComponent), tensorField.getBoundingBox(),
        component, tensorField);
}

template <typename T, int nDim>
std::unique_ptr<ScalarField3D<T> > extractComponent(
    TensorField3D<T, nDim> &tensorField, int iComponent)
{
    ScalarField3D<T> *component =
        new ScalarField3D<T>(tensorField.getNx(), tensorField.getNy(), tensorField.getNz());
    extractComponent(tensorField, *component, iComponent);
    return std::unique_ptr<ScalarField3D<T> >(component);
}

/* *************** Vector-norm of each cell in the field *************** */

template <typename T, int nDim>
void computeNorm(TensorField3D<T, nDim> &tensorField, ScalarField3D<T> &component)
{
    applyProcessingFunctional(
        new ComputeNormFunctional3D<T, nDim>, tensorField.getBoundingBox(), component, tensorField);
}

template <typename T, int nDim>
std::unique_ptr<ScalarField3D<T> > computeNorm(TensorField3D<T, nDim> &tensorField)
{
    ScalarField3D<T> *component =
        new ScalarField3D<T>(tensorField.getNx(), tensorField.getNy(), tensorField.getNz());
    computeNorm(tensorField, *component);
    return std::unique_ptr<ScalarField3D<T> >(component);
}

/* *************** Sqrt operation on each component of each cell *************** */

template <typename T, int nDim>
void computeSqrt(TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &result, Box3D domain)
{
    applyProcessingFunctional(new ComputeTensorSqrtFunctional3D<T, nDim>, domain, A, result);
}

template <typename T, int nDim>
std::unique_ptr<TensorField3D<T, nDim> > computeSqrt(TensorField3D<T, nDim> &A, Box3D domain)
{
    TensorField3D<T, nDim> *result = new TensorField3D<T, nDim>(A.getNx(), A.getNy(), A.getNz());
    computeSqrt(A, *result, domain);
    return result;
}

template <typename T, int nDim>
std::unique_ptr<TensorField3D<T, nDim> > computeSqrt(TensorField3D<T, nDim> &A)
{
    return computeSqrt(A, A.getBoundingBox());
}

/* *************** Squared vector-norm of each cell in the field ******** */

template <typename T, int nDim>
void computeNormSqr(TensorField3D<T, nDim> &tensorField, ScalarField3D<T> &component)
{
    applyProcessingFunctional(
        new ComputeNormSqrFunctional3D<T, nDim>, tensorField.getBoundingBox(), component,
        tensorField);
}

template <typename T, int nDim>
std::unique_ptr<ScalarField3D<T> > computeNormSqr(TensorField3D<T, nDim> &tensorField)
{
    ScalarField3D<T> *component =
        new ScalarField3D<T>(tensorField.getNx(), tensorField.getNy(), tensorField.getNz());
    computeNormSqr(tensorField, *component);
    return std::unique_ptr<ScalarField3D<T> >(component);
}

/* *************** Max element of each array of each cell *************** */

template <typename T, int nDim>
void computeMaximumElement(
    MultiTensorField3D<T, nDim> &A, MultiScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(
        new BoxLocalMaximumPerComponentFunctional3D<T, nDim>(), domain, result, A);
}

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T> > computeMaximumElement(
    MultiTensorField3D<T, nDim> &A, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > result = generateMultiScalarField<T>(A, domain);

    computeMaximumElement(A, *result, domain);
    return result;
}

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T> > computeMaximumElement(MultiTensorField3D<T, nDim> &A)
{
    return computeMaximumElement(A, A.getBoundingBox());
}

/* *************** Tensor-norm of each symmetric tensor of a field ********************** */

template <typename T>
void computeSymmetricTensorNorm(TensorField3D<T, 6> &tensorField, ScalarField3D<T> &norm)
{
    applyProcessingFunctional(
        new ComputeSymmetricTensorNormFunctional3D<T>, tensorField.getBoundingBox(), norm,
        tensorField);
}

template <typename T>
std::unique_ptr<ScalarField3D<T> > computeSymmetricTensorNorm(TensorField3D<T, 6> &tensorField)
{
    ScalarField3D<T> *norm =
        new ScalarField3D<T>(tensorField.getNx(), tensorField.getNy(), tensorField.getNz());
    computeSymmetricTensorNorm(tensorField, *norm);
    return std::unique_ptr<ScalarField3D<T> >(norm);
}

/* *************** Squared Tensor-norm of each symmetric tensor of a field *************** */

template <typename T>
void computeSymmetricTensorNormSqr(TensorField3D<T, 6> &tensorField, ScalarField3D<T> &normSqr)
{
    applyProcessingFunctional(
        new ComputeSymmetricTensorNormSqrFunctional3D<T>, tensorField.getBoundingBox(), normSqr,
        tensorField);
}

template <typename T>
std::unique_ptr<ScalarField3D<T> > computeSymmetricTensorNormSqr(TensorField3D<T, 6> &tensorField)
{
    ScalarField3D<T> *normSqr =
        new ScalarField3D<T>(tensorField.getNx(), tensorField.getNy(), tensorField.getNz());
    computeSymmetricTensorNormSqr(tensorField, *normSqr);
    return std::unique_ptr<ScalarField3D<T> >(normSqr);
}

/* *************** Trace of each symmetric tensor of a field ************ */

template <typename T>
void computeSymmetricTensorTrace(TensorField3D<T, 6> &tensorField, ScalarField3D<T> &trace)
{
    applyProcessingFunctional(
        new ComputeSymmetricTensorTraceFunctional3D<T>, tensorField.getBoundingBox(), trace,
        tensorField);
}

template <typename T>
std::unique_ptr<ScalarField3D<T> > computeSymmetricTensorTrace(TensorField3D<T, 6> &tensorField)
{
    ScalarField3D<T> *trace =
        new ScalarField3D<T>(tensorField.getNx(), tensorField.getNy(), tensorField.getNz());
    computeSymmetricTensorTrace(tensorField, *trace);
    return std::unique_ptr<ScalarField3D<T> >(trace);
}

/* *************** Vorticity from Velocity field *********************** */

template <typename T>
void computeVorticity(TensorField3D<T, 3> &velocity, TensorField3D<T, 3> &vorticity)
{
    plint envelopeWidth = 1;
    applyProcessingFunctional(
        new BoxVorticityFunctional3D<T, 3>, velocity.getBoundingBox(), velocity, vorticity,
        envelopeWidth);
}

template <typename T>
std::unique_ptr<TensorField3D<T, 3> > computeVorticity(TensorField3D<T, 3> &velocity)
{
    TensorField3D<T, 3> *vorticity =
        new TensorField3D<T, 3>(velocity.getNx(), velocity.getNy(), velocity.getNz());
    computeVorticity(velocity, *vorticity);
    return std::unique_ptr<TensorField3D<T, 3> >(vorticity);
}

/* *************** Vorticity, witout boundary treatment, from Velocity field **************** */

template <typename T>
void computeBulkVorticity(TensorField3D<T, 3> &velocity, TensorField3D<T, 3> &vorticity)
{
    applyProcessingFunctional(
        new BoxBulkVorticityFunctional3D<T, 3>, velocity.getBoundingBox(), velocity, vorticity);
}

template <typename T>
std::unique_ptr<TensorField3D<T, 3> > computeBulkVorticity(TensorField3D<T, 3> &velocity)
{
    TensorField3D<T, 3> *vorticity =
        new TensorField3D<T, 3>(velocity.getNx(), velocity.getNy(), velocity.getNz());
    computeBulkVorticity(velocity, *vorticity);
    return std::unique_ptr<TensorField3D<T, 3> >(vorticity);
}

/* *************** Helicity from Velocity field *********************** */

template <typename T>
void computeHelicity(TensorField3D<T, 3> &velocity, ScalarField3D<T> &helicity)
{
    plint envelopeWidth = 1;
    applyProcessingFunctional(
        new BoxHelicityFunctional3D<T, 3>, velocity.getBoundingBox(), helicity, velocity,
        envelopeWidth);
}

template <typename T>
std::unique_ptr<ScalarField3D<T> > computeHelicity(TensorField3D<T, 3> &velocity)
{
    ScalarField3D<T> *helicity =
        new ScalarField3D<T>(velocity.getNx(), velocity.getNy(), velocity.getNz());
    computeHelicity(velocity, *helicity);
    return std::unique_ptr<ScalarField3D<T> >(helicity);
}

/* *************** Helicity, witout boundary treatment, from Velocity field **************** */

template <typename T>
void computeBulkHelicity(TensorField3D<T, 3> &velocity, ScalarField3D<T> &helicity)
{
    applyProcessingFunctional(
        new BoxBulkHelicityFunctional3D<T, 3>, velocity.getBoundingBox(), helicity, velocity);
}

template <typename T>
std::unique_ptr<ScalarField3D<T> > computeBulkHelicity(TensorField3D<T, 3> &velocity)
{
    ScalarField3D<T> *helicity =
        new ScalarField3D<T>(velocity.getNx(), velocity.getNy(), velocity.getNz());
    computeBulkHelicity(velocity, *helicity);
    return std::unique_ptr<ScalarField3D<T> >(helicity);
}

/* *************** Divergence, witout boundary treatment, from Velocity field *************** */

template <typename T>
void computeBulkDivergence(TensorField3D<T, 3> &velocity, ScalarField3D<T> &divergence)
{
    applyProcessingFunctional(
        new BoxBulkDivergenceFunctional3D<T, 3>, velocity.getBoundingBox(), divergence, velocity);
}

template <typename T>
std::unique_ptr<ScalarField3D<T> > computeBulkDivergence(TensorField3D<T, 3> &velocity)
{
    ScalarField3D<T> *divergence =
        new ScalarField3D<T>(velocity.getNx(), velocity.getNy(), velocity.getNz());
    computeBulkDivergence(velocity, *divergence);
    return std::unique_ptr<ScalarField3D<T> >(divergence);
}

/* *************** Strain rate from Velocity field ********************* */

template <typename T>
void computeStrainRate(TensorField3D<T, 3> &velocity, TensorField3D<T, 6> &S)
{
    plint envelopeWidth = 1;
    applyProcessingFunctional(
        new BoxStrainRateFunctional3D<T, 3>, velocity.getBoundingBox(), velocity, S, envelopeWidth);
}

template <typename T>
std::unique_ptr<TensorField3D<T, 3> > computeStrainRate(TensorField3D<T, 3> &velocity)
{
    TensorField3D<T, 6> *S =
        new TensorField3D<T, 6>(velocity.getNx(), velocity.getNy(), velocity.getNz());
    computeStrainRate(velocity, *S);
    return std::unique_ptr<TensorField3D<T, 6> >(S);
}

/* *************** Str. rate, witout boundary treatment, from Velocity field ***************** */

template <typename T>
void computeBulkStrainRate(TensorField3D<T, 3> &velocity, TensorField3D<T, 6> &S)
{
    applyProcessingFunctional(
        new BoxBulkStrainRateFunctional3D<T, 3>, velocity.getBoundingBox(), velocity, S);
}

template <typename T>
std::unique_ptr<TensorField3D<T, 6> > computeBulkStrainRate(TensorField3D<T, 3> &velocity)
{
    TensorField3D<T, 6> *S =
        new TensorField3D<T, 6>(velocity.getNx(), velocity.getNy(), velocity.getNz());
    computeBulkStrainRate(velocity, *S);
    return std::unique_ptr<TensorField3D<T, 6> >(S);
}

/* *************** TensorField - TensorField operations *************** */

template <typename T1, typename T2, int nDim>
void copy(TensorField3D<T1, nDim> &field, TensorField3D<T2, nDim> &convertedField)
{
    applyProcessingFunctional(
        new CopyConvertTensorFunctional3D<T1, T2, nDim>, field.getBoundingBox(), field,
        convertedField);
}

template <typename T1, typename T2, int nDim>
std::unique_ptr<TensorField3D<T2, nDim> > copyConvert(TensorField3D<T1, nDim> &field)
{
    TensorField3D<T2, nDim> *convertedField =
        new TensorField3D<T2, nDim>(field.getNx(), field.getNy(), field.getNz());
    plb::copy(field, *convertedField);
    return std::unique_ptr<TensorField3D<T2, nDim> >(convertedField);
}

template <typename T, int nDim>
void add(TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B, TensorField3D<T, nDim> &result)
{
    std::vector<TensorField3D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(
        new Tensor_A_plus_B_functional3D<T, nDim>, A.getBoundingBox(), fields);
}

template <typename T, int nDim>
std::unique_ptr<TensorField3D<T, nDim> > add(TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B)
{
    TensorField3D<T, nDim> *result = new TensorField3D<T, nDim>(A.getNx(), A.getNy(), A.getNz());
    add(A, B, *result);
    return std::unique_ptr<TensorField3D<T, nDim> >(result);
}

template <typename T, int nDim>
void subtract(TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B, TensorField3D<T, nDim> &result)
{
    std::vector<TensorField3D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(
        new Tensor_A_minus_B_functional3D<T, nDim>, A.getBoundingBox(), fields);
}

template <typename T, int nDim>
std::unique_ptr<TensorField3D<T, nDim> > subtract(
    TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B)
{
    TensorField3D<T, nDim> *result = new TensorField3D<T, nDim>(A.getNx(), A.getNy(), A.getNz());
    subtract(A, B, *result);
    return std::unique_ptr<TensorField3D<T, nDim> >(result);
}

template <typename T, int nDim>
void multiply(TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B, TensorField3D<T, nDim> &result)
{
    std::vector<TensorField3D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(
        new Tensor_A_times_B_functional3D<T, nDim>, A.getBoundingBox(), fields);
}

template <typename T, int nDim>
std::unique_ptr<TensorField3D<T, nDim> > multiply(
    TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B)
{
    TensorField3D<T, nDim> *result = new TensorField3D<T, nDim>(A.getNx(), A.getNy(), A.getNz());
    multiply(A, B, *result);
    return std::unique_ptr<TensorField3D<T, nDim> >(result);
}

template <typename T, int nDim>
void multiply(TensorField3D<T, nDim> &field, T scalar, TensorField3D<T, nDim> &result, Box3D domain)
{
    applyProcessingFunctional(
        new Tensor_A_times_alpha_functional3D<T, nDim>(scalar), domain, field, result);
}

template <typename T, int nDim>
std::unique_ptr<TensorField3D<T, nDim> > multiply(TensorField3D<T, nDim> &field, T scalar)
{
    TensorField3D<T, nDim> *result =
        new TensorField3D<T, nDim>(field.getNx(), field.getNy(), field.getNz());
    multiply(field, scalar, *result);
    return std::unique_ptr<TensorField3D<T, nDim> >(result);
}

template <typename T, int nDim>
void multiply(T scalar, TensorField3D<T, nDim> &field, TensorField3D<T, nDim> &result, Box3D domain)
{
    applyProcessingFunctional(
        new Tensor_A_times_alpha_functional3D<T, nDim>(scalar), domain, field, result);
}

template <typename T, int nDim>
std::unique_ptr<TensorField3D<T, nDim> > multiply(T scalar, TensorField3D<T, nDim> &field)
{
    TensorField3D<T, nDim> *result =
        new TensorField3D<T, nDim>(field.getNx(), field.getNy(), field.getNz());
    multiply(scalar, field, *result);
    return std::unique_ptr<TensorField3D<T, nDim> >(result);
}

template <typename T, int nDim>
void divide(TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B, TensorField3D<T, nDim> &result)
{
    std::vector<TensorField3D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(
        new Tensor_A_dividedBy_B_functional3D<T, nDim>, A.getBoundingBox(), fields);
}

template <typename T, int nDim>
std::unique_ptr<TensorField3D<T, nDim> > divide(
    TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B)
{
    TensorField3D<T, nDim> *result = new TensorField3D<T, nDim>(A.getNx(), A.getNy(), A.getNz());
    divide(A, B, *result);
    return std::unique_ptr<TensorField3D<T, nDim> >(result);
}

/* *************** TensorField - TensorField inplace operations *************** */

template <typename T, int nDim>
void addInPlace(TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B)
{
    applyProcessingFunctional(
        new Tensor_A_plus_B_inplace_functional3D<T, nDim>, A.getBoundingBox(), A, B);
}

template <typename T, int nDim>
void subtractInPlace(TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B)
{
    applyProcessingFunctional(
        new Tensor_A_minus_B_inplace_functional3D<T, nDim>, A.getBoundingBox(), A, B);
}

template <typename T, int nDim>
void multiplyInPlace(TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B)
{
    applyProcessingFunctional(
        new Tensor_A_times_B_inplace_functional3D<T, nDim>, A.getBoundingBox(), A, B);
}

template <typename T, int nDim>
void multiplyInPlace(TensorField3D<T, nDim> &A, T alpha)
{
    applyProcessingFunctional(
        new Tensor_A_times_alpha_inplace_functional3D<T, nDim>(alpha), A.getBoundingBox(), A);
}

template <typename T, int nDim>
void divideInPlace(TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B)
{
    applyProcessingFunctional(
        new Tensor_A_dividedBy_B_inplace_functional3D<T, nDim>, A.getBoundingBox(), A, B);
}

/* ******************************************************************* */
/* *************** PART IV. Multi-block wrappers: Block-Lattice ****** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template <typename T, template <typename U> class Descriptor>
T computeAverageDensity(MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    BoxSumRhoBarFunctional3D<T, Descriptor> functional;
    applyProcessingFunctional(functional, domain, lattice);
    return Descriptor<T>::fullRho(functional.getSumRhoBar() / (T)domain.nCells());
}

template <typename T, template <typename U> class Descriptor>
T computeAverageDensity(MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeAverageDensity(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    BoxSumRhoBarFunctional3D<T, Descriptor> functional;
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getSumRhoBar() / (T)domain.nCells();
}

template <typename T, template <typename U> class Descriptor>
T computeAverageRhoBar(MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeAverageRhoBar(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
T computeAverageEnergy(MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    BoxSumEnergyFunctional3D<T, Descriptor> functional;
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getSumEnergy() / (T)domain.nCells();
}

template <typename T, template <typename U> class Descriptor>
T computeAverageEnergy(MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeAverageEnergy(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
T computeAverageForcedEnergy(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    Box3D domain)
{
    BoxSumForcedEnergyFunctional3D<T, Descriptor> functional;
    applyProcessingFunctional(functional, domain, lattice, force);
    return functional.getSumEnergy() / (T)domain.nCells();
}

template <typename T, template <typename U> class Descriptor>
T computeAverageForcedEnergy(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force)
{
    return computeAverageForcedEnergy(lattice, force, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
T computeAverageForcedEnergy(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force, Box3D domain)
{
    BoxSumConstForcedEnergyFunctional3D<T, Descriptor> functional(force);
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getSumEnergy() / (T)domain.nCells();
}

template <typename T, template <typename U> class Descriptor>
T computeAverageForcedEnergy(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force)
{
    return computeAverageForcedEnergy(lattice, force, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
T computeAverageForcedEnergy(
    MultiBlockLattice3D<T, Descriptor> &lattice, ForceFunction f, Box3D domain)
{
    BoxSumCustomForcedEnergyFunctional3D<T, Descriptor, ForceFunction> functional(f);
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getSumEnergy() / (T)domain.nCells();
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
T computeAverageForcedEnergy(MultiBlockLattice3D<T, Descriptor> &lattice, ForceFunction f)
{
    return computeAverageForcedEnergy(lattice, f, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor, class BoolMask>
plint count(MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, BoolMask boolMask)
{
    CountLatticeElementsFunctional3D<T, Descriptor, BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getCount();
}

template <typename T, template <typename U> class Descriptor, class BoolMask>
plint count(MultiBlockLattice3D<T, Descriptor> &lattice, BoolMask boolMask)
{
    return count(lattice, lattice.getBoundingBox(), boolMask);
}

template <typename T>
std::vector<T> scalarSingleProbes(
    MultiScalarField3D<T> &scalarField, Box3D domain, std::vector<Array<T, 3> > const &positions)
{
    ScalarFieldSingleProbe3D<T> functional(positions);
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getScalars();
}

template <typename T>
std::vector<T> scalarSingleProbes(
    MultiScalarField3D<T> &scalarField, std::vector<Array<T, 3> > const &positions)
{
    return scalarSingleProbes(scalarField, scalarField.getBoundingBox(), positions);
}

template <typename T, template <typename U> class Descriptor>
std::vector<T> densitySingleProbes(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain,
    std::vector<Array<T, 3> > const &positions)
{
    DensitySingleProbe3D<T, Descriptor> functional(positions);
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getDensities();
}

template <typename T, template <typename U> class Descriptor>
std::vector<T> densitySingleProbes(
    MultiBlockLattice3D<T, Descriptor> &lattice, std::vector<Array<T, 3> > const &positions)
{
    return densitySingleProbes(lattice, lattice.getBoundingBox(), positions);
}

template <typename T, template <typename U> class Descriptor>
std::vector<Array<T, 3> > velocitySingleProbes(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain,
    std::vector<Array<T, 3> > const &positions)
{
    VelocitySingleProbe3D<T, Descriptor> functional(positions);
    ;
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getVelocities();
}

template <typename T, template <typename U> class Descriptor>
std::vector<Array<T, 3> > velocitySingleProbes(
    MultiBlockLattice3D<T, Descriptor> &lattice, std::vector<Array<T, 3> > const &positions)
{
    return velocitySingleProbes(lattice, lattice.getBoundingBox(), positions);
}

template <typename T, template <typename U> class Descriptor>
std::vector<Array<T, 3> > vorticitySingleProbes(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain,
    std::vector<Array<T, 3> > const &positions)
{
    VorticitySingleProbe3D<T, Descriptor> functional(positions);
    ;
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getVorticities();
}

template <typename T, template <typename U> class Descriptor>
std::vector<Array<T, 3> > vorticitySingleProbes(
    MultiBlockLattice3D<T, Descriptor> &lattice, std::vector<Array<T, 3> > const &positions)
{
    return vorticitySingleProbes(lattice, lattice.getBoundingBox(), positions);
}

/* *************** Extract Sub-Lattice ******************************* */

template <typename T, template <typename U> class Descriptor>
void extractSubDomain(
    MultiBlockLattice3D<T, Descriptor> &lattice,
    MultiBlockLattice3D<T, Descriptor> &extractedLattice, Box3D domain)
{
    applyProcessingFunctional(
        new LatticeRegenerateFunctional3D<T, Descriptor>, domain, lattice, extractedLattice);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > extractSubDomain(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > extractedLattice =
        generateMultiBlockLattice<T, Descriptor>(lattice, domain);
    extractSubDomain(lattice, *extractedLattice, domain);
    return extractedLattice;
}

/* *************** Density ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeDensity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &density, Box3D domain)
{
    applyProcessingFunctional(new BoxDensityFunctional3D<T, Descriptor>, domain, lattice, density);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeDensity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > density = generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the BoxDensityFunctional3D()
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
std::unique_ptr<MultiScalarField3D<T> > computeDensity(MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeDensity(lattice, lattice.getBoundingBox());
}

/* *************** RhoBar ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeRhoBar(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &rhoBar, Box3D domain)
{
    applyProcessingFunctional(new BoxRhoBarFunctional3D<T, Descriptor>, domain, lattice, rhoBar);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeRhoBar(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > rhoBar = generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the BoxRhoBarFunctional3D() acts
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
std::unique_ptr<MultiScalarField3D<T> > computeRhoBar(MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeRhoBar(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void computeRhoBarJ(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &rhoBar,
    MultiTensorField3D<T, 3> &j, Box3D domain)
{
    std::vector<MultiBlock3D *> fields;
    fields.push_back(&lattice);
    fields.push_back(&rhoBar);
    fields.push_back(&j);
    applyProcessingFunctional(new BoxRhoBarJfunctional3D<T, Descriptor>, domain, fields);
}

template <typename T, template <typename U> class Descriptor>
void maskedComputeRhoBarJ(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &rhoBar,
    MultiTensorField3D<T, 3> &j, MultiScalarField3D<int> &mask, int whichFlag, Box3D domain)
{
    std::vector<MultiBlock3D *> fields;
    fields.push_back(&lattice);
    fields.push_back(&rhoBar);
    fields.push_back(&j);
    fields.push_back(&mask);
    applyProcessingFunctional(
        new MaskedBoxRhoBarJfunctional3D<T, Descriptor>(whichFlag), domain, fields);
}

template <typename T, template <typename U> class Descriptor>
void computeRhoBarJPiNeq(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &rhoBar,
    MultiTensorField3D<T, 3> &j, TensorField3D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq,
    Box3D domain)
{
    std::vector<MultiBlock3D *> fields;
    fields.push_back(&lattice);
    fields.push_back(&rhoBar);
    fields.push_back(&j);
    fields.push_back(&PiNeq);
    applyProcessingFunctional(new BoxRhoBarJPiNeqfunctional3D<T, Descriptor>, domain, fields);
}

/* *************** Kinetic Energy ************************************ */

template <typename T, template <typename U> class Descriptor>
void computeKineticEnergy(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &energy, Box3D domain)
{
    applyProcessingFunctional(
        new BoxKineticEnergyFunctional3D<T, Descriptor>, domain, lattice, energy);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeKineticEnergy(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > energy = generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxKineticEnergyFunctional3D() acts on both bulk and envelope, you would expect the envelope
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
std::unique_ptr<MultiScalarField3D<T> > computeKineticEnergy(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeKineticEnergy(lattice, lattice.getBoundingBox());
}

/* *************** Packed RhoBar J *********************************** */

template <typename T, template <typename U> class Descriptor>
void computePackedRhoBarJ(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiNTensorField3D<T> &rhoBarJ, Box3D domain)
{
    applyProcessingFunctional(
        new PackedRhoBarJfunctional3D<T, Descriptor>, domain, lattice, rhoBarJ);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiNTensorField3D<T> > computePackedRhoBarJ(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    plint nDim = 1 + Descriptor<T>::d;
    std::unique_ptr<MultiNTensorField3D<T> > rhoBarJ(
        generateMultiNTensorField<T>(lattice, domain, nDim));

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the PackedRhoBarJFunctional3D()
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
std::unique_ptr<MultiScalarField3D<T> > computePackedRhoBarJ(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computePackedRhoBarJ(lattice, lattice.getBoundingBox());
}

template <typename T>
void computeDensityFromRhoBarJ(
    MultiNTensorField3D<T> &rhoBarJ, MultiScalarField3D<T> &density, Box3D domain)
{
    applyProcessingFunctional(new DensityFromRhoBarJfunctional3D<T>, domain, density, rhoBarJ);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeDensityFromRhoBarJ(
    MultiNTensorField3D<T> &rhoBarJ, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > density = generateMultiScalarField<T>(rhoBarJ, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // DensityFromRhoBarJFunctional3D() acts on both bulk and envelope, you would expect the
    // envelope layer around the domain, on the density multi-block, to be assigned some proper
    // values too. By default, this is however not what happens, because the physical space occupied
    // by these envelopes does not intersect with the domain "domain". We work around this issue by
    // extending the domain. There's no problem if the enlarged domain gets beyond the actual extent
    // of the lattice, because Palabos handles these situations properly.

    computeDensityFromRhoBarJ(
        rhoBarJ, *density, domain.enlarge(rhoBarJ.getMultiBlockManagement().getEnvelopeWidth()));
    return density;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeDensityFromRhoBarJ(MultiNTensorField3D<T> &rhoBarJ)
{
    return computeDensityFromRhoBarJ(rhoBarJ, rhoBarJ.getBoundingBox());
}

template <typename T>
void computeVelocityFromRhoBarJ(
    MultiNTensorField3D<T> &rhoBarJ, MultiTensorField3D<T, 3> &velocity, Box3D domain, bool velIsJ)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&velocity);
    args.push_back(&rhoBarJ);
    applyProcessingFunctional(new VelocityFromRhoBarJfunctional3D<T>(velIsJ), domain, args);
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeVelocityFromRhoBarJ(
    MultiNTensorField3D<T> &rhoBarJ, Box3D domain, bool velIsJ)
{
    std::unique_ptr<MultiTensorField3D<T, 3> > velocity =
        generateMultiTensorField<T, 3>(rhoBarJ, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // VelocityFromRhoBarJFunctional3D() acts on both bulk and envelope, you would expect the
    // envelope layer around the domain, on the velocity multi-block, to be assigned some proper
    // values too. By default, this is however not what happens, because the physical space occupied
    // by these envelopes does not intersect with the domain "domain". We work around this issue by
    // extending the domain. There's no problem if the enlarged domain gets beyond the actual extent
    // of the lattice, because Palabos handles these situations properly.

    computeVelocityFromRhoBarJ(
        rhoBarJ, *velocity, domain.enlarge(rhoBarJ.getMultiBlockManagement().getEnvelopeWidth()),
        velIsJ);
    return velocity;
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeVelocityFromRhoBarJ(
    MultiNTensorField3D<T> &rhoBarJ, bool velIsJ)
{
    return computeVelocityFromRhoBarJ(rhoBarJ, rhoBarJ.getBoundingBox(), velIsJ);
}

template <typename T>
void computeVelocityFromRhoBarAndJ(
    MultiScalarField3D<T> &rhoBar, MultiTensorField3D<T, 3> &j, MultiTensorField3D<T, 3> &velocity,
    Box3D domain, bool velIsJ)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&velocity);
    args.push_back(&rhoBar);
    args.push_back(&j);
    applyProcessingFunctional(new VelocityFromRhoBarAndJfunctional3D<T>(velIsJ), domain, args);
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeVelocityFromRhoBarAndJ(
    MultiScalarField3D<T> &rhoBar, MultiTensorField3D<T, 3> &j, Box3D domain, bool velIsJ)
{
    std::unique_ptr<MultiTensorField3D<T, 3> > velocity = generateMultiTensorField<T, 3>(j, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. If the VelocityFromRhoBarAndJFunctional3D()
    // acts on both bulk and envelope, you would expect the envelope layer around the domain, on the
    // velocity multi-block, to be assigned some proper values too. By default, this is however not
    // what happens, because the physical space occupied by these envelopes does not intersect with
    // the domain "domain". We work around this issue by extending the domain. There's no problem if
    // the enlarged domain gets beyond the actual extent of the lattice, because Palabos handles
    // these situations properly.

    computeVelocityFromRhoBarAndJ(
        rhoBar, j, *velocity, domain.enlarge(j.getMultiBlockManagement().getEnvelopeWidth()),
        velIsJ);
    return velocity;
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeVelocityFromRhoBarAndJ(
    MultiScalarField3D<T> &rhoBar, MultiTensorField3D<T, 3> &j, bool velIsJ)
{
    return computeVelocityFromRhoBarAndJ(rhoBar, j, rhoBar.getBoundingBox(), velIsJ);
}

/* *************** Velocity Norm ************************************* */

template <typename T, template <typename U> class Descriptor>
void computeVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &velocityNorm, Box3D domain)
{
    applyProcessingFunctional(
        new BoxVelocityNormFunctional3D<T, Descriptor>, domain, lattice, velocityNorm);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > velocityNorm =
        generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxVelocityNormFunctional3D() acts on both bulk and envelope, you would expect the envelope
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
std::unique_ptr<MultiScalarField3D<T> > computeVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeVelocityNorm(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void computeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<T> &velocityNorm, Box3D domain)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&lattice);
    args.push_back(&force);
    args.push_back(&velocityNorm);
    applyProcessingFunctional(new BoxForcedVelocityNormFunctional3D<T, Descriptor>, domain, args);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > velocityNorm =
        generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxForcedVelocityNormFunctional3D() acts on both bulk and envelope, you would expect the
    // envelope layer around the domain, on the velocityNorm multi-block, to be assigned some proper
    // values too. By default, this is however not what happens, because the physical space occupied
    // by these envelopes does not intersect with the domain "domain". We work around this issue by
    // extending the domain. There's no problem if the enlarged domain gets beyond the actual extent
    // of the lattice, because Palabos handles these situations properly.

    computeForcedVelocityNorm(
        lattice, force, *velocityNorm,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return velocityNorm;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force)
{
    return computeForcedVelocityNorm(lattice, force, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void computeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force,
    MultiScalarField3D<T> &velocityNorm, Box3D domain)
{
    applyProcessingFunctional(
        new BoxConstForcedVelocityNormFunctional3D<T, Descriptor>(force), domain, lattice,
        velocityNorm);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > velocityNorm =
        generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxConstForcedVelocityNormFunctional3D() acts on both bulk and envelope, you would expect the
    // envelope layer around the domain, on the velocityNorm multi-block, to be assigned some proper
    // values too. By default, this is however not what happens, because the physical space occupied
    // by these envelopes does not intersect with the domain "domain". We work around this issue by
    // extending the domain. There's no problem if the enlarged domain gets beyond the actual extent
    // of the lattice, because Palabos handles these situations properly.

    computeForcedVelocityNorm(
        lattice, force, *velocityNorm,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return velocityNorm;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force)
{
    return computeForcedVelocityNorm(lattice, force, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
void computeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, ForceFunction f,
    MultiScalarField3D<T> &velocityNorm, Box3D domain)
{
    applyProcessingFunctional(
        new BoxCustomForcedVelocityNormFunctional3D<T, Descriptor, ForceFunction>(f), domain,
        lattice, velocityNorm);
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, ForceFunction f, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > velocityNorm =
        generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxCustomForcedVelocityNormFunctional3D() acts on both bulk and envelope, you would expect
    // the envelope layer around the domain, on the velocityNorm multi-block, to be assigned some
    // proper values too. By default, this is however not what happens, because the physical space
    // occupied by these envelopes does not intersect with the domain "domain". We work around this
    // issue by extending the domain. There's no problem if the enlarged domain gets beyond the
    // actual extent of the lattice, because Palabos handles these situations properly.

    computeForcedVelocityNorm(
        lattice, f, *velocityNorm,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return velocityNorm;
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, ForceFunction f)
{
    return computeForcedVelocityNorm(lattice, f, lattice.getBoundingBox());
}

/* *************** Velocity Component ******************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &velocityComponent,
    Box3D domain, plint iComponent)
{
    applyProcessingFunctional(
        new BoxVelocityComponentFunctional3D<T, Descriptor>(iComponent), domain, lattice,
        velocityComponent);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, plint iComponent)
{
    std::unique_ptr<MultiScalarField3D<T> > velocityComponent =
        generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxVelocityComponentFunctional3D() acts on both bulk and envelope, you would expect the
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
std::unique_ptr<MultiScalarField3D<T> > computeVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, plint iComponent)
{
    return computeVelocityComponent(lattice, lattice.getBoundingBox(), iComponent);
}

template <typename T, template <typename U> class Descriptor>
void computeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<T> &velocityComponent, Box3D domain, plint iComponent)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&lattice);
    args.push_back(&force);
    args.push_back(&velocityComponent);
    applyProcessingFunctional(
        new BoxForcedVelocityComponentFunctional3D<T, Descriptor>(iComponent), domain, args);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    Box3D domain, plint iComponent)
{
    std::unique_ptr<MultiScalarField3D<T> > velocityComponent =
        generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxForcedVelocityComponentFunctional3D() acts on both bulk and envelope, you would expect the
    // envelope layer around the domain, on the velocityComponent multi-block, to be assigned some
    // proper values too. By default, this is however not what happens, because the physical space
    // occupied by these envelopes does not intersect with the domain "domain". We work around this
    // issue by extending the domain. There's no problem if the enlarged domain gets beyond the
    // actual extent of the lattice, because Palabos handles these situations properly.

    computeForcedVelocityComponent(
        lattice, force, *velocityComponent,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()), iComponent);
    return velocityComponent;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    plint iComponent)
{
    return computeForcedVelocityComponent(lattice, force, lattice.getBoundingBox(), iComponent);
}

template <typename T, template <typename U> class Descriptor>
void computeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force,
    MultiScalarField3D<T> &velocityComponent, Box3D domain, plint iComponent)
{
    applyProcessingFunctional(
        new BoxConstForcedVelocityComponentFunctional3D<T, Descriptor>(force, iComponent), domain,
        lattice, velocityComponent);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force, Box3D domain,
    plint iComponent)
{
    std::unique_ptr<MultiScalarField3D<T> > velocityComponent =
        generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxConstForcedVelocityComponentFunctional3D() acts on both bulk and envelope, you would
    // expect the envelope layer around the domain, on the velocityComponent multi-block, to be
    // assigned some proper values too. By default, this is however not what happens, because the
    // physical space occupied by these envelopes does not intersect with the domain "domain". We
    // work around this issue by extending the domain. There's no problem if the enlarged domain
    // gets beyond the actual extent of the lattice, because Palabos handles these situations
    // properly.

    computeForcedVelocityComponent(
        lattice, force, *velocityComponent,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()), iComponent);
    return velocityComponent;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force, plint iComponent)
{
    return computeForcedVelocityComponent(lattice, force, lattice.getBoundingBox(), iComponent);
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
void computeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, ForceFunction f,
    MultiScalarField3D<T> &velocityComponent, Box3D domain, plint iComponent)
{
    applyProcessingFunctional(
        new BoxCustomForcedVelocityComponentFunctional3D<T, Descriptor, ForceFunction>(
            f, iComponent),
        domain, lattice, velocityComponent);
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, ForceFunction f, Box3D domain, plint iComponent)
{
    std::unique_ptr<MultiScalarField3D<T> > velocityComponent =
        generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxCustomForcedVelocityComponentFunctional3D() acts on both bulk and envelope, you would
    // expect the envelope layer around the domain, on the velocityComponent multi-block, to be
    // assigned some proper values too. By default, this is however not what happens, because the
    // physical space occupied by these envelopes does not intersect with the domain "domain". We
    // work around this issue by extending the domain. There's no problem if the enlarged domain
    // gets beyond the actual extent of the lattice, because Palabos handles these situations
    // properly.

    computeForcedVelocityComponent(
        lattice, f, *velocityComponent,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()), iComponent);
    return velocityComponent;
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
std::unique_ptr<MultiScalarField3D<T> > computeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, ForceFunction f, plint iComponent)
{
    return computeForcedVelocityComponent(lattice, f, lattice.getBoundingBox(), iComponent);
}

/* *************** Velocity ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &velocity,
    Box3D domain)
{
    applyProcessingFunctional(
        new BoxVelocityFunctional3D<T, Descriptor>, domain, lattice, velocity);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > velocity =
        generateMultiTensorField<T, Descriptor<T>::d>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the BoxVelocityFunctional3D()
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
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeVelocity(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void computeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiTensorField3D<T, Descriptor<T>::d> &velocity, Box3D domain)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&lattice);
    args.push_back(&force);
    args.push_back(&velocity);
    applyProcessingFunctional(new BoxForcedVelocityFunctional3D<T, Descriptor>, domain, args);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > velocity =
        generateMultiTensorField<T, Descriptor<T>::d>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxForcedVelocityFunctional3D() acts on both bulk and envelope, you would expect the envelope
    // layer around the domain, on the velocity multi-block, to be assigned some proper values too.
    // By default, this is however not what happens, because the physical space occupied by these
    // envelopes does not intersect with the domain "domain". We work around this issue by extending
    // the domain. There's no problem if the enlarged domain gets beyond the actual extent of the
    // lattice, because Palabos handles these situations properly.

    computeForcedVelocity(
        lattice, force, *velocity,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return velocity;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force)
{
    return computeForcedVelocity(lattice, force, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void computeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force,
    MultiTensorField3D<T, Descriptor<T>::d> &velocity, Box3D domain)
{
    applyProcessingFunctional(
        new BoxConstForcedVelocityFunctional3D<T, Descriptor>(force), domain, lattice, velocity);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > velocity =
        generateMultiTensorField<T, Descriptor<T>::d>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxConstForcedVelocityFunctional3D() acts on both bulk and envelope, you would expect the
    // envelope layer around the domain, on the velocity multi-block, to be assigned some proper
    // values too. By default, this is however not what happens, because the physical space occupied
    // by these envelopes does not intersect with the domain "domain". We work around this issue by
    // extending the domain. There's no problem if the enlarged domain gets beyond the actual extent
    // of the lattice, because Palabos handles these situations properly.

    computeForcedVelocity(
        lattice, force, *velocity,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return velocity;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> force)
{
    return computeForcedVelocity(lattice, force, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
void computeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, ForceFunction f,
    MultiTensorField3D<T, Descriptor<T>::d> &velocity, Box3D domain)
{
    applyProcessingFunctional(
        new BoxCustomForcedVelocityFunctional3D<T, Descriptor, ForceFunction>(f), domain, lattice,
        velocity);
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, ForceFunction f, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > velocity =
        generateMultiTensorField<T, Descriptor<T>::d>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxCustomForcedVelocityFunctional3D() acts on both bulk and envelope, you would expect the
    // envelope layer around the domain, on the velocity multi-block, to be assigned some proper
    // values too. By default, this is however not what happens, because the physical space occupied
    // by these envelopes does not intersect with the domain "domain". We work around this issue by
    // extending the domain. There's no problem if the enlarged domain gets beyond the actual extent
    // of the lattice, because Palabos handles these situations properly.

    computeForcedVelocity(
        lattice, f, *velocity,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return velocity;
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, ForceFunction f)
{
    return computeForcedVelocity(lattice, f, lattice.getBoundingBox());
}

/* *************** Temperature ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeTemperature(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &temperature, Box3D domain)
{
    applyProcessingFunctional(
        new BoxTemperatureFunctional3D<T, Descriptor>, domain, lattice, temperature);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeTemperature(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > temperature =
        generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the BoxTemperatureFunctional3D()
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
std::unique_ptr<MultiScalarField3D<T> > computeTemperature(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeTemperature(lattice, lattice.getBoundingBox());
}

/* *************** Pi Neq ********************************* */

template <typename T, template <typename U> class Descriptor>
void computePiNeq(
    MultiBlockLattice3D<T, Descriptor> &lattice,
    MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq, Box3D domain)
{
    applyProcessingFunctional(new BoxPiNeqFunctional3D<T, Descriptor>, domain, lattice, PiNeq);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computePiNeq(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > PiNeq =
        generateMultiTensorField<T, SymmetricTensor<T, Descriptor>::n>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the BoxPiNeqFunctional3D() acts
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
std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computePiNeq(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computePiNeq(lattice, lattice.getBoundingBox());
}

/* *************** Shear Stress ********************************* */

template <typename T, template <typename U> class Descriptor>
void computeShearStress(
    MultiBlockLattice3D<T, Descriptor> &lattice,
    MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq, Box3D domain)
{
    applyProcessingFunctional(
        new BoxShearStressFunctional3D<T, Descriptor>, domain, lattice, PiNeq);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computeShearStress(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > PiNeq =
        generateMultiTensorField<T, SymmetricTensor<T, Descriptor>::n>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the BoxShearStressFunctional3D()
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
std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computeShearStress(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeShearStress(lattice, lattice.getBoundingBox());
}

/* *************** Stress ********************************* */

template <typename T, template <typename U> class Descriptor>
void computeStress(
    MultiBlockLattice3D<T, Descriptor> &lattice,
    MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> &stress, T rho0, bool isCompressible,
    Box3D domain)
{
    applyProcessingFunctional(
        new BoxStressFunctional3D<T, Descriptor>(rho0, isCompressible), domain, lattice, stress);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computeStress(
    MultiBlockLattice3D<T, Descriptor> &lattice, T rho0, bool isCompressible, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > stress =
        generateMultiTensorField<T, SymmetricTensor<T, Descriptor>::n>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the BoxShearStressFunctional3D()
    // acts on both bulk and envelope, you would expect the envelope layer around the domain, on the
    // PiNeq multi-block, to be assigned some proper values too. By default, this is however not
    // what happens, because the physical space occupied by these envelopes does not intersect with
    // the domain "domain". We work around this issue by extending the domain. There's no problem if
    // the enlarged domain gets beyond the actual extent of the lattice, because Palabos handles
    // these situations properly.

    computeStress(
        lattice, *stress, rho0, isCompressible,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return stress;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computeStress(
    MultiBlockLattice3D<T, Descriptor> &lattice, T rho0, bool isCompressible)
{
    return computeStress(lattice, rho0, isCompressible, lattice.getBoundingBox());
}

/* *************** Strain Rate from Stress *************************** */

template <typename T, template <typename U> class Descriptor>
void computeStrainRateFromStress(
    MultiBlockLattice3D<T, Descriptor> &lattice,
    MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> &S, Box3D domain)
{
    applyProcessingFunctional(
        new BoxStrainRateFromStressFunctional3D<T, Descriptor>, domain, lattice, S);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> >
    computeStrainRateFromStress(MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > S =
        generateMultiTensorField<T, SymmetricTensor<T, Descriptor>::n>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxStrainRateFromStressFunctional3D() acts on both bulk and envelope, you would expect the
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
std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> >
    computeStrainRateFromStress(MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeStrainRateFromStress(lattice, lattice.getBoundingBox());
}

/* *************** Shear Rate ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeShearRate(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &shearRate, Box3D domain)
{
    applyProcessingFunctional(
        new BoxShearRateFunctional3D<T, Descriptor>, domain, lattice, shearRate);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeShearRate(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > shearRate =
        generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the BoxShearRateFunctional3D()
    // acts on both bulk and envelope, you would expect the envelope layer around the domain, on the
    // shearRate multi-block, to be assigned some proper values too. By default, this is however not
    // what happens, because the physical space occupied by these envelopes does not intersect with
    // the domain "domain". We work around this issue by extending the domain. There's no problem if
    // the enlarged domain gets beyond the actual extent of the lattice, because Palabos handles
    // these situations properly.

    computeShearRate(
        lattice, *shearRate, domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return shearRate;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeShearRate(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeShearRate(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void computeShearRate_N(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiNTensorField3D<T> &shearRate, Box3D domain)
{
    applyProcessingFunctional(
        new BoxNTensorShearRateFunctional3D<T, Descriptor>, domain, lattice, shearRate);
}

template <typename T, template <typename U> class Descriptor>
MultiNTensorField3D<T> *computeShearRate_N(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    MultiNTensorField3D<T> *shearRate = generateMultiNTensorField<T>(lattice, domain, 1);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxNTensorShearRateFunctional3D() acts on both bulk and envelope, you would expect the
    // envelope layer around the domain, on the shearRate multi-block, to be assigned some proper
    // values too. By default, this is however not what happens, because the physical space occupied
    // by these envelopes does not intersect with the domain "domain". We work around this issue by
    // extending the domain. There's no problem if the enlarged domain gets beyond the actual extent
    // of the lattice, because Palabos handles these situations properly.

    computeShearRate_N(
        lattice, *shearRate, domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return shearRate;
}

/* *************** Population **************************************** */

template <typename T, template <typename U> class Descriptor>
void computePopulation(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &population, Box3D domain,
    plint iPop)
{
    applyProcessingFunctional(
        new BoxPopulationFunctional3D<T, Descriptor>(iPop), domain, lattice, population);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computePopulation(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, plint iPop)
{
    std::unique_ptr<MultiScalarField3D<T> > population =
        generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the BoxPopulationFunctional3D()
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
std::unique_ptr<MultiScalarField3D<T> > computePopulation(
    MultiBlockLattice3D<T, Descriptor> &lattice, plint iPop)
{
    return computePopulation(lattice, lattice.getBoundingBox(), iPop);
}

template <typename T, template <typename U> class Descriptor>
void computeEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &equilibrium, Box3D domain,
    plint iPop)
{
    applyProcessingFunctional(
        new BoxEquilibriumFunctional3D<T, Descriptor>(iPop), domain, lattice, equilibrium);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, plint iPop)
{
    std::unique_ptr<MultiScalarField3D<T> > equilibrium =
        generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the BoxEquilibriumFunctional3D()
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
std::unique_ptr<MultiScalarField3D<T> > computeEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice, plint iPop)
{
    return computeEquilibrium(lattice, lattice.getBoundingBox(), iPop);
}

template <typename T, template <typename U> class Descriptor>
void computeEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice,
    MultiTensorField3D<T, Descriptor<T>::q> &equilibrium, Box3D domain)
{
    applyProcessingFunctional(
        new BoxAllEquilibriumFunctional3D<T, Descriptor>(), domain, lattice, equilibrium);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::q> > computeEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::q> > equilibrium =
        generateMultiTensorField<T, Descriptor<T>::q>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxAllEquilibriumFunctional3D() acts on both bulk and envelope, you would expect the envelope
    // layer around the domain, on the equilibrium multi-block, to be assigned some proper values
    // too. By default, this is however not what happens, because the physical space occupied by
    // these envelopes does not intersect with the domain "domain". We work around this issue by
    // extending the domain. There's no problem if the enlarged domain gets beyond the actual extent
    // of the lattice, because Palabos handles these situations properly.

    computeEquilibrium(
        lattice, *equilibrium,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return equilibrium;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::q> > computeEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeEquilibrium(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void computeNonEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice,
    MultiTensorField3D<T, Descriptor<T>::q> &nonEquilibrium, Box3D domain)
{
    applyProcessingFunctional(
        new BoxAllNonEquilibriumFunctional3D<T, Descriptor>(), domain, lattice, nonEquilibrium);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::q> > computeNonEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::q> > nonEquilibrium =
        generateMultiTensorField<T, Descriptor<T>::q>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxAllNonEquilibriumFunctional3D() acts on both bulk and envelope, you would expect the
    // envelope layer around the domain, on the nonEquilibrium multi-block, to be assigned some
    // proper values too. By default, this is however not what happens, because the physical space
    // occupied by these envelopes does not intersect with the domain "domain". We work around this
    // issue by extending the domain. There's no problem if the enlarged domain gets beyond the
    // actual extent of the lattice, because Palabos handles these situations properly.

    computeNonEquilibrium(
        lattice, *nonEquilibrium,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return nonEquilibrium;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::q> > computeNonEquilibrium(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeNonEquilibrium(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void computeAllPopulations(
    MultiBlockLattice3D<T, Descriptor> &lattice,
    MultiTensorField3D<T, Descriptor<T>::q> &populations, Box3D domain)
{
    applyProcessingFunctional(
        new BoxAllPopulationsFunctional3D<T, Descriptor>(), domain, lattice, populations);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::q> > computeAllPopulations(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::q> > populations =
        generateMultiTensorField<T, Descriptor<T>::q>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxAllPopulationsFunctional3D() acts on both bulk and envelope, you would expect the envelope
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
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::q> > computeAllPopulations(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeAllPopulations(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void computeAllPopulationsFromTensorField(
    MultiBlockLattice3D<T, Descriptor> &lattice,
    MultiTensorField3D<T, Descriptor<T>::q> &populations, Box3D domain)
{
    applyProcessingFunctional(
        new BoxAllPopulationsToLatticeFunctional3D<T, Descriptor>(), domain, lattice, populations);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::q> > computeAllPopulationsFromTensorField(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::q> > populations =
        generateMultiTensorField<T, Descriptor<T>::q>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxAllPopulationsToLatticeFunctional3D() acts on both bulk and envelope, you would expect the
    // envelope layer around the domain, on the populations multi-block, to be assigned some proper
    // values too. By default, this is however not what happens, because the physical space occupied
    // by these envelopes does not intersect with the domain "domain". We work around this issue by
    // extending the domain. There's no problem if the enlarged domain gets beyond the actual extent
    // of the lattice, because Palabos handles these situations properly.

    computeAllPopulationsFromTensorField(
        lattice, *populations,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return populations;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::q> > computeAllPopulationsFromTensorField(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeAllPopulationsFromTensorField(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void copyPopulations(
    MultiBlockLattice3D<T, Descriptor> &latticeFrom, MultiBlockLattice3D<T, Descriptor> &latticeTo,
    Box3D domain)
{
    applyProcessingFunctional(
        new CopyPopulationsFunctional3D<T, Descriptor>(), domain, latticeFrom, latticeTo);
}

template <typename T, template <typename U> class Descriptor>
void copyAll(
    MultiBlockLattice3D<T, Descriptor> &latticeFrom, MultiBlockLattice3D<T, Descriptor> &latticeTo,
    Box3D domain)
{
    applyProcessingFunctional(
        new LatticeCopyAllFunctional3D<T, Descriptor>(), domain, latticeFrom, latticeTo);
}

template <typename T, template <typename U> class Descriptor>
void copyRegenerate(
    MultiBlockLattice3D<T, Descriptor> &latticeFrom, MultiBlockLattice3D<T, Descriptor> &latticeTo,
    Box3D domain)
{
    applyProcessingFunctional(
        new LatticeRegenerateFunctional3D<T, Descriptor>(), domain, latticeFrom, latticeTo);
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void copyConvertPopulations(
    MultiBlockLattice3D<T1, Descriptor1> &latticeFrom,
    MultiBlockLattice3D<T2, Descriptor2> &latticeTo, Box3D domain)
{
    applyProcessingFunctional(
        new CopyConvertPopulationsFunctional3D<T1, Descriptor1, T2, Descriptor2>(), domain,
        latticeFrom, latticeTo);
}

/* *************** Omega ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeOmega(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &omega, Box3D domain)
{
    applyProcessingFunctional(new BoxOmegaFunctional3D<T, Descriptor>, domain, lattice, omega);
}

template <typename T, template <typename U> class Descriptor>
void computeOmega(MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &omega)
{
    computeOmega(lattice, omega, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeOmega(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > omega = generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the BoxOmegaFunctional3D() acts
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
std::unique_ptr<MultiScalarField3D<T> > computeOmega(MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeOmega(lattice, lattice.getBoundingBox());
}

/* *************** KinematicViscosity ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeKinematicViscosity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &omega, Box3D domain)
{
    applyProcessingFunctional(
        new BoxKinematicViscosityFunctional3D<T, Descriptor>, domain, lattice, omega);
}

template <typename T, template <typename U> class Descriptor>
void computeKinematicViscosity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &omega)
{
    computeKinematicViscosity(lattice, omega, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeKinematicViscosity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > omega = generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxKinematicViscosityFunctional3D() acts on both bulk and envelope, you would expect the
    // envelope layer around the domain, on the omega multi-block, to be assigned some proper values
    // too. By default, this is however not what happens, because the physical space occupied by
    // these envelopes does not intersect with the domain "domain". We work around this issue by
    // extending the domain. There's no problem if the enlarged domain gets beyond the actual extent
    // of the lattice, because Palabos handles these situations properly.

    computeKinematicViscosity(
        lattice, *omega, domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return omega;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeKinematicViscosity(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeKinematicViscosity(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void computeKinematicViscosity_N(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiNTensorField3D<T> &nu, Box3D domain)
{
    applyProcessingFunctional(
        new BoxNTensorKinematicViscosityFunctional3D<T, Descriptor>, domain, lattice, nu);
}

template <typename T, template <typename U> class Descriptor>
MultiNTensorField3D<T> *computeKinematicViscosity_N(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    MultiNTensorField3D<T> *nu = generateMultiNTensorField<T>(lattice, domain, 1);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxNTensorKinematicViscosityFunctional3D() acts on both bulk and envelope, you would expect
    // the envelope layer around the domain, on the omega multi-block, to be assigned some proper
    // values too. By default, this is however not what happens, because the physical space occupied
    // by these envelopes does not intersect with the domain "domain". We work around this issue by
    // extending the domain. There's no problem if the enlarged domain gets beyond the actual extent
    // of the lattice, because Palabos handles these situations properly.

    computeKinematicViscosity_N(
        lattice, *nu, domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return nu;
}

template <typename T, template <typename U> class Descriptor>
MultiNTensorField3D<T> *computeKinematicViscosity_N(MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeKinematicViscosity_N(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void computeKinematicEddyViscosity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &omega, Box3D domain)
{
    applyProcessingFunctional(
        new BoxKinematicEddyViscosityFunctional3D<T, Descriptor>, domain, lattice, omega);
}

template <typename T, template <typename U> class Descriptor>
void computeKinematicEddyViscosity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &omega)
{
    computeKinematicEddyViscosity(lattice, omega, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeKinematicEddyViscosity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > omega = generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxKinematicEddyViscosityFunctional3D() acts on both bulk and envelope, you would expect the
    // envelope layer around the domain, on the omega multi-block, to be assigned some proper values
    // too. By default, this is however not what happens, because the physical space occupied by
    // these envelopes does not intersect with the domain "domain". We work around this issue by
    // extending the domain. There's no problem if the enlarged domain gets beyond the actual extent
    // of the lattice, because Palabos handles these situations properly.

    computeKinematicEddyViscosity(
        lattice, *omega, domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return omega;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeKinematicEddyViscosity(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeKinematicEddyViscosity(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void computeKinematicEddyViscosity_N(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiNTensorField3D<T> &nu, Box3D domain)
{
    applyProcessingFunctional(
        new BoxNTensorKinematicEddyViscosityFunctional3D<T, Descriptor>, domain, lattice, nu);
}

template <typename T, template <typename U> class Descriptor>
MultiNTensorField3D<T> *computeKinematicEddyViscosity_N(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    MultiNTensorField3D<T> *nu = generateMultiNTensorField<T>(lattice, domain, 1);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxNTensorKinematicEddyViscosityFunctional3D() acts on both bulk and envelope, you would
    // expect the envelope layer around the domain, on the omega multi-block, to be assigned some
    // proper values too. By default, this is however not what happens, because the physical space
    // occupied by these envelopes does not intersect with the domain "domain". We work around this
    // issue by extending the domain. There's no problem if the enlarged domain gets beyond the
    // actual extent of the lattice, because Palabos handles these situations properly.

    computeKinematicEddyViscosity_N(
        lattice, *nu, domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return nu;
}

template <typename T, template <typename U> class Descriptor>
MultiNTensorField3D<T> *computeKinematicEddyViscosity_N(MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeKinematicEddyViscosity_N(lattice, lattice.getBoundingBox());
}

/* *************** ExternalForce ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeExternalForce(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    Box3D domain)
{
    applyProcessingFunctional(
        new BoxExternalForceFunctional3D<T, Descriptor>, domain, lattice, force);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeExternalForce(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > force =
        generateMultiTensorField<T, Descriptor<T>::d>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxExternalForceFunctional3D() acts on both bulk and envelope, you would expect the envelope
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
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeExternalForce(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeExternalForce(lattice, lattice.getBoundingBox());
}

/* *************** ExternalScalar ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeExternalScalar(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &scalar, int whichScalar,
    Box3D domain)
{
    applyProcessingFunctional(
        new BoxExternalScalarFunctional3D<T, Descriptor>(whichScalar), domain, lattice, scalar);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeExternalScalar(
    MultiBlockLattice3D<T, Descriptor> &lattice, int whichScalar, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > scalar = generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxExternalScalarFunctional3D() acts on both bulk and envelope, you would expect the envelope
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
std::unique_ptr<MultiScalarField3D<T> > computeExternalScalar(
    MultiBlockLattice3D<T, Descriptor> &lattice, int whichScalar)
{
    return computeExternalScalar(lattice, whichScalar, lattice.getBoundingBox());
}

/* *************** ExternalVector ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeExternalVector(
    MultiBlockLattice3D<T, Descriptor> &lattice,
    MultiTensorField3D<T, Descriptor<T>::d> &tensorField, int vectorBeginsAt, Box3D domain)
{
    applyProcessingFunctional(
        new BoxExternalVectorFunctional3D<T, Descriptor>(vectorBeginsAt), domain, lattice,
        tensorField);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeExternalVector(
    MultiBlockLattice3D<T, Descriptor> &lattice, int vectorBeginsAt, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > tensorField =
        generateMultiTensorField<T, Descriptor<T>::d>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxExternalVectorFunctional3D() acts on both bulk and envelope, you would expect the envelope
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
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeExternalVector(
    MultiBlockLattice3D<T, Descriptor> &lattice, int vectorBeginsAt)
{
    return computeExternalVector(lattice, vectorBeginsAt, lattice.getBoundingBox());
}

/* *************** DynamicParameter ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeDynamicParameter(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &scalar,
    plint whichParameter, Box3D domain)
{
    applyProcessingFunctional(
        new BoxDynamicParameterFunctional3D<T, Descriptor>(whichParameter), domain, lattice,
        scalar);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeDynamicParameter(
    MultiBlockLattice3D<T, Descriptor> &lattice, plint whichParameter, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > scalar = generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxDynamicParameterFunctional3D() acts on both bulk and envelope, you would expect the
    // envelope layer around the domain, on the scalar multi-block, to be assigned some proper
    // values too. By default, this is however not what happens, because the physical space occupied
    // by these envelopes does not intersect with the domain "domain". We work around this issue by
    // extending the domain. There's no problem if the enlarged domain gets beyond the actual extent
    // of the lattice, because Palabos handles these situations properly.

    computeDynamicParameter(
        lattice, *scalar, whichParameter,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return scalar;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeDynamicParameter(
    MultiBlockLattice3D<T, Descriptor> &lattice, plint whichParameter)
{
    return computeDynamicParameter(lattice, whichParameter, lattice.getBoundingBox());
}

/* *************** DynamicViscosity ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeDynamicViscosity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &scalar, Box3D domain)
{
    applyProcessingFunctional(
        new BoxDynamicViscosityFunctional3D<T, Descriptor>(), domain, lattice, scalar);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeDynamicViscosity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > scalar = generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxDynamicViscosityFunctional3D() acts on both bulk and envelope, you would expect the
    // envelope layer around the domain, on the scalar multi-block, to be assigned some proper
    // values too. By default, this is however not what happens, because the physical space occupied
    // by these envelopes does not intersect with the domain "domain". We work around this issue by
    // extending the domain. There's no problem if the enlarged domain gets beyond the actual extent
    // of the lattice, because Palabos handles these situations properly.

    computeDynamicViscosity(
        lattice, *scalar, domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return scalar;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeDynamicViscosity(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeDynamicViscosity(lattice, lattice.getBoundingBox());
}

/* ******************************************************************* */
/* *************** PART V. Multi-block wrappers: Scalar-Field ******** */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */
template <typename T>
T computeSum(MultiScalarField3D<T> &scalarField, Box3D domain)
{
    BoxScalarSumFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getSumScalar();
}

template <typename T>
T computeSum(MultiScalarField3D<T> &scalarField)
{
    return computeSum(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeSum(
    MultiScalarField3D<T> &scalarField, MultiScalarField3D<int> &mask, int flag, Box3D domain)
{
    MaskedBoxScalarSumFunctional3D<T> functional(flag);
    applyProcessingFunctional(functional, domain, scalarField, mask);
    return functional.getSumScalar();
}

template <typename T>
T computeSum(MultiScalarField3D<T> &scalarField, MultiScalarField3D<int> &mask, int flag)
{
    return computeSum(scalarField, mask, flag, scalarField.getBoundingBox());
}

template <typename T>
plint computeIntSum(MultiScalarField3D<T> &scalarField, Box3D domain)
{
    BoxScalarIntSumFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getSumScalar();
}

template <typename T>
T computeBoundedSum(MultiScalarField3D<T> &scalarField, Box3D domain)
{
    BoundedBoxScalarSumFunctional3D<T> functional;
    plint envelopeWidth = 1;
    applyProcessingFunctional(functional, domain, scalarField, envelopeWidth);
    return functional.getSumScalar();
}

template <typename T>
T computeBoundedSum(MultiScalarField3D<T> &scalarField)
{
    return computeBoundedSum(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeAverage(MultiScalarField3D<T> &scalarField, Box3D domain)
{
    BoxScalarSumFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getSumScalar() / (T)domain.nCells();
}

template <typename T>
T computeAverage(MultiScalarField3D<T> &scalarField)
{
    return computeAverage(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeAverage(
    MultiScalarField3D<T> &scalarField, MultiScalarField3D<int> &mask, int flag, Box3D domain)
{
    MaskedBoxScalarAverageFunctional3D<T> functional(flag);
    applyProcessingFunctional(functional, domain, scalarField, mask);
    return functional.getAverageScalar();
}

template <typename T>
T computeAverage(MultiScalarField3D<T> &scalarField, MultiScalarField3D<int> &mask, int flag)
{
    return computeAverage(scalarField, mask, flag, scalarField.getBoundingBox());
}

template <typename T>
T computeMin(MultiScalarField3D<T> &scalarField, Box3D domain)
{
    BoxScalarMinFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getMinScalar();
}

template <typename T>
T computeMin(MultiScalarField3D<T> &scalarField)
{
    return computeMin(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeMin(
    MultiScalarField3D<T> &scalarField, MultiScalarField3D<int> &mask, int flag, Box3D domain)
{
    MaskedBoxScalarMinFunctional3D<T> functional(flag);
    applyProcessingFunctional(functional, domain, scalarField, mask);
    return functional.getMinScalar();
}

template <typename T>
T computeMin(MultiScalarField3D<T> &scalarField, MultiScalarField3D<int> &mask, int flag)
{
    return computeMin(scalarField, mask, flag, scalarField.getBoundingBox());
}

template <typename T>
T computeMax(MultiScalarField3D<T> &scalarField, Box3D domain)
{
    BoxScalarMaxFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getMaxScalar();
}

template <typename T>
T computeMax(MultiScalarField3D<T> &scalarField)
{
    return computeMax(scalarField, scalarField.getBoundingBox());
}

template <typename T>
T computeMax(
    MultiScalarField3D<T> &scalarField, MultiScalarField3D<int> &mask, int flag, Box3D domain)
{
    MaskedBoxScalarMaxFunctional3D<T> functional(flag);
    applyProcessingFunctional(functional, domain, scalarField, mask);
    return functional.getMaxScalar();
}

template <typename T>
T computeMax(MultiScalarField3D<T> &scalarField, MultiScalarField3D<int> &mask, int flag)
{
    return computeMax(scalarField, mask, flag, scalarField.getBoundingBox());
}

template <typename T>
T computeBoundedAverage(MultiScalarField3D<T> &scalarField, Box3D domain)
{
    BoundedBoxScalarSumFunctional3D<T> functional;
    plint envelopeWidth = 1;
    applyProcessingFunctional(functional, domain, scalarField, envelopeWidth);
    return functional.getSumScalar()
           / (T)((domain.getNx() - 1) * (domain.getNy() - 1) * (domain.getNz() - 1));
}

template <typename T>
T computeBoundedAverage(MultiScalarField3D<T> &scalarField)
{
    return computeBoundedAverage(scalarField, scalarField.getBoundingBox());
}

template <typename T, class BoolMask>
plint count(MultiScalarField3D<T> &field, Box3D domain, BoolMask boolMask)
{
    CountScalarElementsFunctional3D<T, BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, field);
    return functional.getCount();
}

template <typename T, class BoolMask>
plint count(MultiScalarField3D<T> &field, BoolMask boolMask)
{
    return count(field, field.getBoundingBox(), boolMask);
}

/* *************** Extract Sub-ScalarField *************************** */

template <typename T>
void extractSubDomain(
    MultiScalarField3D<T> &field, MultiScalarField3D<T> &extractedField, Box3D domain)
{
    applyProcessingFunctional(
        new ExtractScalarSubDomainFunctional3D<T>, domain, field, extractedField);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > extractSubDomain(MultiScalarField3D<T> &field, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > extractedField =
        generateMultiScalarField<T>(field, domain);
    extractSubDomain(field, *extractedField, domain);
    return extractedField;
}

/* *************** MultiScalarField - Scalar operations *************** */

template <typename T>
void lessThan(MultiScalarField3D<T> &field, T scalar, MultiScalarField3D<int> &result, Box3D domain)
{
    applyProcessingFunctional(new A_lt_alpha_functional3D<T>(scalar), domain, field, result);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<int> > lessThan(
    MultiScalarField3D<T> &field, T scalar, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<int> > result = generateMultiScalarField<int>(field, domain);
    lessThan(field, scalar, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<int> > lessThan(MultiScalarField3D<T> &field, T scalar)
{
    return lessThan(field, scalar, field.getBoundingBox());
}

template <typename T>
void greaterThan(
    MultiScalarField3D<T> &field, T scalar, MultiScalarField3D<int> &result, Box3D domain)
{
    applyProcessingFunctional(new A_gt_alpha_functional3D<T>(scalar), domain, field, result);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<int> > greaterThan(
    MultiScalarField3D<T> &field, T scalar, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<int> > result = generateMultiScalarField<int>(field, domain);
    greaterThan(field, scalar, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<int> > greaterThan(MultiScalarField3D<T> &field, T scalar)
{
    return greaterThan(field, scalar, field.getBoundingBox());
}

template <typename T>
void add(MultiScalarField3D<T> &field, T scalar, MultiScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(new A_plus_alpha_functional3D<T>(scalar), domain, field, result);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > add(MultiScalarField3D<T> &field, T scalar, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > result = generateMultiScalarField<T>(field, domain);
    add(field, scalar, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > add(MultiScalarField3D<T> &field, T scalar)
{
    return add(field, scalar, field.getBoundingBox());
}

template <typename T>
void add(T scalar, MultiScalarField3D<T> &field, MultiScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(new A_plus_alpha_functional3D<T>(scalar), domain, field, result);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > add(T scalar, MultiScalarField3D<T> &field, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > result = generateMultiScalarField<T>(field, domain);
    add(scalar, field, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > add(T scalar, MultiScalarField3D<T> &field)
{
    return add(scalar, field, field.getBoundingBox());
}

template <typename T>
void subtract(MultiScalarField3D<T> &field, T scalar, MultiScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(new A_minus_alpha_functional3D<T>(scalar), domain, field, result);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > subtract(
    MultiScalarField3D<T> &field, T scalar, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > result = generateMultiScalarField<T>(field, domain);
    subtract(field, scalar, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > subtract(MultiScalarField3D<T> &field, T scalar)
{
    return subtract(field, scalar, field.getBoundingBox());
}

template <typename T>
void subtract(T scalar, MultiScalarField3D<T> &field, MultiScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(new Alpha_minus_A_functional3D<T>(scalar), domain, field, result);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > subtract(
    T scalar, MultiScalarField3D<T> &field, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > result = generateMultiScalarField<T>(field, domain);
    subtract(scalar, field, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > subtract(T scalar, MultiScalarField3D<T> &field)
{
    return subtract(scalar, field, field.getBoundingBox());
}

template <typename T>
void multiply(MultiScalarField3D<T> &field, T scalar, MultiScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(new A_times_alpha_functional3D<T>(scalar), domain, field, result);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > multiply(
    MultiScalarField3D<T> &field, T scalar, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > result = generateMultiScalarField<T>(field, domain);
    multiply(field, scalar, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > multiply(MultiScalarField3D<T> &field, T scalar)
{
    return multiply(field, scalar, field.getBoundingBox());
}

template <typename T>
void multiply(T scalar, MultiScalarField3D<T> &field, MultiScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(new A_times_alpha_functional3D<T>(scalar), domain, field, result);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > multiply(
    T scalar, MultiScalarField3D<T> &field, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > result = generateMultiScalarField<T>(field, domain);
    multiply(scalar, field, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > multiply(T scalar, MultiScalarField3D<T> &field)
{
    return multiply(scalar, field, field.getBoundingBox());
}

template <typename T>
void divide(MultiScalarField3D<T> &field, T scalar, MultiScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(new A_dividedBy_alpha_functional3D<T>(scalar), domain, field, result);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > divide(MultiScalarField3D<T> &field, T scalar, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > result = generateMultiScalarField<T>(field, domain);
    divide(field, scalar, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > divide(MultiScalarField3D<T> &field, T scalar)
{
    return divide(field, scalar, field.getBoundingBox());
}

template <typename T>
void divide(T scalar, MultiScalarField3D<T> &field, MultiScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(new Alpha_dividedBy_A_functional3D<T>(scalar), domain, field, result);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > divide(T scalar, MultiScalarField3D<T> &field, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > result = generateMultiScalarField<T>(field, domain);
    divide(scalar, field, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > divide(T scalar, MultiScalarField3D<T> &field)
{
    return divide(scalar, field, field.getBoundingBox());
}

/* *************** MultiScalarField - Scalar inplace operations *************** */

template <typename T>
void addInPlace(MultiScalarField3D<T> &field, T scalar, Box3D domain)
{
    applyProcessingFunctional(new A_plus_alpha_inplace_functional3D<T>(scalar), domain, field);
}

template <typename T>
void addInPlace(MultiScalarField3D<T> &field, T scalar)
{
    addInPlace(field, scalar, field.getBoundingBox());
}

template <typename T>
void subtractInPlace(MultiScalarField3D<T> &field, T scalar, Box3D domain)
{
    applyProcessingFunctional(new A_minus_alpha_inplace_functional3D<T>(scalar), domain, field);
}

template <typename T>
void subtractInPlace(MultiScalarField3D<T> &field, T scalar)
{
    subtractInPlace(field, scalar, field.getBoundingBox());
}

template <typename T>
void multiplyInPlace(MultiScalarField3D<T> &field, T scalar, Box3D domain)
{
    applyProcessingFunctional(new A_times_alpha_inplace_functional3D<T>(scalar), domain, field);
}

template <typename T>
void multiplyInPlace(MultiScalarField3D<T> &field, T scalar)
{
    multiplyInPlace(field, scalar, field.getBoundingBox());
}

template <typename T>
void divideInPlace(MultiScalarField3D<T> &field, T scalar, Box3D domain)
{
    applyProcessingFunctional(new A_dividedBy_alpha_inplace_functional3D<T>(scalar), domain, field);
}

template <typename T>
void divideInPlace(MultiScalarField3D<T> &field, T scalar)
{
    divideInPlace(field, scalar, field.getBoundingBox());
}

/* *************** MultiScalarField - MultiScalarField operations *************** */

template <typename T1, typename T2>
void copy(MultiScalarField3D<T1> &field, MultiScalarField3D<T2> &convertedField, Box3D domain)
{
    applyProcessingFunctional(
        new CopyConvertScalarFunctional3D<T1, T2>, domain, field, convertedField);
}

template <typename T1, typename T2>
std::unique_ptr<MultiScalarField3D<T2> > copyConvert(MultiScalarField3D<T1> &field, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T2> > convertedField =
        generateMultiScalarField<T2>(field, domain);
    plb::copy(field, *convertedField, domain);
    return convertedField;
}

template <typename T1, typename T2>
std::unique_ptr<MultiScalarField3D<T2> > copyConvert(MultiScalarField3D<T1> &field)
{
    return copyConvert<T1, T2>(field, field.getBoundingBox());
}

template <typename T>
void lessThan(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, MultiScalarField3D<int> &result,
    Box3D domain)
{
    std::vector<MultiBlock3D *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_lt_B_functional3D<T>, domain, fields);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<int> > lessThan(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<int> > result =
        generateIntersectMultiScalarField<int>(A, B, domain);
    lessThan(A, B, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<int> > lessThan(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B)
{
    return lessThan(A, B, A.getBoundingBox());
}

template <typename T>
void greaterThan(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, MultiScalarField3D<int> &result,
    Box3D domain)
{
    std::vector<MultiBlock3D *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_gt_B_functional3D<T>, domain, fields);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<int> > greaterThan(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<int> > result =
        generateIntersectMultiScalarField<int>(A, B, domain);
    greaterThan(A, B, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<int> > greaterThan(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B)
{
    return greaterThan(A, B, A.getBoundingBox());
}

template <typename T>
void add(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, MultiScalarField3D<T> &result, Box3D domain)
{
    std::vector<MultiScalarField3D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_plus_B_functional3D<T>, domain, fields);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > add(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > result =
        generateIntersectMultiScalarField<T>(A, B, domain);
    add(A, B, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > add(MultiScalarField3D<T> &A, MultiScalarField3D<T> &B)
{
    return add(A, B, A.getBoundingBox());
}

template <typename T>
void subtract(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, MultiScalarField3D<T> &result, Box3D domain)
{
    std::vector<MultiScalarField3D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_minus_B_functional3D<T>, domain, fields);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > subtract(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > result =
        generateIntersectMultiScalarField<T>(A, B, domain);
    subtract(A, B, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > subtract(MultiScalarField3D<T> &A, MultiScalarField3D<T> &B)
{
    return subtract(A, B, A.getBoundingBox());
}

template <typename T>
void multiply(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, MultiScalarField3D<T> &result, Box3D domain)
{
    std::vector<MultiScalarField3D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_times_B_functional3D<T>, domain, fields);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > multiply(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > result =
        generateIntersectMultiScalarField<T>(A, B, domain);
    multiply(A, B, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > multiply(MultiScalarField3D<T> &A, MultiScalarField3D<T> &B)
{
    return multiply(A, B, A.getBoundingBox());
}

template <typename T>
void divide(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, MultiScalarField3D<T> &result, Box3D domain)
{
    std::vector<MultiScalarField3D<T> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new A_dividedBy_B_functional3D<T>, domain, fields);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > divide(
    MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > result =
        generateIntersectMultiScalarField<T>(A, B, domain);
    divide(A, B, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > divide(MultiScalarField3D<T> &A, MultiScalarField3D<T> &B)
{
    return divide(A, B, A.getBoundingBox());
}

/* *************** MultiScalarField - Scalar inplace operations *************** */

template <typename T>
void addInPlace(MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, Box3D domain)
{
    applyProcessingFunctional(new A_plus_B_inplace_functional3D<T>, domain, A, B);
}

template <typename T>
void addInPlace(MultiScalarField3D<T> &A, MultiScalarField3D<T> &B)
{
    addInPlace(A, B, A.getBoundingBox());
}

template <typename T>
void subtractInPlace(MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, Box3D domain)
{
    applyProcessingFunctional(new A_minus_B_inplace_functional3D<T>, domain, A, B);
}

template <typename T>
void subtractInPlace(MultiScalarField3D<T> &A, MultiScalarField3D<T> &B)
{
    subtractInPlace(A, B, A.getBoundingBox());
}

template <typename T>
void multiplyInPlace(MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, Box3D domain)
{
    applyProcessingFunctional(new A_times_B_inplace_functional3D<T>, domain, A, B);
}

template <typename T>
void multiplyInPlace(MultiScalarField3D<T> &A, MultiScalarField3D<T> &B)
{
    multiplyInPlace(A, B, A.getBoundingBox());
}

template <typename T>
void divideInPlace(MultiScalarField3D<T> &A, MultiScalarField3D<T> &B, Box3D domain)
{
    applyProcessingFunctional(new A_dividedBy_B_inplace_functional3D<T>, domain, A, B);
}

template <typename T>
void divideInPlace(MultiScalarField3D<T> &A, MultiScalarField3D<T> &B)
{
    divideInPlace(A, B, A.getBoundingBox());
}

template <typename T>
void uniformlyBoundScalarField(MultiScalarField3D<T> &data, Box3D domain, T bound)
{
    applyProcessingFunctional(new UniformlyBoundScalarField3D<T>(bound), domain, data);
}

template <typename T>
void uniformlyBoundScalarField(MultiScalarField3D<T> &data, T bound)
{
    uniformlyBoundScalarField(data, data.getBoundingBox(), bound);
}

template <typename T>
void boundScalarField(MultiScalarField3D<T> &data, Box3D domain, T lowerBound, T upperBound)
{
    applyProcessingFunctional(new BoundScalarField3D<T>(lowerBound, upperBound), domain, data);
}

template <typename T>
void boundScalarField(MultiScalarField3D<T> &data, T lowerBound, T upperBound)
{
    boundScalarField(data, data.getBoundingBox(), lowerBound, upperBound);
}

/* *************** MultiScalarField operations *************** */

template <typename T>
void computeSqrt(MultiScalarField3D<T> &A, MultiScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(new ComputeScalarSqrtFunctional3D<T>, domain, A, result);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeSqrt(MultiScalarField3D<T> &A, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > result = generateMultiScalarField<T>(A, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the ComputeScalarFunctional3D()
    // acts on both bulk and envelope, you would expect the envelope layer around the domain, on the
    // result multi-block, to be assigned some proper values too. By default, this is however not
    // what happens, because the physical space occupied by these envelopes does not intersect with
    // the domain "domain". We work around this issue by extending the domain. There's no problem if
    // the enlarged domain gets beyond the actual extent of the lattice, because Palabos handles
    // these situations properly.

    computeSqrt(A, *result, domain.enlarge(A.getMultiBlockManagement().getEnvelopeWidth()));
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeSqrt(MultiScalarField3D<T> &A)
{
    return computeSqrt(A, A.getBoundingBox());
}

template <typename T>
void computeAbsoluteValue(MultiScalarField3D<T> &A, MultiScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(new ComputeAbsoluteValueFunctional3D<T>, domain, A, result);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeAbsoluteValue(MultiScalarField3D<T> &A, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > result = generateMultiScalarField<T>(A, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // ComputeAbsoluteValueFunctional3D() acts on both bulk and envelope, you would expect the
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
std::unique_ptr<MultiScalarField3D<T> > computeAbsoluteValue(MultiScalarField3D<T> &A)
{
    return computeAbsoluteValue(A, A.getBoundingBox());
}

/* *************** LBMsmoothen3D ******************************************* */

template <typename T, template <typename U> class Descriptor>
void lbmSmoothen(MultiScalarField3D<T> &data, MultiScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(new LBMsmoothen3D<T, Descriptor>, domain, data, result);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > lbmSmoothen(MultiScalarField3D<T> &data, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > result = generateMultiScalarField<T>(data, domain);
    lbmSmoothen<T, Descriptor>(data, *result, domain);
    return result;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > lbmSmoothen(MultiScalarField3D<T> &data)
{
    return lbmSmoothen<T, Descriptor>(data, data.getBoundingBox());
}

/* *************** LBMsmoothenInPlace3D ******************************************* */

template <typename T, template <typename U> class Descriptor>
void lbmSmoothenInPlace(MultiScalarField3D<T> &data, Box3D domain)
{
    applyProcessingFunctional(new LBMsmoothenInPlace3D<T, Descriptor>, domain, data);
}

template <typename T, template <typename U> class Descriptor>
void lbmSmoothenInPlace(MultiScalarField3D<T> &data)
{
    lbmSmoothenInPlace<T, Descriptor>(data, data.getBoundingBox());
}

/* *************** Smoothen3D ******************************************* */

template <typename T>
void smoothen(MultiScalarField3D<T> &data, MultiScalarField3D<T> &result, Box3D domain)
{
    applyProcessingFunctional(new Smoothen3D<T>, domain, data, result);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > smoothen(MultiScalarField3D<T> &data, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > result = generateMultiScalarField<T>(data, domain);
    smoothen<T>(data, *result, domain);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > smoothen(MultiScalarField3D<T> &data)
{
    return smoothen<T>(data, data.getBoundingBox());
}

/* *************** MollifyScalar3D ******************************************* */

template <typename T>
void mollifyScalar(
    T l, plint d, Box3D globalDomain, int exclusionFlag, MultiScalarField3D<T> &data,
    MultiScalarField3D<int> &flag, MultiScalarField3D<T> &result, Box3D domain)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&data);
    args.push_back(&flag);
    args.push_back(&result);
    applyProcessingFunctional(
        new MollifyScalar3D<T>(l, d, globalDomain, exclusionFlag), domain, args);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > mollifyScalar(
    T l, plint d, Box3D globalDomain, int exclusionFlag, MultiScalarField3D<T> &data,
    MultiScalarField3D<int> &flag, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > result = generateMultiScalarField<T>(data, domain);

    std::vector<MultiBlock3D *> args;
    args.push_back(&data);
    args.push_back(&flag);
    args.push_back(result.get());
    applyProcessingFunctional(
        new MollifyScalar3D<T>(l, d, globalDomain, exclusionFlag), domain, args);

    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > mollifyScalar(
    T l, plint d, Box3D globalDomain, int exclusionFlag, MultiScalarField3D<T> &data,
    MultiScalarField3D<int> &flag)
{
    std::unique_ptr<MultiScalarField3D<T> > result =
        generateMultiScalarField<T>(data, data.getBoundingBox());

    std::vector<MultiBlock3D *> args;
    args.push_back(&data);
    args.push_back(&flag);
    args.push_back(result.get());
    applyProcessingFunctional(
        new MollifyScalar3D<T>(l, d, globalDomain, exclusionFlag), data.getBoundingBox(), args);

    return result;
}

/* *************** LBMcomputeGradient3D ******************************************* */

template <typename T, template <typename U> class Descriptor>
void lbmComputeGradient(
    MultiScalarField3D<T> &scalarField, MultiTensorField3D<T, 3> &gradient, Box3D domain)
{
    applyProcessingFunctional(
        new LBMcomputeGradient3D<T, Descriptor>, domain, scalarField, gradient);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, 3> > lbmComputeGradient(
    MultiScalarField3D<T> &scalarField, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, 3> > gradient =
        generateMultiTensorField<T, 3>(scalarField, domain);
    lbmComputeGradient<T, Descriptor>(scalarField, *gradient, domain);
    return gradient;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, 3> > lbmComputeGradient(MultiScalarField3D<T> &scalarField)
{
    return lbmComputeGradient<T, Descriptor>(scalarField, scalarField.getBoundingBox());
}

/* ******************************************************************* */
/* *************** PART VI. Multi-block wrappers: Tensor-field ******* */
/* ******************************************************************* */

/* *************** Reductive functions ******************************* */

template <typename T, plint nDim, class BoolMask>
plint count(MultiTensorField3D<T, nDim> &field, Box3D domain, BoolMask boolMask)
{
    CountTensorElementsFunctional3D<T, nDim, BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, field);
    return functional.getCount();
}

template <typename T, plint nDim, class BoolMask>
plint count(MultiTensorField3D<T, nDim> &field, BoolMask boolMask)
{
    return count(field, field.getBoundingBox(), boolMask);
}

template <typename T, int nDim>
Array<T, nDim> computeSum(MultiTensorField3D<T, nDim> &tensorField, Box3D domain)
{
    BoxTensorSumFunctional3D<T, nDim> functional;
    applyProcessingFunctional(functional, domain, tensorField);
    return functional.getSumTensor();
}

template <typename T, int nDim>
Array<T, nDim> computeSum(MultiTensorField3D<T, nDim> &tensorField)
{
    return computeSum(tensorField, tensorField.getBoundingBox());
}

template <typename T, int nDim>
Array<T, nDim> computeSum(
    MultiTensorField3D<T, nDim> &tensorField, MultiScalarField3D<int> &mask, int flag, Box3D domain)
{
    MaskedBoxTensorSumFunctional3D<T, nDim> functional(flag);
    applyProcessingFunctional(functional, domain, mask, tensorField);
    return functional.getSumTensor();
}

template <typename T, int nDim>
Array<T, nDim> computeSum(
    MultiTensorField3D<T, nDim> &tensorField, MultiScalarField3D<int> &mask, int flag)
{
    return computeSum(tensorField, mask, flag, tensorField.getBoundingBox());
}

template <typename T, int nDim>
Array<T, nDim> computeAverage(MultiTensorField3D<T, nDim> &tensorField, Box3D domain)
{
    BoxTensorSumFunctional3D<T, nDim> functional;
    applyProcessingFunctional(functional, domain, tensorField);
    return functional.getSumTensor() / (T)domain.nCells();
}

template <typename T, int nDim>
Array<T, nDim> computeAverage(MultiTensorField3D<T, nDim> &tensorField)
{
    return computeAverage(tensorField, tensorField.getBoundingBox());
}

template <typename T, int nDim>
Array<T, nDim> computeAverage(
    MultiTensorField3D<T, nDim> &tensorField, MultiScalarField3D<int> &mask, int flag, Box3D domain)
{
    MaskedBoxTensorAverageFunctional3D<T, nDim> functional(flag);
    applyProcessingFunctional(functional, domain, mask, tensorField);
    return functional.getAverageTensor();
}

template <typename T, int nDim>
Array<T, nDim> computeAverage(
    MultiTensorField3D<T, nDim> &tensorField, MultiScalarField3D<int> &mask, int flag)
{
    return computeAverage(tensorField, mask, flag, tensorField.getBoundingBox());
}

/* *************** Extract Sub-TensorField *************************** */

template <typename T, int nDim>
void extractSubDomain(
    MultiTensorField3D<T, nDim> &field, MultiTensorField3D<T, nDim> &extractedField, Box3D domain)
{
    applyProcessingFunctional(
        new ExtractTensorSubDomainFunctional3D<T, nDim>, domain, field, extractedField);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > extractSubDomain(
    MultiTensorField3D<T, nDim> &field, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, nDim> > extractedField =
        generateMultiTensorField<T, nDim>(field, domain);
    extractSubDomain(field, *extractedField, domain);
    return extractedField;
}

/* *************** Component (scalar-field) out of a tensor-field ****** */

template <typename T, int nDim>
void extractComponent(
    MultiTensorField3D<T, nDim> &tensorField, MultiScalarField3D<T> &component, Box3D domain,
    int iComponent)
{
    applyProcessingFunctional(
        new ExtractTensorComponentFunctional3D<T, nDim>(iComponent), domain, component,
        tensorField);
}

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T> > extractComponent(
    MultiTensorField3D<T, nDim> &tensorField, Box3D domain, int iComponent)
{
    std::unique_ptr<MultiScalarField3D<T> > component =
        generateMultiScalarField<T>(tensorField, domain);
    extractComponent(tensorField, *component, domain, iComponent);
    return component;
}

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T> > extractComponent(
    MultiTensorField3D<T, nDim> &tensorField, int iComponent)
{
    return extractComponent(tensorField, tensorField.getBoundingBox(), iComponent);
}

/* *************** Vector-norm of each cell in the field *************** */

template <typename T, int nDim>
void computeNorm(
    MultiTensorField3D<T, nDim> &tensorField, MultiScalarField3D<T> &norm, Box3D domain)
{
    applyProcessingFunctional(new ComputeNormFunctional3D<T, nDim>, domain, norm, tensorField);
}

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T> > computeNorm(
    MultiTensorField3D<T, nDim> &tensorField, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > norm = generateMultiScalarField<T>(tensorField, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the ComputeNormFunctional3D()
    // acts on both bulk and envelope, you would expect the envelope layer around the domain, on the
    // norm multi-block, to be assigned some proper values too. By default, this is however not what
    // happens, because the physical space occupied by these envelopes does not intersect with the
    // domain "domain". We work around this issue by extending the domain. There's no problem if the
    // enlarged domain gets beyond the actual extent of the lattice, because Palabos handles these
    // situations properly.

    computeNorm(
        tensorField, *norm,
        domain.enlarge(tensorField.getMultiBlockManagement().getEnvelopeWidth()));
    return norm;
}

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T> > computeNorm(MultiTensorField3D<T, nDim> &tensorField)
{
    return computeNorm(tensorField, tensorField.getBoundingBox());
}

/* *************** Sqrt operation on each component of each cell *************** */

template <typename T, int nDim>
void computeSqrt(MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &result, Box3D domain)
{
    applyProcessingFunctional(new ComputeTensorSqrtFunctional3D<T, nDim>, domain, A, result);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > computeSqrt(
    MultiTensorField3D<T, nDim> &A, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, nDim> > result =
        generateMultiTensorField<T, nDim>(A, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // ComputeTensorSqrtFunctional3D() acts on both bulk and envelope, you would expect the envelope
    // layer around the domain, on the result multi-block, to be assigned some proper values too. By
    // default, this is however not what happens, because the physical space occupied by these
    // envelopes does not intersect with the domain "domain". We work around this issue by extending
    // the domain. There's no problem if the enlarged domain gets beyond the actual extent of the
    // lattice, because Palabos handles these situations properly.

    computeSqrt(A, *result, domain.enlarge(A.getMultiBlockManagement().getEnvelopeWidth()));
    return result;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > computeSqrt(MultiTensorField3D<T, nDim> &A)
{
    return computeSqrt(A, A.getBoundingBox());
}

/* *************** Squared vector-norm of each cell in the field ******** */

template <typename T, int nDim>
void computeNormSqr(
    MultiTensorField3D<T, nDim> &tensorField, MultiScalarField3D<T> &normSqr, Box3D domain)
{
    applyProcessingFunctional(
        new ComputeNormSqrFunctional3D<T, nDim>, domain, normSqr, tensorField);
}

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T> > computeNormSqr(
    MultiTensorField3D<T, nDim> &tensorField, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > normSqr =
        generateMultiScalarField<T>(tensorField, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the ComputeNormSqrFunctional3D()
    // acts on both bulk and envelope, you would expect the envelope layer around the domain, on the
    // normSqr multi-block, to be assigned some proper values too. By default, this is however not
    // what happens, because the physical space occupied by these envelopes does not intersect with
    // the domain "domain". We work around this issue by extending the domain. There's no problem if
    // the enlarged domain gets beyond the actual extent of the lattice, because Palabos handles
    // these situations properly.

    computeNormSqr(
        tensorField, *normSqr,
        domain.enlarge(tensorField.getMultiBlockManagement().getEnvelopeWidth()));
    return normSqr;
}

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T> > computeNormSqr(MultiTensorField3D<T, nDim> &tensorField)
{
    return computeNormSqr(tensorField, tensorField.getBoundingBox());
}

/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template <typename T>
void computeSymmetricTensorNorm(
    MultiTensorField3D<T, 6> &tensorField, MultiScalarField3D<T> &norm, Box3D domain)
{
    applyProcessingFunctional(
        new ComputeSymmetricTensorNormFunctional3D<T>, domain, norm, tensorField);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeSymmetricTensorNorm(
    MultiTensorField3D<T, 6> &tensorField, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > norm = generateMultiScalarField<T>(tensorField, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // ComputeSymmetricTensorNormFunctional3D() acts on both bulk and envelope, you would expect the
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
std::unique_ptr<MultiScalarField3D<T> > computeSymmetricTensorNorm(
    MultiTensorField3D<T, 6> &tensorField)
{
    return computeSymmetricTensorNorm(tensorField, tensorField.getBoundingBox());
}

/* *************** Squared Tensor-norm of each symmetric tensor of a field ******************** */

template <typename T>
void computeSymmetricTensorNormSqr(
    MultiTensorField3D<T, 6> &tensorField, MultiScalarField3D<T> &normSqr, Box3D domain)
{
    applyProcessingFunctional(
        new ComputeSymmetricTensorNormSqrFunctional3D<T>, domain, normSqr, tensorField);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeSymmetricTensorNormSqr(
    MultiTensorField3D<T, 6> &tensorField, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > normSqr =
        generateMultiScalarField<T>(tensorField, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // ComputeSymmetricTensorNormSqrFunctional3D() acts on both bulk and envelope, you would expect
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
std::unique_ptr<MultiScalarField3D<T> > computeSymmetricTensorNormSqr(
    MultiTensorField3D<T, 6> &tensorField)
{
    return computeSymmetricTensorNormSqr(tensorField, tensorField.getBoundingBox());
}

/* *************** Trace of each symmetric tensor of a field ************ */

template <typename T>
void computeSymmetricTensorTrace(
    MultiTensorField3D<T, 6> &tensorField, MultiScalarField3D<T> &trace, Box3D domain)
{
    applyProcessingFunctional(
        new ComputeSymmetricTensorTraceFunctional3D<T>, domain, trace, tensorField);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeSymmetricTensorTrace(
    MultiTensorField3D<T, 6> &tensorField, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > trace =
        generateMultiScalarField<T>(tensorField, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // ComputeSymmetricTensorTraceFunctional3D() acts on both bulk and envelope, you would expect
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
std::unique_ptr<MultiScalarField3D<T> > computeSymmetricTensorTrace(
    MultiTensorField3D<T, 6> &tensorField)
{
    return computeSymmetricTensorTrace(tensorField, tensorField.getBoundingBox());
}

/* *************** Gradient from scalar field *********************** */

template <typename T>
void computeGradient(MultiScalarField3D<T> &phi, MultiTensorField3D<T, 3> &gradient, Box3D domain)
{
    plint envelopeWidth = 1;
    applyProcessingFunctional(new BoxGradientFunctional3D<T>, domain, phi, gradient, envelopeWidth);
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeGradient(MultiScalarField3D<T> &phi, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, 3> > gradient =
        generateMultiTensorField<T, 3>(phi, domain);
    computeGradient(phi, *gradient, domain);
    return gradient;
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeGradient(MultiScalarField3D<T> &phi)
{
    return computeGradient(phi, phi.getBoundingBox());
}

/* *************** Gradient, without boundary treatment, from scalar field ***************** */

template <typename T>
void computeBulkGradient(
    MultiScalarField3D<T> &phi, MultiTensorField3D<T, 3> &gradient, Box3D domain)
{
    applyProcessingFunctional(new BoxBulkGradientFunctional3D<T>, domain, phi, gradient);
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeBulkGradient(
    MultiScalarField3D<T> &phi, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, 3> > gradient =
        generateMultiTensorField<T, 3>(phi, domain);
    computeBulkGradient(phi, *gradient, domain);
    return gradient;
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeBulkGradient(MultiScalarField3D<T> &phi)
{
    return computeBulkGradient(phi, phi.getBoundingBox());
}

/* *************** Vorticity from Velocity field *********************** */

template <typename T>
void computeVorticity(
    MultiTensorField3D<T, 3> &velocity, MultiTensorField3D<T, 3> &vorticity, Box3D domain)
{
    plint envelopeWidth = 1;
    applyProcessingFunctional(
        new BoxVorticityFunctional3D<T, 3>, domain, velocity, vorticity, envelopeWidth);
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeVorticity(
    MultiTensorField3D<T, 3> &velocity, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, 3> > vorticity =
        generateMultiTensorField<T, 3>(velocity, domain);
    computeVorticity(velocity, *vorticity, domain);
    return vorticity;
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeVorticity(MultiTensorField3D<T, 3> &velocity)
{
    return computeVorticity(velocity, velocity.getBoundingBox());
}

/* *************** Vorticity, without boundary treatment, from Velocity field ***************** */

template <typename T>
void computeBulkVorticity(
    MultiTensorField3D<T, 3> &velocity, MultiTensorField3D<T, 3> &vorticity, Box3D domain)
{
    applyProcessingFunctional(new BoxBulkVorticityFunctional3D<T, 3>, domain, velocity, vorticity);
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeBulkVorticity(
    MultiTensorField3D<T, 3> &velocity, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, 3> > vorticity =
        generateMultiTensorField<T, 3>(velocity, domain);
    computeBulkVorticity(velocity, *vorticity, domain);
    return vorticity;
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeBulkVorticity(MultiTensorField3D<T, 3> &velocity)
{
    return computeBulkVorticity(velocity, velocity.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeBulkVorticity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > vel =
        generateMultiTensorField<T, Descriptor<T>::d>(lattice, domain);
    vel->periodicity().toggle(0, lattice.periodicity().get(0));
    vel->periodicity().toggle(1, lattice.periodicity().get(1));
    vel->periodicity().toggle(2, lattice.periodicity().get(2));
    std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > vorticity =
        generateMultiTensorField<T, Descriptor<T>::d>(*vel, 1);

    computeVelocity(lattice, *vel, domain.enlarge(1));
    computeBulkVorticity(*vel, *vorticity, domain);
    return vorticity;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeBulkVorticity(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeBulkVorticity<T, Descriptor>(lattice, lattice.getBoundingBox());
}

/* *************** Order 4 Vorticity, without boundary treatment, from Velocity field
 * ***************** */

template <typename T>
void computeBulkVorticityOrderFour(
    MultiTensorField3D<T, 3> &velocity, MultiTensorField3D<T, 3> &vorticity, Box3D domain)
{
    applyProcessingFunctional(
        new BoxBulkVorticityOrderFourFunctional3D<T, 3>, domain, velocity, vorticity);
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeBulkVorticityOrderFour(
    MultiTensorField3D<T, 3> &velocity, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, 3> > vorticity =
        generateMultiTensorField<T, 3>(velocity, domain);
    computeBulkVorticityOrderFour(velocity, *vorticity, domain);
    return vorticity;
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeBulkVorticityOrderFour(
    MultiTensorField3D<T, 3> &velocity)
{
    return computeBulkVorticityOrderFour(velocity, velocity.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeBulkVorticityOrderFour(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > vel =
        generateMultiTensorField<T, Descriptor<T>::d>(domain, 3);

    vel->periodicity().toggle(0, lattice.periodicity().get(0));
    vel->periodicity().toggle(1, lattice.periodicity().get(1));
    vel->periodicity().toggle(2, lattice.periodicity().get(2));

    computeVelocity(lattice, *vel, domain);
    return computeBulkVorticityOrderFour(*vel, domain);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeBulkVorticityOrderFour(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeBulkVorticityOrderFour<T, Descriptor>(lattice, lattice.getBoundingBox());
}

/* *************** Order 6 Vorticity, without boundary treatment, from Velocity field
 * ***************** */

template <typename T>
void computeBulkVorticityOrderSix(
    MultiTensorField3D<T, 3> &velocity, MultiTensorField3D<T, 3> &vorticity, Box3D domain)
{
    applyProcessingFunctional(
        new BoxBulkVorticityOrderSixFunctional3D<T, 3>, domain, velocity, vorticity);
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeBulkVorticityOrderSix(
    MultiTensorField3D<T, 3> &velocity, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, 3> > vorticity =
        generateMultiTensorField<T, 3>(velocity, domain);
    computeBulkVorticityOrderSix(velocity, *vorticity, domain);
    return vorticity;
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeBulkVorticityOrderSix(
    MultiTensorField3D<T, 3> &velocity)
{
    return computeBulkVorticityOrderSix(velocity, velocity.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeBulkVorticityOrderSix(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > vel =
        generateMultiTensorField<T, Descriptor<T>::d>(domain, 3);

    vel->periodicity().toggle(0, lattice.periodicity().get(0));
    vel->periodicity().toggle(1, lattice.periodicity().get(1));
    vel->periodicity().toggle(2, lattice.periodicity().get(2));

    computeVelocity(lattice, *vel, domain);
    return computeBulkVorticityOrderSix(*vel, domain);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeBulkVorticityOrderSix(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeBulkVorticityOrderSix<T, Descriptor>(lattice, lattice.getBoundingBox());
}

/* *************** Order 8 Vorticity, without boundary treatment, from Velocity field
 * ***************** */

template <typename T>
void computeBulkVorticityOrderEight(
    MultiTensorField3D<T, 3> &velocity, MultiTensorField3D<T, 3> &vorticity, Box3D domain)
{
    applyProcessingFunctional(
        new BoxBulkVorticityOrderEightFunctional3D<T, 3>, domain, velocity, vorticity);
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeBulkVorticityOrderEight(
    MultiTensorField3D<T, 3> &velocity, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, 3> > vorticity =
        generateMultiTensorField<T, 3>(velocity, domain);
    computeBulkVorticityOrderEight(velocity, *vorticity, domain);
    return vorticity;
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeBulkVorticityOrderEight(
    MultiTensorField3D<T, 3> &velocity)
{
    return computeBulkVorticityOrderEight(velocity, velocity.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeBulkVorticityOrderEight(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > vel =
        generateMultiTensorField<T, Descriptor<T>::d>(domain, 4);

    vel->periodicity().toggle(0, lattice.periodicity().get(0));
    vel->periodicity().toggle(1, lattice.periodicity().get(1));
    vel->periodicity().toggle(2, lattice.periodicity().get(2));

    computeVelocity(lattice, *vel, domain);
    return computeBulkVorticityOrderEight(*vel, domain);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeBulkVorticityOrderEight(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeBulkVorticityOrderEight<T, Descriptor>(lattice, lattice.getBoundingBox());
}

/* *************** Helicity from Velocity field *********************** */

template <typename T>
void computeHelicity(
    MultiTensorField3D<T, 3> &velocity, MultiScalarField3D<T> &helicity, Box3D domain)
{
    plint envelopeWidth = 1;
    applyProcessingFunctional(
        new BoxHelicityFunctional3D<T, 3>, domain, helicity, velocity, envelopeWidth);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeHelicity(
    MultiTensorField3D<T, 3> &velocity, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > helicity =
        generateMultiScalarField<T>(velocity, domain);
    computeHelicity(velocity, *helicity, domain);
    return helicity;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeHelicity(MultiTensorField3D<T, 3> &velocity)
{
    return computeHelicity(velocity, velocity.getBoundingBox());
}

/* *************** Helicity, without boundary treatment, from Velocity field ***************** */

template <typename T>
void computeBulkHelicity(
    MultiTensorField3D<T, 3> &velocity, MultiScalarField3D<T> &helicity, Box3D domain)
{
    applyProcessingFunctional(new BoxBulkHelicityFunctional3D<T, 3>, domain, helicity, velocity);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeBulkHelicity(
    MultiTensorField3D<T, 3> &velocity, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > helicity =
        generateMultiScalarField<T>(velocity, domain);
    computeBulkHelicity(velocity, *helicity, domain);
    return helicity;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeBulkHelicity(MultiTensorField3D<T, 3> &velocity)
{
    return computeBulkHelicity(velocity, velocity.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeBulkHelicity(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > vel =
        generateMultiTensorField<T, Descriptor<T>::d>(lattice, domain);
    vel->periodicity().toggle(0, lattice.periodicity().get(0));
    vel->periodicity().toggle(1, lattice.periodicity().get(1));
    vel->periodicity().toggle(2, lattice.periodicity().get(2));
    std::unique_ptr<MultiScalarField3D<T> > helicity = generateMultiScalarField<T>(*vel, 1);

    computeVelocity(lattice, *vel, domain.enlarge(1));
    computeBulkHelicity(*vel, *helicity, domain);
    return helicity;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeBulkHelicity(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return computeBulkHelicity<T, Descriptor>(lattice, lattice.getBoundingBox());
}

/* *************** Divergence, without boundary treatment, from Velocity field **************** */

template <typename T>
void computeBulkDivergence(
    MultiTensorField3D<T, 3> &velocity, MultiScalarField3D<T> &divergence, Box3D domain)
{
    applyProcessingFunctional(
        new BoxBulkDivergenceFunctional3D<T, 3>, domain, divergence, velocity);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeBulkDivergence(
    MultiTensorField3D<T, 3> &velocity, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > divergence =
        generateMultiScalarField<T>(velocity, domain);
    computeBulkDivergence(velocity, *divergence, domain);
    return divergence;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeBulkDivergence(MultiTensorField3D<T, 3> &velocity)
{
    return computeBulkDivergence(velocity, velocity.getBoundingBox());
}

/* *************** Strain rate from Velocity field ********************* */

template <typename T>
void computeStrainRate(
    MultiTensorField3D<T, 3> &velocity, MultiTensorField3D<T, 6> &S, Box3D domain)
{
    plint envelopeWidth = 1;
    applyProcessingFunctional(
        new BoxStrainRateFunctional3D<T, 3>, domain, velocity, S, envelopeWidth);
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 6> > computeStrainRate(
    MultiTensorField3D<T, 3> &velocity, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, 6> > S = generateMultiTensorField<T, 6>(velocity, domain);
    computeStrainRate(velocity, *S, domain);
    return S;
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 6> > computeStrainRate(MultiTensorField3D<T, 3> &velocity)
{
    return computeStrainRate(velocity, velocity.getBoundingBox());
}

/* *************** Str. rate, without boundary treatment, from Velocity field ****************** */

template <typename T>
void computeBulkStrainRate(
    MultiTensorField3D<T, 3> &velocity, MultiTensorField3D<T, 6> &S, Box3D domain)
{
    applyProcessingFunctional(new BoxBulkStrainRateFunctional3D<T, 3>, domain, velocity, S);
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 6> > computeBulkStrainRate(
    MultiTensorField3D<T, 3> &velocity, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, 6> > S = generateMultiTensorField<T, 6>(velocity, domain);
    computeBulkStrainRate(velocity, *S, domain);
    return S;
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 6> > computeBulkStrainRate(MultiTensorField3D<T, 3> &velocity)
{
    return computeBulkStrainRate(velocity, velocity.getBoundingBox());
}

/* *************** Computes instantaneous Reynolds stresses from average velocity and velocity
 * ******************** */

template <typename T>
void computeInstantaneousReynoldsStress(
    MultiTensorField3D<T, 3> &vel, MultiTensorField3D<T, 3> &avgVel, MultiTensorField3D<T, 6> &tau,
    Box3D domain)
{
    std::vector<MultiBlock3D *> fields;
    fields.push_back(&vel);
    fields.push_back(&avgVel);
    fields.push_back(&tau);

    applyProcessingFunctional(
        new BoxComputeInstantaneousReynoldsStressFunctional3D<T>, domain, fields);
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 6> > computeInstantaneousReynoldsStress(
    MultiTensorField3D<T, 3> &vel, MultiTensorField3D<T, 3> &avgVel, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, 6> > tau = generateMultiTensorField<T, 6>(vel, domain);
    computeInstantaneousReynoldsStress(vel, avgVel, *tau, domain);
    return tau;
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 6> > computeInstantaneousReynoldsStress(
    MultiTensorField3D<T, 3> &vel, MultiTensorField3D<T, 3> &avgVel)
{
    return computeInstantaneousReynoldsStress(vel, avgVel, vel.getBoundingBox());
}

/* *************** Q-criterion from vorticity and strain rate fields ******************** */

template <typename T>
void computeQcriterion(
    MultiTensorField3D<T, 3> &vorticity, MultiTensorField3D<T, 6> &S,
    MultiScalarField3D<T> &qCriterion, Box3D domain)
{
    std::vector<MultiBlock3D *> fields;
    fields.push_back(&vorticity);
    fields.push_back(&S);
    fields.push_back(&qCriterion);

    applyProcessingFunctional(new BoxQcriterionFunctional3D<T>, domain, fields);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeQcriterion(
    MultiTensorField3D<T, 3> &vorticity, MultiTensorField3D<T, 6> &S, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > qCriterion =
        generateMultiScalarField<T>(vorticity, domain);
    computeQcriterion(vorticity, S, *qCriterion, domain);
    return qCriterion;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeQcriterion(
    MultiTensorField3D<T, 3> &vorticity, MultiTensorField3D<T, 6> &S)
{
    return computeQcriterion(vorticity, S, vorticity.getBoundingBox());
}

/* *************** lambda2-criterion from vorticity and strain rate fields ******************** */
#ifndef PLB_BGP
#ifdef PLB_USE_EIGEN
template <typename T>
void computeLambda2(
    MultiTensorField3D<T, 3> &vorticity, MultiTensorField3D<T, 6> &S,
    MultiScalarField3D<T> &lambda2, Box3D domain)
{
    std::vector<MultiBlock3D *> fields;
    fields.push_back(&vorticity);
    fields.push_back(&S);
    fields.push_back(&lambda2);

    applyProcessingFunctional(new BoxLambda2Functional3D<T>, domain, fields);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeLambda2(
    MultiTensorField3D<T, 3> &vorticity, MultiTensorField3D<T, 6> &S, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > lambda2 =
        generateMultiScalarField<T>(vorticity, domain);
    computeLambda2(vorticity, S, *lambda2, domain);
    return lambda2;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeLambda2(
    MultiTensorField3D<T, 3> &vorticity, MultiTensorField3D<T, 6> &S)
{
    return computeLambda2(vorticity, S, vorticity.getBoundingBox());
}
#endif
#endif

/* *************** MultiTensorField - MultiTensorField operations *************** */

template <typename T1, typename T2, int nDim>
void copy(
    MultiTensorField3D<T1, nDim> &field, MultiTensorField3D<T2, nDim> &convertedField, Box3D domain)
{
    applyProcessingFunctional(
        new CopyConvertTensorFunctional3D<T1, T2, nDim>, domain, field, convertedField);
}

template <typename T1, typename T2, int nDim>
std::unique_ptr<MultiTensorField3D<T2, nDim> > copyConvert(
    MultiTensorField3D<T1, nDim> &field, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T2, nDim> > convertedField =
        generateMultiTensorField<T2, nDim>(field, domain);
    plb::copy<T1, T2, nDim>(field, *convertedField, domain);
    return convertedField;
}

template <typename T1, typename T2, int nDim>
std::unique_ptr<MultiTensorField3D<T2, nDim> > copyConvert(MultiTensorField3D<T1, nDim> &field)
{
    return copyConvert<T1, T2, nDim>(field, field.getBoundingBox());
}

template <typename T, int nDim>
void add(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B,
    MultiTensorField3D<T, nDim> &result, Box3D domain)
{
    std::vector<MultiTensorField3D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new Tensor_A_plus_B_functional3D<T, nDim>, domain, fields);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > add(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, nDim> > result =
        generateIntersectMultiTensorField<T, nDim>(A, B, domain);
    add(A, B, *result, domain);
    return result;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > add(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B)
{
    return add(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void subtract(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B,
    MultiTensorField3D<T, nDim> &result, Box3D domain)
{
    std::vector<MultiTensorField3D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new Tensor_A_minus_B_functional3D<T, nDim>, domain, fields);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > subtract(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, nDim> > result =
        generateIntersectMultiTensorField<T, nDim>(A, B, domain);
    subtract(A, B, *result, domain);
    return result;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > subtract(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B)
{
    return subtract(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void multiply(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B,
    MultiTensorField3D<T, nDim> &result, Box3D domain)
{
    std::vector<MultiTensorField3D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new Tensor_A_times_B_functional3D<T, nDim>, domain, fields);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > multiply(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, nDim> > result =
        generateIntersectMultiTensorField<T, nDim>(A, B, domain);
    multiply(A, B, *result, domain);
    return result;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > multiply(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B)
{
    return multiply(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void multiply(
    MultiTensorField3D<T, nDim> &field, T scalar, MultiTensorField3D<T, nDim> &result, Box3D domain)
{
    applyProcessingFunctional(
        new Tensor_A_times_alpha_functional3D<T, nDim>(scalar), domain, field, result);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > multiply(
    MultiTensorField3D<T, nDim> &field, T scalar, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, nDim> > result =
        generateMultiTensorField<T, nDim>(field, domain);
    multiply(field, scalar, *result, domain);
    return result;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > multiply(MultiTensorField3D<T, nDim> &field, T scalar)
{
    return multiply(field, scalar, field.getBoundingBox());
}

template <typename T, int nDim>
void multiply(
    T scalar, MultiTensorField3D<T, nDim> &field, MultiTensorField3D<T, nDim> &result, Box3D domain)
{
    applyProcessingFunctional(
        new Tensor_A_times_alpha_functional3D<T, nDim>(scalar), domain, field, result);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > multiply(
    T scalar, MultiTensorField3D<T, nDim> &field, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, nDim> > result =
        generateMultiTensorField<T, nDim>(field, domain);
    multiply(scalar, field, *result, domain);
    return result;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > multiply(T scalar, MultiTensorField3D<T, nDim> &field)
{
    return multiply(scalar, field, field.getBoundingBox());
}

template <typename T, int nDim>
void divide(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B,
    MultiTensorField3D<T, nDim> &result, Box3D domain)
{
    std::vector<MultiTensorField3D<T, nDim> *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new Tensor_A_dividedBy_B_functional3D<T, nDim>, domain, fields);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > divide(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, nDim> > result =
        generateIntersectMultiTensorField<T, nDim>(A, B, domain);
    divide(A, B, *result, domain);
    return result;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > divide(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B)
{
    return divide(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void normalize(MultiTensorField3D<T, nDim> &data, MultiTensorField3D<T, nDim> &result, Box3D domain)
{
    applyProcessingFunctional(new Normalize_Tensor_functional3D<T, nDim>(), domain, data, result);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > normalize(
    MultiTensorField3D<T, nDim> &data, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, nDim> > result =
        generateMultiTensorField<T, nDim>(data, domain);
    normalize(data, *result, domain);
    return result;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > normalize(MultiTensorField3D<T, nDim> &data)
{
    return normalize(data, data.getBoundingBox());
}

// tensor product of a vector with itself
template <typename T, int nDim>
void symmetricTensorProduct(
    MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, SymmetricTensorImpl<T, nDim>::n> &result,
    Box3D domain)
{
    std::vector<MultiBlock3D *> fields;
    fields.push_back(&A);
    fields.push_back(&result);
    applyProcessingFunctional(new TensorProduct_A_A_functional3D<T, nDim>, domain, fields);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, SymmetricTensorImpl<T, nDim>::n> > symmetricTensorProduct(
    MultiTensorField3D<T, nDim> &A, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, SymmetricTensorImpl<T, nDim>::n> > result =
        generateMultiTensorField<T, SymmetricTensorImpl<T, nDim>::n>(A, domain);
    symmetricTensorProduct(A, *result, domain);
    return result;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, SymmetricTensorImpl<T, nDim>::n> > symmetricTensorProduct(
    MultiTensorField3D<T, nDim> &A)
{
    return symmetricTensorProduct(A, A.getBoundingBox());
}

/* *************** MultiScalarField - MultiTensorField operations *************** */

template <typename T, int nDim>
void multiply(
    MultiScalarField3D<T> &A, MultiTensorField3D<T, nDim> &B, MultiTensorField3D<T, nDim> &result,
    Box3D domain)
{
    std::vector<MultiBlock3D *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(new Scalar_A_times_Tensor_B_functional3D<T, nDim>, domain, fields);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > multiply(
    MultiScalarField3D<T> &A, MultiTensorField3D<T, nDim> &B, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, nDim> > result =
        generateMultiTensorField<T, nDim>(B, domain);
    multiply(A, B, *result, domain);
    return result;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > multiply(
    MultiScalarField3D<T> &A, MultiTensorField3D<T, nDim> &B)
{
    return multiply(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void fullIndexContractionOfSymmetricTensors(
    MultiTensorField3D<T, SymmetricTensorImpl<T, nDim>::n> &A,
    MultiTensorField3D<T, SymmetricTensorImpl<T, nDim>::n> &B, MultiScalarField3D<T> &result,
    Box3D domain)
{
    std::vector<MultiBlock3D *> fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional(
        new IndexContraction_SymmetricTensor_A_SymmetricTensor_B_functional3D<T, nDim>, domain,
        fields);
}

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T> > fullIndexContractionOfSymmetricTensors(
    MultiTensorField3D<T, SymmetricTensorImpl<T, nDim>::n> &A,
    MultiTensorField3D<T, SymmetricTensorImpl<T, nDim>::n> &B, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > result = generateMultiScalarField<T>(B, domain);
    fullIndexContractionOfSymmetricTensors(A, B, *result, domain);
    return result;
}

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T> > fullIndexContractionOfSymmetricTensors(
    MultiTensorField3D<T, SymmetricTensorImpl<T, nDim>::n> &A,
    MultiTensorField3D<T, SymmetricTensorImpl<T, nDim>::n> &B)
{
    return fullIndexContractionOfSymmetricTensors(A, B, A.getBoundingBox());
}

/* *************** TensorField - TensorField inplace operations *************** */

template <typename T, int nDim>
void addInPlace(MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B, Box3D domain)
{
    applyProcessingFunctional(new Tensor_A_plus_B_inplace_functional3D<T, nDim>, domain, A, B);
}

template <typename T, int nDim>
void addInPlace(MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B)
{
    addInPlace(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void subtractInPlace(MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B, Box3D domain)
{
    applyProcessingFunctional(new Tensor_A_minus_B_inplace_functional3D<T, nDim>, domain, A, B);
}

template <typename T, int nDim>
void subtractInPlace(MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B)
{
    subtractInPlace(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void multiplyInPlace(MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B, Box3D domain)
{
    applyProcessingFunctional(new Tensor_A_times_B_inplace_functional3D<T, nDim>, domain, A, B);
}

template <typename T, int nDim>
void multiplyInPlace(MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B)
{
    multiplyInPlace(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void multiplyInPlace(MultiTensorField3D<T, nDim> &A, T alpha, Box3D domain)
{
    applyProcessingFunctional(
        new Tensor_A_times_alpha_inplace_functional3D<T, nDim>(alpha), domain, A);
}

template <typename T, int nDim>
void multiplyInPlace(MultiTensorField3D<T, nDim> &A, T alpha)
{
    multiplyInPlace(A, alpha, A.getBoundingBox());
}

template <typename T, int nDim>
void divideInPlace(MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B, Box3D domain)
{
    applyProcessingFunctional(new Tensor_A_dividedBy_B_inplace_functional3D<T, nDim>, domain, A, B);
}

template <typename T, int nDim>
void divideInPlace(MultiTensorField3D<T, nDim> &A, MultiTensorField3D<T, nDim> &B)
{
    divideInPlace(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void divideInPlace(MultiTensorField3D<T, nDim> &A, MultiScalarField3D<T> &B, Box3D domain)
{
    applyProcessingFunctional(
        new Tensor_A_dividedBy_Scalar_B_inplace_functional3D<T, nDim>, domain, B, A);
}

template <typename T, int nDim>
void divideInPlace(MultiTensorField3D<T, nDim> &A, MultiScalarField3D<T> &B)
{
    divideInPlace(A, B, A.getBoundingBox());
}

template <typename T, int nDim>
void divideInPlace(
    MultiTensorField3D<T, nDim> &A, MultiScalarField3D<T> &B, MultiScalarField3D<int> &mask,
    int flag, Box3D domain)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&A);
    args.push_back(&B);
    args.push_back(&mask);
    applyProcessingFunctional(
        new Masked_Tensor_A_dividedBy_Scalar_B_inplace_functional3D<T, nDim>(flag), domain, args);
}

template <typename T, int nDim>
void divideInPlace(
    MultiTensorField3D<T, nDim> &A, MultiScalarField3D<T> &B, MultiScalarField3D<int> &mask,
    int flag)
{
    divideInPlace(A, B, mask, flag, A.getBoundingBox());
}

template <typename T, int nDim>
void normalizeInPlace(MultiTensorField3D<T, nDim> &data, Box3D domain)
{
    applyProcessingFunctional(new Normalize_Tensor_inplace_functional3D<T, nDim>(), domain, data);
}

template <typename T, int nDim>
void normalizeInPlace(MultiTensorField3D<T, nDim> &data)
{
    normalizeInPlace(data, data.getBoundingBox());
}

/* *************** LBMsmoothenTensor3D ******************************************* */

template <typename T, int nDim, template <typename U> class Descriptor>
void lbmSmoothenTensor(
    MultiTensorField3D<T, nDim> &data, MultiTensorField3D<T, nDim> &result, Box3D domain)
{
    applyProcessingFunctional(new LBMsmoothenTensor3D<T, nDim, Descriptor>, domain, data, result);
}

template <typename T, int nDim, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, nDim> > lbmSmoothenTensor(
    MultiTensorField3D<T, nDim> &data, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, nDim> > result =
        generateMultiTensorField<T, nDim>(data, domain);
    lbmSmoothenTensor<T, nDim, Descriptor>(data, *result, domain);
    return result;
}

template <typename T, int nDim, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, nDim> > lbmSmoothenTensor(MultiTensorField3D<T, nDim> &data)
{
    return lbmSmoothenTensor<T, nDim, Descriptor>(data, data.getBoundingBox());
}

/* *************** LBMsmoothenTensorInPlace3D ******************************************* */

template <typename T, int nDim, template <typename U> class Descriptor>
void lbmSmoothenTensorInPlace(MultiTensorField3D<T, nDim> &data, Box3D domain)
{
    applyProcessingFunctional(new LBMsmoothenTensorInPlace3D<T, nDim, Descriptor>, domain, data);
}

template <typename T, int nDim, template <typename U> class Descriptor>
void lbmSmoothenTensorInPlace(MultiTensorField3D<T, nDim> &data)
{
    lbmSmoothenTensorInPlace<T, nDim, Descriptor>(data, data.getBoundingBox());
}

/* *************** SmoothenTensor3D ******************************************* */

template <typename T, int nDim>
void smoothenTensor(
    MultiTensorField3D<T, nDim> &data, MultiTensorField3D<T, nDim> &result, Box3D domain)
{
    applyProcessingFunctional(new SmoothenTensor3D<T, nDim>, domain, data, result);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > smoothenTensor(
    MultiTensorField3D<T, nDim> &data, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, nDim> > result =
        generateMultiTensorField<T, nDim>(data, domain);
    smoothenTensor<T, nDim>(data, *result, domain);
    return result;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > smoothenTensor(MultiTensorField3D<T, nDim> &data)
{
    return smoothenTensor<T, nDim>(data, data.getBoundingBox());
}

/* *************** MollifyTensor3D ******************************************* */

template <typename T, int nDim>
void mollifyTensor(
    T l, plint d, Box3D globalDomain, int exclusionFlag, MultiTensorField3D<T, nDim> &data,
    MultiScalarField3D<int> &flag, MultiTensorField3D<T, nDim> &result, Box3D domain)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&data);
    args.push_back(&flag);
    args.push_back(&result);
    applyProcessingFunctional(
        new MollifyTensor3D<T, nDim>(l, d, globalDomain, exclusionFlag), domain, args);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > mollifyTensor(
    T l, plint d, Box3D globalDomain, int exclusionFlag, MultiTensorField3D<T, nDim> &data,
    MultiScalarField3D<int> &flag, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, nDim> > result =
        generateMultiTensorField<T, nDim>(data, domain);

    std::vector<MultiBlock3D *> args;
    args.push_back(&data);
    args.push_back(&flag);
    args.push_back(result.get());
    applyProcessingFunctional(
        new MollifyTensor3D<T, nDim>(l, d, globalDomain, exclusionFlag), domain, args);

    return result;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > mollifyTensor(
    T l, plint d, Box3D globalDomain, int exclusionFlag, MultiTensorField3D<T, nDim> &data,
    MultiScalarField3D<int> &flag)
{
    std::unique_ptr<MultiTensorField3D<T, nDim> > result =
        generateMultiTensorField<T, nDim>(data, data.getBoundingBox());

    std::vector<MultiBlock3D *> args;
    args.push_back(&data);
    args.push_back(&flag);
    args.push_back(result.get());
    applyProcessingFunctional(
        new MollifyTensor3D<T, nDim>(l, d, globalDomain, exclusionFlag), data.getBoundingBox(),
        args);

    return result;
}

/* *************** LBMcomputeDivergence3D ******************************************* */

template <typename T, template <typename U> class Descriptor>
void lbmComputeDivergence(
    MultiScalarField3D<T> &divergence, MultiTensorField3D<T, 3> &vectorField, Box3D domain)
{
    applyProcessingFunctional(
        new LBMcomputeDivergence3D<T, Descriptor>, domain, divergence, vectorField);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > lbmComputeDivergence(
    MultiTensorField3D<T, 3> &vectorField, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > divergence =
        generateMultiScalarField<T>(vectorField, domain);
    lbmComputeDivergence<T, Descriptor>(*divergence, vectorField, domain);
    return divergence;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > lbmComputeDivergence(MultiTensorField3D<T, 3> &vectorField)
{
    return lbmComputeDivergence<T, Descriptor>(vectorField, vectorField.getBoundingBox());
}

/* ********************************************************************* */
/* *************** PART VII. Multi-block wrappers: NTensor-field ******* */
/* ********************************************************************* */

/* *************** Extract Sub-NTensorField **************************** */

template <typename T>
void extractSubDomain(
    MultiNTensorField3D<T> &field, MultiNTensorField3D<T> &extractedField, Box3D domain)
{
    applyProcessingFunctional(
        new ExtractNTensorSubDomainFunctional3D<T>, domain, field, extractedField);
}

template <typename T>
std::unique_ptr<MultiNTensorField3D<T> > extractSubDomain(
    MultiNTensorField3D<T> &field, Box3D domain)
{
    MultiNTensorField3D<T> *extractedField =
        generateMultiNTensorField<T>(field, domain, field.getNdim());
    extractSubDomain(field, *extractedField, domain);
    return std::unique_ptr<MultiNTensorField3D<T> >(extractedField);
}

}  // namespace plb

#endif  // DATA_ANALYSIS_WRAPPER_3D_HH
