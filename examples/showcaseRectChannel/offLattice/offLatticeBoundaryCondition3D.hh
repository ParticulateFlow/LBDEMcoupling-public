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

#ifndef OFF_LATTICE_BOUNDARY_CONDITION_3D_HH
#define OFF_LATTICE_BOUNDARY_CONDITION_3D_HH

#include "core/globalDefs.h"
#include "offLattice/offLatticeBoundaryCondition3D.h"
#include "offLattice/offLatticeBoundaryProfiles3D.h"
#include "offLattice/triangularSurfaceMesh.h"
#include "triangleToDef.h"

namespace plb {

/* ********** OffLatticeBoundaryCondition3D ********************************** */

template <typename T, template <typename U> class Descriptor, class BoundaryType>
OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::OffLatticeBoundaryCondition3D(
    OffLatticeModel3D<T, BoundaryType> *offLatticeModel_, VoxelizedDomain3D<T> &voxelizedDomain_,
    MultiBlockLattice3D<T, Descriptor> &lattice_) :
    voxelizedDomain(voxelizedDomain_),
    lattice(lattice_),
    boundaryShapeArg(lattice_),
    offLatticeModel(offLatticeModel_),
    offLatticePattern(lattice)
{
    // It is very important that the "offLatticePattern" container block
    // has the same multi-block management as the lattice used in the
    // simulation.
    std::vector<MultiBlock3D *> offLatticeIniArg;
    // First argument for compute-off-lattice-pattern.
    offLatticeIniArg.push_back(&offLatticePattern);
    // Remaining arguments for inner-flow-shape.
    offLatticeIniArg.push_back(&voxelizedDomain.getVoxelMatrix());
    offLatticeIniArg.push_back(&voxelizedDomain.getTriangleHash());
    offLatticeIniArg.push_back(&boundaryShapeArg);
    applyProcessingFunctional(
        new OffLatticePatternFunctional3D<T, BoundaryType>(offLatticeModel->clone()),
        offLatticePattern.getBoundingBox(), offLatticeIniArg);
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::OffLatticeBoundaryCondition3D(
    OffLatticeModel3D<T, BoundaryType> *offLatticeModel_, VoxelizedDomain3D<T> &voxelizedDomain_,
    MultiBlockLattice3D<T, Descriptor> &lattice_,
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particleField_) :
    offLatticeModel(offLatticeModel_),
    voxelizedDomain(voxelizedDomain_),
    lattice(lattice_),
    boundaryShapeArg(particleField_),
    offLatticePattern(lattice)
{
    // It is very important that the "offLatticePattern" container block
    // has the same multi-block management as the lattice used in the
    // simulation.
    std::vector<MultiBlock3D *> offLatticeIniArg;
    // First argument for compute-off-lattice-pattern.
    offLatticeIniArg.push_back(&offLatticePattern);
    // Remaining arguments for inner-flow-shape.
    offLatticeIniArg.push_back(&voxelizedDomain.getVoxelMatrix());
    offLatticeIniArg.push_back(&voxelizedDomain.getTriangleHash());
    offLatticeIniArg.push_back(&boundaryShapeArg);
    applyProcessingFunctional(
        new OffLatticePatternFunctional3D<T, BoundaryType>(offLatticeModel->clone()),
        offLatticePattern.getBoundingBox(), offLatticeIniArg);
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::OffLatticeBoundaryCondition3D(
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType> const &rhs) :
    offLatticeModel(rhs.offLatticeModel.clone()),
    voxelizedDomain(rhs.voxelizedDomain),
    lattice(rhs.lattice),
    boundaryShapeArg(rhs.boundaryShapeArg),
    offLatticePattern(rhs.offLatticePattern)
{ }

template <typename T, template <typename U> class Descriptor, class BoundaryType>
OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::~OffLatticeBoundaryCondition3D()
{
    delete offLatticeModel;
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
void OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::insert(plint processorLevel)
{
    // It is very important that the "offLatticePattern" container block
    // has the same multi-block management as the lattice used in the
    // simulation.
    std::vector<MultiBlock3D *> offLatticeArg;
    // First two arguments for Guo algorithm.
    offLatticeArg.push_back(&lattice);
    offLatticeArg.push_back(&offLatticePattern);
    // Remaining arguments for inner-flow-shape.
    offLatticeArg.push_back(&voxelizedDomain.getVoxelMatrix());
    offLatticeArg.push_back(&voxelizedDomain.getTriangleHash());
    offLatticeArg.push_back(&boundaryShapeArg);
    plint numShapeArgs = 3;
    plint numCompletionArgs = 0;
    integrateProcessingFunctional(
        new OffLatticeCompletionFunctional3D<T, Descriptor, BoundaryType>(
            offLatticeModel->clone(), numShapeArgs, numCompletionArgs),
        boundaryShapeArg.getBoundingBox(), offLatticeArg, processorLevel);
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
void OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::insert(
    std::vector<MultiBlock3D *> const &completionArg, plint processorLevel)
{
    // It is very important that the "offLatticePattern" container block
    // has the same multi-block management as the lattice used in the
    // simulation.
    std::vector<MultiBlock3D *> offLatticeArg;
    // First three arguments for Guo algorithm.
    offLatticeArg.push_back(&lattice);
    offLatticeArg.push_back(&offLatticePattern);
    // Next arguments for inner-flow-shape.
    offLatticeArg.push_back(&voxelizedDomain.getVoxelMatrix());
    offLatticeArg.push_back(&voxelizedDomain.getTriangleHash());
    offLatticeArg.push_back(&boundaryShapeArg);
    // Remaining are optional arguments for completion algorithm.
    plint numCompletionArgs = (plint)completionArg.size();
    for (plint i = 0; i < numCompletionArgs; ++i) {
        offLatticeArg.push_back(completionArg[i]);
    }
    plint numShapeArgs = 3;
    integrateProcessingFunctional(
        new OffLatticeCompletionFunctional3D<T, Descriptor, BoundaryType>(
            offLatticeModel->clone(), numShapeArgs, numCompletionArgs),
        boundaryShapeArg.getBoundingBox(), offLatticeArg, processorLevel);
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
void OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::apply()
{
    // It is very important that the "offLatticePattern" container block
    // has the same multi-block management as the lattice used in the
    // simulation.
    std::vector<MultiBlock3D *> offLatticeArg;
    // First three arguments for Guo algorithm.
    offLatticeArg.push_back(&lattice);
    offLatticeArg.push_back(&offLatticePattern);
    // Remaining arguments for inner-flow-shape.
    offLatticeArg.push_back(&voxelizedDomain.getVoxelMatrix());
    offLatticeArg.push_back(&voxelizedDomain.getTriangleHash());
    offLatticeArg.push_back(&boundaryShapeArg);
    plint numShapeArgs = 3;
    plint numCompletionArgs = 0;
    applyProcessingFunctional(
        new OffLatticeCompletionFunctional3D<T, Descriptor, BoundaryType>(
            offLatticeModel->clone(), numShapeArgs, numCompletionArgs),
        boundaryShapeArg.getBoundingBox(), offLatticeArg);
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
void OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::apply(
    std::vector<MultiBlock3D *> const &completionArg)
{
    // It is very important that the "offLatticePattern" container block
    // has the same multi-block management as the lattice used in the
    // simulation.
    std::vector<MultiBlock3D *> offLatticeArg;
    // First three arguments for Guo algorithm.
    offLatticeArg.push_back(&lattice);
    offLatticeArg.push_back(&offLatticePattern);
    // Next arguments for inner-flow-shape.
    offLatticeArg.push_back(&voxelizedDomain.getVoxelMatrix());
    offLatticeArg.push_back(&voxelizedDomain.getTriangleHash());
    offLatticeArg.push_back(&boundaryShapeArg);
    // Remaining are optional arguments for completion algorithm.
    plint numCompletionArgs = (plint)completionArg.size();
    for (plint i = 0; i < numCompletionArgs; ++i) {
        offLatticeArg.push_back(completionArg[i]);
    }
    plint numShapeArgs = 3;
    applyProcessingFunctional(
        new OffLatticeCompletionFunctional3D<T, Descriptor, BoundaryType>(
            offLatticeModel->clone(), numShapeArgs, numCompletionArgs),
        boundaryShapeArg.getBoundingBox(), offLatticeArg);
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
Array<T, 3> OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::getForceOnObject()
{
    std::vector<MultiBlock3D *> arg;
    arg.push_back(&offLatticePattern);
    GetForceOnObjectFunctional3D<T, BoundaryType> functional(offLatticeModel->clone());
    applyProcessingFunctional(functional, boundaryShapeArg.getBoundingBox(), arg);
    return functional.getForce();
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
std::unique_ptr<MultiTensorField3D<T, 3> >
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeVelocity(Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, 3> > velocity(plb::computeVelocity(lattice, domain));
    int flowType = voxelizedDomain.getFlowType();
    int solidFlag = voxelFlag::invert(flowType);
    int solidBorderFlag = voxelFlag::borderFlag(solidFlag);
    setToConstant<T, 3>(
        *velocity, voxelizedDomain.getVoxelMatrix(), solidFlag, domain, Array<T, 3>(T(), T(), T()));
    setToConstant<T, 3>(
        *velocity, voxelizedDomain.getVoxelMatrix(), solidBorderFlag, domain,
        Array<T, 3>(T(), T(), T()));
    return velocity;
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
std::unique_ptr<MultiTensorField3D<T, 3> >
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeVelocity()
{
    return computeVelocity(lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
std::unique_ptr<MultiTensorField3D<T, 3> >
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeVorticity(Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, 3> > vorticity(
        plb::computeBulkVorticity(*plb::computeVelocity(lattice, domain), domain));
    int flowType = voxelizedDomain.getFlowType();
    int solidFlag = voxelFlag::invert(flowType);
    int solidBorderFlag = voxelFlag::borderFlag(solidFlag);
    setToConstant<T, 3>(
        *vorticity, voxelizedDomain.getVoxelMatrix(), solidFlag, domain,
        Array<T, 3>(T(), T(), T()));
    setToConstant<T, 3>(
        *vorticity, voxelizedDomain.getVoxelMatrix(), solidBorderFlag, domain,
        Array<T, 3>(T(), T(), T()));
    return vorticity;
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
std::unique_ptr<MultiTensorField3D<T, 3> >
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeVorticity()
{
    return computeVorticity(lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
std::unique_ptr<MultiScalarField3D<T> >
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeVelocityNorm(Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > velocityNorm(plb::computeVelocityNorm(lattice, domain));
    int flowType = voxelizedDomain.getFlowType();
    int solidFlag = voxelFlag::invert(flowType);
    int solidBorderFlag = voxelFlag::borderFlag(solidFlag);
    setToConstant(*velocityNorm, voxelizedDomain.getVoxelMatrix(), solidFlag, domain, (T)0);
    setToConstant(*velocityNorm, voxelizedDomain.getVoxelMatrix(), solidBorderFlag, domain, (T)0);
    return velocityNorm;
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
std::unique_ptr<MultiScalarField3D<T> >
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeVelocityNorm()
{
    return computeVelocityNorm(lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
std::unique_ptr<MultiScalarField3D<T> >
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeVorticityNorm(Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > vorticityNorm(plb::computeNorm(
        *plb::computeBulkVorticity(*plb::computeVelocity(lattice, domain), domain), domain));
    int flowType = voxelizedDomain.getFlowType();
    int solidFlag = voxelFlag::invert(flowType);
    int solidBorderFlag = voxelFlag::borderFlag(solidFlag);
    setToConstant<T>(*vorticityNorm, voxelizedDomain.getVoxelMatrix(), solidFlag, domain, T());
    setToConstant<T>(
        *vorticityNorm, voxelizedDomain.getVoxelMatrix(), solidBorderFlag, domain, T());
    return vorticityNorm;
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
std::unique_ptr<MultiScalarField3D<T> >
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeVorticityNorm()
{
    return computeVorticityNorm(lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
std::unique_ptr<MultiScalarField3D<T> >
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeVelocityComponent(
        Box3D domain, plint iComp)
{
    std::unique_ptr<MultiScalarField3D<T> > velocityComponent(
        plb::computeVelocityComponent(lattice, domain, iComp));
    int flowType = voxelizedDomain.getFlowType();
    int solidFlag = voxelFlag::invert(flowType);
    int solidBorderFlag = voxelFlag::borderFlag(solidFlag);
    setToConstant(*velocityComponent, voxelizedDomain.getVoxelMatrix(), solidFlag, domain, (T)0);
    setToConstant(
        *velocityComponent, voxelizedDomain.getVoxelMatrix(), solidBorderFlag, domain, (T)0);
    return velocityComponent;
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
std::unique_ptr<MultiScalarField3D<T> >
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeVelocityComponent(
        plint iComp)
{
    return computeVelocityComponent(lattice.getBoundingBox(), iComp);
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
std::unique_ptr<MultiScalarField3D<T> >
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computePressure(Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > pressure(plb::computeDensity(lattice, domain));
    T averageDensity = computeAverageDensity(domain);
    subtractInPlace(*pressure, averageDensity, domain);
    multiplyInPlace(*pressure, Descriptor<T>::cs2, domain);
    int flowType = voxelizedDomain.getFlowType();
    int solidFlag = voxelFlag::invert(flowType);
    int solidBorderFlag = voxelFlag::borderFlag(solidFlag);
    setToConstant(*pressure, voxelizedDomain.getVoxelMatrix(), solidFlag, domain, (T)0);
    setToConstant(*pressure, voxelizedDomain.getVoxelMatrix(), solidBorderFlag, domain, (T)0);
    return pressure;
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
std::unique_ptr<MultiScalarField3D<T> >
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computePressure()
{
    return computePressure(lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
std::unique_ptr<MultiScalarField3D<T> >
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeDensity(
        Box3D domain, T solidDensity)
{
    std::unique_ptr<MultiScalarField3D<T> > density(plb::computeDensity(lattice, domain));
    int flowType = voxelizedDomain.getFlowType();
    int solidFlag = voxelFlag::invert(flowType);
    int solidBorderFlag = voxelFlag::borderFlag(solidFlag);
    setToConstant(*density, voxelizedDomain.getVoxelMatrix(), solidFlag, domain, solidDensity);
    setToConstant(
        *density, voxelizedDomain.getVoxelMatrix(), solidBorderFlag, domain, solidDensity);
    return density;
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
std::unique_ptr<MultiScalarField3D<T> >
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeDensity(T solidDensity)
{
    return computeDensity(lattice.getBoundingBox(), solidDensity);
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
std::unique_ptr<MultiScalarField3D<T> >
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeStrainRateNorm(Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > strainRateNorm(
        computeSymmetricTensorNorm(*plb::computeStrainRateFromStress(lattice, domain)));
    int flowType = voxelizedDomain.getFlowType();
    int solidFlag = voxelFlag::invert(flowType);
    int solidBorderFlag = voxelFlag::borderFlag(solidFlag);
    setToConstant(*strainRateNorm, voxelizedDomain.getVoxelMatrix(), solidFlag, domain, (T)0);
    setToConstant(*strainRateNorm, voxelizedDomain.getVoxelMatrix(), solidBorderFlag, domain, (T)0);
    return strainRateNorm;
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
std::unique_ptr<MultiScalarField3D<T> >
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeStrainRateNorm()
{
    return computeStrainRateNorm(lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> >
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeStrainRate(Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > strainRate(
        plb::computeStrainRateFromStress(lattice, domain));
    int flowType = voxelizedDomain.getFlowType();
    int solidFlag = voxelFlag::invert(flowType);
    int solidBorderFlag = voxelFlag::borderFlag(solidFlag);
    Array<T, SymmetricTensor<T, Descriptor>::n> zeros;
    zeros.resetToZero();
    setToConstant<T, SymmetricTensor<T, Descriptor>::n>(
        *strainRate, voxelizedDomain.getVoxelMatrix(), solidFlag, domain, zeros);
    setToConstant<T, SymmetricTensor<T, Descriptor>::n>(
        *strainRate, voxelizedDomain.getVoxelMatrix(), solidBorderFlag, domain, zeros);
    return strainRate;
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> >
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeStrainRate()
{
    return computeStrainRate(lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
std::unique_ptr<MultiScalarField3D<T> >
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeShearStressNorm(Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > shearStressNorm(
        computeSymmetricTensorNorm(*plb::computeShearStress(lattice, domain)));
    int flowType = voxelizedDomain.getFlowType();
    int solidFlag = voxelFlag::invert(flowType);
    int solidBorderFlag = voxelFlag::borderFlag(solidFlag);
    setToConstant(*shearStressNorm, voxelizedDomain.getVoxelMatrix(), solidFlag, domain, (T)0);
    setToConstant(
        *shearStressNorm, voxelizedDomain.getVoxelMatrix(), solidBorderFlag, domain, (T)0);
    return shearStressNorm;
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
std::unique_ptr<MultiScalarField3D<T> >
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeShearStressNorm()
{
    return computeShearStressNorm(lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> >
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeShearStress(Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > shearStress(
        plb::computeShearStress(lattice, domain));
    int flowType = voxelizedDomain.getFlowType();
    int solidFlag = voxelFlag::invert(flowType);
    int solidBorderFlag = voxelFlag::borderFlag(solidFlag);
    Array<T, SymmetricTensor<T, Descriptor>::n> zeros;
    zeros.resetToZero();
    setToConstant<T, SymmetricTensor<T, Descriptor>::n>(
        *shearStress, voxelizedDomain.getVoxelMatrix(), solidFlag, domain, zeros);
    setToConstant<T, SymmetricTensor<T, Descriptor>::n>(
        *shearStress, voxelizedDomain.getVoxelMatrix(), solidBorderFlag, domain, zeros);
    return shearStress;
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> >
    OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeShearStress()
{
    return computeShearStress(lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
T OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeAverageVelocityComponent(
    Box3D domain, plint iComponent)
{
    std::unique_ptr<MultiScalarField3D<T> > density(
        plb::computeVelocityComponent(lattice, domain, iComponent));
    MultiScalarField3D<int> flagMatrix((MultiBlock3D &)voxelizedDomain.getVoxelMatrix());
    int flowType = voxelizedDomain.getFlowType();
    int fluidBorderFlag = voxelFlag::borderFlag(flowType);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(), flowType, domain, 1);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(), fluidBorderFlag, domain, 1);
    return computeAverage(*density, flagMatrix, 1, domain);
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
Array<T, 3> OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeAverageVelocity(
    Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, 3> > velocity(plb::computeVelocity(lattice, domain));
    MultiScalarField3D<int> flagMatrix((MultiBlock3D &)voxelizedDomain.getVoxelMatrix());
    int flowType = voxelizedDomain.getFlowType();
    int fluidBorderFlag = voxelFlag::borderFlag(flowType);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(), flowType, domain, 1);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(), fluidBorderFlag, domain, 1);
    return computeAverage<T, 3>(*velocity, flagMatrix, 1, domain);
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
T OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeAverageDensity()
{
    return computeAverageDensity(voxelizedDomain.getVoxelMatrix().getBoundingBox());
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
T OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeAverageDensity(Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > density(plb::computeDensity(lattice, domain));
    std::unique_ptr<MultiScalarField3D<T> > density2(plb::computeDensity(lattice, domain));
    MultiScalarField3D<int> flagMatrix((MultiBlock3D &)voxelizedDomain.getVoxelMatrix());
    int flowType = voxelizedDomain.getFlowType();
    int fluidBorderFlag = voxelFlag::borderFlag(flowType);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(), flowType, domain, 1);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(), fluidBorderFlag, domain, 1);
    return computeAverage(*density, flagMatrix, 1, domain);
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
T OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeAverageEnergy()
{
    return computeAverageEnergy(voxelizedDomain.getVoxelMatrix().getBoundingBox());
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
T OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeAverageEnergy(Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > energy(plb::computeKineticEnergy(lattice, domain));
    MultiScalarField3D<int> flagMatrix((MultiBlock3D &)voxelizedDomain.getVoxelMatrix());
    int flowType = voxelizedDomain.getFlowType();
    int fluidBorderFlag = voxelFlag::borderFlag(flowType);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(), flowType, domain, 1);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(), fluidBorderFlag, domain, 1);
    return computeAverage(*energy, flagMatrix, 1, domain);
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
T OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeRMSvorticity()
{
    return computeRMSvorticity(voxelizedDomain.getVoxelMatrix().getBoundingBox());
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
T OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeRMSvorticity(Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > vorticityNormSqr(
        plb::computeNormSqr(*plb::computeBulkVorticity(*plb::computeVelocity(lattice, domain))));
    MultiScalarField3D<int> flagMatrix((MultiBlock3D &)voxelizedDomain.getVoxelMatrix());
    int flowType = voxelizedDomain.getFlowType();
    int fluidBorderFlag = voxelFlag::borderFlag(flowType);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(), flowType, domain, 1);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(), fluidBorderFlag, domain, 1);
    return std::sqrt(computeAverage(*vorticityNormSqr, flagMatrix, 1, domain));
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
T OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeAverageShearStressNorm()
{
    return computeAverageShearStressNorm(voxelizedDomain.getVoxelMatrix().getBoundingBox());
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
T OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeAverageShearStressNorm(
    Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > shearStress(
        plb::computeSymmetricTensorNorm(*plb::computeShearStress(lattice, domain)));
    MultiScalarField3D<int> flagMatrix((MultiBlock3D &)voxelizedDomain.getVoxelMatrix());
    int flowType = voxelizedDomain.getFlowType();
    int fluidBorderFlag = voxelFlag::borderFlag(flowType);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(), flowType, domain, 1);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(), fluidBorderFlag, domain, 1);
    return computeAverage(*shearStress, flagMatrix, 1, domain);
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
T OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeRMSshearStressNorm()
{
    return computeRMSshearStressNorm(voxelizedDomain.getVoxelMatrix().getBoundingBox());
}

template <typename T, template <typename U> class Descriptor, class BoundaryType>
T OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType>::computeRMSshearStressNorm(
    Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > shearStressNorm(
        plb::computeSymmetricTensorNorm(*plb::computeShearStress(lattice, domain)));
    T avgShearStress = computeAverageShearStressNorm(domain);

    MultiScalarField3D<int> flagMatrix((MultiBlock3D &)voxelizedDomain.getVoxelMatrix());
    int flowType = voxelizedDomain.getFlowType();
    int fluidBorderFlag = voxelFlag::borderFlag(flowType);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(), flowType, domain, 1);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(), fluidBorderFlag, domain, 1);
    shearStressNorm = subtract(*shearStressNorm, avgShearStress);
    return std::sqrt(
        computeAverage(*multiply(*shearStressNorm, *shearStressNorm), flagMatrix, 1, domain));
}

}  // namespace plb

#endif  // OFF_LATTICE_BOUNDARY_CONDITION_3D_HH
