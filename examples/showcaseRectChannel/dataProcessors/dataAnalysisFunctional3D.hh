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
 * Data processors for data analysis -- header file.
 */

#ifndef DATA_ANALYSIS_FUNCTIONAL_3D_HH
#define DATA_ANALYSIS_FUNCTIONAL_3D_HH

#include <cmath>
#include <limits>

#include "algorithm/linearAlgebra.h"
#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "core/blockStatistics.h"
#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "core/util.h"
#include "dataProcessors/dataAnalysisFunctional3D.h"
#include "finiteDifference/fdStencils1D.h"
#include "finiteDifference/fdWeights.h"
#include "finiteDifference/fdWeights.hh"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#ifndef PLB_BGP
#ifdef PLB_USE_EIGEN
#include <Eigen3/Eigenvalues>
#endif
#endif

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

/* *************** Reductive Data Functionals for BlockLattice ******* */

template <typename T, template <typename U> class Descriptor>
BoxSumRhoBarFunctional3D<T, Descriptor>::BoxSumRhoBarFunctional3D() :
    sumRhoBarId(this->getStatistics().subscribeSum())
{ }

template <typename T, template <typename U> class Descriptor>
void BoxSumRhoBarFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> const &cell = lattice.get(iX, iY, iZ);
                statistics.gatherSum(sumRhoBarId, cell.getDynamics().computeRhoBar(cell));
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxSumRhoBarFunctional3D<T, Descriptor> *BoxSumRhoBarFunctional3D<T, Descriptor>::clone() const
{
    return new BoxSumRhoBarFunctional3D(*this);
}

template <typename T, template <typename U> class Descriptor>
T BoxSumRhoBarFunctional3D<T, Descriptor>::getSumRhoBar() const
{
    return this->getStatistics().getSum(sumRhoBarId);
}

template <typename T, template <typename U> class Descriptor>
BoxSumEnergyFunctional3D<T, Descriptor>::BoxSumEnergyFunctional3D() :
    sumEnergyId(this->getStatistics().subscribeSum())
{ }

template <typename T, template <typename U> class Descriptor>
void BoxSumEnergyFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, Descriptor<T>::d> velocity;
                lattice.get(iX, iY, iZ).computeVelocity(velocity);
                T uNormSqr = VectorTemplate<T, Descriptor>::normSqr(velocity);
                statistics.gatherSum(sumEnergyId, uNormSqr);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxSumEnergyFunctional3D<T, Descriptor> *BoxSumEnergyFunctional3D<T, Descriptor>::clone() const
{
    return new BoxSumEnergyFunctional3D(*this);
}

template <typename T, template <typename U> class Descriptor>
T BoxSumEnergyFunctional3D<T, Descriptor>::getSumEnergy() const
{
    return this->getStatistics().getSum(sumEnergyId) / (T)2;
}

template <typename T, template <typename U> class Descriptor>
void BoxSumEnergyFunctional3D<T, Descriptor>::getDimensionsX(std::vector<int> &dimensions) const
{
    dimensions.resize(1);
    dimensions[0] = 2;
}

template <typename T, template <typename U> class Descriptor>
void BoxSumEnergyFunctional3D<T, Descriptor>::getDimensionsT(std::vector<int> &dimensions) const
{
    dimensions.resize(1);
    dimensions[0] = -2;
}

template <typename T, template <typename U> class Descriptor>
BoxSumForcedEnergyFunctional3D<T, Descriptor>::BoxSumForcedEnergyFunctional3D() :
    sumEnergyId(this->getStatistics().subscribeSum())
{ }

template <typename T, template <typename U> class Descriptor>
void BoxSumForcedEnergyFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, TensorField3D<T, Descriptor<T>::d> &force)
{
    Dot3D offset = computeRelativeDisplacement(lattice, force);
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, Descriptor<T>::d> velocity;
                lattice.get(iX, iY, iZ).computeVelocity(velocity);
                if (lattice.get(iX, iY, iZ).getDynamics().hasMoments()) {
                    Array<T, Descriptor<T>::d> const &f =
                        force.get(iX + offset.x, iY + offset.y, iZ + offset.z);
                    velocity[0] += (T)0.5 * f[0];
                    velocity[1] += (T)0.5 * f[1];
                    velocity[2] += (T)0.5 * f[2];
                }
                T uNormSqr = VectorTemplate<T, Descriptor>::normSqr(velocity);
                statistics.gatherSum(sumEnergyId, uNormSqr);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxSumForcedEnergyFunctional3D<T, Descriptor>
    *BoxSumForcedEnergyFunctional3D<T, Descriptor>::clone() const
{
    return new BoxSumForcedEnergyFunctional3D(*this);
}

template <typename T, template <typename U> class Descriptor>
T BoxSumForcedEnergyFunctional3D<T, Descriptor>::getSumEnergy() const
{
    return this->getStatistics().getSum(sumEnergyId) / (T)2;
}

template <typename T, template <typename U> class Descriptor>
void BoxSumForcedEnergyFunctional3D<T, Descriptor>::getDimensionsX(
    std::vector<int> &dimensions) const
{
    dimensions.resize(1);
    dimensions[0] = 2;
}

template <typename T, template <typename U> class Descriptor>
void BoxSumForcedEnergyFunctional3D<T, Descriptor>::getDimensionsT(
    std::vector<int> &dimensions) const
{
    dimensions.resize(1);
    dimensions[0] = -2;
}

template <typename T, template <typename U> class Descriptor>
BoxSumConstForcedEnergyFunctional3D<T, Descriptor>::BoxSumConstForcedEnergyFunctional3D(
    Array<T, Descriptor<T>::d> force_) :
    force(force_), sumEnergyId(this->getStatistics().subscribeSum())
{ }

template <typename T, template <typename U> class Descriptor>
void BoxSumConstForcedEnergyFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, Descriptor<T>::d> velocity;
                lattice.get(iX, iY, iZ).computeVelocity(velocity);
                if (lattice.get(iX, iY, iZ).getDynamics().hasMoments()) {
                    velocity[0] += (T)0.5 * force[0];
                    velocity[1] += (T)0.5 * force[1];
                    velocity[2] += (T)0.5 * force[2];
                }
                T uNormSqr = VectorTemplate<T, Descriptor>::normSqr(velocity);
                statistics.gatherSum(sumEnergyId, uNormSqr);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxSumConstForcedEnergyFunctional3D<T, Descriptor>
    *BoxSumConstForcedEnergyFunctional3D<T, Descriptor>::clone() const
{
    return new BoxSumConstForcedEnergyFunctional3D(*this);
}

template <typename T, template <typename U> class Descriptor>
T BoxSumConstForcedEnergyFunctional3D<T, Descriptor>::getSumEnergy() const
{
    return this->getStatistics().getSum(sumEnergyId) / (T)2;
}

template <typename T, template <typename U> class Descriptor>
void BoxSumConstForcedEnergyFunctional3D<T, Descriptor>::getDimensionsX(
    std::vector<int> &dimensions) const
{
    dimensions.resize(1);
    dimensions[0] = 2;
}

template <typename T, template <typename U> class Descriptor>
void BoxSumConstForcedEnergyFunctional3D<T, Descriptor>::getDimensionsT(
    std::vector<int> &dimensions) const
{
    dimensions.resize(1);
    dimensions[0] = -2;
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
BoxSumCustomForcedEnergyFunctional3D<
    T, Descriptor, ForceFunction>::BoxSumCustomForcedEnergyFunctional3D(ForceFunction f_) :
    f(f_), sumEnergyId(this->getStatistics().subscribeSum())
{ }

template <typename T, template <typename U> class Descriptor, class ForceFunction>
void BoxSumCustomForcedEnergyFunctional3D<T, Descriptor, ForceFunction>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    BlockStatistics &statistics = this->getStatistics();
    Dot3D location = lattice.getLocation();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint x = iX + location.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint y = iY + location.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint z = iZ + location.z;
                Array<T, Descriptor<T>::d> velocity;
                lattice.get(iX, iY, iZ).computeVelocity(velocity);
                if (lattice.get(iX, iY, iZ).getDynamics().hasMoments()) {
                    Array<T, Descriptor<T>::d> force;
                    f(x, y, z, force);
                    velocity[0] += (T)0.5 * force[0];
                    velocity[1] += (T)0.5 * force[1];
                    velocity[2] += (T)0.5 * force[2];
                }
                T uNormSqr = VectorTemplate<T, Descriptor>::normSqr(velocity);
                statistics.gatherSum(sumEnergyId, uNormSqr);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
BoxSumCustomForcedEnergyFunctional3D<T, Descriptor, ForceFunction>
    *BoxSumCustomForcedEnergyFunctional3D<T, Descriptor, ForceFunction>::clone() const
{
    return new BoxSumCustomForcedEnergyFunctional3D(*this);
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
T BoxSumCustomForcedEnergyFunctional3D<T, Descriptor, ForceFunction>::getSumEnergy() const
{
    return this->getStatistics().getSum(sumEnergyId) / (T)2;
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
void BoxSumCustomForcedEnergyFunctional3D<T, Descriptor, ForceFunction>::getDimensionsX(
    std::vector<int> &dimensions) const
{
    dimensions.resize(1);
    dimensions[0] = 2;
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
void BoxSumCustomForcedEnergyFunctional3D<T, Descriptor, ForceFunction>::getDimensionsT(
    std::vector<int> &dimensions) const
{
    dimensions.resize(1);
    dimensions[0] = -2;
}

/* ******** ScalarFieldSingleProbe3D *********************************** */

template <typename T>
ScalarFieldSingleProbe3D<T>::ScalarFieldSingleProbe3D(std::vector<Array<T, 3> > const &positions_) :
    positions(positions_)
{
    scalarIds.resize(positions.size());
    for (pluint i = 0; i < positions.size(); ++i) {
        scalarIds[i] = this->getStatistics().subscribeSum();
    }
}

template <typename T>
void ScalarFieldSingleProbe3D<T>::process(Box3D domain, ScalarField3D<T> &scalarField)
{
    std::vector<Dot3D> cellPos(8);
    std::vector<T> weights(8);
    T scalar;
    for (pluint i = 0; i < positions.size(); ++i) {
        scalar = (T)0;
        Array<T, 3> position(positions[i]);
        Dot3D referenceCellPos((plint)position[0], (plint)position[1], (plint)position[2]);
        referenceCellPos -= scalarField.getLocation();
        if (contained(referenceCellPos, domain)) {
            linearInterpolationCoefficients(scalarField, position, cellPos, weights);
            for (plint iCell = 0; iCell < 8; ++iCell) {
                T cellScalar =
                    scalarField.get(cellPos[iCell].x, cellPos[iCell].y, cellPos[iCell].z);
                scalar += weights[iCell] * cellScalar;
            }
        }
        this->getStatistics().gatherSum(scalarIds[i], scalar);
    }
}

template <typename T>
ScalarFieldSingleProbe3D<T> *ScalarFieldSingleProbe3D<T>::clone() const
{
    return new ScalarFieldSingleProbe3D<T>(*this);
}

template <typename T>
std::vector<T> ScalarFieldSingleProbe3D<T>::getScalars() const
{
    std::vector<T> scalars(positions.size());
    for (pluint i = 0; i < positions.size(); ++i) {
        scalars[i] = this->getStatistics().getSum(scalarIds[i]);
    }
    return scalars;
}

/* ******** DensitySingleProbe3D *********************************** */

template <typename T, template <typename U> class Descriptor>
DensitySingleProbe3D<T, Descriptor>::DensitySingleProbe3D(
    std::vector<Array<T, 3> > const &positions_) :
    positions(positions_)
{
    densityIds.resize(positions.size());
    for (pluint i = 0; i < positions.size(); ++i) {
        densityIds[i] = this->getStatistics().subscribeSum();
    }
}

template <typename T, template <typename U> class Descriptor>
void DensitySingleProbe3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    std::vector<Dot3D> cellPos(8);
    std::vector<T> weights(8);
    T density;
    for (pluint i = 0; i < positions.size(); ++i) {
        density = (T)0;
        Array<T, 3> position(positions[i]);
        Dot3D referenceCellPos((plint)position[0], (plint)position[1], (plint)position[2]);
        referenceCellPos -= lattice.getLocation();
        if (contained(referenceCellPos, domain)) {
            linearInterpolationCoefficients(lattice, position, cellPos, weights);
            for (plint iCell = 0; iCell < 8; ++iCell) {
                Cell<T, Descriptor> const &cell =
                    lattice.get(cellPos[iCell].x, cellPos[iCell].y, cellPos[iCell].z);
                T cellDensity = cell.computeDensity();
                density += weights[iCell] * cellDensity;
            }
        }
        this->getStatistics().gatherSum(densityIds[i], density);
    }
}

template <typename T, template <typename U> class Descriptor>
DensitySingleProbe3D<T, Descriptor> *DensitySingleProbe3D<T, Descriptor>::clone() const
{
    return new DensitySingleProbe3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
std::vector<T> DensitySingleProbe3D<T, Descriptor>::getDensities() const
{
    std::vector<T> densities(positions.size());
    for (pluint i = 0; i < positions.size(); ++i) {
        densities[i] = this->getStatistics().getSum(densityIds[i]);
    }
    return densities;
}

/* ******** InternalDensitySingleProbe3D *********************************** */

template <typename T, template <typename U> class Descriptor>
InternalDensitySingleProbe3D<T, Descriptor>::InternalDensitySingleProbe3D(
    std::vector<Array<T, 3> > const &positions_, std::vector<plint> const &ids_) :
    positions(positions_), ids(ids_)
{
    PLB_ASSERT(ids.size() == positions.size());
}

template <typename T, template <typename U> class Descriptor>
void InternalDensitySingleProbe3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    std::vector<Dot3D> cellPos(8);
    std::vector<T> weights(8);
    T density;
    for (pluint i = 0; i < positions.size(); ++i) {
        density = (T)0;
        Array<T, 3> position(positions[i]);
        Dot3D referenceCellPos((plint)position[0], (plint)position[1], (plint)position[2]);
        referenceCellPos -= lattice.getLocation();
        if (contained(referenceCellPos, domain)) {
            linearInterpolationCoefficients(lattice, position, cellPos, weights);
            for (plint iCell = 0; iCell < 8; ++iCell) {
                Cell<T, Descriptor> const &cell =
                    lattice.get(cellPos[iCell].x, cellPos[iCell].y, cellPos[iCell].z);
                T cellDensity = cell.computeDensity();
                density += weights[iCell] * cellDensity;
            }
        }
        lattice.getInternalStatistics().gatherSum(ids[i], density);
    }
}

template <typename T, template <typename U> class Descriptor>
InternalDensitySingleProbe3D<T, Descriptor> *InternalDensitySingleProbe3D<T, Descriptor>::clone()
    const
{
    return new InternalDensitySingleProbe3D<T, Descriptor>(*this);
}

/* ******** VelocitySingleProbe3D *********************************** */

template <typename T, template <typename U> class Descriptor>
VelocitySingleProbe3D<T, Descriptor>::VelocitySingleProbe3D(
    std::vector<Array<T, 3> > const &positions_) :
    positions(positions_)
{
    velIds.resize(positions.size());
    for (pluint iVel = 0; iVel < positions.size(); ++iVel) {
        velIds[iVel][0] = this->getStatistics().subscribeSum();
        velIds[iVel][1] = this->getStatistics().subscribeSum();
        velIds[iVel][2] = this->getStatistics().subscribeSum();
    }
}

template <typename T, template <typename U> class Descriptor>
void VelocitySingleProbe3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    std::vector<Dot3D> cellPos(8);
    std::vector<T> weights(8);
    Array<T, 3> velocity;
    for (pluint iVel = 0; iVel < positions.size(); ++iVel) {
        velocity.resetToZero();
        Array<T, 3> position(positions[iVel]);
        Dot3D referenceCellPos((plint)position[0], (plint)position[1], (plint)position[2]);
        referenceCellPos -= lattice.getLocation();
        if (contained(referenceCellPos, domain)) {
            linearInterpolationCoefficients(lattice, position, cellPos, weights);
            for (plint iCell = 0; iCell < 8; ++iCell) {
                Cell<T, Descriptor> const &cell =
                    lattice.get(cellPos[iCell].x, cellPos[iCell].y, cellPos[iCell].z);
                Array<T, 3> cellVelocity;
                cell.computeVelocity(cellVelocity);
                velocity += weights[iCell] * cellVelocity;
            }
        }
        this->getStatistics().gatherSum(velIds[iVel][0], velocity[0]);
        this->getStatistics().gatherSum(velIds[iVel][1], velocity[1]);
        this->getStatistics().gatherSum(velIds[iVel][2], velocity[2]);
    }
}

template <typename T, template <typename U> class Descriptor>
VelocitySingleProbe3D<T, Descriptor> *VelocitySingleProbe3D<T, Descriptor>::clone() const
{
    return new VelocitySingleProbe3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
std::vector<Array<T, 3> > VelocitySingleProbe3D<T, Descriptor>::getVelocities() const
{
    std::vector<Array<T, 3> > velocities(positions.size());
    for (pluint iVel = 0; iVel < positions.size(); ++iVel) {
        velocities[iVel][0] = this->getStatistics().getSum(velIds[iVel][0]);
        velocities[iVel][1] = this->getStatistics().getSum(velIds[iVel][1]);
        velocities[iVel][2] = this->getStatistics().getSum(velIds[iVel][2]);
    }
    return velocities;
}

/* ******** InternalVelocitySingleProbe3D *********************************** */

template <typename T, template <typename U> class Descriptor>
InternalVelocitySingleProbe3D<T, Descriptor>::InternalVelocitySingleProbe3D(
    std::vector<Array<T, 3> > const &positions_, std::vector<plint> const &velIds_) :
    positions(positions_), velIds(velIds_)
{
    PLB_ASSERT(positions.size() == velIds.size() / Descriptor<T>::d);
    PLB_ASSERT(velIds.size() % Descriptor<T>::d == 0);
}

template <typename T, template <typename U> class Descriptor>
void InternalVelocitySingleProbe3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    std::vector<Dot3D> cellPos(8);
    std::vector<T> weights(8);
    Array<T, 3> velocity;
    for (pluint iVel = 0; iVel < positions.size(); ++iVel) {
        velocity.resetToZero();
        Array<T, 3> position(positions[iVel]);
        Dot3D referenceCellPos((plint)position[0], (plint)position[1], (plint)position[2]);
        referenceCellPos -= lattice.getLocation();
        if (contained(referenceCellPos, domain)) {
            linearInterpolationCoefficients(lattice, position, cellPos, weights);
            for (plint iCell = 0; iCell < 8; ++iCell) {
                Cell<T, Descriptor> const &cell =
                    lattice.get(cellPos[iCell].x, cellPos[iCell].y, cellPos[iCell].z);
                Array<T, 3> cellVelocity;
                cell.computeVelocity(cellVelocity);
                velocity += weights[iCell] * cellVelocity;
            }
        }
        for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
            lattice.getInternalStatistics().gatherSum(velIds[3 * iVel + iD], velocity[iD]);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
InternalVelocitySingleProbe3D<T, Descriptor> *InternalVelocitySingleProbe3D<T, Descriptor>::clone()
    const
{
    return new InternalVelocitySingleProbe3D<T, Descriptor>(*this);
}

/* ******** VorticitySingleProbe3D *********************************** */

template <typename T, template <typename U> class Descriptor>
VorticitySingleProbe3D<T, Descriptor>::VorticitySingleProbe3D(
    std::vector<Array<T, 3> > const &positions_) :
    positions(positions_)
{
    vorticityIds.resize(positions.size());
    for (pluint iVel = 0; iVel < positions.size(); ++iVel) {
        vorticityIds[iVel][0] = this->getStatistics().subscribeSum();
        vorticityIds[iVel][1] = this->getStatistics().subscribeSum();
        vorticityIds[iVel][2] = this->getStatistics().subscribeSum();
    }
}

template <typename T, template <typename U> class Descriptor>
void VorticitySingleProbe3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    Array<T, 3> vorticity;
    for (pluint iVort = 0; iVort < positions.size(); ++iVort) {
        vorticity.resetToZero();
        Array<T, 3> position(positions[iVort]);
        Dot3D referenceCellPos((plint)position[0], (plint)position[1], (plint)position[2]);
        referenceCellPos -= lattice.getLocation();
        if (contained(referenceCellPos, domain)) {
            vorticity = fdLattice::firstOrderBulkVorticity(
                lattice, referenceCellPos.x, referenceCellPos.y, referenceCellPos.z);
        }
        this->getStatistics().gatherSum(vorticityIds[iVort][0], vorticity[0]);
        this->getStatistics().gatherSum(vorticityIds[iVort][1], vorticity[1]);
        this->getStatistics().gatherSum(vorticityIds[iVort][2], vorticity[2]);
    }
}

template <typename T, template <typename U> class Descriptor>
VorticitySingleProbe3D<T, Descriptor> *VorticitySingleProbe3D<T, Descriptor>::clone() const
{
    return new VorticitySingleProbe3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
std::vector<Array<T, 3> > VorticitySingleProbe3D<T, Descriptor>::getVorticities() const
{
    std::vector<Array<T, 3> > vorticities(positions.size());
    for (pluint iVort = 0; iVort < positions.size(); ++iVort) {
        vorticities[iVort][0] = this->getStatistics().getSum(vorticityIds[iVort][0]);
        vorticities[iVort][1] = this->getStatistics().getSum(vorticityIds[iVort][1]);
        vorticities[iVort][2] = this->getStatistics().getSum(vorticityIds[iVort][2]);
    }
    return vorticities;
}

/* ******** InternalVorticitySingleProbe3D *********************************** */

template <typename T, template <typename U> class Descriptor>
InternalVorticitySingleProbe3D<T, Descriptor>::InternalVorticitySingleProbe3D(
    std::vector<Array<T, 3> > const &positions_, std::vector<plint> const &vorticityIds_) :
    positions(positions_), vorticityIds(vorticityIds_)
{
    PLB_ASSERT(positions.size() == vorticityIds.size() / Descriptor<T>::d);
    PLB_ASSERT(vorticityIds.size() % Descriptor<T>::d == 0);
}

template <typename T, template <typename U> class Descriptor>
void InternalVorticitySingleProbe3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    Array<T, 3> vorticity;
    for (pluint iVort = 0; iVort < positions.size(); ++iVort) {
        vorticity.resetToZero();
        Array<T, 3> position(positions[iVort]);
        Dot3D referenceCellPos((plint)position[0], (plint)position[1], (plint)position[2]);
        referenceCellPos -= lattice.getLocation();
        if (contained(referenceCellPos, domain)) {
            vorticity = fdLattice::firstOrderBulkVorticity(
                lattice, referenceCellPos.x, referenceCellPos.y, referenceCellPos.z);
        }

        for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
            lattice.getInternalStatistics().gatherSum(vorticityIds[3 * iVort + iD], vorticity[iD]);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
InternalVorticitySingleProbe3D<T, Descriptor>
    *InternalVorticitySingleProbe3D<T, Descriptor>::clone() const
{
    return new InternalVorticitySingleProbe3D<T, Descriptor>(*this);
}

/* *************** Data Functionals for BlockLattice ***************** */

template <typename T, template <typename U> class Descriptor>
void CopyPopulationsFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &latticeFrom,
    BlockLattice3D<T, Descriptor> &latticeTo)
{
    Dot3D offset = computeRelativeDisplacement(latticeFrom, latticeTo);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint ofX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint ofY = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                latticeTo.get(ofX, ofY, iZ + offset.z).attributeValues(latticeFrom.get(iX, iY, iZ));
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
CopyPopulationsFunctional3D<T, Descriptor> *CopyPopulationsFunctional3D<T, Descriptor>::clone()
    const
{
    return new CopyPopulationsFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void CopyPopulationsFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT CopyPopulationsFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void CopyConvertPopulationsFunctional3D<T1, Descriptor1, T2, Descriptor2>::process(
    Box3D domain, BlockLattice3D<T1, Descriptor1> &latticeFrom,
    BlockLattice3D<T2, Descriptor2> &latticeTo)
{
    Dot3D offset = computeRelativeDisplacement(latticeFrom, latticeTo);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint ofX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint ofY = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                latticeTo.get(ofX, ofY, iZ + offset.z).attributeValues(latticeFrom.get(iX, iY, iZ));
            }
        }
    }
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
CopyConvertPopulationsFunctional3D<T1, Descriptor1, T2, Descriptor2>
    *CopyConvertPopulationsFunctional3D<T1, Descriptor1, T2, Descriptor2>::clone() const
{
    return new CopyConvertPopulationsFunctional3D<T1, Descriptor1, T2, Descriptor2>(*this);
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void CopyConvertPopulationsFunctional3D<T1, Descriptor1, T2, Descriptor2>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
BlockDomain::DomainT
    CopyConvertPopulationsFunctional3D<T1, Descriptor1, T2, Descriptor2>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void LatticeCopyAllFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &latticeFrom,
    BlockLattice3D<T, Descriptor> &latticeTo)
{
    std::vector<char> data;
    Dot3D offset = computeRelativeDisplacement(latticeFrom, latticeTo);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> const &fromCell = latticeFrom.get(iX, iY, iZ);
                Cell<T, Descriptor> &toCell =
                    latticeTo.get(iX + offset.x, iY + offset.y, iZ + offset.z);
                toCell.attributeValues(fromCell);

                data.clear();
                HierarchicSerializer serializer(data, fromCell.getDynamics().getId());
                fromCell.getDynamics().serialize(serializer);

                HierarchicUnserializer unSerializer(data, 0);
                toCell.getDynamics().unserialize(unSerializer);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
LatticeCopyAllFunctional3D<T, Descriptor> *LatticeCopyAllFunctional3D<T, Descriptor>::clone() const
{
    return new LatticeCopyAllFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void LatticeCopyAllFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::allVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT LatticeCopyAllFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void LatticeRegenerateFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &latticeFrom,
    BlockLattice3D<T, Descriptor> &latticeTo)
{
    Dot3D offset = computeRelativeDisplacement(latticeFrom, latticeTo);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                latticeTo.get(iX + offset.x, iY + offset.y, iZ + offset.z)
                    .attributeValues(latticeFrom.get(iX, iY, iZ));
                latticeTo.attributeDynamics(
                    iX + offset.x, iY + offset.y, iZ + offset.z,
                    latticeFrom.get(iX, iY, iZ).getDynamics().clone());
                latticeTo.get(iX + offset.x, iY + offset.y, iZ + offset.z)
                    .attributeValues(latticeFrom.get(iX, iY, iZ));
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
LatticeRegenerateFunctional3D<T, Descriptor> *LatticeRegenerateFunctional3D<T, Descriptor>::clone()
    const
{
    return new LatticeRegenerateFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void LatticeRegenerateFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    // Full dynamics object must be recreated, because this data processor
    //   re-attributes a new dynamics and acts on the bulk only.
    modified[1] = modif::dataStructure;
}

template <typename T, template <typename U> class Descriptor>
void BoxDensityFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    lattice.get(iX, iY, iZ).computeDensity();
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxDensityFunctional3D<T, Descriptor> *BoxDensityFunctional3D<T, Descriptor>::clone() const
{
    return new BoxDensityFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxDensityFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
void BoxRhoBarFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> const &cell = lattice.get(iX, iY, iZ);
                scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    cell.getDynamics().computeRhoBar(cell);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxRhoBarFunctional3D<T, Descriptor> *BoxRhoBarFunctional3D<T, Descriptor>::clone() const
{
    return new BoxRhoBarFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxRhoBarFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
void BoxRhoBarJfunctional3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_ASSERT(fields.size() == 3);
    BlockLattice3D<T, Descriptor> &lattice =
        *dynamic_cast<BlockLattice3D<T, Descriptor> *>(fields[0]);
    ScalarField3D<T> &rhoBarField = *dynamic_cast<ScalarField3D<T> *>(fields[1]);
    TensorField3D<T, 3> &jField = *dynamic_cast<TensorField3D<T, 3> *>(fields[2]);
    Dot3D offset1 = computeRelativeDisplacement(lattice, rhoBarField);
    Dot3D offset2 = computeRelativeDisplacement(lattice, jField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> const &cell = lattice.get(iX, iY, iZ);
                cell.getDynamics().computeRhoBarJ(
                    cell, rhoBarField.get(iX + offset1.x, iY + offset1.y, iZ + offset1.z),
                    jField.get(iX + offset2.x, iY + offset2.y, iZ + offset2.z));
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxRhoBarJfunctional3D<T, Descriptor> *BoxRhoBarJfunctional3D<T, Descriptor>::clone() const
{
    return new BoxRhoBarJfunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxRhoBarJfunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // lattice
    modified[1] = modif::staticVariables;  // rhoBar
    modified[2] = modif::staticVariables;  // j
}

template <typename T, template <typename U> class Descriptor>
void MaskedBoxRhoBarJfunctional3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_ASSERT(fields.size() == 4);
    BlockLattice3D<T, Descriptor> &lattice =
        *dynamic_cast<BlockLattice3D<T, Descriptor> *>(fields[0]);
    ScalarField3D<T> &rhoBarField = *dynamic_cast<ScalarField3D<T> *>(fields[1]);
    TensorField3D<T, 3> &jField = *dynamic_cast<TensorField3D<T, 3> *>(fields[2]);
    ScalarField3D<int> &maskField = *dynamic_cast<ScalarField3D<int> *>(fields[3]);
    Dot3D offset1 = computeRelativeDisplacement(lattice, rhoBarField);
    Dot3D offset2 = computeRelativeDisplacement(lattice, jField);
    Dot3D offset3 = computeRelativeDisplacement(lattice, maskField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (maskField.get(iX + offset3.x, iY + offset3.y, iZ + offset3.z) == flag) {
                    Cell<T, Descriptor> const &cell = lattice.get(iX, iY, iZ);
                    cell.getDynamics().computeRhoBarJ(
                        cell, rhoBarField.get(iX + offset1.x, iY + offset1.y, iZ + offset1.z),
                        jField.get(iX + offset2.x, iY + offset2.y, iZ + offset2.z));
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
MaskedBoxRhoBarJfunctional3D<T, Descriptor> *MaskedBoxRhoBarJfunctional3D<T, Descriptor>::clone()
    const
{
    return new MaskedBoxRhoBarJfunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void MaskedBoxRhoBarJfunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // lattice
    modified[1] = modif::staticVariables;  // rhoBar
    modified[2] = modif::staticVariables;  // j
    modified[3] = modif::nothing;          // mask
}

template <typename T, template <typename U> class Descriptor>
void BoxJfunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
    TensorField3D<T, Descriptor<T>::d> &jField)
{
    Dot3D offset2 = computeRelativeDisplacement(lattice, jField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> const &cell = lattice.get(iX, iY, iZ);
                momentTemplates<T, Descriptor>::get_j(
                    cell, jField.get(iX + offset2.x, iY + offset2.y, iZ + offset2.z));
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxJfunctional3D<T, Descriptor> *BoxJfunctional3D<T, Descriptor>::clone() const
{
    return new BoxJfunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxJfunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // lattice
    modified[1] = modif::staticVariables;  // j
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxJfunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxRhoBarJPiNeqfunctional3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_ASSERT(fields.size() == 4);
    BlockLattice3D<T, Descriptor> &lattice =
        *dynamic_cast<BlockLattice3D<T, Descriptor> *>(fields[0]);
    ScalarField3D<T> &rhoBarField = *dynamic_cast<ScalarField3D<T> *>(fields[1]);
    TensorField3D<T, 3> &jField = *dynamic_cast<TensorField3D<T, 3> *>(fields[2]);
    TensorField3D<T, 6> &piNeqField = *dynamic_cast<TensorField3D<T, 6> *>(fields[3]);
    Dot3D offset1 = computeRelativeDisplacement(lattice, rhoBarField);
    Dot3D offset2 = computeRelativeDisplacement(lattice, jField);
    Dot3D offset3 = computeRelativeDisplacement(lattice, piNeqField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> const &cell = lattice.get(iX, iY, iZ);
                cell.getDynamics().computeRhoBarJPiNeq(
                    cell, rhoBarField.get(iX + offset1.x, iY + offset1.y, iZ + offset1.z),
                    jField.get(iX + offset2.x, iY + offset2.y, iZ + offset2.z),
                    piNeqField.get(iX + offset3.x, iY + offset3.y, iZ + offset3.z));
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxRhoBarJPiNeqfunctional3D<T, Descriptor> *BoxRhoBarJPiNeqfunctional3D<T, Descriptor>::clone()
    const
{
    return new BoxRhoBarJPiNeqfunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxRhoBarJPiNeqfunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // lattice
    modified[1] = modif::staticVariables;  // rhoBar
    modified[2] = modif::staticVariables;  // j
    modified[3] = modif::staticVariables;  // piNeq
}

template <typename T, template <typename U> class Descriptor>
void PackedRhoBarJfunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, NTensorField3D<T> &rhoBarJField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, rhoBarJField);
    T rhoBar;
    Array<T, 3> j;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> const &cell = lattice.get(iX, iY, iZ);
                momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
                T *rhoBarJ = rhoBarJField.get(iX + offset.x, iY + offset.y, iZ + offset.z);
                *rhoBarJ = rhoBar;
                j.to_cArray(rhoBarJ + 1);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
PackedRhoBarJfunctional3D<T, Descriptor> *PackedRhoBarJfunctional3D<T, Descriptor>::clone() const
{
    return new PackedRhoBarJfunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void PackedRhoBarJfunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // lattice
    modified[1] = modif::staticVariables;  // rhoBarJ
}

template <typename T>
void DensityFromRhoBarJfunctional3D<T>::process(
    Box3D domain, ScalarField3D<T> &density, NTensorField3D<T> &rhoBarJField)
{
    Dot3D offset = computeRelativeDisplacement(density, rhoBarJField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T rhoBar = *rhoBarJField.get(iX + offset.x, iY + offset.y, iZ + offset.z);
                density.get(iX, iY, iZ) = rhoBar + (T)1;
            }
        }
    }
}

template <typename T>
DensityFromRhoBarJfunctional3D<T> *DensityFromRhoBarJfunctional3D<T>::clone() const
{
    return new DensityFromRhoBarJfunctional3D<T>(*this);
}

template <typename T>
void DensityFromRhoBarJfunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // density
    modified[1] = modif::nothing;          // rhoBarJ
}

template <typename T>
VelocityFromRhoBarJfunctional3D<T>::VelocityFromRhoBarJfunctional3D(bool velIsJ_) : velIsJ(velIsJ_)
{ }

template <typename T>
void VelocityFromRhoBarJfunctional3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_ASSERT(fields.size() == 2);
    TensorField3D<T, 3> *velocity = dynamic_cast<TensorField3D<T, 3> *>(fields[0]);
    NTensorField3D<T> *rhoBarJfield = dynamic_cast<NTensorField3D<T> *>(fields[1]);
    PLB_ASSERT(velocity);
    PLB_ASSERT(rhoBarJfield);
    Dot3D offset = computeRelativeDisplacement(*velocity, *rhoBarJfield);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, 3> nextVelocity;
                nextVelocity.from_cArray(
                    (*rhoBarJfield).get(iX + offset.x, iY + offset.y, iZ + offset.z) + 1);
                if (!velIsJ) {
                    nextVelocity /=
                        (*(*rhoBarJfield).get(iX + offset.x, iY + offset.y, iZ + offset.z)) + (T)1;
                }
                (*velocity).get(iX, iY, iZ) = nextVelocity;
            }
        }
    }
}

template <typename T>
VelocityFromRhoBarJfunctional3D<T> *VelocityFromRhoBarJfunctional3D<T>::clone() const
{
    return new VelocityFromRhoBarJfunctional3D<T>(*this);
}

template <typename T>
void VelocityFromRhoBarJfunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // velocity
    modified[1] = modif::nothing;          // rhoBarJ
}

template <typename T>
VelocityFromRhoBarAndJfunctional3D<T>::VelocityFromRhoBarAndJfunctional3D(bool velIsJ_) :
    velIsJ(velIsJ_)
{ }

template <typename T>
void VelocityFromRhoBarAndJfunctional3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_ASSERT(fields.size() == 3);
    TensorField3D<T, 3> *velocity = dynamic_cast<TensorField3D<T, 3> *>(fields[0]);
    ScalarField3D<T> *rhoBar = dynamic_cast<ScalarField3D<T> *>(fields[1]);
    TensorField3D<T, 3> *j = dynamic_cast<TensorField3D<T, 3> *>(fields[2]);
    PLB_ASSERT(velocity);
    PLB_ASSERT(rhoBar);
    PLB_ASSERT(j);

    Dot3D ofsR = computeRelativeDisplacement(*velocity, *rhoBar);
    Dot3D ofsJ = computeRelativeDisplacement(*velocity, *j);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, 3> nextVelocity = j->get(iX + ofsJ.x, iY + ofsJ.y, iZ + ofsJ.z);
                if (!velIsJ) {
                    nextVelocity /= (rhoBar->get(iX + ofsR.x, iY + ofsR.y, iZ + ofsR.z) + (T)1);
                }
                velocity->get(iX, iY, iZ) = nextVelocity;
            }
        }
    }
}

template <typename T>
VelocityFromRhoBarAndJfunctional3D<T> *VelocityFromRhoBarAndJfunctional3D<T>::clone() const
{
    return new VelocityFromRhoBarAndJfunctional3D<T>(*this);
}

template <typename T>
void VelocityFromRhoBarAndJfunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // velocity
    modified[1] = modif::nothing;          // rhoBar
    modified[2] = modif::nothing;          // j
}

template <typename T, template <typename U> class Descriptor>
void BoxKineticEnergyFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, Descriptor<T>::d> velocity;
                lattice.get(iX, iY, iZ).computeVelocity(velocity);
                scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    VectorTemplate<T, Descriptor>::normSqr(velocity) / (T)2;
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxKineticEnergyFunctional3D<T, Descriptor> *BoxKineticEnergyFunctional3D<T, Descriptor>::clone()
    const
{
    return new BoxKineticEnergyFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxKineticEnergyFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxKineticEnergyFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxVelocityNormFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, Descriptor<T>::d> velocity;
                lattice.get(iX, iY, iZ).computeVelocity(velocity);
                // The type cast converts the result of normSqr to type U in case T is of type
                // Complex<U>. Otherwise, the call to std::sqrt would fail, because std::sqrt is
                // overloaded, but not for Palabos' Complex type.
                scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) = std::sqrt(
                    (typename PlbTraits<T>::BaseType)VectorTemplate<T, Descriptor>::normSqr(
                        velocity));
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxVelocityNormFunctional3D<T, Descriptor> *BoxVelocityNormFunctional3D<T, Descriptor>::clone()
    const
{
    return new BoxVelocityNormFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxVelocityNormFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxVelocityNormFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxForcedVelocityNormFunctional3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 3);
    BlockLattice3D<T, Descriptor> *lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> *>(blocks[0]);
    PLB_ASSERT(lattice);
    TensorField3D<T, Descriptor<T>::d> *force =
        dynamic_cast<TensorField3D<T, Descriptor<T>::d> *>(blocks[1]);
    PLB_ASSERT(force);
    ScalarField3D<T> *velocityNorm = dynamic_cast<ScalarField3D<T> *>(blocks[2]);
    PLB_ASSERT(velocityNorm);

    Dot3D ofsF = computeRelativeDisplacement(*lattice, *force);
    Dot3D ofsV = computeRelativeDisplacement(*lattice, *velocityNorm);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, Descriptor<T>::d> velocity;
                lattice->get(iX, iY, iZ).computeVelocity(velocity);
                if (lattice->get(iX, iY, iZ).getDynamics().hasMoments()) {
                    Array<T, Descriptor<T>::d> const &f =
                        force->get(iX + ofsF.x, iY + ofsF.y, iZ + ofsF.z);
                    velocity[0] += (T)0.5 * f[0];
                    velocity[1] += (T)0.5 * f[1];
                    velocity[2] += (T)0.5 * f[2];
                }

                // The type cast converts the result of normSqr to type U in case T is of type
                // Complex<U>. Otherwise, the call to std::sqrt would fail, because std::sqrt is
                // overloaded, but not for Palabos' Complex type.
                velocityNorm->get(iX + ofsV.x, iY + ofsV.y, iZ + ofsV.z) = std::sqrt(
                    (typename PlbTraits<T>::BaseType)VectorTemplate<T, Descriptor>::normSqr(
                        velocity));
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxForcedVelocityNormFunctional3D<T, Descriptor>
    *BoxForcedVelocityNormFunctional3D<T, Descriptor>::clone() const
{
    return new BoxForcedVelocityNormFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxForcedVelocityNormFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxForcedVelocityNormFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
BoxConstForcedVelocityNormFunctional3D<T, Descriptor>::BoxConstForcedVelocityNormFunctional3D(
    Array<T, Descriptor<T>::d> force_) :
    force(force_)
{ }

template <typename T, template <typename U> class Descriptor>
void BoxConstForcedVelocityNormFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, Descriptor<T>::d> velocity;
                lattice.get(iX, iY, iZ).computeVelocity(velocity);
                if (lattice.get(iX, iY, iZ).getDynamics().hasMoments()) {
                    velocity[0] += (T)0.5 * force[0];
                    velocity[1] += (T)0.5 * force[1];
                    velocity[2] += (T)0.5 * force[2];
                }

                // The type cast converts the result of normSqr to type U in case T is of type
                // Complex<U>. Otherwise, the call to std::sqrt would fail, because std::sqrt is
                // overloaded, but not for Palabos' Complex type.
                scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) = std::sqrt(
                    (typename PlbTraits<T>::BaseType)VectorTemplate<T, Descriptor>::normSqr(
                        velocity));
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxConstForcedVelocityNormFunctional3D<T, Descriptor>
    *BoxConstForcedVelocityNormFunctional3D<T, Descriptor>::clone() const
{
    return new BoxConstForcedVelocityNormFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxConstForcedVelocityNormFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxConstForcedVelocityNormFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
BoxCustomForcedVelocityNormFunctional3D<
    T, Descriptor, ForceFunction>::BoxCustomForcedVelocityNormFunctional3D(ForceFunction f_) :
    f(f_)
{ }

template <typename T, template <typename U> class Descriptor, class ForceFunction>
void BoxCustomForcedVelocityNormFunctional3D<T, Descriptor, ForceFunction>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    Dot3D location = lattice.getLocation();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint x = iX + location.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint y = iY + location.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint z = iZ + location.z;
                Array<T, Descriptor<T>::d> velocity;
                lattice.get(iX, iY, iZ).computeVelocity(velocity);
                if (lattice.get(iX, iY, iZ).getDynamics().hasMoments()) {
                    Array<T, Descriptor<T>::d> force;
                    f(x, y, z, force);
                    velocity[0] += (T)0.5 * force[0];
                    velocity[1] += (T)0.5 * force[1];
                    velocity[2] += (T)0.5 * force[2];
                }

                // The type cast converts the result of normSqr to type U in case T is of type
                // Complex<U>. Otherwise, the call to std::sqrt would fail, because std::sqrt is
                // overloaded, but not for Palabos' Complex type.
                scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) = std::sqrt(
                    (typename PlbTraits<T>::BaseType)VectorTemplate<T, Descriptor>::normSqr(
                        velocity));
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
BoxCustomForcedVelocityNormFunctional3D<T, Descriptor, ForceFunction>
    *BoxCustomForcedVelocityNormFunctional3D<T, Descriptor, ForceFunction>::clone() const
{
    return new BoxCustomForcedVelocityNormFunctional3D<T, Descriptor, ForceFunction>(*this);
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
void BoxCustomForcedVelocityNormFunctional3D<T, Descriptor, ForceFunction>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
BlockDomain::DomainT
    BoxCustomForcedVelocityNormFunctional3D<T, Descriptor, ForceFunction>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
BoxVelocityComponentFunctional3D<T, Descriptor>::BoxVelocityComponentFunctional3D(int iComponent_) :
    iComponent(iComponent_)
{ }

template <typename T, template <typename U> class Descriptor>
void BoxVelocityComponentFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, Descriptor<T>::d> velocity;
                lattice.get(iX, iY, iZ).computeVelocity(velocity);
                scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) = velocity[iComponent];
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxVelocityComponentFunctional3D<T, Descriptor>
    *BoxVelocityComponentFunctional3D<T, Descriptor>::clone() const
{
    return new BoxVelocityComponentFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxVelocityComponentFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxVelocityComponentFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
BoxForcedVelocityComponentFunctional3D<T, Descriptor>::BoxForcedVelocityComponentFunctional3D(
    int iComponent_) :
    iComponent(iComponent_)
{ }

template <typename T, template <typename U> class Descriptor>
void BoxForcedVelocityComponentFunctional3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 3);
    BlockLattice3D<T, Descriptor> *lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> *>(blocks[0]);
    PLB_ASSERT(lattice);
    TensorField3D<T, Descriptor<T>::d> *force =
        dynamic_cast<TensorField3D<T, Descriptor<T>::d> *>(blocks[1]);
    PLB_ASSERT(force);
    ScalarField3D<T> *velocityComponent = dynamic_cast<ScalarField3D<T> *>(blocks[2]);
    PLB_ASSERT(velocityComponent);

    Dot3D ofsF = computeRelativeDisplacement(*lattice, *force);
    Dot3D ofsV = computeRelativeDisplacement(*lattice, *velocityComponent);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, Descriptor<T>::d> velocity;
                lattice->get(iX, iY, iZ).computeVelocity(velocity);
                if (lattice->get(iX, iY, iZ).getDynamics().hasMoments()) {
                    Array<T, Descriptor<T>::d> const &f =
                        force->get(iX + ofsF.x, iY + ofsF.y, iZ + ofsF.z);
                    velocityComponent->get(iX + ofsV.x, iY + ofsV.y, iZ + ofsV.z) =
                        velocity[iComponent] + (T)0.5 * f[iComponent];
                } else {
                    velocityComponent->get(iX + ofsV.x, iY + ofsV.y, iZ + ofsV.z) =
                        velocity[iComponent];
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxForcedVelocityComponentFunctional3D<T, Descriptor>
    *BoxForcedVelocityComponentFunctional3D<T, Descriptor>::clone() const
{
    return new BoxForcedVelocityComponentFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxForcedVelocityComponentFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxForcedVelocityComponentFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
BoxConstForcedVelocityComponentFunctional3D<T, Descriptor>::
    BoxConstForcedVelocityComponentFunctional3D(
        Array<T, Descriptor<T>::d> force_, int iComponent_) :
    force(force_), iComponent(iComponent_)
{ }

template <typename T, template <typename U> class Descriptor>
void BoxConstForcedVelocityComponentFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, Descriptor<T>::d> velocity;
                lattice.get(iX, iY, iZ).computeVelocity(velocity);
                if (lattice.get(iX, iY, iZ).getDynamics().hasMoments()) {
                    scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                        velocity[iComponent] + (T)0.5 * force[iComponent];
                } else {
                    scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                        velocity[iComponent];
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxConstForcedVelocityComponentFunctional3D<T, Descriptor>
    *BoxConstForcedVelocityComponentFunctional3D<T, Descriptor>::clone() const
{
    return new BoxConstForcedVelocityComponentFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxConstForcedVelocityComponentFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxConstForcedVelocityComponentFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
BoxCustomForcedVelocityComponentFunctional3D<T, Descriptor, ForceFunction>::
    BoxCustomForcedVelocityComponentFunctional3D(ForceFunction f_, int iComponent_) :
    f(f_), iComponent(iComponent_)
{ }

template <typename T, template <typename U> class Descriptor, class ForceFunction>
void BoxCustomForcedVelocityComponentFunctional3D<T, Descriptor, ForceFunction>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    Dot3D location = lattice.getLocation();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint x = iX + location.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint y = iY + location.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint z = iZ + location.z;
                Array<T, Descriptor<T>::d> velocity;
                lattice.get(iX, iY, iZ).computeVelocity(velocity);
                if (lattice.get(iX, iY, iZ).getDynamics().hasMoments()) {
                    Array<T, Descriptor<T>::d> force;
                    f(x, y, z, force);
                    scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                        velocity[iComponent] + (T)0.5 * force[iComponent];
                } else {
                    scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                        velocity[iComponent];
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
BoxCustomForcedVelocityComponentFunctional3D<T, Descriptor, ForceFunction>
    *BoxCustomForcedVelocityComponentFunctional3D<T, Descriptor, ForceFunction>::clone() const
{
    return new BoxCustomForcedVelocityComponentFunctional3D<T, Descriptor, ForceFunction>(*this);
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
void BoxCustomForcedVelocityComponentFunctional3D<
    T, Descriptor, ForceFunction>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
BlockDomain::DomainT
    BoxCustomForcedVelocityComponentFunctional3D<T, Descriptor, ForceFunction>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxVelocityFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
    TensorField3D<T, Descriptor<T>::d> &tensorField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                lattice.get(iX, iY, iZ)
                    .computeVelocity(tensorField.get(iX + offset.x, iY + offset.y, iZ + offset.z));
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxVelocityFunctional3D<T, Descriptor> *BoxVelocityFunctional3D<T, Descriptor>::clone() const
{
    return new BoxVelocityFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxVelocityFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxVelocityFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxForcedVelocityFunctional3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 3);
    BlockLattice3D<T, Descriptor> *lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> *>(blocks[0]);
    PLB_ASSERT(lattice);
    TensorField3D<T, Descriptor<T>::d> *force =
        dynamic_cast<TensorField3D<T, Descriptor<T>::d> *>(blocks[1]);
    PLB_ASSERT(force);
    TensorField3D<T, Descriptor<T>::d> *velocity =
        dynamic_cast<TensorField3D<T, Descriptor<T>::d> *>(blocks[2]);
    PLB_ASSERT(velocity);

    Dot3D ofsF = computeRelativeDisplacement(*lattice, *force);
    Dot3D ofsV = computeRelativeDisplacement(*lattice, *velocity);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, Descriptor<T>::d> u;
                lattice->get(iX, iY, iZ).computeVelocity(u);
                Array<T, Descriptor<T>::d> &v =
                    velocity->get(iX + ofsV.x, iY + ofsV.y, iZ + ofsV.z);
                if (lattice->get(iX, iY, iZ).getDynamics().hasMoments()) {
                    Array<T, Descriptor<T>::d> const &f =
                        force->get(iX + ofsF.x, iY + ofsF.y, iZ + ofsF.z);
                    v[0] = u[0] + (T)0.5 * f[0];
                    v[1] = u[1] + (T)0.5 * f[1];
                    v[2] = u[2] + (T)0.5 * f[2];
                } else {
                    v[0] = u[0];
                    v[1] = u[1];
                    v[2] = u[2];
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxForcedVelocityFunctional3D<T, Descriptor> *BoxForcedVelocityFunctional3D<T, Descriptor>::clone()
    const
{
    return new BoxForcedVelocityFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxForcedVelocityFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxForcedVelocityFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
BoxConstForcedVelocityFunctional3D<T, Descriptor>::BoxConstForcedVelocityFunctional3D(
    Array<T, Descriptor<T>::d> force_) :
    force(force_)
{ }

template <typename T, template <typename U> class Descriptor>
void BoxConstForcedVelocityFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
    TensorField3D<T, Descriptor<T>::d> &tensorField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, Descriptor<T>::d> u;
                lattice.get(iX, iY, iZ).computeVelocity(u);
                Array<T, Descriptor<T>::d> &v =
                    tensorField.get(iX + offset.x, iY + offset.y, iZ + offset.z);
                if (lattice.get(iX, iY, iZ).getDynamics().hasMoments()) {
                    v[0] = u[0] + (T)0.5 * force[0];
                    v[1] = u[1] + (T)0.5 * force[1];
                    v[2] = u[2] + (T)0.5 * force[2];
                } else {
                    v[0] = u[0];
                    v[1] = u[1];
                    v[2] = u[2];
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxConstForcedVelocityFunctional3D<T, Descriptor>
    *BoxConstForcedVelocityFunctional3D<T, Descriptor>::clone() const
{
    return new BoxConstForcedVelocityFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxConstForcedVelocityFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxConstForcedVelocityFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
BoxCustomForcedVelocityFunctional3D<
    T, Descriptor, ForceFunction>::BoxCustomForcedVelocityFunctional3D(ForceFunction f_) :
    f(f_)
{ }

template <typename T, template <typename U> class Descriptor, class ForceFunction>
void BoxCustomForcedVelocityFunctional3D<T, Descriptor, ForceFunction>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
    TensorField3D<T, Descriptor<T>::d> &tensorField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, tensorField);
    Dot3D location = lattice.getLocation();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint x = iX + location.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint y = iY + location.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint z = iZ + location.z;
                Array<T, Descriptor<T>::d> u;
                lattice.get(iX, iY, iZ).computeVelocity(u);
                Array<T, Descriptor<T>::d> &v =
                    tensorField.get(iX + offset.x, iY + offset.y, iZ + offset.z);
                if (lattice.get(iX, iY, iZ).getDynamics().hasMoments()) {
                    Array<T, Descriptor<T>::d> force;
                    f(x, y, z, force);
                    v[0] = u[0] + (T)0.5 * force[0];
                    v[1] = u[1] + (T)0.5 * force[1];
                    v[2] = u[2] + (T)0.5 * force[2];
                } else {
                    v[0] = u[0];
                    v[1] = u[1];
                    v[2] = u[2];
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
BoxCustomForcedVelocityFunctional3D<T, Descriptor, ForceFunction>
    *BoxCustomForcedVelocityFunctional3D<T, Descriptor, ForceFunction>::clone() const
{
    return new BoxCustomForcedVelocityFunctional3D<T, Descriptor, ForceFunction>(*this);
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
void BoxCustomForcedVelocityFunctional3D<T, Descriptor, ForceFunction>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor, class ForceFunction>
BlockDomain::DomainT BoxCustomForcedVelocityFunctional3D<T, Descriptor, ForceFunction>::appliesTo()
    const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxTemperatureFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    lattice.get(iX, iY, iZ).computeTemperature();
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxTemperatureFunctional3D<T, Descriptor> *BoxTemperatureFunctional3D<T, Descriptor>::clone() const
{
    return new BoxTemperatureFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxTemperatureFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxTemperatureFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxPiNeqFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
    TensorField3D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq)
{
    Dot3D offset = computeRelativeDisplacement(lattice, PiNeq);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                lattice.get(iX, iY, iZ)
                    .computePiNeq(PiNeq.get(iX + offset.x, iY + offset.y, iZ + offset.z));
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxPiNeqFunctional3D<T, Descriptor> *BoxPiNeqFunctional3D<T, Descriptor>::clone() const
{
    return new BoxPiNeqFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxPiNeqFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
void BoxShearStressFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
    TensorField3D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq)
{
    Dot3D offset = computeRelativeDisplacement(lattice, PiNeq);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                lattice.get(iX, iY, iZ)
                    .computeShearStress(PiNeq.get(iX + offset.x, iY + offset.y, iZ + offset.z));
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxShearStressFunctional3D<T, Descriptor> *BoxShearStressFunctional3D<T, Descriptor>::clone() const
{
    return new BoxShearStressFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxShearStressFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxShearStressFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
BoxStressFunctional3D<T, Descriptor>::BoxStressFunctional3D(T rho0_, bool isCompressible_) :
    rho0(rho0_), isCompressible(isCompressible_)
{ }

template <typename T, template <typename U> class Descriptor>
void BoxStressFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
    TensorField3D<T, SymmetricTensor<T, Descriptor>::n> &stress)
{
    typedef SymmetricTensorImpl<T, 3> tensor;
    Dot3D offset = computeRelativeDisplacement(lattice, stress);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> const &cell = lattice.get(iX, iY, iZ);
                Array<T, tensor::n> &s = stress.get(iX + offset.x, iY + offset.y, iZ + offset.z);
                cell.computeShearStress(s);
                T rho = cell.computeDensity();
                if (isCompressible) {
                    s /= rho;
                }
                T pressure = (rho - rho0) * Descriptor<T>::cs2;
                s[tensor::xx] -= pressure;
                s[tensor::yy] -= pressure;
                s[tensor::zz] -= pressure;
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxStressFunctional3D<T, Descriptor> *BoxStressFunctional3D<T, Descriptor>::clone() const
{
    return new BoxStressFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxStressFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxStressFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxStrainRateFromStressFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
    TensorField3D<T, SymmetricTensor<T, Descriptor>::n> &S)
{
    Dot3D offset = computeRelativeDisplacement(lattice, S);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> const &cell = lattice.get(iX, iY, iZ);
                Array<T, SymmetricTensor<T, Descriptor>::n> &element =
                    S.get(iX + offset.x, iY + offset.y, iZ + offset.z);
                cell.computePiNeq(element);
                T omega = cell.getDynamics().getOmega();
                if (cell.getDynamics().hasMoments()) {
                    T dynamicOmega =
                        cell.getDynamics().getDynamicParameter(dynamicParams::dynamicOmega, cell);
                    if (!util::isZero(dynamicOmega)) {
                        omega = dynamicOmega;
                    }
                }
                T rhoBar = cell.getDynamics().computeRhoBar(cell);
                T prefactor = -omega * Descriptor<T>::invCs2 * Descriptor<T>::invRho(rhoBar) / (T)2;
                for (int iTensor = 0; iTensor < SymmetricTensor<T, Descriptor>::n; ++iTensor) {
                    element[iTensor] *= prefactor;
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxStrainRateFromStressFunctional3D<T, Descriptor>
    *BoxStrainRateFromStressFunctional3D<T, Descriptor>::clone() const
{
    return new BoxStrainRateFromStressFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxStrainRateFromStressFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxStrainRateFromStressFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxShearRateFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &shearRate)
{
    Dot3D offset = computeRelativeDisplacement(lattice, shearRate);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> const &cell = lattice.get(iX, iY, iZ);

                T rhoBar = (T)0;
                Array<T, 3> j((T)0, (T)0, (T)0);
                Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
                momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
                T normPiNeq = std::sqrt(SymmetricTensor<T, Descriptor>::tensorNormSqr(PiNeq));
                T rho = Descriptor<T>::fullRho(rhoBar);

                T dynamicOmega =
                    cell.getDynamics().getDynamicParameter(dynamicParams::dynamicOmega, cell);
                T tau = 0.0;
                if (!util::isZero(dynamicOmega)) {
                    tau = (T)1 / dynamicOmega;
                } else {
                    tau = (T)1 / cell.getDynamics().getOmega();
                }

                shearRate.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    normPiNeq / (tau * Descriptor<T>::cs2 * rho);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxShearRateFunctional3D<T, Descriptor> *BoxShearRateFunctional3D<T, Descriptor>::clone() const
{
    return new BoxShearRateFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxShearRateFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
void BoxNTensorShearRateFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, NTensorField3D<T> &shearRate)
{
    PLB_PRECONDITION(shearRate.getNdim() == 1);
    Dot3D offset = computeRelativeDisplacement(lattice, shearRate);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> const &cell = lattice.get(iX, iY, iZ);

                T rhoBar = (T)0;
                Array<T, 3> j((T)0, (T)0, (T)0);
                Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
                momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
                T normPiNeq = std::sqrt(SymmetricTensor<T, Descriptor>::tensorNormSqr(PiNeq));
                T rho = Descriptor<T>::fullRho(rhoBar);

                T dynamicOmega =
                    cell.getDynamics().getDynamicParameter(dynamicParams::dynamicOmega, cell);
                T tau = 0.0;
                if (!util::isZero(dynamicOmega)) {
                    tau = (T)1 / dynamicOmega;
                } else {
                    tau = (T)1 / cell.getDynamics().getOmega();
                }

                *shearRate.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    normPiNeq / (tau * Descriptor<T>::cs2 * rho);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxNTensorShearRateFunctional3D<T, Descriptor>
    *BoxNTensorShearRateFunctional3D<T, Descriptor>::clone() const
{
    return new BoxNTensorShearRateFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxNTensorShearRateFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
void BoxQcriterionFunctional3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_ASSERT(fields.size() == 3);
    TensorField3D<T, 3> &vorticity = *dynamic_cast<TensorField3D<T, 3> *>(fields[0]);
    TensorField3D<T, 6> &strain = *dynamic_cast<TensorField3D<T, 6> *>(fields[1]);
    ScalarField3D<T> &qCriterion = *dynamic_cast<ScalarField3D<T> *>(fields[2]);
    Dot3D offset1 = computeRelativeDisplacement(vorticity, strain);
    Dot3D offset2 = computeRelativeDisplacement(vorticity, qCriterion);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX1 = iX + offset1.x;
        plint oX2 = iX + offset2.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY1 = iY + offset1.y;
            plint oY2 = iY + offset2.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint oZ1 = iZ + offset1.z;
                plint oZ2 = iZ + offset2.z;

                T vortNorm = VectorTemplateImpl<T, 3>::normSqr(vorticity.get(iX, iY, iZ));
                T normStrain = SymmetricTensorImpl<T, 3>::tensorNormSqr(strain.get(oX1, oY1, oZ1));

                qCriterion.get(oX2, oY2, oZ2) = (vortNorm - (T)2 * normStrain) / (T)4;
            }
        }
    }
}

template <typename T>
BoxQcriterionFunctional3D<T> *BoxQcriterionFunctional3D<T>::clone() const
{
    return new BoxQcriterionFunctional3D<T>(*this);
}

template <typename T>
void BoxQcriterionFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT BoxQcriterionFunctional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T>
void BoxComputeInstantaneousReynoldsStressFunctional3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_ASSERT(fields.size() == 3);
    TensorField3D<T, 3> &vel = *dynamic_cast<TensorField3D<T, 3> *>(fields[0]);
    TensorField3D<T, 3> &avgVel = *dynamic_cast<TensorField3D<T, 3> *>(fields[1]);
    TensorField3D<T, 6> &reynoldsStress = *dynamic_cast<TensorField3D<T, 6> *>(fields[2]);

    Dot3D offset1 = computeRelativeDisplacement(vel, avgVel);
    Dot3D offset2 = computeRelativeDisplacement(vel, reynoldsStress);

    typedef SymmetricTensorImpl<T, 3> S;

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX1 = iX + offset1.x;
        plint oX2 = iX + offset2.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY1 = iY + offset1.y;
            plint oY2 = iY + offset2.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint oZ1 = iZ + offset1.z;
                plint oZ2 = iZ + offset2.z;

                Array<T, 3> diffVel = vel.get(iX, iY, iZ) - avgVel.get(oX1, oY1, oZ1);

                reynoldsStress.get(oX2, oY2, oZ2)[S::xx] = diffVel[0] * diffVel[0];
                reynoldsStress.get(oX2, oY2, oZ2)[S::xy] = diffVel[0] * diffVel[1];
                reynoldsStress.get(oX2, oY2, oZ2)[S::xz] = diffVel[0] * diffVel[2];
                reynoldsStress.get(oX2, oY2, oZ2)[S::yy] = diffVel[1] * diffVel[1];
                reynoldsStress.get(oX2, oY2, oZ2)[S::yz] = diffVel[1] * diffVel[2];
                reynoldsStress.get(oX2, oY2, oZ2)[S::zz] = diffVel[2] * diffVel[2];
            }
        }
    }
}

template <typename T>
BoxComputeInstantaneousReynoldsStressFunctional3D<T>
    *BoxComputeInstantaneousReynoldsStressFunctional3D<T>::clone() const
{
    return new BoxComputeInstantaneousReynoldsStressFunctional3D<T>(*this);
}

template <typename T>
void BoxComputeInstantaneousReynoldsStressFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT BoxComputeInstantaneousReynoldsStressFunctional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

#ifndef PLB_BGP
#ifdef PLB_USE_EIGEN

template <typename T>
void BoxLambda2Functional3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);

    TensorField3D<T, 3> *vorticity = dynamic_cast<TensorField3D<T, 3> *>(fields[0]);
    TensorField3D<T, 6> *strain = dynamic_cast<TensorField3D<T, 6> *>(fields[1]);
    ScalarField3D<T> *lambda2 = dynamic_cast<ScalarField3D<T> *>(fields[2]);
    PLB_ASSERT(vorticity);
    PLB_ASSERT(strain);
    PLB_ASSERT(lambda2);

    Dot3D offset1 = computeRelativeDisplacement(*vorticity, *strain);
    Dot3D offset2 = computeRelativeDisplacement(*vorticity, *lambda2);

    static plint xx = SymmetricTensorImpl<T, 3>::xx;
    static plint xy = SymmetricTensorImpl<T, 3>::xy;
    static plint xz = SymmetricTensorImpl<T, 3>::xz;
    static plint yy = SymmetricTensorImpl<T, 3>::yy;
    static plint yz = SymmetricTensorImpl<T, 3>::yz;
    static plint zz = SymmetricTensorImpl<T, 3>::zz;
    static plint yx = xy;
    static plint zx = xz;
    static plint zy = yz;

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX1 = iX + offset1.x;
        plint oX2 = iX + offset2.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY1 = iY + offset1.y;
            plint oY2 = iY + offset2.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint oZ1 = iZ + offset1.z;
                plint oZ2 = iZ + offset2.z;

                Array<T, 3> const &v = vorticity->get(iX, iY, iZ);
                Array<T, 6> const &s = strain->get(oX1, oY1, oZ1);

                Array<Array<T, 3>, 3> S;  // Strain-rate tensor (symmetric).
                S[0][0] = s[xx];
                S[0][1] = s[xy];
                S[0][2] = s[xz];

                S[1][0] = s[yx];
                S[1][1] = s[yy];
                S[1][2] = s[yz];

                S[2][0] = s[zx];
                S[2][1] = s[zy];
                S[2][2] = s[zz];

                Array<Array<T, 3>, 3> V;  // Vorticity tensor (anti-symmetric).
                V[0][0] = (T)0;
                V[0][1] = -v[2];
                V[0][2] = v[1];

                V[1][0] = v[2];
                V[1][1] = (T)0;
                V[1][2] = -v[0];

                V[2][0] = -v[1];
                V[2][1] = v[0];
                V[2][2] = (T)0;

                Array<Array<T, 3>, 3> S2V2;  // S^2 + V^2 matrix (symmetric).
                for (plint i = 0; i < 3; i++) {
                    for (plint j = 0; j < 3; j++) {
                        T S2 = (T)0;
                        T V2 = (T)0;
                        for (plint k = 0; k < 3; k++) {
                            S2 += S[i][k] * S[k][j];
                            V2 += V[i][k] * V[k][j];
                        }
                        S2V2[i][j] = S2 + V2;
                    }
                }

                /*
                Array<Array<T,3>,3> x;  // Eigenvectors of S2V2.
                Array<T,3> d;           // Eigenvalues of S2V2.
                eigenDecomposition(S2V2, x, d);
                std::vector<T> lambda(3);
                lambda[0] = d[0];
                lambda[1] = d[1];
                lambda[2] = d[2];
                std::sort(lambda.begin(), lambda.end());

                lambda2->get(oX2, oY2, oZ2) = lambda[1];
                */

                Eigen::Matrix<T, 3, 3> A;
                for (plint i = 0; i < 3; i++) {
                    for (plint j = 0; j < 3; j++) {
                        A(i, j) = S2V2[i][j];
                    }
                }

                bool computeEigenvectors = false;
                Eigen::EigenSolver<Eigen::Matrix<T, 3, 3> > es(A, computeEigenvectors);
                std::vector<T> lambda(3);
                lambda[0] = std::real(es.eigenvalues()[0]);
                lambda[1] = std::real(es.eigenvalues()[1]);
                lambda[2] = std::real(es.eigenvalues()[2]);
                std::sort(lambda.begin(), lambda.end());

                lambda2->get(oX2, oY2, oZ2) = lambda[1];
            }
        }
    }
}

template <typename T>
BoxLambda2Functional3D<T> *BoxLambda2Functional3D<T>::clone() const
{
    return new BoxLambda2Functional3D<T>(*this);
}

template <typename T>
void BoxLambda2Functional3D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Vorticity.
    modified[1] = modif::nothing;          // Strain-rate.
    modified[2] = modif::staticVariables;  // lambda2.
}

template <typename T>
BlockDomain::DomainT BoxLambda2Functional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

#endif
#endif

template <typename T, template <typename U> class Descriptor>
BoxPopulationFunctional3D<T, Descriptor>::BoxPopulationFunctional3D(plint iComponent_) :
    iComponent(iComponent_)
{ }

template <typename T, template <typename U> class Descriptor>
void BoxPopulationFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    lattice.get(iX, iY, iZ)[iComponent];
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxPopulationFunctional3D<T, Descriptor> *BoxPopulationFunctional3D<T, Descriptor>::clone() const
{
    return new BoxPopulationFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxPopulationFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxPopulationFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
BoxEquilibriumFunctional3D<T, Descriptor>::BoxEquilibriumFunctional3D(plint iComponent_) :
    iComponent(iComponent_)
{ }

template <typename T, template <typename U> class Descriptor>
void BoxEquilibriumFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &equilibrium)
{
    Dot3D offset = computeRelativeDisplacement(lattice, equilibrium);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T rhoBar;
                Array<T, Descriptor<T>::d> j;
                Cell<T, Descriptor> const &cell = lattice.get(iX, iY, iZ);
                cell.getDynamics().computeRhoBarJ(cell, rhoBar, j);
                T jSqr = normSqr(j);
                equilibrium.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    cell.computeEquilibrium(iComponent, rhoBar, j, jSqr);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxEquilibriumFunctional3D<T, Descriptor> *BoxEquilibriumFunctional3D<T, Descriptor>::clone() const
{
    return new BoxEquilibriumFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxEquilibriumFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxEquilibriumFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
BoxAllPopulationsFunctional3D<T, Descriptor>::BoxAllPopulationsFunctional3D()
{ }

template <typename T, template <typename U> class Descriptor>
void BoxAllPopulationsFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
    TensorField3D<T, Descriptor<T>::q> &tensorField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                    tensorField.get(oX, oY, iZ + offset.z)[iPop] = lattice.get(iX, iY, iZ)[iPop];
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxAllPopulationsFunctional3D<T, Descriptor> *BoxAllPopulationsFunctional3D<T, Descriptor>::clone()
    const
{
    return new BoxAllPopulationsFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxAllPopulationsFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxAllPopulationsFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxAllEquilibriumFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
    TensorField3D<T, Descriptor<T>::q> &equilibrium)
{
    Dot3D offset = computeRelativeDisplacement(lattice, equilibrium);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T rhoBar;
                Array<T, Descriptor<T>::d> j;
                Cell<T, Descriptor> const &cell = lattice.get(iX, iY, iZ);
                cell.getDynamics().computeRhoBarJ(cell, rhoBar, j);
                T jSqr = normSqr(j);
                for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                    equilibrium.get(iX + offset.x, iY + offset.y, iZ + offset.z)[iPop] =
                        cell.computeEquilibrium(iPop, rhoBar, j, jSqr);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxAllEquilibriumFunctional3D<T, Descriptor> *BoxAllEquilibriumFunctional3D<T, Descriptor>::clone()
    const
{
    return new BoxAllEquilibriumFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxAllEquilibriumFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxAllEquilibriumFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

// ========================== BoxAllNonEquilibriumFunctional3D ====================

template <typename T, template <typename U> class Descriptor>
void BoxAllNonEquilibriumFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
    TensorField3D<T, Descriptor<T>::q> &nonEquilibrium)
{
    Dot3D offset = computeRelativeDisplacement(lattice, nonEquilibrium);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T rhoBar;
                Array<T, Descriptor<T>::d> j;
                Cell<T, Descriptor> const &cell = lattice.get(iX, iY, iZ);
                cell.getDynamics().computeRhoBarJ(cell, rhoBar, j);
                T jSqr = normSqr(j);
                for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                    nonEquilibrium.get(iX + offset.x, iY + offset.y, iZ + offset.z)[iPop] =
                        cell[iPop] - cell.computeEquilibrium(iPop, rhoBar, j, jSqr);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxAllNonEquilibriumFunctional3D<T, Descriptor>
    *BoxAllNonEquilibriumFunctional3D<T, Descriptor>::clone() const
{
    return new BoxAllNonEquilibriumFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxAllNonEquilibriumFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxAllNonEquilibriumFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
BoxAllPopulationsToLatticeFunctional3D<T, Descriptor>::BoxAllPopulationsToLatticeFunctional3D()
{ }

template <typename T, template <typename U> class Descriptor>
void BoxAllPopulationsToLatticeFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
    TensorField3D<T, Descriptor<T>::q> &tensorField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                    lattice.get(iX, iY, iZ)[iPop] =
                        tensorField.get(iX + offset.x, iY + offset.y, iZ + offset.z)[iPop];
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxAllPopulationsToLatticeFunctional3D<T, Descriptor>
    *BoxAllPopulationsToLatticeFunctional3D<T, Descriptor>::clone() const
{
    return new BoxAllPopulationsToLatticeFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxAllPopulationsToLatticeFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxAllPopulationsToLatticeFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxOmegaFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> const &cell = lattice.get(iX, iY, iZ);
                scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    cell.getDynamics().getOmega();
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxOmegaFunctional3D<T, Descriptor> *BoxOmegaFunctional3D<T, Descriptor>::clone() const
{
    return new BoxOmegaFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxOmegaFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxOmegaFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxKinematicViscosityFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T omega = lattice.get(iX, iY, iZ).getDynamics().getOmega();
                if (!util::isZero(omega)) {
                    scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                        Descriptor<T>::cs2 * ((T)1 / omega - (T)0.5);
                } else {
                    scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) = (T)0;
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxKinematicViscosityFunctional3D<T, Descriptor>
    *BoxKinematicViscosityFunctional3D<T, Descriptor>::clone() const
{
    return new BoxKinematicViscosityFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxKinematicViscosityFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxKinematicViscosityFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxNTensorKinematicViscosityFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, NTensorField3D<T> &nu)
{
    PLB_PRECONDITION(nu.getNdim() == 1);
    Dot3D offset = computeRelativeDisplacement(lattice, nu);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T omega = lattice.get(iX, iY, iZ).getDynamics().getOmega();
                if (!util::isZero(omega)) {
                    *nu.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                        Descriptor<T>::cs2 * ((T)1 / omega - (T)0.5);
                } else {
                    *nu.get(iX + offset.x, iY + offset.y, iZ + offset.z) = (T)0;
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxNTensorKinematicViscosityFunctional3D<T, Descriptor>
    *BoxNTensorKinematicViscosityFunctional3D<T, Descriptor>::clone() const
{
    return new BoxNTensorKinematicViscosityFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxNTensorKinematicViscosityFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxNTensorKinematicViscosityFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxKinematicEddyViscosityFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
                if (cell.getDynamics().hasMoments()) {
                    T dynamicOmega =
                        cell.getDynamics().getDynamicParameter(dynamicParams::dynamicOmega, cell);
                    T omega = cell.getDynamics().getOmega();
                    if (!util::isZero(dynamicOmega) && !util::isZero(omega)) {
                        scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                            Descriptor<T>::cs2 * ((T)1 / dynamicOmega - (T)1 / omega);
                    } else {
                        scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) = (T)0;
                    }
                } else {
                    scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) = (T)0;
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxKinematicEddyViscosityFunctional3D<T, Descriptor>
    *BoxKinematicEddyViscosityFunctional3D<T, Descriptor>::clone() const
{
    return new BoxKinematicEddyViscosityFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxKinematicEddyViscosityFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxKinematicEddyViscosityFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxNTensorKinematicEddyViscosityFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, NTensorField3D<T> &nu)
{
    PLB_PRECONDITION(nu.getNdim() == 1);
    Dot3D offset = computeRelativeDisplacement(lattice, nu);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
                if (cell.getDynamics().hasMoments()) {
                    T dynamicOmega =
                        cell.getDynamics().getDynamicParameter(dynamicParams::dynamicOmega, cell);
                    T omega = cell.getDynamics().getOmega();
                    if (!util::isZero(dynamicOmega) && !util::isZero(omega)) {
                        *nu.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                            Descriptor<T>::cs2 * ((T)1 / dynamicOmega - (T)1 / omega);
                    } else {
                        *nu.get(iX + offset.x, iY + offset.y, iZ + offset.z) = (T)0;
                    }
                } else {
                    *nu.get(iX + offset.x, iY + offset.y, iZ + offset.z) = (T)0;
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxNTensorKinematicEddyViscosityFunctional3D<T, Descriptor>
    *BoxNTensorKinematicEddyViscosityFunctional3D<T, Descriptor>::clone() const
{
    return new BoxNTensorKinematicEddyViscosityFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxNTensorKinematicEddyViscosityFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxNTensorKinematicEddyViscosityFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxExternalForceFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
    TensorField3D<T, Descriptor<T>::d> &tensorField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T *force = lattice.get(iX, iY, iZ)
                               .getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
                tensorField.get(oX, oY, iZ + offset.z).from_cArray(force);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxExternalForceFunctional3D<T, Descriptor> *BoxExternalForceFunctional3D<T, Descriptor>::clone()
    const
{
    return new BoxExternalForceFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxExternalForceFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BoxExternalScalarFunctional3D<T, Descriptor>::BoxExternalScalarFunctional3D(int whichScalar_) :
    whichScalar(whichScalar_)
{
    PLB_ASSERT(whichScalar < Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor>
void BoxExternalScalarFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                scalarField.get(oX, oY, iZ + offset.z) =
                    *lattice.get(iX, iY, iZ).getExternal(whichScalar);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxExternalScalarFunctional3D<T, Descriptor> *BoxExternalScalarFunctional3D<T, Descriptor>::clone()
    const
{
    return new BoxExternalScalarFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxExternalScalarFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BoxExternalVectorFunctional3D<T, Descriptor>::BoxExternalVectorFunctional3D(int vectorBeginsAt_) :
    vectorBeginsAt(vectorBeginsAt_)
{
    PLB_ASSERT(vectorBeginsAt + Descriptor<T>::d <= Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor>
void BoxExternalVectorFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
    TensorField3D<T, Descriptor<T>::d> &tensorField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T *tensor = lattice.get(iX, iY, iZ).getExternal(vectorBeginsAt);
                tensorField.get(oX, oY, iZ + offset.z).from_cArray(tensor);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxExternalVectorFunctional3D<T, Descriptor> *BoxExternalVectorFunctional3D<T, Descriptor>::clone()
    const
{
    return new BoxExternalVectorFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxExternalVectorFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxExternalVectorFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
BoxDynamicParameterFunctional3D<T, Descriptor>::BoxDynamicParameterFunctional3D(
    plint whichParameter_) :
    whichParameter(whichParameter_)
{ }

template <typename T, template <typename U> class Descriptor>
void BoxDynamicParameterFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
                scalarField.get(oX, oY, iZ + offset.z) =
                    cell.getDynamics().getDynamicParameter(whichParameter, cell);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxDynamicParameterFunctional3D<T, Descriptor>
    *BoxDynamicParameterFunctional3D<T, Descriptor>::clone() const
{
    return new BoxDynamicParameterFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxDynamicParameterFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
void BoxDynamicViscosityFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
                T omega = cell.getDynamics().getDynamicParameter(dynamicParams::dynamicOmega, cell);
                T nu = Descriptor<T>::cs2 * ((T)1. / omega - (T)0.5);
                scalarField.get(oX, oY, iZ + offset.z) = nu;
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxDynamicViscosityFunctional3D<T, Descriptor>
    *BoxDynamicViscosityFunctional3D<T, Descriptor>::clone() const
{
    return new BoxDynamicViscosityFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxDynamicViscosityFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
TagLocalDynamicsFunctional3D<T, Descriptor>::TagLocalDynamicsFunctional3D(
    int dynamicsId_, int tag_) :
    dynamicsId(dynamicsId_), tag(tag_)
{ }

template <typename T, template <typename U> class Descriptor>
void TagLocalDynamicsFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &tags)
{
    Dot3D offset = computeRelativeDisplacement(lattice, tags);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (lattice.get(iX, iY, iZ).getDynamics().getId() == dynamicsId) {
                    tags.get(iX + offset.x, iY + offset.y, iZ + offset.z) = tag;
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
TagLocalDynamicsFunctional3D<T, Descriptor> *TagLocalDynamicsFunctional3D<T, Descriptor>::clone()
    const
{
    return new TagLocalDynamicsFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT TagLocalDynamicsFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void TagLocalDynamicsFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Lattice.
    modified[1] = modif::staticVariables;  // Tags.
}

/* *************** PART II ******************************************* */
/* *************** Analysis of the scalar-field ********************** */
/* ******************************************************************* */

/* *************** Reductive Data Functionals for scalar-field ******* */

template <typename T>
BoxScalarSumFunctional3D<T>::BoxScalarSumFunctional3D() :
    sumScalarId(this->getStatistics().subscribeSum())
{ }

template <typename T>
void BoxScalarSumFunctional3D<T>::process(Box3D domain, ScalarField3D<T> &scalarField)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                statistics.gatherSum(sumScalarId, (double)scalarField.get(iX, iY, iZ));
            }
        }
    }
}

template <typename T>
BoxScalarSumFunctional3D<T> *BoxScalarSumFunctional3D<T>::clone() const
{
    return new BoxScalarSumFunctional3D<T>(*this);
}

template <typename T>
T BoxScalarSumFunctional3D<T>::getSumScalar() const
{
    double doubleSum = this->getStatistics().getSum(sumScalarId);
    // The sum is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    if (std::numeric_limits<T>::is_integer) {
        return (T)util::roundToInt(doubleSum);
    }
    return (T)doubleSum;
}

template <typename T>
MaskedBoxScalarSumFunctional3D<T>::MaskedBoxScalarSumFunctional3D(int flag_) :
    sumScalarId(this->getStatistics().subscribeSum()), flag(flag_)
{ }

template <typename T>
void MaskedBoxScalarSumFunctional3D<T>::process(
    Box3D domain, ScalarField3D<T> &scalarField, ScalarField3D<int> &mask)
{
    Dot3D offset = computeRelativeDisplacement(scalarField, mask);
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (mask.get(iX + offset.x, iY + offset.y, iZ + offset.z) == flag) {
                    statistics.gatherSum(sumScalarId, (double)scalarField.get(iX, iY, iZ));
                }
            }
        }
    }
}

template <typename T>
MaskedBoxScalarSumFunctional3D<T> *MaskedBoxScalarSumFunctional3D<T>::clone() const
{
    return new MaskedBoxScalarSumFunctional3D<T>(*this);
}

template <typename T>
T MaskedBoxScalarSumFunctional3D<T>::getSumScalar() const
{
    double doubleSum = this->getStatistics().getSum(sumScalarId);
    // The sum is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    if (std::numeric_limits<T>::is_integer) {
        return (T)util::roundToInt(doubleSum);
    }
    return (T)doubleSum;
}

template <typename T>
BoxScalarIntSumFunctional3D<T>::BoxScalarIntSumFunctional3D() :
    sumScalarId(this->getStatistics().subscribeIntSum())
{ }

template <typename T>
void BoxScalarIntSumFunctional3D<T>::process(Box3D domain, ScalarField3D<T> &scalarField)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                statistics.gatherIntSum(sumScalarId, (plint)scalarField.get(iX, iY, iZ));
            }
        }
    }
}

template <typename T>
BoxScalarIntSumFunctional3D<T> *BoxScalarIntSumFunctional3D<T>::clone() const
{
    return new BoxScalarIntSumFunctional3D<T>(*this);
}

template <typename T>
plint BoxScalarIntSumFunctional3D<T>::getSumScalar() const
{
    return this->getStatistics().getIntSum(sumScalarId);
}

template <typename T>
MaskedBoxScalarAverageFunctional3D<T>::MaskedBoxScalarAverageFunctional3D(int flag_) :
    averageScalarId(this->getStatistics().subscribeAverage()), flag(flag_)
{ }

template <typename T>
void MaskedBoxScalarAverageFunctional3D<T>::process(
    Box3D domain, ScalarField3D<T> &scalarField, ScalarField3D<int> &mask)
{
    Dot3D offset = computeRelativeDisplacement(scalarField, mask);
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (mask.get(iX + offset.x, iY + offset.y, iZ + offset.z) == flag) {
                    statistics.gatherAverage(averageScalarId, (double)scalarField.get(iX, iY, iZ));
                    statistics.incrementStats();
                }
            }
        }
    }
}

template <typename T>
MaskedBoxScalarAverageFunctional3D<T> *MaskedBoxScalarAverageFunctional3D<T>::clone() const
{
    return new MaskedBoxScalarAverageFunctional3D<T>(*this);
}

template <typename T>
T MaskedBoxScalarAverageFunctional3D<T>::getAverageScalar() const
{
    double doubleAverage = this->getStatistics().getAverage(averageScalarId);
    // The average is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    if (std::numeric_limits<T>::is_integer) {
        return (T)util::roundToInt(doubleAverage);
    }
    return (T)doubleAverage;
}

template <typename T>
BoxScalarMinFunctional3D<T>::BoxScalarMinFunctional3D() :
    maxScalarId(this->getStatistics().subscribeMax())
{ }

template <typename T>
void BoxScalarMinFunctional3D<T>::process(Box3D domain, ScalarField3D<T> &scalarField)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                // BlockStatistics computes only maximum, no minimum. Therefore,
                //   the relation min(x) = -max(-x) is used.
                statistics.gatherMax(maxScalarId, -(double)scalarField.get(iX, iY, iZ));
            }
        }
    }
}

template <typename T>
BoxScalarMinFunctional3D<T> *BoxScalarMinFunctional3D<T>::clone() const
{
    return new BoxScalarMinFunctional3D<T>(*this);
}

template <typename T>
T BoxScalarMinFunctional3D<T>::getMinScalar() const
{
    // The minus sign accounts for the relation min(x) = -max(-x).
    // The sum is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    double doubleMin = -this->getStatistics().getMax(maxScalarId);
    if (std::numeric_limits<T>::is_integer) {
        return (T)util::roundToInt(doubleMin);
    }
    return (T)doubleMin;
}

template <typename T>
MaskedBoxScalarMinFunctional3D<T>::MaskedBoxScalarMinFunctional3D(int flag_) :
    maxScalarId(this->getStatistics().subscribeMax()), flag(flag_)
{ }

template <typename T>
void MaskedBoxScalarMinFunctional3D<T>::process(
    Box3D domain, ScalarField3D<T> &scalarField, ScalarField3D<int> &mask)
{
    BlockStatistics &statistics = this->getStatistics();
    Dot3D offset = computeRelativeDisplacement(scalarField, mask);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (mask.get(iX + offset.x, iY + offset.y, iZ + offset.z) == flag) {
                    // BlockStatistics computes only maximum, no minimum. Therefore,
                    //   the relation min(x) = -max(-x) is used.
                    statistics.gatherMax(maxScalarId, -(double)scalarField.get(iX, iY, iZ));
                }
            }
        }
    }
}

template <typename T>
MaskedBoxScalarMinFunctional3D<T> *MaskedBoxScalarMinFunctional3D<T>::clone() const
{
    return new MaskedBoxScalarMinFunctional3D<T>(*this);
}

template <typename T>
T MaskedBoxScalarMinFunctional3D<T>::getMinScalar() const
{
    // The minus sign accounts for the relation min(x) = -max(-x).
    // The sum is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    double doubleMin = -this->getStatistics().getMax(maxScalarId);
    if (std::numeric_limits<T>::is_integer) {
        return (T)util::roundToInt(doubleMin);
    }
    return (T)doubleMin;
}

template <typename T>
BoxScalarMaxFunctional3D<T>::BoxScalarMaxFunctional3D() :
    maxScalarId(this->getStatistics().subscribeMax())
{ }

template <typename T>
void BoxScalarMaxFunctional3D<T>::process(Box3D domain, ScalarField3D<T> &scalarField)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                statistics.gatherMax(maxScalarId, (double)scalarField.get(iX, iY, iZ));
            }
        }
    }
}

template <typename T>
BoxScalarMaxFunctional3D<T> *BoxScalarMaxFunctional3D<T>::clone() const
{
    return new BoxScalarMaxFunctional3D<T>(*this);
}

template <typename T>
T BoxScalarMaxFunctional3D<T>::getMaxScalar() const
{
    // The sum is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    double doubleMax = this->getStatistics().getMax(maxScalarId);
    if (std::numeric_limits<T>::is_integer) {
        return (T)util::roundToInt(doubleMax);
    }
    return (T)doubleMax;
}

template <typename T>
MaskedBoxScalarMaxFunctional3D<T>::MaskedBoxScalarMaxFunctional3D(int flag_) :
    maxScalarId(this->getStatistics().subscribeMax()), flag(flag_)
{ }

template <typename T>
void MaskedBoxScalarMaxFunctional3D<T>::process(
    Box3D domain, ScalarField3D<T> &scalarField, ScalarField3D<int> &mask)
{
    BlockStatistics &statistics = this->getStatistics();
    Dot3D offset = computeRelativeDisplacement(scalarField, mask);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (mask.get(iX + offset.x, iY + offset.y, iZ + offset.z) == flag) {
                    statistics.gatherMax(maxScalarId, (double)scalarField.get(iX, iY, iZ));
                }
            }
        }
    }
}

template <typename T>
MaskedBoxScalarMaxFunctional3D<T> *MaskedBoxScalarMaxFunctional3D<T>::clone() const
{
    return new MaskedBoxScalarMaxFunctional3D<T>(*this);
}

template <typename T>
T MaskedBoxScalarMaxFunctional3D<T>::getMaxScalar() const
{
    // The sum is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    double doubleMax = this->getStatistics().getMax(maxScalarId);
    if (std::numeric_limits<T>::is_integer) {
        return (T)util::roundToInt(doubleMax);
    }
    return (T)doubleMax;
}

template <typename T>
BoundedBoxScalarSumFunctional3D<T>::BoundedBoxScalarSumFunctional3D() :
    sumScalarId(this->getStatistics().subscribeSum())
{ }

template <typename T>
void BoundedBoxScalarSumFunctional3D<T>::processBulk(Box3D domain, ScalarField3D<T> &scalarField)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                statistics.gatherSum(sumScalarId, (double)scalarField.get(iX, iY, iZ));
            }
        }
    }
}

template <typename T>
void BoundedBoxScalarSumFunctional3D<T>::processPlane(
    int direction, int orientation, Box3D domain, ScalarField3D<T> &scalarField)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                // Plane boundary nodes have a weight of 0.5, because only 50% of the
                //   cell centered at the node is inside the computational domain.
                statistics.gatherSum(sumScalarId, (double)scalarField.get(iX, iY, iZ) / 2.);
            }
        }
    }
}

template <typename T>
void BoundedBoxScalarSumFunctional3D<T>::processEdge(
    int plane, int normal1, int normal2, Box3D domain, ScalarField3D<T> &scalarField)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                // Edge nodes have a weight of 0.25, because only 25% of the
                //   cell centered at the node is inside the computational domain.
                statistics.gatherSum(sumScalarId, (double)scalarField.get(iX, iY, iZ) / 4.);
            }
        }
    }
}

template <typename T>
void BoundedBoxScalarSumFunctional3D<T>::processCorner(
    int normalX, int normalY, int normalZ, Box3D domain, ScalarField3D<T> &scalarField)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                // Corner nodes have a weight of 0.125, because only 1/8 of the
                //   cell centered at the node is inside the computational domain.
                statistics.gatherSum(sumScalarId, (double)scalarField.get(iX, iY, iZ) / 8.);
            }
        }
    }
}

template <typename T>
BoundedBoxScalarSumFunctional3D<T> *BoundedBoxScalarSumFunctional3D<T>::clone() const
{
    return new BoundedBoxScalarSumFunctional3D<T>(*this);
}

template <typename T>
T BoundedBoxScalarSumFunctional3D<T>::getSumScalar() const
{
    double doubleSum = this->getStatistics().getSum(sumScalarId);
    // The sum is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    if (std::numeric_limits<T>::is_integer) {
        return (T)util::roundToInt(doubleSum);
    }
    return (T)doubleSum;
}

template <typename T1, typename T2>
void CopyConvertScalarFunctional3D<T1, T2>::process(
    Box3D domain, ScalarField3D<T1> &field1, ScalarField3D<T2> &field2)
{
    Dot3D offset = computeRelativeDisplacement(field1, field2);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                field2.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    (T2)field1.get(iX, iY, iZ);
            }
        }
    }
}

template <typename T1, typename T2>
CopyConvertScalarFunctional3D<T1, T2> *CopyConvertScalarFunctional3D<T1, T2>::clone() const
{
    return new CopyConvertScalarFunctional3D<T1, T2>(*this);
}

template <typename T1, typename T2>
void CopyConvertScalarFunctional3D<T1, T2>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T1, typename T2>
BlockDomain::DomainT CopyConvertScalarFunctional3D<T1, T2>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* *************** Data Functionals for scalar-fields **************** */

template <typename T>
void ExtractScalarSubDomainFunctional3D<T>::process(
    Box3D domain, ScalarField3D<T> &field1, ScalarField3D<T> &field2)
{
    Dot3D offset = computeRelativeDisplacement(field1, field2);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                field2.get(iX + offset.x, iY + offset.y, iZ + offset.z) = field1.get(iX, iY, iZ);
            }
        }
    }
}

template <typename T>
ExtractScalarSubDomainFunctional3D<T> *ExtractScalarSubDomainFunctional3D<T>::clone() const
{
    return new ExtractScalarSubDomainFunctional3D<T>(*this);
}

template <typename T>
void ExtractScalarSubDomainFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT ExtractScalarSubDomainFunctional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** compute sqrt functional 3D ************************************* */

template <typename T>
void ComputeScalarSqrtFunctional3D<T>::process(
    Box3D domain, ScalarField3D<T> &A, ScalarField3D<T> &B)
{
    Dot3D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                B.get(iX + offset.x, iY + offset.y, iZ + offset.z) = std::sqrt(A.get(iX, iY, iZ));
            }
        }
    }
}

template <typename T>
ComputeScalarSqrtFunctional3D<T> *ComputeScalarSqrtFunctional3D<T>::clone() const
{
    return new ComputeScalarSqrtFunctional3D<T>(*this);
}

template <typename T>
void ComputeScalarSqrtFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT ComputeScalarSqrtFunctional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** compute fabs functional 3D ************************************* */

template <typename T>
void ComputeAbsoluteValueFunctional3D<T>::process(
    Box3D domain, ScalarField3D<T> &A, ScalarField3D<T> &B)
{
    Dot3D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                B.get(iX + offset.x, iY + offset.y, iZ + offset.z) = std::fabs(A.get(iX, iY, iZ));
            }
        }
    }
}

template <typename T>
ComputeAbsoluteValueFunctional3D<T> *ComputeAbsoluteValueFunctional3D<T>::clone() const
{
    return new ComputeAbsoluteValueFunctional3D<T>(*this);
}

template <typename T>
void ComputeAbsoluteValueFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT ComputeAbsoluteValueFunctional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** compute sqrt functional 3D ************************************* */

template <typename T, int nDim>
void ComputeTensorSqrtFunctional3D<T, nDim>::process(
    Box3D domain, TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B)
{
    Dot3D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                for (plint iD = 0; iD < nDim; ++iD) {
                    B.get(iX + offset.x, iY + offset.y, iZ + offset.z)[iD] =
                        std::sqrt(A.get(iX, iY, iZ)[iD]);
                }
            }
        }
    }
}

template <typename T, int nDim>
ComputeTensorSqrtFunctional3D<T, nDim> *ComputeTensorSqrtFunctional3D<T, nDim>::clone() const
{
    return new ComputeTensorSqrtFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void ComputeTensorSqrtFunctional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT ComputeTensorSqrtFunctional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_lt_alpha_functional3D ************************************* */

template <typename T>
A_lt_alpha_functional3D<T>::A_lt_alpha_functional3D(T alpha_) : alpha(alpha_)
{ }

template <typename T>
void A_lt_alpha_functional3D<T>::process(
    Box3D domain, ScalarField3D<T> &A, ScalarField3D<int> &result)
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    A.get(iX, iY, iZ) < alpha ? 1 : 0;
            }
        }
    }
}

template <typename T>
A_lt_alpha_functional3D<T> *A_lt_alpha_functional3D<T>::clone() const
{
    return new A_lt_alpha_functional3D<T>(*this);
}

template <typename T>
void A_lt_alpha_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_lt_alpha_functional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_gt_alpha_functional3D ************************************* */

template <typename T>
A_gt_alpha_functional3D<T>::A_gt_alpha_functional3D(T alpha_) : alpha(alpha_)
{ }

template <typename T>
void A_gt_alpha_functional3D<T>::process(
    Box3D domain, ScalarField3D<T> &A, ScalarField3D<int> &result)
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    A.get(iX, iY, iZ) > alpha ? 1 : 0;
            }
        }
    }
}

template <typename T>
A_gt_alpha_functional3D<T> *A_gt_alpha_functional3D<T>::clone() const
{
    return new A_gt_alpha_functional3D<T>(*this);
}

template <typename T>
void A_gt_alpha_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_gt_alpha_functional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_plus_alpha_functional3D ************************************* */

template <typename T>
A_plus_alpha_functional3D<T>::A_plus_alpha_functional3D(T alpha_) : alpha(alpha_)
{ }

template <typename T>
void A_plus_alpha_functional3D<T>::process(
    Box3D domain, ScalarField3D<T> &A, ScalarField3D<T> &result)
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offset.x, iY + offset.y, iZ + offset.z) = A.get(iX, iY, iZ) + alpha;
            }
        }
    }
}

template <typename T>
A_plus_alpha_functional3D<T> *A_plus_alpha_functional3D<T>::clone() const
{
    return new A_plus_alpha_functional3D<T>(*this);
}

template <typename T>
void A_plus_alpha_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_plus_alpha_functional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_minus_alpha_functional3D ************************************** */

template <typename T>
A_minus_alpha_functional3D<T>::A_minus_alpha_functional3D(T alpha_) : alpha(alpha_)
{ }

template <typename T>
void A_minus_alpha_functional3D<T>::process(
    Box3D domain, ScalarField3D<T> &A, ScalarField3D<T> &result)
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offset.x, iY + offset.y, iZ + offset.z) = A.get(iX, iY, iZ) - alpha;
            }
        }
    }
}

template <typename T>
A_minus_alpha_functional3D<T> *A_minus_alpha_functional3D<T>::clone() const
{
    return new A_minus_alpha_functional3D<T>(*this);
}

template <typename T>
void A_minus_alpha_functional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_minus_alpha_functional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Alpha_minus_A_functional3D ************************************* */

template <typename T>
Alpha_minus_A_functional3D<T>::Alpha_minus_A_functional3D(T alpha_) : alpha(alpha_)
{ }

template <typename T>
void Alpha_minus_A_functional3D<T>::process(
    Box3D domain, ScalarField3D<T> &A, ScalarField3D<T> &result)
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offset.x, iY + offset.y, iZ + offset.z) = alpha - A.get(iX, iY, iZ);
            }
        }
    }
}

template <typename T>
Alpha_minus_A_functional3D<T> *Alpha_minus_A_functional3D<T>::clone() const
{
    return new Alpha_minus_A_functional3D<T>(*this);
}

template <typename T>
void Alpha_minus_A_functional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT Alpha_minus_A_functional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_times_alpha_functional3D ************************************* */

template <typename T>
A_times_alpha_functional3D<T>::A_times_alpha_functional3D(T alpha_) : alpha(alpha_)
{ }

template <typename T>
void A_times_alpha_functional3D<T>::process(
    Box3D domain, ScalarField3D<T> &A, ScalarField3D<T> &result)
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offset.x, iY + offset.y, iZ + offset.z) = A.get(iX, iY, iZ) * alpha;
            }
        }
    }
}

template <typename T>
A_times_alpha_functional3D<T> *A_times_alpha_functional3D<T>::clone() const
{
    return new A_times_alpha_functional3D<T>(*this);
}

template <typename T>
void A_times_alpha_functional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_times_alpha_functional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_dividedBy_alpha_functional3D ************************************* */

template <typename T>
A_dividedBy_alpha_functional3D<T>::A_dividedBy_alpha_functional3D(T alpha_) : alpha(alpha_)
{ }

template <typename T>
void A_dividedBy_alpha_functional3D<T>::process(
    Box3D domain, ScalarField3D<T> &A, ScalarField3D<T> &result)
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offset.x, iY + offset.y, iZ + offset.z) = A.get(iX, iY, iZ) / alpha;
            }
        }
    }
}

template <typename T>
A_dividedBy_alpha_functional3D<T> *A_dividedBy_alpha_functional3D<T>::clone() const
{
    return new A_dividedBy_alpha_functional3D<T>(*this);
}

template <typename T>
void A_dividedBy_alpha_functional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_dividedBy_alpha_functional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Alpha_dividedBy_A_functional3D ************************************* */

template <typename T>
Alpha_dividedBy_A_functional3D<T>::Alpha_dividedBy_A_functional3D(T alpha_) : alpha(alpha_)
{ }

template <typename T>
void Alpha_dividedBy_A_functional3D<T>::process(
    Box3D domain, ScalarField3D<T> &A, ScalarField3D<T> &result)
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offset.x, iY + offset.y, iZ + offset.z) = alpha / A.get(iX, iY, iZ);
            }
        }
    }
}

template <typename T>
Alpha_dividedBy_A_functional3D<T> *Alpha_dividedBy_A_functional3D<T>::clone() const
{
    return new Alpha_dividedBy_A_functional3D<T>(*this);
}

template <typename T>
void Alpha_dividedBy_A_functional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT Alpha_dividedBy_A_functional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_plus_alpha_inplace_functional3D ************************************* */

template <typename T>
A_plus_alpha_inplace_functional3D<T>::A_plus_alpha_inplace_functional3D(T alpha_) : alpha(alpha_)
{ }

template <typename T>
void A_plus_alpha_inplace_functional3D<T>::process(Box3D domain, ScalarField3D<T> &A)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                A.get(iX, iY, iZ) += alpha;
            }
        }
    }
}

template <typename T>
A_plus_alpha_inplace_functional3D<T> *A_plus_alpha_inplace_functional3D<T>::clone() const
{
    return new A_plus_alpha_inplace_functional3D<T>(*this);
}

template <typename T>
void A_plus_alpha_inplace_functional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_plus_alpha_inplace_functional3D<T>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

/* ******** A_minus_alpha_inplace_functional3D ************************************** */

template <typename T>
A_minus_alpha_inplace_functional3D<T>::A_minus_alpha_inplace_functional3D(T alpha_) : alpha(alpha_)
{ }

template <typename T>
void A_minus_alpha_inplace_functional3D<T>::process(Box3D domain, ScalarField3D<T> &A)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                A.get(iX, iY, iZ) -= alpha;
            }
        }
    }
}

template <typename T>
A_minus_alpha_inplace_functional3D<T> *A_minus_alpha_inplace_functional3D<T>::clone() const
{
    return new A_minus_alpha_inplace_functional3D<T>(*this);
}

template <typename T>
void A_minus_alpha_inplace_functional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_minus_alpha_inplace_functional3D<T>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

/* ******** A_times_alpha_inplace_functional3D ************************************* */

template <typename T>
A_times_alpha_inplace_functional3D<T>::A_times_alpha_inplace_functional3D(T alpha_) : alpha(alpha_)
{ }

template <typename T>
void A_times_alpha_inplace_functional3D<T>::process(Box3D domain, ScalarField3D<T> &A)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                A.get(iX, iY, iZ) *= alpha;
            }
        }
    }
}

template <typename T>
A_times_alpha_inplace_functional3D<T> *A_times_alpha_inplace_functional3D<T>::clone() const
{
    return new A_times_alpha_inplace_functional3D<T>(*this);
}

template <typename T>
void A_times_alpha_inplace_functional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_times_alpha_inplace_functional3D<T>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

/* ******** A_dividedBy_alpha_inplace_functional3D ************************************* */

template <typename T>
A_dividedBy_alpha_inplace_functional3D<T>::A_dividedBy_alpha_inplace_functional3D(T alpha_) :
    alpha(alpha_)
{ }

template <typename T>
void A_dividedBy_alpha_inplace_functional3D<T>::process(Box3D domain, ScalarField3D<T> &A)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                A.get(iX, iY, iZ) /= alpha;
            }
        }
    }
}

template <typename T>
A_dividedBy_alpha_inplace_functional3D<T> *A_dividedBy_alpha_inplace_functional3D<T>::clone() const
{
    return new A_dividedBy_alpha_inplace_functional3D<T>(*this);
}

template <typename T>
void A_dividedBy_alpha_inplace_functional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_dividedBy_alpha_inplace_functional3D<T>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

/* ******** A_lt_B_functional3D ****************************************** */

template <typename T>
void A_lt_B_functional3D<T>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);
    ScalarField3D<T> &A = *dynamic_cast<ScalarField3D<T> *>(fields[0]);
    ScalarField3D<T> &B = *dynamic_cast<ScalarField3D<T> *>(fields[1]);
    ScalarField3D<int> &result = *dynamic_cast<ScalarField3D<int> *>(fields[2]);
    Dot3D offsetB = computeRelativeDisplacement(A, B);
    Dot3D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offsetResult.x, iY + offsetResult.y, iZ + offsetResult.z) =
                    A.get(iX, iY, iZ) < B.get(iX + offsetB.x, iY + offsetB.y, iZ + offsetB.z) ? 1
                                                                                              : 0;
            }
        }
    }
}

template <typename T>
A_lt_B_functional3D<T> *A_lt_B_functional3D<T>::clone() const
{
    return new A_lt_B_functional3D<T>(*this);
}

template <typename T>
void A_lt_B_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_lt_B_functional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_gt_B_functional3D ****************************************** */

template <typename T>
void A_gt_B_functional3D<T>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);
    ScalarField3D<T> &A = *dynamic_cast<ScalarField3D<T> *>(fields[0]);
    ScalarField3D<T> &B = *dynamic_cast<ScalarField3D<T> *>(fields[1]);
    ScalarField3D<int> &result = *dynamic_cast<ScalarField3D<int> *>(fields[2]);
    Dot3D offsetB = computeRelativeDisplacement(A, B);
    Dot3D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offsetResult.x, iY + offsetResult.y, iZ + offsetResult.z) =
                    A.get(iX, iY, iZ) > B.get(iX + offsetB.x, iY + offsetB.y, iZ + offsetB.z) ? 1
                                                                                              : 0;
            }
        }
    }
}

template <typename T>
A_gt_B_functional3D<T> *A_gt_B_functional3D<T>::clone() const
{
    return new A_gt_B_functional3D<T>(*this);
}

template <typename T>
void A_gt_B_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_gt_B_functional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_plus_B_functional3D ****************************************** */

template <typename T>
void A_plus_B_functional3D<T>::process(Box3D domain, std::vector<ScalarField3D<T> *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);
    ScalarField3D<T> &A = *fields[0];
    ScalarField3D<T> &B = *fields[1];
    ScalarField3D<T> &result = *fields[2];
    Dot3D offsetB = computeRelativeDisplacement(A, B);
    Dot3D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offsetResult.x, iY + offsetResult.y, iZ + offsetResult.z) =
                    A.get(iX, iY, iZ) + B.get(iX + offsetB.x, iY + offsetB.y, iZ + offsetB.z);
            }
        }
    }
}

template <typename T>
A_plus_B_functional3D<T> *A_plus_B_functional3D<T>::clone() const
{
    return new A_plus_B_functional3D<T>(*this);
}

template <typename T>
void A_plus_B_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_plus_B_functional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_minus_B_functional3D ****************************************** */

template <typename T>
void A_minus_B_functional3D<T>::process(Box3D domain, std::vector<ScalarField3D<T> *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);
    ScalarField3D<T> &A = *fields[0];
    ScalarField3D<T> &B = *fields[1];
    ScalarField3D<T> &result = *fields[2];
    Dot3D offsetB = computeRelativeDisplacement(A, B);
    Dot3D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offsetResult.x, iY + offsetResult.y, iZ + offsetResult.z) =
                    A.get(iX, iY, iZ) - B.get(iX + offsetB.x, iY + offsetB.y, iZ + offsetB.z);
            }
        }
    }
}

template <typename T>
A_minus_B_functional3D<T> *A_minus_B_functional3D<T>::clone() const
{
    return new A_minus_B_functional3D<T>(*this);
}

template <typename T>
void A_minus_B_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_minus_B_functional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_times_B_functional3D ****************************************** */

template <typename T>
void A_times_B_functional3D<T>::process(Box3D domain, std::vector<ScalarField3D<T> *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);
    ScalarField3D<T> &A = *fields[0];
    ScalarField3D<T> &B = *fields[1];
    ScalarField3D<T> &result = *fields[2];
    Dot3D offsetB = computeRelativeDisplacement(A, B);
    Dot3D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offsetResult.x, iY + offsetResult.y, iZ + offsetResult.z) =
                    A.get(iX, iY, iZ) * B.get(iX + offsetB.x, iY + offsetB.y, iZ + offsetB.z);
            }
        }
    }
}

template <typename T>
A_times_B_functional3D<T> *A_times_B_functional3D<T>::clone() const
{
    return new A_times_B_functional3D<T>(*this);
}

template <typename T>
void A_times_B_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_times_B_functional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_dividedBy_B_functional3D ****************************************** */

template <typename T>
void A_dividedBy_B_functional3D<T>::process(Box3D domain, std::vector<ScalarField3D<T> *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);
    ScalarField3D<T> &A = *fields[0];
    ScalarField3D<T> &B = *fields[1];
    ScalarField3D<T> &result = *fields[2];
    Dot3D offsetB = computeRelativeDisplacement(A, B);
    Dot3D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offsetResult.x, iY + offsetResult.y, iZ + offsetResult.z) =
                    A.get(iX, iY, iZ) / B.get(iX + offsetB.x, iY + offsetB.y, iZ + offsetB.z);
            }
        }
    }
}

template <typename T>
A_dividedBy_B_functional3D<T> *A_dividedBy_B_functional3D<T>::clone() const
{
    return new A_dividedBy_B_functional3D<T>(*this);
}

template <typename T>
void A_dividedBy_B_functional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_dividedBy_B_functional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_plus_B_inplace_functional3D ****************************************** */

template <typename T>
void A_plus_B_inplace_functional3D<T>::process(
    Box3D domain, ScalarField3D<T> &A, ScalarField3D<T> &B)
{
    Dot3D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                A.get(iX, iY, iZ) += B.get(iX + offset.x, iY + offset.y, iZ + offset.z);
            }
        }
    }
}

template <typename T>
A_plus_B_inplace_functional3D<T> *A_plus_B_inplace_functional3D<T>::clone() const
{
    return new A_plus_B_inplace_functional3D<T>(*this);
}

template <typename T>
void A_plus_B_inplace_functional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T>
BlockDomain::DomainT A_plus_B_inplace_functional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_minus_B_inplace_functional3D ****************************************** */

template <typename T>
void A_minus_B_inplace_functional3D<T>::process(
    Box3D domain, ScalarField3D<T> &A, ScalarField3D<T> &B)
{
    Dot3D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                A.get(iX, iY, iZ) -= B.get(iX + offset.x, iY + offset.y, iZ + offset.z);
            }
        }
    }
}

template <typename T>
A_minus_B_inplace_functional3D<T> *A_minus_B_inplace_functional3D<T>::clone() const
{
    return new A_minus_B_inplace_functional3D<T>(*this);
}

template <typename T>
void A_minus_B_inplace_functional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T>
BlockDomain::DomainT A_minus_B_inplace_functional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_times_B_inplace_functional3D ****************************************** */

template <typename T>
void A_times_B_inplace_functional3D<T>::process(
    Box3D domain, ScalarField3D<T> &A, ScalarField3D<T> &B)
{
    Dot3D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                A.get(iX, iY, iZ) *= B.get(iX + offset.x, iY + offset.y, iZ + offset.z);
            }
        }
    }
}

template <typename T>
A_times_B_inplace_functional3D<T> *A_times_B_inplace_functional3D<T>::clone() const
{
    return new A_times_B_inplace_functional3D<T>(*this);
}

template <typename T>
void A_times_B_inplace_functional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T>
BlockDomain::DomainT A_times_B_inplace_functional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_dividedBy_B_inplace_functional3D ****************************************** */

template <typename T>
void A_dividedBy_B_inplace_functional3D<T>::process(
    Box3D domain, ScalarField3D<T> &A, ScalarField3D<T> &B)
{
    Dot3D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                A.get(iX, iY, iZ) /= B.get(iX + offset.x, iY + offset.y, iZ + offset.z);
            }
        }
    }
}

template <typename T>
A_dividedBy_B_inplace_functional3D<T> *A_dividedBy_B_inplace_functional3D<T>::clone() const
{
    return new A_dividedBy_B_inplace_functional3D<T>(*this);
}

template <typename T>
void A_dividedBy_B_inplace_functional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T>
BlockDomain::DomainT A_dividedBy_B_inplace_functional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** MultiScalarField inplace bounding operations ******************************************
 */

template <typename T>
void UniformlyBoundScalarField3D<T>::process(Box3D domain, ScalarField3D<T> &data)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T &val = data.get(iX, iY, iZ);
                if (std::fabs(val) > bound)
                    val = val > T() ? bound : -bound;
            }
        }
    }
}

template <typename T>
void BoundScalarField3D<T>::process(Box3D domain, ScalarField3D<T> &data)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T &val = data.get(iX, iY, iZ);
                if (val < lowerBound) {
                    val = lowerBound;
                } else if (val > upperBound) {
                    val = upperBound;
                }
            }
        }
    }
}

/* ************* Class LBMsmoothen3D ******************* */

template <typename T, template <typename U> class Descriptor>
void LBMsmoothen3D<T, Descriptor>::process(
    Box3D domain, ScalarField3D<T> &data, ScalarField3D<T> &result)
{
    typedef Descriptor<T> D;
    Dot3D offset = computeRelativeDisplacement(data, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offset.x, iY + offset.y, iZ + offset.z) = (T)0;
                T sum = (T)0;
                for (plint iPop = 1; iPop < D::q; ++iPop) {
                    plint nextX = iX + D::c[iPop][0];
                    plint nextY = iY + D::c[iPop][1];
                    plint nextZ = iZ + D::c[iPop][2];
                    sum += D::t[iPop];
                    result.get(iX + offset.x, iY + offset.y, iZ + offset.z) +=
                        D::t[iPop] * data.get(nextX, nextY, nextZ);
                }
                result.get(iX + offset.x, iY + offset.y, iZ + offset.z) /= sum;
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
LBMsmoothen3D<T, Descriptor> *LBMsmoothen3D<T, Descriptor>::clone() const
{
    return new LBMsmoothen3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void LBMsmoothen3D<T, Descriptor>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

/* ************* Class LBMsmoothenInPlace3D ******************* */

template <typename T, template <typename U> class Descriptor>
void LBMsmoothenInPlace3D<T, Descriptor>::process(Box3D domain, ScalarField3D<T> &data)
{
    typedef Descriptor<T> D;

    ScalarField3D<T> smoothData(domain.getNx(), domain.getNy(), domain.getNz());
    Dot3D offset(-domain.x0, -domain.y0, -domain.z0);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                smoothData.get(iX + offset.x, iY + offset.y, iZ + offset.z) = (T)0;
                T sum = (T)0;
                for (plint iPop = 1; iPop < D::q; ++iPop) {
                    plint nextX = iX + D::c[iPop][0];
                    plint nextY = iY + D::c[iPop][1];
                    plint nextZ = iZ + D::c[iPop][2];
                    sum += D::t[iPop];
                    smoothData.get(iX + offset.x, iY + offset.y, iZ + offset.z) +=
                        D::t[iPop] * data.get(nextX, nextY, nextZ);
                }
                smoothData.get(iX + offset.x, iY + offset.y, iZ + offset.z) /= sum;
            }
        }
    }

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                data.get(iX, iY, iZ) = smoothData.get(iX + offset.x, iY + offset.y, iZ + offset.z);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
LBMsmoothenInPlace3D<T, Descriptor> *LBMsmoothenInPlace3D<T, Descriptor>::clone() const
{
    return new LBMsmoothenInPlace3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void LBMsmoothenInPlace3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************* Class Smoothen3D ******************* */

template <typename T>
void Smoothen3D<T>::process(Box3D domain, ScalarField3D<T> &data, ScalarField3D<T> &result)
{
    Dot3D offset = computeRelativeDisplacement(data, result);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T *res = &result.get(iX + offset.x, iY + offset.y, iZ + offset.z);
                *res = (T)0;
                int n = 0;
                for (int i = -1; i < 2; i++) {
                    plint nextX = iX + i;
                    for (int j = -1; j < 2; j++) {
                        plint nextY = iY + j;
                        for (int k = -1; k < 2; k++) {
                            plint nextZ = iZ + k;
                            if (!(i == 0 && j == 0 && k == 0)) {
                                n++;
                                *res += data.get(nextX, nextY, nextZ);
                            }
                        }
                    }
                }
                *res /= (T)n;
            }
        }
    }
}

template <typename T>
Smoothen3D<T> *Smoothen3D<T>::clone() const
{
    return new Smoothen3D<T>(*this);
}

template <typename T>
void Smoothen3D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

/* ************* Class MollifyScalar3D ******************* */

template <typename T>
MollifyScalar3D<T>::MollifyScalar3D(T l_, plint d_, Box3D globalDomain_, int exclusionFlag_) :
    l(l_), d(d_), globalDomain(globalDomain_), exclusionFlag(exclusionFlag_)
{ }

template <typename T>
void MollifyScalar3D<T>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_ASSERT(fields.size() == 3);
    ScalarField3D<T> *scalar = dynamic_cast<ScalarField3D<T> *>(fields[0]);
    PLB_ASSERT(scalar);
    ScalarField3D<int> *flag = dynamic_cast<ScalarField3D<int> *>(fields[1]);
    PLB_ASSERT(flag);
    ScalarField3D<T> *result = dynamic_cast<ScalarField3D<T> *>(fields[2]);
    PLB_ASSERT(result);

    Dot3D absOfs = scalar->getLocation();

    Dot3D ofsF = computeRelativeDisplacement(*scalar, *flag);
    Dot3D ofsR = computeRelativeDisplacement(*scalar, *result);

    static T pi = std::acos((T)-1);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (flag->get(iX + ofsF.x, iY + ofsF.y, iZ + ofsF.z) == exclusionFlag) {
                    result->get(iX + ofsR.x, iY + ofsR.y, iZ + ofsR.z) = scalar->get(iX, iY, iZ);
                    continue;
                }

                T sum = 0.0;
                for (plint dx = -d; dx <= d; dx++) {
                    plint i = iX + dx;
                    for (plint dy = -d; dy <= d; dy++) {
                        plint j = iY + dy;
                        for (plint dz = -d; dz <= d; dz++) {
                            plint k = iZ + dz;
                            if (contained(i + absOfs.x, j + absOfs.y, k + absOfs.z, globalDomain)) {
                                if (flag->get(i + ofsF.x, j + ofsF.y, k + ofsF.z) != exclusionFlag)
                                {
                                    T r = std::sqrt((T)dx * dx + (T)dy * dy + (T)dz * dz);
                                    T integrand = scalar->get(i, j, k);
                                    T mollifier = 0.0;
                                    if (r < l) {
                                        mollifier = (1.0 + std::cos(pi * r / l));
                                    }
                                    sum += integrand * mollifier;
                                }
                            }
                        }
                    }
                }
                sum *= (1.0 / (2.0 * l));
                result->get(iX + ofsR.x, iY + ofsR.y, iZ + ofsR.z) = sum;
            }
        }
    }
}

template <typename T>
MollifyScalar3D<T> *MollifyScalar3D<T>::clone() const
{
    return new MollifyScalar3D<T>(*this);
}

template <typename T>
void MollifyScalar3D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Data.
    modified[1] = modif::nothing;          // Flags.
    modified[2] = modif::staticVariables;  // Result.
}

/* ************* Class LBMcomputeGradient3D ******************* */

template <typename T, template <typename U> class Descriptor>
void LBMcomputeGradient3D<T, Descriptor>::process(
    Box3D domain, ScalarField3D<T> &scalarField, TensorField3D<T, 3> &gradient)
{
    typedef Descriptor<T> D;
    Dot3D ofs = computeRelativeDisplacement(scalarField, gradient);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, 3> &gradientOfScalar = gradient.get(iX + ofs.x, iY + ofs.y, iZ + ofs.z);
                gradientOfScalar.resetToZero();
                // Compute the gradient of a scalar function "the lattice Boltzmann way".
                for (plint iPop = 1; iPop < D::q; ++iPop) {
                    plint nextX = iX + D::c[iPop][0];
                    plint nextY = iY + D::c[iPop][1];
                    plint nextZ = iZ + D::c[iPop][2];

                    gradientOfScalar[0] +=
                        D::t[iPop] * D::c[iPop][0] * scalarField.get(nextX, nextY, nextZ);
                    gradientOfScalar[1] +=
                        D::t[iPop] * D::c[iPop][1] * scalarField.get(nextX, nextY, nextZ);
                    gradientOfScalar[2] +=
                        D::t[iPop] * D::c[iPop][2] * scalarField.get(nextX, nextY, nextZ);
                }
                gradientOfScalar *= D::invCs2;
            }
        }
    }
}

/* ************* Class UpdateMinScalarTransientStatistics3D ******************* */

template <typename T>
void UpdateMinScalarTransientStatistics3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 2);
    ScalarField3D<T> *scalar = dynamic_cast<ScalarField3D<T> *>(blocks[0]);
    ScalarField3D<T> *statistics = dynamic_cast<ScalarField3D<T> *>(blocks[1]);
    PLB_ASSERT(scalar);
    PLB_ASSERT(statistics);

    Dot3D offset = computeRelativeDisplacement(*scalar, *statistics);

    for (plint iX = domain.x0; iX <= domain.x1; iX++) {
        for (plint iY = domain.y0; iY <= domain.y1; iY++) {
            for (plint iZ = domain.z0; iZ <= domain.z1; iZ++) {
                statistics->get(iX + offset.x, iY + offset.y, iZ + offset.z) = std::min(
                    statistics->get(iX + offset.x, iY + offset.y, iZ + offset.z),
                    scalar->get(iX, iY, iZ));
            }
        }
    }
}

template <typename T>
UpdateMinScalarTransientStatistics3D<T> *UpdateMinScalarTransientStatistics3D<T>::clone() const
{
    return new UpdateMinScalarTransientStatistics3D<T>(*this);
}

template <typename T>
void UpdateMinScalarTransientStatistics3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Scalar field.
    modified[1] = modif::staticVariables;  // Statistics field.
}

template <typename T>
BlockDomain::DomainT UpdateMinScalarTransientStatistics3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ************* Class UpdateMaxScalarTransientStatistics3D ******************* */

template <typename T>
void UpdateMaxScalarTransientStatistics3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 2);
    ScalarField3D<T> *scalar = dynamic_cast<ScalarField3D<T> *>(blocks[0]);
    ScalarField3D<T> *statistics = dynamic_cast<ScalarField3D<T> *>(blocks[1]);
    PLB_ASSERT(scalar);
    PLB_ASSERT(statistics);

    Dot3D offset = computeRelativeDisplacement(*scalar, *statistics);

    for (plint iX = domain.x0; iX <= domain.x1; iX++) {
        for (plint iY = domain.y0; iY <= domain.y1; iY++) {
            for (plint iZ = domain.z0; iZ <= domain.z1; iZ++) {
                statistics->get(iX + offset.x, iY + offset.y, iZ + offset.z) = std::max(
                    statistics->get(iX + offset.x, iY + offset.y, iZ + offset.z),
                    scalar->get(iX, iY, iZ));
            }
        }
    }
}

template <typename T>
UpdateMaxScalarTransientStatistics3D<T> *UpdateMaxScalarTransientStatistics3D<T>::clone() const
{
    return new UpdateMaxScalarTransientStatistics3D<T>(*this);
}

template <typename T>
void UpdateMaxScalarTransientStatistics3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Scalar field.
    modified[1] = modif::staticVariables;  // Statistics field.
}

template <typename T>
BlockDomain::DomainT UpdateMaxScalarTransientStatistics3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ************* Class UpdateAveScalarTransientStatistics3D ******************* */

template <typename T>
UpdateAveScalarTransientStatistics3D<T>::UpdateAveScalarTransientStatistics3D(plint n_) : n(n_)
{ }

template <typename T>
void UpdateAveScalarTransientStatistics3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 2);
    ScalarField3D<T> *scalar = dynamic_cast<ScalarField3D<T> *>(blocks[0]);
    ScalarField3D<T> *statistics = dynamic_cast<ScalarField3D<T> *>(blocks[1]);
    PLB_ASSERT(scalar);
    PLB_ASSERT(statistics);

    Dot3D offset = computeRelativeDisplacement(*scalar, *statistics);

    T nMinusOne = (T)n - (T)1;
    T oneOverN = (T)1 / (T)n;

    for (plint iX = domain.x0; iX <= domain.x1; iX++) {
        for (plint iY = domain.y0; iY <= domain.y1; iY++) {
            for (plint iZ = domain.z0; iZ <= domain.z1; iZ++) {
                statistics->get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    oneOverN
                    * (nMinusOne * statistics->get(iX + offset.x, iY + offset.y, iZ + offset.z)
                       + scalar->get(iX, iY, iZ));
            }
        }
    }
    ++n;
}

template <typename T>
UpdateAveScalarTransientStatistics3D<T> *UpdateAveScalarTransientStatistics3D<T>::clone() const
{
    return new UpdateAveScalarTransientStatistics3D<T>(*this);
}

template <typename T>
void UpdateAveScalarTransientStatistics3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Scalar field.
    modified[1] = modif::staticVariables;  // Statistics field.
}

template <typename T>
BlockDomain::DomainT UpdateAveScalarTransientStatistics3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ************* Class UpdateRmsScalarTransientStatistics3D ******************* */

template <typename T>
UpdateRmsScalarTransientStatistics3D<T>::UpdateRmsScalarTransientStatistics3D(plint n_) : n(n_)
{ }

template <typename T>
void UpdateRmsScalarTransientStatistics3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 2);
    ScalarField3D<T> *scalar = dynamic_cast<ScalarField3D<T> *>(blocks[0]);
    ScalarField3D<T> *statistics = dynamic_cast<ScalarField3D<T> *>(blocks[1]);
    PLB_ASSERT(scalar);
    PLB_ASSERT(statistics);

    Dot3D offset = computeRelativeDisplacement(*scalar, *statistics);

    T nMinusOne = (T)n - (T)1;
    T oneOverN = (T)1 / (T)n;

    for (plint iX = domain.x0; iX <= domain.x1; iX++) {
        for (plint iY = domain.y0; iY <= domain.y1; iY++) {
            for (plint iZ = domain.z0; iZ <= domain.z1; iZ++) {
                T oldRms = statistics->get(iX + offset.x, iY + offset.y, iZ + offset.z);
                T newData = scalar->get(iX, iY, iZ);
                statistics->get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    std::sqrt(oneOverN * (nMinusOne * oldRms * oldRms + newData * newData));
            }
        }
    }
    ++n;
}

template <typename T>
UpdateRmsScalarTransientStatistics3D<T> *UpdateRmsScalarTransientStatistics3D<T>::clone() const
{
    return new UpdateRmsScalarTransientStatistics3D<T>(*this);
}

template <typename T>
void UpdateRmsScalarTransientStatistics3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Scalar field.
    modified[1] = modif::staticVariables;  // Statistics field.
}

template <typename T>
BlockDomain::DomainT UpdateRmsScalarTransientStatistics3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ************* Class UpdateDevScalarTransientStatistics3D ******************* */

template <typename T>
UpdateDevScalarTransientStatistics3D<T>::UpdateDevScalarTransientStatistics3D(plint n_) : n(n_)
{ }

template <typename T>
void UpdateDevScalarTransientStatistics3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 3);
    ScalarField3D<T> *scalar = dynamic_cast<ScalarField3D<T> *>(blocks[0]);
    ScalarField3D<T> *average = dynamic_cast<ScalarField3D<T> *>(blocks[1]);
    ScalarField3D<T> *statistics = dynamic_cast<ScalarField3D<T> *>(blocks[2]);
    PLB_ASSERT(scalar);
    PLB_ASSERT(average);
    PLB_ASSERT(statistics);

    Dot3D ofsA = computeRelativeDisplacement(*scalar, *average);
    Dot3D ofsS = computeRelativeDisplacement(*scalar, *statistics);

    T nMinusOne = (T)n - (T)1;
    T oneOverN = (T)1 / (T)n;

    for (plint iX = domain.x0; iX <= domain.x1; iX++) {
        for (plint iY = domain.y0; iY <= domain.y1; iY++) {
            for (plint iZ = domain.z0; iZ <= domain.z1; iZ++) {
                T oldDev = statistics->get(iX + ofsS.x, iY + ofsS.y, iZ + ofsS.z);
                T newData =
                    scalar->get(iX, iY, iZ) - average->get(iX + ofsA.x, iY + ofsA.y, iZ + ofsA.z);
                statistics->get(iX + ofsS.x, iY + ofsS.y, iZ + ofsS.z) =
                    std::sqrt(oneOverN * (nMinusOne * oldDev * oldDev + newData * newData));
            }
        }
    }
    ++n;
}

template <typename T>
UpdateDevScalarTransientStatistics3D<T> *UpdateDevScalarTransientStatistics3D<T>::clone() const
{
    return new UpdateDevScalarTransientStatistics3D<T>(*this);
}

template <typename T>
void UpdateDevScalarTransientStatistics3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Scalar field.
    modified[1] = modif::nothing;          // Mean value field.
    modified[2] = modif::staticVariables;  // Statistics field.
}

template <typename T>
BlockDomain::DomainT UpdateDevScalarTransientStatistics3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* *************** PART III ****************************************** */
/* *************** Analysis of the tensor-field ********************** */
/* ******************************************************************* */

template <typename T, int nDim>
BoxTensorSumFunctional3D<T, nDim>::BoxTensorSumFunctional3D()
{
    for (plint i = 0; i < nDim; i++) {
        sumTensorId[i] = this->getStatistics().subscribeSum();
    }
}

template <typename T, int nDim>
void BoxTensorSumFunctional3D<T, nDim>::process(Box3D domain, TensorField3D<T, nDim> &tensorField)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                for (plint i = 0; i < nDim; i++) {
                    statistics.gatherSum(sumTensorId[i], (double)tensorField.get(iX, iY, iZ)[i]);
                }
            }
        }
    }
}

template <typename T, int nDim>
BoxTensorSumFunctional3D<T, nDim> *BoxTensorSumFunctional3D<T, nDim>::clone() const
{
    return new BoxTensorSumFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
Array<T, nDim> BoxTensorSumFunctional3D<T, nDim>::getSumTensor() const
{
    Array<T, nDim> sum;
    for (plint i = 0; i < nDim; i++) {
        double doubleSum = this->getStatistics().getSum(sumTensorId[i]);
        // The sum is internally computed on floating-point values. If T is
        //   integer, the value must be rounded at the end.
        if (std::numeric_limits<T>::is_integer) {
            doubleSum = util::roundToInt(doubleSum);
        }
        sum[i] = doubleSum;
    }
    return sum;
}

template <typename T, int nDim>
MaskedBoxTensorSumFunctional3D<T, nDim>::MaskedBoxTensorSumFunctional3D(int flag_) : flag(flag_)
{
    for (plint i = 0; i < nDim; i++) {
        sumTensorId[i] = this->getStatistics().subscribeSum();
    }
}

template <typename T, int nDim>
void MaskedBoxTensorSumFunctional3D<T, nDim>::process(
    Box3D domain, ScalarField3D<int> &mask, TensorField3D<T, nDim> &tensorField)
{
    Dot3D offset = computeRelativeDisplacement(mask, tensorField);
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (mask.get(iX, iY, iZ) == flag) {
                    for (plint i = 0; i < nDim; i++) {
                        statistics.gatherSum(
                            sumTensorId[i], (double)tensorField.get(
                                                iX + offset.x, iY + offset.y, iZ + offset.z)[i]);
                    }
                }
            }
        }
    }
}

template <typename T, int nDim>
MaskedBoxTensorSumFunctional3D<T, nDim> *MaskedBoxTensorSumFunctional3D<T, nDim>::clone() const
{
    return new MaskedBoxTensorSumFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
Array<T, nDim> MaskedBoxTensorSumFunctional3D<T, nDim>::getSumTensor() const
{
    Array<T, nDim> sum;
    for (plint i = 0; i < nDim; i++) {
        double doubleSum = this->getStatistics().getSum(sumTensorId[i]);
        // The sum is internally computed on floating-point values. If T is
        //   integer, the value must be rounded at the end.
        if (std::numeric_limits<T>::is_integer) {
            doubleSum = util::roundToInt(doubleSum);
        }
        sum[i] = doubleSum;
    }
    return sum;
}

template <typename T, int nDim>
MaskedBoxTensorAverageFunctional3D<T, nDim>::MaskedBoxTensorAverageFunctional3D(int flag_) :
    flag(flag_)
{
    for (plint i = 0; i < nDim; i++) {
        averageTensorId[i] = this->getStatistics().subscribeAverage();
    }
}

template <typename T, int nDim>
void MaskedBoxTensorAverageFunctional3D<T, nDim>::process(
    Box3D domain, ScalarField3D<int> &mask, TensorField3D<T, nDim> &tensorField)
{
    Dot3D offset = computeRelativeDisplacement(mask, tensorField);
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (mask.get(iX, iY, iZ) == flag) {
                    for (plint i = 0; i < nDim; i++) {
                        statistics.gatherAverage(
                            averageTensorId[i],
                            (double)tensorField.get(
                                iX + offset.x, iY + offset.y, iZ + offset.z)[i]);
                    }
                    statistics.incrementStats();
                }
            }
        }
    }
}

template <typename T, int nDim>
MaskedBoxTensorAverageFunctional3D<T, nDim> *MaskedBoxTensorAverageFunctional3D<T, nDim>::clone()
    const
{
    return new MaskedBoxTensorAverageFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
Array<T, nDim> MaskedBoxTensorAverageFunctional3D<T, nDim>::getAverageTensor() const
{
    Array<T, nDim> average;
    for (plint i = 0; i < nDim; i++) {
        double doubleAverage = this->getStatistics().getAverage(averageTensorId[i]);
        // The average is internally computed on floating-point values. If T is
        //   integer, the value must be rounded at the end.
        if (std::numeric_limits<T>::is_integer) {
            doubleAverage = util::roundToInt(doubleAverage);
        }
        average[i] = doubleAverage;
    }
    return average;
}

template <typename T1, typename T2, int nDim>
void CopyConvertTensorFunctional3D<T1, T2, nDim>::process(
    Box3D domain, TensorField3D<T1, nDim> &field1, TensorField3D<T2, nDim> &field2)
{
    Dot3D offset = computeRelativeDisplacement(field1, field2);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                for (int iDim = 0; iDim < nDim; ++iDim) {
                    field2.get(iX + offset.x, iY + offset.y, iZ + offset.z)[iDim] =
                        (T2)field1.get(iX, iY, iZ)[iDim];
                }
            }
        }
    }
}

template <typename T1, typename T2, int nDim>
CopyConvertTensorFunctional3D<T1, T2, nDim> *CopyConvertTensorFunctional3D<T1, T2, nDim>::clone()
    const
{
    return new CopyConvertTensorFunctional3D<T1, T2, nDim>(*this);
}

template <typename T1, typename T2, int nDim>
void CopyConvertTensorFunctional3D<T1, T2, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T1, typename T2, int nDim>
BlockDomain::DomainT CopyConvertTensorFunctional3D<T1, T2, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void ExtractTensorSubDomainFunctional3D<T, nDim>::process(
    Box3D domain, TensorField3D<T, nDim> &field1, TensorField3D<T, nDim> &field2)
{
    Dot3D offset = computeRelativeDisplacement(field1, field2);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                for (int iDim = 0; iDim < nDim; ++iDim) {
                    field2.get(iX + offset.x, iY + offset.y, iZ + offset.z)[iDim] =
                        field1.get(iX, iY, iZ)[iDim];
                }
            }
        }
    }
}

template <typename T, int nDim>
ExtractTensorSubDomainFunctional3D<T, nDim> *ExtractTensorSubDomainFunctional3D<T, nDim>::clone()
    const
{
    return new ExtractTensorSubDomainFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void ExtractTensorSubDomainFunctional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT ExtractTensorSubDomainFunctional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, int nDim>
ExtractTensorComponentFunctional3D<T, nDim>::ExtractTensorComponentFunctional3D(int iComponent_) :
    iComponent(iComponent_)
{
    PLB_ASSERT(iComponent < nDim);
}

template <typename T, int nDim>
void ExtractTensorComponentFunctional3D<T, nDim>::process(
    Box3D domain, ScalarField3D<T> &scalarField, TensorField3D<T, nDim> &tensorField)
{
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                scalarField.get(iX, iY, iZ) =
                    tensorField.get(iX + offset.x, iY + offset.y, iZ + offset.z)[iComponent];
            }
        }
    }
}

template <typename T, int nDim>
ExtractTensorComponentFunctional3D<T, nDim> *ExtractTensorComponentFunctional3D<T, nDim>::clone()
    const
{
    return new ExtractTensorComponentFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void ExtractTensorComponentFunctional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, int nDim>
BlockDomain::DomainT ExtractTensorComponentFunctional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void ComputeNormFunctional3D<T, nDim>::process(
    Box3D domain, ScalarField3D<T> &scalarField, TensorField3D<T, nDim> &tensorField)
{
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                scalarField.get(iX, iY, iZ) = std::sqrt(VectorTemplateImpl<T, nDim>::normSqr(
                    tensorField.get(iX + offset.x, iY + offset.y, iZ + offset.z)));
            }
        }
    }
}

template <typename T, int nDim>
ComputeNormFunctional3D<T, nDim> *ComputeNormFunctional3D<T, nDim>::clone() const
{
    return new ComputeNormFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void ComputeNormFunctional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, int nDim>
BlockDomain::DomainT ComputeNormFunctional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void ComputeNormSqrFunctional3D<T, nDim>::process(
    Box3D domain, ScalarField3D<T> &scalarField, TensorField3D<T, nDim> &tensorField)
{
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                scalarField.get(iX, iY, iZ) = VectorTemplateImpl<T, nDim>::normSqr(
                    tensorField.get(iX + offset.x, iY + offset.y, iZ + offset.z));
            }
        }
    }
}

template <typename T, int nDim>
ComputeNormSqrFunctional3D<T, nDim> *ComputeNormSqrFunctional3D<T, nDim>::clone() const
{
    return new ComputeNormSqrFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void ComputeNormSqrFunctional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, int nDim>
BlockDomain::DomainT ComputeNormSqrFunctional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

// ============= BoxLocalMaximumPerComponentFunctional3D ============== //

template <typename T, int nDim>
void BoxLocalMaximumPerComponentFunctional3D<T, nDim>::process(
    Box3D domain, ScalarField3D<T> &scalarField, TensorField3D<T, nDim> &tensorField)
{
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                scalarField.get(iX, iY, iZ) =
                    maxElement(tensorField.get(iX + offset.x, iY + offset.y, iZ + offset.z));
            }
        }
    }
}

template <typename T, int nDim>
BoxLocalMaximumPerComponentFunctional3D<T, nDim>
    *BoxLocalMaximumPerComponentFunctional3D<T, nDim>::clone() const
{
    return new BoxLocalMaximumPerComponentFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void BoxLocalMaximumPerComponentFunctional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T>
void ComputeSymmetricTensorNormFunctional3D<T>::process(
    Box3D domain, ScalarField3D<T> &scalarField, TensorField3D<T, 6> &tensorField)
{
    typedef SymmetricTensorImpl<T, 3> tensor;
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, 6> &el = tensorField.get(iX + offset.x, iY + offset.y, iZ + offset.z);
                scalarField.get(iX, iY, iZ) = std::sqrt(
                    // Count diagonal components once ...
                    util::sqr(el[tensor::xx]) + util::sqr(el[tensor::yy])
                    + util::sqr(el[tensor::zz]) +
                    // .. and off-diagonal components twice, due to symmetry.
                    (T)2
                        * (util::sqr(el[tensor::xy]) + util::sqr(el[tensor::xz])
                           + util::sqr(el[tensor::yz])));
            }
        }
    }
}

template <typename T>
ComputeSymmetricTensorNormFunctional3D<T> *ComputeSymmetricTensorNormFunctional3D<T>::clone() const
{
    return new ComputeSymmetricTensorNormFunctional3D<T>(*this);
}

template <typename T>
void ComputeSymmetricTensorNormFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T>
BlockDomain::DomainT ComputeSymmetricTensorNormFunctional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T>
void ComputeSymmetricTensorNormSqrFunctional3D<T>::process(
    Box3D domain, ScalarField3D<T> &scalarField, TensorField3D<T, 6> &tensorField)
{
    typedef SymmetricTensorImpl<T, 3> tensor;
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, 6> &el = tensorField.get(iX + offset.x, iY + offset.y, iZ + offset.z);
                scalarField.get(iX, iY, iZ) =
                    // Count diagonal components once ...
                    util::sqr(el[tensor::xx]) + util::sqr(el[tensor::yy])
                    + util::sqr(el[tensor::zz]) +
                    // .. and off-diagonal components twice, due to symmetry.
                    (T)2
                        * (util::sqr(el[tensor::xy]) + util::sqr(el[tensor::xz])
                           + util::sqr(el[tensor::yz]));
            }
        }
    }
}

template <typename T>
ComputeSymmetricTensorNormSqrFunctional3D<T> *ComputeSymmetricTensorNormSqrFunctional3D<T>::clone()
    const
{
    return new ComputeSymmetricTensorNormSqrFunctional3D<T>(*this);
}

template <typename T>
void ComputeSymmetricTensorNormSqrFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T>
BlockDomain::DomainT ComputeSymmetricTensorNormSqrFunctional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T>
void ComputeSymmetricTensorTraceFunctional3D<T>::process(
    Box3D domain, ScalarField3D<T> &scalarField, TensorField3D<T, 6> &tensorField)
{
    typedef SymmetricTensorImpl<T, 3> tensor;
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, 6> &el = tensorField.get(iX + offset.x, iY + offset.y, iZ + offset.z);
                scalarField.get(iX, iY, iZ) = el[tensor::xx] + el[tensor::yy] + el[tensor::zz];
            }
        }
    }
}

template <typename T>
ComputeSymmetricTensorTraceFunctional3D<T> *ComputeSymmetricTensorTraceFunctional3D<T>::clone()
    const
{
    return new ComputeSymmetricTensorTraceFunctional3D<T>(*this);
}

template <typename T>
void ComputeSymmetricTensorTraceFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T>
BlockDomain::DomainT ComputeSymmetricTensorTraceFunctional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T>
void BoxBulkGradientFunctional3D<T>::process(
    Box3D domain, ScalarField3D<T> &phi, TensorField3D<T, 3> &gradient)
{
    Dot3D offset = computeRelativeDisplacement(phi, gradient);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iX2 = iX + offset.x;
                plint iY2 = iY + offset.y;
                plint iZ2 = iZ + offset.z;
                gradient.get(iX2, iY2, iZ2)[0] = fdDataField::bulkXderiv(phi, iX, iY, iZ);
                gradient.get(iX2, iY2, iZ2)[1] = fdDataField::bulkYderiv(phi, iX, iY, iZ);
                gradient.get(iX2, iY2, iZ2)[2] = fdDataField::bulkZderiv(phi, iX, iY, iZ);
            }
        }
    }
}

template <typename T>
BoxBulkGradientFunctional3D<T> *BoxBulkGradientFunctional3D<T>::clone() const
{
    return new BoxBulkGradientFunctional3D<T>(*this);
}

template <typename T>
void BoxBulkGradientFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT BoxBulkGradientFunctional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T>
void BoxGradientFunctional3D<T>::processBulk(
    Box3D domain, ScalarField3D<T> &phi, TensorField3D<T, 3> &gradient)
{
    Dot3D offset = computeRelativeDisplacement(phi, gradient);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iX2 = iX + offset.x;
                plint iY2 = iY + offset.y;
                plint iZ2 = iZ + offset.z;
                gradient.get(iX2, iY2, iZ2)[0] = fdDataField::bulkXderiv(phi, iX, iY, iZ);
                gradient.get(iX2, iY2, iZ2)[1] = fdDataField::bulkYderiv(phi, iX, iY, iZ);
                gradient.get(iX2, iY2, iZ2)[2] = fdDataField::bulkZderiv(phi, iX, iY, iZ);
            }
        }
    }
}

template <typename T>
void BoxGradientFunctional3D<T>::processPlane(
    int direction, int orientation, Box3D domain, ScalarField3D<T> &phi,
    TensorField3D<T, 3> &gradient)
{
    Dot3D offset = computeRelativeDisplacement(phi, gradient);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iX2 = iX + offset.x;
                plint iY2 = iY + offset.y;
                plint iZ2 = iZ + offset.z;
                gradient.get(iX2, iY2, iZ2)[0] =
                    fdDataField::planeXderiv(phi, direction, orientation, iX, iY, iZ);
                gradient.get(iX2, iY2, iZ2)[1] =
                    fdDataField::planeYderiv(phi, direction, orientation, iX, iY, iZ);
                gradient.get(iX2, iY2, iZ2)[2] =
                    fdDataField::planeZderiv(phi, direction, orientation, iX, iY, iZ);
            }
        }
    }
}

template <typename T>
void BoxGradientFunctional3D<T>::processEdge(
    int plane, int normal1, int normal2, Box3D domain, ScalarField3D<T> &phi,
    TensorField3D<T, 3> &gradient)
{
    Dot3D offset = computeRelativeDisplacement(phi, gradient);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iX2 = iX + offset.x;
                plint iY2 = iY + offset.y;
                plint iZ2 = iZ + offset.z;
                gradient.get(iX2, iY2, iZ2)[0] =
                    fdDataField::edgeXderiv(phi, plane, normal1, normal2, iX, iY, iZ);
                gradient.get(iX2, iY2, iZ2)[1] =
                    fdDataField::edgeYderiv(phi, plane, normal1, normal2, iX, iY, iZ);
                gradient.get(iX2, iY2, iZ2)[2] =
                    fdDataField::edgeZderiv(phi, plane, normal1, normal2, iX, iY, iZ);
            }
        }
    }
}

template <typename T>
void BoxGradientFunctional3D<T>::processCorner(
    int normalX, int normalY, int normalZ, Box3D domain, ScalarField3D<T> &phi,
    TensorField3D<T, 3> &gradient)
{
    Dot3D offset = computeRelativeDisplacement(phi, gradient);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iX2 = iX + offset.x;
                plint iY2 = iY + offset.y;
                plint iZ2 = iZ + offset.z;
                gradient.get(iX2, iY2, iZ2)[0] =
                    fdDataField::cornerXderiv(phi, normalX, normalY, normalZ, iX, iY, iZ);
                gradient.get(iX2, iY2, iZ2)[1] =
                    fdDataField::cornerYderiv(phi, normalX, normalY, normalZ, iX, iY, iZ);
                gradient.get(iX2, iY2, iZ2)[2] =
                    fdDataField::cornerZderiv(phi, normalX, normalY, normalZ, iX, iY, iZ);
            }
        }
    }
}

template <typename T>
BoxGradientFunctional3D<T> *BoxGradientFunctional3D<T>::clone() const
{
    return new BoxGradientFunctional3D<T>(*this);
}

template <typename T>
void BoxGradientFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT BoxGradientFunctional3D<T>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void BoxBulkVorticityFunctional3D<T, nDim>::process(
    Box3D domain, TensorField3D<T, nDim> &velocity, TensorField3D<T, nDim> &vorticity)
{
    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint iX2 = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint iY2 = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iZ2 = iZ + offset.z;
                vorticity.get(iX2, iY2, iZ2)[0] = fdDataField::bulkVorticityX(velocity, iX, iY, iZ);
                vorticity.get(iX2, iY2, iZ2)[1] = fdDataField::bulkVorticityY(velocity, iX, iY, iZ);
                vorticity.get(iX2, iY2, iZ2)[2] = fdDataField::bulkVorticityZ(velocity, iX, iY, iZ);
            }
        }
    }
}

template <typename T, int nDim>
BoxBulkVorticityFunctional3D<T, nDim> *BoxBulkVorticityFunctional3D<T, nDim>::clone() const
{
    return new BoxBulkVorticityFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void BoxBulkVorticityFunctional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT BoxBulkVorticityFunctional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

//===================  BoxBulkVorticityOrderFourFunctional3D ==================//

template <typename T, int nDim>
void BoxBulkVorticityOrderFourFunctional3D<T, nDim>::process(
    Box3D domain, TensorField3D<T, nDim> &velocity, TensorField3D<T, nDim> &vorticity)
{
    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint iX2 = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint iY2 = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iZ2 = iZ + offset.z;
                vorticity.get(iX2, iY2, iZ2)[0] =
                    fdDataField::bulkVorticityXOrderFour(velocity, iX, iY, iZ);
                vorticity.get(iX2, iY2, iZ2)[1] =
                    fdDataField::bulkVorticityYOrderFour(velocity, iX, iY, iZ);
                vorticity.get(iX2, iY2, iZ2)[2] =
                    fdDataField::bulkVorticityZOrderFour(velocity, iX, iY, iZ);
            }
        }
    }
}

template <typename T, int nDim>
BoxBulkVorticityOrderFourFunctional3D<T, nDim>
    *BoxBulkVorticityOrderFourFunctional3D<T, nDim>::clone() const
{
    return new BoxBulkVorticityOrderFourFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void BoxBulkVorticityOrderFourFunctional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT BoxBulkVorticityOrderFourFunctional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

//===================  BoxBulkVorticityOrderSixFunctional3D ==================//

template <typename T, int nDim>
void BoxBulkVorticityOrderSixFunctional3D<T, nDim>::process(
    Box3D domain, TensorField3D<T, nDim> &velocity, TensorField3D<T, nDim> &vorticity)
{
    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint iX2 = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint iY2 = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iZ2 = iZ + offset.z;
                vorticity.get(iX2, iY2, iZ2)[0] =
                    fdDataField::bulkVorticityXOrderSix(velocity, iX, iY, iZ);
                vorticity.get(iX2, iY2, iZ2)[1] =
                    fdDataField::bulkVorticityYOrderSix(velocity, iX, iY, iZ);
                vorticity.get(iX2, iY2, iZ2)[2] =
                    fdDataField::bulkVorticityZOrderSix(velocity, iX, iY, iZ);
            }
        }
    }
}

template <typename T, int nDim>
BoxBulkVorticityOrderSixFunctional3D<T, nDim>
    *BoxBulkVorticityOrderSixFunctional3D<T, nDim>::clone() const
{
    return new BoxBulkVorticityOrderSixFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void BoxBulkVorticityOrderSixFunctional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT BoxBulkVorticityOrderSixFunctional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

//===================  BoxBulkVorticityOrderEightFunctional3D ==================//

template <typename T, int nDim>
void BoxBulkVorticityOrderEightFunctional3D<T, nDim>::process(
    Box3D domain, TensorField3D<T, nDim> &velocity, TensorField3D<T, nDim> &vorticity)
{
    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint iX2 = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint iY2 = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iZ2 = iZ + offset.z;
                vorticity.get(iX2, iY2, iZ2)[0] =
                    fdDataField::bulkVorticityXOrderEight(velocity, iX, iY, iZ);
                vorticity.get(iX2, iY2, iZ2)[1] =
                    fdDataField::bulkVorticityYOrderEight(velocity, iX, iY, iZ);
                vorticity.get(iX2, iY2, iZ2)[2] =
                    fdDataField::bulkVorticityZOrderEight(velocity, iX, iY, iZ);
            }
        }
    }
}

template <typename T, int nDim>
BoxBulkVorticityOrderEightFunctional3D<T, nDim>
    *BoxBulkVorticityOrderEightFunctional3D<T, nDim>::clone() const
{
    return new BoxBulkVorticityOrderEightFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void BoxBulkVorticityOrderEightFunctional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT BoxBulkVorticityOrderEightFunctional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void BoxVorticityFunctional3D<T, nDim>::processBulk(
    Box3D domain, TensorField3D<T, nDim> &velocity, TensorField3D<T, nDim> &vorticity)
{
    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iX2 = iX + offset.x;
                plint iY2 = iY + offset.y;
                plint iZ2 = iZ + offset.z;
                vorticity.get(iX2, iY2, iZ2)[0] = fdDataField::bulkVorticityX(velocity, iX, iY, iZ);
                vorticity.get(iX2, iY2, iZ2)[1] = fdDataField::bulkVorticityY(velocity, iX, iY, iZ);
                vorticity.get(iX2, iY2, iZ2)[2] = fdDataField::bulkVorticityZ(velocity, iX, iY, iZ);
            }
        }
    }
}

template <typename T, int nDim>
void BoxVorticityFunctional3D<T, nDim>::processPlane(
    int direction, int orientation, Box3D domain, TensorField3D<T, nDim> &velocity,
    TensorField3D<T, nDim> &vorticity)
{
    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iX2 = iX + offset.x;
                plint iY2 = iY + offset.y;
                plint iZ2 = iZ + offset.z;
                vorticity.get(iX2, iY2, iZ2)[0] =
                    fdDataField::planeVorticityX(velocity, direction, orientation, iX, iY, iZ);
                vorticity.get(iX2, iY2, iZ2)[1] =
                    fdDataField::planeVorticityY(velocity, direction, orientation, iX, iY, iZ);
                vorticity.get(iX2, iY2, iZ2)[2] =
                    fdDataField::planeVorticityZ(velocity, direction, orientation, iX, iY, iZ);
            }
        }
    }
}

template <typename T, int nDim>
void BoxVorticityFunctional3D<T, nDim>::processEdge(
    int plane, int normal1, int normal2, Box3D domain, TensorField3D<T, nDim> &velocity,
    TensorField3D<T, nDim> &vorticity)
{
    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iX2 = iX + offset.x;
                plint iY2 = iY + offset.y;
                plint iZ2 = iZ + offset.z;
                vorticity.get(iX2, iY2, iZ2)[0] =
                    fdDataField::edgeVorticityX(velocity, plane, normal1, normal2, iX, iY, iZ);
                vorticity.get(iX2, iY2, iZ2)[1] =
                    fdDataField::edgeVorticityY(velocity, plane, normal1, normal2, iX, iY, iZ);
                vorticity.get(iX2, iY2, iZ2)[2] =
                    fdDataField::edgeVorticityZ(velocity, plane, normal1, normal2, iX, iY, iZ);
            }
        }
    }
}

template <typename T, int nDim>
void BoxVorticityFunctional3D<T, nDim>::processCorner(
    int normalX, int normalY, int normalZ, Box3D domain, TensorField3D<T, nDim> &velocity,
    TensorField3D<T, nDim> &vorticity)
{
    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iX2 = iX + offset.x;
                plint iY2 = iY + offset.y;
                plint iZ2 = iZ + offset.z;
                vorticity.get(iX2, iY2, iZ2)[0] =
                    fdDataField::cornerVorticityX(velocity, normalX, normalY, normalZ, iX, iY, iZ);
                vorticity.get(iX2, iY2, iZ2)[1] =
                    fdDataField::cornerVorticityY(velocity, normalX, normalY, normalZ, iX, iY, iZ);
                vorticity.get(iX2, iY2, iZ2)[2] =
                    fdDataField::cornerVorticityZ(velocity, normalX, normalY, normalZ, iX, iY, iZ);
            }
        }
    }
}

template <typename T, int nDim>
BoxVorticityFunctional3D<T, nDim> *BoxVorticityFunctional3D<T, nDim>::clone() const
{
    return new BoxVorticityFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void BoxVorticityFunctional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT BoxVorticityFunctional3D<T, nDim>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void BoxBulkHelicityFunctional3D<T, nDim>::process(
    Box3D domain, ScalarField3D<T> &helicity, TensorField3D<T, nDim> &velocity)
{
    Dot3D offset = computeRelativeDisplacement(helicity, velocity);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint iX2 = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint iY2 = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iZ2 = iZ + offset.z;
                helicity.get(iX, iY, iZ) =
                    velocity.get(iX2, iY2, iZ2)[0]
                        * fdDataField::bulkVorticityX(velocity, iX2, iY2, iZ2)
                    + velocity.get(iX2, iY2, iZ2)[1]
                          * fdDataField::bulkVorticityY(velocity, iX2, iY2, iZ2)
                    + velocity.get(iX2, iY2, iZ2)[2]
                          * fdDataField::bulkVorticityZ(velocity, iX2, iY2, iZ2);
            }
        }
    }
}

template <typename T, int nDim>
BoxBulkHelicityFunctional3D<T, nDim> *BoxBulkHelicityFunctional3D<T, nDim>::clone() const
{
    return new BoxBulkHelicityFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void BoxBulkHelicityFunctional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, int nDim>
BlockDomain::DomainT BoxBulkHelicityFunctional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void BoxHelicityFunctional3D<T, nDim>::processBulk(
    Box3D domain, ScalarField3D<T> &helicity, TensorField3D<T, nDim> &velocity)
{
    Dot3D offset = computeRelativeDisplacement(helicity, velocity);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iX2 = iX + offset.x;
                plint iY2 = iY + offset.y;
                plint iZ2 = iZ + offset.z;
                helicity.get(iX, iY, iZ) =
                    velocity.get(iX2, iY2, iZ2)[0]
                        * fdDataField::bulkVorticityX(velocity, iX2, iY2, iZ2)
                    + velocity.get(iX2, iY2, iZ2)[1]
                          * fdDataField::bulkVorticityY(velocity, iX2, iY2, iZ2)
                    + velocity.get(iX2, iY2, iZ2)[2]
                          * fdDataField::bulkVorticityZ(velocity, iX2, iY2, iZ2);
            }
        }
    }
}

template <typename T, int nDim>
void BoxHelicityFunctional3D<T, nDim>::processPlane(
    int direction, int orientation, Box3D domain, ScalarField3D<T> &helicity,
    TensorField3D<T, nDim> &velocity)
{
    Dot3D offset = computeRelativeDisplacement(helicity, velocity);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iX2 = iX + offset.x;
                plint iY2 = iY + offset.y;
                plint iZ2 = iZ + offset.z;
                helicity.get(iX, iY, iZ) =
                    velocity.get(iX2, iY2, iZ2)[0]
                        * fdDataField::planeVorticityX(
                            velocity, direction, orientation, iX2, iY2, iZ2)
                    + velocity.get(iX2, iY2, iZ2)[1]
                          * fdDataField::planeVorticityY(
                              velocity, direction, orientation, iX2, iY2, iZ2)
                    + velocity.get(iX2, iY2, iZ2)[2]
                          * fdDataField::planeVorticityZ(
                              velocity, direction, orientation, iX2, iY2, iZ2);
            }
        }
    }
}

template <typename T, int nDim>
void BoxHelicityFunctional3D<T, nDim>::processEdge(
    int plane, int normal1, int normal2, Box3D domain, ScalarField3D<T> &helicity,
    TensorField3D<T, nDim> &velocity)
{
    Dot3D offset = computeRelativeDisplacement(helicity, velocity);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iX2 = iX + offset.x;
                plint iY2 = iY + offset.y;
                plint iZ2 = iZ + offset.z;
                helicity.get(iX, iY, iZ) =
                    velocity.get(iX2, iY2, iZ2)[0]
                        * fdDataField::edgeVorticityX(
                            velocity, plane, normal1, normal2, iX2, iY2, iZ2)
                    + velocity.get(iX2, iY2, iZ2)[1]
                          * fdDataField::edgeVorticityY(
                              velocity, plane, normal1, normal2, iX2, iY2, iZ2)
                    + velocity.get(iX2, iY2, iZ2)[2]
                          * fdDataField::edgeVorticityZ(
                              velocity, plane, normal1, normal2, iX2, iY2, iZ2);
            }
        }
    }
}

template <typename T, int nDim>
void BoxHelicityFunctional3D<T, nDim>::processCorner(
    int normalX, int normalY, int normalZ, Box3D domain, ScalarField3D<T> &helicity,
    TensorField3D<T, nDim> &velocity)
{
    Dot3D offset = computeRelativeDisplacement(helicity, velocity);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iX2 = iX + offset.x;
                plint iY2 = iY + offset.y;
                plint iZ2 = iZ + offset.z;
                helicity.get(iX, iY, iZ) =
                    velocity.get(iX2, iY2, iZ2)[0]
                        * fdDataField::cornerVorticityX(
                            velocity, normalX, normalY, normalZ, iX2, iY2, iZ2)
                    + velocity.get(iX2, iY2, iZ2)[1]
                          * fdDataField::cornerVorticityY(
                              velocity, normalX, normalY, normalZ, iX2, iY2, iZ2)
                    + velocity.get(iX2, iY2, iZ2)[2]
                          * fdDataField::cornerVorticityZ(
                              velocity, normalX, normalY, normalZ, iX2, iY2, iZ2);
            }
        }
    }
}

template <typename T, int nDim>
BoxHelicityFunctional3D<T, nDim> *BoxHelicityFunctional3D<T, nDim>::clone() const
{
    return new BoxHelicityFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void BoxHelicityFunctional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, int nDim>
BlockDomain::DomainT BoxHelicityFunctional3D<T, nDim>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void BoxBulkStrainRateFunctional3D<T, nDim>::process(
    Box3D domain, TensorField3D<T, nDim> &velocity,
    TensorField3D<T, SymmetricTensorImpl<T, nDim>::n> &S)
{
    typedef SymmetricTensorImpl<T, nDim> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iX2 = iX + offset.x;
                plint iY2 = iY + offset.y;
                plint iZ2 = iZ + offset.z;
                Array<T, SymmetricTensorImpl<T, nDim>::n> &el = S.get(iX2, iY2, iZ2);
                el[tensor::xx] = fdDataField::bulkXderiv(velocity, iX, iY, iZ, 0);
                el[tensor::xy] = (fdDataField::bulkXderiv(velocity, iX, iY, iZ, 1)
                                  + fdDataField::bulkYderiv(velocity, iX, iY, iZ, 0))
                                 / (T)2;
                el[tensor::xz] = (fdDataField::bulkXderiv(velocity, iX, iY, iZ, 2)
                                  + fdDataField::bulkZderiv(velocity, iX, iY, iZ, 0))
                                 / (T)2;
                el[tensor::yy] = fdDataField::bulkYderiv(velocity, iX, iY, iZ, 1);
                el[tensor::yz] = (fdDataField::bulkYderiv(velocity, iX, iY, iZ, 2)
                                  + fdDataField::bulkZderiv(velocity, iX, iY, iZ, 1))
                                 / (T)2;
                el[tensor::zz] = fdDataField::bulkZderiv(velocity, iX, iY, iZ, 2);
            }
        }
    }
}

template <typename T, int nDim>
BoxBulkStrainRateFunctional3D<T, nDim> *BoxBulkStrainRateFunctional3D<T, nDim>::clone() const
{
    return new BoxBulkStrainRateFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void BoxBulkStrainRateFunctional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT BoxBulkStrainRateFunctional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void BoxStrainRateFunctional3D<T, nDim>::processBulk(
    Box3D domain, TensorField3D<T, nDim> &velocity,
    TensorField3D<T, SymmetricTensorImpl<T, nDim>::n> &S)
{
    typedef SymmetricTensorImpl<T, nDim> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iX2 = iX + offset.x;
                plint iY2 = iY + offset.y;
                plint iZ2 = iZ + offset.z;
                Array<T, SymmetricTensorImpl<T, nDim>::n> &el = S.get(iX2, iY2, iZ2);
                el[tensor::xx] = fdDataField::bulkXderiv(velocity, iX, iY, iZ, 0);
                el[tensor::xy] = (fdDataField::bulkXderiv(velocity, iX, iY, iZ, 1)
                                  + fdDataField::bulkYderiv(velocity, iX, iY, iZ, 0))
                                 / (T)2;
                el[tensor::xz] = (fdDataField::bulkXderiv(velocity, iX, iY, iZ, 2)
                                  + fdDataField::bulkZderiv(velocity, iX, iY, iZ, 0))
                                 / (T)2;
                el[tensor::yy] = fdDataField::bulkYderiv(velocity, iX, iY, iZ, 1);
                el[tensor::yz] = (fdDataField::bulkYderiv(velocity, iX, iY, iZ, 2)
                                  + fdDataField::bulkZderiv(velocity, iX, iY, iZ, 1))
                                 / (T)2;
                el[tensor::zz] = fdDataField::bulkZderiv(velocity, iX, iY, iZ, 2);
            }
        }
    }
}

template <typename T, int nDim>
void BoxStrainRateFunctional3D<T, nDim>::processPlane(
    int direction, int orientation, Box3D domain, TensorField3D<T, nDim> &velocity,
    TensorField3D<T, SymmetricTensorImpl<T, nDim>::n> &S)
{
    typedef SymmetricTensorImpl<T, nDim> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iX2 = iX + offset.x;
                plint iY2 = iY + offset.y;
                plint iZ2 = iZ + offset.z;
                Array<T, SymmetricTensorImpl<T, nDim>::n> &el = S.get(iX2, iY2, iZ2);
                el[tensor::xx] =
                    fdDataField::planeXderiv(velocity, direction, orientation, iX, iY, iZ, 0);
                el[tensor::xy] =
                    (fdDataField::planeXderiv(velocity, direction, orientation, iX, iY, iZ, 1)
                     + fdDataField::planeYderiv(velocity, direction, orientation, iX, iY, iZ, 0))
                    / (T)2;
                el[tensor::xz] =
                    (fdDataField::planeXderiv(velocity, direction, orientation, iX, iY, iZ, 2)
                     + fdDataField::planeZderiv(velocity, direction, orientation, iX, iY, iZ, 0))
                    / (T)2;
                el[tensor::yy] =
                    fdDataField::planeYderiv(velocity, direction, orientation, iX, iY, iZ, 1);
                el[tensor::yz] =
                    (fdDataField::planeYderiv(velocity, direction, orientation, iX, iY, iZ, 2)
                     + fdDataField::planeZderiv(velocity, direction, orientation, iX, iY, iZ, 1))
                    / (T)2;
                el[tensor::zz] =
                    fdDataField::planeZderiv(velocity, direction, orientation, iX, iY, iZ, 2);
            }
        }
    }
}

template <typename T, int nDim>
void BoxStrainRateFunctional3D<T, nDim>::processEdge(
    int plane, int normal1, int normal2, Box3D domain, TensorField3D<T, nDim> &velocity,
    TensorField3D<T, SymmetricTensorImpl<T, nDim>::n> &S)
{
    typedef SymmetricTensorImpl<T, nDim> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iX2 = iX + offset.x;
                plint iY2 = iY + offset.y;
                plint iZ2 = iZ + offset.z;
                Array<T, SymmetricTensorImpl<T, nDim>::n> &el = S.get(iX2, iY2, iZ2);
                el[tensor::xx] =
                    fdDataField::edgeXderiv(velocity, plane, normal1, normal2, iX, iY, iZ, 0);
                el[tensor::xy] =
                    (fdDataField::edgeXderiv(velocity, plane, normal1, normal2, iX, iY, iZ, 1)
                     + fdDataField::edgeYderiv(velocity, plane, normal1, normal2, iX, iY, iZ, 0))
                    / (T)2;
                el[tensor::xz] =
                    (fdDataField::edgeXderiv(velocity, plane, normal1, normal2, iX, iY, iZ, 2)
                     + fdDataField::edgeZderiv(velocity, plane, normal1, normal2, iX, iY, iZ, 0))
                    / (T)2;
                el[tensor::yy] =
                    fdDataField::edgeYderiv(velocity, plane, normal1, normal2, iX, iY, iZ, 1);
                el[tensor::yz] =
                    (fdDataField::edgeYderiv(velocity, plane, normal1, normal2, iX, iY, iZ, 2)
                     + fdDataField::edgeZderiv(velocity, plane, normal1, normal2, iX, iY, iZ, 1))
                    / (T)2;
                el[tensor::zz] =
                    fdDataField::edgeZderiv(velocity, plane, normal1, normal2, iX, iY, iZ, 2);
            }
        }
    }
}

template <typename T, int nDim>
void BoxStrainRateFunctional3D<T, nDim>::processCorner(
    int normalX, int normalY, int normalZ, Box3D domain, TensorField3D<T, nDim> &velocity,
    TensorField3D<T, SymmetricTensorImpl<T, nDim>::n> &S)
{
    typedef SymmetricTensorImpl<T, nDim> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iX2 = iX + offset.x;
                plint iY2 = iY + offset.y;
                plint iZ2 = iZ + offset.z;
                Array<T, SymmetricTensorImpl<T, nDim>::n> &el = S.get(iX2, iY2, iZ2);
                el[tensor::xx] =
                    fdDataField::cornerXderiv(velocity, normalX, normalY, normalZ, iX, iY, iZ, 0);
                el[tensor::xy] =
                    (fdDataField::cornerXderiv(velocity, normalX, normalY, normalZ, iX, iY, iZ, 1)
                     + fdDataField::cornerYderiv(
                         velocity, normalX, normalY, normalZ, iX, iY, iZ, 0))
                    / (T)2;
                el[tensor::xz] =
                    (fdDataField::cornerXderiv(velocity, normalX, normalY, normalZ, iX, iY, iZ, 2)
                     + fdDataField::cornerZderiv(
                         velocity, normalX, normalY, normalZ, iX, iY, iZ, 0))
                    / (T)2;
                el[tensor::yy] =
                    fdDataField::cornerYderiv(velocity, normalX, normalY, normalZ, iX, iY, iZ, 1);
                el[tensor::yz] =
                    (fdDataField::cornerYderiv(velocity, normalX, normalY, normalZ, iX, iY, iZ, 2)
                     + fdDataField::cornerZderiv(
                         velocity, normalX, normalY, normalZ, iX, iY, iZ, 1))
                    / (T)2;
                el[tensor::zz] =
                    fdDataField::cornerZderiv(velocity, normalX, normalY, normalZ, iX, iY, iZ, 2);
            }
        }
    }
}

template <typename T, int nDim>
BoxStrainRateFunctional3D<T, nDim> *BoxStrainRateFunctional3D<T, nDim>::clone() const
{
    return new BoxStrainRateFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void BoxStrainRateFunctional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT BoxStrainRateFunctional3D<T, nDim>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

/* ******** BoxBulkDivergenceFunctional3D ****************************************** */

template <typename T, int nDim>
void BoxBulkDivergenceFunctional3D<T, nDim>::process(
    Box3D domain, ScalarField3D<T> &divergence, TensorField3D<T, nDim> &velocity)
{
    Dot3D offset = computeRelativeDisplacement(divergence, velocity);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iX2 = iX + offset.x;
                plint iY2 = iY + offset.y;
                plint iZ2 = iZ + offset.z;
                divergence.get(iX, iY, iZ) = fdDataField::bulkXderiv(velocity, iX2, iY2, iZ2, 0)
                                             + fdDataField::bulkYderiv(velocity, iX2, iY2, iZ2, 1)
                                             + fdDataField::bulkZderiv(velocity, iX2, iY2, iZ2, 2);
            }
        }
    }
}

template <typename T, int nDim>
BoxBulkDivergenceFunctional3D<T, nDim> *BoxBulkDivergenceFunctional3D<T, nDim>::clone() const
{
    return new BoxBulkDivergenceFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void BoxBulkDivergenceFunctional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Divergence.
    modified[1] = modif::nothing;          // Velocity.
}

template <typename T, int nDim>
BlockDomain::DomainT BoxBulkDivergenceFunctional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Tensor_A_plus_B_functional3D ****************************************** */

template <typename T, int nDim>
void Tensor_A_plus_B_functional3D<T, nDim>::process(
    Box3D domain, std::vector<TensorField3D<T, nDim> *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);
    TensorField3D<T, nDim> &A = *fields[0];
    TensorField3D<T, nDim> &B = *fields[1];
    TensorField3D<T, nDim> &result = *fields[2];
    Dot3D offsetB = computeRelativeDisplacement(A, B);
    Dot3D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offsetResult.x, iY + offsetResult.y, iZ + offsetResult.z) =
                    A.get(iX, iY, iZ) + B.get(iX + offsetB.x, iY + offsetB.y, iZ + offsetB.z);
            }
        }
    }
}

template <typename T, int nDim>
Tensor_A_plus_B_functional3D<T, nDim> *Tensor_A_plus_B_functional3D<T, nDim>::clone() const
{
    return new Tensor_A_plus_B_functional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void Tensor_A_plus_B_functional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT Tensor_A_plus_B_functional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Tensor_A_minus_B_functional3D ************************************ */

template <typename T, int nDim>
void Tensor_A_minus_B_functional3D<T, nDim>::process(
    Box3D domain, std::vector<TensorField3D<T, nDim> *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);
    TensorField3D<T, nDim> &A = *fields[0];
    TensorField3D<T, nDim> &B = *fields[1];
    TensorField3D<T, nDim> &result = *fields[2];
    Dot3D offsetB = computeRelativeDisplacement(A, B);
    Dot3D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offsetResult.x, iY + offsetResult.y, iZ + offsetResult.z) =
                    A.get(iX, iY, iZ) - B.get(iX + offsetB.x, iY + offsetB.y, iZ + offsetB.z);
            }
        }
    }
}

template <typename T, int nDim>
Tensor_A_minus_B_functional3D<T, nDim> *Tensor_A_minus_B_functional3D<T, nDim>::clone() const
{
    return new Tensor_A_minus_B_functional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void Tensor_A_minus_B_functional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT Tensor_A_minus_B_functional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Tensor_A_times_B_functional3D ************************************ */

template <typename T, int nDim>
void Tensor_A_times_B_functional3D<T, nDim>::process(
    Box3D domain, std::vector<TensorField3D<T, nDim> *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);
    TensorField3D<T, nDim> &A = *fields[0];
    TensorField3D<T, nDim> &B = *fields[1];
    TensorField3D<T, nDim> &result = *fields[2];
    Dot3D offsetB = computeRelativeDisplacement(A, B);
    Dot3D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offsetResult.x, iY + offsetResult.y, iZ + offsetResult.z) =
                    A.get(iX, iY, iZ) * B.get(iX + offsetB.x, iY + offsetB.y, iZ + offsetB.z);
            }
        }
    }
}

template <typename T, int nDim>
Tensor_A_times_B_functional3D<T, nDim> *Tensor_A_times_B_functional3D<T, nDim>::clone() const
{
    return new Tensor_A_times_B_functional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void Tensor_A_times_B_functional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT Tensor_A_times_B_functional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** IndexContraction_SymmetricTensor_A_SymmetricTensor_B_functional3D
 * ************************************ */

template <typename T, int nDim>
void IndexContraction_SymmetricTensor_A_SymmetricTensor_B_functional3D<
    T, nDim>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);
    TensorField3D<T, SymmetricTensorImpl<T, nDim>::n> &A =
        *dynamic_cast<TensorField3D<T, SymmetricTensorImpl<T, nDim>::n> *>(fields[0]);
    TensorField3D<T, SymmetricTensorImpl<T, nDim>::n> &B =
        *dynamic_cast<TensorField3D<T, SymmetricTensorImpl<T, nDim>::n> *>(fields[1]);
    ScalarField3D<T> &result = *dynamic_cast<ScalarField3D<T> *>(fields[2]);

    Dot3D offsetB = computeRelativeDisplacement(A, B);
    Dot3D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offsetResult.x, iY + offsetResult.y, iZ + offsetResult.z) =
                    SymmetricTensorImpl<T, nDim>::contractIndexes(
                        A.get(iX, iY, iZ), B.get(iX + offsetB.x, iY + offsetB.y, iZ + offsetB.z));
            }
        }
    }
}

template <typename T, int nDim>
IndexContraction_SymmetricTensor_A_SymmetricTensor_B_functional3D<T, nDim>
    *IndexContraction_SymmetricTensor_A_SymmetricTensor_B_functional3D<T, nDim>::clone() const
{
    return new IndexContraction_SymmetricTensor_A_SymmetricTensor_B_functional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void IndexContraction_SymmetricTensor_A_SymmetricTensor_B_functional3D<
    T, nDim>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT
    IndexContraction_SymmetricTensor_A_SymmetricTensor_B_functional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** TensorProduct_A_A_functional3D ************************************ */

template <typename T, int nDim>
void TensorProduct_A_A_functional3D<T, nDim>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_PRECONDITION(fields.size() == 2);
    TensorField3D<T, nDim> &A = *dynamic_cast<TensorField3D<T, nDim> *>(fields[0]);
    TensorField3D<T, SymmetricTensorImpl<T, nDim>::n> &result =
        *dynamic_cast<TensorField3D<T, SymmetricTensorImpl<T, nDim>::n> *>(fields[1]);

    Dot3D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iPi = 0;
                for (plint iA = 0; iA < nDim; ++iA) {
                    for (plint iB = iA; iB < nDim; ++iB) {
                        result.get(
                            iX + offsetResult.x, iY + offsetResult.y, iZ + offsetResult.z)[iPi] =
                            A.get(iX, iY, iZ)[iA] * A.get(iX, iY, iZ)[iB];
                        ++iPi;
                    }
                }
            }
        }
    }
}

template <typename T, int nDim>
TensorProduct_A_A_functional3D<T, nDim> *TensorProduct_A_A_functional3D<T, nDim>::clone() const
{
    return new TensorProduct_A_A_functional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void TensorProduct_A_A_functional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT TensorProduct_A_A_functional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** InterpolateTensorFieldFunctional3D ************************************ */

template <typename T, int nDim>
InterpolateTensorFieldFunctional3D<T, nDim>::InterpolateTensorFieldFunctional3D(plint N_, T t_) :
    N(N_), t(t_)
{
    PLB_ASSERT(N >= 2);
    w.resize(N);
    plb::fdWeights<T, 0, 1>().getWeights(N, t, &w[0]);
}

template <typename T, int nDim>
void InterpolateTensorFieldFunctional3D<T, nDim>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT((plint)blocks.size() == N + 1);
    std::vector<TensorField3D<T, nDim> *> fields(N);
    for (plint iField = 0; iField < N; iField++) {
        fields[iField] = dynamic_cast<TensorField3D<T, nDim> *>(blocks[iField]);
        PLB_ASSERT(fields[iField]);
    }
    TensorField3D<T, nDim> *result = dynamic_cast<TensorField3D<T, nDim> *>(blocks[N]);
    PLB_ASSERT(result);

    std::vector<Dot3D> ofsF(N);
    for (plint iField = 0; iField < N; iField++) {
        ofsF[iField] = computeRelativeDisplacement(*fields[0], *fields[iField]);
    }
    Dot3D ofsR = computeRelativeDisplacement(*fields[0], *result);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, nDim> &res = result->get(iX + ofsR.x, iY + ofsR.y, iZ + ofsR.z);
                res.resetToZero();
                for (plint iField = 0; iField < N; iField++) {
                    res += w[iField]
                           * fields[iField]->get(
                               iX + ofsF[iField].x, iY + ofsF[iField].y, iZ + ofsF[iField].z);
                }
            }
        }
    }
}

template <typename T, int nDim>
InterpolateTensorFieldFunctional3D<T, nDim> *InterpolateTensorFieldFunctional3D<T, nDim>::clone()
    const
{
    return new InterpolateTensorFieldFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void InterpolateTensorFieldFunctional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    std::fill(
        modified.begin(), modified.end(),
        modif::nothing);  // All tensor fields that participate in the interpolation.
    modified.back() = modif::staticVariables;  // The result of the interpolation.
}

template <typename T, int nDim>
BlockDomain::DomainT InterpolateTensorFieldFunctional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Scalar_A_times_Tensor_B_functional3D ************************************ */

template <typename T, int nDim>
void Scalar_A_times_Tensor_B_functional3D<T, nDim>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);
    ScalarField3D<T> &A = *dynamic_cast<ScalarField3D<T> *>(fields[0]);
    TensorField3D<T, nDim> &B = *dynamic_cast<TensorField3D<T, nDim> *>(fields[1]);
    TensorField3D<T, nDim> &result = *dynamic_cast<TensorField3D<T, nDim> *>(fields[2]);

    Dot3D offsetB = computeRelativeDisplacement(A, B);
    Dot3D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offsetResult.x, iY + offsetResult.y, iZ + offsetResult.z) =
                    A.get(iX, iY, iZ) * B.get(iX + offsetB.x, iY + offsetB.y, iZ + offsetB.z);
            }
        }
    }
}

template <typename T, int nDim>
Scalar_A_times_Tensor_B_functional3D<T, nDim>
    *Scalar_A_times_Tensor_B_functional3D<T, nDim>::clone() const
{
    return new Scalar_A_times_Tensor_B_functional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void Scalar_A_times_Tensor_B_functional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT Scalar_A_times_Tensor_B_functional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Tensor_A_dividedBy_B_functional3D ************************************ */

template <typename T, int nDim>
void Tensor_A_dividedBy_B_functional3D<T, nDim>::process(
    Box3D domain, std::vector<TensorField3D<T, nDim> *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);
    TensorField3D<T, nDim> &A = *fields[0];
    TensorField3D<T, nDim> &B = *fields[1];
    TensorField3D<T, nDim> &result = *fields[2];
    Dot3D offsetB = computeRelativeDisplacement(A, B);
    Dot3D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offsetResult.x, iY + offsetResult.y, iZ + offsetResult.z) =
                    A.get(iX, iY, iZ) / B.get(iX + offsetB.x, iY + offsetB.y, iZ + offsetB.z);
            }
        }
    }
}

template <typename T, int nDim>
Tensor_A_dividedBy_B_functional3D<T, nDim> *Tensor_A_dividedBy_B_functional3D<T, nDim>::clone()
    const
{
    return new Tensor_A_dividedBy_B_functional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void Tensor_A_dividedBy_B_functional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT Tensor_A_dividedBy_B_functional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Tensor_A_plus_B_inplace_functional3D ************************************ */

template <typename T, int nDim>
void Tensor_A_plus_B_inplace_functional3D<T, nDim>::process(
    Box3D domain, TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B)
{
    Dot3D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                A.get(iX, iY, iZ) += B.get(iX + offset.x, iY + offset.y, iZ + offset.z);
            }
        }
    }
}

template <typename T, int nDim>
Tensor_A_plus_B_inplace_functional3D<T, nDim>
    *Tensor_A_plus_B_inplace_functional3D<T, nDim>::clone() const
{
    return new Tensor_A_plus_B_inplace_functional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void Tensor_A_plus_B_inplace_functional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, int nDim>
BlockDomain::DomainT Tensor_A_plus_B_inplace_functional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Tensor_A_minus_B_inplace_functional3D ************************************ */

template <typename T, int nDim>
void Tensor_A_minus_B_inplace_functional3D<T, nDim>::process(
    Box3D domain, TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B)
{
    Dot3D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                A.get(iX, iY, iZ) -= B.get(iX + offset.x, iY + offset.y, iZ + offset.z);
            }
        }
    }
}

template <typename T, int nDim>
Tensor_A_minus_B_inplace_functional3D<T, nDim>
    *Tensor_A_minus_B_inplace_functional3D<T, nDim>::clone() const
{
    return new Tensor_A_minus_B_inplace_functional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void Tensor_A_minus_B_inplace_functional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, int nDim>
BlockDomain::DomainT Tensor_A_minus_B_inplace_functional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Tensor_A_times_B_inplace_functional3D ************************************ */

template <typename T, int nDim>
void Tensor_A_times_B_inplace_functional3D<T, nDim>::process(
    Box3D domain, TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B)
{
    Dot3D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                A.get(iX, iY, iZ) *= B.get(iX + offset.x, iY + offset.y, iZ + offset.z);
            }
        }
    }
}

template <typename T, int nDim>
Tensor_A_times_B_inplace_functional3D<T, nDim>
    *Tensor_A_times_B_inplace_functional3D<T, nDim>::clone() const
{
    return new Tensor_A_times_B_inplace_functional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void Tensor_A_times_B_inplace_functional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, int nDim>
BlockDomain::DomainT Tensor_A_times_B_inplace_functional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Tensor_A_times_alpha_inplace_functional3D ************************************ */

template <typename T, int nDim>
Tensor_A_times_alpha_inplace_functional3D<T, nDim>::Tensor_A_times_alpha_inplace_functional3D(
    T alpha_) :
    alpha(alpha_)
{ }

template <typename T, int nDim>
void Tensor_A_times_alpha_inplace_functional3D<T, nDim>::process(
    Box3D domain, TensorField3D<T, nDim> &A)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                A.get(iX, iY, iZ) *= alpha;
            }
        }
    }
}

template <typename T, int nDim>
Tensor_A_times_alpha_inplace_functional3D<T, nDim>
    *Tensor_A_times_alpha_inplace_functional3D<T, nDim>::clone() const
{
    return new Tensor_A_times_alpha_inplace_functional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void Tensor_A_times_alpha_inplace_functional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT Tensor_A_times_alpha_inplace_functional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

/* ******** Tensor_A_times_alpha_functional3D ************************************ */

template <typename T, int nDim>
Tensor_A_times_alpha_functional3D<T, nDim>::Tensor_A_times_alpha_functional3D(T alpha_) :
    alpha(alpha_)
{ }

template <typename T, int nDim>
void Tensor_A_times_alpha_functional3D<T, nDim>::process(
    Box3D domain, TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &result)
{
    Dot3D offset = computeRelativeDisplacement(A, result);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint bX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint bY = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint bZ = iZ + offset.z;
                result.get(bX, bY, bZ) = alpha * A.get(iX, iY, iZ);
            }
        }
    }
}

template <typename T, int nDim>
Tensor_A_times_alpha_functional3D<T, nDim> *Tensor_A_times_alpha_functional3D<T, nDim>::clone()
    const
{
    return new Tensor_A_times_alpha_functional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void Tensor_A_times_alpha_functional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT Tensor_A_times_alpha_functional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Tensor_A_dividedBy_B_inplace_functional3D ************************************ */

template <typename T, int nDim>
void Tensor_A_dividedBy_B_inplace_functional3D<T, nDim>::process(
    Box3D domain, TensorField3D<T, nDim> &A, TensorField3D<T, nDim> &B)
{
    Dot3D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                A.get(iX, iY, iZ) /= B.get(iX + offset.x, iY + offset.y, iZ + offset.z);
            }
        }
    }
}

template <typename T, int nDim>
Tensor_A_dividedBy_B_inplace_functional3D<T, nDim>
    *Tensor_A_dividedBy_B_inplace_functional3D<T, nDim>::clone() const
{
    return new Tensor_A_dividedBy_B_inplace_functional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void Tensor_A_dividedBy_B_inplace_functional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, int nDim>
BlockDomain::DomainT Tensor_A_dividedBy_B_inplace_functional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Tensor_A_dividedBy_Scalar_B_inplace_functional3D ************************************ */

template <typename T, int nDim>
void Tensor_A_dividedBy_Scalar_B_inplace_functional3D<T, nDim>::process(
    Box3D domain, ScalarField3D<T> &B, TensorField3D<T, nDim> &A)
{
    Dot3D offset = computeRelativeDisplacement(B, A);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                A.get(iX + offset.x, iY + offset.y, iZ + offset.z) /= B.get(iX, iY, iZ);
            }
        }
    }
}

template <typename T, int nDim>
Tensor_A_dividedBy_Scalar_B_inplace_functional3D<T, nDim>
    *Tensor_A_dividedBy_Scalar_B_inplace_functional3D<T, nDim>::clone() const
{
    return new Tensor_A_dividedBy_Scalar_B_inplace_functional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void Tensor_A_dividedBy_Scalar_B_inplace_functional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, int nDim>
BlockDomain::DomainT Tensor_A_dividedBy_Scalar_B_inplace_functional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Masked_Tensor_A_dividedBy_Scalar_B_inplace_functional3D
 * ************************************ */

template <typename T, int nDim>
Masked_Tensor_A_dividedBy_Scalar_B_inplace_functional3D<
    T, nDim>::Masked_Tensor_A_dividedBy_Scalar_B_inplace_functional3D(int flag_) :
    flag(flag_)
{ }

template <typename T, int nDim>
void Masked_Tensor_A_dividedBy_Scalar_B_inplace_functional3D<T, nDim>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_ASSERT(fields.size() == 3);
    TensorField3D<T, nDim> *tensor = dynamic_cast<TensorField3D<T, nDim> *>(fields[0]);
    PLB_ASSERT(tensor);
    ScalarField3D<T> *scalar = dynamic_cast<ScalarField3D<T> *>(fields[1]);
    PLB_ASSERT(scalar);
    ScalarField3D<int> *mask = dynamic_cast<ScalarField3D<int> *>(fields[2]);
    PLB_ASSERT(mask);

    Dot3D ofsS = computeRelativeDisplacement(*tensor, *scalar);
    Dot3D ofsM = computeRelativeDisplacement(*tensor, *mask);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (mask->get(iX + ofsM.x, iY + ofsM.y, iZ + ofsM.z) == flag) {
                    tensor->get(iX, iY, iZ) /= scalar->get(iX + ofsS.x, iY + ofsS.y, iZ + ofsS.z);
                }
            }
        }
    }
}

template <typename T, int nDim>
Masked_Tensor_A_dividedBy_Scalar_B_inplace_functional3D<T, nDim>
    *Masked_Tensor_A_dividedBy_Scalar_B_inplace_functional3D<T, nDim>::clone() const
{
    return new Masked_Tensor_A_dividedBy_Scalar_B_inplace_functional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void Masked_Tensor_A_dividedBy_Scalar_B_inplace_functional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template <typename T, int nDim>
BlockDomain::DomainT Masked_Tensor_A_dividedBy_Scalar_B_inplace_functional3D<T, nDim>::appliesTo()
    const
{
    return BlockDomain::bulk;
}

/* ******** Normalize_Tensor_functional3D ****************************************** */

template <typename T, int nDim>
void Normalize_Tensor_functional3D<T, nDim>::process(
    Box3D domain, TensorField3D<T, nDim> &data, TensorField3D<T, nDim> &result)
{
    Dot3D offset = computeRelativeDisplacement(data, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, nDim> &a = data.get(iX, iY, iZ);
                T normA = norm(a);
                if (util::isZero(normA)) {
                    result.get(iX + offset.x, iY + offset.y, iZ + offset.z).resetToZero();
                } else {
                    result.get(iX + offset.x, iY + offset.y, iZ + offset.z) = a / normA;
                }
            }
        }
    }
}

template <typename T, int nDim>
Normalize_Tensor_functional3D<T, nDim> *Normalize_Tensor_functional3D<T, nDim>::clone() const
{
    return new Normalize_Tensor_functional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void Normalize_Tensor_functional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

/* ******** Normalize_Tensor_inplace_functional3D ************************************ */

template <typename T, int nDim>
void Normalize_Tensor_inplace_functional3D<T, nDim>::process(
    Box3D domain, TensorField3D<T, nDim> &data)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, nDim> &a = data.get(iX, iY, iZ);
                T normA = norm(a);
                if (util::isZero(normA)) {
                    a.resetToZero();
                } else {
                    a /= normA;
                }
            }
        }
    }
}

template <typename T, int nDim>
Normalize_Tensor_inplace_functional3D<T, nDim>
    *Normalize_Tensor_inplace_functional3D<T, nDim>::clone() const
{
    return new Normalize_Tensor_inplace_functional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void Normalize_Tensor_inplace_functional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************* Class LBMsmoothenTensor3D ******************* */

template <typename T, int nDim, template <typename U> class Descriptor>
void LBMsmoothenTensor3D<T, nDim, Descriptor>::process(
    Box3D domain, TensorField3D<T, nDim> &data, TensorField3D<T, nDim> &result)
{
    typedef Descriptor<T> D;
    Dot3D offset = computeRelativeDisplacement(data, result);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offset.x, iY + offset.y, iZ + offset.z).resetToZero();
                T sum = (T)0;
                for (plint iPop = 1; iPop < D::q; ++iPop) {
                    plint nextX = iX + D::c[iPop][0];
                    plint nextY = iY + D::c[iPop][1];
                    plint nextZ = iZ + D::c[iPop][2];
                    sum += D::t[iPop];
                    result.get(iX + offset.x, iY + offset.y, iZ + offset.z) +=
                        D::t[iPop] * data.get(nextX, nextY, nextZ);
                }
                result.get(iX + offset.x, iY + offset.y, iZ + offset.z) /= sum;
            }
        }
    }
}

template <typename T, int nDim, template <typename U> class Descriptor>
LBMsmoothenTensor3D<T, nDim, Descriptor> *LBMsmoothenTensor3D<T, nDim, Descriptor>::clone() const
{
    return new LBMsmoothenTensor3D<T, nDim, Descriptor>(*this);
}

template <typename T, int nDim, template <typename U> class Descriptor>
void LBMsmoothenTensor3D<T, nDim, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

/* ************* Class LBMsmoothenTensorInPlace3D ******************* */

template <typename T, int nDim, template <typename U> class Descriptor>
void LBMsmoothenTensorInPlace3D<T, nDim, Descriptor>::process(
    Box3D domain, TensorField3D<T, nDim> &data)
{
    typedef Descriptor<T> D;

    TensorField3D<T, nDim> smoothData(domain.getNx(), domain.getNy(), domain.getNz());
    Dot3D offset(-domain.x0, -domain.y0, -domain.z0);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                smoothData.get(iX + offset.x, iY + offset.y, iZ + offset.z).resetToZero();
                T sum = (T)0;
                for (plint iPop = 1; iPop < D::q; ++iPop) {
                    plint nextX = iX + D::c[iPop][0];
                    plint nextY = iY + D::c[iPop][1];
                    plint nextZ = iZ + D::c[iPop][2];
                    sum += D::t[iPop];
                    smoothData.get(iX + offset.x, iY + offset.y, iZ + offset.z) +=
                        D::t[iPop] * data.get(nextX, nextY, nextZ);
                }
                smoothData.get(iX + offset.x, iY + offset.y, iZ + offset.z) /= sum;
            }
        }
    }

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                data.get(iX, iY, iZ) = smoothData.get(iX + offset.x, iY + offset.y, iZ + offset.z);
            }
        }
    }
}

template <typename T, int nDim, template <typename U> class Descriptor>
LBMsmoothenTensorInPlace3D<T, nDim, Descriptor>
    *LBMsmoothenTensorInPlace3D<T, nDim, Descriptor>::clone() const
{
    return new LBMsmoothenTensorInPlace3D<T, nDim, Descriptor>(*this);
}

template <typename T, int nDim, template <typename U> class Descriptor>
void LBMsmoothenTensorInPlace3D<T, nDim, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************* Class SmoothenTensor3D ******************* */

template <typename T, int nDim>
void SmoothenTensor3D<T, nDim>::process(
    Box3D domain, TensorField3D<T, nDim> &data, TensorField3D<T, nDim> &result)
{
    Dot3D offset = computeRelativeDisplacement(data, result);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, 3> *res = &result.get(iX + offset.x, iY + offset.y, iZ + offset.z);
                res->resetToZero();
                int n = 0;
                for (int i = -1; i < 2; i++) {
                    plint nextX = iX + i;
                    for (int j = -1; j < 2; j++) {
                        plint nextY = iY + j;
                        for (int k = -1; k < 2; k++) {
                            plint nextZ = iZ + k;
                            if (!(i == 0 && j == 0 && k == 0)) {
                                n++;
                                *res += data.get(nextX, nextY, nextZ);
                            }
                        }
                    }
                }
                *res /= (T)n;
            }
        }
    }
}

template <typename T, int nDim>
SmoothenTensor3D<T, nDim> *SmoothenTensor3D<T, nDim>::clone() const
{
    return new SmoothenTensor3D<T, nDim>(*this);
}

template <typename T, int nDim>
void SmoothenTensor3D<T, nDim>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

/* ************* Class MollifyTensor3D ******************* */

template <typename T, int nDim>
MollifyTensor3D<T, nDim>::MollifyTensor3D(T l_, plint d_, Box3D globalDomain_, int exclusionFlag_) :
    l(l_), d(d_), globalDomain(globalDomain_), exclusionFlag(exclusionFlag_)
{ }

template <typename T, int nDim>
void MollifyTensor3D<T, nDim>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_ASSERT(fields.size() == 3);
    TensorField3D<T, nDim> *tensor = dynamic_cast<TensorField3D<T, nDim> *>(fields[0]);
    PLB_ASSERT(tensor);
    ScalarField3D<int> *flag = dynamic_cast<ScalarField3D<int> *>(fields[1]);
    PLB_ASSERT(flag);
    TensorField3D<T, nDim> *result = dynamic_cast<TensorField3D<T, nDim> *>(fields[2]);
    PLB_ASSERT(result);

    Dot3D absOfs = tensor->getLocation();

    Dot3D ofsF = computeRelativeDisplacement(*tensor, *flag);
    Dot3D ofsR = computeRelativeDisplacement(*tensor, *result);

    static T pi = std::acos((T)-1);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (flag->get(iX + ofsF.x, iY + ofsF.y, iZ + ofsF.z) == exclusionFlag) {
                    result->get(iX + ofsR.x, iY + ofsR.y, iZ + ofsR.z) = tensor->get(iX, iY, iZ);
                    continue;
                }

                Array<T, nDim> sum;
                sum.resetToZero();
                for (plint dx = -d; dx <= d; dx++) {
                    plint i = iX + dx;
                    for (plint dy = -d; dy <= d; dy++) {
                        plint j = iY + dy;
                        for (plint dz = -d; dz <= d; dz++) {
                            plint k = iZ + dz;
                            if (contained(i + absOfs.x, j + absOfs.y, k + absOfs.z, globalDomain)) {
                                if (flag->get(i + ofsF.x, j + ofsF.y, k + ofsF.z) != exclusionFlag)
                                {
                                    T r = std::sqrt((T)dx * dx + (T)dy * dy + (T)dz * dz);
                                    Array<T, nDim> integrand = tensor->get(i, j, k);
                                    T mollifier = 0.0;
                                    if (r < l) {
                                        mollifier = (1.0 + std::cos(pi * r / l));
                                    }
                                    sum += integrand * mollifier;
                                }
                            }
                        }
                    }
                }
                sum *= (1.0 / (2.0 * l));
                result->get(iX + ofsR.x, iY + ofsR.y, iZ + ofsR.z) = sum;
            }
        }
    }
}

template <typename T, int nDim>
MollifyTensor3D<T, nDim> *MollifyTensor3D<T, nDim>::clone() const
{
    return new MollifyTensor3D<T, nDim>(*this);
}

template <typename T, int nDim>
void MollifyTensor3D<T, nDim>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Data.
    modified[1] = modif::nothing;          // Flags.
    modified[2] = modif::staticVariables;  // Result.
}

/* ************* Class LBMcomputeDivergence3D ******************* */

template <typename T, template <typename U> class Descriptor>
void LBMcomputeDivergence3D<T, Descriptor>::process(
    Box3D domain, ScalarField3D<T> &divergence, TensorField3D<T, 3> &vectorField)
{
    typedef Descriptor<T> D;
    Dot3D ofs = computeRelativeDisplacement(divergence, vectorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                // Compute the divergence of a vector function "the lattice Boltzmann way".
                T sum = T();
                for (plint iPop = 1; iPop < D::q; ++iPop) {
                    plint nextX = iX + D::c[iPop][0];
                    plint nextY = iY + D::c[iPop][1];
                    plint nextZ = iZ + D::c[iPop][2];

                    Array<T, 3> &vectorFunction =
                        vectorField.get(nextX + ofs.x, nextY + ofs.y, nextZ + ofs.z);

                    T tmp = D::c[iPop][0] * vectorFunction[0] + D::c[iPop][1] * vectorFunction[1]
                            + D::c[iPop][2] * vectorFunction[2];

                    sum += D::t[iPop] * tmp;
                }
                divergence.get(iX, iY, iZ) = D::invCs2 * sum;
            }
        }
    }
}

/* ************* Class UpdateAveTensorTransientStatistics3D ******************* */

template <typename T, int nDim>
UpdateAveTensorTransientStatistics3D<T, nDim>::UpdateAveTensorTransientStatistics3D(plint n_) :
    n(n_)
{ }

template <typename T, int nDim>
void UpdateAveTensorTransientStatistics3D<T, nDim>::process(
    Box3D domain, TensorField3D<T, nDim> &tensor, TensorField3D<T, nDim> &avg)
{
    Dot3D offset = computeRelativeDisplacement(tensor, avg);

    T nMinusOne = (T)n - (T)1;
    T oneOverN = (T)1 / (T)n;

    for (plint iX = domain.x0; iX <= domain.x1; iX++) {
        for (plint iY = domain.y0; iY <= domain.y1; iY++) {
            for (plint iZ = domain.z0; iZ <= domain.z1; iZ++) {
                for (plint iD = 0; iD < nDim; ++iD) {
                    avg.get(iX + offset.x, iY + offset.y, iZ + offset.z)[iD] =
                        oneOverN
                        * (nMinusOne * avg.get(iX + offset.x, iY + offset.y, iZ + offset.z)[iD]
                           + tensor.get(iX, iY, iZ)[iD]);
                }
            }
        }
    }
    ++n;
}

template <typename T, int nDim>
UpdateAveTensorTransientStatistics3D<T, nDim>
    *UpdateAveTensorTransientStatistics3D<T, nDim>::clone() const
{
    return new UpdateAveTensorTransientStatistics3D<T, nDim>(*this);
}

template <typename T, int nDim>
void UpdateAveTensorTransientStatistics3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Tensor field.
    modified[1] = modif::staticVariables;  // Statistics field.
}

template <typename T, int nDim>
BlockDomain::DomainT UpdateAveTensorTransientStatistics3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ************* Class UpdateDevTensorTransientStatistics3D ******************* */

template <typename T, int nDim>
UpdateDevTensorTransientStatistics3D<T, nDim>::UpdateDevTensorTransientStatistics3D(plint n_) :
    n(n_)
{ }

template <typename T, int nDim>
void UpdateDevTensorTransientStatistics3D<T, nDim>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 3);
    TensorField3D<T, nDim> *tensor = dynamic_cast<TensorField3D<T, nDim> *>(blocks[0]);
    TensorField3D<T, nDim> *average = dynamic_cast<TensorField3D<T, nDim> *>(blocks[1]);
    TensorField3D<T, (nDim * (nDim + 1)) / 2> *statistics =
        dynamic_cast<TensorField3D<T, (nDim * (nDim + 1)) / 2> *>(blocks[2]);
    PLB_ASSERT(tensor);
    PLB_ASSERT(average);
    PLB_ASSERT(statistics);

    Dot3D ofsA = computeRelativeDisplacement(*tensor, *average);
    Dot3D ofsS = computeRelativeDisplacement(*tensor, *statistics);

    T nMinusOne = (T)n - (T)1;
    T oneOverN = (T)1 / (T)n;

    for (plint iX = domain.x0; iX <= domain.x1; iX++) {
        for (plint iY = domain.y0; iY <= domain.y1; iY++) {
            for (plint iZ = domain.z0; iZ <= domain.z1; iZ++) {
                plint iPi = 0;
                for (plint iA = 0; iA < nDim; ++iA) {
                    T newDataA = tensor->get(iX, iY, iZ)[iA]
                                 - average->get(iX + ofsA.x, iY + ofsA.y, iZ + ofsA.z)[iA];

                    for (plint iB = iA; iB < nDim; ++iB) {
                        T newDataB = tensor->get(iX, iY, iZ)[iB]
                                     - average->get(iX + ofsA.x, iY + ofsA.y, iZ + ofsA.z)[iB];

                        T oldDev = statistics->get(iX + ofsS.x, iY + ofsS.y, iZ + ofsS.z)[iPi];
                        statistics->get(iX + ofsS.x, iY + ofsS.y, iZ + ofsS.z)[iPi] =
                            oneOverN * (nMinusOne * oldDev + newDataA * newDataB);

                        ++iPi;
                    }
                }
            }
        }
    }
    ++n;
}

template <typename T, int nDim>
UpdateDevTensorTransientStatistics3D<T, nDim>
    *UpdateDevTensorTransientStatistics3D<T, nDim>::clone() const
{
    return new UpdateDevTensorTransientStatistics3D<T, nDim>(*this);
}

template <typename T, int nDim>
void UpdateDevTensorTransientStatistics3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Tensor field.
    modified[1] = modif::nothing;          // Mean value field.
    modified[2] = modif::staticVariables;  // Statistics field.
}

template <typename T, int nDim>
BlockDomain::DomainT UpdateDevTensorTransientStatistics3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* *************** PART IV ******************************************** */
/* *************** Analysis of the NTensor-field ********************** */
/* ******************************************************************** */

template <typename T>
void ExtractNTensorSubDomainFunctional3D<T>::process(
    Box3D domain, NTensorField3D<T> &field1, NTensorField3D<T> &field2)
{
    PLB_ASSERT(field1.getNdim() == field2.getNdim());

    Dot3D offset = computeRelativeDisplacement(field1, field2);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                for (plint iDim = 0; iDim < field1.getNdim(); ++iDim) {
                    field2.get(iX + offset.x, iY + offset.y, iZ + offset.z)[iDim] =
                        field1.get(iX, iY, iZ)[iDim];
                }
            }
        }
    }
}

template <typename T>
ExtractNTensorSubDomainFunctional3D<T> *ExtractNTensorSubDomainFunctional3D<T>::clone() const
{
    return new ExtractNTensorSubDomainFunctional3D<T>(*this);
}

template <typename T>
void ExtractNTensorSubDomainFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT ExtractNTensorSubDomainFunctional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T>
MaskedNTensorNeumannInLayersFunctional3D<T>::MaskedNTensorNeumannInLayersFunctional3D(
    int toFlag_, int fromFlag_, int noFlag_) :
    toFlag(toFlag_), fromFlag(fromFlag_), noFlag(noFlag_)
{
    PLB_ASSERT(toFlag != fromFlag);
    PLB_ASSERT(noFlag != toFlag && noFlag != fromFlag);
}

template <typename T>
void MaskedNTensorNeumannInLayersFunctional3D<T>::process(
    Box3D domain, NTensorField3D<T> &data, NTensorField3D<int> &mask)
{
    PLB_ASSERT(mask.getNdim() == 1);
    plint nDim = data.getNdim();
    Dot3D offset = computeRelativeDisplacement(data, mask);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (*mask.get(iX + offset.x, iY + offset.y, iZ + offset.z) == toFlag) {
                    std::vector<T> sum(nDim, (T)0);
                    plint n = 0;
                    for (plint dx = -1; dx <= +1; dx++) {
                        plint nextX = iX + dx;
                        for (plint dy = -1; dy <= +1; dy++) {
                            plint nextY = iY + dy;
                            for (plint dz = -1; dz <= +1; dz++) {
                                plint nextZ = iZ + dz;
                                if (*mask.get(nextX + offset.x, nextY + offset.y, nextZ + offset.z)
                                    == fromFlag) {
                                    T *fromData = data.get(nextX, nextY, nextZ);
                                    for (plint iDim = 0; iDim < nDim; iDim++) {
                                        sum[iDim] += fromData[iDim];
                                    }
                                    n++;
                                }
                            }
                        }
                    }
                    if (n != 0) {
                        T *toData = data.get(iX, iY, iZ);
                        for (plint iDim = 0; iDim < nDim; iDim++) {
                            toData[iDim] = sum[iDim] / (T)n;
                        }
                        *mask.get(iX + offset.x, iY + offset.y, iZ + offset.z) = noFlag;
                    }  // If no admissible neighborhood is found, then the data is not changed.
                }
            }
        }
    }
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (*mask.get(iX + offset.x, iY + offset.y, iZ + offset.z) == noFlag) {
                    *mask.get(iX + offset.x, iY + offset.y, iZ + offset.z) = fromFlag;
                }
            }
        }
    }
}

template <typename T>
MaskedNTensorNeumannInLayersFunctional3D<T> *MaskedNTensorNeumannInLayersFunctional3D<T>::clone()
    const
{
    return new MaskedNTensorNeumannInLayersFunctional3D<T>(*this);
}

template <typename T>
void MaskedNTensorNeumannInLayersFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::staticVariables;
}

}  // namespace plb

#endif  // DATA_ANALYSIS_FUNCTIONAL_3D_HH
