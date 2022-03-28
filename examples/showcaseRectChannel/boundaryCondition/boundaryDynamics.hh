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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- generic implementation.
 */
#ifndef BOUNDARY_DYNAMICS_HH
#define BOUNDARY_DYNAMICS_HH

#include "boundaryCondition/boundaryDynamics.h"
#include "core/cell.h"
#include "core/dynamicsIdentifiers.h"
#include "latticeBoltzmann/externalFieldAccess.h"
#include "latticeBoltzmann/indexTemplates.h"

namespace plb {

/* *************** Class BoundaryCompositeDynamics *********************** */

template <typename T, template <typename U> class Descriptor>
int BoundaryCompositeDynamics<T, Descriptor>::id =
    meta::registerCompositeDynamics<T, Descriptor, BoundaryCompositeDynamics<T, Descriptor> >(
        "Boundary_Composite");

template <typename T, template <typename U> class Descriptor>
BoundaryCompositeDynamics<T, Descriptor>::BoundaryCompositeDynamics(
    Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_) :
    PreparePopulationsDynamics<T, Descriptor>(baseDynamics_, automaticPrepareCollision_)
{ }

template <typename T, template <typename U> class Descriptor>
BoundaryCompositeDynamics<T, Descriptor> *BoundaryCompositeDynamics<T, Descriptor>::clone() const
{
    return new BoundaryCompositeDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int BoundaryCompositeDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
bool BoundaryCompositeDynamics<T, Descriptor>::isBoundary() const
{
    return true;
}

template <typename T, template <typename U> class Descriptor>
T BoundaryCompositeDynamics<T, Descriptor>::computeDensity(Cell<T, Descriptor> const &cell) const
{
    Cell<T, Descriptor> tmpCell(cell);
    this->completePopulations(tmpCell);
    return this->getBaseDynamics().computeDensity(tmpCell);
}

template <typename T, template <typename U> class Descriptor>
T BoundaryCompositeDynamics<T, Descriptor>::computePressure(Cell<T, Descriptor> const &cell) const
{
    Cell<T, Descriptor> tmpCell(cell);
    this->completePopulations(tmpCell);
    return this->getBaseDynamics().computePressure(tmpCell);
}

template <typename T, template <typename U> class Descriptor>
void BoundaryCompositeDynamics<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &velocity) const
{
    Cell<T, Descriptor> tmpCell(cell);
    this->completePopulations(tmpCell);
    this->getBaseDynamics().computeVelocity(tmpCell, velocity);
}

template <typename T, template <typename U> class Descriptor>
T BoundaryCompositeDynamics<T, Descriptor>::computeTemperature(
    Cell<T, Descriptor> const &cell) const
{
    Cell<T, Descriptor> tmpCell(cell);
    this->completePopulations(tmpCell);
    return this->getBaseDynamics().computeTemperature(tmpCell);
}

template <typename T, template <typename U> class Descriptor>
void BoundaryCompositeDynamics<T, Descriptor>::computePiNeq(
    Cell<T, Descriptor> const &cell, Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const
{
    Cell<T, Descriptor> tmpCell(cell);
    this->completePopulations(tmpCell);
    return this->getBaseDynamics().computePiNeq(tmpCell, PiNeq);
}

template <typename T, template <typename U> class Descriptor>
void BoundaryCompositeDynamics<T, Descriptor>::computeShearStress(
    Cell<T, Descriptor> const &cell, Array<T, SymmetricTensor<T, Descriptor>::n> &stress) const
{
    Cell<T, Descriptor> tmpCell(cell);
    this->completePopulations(tmpCell);
    return this->getBaseDynamics().computeShearStress(tmpCell, stress);
}

template <typename T, template <typename U> class Descriptor>
void BoundaryCompositeDynamics<T, Descriptor>::computeHeatFlux(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &q) const
{
    Cell<T, Descriptor> tmpCell(cell);
    this->completePopulations(tmpCell);
    this->getBaseDynamics().computeHeatFlux(tmpCell, q);
}

template <typename T, template <typename U> class Descriptor>
void BoundaryCompositeDynamics<T, Descriptor>::computeMoment(
    Cell<T, Descriptor> const &cell, plint momentId, T *moment) const
{
    Cell<T, Descriptor> tmpCell(cell);
    this->completePopulations(tmpCell);
    this->getBaseDynamics().computeMoment(tmpCell, momentId, moment);
}

template <typename T, template <typename U> class Descriptor>
T BoundaryCompositeDynamics<T, Descriptor>::computeRhoBar(Cell<T, Descriptor> const &cell) const
{
    Cell<T, Descriptor> tmpCell(cell);
    this->completePopulations(tmpCell);
    return this->getBaseDynamics().computeRhoBar(tmpCell);
}

template <typename T, template <typename U> class Descriptor>
void BoundaryCompositeDynamics<T, Descriptor>::computeRhoBarJ(
    Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j) const
{
    Cell<T, Descriptor> tmpCell(cell);
    this->completePopulations(tmpCell);
    this->getBaseDynamics().computeRhoBarJ(tmpCell, rhoBar, j);
}

template <typename T, template <typename U> class Descriptor>
void BoundaryCompositeDynamics<T, Descriptor>::computeRhoBarJPiNeq(
    Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
    Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const
{
    Cell<T, Descriptor> tmpCell(cell);
    this->completePopulations(tmpCell);
    this->getBaseDynamics().computeRhoBarJPiNeq(tmpCell, rhoBar, j, PiNeq);
}

template <typename T, template <typename U> class Descriptor>
T BoundaryCompositeDynamics<T, Descriptor>::computeEbar(Cell<T, Descriptor> const &cell) const
{
    Cell<T, Descriptor> tmpCell(cell);
    this->completePopulations(tmpCell);
    return this->getBaseDynamics().computeEbar(tmpCell);
}

/** Do nothing inside this functions. This defaults to a behavior where
 *  the dynamics of BoundaryCompositeDynamics is identical to baseDynamics.
 *  More interesting behavior is achieved in derived classes which overload
 *  method completePopulations().
 */
template <typename T, template <typename U> class Descriptor>
void BoundaryCompositeDynamics<T, Descriptor>::completePopulations(Cell<T, Descriptor> &cell) const
{ }

template <typename T, template <typename U> class Descriptor>
void BoundaryCompositeDynamics<T, Descriptor>::decompose(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData, plint order) const
{
    Cell<T, Descriptor> tmpCell(cell);
    this->completePopulations(tmpCell);
    this->getBaseDynamics().decompose(tmpCell, rawData, order);
}

/* *************** Class StoreDensityDynamics *********************** */

template <typename T, template <typename U> class Descriptor>
int StoreDensityDynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, StoreDensityDynamics<T, Descriptor> >(
        "Boundary_StoreDensity");

template <typename T, template <typename U> class Descriptor>
StoreDensityDynamics<T, Descriptor>::StoreDensityDynamics(
    Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_) :
    BoundaryCompositeDynamics<T, Descriptor>(baseDynamics_, automaticPrepareCollision_)
{
    rhoBar = Descriptor<T>::rhoBar((T)1);
}

template <typename T, template <typename U> class Descriptor>
StoreDensityDynamics<T, Descriptor>::StoreDensityDynamics(HierarchicUnserializer &unserializer) :
    BoundaryCompositeDynamics<T, Descriptor>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
StoreDensityDynamics<T, Descriptor> *StoreDensityDynamics<T, Descriptor>::clone() const
{
    return new StoreDensityDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int StoreDensityDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void StoreDensityDynamics<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    BoundaryCompositeDynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(rhoBar);
}

template <typename T, template <typename U> class Descriptor>
void StoreDensityDynamics<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    BoundaryCompositeDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(rhoBar);
}

template <typename T, template <typename U> class Descriptor>
T StoreDensityDynamics<T, Descriptor>::computeDensity(Cell<T, Descriptor> const &cell) const
{
    return Descriptor<T>::fullRho(rhoBar);
}

template <typename T, template <typename U> class Descriptor>
void StoreDensityDynamics<T, Descriptor>::defineDensity(Cell<T, Descriptor> &cell, T rho_)
{
    rhoBar = Descriptor<T>::rhoBar(rho_);
}

template <typename T, template <typename U> class Descriptor>
T StoreDensityDynamics<T, Descriptor>::computeRhoBar(Cell<T, Descriptor> const &cell) const
{
    return rhoBar;
}

template <typename T, template <typename U> class Descriptor>
void StoreDensityDynamics<T, Descriptor>::computeRhoBarJ(
    Cell<T, Descriptor> const &cell, T &rhoBar_, Array<T, Descriptor<T>::d> &j) const
{
    BoundaryCompositeDynamics<T, Descriptor>::computeRhoBarJ(cell, rhoBar_, j);
    rhoBar_ = rhoBar;
}

/* *************** Class StoreVelocityDynamics *********************** */

template <typename T, template <typename U> class Descriptor>
int StoreVelocityDynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, StoreVelocityDynamics<T, Descriptor> >(
        "Boundary_StoreVelocity");

template <typename T, template <typename U> class Descriptor>
StoreVelocityDynamics<T, Descriptor>::StoreVelocityDynamics(
    Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_) :
    BoundaryCompositeDynamics<T, Descriptor>(baseDynamics_, automaticPrepareCollision_)
{
    velocity.resetToZero();
}

template <typename T, template <typename U> class Descriptor>
StoreVelocityDynamics<T, Descriptor>::StoreVelocityDynamics(HierarchicUnserializer &unserializer) :
    BoundaryCompositeDynamics<T, Descriptor>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
StoreVelocityDynamics<T, Descriptor> *StoreVelocityDynamics<T, Descriptor>::clone() const
{
    return new StoreVelocityDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int StoreVelocityDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void StoreVelocityDynamics<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    PLB_PRECONDITION(unserializer.getId() == this->getId());
    BoundaryCompositeDynamics<T, Descriptor>::unserialize(unserializer);
    for (plint iDim = 0; iDim < Descriptor<T>::d; ++iDim) {
        unserializer.readValue(velocity[iDim]);
    }
}

template <typename T, template <typename U> class Descriptor>
void StoreVelocityDynamics<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    BoundaryCompositeDynamics<T, Descriptor>::serialize(serializer);
    for (plint iDim = 0; iDim < Descriptor<T>::d; ++iDim) {
        serializer.addValue(velocity[iDim]);
    }
}

template <typename T, template <typename U> class Descriptor>
void StoreVelocityDynamics<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &velocity_) const
{
    for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
        velocity_[iD] = velocity[iD];
    }
}

template <typename T, template <typename U> class Descriptor>
void StoreVelocityDynamics<T, Descriptor>::defineVelocity(
    Cell<T, Descriptor> &cell, Array<T, Descriptor<T>::d> const &velocity_)
{
    for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
        velocity[iD] = velocity_[iD];
    }
}

template <typename T, template <typename U> class Descriptor>
T StoreVelocityDynamics<T, Descriptor>::computeRhoBar(Cell<T, Descriptor> const &cell) const
{
    T rho = this->computeDensity(cell);
    return Descriptor<T>::rhoBar(rho);
}

template <typename T, template <typename U> class Descriptor>
void StoreVelocityDynamics<T, Descriptor>::computeRhoBarJ(
    Cell<T, Descriptor> const &cell, T &rhoBar_, Array<T, Descriptor<T>::d> &j) const
{
    T rho = this->computeDensity(cell);
    rhoBar_ = Descriptor<T>::rhoBar(rho);
    for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
        // Use the formula uLB = uP - 1/2 g. If there is no external force,
        //   the force term automatically evaluates to zero.
        j[iD] = rho * (velocity[iD] - 0.5 * getExternalForceComponent(cell, iD));
    }
}

/* *************** Class StoreDensityAndVelocityDynamics ************* */

template <typename T, template <typename U> class Descriptor>
int StoreDensityAndVelocityDynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, StoreDensityAndVelocityDynamics<T, Descriptor> >(
        "Boundary_StoreDensityAndVelocity");

template <typename T, template <typename U> class Descriptor>
StoreDensityAndVelocityDynamics<T, Descriptor>::StoreDensityAndVelocityDynamics(
    Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_) :
    BoundaryCompositeDynamics<T, Descriptor>(baseDynamics_, automaticPrepareCollision_)
{
    rhoBar = Descriptor<T>::rhoBar((T)1);
    velocity.resetToZero();
}

template <typename T, template <typename U> class Descriptor>
StoreDensityAndVelocityDynamics<T, Descriptor>::StoreDensityAndVelocityDynamics(
    HierarchicUnserializer &unserializer) :
    BoundaryCompositeDynamics<T, Descriptor>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
StoreDensityAndVelocityDynamics<T, Descriptor>
    *StoreDensityAndVelocityDynamics<T, Descriptor>::clone() const
{
    return new StoreDensityAndVelocityDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int StoreDensityAndVelocityDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void StoreDensityAndVelocityDynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    PLB_PRECONDITION(unserializer.getId() == this->getId());
    BoundaryCompositeDynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(rhoBar);
    for (plint iDim = 0; iDim < Descriptor<T>::d; ++iDim) {
        unserializer.readValue(velocity[iDim]);
    }
}

template <typename T, template <typename U> class Descriptor>
void StoreDensityAndVelocityDynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    BoundaryCompositeDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(rhoBar);
    for (plint iDim = 0; iDim < Descriptor<T>::d; ++iDim) {
        serializer.addValue(velocity[iDim]);
    }
}

template <typename T, template <typename U> class Descriptor>
T StoreDensityAndVelocityDynamics<T, Descriptor>::computeDensity(
    Cell<T, Descriptor> const &cell) const
{
    return Descriptor<T>::fullRho(rhoBar);
}

template <typename T, template <typename U> class Descriptor>
void StoreDensityAndVelocityDynamics<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &velocity_) const
{
    if (this->velIsJ()) {
        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            velocity_[iD] = Descriptor<T>::fullRho(rhoBar) * velocity[iD];
        }
    } else {
        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            velocity_[iD] = velocity[iD];
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void StoreDensityAndVelocityDynamics<T, Descriptor>::defineDensity(
    Cell<T, Descriptor> &cell, T rho_)
{
    rhoBar = Descriptor<T>::rhoBar(rho_);
}

template <typename T, template <typename U> class Descriptor>
void StoreDensityAndVelocityDynamics<T, Descriptor>::defineVelocity(
    Cell<T, Descriptor> &cell, Array<T, Descriptor<T>::d> const &velocity_)
{
    for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
        velocity[iD] = velocity_[iD];
    }
}

template <typename T, template <typename U> class Descriptor>
T StoreDensityAndVelocityDynamics<T, Descriptor>::computeRhoBar(
    Cell<T, Descriptor> const &cell) const
{
    return rhoBar;
}

template <typename T, template <typename U> class Descriptor>
void StoreDensityAndVelocityDynamics<T, Descriptor>::computeRhoBarJ(
    Cell<T, Descriptor> const &cell, T &rhoBar_, Array<T, Descriptor<T>::d> &j) const
{
    rhoBar_ = rhoBar;
    T rho = Descriptor<T>::fullRho(rhoBar);
    for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
        // Use the formula uLB = uP - 1/2 g. If there is no external force,
        //   the force term automatically evaluates to zero.
        j[iD] = rho * (velocity[iD] - 0.5 * getExternalForceComponent(cell, iD));
    }
}

/* *************** Class StoreTemperatureAndVelocityDynamics ************* */

template <typename T, template <typename U> class Descriptor>
int StoreTemperatureAndVelocityDynamics<T, Descriptor>::id = meta::registerGeneralDynamics<
    T, Descriptor, StoreTemperatureAndVelocityDynamics<T, Descriptor> >(
    "Boundary_StoreTemperatureAndVelocity");

template <typename T, template <typename U> class Descriptor>
StoreTemperatureAndVelocityDynamics<T, Descriptor>::StoreTemperatureAndVelocityDynamics(
    Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_) :
    BoundaryCompositeDynamics<T, Descriptor>(baseDynamics_, automaticPrepareCollision_)
{
    thetaBar = T();
    velocity.resetToZero();
}

template <typename T, template <typename U> class Descriptor>
StoreTemperatureAndVelocityDynamics<T, Descriptor>::StoreTemperatureAndVelocityDynamics(
    HierarchicUnserializer &unserializer) :
    BoundaryCompositeDynamics<T, Descriptor>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
StoreTemperatureAndVelocityDynamics<T, Descriptor>
    *StoreTemperatureAndVelocityDynamics<T, Descriptor>::clone() const
{
    return new StoreTemperatureAndVelocityDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int StoreTemperatureAndVelocityDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void StoreTemperatureAndVelocityDynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    PLB_PRECONDITION(unserializer.getId() == this->getId());
    BoundaryCompositeDynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(thetaBar);
    for (plint iDim = 0; iDim < Descriptor<T>::d; ++iDim) {
        unserializer.readValue(velocity[iDim]);
    }
}

template <typename T, template <typename U> class Descriptor>
void StoreTemperatureAndVelocityDynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    BoundaryCompositeDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(thetaBar);
    for (plint iDim = 0; iDim < Descriptor<T>::d; ++iDim) {
        serializer.addValue(velocity[iDim]);
    }
}

template <typename T, template <typename U> class Descriptor>
void StoreTemperatureAndVelocityDynamics<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &velocity_) const
{
    for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
        velocity_[iD] = velocity[iD];
    }
}

template <typename T, template <typename U> class Descriptor>
T StoreTemperatureAndVelocityDynamics<T, Descriptor>::computeTemperature(
    Cell<T, Descriptor> const &cell) const
{
    return (thetaBar + (T)1) * Descriptor<T>::cs2;
}

template <typename T, template <typename U> class Descriptor>
void StoreTemperatureAndVelocityDynamics<T, Descriptor>::defineVelocity(
    Cell<T, Descriptor> &cell, Array<T, Descriptor<T>::d> const &velocity_)
{
    for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
        velocity[iD] = velocity_[iD];
    }
}

template <typename T, template <typename U> class Descriptor>
void StoreTemperatureAndVelocityDynamics<T, Descriptor>::defineTemperature(
    Cell<T, Descriptor> &cell, T theta_)
{
    thetaBar = theta_ * Descriptor<T>::invCs2 - (T)1;
}

template <typename T, template <typename U> class Descriptor>
T StoreTemperatureAndVelocityDynamics<T, Descriptor>::computeRhoBar(
    Cell<T, Descriptor> const &cell) const
{
    T rho = this->computeDensity(cell);
    return Descriptor<T>::rhoBar(rho);
}

template <typename T, template <typename U> class Descriptor>
void StoreTemperatureAndVelocityDynamics<T, Descriptor>::computeRhoBarJ(
    Cell<T, Descriptor> const &cell, T &rhoBar_, Array<T, Descriptor<T>::d> &j) const
{
    T rho = this->computeDensity(cell);
    rhoBar_ = Descriptor<T>::rhoBar(rho);
    for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
        // Use the formula uLB = uP - 1/2 g. If there is no external force,
        //   the force term automatically evaluates to zero.
        j[iD] = rho * (velocity[iD] - 0.5 * getExternalForceComponent(cell, iD));
    }
}

/* *************** Class VelocityDirichletBoundaryDynamics ************* */

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int VelocityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>::id =
    meta::registerGeneralDynamics<
        T, Descriptor, VelocityDirichletBoundaryDynamics<T, Descriptor, direction, orientation> >(
        std::string("Boundary_VelocityDirichlet_") + util::val2str(direction) + std::string("_")
        + util::val2str(orientation));

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
VelocityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>::
    VelocityDirichletBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_) :
    StoreVelocityDynamics<T, Descriptor>(baseDynamics_, automaticPrepareCollision_)
{ }

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
VelocityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>::
    VelocityDirichletBoundaryDynamics(HierarchicUnserializer &unserializer) :
    StoreVelocityDynamics<T, Descriptor>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
VelocityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>
    *VelocityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>::clone() const
{
    return new VelocityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>(*this);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int VelocityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
T VelocityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>::computeRhoBar(
    Cell<T, Descriptor> const &cell) const
{
    std::vector<plint> const &onWallIndices =
        indexTemplates::subIndex<Descriptor<T>, direction, 0>();

    std::vector<plint> const &normalIndices =
        indexTemplates::subIndex<Descriptor<T>, direction, orientation>();

    T rhoOnWall = T();
    for (pluint fIndex = 0; fIndex < onWallIndices.size(); ++fIndex) {
        rhoOnWall += cell[onWallIndices[fIndex]];
    }

    T rhoNormal = T();
    for (pluint fIndex = 0; fIndex < normalIndices.size(); ++fIndex) {
        rhoNormal += cell[normalIndices[fIndex]];
    }

    T velNormal = (T)orientation
                  * (this->velocity[direction] - 0.5 * getExternalForceComponent(cell, direction));
    if (this->velIsJ()) {
        T rhoBar = (T)2 * rhoNormal + rhoOnWall - velNormal;
        return rhoBar;
    } else {
        T rhoBar = ((T)2 * rhoNormal + rhoOnWall - Descriptor<T>::SkordosFactor() * velNormal)
                   / ((T)1 + velNormal);
        return rhoBar;
    }
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
T VelocityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>::computeDensity(
    Cell<T, Descriptor> const &cell) const
{
    return Descriptor<T>::fullRho(computeRhoBar(cell));
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void VelocityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>::computeRhoBarJ(
    Cell<T, Descriptor> const &cell, T &rhoBar_, Array<T, Descriptor<T>::d> &j) const
{
    rhoBar_ = computeRhoBar(cell);
    T rho = Descriptor<T>::fullRho(rhoBar_);
    if (this->velIsJ()) {
        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            // Use the formula uLB = uP - 1/2 g. If there is no external force,
            //   the force term automatically evaluates to zero.
            j[iD] = this->velocity[iD] - 0.5 * getExternalForceComponent(cell, iD);
        }
    } else {
        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            // Use the formula uLB = uP - 1/2 g. If there is no external force,
            //   the force term automatically evaluates to zero.
            j[iD] = rho * (this->velocity[iD] - 0.5 * getExternalForceComponent(cell, iD));
        }
    }
}

/* *************** Class VelocityDirichletConstRhoBoundaryDynamics ************* */

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int VelocityDirichletConstRhoBoundaryDynamics<T, Descriptor, direction, orientation>::id =
    meta::registerGeneralDynamics<
        T, Descriptor,
        VelocityDirichletConstRhoBoundaryDynamics<T, Descriptor, direction, orientation> >(
        std::string("Boundary_VelocityDirichletConstRho_") + util::val2str(direction)
        + std::string("_") + util::val2str(orientation));

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
VelocityDirichletConstRhoBoundaryDynamics<T, Descriptor, direction, orientation>::
    VelocityDirichletConstRhoBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_) :
    StoreVelocityDynamics<T, Descriptor>(baseDynamics_, automaticPrepareCollision_)
{ }

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
VelocityDirichletConstRhoBoundaryDynamics<T, Descriptor, direction, orientation>::
    VelocityDirichletConstRhoBoundaryDynamics(HierarchicUnserializer &unserializer) :
    StoreVelocityDynamics<T, Descriptor>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
VelocityDirichletConstRhoBoundaryDynamics<T, Descriptor, direction, orientation>
    *VelocityDirichletConstRhoBoundaryDynamics<T, Descriptor, direction, orientation>::clone() const
{
    return new VelocityDirichletConstRhoBoundaryDynamics<T, Descriptor, direction, orientation>(
        *this);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int VelocityDirichletConstRhoBoundaryDynamics<T, Descriptor, direction, orientation>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
T VelocityDirichletConstRhoBoundaryDynamics<T, Descriptor, direction, orientation>::computeRhoBar(
    Cell<T, Descriptor> const &cell) const
{
    return Descriptor<T>::rhoBar((T)1.0);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
T VelocityDirichletConstRhoBoundaryDynamics<T, Descriptor, direction, orientation>::computeDensity(
    Cell<T, Descriptor> const &cell) const
{
    return (T)1.0;
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void VelocityDirichletConstRhoBoundaryDynamics<T, Descriptor, direction, orientation>::
    computeRhoBarJ(Cell<T, Descriptor> const &cell, T &rhoBar_, Array<T, Descriptor<T>::d> &j) const
{
    rhoBar_ = computeRhoBar(cell);
    T rho = Descriptor<T>::fullRho(rhoBar_);
    if (this->velIsJ()) {
        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            // Use the formula uLB = uP - 1/2 g. If there is no external force,
            //   the force term automatically evaluates to zero.
            j[iD] = this->velocity[iD] - 0.5 * getExternalForceComponent(cell, iD);
        }
    } else {
        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            // Use the formula uLB = uP - 1/2 g. If there is no external force,
            //   the force term automatically evaluates to zero.
            j[iD] = rho * (this->velocity[iD] - 0.5 * getExternalForceComponent(cell, iD));
        }
    }
}

/* *************** Class DensityDirichletBoundaryDynamics ************* */

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>::id =
    meta::registerGeneralDynamics<
        T, Descriptor, DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation> >(
        std::string("Boundary_DensityDirichlet_") + util::val2str(direction) + std::string("_")
        + util::val2str(orientation));

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>::
    DensityDirichletBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_) :
    StoreDensityDynamics<T, Descriptor>(baseDynamics_, automaticPrepareCollision_)
{ }

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>::
    DensityDirichletBoundaryDynamics(HierarchicUnserializer &unserializer) :
    StoreDensityDynamics<T, Descriptor>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>
    *DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>::clone() const
{
    return new DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>(*this);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>::computeJ(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &j_) const
{
    // All velocity components parallel to the wall are zero by definition.
    for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
        j_[iD] = T();
    }
    T rhoBar = this->rhoBar;

    std::vector<plint> const &onWallIndices =
        indexTemplates::subIndex<Descriptor<T>, direction, 0>();

    std::vector<plint> const &normalIndices =
        indexTemplates::subIndex<Descriptor<T>, direction, orientation>();

    T rhoOnWall = T();
    for (pluint fIndex = 0; fIndex < onWallIndices.size(); ++fIndex) {
        rhoOnWall += cell[onWallIndices[fIndex]];
    }

    T rhoNormal = T();
    for (pluint fIndex = 0; fIndex < normalIndices.size(); ++fIndex) {
        rhoNormal += cell[normalIndices[fIndex]];
    }

    j_[direction] = (T)orientation * ((T)2 * rhoNormal + rhoOnWall - rhoBar);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &velocity_) const
{
    this->computeJ(cell, velocity_);
    if (this->velIsJ()) {
        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            // Use the formula uLB = uP - 1/2 g. If there is no external force,
            //   the force term automatically evaluates to zero.
            velocity_[iD] += 0.5 * getExternalForceComponent(cell, iD);
        }
    } else {
        T invRho = Descriptor<T>::invRho(this->rhoBar);
        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            // Use the formula uLB = uP - 1/2 g. If there is no external force,
            //   the force term automatically evaluates to zero.
            velocity_[iD] *= invRho;
            velocity_[iD] += 0.5 * getExternalForceComponent(cell, iD);
        }
    }
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>::computeRhoBarJ(
    Cell<T, Descriptor> const &cell, T &rhoBar_, Array<T, Descriptor<T>::d> &j) const
{
    rhoBar_ = this->rhoBar;
    this->computeJ(cell, j);
}

}  // namespace plb

#endif  // BOUNDARY_DYNAMICS_HH
