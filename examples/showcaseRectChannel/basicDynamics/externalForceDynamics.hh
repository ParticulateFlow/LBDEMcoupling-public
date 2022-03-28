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
 * Collision terms with external force -- generic code.
 */
#ifndef EXTERNAL_FORCE_DYNAMICS_HH
#define EXTERNAL_FORCE_DYNAMICS_HH

#include <algorithm>
#include <limits>

#include "basicDynamics/externalForceDynamics.h"
#include "basicDynamics/isoThermalDynamics.h"
#include "complexDynamics/smagorinskyDynamics.hh"
#include "core/cell.h"
#include "core/dynamicsIdentifiers.h"
#include "core/latticeStatistics.h"
#include "latticeBoltzmann/d3q13Templates.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/externalForceTemplates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/offEquilibriumTemplates.h"

namespace plb {

/* *************** Class ExternalForceDynamics *********************************************** */

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
ExternalForceDynamics<T, Descriptor>::ExternalForceDynamics(T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
void ExternalForceDynamics<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> force, j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));

    T invRho = Descriptor<T>::invRho(rhoBar);
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        u[iD] = j[iD] * invRho + force[iD] / (T)2;
    }
}

template <typename T, template <typename U> class Descriptor>
void ExternalForceDynamics<T, Descriptor>::computeVelocityExternal(
    Cell<T, Descriptor> const &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j,
    Array<T, Descriptor<T>::d> &u) const
{
    Array<T, Descriptor<T>::d> force;
    force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));

    T invRho = Descriptor<T>::invRho(rhoBar);
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        u[iD] = j[iD] * invRho + force[iD] / (T)2;
    }
}

/* *************** Class IncExternalForceDynamics *********************************************** */

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
IncExternalForceDynamics<T, Descriptor>::IncExternalForceDynamics(T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
void IncExternalForceDynamics<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
{
    Array<T, Descriptor<T>::d> force, j;
    momentTemplates<T, Descriptor>::get_j(cell, j);
    force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));
    u = j + (T)0.5 * force;
}

template <typename T, template <typename U> class Descriptor>
void IncExternalForceDynamics<T, Descriptor>::computeVelocityExternal(
    Cell<T, Descriptor> const &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j,
    Array<T, Descriptor<T>::d> &u) const
{
    Array<T, Descriptor<T>::d> force;
    force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));
    u = j + (T)0.5 * force;
}

template <typename T, template <typename U> class Descriptor>
void IncExternalForceDynamics<T, Descriptor>::computeRhoBarJPiNeq(
    Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
    Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const
{
    T invRho0 = (T)1.;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq, invRho0);
}

/* *************** Class NaiveExternalForceBGKdynamics ********************************** */

template <typename T, template <typename U> class Descriptor>
int NaiveExternalForceBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, NaiveExternalForceBGKdynamics<T, Descriptor> >(
        "BGK_ExternalForce_Naive");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
NaiveExternalForceBGKdynamics<T, Descriptor>::NaiveExternalForceBGKdynamics(T omega_) :
    ExternalForceDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
NaiveExternalForceBGKdynamics<T, Descriptor>::NaiveExternalForceBGKdynamics(
    HierarchicUnserializer &unserializer) :
    ExternalForceDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
NaiveExternalForceBGKdynamics<T, Descriptor> *NaiveExternalForceBGKdynamics<T, Descriptor>::clone()
    const
{
    return new NaiveExternalForceBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int NaiveExternalForceBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void NaiveExternalForceBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, j, this->getOmega());
    externalForceTemplates<T, Descriptor>::addNaiveForce(cell);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void NaiveExternalForceBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, j, this->getOmega());
    externalForceTemplates<T, Descriptor>::addNaiveForce(cell);

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T NaiveExternalForceBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class NaiveExternalForcePrecondBGKdynamics ********************************** */

template <typename T, template <typename U> class Descriptor>
int NaiveExternalForcePrecondBGKdynamics<T, Descriptor>::id = meta::registerGeneralDynamics<
    T, Descriptor, NaiveExternalForcePrecondBGKdynamics<T, Descriptor> >(
    "Precond_BGK_ExternalForce_Naive");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
NaiveExternalForcePrecondBGKdynamics<T, Descriptor>::NaiveExternalForcePrecondBGKdynamics(
    T omega_, T invGamma_) :
    ExternalForceDynamics<T, Descriptor>(omega_), invGamma(invGamma_)
{ }

template <typename T, template <typename U> class Descriptor>
NaiveExternalForcePrecondBGKdynamics<T, Descriptor>::NaiveExternalForcePrecondBGKdynamics(
    HierarchicUnserializer &unserializer) :
    ExternalForceDynamics<T, Descriptor>(T()), invGamma(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
NaiveExternalForcePrecondBGKdynamics<T, Descriptor>
    *NaiveExternalForcePrecondBGKdynamics<T, Descriptor>::clone() const
{
    return new NaiveExternalForcePrecondBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int NaiveExternalForcePrecondBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void NaiveExternalForcePrecondBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T uSqr = dynamicsTemplates<T, Descriptor>::precond_bgk_ma2_collision(
        cell, rhoBar, j, this->getOmega(), invGamma);
    externalForceTemplates<T, Descriptor>::addNaiveForce(cell);
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void NaiveExternalForcePrecondBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T uSqr = dynamicsTemplates<T, Descriptor>::precond_bgk_ma2_collision(
        cell, rhoBar, j, this->getOmega(), invGamma);
    externalForceTemplates<T, Descriptor>::addNaiveForce(cell);

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T NaiveExternalForcePrecondBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::precond_bgk_ma2_equilibrium(
        iPop, rhoBar, invRho, j, jSqr, invGamma);
}

template <typename T, template <typename U> class Descriptor>
void NaiveExternalForcePrecondBGKdynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    ExternalForceDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(invGamma);
}

template <typename T, template <typename U> class Descriptor>
void NaiveExternalForcePrecondBGKdynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    ExternalForceDynamics<T, Descriptor>::unserialize(unserializer);
    invGamma = unserializer.readValue<T>();
}

/* *************** Class GuoExternalForceBGKdynamics ********************************** */

template <typename T, template <typename U> class Descriptor>
int GuoExternalForceBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, GuoExternalForceBGKdynamics<T, Descriptor> >(
        "BGK_ExternalForce_Guo");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
GuoExternalForceBGKdynamics<T, Descriptor>::GuoExternalForceBGKdynamics(T omega_) :
    ExternalForceDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
GuoExternalForceBGKdynamics<T, Descriptor>::GuoExternalForceBGKdynamics(
    HierarchicUnserializer &unserializer) :
    ExternalForceDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
GuoExternalForceBGKdynamics<T, Descriptor> *GuoExternalForceBGKdynamics<T, Descriptor>::clone()
    const
{
    return new GuoExternalForceBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int GuoExternalForceBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar = this->computeRhoBar(cell);
    Array<T, Descriptor<T>::d> u, j;
    this->computeVelocity(cell, u);
    T rho = Descriptor<T>::fullRho(rhoBar);
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        j[iD] = rho * u[iD];
    }

    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, j, this->getOmega());
    externalForceTemplates<T, Descriptor>::addGuoForce(cell, u, this->getOmega(), (T)1);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    Array<T, Descriptor<T>::d> u;
    this->computeVelocityExternal(cell, rhoBar, j, u);
    T rho = Descriptor<T>::fullRho(rhoBar);
    Array<T, Descriptor<T>::d> newJ;
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        newJ[iD] = rho * u[iD];
    }

    T uSqr =
        dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, newJ, this->getOmega());
    externalForceTemplates<T, Descriptor>::addGuoForce(cell, u, this->getOmega(), (T)1);

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T GuoExternalForceBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class GuoExternalForceCompleteRegularizedBGKdynamics
 * ********************************** */

template <typename T, template <typename U> class Descriptor>
int GuoExternalForceCompleteRegularizedBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<
        T, Descriptor, GuoExternalForceCompleteRegularizedBGKdynamics<T, Descriptor> >(
        "RegularizedCompleteBGK_ExternalForce_Guo");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
GuoExternalForceCompleteRegularizedBGKdynamics<
    T, Descriptor>::GuoExternalForceCompleteRegularizedBGKdynamics(T omega_) :
    ExternalForceDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
GuoExternalForceCompleteRegularizedBGKdynamics<T, Descriptor>::
    GuoExternalForceCompleteRegularizedBGKdynamics(HierarchicUnserializer &unserializer) :
    ExternalForceDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
GuoExternalForceCompleteRegularizedBGKdynamics<T, Descriptor>
    *GuoExternalForceCompleteRegularizedBGKdynamics<T, Descriptor>::clone() const
{
    return new GuoExternalForceCompleteRegularizedBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int GuoExternalForceCompleteRegularizedBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceCompleteRegularizedBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar = this->computeRhoBar(cell);
    Array<T, Descriptor<T>::d> u, j;
    this->computeVelocity(cell, u);
    T rho = Descriptor<T>::fullRho(rhoBar);
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        j[iD] = rho * u[iD];
    }

    T rhoBarLb;
    Array<T, Descriptor<T>::d> jLb;
    Array<T, SymmetricTensor<T, Descriptor>::n> piNeqLb;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBarLb, jLb, piNeqLb);
    T jSqrLb = VectorTemplate<T, Descriptor>::normSqr(jLb);
    // T invRhoLb = Descriptor<T>::invRho(rhoBarLb);
    regularize(cell, rhoBarLb, jLb, jSqrLb, piNeqLb);

    T uSqr = dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_collision(
        cell, rhoBar, Descriptor<T>::invRho(rhoBar), j, this->getOmega());

    externalForceTemplates<T, Descriptor>::addGuoForce(cell, u, this->getOmega(), (T)1);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceCompleteRegularizedBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    Array<T, Descriptor<T>::d> u;
    this->computeVelocityExternal(cell, rhoBar, j, u);
    T rho = Descriptor<T>::fullRho(rhoBar);
    Array<T, Descriptor<T>::d> newJ;
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        newJ[iD] = rho * u[iD];
    }

    T rhoBarLb;
    Array<T, Descriptor<T>::d> jLb;
    Array<T, SymmetricTensor<T, Descriptor>::n> piNeqLb;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBarLb, jLb, piNeqLb);
    T jSqrLb = VectorTemplate<T, Descriptor>::normSqr(jLb);
    // T invRhoLb = Descriptor<T>::invRho(rhoBarLb);
    regularize(cell, rhoBarLb, jLb, jSqrLb, piNeqLb);

    T uSqr = dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_collision(
        cell, rhoBar, Descriptor<T>::invRho(rhoBar), newJ, this->getOmega());

    externalForceTemplates<T, Descriptor>::addGuoForce(cell, u, this->getOmega(), (T)1);

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceCompleteRegularizedBGKdynamics<T, Descriptor>::computeEquilibria(
    Array<T, Descriptor<T>::q> &fEq, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_equilibria(rhoBar, invRho, j, jSqr, fEq);
}

template <typename T, template <typename U> class Descriptor>
T GuoExternalForceCompleteRegularizedBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_equilibrium(
        iPop, rhoBar, invRho, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceCompleteRegularizedBGKdynamics<T, Descriptor>::regularize(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_regularize(
        cell, rhoBar, invRho, j, jSqr, PiNeq, this->getOmega(), this->getOmega());
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceCompleteRegularizedBGKdynamics<T, Descriptor>::decomposeOrder0(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    rawData[0] = rhoBar;
    j.to_cArray(&rawData[1]);

    Array<T, Descriptor<T>::q> fEq;
    dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_equilibria(
        rhoBar, Descriptor<T>::invRho(rhoBar), j, jSqr, fEq);

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        rawData[1 + Descriptor<T>::d + iPop] = cell[iPop] - fEq[iPop];
    }

    int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        rawData[offset + iExt] = *cell.getExternal(iExt);
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceCompleteRegularizedBGKdynamics<T, Descriptor>::recomposeOrder0(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    T rhoBar = rawData[0];
    Array<T, Descriptor<T>::d> j;
    j.from_cArray(&rawData[1]);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    Array<T, Descriptor<T>::q> fEq;
    dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_equilibria(
        rhoBar, Descriptor<T>::invRho(rhoBar), j, jSqr, fEq);

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] = fEq[iPop] + rawData[1 + Descriptor<T>::d + iPop];
    }

    int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *cell.getExternal(iExt) = rawData[offset + iExt];
    }
}

/* *************** Class GuoExternalForceConsistentSmagorinskyCompleteRegularizedBGKdynamics
 * ********************************** */

template <typename T, template <typename U> class Descriptor>
int GuoExternalForceConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<
        T, Descriptor,
        GuoExternalForceConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor> >(
        "RegularizedCompleteBGK_Consitent_Smagorinsky_ExternalForce_Guo");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor>::
    GuoExternalForceConsistentSmagorinskyCompleteRegularizedBGKdynamics(T omega_, T cSmago_) :
    GuoExternalForceCompleteRegularizedBGKdynamics<T, Descriptor>(omega_), cSmago(cSmago_)
{ }

template <typename T, template <typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor>::
    GuoExternalForceConsistentSmagorinskyCompleteRegularizedBGKdynamics(
        HierarchicUnserializer &unserializer) :
    GuoExternalForceCompleteRegularizedBGKdynamics<T, Descriptor>(T()), cSmago(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor>
    *GuoExternalForceConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor>::clone()
        const
{
    return new GuoExternalForceConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor>(
        *this);
}

template <typename T, template <typename U> class Descriptor>
int GuoExternalForceConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor>::getId()
    const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    GuoExternalForceCompleteRegularizedBGKdynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceConsistentSmagorinskyCompleteRegularizedBGKdynamics<
    T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    GuoExternalForceCompleteRegularizedBGKdynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
T GuoExternalForceConsistentSmagorinskyCompleteRegularizedBGKdynamics<
    T, Descriptor>::getDynamicParameter(plint whichParameter, Cell<T, Descriptor> const &cell) const
{
    if (whichParameter == dynamicParams::dynamicOmega) {
        T rhoBar;
        Array<T, Descriptor<T>::d> j;
        Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
        momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
        T preFactor = SmagoOperations<T, Descriptor>::computePrefactor(this->getOmega(), cSmago);
        T omega = SmagoOperations<T, Descriptor>::computeOmega(
            this->getOmega(), preFactor, rhoBar, PiNeq);
        return omega;
    } else if (whichParameter == dynamicParams::smagorinskyConstant) {
        return cSmago;
    } else {
        return IsoThermalBulkDynamics<T, Descriptor>::getParameter(whichParameter);
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar = this->computeRhoBar(cell);
    Array<T, Descriptor<T>::d> u, j;
    this->computeVelocity(cell, u);
    T rho = Descriptor<T>::fullRho(rhoBar);
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        j[iD] = rho * u[iD];
    }

    T rhoBarLb;
    Array<T, Descriptor<T>::d> jLb;
    Array<T, SymmetricTensor<T, Descriptor>::n> piNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBarLb, jLb, piNeq);
    T jSqrLb = VectorTemplate<T, Descriptor>::normSqr(jLb);
    // T invRhoLb = Descriptor<T>::invRho(rhoBarLb);
    this->regularize(cell, rhoBarLb, jLb, jSqrLb, piNeq);

    T uSqr = dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_collision(
        cell, rhoBar, Descriptor<T>::invRho(rhoBar), j, this->getOmega());

    externalForceTemplates<T, Descriptor>::addGuoForce(cell, u, this->getOmega(), (T)1);

    T tau = (T)1 / this->getOmega();

    T normPiNeq = std::sqrt((T)2 * SymmetricTensor<T, Descriptor>::tensorNormSqr(piNeq));
    normPiNeq *= (T)2 * rho * util::sqr(Descriptor<T>::cs2 * cSmago * cSmago);

    Array<T, SymmetricTensor<T, Descriptor>::n> S;
    S.resetToZero();
    if (normPiNeq != T()) {  // test to avoid division per 0
        S = (-rho * tau * Descriptor<T>::cs2
             + std::sqrt(util::sqr(rho * tau * Descriptor<T>::cs2) + normPiNeq))
            / normPiNeq * piNeq;
    }
    T sNorm = std::sqrt((T)2 * SymmetricTensor<T, Descriptor>::tensorNormSqr(S));
    T preFactor = (T)0.5 * Descriptor<T>::invCs2 * Descriptor<T>::invCs2 * this->getOmega() * cSmago
                  * cSmago * sNorm;

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        Array<T, SymmetricTensor<T, Descriptor>::n> H2 =
            HermiteTemplate<T, Descriptor>::order2(iPop);
        T H_S = SymmetricTensor<T, Descriptor>::contractIndexes(H2, S);
        cell[iPop] += Descriptor<T>::t[iPop] * preFactor * H_S;
    }

    // T rhoBar = this->computeRhoBar(cell);
    // Array<T,Descriptor<T>::d> u;
    // this->computeVelocity(cell, u);

    // T uSqr = dynamicsTemplates<T,Descriptor>::
    //     consistent_smagorinsky_complete_regularized_mrt_ma2_collision(cell, 6,  cSmago,
    //     this->getOmega(), this->getOmega());

    externalForceTemplates<T, Descriptor>::addGuoForce(cell, u, this->getOmega(), (T)1);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor>::
    collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat)
{
    PLB_ASSERT(false);
    // Array<T,Descriptor<T>::d> u;
    // this->computeVelocityExternal(cell, rhoBar, j, u);
    // T rho = Descriptor<T>::fullRho(rhoBar);
    // Array<T,Descriptor<T>::d> newJ;
    // for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
    // {
    //     newJ[iD] = rho * u[iD];
    // }

    // T rhoBarLb;
    // Array<T,Descriptor<T>::d> jLb;
    // Array<T,SymmetricTensor<T,Descriptor>::n> piNeqLb;
    // momentTemplates<T,Descriptor>::compute_rhoBar_j_PiNeq(cell,rhoBarLb,jLb,piNeqLb);
    // T jSqrLb = VectorTemplate<T,Descriptor>::normSqr(jLb);
    // // T invRhoLb = Descriptor<T>::invRho(rhoBarLb);
    // this->regularize(cell,rhoBarLb,jLb,jSqrLb,piNeqLb);

    // T uSqr = dynamicsTemplates<T,Descriptor>::complete_bgk_ma2_collision(cell, rhoBar,
    // Descriptor<T>::invRho(rhoBar),
    //         newJ, this->getOmega());

    // externalForceTemplates<T,Descriptor>::addGuoForce(cell, u, this->getOmega(), (T)1);

    // if (cell.takesStatistics()) {
    //     gatherStatistics(stat, rhoBar, uSqr);
    // }
}

/* *************** Class ShanChenExternalForceBGKdynamics ********************************** */

template <typename T, template <typename U> class Descriptor>
int ShanChenExternalForceBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, ShanChenExternalForceBGKdynamics<T, Descriptor> >(
        "BGK_ExternalForce_ShanChen");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
ShanChenExternalForceBGKdynamics<T, Descriptor>::ShanChenExternalForceBGKdynamics(T omega_) :
    ExternalForceDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
ShanChenExternalForceBGKdynamics<T, Descriptor>::ShanChenExternalForceBGKdynamics(
    HierarchicUnserializer &unserializer) :
    ExternalForceDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
ShanChenExternalForceBGKdynamics<T, Descriptor>
    *ShanChenExternalForceBGKdynamics<T, Descriptor>::clone() const
{
    return new ShanChenExternalForceBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int ShanChenExternalForceBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void ShanChenExternalForceBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);

    T invOmega = 1. / this->getOmega();
    Array<T, Descriptor<T>::d> force;
    force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));
    Array<T, Descriptor<T>::d> jCorrected(j + invOmega * Descriptor<T>::fullRho(rhoBar) * force);

    dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, jCorrected, this->getOmega());
    if (cell.takesStatistics()) {
        T uSqr = T();
        T invRho = Descriptor<T>::invRho(rhoBar);
        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            uSqr += util::sqr(j[iD] * invRho + 0.5 * force[iD]);
        }
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void ShanChenExternalForceBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T invOmega = 1. / this->getOmega();
    Array<T, Descriptor<T>::d> force;
    force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));
    Array<T, Descriptor<T>::d> jCorrected(j + invOmega * Descriptor<T>::fullRho(rhoBar) * force);

    dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, jCorrected, this->getOmega());
    if (cell.takesStatistics()) {
        T uSqr = T();
        T invRho = Descriptor<T>::invRho(rhoBar);
        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            uSqr += util::sqr(j[iD] * invRho + 0.5 * force[iD]);
        }
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T ShanChenExternalForceBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class HeExternalForceBGKdynamics ********************************** */

template <typename T, template <typename U> class Descriptor>
int HeExternalForceBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, HeExternalForceBGKdynamics<T, Descriptor> >(
        "BGK_ExternalForce_He");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
HeExternalForceBGKdynamics<T, Descriptor>::HeExternalForceBGKdynamics(T omega_) :
    ExternalForceDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
HeExternalForceBGKdynamics<T, Descriptor>::HeExternalForceBGKdynamics(
    HierarchicUnserializer &unserializer) :
    ExternalForceDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
HeExternalForceBGKdynamics<T, Descriptor> *HeExternalForceBGKdynamics<T, Descriptor>::clone() const
{
    return new HeExternalForceBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int HeExternalForceBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void HeExternalForceBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar = this->computeRhoBar(cell);
    Array<T, Descriptor<T>::d> u, j;
    this->computeVelocity(cell, u);
    T rho = Descriptor<T>::fullRho(rhoBar);
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        j[iD] = rho * u[iD];
    }

    T uSqr = externalForceTemplates<T, Descriptor>::heForcedBGKCollision(
        cell, rhoBar, j, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void HeExternalForceBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    Array<T, Descriptor<T>::d> u;
    this->computeVelocityExternal(cell, rhoBar, j, u);
    T rho = Descriptor<T>::fullRho(rhoBar);
    Array<T, Descriptor<T>::d> newJ;
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        newJ[iD] = rho * u[iD];
    }

    T uSqr = externalForceTemplates<T, Descriptor>::heForcedBGKCollision(
        cell, rhoBar, newJ, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T HeExternalForceBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class IncGuoExternalForceBGKdynamics ********************************** */

template <typename T, template <typename U> class Descriptor>
int IncGuoExternalForceBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, IncGuoExternalForceBGKdynamics<T, Descriptor> >(
        "IncBGK_ExternalForce_Guo");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
IncGuoExternalForceBGKdynamics<T, Descriptor>::IncGuoExternalForceBGKdynamics(T omega_) :
    IncExternalForceDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
IncGuoExternalForceBGKdynamics<T, Descriptor>::IncGuoExternalForceBGKdynamics(
    HierarchicUnserializer &unserializer) :
    IncExternalForceDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
IncGuoExternalForceBGKdynamics<T, Descriptor>
    *IncGuoExternalForceBGKdynamics<T, Descriptor>::clone() const
{
    return new IncGuoExternalForceBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int IncGuoExternalForceBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void IncGuoExternalForceBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar = this->computeRhoBar(cell);
    Array<T, Descriptor<T>::d> u, j;
    this->computeVelocity(cell, u);
    T rho = Descriptor<T>::fullRho(rhoBar);
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        j[iD] = rho * u[iD];
    }

    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_inc_collision(cell, rhoBar, j, this->getOmega());
    externalForceTemplates<T, Descriptor>::addGuoForce(cell, u, this->getOmega(), (T)1);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void IncGuoExternalForceBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    Array<T, Descriptor<T>::d> u;
    this->computeVelocityExternal(cell, rhoBar, j, u);
    T rho = Descriptor<T>::fullRho(rhoBar);
    Array<T, Descriptor<T>::d> newJ;
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        newJ[iD] = rho * u[iD];
    }

    T uSqr =
        dynamicsTemplates<T, Descriptor>::bgk_inc_collision(cell, rhoBar, newJ, this->getOmega());
    externalForceTemplates<T, Descriptor>::addGuoForce(cell, u, this->getOmega(), (T)1);

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T IncGuoExternalForceBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    // Inc => invRho = 1.0;
    T invRho = (T)1;
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
bool IncGuoExternalForceBGKdynamics<T, Descriptor>::velIsJ() const
{
    return true;
}

/* *************** Class ShanChenExternalForceRegularizedBGKdynamics
 * ************************************ */

template <typename T, template <typename U> class Descriptor>
int ShanChenExternalForceRegularizedBGKdynamics<T, Descriptor>::id = meta::registerGeneralDynamics<
    T, Descriptor, ShanChenExternalForceRegularizedBGKdynamics<T, Descriptor> >(
    "Regularized_BGK_ExternalForce_ShanChen");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
ShanChenExternalForceRegularizedBGKdynamics<
    T, Descriptor>::ShanChenExternalForceRegularizedBGKdynamics(T omega_) :
    ExternalForceDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
ShanChenExternalForceRegularizedBGKdynamics<T, Descriptor>::
    ShanChenExternalForceRegularizedBGKdynamics(HierarchicUnserializer &unserializer) :
    ExternalForceDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
ShanChenExternalForceRegularizedBGKdynamics<T, Descriptor>
    *ShanChenExternalForceRegularizedBGKdynamics<T, Descriptor>::clone() const
{
    return new ShanChenExternalForceRegularizedBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int ShanChenExternalForceRegularizedBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void ShanChenExternalForceRegularizedBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);

    T invOmega = 1. / this->getOmega();
    Array<T, Descriptor<T>::d> force;
    force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));
    Array<T, Descriptor<T>::d> jCorrected(j + invOmega * Descriptor<T>::fullRho(rhoBar) * force);

    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T, Descriptor>::rlb_collision(
        cell, rhoBar, invRho, jCorrected, PiNeq, this->getOmega());

    if (cell.takesStatistics()) {
        T uSqr = T();
        T invRho = Descriptor<T>::invRho(rhoBar);
        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            uSqr += util::sqr(j[iD] * invRho + 0.5 * force[iD]);
        }
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void ShanChenExternalForceRegularizedBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);

    T invOmega = 1. / this->getOmega();
    Array<T, Descriptor<T>::d> force;
    force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));
    Array<T, Descriptor<T>::d> jCorrected(j + invOmega * Descriptor<T>::fullRho(rhoBar) * force);

    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T, Descriptor>::rlb_collision(
        cell, rhoBar, invRho, jCorrected, PiNeq, this->getOmega());

    if (cell.takesStatistics()) {
        T uSqr = T();
        T invRho = Descriptor<T>::invRho(rhoBar);
        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            uSqr += util::sqr(j[iD] * invRho + 0.5 * force[iD]);
        }
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T ShanChenExternalForceRegularizedBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

}  // namespace plb

#endif  // EXTERNAL_FORCE_DYNAMICS_HH
