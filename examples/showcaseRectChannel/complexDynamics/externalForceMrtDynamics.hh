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

/* Orestis Malaspinas contributed this code. */

/** \file
 * MRT dynamics -- generic implementation.
 */
#ifndef EXTERNAL_FORCE_MRT_DYNAMICS_HH
#define EXTERNAL_FORCE_MRT_DYNAMICS_HH

#include <algorithm>
#include <limits>

#include "core/latticeStatistics.h"
#include "externalForceMrtDynamics.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/mrtTemplates.h"

namespace plb {

/* *************** Class GuoExternalForceMRTdynamics ***********************************************
 */
template <typename T, template <typename U> class Descriptor>
int GuoExternalForceMRTdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, GuoExternalForceMRTdynamics<T, Descriptor> >(
        "MRT_ExternalForce_Guo");

template <typename T, template <typename U> class Descriptor>
GuoExternalForceMRTdynamics<T, Descriptor>::GuoExternalForceMRTdynamics(T omega_) :
    ExternalForceDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
GuoExternalForceMRTdynamics<T, Descriptor>::GuoExternalForceMRTdynamics(
    HierarchicUnserializer &unserializer) :
    ExternalForceDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
GuoExternalForceMRTdynamics<T, Descriptor> *GuoExternalForceMRTdynamics<T, Descriptor>::clone()
    const
{
    return new GuoExternalForceMRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int GuoExternalForceMRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceMRTdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    typedef mrtTemplates<T, Descriptor> mrtTemp;

    T rhoBar = momentTemplates<T, Descriptor>::get_rhoBar(cell);

    Array<T, Descriptor<T>::d> u;
    this->computeVelocity(cell, u);

    T jSqr = mrtTemp::mrtCollisionWithForce(cell, rhoBar, u, this->getOmega(), (T)1);

    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(statistics, rhoBar, jSqr * invRho * invRho);
    }
}

template <typename T, template <typename U> class Descriptor>
T GuoExternalForceMRTdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class GuoExternalForceSmagorinskyMRTdynamics
 * *********************************************** */
template <typename T, template <typename U> class Descriptor>
int GuoExternalForceSmagorinskyMRTdynamics<T, Descriptor>::id = meta::registerGeneralDynamics<
    T, Descriptor, GuoExternalForceSmagorinskyMRTdynamics<T, Descriptor> >(
    "MRT_ExternalForce_Guo_Smagorinsky");

template <typename T, template <typename U> class Descriptor>
GuoExternalForceSmagorinskyMRTdynamics<T, Descriptor>::GuoExternalForceSmagorinskyMRTdynamics(
    T omega_, T cSmago_) :
    ExternalForceDynamics<T, Descriptor>(omega_), cSmago(cSmago_)
{ }

template <typename T, template <typename U> class Descriptor>
GuoExternalForceSmagorinskyMRTdynamics<T, Descriptor>::GuoExternalForceSmagorinskyMRTdynamics(
    HierarchicUnserializer &unserializer) :
    ExternalForceDynamics<T, Descriptor>(T()), cSmago(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
GuoExternalForceSmagorinskyMRTdynamics<T, Descriptor>
    *GuoExternalForceSmagorinskyMRTdynamics<T, Descriptor>::clone() const
{
    return new GuoExternalForceSmagorinskyMRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int GuoExternalForceSmagorinskyMRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceSmagorinskyMRTdynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    ExternalForceDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceSmagorinskyMRTdynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    ExternalForceDynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceSmagorinskyMRTdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    typedef mrtTemplates<T, Descriptor> mrtTemp;

    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T rho = Descriptor<T>::fullRho(rhoBar);

    T tau = (T)1 / this->getOmega();

    T normPiNeq = std::sqrt((T)2 * SymmetricTensor<T, Descriptor>::tensorNormSqr(PiNeq));
    normPiNeq *= (T)2 * rho * util::sqr(Descriptor<T>::cs2 * cSmago * cSmago);

    T tauSGS = T();
    if (normPiNeq != T()) {  // test to avoid division by 0
        Array<T, SymmetricTensor<T, Descriptor>::n> S =
            (-rho * tau * Descriptor<T>::cs2
             + std::sqrt(util::sqr(rho * tau * Descriptor<T>::cs2) + normPiNeq))
            / normPiNeq * PiNeq;

        T sNorm = std::sqrt((T)2 * SymmetricTensor<T, Descriptor>::tensorNormSqr(S));

        tauSGS = cSmago * cSmago * sNorm * Descriptor<T>::invCs2;
    }
    T omega = (T)1 / (tau + tauSGS);

    Array<T, Descriptor<T>::d> u;
    this->computeVelocity(cell, u);

    T jSqr = mrtTemp::mrtCollisionWithForce(cell, rhoBar, u, omega, (T)1);

    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(statistics, rhoBar, jSqr * invRho * invRho);
    }
}

template <typename T, template <typename U> class Descriptor>
T GuoExternalForceSmagorinskyMRTdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class GuoExternalForceSmagorinskyIncMRTdynamics
 * *********************************************** */
template <typename T, template <typename U> class Descriptor>
int GuoExternalForceSmagorinskyIncMRTdynamics<T, Descriptor>::id = meta::registerGeneralDynamics<
    T, Descriptor, GuoExternalForceSmagorinskyIncMRTdynamics<T, Descriptor> >(
    "IncMRT_ExternalForce_Guo_Smagorinsky");

template <typename T, template <typename U> class Descriptor>
GuoExternalForceSmagorinskyIncMRTdynamics<T, Descriptor>::GuoExternalForceSmagorinskyIncMRTdynamics(
    T omega_, T cSmago_) :
    IncExternalForceDynamics<T, Descriptor>(omega_), cSmago(cSmago_)
{ }

template <typename T, template <typename U> class Descriptor>
GuoExternalForceSmagorinskyIncMRTdynamics<T, Descriptor>::GuoExternalForceSmagorinskyIncMRTdynamics(
    HierarchicUnserializer &unserializer) :
    IncExternalForceDynamics<T, Descriptor>(T()), cSmago(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
GuoExternalForceSmagorinskyIncMRTdynamics<T, Descriptor>
    *GuoExternalForceSmagorinskyIncMRTdynamics<T, Descriptor>::clone() const
{
    return new GuoExternalForceSmagorinskyIncMRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int GuoExternalForceSmagorinskyIncMRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceSmagorinskyIncMRTdynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    IncExternalForceDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceSmagorinskyIncMRTdynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    IncExternalForceDynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
bool GuoExternalForceSmagorinskyIncMRTdynamics<T, Descriptor>::velIsJ() const
{
    return true;
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceSmagorinskyIncMRTdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    typedef mrtTemplates<T, Descriptor> mrtTemp;

    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T rho = Descriptor<T>::fullRho(rhoBar);

    T tau = (T)1 / this->getOmega();

    T normPiNeq = std::sqrt((T)2 * SymmetricTensor<T, Descriptor>::tensorNormSqr(PiNeq));
    normPiNeq *= (T)2 * rho * util::sqr(Descriptor<T>::cs2 * cSmago * cSmago);

    T tauSGS = T();
    if (normPiNeq != T()) {  // test to avoid division by 0
        Array<T, SymmetricTensor<T, Descriptor>::n> S =
            (-rho * tau * Descriptor<T>::cs2
             + std::sqrt(util::sqr(rho * tau * Descriptor<T>::cs2) + normPiNeq))
            / normPiNeq * PiNeq;

        T sNorm = std::sqrt((T)2 * SymmetricTensor<T, Descriptor>::tensorNormSqr(S));

        tauSGS = cSmago * cSmago * sNorm * Descriptor<T>::invCs2;
    }
    T omega = (T)1 / (tau + tauSGS);

    Array<T, Descriptor<T>::d> u;
    this->computeVelocity(cell, u);

    T jSqr = mrtTemp::incMrtCollisionWithForce(cell, rhoBar, u, omega, (T)1);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T GuoExternalForceSmagorinskyIncMRTdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, (T)1, j, jSqr);
}

/* *************** Class GuoExternalForceConsistentSmagorinskyMRTdynamics
 * *********************************************** */
template <typename T, template <typename U> class Descriptor>
int GuoExternalForceConsistentSmagorinskyMRTdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<
        T, Descriptor, GuoExternalForceConsistentSmagorinskyMRTdynamics<T, Descriptor> >(
        "MRT_ExternalForce_Guo_ConsistentSmagorinsky");

template <typename T, template <typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyMRTdynamics<
    T, Descriptor>::GuoExternalForceConsistentSmagorinskyMRTdynamics(T omega_, T cSmago_) :
    ExternalForceDynamics<T, Descriptor>(omega_), cSmago(cSmago_)
{ }

template <typename T, template <typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyMRTdynamics<T, Descriptor>::
    GuoExternalForceConsistentSmagorinskyMRTdynamics(HierarchicUnserializer &unserializer) :
    ExternalForceDynamics<T, Descriptor>(T()), cSmago(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyMRTdynamics<T, Descriptor>
    *GuoExternalForceConsistentSmagorinskyMRTdynamics<T, Descriptor>::clone() const
{
    return new GuoExternalForceConsistentSmagorinskyMRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int GuoExternalForceConsistentSmagorinskyMRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceConsistentSmagorinskyMRTdynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    ExternalForceDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceConsistentSmagorinskyMRTdynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    ExternalForceDynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceConsistentSmagorinskyMRTdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    typedef mrtTemplates<T, Descriptor> mrtTemp;

    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T rho = Descriptor<T>::fullRho(rhoBar);

    T tau = (T)1 / this->getOmega();

    T normPiNeq = std::sqrt((T)2 * SymmetricTensor<T, Descriptor>::tensorNormSqr(PiNeq));
    normPiNeq *= (T)2 * rho * util::sqr(Descriptor<T>::cs2 * cSmago * cSmago);

    Array<T, SymmetricTensor<T, Descriptor>::n> S;
    S.resetToZero();
    if (normPiNeq != T()) {  // test to avoid division by 0
        Array<T, SymmetricTensor<T, Descriptor>::n> S =
            (-rho * tau * Descriptor<T>::cs2
             + std::sqrt(util::sqr(rho * tau * Descriptor<T>::cs2) + normPiNeq))
            / normPiNeq * PiNeq;
    }

    Array<T, Descriptor<T>::d> u;
    this->computeVelocity(cell, u);

    T jSqr = mrtTemp::smagorinskyMrtCollisionWithForce(
        cell, rhoBar, u, S, cSmago, this->getOmega(), (T)1);

    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(statistics, rhoBar, jSqr * invRho * invRho);
    }
}

template <typename T, template <typename U> class Descriptor>
T GuoExternalForceConsistentSmagorinskyMRTdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class GuoExternalForceConsistentSmagorinskyIncMRTdynamics
 * *********************************************** */
template <typename T, template <typename U> class Descriptor>
int GuoExternalForceConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<
        T, Descriptor, GuoExternalForceConsistentSmagorinskyIncMRTdynamics<T, Descriptor> >(
        "IncMRT_ExternalForce_Guo_ConsistentSmagorinsky");

template <typename T, template <typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyIncMRTdynamics<
    T, Descriptor>::GuoExternalForceConsistentSmagorinskyIncMRTdynamics(T omega_, T cSmago_) :
    IncExternalForceDynamics<T, Descriptor>(omega_), cSmago(cSmago_)
{ }

template <typename T, template <typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::
    GuoExternalForceConsistentSmagorinskyIncMRTdynamics(HierarchicUnserializer &unserializer) :
    IncExternalForceDynamics<T, Descriptor>(T()), cSmago(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
GuoExternalForceConsistentSmagorinskyIncMRTdynamics<T, Descriptor>
    *GuoExternalForceConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::clone() const
{
    return new GuoExternalForceConsistentSmagorinskyIncMRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int GuoExternalForceConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    IncExternalForceDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    IncExternalForceDynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
bool GuoExternalForceConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::velIsJ() const
{
    return true;
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    typedef mrtTemplates<T, Descriptor> mrtTemp;

    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T rho = Descriptor<T>::fullRho(rhoBar);

    T tau = (T)1 / this->getOmega();

    T normPiNeq = std::sqrt((T)2 * SymmetricTensor<T, Descriptor>::tensorNormSqr(PiNeq));
    normPiNeq *= (T)2 * rho * util::sqr(Descriptor<T>::cs2 * cSmago * cSmago);

    Array<T, SymmetricTensor<T, Descriptor>::n> S;
    S.resetToZero();
    if (normPiNeq != T()) {  // test to avoid division by 0
        S = (-rho * tau * Descriptor<T>::cs2
             + std::sqrt(util::sqr(rho * tau * Descriptor<T>::cs2) + normPiNeq))
            / normPiNeq * PiNeq;
    }

    Array<T, Descriptor<T>::d> u;
    this->computeVelocity(cell, u);

    T jSqr = mrtTemp::incSmagorinskyMrtCollisionWithForce(
        cell, rhoBar, u, S, cSmago, this->getOmega(), (T)1);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T GuoExternalForceConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, (T)1, j, jSqr);
}

/* *************** Class GuoExternalForceAndMomentMRTdynamics
 * *********************************************** */

template <typename T, template <typename U> class Descriptor>
int GuoExternalForceAndMomentMRTdynamics<T, Descriptor>::id = meta::registerGeneralDynamics<
    T, Descriptor, GuoExternalForceAndMomentMRTdynamics<T, Descriptor> >(
    "MRT_ExternalForceAndMoment_Guo");

template <typename T, template <typename U> class Descriptor>
GuoExternalForceAndMomentMRTdynamics<T, Descriptor>::GuoExternalForceAndMomentMRTdynamics(
    T omega_) :
    ExternalForceDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
GuoExternalForceAndMomentMRTdynamics<T, Descriptor>::GuoExternalForceAndMomentMRTdynamics(
    HierarchicUnserializer &unserializer) :
    ExternalForceDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
GuoExternalForceAndMomentMRTdynamics<T, Descriptor>
    *GuoExternalForceAndMomentMRTdynamics<T, Descriptor>::clone() const
{
    return new GuoExternalForceAndMomentMRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int GuoExternalForceAndMomentMRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceAndMomentMRTdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    typedef mrtTemplates<T, Descriptor> mrtTemp;

    T rho = *cell.getExternal(Descriptor<T>::ExternalField::densityBeginsAt);
    T rhoBar = Descriptor<T>::rhoBar(rho);
    T invRho = Descriptor<T>::invRho(rhoBar);
    Array<T, Descriptor<T>::d> u;
    u.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::momentumBeginsAt));
    u = u * invRho;

    T jSqr = mrtTemp::mrtCollisionWithForce(cell, rhoBar, u, this->getOmega(), (T)1);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr * invRho * invRho);
    }
}

template <typename T, template <typename U> class Descriptor>
T GuoExternalForceAndMomentMRTdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class HeExternalForceMRTdynamics ***********************************************
 */
template <typename T, template <typename U> class Descriptor>
int HeExternalForceMRTdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, HeExternalForceMRTdynamics<T, Descriptor> >(
        "MRT_ExternalForce_He");

template <typename T, template <typename U> class Descriptor>
HeExternalForceMRTdynamics<T, Descriptor>::HeExternalForceMRTdynamics(T omega_) :
    ExternalForceDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
HeExternalForceMRTdynamics<T, Descriptor>::HeExternalForceMRTdynamics(
    HierarchicUnserializer &unserializer) :
    ExternalForceDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
HeExternalForceMRTdynamics<T, Descriptor> *HeExternalForceMRTdynamics<T, Descriptor>::clone() const
{
    return new HeExternalForceMRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int HeExternalForceMRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void HeExternalForceMRTdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    typedef mrtTemplates<T, Descriptor> mrtTemp;

    // get rho and j from external mometums (calculated during dataProcessors step in ShanChen
    // processor)

    T rho = *cell.getExternal(Descriptor<T>::ExternalField::densityBeginsAt);
    T rhoBar = Descriptor<T>::rhoBar(rho);
    T invRho = Descriptor<T>::invRho(rhoBar);
    Array<T, Descriptor<T>::d> j;
    j.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::momentumBeginsAt));

    T jSqr = mrtTemp::mrtCollisionWithHeForce(cell, rhoBar, j * invRho, this->getOmega(), (T)1);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr * invRho * invRho);
    }
}

template <typename T, template <typename U> class Descriptor>
T HeExternalForceMRTdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

}  // namespace plb
#endif
