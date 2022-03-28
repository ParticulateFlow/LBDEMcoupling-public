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

/* Orestis Malaspinas designed some of the classes and concepts contained
 * in this file. */

#ifndef SMAGORINSKY_DYNAMICS_HH
#define SMAGORINSKY_DYNAMICS_HH

#include <cmath>

#include "complexDynamics/smagorinskyDynamics.h"
#include "core/dynamicsIdentifiers.h"
#include "core/latticeStatistics.h"
#include "core/util.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/externalForceTemplates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/mrtTemplates.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
struct SmagoOperations {
    static T computePrefactor(T omega0, T cSmago)
    {
        return (T)0.5 * util::sqr(cSmago * omega0 * Descriptor<T>::invCs2);
    }
    static T computeOmega(
        T omega0, T preFactor, T rhoBar, Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq)
    {
        T PiNeqNormSqr = SymmetricTensor<T, Descriptor>::tensorNormSqr(PiNeq);
        T PiNeqNorm = std::sqrt(PiNeqNormSqr);
        T alpha = preFactor * Descriptor<T>::invRho(rhoBar);
        T linearTerm = alpha * PiNeqNorm;
        T squareTerm = (T)2 * alpha * alpha * PiNeqNormSqr;
        // In the following formula, the square-root appearing in the explicit form of
        //   omega is developed to second-order.
        return omega0 * (1 - linearTerm + squareTerm);
    }

    // In the consistent formulation for SGS terms
    // the turbulence modelling is done by adding a term
    // of the form t_i/(2*tau_0*cs^4)*H_i^2:T to the collision, where
    // T = \overline{rho u u}-\bar{rho} \bar{u} \bar{u}.
    static Array<T, Descriptor<T>::q> computeSgsTensorTerm(
        T rho, const Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq, T cSmago, T omega0)
    {
        T tau = (T)1 / omega0;

        //     Array<T,SymmetricTensor<T,Descriptor>::n> S = -PiNeq;

        T normPiNeq = std::sqrt((T)2 * SymmetricTensor<T, Descriptor>::tensorNormSqr(PiNeq));
        normPiNeq *= (T)2 * rho * util::sqr(Descriptor<T>::cs2 * cSmago * cSmago);
        Array<T, SymmetricTensor<T, Descriptor>::n> S;
        S.resetToZero();
        if (normPiNeq != T()) {  // test to avoid division per 0
            S = -(-rho * tau * Descriptor<T>::cs2
                  + std::sqrt(util::sqr(rho * tau * Descriptor<T>::cs2) + normPiNeq))
                / normPiNeq * PiNeq;
        }
        T sNorm = std::sqrt((T)2 * SymmetricTensor<T, Descriptor>::tensorNormSqr(S));
        T preFactor = (T)0.5 * Descriptor<T>::invCs2 * Descriptor<T>::invCs2 * omega0 * cSmago
                      * cSmago * sNorm;

        Array<T, Descriptor<T>::q> sgsTerm;
        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            Array<T, SymmetricTensor<T, Descriptor>::n> H2 =
                HermiteTemplate<T, Descriptor>::order2(iPop);
            T H_S = SymmetricTensor<T, Descriptor>::contractIndexes(H2, S);
            sgsTerm[iPop] = -Descriptor<T>::t[iPop] * preFactor * H_S;
        }
        return sgsTerm;
    }

    // In the consistent formulation for SGS terms
    // the turbulence modelling is done by adding a term
    // of the form t_i/(2*tau_0*cs^4)*H_i^2:T to the collision, where
    // T = \overline{rho u u}-\bar{rho} \bar{u} \bar{u}.
    static Array<T, SymmetricTensor<T, Descriptor>::n> computeSmagorinskyTensorTerm(
        T rho, const Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq, T cSmago, T omega0)
    {
        T tau = (T)1 / omega0;

        T normPiNeq = std::sqrt((T)2 * SymmetricTensor<T, Descriptor>::tensorNormSqr(PiNeq));
        normPiNeq *= (T)2 * rho * util::sqr(Descriptor<T>::cs2 * cSmago * cSmago);
        Array<T, SymmetricTensor<T, Descriptor>::n> S;
        S.resetToZero();
        if (normPiNeq != T()) {  // test to avoid division per 0
            S = -(-rho * tau * Descriptor<T>::cs2
                  + std::sqrt(util::sqr(rho * tau * Descriptor<T>::cs2) + normPiNeq))
                / normPiNeq * PiNeq;
        }
        T sNorm = std::sqrt((T)2 * SymmetricTensor<T, Descriptor>::tensorNormSqr(S));
        return omega0 * cSmago * cSmago * sNorm * S;
    }

    // In the consistent formulation for SGS terms
    // the turbulence modelling is done by adding a term
    // of the form t_i/(2*tau_0*cs^4)*H_i^2:T to the collision, where
    // T = \overline{rho u u}-\bar{rho} \bar{u} \bar{u}.
    static Array<T, Descriptor<T>::q> computeSgsTensorTerm(
        const Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq, T omega0)
    {
        T preFactor = (T)0.5 * Descriptor<T>::invCs2 * Descriptor<T>::invCs2 * omega0;

        Array<T, Descriptor<T>::q> sgsTerm;
        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            Array<T, SymmetricTensor<T, Descriptor>::n> H2 =
                HermiteTemplate<T, Descriptor>::order2(iPop);
            T H_PiNeq = SymmetricTensor<T, Descriptor>::contractIndexes(H2, PiNeq);
            sgsTerm[iPop] = Descriptor<T>::t[iPop] * preFactor * H_PiNeq;
        }
        return sgsTerm;
    }
};

template <typename T, template <typename U> class Descriptor>
int SmagorinskyDynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, SmagorinskyDynamics<T, Descriptor> >(
        "Smagorinsky_Generic");

template <typename T, template <typename U> class Descriptor>
SmagorinskyDynamics<T, Descriptor>::SmagorinskyDynamics(
    Dynamics<T, Descriptor> *baseDynamics_, T omega0_, T cSmago_, bool automaticPrepareCollision) :
    OmegaFromPiDynamics<T, Descriptor>(baseDynamics_, automaticPrepareCollision),
    omega0(omega0_),
    cSmago(cSmago_),
    preFactor(SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago))
{ }

template <typename T, template <typename U> class Descriptor>
SmagorinskyDynamics<T, Descriptor>::SmagorinskyDynamics(HierarchicUnserializer &unserializer) :
    OmegaFromPiDynamics<T, Descriptor>(0, false), omega0(T()), cSmago(T()), preFactor(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
T SmagorinskyDynamics<T, Descriptor>::getOmegaFromPiAndRhoBar(
    Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T rhoBar) const
{
    return SmagoOperations<T, Descriptor>::computeOmega(omega0, preFactor, rhoBar, PiNeq);
}

template <typename T, template <typename U> class Descriptor>
SmagorinskyDynamics<T, Descriptor> *SmagorinskyDynamics<T, Descriptor>::clone() const
{
    return new SmagorinskyDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int SmagorinskyDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyDynamics<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    OmegaFromPiDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(omega0);
    serializer.addValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyDynamics<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    OmegaFromPiDynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(omega0);
    unserializer.readValue(cSmago);
    preFactor = SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago);
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyDynamics<T, Descriptor>::setOmega(T omega0_)
{
    // Just to be sure to avoid an undefined state, also reset the value of
    // omega in the baseDynamics, although the actual value of omega in the
    // baseDynamics is omega=omega0+deltaOmega and is recomputed before each
    // collision.
    omega0 = omega0_;
    preFactor = SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago);
    this->getBaseDynamics().setOmega(omega0);
}

template <typename T, template <typename U> class Descriptor>
T SmagorinskyDynamics<T, Descriptor>::getOmega() const
{
    return omega0;
}

template <typename T, template <typename U> class Descriptor>
T SmagorinskyDynamics<T, Descriptor>::getDynamicParameter(
    plint whichParameter, Cell<T, Descriptor> const &cell) const
{
    if (whichParameter == dynamicParams::smagorinskyConstant) {
        return cSmago;
    } else if (whichParameter == dynamicParams::dynamicOmega) {
        T rhoBar;
        Array<T, Descriptor<T>::d> j;
        Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
        cell.getDynamics().computeRhoBarJPiNeq(cell, rhoBar, j, PiNeq);
        T preFactor = SmagoOperations<T, Descriptor>::computePrefactor(this->getOmega(), cSmago);
        T omega = SmagoOperations<T, Descriptor>::computeOmega(
            this->getOmega(), preFactor, rhoBar, PiNeq);
        return omega;
    } else {
        return OmegaFromPiDynamics<T, Descriptor>::getParameter(whichParameter);
    }
}

/* *************** Class SmagorinskyBGKdynamics ************************************** */

template <typename T, template <typename U> class Descriptor>
int SmagorinskyBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, SmagorinskyBGKdynamics<T, Descriptor> >(
        "BGK_Smagorinsky");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
SmagorinskyBGKdynamics<T, Descriptor>::SmagorinskyBGKdynamics(T omega0_, T cSmago_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega0_),
    omega0(omega0_),
    cSmago(cSmago_),
    preFactor(SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago))
{ }

template <typename T, template <typename U> class Descriptor>
SmagorinskyBGKdynamics<T, Descriptor>::SmagorinskyBGKdynamics(
    HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T()), omega0(T()), cSmago(T()), preFactor(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyBGKdynamics<T, Descriptor>::setOmega(T omega0_)
{
    omega0 = omega0_;
    preFactor = SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago);
}

template <typename T, template <typename U> class Descriptor>
T SmagorinskyBGKdynamics<T, Descriptor>::getOmega() const
{
    return omega0;
}

template <typename T, template <typename U> class Descriptor>
T SmagorinskyBGKdynamics<T, Descriptor>::getDynamicParameter(
    plint whichParameter, Cell<T, Descriptor> const &cell) const
{
    if (whichParameter == dynamicParams::dynamicOmega) {
        T rhoBar;
        Array<T, Descriptor<T>::d> j;
        Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
        momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
        T omega = SmagoOperations<T, Descriptor>::computeOmega(omega0, preFactor, rhoBar, PiNeq);
        return omega;
    } else if (whichParameter == dynamicParams::smagorinskyConstant) {
        return cSmago;
    } else {
        return IsoThermalBulkDynamics<T, Descriptor>::getParameter(whichParameter);
    }
}

template <typename T, template <typename U> class Descriptor>
SmagorinskyBGKdynamics<T, Descriptor> *SmagorinskyBGKdynamics<T, Descriptor>::clone() const
{
    return new SmagorinskyBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int SmagorinskyBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyBGKdynamics<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    IsoThermalBulkDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(omega0);
    serializer.addValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyBGKdynamics<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    IsoThermalBulkDynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(omega0);
    unserializer.readValue(cSmago);
    preFactor = SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago);
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T omega = SmagoOperations<T, Descriptor>::computeOmega(omega0, preFactor, rhoBar, PiNeq);
    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, j, omega);
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);
    T omega = SmagoOperations<T, Descriptor>::computeOmega(omega0, preFactor, rhoBar, PiNeq);
    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, j, omega);
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T SmagorinskyBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class ConsistentSmagorinskyBGKdynamics ************************************** */

template <typename T, template <typename U> class Descriptor>
int ConsistentSmagorinskyBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, ConsistentSmagorinskyBGKdynamics<T, Descriptor> >(
        "BGK_Consistent_Smagorinsky");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyBGKdynamics<T, Descriptor>::ConsistentSmagorinskyBGKdynamics(
    T omega0_, T cSmago_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega0_), omega0(omega0_), cSmago(cSmago_)
{ }

template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyBGKdynamics<T, Descriptor>::ConsistentSmagorinskyBGKdynamics(
    HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T()), omega0(T()), cSmago(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
T ConsistentSmagorinskyBGKdynamics<T, Descriptor>::getOmega() const
{
    return omega0;
}

template <typename T, template <typename U> class Descriptor>
T ConsistentSmagorinskyBGKdynamics<T, Descriptor>::getDynamicParameter(
    plint whichParameter, Cell<T, Descriptor> const &cell) const
{
    if (whichParameter == dynamicParams::smagorinskyConstant) {
        return cSmago;
    } else if (whichParameter == dynamicParams::dynamicOmega) {
        T rhoBar;
        Array<T, Descriptor<T>::d> j;
        Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
        momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
        T preFactor = SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago);
        T omega = SmagoOperations<T, Descriptor>::computeOmega(omega0, preFactor, rhoBar, PiNeq);
        return omega;
    } else {
        return IsoThermalBulkDynamics<T, Descriptor>::getParameter(whichParameter);
    }
}

template <typename T, template <typename U> class Descriptor>
T ConsistentSmagorinskyBGKdynamics<T, Descriptor>::getParameter(plint whichParameter) const
{
    switch (whichParameter) {
    case dynamicParams::omega_shear:
        return omega0;
    case dynamicParams::smagorinskyConstant:
        return cSmago;
    default:
        return IsoThermalBulkDynamics<T, Descriptor>::getParameter(whichParameter);
    }
}

template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyBGKdynamics<T, Descriptor>
    *ConsistentSmagorinskyBGKdynamics<T, Descriptor>::clone() const
{
    return new ConsistentSmagorinskyBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int ConsistentSmagorinskyBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyBGKdynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    IsoThermalBulkDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(omega0);
    serializer.addValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyBGKdynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    IsoThermalBulkDynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(omega0);
    unserializer.readValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyBGKdynamics<T, Descriptor>::setOmega(T omega0_)
{
    omega0 = omega0_;
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T rho = Descriptor<T>::fullRho(rhoBar);

    Array<T, Descriptor<T>::q> sgsTerm =
        SmagoOperations<T, Descriptor>::computeSgsTensorTerm(rho, PiNeq, cSmago, this->getOmega());

    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, j, omega0);
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] += sgsTerm[iPop];
    }

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);
    T rho = Descriptor<T>::fullRho(rhoBar);

    T tau = (T)1 / this->getOmega();

    T normPiNeq = std::sqrt((T)2 * SymmetricTensor<T, Descriptor>::tensorNormSqr(PiNeq));
    normPiNeq *= (T)2 * rho * util::sqr(Descriptor<T>::cs2 * cSmago * cSmago);

    Array<T, SymmetricTensor<T, Descriptor>::n> S;
    S.resetToZero();
    if (normPiNeq != T()) {  // test to avoid division per 0
        S = (-rho * tau * Descriptor<T>::cs2
             + std::sqrt(util::sqr(rho * tau * Descriptor<T>::cs2) + normPiNeq))
            / normPiNeq * PiNeq;
    }
    T sNorm = std::sqrt((T)2 * SymmetricTensor<T, Descriptor>::tensorNormSqr(S));
    T preFactor =
        (T)0.5 * Descriptor<T>::invCs2 * Descriptor<T>::invCs2 * omega0 * cSmago * cSmago * sNorm;

    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, j, omega0);
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        Array<T, SymmetricTensor<T, Descriptor>::n> H2 =
            HermiteTemplate<T, Descriptor>::order2(iPop);
        T H_S = SymmetricTensor<T, Descriptor>::contractIndexes(H2, S);
        cell[iPop] += Descriptor<T>::t[iPop] * preFactor * H_S;
    }

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T ConsistentSmagorinskyBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class ConsistentSmagorinskyCompleteTRTdynamics
 * ************************************** */

template <typename T, template <typename U> class Descriptor>
int ConsistentSmagorinskyCompleteTRTdynamics<T, Descriptor>::id = meta::registerGeneralDynamics<
    T, Descriptor, ConsistentSmagorinskyCompleteTRTdynamics<T, Descriptor> >(
    "CompleteTRT_Consistent_Smagorinsky");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyCompleteTRTdynamics<T, Descriptor>::ConsistentSmagorinskyCompleteTRTdynamics(
    T omega_, T psi_, T cSmago_) :
    CompleteTRTdynamics<T, Descriptor>(omega_, psi_), cSmago(cSmago_)
{ }

template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyCompleteTRTdynamics<T, Descriptor>::ConsistentSmagorinskyCompleteTRTdynamics(
    T omega_, T cSmago_) :
    CompleteTRTdynamics<T, Descriptor>(omega_), cSmago(cSmago_)
{ }

template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyCompleteTRTdynamics<T, Descriptor>::ConsistentSmagorinskyCompleteTRTdynamics(
    HierarchicUnserializer &unserializer) :
    CompleteTRTdynamics<T, Descriptor>(T(), T()), cSmago(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
T ConsistentSmagorinskyCompleteTRTdynamics<T, Descriptor>::getDynamicParameter(
    plint whichParameter, Cell<T, Descriptor> const &cell) const
{
    if (whichParameter == dynamicParams::smagorinskyConstant) {
        return cSmago;
    } else if (whichParameter == dynamicParams::dynamicOmega) {
        T rhoBar;
        Array<T, Descriptor<T>::d> j;
        Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
        momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
        T preFactor = SmagoOperations<T, Descriptor>::computePrefactor(this->getOmega(), cSmago);
        T omega = SmagoOperations<T, Descriptor>::computeOmega(
            this->getOmega(), preFactor, rhoBar, PiNeq);
        return omega;
    } else {
        return CompleteTRTdynamics<T, Descriptor>::getParameter(whichParameter);
    }
}

template <typename T, template <typename U> class Descriptor>
int ConsistentSmagorinskyCompleteTRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyCompleteTRTdynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    CompleteTRTdynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyCompleteTRTdynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    CompleteTRTdynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyCompleteTRTdynamics<T, Descriptor>
    *ConsistentSmagorinskyCompleteTRTdynamics<T, Descriptor>::clone() const
{
    return new ConsistentSmagorinskyCompleteTRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyCompleteTRTdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T uSqr = dynamicsTemplates<T, Descriptor>::complete_mrt_smagorinsky_ma2_collision(
        cell, cSmago, this->getOmega(), this->getPsi());

    if (cell.takesStatistics()) {
        T rhoBar = momentTemplates<T, Descriptor>::get_rhoBar(cell);
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyCompleteTRTdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T uSqr = dynamicsTemplates<T, Descriptor>::complete_mrt_smagorinsky_ma2_ext_rhoBar_j_collision(
        cell, rhoBar, j, 6, cSmago, this->getOmega(), this->getPsi());

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

/* *************** Class ConsistentSmagorinskyCompleteRegularizedBGKdynamics
 * ************************************** */

template <typename T, template <typename U> class Descriptor>
int ConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<
        T, Descriptor, ConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor> >(
        "RegularizedCompleteBGK_Consistent_Smagorinsky");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyCompleteRegularizedBGKdynamics<
    T, Descriptor>::ConsistentSmagorinskyCompleteRegularizedBGKdynamics(T omega_, T cSmago_) :
    CompleteRegularizedBGKdynamics<T, Descriptor>(omega_), cSmago(cSmago_)
{ }

template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor>::
    ConsistentSmagorinskyCompleteRegularizedBGKdynamics(HierarchicUnserializer &unserializer) :
    CompleteRegularizedBGKdynamics<T, Descriptor>(T()), cSmago(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
T ConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor>::getDynamicParameter(
    plint whichParameter, Cell<T, Descriptor> const &cell) const
{
    if (whichParameter == dynamicParams::smagorinskyConstant) {
        return cSmago;
    } else if (whichParameter == dynamicParams::dynamicOmega) {
        T rhoBar;
        Array<T, Descriptor<T>::d> j;
        Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
        momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
        T preFactor = SmagoOperations<T, Descriptor>::computePrefactor(this->getOmega(), cSmago);
        T omega = SmagoOperations<T, Descriptor>::computeOmega(
            this->getOmega(), preFactor, rhoBar, PiNeq);
        return omega;
    } else {
        return CompleteRegularizedBGKdynamics<T, Descriptor>::getParameter(whichParameter);
    }
}

template <typename T, template <typename U> class Descriptor>
int ConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    CompleteRegularizedBGKdynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    CompleteRegularizedBGKdynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor>
    *ConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor>::clone() const
{
    return new ConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> piNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, piNeq);
    T uSqr = dynamicsTemplates<T, Descriptor>::complete_regularized_bgk_ma2_collision(
        cell, rhoBar, j, piNeq, this->getOmega());

    T rho = Descriptor<T>::fullRho(rhoBar);

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

    if (cell.takesStatistics()) {
        T rhoBar = momentTemplates<T, Descriptor>::get_rhoBar(cell);
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    Array<T, SymmetricTensor<T, Descriptor>::n> piNeq;
    momentTemplates<T, Descriptor>::compute_PiNeq(cell, rhoBar, j, piNeq, invRho);
    T uSqr = dynamicsTemplates<T, Descriptor>::complete_regularized_bgk_ma2_collision(
        cell, rhoBar, j, piNeq, this->getOmega());

    T rho = Descriptor<T>::fullRho(rhoBar);

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

    if (cell.takesStatistics()) {
        T rhoBar = momentTemplates<T, Descriptor>::get_rhoBar(cell);
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

/* *************** Class ConsistentSmagorinskyCompleteRegularizedTRTdynamics
 * ************************************** */

template <typename T, template <typename U> class Descriptor>
int ConsistentSmagorinskyCompleteRegularizedTRTdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<
        T, Descriptor, ConsistentSmagorinskyCompleteRegularizedTRTdynamics<T, Descriptor> >(
        "RegularizedCompleteTRT_Consistent_Smagorinsky");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyCompleteRegularizedTRTdynamics<T, Descriptor>::
    ConsistentSmagorinskyCompleteRegularizedTRTdynamics(T omega_, T psi_, int order_, T cSmago_) :
    CompleteRegularizedTRTdynamics<T, Descriptor>(omega_, psi_, order_),
    order(order_),
    cSmago(cSmago_)
{ }

template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyCompleteRegularizedTRTdynamics<T, Descriptor>::
    ConsistentSmagorinskyCompleteRegularizedTRTdynamics(T omega_, int order_, T cSmago_) :
    CompleteRegularizedTRTdynamics<T, Descriptor>(omega_, order_), order(order_), cSmago(cSmago_)
{ }

template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyCompleteRegularizedTRTdynamics<T, Descriptor>::
    ConsistentSmagorinskyCompleteRegularizedTRTdynamics(HierarchicUnserializer &unserializer) :
    CompleteRegularizedTRTdynamics<T, Descriptor>(T(), T(), 0), order(0), cSmago(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
T ConsistentSmagorinskyCompleteRegularizedTRTdynamics<T, Descriptor>::getDynamicParameter(
    plint whichParameter, Cell<T, Descriptor> const &cell) const
{
    if (whichParameter == dynamicParams::smagorinskyConstant) {
        return cSmago;
    } else if (whichParameter == dynamicParams::dynamicOmega) {
        T rhoBar;
        Array<T, Descriptor<T>::d> j;
        Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
        momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
        T preFactor = SmagoOperations<T, Descriptor>::computePrefactor(this->getOmega(), cSmago);
        T omega = SmagoOperations<T, Descriptor>::computeOmega(
            this->getOmega(), preFactor, rhoBar, PiNeq);
        return omega;
    } else {
        return CompleteRegularizedTRTdynamics<T, Descriptor>::getParameter(whichParameter);
    }
}

template <typename T, template <typename U> class Descriptor>
int ConsistentSmagorinskyCompleteRegularizedTRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyCompleteRegularizedTRTdynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    CompleteRegularizedTRTdynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(cSmago);
    serializer.addValue(order);
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyCompleteRegularizedTRTdynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    CompleteRegularizedTRTdynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(cSmago);
    unserializer.readValue(order);
}

template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyCompleteRegularizedTRTdynamics<T, Descriptor>
    *ConsistentSmagorinskyCompleteRegularizedTRTdynamics<T, Descriptor>::clone() const
{
    return new ConsistentSmagorinskyCompleteRegularizedTRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyCompleteRegularizedTRTdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T uSqr = dynamicsTemplates<T, Descriptor>::
        consistent_smagorinsky_complete_regularized_mrt_ma2_collision(
            cell, order, cSmago, this->getOmega(), this->getPsi());

    if (cell.takesStatistics()) {
        T rhoBar = momentTemplates<T, Descriptor>::get_rhoBar(cell);
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyCompleteRegularizedTRTdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T rhoBarLb;
    Array<T, Descriptor<T>::d> jLb;
    Array<T, SymmetricTensor<T, Descriptor>::n> piNeqLb;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBarLb, jLb, piNeqLb);
    T jSqrLb = VectorTemplate<T, Descriptor>::normSqr(jLb);
    this->regularize(cell, rhoBarLb, jLb, jSqrLb, piNeqLb);

    T uSqr = dynamicsTemplates<T, Descriptor>::complete_mrt_smagorinsky_ma2_ext_rhoBar_j_collision(
        cell, rhoBar, j, order, cSmago, this->getOmega(), this->getPsi());

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

/* *************** Class ConsistentSmagorinskyTruncatedTRTdynamics
 * ************************************** */

template <typename T, template <typename U> class Descriptor>
int ConsistentSmagorinskyTruncatedTRTdynamics<T, Descriptor>::id = meta::registerGeneralDynamics<
    T, Descriptor, ConsistentSmagorinskyTruncatedTRTdynamics<T, Descriptor> >(
    "TruncatedTRT_Consistent_Smagorinsky");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyTruncatedTRTdynamics<T, Descriptor>::ConsistentSmagorinskyTruncatedTRTdynamics(
    T omega_, T psi_, T cSmago_) :
    TruncatedTRTdynamics<T, Descriptor>(omega_, psi_), cSmago(cSmago_)
{ }

template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyTruncatedTRTdynamics<T, Descriptor>::ConsistentSmagorinskyTruncatedTRTdynamics(
    T omega_, T cSmago_) :
    TruncatedTRTdynamics<T, Descriptor>(omega_), cSmago(cSmago_)
{ }

template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyTruncatedTRTdynamics<T, Descriptor>::ConsistentSmagorinskyTruncatedTRTdynamics(
    HierarchicUnserializer &unserializer) :
    TruncatedTRTdynamics<T, Descriptor>(T(), T()), cSmago(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
T ConsistentSmagorinskyTruncatedTRTdynamics<T, Descriptor>::getDynamicParameter(
    plint whichParameter, Cell<T, Descriptor> const &cell) const
{
    if (whichParameter == dynamicParams::smagorinskyConstant) {
        return cSmago;
    } else if (whichParameter == dynamicParams::dynamicOmega) {
        T rhoBar;
        Array<T, Descriptor<T>::d> j;
        Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
        momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
        T preFactor = SmagoOperations<T, Descriptor>::computePrefactor(this->getOmega(), cSmago);
        T omega = SmagoOperations<T, Descriptor>::computeOmega(
            this->getOmega(), preFactor, rhoBar, PiNeq);
        return omega;
    } else {
        return TruncatedTRTdynamics<T, Descriptor>::getParameter(whichParameter);
    }
}

template <typename T, template <typename U> class Descriptor>
int ConsistentSmagorinskyTruncatedTRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyTruncatedTRTdynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    TruncatedTRTdynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyTruncatedTRTdynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    TruncatedTRTdynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyTruncatedTRTdynamics<T, Descriptor>
    *ConsistentSmagorinskyTruncatedTRTdynamics<T, Descriptor>::clone() const
{
    return new ConsistentSmagorinskyTruncatedTRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyTruncatedTRTdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T uSqr = dynamicsTemplates<T, Descriptor>::truncated_mrt_smagorinsky_ma2_collision(
        cell, cSmago, this->getOmega(), this->getPsi());

    if (cell.takesStatistics()) {
        T rhoBar = momentTemplates<T, Descriptor>::get_rhoBar(cell);
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyTruncatedTRTdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T uSqr = dynamicsTemplates<T, Descriptor>::truncated_mrt_smagorinsky_ma2_ext_rhoBar_j_collision(
        cell, rhoBar, j, cSmago, this->getOmega(), this->getPsi());

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

/* *************** Class ConsistentSgsBGKdynamics ************************************** */

template <typename T, template <typename U> class Descriptor>
int ConsistentSgsBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, ConsistentSgsBGKdynamics<T, Descriptor> >(
        "BGK_Consistent_Sgs");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
ConsistentSgsBGKdynamics<T, Descriptor>::ConsistentSgsBGKdynamics(T omega0_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega0_), omega0(omega0_)
{ }

template <typename T, template <typename U> class Descriptor>
ConsistentSgsBGKdynamics<T, Descriptor>::ConsistentSgsBGKdynamics(
    HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T()), omega0(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
T ConsistentSgsBGKdynamics<T, Descriptor>::getOmega() const
{
    return omega0;
}

template <typename T, template <typename U> class Descriptor>
T ConsistentSgsBGKdynamics<T, Descriptor>::getDynamicParameter(
    plint whichParameter, Cell<T, Descriptor> const &cell) const
{
    //     if (whichParameter==dynamicParams::dynamicOmega) {
    //         T rhoBar;
    //         Array<T,Descriptor<T>::d> j;
    //         Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    //         momentTemplates<T,Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    //         T preFactor = (SmagoOperations<T,Descriptor>::computePrefactor(omega0,cSmago));
    //         T omega = SmagoOperations<T,Descriptor>::computeOmega (
    //                 omega0, preFactor, rhoBar, PiNeq );
    //         return omega;
    //     }
    //     else {
    return IsoThermalBulkDynamics<T, Descriptor>::getParameter(whichParameter);
    //     }
}

template <typename T, template <typename U> class Descriptor>
T ConsistentSgsBGKdynamics<T, Descriptor>::getParameter(plint whichParameter) const
{
    switch (whichParameter) {
    case dynamicParams::omega_shear:
        return omega0;
    default:
        return IsoThermalBulkDynamics<T, Descriptor>::getParameter(whichParameter);
    }
}

template <typename T, template <typename U> class Descriptor>
ConsistentSgsBGKdynamics<T, Descriptor> *ConsistentSgsBGKdynamics<T, Descriptor>::clone() const
{
    return new ConsistentSgsBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int ConsistentSgsBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSgsBGKdynamics<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    IsoThermalBulkDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(omega0);
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSgsBGKdynamics<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    IsoThermalBulkDynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(omega0);
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSgsBGKdynamics<T, Descriptor>::setOmega(T omega0_)
{
    omega0 = omega0_;
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSgsBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);

    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    PiNeq.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::tensorBeginsAt));
    Array<T, Descriptor<T>::q> sgsTerm =
        SmagoOperations<T, Descriptor>::computeSgsTensorTerm(PiNeq, omega0);

    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, j, omega0);
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] += sgsTerm[iPop];
    }

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSgsBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    PiNeq.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::tensorBeginsAt));
    Array<T, Descriptor<T>::q> sgsTerm =
        SmagoOperations<T, Descriptor>::computeSgsTensorTerm(PiNeq, omega0);

    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, j, omega0);
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] += sgsTerm[iPop];
    }

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T ConsistentSgsBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class GuoExternalForceSmagorinskyBGKdynamics
 * ************************************** */

template <typename T, template <typename U> class Descriptor>
int GuoExternalForceSmagorinskyBGKdynamics<T, Descriptor>::id = meta::registerGeneralDynamics<
    T, Descriptor, GuoExternalForceSmagorinskyBGKdynamics<T, Descriptor> >(
    "BGK_Smagorinsky_ExternalForce_Guo");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
GuoExternalForceSmagorinskyBGKdynamics<T, Descriptor>::GuoExternalForceSmagorinskyBGKdynamics(
    T omega0_, T cSmago_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega0_),
    omega0(omega0_),
    cSmago(cSmago_),
    preFactor(SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago))
{ }

template <typename T, template <typename U> class Descriptor>
GuoExternalForceSmagorinskyBGKdynamics<T, Descriptor>::GuoExternalForceSmagorinskyBGKdynamics(
    HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T()), omega0(T()), cSmago(T()), preFactor(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceSmagorinskyBGKdynamics<T, Descriptor>::setOmega(T omega0_)
{
    omega0 = omega0_;
    preFactor = SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago);
}

template <typename T, template <typename U> class Descriptor>
T GuoExternalForceSmagorinskyBGKdynamics<T, Descriptor>::getOmega() const
{
    return omega0;
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceSmagorinskyBGKdynamics<T, Descriptor>::setParameter(
    plint whichParameter, T value)
{
    switch (whichParameter) {
    case dynamicParams::omega_shear:
        omega0 = value;
        break;
    case dynamicParams::smagorinskyConstant:
        cSmago = value;
        break;
    default:
        IsoThermalBulkDynamics<T, Descriptor>::setParameter(whichParameter, value);
    }

    preFactor = SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago);
}

template <typename T, template <typename U> class Descriptor>
GuoExternalForceSmagorinskyBGKdynamics<T, Descriptor>
    *GuoExternalForceSmagorinskyBGKdynamics<T, Descriptor>::clone() const
{
    return new GuoExternalForceSmagorinskyBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int GuoExternalForceSmagorinskyBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceSmagorinskyBGKdynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    IsoThermalBulkDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(omega0);
    serializer.addValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceSmagorinskyBGKdynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    IsoThermalBulkDynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(omega0);
    unserializer.readValue(cSmago);
    preFactor = SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago);
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceSmagorinskyBGKdynamics<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> force, j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));

    T invRho = Descriptor<T>::invRho(rhoBar);
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        u[iD] = j[iD] * invRho + (T)0.5 * force[iD];
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceSmagorinskyBGKdynamics<T, Descriptor>::computeVelocityExternal(
    Cell<T, Descriptor> const &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j,
    Array<T, Descriptor<T>::d> &u) const
{
    Array<T, Descriptor<T>::d> force;
    force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));

    T invRho = Descriptor<T>::invRho(rhoBar);
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        u[iD] = j[iD] * invRho + (T)0.5 * force[iD];
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceSmagorinskyBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T omega = SmagoOperations<T, Descriptor>::computeOmega(omega0, preFactor, rhoBar, PiNeq);

    Array<T, Descriptor<T>::d> u;
    computeVelocityExternal(cell, rhoBar, j, u);
    T rho = Descriptor<T>::fullRho(rhoBar);
    j = rho * u;

    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, j, omega);
    externalForceTemplates<T, Descriptor>::addGuoForce(cell, u, omega, (T)1);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceSmagorinskyBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);
    T omega = SmagoOperations<T, Descriptor>::computeOmega(omega0, preFactor, rhoBar, PiNeq);

    Array<T, Descriptor<T>::d> u;
    computeVelocityExternal(cell, rhoBar, j, u);
    T rho = Descriptor<T>::fullRho(rhoBar);
    Array<T, Descriptor<T>::d> jP = rho * u;

    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, jP, omega);
    externalForceTemplates<T, Descriptor>::addGuoForce(cell, u, omega, (T)1);

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T GuoExternalForceSmagorinskyBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class SmagorinskyIncBGKdynamics ************************************** */

template <typename T, template <typename U> class Descriptor>
int SmagorinskyIncBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, SmagorinskyIncBGKdynamics<T, Descriptor> >(
        "IncBGK_Smagorinsky");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
SmagorinskyIncBGKdynamics<T, Descriptor>::SmagorinskyIncBGKdynamics(T omega0_, T cSmago_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega0_),
    omega0(omega0_),
    cSmago(cSmago_),
    preFactor(SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago))
{ }

template <typename T, template <typename U> class Descriptor>
SmagorinskyIncBGKdynamics<T, Descriptor>::SmagorinskyIncBGKdynamics(
    HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T()), omega0(T()), cSmago(T()), preFactor(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
SmagorinskyIncBGKdynamics<T, Descriptor> *SmagorinskyIncBGKdynamics<T, Descriptor>::clone() const
{
    return new SmagorinskyIncBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int SmagorinskyIncBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyIncBGKdynamics<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    IsoThermalBulkDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(omega0);
    serializer.addValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyIncBGKdynamics<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    IsoThermalBulkDynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(omega0);
    unserializer.readValue(cSmago);
    preFactor = SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago);
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyIncBGKdynamics<T, Descriptor>::setOmega(T omega0_)
{
    omega0 = omega0_;
    preFactor = SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago);
}

template <typename T, template <typename U> class Descriptor>
T SmagorinskyIncBGKdynamics<T, Descriptor>::getOmega() const
{
    return omega0;
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyIncBGKdynamics<T, Descriptor>::computeRhoBarJPiNeq(
    Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
    Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const
{
    // Incompressible: rho0=1
    T invRho0 = (T)1.;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq, invRho0);
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyIncBGKdynamics<T, Descriptor>::setParameter(plint whichParameter, T value)
{
    switch (whichParameter) {
    case dynamicParams::omega_shear:
        omega0 = value;
        break;
    case dynamicParams::smagorinskyConstant:
        cSmago = value;
        break;
    default:
        IsoThermalBulkDynamics<T, Descriptor>::setParameter(whichParameter, value);
    }
}

template <typename T, template <typename U> class Descriptor>
T SmagorinskyIncBGKdynamics<T, Descriptor>::getParameter(plint whichParameter) const
{
    switch (whichParameter) {
    case dynamicParams::omega_shear:
        return omega0;
    case dynamicParams::smagorinskyConstant:
        return cSmago;
    default:
        return IsoThermalBulkDynamics<T, Descriptor>::getParameter(whichParameter);
    }
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyIncBGKdynamics<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
{
    // Incompressible: rho0=1 => u=j
    T rhoBar;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, u);
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyIncBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T omega = SmagoOperations<T, Descriptor>::computeOmega(omega0, preFactor, rhoBar, PiNeq);
    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_inc_collision(cell, rhoBar, j, omega);
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyIncBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);
    T omega = SmagoOperations<T, Descriptor>::computeOmega(omega0, preFactor, rhoBar, PiNeq);
    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_inc_collision(cell, rhoBar, j, omega);
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T SmagorinskyIncBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    // Incompressible: rho0=1
    T invRho0 = (T)1;
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho0, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
bool SmagorinskyIncBGKdynamics<T, Descriptor>::velIsJ() const
{
    return true;
}

/* *************** Class SmagorinskyRegularizedDynamics ************************************ */

template <typename T, template <typename U> class Descriptor>
int SmagorinskyRegularizedDynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, SmagorinskyRegularizedDynamics<T, Descriptor> >(
        "Smagorinsky_Regularized");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
SmagorinskyRegularizedDynamics<T, Descriptor>::SmagorinskyRegularizedDynamics(
    T omega0_, T cSmago_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega0_),
    omega0(omega0_),
    cSmago(cSmago_),
    preFactor(SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago))
{ }

template <typename T, template <typename U> class Descriptor>
SmagorinskyRegularizedDynamics<T, Descriptor>::SmagorinskyRegularizedDynamics(
    HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T()), omega0(T()), cSmago(T()), preFactor(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyRegularizedDynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    IsoThermalBulkDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(omega0);
    serializer.addValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyRegularizedDynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    IsoThermalBulkDynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(omega0);
    unserializer.readValue(cSmago);
    preFactor = SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago);
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyRegularizedDynamics<T, Descriptor>::setOmega(T omega0_)
{
    omega0 = omega0_;
    preFactor = SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago);
}

template <typename T, template <typename U> class Descriptor>
T SmagorinskyRegularizedDynamics<T, Descriptor>::getOmega() const
{
    return omega0;
}

template <typename T, template <typename U> class Descriptor>
SmagorinskyRegularizedDynamics<T, Descriptor>
    *SmagorinskyRegularizedDynamics<T, Descriptor>::clone() const
{
    return new SmagorinskyRegularizedDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int SmagorinskyRegularizedDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyRegularizedDynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T omega = SmagoOperations<T, Descriptor>::computeOmega(omega0, preFactor, rhoBar, PiNeq);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T uSqr = dynamicsTemplates<T, Descriptor>::rlb_collision(cell, rhoBar, invRho, j, PiNeq, omega);
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyRegularizedDynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);
    T omega = SmagoOperations<T, Descriptor>::computeOmega(omega0, preFactor, rhoBar, PiNeq);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T uSqr = dynamicsTemplates<T, Descriptor>::rlb_collision(cell, rhoBar, invRho, j, PiNeq, omega);
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T SmagorinskyRegularizedDynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class SecuredSmagorinskyRegularizedDynamics ************************************
 */

template <typename T, template <typename U> class Descriptor>
int SecuredSmagorinskyRegularizedDynamics<T, Descriptor>::id = meta::registerGeneralDynamics<
    T, Descriptor, SecuredSmagorinskyRegularizedDynamics<T, Descriptor> >(
    "Secured_Smagorinsky_Regularized");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
SecuredSmagorinskyRegularizedDynamics<T, Descriptor>::SecuredSmagorinskyRegularizedDynamics(
    T omega0_, T cSmago_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega0_),
    omega0(omega0_),
    cSmago(cSmago_),
    preFactor(SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago))
{ }

template <typename T, template <typename U> class Descriptor>
SecuredSmagorinskyRegularizedDynamics<T, Descriptor>::SecuredSmagorinskyRegularizedDynamics(
    HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T()), omega0(T()), cSmago(T()), preFactor(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
void SecuredSmagorinskyRegularizedDynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    IsoThermalBulkDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(omega0);
    serializer.addValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void SecuredSmagorinskyRegularizedDynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    IsoThermalBulkDynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(omega0);
    unserializer.readValue(cSmago);
    preFactor = SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago);
}

template <typename T, template <typename U> class Descriptor>
void SecuredSmagorinskyRegularizedDynamics<T, Descriptor>::setOmega(T omega0_)
{
    omega0 = omega0_;
    preFactor = SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago);
}

template <typename T, template <typename U> class Descriptor>
T SecuredSmagorinskyRegularizedDynamics<T, Descriptor>::getOmega() const
{
    return omega0;
}

template <typename T, template <typename U> class Descriptor>
SecuredSmagorinskyRegularizedDynamics<T, Descriptor>
    *SecuredSmagorinskyRegularizedDynamics<T, Descriptor>::clone() const
{
    return new SecuredSmagorinskyRegularizedDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int SecuredSmagorinskyRegularizedDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void SecuredSmagorinskyRegularizedDynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);

    T softLimit = 0.2;  // If you change this, please also change it
                        // in collideExternal()
    T hardLimit = 0.6;
    for (plint i = 0; i < Descriptor<T>::d; ++i) {
        constrainValue(j[i], softLimit, hardLimit);
    }
    softLimit *= 0.1;
    hardLimit *= 0.1;
    for (plint i = 0; i < SymmetricTensor<T, Descriptor>::n; ++i) {
        constrainValue(PiNeq[i], softLimit, hardLimit);
    }

    T omega = SmagoOperations<T, Descriptor>::computeOmega(omega0, preFactor, rhoBar, PiNeq);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T uSqr = dynamicsTemplates<T, Descriptor>::rlb_collision(cell, rhoBar, invRho, j, PiNeq, omega);
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void SecuredSmagorinskyRegularizedDynamics<T, Descriptor>::constrainValue(
    T &value, T softLimit, T hardLimit)
{
    T fvalue = std::fabs(value);
    plint sign = fvalue > 0 ? +1 : -1;
    if (fvalue > softLimit) {
        if (fvalue > hardLimit) {
            value = softLimit * sign;
        } else {
            T diff = fvalue - softLimit;
            T correction = diff * diff / (hardLimit - softLimit);
            value -= correction * sign;
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void SecuredSmagorinskyRegularizedDynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);

    T rho = Descriptor<T>::fullRho(rhoBar);
    T softLimit = 0.2;  // If you change this, please also change it
                        // in collide()
    T hardLimit = 0.6;
    Array<T, 3> newJ(j);
    for (plint i = 0; i < Descriptor<T>::d; ++i) {
        constrainValue(newJ[i], softLimit, hardLimit);
    }
    softLimit *= 0.1;
    hardLimit *= 0.1;
    for (plint i = 0; i < SymmetricTensor<T, Descriptor>::n; ++i) {
        constrainValue(PiNeq[i], softLimit, hardLimit);
    }

    T omega = SmagoOperations<T, Descriptor>::computeOmega(omega0, preFactor, rhoBar, PiNeq);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T uSqr =
        dynamicsTemplates<T, Descriptor>::rlb_collision(cell, rhoBar, invRho, newJ, PiNeq, omega);
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T SecuredSmagorinskyRegularizedDynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class SmagorinskyMRTdynamics ************************************** */

template <typename T, template <typename U> class Descriptor>
int SmagorinskyMRTdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, SmagorinskyMRTdynamics<T, Descriptor> >(
        "MRT_Smagorinsky");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
SmagorinskyMRTdynamics<T, Descriptor>::SmagorinskyMRTdynamics(T omega0_, T cSmago_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega0_),
    omega0(omega0_),
    cSmago(cSmago_),
    preFactor(SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago))
{ }

template <typename T, template <typename U> class Descriptor>
SmagorinskyMRTdynamics<T, Descriptor>::SmagorinskyMRTdynamics(
    HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T()), omega0(T()), cSmago(T()), preFactor(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyMRTdynamics<T, Descriptor>::setOmega(T omega0_)
{
    omega0 = omega0_;
    preFactor = SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago);
}

template <typename T, template <typename U> class Descriptor>
T SmagorinskyMRTdynamics<T, Descriptor>::getOmega() const
{
    return omega0;
}

template <typename T, template <typename U> class Descriptor>
T SmagorinskyMRTdynamics<T, Descriptor>::getDynamicParameter(
    plint whichParameter, Cell<T, Descriptor> const &cell) const
{
    if (whichParameter == dynamicParams::dynamicOmega) {
        T rhoBar;
        Array<T, Descriptor<T>::d> j;
        Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
        momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
        T omega = SmagoOperations<T, Descriptor>::computeOmega(omega0, preFactor, rhoBar, PiNeq);
        return omega;
    } else if (whichParameter == dynamicParams::smagorinskyConstant) {
        return cSmago;
    } else {
        return IsoThermalBulkDynamics<T, Descriptor>::getParameter(whichParameter);
    }
}

template <typename T, template <typename U> class Descriptor>
SmagorinskyMRTdynamics<T, Descriptor> *SmagorinskyMRTdynamics<T, Descriptor>::clone() const
{
    return new SmagorinskyMRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int SmagorinskyMRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyMRTdynamics<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    IsoThermalBulkDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(omega0);
    serializer.addValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyMRTdynamics<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    IsoThermalBulkDynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(omega0);
    unserializer.readValue(cSmago);
    preFactor = SmagoOperations<T, Descriptor>::computePrefactor(omega0, cSmago);
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyMRTdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    typedef mrtTemplates<T, Descriptor> mrtTemp;

    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T omega = SmagoOperations<T, Descriptor>::computeOmega(omega0, preFactor, rhoBar, PiNeq);
    T jSqr = mrtTemp::mrtCollision(cell, omega);
    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(statistics, rhoBar, jSqr * invRho * invRho);
    }
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyMRTdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    typedef mrtTemplates<T, Descriptor> mrtTemp;

    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);
    T omega = SmagoOperations<T, Descriptor>::computeOmega(omega0, preFactor, rhoBar, PiNeq);
    T jSqr = mrtTemp::mrtCollision(cell, rhoBar, j, omega);
    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(stat, rhoBar, jSqr * invRho * invRho);
    }
}

template <typename T, template <typename U> class Descriptor>
T SmagorinskyMRTdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class ConsistentSmagorinskyMRTdynamics
 * *********************************************** */

template <typename T, template <typename U> class Descriptor>
int ConsistentSmagorinskyMRTdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, ConsistentSmagorinskyMRTdynamics<T, Descriptor> >(
        "Consistent_MRT_Smagorinsky");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyMRTdynamics<T, Descriptor>::ConsistentSmagorinskyMRTdynamics(
    T omega0_, T cSmago_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega0_), cSmago(cSmago_)
{ }

template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyMRTdynamics<T, Descriptor>::ConsistentSmagorinskyMRTdynamics(
    HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T()), cSmago(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
T ConsistentSmagorinskyMRTdynamics<T, Descriptor>::getDynamicParameter(
    plint whichParameter, Cell<T, Descriptor> const &cell) const
{
    if (whichParameter == dynamicParams::smagorinskyConstant) {
        return cSmago;
    } else if (whichParameter == dynamicParams::dynamicOmega) {
        T rhoBar;
        Array<T, Descriptor<T>::d> j;
        Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
        momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
        T preFactor = (SmagoOperations<T, Descriptor>::computePrefactor(this->getOmega(), cSmago));
        T omega = SmagoOperations<T, Descriptor>::computeOmega(
            this->getOmega(), preFactor, rhoBar, PiNeq);
        return omega;
    } else {
        return IsoThermalBulkDynamics<T, Descriptor>::getParameter(whichParameter);
    }
}

template <typename T, template <typename U> class Descriptor>
T ConsistentSmagorinskyMRTdynamics<T, Descriptor>::getParameter(plint whichParameter) const
{
    switch (whichParameter) {
    case dynamicParams::smagorinskyConstant:
        return cSmago;
    default:
        return IsoThermalBulkDynamics<T, Descriptor>::getParameter(whichParameter);
    }
}

template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyMRTdynamics<T, Descriptor>
    *ConsistentSmagorinskyMRTdynamics<T, Descriptor>::clone() const
{
    return new ConsistentSmagorinskyMRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int ConsistentSmagorinskyMRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyMRTdynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    IsoThermalBulkDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyMRTdynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    IsoThermalBulkDynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyMRTdynamics<T, Descriptor>::collide(
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
    if (normPiNeq != T()) {  // test to avoid division per 0
        S = (-rho * tau * Descriptor<T>::cs2
             + std::sqrt(util::sqr(rho * tau * Descriptor<T>::cs2) + normPiNeq))
            / normPiNeq * PiNeq;
    }

    T jSqr = mrtTemp::smagorinskyMrtCollision(cell, rhoBar, j, S, cSmago, this->getOmega());

    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(statistics, rhoBar, jSqr * invRho * invRho);
    }
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyMRTdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    typedef mrtTemplates<T, Descriptor> mrtTemp;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);
    T rho = Descriptor<T>::fullRho(rhoBar);

    T tau = (T)1 / this->getOmega();

    T normPiNeq = std::sqrt((T)2 * SymmetricTensor<T, Descriptor>::tensorNormSqr(PiNeq));
    normPiNeq *= (T)2 * rho * util::sqr(Descriptor<T>::cs2 * cSmago * cSmago);

    Array<T, SymmetricTensor<T, Descriptor>::n> S;
    S.resetToZero();
    if (normPiNeq != T()) {  // test to avoid division per 0
        S = (-rho * tau * Descriptor<T>::cs2
             + std::sqrt(util::sqr(rho * tau * Descriptor<T>::cs2) + normPiNeq))
            / normPiNeq * PiNeq;
    }

    T jSqr = mrtTemp::smagorinskyMrtCollision(cell, rhoBar, j, S, cSmago, this->getOmega());

    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(stat, rhoBar, jSqr * invRho * invRho);
    }
}

template <typename T, template <typename U> class Descriptor>
T ConsistentSmagorinskyMRTdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class ConsistentSgsMRTdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
int ConsistentSgsMRTdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, ConsistentSgsMRTdynamics<T, Descriptor> >(
        "Consistent_MRT_Sgs");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
ConsistentSgsMRTdynamics<T, Descriptor>::ConsistentSgsMRTdynamics(T omega0_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega0_)
{ }

template <typename T, template <typename U> class Descriptor>
ConsistentSgsMRTdynamics<T, Descriptor>::ConsistentSgsMRTdynamics(
    HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
T ConsistentSgsMRTdynamics<T, Descriptor>::getDynamicParameter(
    plint whichParameter, Cell<T, Descriptor> const &cell) const
{
    //     if (whichParameter==dynamicParams::dynamicOmega) {
    //         T rhoBar;
    //         Array<T,Descriptor<T>::d> j;
    //         Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    //         momentTemplates<T,Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    //         T preFactor = SmagoOperations<T,Descriptor>::computePrefactor(omega0,cSmago);
    //         T omega = SmagoOperations<T,Descriptor>::computeOmega (
    //                 this->getOmega(), preFactor, rhoBar, PiNeq );
    //         return omega;
    //     }
    //     else {
    return IsoThermalBulkDynamics<T, Descriptor>::getParameter(whichParameter);
    //     }
}

template <typename T, template <typename U> class Descriptor>
T ConsistentSgsMRTdynamics<T, Descriptor>::getParameter(plint whichParameter) const
{
    return IsoThermalBulkDynamics<T, Descriptor>::getParameter(whichParameter);
}

template <typename T, template <typename U> class Descriptor>
ConsistentSgsMRTdynamics<T, Descriptor> *ConsistentSgsMRTdynamics<T, Descriptor>::clone() const
{
    return new ConsistentSgsMRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int ConsistentSgsMRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSgsMRTdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    typedef mrtTemplates<T, Descriptor> mrtTemp;
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);

    T jSqr = mrtTemp::mrtCollision(cell, rhoBar, j, this->getOmega());
    Array<T, Descriptor<T>::q> sgs, mNeq;
    sgs.resetToZero();
    mNeq.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::tensorBeginsAt));
    mrtTemplatesImpl<T, typename Descriptor<T>::SecondBaseDescriptor>::computef_InvM_Smoments(
        sgs, mNeq, this->getOmega());

    for (pluint iPop = 0; iPop < Descriptor<T>::q; ++iPop)
        cell[iPop] += sgs[iPop];

    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(statistics, rhoBar, jSqr * invRho * invRho);
    }
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSgsMRTdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    typedef mrtTemplates<T, Descriptor> mrtTemp;

    T jSqr = mrtTemp::mrtCollision(cell, rhoBar, j, this->getOmega());

    Array<T, Descriptor<T>::q> sgs, mNeq;
    sgs.resetToZero();
    mNeq.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::tensorBeginsAt));
    mrtTemplatesImpl<T, typename Descriptor<T>::SecondBaseDescriptor>::computef_InvM_Smoments(
        sgs, mNeq, this->getOmega());

    for (pluint iPop = 0; iPop < Descriptor<T>::q; ++iPop)
        cell[iPop] += sgs[iPop];

    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(stat, rhoBar, jSqr * invRho * invRho);
    }
}

template <typename T, template <typename U> class Descriptor>
T ConsistentSgsMRTdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class GuoExternalForceConsistentSgsMRTdynamics
 * *********************************************** */

template <typename T, template <typename U> class Descriptor>
int GuoExternalForceConsistentSgsMRTdynamics<T, Descriptor>::id = meta::registerGeneralDynamics<
    T, Descriptor, GuoExternalForceConsistentSgsMRTdynamics<T, Descriptor> >(
    "Guo_Consistent_MRT_Sgs");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
GuoExternalForceConsistentSgsMRTdynamics<T, Descriptor>::GuoExternalForceConsistentSgsMRTdynamics(
    T omega0_) :
    ExternalForceDynamics<T, Descriptor>(omega0_)
{ }

template <typename T, template <typename U> class Descriptor>
GuoExternalForceConsistentSgsMRTdynamics<T, Descriptor>::GuoExternalForceConsistentSgsMRTdynamics(
    HierarchicUnserializer &unserializer) :
    ExternalForceDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
GuoExternalForceConsistentSgsMRTdynamics<T, Descriptor>
    *GuoExternalForceConsistentSgsMRTdynamics<T, Descriptor>::clone() const
{
    return new GuoExternalForceConsistentSgsMRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int GuoExternalForceConsistentSgsMRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void GuoExternalForceConsistentSgsMRTdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    typedef mrtTemplates<T, Descriptor> mrtTemp;
    T rhoBar = momentTemplates<T, Descriptor>::get_rhoBar(cell);

    Array<T, Descriptor<T>::d> u;
    this->computeVelocity(cell, u);

    T jSqr = mrtTemp::mrtCollisionWithForce(cell, rhoBar, u, this->getOmega(), (T)1);

    Array<T, Descriptor<T>::q> sgs, mNeq;
    sgs.resetToZero();
    mNeq.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::tensorBeginsAt));
    mrtTemplatesImpl<T, typename Descriptor<T>::SecondBaseDescriptor>::computef_InvM_Smoments(
        sgs, mNeq, this->getOmega());

    for (pluint iPop = 0; iPop < Descriptor<T>::q; ++iPop)
        cell[iPop] += sgs[iPop];

    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(statistics, rhoBar, jSqr * invRho * invRho);
    }
}

template <typename T, template <typename U> class Descriptor>
T GuoExternalForceConsistentSgsMRTdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class ConsistentSmagorinskyIncMRTdynamics
 * *********************************************** */

template <typename T, template <typename U> class Descriptor>
int ConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::id = meta::registerGeneralDynamics<
    T, Descriptor, ConsistentSmagorinskyIncMRTdynamics<T, Descriptor> >(
    "Consistent_IncMRT_Smagorinsky");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::ConsistentSmagorinskyIncMRTdynamics(
    T omega0_, T cSmago_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega0_), cSmago(cSmago_)
{ }

template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::ConsistentSmagorinskyIncMRTdynamics(
    HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T()), cSmago(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
T ConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::getDynamicParameter(
    plint whichParameter, Cell<T, Descriptor> const &cell) const
{
    if (whichParameter == dynamicParams::smagorinskyConstant) {
        return cSmago;
    } else if (whichParameter == dynamicParams::dynamicOmega) {
        return this->getOmega();
    } else {
        return IsoThermalBulkDynamics<T, Descriptor>::getParameter(whichParameter);
    }
}

template <typename T, template <typename U> class Descriptor>
T ConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::getParameter(plint whichParameter) const
{
    switch (whichParameter) {
    case dynamicParams::smagorinskyConstant:
        return cSmago;
    default:
        return IsoThermalBulkDynamics<T, Descriptor>::getParameter(whichParameter);
    }
}

template <typename T, template <typename U> class Descriptor>
ConsistentSmagorinskyIncMRTdynamics<T, Descriptor>
    *ConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::clone() const
{
    return new ConsistentSmagorinskyIncMRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int ConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    IsoThermalBulkDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    IsoThermalBulkDynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
{
    T dummyRhoBar;
    this->computeRhoBarJ(cell, dummyRhoBar, u);
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::computeRhoBarJPiNeq(
    Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
    Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const
{
    // Incompressible: rho0=1
    T invRho0 = (T)1.;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq, invRho0);
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::collide(
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
    if (normPiNeq != T()) {  // test to avoid division per 0
        S = (-rho * tau * Descriptor<T>::cs2
             + std::sqrt(util::sqr(rho * tau * Descriptor<T>::cs2) + normPiNeq))
            / normPiNeq * PiNeq;
    }

    T jSqr = mrtTemp::incSmagorinskyMrtCollision(cell, rhoBar, j, S, cSmago, this->getOmega());

    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(statistics, rhoBar, jSqr * invRho * invRho);
    }
}

template <typename T, template <typename U> class Descriptor>
void ConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    typedef mrtTemplates<T, Descriptor> mrtTemp;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);
    T rho = Descriptor<T>::fullRho(rhoBar);

    T tau = (T)1 / this->getOmega();

    T normPiNeq = std::sqrt((T)2 * SymmetricTensor<T, Descriptor>::tensorNormSqr(PiNeq));
    normPiNeq *= (T)2 * rho * util::sqr(Descriptor<T>::cs2 * cSmago * cSmago);

    Array<T, SymmetricTensor<T, Descriptor>::n> S;
    S.resetToZero();
    if (normPiNeq != T()) {  // test to avoid division per 0
        S = (-rho * tau * Descriptor<T>::cs2
             + std::sqrt(util::sqr(rho * tau * Descriptor<T>::cs2) + normPiNeq))
            / normPiNeq * PiNeq;
    }

    T jSqr = mrtTemp::incSmagorinskyMrtCollision(cell, rhoBar, j, S, cSmago, this->getOmega());

    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(stat, rhoBar, jSqr * invRho * invRho);
    }
}

template <typename T, template <typename U> class Descriptor>
T ConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, (T)1, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
bool ConsistentSmagorinskyIncMRTdynamics<T, Descriptor>::velIsJ() const
{
    return true;
}

}  // namespace plb

#endif  // SMAGORINSKY_DYNAMICS_HH
