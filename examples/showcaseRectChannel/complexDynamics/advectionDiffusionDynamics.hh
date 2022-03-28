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
#ifndef ADVECTION_DIFFUSION_DYNAMICS_HH
#define ADVECTION_DIFFUSION_DYNAMICS_HH

#include "complexDynamics/advectionDiffusionDynamics.h"
#include "core/dynamicsIdentifiers.h"
#include "core/latticeStatistics.h"
#include "latticeBoltzmann/advectionDiffusionDynamicsTemplates.h"
#include "latticeBoltzmann/advectionDiffusionMomentTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/offEquilibriumAdvectionDiffusionTemplates.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
struct AD_SmagoOperations {
    static T computePrefactor(T omega0, T cSmago)
    {
        return util::sqr(cSmago * omega0 * Descriptor<T>::invCs2);
    }
    static T computeOmega(T omega0, T alpha, Array<T, Descriptor<T>::d> const &j1)
    {
        T j1Norm = norm(j1);
        T linearTerm = alpha * j1Norm;
        T squareTerm = linearTerm * linearTerm;
        // In the following formula, the square-root appearing in the explicit form of
        //   omega is developed to second-order.
        return omega0 * (1 - linearTerm + squareTerm);
    }
};

/* *************** Class AdvectionDiffusionDynamics ************************************ */

template <typename T, template <typename U> class Descriptor>
AdvectionDiffusionDynamics<T, Descriptor>::AdvectionDiffusionDynamics(T omega_) :
    BasicBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
void AdvectionDiffusionDynamics<T, Descriptor>::regularize(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar) const
{
    // jAdvDiff is the first order moment of
    Array<T, Descriptor<T>::d> jEq;

    advectionDiffusionMomentTemplates<T, Descriptor>::get_jEq(cell, rhoBar, jEq);

    advectionDiffusionDynamicsTemplates<T, Descriptor>::regularize(cell, rhoBar, j, jEq);
}

template <typename T, template <typename U> class Descriptor>
T AdvectionDiffusionDynamics<T, Descriptor>::computeEbar(Cell<T, Descriptor> const &cell) const
{
    return T();
}

template <typename T, template <typename U> class Descriptor>
void AdvectionDiffusionDynamics<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
{
    const T *u_ = cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt);
    for (plint i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = u_[i];
    }
}

template <typename T, template <typename U> class Descriptor>
plint AdvectionDiffusionDynamics<T, Descriptor>::numDecomposedVariables(plint order) const
{
    // Start with the decomposed version of the populations.
    plint numVariables = 1 + Descriptor<T>::d + Descriptor<T>::q;

    // Add the variables in the external scalars.
    numVariables += Descriptor<T>::ExternalField::numScalars;
    return numVariables;
}

template <typename T, template <typename U> class Descriptor>
void AdvectionDiffusionDynamics<T, Descriptor>::decompose(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData, plint order) const
{
    rawData.resize(numDecomposedVariables(order));

    if (order == 0) {
        T rhoBar;
        Array<T, Descriptor<T>::d> j;
        advectionDiffusionMomentTemplates<T, Descriptor>::get_rhoBar_jEq(cell, rhoBar, j);
        T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
        rawData[0] = rhoBar;
        j.to_cArray(&rawData[1]);

        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            rawData[1 + Descriptor<T>::d + iPop] =
                cell[iPop] - this->computeEquilibrium(iPop, rhoBar, j, jSqr);
        }

        int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
        for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
            rawData[offset + iExt] = *cell.getExternal(iExt);
        }
    } else {
        PLB_ASSERT(false);
        // there is no PiNeq in Advection/Diffusion - let's hope, this is never needed!
    }
}

template <typename T, template <typename U> class Descriptor>
void AdvectionDiffusionDynamics<T, Descriptor>::recompose(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData, plint order) const
{
    PLB_PRECONDITION((plint)rawData.size() == numDecomposedVariables(order));

    if (order == 0) {
        T rhoBar = rawData[0];
        Array<T, Descriptor<T>::d> j;
        j.from_cArray(&rawData[1]);
        T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            cell[iPop] = this->computeEquilibrium(iPop, rhoBar, j, jSqr)
                         + rawData[1 + Descriptor<T>::d + iPop];
        }

        int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
        for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
            *cell.getExternal(iExt) = rawData[offset + iExt];
        }
    } else {
        PLB_ASSERT(false);
        // there is no PiNeq in Advection/Diffusion - let's hope, this is never needed!
    }
}

template <typename T, template <typename U> class Descriptor>
void AdvectionDiffusionDynamics<T, Descriptor>::rescale(
    std::vector<T> &rawData, T xDxInv, T xDt, plint order) const
{
    // TODO: rescale velcotiy in the external scalars.
}

/* *************** Class SmagorinskyAdvectionDiffusionRLBdynamics *************** */

template <typename T, template <typename U> class Descriptor>
int SmagorinskyAdvectionDiffusionRLBdynamics<T, Descriptor>::id = meta::registerGeneralDynamics<
    T, Descriptor, SmagorinskyAdvectionDiffusionRLBdynamics<T, Descriptor> >(
    "SmagoAdvectionDiffusion_RLB");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
SmagorinskyAdvectionDiffusionRLBdynamics<T, Descriptor>::SmagorinskyAdvectionDiffusionRLBdynamics(
    T omega_, T T0_, T cSmago_) :
    AdvectionDiffusionDynamics<T, Descriptor>(omega_), invT0((T)1 / T0_), cSmago(cSmago_)
{ }

template <typename T, template <typename U> class Descriptor>
SmagorinskyAdvectionDiffusionRLBdynamics<T, Descriptor>::SmagorinskyAdvectionDiffusionRLBdynamics(
    HierarchicUnserializer &unserializer) :
    AdvectionDiffusionDynamics<T, Descriptor>((T)1)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
SmagorinskyAdvectionDiffusionRLBdynamics<T, Descriptor>
    *SmagorinskyAdvectionDiffusionRLBdynamics<T, Descriptor>::clone() const
{
    return new SmagorinskyAdvectionDiffusionRLBdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int SmagorinskyAdvectionDiffusionRLBdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyAdvectionDiffusionRLBdynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    AdvectionDiffusionDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(invT0);
    serializer.addValue(cSmago);
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyAdvectionDiffusionRLBdynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    AdvectionDiffusionDynamics<T, Descriptor>::unserialize(unserializer);
    invT0 = unserializer.readValue<T>();
    cSmago = unserializer.readValue<T>();
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyAdvectionDiffusionRLBdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> jEq, jNeq;
    advectionDiffusionMomentTemplates<T, Descriptor>::get_rhoBar_jEq_jNeq(cell, rhoBar, jEq, jNeq);

    T omega0 = this->getOmega();
    T omega = omega0
              / (1.
                 + (1. - omega0 / 2.0) * cSmago * cSmago * omega0 * Descriptor<T>::invCs2 * invT0
                       * norm(jNeq));

    T uSqr = advectionDiffusionDynamicsTemplates<T, Descriptor>::no_corr_rlb_collision(
        cell, rhoBar, jEq, jNeq, omega);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void SmagorinskyAdvectionDiffusionRLBdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T thetaBar,
    BlockStatistics &statistics)
{
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_j(cell, j);
    Array<T, Descriptor<T>::d> jNeq(j - jEq);

    T omega0 = this->getOmega();
    T omega = omega0
              / (1.
                 + (1. - omega0 / 2.0) * cSmago * cSmago * omega0 * Descriptor<T>::invCs2 * invT0
                       * norm(jNeq));

    T uSqr = advectionDiffusionDynamicsTemplates<T, Descriptor>::no_corr_rlb_collision(
        cell, rhoBar, jEq, jNeq, omega);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

/** \param jEq The parameter jEq is defined as jEq = j_advDiff = rho_advDiff*u_fluid
 */
template <typename T, template <typename U> class Descriptor>
T SmagorinskyAdvectionDiffusionRLBdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T jSqr, T thetaBar) const
{
    return advectionDiffusionDynamicsTemplates<T, Descriptor>::bgk_ma1_equilibrium(
        iPop, rhoBar, jEq);
}

/* *************** Class AdvectionDiffusionPerkoDynamics *************** */

template <typename T, template <typename U> class Descriptor>
int AdvectionDiffusionPerkoDynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, AdvectionDiffusionPerkoDynamics<T, Descriptor> >(
        "AdvectionDiffusion_Perko");

template <typename T, template <typename U> class Descriptor>
AdvectionDiffusionPerkoDynamics<T, Descriptor>::AdvectionDiffusionPerkoDynamics(
    T omega_, T omegaRef_) :
    AdvectionDiffusionDynamics<T, Descriptor>(omega_), omegaRef(omegaRef_)
{ }

template <typename T, template <typename U> class Descriptor>
AdvectionDiffusionPerkoDynamics<T, Descriptor>::AdvectionDiffusionPerkoDynamics(
    HierarchicUnserializer &unserializer) :
    AdvectionDiffusionDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
void AdvectionDiffusionPerkoDynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    AdvectionDiffusionDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(omegaRef);
}

template <typename T, template <typename U> class Descriptor>
void AdvectionDiffusionPerkoDynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    AdvectionDiffusionDynamics<T, Descriptor>::unserialize(unserializer);
    omegaRef = unserializer.readValue<T>();
}

template <typename T, template <typename U> class Descriptor>
AdvectionDiffusionPerkoDynamics<T, Descriptor>
    *AdvectionDiffusionPerkoDynamics<T, Descriptor>::clone() const
{
    return new AdvectionDiffusionPerkoDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int AdvectionDiffusionPerkoDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void AdvectionDiffusionPerkoDynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> jEq, jNeq;
    advectionDiffusionMomentTemplates<T, Descriptor>::get_rhoBar_jEq_jNeq(cell, rhoBar, jEq, jNeq);

    T uSqr = advectionDiffusionDynamicsTemplates<T, Descriptor>::no_corr_bgk_collision(
        cell, rhoBar, jEq, omegaRef);

    T kappa = (T)1 - this->getOmega() / omegaRef;

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        T correction = omegaRef * Descriptor<T>::t[iPop] * kappa * Descriptor<T>::invCs2
                       * intDot<T, Descriptor<T>::d>(Descriptor<T>::c[iPop], jNeq);
        cell[iPop] += correction;
    }

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void AdvectionDiffusionPerkoDynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T thetaBar,
    BlockStatistics &statistics)
{
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_j(cell, j);
    Array<T, Descriptor<T>::d> jNeq(j - jEq);

    T uSqr = advectionDiffusionDynamicsTemplates<T, Descriptor>::no_corr_bgk_collision(
        cell, rhoBar, jEq, omegaRef);

    T kappa = (T)1 - this->getOmega() / omegaRef;

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        T correction = omegaRef * Descriptor<T>::t[iPop] * kappa * Descriptor<T>::invCs2
                       * intDot(Descriptor<T>::c[iPop], jNeq);
        cell[iPop] += correction;
    }

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

/** \param jEq The parameter jEq is defined as jEq = j_advDiff = rho_advDiff*u_fluid
 */
template <typename T, template <typename U> class Descriptor>
T AdvectionDiffusionPerkoDynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T jSqr, T thetaBar) const
{
    return advectionDiffusionDynamicsTemplates<T, Descriptor>::bgk_ma1_equilibrium(
        iPop, rhoBar, jEq);
}

/* *************** Class AdvectionDiffusionRLBdynamics *************** */

template <typename T, template <typename U> class Descriptor>
int AdvectionDiffusionRLBdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, AdvectionDiffusionRLBdynamics<T, Descriptor> >(
        "AdvectionDiffusion_RLB");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
AdvectionDiffusionRLBdynamics<T, Descriptor>::AdvectionDiffusionRLBdynamics(T omega_) :
    AdvectionDiffusionDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
AdvectionDiffusionRLBdynamics<T, Descriptor>::AdvectionDiffusionRLBdynamics(
    HierarchicUnserializer &unserializer) :
    AdvectionDiffusionDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
AdvectionDiffusionRLBdynamics<T, Descriptor> *AdvectionDiffusionRLBdynamics<T, Descriptor>::clone()
    const
{
    return new AdvectionDiffusionRLBdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int AdvectionDiffusionRLBdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void AdvectionDiffusionRLBdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> jEq, jNeq;
    advectionDiffusionMomentTemplates<T, Descriptor>::get_rhoBar_jEq_jNeq(cell, rhoBar, jEq, jNeq);

    T uSqr = advectionDiffusionDynamicsTemplates<T, Descriptor>::no_corr_rlb_collision(
        cell, rhoBar, jEq, jNeq, this->getOmega());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void AdvectionDiffusionRLBdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T thetaBar,
    BlockStatistics &statistics)
{
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_j(cell, j);
    Array<T, Descriptor<T>::d> jNeq(j - jEq);

    T uSqr = advectionDiffusionDynamicsTemplates<T, Descriptor>::no_corr_rlb_collision(
        cell, rhoBar, jEq, jNeq, this->getOmega());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

/** \param jEq The parameter jEq is defined as jEq = j_advDiff = rho_advDiff*u_fluid
 */
template <typename T, template <typename U> class Descriptor>
T AdvectionDiffusionRLBdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T jSqr, T thetaBar) const
{
    return advectionDiffusionDynamicsTemplates<T, Descriptor>::bgk_ma1_equilibrium(
        iPop, rhoBar, jEq);
}

/* *************** Class AdvectionDiffusionWithSourceLinearRLBdynamics *************** */

template <typename T, template <typename U> class Descriptor>
int AdvectionDiffusionWithSourceLinearRLBdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<
        T, Descriptor, AdvectionDiffusionWithSourceLinearRLBdynamics<T, Descriptor> >(
        "AdvectionDiffusionWithSourceLinear_RLB");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
AdvectionDiffusionWithSourceLinearRLBdynamics<
    T, Descriptor>::AdvectionDiffusionWithSourceLinearRLBdynamics(T omega_) :
    AdvectionDiffusionDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
AdvectionDiffusionWithSourceLinearRLBdynamics<T, Descriptor>::
    AdvectionDiffusionWithSourceLinearRLBdynamics(HierarchicUnserializer &unserializer) :
    AdvectionDiffusionDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
AdvectionDiffusionWithSourceLinearRLBdynamics<T, Descriptor>
    *AdvectionDiffusionWithSourceLinearRLBdynamics<T, Descriptor>::clone() const
{
    return new AdvectionDiffusionWithSourceLinearRLBdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int AdvectionDiffusionWithSourceLinearRLBdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void AdvectionDiffusionWithSourceLinearRLBdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> jEq, jNeq;
    advectionDiffusionMomentTemplates<T, Descriptor>::get_rhoBar_jEq_jNeq_linear(
        cell, rhoBar, jEq, jNeq);

    T sourceTerm = *cell.getExternal(Descriptor<T>::ExternalField::scalarBeginsAt);

    T uSqr = advectionDiffusionDynamicsTemplates<T, Descriptor>::no_corr_rlb_collision(
        cell, rhoBar, jEq, jNeq, this->getOmega(), sourceTerm);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void AdvectionDiffusionWithSourceLinearRLBdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T thetaBar,
    BlockStatistics &statistics)
{
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_j(cell, j);
    Array<T, Descriptor<T>::d> jNeq(j - jEq);

    T sourceTerm = *cell.getExternal(Descriptor<T>::ExternalField::scalarBeginsAt);

    T uSqr = advectionDiffusionDynamicsTemplates<T, Descriptor>::no_corr_rlb_collision(
        cell, rhoBar, jEq, jNeq, this->getOmega(), sourceTerm);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

/** \param jEq The parameter jEq is defined as jEq = j_advDiff = rho_advDiff*u_fluid
 */
template <typename T, template <typename U> class Descriptor>
T AdvectionDiffusionWithSourceLinearRLBdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T jSqr, T thetaBar) const
{
    return advectionDiffusionDynamicsTemplates<T, Descriptor>::bgk_ma1_equilibrium(
        iPop, rhoBar, jEq);
}

/* *************** Class AdvectionDiffusionWithSourceRLBdynamics *************** */

template <typename T, template <typename U> class Descriptor>
int AdvectionDiffusionWithSourceRLBdynamics<T, Descriptor>::id = meta::registerGeneralDynamics<
    T, Descriptor, AdvectionDiffusionWithSourceRLBdynamics<T, Descriptor> >(
    "AdvectionDiffusionWithSource_RLB");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
AdvectionDiffusionWithSourceRLBdynamics<T, Descriptor>::AdvectionDiffusionWithSourceRLBdynamics(
    T omega_) :
    AdvectionDiffusionDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
AdvectionDiffusionWithSourceRLBdynamics<T, Descriptor>::AdvectionDiffusionWithSourceRLBdynamics(
    HierarchicUnserializer &unserializer) :
    AdvectionDiffusionDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
AdvectionDiffusionWithSourceRLBdynamics<T, Descriptor>
    *AdvectionDiffusionWithSourceRLBdynamics<T, Descriptor>::clone() const
{
    return new AdvectionDiffusionWithSourceRLBdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int AdvectionDiffusionWithSourceRLBdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void AdvectionDiffusionWithSourceRLBdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> jEq, jNeq;
    advectionDiffusionMomentTemplates<T, Descriptor>::get_rhoBar_jEq_jNeq(cell, rhoBar, jEq, jNeq);

    T sourceTerm = *cell.getExternal(Descriptor<T>::ExternalField::scalarBeginsAt);

    T uSqr = advectionDiffusionDynamicsTemplates<T, Descriptor>::no_corr_rlb_collision(
        cell, rhoBar, jEq, jNeq, this->getOmega(), sourceTerm);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void AdvectionDiffusionWithSourceRLBdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T thetaBar,
    BlockStatistics &statistics)
{
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_j(cell, j);
    Array<T, Descriptor<T>::d> jNeq(j - jEq);

    T sourceTerm = *cell.getExternal(Descriptor<T>::ExternalField::scalarBeginsAt);

    T uSqr = advectionDiffusionDynamicsTemplates<T, Descriptor>::no_corr_rlb_collision(
        cell, rhoBar, jEq, jNeq, this->getOmega(), sourceTerm);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

/** \param jEq The parameter jEq is defined as jEq = j_advDiff = rho_advDiff*u_fluid
 */
template <typename T, template <typename U> class Descriptor>
T AdvectionDiffusionWithSourceRLBdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T jSqr, T thetaBar) const
{
    return advectionDiffusionDynamicsTemplates<T, Descriptor>::bgk_ma1_equilibrium(
        iPop, rhoBar, jEq);
}

/* *************** Class AdvectionDiffusionWithSourcePerkoDynamics *************** */

template <typename T, template <typename U> class Descriptor>
int AdvectionDiffusionWithSourcePerkoDynamics<T, Descriptor>::id = meta::registerGeneralDynamics<
    T, Descriptor, AdvectionDiffusionWithSourcePerkoDynamics<T, Descriptor> >(
    "AdvectionDiffusionWithSource_Perko");

template <typename T, template <typename U> class Descriptor>
AdvectionDiffusionWithSourcePerkoDynamics<T, Descriptor>::AdvectionDiffusionWithSourcePerkoDynamics(
    T omega_, T omegaRef_) :
    AdvectionDiffusionDynamics<T, Descriptor>(omega_), omegaRef(omegaRef_)
{ }

template <typename T, template <typename U> class Descriptor>
AdvectionDiffusionWithSourcePerkoDynamics<T, Descriptor>::AdvectionDiffusionWithSourcePerkoDynamics(
    HierarchicUnserializer &unserializer) :
    AdvectionDiffusionDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
void AdvectionDiffusionWithSourcePerkoDynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    AdvectionDiffusionDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(omegaRef);
}

template <typename T, template <typename U> class Descriptor>
void AdvectionDiffusionWithSourcePerkoDynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    AdvectionDiffusionDynamics<T, Descriptor>::unserialize(unserializer);
    omegaRef = unserializer.readValue<T>();
}

template <typename T, template <typename U> class Descriptor>
AdvectionDiffusionWithSourcePerkoDynamics<T, Descriptor>
    *AdvectionDiffusionWithSourcePerkoDynamics<T, Descriptor>::clone() const
{
    return new AdvectionDiffusionWithSourcePerkoDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int AdvectionDiffusionWithSourcePerkoDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void AdvectionDiffusionWithSourcePerkoDynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> jEq, jNeq;
    advectionDiffusionMomentTemplates<T, Descriptor>::get_rhoBar_jEq_jNeq(cell, rhoBar, jEq, jNeq);

    T sourceTerm = *cell.getExternal(Descriptor<T>::ExternalField::scalarBeginsAt);

    T uSqr = advectionDiffusionDynamicsTemplates<T, Descriptor>::no_corr_bgk_collision(
        cell, rhoBar, jEq, omegaRef, sourceTerm);

    T kappa = (T)1 - this->getOmega() / omegaRef;

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        T correction = omegaRef * Descriptor<T>::t[iPop] * kappa * Descriptor<T>::invCs2
                       * intDot<T, Descriptor<T>::d>(Descriptor<T>::c[iPop], jNeq);
        cell[iPop] += correction;
    }

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void AdvectionDiffusionWithSourcePerkoDynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T thetaBar,
    BlockStatistics &statistics)
{
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_j(cell, j);
    Array<T, Descriptor<T>::d> jNeq(j - jEq);

    T sourceTerm = *cell.getExternal(Descriptor<T>::ExternalField::scalarBeginsAt);

    T uSqr = advectionDiffusionDynamicsTemplates<T, Descriptor>::no_corr_bgk_collision(
        cell, rhoBar, jEq, omegaRef, sourceTerm);

    T kappa = (T)1 - this->getOmega() / omegaRef;

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        T correction = omegaRef * Descriptor<T>::t[iPop] * kappa * Descriptor<T>::invCs2
                       * intDot(Descriptor<T>::c[iPop], jNeq);
        cell[iPop] += correction;
    }

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

/** \param jEq The parameter jEq is defined as jEq = j_advDiff = rho_advDiff*u_fluid
 */
template <typename T, template <typename U> class Descriptor>
T AdvectionDiffusionWithSourcePerkoDynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T jSqr, T thetaBar) const
{
    return advectionDiffusionDynamicsTemplates<T, Descriptor>::bgk_ma1_equilibrium(
        iPop, rhoBar, jEq);
}

/* *************** Class AdvectionDiffusionBGKdynamics *************** */

template <typename T, template <typename U> class Descriptor>
int AdvectionDiffusionBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, AdvectionDiffusionBGKdynamics<T, Descriptor> >(
        "AdvectionDiffusion_BGK");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
AdvectionDiffusionBGKdynamics<T, Descriptor>::AdvectionDiffusionBGKdynamics(T omega_) :
    AdvectionDiffusionDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
AdvectionDiffusionBGKdynamics<T, Descriptor>::AdvectionDiffusionBGKdynamics(
    HierarchicUnserializer &unserializer) :
    AdvectionDiffusionDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
AdvectionDiffusionBGKdynamics<T, Descriptor> *AdvectionDiffusionBGKdynamics<T, Descriptor>::clone()
    const
{
    return new AdvectionDiffusionBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int AdvectionDiffusionBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void AdvectionDiffusionBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> jEq;
    advectionDiffusionMomentTemplates<T, Descriptor>::get_rhoBar_jEq(cell, rhoBar, jEq);

    T uSqr = advectionDiffusionDynamicsTemplates<T, Descriptor>::no_corr_bgk_collision(
        cell, rhoBar, jEq, this->getOmega());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void AdvectionDiffusionBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T thetaBar,
    BlockStatistics &statistics)
{
    T uSqr = advectionDiffusionDynamicsTemplates<T, Descriptor>::no_corr_bgk_collision(
        cell, rhoBar, jEq, this->getOmega());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

/** \param jEq The parameter jEq is defined as jEq = j_advDiff = rho_advDiff*u_fluid
 */
template <typename T, template <typename U> class Descriptor>
T AdvectionDiffusionBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T jSqr, T thetaBar) const
{
    return advectionDiffusionDynamicsTemplates<T, Descriptor>::bgk_ma1_equilibrium(
        iPop, rhoBar, jEq);
}

// TODO implement decompose and recompose.
/* *************** Class CompleteAdvectionDiffusionBGKdynamics *************** */

template <typename T, template <typename U> class Descriptor>
int CompleteAdvectionDiffusionBGKdynamics<T, Descriptor>::id = meta::registerGeneralDynamics<
    T, Descriptor, CompleteAdvectionDiffusionBGKdynamics<T, Descriptor> >(
    "CompleteAdvectionDiffusion_BGK");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
CompleteAdvectionDiffusionBGKdynamics<T, Descriptor>::CompleteAdvectionDiffusionBGKdynamics(
    T omega_) :
    AdvectionDiffusionDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
CompleteAdvectionDiffusionBGKdynamics<T, Descriptor>::CompleteAdvectionDiffusionBGKdynamics(
    HierarchicUnserializer &unserializer) :
    AdvectionDiffusionDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
CompleteAdvectionDiffusionBGKdynamics<T, Descriptor>
    *CompleteAdvectionDiffusionBGKdynamics<T, Descriptor>::clone() const
{
    return new CompleteAdvectionDiffusionBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int CompleteAdvectionDiffusionBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
T CompleteAdvectionDiffusionBGKdynamics<T, Descriptor>::computeDensity(
    Cell<T, Descriptor> const &cell) const
{
    T rhoBar = *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt);
    return momentTemplates<T, Descriptor>::compute_rho(cell) * Descriptor<T>::invRho(rhoBar);
}

template <typename T, template <typename U> class Descriptor>
void CompleteAdvectionDiffusionBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoPhiBar = momentTemplates<T, Descriptor>::get_rhoBar(cell);
    T rhoBar = *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt);

    Array<T, Descriptor<T>::d> j;
    j.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));
    j *= Descriptor<T>::fullRho(rhoBar);
    //     Array<T,SymmetricTensor<T,Descriptor>::n> piNeq;
    //     piNeq.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::piNeqBeginsAt));

    T rhoPhi = Descriptor<T>::fullRho(rhoPhiBar);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T phi = rhoPhi * invRho;

    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(
        cell, rhoPhiBar, phi * j, this->getOmega());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, Descriptor<T>::rhoBar(rhoPhi), uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void CompleteAdvectionDiffusionBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &statistics)
{
    // TODO: IMPLEMENT
    PLB_ASSERT(false);
}

/** \param j The parameter j is defined as j = j_advDiff = phi*rho_fluid*u_fluid
 */
template <typename T, template <typename U> class Descriptor>
T CompleteAdvectionDiffusionBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

// j = sum_i c_i*f_i, piNeq is nothing, jSqr also
template <typename T, template <typename U> class Descriptor>
void CompleteAdvectionDiffusionBGKdynamics<T, Descriptor>::regularize(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar) const
{
    T rhoPhiBar = rhoBar;
    rhoBar = *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt);

    T phi = Descriptor<T>::fullRho(rhoPhiBar) * Descriptor<T>::invRho(rhoBar);

    Array<T, Descriptor<T>::d> jEq;
    jEq.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));

    jEq *= phi * Descriptor<T>::fullRho(rhoBar);

    Array<T, Descriptor<T>::d> jNeq = j - jEq;

    Array<T, SymmetricTensor<T, Descriptor>::n> pi;
    pi.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::piNeqBeginsAt));

    T omegaFluid = *cell.getExternal(Descriptor<T>::ExternalField::omegaBeginsAt);

    advectionDiffusionDynamicsTemplates<T, Descriptor>::complete_bgk_ma2_regularize(
        cell, rhoPhiBar, rhoBar, jEq, jNeq, pi, this->getOmega(), this->getOmega(), omegaFluid,
        omegaFluid);
}

// TODO implement decompose and recompose.
/* *************** Class CompleteRegularizedAdvectionDiffusionBGKdynamics *************** */

template <typename T, template <typename U> class Descriptor>
int CompleteRegularizedAdvectionDiffusionBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<
        T, Descriptor, CompleteRegularizedAdvectionDiffusionBGKdynamics<T, Descriptor> >(
        "CompleteRegularizedAdvectionDiffusion_BGK");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
CompleteRegularizedAdvectionDiffusionBGKdynamics<
    T, Descriptor>::CompleteRegularizedAdvectionDiffusionBGKdynamics(T omega_) :
    CompleteAdvectionDiffusionBGKdynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
CompleteRegularizedAdvectionDiffusionBGKdynamics<T, Descriptor>::
    CompleteRegularizedAdvectionDiffusionBGKdynamics(HierarchicUnserializer &unserializer) :
    CompleteAdvectionDiffusionBGKdynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
CompleteRegularizedAdvectionDiffusionBGKdynamics<T, Descriptor>
    *CompleteRegularizedAdvectionDiffusionBGKdynamics<T, Descriptor>::clone() const
{
    return new CompleteRegularizedAdvectionDiffusionBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int CompleteRegularizedAdvectionDiffusionBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
T CompleteRegularizedAdvectionDiffusionBGKdynamics<T, Descriptor>::computeDensity(
    Cell<T, Descriptor> const &cell) const
{
    T rhoBar = *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt);
    return momentTemplates<T, Descriptor>::compute_rho(cell) * Descriptor<T>::invRho(rhoBar);
}

template <typename T, template <typename U> class Descriptor>
void CompleteRegularizedAdvectionDiffusionBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoPhiBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoPhiBar, j);
    T rhoBar = *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt);

    //     Array<T,SymmetricTensor<T,Descriptor>::n> piNeq;
    //     piNeq.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::piNeqBeginsAt));

    T rhoPhi = Descriptor<T>::fullRho(rhoPhiBar);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T phi = rhoPhi * invRho;

    rhoBar = *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt);

    Array<T, Descriptor<T>::d> jEq;
    jEq.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));

    jEq *= phi * Descriptor<T>::fullRho(rhoBar);

    Array<T, Descriptor<T>::d> jNeq = j - jEq;

    Array<T, SymmetricTensor<T, Descriptor>::n> pi;
    pi.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::piNeqBeginsAt));

    T omegaFluid = *cell.getExternal(Descriptor<T>::ExternalField::omegaBeginsAt);

    T uSqr =
        advectionDiffusionDynamicsTemplates<T, Descriptor>::complete_bgk_ma2_regularized_collision(
            cell, rhoPhiBar, rhoBar, jEq, jNeq, pi, this->getOmega(), this->getOmega(), omegaFluid,
            omegaFluid);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, Descriptor<T>::rhoBar(rhoPhi), uSqr);
    }
}

// TODO implement decompose and recompose.
/* *************** Class CompleteRegularizedAdvectionDiffusionWithSourceBGKdynamics ***************
 */

template <typename T, template <typename U> class Descriptor>
int CompleteRegularizedAdvectionDiffusionWithSourceBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<
        T, Descriptor, CompleteRegularizedAdvectionDiffusionWithSourceBGKdynamics<T, Descriptor> >(
        "CompleteRegularizedAdvectionDiffusion_BGK");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
CompleteRegularizedAdvectionDiffusionWithSourceBGKdynamics<
    T, Descriptor>::CompleteRegularizedAdvectionDiffusionWithSourceBGKdynamics(T omega_) :
    CompleteAdvectionDiffusionBGKdynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
CompleteRegularizedAdvectionDiffusionWithSourceBGKdynamics<T, Descriptor>::
    CompleteRegularizedAdvectionDiffusionWithSourceBGKdynamics(
        HierarchicUnserializer &unserializer) :
    CompleteAdvectionDiffusionBGKdynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
CompleteRegularizedAdvectionDiffusionWithSourceBGKdynamics<T, Descriptor>
    *CompleteRegularizedAdvectionDiffusionWithSourceBGKdynamics<T, Descriptor>::clone() const
{
    return new CompleteRegularizedAdvectionDiffusionWithSourceBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int CompleteRegularizedAdvectionDiffusionWithSourceBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
T CompleteRegularizedAdvectionDiffusionWithSourceBGKdynamics<T, Descriptor>::computeDensity(
    Cell<T, Descriptor> const &cell) const
{
    T rhoBar = *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt);
    return momentTemplates<T, Descriptor>::compute_rho(cell) * Descriptor<T>::invRho(rhoBar);
}

template <typename T, template <typename U> class Descriptor>
void CompleteRegularizedAdvectionDiffusionWithSourceBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoPhiBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoPhiBar, j);
    T rhoBar = *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt);

    //     Array<T,SymmetricTensor<T,Descriptor>::n> piNeq;
    //     piNeq.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::piNeqBeginsAt));

    T rhoPhi = Descriptor<T>::fullRho(rhoPhiBar);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T phi = rhoPhi * invRho;

    rhoBar = *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt);

    Array<T, Descriptor<T>::d> jEq;
    jEq.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));

    jEq *= phi * Descriptor<T>::fullRho(rhoBar);

    Array<T, Descriptor<T>::d> jNeq = j - jEq;

    Array<T, SymmetricTensor<T, Descriptor>::n> pi;
    pi.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::piNeqBeginsAt));

    T omegaFluid = *cell.getExternal(Descriptor<T>::ExternalField::omegaBeginsAt);

    T uSqr =
        advectionDiffusionDynamicsTemplates<T, Descriptor>::complete_bgk_ma2_regularized_collision(
            cell, rhoPhiBar, rhoBar, jEq, jNeq, pi, this->getOmega(), this->getOmega(), omegaFluid,
            omegaFluid);

    T sourceTerm = *cell.getExternal(Descriptor<T>::ExternalField::scalarBeginsAt);
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] += Descriptor<T>::t[iPop] * sourceTerm;
    }

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, Descriptor<T>::rhoBar(rhoPhi), uSqr);
    }
}

// TODO implement decompose and recompose.
/* *************** Class CompleteAdvectionDiffusionTRTdynamics *************** */

template <typename T, template <typename U> class Descriptor>
int CompleteAdvectionDiffusionTRTdynamics<T, Descriptor>::id = meta::registerGeneralDynamics<
    T, Descriptor, CompleteAdvectionDiffusionTRTdynamics<T, Descriptor> >(
    "CompleteAdvectionDiffusion_TRT");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
CompleteAdvectionDiffusionTRTdynamics<T, Descriptor>::CompleteAdvectionDiffusionTRTdynamics(
    T omega_, T psi_) :
    AdvectionDiffusionDynamics<T, Descriptor>(omega_), psi(psi_)
{ }

template <typename T, template <typename U> class Descriptor>
CompleteAdvectionDiffusionTRTdynamics<T, Descriptor>::CompleteAdvectionDiffusionTRTdynamics(
    HierarchicUnserializer &unserializer) :
    AdvectionDiffusionDynamics<T, Descriptor>(T()), psi(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
CompleteAdvectionDiffusionTRTdynamics<T, Descriptor>::CompleteAdvectionDiffusionTRTdynamics(
    T omega_) :
    AdvectionDiffusionDynamics<T, Descriptor>(omega_)
{
    psi = advectionDiffusionDynamicsTemplatesImpl<
        T, typename Descriptor<T>::BaseDescriptor>::computePsiComplete(omega_);
}

template <typename T, template <typename U> class Descriptor>
CompleteAdvectionDiffusionTRTdynamics<T, Descriptor>
    *CompleteAdvectionDiffusionTRTdynamics<T, Descriptor>::clone() const
{
    return new CompleteAdvectionDiffusionTRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int CompleteAdvectionDiffusionTRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void CompleteAdvectionDiffusionTRTdynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    AdvectionDiffusionDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(psi);
}

template <typename T, template <typename U> class Descriptor>
void CompleteAdvectionDiffusionTRTdynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    AdvectionDiffusionDynamics<T, Descriptor>::unserialize(unserializer);
    psi = unserializer.readValue<T>();
}

template <typename T, template <typename U> class Descriptor>
T CompleteAdvectionDiffusionTRTdynamics<T, Descriptor>::computeDensity(
    Cell<T, Descriptor> const &cell) const
{
    T rhoBar = *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt);
    return momentTemplates<T, Descriptor>::compute_rho(cell) * Descriptor<T>::invRho(rhoBar);
}

template <typename T, template <typename U> class Descriptor>
void CompleteAdvectionDiffusionTRTdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoPhiBar = momentTemplates<T, Descriptor>::get_rhoBar(cell);
    T rhoBar = *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt);

    Array<T, Descriptor<T>::d> j;
    j.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));

    T rhoPhi = Descriptor<T>::fullRho(rhoPhiBar);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T phi = rhoPhi * invRho;

    T uSqr =
        advectionDiffusionDynamicsTemplates<T, Descriptor>::complete_mrt_ma2_ext_rhoBar_j_collision(
            cell, rhoPhiBar, phi * j, this->getOmega(), psi);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, Descriptor<T>::rhoBar(rhoPhi), uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void CompleteAdvectionDiffusionTRTdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &statistics)
{
    T uSqr =
        advectionDiffusionDynamicsTemplates<T, Descriptor>::complete_mrt_ma2_ext_rhoBar_j_collision(
            cell, rhoBar, j, this->getOmega(), psi);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

/** \param j The parameter j is defined as j = j_advDiff = phi*rho_fluid*u_fluid
 */
template <typename T, template <typename U> class Descriptor>
T CompleteAdvectionDiffusionTRTdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_equilibrium(
        iPop, rhoBar, invRho, j, jSqr);
}

// j = sum_i c_i*f_i, piNeq is nothing, jSqr also
template <typename T, template <typename U> class Descriptor>
void CompleteAdvectionDiffusionTRTdynamics<T, Descriptor>::regularize(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar) const
{
    T rhoPhiBar = rhoBar;
    rhoBar = *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt);

    T phi = Descriptor<T>::fullRho(rhoPhiBar) * Descriptor<T>::invRho(rhoBar);

    Array<T, Descriptor<T>::d> jEq;
    jEq.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));

    jEq *= phi;

    Array<T, Descriptor<T>::d> jNeq = j - jEq;

    Array<T, SymmetricTensor<T, Descriptor>::n> pi;
    pi.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::piNeqBeginsAt));

    T omegaFluid = *cell.getExternal(Descriptor<T>::ExternalField::omegaBeginsAt);

    advectionDiffusionDynamicsTemplates<T, Descriptor>::complete_bgk_ma2_regularize(
        cell, rhoPhiBar, rhoBar, jEq, jNeq, pi, this->getOmega(), psi, omegaFluid, omegaFluid);
}

template <typename T, template <typename U> class Descriptor>
T CompleteAdvectionDiffusionTRTdynamics<T, Descriptor>::getParameter(plint whichParameter) const
{
    switch (whichParameter) {
    case dynamicParams::omega_shear:
        return this->getOmega();
    case dynamicParams::psi:
        return getPsi();
    };
    return 0.;
}

template <typename T, template <typename U> class Descriptor>
void CompleteAdvectionDiffusionTRTdynamics<T, Descriptor>::setParameter(
    plint whichParameter, T value)
{
    switch (whichParameter) {
    case dynamicParams::omega_shear:
        this->setOmega(value);
    case dynamicParams::psi:
        setPsi(value);
    };
}

template <typename T, template <typename U> class Descriptor>
T CompleteAdvectionDiffusionTRTdynamics<T, Descriptor>::getPsi() const
{
    return psi;
}

template <typename T, template <typename U> class Descriptor>
void CompleteAdvectionDiffusionTRTdynamics<T, Descriptor>::setPsi(T psi_)
{
    psi = psi_;
}

/* *************** Class AdvectionDiffusionWithSourceBGKdynamics *************** */

template <typename T, template <typename U> class Descriptor>
int AdvectionDiffusionWithSourceBGKdynamics<T, Descriptor>::id = meta::registerGeneralDynamics<
    T, Descriptor, AdvectionDiffusionWithSourceBGKdynamics<T, Descriptor> >(
    "AdvectionDiffusionWithSource_BGK");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
AdvectionDiffusionWithSourceBGKdynamics<T, Descriptor>::AdvectionDiffusionWithSourceBGKdynamics(
    T omega_) :
    AdvectionDiffusionDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
AdvectionDiffusionWithSourceBGKdynamics<T, Descriptor>::AdvectionDiffusionWithSourceBGKdynamics(
    HierarchicUnserializer &unserializer) :
    AdvectionDiffusionDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
AdvectionDiffusionWithSourceBGKdynamics<T, Descriptor>
    *AdvectionDiffusionWithSourceBGKdynamics<T, Descriptor>::clone() const
{
    return new AdvectionDiffusionWithSourceBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int AdvectionDiffusionWithSourceBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void AdvectionDiffusionWithSourceBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> jEq;
    advectionDiffusionMomentTemplates<T, Descriptor>::get_rhoBar_jEq(cell, rhoBar, jEq);

    T sourceTerm = *cell.getExternal(Descriptor<T>::ExternalField::scalarBeginsAt);

    T uSqr = advectionDiffusionDynamicsTemplates<T, Descriptor>::no_corr_bgk_collision(
        cell, rhoBar, jEq, this->getOmega(), sourceTerm);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

/** \param jEq The parameter jEq is defined as jEq = j_advDiff = rho_advDiff*u_fluid
 */
template <typename T, template <typename U> class Descriptor>
T AdvectionDiffusionWithSourceBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T jSqr, T thetaBar) const
{
    return advectionDiffusionDynamicsTemplates<T, Descriptor>::bgk_ma1_equilibrium(
        iPop, rhoBar, jEq);
}

}  // namespace plb

#endif  // ADVECTION_DIFFUSION_DYNAMICS_HH
