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
#ifndef ISO_THERMAL_DYNAMICS_HH
#define ISO_THERMAL_DYNAMICS_HH

#include <algorithm>
#include <cstdlib>
#include <limits>

#include "basicDynamics/isoThermalDynamics.h"
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

/* *************** Class IsoThermalBulkDynamics ************************************ */

template <typename T, template <typename U> class Descriptor>
IsoThermalBulkDynamics<T, Descriptor>::IsoThermalBulkDynamics(T omega_) :
    BasicBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
void IsoThermalBulkDynamics<T, Descriptor>::regularize(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar) const
{
    typedef Descriptor<T> L;
    cell[0] = this->computeEquilibrium(0, rhoBar, j, jSqr)
              + offEquilibriumTemplates<T, Descriptor>::fromPiToFneq(0, PiNeq);
    for (plint iPop = 1; iPop <= L::q / 2; ++iPop) {
        cell[iPop] = this->computeEquilibrium(iPop, rhoBar, j, jSqr);
        cell[iPop + L::q / 2] = this->computeEquilibrium(iPop + L::q / 2, rhoBar, j, jSqr);
        T fNeq = offEquilibriumTemplates<T, Descriptor>::fromPiToFneq(iPop, PiNeq);
        cell[iPop] += fNeq;
        cell[iPop + L::q / 2] += fNeq;
    }
}

template <typename T, template <typename U> class Descriptor>
T IsoThermalBulkDynamics<T, Descriptor>::computeTemperature(Cell<T, Descriptor> const &cell) const
{
    return (T)1;
}

template <typename T, template <typename U> class Descriptor>
void IsoThermalBulkDynamics<T, Descriptor>::computePiNeq(
    Cell<T, Descriptor> const &cell, Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    momentTemplates<T, Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);
}

template <typename T, template <typename U> class Descriptor>
void IsoThermalBulkDynamics<T, Descriptor>::computeShearStress(
    Cell<T, Descriptor> const &cell, Array<T, SymmetricTensor<T, Descriptor>::n> &stress) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    momentTemplates<T, Descriptor>::compute_PiNeq(cell, rhoBar, j, stress);
    T omega = cell.getDynamics().getOmega();
    stress *= ((T)0.5 * omega - (T)1);
}

template <typename T, template <typename U> class Descriptor>
void IsoThermalBulkDynamics<T, Descriptor>::computeHeatFlux(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &q) const
{
    q.resetToZero();
}

template <typename T, template <typename U> class Descriptor>
T IsoThermalBulkDynamics<T, Descriptor>::computeEbar(Cell<T, Descriptor> const &cell) const
{
    return T();
}

template <typename T, template <typename U> class Descriptor>
plint IsoThermalBulkDynamics<T, Descriptor>::numDecomposedVariables(plint order) const
{
    // Start with the decomposed version of the populations.
    plint numVariables =
        // Order 0: density + velocity + fNeq
        (order == 0) ? (1 + Descriptor<T>::d + Descriptor<T>::q)
                     // Order >=1: density + velocity + PiNeq
                     : (1 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n);

    // Add the variables in the external scalars.
    numVariables += Descriptor<T>::ExternalField::numScalars;
    return numVariables;
}

template <typename T, template <typename U> class Descriptor>
void IsoThermalBulkDynamics<T, Descriptor>::decompose(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData, plint order) const
{
    rawData.resize(numDecomposedVariables(order));

    if (order == 0) {
        decomposeOrder0(cell, rawData);
    } else {
        decomposeOrder1(cell, rawData);
    }
}

template <typename T, template <typename U> class Descriptor>
void IsoThermalBulkDynamics<T, Descriptor>::recompose(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData, plint order) const
{
    PLB_PRECONDITION((plint)rawData.size() == numDecomposedVariables(order));

    if (order == 0) {
        recomposeOrder0(cell, rawData);
    } else {
        recomposeOrder1(cell, rawData);
    }
}

template <typename T, template <typename U> class Descriptor>
void IsoThermalBulkDynamics<T, Descriptor>::rescale(
    std::vector<T> &rawData, T xDxInv, T xDt, plint order) const
{
    PLB_PRECONDITION((plint)rawData.size() == numDecomposedVariables(order));

    if (order == 0) {
        rescaleOrder0(rawData, xDxInv, xDt);
    } else {
        rescaleOrder1(rawData, xDxInv, xDt);
    }
}

template <typename T, template <typename U> class Descriptor>
void IsoThermalBulkDynamics<T, Descriptor>::decomposeOrder0(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
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
}

template <typename T, template <typename U> class Descriptor>
void IsoThermalBulkDynamics<T, Descriptor>::decomposeOrder1(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);

    rawData[0] = rhoBar;
    j.to_cArray(&rawData[1]);
    PiNeq.to_cArray(&rawData[1 + Descriptor<T>::d]);

    int offset = 1 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        rawData[offset + iExt] = *cell.getExternal(iExt);
    }
}

template <typename T, template <typename U> class Descriptor>
void IsoThermalBulkDynamics<T, Descriptor>::recomposeOrder0(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    T rhoBar = rawData[0];
    Array<T, Descriptor<T>::d> j;
    j.from_cArray(&rawData[1]);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] =
            this->computeEquilibrium(iPop, rhoBar, j, jSqr) + rawData[1 + Descriptor<T>::d + iPop];
    }

    int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *cell.getExternal(iExt) = rawData[offset + iExt];
    }
}

template <typename T, template <typename U> class Descriptor>
void IsoThermalBulkDynamics<T, Descriptor>::recomposeOrder1(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;

    rhoBar = rawData[0];
    j.from_cArray(&rawData[1]);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    PiNeq.from_cArray(&rawData[1 + Descriptor<T>::d]);

    this->regularize(cell, rhoBar, j, jSqr, PiNeq);

    int offset = 1 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *cell.getExternal(iExt) = rawData[offset + iExt];
    }
}

template <typename T, template <typename U> class Descriptor>
void IsoThermalBulkDynamics<T, Descriptor>::rescaleOrder0(
    std::vector<T> &rawData, T xDxInv, T xDt) const
{
    // Don't change rho (rawData[0]), because it is invariant

    // Change velocity, according to its units dx/dt
    T velScale = xDt * xDxInv;
    for (plint iVel = 0; iVel < Descriptor<T>::d; ++iVel) {
        rawData[1 + iVel] *= velScale;
    }

    // Change off-equilibrium, according to its units 1/dt
    T fNeqScale = xDt;
    for (plint iFneq = 0; iFneq < Descriptor<T>::q; ++iFneq) {
        rawData[1 + Descriptor<T>::d + iFneq] *= fNeqScale;
    }

    // Don't change external fields; their scaling must be taken care of
    //   in specialized versions of this class.
}

template <typename T, template <typename U> class Descriptor>
void IsoThermalBulkDynamics<T, Descriptor>::rescaleOrder1(
    std::vector<T> &rawData, T xDxInv, T xDt) const
{
    // Don't change rho (rawData[0]), because it is invariant

    // Change velocity, according to its units dx/dt
    T velScale = xDt * xDxInv;
    for (plint iVel = 0; iVel < Descriptor<T>::d; ++iVel) {
        rawData[1 + iVel] *= velScale;
    }

    // Change off-equilibrium stress, according to its units 1/dt
    T PiNeqScale = xDt;
    for (plint iPi = 0; iPi < SymmetricTensor<T, Descriptor>::n; ++iPi) {
        rawData[1 + Descriptor<T>::d + iPi] *= PiNeqScale;
    }

    // Don't change external fields; their scaling must be taken care of
    //   in specialized versions of this class.
}

/* *************** Class BGKdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
int BGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, BGKdynamics<T, Descriptor> >("BGK");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
BGKdynamics<T, Descriptor>::BGKdynamics(T omega_) : IsoThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
BGKdynamics<T, Descriptor>::BGKdynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
BGKdynamics<T, Descriptor> *BGKdynamics<T, Descriptor>::clone() const
{
    return new BGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int BGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void BGKdynamics<T, Descriptor>::collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, j, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void BGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, j, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T BGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
void BGKdynamics<T, Descriptor>::decomposeOrder0(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    rawData[0] = rhoBar;
    j.to_cArray(&rawData[1]);

    Array<T, Descriptor<T>::q> fEq;
    dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibria(
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
void BGKdynamics<T, Descriptor>::recomposeOrder0(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    T rhoBar = rawData[0];
    Array<T, Descriptor<T>::d> j;
    j.from_cArray(&rawData[1]);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    Array<T, Descriptor<T>::q> fEq;
    dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibria(
        rhoBar, Descriptor<T>::invRho(rhoBar), j, jSqr, fEq);

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] = fEq[iPop] + rawData[1 + Descriptor<T>::d + iPop];
    }

    int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *cell.getExternal(iExt) = rawData[offset + iExt];
    }
}

/* *************** Class StochasticDynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
int StochasticDynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, StochasticDynamics<T, Descriptor> >("Stochastic");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
StochasticDynamics<T, Descriptor>::StochasticDynamics(T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
StochasticDynamics<T, Descriptor>::StochasticDynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
StochasticDynamics<T, Descriptor> *StochasticDynamics<T, Descriptor>::clone() const
{
    return new StochasticDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int StochasticDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void StochasticDynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T uSqr = 0.;
    if (rand() < RAND_MAX / 4) {
        uSqr =
            dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, j, this->getOmega());
    } else {
        uSqr = normSqr(j * Descriptor<T>::invRho(rhoBar));
    }
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void StochasticDynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T uSqr = 0.;
    if (rand() < RAND_MAX / 4) {
        uSqr =
            dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, j, this->getOmega());
    } else {
        uSqr = normSqr(j * Descriptor<T>::invRho(rhoBar));
    }
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T StochasticDynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
void StochasticDynamics<T, Descriptor>::decomposeOrder0(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    rawData[0] = rhoBar;
    j.to_cArray(&rawData[1]);

    Array<T, Descriptor<T>::q> fEq;
    dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibria(
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
void StochasticDynamics<T, Descriptor>::recomposeOrder0(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    T rhoBar = rawData[0];
    Array<T, Descriptor<T>::d> j;
    j.from_cArray(&rawData[1]);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    Array<T, Descriptor<T>::q> fEq;
    dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibria(
        rhoBar, Descriptor<T>::invRho(rhoBar), j, jSqr, fEq);

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] = fEq[iPop] + rawData[1 + Descriptor<T>::d + iPop];
    }

    int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *cell.getExternal(iExt) = rawData[offset + iExt];
    }
}

/* *************** Class CompleteBGKdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
int CompleteBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, CompleteBGKdynamics<T, Descriptor> >(
        "Complete_BGK");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
CompleteBGKdynamics<T, Descriptor>::CompleteBGKdynamics(T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
CompleteBGKdynamics<T, Descriptor>::CompleteBGKdynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
CompleteBGKdynamics<T, Descriptor> *CompleteBGKdynamics<T, Descriptor>::clone() const
{
    return new CompleteBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int CompleteBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void CompleteBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T uSqr = dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_collision(
        cell, rhoBar, Descriptor<T>::invRho(rhoBar), j, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void CompleteBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T uSqr = dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_collision(
        cell, rhoBar, Descriptor<T>::invRho(rhoBar), j, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void CompleteBGKdynamics<T, Descriptor>::computeEquilibria(
    Array<T, Descriptor<T>::q> &fEq, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_equilibria(rhoBar, invRho, j, jSqr, fEq);
}

template <typename T, template <typename U> class Descriptor>
T CompleteBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_equilibrium(
        iPop, rhoBar, invRho, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
void CompleteBGKdynamics<T, Descriptor>::regularize(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_regularize(
        cell, rhoBar, invRho, j, jSqr, PiNeq);
}

template <typename T, template <typename U> class Descriptor>
void CompleteBGKdynamics<T, Descriptor>::decomposeOrder0(
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
void CompleteBGKdynamics<T, Descriptor>::recomposeOrder0(
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

/* *************** Class CompleteRegularizedBGKdynamics
 * *********************************************** */

template <typename T, template <typename U> class Descriptor>
int CompleteRegularizedBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, CompleteRegularizedBGKdynamics<T, Descriptor> >(
        "Complete_Regularized_BGK");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
CompleteRegularizedBGKdynamics<T, Descriptor>::CompleteRegularizedBGKdynamics(T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
CompleteRegularizedBGKdynamics<T, Descriptor>::CompleteRegularizedBGKdynamics(
    HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
CompleteRegularizedBGKdynamics<T, Descriptor>
    *CompleteRegularizedBGKdynamics<T, Descriptor>::clone() const
{
    return new CompleteRegularizedBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int CompleteRegularizedBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void CompleteRegularizedBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> piNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, piNeq);
    T uSqr = dynamicsTemplates<T, Descriptor>::complete_regularized_bgk_ma2_collision(
        cell, rhoBar, j, piNeq, this->getOmega());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void CompleteRegularizedBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T rhoBarLb;
    Array<T, Descriptor<T>::d> jLb;
    Array<T, SymmetricTensor<T, Descriptor>::n> piNeqLb;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBarLb, jLb, piNeqLb);
    T jSqrLb = VectorTemplate<T, Descriptor>::normSqr(jLb);
    regularize(cell, rhoBarLb, jLb, jSqrLb, piNeqLb);

    T uSqr = dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_collision(
        cell, rhoBar, Descriptor<T>::invRho(rhoBar), j, this->getOmega());

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void CompleteRegularizedBGKdynamics<T, Descriptor>::computeEquilibria(
    Array<T, Descriptor<T>::q> &fEq, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_equilibria(rhoBar, invRho, j, jSqr, fEq);
}

template <typename T, template <typename U> class Descriptor>
T CompleteRegularizedBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_equilibrium(
        iPop, rhoBar, invRho, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
void CompleteRegularizedBGKdynamics<T, Descriptor>::regularize(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_regularize(
        cell, rhoBar, invRho, j, jSqr, PiNeq, this->getOmega(), this->getOmega());
}

template <typename T, template <typename U> class Descriptor>
void CompleteRegularizedBGKdynamics<T, Descriptor>::decomposeOrder0(
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
void CompleteRegularizedBGKdynamics<T, Descriptor>::recomposeOrder0(
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

/* *************** Class CompleteTRTdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
int CompleteTRTdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, CompleteTRTdynamics<T, Descriptor> >(
        "Complete_TRT");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
CompleteTRTdynamics<T, Descriptor>::CompleteTRTdynamics(T omega_, T psi_, int order_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_), psi(psi_), order(order_)
{
    PLB_ASSERT(order >= 2 && "Order must bne greater than 2.");
}

template <typename T, template <typename U> class Descriptor>
CompleteTRTdynamics<T, Descriptor>::CompleteTRTdynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T()), psi(T()), order(0)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
CompleteTRTdynamics<T, Descriptor>::CompleteTRTdynamics(T omega_, int order_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{
    psi = dynamicsTemplates<T, Descriptor>::computePsiComplete(this->getOmega());
    order = order_;
    PLB_ASSERT(order >= 2 && "Order must bne greater than 2.");
}

template <typename T, template <typename U> class Descriptor>
void CompleteTRTdynamics<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    IsoThermalBulkDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(psi);
    serializer.addValue(order);
}

template <typename T, template <typename U> class Descriptor>
void CompleteTRTdynamics<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    IsoThermalBulkDynamics<T, Descriptor>::unserialize(unserializer);
    psi = unserializer.readValue<T>();
    order = unserializer.readValue<int>();
}

template <typename T, template <typename U> class Descriptor>
CompleteTRTdynamics<T, Descriptor> *CompleteTRTdynamics<T, Descriptor>::clone() const
{
    return new CompleteTRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int CompleteTRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void CompleteTRTdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T uSqr = dynamicsTemplates<T, Descriptor>::complete_mrt_ma2_collision(
        cell, order, this->getOmega(), psi);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, momentTemplates<T, Descriptor>::get_rhoBar(cell), uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void CompleteTRTdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T uSqr = dynamicsTemplates<T, Descriptor>::complete_mrt_ma2_ext_rhoBar_j_collision(
        cell, rhoBar, j, order, this->getOmega(), psi);

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void CompleteTRTdynamics<T, Descriptor>::computeEquilibria(
    Array<T, Descriptor<T>::q> &fEq, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_equilibria(rhoBar, invRho, j, jSqr, fEq);
}

template <typename T, template <typename U> class Descriptor>
T CompleteTRTdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_equilibrium(
        iPop, rhoBar, invRho, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
void CompleteTRTdynamics<T, Descriptor>::regularize(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_regularize(
        cell, rhoBar, invRho, j, jSqr, PiNeq, this->getOmega(), psi);
}

template <typename T, template <typename U> class Descriptor>
Array<T, Descriptor<T>::q> CompleteTRTdynamics<T, Descriptor>::getRelaxationFrequencies() const
{
    Array<T, Descriptor<T>::q> frequencies;
    for (plint i = 0; i <= Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n; ++i) {
        frequencies[i] = this->getOmega();
    }
    for (plint i = Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n + 1; i < Descriptor<T>::q;
         ++i) {
        frequencies[i] = psi;
    }
    return frequencies;
}

template <typename T, template <typename U> class Descriptor>
void CompleteTRTdynamics<T, Descriptor>::setRelaxationFrequencies(
    Array<T, Descriptor<T>::q> const &frequencies)
{
    this->setOmega(frequencies[0]);
    setPsi(frequencies[1 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n]);
}

template <typename T, template <typename U> class Descriptor>
T CompleteTRTdynamics<T, Descriptor>::getParameter(plint whichParameter) const
{
    switch (whichParameter) {
    case dynamicParams::omega_shear:
        return this->getOmega();
    case dynamicParams::psi:
        return this->getPsi();
    };
    return 0.;
}

template <typename T, template <typename U> class Descriptor>
void CompleteTRTdynamics<T, Descriptor>::setParameter(plint whichParameter, T value)
{
    switch (whichParameter) {
    case dynamicParams::omega_shear:
        this->setOmega(value);
    case dynamicParams::psi:
        setPsi(value);
    };
}

template <typename T, template <typename U> class Descriptor>
T CompleteTRTdynamics<T, Descriptor>::getPsi() const
{
    return psi;
}

template <typename T, template <typename U> class Descriptor>
void CompleteTRTdynamics<T, Descriptor>::setPsi(T psi_)
{
    psi = psi_;
}

template <typename T, template <typename U> class Descriptor>
void CompleteTRTdynamics<T, Descriptor>::decomposeOrder0(
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
void CompleteTRTdynamics<T, Descriptor>::recomposeOrder0(
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

/* *************** Class CompleteRegularizedTRTdynamics
 * *********************************************** */

template <typename T, template <typename U> class Descriptor>
int CompleteRegularizedTRTdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, CompleteRegularizedTRTdynamics<T, Descriptor> >(
        "Complete_Regularized_TRT");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
CompleteRegularizedTRTdynamics<T, Descriptor>::CompleteRegularizedTRTdynamics(
    T omega_, T psi_, int order_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_), psi(psi_), order(order_)
{ }

template <typename T, template <typename U> class Descriptor>
CompleteRegularizedTRTdynamics<T, Descriptor>::CompleteRegularizedTRTdynamics(
    HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T()), psi(T()), order(0)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
CompleteRegularizedTRTdynamics<T, Descriptor>::CompleteRegularizedTRTdynamics(
    T omega_, int order_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{
    psi = dynamicsTemplates<T, Descriptor>::computePsiComplete(this->getOmega());
    order = order_;
    PLB_ASSERT(order >= 2 && "Order must be greater than 2.");
}

template <typename T, template <typename U> class Descriptor>
void CompleteRegularizedTRTdynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    IsoThermalBulkDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(psi);
    serializer.addValue(order);
}

template <typename T, template <typename U> class Descriptor>
void CompleteRegularizedTRTdynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    IsoThermalBulkDynamics<T, Descriptor>::unserialize(unserializer);
    psi = unserializer.readValue<T>();
    order = unserializer.readValue<int>();
}

template <typename T, template <typename U> class Descriptor>
CompleteRegularizedTRTdynamics<T, Descriptor>
    *CompleteRegularizedTRTdynamics<T, Descriptor>::clone() const
{
    return new CompleteRegularizedTRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int CompleteRegularizedTRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void CompleteRegularizedTRTdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> piNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, piNeq);
    T uSqr = dynamicsTemplates<T, Descriptor>::complete_regularized_mrt_ma2_collision(
        cell, rhoBar, j, piNeq, order, this->getOmega(), psi);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void CompleteRegularizedTRTdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T rhoBarLb;
    Array<T, Descriptor<T>::d> jLb;
    Array<T, SymmetricTensor<T, Descriptor>::n> piNeqLb;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBarLb, jLb, piNeqLb);
    T jSqrLb = VectorTemplate<T, Descriptor>::normSqr(jLb);
    regularize(cell, rhoBarLb, jLb, jSqrLb, piNeqLb);

    T uSqr = dynamicsTemplates<T, Descriptor>::complete_mrt_ma2_ext_rhoBar_j_collision(
        cell, rhoBar, j, order, this->getOmega(), psi);

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void CompleteRegularizedTRTdynamics<T, Descriptor>::computeEquilibria(
    Array<T, Descriptor<T>::q> &fEq, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_equilibria(rhoBar, invRho, j, jSqr, fEq);
}

template <typename T, template <typename U> class Descriptor>
T CompleteRegularizedTRTdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_equilibrium(
        iPop, rhoBar, invRho, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
void CompleteRegularizedTRTdynamics<T, Descriptor>::regularize(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_regularize(
        cell, rhoBar, invRho, j, jSqr, PiNeq, this->getOmega(), psi);
}

template <typename T, template <typename U> class Descriptor>
Array<T, Descriptor<T>::q> CompleteRegularizedTRTdynamics<T, Descriptor>::getRelaxationFrequencies()
    const
{
    Array<T, Descriptor<T>::q> frequencies;
    for (plint i = 0; i <= Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n; ++i) {
        frequencies[i] = this->getOmega();
    }
    for (plint i = Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n + 1; i < Descriptor<T>::q;
         ++i) {
        frequencies[i] = psi;
    }
    return frequencies;
}

template <typename T, template <typename U> class Descriptor>
void CompleteRegularizedTRTdynamics<T, Descriptor>::setRelaxationFrequencies(
    Array<T, Descriptor<T>::q> const &frequencies)
{
    this->setOmega(frequencies[0]);
    setPsi(frequencies[1 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n]);
}

template <typename T, template <typename U> class Descriptor>
T CompleteRegularizedTRTdynamics<T, Descriptor>::getParameter(plint whichParameter) const
{
    switch (whichParameter) {
    case dynamicParams::omega_shear:
        return this->getOmega();
    case dynamicParams::psi:
        return this->getPsi();
    };
    return 0.;
}

template <typename T, template <typename U> class Descriptor>
void CompleteRegularizedTRTdynamics<T, Descriptor>::setParameter(plint whichParameter, T value)
{
    switch (whichParameter) {
    case dynamicParams::omega_shear:
        this->setOmega(value);
        break;
    case dynamicParams::psi:
        setPsi(value);
    };
}

template <typename T, template <typename U> class Descriptor>
T CompleteRegularizedTRTdynamics<T, Descriptor>::getPsi() const
{
    return psi;
}

template <typename T, template <typename U> class Descriptor>
void CompleteRegularizedTRTdynamics<T, Descriptor>::setPsi(T psi_)
{
    psi = psi_;
}

template <typename T, template <typename U> class Descriptor>
void CompleteRegularizedTRTdynamics<T, Descriptor>::decomposeOrder0(
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
void CompleteRegularizedTRTdynamics<T, Descriptor>::recomposeOrder0(
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

/* *************** Class TruncatedTRTdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
int TruncatedTRTdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, TruncatedTRTdynamics<T, Descriptor> >(
        "Truncated_TRT");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
TruncatedTRTdynamics<T, Descriptor>::TruncatedTRTdynamics(T omega_, T psi_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_), psi(psi_)
{ }

template <typename T, template <typename U> class Descriptor>
TruncatedTRTdynamics<T, Descriptor>::TruncatedTRTdynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T()), psi(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
TruncatedTRTdynamics<T, Descriptor>::TruncatedTRTdynamics(T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{
    psi = dynamicsTemplates<T, Descriptor>::computePsiTruncated(this->getOmega());
}

template <typename T, template <typename U> class Descriptor>
void TruncatedTRTdynamics<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    IsoThermalBulkDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(psi);
}

template <typename T, template <typename U> class Descriptor>
void TruncatedTRTdynamics<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    IsoThermalBulkDynamics<T, Descriptor>::unserialize(unserializer);
    psi = unserializer.readValue<T>();
}

template <typename T, template <typename U> class Descriptor>
TruncatedTRTdynamics<T, Descriptor> *TruncatedTRTdynamics<T, Descriptor>::clone() const
{
    return new TruncatedTRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int TruncatedTRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void TruncatedTRTdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T uSqr =
        dynamicsTemplates<T, Descriptor>::truncated_mrt_ma2_collision(cell, this->getOmega(), psi);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, momentTemplates<T, Descriptor>::get_rhoBar(cell), uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void TruncatedTRTdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T uSqr = dynamicsTemplates<T, Descriptor>::truncated_mrt_ma2_ext_rhoBar_j_collision(
        cell, rhoBar, j, this->getOmega(), psi);

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void TruncatedTRTdynamics<T, Descriptor>::computeEquilibria(
    Array<T, Descriptor<T>::q> &fEq, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibria(rhoBar, invRho, j, jSqr, fEq);
}

template <typename T, template <typename U> class Descriptor>
T TruncatedTRTdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
Array<T, Descriptor<T>::q> TruncatedTRTdynamics<T, Descriptor>::getRelaxationFrequencies() const
{
    Array<T, Descriptor<T>::q> frequencies;
    for (plint i = 0; i <= Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n; ++i) {
        frequencies[i] = this->getOmega();
    }
    for (plint i = Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n + 1; i < Descriptor<T>::q;
         ++i) {
        frequencies[i] = psi;
    }
    return frequencies;
}

template <typename T, template <typename U> class Descriptor>
void TruncatedTRTdynamics<T, Descriptor>::setRelaxationFrequencies(
    Array<T, Descriptor<T>::q> const &frequencies)
{
    this->setOmega(frequencies[0]);
    setPsi(frequencies[1 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n]);
}

template <typename T, template <typename U> class Descriptor>
T TruncatedTRTdynamics<T, Descriptor>::getParameter(plint whichParameter) const
{
    switch (whichParameter) {
    case dynamicParams::omega_shear:
        return this->getOmega();
    case dynamicParams::psi:
        return this->getPsi();
    };
    return 0.;
}

template <typename T, template <typename U> class Descriptor>
void TruncatedTRTdynamics<T, Descriptor>::setParameter(plint whichParameter, T value)
{
    switch (whichParameter) {
    case dynamicParams::omega_shear:
        this->setOmega(value);
        break;
    case dynamicParams::psi:
        setPsi(value);
    };
}

template <typename T, template <typename U> class Descriptor>
T TruncatedTRTdynamics<T, Descriptor>::getPsi() const
{
    return psi;
}

template <typename T, template <typename U> class Descriptor>
void TruncatedTRTdynamics<T, Descriptor>::setPsi(T psi_)
{
    psi = psi_;
}

template <typename T, template <typename U> class Descriptor>
void TruncatedTRTdynamics<T, Descriptor>::decomposeOrder0(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    rawData[0] = rhoBar;
    j.to_cArray(&rawData[1]);

    Array<T, Descriptor<T>::q> fEq;
    dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibria(
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
void TruncatedTRTdynamics<T, Descriptor>::recomposeOrder0(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    T rhoBar = rawData[0];
    Array<T, Descriptor<T>::d> j;
    j.from_cArray(&rawData[1]);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    Array<T, Descriptor<T>::q> fEq;
    dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibria(
        rhoBar, Descriptor<T>::invRho(rhoBar), j, jSqr, fEq);

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] = fEq[iPop] + rawData[1 + Descriptor<T>::d + iPop];
    }

    int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *cell.getExternal(iExt) = rawData[offset + iExt];
    }
}

/* *************** Class StoreRhoBarJBGKdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
int StoreRhoBarJBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, StoreRhoBarJBGKdynamics<T, Descriptor> >(
        "StoreRhoBarJBGK");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
StoreRhoBarJBGKdynamics<T, Descriptor>::StoreRhoBarJBGKdynamics(T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
StoreRhoBarJBGKdynamics<T, Descriptor>::StoreRhoBarJBGKdynamics(
    HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
StoreRhoBarJBGKdynamics<T, Descriptor> *StoreRhoBarJBGKdynamics<T, Descriptor>::clone() const
{
    return new StoreRhoBarJBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int StoreRhoBarJBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void StoreRhoBarJBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, j, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
    *(cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt)) = rhoBar;
    j.to_cArray(cell.getExternal(Descriptor<T>::ExternalField::jBeginsAt));
}

template <typename T, template <typename U> class Descriptor>
void StoreRhoBarJBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, j, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
    *(cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt)) = rhoBar;
    j.to_cArray(cell.getExternal(Descriptor<T>::ExternalField::jBeginsAt));
}

template <typename T, template <typename U> class Descriptor>
T StoreRhoBarJBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class ExternalMomentBGKdynamics ********************************** */

template <typename T, template <typename U> class Descriptor>
int ExternalMomentBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, ExternalMomentBGKdynamics<T, Descriptor> >(
        "BGK_ExternalMoment");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
ExternalMomentBGKdynamics<T, Descriptor>::ExternalMomentBGKdynamics(T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
ExternalMomentBGKdynamics<T, Descriptor>::ExternalMomentBGKdynamics(
    HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
ExternalMomentBGKdynamics<T, Descriptor> *ExternalMomentBGKdynamics<T, Descriptor>::clone() const
{
    return new ExternalMomentBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int ExternalMomentBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void ExternalMomentBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T &rho = *cell.getExternal(Descriptor<T>::ExternalField::densityBeginsAt);
    T rhoBar = Descriptor<T>::rhoBar(rho);
    Array<T, Descriptor<T>::d> j;
    j.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::momentumBeginsAt));
    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, j, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void ExternalMomentBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, j, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T ExternalMomentBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
void ExternalMomentBGKdynamics<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
{
    T const &rho = *cell.getExternal(Descriptor<T>::ExternalField::densityBeginsAt);
    u.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::momentumBeginsAt));
    u /= rho;
}

template <typename T, template <typename U> class Descriptor>
T ExternalMomentBGKdynamics<T, Descriptor>::computeRhoBar(Cell<T, Descriptor> const &cell) const
{
    return momentTemplates<T, Descriptor>::get_rhoBar(cell);
}

template <typename T, template <typename U> class Descriptor>
void ExternalMomentBGKdynamics<T, Descriptor>::computeRhoBarJ(
    Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j) const
{
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
}

/* *************** Class ExternalVelocityBGKdynamics ********************************** */

template <typename T, template <typename U> class Descriptor>
int ExternalVelocityBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, ExternalVelocityBGKdynamics<T, Descriptor> >(
        "BGK_ExternalVelocity");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
ExternalVelocityBGKdynamics<T, Descriptor>::ExternalVelocityBGKdynamics(T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
ExternalVelocityBGKdynamics<T, Descriptor>::ExternalVelocityBGKdynamics(
    HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
ExternalVelocityBGKdynamics<T, Descriptor> *ExternalVelocityBGKdynamics<T, Descriptor>::clone()
    const
{
    return new ExternalVelocityBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int ExternalVelocityBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void ExternalVelocityBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar = momentTemplates<T, Descriptor>::get_rhoBar(cell);
    T rho = Descriptor<T>::fullRho(rhoBar);
    Array<T, Descriptor<T>::d> j;
    j.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));
    j *= rho;
    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, j, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void ExternalVelocityBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, j, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T ExternalVelocityBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
void ExternalVelocityBGKdynamics<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
{
    u.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));
}

template <typename T, template <typename U> class Descriptor>
T ExternalVelocityBGKdynamics<T, Descriptor>::computeRhoBar(Cell<T, Descriptor> const &cell) const
{
    return momentTemplates<T, Descriptor>::get_rhoBar(cell);
}

template <typename T, template <typename U> class Descriptor>
void ExternalVelocityBGKdynamics<T, Descriptor>::computeRhoBarJ(
    Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j) const
{
    rhoBar = momentTemplates<T, Descriptor>::get_rhoBar(cell);
    j.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));
    j *= Descriptor<T>::fullRho(rhoBar);
}

/* *************** Class QuasiIncBGKdynamics ************************************ */

template <typename T, template <typename U> class Descriptor>
int QuasiIncBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, QuasiIncBGKdynamics<T, Descriptor> >(
        "QuasiInc_BGK");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
QuasiIncBGKdynamics<T, Descriptor>::QuasiIncBGKdynamics(T omega_) :
    BGKdynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
QuasiIncBGKdynamics<T, Descriptor>::QuasiIncBGKdynamics(HierarchicUnserializer &unserializer) :
    BGKdynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
QuasiIncBGKdynamics<T, Descriptor> *QuasiIncBGKdynamics<T, Descriptor>::clone() const
{
    return new QuasiIncBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int QuasiIncBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void QuasiIncBGKdynamics<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
{
    T dummyRhoBar;
    this->computeRhoBarJ(cell, dummyRhoBar, u);
}

template <typename T, template <typename U> class Descriptor>
void QuasiIncBGKdynamics<T, Descriptor>::computeRhoBarJPiNeq(
    Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
    Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const
{
    T invRho0 = (T)1.;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq, invRho0);
}

template <typename T, template <typename U> class Descriptor>
bool QuasiIncBGKdynamics<T, Descriptor>::velIsJ() const
{
    return true;
}

/* *************** Class IncBGKdynamics ******************************************** */

template <typename T, template <typename U> class Descriptor>
int IncBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, IncBGKdynamics<T, Descriptor> >(
        "BGK_Incompressible");

/** \param omega_ relaxation parameter, related to the dynamic viscosity */
template <typename T, template <typename U> class Descriptor>
IncBGKdynamics<T, Descriptor>::IncBGKdynamics(T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
IncBGKdynamics<T, Descriptor>::IncBGKdynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>((T)1)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
IncBGKdynamics<T, Descriptor> *IncBGKdynamics<T, Descriptor>::clone() const
{
    return new IncBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int IncBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void IncBGKdynamics<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
{
    T dummyRhoBar;
    this->computeRhoBarJ(cell, dummyRhoBar, u);
}

template <typename T, template <typename U> class Descriptor>
void IncBGKdynamics<T, Descriptor>::computeRhoBarJPiNeq(
    Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
    Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const
{
    // Incompressible: rho0=1
    T invRho0 = (T)1.;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq, invRho0);
}

template <typename T, template <typename U> class Descriptor>
void IncBGKdynamics<T, Descriptor>::collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_inc_collision(cell, rhoBar, j, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void IncBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_inc_collision(cell, rhoBar, j, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T IncBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    // For the incompressible BGK dynamics, the "1/rho" pre-factor of
    // the O(Ma^2) term is equal to 1/rho0.
    // Incompressible: rho0=1
    T invRho0 = (T)1;
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho0, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
bool IncBGKdynamics<T, Descriptor>::velIsJ() const
{
    return true;
}

/* *************** Class ConstRhoBGKdynamics *************************************** */

template <typename T, template <typename U> class Descriptor>
int ConstRhoBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, ConstRhoBGKdynamics<T, Descriptor> >(
        "BGK_ConstRho");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
ConstRhoBGKdynamics<T, Descriptor>::ConstRhoBGKdynamics(T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
ConstRhoBGKdynamics<T, Descriptor>::ConstRhoBGKdynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
ConstRhoBGKdynamics<T, Descriptor> *ConstRhoBGKdynamics<T, Descriptor>::clone() const
{
    return new ConstRhoBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int ConstRhoBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void ConstRhoBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T rho = Descriptor<T>::fullRho(rhoBar);

    T deltaRho =
        -statistics.getAverage(LatticeStatistics::avRhoBar) + (1 - Descriptor<T>::SkordosFactor());
    T ratioRho = (T)1 + deltaRho / rho;

    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma2_constRho_collision(
        cell, rhoBar, j, ratioRho, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar + deltaRho, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T ConstRhoBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class RLBdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
int RLBdynamics<T, Descriptor>::id =
    meta::registerCompositeDynamics<T, Descriptor, RLBdynamics<T, Descriptor> >(
        "Regularized_generic");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
RLBdynamics<T, Descriptor>::RLBdynamics(
    Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision) :
    BulkCompositeDynamics<T, Descriptor>(baseDynamics, automaticPrepareCollision)
{ }

template <typename T, template <typename U> class Descriptor>
RLBdynamics<T, Descriptor> *RLBdynamics<T, Descriptor>::clone() const
{
    return new RLBdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int RLBdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void RLBdynamics<T, Descriptor>::completePopulations(Cell<T, Descriptor> &cell) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] = this->computeEquilibrium(iPop, rhoBar, j, jSqr)
                     + offEquilibriumTemplates<T, Descriptor>::fromPiToFneq(iPop, PiNeq);
    }
}

/* *************** Class RegularizedBGKdynamics ************************************ */

template <typename T, template <typename U> class Descriptor>
int RegularizedBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, RegularizedBGKdynamics<T, Descriptor> >(
        "Regularized_BGK");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
RegularizedBGKdynamics<T, Descriptor>::RegularizedBGKdynamics(T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
RegularizedBGKdynamics<T, Descriptor>::RegularizedBGKdynamics(
    HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
RegularizedBGKdynamics<T, Descriptor> *RegularizedBGKdynamics<T, Descriptor>::clone() const
{
    return new RegularizedBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int RegularizedBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void RegularizedBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T uSqr = dynamicsTemplates<T, Descriptor>::rlb_collision(
        cell, rhoBar, invRho, j, PiNeq, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void RegularizedBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T uSqr = dynamicsTemplates<T, Descriptor>::rlb_collision(
        cell, rhoBar, invRho, j, PiNeq, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T RegularizedBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class SecuredRegularizedBGKdynamics ************************************ */

template <typename T, template <typename U> class Descriptor>
int SecuredRegularizedBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, SecuredRegularizedBGKdynamics<T, Descriptor> >(
        "SecuredRegularized_BGK");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
SecuredRegularizedBGKdynamics<T, Descriptor>::SecuredRegularizedBGKdynamics(T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
SecuredRegularizedBGKdynamics<T, Descriptor>::SecuredRegularizedBGKdynamics(
    HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
SecuredRegularizedBGKdynamics<T, Descriptor> *SecuredRegularizedBGKdynamics<T, Descriptor>::clone()
    const
{
    return new SecuredRegularizedBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int SecuredRegularizedBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void SecuredRegularizedBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T rho = Descriptor<T>::fullRho(rhoBar);
    T softLimit = rho * 0.1;
    T hardLimit = rho * 0.2;
    for (plint i = 0; i < Descriptor<T>::d; ++i) {
        constrainValue(j[i], softLimit, hardLimit);
    }
    softLimit *= 0.05;
    hardLimit *= 0.05;
    for (plint i = 0; i < SymmetricTensor<T, Descriptor>::n; ++i) {
        constrainValue(PiNeq[i], softLimit, hardLimit);
    }
    T invRho = Descriptor<T>::invRho(rhoBar);
    T uSqr = dynamicsTemplates<T, Descriptor>::rlb_collision(
        cell, rhoBar, invRho, j, PiNeq, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void SecuredRegularizedBGKdynamics<T, Descriptor>::constrainValue(
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
void SecuredRegularizedBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T uSqr = dynamicsTemplates<T, Descriptor>::rlb_collision(
        cell, rhoBar, invRho, j, PiNeq, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T SecuredRegularizedBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class IncRegularizedBGKdynamics ************************************ */

template <typename T, template <typename U> class Descriptor>
int IncRegularizedBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, IncRegularizedBGKdynamics<T, Descriptor> >(
        "IncRegularized_BGK");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
IncRegularizedBGKdynamics<T, Descriptor>::IncRegularizedBGKdynamics(T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
IncRegularizedBGKdynamics<T, Descriptor>::IncRegularizedBGKdynamics(
    HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
IncRegularizedBGKdynamics<T, Descriptor> *IncRegularizedBGKdynamics<T, Descriptor>::clone() const
{
    return new IncRegularizedBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int IncRegularizedBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void IncRegularizedBGKdynamics<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
{
    T dummyRhoBar;
    this->computeRhoBarJ(cell, dummyRhoBar, u);
}

template <typename T, template <typename U> class Descriptor>
void IncRegularizedBGKdynamics<T, Descriptor>::computeRhoBarJPiNeq(
    Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
    Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const
{
    T invRho0 = (T)1.;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq, invRho0);
}

template <typename T, template <typename U> class Descriptor>
void IncRegularizedBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T invRho = (T)1.;
    T uSqr = dynamicsTemplates<T, Descriptor>::rlb_collision(
        cell, rhoBar, invRho, j, PiNeq, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void IncRegularizedBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);
    T invRho = (T)1.;
    T uSqr = dynamicsTemplates<T, Descriptor>::rlb_collision(
        cell, rhoBar, invRho, j, PiNeq, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T IncRegularizedBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = (T)1.;
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
bool IncRegularizedBGKdynamics<T, Descriptor>::velIsJ() const
{
    return true;
}

/* *************** Class ExternalMomentRegularizedBGKdynamics *********************** */

template <typename T, template <typename U> class Descriptor>
int ExternalMomentRegularizedBGKdynamics<T, Descriptor>::id = meta::registerGeneralDynamics<
    T, Descriptor, ExternalMomentRegularizedBGKdynamics<T, Descriptor> >(
    "Regularized_BGK_ExternalMoment");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
ExternalMomentRegularizedBGKdynamics<T, Descriptor>::ExternalMomentRegularizedBGKdynamics(
    T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
ExternalMomentRegularizedBGKdynamics<T, Descriptor>::ExternalMomentRegularizedBGKdynamics(
    HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
ExternalMomentRegularizedBGKdynamics<T, Descriptor>
    *ExternalMomentRegularizedBGKdynamics<T, Descriptor>::clone() const
{
    return new ExternalMomentRegularizedBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int ExternalMomentRegularizedBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void ExternalMomentRegularizedBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T &rho = *cell.getExternal(Descriptor<T>::ExternalField::densityBeginsAt);
    T rhoBar = Descriptor<T>::rhoBar(rho);
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    j.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::momentumBeginsAt));
    momentTemplates<T, Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T uSqr = dynamicsTemplates<T, Descriptor>::rlb_collision(
        cell, rhoBar, invRho, j, PiNeq, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void ExternalMomentRegularizedBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_PiNeq(cell, rhoBar, j, PiNeq);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T uSqr = dynamicsTemplates<T, Descriptor>::rlb_collision(
        cell, rhoBar, invRho, j, PiNeq, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T ExternalMomentRegularizedBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
void ExternalMomentRegularizedBGKdynamics<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
{
    T const &rho = *cell.getExternal(Descriptor<T>::ExternalField::densityBeginsAt);
    u.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::momentumBeginsAt));
    u /= rho;
}

template <typename T, template <typename U> class Descriptor>
T ExternalMomentRegularizedBGKdynamics<T, Descriptor>::computeRhoBar(
    Cell<T, Descriptor> const &cell) const
{
    return momentTemplates<T, Descriptor>::get_rhoBar(cell);
}

template <typename T, template <typename U> class Descriptor>
void ExternalMomentRegularizedBGKdynamics<T, Descriptor>::computeRhoBarJ(
    Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j) const
{
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
}

/* *************** Class ChopardDynamics ************************************ */

template <typename T, template <typename U> class Descriptor>
int ChopardDynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, ChopardDynamics<T, Descriptor> >(
        "Chopard_VariableCs2");

/** \param vs2_ speed of sound
 *  \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
ChopardDynamics<T, Descriptor>::ChopardDynamics(T vs2_, T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_), vs2(vs2_)
{ }

template <typename T, template <typename U> class Descriptor>
ChopardDynamics<T, Descriptor>::ChopardDynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T()), vs2(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
void ChopardDynamics<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    IsoThermalBulkDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(vs2);
}

template <typename T, template <typename U> class Descriptor>
void ChopardDynamics<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    IsoThermalBulkDynamics<T, Descriptor>::unserialize(unserializer);
    vs2 = unserializer.readValue<T>();
}

template <typename T, template <typename U> class Descriptor>
ChopardDynamics<T, Descriptor> *ChopardDynamics<T, Descriptor>::clone() const
{
    return new ChopardDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int ChopardDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void ChopardDynamics<T, Descriptor>::collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T uSqr = chopardBgkCollision(cell, rhoBar, j, vs2, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void ChopardDynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T uSqr = chopardBgkCollision(cell, rhoBar, j, vs2, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T ChopardDynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return chopardEquilibrium(iPop, rhoBar, invRho, j, jSqr, vs2);
}

template <typename T, template <typename U> class Descriptor>
T ChopardDynamics<T, Descriptor>::getParameter(plint whichParameter) const
{
    switch (whichParameter) {
    case dynamicParams::omega_shear:
        return this->getOmega();
    case dynamicParams::sqrSpeedOfSound:
        return this->getVs2();
    };
    return 0.;
}

template <typename T, template <typename U> class Descriptor>
void ChopardDynamics<T, Descriptor>::setParameter(plint whichParameter, T value)
{
    switch (whichParameter) {
    case dynamicParams::omega_shear:
        this->setOmega(value);
    case dynamicParams::sqrSpeedOfSound:
        setVs2(value);
    };
}

template <typename T, template <typename U> class Descriptor>
T ChopardDynamics<T, Descriptor>::getVs2() const
{
    return vs2;
}

template <typename T, template <typename U> class Descriptor>
void ChopardDynamics<T, Descriptor>::setVs2(T vs2_)
{
    vs2 = vs2_;
}

template <typename T, template <typename U> class Descriptor>
T ChopardDynamics<T, Descriptor>::chopardBgkCollision(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T vs2, T omega)
{
    const T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    T invRho = Descriptor<T>::invRho(rhoBar);
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] *= (T)1 - omega;
        cell[iPop] += omega * chopardEquilibrium(iPop, rhoBar, invRho, j, jSqr, vs2);
    }
    return invRho * invRho * jSqr;
}

template <typename T, template <typename U> class Descriptor>
T ChopardDynamics<T, Descriptor>::chopardEquilibrium(
    plint iPop, T rhoBar, T invRho, Array<T, Descriptor<T>::d> const &j, T jSqr, T vs2)
{
    T kappa = vs2 - Descriptor<T>::cs2;
    if (iPop == 0) {
        return Descriptor<T>::invCs2
               * (kappa * (Descriptor<T>::t[0] - (T)1)
                  + rhoBar * (Descriptor<T>::t[0] * vs2 - kappa)
                  - invRho * jSqr * Descriptor<T>::t[0] / (T)2 * Descriptor<T>::invCs2);
    } else {
        T c_j = T();
        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            c_j += Descriptor<T>::c[iPop][iD] * j[iD];
        }
        return Descriptor<T>::invCs2 * Descriptor<T>::t[iPop]
               * (kappa + rhoBar * vs2 + c_j
                  + invRho / (T)2 * (Descriptor<T>::invCs2 * c_j * c_j - jSqr));
    }
}

/* *************** Class PrecondBGKdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
int PrecondBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, PrecondBGKdynamics<T, Descriptor> >("PrecondBGK");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
PrecondBGKdynamics<T, Descriptor>::PrecondBGKdynamics(T omega_, T invGamma_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_), invGamma(invGamma_)
{ }

template <typename T, template <typename U> class Descriptor>
PrecondBGKdynamics<T, Descriptor>::PrecondBGKdynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T()), invGamma(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
PrecondBGKdynamics<T, Descriptor> *PrecondBGKdynamics<T, Descriptor>::clone() const
{
    return new PrecondBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int PrecondBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void PrecondBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T uSqr = dynamicsTemplates<T, Descriptor>::precond_bgk_ma2_collision(
        cell, rhoBar, j, this->getOmega(), invGamma);
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void PrecondBGKdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T uSqr = dynamicsTemplates<T, Descriptor>::precond_bgk_ma2_collision(
        cell, rhoBar, j, this->getOmega(), invGamma);
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T PrecondBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::precond_bgk_ma2_equilibrium(
        iPop, rhoBar, invRho, j, jSqr, invGamma);
}

template <typename T, template <typename U> class Descriptor>
void PrecondBGKdynamics<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    IsoThermalBulkDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(invGamma);
}

template <typename T, template <typename U> class Descriptor>
void PrecondBGKdynamics<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    IsoThermalBulkDynamics<T, Descriptor>::unserialize(unserializer);
    invGamma = unserializer.readValue<T>();
}

}  // namespace plb

#endif  // ISO_THERMAL_DYNAMICS_HH
