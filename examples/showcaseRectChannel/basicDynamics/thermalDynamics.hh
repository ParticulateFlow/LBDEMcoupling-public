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
#ifndef THERMAL_DYNAMICS_HH
#define THERMAL_DYNAMICS_HH

#include <algorithm>
#include <limits>

#include "basicDynamics/thermalDynamics.h"
#include "core/cell.h"
#include "core/latticeStatistics.h"
#include "latticeBoltzmann/d3q13Templates.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/externalForceTemplates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/offEquilibriumTemplates.h"

namespace plb {

/* *************** Class ThermalBulkDynamics ************************************ */

template <typename T, template <typename U> class Descriptor>
ThermalBulkDynamics<T, Descriptor>::ThermalBulkDynamics(T omega_) :
    BasicBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
void ThermalBulkDynamics<T, Descriptor>::regularize(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar) const
{
    Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> qNeq;
    momentTemplates<T, Descriptor>::compute_Qneq(cell, rhoBar, j, thetaBar, qNeq);

    typedef Descriptor<T> L;
    for (plint iPop = 0; iPop < L::q; ++iPop) {
        cell[iPop] = this->computeEquilibrium(iPop, rhoBar, j, jSqr, thetaBar);
        T fNeq = offEquilibriumTemplatesImpl<T, L>::fromPiAndQtoFneq(iPop, PiNeq, qNeq);
        cell[iPop] += fNeq;
    }
}

template <typename T, template <typename U> class Descriptor>
T ThermalBulkDynamics<T, Descriptor>::computeTemperature(Cell<T, Descriptor> const &cell) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    this->computeRhoBarJ(cell, rhoBar, j);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    return momentTemplates<T, Descriptor>::compute_theta(cell, rhoBar, jSqr);
}

template <typename T, template <typename U> class Descriptor>
void ThermalBulkDynamics<T, Descriptor>::computePiNeq(
    Cell<T, Descriptor> const &cell, Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    this->computeRhoBarJ(cell, rhoBar, j);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    T thetaBar = momentTemplates<T, Descriptor>::compute_theta(cell, rhoBar, jSqr) - (T)1;
    momentTemplates<T, Descriptor>::compute_thermal_PiNeq(cell, rhoBar, thetaBar, j, PiNeq);
}

template <typename T, template <typename U> class Descriptor>
void ThermalBulkDynamics<T, Descriptor>::computeShearStress(
    Cell<T, Descriptor> const &cell, Array<T, SymmetricTensor<T, Descriptor>::n> &stress) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    this->computeRhoBarJ(cell, rhoBar, j);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    T thetaBar = momentTemplates<T, Descriptor>::compute_theta(cell, rhoBar, jSqr) - (T)1;
    momentTemplates<T, Descriptor>::compute_thermal_PiNeq(cell, rhoBar, thetaBar, j, stress);
    T omega = cell.getDynamics().getOmega();
    stress *= ((T)0.5 * omega - (T)1);
}

template <typename T, template <typename U> class Descriptor>
void ThermalBulkDynamics<T, Descriptor>::computeHeatFlux(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &q) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    this->computeRhoBarJ(cell, rhoBar, j);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    T thetaBar = momentTemplates<T, Descriptor>::compute_theta(cell, rhoBar, jSqr) - (T)1;
    momentTemplates<T, Descriptor>::compute_heat_flux(cell, rhoBar, j, thetaBar, q);
}

template <typename T, template <typename U> class Descriptor>
T ThermalBulkDynamics<T, Descriptor>::computeEbar(Cell<T, Descriptor> const &cell) const
{
    return momentTemplates<T, Descriptor>::get_eBar(cell);
}

// TODO needs to be extended from iso-thermal to thermal case
template <typename T, template <typename U> class Descriptor>
plint ThermalBulkDynamics<T, Descriptor>::numDecomposedVariables(plint order) const
{
    plint numVariables =
        // Order 0: density + velocity + fNeq
        (order == 0)
            ? (1 + Descriptor<T>::d + 1
               + Descriptor<T>::q)  // rhoBar, j, thetaBar, fNeq
                                    // Order >=1: density + velocity + PiNeq
            : (1 + Descriptor<T>::d + 1 + SymmetricTensor<T, Descriptor>::n
               + SymmetricRankThreeTensor<T, Descriptor>::n);  // rhoBar, j, thetaBar, piNeq, qNeq

    numVariables += Descriptor<T>::ExternalField::numScalars;
    return numVariables;
}

// TODO needs to be extended from iso-thermal to thermal case
template <typename T, template <typename U> class Descriptor>
void ThermalBulkDynamics<T, Descriptor>::decompose(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData, plint order) const
{
    rawData.resize(numDecomposedVariables(order));

    if (order == 0) {
        decomposeOrder0(cell, rawData);
    } else {
        decomposeOrder1(cell, rawData);
    }
}

// TODO needs to be extended from iso-thermal to thermal case
template <typename T, template <typename U> class Descriptor>
void ThermalBulkDynamics<T, Descriptor>::recompose(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData, plint order) const
{
    PLB_PRECONDITION((plint)rawData.size() == numDecomposedVariables(order));

    if (order == 0) {
        recomposeOrder0(cell, rawData);
    } else {
        recomposeOrder1(cell, rawData);
    }
}

// TODO needs to be extended from iso-thermal to thermal case
template <typename T, template <typename U> class Descriptor>
void ThermalBulkDynamics<T, Descriptor>::rescale(
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
void ThermalBulkDynamics<T, Descriptor>::decomposeOrder0(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    T rhoBar, thetaBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j_thetaBar(cell, rhoBar, j, thetaBar);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    rawData[0] = rhoBar;
    j.to_cArray(&rawData[1]);
    rawData[1 + Descriptor<T>::d] = thetaBar;

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        rawData[2 + Descriptor<T>::d + iPop] =
            cell[iPop] - this->computeEquilibrium(iPop, rhoBar, j, jSqr);
    }

    int offset = 2 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        rawData[offset + iExt] = *cell.getExternal(iExt);
    }
}

template <typename T, template <typename U> class Descriptor>
void ThermalBulkDynamics<T, Descriptor>::decomposeOrder1(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    T rhoBar, thetaBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> qNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_thetaBar_j_PiNeq_qNeq(
        cell, rhoBar, thetaBar, j, PiNeq, qNeq);

    rawData[0] = rhoBar;
    j.to_cArray(&rawData[1]);
    rawData[1 + Descriptor<T>::d] = thetaBar;
    PiNeq.to_cArray(&rawData[2 + Descriptor<T>::d]);
    qNeq.to_cArray(&rawData
                       [2 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n
                        + SymmetricRankThreeTensor<T, Descriptor>::n]);

    int offset = 2 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n
                 + SymmetricRankThreeTensor<T, Descriptor>::n;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        rawData[offset + iExt] = *cell.getExternal(iExt);
    }
}

template <typename T, template <typename U> class Descriptor>
void ThermalBulkDynamics<T, Descriptor>::recomposeOrder0(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    T rhoBar = rawData[0];
    Array<T, Descriptor<T>::d> j;
    j.from_cArray(&rawData[1]);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    T thetaBar = rawData[1 + Descriptor<T>::d];

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] = this->computeEquilibrium(iPop, rhoBar, j, jSqr, thetaBar)
                     + rawData[2 + Descriptor<T>::d + iPop];
    }

    int offset = 2 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *cell.getExternal(iExt) = rawData[offset + iExt];
    }
}

template <typename T, template <typename U> class Descriptor>
void ThermalBulkDynamics<T, Descriptor>::recomposeOrder1(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    typedef Descriptor<T> L;

    T rhoBar, thetaBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> qNeq;

    rhoBar = rawData[0];
    j.from_cArray(&rawData[1]);
    thetaBar = rawData[1 + Descriptor<T>::d];
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    PiNeq.from_cArray(&rawData[2 + Descriptor<T>::d]);
    qNeq.from_cArray(&rawData[2 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n]);

    for (plint iPop = 0; iPop < L::q; ++iPop) {
        cell[iPop] = this->computeEquilibrium(0, rhoBar, j, jSqr, thetaBar)
                     + offEquilibriumTemplates<T, Descriptor>::fromPiAndQtoFneq(0, PiNeq, qNeq);
    }

    int offset = 2 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *cell.getExternal(iExt) = rawData[offset + iExt];
    }
}

template <typename T, template <typename U> class Descriptor>
void ThermalBulkDynamics<T, Descriptor>::rescaleOrder0(
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
void ThermalBulkDynamics<T, Descriptor>::rescaleOrder1(
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

/* *************** Class IsoThermalBGKdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
int IsoThermalBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, IsoThermalBGKdynamics<T, Descriptor> >(
        "BGK_IsoThermal");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
IsoThermalBGKdynamics<T, Descriptor>::IsoThermalBGKdynamics(T omega_) :
    ThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
IsoThermalBGKdynamics<T, Descriptor>::IsoThermalBGKdynamics(HierarchicUnserializer &unserializer) :
    ThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
IsoThermalBGKdynamics<T, Descriptor> *IsoThermalBGKdynamics<T, Descriptor>::clone() const
{
    return new IsoThermalBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int IsoThermalBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
T IsoThermalBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    T rhoThetaBar = Descriptor<T>::fullRho(rhoBar) * thetaBar;  // thetaBar = theta - 1
    return dynamicsTemplates<T, Descriptor>::bgk_ma3_equilibrium(
        iPop, rhoBar, invRho, j, jSqr, rhoThetaBar);
}

template <typename T, template <typename U> class Descriptor>
void IsoThermalBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    this->computeRhoBarJ(cell, rhoBar, j);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    T rhoThetaBar = momentTemplates<T, Descriptor>::compute_rhoThetaBar(cell, rhoBar, jSqr);

    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma3_collision(
        cell, rhoBar, j, rhoThetaBar, this->getOmega());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

/* *************** Class ThermalBGKdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
int ThermalBGKdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, ThermalBGKdynamics<T, Descriptor> >("BGK_Thermal");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
ThermalBGKdynamics<T, Descriptor>::ThermalBGKdynamics(T omega_) :
    ThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
ThermalBGKdynamics<T, Descriptor>::ThermalBGKdynamics(HierarchicUnserializer &unserializer) :
    ThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
ThermalBGKdynamics<T, Descriptor> *ThermalBGKdynamics<T, Descriptor>::clone() const
{
    return new ThermalBGKdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int ThermalBGKdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
T ThermalBGKdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    T rhoThetaBar = Descriptor<T>::fullRho(rhoBar) * thetaBar;  // thetaBar = theta - 1
    return dynamicsTemplates<T, Descriptor>::bgk_ma4_equilibrium(
        iPop, rhoBar, invRho, j, jSqr, rhoThetaBar);
}

template <typename T, template <typename U> class Descriptor>
void ThermalBGKdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    this->computeRhoBarJ(cell, rhoBar, j);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    T rhoThetaBar = momentTemplates<T, Descriptor>::compute_rhoThetaBar(cell, rhoBar, jSqr);

    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma4_collision(
        cell, rhoBar, j, rhoThetaBar, this->getOmega());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

/* *************** Class ThermalRLBdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
int ThermalRLBdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, ThermalRLBdynamics<T, Descriptor> >("RLB_Thermal");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
ThermalRLBdynamics<T, Descriptor>::ThermalRLBdynamics(T omega_) :
    ThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
ThermalRLBdynamics<T, Descriptor>::ThermalRLBdynamics(HierarchicUnserializer &unserializer) :
    ThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
ThermalRLBdynamics<T, Descriptor> *ThermalRLBdynamics<T, Descriptor>::clone() const
{
    return new ThermalRLBdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int ThermalRLBdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
T ThermalRLBdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    T rhoThetaBar = Descriptor<T>::fullRho(rhoBar) * thetaBar;  // thetaBar = theta - 1
    return dynamicsTemplates<T, Descriptor>::bgk_ma4_equilibrium(
        iPop, rhoBar, invRho, j, jSqr, rhoThetaBar);
}

template <typename T, template <typename U> class Descriptor>
void ThermalRLBdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar, thetaBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> piNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_thetaBar_j_PiNeq(
        cell, rhoBar, thetaBar, j, piNeq);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    this->regularize(cell, rhoBar, j, jSqr, piNeq, thetaBar);
    T rhoThetaBar = Descriptor<T>::fullRho(rhoBar) * thetaBar;

    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma4_collision(
        cell, rhoBar, j, rhoThetaBar, this->getOmega());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

}  // namespace plb

#endif
