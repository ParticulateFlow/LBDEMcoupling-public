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
 * A collection of dynamics classes with which a Cell object
 * can be instantiated -- generic implementation.
 * Theoretical background about these collision models can be found in
 * Coreixas et al. 'Comprehensive comparison of collision models in the
 * lattice Boltzmann framework: Theoretical investigations', PRE, 2019.
 */
#ifndef COMPREHENSIVE_ISO_THERMAL_DYNAMICS_HH
#define COMPREHENSIVE_ISO_THERMAL_DYNAMICS_HH

#include <algorithm>
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

/* *************** Class RMdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
Array<T, Descriptor<T>::numRelaxationTimes> RMdynamics<T, Descriptor>::allOmega =
    Array<T, Descriptor<T>::numRelaxationTimes>::ones();

template <typename T, template <typename U> class Descriptor>
int RMdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, RMdynamics<T, Descriptor> >("RM");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
RMdynamics<T, Descriptor>::RMdynamics(T omegaVisc_) :
    IsoThermalBulkDynamics<T, Descriptor>(omegaVisc_)
{ }

template <typename T, template <typename U> class Descriptor>
RMdynamics<T, Descriptor>::RMdynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
RMdynamics<T, Descriptor> *RMdynamics<T, Descriptor>::clone() const
{
    return new RMdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int RMdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void RMdynamics<T, Descriptor>::collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rho;
    Array<T, Descriptor<T>::q> RM, RMeq;
    comprehensiveDynamicsTemplates<T, Descriptor>::RMcomputeMoments(cell, RM, rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = RM[i + 1];
    }
    comprehensiveDynamicsTemplates<T, Descriptor>::RMcomputeEquilibriumMoments(u, RMeq);
    // All relaxation times are static class members, except for the first one
    // (the viscous relaxation time), which is dynamic.
    Array<T, Descriptor<T>::numRelaxationTimes> allOmegaTmp(allOmega);
    allOmegaTmp[0] = this->getOmega();
    comprehensiveDynamicsTemplates<T, Descriptor>::RMcollide(cell, rho, u, RM, RMeq, allOmegaTmp);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rho - 1.0, normSqr(u));
    }
}

template <typename T, template <typename U> class Descriptor>
void RMdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T rho;
    Array<T, Descriptor<T>::q> RM, RMeq;
    comprehensiveDynamicsTemplates<T, Descriptor>::RMcomputeMoments(cell, RM, rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    comprehensiveDynamicsTemplates<T, Descriptor>::RMcomputeEquilibriumMoments(u, RMeq);
    // All relaxation times are static class members, except for the first one
    // (the viscous relaxation time), which is dynamic.
    Array<T, Descriptor<T>::numRelaxationTimes> allOmegaTmp(allOmega);
    allOmegaTmp[0] = this->getOmega();
    comprehensiveDynamicsTemplates<T, Descriptor>::RMcollide(cell, rho, u, RM, RMeq, allOmegaTmp);

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rho - 1.0, normSqr(u));
    }
}

template <typename T, template <typename U> class Descriptor>
void RMdynamics<T, Descriptor>::computeEquilibria(
    Array<T, Descriptor<T>::q> &fEq, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    T thetaBar) const
{
    T rho = rhoBar + 1.0;
    // Array<T, 3> u(j[0] / rho, j[1] / rho, j[2] / rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    Array<T, Descriptor<T>::q> RMeq;
    comprehensiveDynamicsTemplates<T, Descriptor>::RMcomputeEquilibriumMoments(u, RMeq);
    Array<T, Descriptor<T>::q> eq;
    comprehensiveDynamicsTemplates<T, Descriptor>::RMcomputeEquilibrium(rho, RMeq, eq);
    for (int i = 0; i < Descriptor<T>::q; ++i) {
        fEq[i] = eq[i] - Descriptor<T>::t[i];
    }
}

template <typename T, template <typename U> class Descriptor>
T RMdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    PLB_ASSERT(iPop < Descriptor<T>::q);
    Array<T, Descriptor<T>::q> eq;
    computeEquilibria(eq, rhoBar, j, jSqr, thetaBar);
    return eq[iPop];
}

template <typename T, template <typename U> class Descriptor>
void RMdynamics<T, Descriptor>::regularize(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_regularize(
        cell, rhoBar, invRho, j, jSqr, PiNeq, this->getOmega(), this->getOmega());
}

template <typename T, template <typename U> class Descriptor>
void RMdynamics<T, Descriptor>::decomposeOrder0(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    rawData[0] = rhoBar;
    j.to_cArray(&rawData[1]);

    Array<T, Descriptor<T>::q> eq, fEq, RMeq;
    T rho = rhoBar + 1.0;
    // Array<T, 3> u(j[0]/rho, j[1]/rho, j[2]/rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    comprehensiveDynamicsTemplates<T, Descriptor>::RMcomputeEquilibriumMoments(u, RMeq);
    comprehensiveDynamicsTemplates<T, Descriptor>::RMcomputeEquilibrium(rho, RMeq, eq);
    for (int i = 0; i < Descriptor<T>::q; ++i) {
        fEq[i] = eq[i] - Descriptor<T>::t[i];
    }

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        rawData[1 + Descriptor<T>::d + iPop] = cell[iPop] - fEq[iPop];
    }

    int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        rawData[offset + iExt] = *cell.getExternal(iExt);
    }
}

template <typename T, template <typename U> class Descriptor>
void RMdynamics<T, Descriptor>::decomposeOrder1(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    PLB_ASSERT(false);
}

template <typename T, template <typename U> class Descriptor>
void RMdynamics<T, Descriptor>::recomposeOrder0(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    T rhoBar = rawData[0];
    Array<T, Descriptor<T>::d> j;
    j.from_cArray(&rawData[1]);

    Array<T, Descriptor<T>::q> eq, fEq, RMeq;
    T rho = rhoBar + 1.0;
    // Array<T, 3> u(j[0]/rho, j[1]/rho, j[2]/rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    comprehensiveDynamicsTemplates<T, Descriptor>::RMcomputeEquilibriumMoments(u, RMeq);
    comprehensiveDynamicsTemplates<T, Descriptor>::RMcomputeEquilibrium(rho, RMeq, eq);
    for (int i = 0; i < Descriptor<T>::q; ++i) {
        fEq[i] = eq[i] - Descriptor<T>::t[i];
    }

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] = fEq[iPop] + rawData[1 + Descriptor<T>::d + iPop];
    }

    int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *cell.getExternal(iExt) = rawData[offset + iExt];
    }
}

template <typename T, template <typename U> class Descriptor>
void RMdynamics<T, Descriptor>::recomposeOrder1(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    PLB_ASSERT(false);
}

/* *************** Class HMdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
Array<T, Descriptor<T>::numRelaxationTimes> HMdynamics<T, Descriptor>::allOmega =
    Array<T, Descriptor<T>::numRelaxationTimes>::ones();

template <typename T, template <typename U> class Descriptor>
int HMdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, HMdynamics<T, Descriptor> >("HM");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
HMdynamics<T, Descriptor>::HMdynamics(T omegaVisc_) :
    IsoThermalBulkDynamics<T, Descriptor>(omegaVisc_)
{ }

template <typename T, template <typename U> class Descriptor>
HMdynamics<T, Descriptor>::HMdynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
HMdynamics<T, Descriptor> *HMdynamics<T, Descriptor>::clone() const
{
    return new HMdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int HMdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void HMdynamics<T, Descriptor>::collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rho;
    Array<T, Descriptor<T>::q> HM, HMeq;
    comprehensiveDynamicsTemplates<T, Descriptor>::HMcomputeMoments(cell, HM, rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = HM[i + 1];
    }
    comprehensiveDynamicsTemplates<T, Descriptor>::HMcomputeEquilibriumMoments(u, HMeq);
    // All relaxation times are static class members, except for the first one
    // (the viscous relaxation time), which is dynamic.
    Array<T, Descriptor<T>::numRelaxationTimes> allOmegaTmp(allOmega);
    allOmegaTmp[0] = this->getOmega();
    comprehensiveDynamicsTemplates<T, Descriptor>::HMcollide(cell, rho, u, HM, HMeq, allOmegaTmp);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rho - 1.0, normSqr(u));
    }
}

template <typename T, template <typename U> class Descriptor>
void HMdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T falseRho;
    Array<T, Descriptor<T>::q> HM, HMeq;
    comprehensiveDynamicsTemplates<T, Descriptor>::HMcomputeMoments(cell, HM, falseRho);
    T rho = rhoBar + 1.0;
    for (int i = (Descriptor<T>::q + 1); i < Descriptor<T>::q; ++i) {
        HM[i] *= falseRho / rho;
    }
    // HM[0] = 1.0;
    // for (int i=0; i<Descriptor<T>::d; ++i) {
    //     HM[i+1] = j[i] / rho;
    // }
    // Array<T, Descriptor<T>::d> u;
    // for (int i=0; i<Descriptor<T>::d; ++i) {
    //     u[i] = HM[i+1];
    // }
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    comprehensiveDynamicsTemplates<T, Descriptor>::HMcomputeEquilibriumMoments(u, HMeq);
    // All relaxation times are static class members, except for the first one
    // (the viscous relaxation time), which is dynamic.
    Array<T, Descriptor<T>::numRelaxationTimes> allOmegaTmp(allOmega);
    allOmegaTmp[0] = this->getOmega();
    comprehensiveDynamicsTemplates<T, Descriptor>::HMcollide(cell, rho, u, HM, HMeq, allOmegaTmp);

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rho - 1.0, normSqr(u));
    }
}

template <typename T, template <typename U> class Descriptor>
void HMdynamics<T, Descriptor>::computeEquilibria(
    Array<T, Descriptor<T>::q> &fEq, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    T thetaBar) const
{
    T rho = rhoBar + 1.0;
    // Array<T, 3> u(j[0] / rho, j[1] / rho, j[2] / rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    Array<T, Descriptor<T>::q> HMeq;
    comprehensiveDynamicsTemplates<T, Descriptor>::HMcomputeEquilibriumMoments(u, HMeq);
    Array<T, Descriptor<T>::q> eq;
    comprehensiveDynamicsTemplates<T, Descriptor>::HMcomputeEquilibrium(rho, HMeq, eq);
    for (int i = 0; i < Descriptor<T>::q; ++i) {
        fEq[i] = eq[i] - Descriptor<T>::t[i];
    }
}

template <typename T, template <typename U> class Descriptor>
T HMdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    PLB_ASSERT(iPop < Descriptor<T>::q);
    Array<T, Descriptor<T>::q> eq;
    computeEquilibria(eq, rhoBar, j, jSqr, thetaBar);
    return eq[iPop];
}

template <typename T, template <typename U> class Descriptor>
void HMdynamics<T, Descriptor>::regularize(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_regularize(
        cell, rhoBar, invRho, j, jSqr, PiNeq, this->getOmega(), this->getOmega());
}

template <typename T, template <typename U> class Descriptor>
void HMdynamics<T, Descriptor>::decomposeOrder0(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    rawData[0] = rhoBar;
    j.to_cArray(&rawData[1]);

    Array<T, Descriptor<T>::q> eq, fEq, HMeq;
    T rho = rhoBar + 1.0;
    // Array<T, 3> u(j[0]/rho, j[1]/rho, j[2]/rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    comprehensiveDynamicsTemplates<T, Descriptor>::HMcomputeEquilibriumMoments(u, HMeq);
    comprehensiveDynamicsTemplates<T, Descriptor>::HMcomputeEquilibrium(rho, HMeq, eq);
    for (int i = 0; i < Descriptor<T>::q; ++i) {
        fEq[i] = eq[i] - Descriptor<T>::t[i];
    }

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        rawData[1 + Descriptor<T>::d + iPop] = cell[iPop] - fEq[iPop];
    }

    int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        rawData[offset + iExt] = *cell.getExternal(iExt);
    }
}

template <typename T, template <typename U> class Descriptor>
void HMdynamics<T, Descriptor>::decomposeOrder1(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    PLB_ASSERT(false);
}

template <typename T, template <typename U> class Descriptor>
void HMdynamics<T, Descriptor>::recomposeOrder0(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    T rhoBar = rawData[0];
    Array<T, Descriptor<T>::d> j;
    j.from_cArray(&rawData[1]);

    Array<T, Descriptor<T>::q> eq, fEq, HMeq;
    T rho = rhoBar + 1.0;
    // Array<T, 3> u(j[0]/rho, j[1]/rho, j[2]/rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    comprehensiveDynamicsTemplates<T, Descriptor>::HMcomputeEquilibriumMoments(u, HMeq);
    comprehensiveDynamicsTemplates<T, Descriptor>::HMcomputeEquilibrium(rho, HMeq, eq);
    for (int i = 0; i < Descriptor<T>::q; ++i) {
        fEq[i] = eq[i] - Descriptor<T>::t[i];
    }

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] = fEq[iPop] + rawData[1 + Descriptor<T>::d + iPop];
    }

    int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *cell.getExternal(iExt) = rawData[offset + iExt];
    }
}

template <typename T, template <typename U> class Descriptor>
void HMdynamics<T, Descriptor>::recomposeOrder1(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    PLB_ASSERT(false);
}

/* *************** Class CMdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
Array<T, Descriptor<T>::numRelaxationTimes> CMdynamics<T, Descriptor>::allOmega =
    Array<T, Descriptor<T>::numRelaxationTimes>::ones();

template <typename T, template <typename U> class Descriptor>
int CMdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, CMdynamics<T, Descriptor> >("CM");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
CMdynamics<T, Descriptor>::CMdynamics(T omegaVisc_) :
    IsoThermalBulkDynamics<T, Descriptor>(omegaVisc_)
{ }

template <typename T, template <typename U> class Descriptor>
CMdynamics<T, Descriptor>::CMdynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
CMdynamics<T, Descriptor> *CMdynamics<T, Descriptor>::clone() const
{
    return new CMdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int CMdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void CMdynamics<T, Descriptor>::collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rho;
    Array<T, Descriptor<T>::d> u;
    Array<T, Descriptor<T>::q> CM, CMeq;
    // Here we cannot use CMs to compute u... That is why it is done during CMcomputeMoments
    comprehensiveDynamicsTemplates<T, Descriptor>::CMcomputeMoments(cell, CM, rho, u);
    comprehensiveDynamicsTemplates<T, Descriptor>::CMcomputeEquilibriumMoments(CMeq);
    // All relaxation times are static class members, except for the first one
    // (the viscous relaxation time), which is dynamic.
    Array<T, Descriptor<T>::numRelaxationTimes> allOmegaTmp(allOmega);
    allOmegaTmp[0] = this->getOmega();
    comprehensiveDynamicsTemplates<T, Descriptor>::CMcollide(cell, rho, u, CM, CMeq, allOmegaTmp);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rho - 1.0, normSqr(u));
    }
}

template <typename T, template <typename U> class Descriptor>
void CMdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T falseRho;
    Array<T, Descriptor<T>::d> u;
    Array<T, Descriptor<T>::q> CM, CMeq;
    comprehensiveDynamicsTemplates<T, Descriptor>::CMcomputeMoments(cell, CM, falseRho, u);
    T rho = rhoBar + 1.0;
    for (int i = (Descriptor<T>::q + 1); i < Descriptor<T>::q; ++i) {
        CM[i] *= falseRho / rho;
    }
    // CM[0] = 1.0;
    // for (int i=0; i<Descriptor<T>::d; ++i) {
    //     CM[i+1] = 0.;
    // }
    // Array<T, Descriptor<T>::d> u;
    // for (int i=0; i<Descriptor<T>::d; ++i) {
    //     u[i] = j[i] / rho;
    // }
    comprehensiveDynamicsTemplates<T, Descriptor>::CMcomputeEquilibriumMoments(CMeq);
    // All relaxation times are static class members, except for the first one
    // (the viscous relaxation time), which is dynamic.
    Array<T, Descriptor<T>::numRelaxationTimes> allOmegaTmp(allOmega);
    allOmegaTmp[0] = this->getOmega();
    comprehensiveDynamicsTemplates<T, Descriptor>::CMcollide(cell, rho, u, CM, CMeq, allOmegaTmp);

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rho - 1.0, normSqr(u));
    }
}

template <typename T, template <typename U> class Descriptor>
void CMdynamics<T, Descriptor>::computeEquilibria(
    Array<T, Descriptor<T>::q> &fEq, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    T thetaBar) const
{
    T rho = rhoBar + 1.0;
    // Array<T, 3> u(j[0] / rho, j[1] / rho, j[2] / rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    Array<T, Descriptor<T>::q> CMeq;
    comprehensiveDynamicsTemplates<T, Descriptor>::CMcomputeEquilibriumMoments(CMeq);
    Array<T, Descriptor<T>::q> eq;
    comprehensiveDynamicsTemplates<T, Descriptor>::CMcomputeEquilibrium(rho, u, CMeq, eq);
    for (int i = 0; i < Descriptor<T>::q; ++i) {
        fEq[i] = eq[i] - Descriptor<T>::t[i];
    }
}

template <typename T, template <typename U> class Descriptor>
T CMdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    PLB_ASSERT(iPop < Descriptor<T>::q);
    Array<T, Descriptor<T>::q> eq;
    computeEquilibria(eq, rhoBar, j, jSqr, thetaBar);
    return eq[iPop];
}

template <typename T, template <typename U> class Descriptor>
void CMdynamics<T, Descriptor>::regularize(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_regularize(
        cell, rhoBar, invRho, j, jSqr, PiNeq, this->getOmega(), this->getOmega());
}

template <typename T, template <typename U> class Descriptor>
void CMdynamics<T, Descriptor>::decomposeOrder0(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    rawData[0] = rhoBar;
    j.to_cArray(&rawData[1]);

    Array<T, Descriptor<T>::q> eq, fEq, CMeq;
    T rho = rhoBar + 1.0;
    // Array<T, 3> u(j[0]/rho, j[1]/rho, j[2]/rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    comprehensiveDynamicsTemplates<T, Descriptor>::CMcomputeEquilibriumMoments(CMeq);
    comprehensiveDynamicsTemplates<T, Descriptor>::CMcomputeEquilibrium(rho, u, CMeq, eq);
    for (int i = 0; i < Descriptor<T>::q; ++i) {
        fEq[i] = eq[i] - Descriptor<T>::t[i];
    }

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        rawData[1 + Descriptor<T>::d + iPop] = cell[iPop] - fEq[iPop];
    }

    int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        rawData[offset + iExt] = *cell.getExternal(iExt);
    }
}

template <typename T, template <typename U> class Descriptor>
void CMdynamics<T, Descriptor>::decomposeOrder1(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    PLB_ASSERT(false);
}

template <typename T, template <typename U> class Descriptor>
void CMdynamics<T, Descriptor>::recomposeOrder0(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    T rhoBar = rawData[0];
    Array<T, Descriptor<T>::d> j;
    j.from_cArray(&rawData[1]);

    Array<T, Descriptor<T>::q> eq, fEq, CMeq;
    T rho = rhoBar + 1.0;
    // Array<T, 3> u(j[0]/rho, j[1]/rho, j[2]/rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    comprehensiveDynamicsTemplates<T, Descriptor>::CMcomputeEquilibriumMoments(CMeq);
    comprehensiveDynamicsTemplates<T, Descriptor>::CMcomputeEquilibrium(rho, u, CMeq, eq);
    for (int i = 0; i < Descriptor<T>::q; ++i) {
        fEq[i] = eq[i] - Descriptor<T>::t[i];
    }

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] = fEq[iPop] + rawData[1 + Descriptor<T>::d + iPop];
    }

    int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *cell.getExternal(iExt) = rawData[offset + iExt];
    }
}

template <typename T, template <typename U> class Descriptor>
void CMdynamics<T, Descriptor>::recomposeOrder1(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    PLB_ASSERT(false);
}

/* *************** Class CHMdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
Array<T, Descriptor<T>::numRelaxationTimes> CHMdynamics<T, Descriptor>::allOmega =
    Array<T, Descriptor<T>::numRelaxationTimes>::ones();

template <typename T, template <typename U> class Descriptor>
int CHMdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, CHMdynamics<T, Descriptor> >("CHM");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
CHMdynamics<T, Descriptor>::CHMdynamics(T omegaVisc_) :
    IsoThermalBulkDynamics<T, Descriptor>(omegaVisc_)
{ }

template <typename T, template <typename U> class Descriptor>
CHMdynamics<T, Descriptor>::CHMdynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
CHMdynamics<T, Descriptor> *CHMdynamics<T, Descriptor>::clone() const
{
    return new CHMdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int CHMdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void CHMdynamics<T, Descriptor>::collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rho;
    Array<T, Descriptor<T>::d> u;
    Array<T, Descriptor<T>::q> CHM, CHMeq;
    comprehensiveDynamicsTemplates<T, Descriptor>::CHMcomputeMoments(cell, CHM, rho, u);
    // Here we cannot use CMs to compute u... That is why it is done during CHMcomputeMoments
    comprehensiveDynamicsTemplates<T, Descriptor>::CHMcomputeEquilibriumMoments(CHMeq);
    // All relaxation times are static class members, except for the first one
    // (the viscous relaxation time), which is dynamic.
    Array<T, Descriptor<T>::numRelaxationTimes> allOmegaTmp(allOmega);
    allOmegaTmp[0] = this->getOmega();
    comprehensiveDynamicsTemplates<T, Descriptor>::CHMcollide(
        cell, rho, u, CHM, CHMeq, allOmegaTmp);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rho - 1.0, normSqr(u));
    }
}

template <typename T, template <typename U> class Descriptor>
void CHMdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T falseRho;
    Array<T, Descriptor<T>::d> u;
    Array<T, Descriptor<T>::q> CHM, CHMeq;
    comprehensiveDynamicsTemplates<T, Descriptor>::CHMcomputeMoments(cell, CHM, falseRho, u);
    T rho = rhoBar + 1.0;
    for (int i = (Descriptor<T>::q + 1); i < Descriptor<T>::q; ++i) {
        CHM[i] *= falseRho / rho;
    }
    // CHM[0] = 1.0;
    // for (int i=0; i<Descriptor<T>::d; ++i) {
    //     CHM[i+1] = 0.;
    // }
    // Array<T, Descriptor<T>::d> u;
    // for (int i=0; i<Descriptor<T>::d; ++i) {
    //     u[i] = j[i] / rho;
    // }
    comprehensiveDynamicsTemplates<T, Descriptor>::CHMcomputeEquilibriumMoments(CHMeq);
    // All relaxation times are static class members, except for the first one
    // (the viscous relaxation time), which is dynamic.
    Array<T, Descriptor<T>::numRelaxationTimes> allOmegaTmp(allOmega);
    allOmegaTmp[0] = this->getOmega();
    comprehensiveDynamicsTemplates<T, Descriptor>::CHMcollide(
        cell, rho, u, CHM, CHMeq, allOmegaTmp);

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rho - 1.0, normSqr(u));
    }
}

template <typename T, template <typename U> class Descriptor>
void CHMdynamics<T, Descriptor>::computeEquilibria(
    Array<T, Descriptor<T>::q> &fEq, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    T thetaBar) const
{
    T rho = rhoBar + 1.0;
    // Array<T, 3> u(j[0] / rho, j[1] / rho, j[2] / rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    Array<T, Descriptor<T>::q> CHMeq;
    comprehensiveDynamicsTemplates<T, Descriptor>::CHMcomputeEquilibriumMoments(CHMeq);
    Array<T, Descriptor<T>::q> eq;
    comprehensiveDynamicsTemplates<T, Descriptor>::CHMcomputeEquilibrium(rho, u, CHMeq, eq);
    for (int i = 0; i < Descriptor<T>::q; ++i) {
        fEq[i] = eq[i] - Descriptor<T>::t[i];
    }
}

template <typename T, template <typename U> class Descriptor>
T CHMdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    PLB_ASSERT(iPop < Descriptor<T>::q);
    Array<T, Descriptor<T>::q> eq;
    computeEquilibria(eq, rhoBar, j, jSqr, thetaBar);
    return eq[iPop];
}

template <typename T, template <typename U> class Descriptor>
void CHMdynamics<T, Descriptor>::regularize(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_regularize(
        cell, rhoBar, invRho, j, jSqr, PiNeq, this->getOmega(), this->getOmega());
}

template <typename T, template <typename U> class Descriptor>
void CHMdynamics<T, Descriptor>::decomposeOrder0(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    rawData[0] = rhoBar;
    j.to_cArray(&rawData[1]);

    Array<T, Descriptor<T>::q> eq, fEq, CHMeq;
    T rho = rhoBar + 1.0;
    // Array<T, 3> u(j[0]/rho, j[1]/rho, j[2]/rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    comprehensiveDynamicsTemplates<T, Descriptor>::CHMcomputeEquilibriumMoments(CHMeq);
    comprehensiveDynamicsTemplates<T, Descriptor>::CHMcomputeEquilibrium(rho, u, CHMeq, eq);
    for (int i = 0; i < Descriptor<T>::q; ++i) {
        fEq[i] = eq[i] - Descriptor<T>::t[i];
    }

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        rawData[1 + Descriptor<T>::d + iPop] = cell[iPop] - fEq[iPop];
    }

    int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        rawData[offset + iExt] = *cell.getExternal(iExt);
    }
}

template <typename T, template <typename U> class Descriptor>
void CHMdynamics<T, Descriptor>::decomposeOrder1(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    PLB_ASSERT(false);
}

template <typename T, template <typename U> class Descriptor>
void CHMdynamics<T, Descriptor>::recomposeOrder0(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    T rhoBar = rawData[0];
    Array<T, Descriptor<T>::d> j;
    j.from_cArray(&rawData[1]);

    Array<T, Descriptor<T>::q> eq, fEq, CHMeq;
    T rho = rhoBar + 1.0;
    // Array<T, 3> u(j[0]/rho, j[1]/rho, j[2]/rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    comprehensiveDynamicsTemplates<T, Descriptor>::CHMcomputeEquilibriumMoments(CHMeq);
    comprehensiveDynamicsTemplates<T, Descriptor>::CHMcomputeEquilibrium(rho, u, CHMeq, eq);
    for (int i = 0; i < Descriptor<T>::q; ++i) {
        fEq[i] = eq[i] - Descriptor<T>::t[i];
    }

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] = fEq[iPop] + rawData[1 + Descriptor<T>::d + iPop];
    }

    int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *cell.getExternal(iExt) = rawData[offset + iExt];
    }
}

template <typename T, template <typename U> class Descriptor>
void CHMdynamics<T, Descriptor>::recomposeOrder1(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    PLB_ASSERT(false);
}

/* *************** Class Kdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
Array<T, Descriptor<T>::numRelaxationTimes> Kdynamics<T, Descriptor>::allOmega =
    Array<T, Descriptor<T>::numRelaxationTimes>::ones();

template <typename T, template <typename U> class Descriptor>
int Kdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, Kdynamics<T, Descriptor> >("K");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
Kdynamics<T, Descriptor>::Kdynamics(T omegaVisc_) :
    IsoThermalBulkDynamics<T, Descriptor>(omegaVisc_)
{ }

template <typename T, template <typename U> class Descriptor>
Kdynamics<T, Descriptor>::Kdynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
Kdynamics<T, Descriptor> *Kdynamics<T, Descriptor>::clone() const
{
    return new Kdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int Kdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void Kdynamics<T, Descriptor>::collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rho;
    Array<T, Descriptor<T>::d> u;
    Array<T, Descriptor<T>::q> K, Keq;
    comprehensiveDynamicsTemplates<T, Descriptor>::KcomputeMoments(cell, K, rho, u);
    // Here we cannot use CMs to compute u... That is why it is done during KcomputeMoments
    comprehensiveDynamicsTemplates<T, Descriptor>::KcomputeEquilibriumMoments(u, Keq);
    // All relaxation times are static class members, except for the first one
    // (the viscous relaxation time), which is dynamic.
    Array<T, Descriptor<T>::numRelaxationTimes> allOmegaTmp(allOmega);
    allOmegaTmp[0] = this->getOmega();
    comprehensiveDynamicsTemplates<T, Descriptor>::Kcollide(cell, rho, u, K, Keq, allOmegaTmp);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rho - 1.0, normSqr(u));
    }
}

template <typename T, template <typename U> class Descriptor>
void Kdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T falseRho;
    Array<T, Descriptor<T>::d> u;
    Array<T, Descriptor<T>::q> K, Keq;
    comprehensiveDynamicsTemplates<T, Descriptor>::KcomputeMoments(cell, K, falseRho, u);
    T rho = rhoBar + 1.0;
    for (int i = (Descriptor<T>::q + 1); i < Descriptor<T>::q; ++i) {
        K[i] *= falseRho / rho;
    }
    // K[0] = 1.0;
    // for (int i=0; i<Descriptor<T>::d; ++i) {
    //     K[i+1] = 0.;
    // }
    // Array<T, Descriptor<T>::d> u;
    // for (int i=0; i<Descriptor<T>::d; ++i) {
    //     u[i] = j[i] / rho;
    // }
    comprehensiveDynamicsTemplates<T, Descriptor>::KcomputeEquilibriumMoments(u, Keq);
    // All relaxation times are static class members, except for the first one
    // (the viscous relaxation time), which is dynamic.
    Array<T, Descriptor<T>::numRelaxationTimes> allOmegaTmp(allOmega);
    allOmegaTmp[0] = this->getOmega();
    comprehensiveDynamicsTemplates<T, Descriptor>::Kcollide(cell, rho, u, K, Keq, allOmegaTmp);

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rho - 1.0, normSqr(u));
    }
}

template <typename T, template <typename U> class Descriptor>
void Kdynamics<T, Descriptor>::computeEquilibria(
    Array<T, Descriptor<T>::q> &fEq, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    T thetaBar) const
{
    T rho = rhoBar + 1.0;
    // Array<T, 3> u(j[0] / rho, j[1] / rho, j[2] / rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    Array<T, Descriptor<T>::q> Keq;
    comprehensiveDynamicsTemplates<T, Descriptor>::KcomputeEquilibriumMoments(u, Keq);
    Array<T, Descriptor<T>::q> eq;
    comprehensiveDynamicsTemplates<T, Descriptor>::KcomputeEquilibrium(rho, u, Keq, eq);
    for (int i = 0; i < Descriptor<T>::q; ++i) {
        fEq[i] = eq[i] - Descriptor<T>::t[i];
    }
}

template <typename T, template <typename U> class Descriptor>
T Kdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    PLB_ASSERT(iPop < Descriptor<T>::q);
    Array<T, Descriptor<T>::q> eq;
    computeEquilibria(eq, rhoBar, j, jSqr, thetaBar);
    return eq[iPop];
}

template <typename T, template <typename U> class Descriptor>
void Kdynamics<T, Descriptor>::regularize(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_regularize(
        cell, rhoBar, invRho, j, jSqr, PiNeq, this->getOmega(), this->getOmega());
}

template <typename T, template <typename U> class Descriptor>
void Kdynamics<T, Descriptor>::decomposeOrder0(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    rawData[0] = rhoBar;
    j.to_cArray(&rawData[1]);

    Array<T, Descriptor<T>::q> eq, fEq, Keq;
    T rho = rhoBar + 1.0;
    // Array<T, 3> u(j[0]/rho, j[1]/rho, j[2]/rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    comprehensiveDynamicsTemplates<T, Descriptor>::KcomputeEquilibriumMoments(u, Keq);
    comprehensiveDynamicsTemplates<T, Descriptor>::KcomputeEquilibrium(rho, u, Keq, eq);
    for (int i = 0; i < Descriptor<T>::q; ++i) {
        fEq[i] = eq[i] - Descriptor<T>::t[i];
    }

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        rawData[1 + Descriptor<T>::d + iPop] = cell[iPop] - fEq[iPop];
    }

    int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        rawData[offset + iExt] = *cell.getExternal(iExt);
    }
}

template <typename T, template <typename U> class Descriptor>
void Kdynamics<T, Descriptor>::decomposeOrder1(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    PLB_ASSERT(false);
}

template <typename T, template <typename U> class Descriptor>
void Kdynamics<T, Descriptor>::recomposeOrder0(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    T rhoBar = rawData[0];
    Array<T, Descriptor<T>::d> j;
    j.from_cArray(&rawData[1]);

    Array<T, Descriptor<T>::q> eq, fEq, Keq;
    T rho = rhoBar + 1.0;
    // Array<T, 3> u(j[0]/rho, j[1]/rho, j[2]/rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    comprehensiveDynamicsTemplates<T, Descriptor>::KcomputeEquilibriumMoments(u, Keq);
    comprehensiveDynamicsTemplates<T, Descriptor>::KcomputeEquilibrium(rho, u, Keq, eq);
    for (int i = 0; i < Descriptor<T>::q; ++i) {
        fEq[i] = eq[i] - Descriptor<T>::t[i];
    }

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] = fEq[iPop] + rawData[1 + Descriptor<T>::d + iPop];
    }

    int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *cell.getExternal(iExt) = rawData[offset + iExt];
    }
}

template <typename T, template <typename U> class Descriptor>
void Kdynamics<T, Descriptor>::recomposeOrder1(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    PLB_ASSERT(false);
}

/* *************** Class GHdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
Array<T, Descriptor<T>::numRelaxationTimes> GHdynamics<T, Descriptor>::allOmega =
    Array<T, Descriptor<T>::numRelaxationTimes>::ones();

template <typename T, template <typename U> class Descriptor>
int GHdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, GHdynamics<T, Descriptor> >("GH");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
GHdynamics<T, Descriptor>::GHdynamics(T omegaVisc_) :
    IsoThermalBulkDynamics<T, Descriptor>(omegaVisc_)
{ }

template <typename T, template <typename U> class Descriptor>
GHdynamics<T, Descriptor>::GHdynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
GHdynamics<T, Descriptor> *GHdynamics<T, Descriptor>::clone() const
{
    return new GHdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int GHdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void GHdynamics<T, Descriptor>::collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rho;
    Array<T, Descriptor<T>::q> GH, GHeq;
    comprehensiveDynamicsTemplates<T, Descriptor>::GHcomputeMoments(cell, GH, rho);
    // Array<T, Descriptor<T>::d> u;
    // for (int i=0; i<Descriptor<T>::d; ++i) {
    //     u[i] = GH[i+1];
    // }
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = GH[i + 1];
    }
    comprehensiveDynamicsTemplates<T, Descriptor>::GHcomputeEquilibriumMoments(u, GHeq);
    // All relaxation times are static class members, except for the first one
    // (the viscous relaxation time), which is dynamic.
    Array<T, Descriptor<T>::numRelaxationTimes> allOmegaTmp(allOmega);
    allOmegaTmp[0] = this->getOmega();
    comprehensiveDynamicsTemplates<T, Descriptor>::GHcollide(cell, rho, u, GH, GHeq, allOmegaTmp);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rho - 1.0, normSqr(u));
    }
}

template <typename T, template <typename U> class Descriptor>
void GHdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T falseRho;
    Array<T, Descriptor<T>::q> GH, GHeq;
    comprehensiveDynamicsTemplates<T, Descriptor>::GHcomputeMoments(cell, GH, falseRho);
    T rho = rhoBar + 1.0;
    for (int i = (Descriptor<T>::q + 1); i < Descriptor<T>::q; ++i) {
        GH[i] *= falseRho / rho;
    }
    // GH[0] = 1.0;
    // for (int i=0; i<Descriptor<T>::d; ++i) {
    //     GH[i+1] = j[i] / rho;
    // }
    // Array<T, Descriptor<T>::d> u;
    // for (int i=0; i<Descriptor<T>::d; ++i) {
    //     u[i] = GH[i+1];
    // }
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }

    comprehensiveDynamicsTemplates<T, Descriptor>::GHcomputeEquilibriumMoments(u, GHeq);
    // All relaxation times are static class members, except for the first one
    // (the viscous relaxation time), which is dynamic.
    Array<T, Descriptor<T>::numRelaxationTimes> allOmegaTmp(allOmega);
    allOmegaTmp[0] = this->getOmega();
    comprehensiveDynamicsTemplates<T, Descriptor>::GHcollide(cell, rho, u, GH, GHeq, allOmegaTmp);

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rho - 1.0, normSqr(u));
    }
}

template <typename T, template <typename U> class Descriptor>
void GHdynamics<T, Descriptor>::computeEquilibria(
    Array<T, Descriptor<T>::q> &fEq, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    T thetaBar) const
{
    T rho = rhoBar + 1.0;
    // Array<T, 3> u(j[0] / rho, j[1] / rho, j[2] / rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    Array<T, Descriptor<T>::q> GHeq;
    comprehensiveDynamicsTemplates<T, Descriptor>::GHcomputeEquilibriumMoments(u, GHeq);
    Array<T, Descriptor<T>::q> eq;
    comprehensiveDynamicsTemplates<T, Descriptor>::GHcomputeEquilibrium(rho, GHeq, eq);
    for (int i = 0; i < Descriptor<T>::q; ++i) {
        fEq[i] = eq[i] - Descriptor<T>::t[i];
    }
}

template <typename T, template <typename U> class Descriptor>
T GHdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    PLB_ASSERT(iPop < Descriptor<T>::q);
    Array<T, Descriptor<T>::q> eq;
    computeEquilibria(eq, rhoBar, j, jSqr, thetaBar);
    return eq[iPop];
}

template <typename T, template <typename U> class Descriptor>
void GHdynamics<T, Descriptor>::regularize(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_regularize(
        cell, rhoBar, invRho, j, jSqr, PiNeq, this->getOmega(), this->getOmega());
}

template <typename T, template <typename U> class Descriptor>
void GHdynamics<T, Descriptor>::decomposeOrder0(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    rawData[0] = rhoBar;
    j.to_cArray(&rawData[1]);

    Array<T, Descriptor<T>::q> eq, fEq, GHeq;
    T rho = rhoBar + 1.0;
    // Array<T, 3> u(j[0]/rho, j[1]/rho, j[2]/rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    comprehensiveDynamicsTemplates<T, Descriptor>::GHcomputeEquilibriumMoments(u, GHeq);
    comprehensiveDynamicsTemplates<T, Descriptor>::GHcomputeEquilibrium(rho, GHeq, eq);
    for (int i = 0; i < Descriptor<T>::q; ++i) {
        fEq[i] = eq[i] - Descriptor<T>::t[i];
    }

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        rawData[1 + Descriptor<T>::d + iPop] = cell[iPop] - fEq[iPop];
    }

    int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        rawData[offset + iExt] = *cell.getExternal(iExt);
    }
}

template <typename T, template <typename U> class Descriptor>
void GHdynamics<T, Descriptor>::decomposeOrder1(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    PLB_ASSERT(false);
}

template <typename T, template <typename U> class Descriptor>
void GHdynamics<T, Descriptor>::recomposeOrder0(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    T rhoBar = rawData[0];
    Array<T, Descriptor<T>::d> j;
    j.from_cArray(&rawData[1]);

    Array<T, Descriptor<T>::q> eq, fEq, GHeq;
    T rho = rhoBar + 1.0;
    // Array<T, 3> u(j[0]/rho, j[1]/rho, j[2]/rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    comprehensiveDynamicsTemplates<T, Descriptor>::GHcomputeEquilibriumMoments(u, GHeq);
    comprehensiveDynamicsTemplates<T, Descriptor>::GHcomputeEquilibrium(rho, GHeq, eq);
    for (int i = 0; i < Descriptor<T>::q; ++i) {
        fEq[i] = eq[i] - Descriptor<T>::t[i];
    }

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] = fEq[iPop] + rawData[1 + Descriptor<T>::d + iPop];
    }

    int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *cell.getExternal(iExt) = rawData[offset + iExt];
    }
}

template <typename T, template <typename U> class Descriptor>
void GHdynamics<T, Descriptor>::recomposeOrder1(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    PLB_ASSERT(false);
}

/* *************** Class RRdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
Array<T, Descriptor<T>::numRelaxationTimes> RRdynamics<T, Descriptor>::allOmega =
    Array<T, Descriptor<T>::numRelaxationTimes>::ones();

template <typename T, template <typename U> class Descriptor>
int RRdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, RRdynamics<T, Descriptor> >("RR");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
RRdynamics<T, Descriptor>::RRdynamics(T omegaVisc_) :
    IsoThermalBulkDynamics<T, Descriptor>(omegaVisc_)
{ }

template <typename T, template <typename U> class Descriptor>
RRdynamics<T, Descriptor>::RRdynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
RRdynamics<T, Descriptor> *RRdynamics<T, Descriptor>::clone() const
{
    return new RRdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int RRdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void RRdynamics<T, Descriptor>::collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rho;
    Array<T, Descriptor<T>::q> RR, RReq;
    comprehensiveDynamicsTemplates<T, Descriptor>::RRcomputeMoments(cell, RR, rho);
    // Array<T, Descriptor<T>::d> u;
    // for (int i=0; i<Descriptor<T>::d; ++i) {
    //     u[i] = RR[i+1];
    // }
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = RR[i + 1];
    }
    comprehensiveDynamicsTemplates<T, Descriptor>::RRcomputeEquilibriumMoments(u, RReq);
    // All relaxation times are static class members, except for the first one
    // (the viscous relaxation time), which is dynamic.
    Array<T, Descriptor<T>::numRelaxationTimes> allOmegaTmp(allOmega);
    allOmegaTmp[0] = this->getOmega();
    comprehensiveDynamicsTemplates<T, Descriptor>::RRcollide(cell, rho, u, RR, RReq, allOmegaTmp);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rho - 1.0, normSqr(u));
    }
}

template <typename T, template <typename U> class Descriptor>
void RRdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T falseRho;
    Array<T, Descriptor<T>::q> RR, RReq;
    comprehensiveDynamicsTemplates<T, Descriptor>::RRcomputeMoments(cell, RR, falseRho);
    T rho = rhoBar + 1.0;
    for (int i = (Descriptor<T>::q + 1); i < Descriptor<T>::q; ++i) {
        RR[i] *= falseRho / rho;
    }
    // RR[0] = 1.0;
    // for (int i=0; i<Descriptor<T>::d; ++i) {
    //     RR[i+1] = j[i] / rho;
    // }
    // Array<T, Descriptor<T>::d> u;
    // for (int i=0; i<Descriptor<T>::d; ++i) {
    //     u[i] = RR[i+1];
    // }
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }

    comprehensiveDynamicsTemplates<T, Descriptor>::RRcomputeEquilibriumMoments(u, RReq);
    // All relaxation times are static class members, except for the first one
    // (the viscous relaxation time), which is dynamic.
    Array<T, Descriptor<T>::numRelaxationTimes> allOmegaTmp(allOmega);
    allOmegaTmp[0] = this->getOmega();
    comprehensiveDynamicsTemplates<T, Descriptor>::RRcollide(cell, rho, u, RR, RReq, allOmegaTmp);

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rho - 1.0, normSqr(u));
    }
}

template <typename T, template <typename U> class Descriptor>
void RRdynamics<T, Descriptor>::computeEquilibria(
    Array<T, Descriptor<T>::q> &fEq, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    T thetaBar) const
{
    T rho = rhoBar + 1.0;
    // Array<T, 3> u(j[0] / rho, j[1] / rho, j[2] / rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    Array<T, Descriptor<T>::q> RReq;
    comprehensiveDynamicsTemplates<T, Descriptor>::RRcomputeEquilibriumMoments(u, RReq);
    Array<T, Descriptor<T>::q> eq;
    comprehensiveDynamicsTemplates<T, Descriptor>::RRcomputeEquilibrium(rho, RReq, eq);
    for (int i = 0; i < Descriptor<T>::q; ++i) {
        fEq[i] = eq[i] - Descriptor<T>::t[i];
    }
}

template <typename T, template <typename U> class Descriptor>
T RRdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    PLB_ASSERT(iPop < Descriptor<T>::q);
    Array<T, Descriptor<T>::q> eq;
    computeEquilibria(eq, rhoBar, j, jSqr, thetaBar);
    return eq[iPop];
}

template <typename T, template <typename U> class Descriptor>
void RRdynamics<T, Descriptor>::regularize(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T, Descriptor>::complete_bgk_ma2_regularize(
        cell, rhoBar, invRho, j, jSqr, PiNeq, this->getOmega(), this->getOmega());
}

template <typename T, template <typename U> class Descriptor>
void RRdynamics<T, Descriptor>::decomposeOrder0(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    rawData[0] = rhoBar;
    j.to_cArray(&rawData[1]);

    Array<T, Descriptor<T>::q> eq, fEq, RReq;
    T rho = rhoBar + 1.0;
    // Array<T, 3> u(j[0]/rho, j[1]/rho, j[2]/rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    comprehensiveDynamicsTemplates<T, Descriptor>::RRcomputeEquilibriumMoments(u, RReq);
    comprehensiveDynamicsTemplates<T, Descriptor>::RRcomputeEquilibrium(rho, RReq, eq);
    for (int i = 0; i < Descriptor<T>::q; ++i) {
        fEq[i] = eq[i] - Descriptor<T>::t[i];
    }

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        rawData[1 + Descriptor<T>::d + iPop] = cell[iPop] - fEq[iPop];
    }

    int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        rawData[offset + iExt] = *cell.getExternal(iExt);
    }
}

template <typename T, template <typename U> class Descriptor>
void RRdynamics<T, Descriptor>::decomposeOrder1(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    PLB_ASSERT(false);
}

template <typename T, template <typename U> class Descriptor>
void RRdynamics<T, Descriptor>::recomposeOrder0(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    T rhoBar = rawData[0];
    Array<T, Descriptor<T>::d> j;
    j.from_cArray(&rawData[1]);

    Array<T, Descriptor<T>::q> eq, fEq, RReq;
    T rho = rhoBar + 1.0;
    // Array<T, 3> u(j[0]/rho, j[1]/rho, j[2]/rho);
    Array<T, Descriptor<T>::d> u;
    for (int i = 0; i < Descriptor<T>::d; ++i) {
        u[i] = j[i] / rho;
    }
    comprehensiveDynamicsTemplates<T, Descriptor>::RRcomputeEquilibriumMoments(u, RReq);
    comprehensiveDynamicsTemplates<T, Descriptor>::RRcomputeEquilibrium(rho, RReq, eq);
    for (int i = 0; i < Descriptor<T>::q; ++i) {
        fEq[i] = eq[i] - Descriptor<T>::t[i];
    }

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] = fEq[iPop] + rawData[1 + Descriptor<T>::d + iPop];
    }

    int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *cell.getExternal(iExt) = rawData[offset + iExt];
    }
}

template <typename T, template <typename U> class Descriptor>
void RRdynamics<T, Descriptor>::recomposeOrder1(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    PLB_ASSERT(false);
}

}  // namespace plb

#endif  // COMPREHENSIVE_ISO_THERMAL_DYNAMICS_HH
