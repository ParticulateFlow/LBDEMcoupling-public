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

#ifndef CARREAU_DYNAMICS_HH
#define CARREAU_DYNAMICS_HH

#include "complexDynamics/carreauDynamics.h"
#include "complexDynamics/carreauDynamicsTemplates.h"
#include "complexDynamics/carreauGlobalDefs.h"
#include "core/dynamicsIdentifiers.h"
#include "core/util.h"
#include "latticeBoltzmann/momentTemplates.h"

namespace plb {

template <typename T, template <typename U> class Descriptor, int N>
int CarreauDynamics<T, Descriptor, N>::id =
    meta::registerCompositeDynamics<T, Descriptor, CarreauDynamics<T, Descriptor, N> >(
        std::string("CarreauDynamics_") + util::val2str(N));

template <typename T, template <typename U> class Descriptor, int N>
CarreauDynamics<T, Descriptor, N>::CarreauDynamics(
    Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision) :
    OmegaFromPiDynamics<T, Descriptor>(baseDynamics_, automaticPrepareCollision)
{ }

template <typename T, template <typename U> class Descriptor, int N>
T CarreauDynamics<T, Descriptor, N>::getOmegaFromPiAndRhoBar(
    Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T rhoBar) const
{
    T nu0_nuInfoverCs2 =
        (global::CarreauParameters().getNu0() - global::CarreauParameters().getNuInf())
        * Descriptor<T>::invCs2;
    T nuInfoverCs2 = global::CarreauParameters().getNuInf() * Descriptor<T>::invCs2;
    T nMinusOneOverTwo = (global::CarreauParameters().getExponent() - (T)1) / (T)2;
    T lambdaOverCs2sqr = global::CarreauParameters().getLambda() * Descriptor<T>::invCs2;
    lambdaOverCs2sqr *= lambdaOverCs2sqr;

    T piNeqNormSqr = SymmetricTensor<T, Descriptor>::tensorNormSqr(PiNeq);
    T alpha = lambdaOverCs2sqr * piNeqNormSqr * (T)0.5 * Descriptor<T>::invRho(rhoBar)
              * Descriptor<T>::invRho(rhoBar);

    T omega = carreauDynamicsTemplates<T, N>::fromPiAndRhoToOmega(
        alpha, nu0_nuInfoverCs2, nuInfoverCs2, nMinusOneOverTwo, this->getOmega());

    return omega;
}

template <typename T, template <typename U> class Descriptor, int N>
CarreauDynamics<T, Descriptor, N> *CarreauDynamics<T, Descriptor, N>::clone() const
{
    return new CarreauDynamics<T, Descriptor, N>(*this);
}

template <typename T, template <typename U> class Descriptor, int N>
int CarreauDynamics<T, Descriptor, N>::getId() const
{
    return id;
}

/* *************** Class BGKCarreauDynamics ************************************ */

template <typename T, template <typename U> class Descriptor, int N>
int BGKCarreauDynamics<T, Descriptor, N>::id =
    meta::registerGeneralDynamics<T, Descriptor, BGKCarreauDynamics<T, Descriptor, N> >(
        std::string("CarreauDynamics_BGK_") + util::val2str(N));

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor, int N>
BGKCarreauDynamics<T, Descriptor, N>::BGKCarreauDynamics(T omega) :
    IsoThermalBulkDynamics<T, Descriptor>(omega)
{ }

template <typename T, template <typename U> class Descriptor, int N>
BGKCarreauDynamics<T, Descriptor, N>::BGKCarreauDynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor, int N>
BGKCarreauDynamics<T, Descriptor, N> *BGKCarreauDynamics<T, Descriptor, N>::clone() const
{
    return new BGKCarreauDynamics<T, Descriptor, N>(*this);
}

template <typename T, template <typename U> class Descriptor, int N>
int BGKCarreauDynamics<T, Descriptor, N>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor, int N>
void BGKCarreauDynamics<T, Descriptor, N>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T nu0_nuInfoverCs2 =
        (global::CarreauParameters().getNu0() - global::CarreauParameters().getNuInf())
        * Descriptor<T>::invCs2;
    T nuInfoverCs2 = global::CarreauParameters().getNuInf() * Descriptor<T>::invCs2;
    T nMinusOneOverTwo = (global::CarreauParameters().getExponent() - (T)1) / (T)2;
    T lambdaOverCs2sqr = global::CarreauParameters().getLambda() * Descriptor<T>::invCs2;
    lambdaOverCs2sqr *= lambdaOverCs2sqr;

    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);

    T piNeqNormSqr = SymmetricTensor<T, Descriptor>::tensorNormSqr(PiNeq);
    T alpha = lambdaOverCs2sqr * piNeqNormSqr * (T)0.5 * Descriptor<T>::invRho(rhoBar)
              * Descriptor<T>::invRho(rhoBar);

    T omegaTmp = (T)1;
    T omega0 = this->getOmega();

    for (int iN = 0; iN < N; ++iN) {
        omegaTmp = carreauDynamicsTemplates<T, 0>::fromPiAndRhoToOmega(
            alpha, nu0_nuInfoverCs2, nuInfoverCs2, nMinusOneOverTwo, omega0);
        if (std::fabs((omegaTmp - omega0) / omega0) < 1.0e-10) {
            omega0 = omegaTmp;
            break;
        }
        omega0 = omegaTmp;
    }

    this->setOmega(omega0);

    T uSqr = dynamicsTemplates<T, Descriptor>::bgk_ma2_collision(cell, rhoBar, j, this->getOmega());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor, int N>
T BGKCarreauDynamics<T, Descriptor, N>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class RegularizedBGKCarreauDynamics ************************************ */

template <typename T, template <typename U> class Descriptor, int N>
int RegularizedBGKCarreauDynamics<T, Descriptor, N>::id =
    meta::registerGeneralDynamics<T, Descriptor, RegularizedBGKCarreauDynamics<T, Descriptor, N> >(
        std::string("CarreauDynamics_RLB_") + util::val2str(N));

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor, int N>
RegularizedBGKCarreauDynamics<T, Descriptor, N>::RegularizedBGKCarreauDynamics(T omega) :
    IsoThermalBulkDynamics<T, Descriptor>(omega)
{ }

template <typename T, template <typename U> class Descriptor, int N>
RegularizedBGKCarreauDynamics<T, Descriptor, N>::RegularizedBGKCarreauDynamics(
    HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor, int N>
RegularizedBGKCarreauDynamics<T, Descriptor, N>
    *RegularizedBGKCarreauDynamics<T, Descriptor, N>::clone() const
{
    return new RegularizedBGKCarreauDynamics<T, Descriptor, N>(*this);
}

template <typename T, template <typename U> class Descriptor, int N>
int RegularizedBGKCarreauDynamics<T, Descriptor, N>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor, int N>
void RegularizedBGKCarreauDynamics<T, Descriptor, N>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T nu0overCs2 = global::CarreauParameters().getNu0() * Descriptor<T>::invCs2;
    T nMinusOneOverTwo = (global::CarreauParameters().getExponent() - (T)1) / (T)2;
    T lambdaOverCs2sqr = global::CarreauParameters().getLambda() * Descriptor<T>::invCs2;
    lambdaOverCs2sqr *= lambdaOverCs2sqr;

    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);

    T piNeqNormSqr = SymmetricTensor<T, Descriptor>::tensorNormSqr(PiNeq);
    T alpha = lambdaOverCs2sqr * piNeqNormSqr * (T)0.5 * Descriptor<T>::invRho(rhoBar)
              * Descriptor<T>::invRho(rhoBar);

    T omegaTmp = (T)1;
    T omega0 = this->getOmega();

    for (int iN = 0; iN < N; ++iN) {
        omegaTmp = carreauDynamicsTemplates<T, 0>::fromPiAndRhoToOmega(
            alpha, nu0overCs2, nMinusOneOverTwo, omega0);
        if (std::fabs((omegaTmp - omega0) / omega0) < 1.0e-3) {
            break;
        }
        omega0 = omegaTmp;
    }

    this->setOmega(omegaTmp);

    T uSqr =
        dynamicsTemplates<T, Descriptor>::rlb_collision(cell, rhoBar, j, PiNeq, this->getOmega());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor, int N>
T RegularizedBGKCarreauDynamics<T, Descriptor, N>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

}  // namespace plb

#endif  // VARIABLE_OMEGA_DYNAMICS_HH
