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

#ifndef KBCDYNAMICS_HH
#define KBCDYNAMICS_HH

#include <algorithm>
#include <limits>

#include "complexDynamics/kbcDynamics.h"
#include "core/dynamicsIdentifiers.h"
#include "core/latticeStatistics.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/externalForceTemplates.h"
#include "latticeBoltzmann/kbcTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"

namespace plb {

//============================================================================//
/////////////////////////// Class KBCDynamics ////////////////////////////////
//============================================================================//
template <typename T, template <typename U> class Descriptor>
int KBCDynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, KBCDynamics<T, Descriptor> >("KBCDynamics");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
KBCDynamics<T, Descriptor>::KBCDynamics(T omega_) : IsoThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
KBCDynamics<T, Descriptor>::KBCDynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
KBCDynamics<T, Descriptor> *KBCDynamics<T, Descriptor>::clone() const
{
    return new KBCDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int KBCDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
T KBCDynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    typedef Descriptor<T> L;
    //    typedef kbcTemplates<T,Descriptor> kbc;
    //    Array<T, L::q > sEq, hEq, kEq, moments;

    //    moments[0] = L::fullRho(rhoBar);
    T invRho = L::invRho(rhoBar);
    //    moments[1] = j[0]*invRho;
    //    moments[2] = j[1]*invRho;
    //    moments[3] = 2*L::cs2;
    //    moments[4] = (T)0.;
    //    moments[5] = (T)0.;
    //    moments[6] = (T)0.;
    //    moments[7] = (T)0.;
    //    moments[8] = L::cs2*L::cs2;

    //    kbc::computeKEq(kEq, moments);
    //    kbc::computeSEq(sEq, moments);
    //    kbc::computeHEq(hEq, moments);

    //    return kEq[iPop] + sEq[iPop] + hEq[iPop] - L::SkordosFactor()*L::t[iPop];
    return dynamicsTemplatesImpl<T, Descriptor<T> >::bgk_ma2_equilibrium(
        iPop, rhoBar, invRho, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
void KBCDynamics<T, Descriptor>::collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    typedef Descriptor<T> L;
    typedef kbcTemplates<T, Descriptor> kbc;

    T rho;
    Array<T, Descriptor<T>::d> u;
    momentTemplates<T, Descriptor>::compute_rho_uLb(cell, rho, u);
    T uSqr = VectorTemplate<T, Descriptor>::normSqr(u);
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    T invRho = L::invRho(rhoBar);

    Array<T, L::q> deltaS, deltaH, fEq, sEq, sFull, moments;

    kbc::computeMomentsForS(moments, cell.getRawPopulations());

    kbc::computeSEq(sEq, moments);
    kbc::computeS(sFull, moments);

    for (plint iPop = 0; iPop < L::q; ++iPop) {
        // fEq[iPop] = kEq[iPop] + sEq[iPop] + hEq[iPop];
        fEq[iPop] = dynamicsTemplatesImpl<T, Descriptor<T> >::bgk_ma2_equilibrium(
            iPop, rhoBar, invRho, j, jSqr);
        fEq[iPop] += L::SkordosFactor() * L::t[iPop];

        deltaS[iPop] = sFull[iPop] - sEq[iPop];
        deltaH[iPop] = (cell[iPop] + L::SkordosFactor() * L::t[iPop]) - fEq[iPop] - deltaS[iPop];
    }

    T gamma = computeGamma(deltaS, deltaH, fEq);

    for (plint iPop = 0; iPop < L::q; ++iPop) {
        cell[iPop] -= this->getOmega() * 0.5 * (2 * deltaS[iPop] + gamma * deltaH[iPop]);
    }

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T KBCDynamics<T, Descriptor>::computeGamma(
    Array<T, Descriptor<T>::q> const &deltaS, Array<T, Descriptor<T>::q> const &deltaH,
    Array<T, Descriptor<T>::q> const &fEq)
{
    typedef Descriptor<T> L;

    T invBeta = 2. / this->getOmega();

    T deltaS_deltaH = T(0.);
    T deltaH_deltaH = T(0.);
    for (plint iPop = 0; iPop < L::q; ++iPop) {
        deltaS_deltaH += deltaS[iPop] * deltaH[iPop] / fEq[iPop];
        deltaH_deltaH += deltaH[iPop] * deltaH[iPop] / fEq[iPop];
    }
    if (deltaS_deltaH == (T)0.)
        return 1.;
    else
        return invBeta - (2 - invBeta) * deltaS_deltaH / deltaH_deltaH;
}

}  // namespace plb

#endif  // KBCDYNAMICS_HH
