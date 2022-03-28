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

/* Orestis Malaspinas contributed this code.
 */

#ifndef INAMURO_ANALYTICAL_DYNAMICS_HH
#define INAMURO_ANALYTICAL_DYNAMICS_HH

#include <cmath>

#include "boundaryCondition/inamuroAnalyticalDynamics.h"
#include "core/cell.h"
#include "core/dynamicsIdentifiers.h"
#include "core/latticeStatistics.h"
#include "core/util.h"
#include "latticeBoltzmann/indexTemplates.h"

namespace plb {

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void InamuroAnalyticalClosure(Cell<T, Descriptor> &cell, Dynamics<T, Descriptor> const &dynamics)
{
    typedef Descriptor<T> L;

    // Along all the commented parts of this code there will be an example based
    // on the situation where the wall's normal vector if (0,1) and the
    // numerotation of the velocites are done according to the D2Q9
    // lattice of the Palabos library.

    // Find all the missing populations
    // (directions 3,4,5)
    std::vector<plint> missInd = indexTemplates::subIndexOutgoing<L, direction, orientation>();

    // Will contain the missing poputations that are not normal to the wall.
    // (directions 3,5)
    std::vector<plint> missDiagInd = missInd;

    for (pluint iPop = 0; iPop < missInd.size(); ++iPop) {
        plint numOfNonNullComp = 0;
        for (int iDim = 0; iDim < L::d; ++iDim)
            numOfNonNullComp += abs(L::c[missInd[iPop]][iDim]);

        if (numOfNonNullComp == 1) {
            missDiagInd.erase(missDiagInd.begin() + iPop);
            break;
        }
    }

    // Will contain the populations normal to the wall's normal vector.
    // (directions 2,6)
    std::vector<plint> perpInd = indexTemplates::subIndex<L, direction, 0>();
    for (pluint iPop = 0; iPop < perpInd.size(); ++iPop) {
        if (L::c[perpInd[iPop]][0] == 0 && L::c[perpInd[iPop]][1] == 0) {
            perpInd.erase(perpInd.begin() + iPop);
            break;
        }
    }

    T rho = dynamics.computeDensity(cell);
    Array<T, L::d> u;
    dynamics.computeVelocity(cell, u);

    T rhoCs = T();
    Array<T, L::d> jCs;
    for (int iDim = 0; iDim < L::d; ++iDim)
        jCs[iDim] = T();

    T fSum = T();
    for (pluint iPop = 0; iPop < missInd.size(); ++iPop) {
        fSum += cell[indexTemplates::opposite<L>(missInd[iPop])];
    }
    // do not forget the "+1" in the rhoCs equation in the numerator (it's
    // here because fEq = usualfEq - t[i])
    rhoCs = ((T)6 * (-orientation * rho * u[direction] + fSum) + (T)1)
            / ((T)3 * u[direction] * u[direction] - orientation * (T)3 * u[direction] + (T)1);

    T fDiffPerp = T();
    for (pluint iPop = 0; iPop < perpInd.size(); ++iPop)
        fDiffPerp += L::c[perpInd[iPop]][(direction + 1) % 2] * cell[perpInd[iPop]];
    fDiffPerp *= orientation;

    T fDiffDiag = T();
    for (pluint iPop = 0; iPop < missDiagInd.size(); ++iPop)
        fDiffDiag += L::c[indexTemplates::opposite<L>(missDiagInd[iPop])][(direction + 1) % 2]
                     * cell[indexTemplates::opposite<L>(missDiagInd[iPop])];
    fDiffDiag *= orientation;

    jCs[(direction + 1) % L::d] =
        (-orientation * (T)6 * rho * u[(direction + 1) % L::d]
         + orientation * rhoCs * u[(direction + 1) % L::d]
         - (T)3 * rhoCs * u[direction] * u[(direction + 1) % L::d] + (T)6 * (fDiffPerp + fDiffDiag))
        / (-orientation + (T)3 * u[direction]);

    for (int iDim = 0; iDim < L::d; ++iDim)
        jCs[iDim] += rhoCs * u[iDim];

    T jSqr = VectorTemplate<T, Descriptor>::normSqr(jCs);

    for (pluint iPop = 0; iPop < missInd.size(); ++iPop) {
        cell[missInd[iPop]] =
            dynamics.computeEquilibrium(missInd[iPop], Descriptor<T>::rhoBar(rhoCs), jCs, jSqr);
    }
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int InamuroAnalyticalVelocityDynamics<T, Descriptor, direction, orientation>::id =
    meta::registerCompositeDynamics<
        T, Descriptor, InamuroAnalyticalVelocityDynamics<T, Descriptor, direction, orientation> >(
        std::string("Boundary_InamuroAnalytical_") + util::val2str(direction) + std::string("_")
        + util::val2str(orientation));

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
InamuroAnalyticalVelocityDynamics<T, Descriptor, direction, orientation>::
    InamuroAnalyticalVelocityDynamics(
        Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision) :
    VelocityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics, automaticPrepareCollision)
{ }

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
InamuroAnalyticalVelocityDynamics<T, Descriptor, direction, orientation>
    *InamuroAnalyticalVelocityDynamics<T, Descriptor, direction, orientation>::clone() const
{
    return new InamuroAnalyticalVelocityDynamics<T, Descriptor, direction, orientation>(*this);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int InamuroAnalyticalVelocityDynamics<T, Descriptor, direction, orientation>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void InamuroAnalyticalVelocityDynamics<T, Descriptor, direction, orientation>::completePopulations(
    Cell<T, Descriptor> &cell) const
{
    InamuroAnalyticalClosure<T, Descriptor, direction, orientation>(cell, *this);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int InamuroAnalyticalPressureDynamics<T, Descriptor, direction, orientation>::id =
    meta::registerCompositeDynamics<
        T, Descriptor, InamuroAnalyticalPressureDynamics<T, Descriptor, direction, orientation> >(
        std::string("Boundary_InamuroAnalyticalPressure_") + util::val2str(direction)
        + std::string("_") + util::val2str(orientation));

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
InamuroAnalyticalPressureDynamics<T, Descriptor, direction, orientation>::
    InamuroAnalyticalPressureDynamics(
        Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision) :
    DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics, automaticPrepareCollision)
{ }

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
InamuroAnalyticalPressureDynamics<T, Descriptor, direction, orientation>
    *InamuroAnalyticalPressureDynamics<T, Descriptor, direction, orientation>::clone() const
{
    return new InamuroAnalyticalPressureDynamics<T, Descriptor, direction, orientation>(*this);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int InamuroAnalyticalPressureDynamics<T, Descriptor, direction, orientation>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void InamuroAnalyticalPressureDynamics<T, Descriptor, direction, orientation>::completePopulations(
    Cell<T, Descriptor> &cell) const
{
    InamuroAnalyticalClosure<T, Descriptor, direction, orientation>(cell, *this);
}

}  // namespace plb

#endif  // INAMURO_ANALYTICAL_DYNAMICS_HH
