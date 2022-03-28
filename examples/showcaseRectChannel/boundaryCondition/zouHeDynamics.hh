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

#ifndef ZOU_HE_DYNAMICS_HH
#define ZOU_HE_DYNAMICS_HH

#include <cmath>

#include "boundaryCondition/zouHeDynamics.h"
#include "core/dynamicsIdentifiers.h"
#include "core/latticeStatistics.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/indexTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"

namespace plb {

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void ZouHeClosure(Cell<T, Descriptor> &cell, Dynamics<T, Descriptor> const &dynamics)
{
    typedef momentTemplates<T, Descriptor> mtl;
    typedef Descriptor<T> L;

    // Along all the commented parts of this code there will be an example based
    // on the situation where the wall's normal vector if (0,1) and the
    // numbering of the velocites are done according to the D2Q9
    // lattice of the Palabos library.

    // Find all the missing populations
    // (directions 3,4,5)
    std::vector<plint> missingIndexes =
        indexTemplates::subIndexOutgoing<L, direction, orientation>();

    // Will contain the missing populations that are not normal to the wall.
    // (directions 3,5)
    std::vector<plint> missingDiagonalIndexes = missingIndexes;
    for (pluint iPop = 0; iPop < missingIndexes.size(); ++iPop) {
        plint numOfNonNullComp = 0;
        for (int iDim = 0; iDim < L::d; ++iDim)
            numOfNonNullComp += abs(L::c[missingIndexes[iPop]][iDim]);

        if (numOfNonNullComp == 1) {
            missingDiagonalIndexes.erase(missingDiagonalIndexes.begin() + iPop);
            break;
        }
    }

    T rhoBar;
    Array<T, L::d> j;
    dynamics.computeRhoBarJ(cell, rhoBar, j);

    T falseRho;
    Array<T, L::d> falseU;

    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    // The unknown non equilibrium populations are bounced back
    // (f[3] = feq[3] + fneq[7], f[4] = feq[4] + fneq[8],
    //  f[5] = feq[5] + fneq[1])
    for (pluint iPop = 0; iPop < missingIndexes.size(); ++iPop) {
        plint pop = missingIndexes[iPop];
        plint opp = indexTemplates::opposite<L>(pop);
        cell[pop] = cell[opp] - dynamics.computeEquilibrium(opp, rhoBar, j, jSqr)
                    + dynamics.computeEquilibrium(pop, rhoBar, j, jSqr);
    }

    // We recompute rho and u in order to have the new momentum and density. Since
    // the momentum is not conserved from this scheme, we will corect it. By adding
    // a contribution to the missingDiagonalVelocities.
    mtl::compute_rho_uLb(cell, falseRho, falseU);

    Array<T, L::d> diff;
    for (int iDim = 0; iDim < L::d; ++iDim) {
        diff[iDim] = (j[iDim] - falseRho * falseU[iDim]) / (T)missingDiagonalIndexes.size();
    }

    for (pluint iPop = 0; iPop < missingDiagonalIndexes.size(); ++iPop) {
        for (int iDim = 1; iDim < L::d; ++iDim) {
            cell[missingDiagonalIndexes[iPop]] +=
                L::c[missingDiagonalIndexes[iPop]][(direction + iDim) % L::d]
                * diff[(direction + iDim) % L::d];
        }
    }
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int ZouHeVelocityDynamics<T, Descriptor, direction, orientation>::id =
    meta::registerGeneralDynamics<
        T, Descriptor, ZouHeVelocityDynamics<T, Descriptor, direction, orientation> >(
        std::string("Boundary_ZouHeVelocity_") + util::val2str(direction) + std::string("_")
        + util::val2str(orientation));

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
ZouHeVelocityDynamics<T, Descriptor, direction, orientation>::ZouHeVelocityDynamics(
    Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision) :
    VelocityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics, automaticPrepareCollision)
{ }

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
ZouHeVelocityDynamics<T, Descriptor, direction, orientation>::ZouHeVelocityDynamics(
    HierarchicUnserializer &unserializer) :
    VelocityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
ZouHeVelocityDynamics<T, Descriptor, direction, orientation>
    *ZouHeVelocityDynamics<T, Descriptor, direction, orientation>::clone() const
{
    return new ZouHeVelocityDynamics<T, Descriptor, direction, orientation>(*this);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int ZouHeVelocityDynamics<T, Descriptor, direction, orientation>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void ZouHeVelocityDynamics<T, Descriptor, direction, orientation>::completePopulations(
    Cell<T, Descriptor> &cell) const
{
    ZouHeClosure<T, Descriptor, direction, orientation>(cell, *this);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int ZouHePressureDynamics<T, Descriptor, direction, orientation>::id =
    meta::registerGeneralDynamics<
        T, Descriptor, ZouHePressureDynamics<T, Descriptor, direction, orientation> >(
        std::string("Boundary_ZouHePressure_") + util::val2str(direction) + std::string("_")
        + util::val2str(orientation));

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
ZouHePressureDynamics<T, Descriptor, direction, orientation>::ZouHePressureDynamics(
    Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision) :
    DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics, automaticPrepareCollision)
{ }

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
ZouHePressureDynamics<T, Descriptor, direction, orientation>::ZouHePressureDynamics(
    HierarchicUnserializer &unserializer) :
    DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
ZouHePressureDynamics<T, Descriptor, direction, orientation>
    *ZouHePressureDynamics<T, Descriptor, direction, orientation>::clone() const
{
    return new ZouHePressureDynamics<T, Descriptor, direction, orientation>(*this);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int ZouHePressureDynamics<T, Descriptor, direction, orientation>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void ZouHePressureDynamics<T, Descriptor, direction, orientation>::completePopulations(
    Cell<T, Descriptor> &cell) const
{
    ZouHeClosure<T, Descriptor, direction, orientation>(cell, *this);
}

}  // namespace plb

#endif  // ZOU_HE_DYNAMICS_HH
