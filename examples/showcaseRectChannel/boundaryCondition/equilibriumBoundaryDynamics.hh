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
 * Dirichlet boundary condition which imposes equilibrium (but computes
 * density properly from velocity, or vice versa)
 */
#ifndef EQUILIBRIUM_BOUNDARY_DYNAMICS_HH
#define EQUILIBRIUM_BOUNDARY_DYNAMICS_HH

#include "boundaryCondition/equilibriumBoundaryDynamics.h"
#include "core/cell.h"
#include "core/dynamicsIdentifiers.h"

namespace plb {

/* *************** Class EquilibriumVelocityBoundaryDynamics ************* */

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int EquilibriumVelocityBoundaryDynamics<T, Descriptor, direction, orientation>::id =
    meta::registerCompositeDynamics<
        T, Descriptor, EquilibriumVelocityBoundaryDynamics<T, Descriptor, direction, orientation> >(
        std::string("Boundary_EquilibriumVelocity_") + util::val2str(direction) + std::string("_")
        + util::val2str(orientation));

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
EquilibriumVelocityBoundaryDynamics<T, Descriptor, direction, orientation>::
    EquilibriumVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision) :
    VelocityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics_, automaticPrepareCollision)
{ }

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
EquilibriumVelocityBoundaryDynamics<T, Descriptor, direction, orientation>
    *EquilibriumVelocityBoundaryDynamics<T, Descriptor, direction, orientation>::clone() const
{
    return new EquilibriumVelocityBoundaryDynamics<T, Descriptor, direction, orientation>(*this);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int EquilibriumVelocityBoundaryDynamics<T, Descriptor, direction, orientation>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void EquilibriumVelocityBoundaryDynamics<
    T, Descriptor, direction, orientation>::completePopulations(Cell<T, Descriptor> &cell) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    this->computeRhoBarJ(cell, rhoBar, j);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] = this->getBaseDynamics().computeEquilibrium(iPop, rhoBar, j, jSqr);
    }
}

/* *************** Class EquilibriumDensityBoundaryDynamics ************* */

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int EquilibriumDensityBoundaryDynamics<T, Descriptor, direction, orientation>::id =
    meta::registerCompositeDynamics<
        T, Descriptor, EquilibriumDensityBoundaryDynamics<T, Descriptor, direction, orientation> >(
        std::string("Boundary_EquilibriumDensity_") + util::val2str(direction) + std::string("_")
        + util::val2str(orientation));

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
EquilibriumDensityBoundaryDynamics<T, Descriptor, direction, orientation>::
    EquilibriumDensityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision) :
    DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics_, automaticPrepareCollision)
{ }

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
EquilibriumDensityBoundaryDynamics<T, Descriptor, direction, orientation>
    *EquilibriumDensityBoundaryDynamics<T, Descriptor, direction, orientation>::clone() const
{
    return new EquilibriumDensityBoundaryDynamics<T, Descriptor, direction, orientation>(*this);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int EquilibriumDensityBoundaryDynamics<T, Descriptor, direction, orientation>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void EquilibriumDensityBoundaryDynamics<T, Descriptor, direction, orientation>::completePopulations(
    Cell<T, Descriptor> &cell) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    this->computeRhoBarJ(cell, rhoBar, j);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] = this->getBaseDynamics().computeEquilibrium(iPop, rhoBar, j, jSqr);
    }
}

/* *************** Class EquilibriumDensityAndVelocityBoundaryDynamics ************* */

template <typename T, template <typename U> class Descriptor>
int EquilibriumDensityAndVelocityBoundaryDynamics<T, Descriptor>::id =
    meta::registerCompositeDynamics<
        T, Descriptor, EquilibriumDensityAndVelocityBoundaryDynamics<T, Descriptor> >(
        "Boundary_EquilibriumDensityAndVelocity");

template <typename T, template <typename U> class Descriptor>
EquilibriumDensityAndVelocityBoundaryDynamics<T, Descriptor>::
    EquilibriumDensityAndVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision) :
    StoreDensityAndVelocityDynamics<T, Descriptor>(baseDynamics_, automaticPrepareCollision)
{ }

template <typename T, template <typename U> class Descriptor>
EquilibriumDensityAndVelocityBoundaryDynamics<T, Descriptor>
    *EquilibriumDensityAndVelocityBoundaryDynamics<T, Descriptor>::clone() const
{
    return new EquilibriumDensityAndVelocityBoundaryDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int EquilibriumDensityAndVelocityBoundaryDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void EquilibriumDensityAndVelocityBoundaryDynamics<T, Descriptor>::completePopulations(
    Cell<T, Descriptor> &cell) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    this->computeRhoBarJ(cell, rhoBar, j);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] = this->getBaseDynamics().computeEquilibrium(iPop, rhoBar, j, jSqr);
    }
}

}  // namespace plb

#endif  // EQUILIBRIUM_BOUNDARY_DYNAMICS_HH
