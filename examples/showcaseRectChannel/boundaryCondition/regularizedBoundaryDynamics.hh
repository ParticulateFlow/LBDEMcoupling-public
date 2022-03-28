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
#ifndef REGULARIZED_BOUNDARY_DYNAMICS_HH
#define REGULARIZED_BOUNDARY_DYNAMICS_HH

#include "boundaryCondition/boundaryTemplates.h"
#include "boundaryCondition/regularizedBoundaryDynamics.h"
#include "core/cell.h"
#include "core/dynamicsIdentifiers.h"
#include "latticeBoltzmann/indexTemplates.h"

namespace plb {

/* *************** Class RegularizedVelocityBoundaryDynamics ************* */

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int RegularizedVelocityBoundaryDynamics<T, Descriptor, direction, orientation>::id =
    meta::registerGeneralDynamics<
        T, Descriptor, RegularizedVelocityBoundaryDynamics<T, Descriptor, direction, orientation> >(
        std::string("Boundary_RegularizedVelocity_") + util::val2str(direction) + std::string("_")
        + util::val2str(orientation));

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
RegularizedVelocityBoundaryDynamics<T, Descriptor, direction, orientation>::
    RegularizedVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_) :
    VelocityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics_, automaticPrepareCollision_)
{ }

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
RegularizedVelocityBoundaryDynamics<T, Descriptor, direction, orientation>::
    RegularizedVelocityBoundaryDynamics(HierarchicUnserializer &unserializer) :
    VelocityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
RegularizedVelocityBoundaryDynamics<T, Descriptor, direction, orientation>
    *RegularizedVelocityBoundaryDynamics<T, Descriptor, direction, orientation>::clone() const
{
    return new RegularizedVelocityBoundaryDynamics<T, Descriptor, direction, orientation>(*this);
}
template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int RegularizedVelocityBoundaryDynamics<T, Descriptor, direction, orientation>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void RegularizedVelocityBoundaryDynamics<
    T, Descriptor, direction, orientation>::completePopulations(Cell<T, Descriptor> &cell) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    this->computeRhoBarJ(cell, rhoBar, j);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    boundaryTemplates<T, Descriptor, direction, orientation>::compute_PiNeq(
        this->getBaseDynamics(), cell, rhoBar, j, jSqr, PiNeq);

    this->getBaseDynamics().regularize(cell, rhoBar, j, jSqr, PiNeq);
}

/* *************** Class RegularizedVelocityConstRhoBoundaryDynamics ************* */

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int RegularizedVelocityConstRhoBoundaryDynamics<T, Descriptor, direction, orientation>::id =
    meta::registerGeneralDynamics<
        T, Descriptor,
        RegularizedVelocityConstRhoBoundaryDynamics<T, Descriptor, direction, orientation> >(
        std::string("Boundary_RegularizedVelocityConstRho_") + util::val2str(direction)
        + std::string("_") + util::val2str(orientation));

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
RegularizedVelocityConstRhoBoundaryDynamics<T, Descriptor, direction, orientation>::
    RegularizedVelocityConstRhoBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_) :
    VelocityDirichletConstRhoBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics_, automaticPrepareCollision_)
{ }

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
RegularizedVelocityConstRhoBoundaryDynamics<T, Descriptor, direction, orientation>::
    RegularizedVelocityConstRhoBoundaryDynamics(HierarchicUnserializer &unserializer) :
    VelocityDirichletConstRhoBoundaryDynamics<T, Descriptor, direction, orientation>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
RegularizedVelocityConstRhoBoundaryDynamics<T, Descriptor, direction, orientation>
    *RegularizedVelocityConstRhoBoundaryDynamics<T, Descriptor, direction, orientation>::clone()
        const
{
    return new RegularizedVelocityConstRhoBoundaryDynamics<T, Descriptor, direction, orientation>(
        *this);
}
template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int RegularizedVelocityConstRhoBoundaryDynamics<T, Descriptor, direction, orientation>::getId()
    const
{
    return id;
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void RegularizedVelocityConstRhoBoundaryDynamics<
    T, Descriptor, direction, orientation>::completePopulations(Cell<T, Descriptor> &cell) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    this->computeRhoBarJ(cell, rhoBar, j);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    boundaryTemplates<T, Descriptor, direction, orientation>::compute_PiNeq(
        this->getBaseDynamics(), cell, rhoBar, j, jSqr, PiNeq);

    this->getBaseDynamics().regularize(cell, rhoBar, j, jSqr, PiNeq);

    T newRhoBar = Descriptor<T>::rhoBar((T)1.0);
    Array<T, 3> newJ(j * (Descriptor<T>::fullRho(newRhoBar) / Descriptor<T>::fullRho(rhoBar)));
    T newJsqr = VectorTemplate<T, Descriptor>::normSqr(j);

    for (pluint iPop = 0; iPop < Descriptor<T>::d; ++iPop) {
        cell[iPop] += this->getBaseDynamics().computeEquilibrium(iPop, newRhoBar, newJ, newJsqr)
                      - this->getBaseDynamics().computeEquilibrium(iPop, rhoBar, j, jSqr);
    }
}

/* *************** Class RegularizedDensityBoundaryDynamics ************* */

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int RegularizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>::id =
    meta::registerGeneralDynamics<
        T, Descriptor, RegularizedDensityBoundaryDynamics<T, Descriptor, direction, orientation> >(
        std::string("Boundary_RegularizedDensity_") + util::val2str(direction) + std::string("_")
        + util::val2str(orientation));

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
RegularizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>::
    RegularizedDensityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_) :
    DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics_, automaticPrepareCollision_)
{ }

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
RegularizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>::
    RegularizedDensityBoundaryDynamics(HierarchicUnserializer &unserializer) :
    DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
RegularizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>
    *RegularizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>::clone() const
{
    return new RegularizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>(*this);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int RegularizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void RegularizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>::completePopulations(
    Cell<T, Descriptor> &cell) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    this->computeRhoBarJ(cell, rhoBar, j);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    boundaryTemplates<T, Descriptor, direction, orientation>::compute_PiNeq(
        this->getBaseDynamics(), cell, rhoBar, j, jSqr, PiNeq);

    this->getBaseDynamics().regularize(cell, rhoBar, j, jSqr, PiNeq);
}

}  // namespace plb

#endif  // REGULARIZED_BOUNDARY_DYNAMICS_HH
