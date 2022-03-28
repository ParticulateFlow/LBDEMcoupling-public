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
#ifndef REGULARIZED_BOUNDARY_DYNAMICS_2D_HH
#define REGULARIZED_BOUNDARY_DYNAMICS_2D_HH

#include "boundaryCondition/regularizedBoundaryDynamics2D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor, int normalX, int normalY>
RegularizedVelocityInnerCornerDynamics2D<T, Descriptor, normalX, normalY>::
    RegularizedVelocityInnerCornerDynamics2D(
        Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_) :
    BoundaryCompositeDynamics<T, Descriptor>(baseDynamics_, automaticPrepareCollision_),
    xDynamics(baseDynamics_->clone()),
    yDynamics(baseDynamics_->clone())

{ }

template <typename T, template <typename U> class Descriptor, int normalX, int normalY>
RegularizedVelocityInnerCornerDynamics2D<T, Descriptor, normalX, normalY>
    *RegularizedVelocityInnerCornerDynamics2D<T, Descriptor, normalX, normalY>::clone() const
{
    return new RegularizedVelocityInnerCornerDynamics2D<T, Descriptor, normalX, normalY>(*this);
}

template <typename T, template <typename U> class Descriptor, int normalX, int normalY>
void RegularizedVelocityInnerCornerDynamics2D<T, Descriptor, normalX, normalY>::replaceBaseDynamics(
    Dynamics<T, Descriptor> *newBaseDynamics)
{
    BoundaryCompositeDynamics<T, Descriptor>::replaceBaseDynamics(newBaseDynamics);
    xDynamics.replaceBaseDynamics(newBaseDynamics->clone());
    yDynamics.replaceBaseDynamics(newBaseDynamics->clone());
}

template <typename T, template <typename U> class Descriptor, int normalX, int normalY>
void RegularizedVelocityInnerCornerDynamics2D<T, Descriptor, normalX, normalY>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u_) const
{
    xDynamics.computeVelocity(cell, u_);
}

template <typename T, template <typename U> class Descriptor, int normalX, int normalY>
void RegularizedVelocityInnerCornerDynamics2D<T, Descriptor, normalX, normalY>::defineVelocity(
    Cell<T, Descriptor> &cell, Array<T, Descriptor<T>::d> const &u_)
{
    xDynamics.defineVelocity(cell, u_);
    yDynamics.defineVelocity(cell, u_);
}

template <typename T, template <typename U> class Descriptor, int normalX, int normalY>
T RegularizedVelocityInnerCornerDynamics2D<T, Descriptor, normalX, normalY>::computeDensity(
    Cell<T, Descriptor> const &cell) const
{
    return (xDynamics.computeDensity(cell) + yDynamics.computeDensity(cell)) / (T)2;
}

template <typename T, template <typename U> class Descriptor, int normalX, int normalY>
T RegularizedVelocityInnerCornerDynamics2D<T, Descriptor, normalX, normalY>::computeRhoBar(
    Cell<T, Descriptor> const &cell) const
{
    return (xDynamics.computeRhoBar(cell) + yDynamics.computeRhoBar(cell)) / (T)2;
}

template <typename T, template <typename U> class Descriptor, int normalX, int normalY>
void RegularizedVelocityInnerCornerDynamics2D<T, Descriptor, normalX, normalY>::computeRhoBarJ(
    Cell<T, Descriptor> const &cell, T &rhoBar_, Array<T, Descriptor<T>::d> &j_) const
{
    rhoBar_ = this->computeRhoBar(cell);
    T rho = Descriptor<T>::fullRho(rhoBar_);
    this->computeVelocity(cell, j_);
    if (!this->velIsJ()) {
        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            j_[iD] *= rho;
        }
    }
}

template <typename T, template <typename U> class Descriptor, int normalX, int normalY>
void RegularizedVelocityInnerCornerDynamics2D<T, Descriptor, normalX, normalY>::computeRhoBarJPiNeq(
    Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
    Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const
{
    T tmpRhoBar;
    Array<T, Descriptor<T>::d> tmpJ;
    Array<T, SymmetricTensor<T, Descriptor>::n> tmpPiNeq;
    xDynamics.computeRhoBarJPiNeq(cell, rhoBar, j, PiNeq);
    yDynamics.computeRhoBarJPiNeq(cell, tmpRhoBar, tmpJ, tmpPiNeq);

    rhoBar = (rhoBar + tmpRhoBar) / (T)2;
    for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
        j[iD] = (j[iD] + tmpJ[iD]) / (T)2;
    }
    for (int iPi = 0; iPi < Descriptor<T>::d; ++iPi) {
        PiNeq[iPi] = (PiNeq[iPi] + tmpPiNeq[iPi]) / (T)2;
    }
}

template <typename T, template <typename U> class Descriptor, int normalX, int normalY>
void RegularizedVelocityInnerCornerDynamics2D<T, Descriptor, normalX, normalY>::completePopulations(
    Cell<T, Descriptor> &cell) const
{
    // 1. Assign "bounce-back of off-equilibrium" to the unknown population
    Array<int, Descriptor<T>::d> v(-normalX, -normalY);
    plint unknownF = indexTemplates::findVelocity<Descriptor<T> >(v);

    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    this->computeRhoBarJ(cell, rhoBar, j);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    // Do nothing if there is no unknown population
    if (unknownF != Descriptor<T>::q) {
        plint oppositeF = indexTemplates::opposite<Descriptor<T> >(unknownF);
        cell[unknownF] = cell[oppositeF] - this->computeEquilibrium(oppositeF, rhoBar, j, jSqr)
                         + this->computeEquilibrium(unknownF, rhoBar, j, jSqr);
    }

    // 2. Regularize all populations
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    this->getBaseDynamics().computePiNeq(cell, PiNeq);
    ;
    this->getBaseDynamics().regularize(cell, rhoBar, j, jSqr, PiNeq);
}

}  // namespace plb

#endif  // REGULARIZED_BOUNDARY_DYNAMICS_2D_HH
