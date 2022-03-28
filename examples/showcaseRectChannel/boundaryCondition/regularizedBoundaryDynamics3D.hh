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
#ifndef REGULARIZED_BOUNDARY_DYNAMICS_3D_HH
#define REGULARIZED_BOUNDARY_DYNAMICS_3D_HH

#include "boundaryCondition/regularizedBoundaryDynamics3D.h"
#include "core/cell.h"
#include "latticeBoltzmann/indexTemplates.h"

namespace plb {

/* *************** Class RegularizedVelocityInnerEdgeDynamics3D ****** */

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
RegularizedVelocityInnerEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>::
    RegularizedVelocityInnerEdgeDynamics3D(
        Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_) :
    BoundaryCompositeDynamics<T, Descriptor>(baseDynamics_, automaticPrepareCollision_),
    dynamics1(baseDynamics_->clone()),
    dynamics2(baseDynamics_->clone())

{ }

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
RegularizedVelocityInnerEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>
    *RegularizedVelocityInnerEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>::clone() const
{
    return new RegularizedVelocityInnerEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>(
        *this);
}

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
void RegularizedVelocityInnerEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>::
    replaceBaseDynamics(Dynamics<T, Descriptor> *newBaseDynamics)
{
    BoundaryCompositeDynamics<T, Descriptor>::replaceBaseDynamics(newBaseDynamics);
    dynamics1.replaceBaseDynamics(newBaseDynamics->clone());
    dynamics2.replaceBaseDynamics(newBaseDynamics->clone());
}

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
void RegularizedVelocityInnerEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>::
    computeVelocity(Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u_) const
{
    dynamics1.computeVelocity(cell, u_);
}

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
void RegularizedVelocityInnerEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>::defineVelocity(
    Cell<T, Descriptor> &cell, Array<T, Descriptor<T>::d> const &u_)
{
    dynamics1.defineVelocity(cell, u_);
    dynamics2.defineVelocity(cell, u_);
}

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
T RegularizedVelocityInnerEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>::computeDensity(
    Cell<T, Descriptor> const &cell) const
{
    return (dynamics1.computeDensity(cell) + dynamics2.computeDensity(cell)) / (T)2;
}

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
T RegularizedVelocityInnerEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>::computeRhoBar(
    Cell<T, Descriptor> const &cell) const
{
    return (dynamics1.computeRhoBar(cell) + dynamics2.computeRhoBar(cell)) / (T)2;
}

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
void RegularizedVelocityInnerEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>::computeRhoBarJ(
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

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
void RegularizedVelocityInnerEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>::
    computeRhoBarJPiNeq(
        Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const
{
    T tmpRhoBar;
    Array<T, Descriptor<T>::d> tmpJ;
    Array<T, SymmetricTensor<T, Descriptor>::n> tmpPiNeq;
    dynamics1.computeRhoBarJPiNeq(cell, rhoBar, j, PiNeq);
    dynamics2.computeRhoBarJPiNeq(cell, tmpRhoBar, tmpJ, tmpPiNeq);

    rhoBar = (rhoBar + tmpRhoBar) / (T)2;
    for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
        j[iD] = (j[iD] + tmpJ[iD]) / (T)2;
    }
    for (int iPi = 0; iPi < Descriptor<T>::d; ++iPi) {
        PiNeq[iPi] = (PiNeq[iPi] + tmpPiNeq[iPi]) / (T)2;
    }
}

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
void RegularizedVelocityInnerEdgeDynamics3D<
    T, Descriptor, plane, normal1, normal2>::completePopulations(Cell<T, Descriptor> &cell) const
{
    // 1. Assign "bounce-back of off-equilibrium" to the unknown population
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    this->computeRhoBarJ(cell, rhoBar, j);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        if ((Descriptor<T>::c[iPop][direction1] == -normal1)
            && (Descriptor<T>::c[iPop][direction2] == -normal2))
        {
            plint opp = indexTemplates::opposite<Descriptor<T> >(iPop);
            cell[iPop] = cell[opp] - this->computeEquilibrium(opp, rhoBar, j, jSqr)
                         + this->computeEquilibrium(iPop, rhoBar, j, jSqr);
        }
    }

    // 2. Regularize all populations
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    this->getBaseDynamics().computePiNeq(cell, PiNeq);
    ;
    this->getBaseDynamics().regularize(cell, rhoBar, j, jSqr, PiNeq);
}

/* *************** Class RegularizedVelocityInnerCornerDynamics3D ****** */

template <typename T, template <typename U> class Descriptor, int normalX, int normalY, int normalZ>
RegularizedVelocityInnerCornerDynamics3D<T, Descriptor, normalX, normalY, normalZ>::
    RegularizedVelocityInnerCornerDynamics3D(
        Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_) :
    BoundaryCompositeDynamics<T, Descriptor>(baseDynamics_, automaticPrepareCollision_),
    xDynamics(baseDynamics_->clone()),
    yDynamics(baseDynamics_->clone()),
    zDynamics(baseDynamics_->clone())

{ }

template <typename T, template <typename U> class Descriptor, int normalX, int normalY, int normalZ>
RegularizedVelocityInnerCornerDynamics3D<T, Descriptor, normalX, normalY, normalZ>
    *RegularizedVelocityInnerCornerDynamics3D<T, Descriptor, normalX, normalY, normalZ>::clone()
        const
{
    return new RegularizedVelocityInnerCornerDynamics3D<T, Descriptor, normalX, normalY, normalZ>(
        *this);
}

template <typename T, template <typename U> class Descriptor, int normalX, int normalY, int normalZ>
void RegularizedVelocityInnerCornerDynamics3D<T, Descriptor, normalX, normalY, normalZ>::
    replaceBaseDynamics(Dynamics<T, Descriptor> *newBaseDynamics)
{
    BoundaryCompositeDynamics<T, Descriptor>::replaceBaseDynamics(newBaseDynamics);
    xDynamics.replaceBaseDynamics(newBaseDynamics->clone());
    yDynamics.replaceBaseDynamics(newBaseDynamics->clone());
    zDynamics.replaceBaseDynamics(newBaseDynamics->clone());
}

template <typename T, template <typename U> class Descriptor, int normalX, int normalY, int normalZ>
void RegularizedVelocityInnerCornerDynamics3D<T, Descriptor, normalX, normalY, normalZ>::
    computeVelocity(Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u_) const
{
    xDynamics.computeVelocity(cell, u_);
}

template <typename T, template <typename U> class Descriptor, int normalX, int normalY, int normalZ>
void RegularizedVelocityInnerCornerDynamics3D<T, Descriptor, normalX, normalY, normalZ>::
    defineVelocity(Cell<T, Descriptor> &cell, Array<T, Descriptor<T>::d> const &u_)
{
    xDynamics.defineVelocity(cell, u_);
    yDynamics.defineVelocity(cell, u_);
    zDynamics.defineVelocity(cell, u_);
}

template <typename T, template <typename U> class Descriptor, int normalX, int normalY, int normalZ>
T RegularizedVelocityInnerCornerDynamics3D<
    T, Descriptor, normalX, normalY, normalZ>::computeDensity(Cell<T, Descriptor> const &cell) const
{
    return (xDynamics.computeDensity(cell) + yDynamics.computeDensity(cell)
            + zDynamics.computeDensity(cell))
           / (T)3;
}

template <typename T, template <typename U> class Descriptor, int normalX, int normalY, int normalZ>
T RegularizedVelocityInnerCornerDynamics3D<T, Descriptor, normalX, normalY, normalZ>::computeRhoBar(
    Cell<T, Descriptor> const &cell) const
{
    return (xDynamics.computeRhoBar(cell) + yDynamics.computeRhoBar(cell)
            + zDynamics.computeRhoBar(cell))
           / (T)3;
}

template <typename T, template <typename U> class Descriptor, int normalX, int normalY, int normalZ>
void RegularizedVelocityInnerCornerDynamics3D<T, Descriptor, normalX, normalY, normalZ>::
    computeRhoBarJ(
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

template <typename T, template <typename U> class Descriptor, int normalX, int normalY, int normalZ>
void RegularizedVelocityInnerCornerDynamics3D<T, Descriptor, normalX, normalY, normalZ>::
    computeRhoBarJPiNeq(
        Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const
{
    T rhoBar2, rhoBar3;
    Array<T, Descriptor<T>::d> j2, j3;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq2, PiNeq3;

    xDynamics.computeRhoBarJPiNeq(cell, rhoBar, j, PiNeq);
    yDynamics.computeRhoBarJPiNeq(cell, rhoBar2, j2, PiNeq2);
    zDynamics.computeRhoBarJPiNeq(cell, rhoBar3, j3, PiNeq3);

    rhoBar = (rhoBar + rhoBar2 + rhoBar3) / (T)3;
    for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
        j[iD] = (j[iD] + j2[iD] + j3[iD]) / (T)3;
    }
    for (int iPi = 0; iPi < Descriptor<T>::d; ++iPi) {
        PiNeq[iPi] = (PiNeq[iPi] + PiNeq2[iPi] + PiNeq3[iPi]) / (T)3;
    }
}

template <typename T, template <typename U> class Descriptor, int normalX, int normalY, int normalZ>
void RegularizedVelocityInnerCornerDynamics3D<
    T, Descriptor, normalX, normalY, normalZ>::completePopulations(Cell<T, Descriptor> &cell) const
{
    // 1. Assign "bounce-back of off-equilibrium" to the unknown population
    Array<int, Descriptor<T>::d> v(-normalX, -normalY, -normalZ);
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

// namespace plb

#endif  // REGULARIZED_BOUNDARY_DYNAMICS_3D_HH
