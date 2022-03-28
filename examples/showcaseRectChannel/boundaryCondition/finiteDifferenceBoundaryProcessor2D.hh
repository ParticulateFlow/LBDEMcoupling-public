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

#ifndef FINITE_DIFFERENCE_BOUNDARY_PROCESSOR_2D_HH
#define FINITE_DIFFERENCE_BOUNDARY_PROCESSOR_2D_HH

#include <typeinfo>

#include "atomicBlock/blockLattice2D.h"
#include "boundaryCondition/finiteDifferenceBoundaryProcessor2D.h"
#include "core/processorIdentifiers2D.h"
#include "finiteDifference/finiteDifference2D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/offEquilibriumTemplates.h"

namespace plb {

///////////  StraightFdBoundaryFunctional2D ///////////////////////////////////

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
const int StraightFdBoundaryFunctional2D<T, Descriptor, direction, orientation>::staticId =
    meta::registerProcessor2D<
        StraightFdBoundaryFunctional2D<T, Descriptor, direction, orientation>, T, Descriptor,
        direction, orientation>(std::string("StraightFdBoundary2D"));

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void StraightFdBoundaryFunctional2D<T, Descriptor, direction, orientation>::processCell(
    plint iX, plint iY, BlockLattice2D<T, Descriptor> &lattice)
{
    typedef SymmetricTensorImpl<T, Descriptor<T>::d> S;
    Array<T, Descriptor<T>::d> dx_u, dy_u;
    Cell<T, Descriptor> &cell = lattice.get(iX, iY);

    T rho = cell.computeDensity();
    Array<T, Descriptor<T>::d> vel;
    cell.computeVelocity(vel);
    interpolateGradients<0>(lattice, dx_u, iX, iY);
    interpolateGradients<1>(lattice, dy_u, iX, iY);
    T dx_ux = dx_u[0];
    T dy_ux = dy_u[0];
    T dx_uy = dx_u[1];
    T dy_uy = dy_u[1];
    T omega = cell.getDynamics().getOmega();
    T sToPi = -rho / Descriptor<T>::invCs2 / omega;
    Array<T, SymmetricTensor<T, Descriptor>::n> pi;
    pi[S::xx] = (T)2 * dx_ux * sToPi;
    pi[S::yy] = (T)2 * dy_uy * sToPi;
    pi[S::xy] = (dx_uy + dy_ux) * sToPi;

    Array<T, Descriptor<T>::d> j;
    if (cell.getDynamics().velIsJ()) {
        j = vel;
    } else {
        j = rho * vel;
    }
    // Computation of the particle distribution functions
    // according to the regularized formula
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    cell.getDynamics().regularize(cell, Descriptor<T>::rhoBar(rho), j, jSqr, pi);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void StraightFdBoundaryFunctional2D<T, Descriptor, direction, orientation>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    PLB_PRECONDITION(domain.x0 == domain.x1 || domain.y0 == domain.y1);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            processCell(iX, iY, lattice);
        }
    }
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
StraightFdBoundaryFunctional2D<T, Descriptor, direction, orientation>
    *StraightFdBoundaryFunctional2D<T, Descriptor, direction, orientation>::clone() const
{
    return new StraightFdBoundaryFunctional2D<T, Descriptor, direction, orientation>(*this);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
template <int deriveDirection>
void StraightFdBoundaryFunctional2D<T, Descriptor, direction, orientation>::interpolateGradients(
    BlockLattice2D<T, Descriptor> const &lattice, Array<T, Descriptor<T>::d> &velDeriv, plint iX,
    plint iY)
{
    fd::DirectedGradients2D<T, Descriptor, direction, orientation, direction == deriveDirection>::
        o1_velocityDerivative(velDeriv, lattice, iX, iY);
}

/////////// OuterVelocityCornerFunctional2D /////////////////////////////////////

template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal>
const int OuterVelocityCornerFunctional2D<T, Descriptor, xNormal, yNormal>::staticId =
    meta::registerProcessor2D<
        OuterVelocityCornerFunctional2D<T, Descriptor, xNormal, yNormal>, T, Descriptor, xNormal,
        yNormal>(std::string("OuterVelocityCorner2D"));

template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal>
void OuterVelocityCornerFunctional2D<T, Descriptor, xNormal, yNormal>::processCell(
    plint iX, plint iY, BlockLattice2D<T, Descriptor> &lattice)
{
    typedef SymmetricTensorImpl<T, Descriptor<T>::d> S;

    T rho10 = lattice.get(iX - 1 * xNormal, iY - 0 * yNormal).computeDensity();
    T rho01 = lattice.get(iX - 0 * xNormal, iY - 1 * yNormal).computeDensity();

    T rho = (T)1 / (T)2 * (rho01 + rho10);

    Array<T, Descriptor<T>::d> dx_u, dy_u;
    fd::DirectedGradients2D<T, Descriptor, 0, xNormal, true>::o1_velocityDerivative(
        dx_u, lattice, iX, iY);
    fd::DirectedGradients2D<T, Descriptor, 1, yNormal, true>::o1_velocityDerivative(
        dy_u, lattice, iX, iY);
    T dx_ux = dx_u[0];
    T dy_ux = dy_u[0];
    T dx_uy = dx_u[1];
    T dy_uy = dy_u[1];

    Cell<T, Descriptor> &cell = lattice.get(iX, iY);
    T omega = cell.getDynamics().getOmega();

    T sToPi = -rho / Descriptor<T>::invCs2 / omega;
    Array<T, SymmetricTensor<T, Descriptor>::n> pi;
    pi[S::xx] = (T)2 * dx_ux * sToPi;
    pi[S::yy] = (T)2 * dy_uy * sToPi;
    pi[S::xy] = (dx_uy + dy_ux) * sToPi;

    // Computation of the particle distribution functions
    // according to the regularized formula
    Array<T, Descriptor<T>::d> vel, j;
    lattice.get(iX, iY).computeVelocity(vel);
    if (cell.getDynamics().velIsJ()) {
        j = vel;
    } else {
        for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
            j[iD] = rho * vel[iD];
        }
    }
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    cell.getDynamics().regularize(cell, Descriptor<T>::rhoBar(rho), j, jSqr, pi);
}

template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal>
void OuterVelocityCornerFunctional2D<T, Descriptor, xNormal, yNormal>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    plint x = domain.x0;
    plint y = domain.y0;
    PLB_ASSERT(x == domain.x1 && y == domain.y1);

    processCell(x, y, lattice);
}

template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal>
OuterVelocityCornerFunctional2D<T, Descriptor, xNormal, yNormal>
    *OuterVelocityCornerFunctional2D<T, Descriptor, xNormal, yNormal>::clone() const
{
    return new OuterVelocityCornerFunctional2D<T, Descriptor, xNormal, yNormal>(*this);
}

}  // namespace plb

#endif  // FINITE_DIFFERENCE_BOUNDARY_PROCESSOR_2D_HH
