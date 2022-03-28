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

#ifndef FINITE_DIFFERENCE_BOUNDARY_PROCESSOR_3D_HH
#define FINITE_DIFFERENCE_BOUNDARY_PROCESSOR_3D_HH

#include <typeinfo>

#include "atomicBlock/blockLattice3D.h"
#include "boundaryCondition/finiteDifferenceBoundaryProcessor3D.h"
#include "core/processorIdentifiers3D.h"
#include "finiteDifference/finiteDifference3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/offEquilibriumTemplates.h"

namespace plb {

namespace fdBoundaryAlgorithm {

}  // namespace fdBoundaryAlgorithm

////////  PlaneFdBoundaryFunctional3D ///////////////////////////////////

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
const int PlaneFdBoundaryFunctional3D<T, Descriptor, direction, orientation>::staticId =
    meta::registerProcessor3D<
        PlaneFdBoundaryFunctional3D<T, Descriptor, direction, orientation>, T, Descriptor,
        direction, orientation>(std::string("PlaneFdBoundary3D"));

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void PlaneFdBoundaryFunctional3D<T, Descriptor, direction, orientation>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    PLB_ASSERT(domain.x0 == domain.x1 || domain.y0 == domain.y1 || domain.z0 == domain.z1);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                processCell(iX, iY, iZ, lattice);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void PlaneFdBoundaryFunctional3D<T, Descriptor, direction, orientation>::processCell(
    plint iX, plint iY, plint iZ, BlockLattice3D<T, Descriptor> &lattice)
{
    typedef SymmetricTensorImpl<T, Descriptor<T>::d> S;
    Array<T, Descriptor<T>::d> dx_u, dy_u, dz_u;
    Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
    Dynamics<T, Descriptor> &dynamics = cell.getDynamics();
    T rho = cell.computeDensity();
    Array<T, Descriptor<T>::d> vel;
    cell.computeVelocity(vel);

    interpolateGradients<0>(lattice, dx_u, iX, iY, iZ);
    interpolateGradients<1>(lattice, dy_u, iX, iY, iZ);
    interpolateGradients<2>(lattice, dz_u, iX, iY, iZ);
    T dx_ux = dx_u[0];
    T dy_ux = dy_u[0];
    T dz_ux = dz_u[0];
    T dx_uy = dx_u[1];
    T dy_uy = dy_u[1];
    T dz_uy = dz_u[1];
    T dx_uz = dx_u[2];
    T dy_uz = dy_u[2];
    T dz_uz = dz_u[2];
    T omega = cell.getDynamics().getOmega();
    T sToPi = -rho / Descriptor<T>::invCs2 / omega;
    Array<T, SymmetricTensor<T, Descriptor>::n> pi;
    pi[S::xx] = (T)2 * dx_ux * sToPi;
    pi[S::yy] = (T)2 * dy_uy * sToPi;
    pi[S::zz] = (T)2 * dz_uz * sToPi;
    pi[S::xy] = (dx_uy + dy_ux) * sToPi;
    pi[S::xz] = (dx_uz + dz_ux) * sToPi;
    pi[S::yz] = (dy_uz + dz_uy) * sToPi;

    Array<T, Descriptor<T>::d> j;
    if (cell.getDynamics().velIsJ()) {
        j = vel;
    } else {
        j = rho * vel;
    }
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    dynamics.regularize(cell, Descriptor<T>::rhoBar(rho), j, jSqr, pi);

    //     // Computation of the particle distribution functions
    //     // according to the regularized formula
    //     for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop)
    //         cell[iPop] = dynamics.computeEquilibrium(iPop,Descriptor<T>::rhoBar(rho),j,jSqr) +
    //                          offEquilibriumTemplates<T,Descriptor>::fromPiToFneq(iPop, pi);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
PlaneFdBoundaryFunctional3D<T, Descriptor, direction, orientation>
    *PlaneFdBoundaryFunctional3D<T, Descriptor, direction, orientation>::clone() const
{
    return new PlaneFdBoundaryFunctional3D<T, Descriptor, direction, orientation>(*this);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
template <int deriveDirection>
void PlaneFdBoundaryFunctional3D<T, Descriptor, direction, orientation>::interpolateGradients(
    BlockLattice3D<T, Descriptor> const &lattice, Array<T, Descriptor<T>::d> &velDeriv, plint iX,
    plint iY, plint iZ)
{
    fd::DirectedGradients3D<
        T, Descriptor, direction, orientation, deriveDirection,
        direction == deriveDirection>::o1_velocityDerivative(velDeriv, lattice, iX, iY, iZ);
}

////////  OuterVelocityEdgeFunctional3D ///////////////////////////////////

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
const int OuterVelocityEdgeFunctional3D<T, Descriptor, plane, normal1, normal2>::staticId =
    meta::registerProcessor3D<
        OuterVelocityEdgeFunctional3D<T, Descriptor, plane, normal1, normal2>, T, Descriptor, plane,
        normal1, normal2>(std::string("OuterVelocityEdge"));

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
void OuterVelocityEdgeFunctional3D<T, Descriptor, plane, normal1, normal2>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    PLB_ASSERT(
        (plane == 2 && domain.x0 == domain.x1 && domain.y0 == domain.y1)
        || (plane == 1 && domain.x0 == domain.x1 && domain.z0 == domain.z1)
        || (plane == 0 && domain.y0 == domain.y1 && domain.z0 == domain.z1));

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                processCell(iX, iY, iZ, lattice);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
void OuterVelocityEdgeFunctional3D<T, Descriptor, plane, normal1, normal2>::processCell(
    plint iX, plint iY, plint iZ, BlockLattice3D<T, Descriptor> &lattice)
{
    typedef SymmetricTensorImpl<T, Descriptor<T>::d> S;
    Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
    Dynamics<T, Descriptor> &dynamics = cell.getDynamics();

    T rho10 = getNeighborRho(iX, iY, iZ, 1, 0, lattice);
    T rho01 = getNeighborRho(iX, iY, iZ, 0, 1, lattice);
    T rho = (T)0.5 * (rho01 + rho10);

    std::vector<Array<T, 3> > dA_uB_(3);
    interpolateGradients<plane, 0>(lattice, dA_uB_[0], iX, iY, iZ);
    interpolateGradients<direction1, normal1>(lattice, dA_uB_[1], iX, iY, iZ);
    interpolateGradients<direction2, normal2>(lattice, dA_uB_[2], iX, iY, iZ);
    std::vector<Array<T, 3> > dA_uB(3);
    for (int iBeta = 0; iBeta < 3; ++iBeta) {
        dA_uB[plane][iBeta] = dA_uB_[0][iBeta];
        dA_uB[direction1][iBeta] = dA_uB_[1][iBeta];
        dA_uB[direction2][iBeta] = dA_uB_[2][iBeta];
    }
    T omega = dynamics.getOmega();
    T sToPi = -rho / Descriptor<T>::invCs2 / omega;
    Array<T, SymmetricTensor<T, Descriptor>::n> pi;
    pi[S::xx] = (T)2 * dA_uB[0][0] * sToPi;
    pi[S::yy] = (T)2 * dA_uB[1][1] * sToPi;
    pi[S::zz] = (T)2 * dA_uB[2][2] * sToPi;
    pi[S::xy] = (dA_uB[0][1] + dA_uB[1][0]) * sToPi;
    pi[S::xz] = (dA_uB[0][2] + dA_uB[2][0]) * sToPi;
    pi[S::yz] = (dA_uB[1][2] + dA_uB[2][1]) * sToPi;

    // Computation of the particle distribution functions
    // according to the regularized formula
    Array<T, Descriptor<T>::d> vel, j;
    cell.computeVelocity(vel);

    if (cell.getDynamics().velIsJ()) {
        j = vel;
    } else {
        j = rho * vel;
    }
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    dynamics.regularize(cell, Descriptor<T>::rhoBar(rho), j, jSqr, pi);

    //     for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
    //         cell[iPop] = dynamics.computeEquilibrium(iPop,Descriptor<T>::rhoBar(rho),j,jSqr) +
    //                          offEquilibriumTemplates<T,Descriptor>::fromPiToFneq(iPop, pi);
    //     }
}

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
OuterVelocityEdgeFunctional3D<T, Descriptor, plane, normal1, normal2>
    *OuterVelocityEdgeFunctional3D<T, Descriptor, plane, normal1, normal2>::clone() const
{
    return new OuterVelocityEdgeFunctional3D<T, Descriptor, plane, normal1, normal2>(*this);
}

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
T OuterVelocityEdgeFunctional3D<T, Descriptor, plane, normal1, normal2>::getNeighborRho(
    plint x, plint y, plint z, plint step1, plint step2,
    BlockLattice3D<T, Descriptor> const &lattice)
{
    Array<int, 3> coords(x, y, z);
    coords[direction1] += -normal1 * step1;
    coords[direction2] += -normal2 * step2;
    return lattice.get(coords[0], coords[1], coords[2]).computeDensity();
}

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
template <int deriveDirection, int orientation>
void OuterVelocityEdgeFunctional3D<T, Descriptor, plane, normal1, normal2>::interpolateGradients(
    BlockLattice3D<T, Descriptor> const &lattice, Array<T, Descriptor<T>::d> &velDeriv, plint iX,
    plint iY, plint iZ)
{
    fd::DirectedGradients3D<
        T, Descriptor, deriveDirection, orientation, deriveDirection,
        deriveDirection != plane>::o1_velocityDerivative(velDeriv, lattice, iX, iY, iZ);
}

/////////// OuterVelocityCornerFunctional3D /////////////////////////////////////

template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
const int OuterVelocityCornerFunctional3D<T, Descriptor, xNormal, yNormal, zNormal>::staticId =
    meta::registerProcessor3D<
        OuterVelocityCornerFunctional3D<T, Descriptor, xNormal, yNormal, zNormal>, T, Descriptor,
        xNormal, yNormal, zNormal>(std::string("OuterVelocityCorner3D"));

template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
void OuterVelocityCornerFunctional3D<T, Descriptor, xNormal, yNormal, zNormal>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    plint x = domain.x0;
    plint y = domain.y0;
    plint z = domain.z0;
    PLB_ASSERT(x == domain.x1 && y == domain.y1 && z == domain.z1);

    processCell(x, y, z, lattice);
}

template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
void OuterVelocityCornerFunctional3D<T, Descriptor, xNormal, yNormal, zNormal>::processCell(
    plint iX, plint iY, plint iZ, BlockLattice3D<T, Descriptor> &lattice)
{
    typedef SymmetricTensorImpl<T, Descriptor<T>::d> S;
    Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
    Dynamics<T, Descriptor> &dynamics = cell.getDynamics();

    T rho100 = lattice.get(iX - 1 * xNormal, iY - 0 * yNormal, iZ - 0 * zNormal).computeDensity();
    T rho010 = lattice.get(iX - 0 * xNormal, iY - 1 * yNormal, iZ - 0 * zNormal).computeDensity();
    T rho001 = lattice.get(iX - 0 * xNormal, iY - 0 * yNormal, iZ - 1 * zNormal).computeDensity();
    T rho = (T)1 / (T)3 * (rho001 + rho010 + rho100);

    Array<T, Descriptor<T>::d> dx_u, dy_u, dz_u;
    fd::DirectedGradients3D<T, Descriptor, 0, xNormal, 0, true>::o1_velocityDerivative(
        dx_u, lattice, iX, iY, iZ);
    fd::DirectedGradients3D<T, Descriptor, 1, yNormal, 0, true>::o1_velocityDerivative(
        dy_u, lattice, iX, iY, iZ);
    fd::DirectedGradients3D<T, Descriptor, 2, zNormal, 0, true>::o1_velocityDerivative(
        dz_u, lattice, iX, iY, iZ);

    T dx_ux = dx_u[0];
    T dy_ux = dy_u[0];
    T dz_ux = dz_u[0];
    T dx_uy = dx_u[1];
    T dy_uy = dy_u[1];
    T dz_uy = dz_u[1];
    T dx_uz = dx_u[2];
    T dy_uz = dy_u[2];
    T dz_uz = dz_u[2];
    T omega = dynamics.getOmega();
    T sToPi = -rho / Descriptor<T>::invCs2 / omega;
    Array<T, SymmetricTensor<T, Descriptor>::n> pi;
    pi[S::xx] = (T)2 * dx_ux * sToPi;
    pi[S::yy] = (T)2 * dy_uy * sToPi;
    pi[S::zz] = (T)2 * dz_uz * sToPi;
    pi[S::xy] = (dx_uy + dy_ux) * sToPi;
    pi[S::xz] = (dx_uz + dz_ux) * sToPi;
    pi[S::yz] = (dy_uz + dz_uy) * sToPi;

    // Computation of the particle distribution functions
    // according to the regularized formula
    Array<T, Descriptor<T>::d> vel, j;
    cell.computeVelocity(vel);
    if (cell.getDynamics().velIsJ()) {
        j = vel;
    } else {
        j = rho * vel;
    }
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    dynamics.regularize(cell, Descriptor<T>::rhoBar(rho), j, jSqr, pi);

    //     for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
    //         cell[iPop] = dynamics.computeEquilibrium(iPop,Descriptor<T>::rhoBar(rho),j,jSqr) +
    //                          offEquilibriumTemplates<T,Descriptor>::fromPiToFneq(iPop, pi);
    //     }
}

template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
OuterVelocityCornerFunctional3D<T, Descriptor, xNormal, yNormal, zNormal>
    *OuterVelocityCornerFunctional3D<T, Descriptor, xNormal, yNormal, zNormal>::clone() const
{
    return new OuterVelocityCornerFunctional3D<T, Descriptor, xNormal, yNormal, zNormal>(*this);
}

}  // namespace plb

#endif  // FINITE_DIFFERENCE_BOUNDARY_PROCESSOR_3D_HH
