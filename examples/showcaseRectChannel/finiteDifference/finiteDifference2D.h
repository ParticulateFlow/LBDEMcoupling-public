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

#ifndef FINITE_DIFFERENCE_2D_H
#define FINITE_DIFFERENCE_2D_H

#include "core/globalDefs.h"
#include "finiteDifference/fdStencils1D.h"

namespace plb {

namespace fd {

template <
    typename T, template <typename U> class Descriptor, int direction, int orientation,
    bool orthogonal>
struct DirectedGradients2D {
    /// Nearest-neighbor evaluation of velocity derivative, first-order accurate
    ///   only along the boundary normal.
    static void o1_velocityDerivative(
        Array<T, Descriptor<T>::d> &velDeriv, BlockLattice2D<T, Descriptor> const &blockLattice,
        plint iX, plint iY);
    /// Nearest-neighbor evaluation of density derivative, first-order accurate
    ///   only along the boundary normal.
    static void o1_densityDerivative(
        T &rhoDeriv, BlockLattice2D<T, Descriptor> const &blockLattice, plint iX, plint iY);
    /// Next-to-nearest-neibhbor, second-order accurate evaluation of velocity
    ///   derivative.
    static void o2_velocityDerivative(
        Array<T, Descriptor<T>::d> &velDeriv, BlockLattice2D<T, Descriptor> const &blockLattice,
        plint iX, plint iY);
    /// Next-to-nearest-neibhbor, second-order accurate evaluation of density
    ///   derivative.
    static void o2_densityDerivative(
        T &rhoDeriv, BlockLattice2D<T, Descriptor> const &blockLattice, plint iX, plint iY);
};

// Implementation for orthogonal==true; i.e. the derivative is along
//   the boundary normal.
template <typename T, template <typename U> class Descriptor, int direction, int orientation>
struct DirectedGradients2D<T, Descriptor, direction, orientation, true> {
    static void o1_velocityDerivative(
        Array<T, Descriptor<T>::d> &velDeriv, BlockLattice2D<T, Descriptor> const &blockLattice,
        plint iX, plint iY)
    {
        Array<T, Descriptor<T>::d> u0, u1;

        blockLattice.get(iX, iY).computeVelocity(u0);
        blockLattice
            .get(
                iX + (direction == 0 ? (-orientation) : 0),
                iY + (direction == 1 ? (-orientation) : 0))
            .computeVelocity(u1);

        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            velDeriv[iD] = -orientation * fd::o1_fwd_diff(u0[iD], u1[iD]);
        }
    }

    static void o1_densityDerivative(
        T &rhoDeriv, BlockLattice2D<T, Descriptor> const &blockLattice, plint iX, plint iY)
    {
        T rho0 = blockLattice.get(iX, iY).computeDensity();
        T rho1 = blockLattice
                     .get(
                         iX + (direction == 0 ? (-orientation) : 0),
                         iY + (direction == 1 ? (-orientation) : 0))
                     .computeDensity();

        rhoDeriv = -orientation * fd::o1_fwd_diff(rho0, rho1);
    }
    static void o2_velocityDerivative(
        Array<T, Descriptor<T>::d> &velDeriv, BlockLattice2D<T, Descriptor> const &blockLattice,
        plint iX, plint iY)
    {
        Array<T, Descriptor<T>::d> u0, u1, u2;

        blockLattice.get(iX, iY).computeVelocity(u0);
        blockLattice
            .get(
                iX + (direction == 0 ? (-orientation) : 0),
                iY + (direction == 1 ? (-orientation) : 0))
            .computeVelocity(u1);
        blockLattice
            .get(
                iX + (direction == 0 ? (-2 * orientation) : 0),
                iY + (direction == 1 ? (-2 * orientation) : 0))
            .computeVelocity(u2);

        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            velDeriv[iD] = -orientation * fd::fwd_diff(u0[iD], u1[iD], u2[iD]);
        }
    }

    static void o2_densityDerivative(
        T &rhoDeriv, BlockLattice2D<T, Descriptor> const &blockLattice, plint iX, plint iY)
    {
        T rho0 = blockLattice.get(iX, iY).computeDensity();
        T rho1 = blockLattice
                     .get(
                         iX + (direction == 0 ? (-orientation) : 0),
                         iY + (direction == 1 ? (-orientation) : 0))
                     .computeDensity();
        T rho2 = blockLattice
                     .get(
                         iX + (direction == 0 ? (-2 * orientation) : 0),
                         iY + (direction == 1 ? (-2 * orientation) : 0))
                     .computeDensity();

        rhoDeriv = -orientation * fd::fwd_diff(rho0, rho1, rho2);
    }
};

// Implementation for orthogonal==false; i.e. the derivative is aligned
// with the boundary.
template <typename T, template <typename U> class Descriptor, int direction, int orientation>
struct DirectedGradients2D<T, Descriptor, direction, orientation, false> {
    static void o1_velocityDerivative(
        Array<T, Descriptor<T>::d> &velDeriv, BlockLattice2D<T, Descriptor> const &blockLattice,
        plint iX, plint iY)
    {
        // Along the boundary, second-order accuracy is achieved with a nearest-
        //   neighbor scheme.
        o2_velocityDerivative(velDeriv, blockLattice, iX, iY);
    }
    static void o1_densityDerivative(
        T &rhoDeriv, BlockLattice2D<T, Descriptor> const &blockLattice, plint iX, plint iY)
    {
        // Along the boundary, second-order accuracy is achieved with a nearest-
        //   neighbor scheme.
        o2_densityDerivative(rhoDeriv, blockLattice, iX, iY);
    }
    static void o2_velocityDerivative(
        Array<T, Descriptor<T>::d> &velDeriv, BlockLattice2D<T, Descriptor> const &blockLattice,
        plint iX, plint iY)
    {
        Array<T, Descriptor<T>::d> u_p1, u_m1;

        int deriveDirection = 1 - direction;
        blockLattice.get(iX + (deriveDirection == 0 ? 1 : 0), iY + (deriveDirection == 1 ? 1 : 0))
            .computeVelocity(u_p1);
        blockLattice
            .get(iX + (deriveDirection == 0 ? (-1) : 0), iY + (deriveDirection == 1 ? (-1) : 0))
            .computeVelocity(u_m1);

        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            velDeriv[iD] = fd::ctl_diff(u_p1[iD], u_m1[iD]);
        }
    }

    static void o2_densityDerivative(
        T &rhoDeriv, BlockLattice2D<T, Descriptor> const &blockLattice, plint iX, plint iY)
    {
        int deriveDirection = 1 - direction;
        T rho_p1 =
            blockLattice
                .get(iX + (deriveDirection == 0 ? 1 : 0), iY + (deriveDirection == 1 ? 1 : 0))
                .computeDensity();
        T rho_m1 =
            blockLattice
                .get(iX + (deriveDirection == 0 ? (-1) : 0), iY + (deriveDirection == 1 ? (-1) : 0))
                .computeDensity();

        rhoDeriv = fd::ctl_diff(rho_p1, rho_m1);
    }
};

}  // namespace fd

}  // namespace plb

#endif
