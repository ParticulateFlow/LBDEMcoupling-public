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

#ifndef KBCTEMPLATES2D_H
#define KBCTEMPLATES2D_H

#include "core/globalDefs.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/momentTemplates2D.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"

namespace plb {

template <typename T>
struct kbcTemplates<T, descriptors::D2Q9Descriptor> {
    typedef descriptors::D2Q9DescriptorBase<T> Descriptor;

    static void computeMoments(Array<T, Descriptor::q> &moments, const Array<T, Descriptor::q> &f)
    {
        typedef momentTemplatesImpl<T, Descriptor> momentTemp;

        T rho;
        Array<T, 2> u;
        momentTemp::compute_rho_uLb(f, rho, u);
        T rhoBar = Descriptor::rhoBar(rho);
        T invRho = Descriptor::invRho(rhoBar);

        T Pi_xy = invRho * (-f[1] + f[3] - f[5] + f[7]) - u[0] * u[1];
        T T_ = invRho * (momentTemp::compute_e(f)) - normSqr(u);
        T N_ = invRho * (f[2] + f[6] - f[4] - f[8]) - u[0] * u[0] + u[1] * u[1];
        T Q_xyy = invRho * (-f[1] - f[3] + f[5] + f[7]) - 2 * u[1] * Pi_xy + 0.5 * u[0] * (N_ - T_)
                  - u[0] * u[1] * u[1];
        T Q_yxx = invRho * (f[1] - f[3] - f[5] + f[7]) - 2 * u[0] * Pi_xy - 0.5 * u[1] * (N_ + T_)
                  - u[1] * u[0] * u[0];
        T A_ = invRho
                   * (f[1] + f[3] + f[5] + f[7]
                      + Descriptor::SkordosFactor() * Descriptor::cs2 * Descriptor::cs2)
               - 2 * (u[0] * Q_xyy + u[1] * Q_yxx) - 4 * u[0] * u[1] * Pi_xy - 0.5 * normSqr(u) * T_
               + 0.5 * (u[0] * u[0] - u[1] * u[1]) * N_ - u[0] * u[0] * u[1] * u[1];

        moments[0] = rho;
        moments[1] = u[0];
        moments[2] = u[1];
        moments[3] = T_;
        moments[4] = N_;
        moments[5] = Pi_xy;
        moments[6] = Q_xyy;
        moments[7] = Q_yxx;
        moments[8] = A_;
    }

    static void computeMomentsForS(
        Array<T, Descriptor::q> &moments, const Array<T, Descriptor::q> &f)
    {
        typedef momentTemplatesImpl<T, Descriptor> momentTemp;

        T rho;
        Array<T, 2> u;
        momentTemp::compute_rho_uLb(f, rho, u);
        T rhoBar = Descriptor::rhoBar(rho);
        T invRho = Descriptor::invRho(rhoBar);

        T Pi_xy = invRho * (-f[1] + f[3] - f[5] + f[7]) - u[0] * u[1];
        T T_ = invRho * (momentTemp::compute_e(f)) - normSqr(u);
        T N_ = invRho * (f[2] + f[6] - f[4] - f[8]) - u[0] * u[0] + u[1] * u[1];

        moments[0] = rho;
        moments[1] = u[0];
        moments[2] = u[1];
        moments[3] = T_;
        moments[4] = N_;
        moments[5] = Pi_xy;
    }

    static void computeEqMoments(
        Array<T, Descriptor::q> &momentsEq, Array<T, Descriptor::q> moments)
    {
        momentsEq = Array<T, 9>(
            moments[0], moments[1], moments[2], 2 * Descriptor::cs2, (T)0., (T)0., (T)0., (T)0.,
            Descriptor::cs2 * Descriptor::cs2);
    }

    static void computeK(Array<T, Descriptor::q> &k, Array<T, Descriptor::q> moments)
    {
        T rho = moments[0];
        T ux = moments[1];
        T uy = moments[2];
        T uSqr = ux * ux + uy * uy;

        k[0] = rho * (1 - uSqr);
        k[1] = -rho * 0.25 * ux * uy;
        k[2] = 0.5 * rho * (ux * ux - ux);
        k[3] = rho * 0.25 * ux * uy;
        k[4] = 0.5 * rho * (uy * uy - uy);
        k[5] = k[1];
        k[6] = 0.5 * rho * (ux * ux + ux);
        k[7] = k[3];
        k[8] = 0.5 * rho * (uy * uy + uy);
    }

    static void computeS(Array<T, Descriptor::q> &s, Array<T, Descriptor::q> moments)
    {
        T rho = moments[0];
        T ux = moments[1];
        T uy = moments[2];
        T uSqr = ux * ux + uy * uy;
        T T_ = moments[3];
        T N_ = moments[4];
        T Pi_xy = moments[5];

        s[0] = rho * (4 * ux * uy * Pi_xy - 0.5 * (ux * ux - uy * uy) * N_ + 0.5 * (uSqr - 2) * T_);

        for (plint iPop = 1; iPop < Descriptor::q; ++iPop) {
            plint sigma = Descriptor::c[iPop][0];
            plint lambda = Descriptor::c[iPop][1];
            if (lambda == 0) {
                s[iPop] = 0.5 * rho
                          * (0.5 * (1 + sigma * ux + ux * ux - uy * uy) * N_
                             - (2 * sigma * uy + 4 * ux * uy) * Pi_xy
                             + 0.5 * (1 - sigma * ux - uSqr) * T_);
            } else if (sigma == 0) {
                s[iPop] = 0.5 * rho
                          * (0.5 * (-1 - lambda * uy + ux * ux - uy * uy) * N_
                             - (2 * lambda * ux + 4 * ux * uy) * Pi_xy
                             + 0.5 * (1 - lambda * uy - uSqr) * T_);
            } else {
                s[iPop] =
                    0.25 * rho
                    * ((4 * ux * uy + sigma * lambda + 2 * sigma * uy + 2 * lambda * ux) * Pi_xy
                       + 0.5 * (-ux * ux + uy * uy - sigma * ux + lambda * uy) * N_
                       + 0.5 * (uSqr + sigma * ux + lambda * uy) * T_);
            }
        }
    }

    static void computeH(Array<T, Descriptor::q> &h, Array<T, Descriptor::q> moments)
    {
        T rho = moments[0];
        T ux = moments[1];
        T uy = moments[2];
        T Q_xyy = moments[6];
        T Q_yxx = moments[7];
        T A_ = moments[8];

        h[0] = rho * (2 * ux * Q_xyy + 2 * uy * Q_yxx + A_);

        for (plint iPop = 1; iPop < Descriptor::q; ++iPop) {
            plint sigma = Descriptor::c[iPop][0];
            plint lambda = Descriptor::c[iPop][1];
            if (lambda == 0) {
                h[iPop] = 0.5 * rho * (-(sigma + 2 * ux) * Q_xyy - 2 * uy * Q_yxx - A_);
            } else if (sigma == 0) {
                h[iPop] = 0.5 * rho * (-(lambda + 2 * uy) * Q_yxx - 2 * ux * Q_xyy - A_);
            } else {
                h[iPop] = 0.25 * rho * ((sigma + 2 * ux) * Q_xyy + (lambda + 2 * uy) * Q_yxx + A_);
            }
        }
    }

    static void computeKEq(Array<T, Descriptor::q> &kEq, Array<T, Descriptor::q> moments)
    {
        computeK(kEq, moments);
    }

    static void computeSEq(Array<T, Descriptor::q> &sEq, Array<T, Descriptor::q> moments)
    {
        T rho = moments[0];
        T ux = moments[1];
        T uy = moments[2];
        T uSqr = ux * ux + uy * uy;

        sEq[0] = rho * (uSqr - 2) * Descriptor::cs2;

        for (plint iPop = 1; iPop < Descriptor::q; ++iPop) {
            plint sigma = Descriptor::c[iPop][0];
            plint lambda = Descriptor::c[iPop][1];
            if (lambda == 0 || sigma == 0) {
                sEq[iPop] = 0.5 * rho * (1 - sigma * ux - lambda * uy - uSqr) * Descriptor::cs2;
            } else {
                sEq[iPop] = 0.25 * rho * (uSqr + sigma * ux + lambda * uy) * Descriptor::cs2;
            }
        }
    }

    static void computeHEq(Array<T, Descriptor::q> &hEq, Array<T, Descriptor::q> moments)
    {
        T rho = moments[0];
        T cs4 = Descriptor::cs2 * Descriptor::cs2;

        hEq[0] = rho * cs4;

        for (plint iPop = 1; iPop < Descriptor::q; ++iPop) {
            plint sigma = Descriptor::c[iPop][0];
            plint lambda = Descriptor::c[iPop][1];
            if (lambda == 0 || sigma == 0) {
                hEq[iPop] = -0.5 * rho * cs4;
            } else {
                hEq[iPop] = 0.25 * rho * cs4;
            }
        }
    }
};

}  // namespace plb

#endif  // KBCTEMPLATES2D_H
