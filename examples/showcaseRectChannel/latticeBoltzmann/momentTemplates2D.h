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
 * 2D specialization of momentTemplates functions.
 */

#ifndef MOMENT_TEMPLATES_2D_H
#define MOMENT_TEMPLATES_2D_H

#include "core/globalDefs.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"

namespace plb {

// Efficient specialization for D2Q9 base lattice
template <typename T>
struct momentTemplatesImpl<T, descriptors::D2Q9DescriptorBase<T> > {
    typedef descriptors::D2Q9DescriptorBase<T> Descriptor;

    static void partial_rho(
        Array<T, Descriptor::q> const &f, T &lineX_P1, T &lineX_0, T &lineX_M1, T &lineY_P1,
        T &lineY_M1)
    {
        lineX_P1 = f[5] + f[6] + f[7];
        lineX_0 = f[0] + f[4] + f[8];
        lineX_M1 = f[1] + f[2] + f[3];

        lineY_P1 = f[7] + f[8] + f[1];
        lineY_M1 = f[3] + f[4] + f[5];
    }

    static T get_rhoBar(Array<T, Descriptor::q> const &f)
    {
        T rhoBar = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
        return rhoBar;
    }

    static void get_j(Array<T, Descriptor::q> const &f, Array<T, 2> &j)
    {
        T lineX_P1, lineX_M1, lineY_P1, lineY_M1;

        lineX_P1 = f[5] + f[6] + f[7];
        lineX_M1 = f[1] + f[2] + f[3];
        lineY_P1 = f[7] + f[8] + f[1];
        lineY_M1 = f[3] + f[4] + f[5];

        j[0] = (lineX_P1 - lineX_M1);
        j[1] = (lineY_P1 - lineY_M1);
    }

    static T get_eBar(Array<T, Descriptor::q> const &f)
    {
        T eBar = (T)2 * (f[1] + f[3] + f[5] + f[7]) + f[2] + f[4] + f[6] + f[8];
        return eBar;
    }

    static void get_rhoBar_j(Array<T, Descriptor::q> const &f, T &rhoBar, Array<T, 2> &j)
    {
        T lineX_P1, lineX_0, lineX_M1, lineY_P1, lineY_M1;
        partial_rho(f, lineX_P1, lineX_0, lineX_M1, lineY_P1, lineY_M1);

        rhoBar = lineX_P1 + lineX_0 + lineX_M1;
        j[0] = (lineX_P1 - lineX_M1);
        j[1] = (lineY_P1 - lineY_M1);
    }

    static T compute_rho(Array<T, Descriptor::q> const &f)
    {
        return Descriptor::fullRho(get_rhoBar(f));
    }

    static void compute_uLb(Array<T, Descriptor::q> const &f, Array<T, 2> &uLb)
    {
        T rhoBar;
        get_rhoBar_j(f, rhoBar, uLb);
        T invRho = Descriptor::invRho(rhoBar);
        uLb[0] *= invRho;
        uLb[1] *= invRho;
    }

    static void compute_rho_uLb(Array<T, Descriptor::q> const &f, T &rho, Array<T, 2> &uLb)
    {
        T rhoBar;
        get_rhoBar_j(f, rhoBar, uLb);
        rho = Descriptor::fullRho(rhoBar);
        T invRho = Descriptor::invRho(rhoBar);
        uLb[0] *= invRho;
        uLb[1] *= invRho;
    }

    static T compute_e(Array<T, Descriptor::q> const &f)
    {
        return get_eBar(f) + (T)Descriptor::SkordosFactor() * Descriptor::d * Descriptor::cs2;
    }

    static T compute_rhoThetaBar(Array<T, Descriptor::q> const &f, T rhoBar, T jSqr)
    {
        T invRho = Descriptor::invRho(rhoBar);
        return Descriptor::invCs2 * Descriptor::invD * (get_eBar(f) - invRho * jSqr) - rhoBar;
    }

    static void compute_rho_rhoThetaBar(Array<T, Descriptor::q> const &f, T &rho, T &rhoThetaBar)
    {
        T j[Descriptor::d], rhoBar;
        get_rhoBar_j(f, rhoBar, j);
        T jSqr = VectorTemplateImpl<T, Descriptor::d>::normSqr(j);
        rho = Descriptor::fullRho(rhoBar);
        rhoThetaBar = compute_rhoThetaBar(f, rhoBar, jSqr);
    }

    static T compute_theta(Array<T, Descriptor::q> const &f, T rhoBar, T jSqr)
    {
        T invRho = Descriptor::invRho(rhoBar);
        T e = compute_e(f);
        return invRho * Descriptor::invD * Descriptor::invCs2 * (e - invRho * jSqr);
    }

    static T compute_rhoEpsilon(Array<T, Descriptor::q> const &f, T rhoBar, T jSqr)
    {
        T invRho = Descriptor::invRho(rhoBar);
        T e = compute_e(f);
        return (e - invRho * jSqr) / (T)2;
    }

    static void compute_PiNeq(
        Array<T, Descriptor::q> const &f, T rhoBar, Array<T, 2> const &j, Array<T, 3> &PiNeq)
    {
        T invRho = Descriptor::invRho(rhoBar);
        compute_PiNeq(f, rhoBar, j, PiNeq, invRho);
    }

    static void compute_PiNeq(
        Array<T, Descriptor::q> const &f, T rhoBar, Array<T, 2> const &j, Array<T, 3> &PiNeq,
        T invRho)
    {
        typedef SymmetricTensorImpl<T, 2> S;

        T lineX_P1, lineX_0, lineX_M1, lineY_P1, lineY_M1;
        partial_rho(f, lineX_P1, lineX_0, lineX_M1, lineY_P1, lineY_M1);

        PiNeq[S::xx] = lineX_P1 + lineX_M1 - Descriptor::cs2 * rhoBar - invRho * j[0] * j[0];
        PiNeq[S::yy] = lineY_P1 + lineY_M1 - Descriptor::cs2 * rhoBar - invRho * j[1] * j[1];
        PiNeq[S::xy] = -f[1] + f[3] - f[5] + f[7] - invRho * j[0] * j[1];
    }

    static void compute_thermal_PiNeq(
        Array<T, Descriptor::q> const &f, T rhoBar, T thetaBar, Array<T, 2> const &j,
        Array<T, 3> &PiNeq)
    {
        typedef SymmetricTensorImpl<T, 2> S;

        T rhoTheta_bar = rhoBar * thetaBar + rhoBar + Descriptor::SkordosFactor * thetaBar;
        T invRho = Descriptor::invRho(rhoBar);
        T lineX_P1, lineX_0, lineX_M1, lineY_P1, lineY_M1;
        partial_rho(f, lineX_P1, lineX_0, lineX_M1, lineY_P1, lineY_M1);

        PiNeq[S::xx] = lineX_P1 + lineX_M1 - Descriptor::cs2 * rhoTheta_bar - invRho * j[0] * j[0];
        PiNeq[S::yy] = lineY_P1 + lineY_M1 - Descriptor::cs2 * rhoTheta_bar - invRho * j[1] * j[1];
        PiNeq[S::xy] = -f[1] + f[3] - f[5] + f[7] - invRho * j[0] * j[1];
    }

    static void compute_rhoBar_j_PiNeq(
        Array<T, Descriptor::q> const &f, T &rhoBar, Array<T, 2> &j, Array<T, 3> &PiNeq)
    {
        typedef SymmetricTensorImpl<T, 2> S;

        T lineX_P1, lineX_0, lineX_M1, lineY_P1, lineY_M1;
        partial_rho(f, lineX_P1, lineX_0, lineX_M1, lineY_P1, lineY_M1);

        rhoBar = lineX_P1 + lineX_0 + lineX_M1;
        j[0] = lineX_P1 - lineX_M1;
        j[1] = lineY_P1 - lineY_M1;
        T invRho = Descriptor::invRho(rhoBar);
        PiNeq[S::xx] = lineX_P1 + lineX_M1 - Descriptor::cs2 * rhoBar - invRho * j[0] * j[0];
        PiNeq[S::yy] = lineY_P1 + lineY_M1 - Descriptor::cs2 * rhoBar - invRho * j[1] * j[1];
        PiNeq[S::xy] = -f[1] + f[3] - f[5] + f[7] - invRho * j[0] * j[1];
    }

    static void compute_rhoBar_j_PiNeq(
        Array<T, Descriptor::q> const &f, T &rhoBar, Array<T, 2> &j, Array<T, 3> &PiNeq, T invRho)
    {
        typedef SymmetricTensorImpl<T, 2> S;

        T lineX_P1, lineX_0, lineX_M1, lineY_P1, lineY_M1;
        partial_rho(f, lineX_P1, lineX_0, lineX_M1, lineY_P1, lineY_M1);

        rhoBar = lineX_P1 + lineX_0 + lineX_M1;
        j[0] = lineX_P1 - lineX_M1;
        j[1] = lineY_P1 - lineY_M1;
        PiNeq[S::xx] = lineX_P1 + lineX_M1 - Descriptor::cs2 * rhoBar - invRho * j[0] * j[0];
        PiNeq[S::yy] = lineY_P1 + lineY_M1 - Descriptor::cs2 * rhoBar - invRho * j[1] * j[1];
        PiNeq[S::xy] = -f[1] + f[3] - f[5] + f[7] - invRho * j[0] * j[1];
    }

    static void compute_rhoBar_thetaBar_j_PiNeq(
        Array<T, Descriptor::q> const &f, T &rhoBar, T &thetaBar, Array<T, 2> const &j,
        Array<T, 3> &PiNeq)
    {
        get_rhoBar_j(f, rhoBar, j);
        compute_PiNeq(f, rhoBar, j, PiNeq);
        T rhoThetaBar = compute_rhoThetaBar(f);
        thetaBar = rhoThetaBar * Descriptor::invRho(rhoBar);
        compute_thermal_PiNeq(f, rhoBar, thetaBar, j, PiNeq);
    }

    static void compute_P(
        Array<T, Descriptor::q> const &f, T rhoBar, Array<T, 2> const &j, Array<T, 3> &P)
    {
        typedef SymmetricTensorImpl<T, 2> S;

        T invRho = Descriptor::invRho(rhoBar);
        T lineX_P1, lineX_0, lineX_M1, lineY_P1, lineY_M1;
        partial_rho(f, lineX_P1, lineX_0, lineX_M1, lineY_P1, lineY_M1);

        P[S::xx] = lineX_P1 + lineX_M1 - invRho * j[0] * j[0];
        P[S::yy] = lineY_P1 + lineY_M1 - invRho * j[1] * j[1];
        P[S::xy] = -f[1] + f[3] - f[5] + f[7] - invRho * j[0] * j[1];
    }

    static void modifyJ(Array<T, Descriptor::q> &f, Array<T, 2> const &newJ)
    {
        T rhoBar;
        Array<T, 2> oldJ;
        get_rhoBar_j(f, rhoBar, oldJ);
        const T oldJSqr = VectorTemplateImpl<T, Descriptor::d>::normSqr(oldJ);
        const T newJSqr = VectorTemplateImpl<T, Descriptor::d>::normSqr(newJ);
        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            f[iPop] = f[iPop] - equilibrium(iPop, rhoBar, oldJ, oldJSqr)
                      + equilibrium(iPop, rhoBar, newJ, newJSqr);
        }
    }

};  // struct momentTemplatesImlp<D2Q9DescriptorBase>

}  // namespace plb

#endif  // MOMENT_TEMPLATES_2D_H
