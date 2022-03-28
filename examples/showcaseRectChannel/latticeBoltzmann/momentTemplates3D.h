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
 * 3D specialization of momentTemplates functions.
 */

#ifndef MOMENT_TEMPLATES_3D_H
#define MOMENT_TEMPLATES_3D_H

#include "core/globalDefs.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"

namespace plb {

// Efficient specialization for D3Q27 lattice
template <typename T>
struct momentTemplatesImpl<T, descriptors::D3Q27DescriptorBase<T> > {
    typedef descriptors::D3Q27DescriptorBase<T> Descriptor;

    static void partial_rho(
        Array<T, Descriptor::q> const &f, T &surfX_M1, T &surfX_P1, T &surfY_M1, T &surfY_P1,
        T &surfZ_M1, T &surfZ_P1)
    {
        surfX_M1 = f[1] + f[4] + f[5] + f[6] + f[7] + f[10] + f[11] + f[12] + f[13];
        surfX_P1 = f[14] + f[17] + f[18] + f[19] + f[20] + f[23] + f[24] + f[25] + f[26];

        surfY_M1 = f[2] + f[4] + f[8] + f[9] + f[10] + f[11] + f[18] + f[25] + f[26];
        surfY_P1 = f[5] + f[12] + f[13] + f[15] + f[17] + f[21] + f[22] + f[23] + f[24];

        surfZ_M1 = f[3] + f[6] + f[8] + f[10] + f[12] + f[20] + f[22] + f[24] + f[26];
        surfZ_P1 = f[7] + f[9] + f[11] + f[13] + f[16] + f[19] + f[21] + f[23] + f[25];
    }

    static void partial_rho(
        Array<T, Descriptor::q> const &f, T &surfX_M1, T &surfX_0, T &surfX_P1, T &surfY_M1,
        T &surfY_P1, T &surfZ_M1, T &surfZ_P1)
    {
        partial_rho(f, surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);
        surfX_0 = f[0] + f[2] + f[3] + f[8] + f[9] + f[15] + f[16] + f[21] + f[22];
    }

    static T get_rhoBar(Array<T, Descriptor::q> const &f)
    {
        T rhoBar = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[10]
                   + f[11] + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18] + f[19] + f[20]
                   + f[21] + f[22] + f[23] + f[24] + f[25] + f[26];
        return rhoBar;
    }

    static void get_j(Array<T, Descriptor::q> const &f, Array<T, 3> &j)
    {
        T surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

        surfX_M1 = f[1] + f[4] + f[5] + f[6] + f[7] + f[10] + f[11] + f[12] + f[13];
        surfX_P1 = f[14] + f[17] + f[18] + f[19] + f[20] + f[23] + f[24] + f[25] + f[26];

        surfY_M1 = f[2] + f[4] + f[8] + f[9] + f[10] + f[11] + f[18] + f[25] + f[26];
        surfY_P1 = f[5] + f[12] + f[13] + f[15] + f[17] + f[21] + f[22] + f[23] + f[24];

        surfZ_M1 = f[3] + f[6] + f[8] + f[10] + f[12] + f[20] + f[22] + f[24] + f[26];
        surfZ_P1 = f[7] + f[9] + f[11] + f[13] + f[16] + f[19] + f[21] + f[23] + f[25];

        j[0] = (surfX_P1 - surfX_M1);
        j[1] = (surfY_P1 - surfY_M1);
        j[2] = (surfZ_P1 - surfZ_M1);
    }

    static T get_eBar(Array<T, Descriptor::q> const &f)
    {
        T eBar = (T)3 * (f[10] + f[11] + f[12] + f[13] + f[23] + f[24] + f[25] + f[26])
                 + (T)2
                       * (f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[17] + f[18] + f[19] + f[20]
                          + f[21] + f[22])
                 + f[1] + f[2] + f[3] + f[14] + f[15] + f[16];
        return eBar;
    }

    static void get_rhoBar_j(Array<T, Descriptor::q> const &f, T &rhoBar, Array<T, 3> &j)
    {
        T surfX_M1, surfX_0, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

        partial_rho(f, surfX_M1, surfX_0, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

        rhoBar = surfX_M1 + surfX_0 + surfX_P1;

        j[0] = (surfX_P1 - surfX_M1);
        j[1] = (surfY_P1 - surfY_M1);
        j[2] = (surfZ_P1 - surfZ_M1);
    }

    static T compute_rho(Array<T, Descriptor::q> const &f)
    {
        return Descriptor::fullRho(get_rhoBar(f));
    }

    static void compute_uLb(Array<T, Descriptor::q> const &f, Array<T, 3> &uLb)
    {
        get_j(f, uLb);
        T invRho = Descriptor::invRho(get_rhoBar(f));
        uLb[0] *= invRho;
        uLb[1] *= invRho;
        uLb[2] *= invRho;
    }

    static void compute_rho_uLb(Array<T, Descriptor::q> const &f, T &rho, Array<T, 3> &uLb)
    {
        T rhoBar;
        get_rhoBar_j(f, rhoBar, uLb);
        T invRho = Descriptor::invRho(rhoBar);
        rho = Descriptor::fullRho(rhoBar);
        uLb[0] *= invRho;
        uLb[1] *= invRho;
        uLb[2] *= invRho;
    }

    static T compute_e(Array<T, Descriptor::q> const &f)
    {
        return get_eBar(f) + Descriptor::SkordosFactor * Descriptor::d * Descriptor::cs2;
    }

    static T compute_rhoThetaBar(Array<T, Descriptor::q> const &f, T rhoBar, T jSqr)
    {
        T invRho = Descriptor::invRho(rhoBar);
        return Descriptor::invCs2 * Descriptor::invD * (get_eBar(f) - invRho * jSqr) - rhoBar;
    }

    static void compute_rho_rhoThetaBar(Array<T, Descriptor::q> const &f, T &rho, T &rhoThetaBar)
    {
        T rhoBar, j[3];
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
        Array<T, Descriptor::q> const &f, T rhoBar, Array<T, 3> const &j, Array<T, 6> &PiNeq)
    {
        T invRho = Descriptor::invRho(rhoBar);
        compute_PiNeq(f, rhoBar, j, PiNeq, invRho);
    }

    static void compute_PiNeq(
        Array<T, Descriptor::q> const &f, T rhoBar, Array<T, 3> const &j, Array<T, 6> &PiNeq,
        T invRho)
    {
        typedef SymmetricTensorImpl<T, 3> S;

        T surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;
        partial_rho(f, surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

        PiNeq[S::xx] = surfX_P1 + surfX_M1 - Descriptor::cs2 * rhoBar - invRho * j[0] * j[0];
        PiNeq[S::yy] = surfY_P1 + surfY_M1 - Descriptor::cs2 * rhoBar - invRho * j[1] * j[1];
        PiNeq[S::zz] = surfZ_P1 + surfZ_M1 - Descriptor::cs2 * rhoBar - invRho * j[2] * j[2];

        PiNeq[S::xy] = f[4] - f[5] + f[10] + f[11] - f[12] - f[13] + f[17] - f[18] + f[23] + f[24]
                       - f[25] - f[26] - invRho * j[0] * j[1];
        PiNeq[S::xz] = f[6] - f[7] + f[10] - f[11] + f[12] - f[13] + f[19] - f[20] + f[23] - f[24]
                       + f[25] - f[26] - invRho * j[0] * j[2];
        PiNeq[S::yz] = f[8] - f[9] + f[10] - f[11] - f[12] + f[13] + f[21] - f[22] + f[23] - f[24]
                       - f[25] + f[26] - invRho * j[1] * j[2];
    }

    static void compute_thermal_PiNeq(
        Array<T, Descriptor::q> const &f, T rhoBar, T thetaBar, Array<T, 3> const &j,
        Array<T, 6> &PiNeq)
    {
        typedef SymmetricTensorImpl<T, 3> S;

        T rhoTheta_bar = rhoBar * thetaBar + rhoBar + Descriptor::SkordosFactor * thetaBar;
        T invRho = Descriptor::invRho(rhoBar);
        T surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;
        partial_rho(f, surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

        PiNeq[S::xx] = surfX_P1 + surfX_M1 - Descriptor::cs2 * rhoTheta_bar - invRho * j[0] * j[0];
        PiNeq[S::yy] = surfY_P1 + surfY_M1 - Descriptor::cs2 * rhoTheta_bar - invRho * j[1] * j[1];
        PiNeq[S::zz] = surfZ_P1 + surfZ_M1 - Descriptor::cs2 * rhoTheta_bar - invRho * j[2] * j[2];

        PiNeq[S::xy] = f[4] - f[5] + f[10] + f[11] - f[12] - f[13] + f[17] - f[18] + f[23] + f[24]
                       - f[25] - f[26] - invRho * j[0] * j[1];
        PiNeq[S::xz] = f[6] - f[7] + f[10] - f[11] + f[12] - f[13] + f[19] - f[20] + f[23] - f[24]
                       + f[25] - f[26] - invRho * j[0] * j[2];
        PiNeq[S::yz] = f[8] - f[9] + f[10] - f[11] - f[12] + f[13] + f[21] - f[22] + f[23] - f[24]
                       - f[25] + f[26] - invRho * j[1] * j[2];
    }

    static void compute_rhoBar_j_PiNeq(
        Array<T, Descriptor::q> const &f, T &rhoBar, Array<T, 3> &j, Array<T, 6> &PiNeq)
    {
        typedef SymmetricTensorImpl<T, 3> S;

        T surfX_M1, surfX_0, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

        partial_rho(f, surfX_M1, surfX_0, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);
        rhoBar = surfX_M1 + surfX_0 + surfX_P1;
        j[0] = (surfX_P1 - surfX_M1);
        j[1] = (surfY_P1 - surfY_M1);
        j[2] = (surfZ_P1 - surfZ_M1);
        T invRho = Descriptor::invRho(rhoBar);

        PiNeq[S::xx] = surfX_P1 + surfX_M1 - Descriptor::cs2 * rhoBar - invRho * j[0] * j[0];
        PiNeq[S::yy] = surfY_P1 + surfY_M1 - Descriptor::cs2 * rhoBar - invRho * j[1] * j[1];
        PiNeq[S::zz] = surfZ_P1 + surfZ_M1 - Descriptor::cs2 * rhoBar - invRho * j[2] * j[2];

        PiNeq[S::xy] = f[4] - f[5] + f[10] + f[11] - f[12] - f[13] + f[17] - f[18] + f[23] + f[24]
                       - f[25] - f[26] - invRho * j[0] * j[1];
        PiNeq[S::xz] = f[6] - f[7] + f[10] - f[11] + f[12] - f[13] + f[19] - f[20] + f[23] - f[24]
                       + f[25] - f[26] - invRho * j[0] * j[2];
        PiNeq[S::yz] = f[8] - f[9] + f[10] - f[11] - f[12] + f[13] + f[21] - f[22] + f[23] - f[24]
                       - f[25] + f[26] - invRho * j[1] * j[2];
    }

    static void compute_rhoBar_j_PiNeq(
        Array<T, Descriptor::q> const &f, T &rhoBar, Array<T, 3> &j, Array<T, 6> &PiNeq, T invRho)
    {
        typedef SymmetricTensorImpl<T, 3> S;

        T surfX_M1, surfX_0, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

        partial_rho(f, surfX_M1, surfX_0, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);
        rhoBar = surfX_M1 + surfX_0 + surfX_P1;
        j[0] = (surfX_P1 - surfX_M1);
        j[1] = (surfY_P1 - surfY_M1);
        j[2] = (surfZ_P1 - surfZ_M1);

        PiNeq[S::xx] = surfX_P1 + surfX_M1 - Descriptor::cs2 * rhoBar - invRho * j[0] * j[0];
        PiNeq[S::yy] = surfY_P1 + surfY_M1 - Descriptor::cs2 * rhoBar - invRho * j[1] * j[1];
        PiNeq[S::zz] = surfZ_P1 + surfZ_M1 - Descriptor::cs2 * rhoBar - invRho * j[2] * j[2];

        PiNeq[S::xy] = f[4] - f[5] + f[10] + f[11] - f[12] - f[13] + f[17] - f[18] + f[23] + f[24]
                       - f[25] - f[26] - invRho * j[0] * j[1];
        PiNeq[S::xz] = f[6] - f[7] + f[10] - f[11] + f[12] - f[13] + f[19] - f[20] + f[23] - f[24]
                       + f[25] - f[26] - invRho * j[0] * j[2];
        PiNeq[S::yz] = f[8] - f[9] + f[10] - f[11] - f[12] + f[13] + f[21] - f[22] + f[23] - f[24]
                       - f[25] + f[26] - invRho * j[1] * j[2];
    }

    static void compute_rhoBar_thetaBar_j_PiNeq(
        Array<T, Descriptor::q> const &f, T &rhoBar, T &thetaBar, Array<T, 3> const &j,
        Array<T, 6> &PiNeq)
    {
        get_rhoBar_j(f, rhoBar, j);
        compute_PiNeq(f, rhoBar, j, PiNeq);
        T rhoThetaBar = compute_rhoThetaBar(f);
        thetaBar = rhoThetaBar * Descriptor::invRho(rhoBar);
        compute_thermal_PiNeq(f, rhoBar, thetaBar, j, PiNeq);
    }

    static void compute_P(
        Array<T, Descriptor::q> const &f, T rhoBar, Array<T, 3> const &j, Array<T, 6> &P)
    {
        typedef SymmetricTensorImpl<T, 3> S;

        T invRho = Descriptor::invRho(rhoBar);
        T surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;
        partial_rho(f, surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

        P[S::xx] = surfX_P1 + surfX_M1 - invRho * j[0] * j[0];
        P[S::yy] = surfY_P1 + surfY_M1 - invRho * j[1] * j[1];
        P[S::zz] = surfZ_P1 + surfZ_M1 - invRho * j[2] * j[2];

        P[S::xy] = f[4] - f[5] + f[10] + f[11] - f[12] - f[13] + f[17] - f[18] + f[23] + f[24]
                   - f[25] - f[26] - invRho * j[0] * j[1];
        P[S::xz] = f[6] - f[7] + f[10] - f[11] + f[12] - f[13] + f[19] - f[20] + f[23] - f[24]
                   + f[25] - f[26] - invRho * j[0] * j[2];
        P[S::yz] = f[8] - f[9] + f[10] - f[11] - f[12] + f[13] + f[21] - f[22] + f[23] - f[24]
                   - f[25] + f[26] - invRho * j[1] * j[2];
    }

    static void modifyJ(Array<T, Descriptor::q> &f, Array<T, 3> const &newJ)
    {
        T rhoBar, oldJ[3];
        get_rhoBar_j(f, rhoBar, oldJ);
        const T oldJSqr = VectorTemplateImpl<T, 3>::normSqr(oldJ);
        const T newJSqr = VectorTemplateImpl<T, 3>::normSqr(newJ);
        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            f[iPop] = f[iPop] - equilibrium(iPop, rhoBar, oldJ, oldJSqr)
                      + equilibrium(iPop, rhoBar, newJ, newJSqr);
        }
    }

};  // struct momentTemplatesImpl<D3Q27DescriptorBase>

// Efficient specialization for D3Q19 lattice
template <typename T>
struct momentTemplatesImpl<T, descriptors::D3Q19DescriptorBase<T> > {
    typedef descriptors::D3Q19DescriptorBase<T> Descriptor;

    static void partial_rho(
        Array<T, Descriptor::q> const &f, T &surfX_M1, T &surfX_P1, T &surfY_M1, T &surfY_P1,
        T &surfZ_M1, T &surfZ_P1)
    {
        surfX_M1 = f[1] + f[4] + f[5] + f[6] + f[7];
        surfX_P1 = f[10] + f[13] + f[14] + f[15] + f[16];

        surfY_M1 = f[2] + f[4] + f[8] + f[9] + f[14];
        surfY_P1 = f[5] + f[11] + f[13] + f[17] + f[18];

        surfZ_M1 = f[3] + f[6] + f[8] + f[16] + f[18];
        surfZ_P1 = f[7] + f[9] + f[12] + f[15] + f[17];
    }

    static void partial_rho(
        Array<T, Descriptor::q> const &f, T &surfX_M1, T &surfX_0, T &surfX_P1, T &surfY_M1,
        T &surfY_P1, T &surfZ_M1, T &surfZ_P1)
    {
        partial_rho(f, surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);
        surfX_0 = f[0] + f[2] + f[3] + f[8] + f[9] + f[11] + f[12] + f[17] + f[18];
    }

    static T get_rhoBar(Array<T, Descriptor::q> const &f)
    {
        T rhoBar = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[10]
                   + f[11] + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18];
        return rhoBar;
    }

    static void get_j(Array<T, Descriptor::q> const &f, Array<T, 3> &j)
    {
        T surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

        surfX_M1 = f[1] + f[4] + f[5] + f[6] + f[7];
        surfX_P1 = f[10] + f[13] + f[14] + f[15] + f[16];

        surfY_M1 = f[2] + f[4] + f[8] + f[9] + f[14];
        surfY_P1 = f[5] + f[11] + f[13] + f[17] + f[18];

        surfZ_M1 = f[3] + f[6] + f[8] + f[16] + f[18];
        surfZ_P1 = f[7] + f[9] + f[12] + f[15] + f[17];

        j[0] = (surfX_P1 - surfX_M1);
        j[1] = (surfY_P1 - surfY_M1);
        j[2] = (surfZ_P1 - surfZ_M1);
    }

    static T get_eBar(Array<T, Descriptor::q> const &f)
    {
        T eBar = (T)2
                     * (f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[13] + f[14] + f[15] + f[16]
                        + f[17] + f[18])
                 + f[1] + f[2] + f[3] + f[10] + f[11] + f[12];
        return eBar;
    }

    static void get_rhoBar_j(Array<T, Descriptor::q> const &f, T &rhoBar, Array<T, 3> &j)
    {
        T surfX_M1, surfX_0, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

        partial_rho(f, surfX_M1, surfX_0, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

        rhoBar = surfX_M1 + surfX_0 + surfX_P1;

        j[0] = (surfX_P1 - surfX_M1);
        j[1] = (surfY_P1 - surfY_M1);
        j[2] = (surfZ_P1 - surfZ_M1);
    }

    static T compute_rho(Array<T, Descriptor::q> const &f)
    {
        return Descriptor::fullRho(get_rhoBar(f));
    }

    static void compute_uLb(Array<T, Descriptor::q> const &f, Array<T, 3> &uLb)
    {
        get_j(f, uLb);
        T invRho = Descriptor::invRho(get_rhoBar(f));
        uLb[0] *= invRho;
        uLb[1] *= invRho;
        uLb[2] *= invRho;
    }

    static void compute_rho_uLb(Array<T, Descriptor::q> const &f, T &rho, Array<T, 3> &uLb)
    {
        T rhoBar;
        get_rhoBar_j(f, rhoBar, uLb);
        T invRho = Descriptor::invRho(rhoBar);
        rho = Descriptor::fullRho(rhoBar);
        uLb[0] *= invRho;
        uLb[1] *= invRho;
        uLb[2] *= invRho;
    }

    static T compute_e(Array<T, Descriptor::q> const &f)
    {
        return get_eBar(f) + Descriptor::SkordosFactor * Descriptor::d * Descriptor::cs2;
    }

    static T compute_rhoThetaBar(Array<T, Descriptor::q> const &f, T rhoBar, T jSqr)
    {
        T invRho = Descriptor::invRho(rhoBar);
        return Descriptor::invCs2 * Descriptor::invD * (get_eBar(f) - invRho * jSqr) - rhoBar;
    }

    static void compute_rho_rhoThetaBar(Array<T, Descriptor::q> const &f, T &rho, T &rhoThetaBar)
    {
        T rhoBar, j[3];
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
        Array<T, Descriptor::q> const &f, T rhoBar, Array<T, 3> const &j, Array<T, 6> &PiNeq)
    {
        T invRho = Descriptor::invRho(rhoBar);
        compute_PiNeq(f, rhoBar, j, PiNeq, invRho);
    }

    static void compute_PiNeq(
        Array<T, Descriptor::q> const &f, T rhoBar, Array<T, 3> const &j, Array<T, 6> &PiNeq,
        T invRho)
    {
        typedef SymmetricTensorImpl<T, 3> S;

        T surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;
        partial_rho(f, surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

        PiNeq[S::xx] = surfX_P1 + surfX_M1 - Descriptor::cs2 * rhoBar - invRho * j[0] * j[0];
        PiNeq[S::yy] = surfY_P1 + surfY_M1 - Descriptor::cs2 * rhoBar - invRho * j[1] * j[1];
        PiNeq[S::zz] = surfZ_P1 + surfZ_M1 - Descriptor::cs2 * rhoBar - invRho * j[2] * j[2];

        PiNeq[S::xy] = f[4] - f[5] + f[13] - f[14] - invRho * j[0] * j[1];
        PiNeq[S::xz] = f[6] - f[7] + f[15] - f[16] - invRho * j[0] * j[2];
        PiNeq[S::yz] = f[8] - f[9] + f[17] - f[18] - invRho * j[1] * j[2];
    }

    static void compute_thermal_PiNeq(
        Array<T, Descriptor::q> const &f, T rhoBar, T thetaBar, Array<T, 3> const &j,
        Array<T, 6> &PiNeq)
    {
        typedef SymmetricTensorImpl<T, 3> S;

        T rhoTheta_bar = rhoBar * thetaBar + rhoBar + Descriptor::SkordosFactor * thetaBar;
        T invRho = Descriptor::invRho(rhoBar);
        T surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;
        partial_rho(f, surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

        PiNeq[S::xx] = surfX_P1 + surfX_M1 - Descriptor::cs2 * rhoTheta_bar - invRho * j[0] * j[0];
        PiNeq[S::yy] = surfY_P1 + surfY_M1 - Descriptor::cs2 * rhoTheta_bar - invRho * j[1] * j[1];
        PiNeq[S::zz] = surfZ_P1 + surfZ_M1 - Descriptor::cs2 * rhoTheta_bar - invRho * j[2] * j[2];

        PiNeq[S::xy] = f[4] - f[5] + f[13] - f[14] - invRho * j[0] * j[1];
        PiNeq[S::xz] = f[6] - f[7] + f[15] - f[16] - invRho * j[0] * j[2];
        PiNeq[S::yz] = f[8] - f[9] + f[17] - f[18] - invRho * j[1] * j[2];
    }

    static void compute_rhoBar_j_PiNeq(
        Array<T, Descriptor::q> const &f, T &rhoBar, Array<T, 3> &j, Array<T, 6> &PiNeq)
    {
        typedef SymmetricTensorImpl<T, 3> S;

        T surfX_M1, surfX_0, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

        partial_rho(f, surfX_M1, surfX_0, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);
        rhoBar = surfX_M1 + surfX_0 + surfX_P1;
        j[0] = (surfX_P1 - surfX_M1);
        j[1] = (surfY_P1 - surfY_M1);
        j[2] = (surfZ_P1 - surfZ_M1);
        T invRho = Descriptor::invRho(rhoBar);

        PiNeq[S::xx] = surfX_P1 + surfX_M1 - Descriptor::cs2 * rhoBar - invRho * j[0] * j[0];
        PiNeq[S::yy] = surfY_P1 + surfY_M1 - Descriptor::cs2 * rhoBar - invRho * j[1] * j[1];
        PiNeq[S::zz] = surfZ_P1 + surfZ_M1 - Descriptor::cs2 * rhoBar - invRho * j[2] * j[2];

        PiNeq[S::xy] = f[4] - f[5] + f[13] - f[14] - invRho * j[0] * j[1];
        PiNeq[S::xz] = f[6] - f[7] + f[15] - f[16] - invRho * j[0] * j[2];
        PiNeq[S::yz] = f[8] - f[9] + f[17] - f[18] - invRho * j[1] * j[2];
    }

    static void compute_rhoBar_j_PiNeq(
        Array<T, Descriptor::q> const &f, T &rhoBar, Array<T, 3> &j, Array<T, 6> &PiNeq, T invRho)
    {
        typedef SymmetricTensorImpl<T, 3> S;

        T surfX_M1, surfX_0, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

        partial_rho(f, surfX_M1, surfX_0, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);
        rhoBar = surfX_M1 + surfX_0 + surfX_P1;
        j[0] = (surfX_P1 - surfX_M1);
        j[1] = (surfY_P1 - surfY_M1);
        j[2] = (surfZ_P1 - surfZ_M1);

        PiNeq[S::xx] = surfX_P1 + surfX_M1 - Descriptor::cs2 * rhoBar - invRho * j[0] * j[0];
        PiNeq[S::yy] = surfY_P1 + surfY_M1 - Descriptor::cs2 * rhoBar - invRho * j[1] * j[1];
        PiNeq[S::zz] = surfZ_P1 + surfZ_M1 - Descriptor::cs2 * rhoBar - invRho * j[2] * j[2];

        PiNeq[S::xy] = f[4] - f[5] + f[13] - f[14] - invRho * j[0] * j[1];
        PiNeq[S::xz] = f[6] - f[7] + f[15] - f[16] - invRho * j[0] * j[2];
        PiNeq[S::yz] = f[8] - f[9] + f[17] - f[18] - invRho * j[1] * j[2];
    }

    static void compute_rhoBar_thetaBar_j_PiNeq(
        Array<T, Descriptor::q> const &f, T &rhoBar, T &thetaBar, Array<T, 3> const &j,
        Array<T, 6> &PiNeq)
    {
        get_rhoBar_j(f, rhoBar, j);
        compute_PiNeq(f, rhoBar, j, PiNeq);
        T rhoThetaBar = compute_rhoThetaBar(f);
        thetaBar = rhoThetaBar * Descriptor::invRho(rhoBar);
        compute_thermal_PiNeq(f, rhoBar, thetaBar, j, PiNeq);
    }

    static void compute_P(
        Array<T, Descriptor::q> const &f, T rhoBar, Array<T, 3> const &j, Array<T, 6> &P)
    {
        typedef SymmetricTensorImpl<T, 3> S;

        T invRho = Descriptor::invRho(rhoBar);
        T surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;
        partial_rho(f, surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

        P[S::xx] = surfX_P1 + surfX_M1 - invRho * j[0] * j[0];
        P[S::yy] = surfY_P1 + surfY_M1 - invRho * j[1] * j[1];
        P[S::zz] = surfZ_P1 + surfZ_M1 - invRho * j[2] * j[2];

        P[S::xy] = f[4] - f[5] + f[13] - f[14] - invRho * j[0] * j[1];
        P[S::xz] = f[6] - f[7] + f[15] - f[16] - invRho * j[0] * j[2];
        P[S::yz] = f[8] - f[9] + f[17] - f[18] - invRho * j[1] * j[2];
    }

    static void modifyJ(Array<T, Descriptor::q> &f, Array<T, 3> const &newJ)
    {
        T rhoBar, oldJ[3];
        get_rhoBar_j(f, rhoBar, oldJ);
        const T oldJSqr = VectorTemplateImpl<T, 3>::normSqr(oldJ);
        const T newJSqr = VectorTemplateImpl<T, 3>::normSqr(newJ);
        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            f[iPop] = f[iPop] - equilibrium(iPop, rhoBar, oldJ, oldJSqr)
                      + equilibrium(iPop, rhoBar, newJ, newJSqr);
        }
    }

};  // struct momentTemplatesImpl<D3Q19DescriptorBase>

// Efficient specialization for D3Q15 lattice
template <typename T>
struct momentTemplatesImpl<T, descriptors::D3Q15DescriptorBase<T> > {
    typedef descriptors::D3Q15DescriptorBase<T> Descriptor;

    static void partial_rho(
        Array<T, Descriptor::q> const &f, T &surfX_M1, T &surfX_P1, T &surfY_M1, T &surfY_P1,
        T &surfZ_M1, T &surfZ_P1)
    {
        surfX_M1 = f[1] + f[4] + f[5] + f[6] + f[7];
        surfX_P1 = f[8] + f[11] + f[12] + f[13] + f[14];

        surfY_M1 = f[2] + f[4] + f[5] + f[13] + f[14];
        surfY_P1 = f[6] + f[7] + f[9] + f[11] + f[12];

        surfZ_M1 = f[3] + f[4] + f[6] + f[12] + f[14];
        surfZ_P1 = f[5] + f[7] + f[10] + f[11] + f[13];
    }

    static void partial_rho(
        Array<T, Descriptor::q> const &f, T &surfX_M1, T &surfX_0, T &surfX_P1, T &surfY_M1,
        T &surfY_P1, T &surfZ_M1, T &surfZ_P1)
    {
        partial_rho(f, surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);
        surfX_0 = f[0] + f[2] + f[3] + f[9] + f[10];
    }

    static T get_rhoBar(Array<T, Descriptor::q> const &f)
    {
        T rhoBar = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[10]
                   + f[11] + f[12] + f[13] + f[14];
        return rhoBar;
    }

    static void get_j(Array<T, Descriptor::q> const &f, Array<T, 3> &j)
    {
        T surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

        surfX_M1 = f[1] + f[4] + f[5] + f[6] + f[7];
        surfX_P1 = f[8] + f[11] + f[12] + f[13] + f[14];

        surfY_M1 = f[2] + f[4] + f[5] + f[13] + f[14];
        surfY_P1 = f[6] + f[7] + f[9] + f[11] + f[12];

        surfZ_M1 = f[3] + f[4] + f[6] + f[12] + f[14];
        surfZ_P1 = f[5] + f[7] + f[10] + f[11] + f[13];

        j[0] = (surfX_P1 - surfX_M1);
        j[1] = (surfY_P1 - surfY_M1);
        j[2] = (surfZ_P1 - surfZ_M1);
    }

    static T get_eBar(Array<T, Descriptor::q> const &f)
    {
        T eBar = (T)3 * (f[4] + f[5] + f[6] + f[7] + f[11] + f[12] + f[13] + f[14]) + f[1] + f[2]
                 + f[3] + f[8] + f[9] + f[10];
        return eBar;
    }

    static void get_rhoBar_j(Array<T, Descriptor::q> const &f, T &rhoBar, Array<T, 3> &j)
    {
        T surfX_M1, surfX_0, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

        partial_rho(f, surfX_M1, surfX_0, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

        rhoBar = surfX_M1 + surfX_0 + surfX_P1;

        j[0] = (surfX_P1 - surfX_M1);
        j[1] = (surfY_P1 - surfY_M1);
        j[2] = (surfZ_P1 - surfZ_M1);
    }

    static T compute_rho(Array<T, Descriptor::q> const &f)
    {
        return Descriptor::fullRho(get_rhoBar(f));
    }

    static void compute_uLb(Array<T, Descriptor::q> const &f, Array<T, 3> &uLb)
    {
        get_j(f, uLb);
        T invRho = Descriptor::invRho(get_rhoBar(f));
        uLb[0] *= invRho;
        uLb[1] *= invRho;
        uLb[2] *= invRho;
    }

    static void compute_rho_uLb(Array<T, Descriptor::q> const &f, T &rho, Array<T, 3> &uLb)
    {
        T rhoBar;
        get_rhoBar_j(f, rhoBar, uLb);
        T invRho = Descriptor::invRho(rhoBar);
        rho = Descriptor::fullRho(rhoBar);
        uLb[0] *= invRho;
        uLb[1] *= invRho;
        uLb[2] *= invRho;
    }

    static T compute_e(Array<T, Descriptor::q> const &f)
    {
        return get_eBar(f) + Descriptor::SkordosFactor * Descriptor::d * Descriptor::cs2;
    }

    static T compute_rhoThetaBar(Array<T, Descriptor::q> const &f, T rhoBar, T jSqr)
    {
        T invRho = Descriptor::invRho(rhoBar);
        return Descriptor::invCs2 * Descriptor::invD * (get_eBar(f) - invRho * jSqr) - rhoBar;
    }

    static void compute_rho_rhoThetaBar(Array<T, Descriptor::q> const &f, T &rho, T &rhoThetaBar)
    {
        T rhoBar;
        Array<T, 3> j;
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
        Array<T, Descriptor::q> const &f, T rhoBar, Array<T, 3> const &j, Array<T, 6> &PiNeq)
    {
        T invRho = Descriptor::invRho(rhoBar);
        compute_PiNeq(f, rhoBar, j, PiNeq, invRho);
    }

    static void compute_PiNeq(
        Array<T, Descriptor::q> const &f, T rhoBar, Array<T, 3> const &j, Array<T, 6> &PiNeq,
        T invRho)
    {
        typedef SymmetricTensorImpl<T, 3> S;

        T surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;
        partial_rho(f, surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

        PiNeq[S::xx] = surfX_P1 + surfX_M1 - Descriptor::cs2 * rhoBar - invRho * j[0] * j[0];
        PiNeq[S::yy] = surfY_P1 + surfY_M1 - Descriptor::cs2 * rhoBar - invRho * j[1] * j[1];
        PiNeq[S::zz] = surfZ_P1 + surfZ_M1 - Descriptor::cs2 * rhoBar - invRho * j[2] * j[2];

        PiNeq[S::xy] =
            f[4] + f[5] - f[6] - f[7] + f[11] + f[12] - f[13] - f[14] - invRho * j[0] * j[1];
        PiNeq[S::xz] =
            f[4] - f[5] + f[6] - f[7] + f[11] - f[12] + f[13] - f[14] - invRho * j[0] * j[2];
        PiNeq[S::yz] =
            f[4] - f[5] - f[6] + f[7] + f[11] - f[12] - f[13] + f[14] - invRho * j[1] * j[2];
    }

    static void compute_thermal_PiNeq(
        Array<T, Descriptor::q> const &f, T rhoBar, T thetaBar, Array<T, 3> const &j,
        Array<T, 6> &PiNeq)
    {
        typedef SymmetricTensorImpl<T, 3> S;

        T rhoTheta_bar = rhoBar * thetaBar + rhoBar + Descriptor::SkordosFactor * thetaBar;
        T invRho = Descriptor::invRho(rhoBar);
        T surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;
        partial_rho(f, surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

        PiNeq[S::xx] = surfX_P1 + surfX_M1 - Descriptor::cs2 * rhoTheta_bar - invRho * j[0] * j[0];
        PiNeq[S::yy] = surfY_P1 + surfY_M1 - Descriptor::cs2 * rhoTheta_bar - invRho * j[1] * j[1];
        PiNeq[S::zz] = surfZ_P1 + surfZ_M1 - Descriptor::cs2 * rhoTheta_bar - invRho * j[2] * j[2];

        PiNeq[S::xy] =
            f[4] + f[5] - f[6] - f[7] + f[11] + f[12] - f[13] - f[14] - invRho * j[0] * j[1];
        PiNeq[S::xz] =
            f[4] - f[5] + f[6] - f[7] + f[11] - f[12] + f[13] - f[14] - invRho * j[0] * j[2];
        PiNeq[S::yz] =
            f[4] - f[5] - f[6] + f[7] + f[11] - f[12] - f[13] + f[14] - invRho * j[1] * j[2];
    }

    static void compute_rhoBar_j_PiNeq(
        Array<T, Descriptor::q> const &f, T &rhoBar, Array<T, 3> &j, Array<T, 6> &PiNeq)
    {
        typedef SymmetricTensorImpl<T, 3> S;

        T surfX_M1, surfX_0, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

        partial_rho(f, surfX_M1, surfX_0, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);
        rhoBar = surfX_M1 + surfX_0 + surfX_P1;
        j[0] = (surfX_P1 - surfX_M1);
        j[1] = (surfY_P1 - surfY_M1);
        j[2] = (surfZ_P1 - surfZ_M1);
        T invRho = Descriptor::invRho(rhoBar);

        PiNeq[S::xx] = surfX_P1 + surfX_M1 - Descriptor::cs2 * rhoBar - invRho * j[0] * j[0];
        PiNeq[S::yy] = surfY_P1 + surfY_M1 - Descriptor::cs2 * rhoBar - invRho * j[1] * j[1];
        PiNeq[S::zz] = surfZ_P1 + surfZ_M1 - Descriptor::cs2 * rhoBar - invRho * j[2] * j[2];

        PiNeq[S::xy] =
            f[4] + f[5] - f[6] - f[7] + f[11] + f[12] - f[13] - f[14] - invRho * j[0] * j[1];
        PiNeq[S::xz] =
            f[4] - f[5] + f[6] - f[7] + f[11] - f[12] + f[13] - f[14] - invRho * j[0] * j[2];
        PiNeq[S::yz] =
            f[4] - f[5] - f[6] + f[7] + f[11] - f[12] - f[13] + f[14] - invRho * j[1] * j[2];
    }

    static void compute_rhoBar_j_PiNeq(
        Array<T, Descriptor::q> const &f, T &rhoBar, Array<T, 3> &j, Array<T, 6> &PiNeq, T invRho)
    {
        typedef SymmetricTensorImpl<T, 3> S;

        T surfX_M1, surfX_0, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

        partial_rho(f, surfX_M1, surfX_0, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);
        rhoBar = surfX_M1 + surfX_0 + surfX_P1;
        j[0] = (surfX_P1 - surfX_M1);
        j[1] = (surfY_P1 - surfY_M1);
        j[2] = (surfZ_P1 - surfZ_M1);

        PiNeq[S::xx] = surfX_P1 + surfX_M1 - Descriptor::cs2 * rhoBar - invRho * j[0] * j[0];
        PiNeq[S::yy] = surfY_P1 + surfY_M1 - Descriptor::cs2 * rhoBar - invRho * j[1] * j[1];
        PiNeq[S::zz] = surfZ_P1 + surfZ_M1 - Descriptor::cs2 * rhoBar - invRho * j[2] * j[2];

        PiNeq[S::xy] =
            f[4] + f[5] - f[6] - f[7] + f[11] + f[12] - f[13] - f[14] - invRho * j[0] * j[1];
        PiNeq[S::xz] =
            f[4] - f[5] + f[6] - f[7] + f[11] - f[12] + f[13] - f[14] - invRho * j[0] * j[2];
        PiNeq[S::yz] =
            f[4] - f[5] - f[6] + f[7] + f[11] - f[12] - f[13] + f[14] - invRho * j[1] * j[2];
    }

    static void compute_rhoBar_thetaBar_j_PiNeq(
        Array<T, Descriptor::q> const &f, T &rhoBar, T &thetaBar, Array<T, 3> const &j,
        Array<T, 6> &PiNeq)
    {
        get_rhoBar_j(f, rhoBar, j);
        compute_PiNeq(f, rhoBar, j, PiNeq);
        T rhoThetaBar = compute_rhoThetaBar(f);
        thetaBar = rhoThetaBar * Descriptor::invRho(rhoBar);
        compute_thermal_PiNeq(f, rhoBar, thetaBar, j, PiNeq);
    }

    static void compute_P(
        Array<T, Descriptor::q> const &f, T rhoBar, Array<T, 3> const &j, Array<T, 6> &P)
    {
        typedef SymmetricTensorImpl<T, 3> S;

        T invRho = Descriptor::invRho(rhoBar);
        T surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;
        partial_rho(f, surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

        P[S::xx] = surfX_P1 + surfX_M1 - invRho * j[0] * j[0];
        P[S::yy] = surfY_P1 + surfY_M1 - invRho * j[1] * j[1];
        P[S::zz] = surfZ_P1 + surfZ_M1 - invRho * j[2] * j[2];

        P[S::xy] = f[4] + f[5] - f[6] - f[7] + f[11] + f[12] - f[13] - f[14] - invRho * j[0] * j[1];
        P[S::xz] = f[4] - f[5] + f[6] - f[7] + f[11] - f[12] + f[13] - f[14] - invRho * j[0] * j[2];
        P[S::yz] = f[4] - f[5] - f[6] + f[7] + f[11] - f[12] - f[13] + f[14] - invRho * j[1] * j[2];
    }

    static void modifyJ(Array<T, Descriptor::q> &f, Array<T, 3> const &newJ)
    {
        T rhoBar, oldJ[3];
        get_rhoBar_j(f, rhoBar, oldJ);
        const T oldJSqr = VectorTemplateImpl<T, 3>::normSqr(oldJ);
        const T newJSqr = VectorTemplateImpl<T, 3>::normSqr(newJ);
        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            f[iPop] = f[iPop] - equilibrium(iPop, rhoBar, oldJ, oldJSqr)
                      + equilibrium(iPop, rhoBar, newJ, newJSqr);
        }
    }

};  // struct momentTemplatesImpl<D3Q15DescriptorBase>

}  // namespace plb

#endif
