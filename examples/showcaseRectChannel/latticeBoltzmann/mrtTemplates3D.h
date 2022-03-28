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

/* Orestis Malaspinas contributed this code.
 */

/** \file
 * Template specializations for some computationally intensive LB
 * functions of the header file mrtTemplates.h, for the D3Q19 grid.
 */

#ifndef MRT_TEMPLATES_3D_H
#define MRT_TEMPLATES_3D_H

#include "core/globalDefs.h"

namespace plb {

// Efficient specialization for D3Q19 lattice
template <typename T>
struct mrtTemplatesImpl<T, descriptors::MRTD3Q19DescriptorBase<T> > {
    typedef descriptors::D3Q19DescriptorBase<T> Descriptor;
    typedef descriptors::MRTD3Q19DescriptorBase<T> MRTDescriptor;

    /// Computation of all equilibrium distribution (in moments space)
    static void computeEquilibriumMoments(
        Array<T, Descriptor::q> &momentsEq, T rhoBar, Array<T, 3> const &j, T jSqr)
    {
        T invRho = Descriptor::invRho(rhoBar);

        momentsEq[0] = rhoBar;
        momentsEq[1] = (T)19 * jSqr * invRho - (T)11 * rhoBar;
        momentsEq[2] = -(T)5.5 * jSqr * invRho + (T)3 * rhoBar;
        momentsEq[3] = j[0];
        momentsEq[4] = -((T)2 / (T)3) * j[0];
        momentsEq[5] = j[1];
        momentsEq[6] = -((T)2 / 3) * j[1];
        momentsEq[7] = j[2];
        momentsEq[8] = -((T)2 / (T)3) * j[2];
        momentsEq[9] = ((T)2 * j[0] * j[0] - j[1] * j[1] - j[2] * j[2]) * invRho;
        momentsEq[10] = (-j[0] * j[0] + (T)0.5 * j[1] * j[1] + (T)0.5 * j[2] * j[2]) * invRho;
        momentsEq[11] = (j[1] * j[1] - j[2] * j[2]) * invRho;
        momentsEq[12] = (-(T)0.5 * j[1] * j[1] + (T)0.5 * j[2] * j[2]) * invRho;
        momentsEq[13] = j[1] * j[0] * invRho;
        momentsEq[14] = j[2] * j[1] * invRho;
        momentsEq[15] = j[2] * j[0] * invRho;
        momentsEq[16] = T();
        momentsEq[17] = T();
        momentsEq[18] = T();
    }

    /// Computation of all equilibrium distribution (in moments space)
    static void computeSmagorinskyEquilibriumMoments(
        Array<T, Descriptor::q> &momentsEq, T rhoBar, Array<T, 3> const &j, T jSqr,
        const Array<T, 6> &strain, T cSmago)
    {
        typedef SymmetricTensorImpl<T, 3> S;

        T invRho = Descriptor::invRho(rhoBar);
        T rho = Descriptor::fullRho(rhoBar);
        T rho2 = rho * rho;

        T sNorm = std::sqrt((T)2 * SymmetricTensorImpl<T, 3>::tensorNormSqr(strain));
        T smagoFactor = (T)2 * cSmago * cSmago * sNorm;

        T ux2 = rho2 * smagoFactor * strain[S::xx];
        T uy2 = rho2 * smagoFactor * strain[S::yy];
        T uz2 = rho2 * smagoFactor * strain[S::zz];

        T uxuy = rho2 * smagoFactor * strain[S::xy];
        T uyuz = rho2 * smagoFactor * strain[S::yz];
        T uxuz = rho2 * smagoFactor * strain[S::xz];

        momentsEq[0] = rhoBar;
        momentsEq[1] = (T)19 * (jSqr + ux2 + uy2 + uz2) * invRho - (T)11 * rhoBar;
        momentsEq[2] = -(T)5.5 * (jSqr + ux2 + uy2 + uz2) * invRho + (T)3 * rhoBar;
        momentsEq[3] = j[0];
        momentsEq[4] = -((T)2 / (T)3) * j[0];
        momentsEq[5] = j[1];
        momentsEq[6] = -((T)2 / 3) * j[1];
        momentsEq[7] = j[2];
        momentsEq[8] = -((T)2 / (T)3) * j[2];
        momentsEq[9] =
            ((T)2 * (j[0] * j[0] + ux2) - (j[1] * j[1] + uy2) - (j[2] * j[2] + uz2)) * invRho;
        momentsEq[10] =
            (-(j[0] * j[0] + ux2) + (T)0.5 * (j[1] * j[1] + uy2) + (T)0.5 * (j[2] * j[2] + uz2))
            * invRho;
        momentsEq[11] = (j[1] * j[1] + uy2 - (j[2] * j[2] + uz2)) * invRho;
        momentsEq[12] = (-(T)0.5 * (j[1] * j[1] + uy2) + (T)0.5 * (j[2] * j[2] + uz2)) * invRho;
        momentsEq[13] = (j[1] * j[0] + uxuy) * invRho;
        momentsEq[14] = (j[2] * j[1] + uyuz) * invRho;
        momentsEq[15] = (j[2] * j[0] + uxuz) * invRho;
        momentsEq[16] = T();
        momentsEq[17] = T();
        momentsEq[18] = T();
    }

    /// Computation of all moments (specialized for d3q19)
    static void computeMoments(Array<T, Descriptor::q> &moments, const Array<T, Descriptor::q> &f)
    {
        moments[0] = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[10]
                     + f[11] + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18];
        moments[1] = -(T)30 * f[0] - (T)11 * (f[1] + f[2] + f[3] + f[10] + f[11] + f[12])
                     + (T)8
                           * (f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[13] + f[14] + f[15]
                              + f[16] + f[17] + f[18]);
        moments[2] = 12 * f[0] - 4 * f[1] - 4 * f[2] - 4 * f[3] + f[4] + f[5] + f[6] + f[7] + f[8]
                     + f[9] - 4 * f[10] - 4 * f[11] - 4 * f[12] + f[13] + f[14] + f[15] + f[16]
                     + f[17] + f[18];
        moments[3] = -f[1] - f[4] - f[5] - f[6] - f[7] + f[10] + f[13] + f[14] + f[15] + f[16];
        moments[4] =
            4 * f[1] - f[4] - f[5] - f[6] - f[7] - 4 * f[10] + f[13] + f[14] + f[15] + f[16];
        moments[5] = -f[2] - f[4] + f[5] - f[8] - f[9] + f[11] + f[13] - f[14] + f[17] + f[18];
        moments[6] =
            4 * f[2] - f[4] + f[5] - f[8] - f[9] - 4 * f[11] + f[13] - f[14] + f[17] + f[18];
        moments[7] = -f[3] - f[6] + f[7] - f[8] + f[9] + f[12] + f[15] - f[16] + f[17] - f[18];
        moments[8] =
            4 * f[3] - f[6] + f[7] - f[8] + f[9] - 4 * f[12] + f[15] - f[16] + f[17] - f[18];
        moments[9] = 2 * f[1] - f[2] - f[3] + f[4] + f[5] + f[6] + f[7] - 2 * f[8] - 2 * f[9]
                     + 2 * f[10] - f[11] - f[12] + f[13] + f[14] + f[15] + f[16] - 2 * f[17]
                     - 2 * f[18];
        moments[10] = -4 * f[1] + 2 * f[2] + 2 * f[3] + f[4] + f[5] + f[6] + f[7] - 2 * f[8]
                      - 2 * f[9] - 4 * f[10] + 2 * f[11] + 2 * f[12] + f[13] + f[14] + f[15] + f[16]
                      - 2 * f[17] - 2 * f[18];
        moments[11] =
            f[2] - f[3] + f[4] + f[5] - f[6] - f[7] + f[11] - f[12] + f[13] + f[14] - f[15] - f[16];
        moments[12] = -2 * f[2] + 2 * f[3] + f[4] + f[5] - f[6] - f[7] - 2 * f[11] + 2 * f[12]
                      + f[13] + f[14] - f[15] - f[16];
        moments[13] = f[4] - f[5] + f[13] - f[14];
        moments[14] = f[8] - f[9] + f[17] - f[18];
        moments[15] = f[6] - f[7] + f[15] - f[16];
        moments[16] = -f[4] - f[5] + f[6] + f[7] + f[13] + f[14] - f[15] - f[16];
        moments[17] = f[4] - f[5] - f[8] - f[9] - f[13] + f[14] + f[17] + f[18];
        moments[18] = -f[6] + f[7] + f[8] - f[9] + f[15] - f[16] - f[17] + f[18];
    }

    static void computef_InvM_Smoments(Array<T, 19> &f, const Array<T, 19> &moments, const T &omega)
    {
        T mom0 = moments[0] * MRTDescriptor::S[0];
        T mom1 = moments[1] * MRTDescriptor::S[1];
        T mom2 = moments[2] * MRTDescriptor::S[2];
        T mom3 = moments[3] * MRTDescriptor::S[3];
        T mom4 = moments[4] * MRTDescriptor::S[4];
        T mom5 = moments[5] * MRTDescriptor::S[5];
        T mom6 = moments[6] * MRTDescriptor::S[6];
        T mom7 = moments[7] * MRTDescriptor::S[7];
        T mom8 = moments[8] * MRTDescriptor::S[8];
        T mom9 = moments[9] * omega;
        T mom10 = moments[10] * MRTDescriptor::S[10];
        T mom11 = moments[11] * omega;
        T mom12 = moments[12] * MRTDescriptor::S[12];
        T mom13 = moments[13] * omega;
        T mom14 = moments[14] * omega;
        T mom15 = moments[15] * omega;
        T mom16 = moments[16] * MRTDescriptor::S[16];
        T mom17 = moments[17] * MRTDescriptor::S[17];
        T mom18 = moments[18] * MRTDescriptor::S[18];

        T mom0tmp = mom0 / (T)19;

        f[0] -= mom0tmp - 5 / (T)399 * mom1 + 1 / (T)21 * mom2;

        T mom1tmp = (T)11 / (T)2394 * mom1;
        T mom2tmp = mom2 / (T)63;
        T mom3_m4 = (T)0.1 * (mom3 - mom4);
        T mom9_m10 = (mom9 - mom10) / (T)18;
        f[1] -= mom0tmp - mom1tmp - mom2tmp - mom3_m4 + mom9_m10;
        f[10] -= mom0tmp - mom1tmp - mom2tmp + mom3_m4 + mom9_m10;

        T mom5_m6 = (T)0.1 * (mom5 - mom6);
        T mom7_m8 = (T)0.1 * (mom7 - mom8);
        mom9_m10 *= (T)0.5;
        T mom11_m12 = (mom11 - mom12) / (T)12;
        f[2] -= mom0tmp - mom1tmp - mom2tmp - mom5_m6 - mom9_m10 + mom11_m12;
        f[3] -= mom0tmp - mom1tmp - mom2tmp - mom7_m8 - mom9_m10 - mom11_m12;
        f[11] -= mom0tmp - mom1tmp - mom2tmp + mom5_m6 - mom9_m10 + mom11_m12;
        f[12] -= mom0tmp - mom1tmp - mom2tmp + mom7_m8 - mom9_m10 - mom11_m12;

        mom1tmp = (T)4 / (T)1197 * mom1;
        mom2tmp *= (T)0.25;
        T mom5_p6 = (T)0.1 * mom5 + (T)0.025 * mom6;
        T mom7_p8 = (T)0.1 * mom7 + (T)0.025 * mom8;
        T mom9_p10 = (mom9 + mom10 * (T)0.5) / (T)18;
        mom14 *= (T)0.25;
        mom17 *= (T)0.125;
        mom18 *= (T)0.125;
        f[8] -= mom0tmp + mom1tmp + mom2tmp - mom5_p6 - mom7_p8 - mom9_p10 + mom14 - mom17 + mom18;
        f[9] -= mom0tmp + mom1tmp + mom2tmp - mom5_p6 + mom7_p8 - mom9_p10 - mom14 - mom17 - mom18;
        f[17] -= mom0tmp + mom1tmp + mom2tmp + mom5_p6 + mom7_p8 - mom9_p10 + mom14 + mom17 - mom18;
        f[18] -= mom0tmp + mom1tmp + mom2tmp + mom5_p6 - mom7_p8 - mom9_p10 - mom14 + mom17 + mom18;

        T mom3_p4 = (T)0.1 * mom3 + (T)0.025 * mom4;
        mom9_p10 *= (T)0.5;
        T mom11_p12 = (mom11 + (T)0.5 * mom12) / (T)12;
        mom13 *= (T)0.25;
        mom16 *= (T)0.125;
        f[4] -= mom0tmp + mom1tmp + mom2tmp - mom3_p4 - mom5_p6 + mom9_p10 + mom11_p12 + mom13
                - mom16 + mom17;
        f[5] -= mom0tmp + mom1tmp + mom2tmp - mom3_p4 + mom5_p6 + mom9_p10 + mom11_p12 - mom13
                - mom16 - mom17;
        f[13] -= mom0tmp + mom1tmp + mom2tmp + mom3_p4 + mom5_p6 + mom9_p10 + mom11_p12 + mom13
                 + mom16 - mom17;
        f[14] -= mom0tmp + mom1tmp + mom2tmp + mom3_p4 - mom5_p6 + mom9_p10 + mom11_p12 - mom13
                 + mom16 + mom17;

        mom15 *= (T)0.25;
        f[15] -= mom0tmp + mom1tmp + mom2tmp + mom3_p4 + mom7_p8 + mom9_p10 - mom11_p12 + mom15
                 - mom16 + mom18;
        f[16] -= mom0tmp + mom1tmp + mom2tmp + mom3_p4 - mom7_p8 + mom9_p10 - mom11_p12 - mom15
                 - mom16 - mom18;
        f[6] -= mom0tmp + mom1tmp + mom2tmp - mom3_p4 - mom7_p8 + mom9_p10 - mom11_p12 + mom15
                + mom16 - mom18;
        f[7] -= mom0tmp + mom1tmp + mom2tmp - mom3_p4 + mom7_p8 + mom9_p10 - mom11_p12 - mom15
                + mom16 + mom18;
    }

    static void computeMneqInPlace(Array<T, 19> &moments, const Array<T, 19> &momentsEq)
    {
        moments[0] -= momentsEq[0];
        moments[1] -= momentsEq[1];
        moments[2] -= momentsEq[2];
        moments[3] -= momentsEq[3];
        moments[4] -= momentsEq[4];
        moments[5] -= momentsEq[5];
        moments[6] -= momentsEq[6];
        moments[7] -= momentsEq[7];
        moments[8] -= momentsEq[8];
        moments[9] -= momentsEq[9];
        moments[10] -= momentsEq[10];
        moments[11] -= momentsEq[11];
        moments[12] -= momentsEq[12];
        moments[13] -= momentsEq[13];
        moments[14] -= momentsEq[14];
        moments[15] -= momentsEq[15];
        moments[16] -= momentsEq[16];
        moments[17] -= momentsEq[17];
        moments[18] -= momentsEq[18];
    }

    /// MRT collision step
    static T mrtCollision(Array<T, Descriptor::q> &f, const T &omega)
    {
        Array<T, 19> moments, momentsEq;

        computeMoments(moments, f);
        T rhoBar = moments[0];
        Array<T, 3> j(
            moments[MRTDescriptor::momentumIndexes[0]], moments[MRTDescriptor::momentumIndexes[1]],
            moments[MRTDescriptor::momentumIndexes[2]]);
        T jSqr = VectorTemplateImpl<T, 3>::normSqr(j);

        computeEquilibriumMoments(momentsEq, rhoBar, j, jSqr);
        computeMneqInPlace(moments, momentsEq);  // moments become mNeq
        computef_InvM_Smoments(f, moments, omega);

        return jSqr;
    }

    /// MRT collision step
    static T mrtCollision(
        Array<T, Descriptor::q> &f, const T &rhoBar, const Array<T, 3> &j, const T omega)
    {
        Array<T, 19> moments, momentsEq;

        computeMoments(moments, f);
        T jSqr = VectorTemplateImpl<T, 3>::normSqr(j);
        computeEquilibriumMoments(momentsEq, rhoBar, j, jSqr);
        computeMneqInPlace(moments, momentsEq);  // moments become mNeq
        computef_InvM_Smoments(f, moments, omega);

        return jSqr;
    }

    /// MRT collision step
    static T smagorinskyMrtCollision(
        Array<T, Descriptor::q> &f, const T &rhoBar, const Array<T, 3> &j, const T &omega,
        const Array<T, 6> &strain, T cSmago)
    {
        Array<T, 19> moments, momentsEq;
        computeMoments(moments, f);
        T jSqr = VectorTemplateImpl<T, 3>::normSqr(j);

        computeSmagorinskyEquilibriumMoments(momentsEq, rhoBar, j, jSqr, strain, cSmago);
        computeMneqInPlace(moments, momentsEq);  // moments become mNeq
        computef_InvM_Smoments(f, moments, omega);

        return jSqr;
    }

    /// MRT collision step
    static T smagorinskyMrtCollision(
        Array<T, Descriptor::q> &f, const T &omega, const Array<T, 6> &strain, T cSmago)
    {
        Array<T, 19> moments, momentsEq;
        computeMoments(moments, f);
        T rhoBar = moments[0];
        Array<T, 3> j(
            moments[MRTDescriptor::momentumIndexes[0]], moments[MRTDescriptor::momentumIndexes[1]],
            moments[MRTDescriptor::momentumIndexes[2]]);
        T jSqr = VectorTemplateImpl<T, 3>::normSqr(j);

        computeSmagorinskyEquilibriumMoments(momentsEq, rhoBar, j, jSqr, strain, cSmago);
        computeMneqInPlace(moments, momentsEq);  // moments become mNeq
        computef_InvM_Smoments(f, moments, omega);

        return jSqr;
    }

    static void addGuoForce(
        Array<T, Descriptor::q> &f, const Array<T, Descriptor::d> &g,
        Array<T, Descriptor::d> const &u, const T &omega, T amplitude)
    {
        Array<T, Descriptor::q> forcing, momForce;
        T gx = amplitude * g[0];
        T gy = amplitude * g[1];
        T gz = amplitude * g[2];

        T g_u = gx * u[0] + gy * u[1] + gz * u[2];

        momForce[0] = 0;
        momForce[1] = 38 * g_u;
        momForce[2] = -11 * g_u;
        momForce[3] = gx;
        momForce[4] = -(T)2 / 3 * gx;
        momForce[5] = gy;
        momForce[6] = -(T)2 / 3 * gy;
        momForce[7] = gz;
        momForce[8] = -(T)2 / 3 * gz;
        momForce[9] = 4 * gx * u[0] - 2 * gy * u[1] - 2 * gz * u[2];
        momForce[10] = -2 * gx * u[0] + gy * u[1] + gz * u[2];
        momForce[11] = 2 * gy * u[1] - 2 * gz * u[2];
        momForce[12] = -gy * u[1] + gz * u[2];
        momForce[13] = gx * u[1] + gy * u[0];
        momForce[14] = gy * u[2] + gz * u[1];
        momForce[15] = gx * u[2] + gz * u[0];
        momForce[16] = 0;
        momForce[17] = 0;
        momForce[18] = 0;

        momForce[0] *= (T)0.5;
        momForce[1] *= (T)0.5;
        momForce[2] *= (T)0.5;
        momForce[3] *= (T)0.5;
        momForce[4] *= (T)0.5;
        momForce[5] *= (T)0.5;
        momForce[6] *= (T)0.5;
        momForce[7] *= (T)0.5;
        momForce[8] *= (T)0.5;
        momForce[9] *= (T)0.5;
        momForce[10] *= (T)0.5;
        momForce[11] *= (T)0.5;
        momForce[12] *= (T)0.5;
        momForce[13] *= (T)0.5;
        momForce[14] *= (T)0.5;
        momForce[15] *= (T)0.5;
        momForce[16] *= (T)0.5;
        momForce[17] *= (T)0.5;
        momForce[18] *= (T)0.5;

        computef_InvM_Smoments(f, momForce, omega);

        static const T oneOver6 = (T)1 / (T)6;
        static const T oneOver12 = (T)1 / (T)12;

        f[0] += -g_u;

        f[1] += oneOver6 * (gx * (-(T)1 + 2 * u[0]) - gy * u[1] - gz * u[2]);
        f[10] += oneOver6 * (gx * ((T)1 + 2 * u[0]) - gy * u[1] - gz * u[2]);

        f[2] += -oneOver6 * (gx * u[0] + gy * ((T)1 - 2 * u[1]) + gz * u[2]);
        f[3] += -oneOver6 * (gx * u[0] + gy * u[1] + gz * ((T)1 - 2 * u[2]));
        f[11] += -oneOver6 * (gx * u[0] + gy * (-(T)1 - 2 * u[1]) + gz * u[2]);
        f[12] += -oneOver6 * (gx * u[0] + gy * u[1] + gz * (-(T)1 - 2 * u[2]));

        f[4] +=
            oneOver12
            * (gx * (-(T)1 + 2 * u[0] + 3 * u[1]) + gy * (-(T)1 + 2 * u[1] + 3 * u[0]) - gz * u[2]);
        f[5] +=
            oneOver12
            * (gx * (-(T)1 + 2 * u[0] - 3 * u[1]) + gy * ((T)1 + 2 * u[1] - 3 * u[0]) - gz * u[2]);
        f[6] +=
            oneOver12
            * (gx * (-(T)1 + 2 * u[0] + 3 * u[2]) - gy * u[1] + gz * (-(T)1 + 2 * u[2] + 3 * u[0]));
        f[7] +=
            oneOver12
            * (gx * (-(T)1 + 2 * u[0] - 3 * u[2]) - gy * u[1] + gz * ((T)1 + 2 * u[2] - 3 * u[0]));
        f[8] +=
            -oneOver12
            * (gx * u[0] + gy * ((T)1 - 2 * u[1] - 3 * u[2]) + gz * ((T)1 - 2 * u[2] - 3 * u[1]));
        f[9] +=
            -oneOver12
            * (gx * u[0] + gy * ((T)1 - 2 * u[1] + 3 * u[2]) + gz * (-(T)1 - 2 * u[2] + 3 * u[1]));
        f[13] +=
            oneOver12
            * (gx * ((T)1 + 2 * u[0] + 3 * u[1]) + gy * ((T)1 + 2 * u[1] + 3 * u[0]) - gz * u[2]);
        f[14] +=
            oneOver12
            * (gx * ((T)1 + 2 * u[0] - 3 * u[1]) + gy * (-(T)1 + 2 * u[1] - 3 * u[0]) - gz * u[2]);
        f[15] +=
            oneOver12
            * (gx * ((T)1 + 2 * u[0] + 3 * u[2]) - gy * u[1] + gz * ((T)1 + 2 * u[2] + 3 * u[0]));
        f[16] +=
            oneOver12
            * (gx * ((T)1 + 2 * u[0] - 3 * u[2]) - gy * u[1] + gz * (-(T)1 + 2 * u[2] - 3 * u[0]));
        f[17] +=
            -oneOver12
            * (gx * u[0] + gy * (-(T)1 - 2 * u[1] - 3 * u[2]) + gz * (-(T)1 - 2 * u[2] - 3 * u[1]));
        f[18] +=
            -oneOver12
            * (gx * u[0] + gy * (-(T)1 - 2 * u[1] + 3 * u[2]) + gz * ((T)1 - 2 * u[2] + 3 * u[1]));
    }

    /// MRT collision step
    static T mrtCollisionWithForce(
        Array<T, Descriptor::q> &f, const T &rhoBar, const Array<T, Descriptor::d> &u,
        const T &omega, const Array<T, Descriptor::d> &force, T amplitude)
    {
        Array<T, Descriptor::d> j = Descriptor::fullRho(rhoBar) * u;
        T jSqr = mrtCollision(f, rhoBar, j, omega);
        addGuoForce(f, force, u, omega, amplitude);

        return jSqr;
    }

    static void addHeForce(
        Array<T, Descriptor::q> &f, const Array<T, Descriptor::d> &force, const T &rhoBar,
        Array<T, Descriptor::d> const &uLB, const T &omega, T amplitude)
    {
        ///////////// new
        // NOW WE CALCULATE THE moments of the forcing (HE forcing)

        T rhoFull = Descriptor::fullRho(rhoBar);
        T invRho = Descriptor::invRho(rhoBar);
        T uSqrLB = VectorTemplateImpl<T, Descriptor::d>::normSqr(uLB);

        T c_u_1 = -uLB[0];
        T c_u_2 = -uLB[1];
        T c_u_3 = -uLB[2];

        T c_u_4 = -uLB[0] - uLB[1];
        T c_u_5 = -uLB[0] + uLB[1];
        T c_u_6 = -uLB[0] - uLB[2];

        T c_u_7 = -uLB[0] + uLB[2];  // {-1, 0, 1}
        T c_u_8 = -uLB[1] - uLB[2];  // { 0,-1,-1}
        T c_u_9 = -uLB[1] + uLB[2];  // { 0,-1, 1}

        T c_u_10 = uLB[0];  // { 1, 0, 0}
        T c_u_11 = uLB[1];  // { 0, 1, 0}
        T c_u_12 = uLB[2];  // { 0, 0, 1}

        T c_u_13 = uLB[0] + uLB[1];  // { 1, 1, 0}
        T c_u_14 = uLB[0] - uLB[1];  //  { 1,-1, 0}
        T c_u_15 = uLB[0] + uLB[2];  // { 1, 0, 1}

        T c_u_16 = uLB[0] - uLB[2];  //  { 1, 0,-1}
        T c_u_17 = uLB[1] + uLB[2];  //  { 0, 1, 1}
        T c_u_18 = uLB[1] - uLB[2];  // { 0, 1,-1}

        // common terms fo ug_i calculations
        T uLB0_force0 = uLB[0] * force[0];
        T uLB1_force1 = uLB[1] * force[1];
        T uLB2_force2 = uLB[2] * force[2];

        T mines_force_0 = (T(-1) - uLB[0]) * force[0];
        T mines_force_1 = (T(-1) - uLB[1]) * force[1];
        T mines_force_2 = (T(-1) - uLB[2]) * force[2];

        T plus_force_0 = (T(1) - uLB[0]) * force[0];
        T plus_force_1 = (T(1) - uLB[1]) * force[1];
        T plus_force_2 = (T(1) - uLB[2]) * force[2];

        T ug_0 = -uLB0_force0 - uLB1_force1 - uLB2_force2;    //{0, 0, 0}
        T ug_1 = mines_force_0 - uLB1_force1 - uLB2_force2;   //{-1, 0, 0}
        T ug_2 = -uLB0_force0 + mines_force_1 - uLB2_force2;  //{ 0,-1, 0}
        T ug_3 = -uLB0_force0 - uLB1_force1 + mines_force_2;  //{ 0, 0,-1}

        T ug_4 = mines_force_0 + mines_force_1 - uLB2_force2;  //{-1,-1, 0}
        T ug_5 = mines_force_0 + plus_force_1 - uLB2_force2;   //{-1, 1, 0}
        T ug_6 = mines_force_0 - uLB1_force1 + mines_force_2;  //{-1, 0,-1}

        T ug_7 = mines_force_0 - uLB1_force1 + plus_force_2;    //{-1, 0, 1}
        T ug_8 = -uLB0_force0 + mines_force_1 + mines_force_2;  // { 0,-1,-1}
        T ug_9 = -uLB0_force0 + mines_force_1 + plus_force_2;   //{ 0,-1, 1}

        T ug_10 = plus_force_0 - uLB1_force1 - uLB2_force2;   //{1, 0, 0}
        T ug_11 = -uLB0_force0 + plus_force_1 - uLB2_force2;  //{ 0,1, 0}
        T ug_12 = -uLB0_force0 - uLB1_force1 + plus_force_2;  //{ 0, 0,1}

        T ug_13 = plus_force_0 + plus_force_1 - uLB2_force2;   //   { 1, 1, 0}
        T ug_14 = plus_force_0 + mines_force_1 - uLB2_force2;  //  { 1,-1, 0}
        T ug_15 = plus_force_0 - uLB1_force1 + plus_force_2;   // { 1, 0, 1}

        T ug_16 = plus_force_0 - uLB1_force1 + mines_force_2;   //   { 1, 0,-1}
        T ug_17 = -uLB0_force0 + plus_force_1 + plus_force_2;   //  { 0, 1, 1}
        T ug_18 = -uLB0_force0 + plus_force_1 + mines_force_2;  // { 0, 1,-1}

        // for equilibrium
        T invCs2_term = rhoFull * Descriptor::invCs2;
        T invCs2_sqr_term = rhoFull * (0.5 * Descriptor::invCs2 * Descriptor::invCs2);
        T common_eq_term = rhoFull * (1 - 0.5 * Descriptor::invCs2 * uSqrLB);

        T eqContribution_0 = Descriptor::t[0] * common_eq_term;
        T eqContribution_1 =
            Descriptor::t[1]
            * (common_eq_term + invCs2_term * c_u_1 + invCs2_sqr_term * c_u_1 * c_u_1);
        T eqContribution_2 =
            Descriptor::t[2]
            * (common_eq_term + invCs2_term * c_u_2 + invCs2_sqr_term * c_u_2 * c_u_2);
        T eqContribution_3 =
            Descriptor::t[3]
            * (common_eq_term + invCs2_term * c_u_3 + invCs2_sqr_term * c_u_3 * c_u_3);
        T eqContribution_4 =
            Descriptor::t[4]
            * (common_eq_term + invCs2_term * c_u_4 + invCs2_sqr_term * c_u_4 * c_u_4);
        T eqContribution_5 =
            Descriptor::t[5]
            * (common_eq_term + invCs2_term * c_u_5 + invCs2_sqr_term * c_u_5 * c_u_5);
        T eqContribution_6 =
            Descriptor::t[6]
            * (common_eq_term + invCs2_term * c_u_6 + invCs2_sqr_term * c_u_6 * c_u_6);
        T eqContribution_7 =
            Descriptor::t[7]
            * (common_eq_term + invCs2_term * c_u_7 + invCs2_sqr_term * c_u_7 * c_u_7);
        T eqContribution_8 =
            Descriptor::t[8]
            * (common_eq_term + invCs2_term * c_u_8 + invCs2_sqr_term * c_u_8 * c_u_8);
        T eqContribution_9 =
            Descriptor::t[9]
            * (common_eq_term + invCs2_term * c_u_9 + invCs2_sqr_term * c_u_9 * c_u_9);
        T eqContribution_10 =
            Descriptor::t[10]
            * (common_eq_term + invCs2_term * c_u_10 + invCs2_sqr_term * c_u_10 * c_u_10);
        T eqContribution_11 =
            Descriptor::t[11]
            * (common_eq_term + invCs2_term * c_u_11 + invCs2_sqr_term * c_u_11 * c_u_11);
        T eqContribution_12 =
            Descriptor::t[12]
            * (common_eq_term + invCs2_term * c_u_12 + invCs2_sqr_term * c_u_12 * c_u_12);
        T eqContribution_13 =
            Descriptor::t[13]
            * (common_eq_term + invCs2_term * c_u_13 + invCs2_sqr_term * c_u_13 * c_u_13);
        T eqContribution_14 =
            Descriptor::t[14]
            * (common_eq_term + invCs2_term * c_u_14 + invCs2_sqr_term * c_u_14 * c_u_14);
        T eqContribution_15 =
            Descriptor::t[15]
            * (common_eq_term + invCs2_term * c_u_15 + invCs2_sqr_term * c_u_15 * c_u_15);
        T eqContribution_16 =
            Descriptor::t[16]
            * (common_eq_term + invCs2_term * c_u_16 + invCs2_sqr_term * c_u_16 * c_u_16);
        T eqContribution_17 =
            Descriptor::t[17]
            * (common_eq_term + invCs2_term * c_u_17 + invCs2_sqr_term * c_u_17 * c_u_17);
        T eqContribution_18 =
            Descriptor::t[18]
            * (common_eq_term + invCs2_term * c_u_18 + invCs2_sqr_term * c_u_18 * c_u_18);

        Array<T, 19> f_full, forceMoments, forcing;

        forcing[0] = invRho * Descriptor::invCs2 * ug_0 * eqContribution_0;
        forcing[1] = invRho * Descriptor::invCs2 * ug_1 * eqContribution_1;
        forcing[2] = invRho * Descriptor::invCs2 * ug_2 * eqContribution_2;
        forcing[3] = invRho * Descriptor::invCs2 * ug_3 * eqContribution_3;
        forcing[4] = invRho * Descriptor::invCs2 * ug_4 * eqContribution_4;
        forcing[5] = invRho * Descriptor::invCs2 * ug_5 * eqContribution_5;
        forcing[6] = invRho * Descriptor::invCs2 * ug_6 * eqContribution_6;
        forcing[7] = invRho * Descriptor::invCs2 * ug_7 * eqContribution_7;
        forcing[8] = invRho * Descriptor::invCs2 * ug_8 * eqContribution_8;
        forcing[9] = invRho * Descriptor::invCs2 * ug_9 * eqContribution_9;
        forcing[10] = invRho * Descriptor::invCs2 * ug_10 * eqContribution_10;
        forcing[11] = invRho * Descriptor::invCs2 * ug_11 * eqContribution_11;
        forcing[12] = invRho * Descriptor::invCs2 * ug_12 * eqContribution_12;
        forcing[13] = invRho * Descriptor::invCs2 * ug_13 * eqContribution_13;
        forcing[14] = invRho * Descriptor::invCs2 * ug_14 * eqContribution_14;
        forcing[15] = invRho * Descriptor::invCs2 * ug_15 * eqContribution_15;
        forcing[16] = invRho * Descriptor::invCs2 * ug_16 * eqContribution_16;
        forcing[17] = invRho * Descriptor::invCs2 * ug_17 * eqContribution_17;
        forcing[18] = invRho * Descriptor::invCs2 * ug_18 * eqContribution_18;

        computeMoments(forceMoments, forcing);

        forceMoments[0] *= 0.5;
        forceMoments[1] *= 0.5;
        forceMoments[2] *= 0.5;
        forceMoments[3] *= 0.5;
        forceMoments[4] *= 0.5;
        forceMoments[5] *= 0.5;
        forceMoments[6] *= 0.5;
        forceMoments[7] *= 0.5;
        forceMoments[8] *= 0.5;
        forceMoments[9] *= 0.5;
        forceMoments[10] *= 0.5;
        forceMoments[11] *= 0.5;
        forceMoments[12] *= 0.5;
        forceMoments[13] *= 0.5;
        forceMoments[14] *= 0.5;
        forceMoments[15] *= 0.5;
        forceMoments[16] *= 0.5;
        forceMoments[17] *= 0.5;
        forceMoments[18] *= 0.5;

        // move from f_bar to full_f
        // 1. let's get full_f first
        f_full[0] = f[0] + Descriptor::SkordosFactor() * Descriptor::t[0];
        f_full[1] = f[1] + Descriptor::SkordosFactor() * Descriptor::t[1];
        f_full[2] = f[2] + Descriptor::SkordosFactor() * Descriptor::t[2];
        f_full[3] = f[3] + Descriptor::SkordosFactor() * Descriptor::t[3];
        f_full[4] = f[4] + Descriptor::SkordosFactor() * Descriptor::t[4];
        f_full[5] = f[5] + Descriptor::SkordosFactor() * Descriptor::t[5];
        f_full[6] = f[6] + Descriptor::SkordosFactor() * Descriptor::t[6];
        f_full[7] = f[7] + Descriptor::SkordosFactor() * Descriptor::t[7];
        f_full[8] = f[8] + Descriptor::SkordosFactor() * Descriptor::t[8];
        f_full[9] = f[9] + Descriptor::SkordosFactor() * Descriptor::t[9];
        f_full[10] = f[10] + Descriptor::SkordosFactor() * Descriptor::t[10];
        f_full[11] = f[11] + Descriptor::SkordosFactor() * Descriptor::t[11];
        f_full[12] = f[12] + Descriptor::SkordosFactor() * Descriptor::t[12];
        f_full[13] = f[13] + Descriptor::SkordosFactor() * Descriptor::t[13];
        f_full[14] = f[14] + Descriptor::SkordosFactor() * Descriptor::t[14];
        f_full[15] = f[15] + Descriptor::SkordosFactor() * Descriptor::t[15];
        f_full[16] = f[16] + Descriptor::SkordosFactor() * Descriptor::t[16];
        f_full[17] = f[17] + Descriptor::SkordosFactor() * Descriptor::t[17];
        f_full[18] = f[18] + Descriptor::SkordosFactor() * Descriptor::t[18];

        // I think that the next step is the tricky one
        computef_InvM_Smoments(f_full, forceMoments, omega);

        f_full[0] += forcing[0];
        f_full[1] += forcing[1];
        f_full[2] += forcing[2];
        f_full[3] += forcing[3];
        f_full[4] += forcing[4];
        f_full[5] += forcing[5];
        f_full[6] += forcing[6];
        f_full[7] += forcing[7];
        f_full[8] += forcing[8];
        f_full[9] += forcing[9];
        f_full[10] += forcing[10];
        f_full[11] += forcing[11];
        f_full[12] += forcing[12];
        f_full[13] += forcing[13];
        f_full[14] += forcing[14];
        f_full[15] += forcing[15];
        f_full[16] += forcing[16];
        f_full[17] += forcing[17];
        f_full[18] += forcing[18];

        f[0] = f_full[0] - Descriptor::SkordosFactor() * Descriptor::t[0];
        f[1] = f_full[1] - Descriptor::SkordosFactor() * Descriptor::t[1];
        f[2] = f_full[2] - Descriptor::SkordosFactor() * Descriptor::t[2];
        f[3] = f_full[3] - Descriptor::SkordosFactor() * Descriptor::t[3];
        f[4] = f_full[4] - Descriptor::SkordosFactor() * Descriptor::t[4];
        f[5] = f_full[5] - Descriptor::SkordosFactor() * Descriptor::t[5];
        f[6] = f_full[6] - Descriptor::SkordosFactor() * Descriptor::t[6];
        f[7] = f_full[7] - Descriptor::SkordosFactor() * Descriptor::t[7];
        f[8] = f_full[8] - Descriptor::SkordosFactor() * Descriptor::t[8];
        f[9] = f_full[9] - Descriptor::SkordosFactor() * Descriptor::t[9];
        f[10] = f_full[10] - Descriptor::SkordosFactor() * Descriptor::t[10];
        f[11] = f_full[11] - Descriptor::SkordosFactor() * Descriptor::t[11];
        f[12] = f_full[12] - Descriptor::SkordosFactor() * Descriptor::t[12];
        f[13] = f_full[13] - Descriptor::SkordosFactor() * Descriptor::t[13];
        f[14] = f_full[14] - Descriptor::SkordosFactor() * Descriptor::t[14];
        f[15] = f_full[15] - Descriptor::SkordosFactor() * Descriptor::t[15];
        f[16] = f_full[16] - Descriptor::SkordosFactor() * Descriptor::t[16];
        f[17] = f_full[17] - Descriptor::SkordosFactor() * Descriptor::t[17];
        f[18] = f_full[18] - Descriptor::SkordosFactor() * Descriptor::t[18];
    }

    /// MRT collision step
    static T mrtCollisionWithHeForce(
        Array<T, Descriptor::q> &f, const T &rhoBar, const Array<T, Descriptor::d> &u,
        const T &omega, const Array<T, Descriptor::d> &force, T amplitude)
    {
        Array<T, Descriptor::d> j = Descriptor::fullRho(rhoBar) * u;
        T jSqr = mrtCollision(f, rhoBar, j, omega);
        addHeForce(f, force, rhoBar, u, omega, amplitude);

        return jSqr;
    }

    /// Smagorinsky MRT collision step
    static T smagorinskyMrtCollisionWithForce(
        Array<T, Descriptor::q> &f, const T &rhoBar, const Array<T, Descriptor::d> &u,
        const Array<T, SymmetricTensorImpl<T, Descriptor::d>::n> &strain, T omega, T cSmago,
        const Array<T, Descriptor::d> &force, T amplitude)
    {
        Array<T, Descriptor::d> j = Descriptor::fullRho(rhoBar) * u;
        T jSqr = smagorinskyMrtCollision(f, rhoBar, j, strain, omega, cSmago);
        addGuoForce(f, force, u, omega, amplitude);

        return jSqr;
    }

    /// MRT collision step
    static T incMrtCollisionWithForce(
        Array<T, Descriptor::q> &f, const T &rhoBar, const Array<T, Descriptor::d> &u, T omega,
        const Array<T, Descriptor::d> &force, T amplitude)
    {
        Array<T, Descriptor::d> j = Descriptor::fullRho(rhoBar) * u;
        T jSqr = incMrtCollision(f, rhoBar, j, omega);
        addGuoForce(f, force, u, omega, amplitude);

        return jSqr;
    }

    /// Smagorinsky MRT collision step
    static T incSmagorinskyMrtCollisionWithForce(
        Array<T, Descriptor::q> &f, const T &rhoBar, const Array<T, Descriptor::d> &u,
        const Array<T, SymmetricTensorImpl<T, Descriptor::d>::n> &strain, T omega, T cSmago,
        const Array<T, Descriptor::d> &force, T amplitude)
    {
        T jSqr = incSmagorinskyMrtCollision(f, omega, strain, cSmago);
        addGuoForce(f, force, u, omega, amplitude);

        return jSqr;
    }

    /// Computation of all equilibrium distribution (in moments space)
    static void computeIncEquilibriumMoments(
        Array<T, Descriptor::q> &momentsEq, T rhoBar, Array<T, 3> const &j, T jSqr)
    {
        momentsEq[0] = rhoBar;
        momentsEq[1] = (T)19 * jSqr - (T)11 * rhoBar;
        momentsEq[2] = -(T)5.5 * jSqr + (T)3 * rhoBar;
        momentsEq[3] = j[0];
        momentsEq[4] = -((T)2 / (T)3) * j[0];
        momentsEq[5] = j[1];
        momentsEq[6] = -((T)2 / 3) * j[1];
        momentsEq[7] = j[2];
        momentsEq[8] = -((T)2 / (T)3) * j[2];
        momentsEq[9] = ((T)2 * j[0] * j[0] - j[1] * j[1] - j[2] * j[2]);
        momentsEq[10] = (-j[0] * j[0] + (T)0.5 * j[1] * j[1] + (T)0.5 * j[2] * j[2]);
        momentsEq[11] = (j[1] * j[1] - j[2] * j[2]);
        momentsEq[12] = (-(T)0.5 * j[1] * j[1] + (T)0.5 * j[2] * j[2]);
        momentsEq[13] = j[1] * j[0];
        momentsEq[14] = j[2] * j[1];
        momentsEq[15] = j[2] * j[0];
        momentsEq[16] = T();
        momentsEq[17] = T();
        momentsEq[18] = T();
    }

    /// Computation of all equilibrium distribution (in moments space)
    static void computeIncSmagorinskyEquilibrium(
        Array<T, Descriptor::q> &momentsEq, T rhoBar, Array<T, 3> const &j, T jSqr,
        const Array<T, 6> &strain, T cSmago)
    {
        typedef SymmetricTensorImpl<T, 3> S;
        T sNorm = std::sqrt((T)2 * SymmetricTensorImpl<T, 3>::tensorNormSqr(strain));
        T smagoFactor = (T)2 * cSmago * cSmago * sNorm;

        T ux2 = smagoFactor * strain[S::xx];
        T uy2 = smagoFactor * strain[S::yy];
        T uz2 = smagoFactor * strain[S::zz];

        T uxuy = smagoFactor * strain[S::xy];
        T uyuz = smagoFactor * strain[S::yz];
        T uxuz = smagoFactor * strain[S::xz];

        momentsEq[0] = rhoBar;
        momentsEq[1] = (T)19 * (jSqr + ux2 + uy2 + uz2) - (T)11 * rhoBar;
        momentsEq[2] = -(T)5.5 * (jSqr + ux2 + uy2 + uz2) + (T)3 * rhoBar;
        momentsEq[3] = j[0];
        momentsEq[4] = -((T)2 / (T)3) * j[0];
        momentsEq[5] = j[1];
        momentsEq[6] = -((T)2 / 3) * j[1];
        momentsEq[7] = j[2];
        momentsEq[8] = -((T)2 / (T)3) * j[2];
        momentsEq[9] = ((T)2 * (j[0] * j[0] + ux2) - (j[1] * j[1] + uy2) - (j[2] * j[2] + uz2));
        momentsEq[10] =
            (-(j[0] * j[0] + ux2) + (T)0.5 * (j[1] * j[1] + uy2) + (T)0.5 * (j[2] * j[2] + uz2));
        momentsEq[11] = (j[1] * j[1] + uy2 - (j[2] * j[2] + uz2));
        momentsEq[12] = (-(T)0.5 * (j[1] * j[1] + uy2) + (T)0.5 * (j[2] * j[2] + uz2));
        momentsEq[13] = j[1] * j[0] + uxuy;
        momentsEq[14] = j[2] * j[1] + uyuz;
        momentsEq[15] = j[2] * j[0] + uxuz;
        momentsEq[16] = T();
        momentsEq[17] = T();
        momentsEq[18] = T();
    }

    /// MRT collision step
    static T incMrtCollision(Array<T, Descriptor::q> &f, const T &omega)
    {
        Array<T, 19> moments, momentsEq;

        computeMoments(moments, f);
        T rhoBar = moments[0];
        Array<T, 3> j(
            moments[MRTDescriptor::momentumIndexes[0]], moments[MRTDescriptor::momentumIndexes[1]],
            moments[MRTDescriptor::momentumIndexes[2]]);
        T jSqr = VectorTemplateImpl<T, 3>::normSqr(j);
        computeIncEquilibriumMoments(momentsEq, rhoBar, j, jSqr);
        computeMneqInPlace(moments, momentsEq);  // moments become mNeq
        computef_InvM_Smoments(f, moments, omega);

        return jSqr;
    }

    /// MRT collision step
    static T incMrtCollision(
        Array<T, Descriptor::q> &f, const T &rhoBar, const Array<T, 3> &j, const T &omega)
    {
        Array<T, 19> moments, momentsEq;

        computeMoments(moments, f);
        T jSqr = VectorTemplateImpl<T, 3>::normSqr(j);
        computeIncEquilibriumMoments(momentsEq, rhoBar, j, jSqr);
        computeMneqInPlace(moments, momentsEq);  // moments become mNeq
        computef_InvM_Smoments(f, moments, omega);

        return jSqr;
    }

    /// MRT collision step
    static T incSmagorinskyMrtCollision(
        Array<T, Descriptor::q> &f, const T &omega, const Array<T, 6> &strain, T cSmago)
    {
        Array<T, 19> moments, momentsEq;

        computeMoments(moments, f);
        T rhoBar = moments[0];
        Array<T, 3> j(
            moments[MRTDescriptor::momentumIndexes[0]], moments[MRTDescriptor::momentumIndexes[1]],
            moments[MRTDescriptor::momentumIndexes[2]]);
        T jSqr = VectorTemplateImpl<T, 3>::normSqr(j);

        computeIncSmagorinskyEquilibrium(momentsEq, rhoBar, j, jSqr, strain, cSmago);
        computeMneqInPlace(moments, momentsEq);  // moments become mNeq
        computef_InvM_Smoments(f, moments, omega);

        return jSqr;
    }

    /// MRT collision step
    static T incSmagorinskyMrtCollision(
        Array<T, Descriptor::q> &f, const T &rhoBar, const Array<T, 3> &j, const T &omega,
        const Array<T, 6> &strain, T cSmago)
    {
        Array<T, 19> moments, momentsEq;

        computeMoments(moments, f);
        T jSqr = VectorTemplateImpl<T, 3>::normSqr(j);

        computeIncSmagorinskyEquilibrium(momentsEq, rhoBar, j, jSqr, strain, cSmago);
        computeMneqInPlace(moments, momentsEq);  // moments become mNeq
        computef_InvM_Smoments(f, moments, omega);

        return jSqr;
    }
};

}  // namespace plb

#endif
