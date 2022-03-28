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
 * 3D specialization of dynamicsTemplates functions.
 * Theoretical background about these collision models can be found in
 * Coreixas et al. 'Comprehensive comparison of collision models in the
 * lattice Boltzmann framework: Theoretical investigations', PRE, 2019.
 */

#ifndef COMPREHENSIVE_MODELS_TEMPLATES_3D_H
#define COMPREHENSIVE_MODELS_TEMPLATES_3D_H

#include "core/globalDefs.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"

namespace plb {
// Efficient specialization for D3Q27 lattice
template <typename T>
struct comprehensiveDynamicsTemplatesImpl<T, descriptors::D3Q27DescriptorBase<T> > {
    typedef descriptors::D3Q27DescriptorBase<T> D;
    // Same order as in E1.
    enum {
        // Order 0
        M000 = 0,

        // Order 1
        M100 = 1,
        M010 = 2,
        M001 = 3,

        // Order 2
        M200 = 4,
        M020 = 5,
        M002 = 6,
        M110 = 7,
        M101 = 8,
        M011 = 9,

        // Order 3
        M210 = 10,
        M201 = 11,
        M021 = 12,
        M120 = 13,
        M102 = 14,
        M012 = 15,
        M111 = 16,

        // Order 4
        M220 = 17,
        M202 = 18,
        M022 = 19,
        M211 = 20,
        M121 = 21,
        M112 = 22,

        // Order 5
        M221 = 23,
        M212 = 24,
        M122 = 25,

        // Order 6
        M222 = 26
    };

    /**
     * // General way to compute RMs
     * static void RMcomputeMoments(Array<T,D::q> const& cell, Array<T, D::q>& RM, T& rho) {
     *
     *     Array<T, D::q> f;
     *     for (int i = 0; i<27; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *         RM[i] = 0.;
     *     }
     *
     *     for (int i = 0; i<27; ++i) {
     *         // Order 0
     *         RM[M000] += f[i];
     *         // Order 1
     *         RM[M100] += D::c[i][0] * f[i];
     *         RM[M010] += D::c[i][1] * f[i];
     *         RM[M001] += D::c[i][2] * f[i];
     *         // Order 2
     *         RM[M200] += D::c[i][0] * D::c[i][0] * f[i];
     *         RM[M020] += D::c[i][1] * D::c[i][1] * f[i];
     *         RM[M002] += D::c[i][2] * D::c[i][2] * f[i];
     *         RM[M110] += D::c[i][0] * D::c[i][1] * f[i];
     *         RM[M101] += D::c[i][0] * D::c[i][2] * f[i];
     *         RM[M011] += D::c[i][1] * D::c[i][2] * f[i];
     *         // Order 3
     *         RM[M210] += D::c[i][0] * D::c[i][0] * D::c[i][1] * f[i];
     *         RM[M201] += D::c[i][0] * D::c[i][0] * D::c[i][2] * f[i];
     *         RM[M021] += D::c[i][1] * D::c[i][1] * D::c[i][2] * f[i];
     *         RM[M120] += D::c[i][0] * D::c[i][1] * D::c[i][1] * f[i];
     *         RM[M102] += D::c[i][0] * D::c[i][2] * D::c[i][2] * f[i];
     *         RM[M012] += D::c[i][1] * D::c[i][2] * D::c[i][2] * f[i];
     *         RM[M111] += D::c[i][0] * D::c[i][1] * D::c[i][2] * f[i];
     *         // Order 4
     *         RM[M220] += D::c[i][0] * D::c[i][0] * D::c[i][1] * D::c[i][1] * f[i];
     *         RM[M202] += D::c[i][0] * D::c[i][0] * D::c[i][2] * D::c[i][2] * f[i];
     *         RM[M022] += D::c[i][1] * D::c[i][1] * D::c[i][2] * D::c[i][2] * f[i];
     *         RM[M211] += D::c[i][0] * D::c[i][0] * D::c[i][1] * D::c[i][2] * f[i];
     *         RM[M121] += D::c[i][0] * D::c[i][1] * D::c[i][1] * D::c[i][2] * f[i];
     *         RM[M112] += D::c[i][0] * D::c[i][1] * D::c[i][2] * D::c[i][2] * f[i];
     *         // Order 5
     *         RM[M221] += D::c[i][0] * D::c[i][0] * D::c[i][1] * D::c[i][1] * D::c[i][2] * f[i];
     *         RM[M212] += D::c[i][0] * D::c[i][0] * D::c[i][1] * D::c[i][2] * D::c[i][2] * f[i];
     *         RM[M122] += D::c[i][0] * D::c[i][1] * D::c[i][1] * D::c[i][2] * D::c[i][2] * f[i];
     *         // Order 6
     *         RM[M222] += D::c[i][0] * D::c[i][0] * D::c[i][1] * D::c[i][1] * D::c[i][2] *
     * D::c[i][2] * f[i];
     *     }
     *
     *     rho = RM[M000];
     *     T invRho = 1. / RM[M000];
     *     for (int i = 0; i<27; ++i) {
     *         RM[i] *= invRho;
     *     }
     * };
     *
     * // Optimized way to compute RMs based on Palabos ordering of discrete velocities
     * static void RMcomputeMoments(Array<T,D::q> const& cell, Array<T, D::q>& RM, T& rho) {
     *
     *     Array<T, D::q> f;
     *     for (int i = 0; i<27; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *     }
     *
     *     // Order 0
     *     rho =   f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[10] +
     * f[11] + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18] + f[19] + f[20] + f[21] + f[22]
     * + f[23] + f[24] + f[25] + f[26]; RM[M000] = 1.0; T invRho = 1./rho;
     *     // Order 1
     *     RM[M100] = invRho * (- f[1] - f[4] - f[5] - f[6] - f[7] - f[10] - f[11] - f[12] - f[13] +
     * f[14] + f[17] + f[18] + f[19] + f[20] + f[23] + f[24] + f[25] + f[26]); RM[M010] = invRho *
     * (- f[2] - f[4] + f[5] - f[8] - f[9] - f[10] - f[11] + f[12] + f[13] + f[15] + f[17] - f[18] +
     * f[21] + f[22] + f[23] + f[24] - f[25] - f[26]); RM[M001] = invRho * (- f[3] - f[6] + f[7] -
     * f[8] + f[9] - f[10] + f[11] - f[12] + f[13] + f[16] + f[19] - f[20] + f[21] - f[22] + f[23] -
     * f[24] + f[25] - f[26]);
     *     // Order 2
     *     RM[M200] = invRho * (  f[1] + f[4] + f[5] + f[6] + f[7] + f[10] + f[11] + f[12] + f[13] +
     * f[14] + f[17] + f[18] + f[19] + f[20] + f[23] + f[24] + f[25] + f[26]); RM[M020] = invRho * (
     * f[2] + f[4] + f[5] + f[8] + f[9] + f[10] + f[11] + f[12] + f[13] + f[15] + f[17] + f[18] +
     * f[21] + f[22] + f[23] + f[24] + f[25] + f[26]); RM[M002] = invRho * (  f[3] + f[6] + f[7] +
     * f[8] + f[9] + f[10] + f[11] + f[12] + f[13] + f[16] + f[19] + f[20] + f[21] + f[22] + f[23] +
     * f[24] + f[25] + f[26]); RM[M110] = invRho * (  f[4] - f[5] + f[10] + f[11] - f[12] - f[13] +
     * f[17] - f[18] + f[23] + f[24] - f[25] - f[26]); RM[M101] = invRho * (  f[6] - f[7] + f[10] -
     * f[11] + f[12] - f[13] + f[19] - f[20] + f[23] - f[24] + f[25] - f[26]); RM[M011] = invRho * (
     * f[8] - f[9] + f[10] - f[11] - f[12] + f[13] + f[21] - f[22] + f[23] - f[24] - f[25] + f[26]);
     *     // Order 3
     *     RM[M210] = invRho * (- f[4] + f[5] - f[10] - f[11] + f[12] + f[13] + f[17] - f[18] +
     * f[23] + f[24] - f[25] - f[26]); RM[M201] = invRho * (- f[6] + f[7] - f[10] + f[11] - f[12] +
     * f[13] + f[19] - f[20] + f[23] - f[24] + f[25] - f[26]); RM[M021] = invRho * (- f[8] + f[9] -
     * f[10] + f[11] - f[12] + f[13] + f[21] - f[22] + f[23] - f[24] + f[25] - f[26]); RM[M120] =
     * invRho * (- f[4] - f[5] - f[10] - f[11] - f[12] - f[13] + f[17] + f[18] + f[23] + f[24] +
     * f[25] + f[26]); RM[M102] = invRho * (- f[6] - f[7] - f[10] - f[11] - f[12] - f[13] + f[19] +
     * f[20] + f[23] + f[24] + f[25] + f[26]); RM[M012] = invRho * (- f[8] - f[9] - f[10] - f[11] +
     * f[12] + f[13] + f[21] + f[22] + f[23] + f[24] - f[25] - f[26]); RM[M111] = invRho * (- f[10]
     * + f[11] + f[12] - f[13] + f[23] - f[24] - f[25] + f[26]);
     *     // Order 4
     *     RM[M220] = invRho * (  f[4] + f[5] + f[10] + f[11] + f[12] + f[13] + f[17] + f[18] +
     * f[23] + f[24] + f[25] + f[26]); RM[M202] = invRho * (  f[6] + f[7] + f[10] + f[11] + f[12] +
     * f[13] + f[19] + f[20] + f[23] + f[24] + f[25] + f[26]); RM[M022] = invRho * (  f[8] + f[9] +
     * f[10] + f[11] + f[12] + f[13] + f[21] + f[22] + f[23] + f[24] + f[25] + f[26]); RM[M211] =
     * invRho * (  f[10] - f[11] - f[12] + f[13] + f[23] - f[24] - f[25] + f[26]); RM[M121] = invRho
     * * (  f[10] - f[11] + f[12] - f[13] + f[23] - f[24] + f[25] - f[26]); RM[M112] = invRho * (
     * f[10] + f[11] - f[12] - f[13] + f[23] + f[24] - f[25] - f[26]);
     *     // Order 5
     *     RM[M221] = invRho * (- f[10] + f[11] - f[12] + f[13] + f[23] - f[24] + f[25] - f[26]);
     *     RM[M212] = invRho * (- f[10] - f[11] + f[12] + f[13] + f[23] + f[24] - f[25] - f[26]);
     *     RM[M122] = invRho * (- f[10] - f[11] - f[12] - f[13] + f[23] + f[24] + f[25] + f[26]);
     *     // Order 6
     *     RM[M222] = invRho * (  f[10] + f[11] + f[12] + f[13] + f[23] + f[24] + f[25] + f[26]);
     * };
     */

    // Optimized way to compute RMs based on the general ordering of discrete velocities
    static void RMcomputeMoments(Array<T, D::q> const &cell, Array<T, D::q> &RM, T &rho)
    {
        Array<T, D::q> f;
        for (int i = 0; i < 27; ++i) {
            f[i] = cell[i] + D::SkordosFactor() * D::t[i];
        }

        T A1 = f[FMMM] + f[FMMP] + f[FMPM] + f[FMPP];
        T A2 = f[FPPP] + f[FPPM] + f[FPMP] + f[FPMM];
        T A3 = f[FMMM] + f[FMMP] - f[FMPM] - f[FMPP];
        T A4 = f[FPPP] + f[FPPM] - f[FPMP] - f[FPMM];
        T A5 = f[FMMM] - f[FMMP] + f[FMPM] - f[FMPP];
        T A6 = f[FPPP] - f[FPPM] + f[FPMP] - f[FPMM];
        T A7 = f[FMMM] - f[FMMP] - f[FMPM] + f[FMPP];
        T A8 = f[FPPP] - f[FPPM] - f[FPMP] + f[FPMM];

        T X_M1 = f[FM00] + f[FMM0] + f[FMP0] + f[FM0M] + f[FM0P] + A1;
        T X_P1 = f[FP00] + f[FPP0] + f[FPM0] + f[FP0P] + f[FP0M] + A2;
        T X_0 =
            f[F0M0] + f[F00M] + f[F0MM] + f[F0MP] + f[F000] + f[F0P0] + f[F00P] + f[F0PP] + f[F0PM];

        T Y_M1 =
            f[F0M0] + f[FMM0] + f[F0MM] + f[F0MP] + f[FMMM] + f[FMMP] + f[FPM0] + f[FPMP] + f[FPMM];
        T Y_P1 =
            f[F0P0] + f[FPP0] + f[F0PP] + f[F0PM] + f[FPPP] + f[FPPM] + f[FMP0] + f[FMPM] + f[FMPP];

        T Z_M1 =
            f[F00M] + f[FM0M] + f[F0MM] + f[FMMM] + f[FMPM] + f[FP0M] + f[F0PM] + f[FPPM] + f[FPMM];
        T Z_P1 =
            f[F00P] + f[FP0P] + f[F0PP] + f[FPPP] + f[FPMP] + f[FM0P] + f[F0MP] + f[FMMP] + f[FMPP];

        rho = X_M1 + X_P1 + X_0;
        T invRho = 1. / rho;

        // Order 6
        RM[M222] = invRho * (A1 + A2);
        // Order 5
        RM[M221] = invRho * (-A5 + A6);
        RM[M212] = invRho * (-A3 + A4);
        RM[M122] = invRho * (-A1 + A2);
        // Order 4
        RM[M220] = invRho * (f[FMM0] + f[FMP0] + f[FPP0] + f[FPM0]) + RM[M222];
        RM[M202] = invRho * (f[FM0M] + f[FM0P] + f[FP0P] + f[FP0M]) + RM[M222];
        RM[M022] = invRho * (f[F0MM] + f[F0MP] + f[F0PP] + f[F0PM]) + RM[M222];
        RM[M211] = invRho * (A7 + A8);
        RM[M121] = invRho * (A5 + A6);
        RM[M112] = invRho * (A3 + A4);
        // Order 3
        RM[M210] = invRho * (-f[FMM0] + f[FMP0] + f[FPP0] - f[FPM0]) + RM[M212];
        RM[M201] = invRho * (-f[FM0M] + f[FM0P] + f[FP0P] - f[FP0M]) + RM[M221];
        RM[M021] = invRho * (-f[F0MM] + f[F0MP] + f[F0PP] - f[F0PM]) + RM[M221];
        RM[M120] = invRho * (-f[FMM0] - f[FMP0] + f[FPP0] + f[FPM0]) + RM[M122];
        RM[M102] = invRho * (-f[FM0M] - f[FM0P] + f[FP0P] + f[FP0M]) + RM[M122];
        RM[M012] = invRho * (-f[F0MM] - f[F0MP] + f[F0PP] + f[F0PM]) + RM[M212];
        RM[M111] = invRho * (-A7 + A8);
        // Order 2
        RM[M200] = invRho * (X_P1 + X_M1);
        RM[M020] = invRho * (Y_P1 + Y_M1);
        RM[M002] = invRho * (Z_P1 + Z_M1);
        RM[M110] = invRho * (f[FMM0] - f[FMP0] + f[FPP0] - f[FPM0]) + RM[M112];
        RM[M101] = invRho * (f[FM0M] - f[FM0P] + f[FP0P] - f[FP0M]) + RM[M121];
        RM[M011] = invRho * (f[F0MM] - f[F0MP] + f[F0PP] - f[F0PM]) + RM[M211];
        // Order 1
        RM[M100] = invRho * (X_P1 - X_M1);
        RM[M010] = invRho * (Y_P1 - Y_M1);
        RM[M001] = invRho * (Z_P1 - Z_M1);

        // Order 0
        RM[M000] = 1.;
    };

    static void RMcomputeEquilibriumMoments(Array<T, D::d> const &u, Array<T, D::q> &RMeq)
    {
        // Order 0
        RMeq[M000] = 1.;
        // Order 1
        RMeq[M100] = u[0];
        RMeq[M010] = u[1];
        RMeq[M001] = u[2];
        // Order 2
        RMeq[M200] = u[0] * u[0] + D::cs2;
        RMeq[M020] = u[1] * u[1] + D::cs2;
        RMeq[M002] = u[2] * u[2] + D::cs2;
        RMeq[M110] = u[0] * u[1];
        RMeq[M101] = u[0] * u[2];
        RMeq[M011] = u[1] * u[2];
        // Order 3
        RMeq[M210] = RMeq[M200] * u[1];
        RMeq[M201] = RMeq[M200] * u[2];
        RMeq[M021] = RMeq[M020] * u[2];
        RMeq[M120] = RMeq[M020] * u[0];
        RMeq[M102] = RMeq[M002] * u[0];
        RMeq[M012] = RMeq[M002] * u[1];
        RMeq[M111] = RMeq[M110] * u[2];
        // Order 4
        RMeq[M220] = RMeq[M200] * RMeq[M020];
        RMeq[M202] = RMeq[M200] * RMeq[M002];
        RMeq[M022] = RMeq[M020] * RMeq[M002];
        RMeq[M211] = RMeq[M200] * RMeq[M011];
        RMeq[M121] = RMeq[M020] * RMeq[M101];
        RMeq[M112] = RMeq[M002] * RMeq[M110];
        // Order 5
        RMeq[M221] = RMeq[M220] * RMeq[M001];
        RMeq[M212] = RMeq[M202] * RMeq[M010];
        RMeq[M122] = RMeq[M022] * RMeq[M100];
        // Order 6
        RMeq[M222] = RMeq[M220] * RMeq[M002];
    };

    enum {
        F000 = 0,
        FM00 = 1,
        F0M0 = 2,
        F00M = 3,
        FMM0 = 4,
        FMP0 = 5,
        FM0M = 6,
        FM0P = 7,
        F0MM = 8,
        F0MP = 9,
        FMMM = 10,
        FMMP = 11,
        FMPM = 12,
        FMPP = 13,

        FP00 = 14,
        F0P0 = 15,
        F00P = 16,
        FPP0 = 17,
        FPM0 = 18,
        FP0P = 19,
        FP0M = 20,
        F0PP = 21,
        F0PM = 22,
        FPPP = 23,
        FPPM = 24,
        FPMP = 25,
        FPMM = 26
    };

    // Equilibrium populations based on 27 moments can be computed using either RM, HM, CM, CHM or
    // Gauss-Hermite formalisms. Here we use raw moments (RMs)
    static void RMcomputeEquilibrium(T rho, Array<T, D::q> const &RMeq, Array<T, D::q> &eq)
    {
        Array<T, D::d> u(RMeq[1], RMeq[2], RMeq[3]);
        // Optimization based on symmetries between populations and their opposite counterpart
        eq[F000] = rho
                   * (1. - RMeq[M200] - RMeq[M020] - RMeq[M002] + RMeq[M220] + RMeq[M202]
                      + RMeq[M022] - RMeq[M222]);

        eq[FP00] = 0.5 * rho
                   * (u[0] + RMeq[M200] - RMeq[M120] - RMeq[M102] - RMeq[M220] - RMeq[M202]
                      + RMeq[M122] + RMeq[M222]);
        eq[FM00] = rho * (-u[0] + RMeq[M120] + RMeq[M102] - RMeq[M122]) + eq[FP00];

        eq[F0P0] = 0.5 * rho
                   * (u[1] + RMeq[M020] - RMeq[M210] - RMeq[M012] - RMeq[M220] - RMeq[M022]
                      + RMeq[M212] + RMeq[M222]);
        eq[F0M0] = rho * (-u[1] + RMeq[M210] + RMeq[M012] - RMeq[M212]) + eq[F0P0];

        eq[F00P] = 0.5 * rho
                   * (u[2] + RMeq[M002] - RMeq[M201] - RMeq[M021] - RMeq[M202] - RMeq[M022]
                      + RMeq[M221] + RMeq[M222]);
        eq[F00M] = rho * (-u[2] + RMeq[M201] + RMeq[M021] - RMeq[M221]) + eq[F00P];

        eq[FPP0] = 0.25 * rho
                   * (RMeq[M110] + RMeq[M210] + RMeq[M120] - RMeq[M112] + RMeq[M220] - RMeq[M212]
                      - RMeq[M122] - RMeq[M222]);
        eq[FMP0] = 0.5 * rho * (-RMeq[M110] - RMeq[M120] + RMeq[M112] + RMeq[M122]) + eq[FPP0];
        eq[FPM0] = 0.5 * rho * (-RMeq[M110] - RMeq[M210] + RMeq[M112] + RMeq[M212]) + eq[FPP0];
        eq[FMM0] = 0.5 * rho * (-RMeq[M210] - RMeq[M120] + RMeq[M212] + RMeq[M122]) + eq[FPP0];

        eq[FP0P] = 0.25 * rho
                   * (RMeq[M101] + RMeq[M201] + RMeq[M102] - RMeq[M121] + RMeq[M202] - RMeq[M221]
                      - RMeq[M122] - RMeq[M222]);
        eq[FM0P] = 0.5 * rho * (-RMeq[M101] - RMeq[M102] + RMeq[M121] + RMeq[M122]) + eq[FP0P];
        eq[FP0M] = 0.5 * rho * (-RMeq[M101] - RMeq[M201] + RMeq[M121] + RMeq[M221]) + eq[FP0P];
        eq[FM0M] = 0.5 * rho * (-RMeq[M201] - RMeq[M102] + RMeq[M221] + RMeq[M122]) + eq[FP0P];

        eq[F0PP] = 0.25 * rho
                   * (RMeq[M011] + RMeq[M021] + RMeq[M012] - RMeq[M211] + RMeq[M022] - RMeq[M221]
                      - RMeq[M212] - RMeq[M222]);
        eq[F0MP] = 0.5 * rho * (-RMeq[M011] - RMeq[M012] + RMeq[M211] + RMeq[M212]) + eq[F0PP];
        eq[F0PM] = 0.5 * rho * (-RMeq[M011] - RMeq[M021] + RMeq[M211] + RMeq[M221]) + eq[F0PP];
        eq[F0MM] = 0.5 * rho * (-RMeq[M021] - RMeq[M012] + RMeq[M221] + RMeq[M212]) + eq[F0PP];

        eq[FPPP] = 0.125 * rho
                   * (RMeq[M111] + RMeq[M211] + RMeq[M121] + RMeq[M112] + RMeq[M221] + RMeq[M212]
                      + RMeq[M122] + RMeq[M222]);
        eq[FMPP] = 0.25 * rho * (-RMeq[M111] - RMeq[M121] - RMeq[M112] - RMeq[M122]) + eq[FPPP];
        eq[FPMP] = 0.25 * rho * (-RMeq[M111] - RMeq[M211] - RMeq[M112] - RMeq[M212]) + eq[FPPP];
        eq[FPPM] = 0.25 * rho * (-RMeq[M111] - RMeq[M211] - RMeq[M121] - RMeq[M221]) + eq[FPPP];
        eq[FMMP] = 0.25 * rho * (-RMeq[M211] - RMeq[M121] - RMeq[M212] - RMeq[M122]) + eq[FPPP];
        eq[FMPM] = 0.25 * rho * (-RMeq[M211] - RMeq[M112] - RMeq[M221] - RMeq[M122]) + eq[FPPP];
        eq[FPMM] = 0.25 * rho * (-RMeq[M121] - RMeq[M112] - RMeq[M221] - RMeq[M212]) + eq[FPPP];
        eq[FMMM] = 0.25 * rho * (-RMeq[M111] - RMeq[M221] - RMeq[M212] - RMeq[M122]) + eq[FPPP];
    };

    static void RMcollide(
        Array<T, D::q> &cell, T rho, Array<T, D::d> const &u,
        Array<T, D::q> const &RM,    // Raw moments
        Array<T, D::q> const &RMeq,  // Equilibrium moments (raw)
        Array<T, D::numRelaxationTimes> const &omega)
    {
        T omega1 = omega[0];
        T omega2 = omega[1];
        T omega3 = omega[2];
        T omega4 = omega[3];
        T omega5 = omega[4];
        T omega6 = omega[5];
        T omega7 = omega[6];
        T omega8 = omega[7];

        T omegaBulk = omega[8];
        T omegaPlus = (omegaBulk + 2. * omega1) / 3.;  // Notation used by Fei
        T omegaMinus = (omegaBulk - omega1) / 3.;      // Notation used by Fei

        // Post-collision moments.
        Array<T, D::q> RMcoll;

        // Order 2 (non-diagonal collision so that we can easily modify the bulk viscosity)
        RMcoll[M200] = RM[M200] - omegaPlus * (RM[M200] - RMeq[M200])
                       - omegaMinus * (RM[M020] - RMeq[M020])
                       - omegaMinus * (RM[M002] - RMeq[M002]);
        RMcoll[M020] = RM[M020] - omegaMinus * (RM[M200] - RMeq[M200])
                       - omegaPlus * (RM[M020] - RMeq[M020]) - omegaMinus * (RM[M002] - RMeq[M002]);
        RMcoll[M002] = RM[M002] - omegaMinus * (RM[M200] - RMeq[M200])
                       - omegaMinus * (RM[M020] - RMeq[M020]) - omegaPlus * (RM[M002] - RMeq[M002]);

        RMcoll[M110] = (1. - omega2) * RM[M110] + omega2 * RMeq[M110];
        RMcoll[M101] = (1. - omega2) * RM[M101] + omega2 * RMeq[M101];
        RMcoll[M011] = (1. - omega2) * RM[M011] + omega2 * RMeq[M011];

        // Order 3
        RMcoll[M210] = (1. - omega3) * RM[M210] + omega3 * RMeq[M210];
        RMcoll[M201] = (1. - omega3) * RM[M201] + omega3 * RMeq[M201];
        RMcoll[M021] = (1. - omega3) * RM[M021] + omega3 * RMeq[M021];
        RMcoll[M120] = (1. - omega3) * RM[M120] + omega3 * RMeq[M120];
        RMcoll[M102] = (1. - omega3) * RM[M102] + omega3 * RMeq[M102];
        RMcoll[M012] = (1. - omega3) * RM[M012] + omega3 * RMeq[M012];

        RMcoll[M111] = (1. - omega4) * RM[M111] + omega4 * RMeq[M111];

        // Order 4
        RMcoll[M220] = (1. - omega5) * RM[M220] + omega5 * RMeq[M220];
        RMcoll[M202] = (1. - omega5) * RM[M202] + omega5 * RMeq[M202];
        RMcoll[M022] = (1. - omega5) * RM[M022] + omega5 * RMeq[M022];

        RMcoll[M211] = (1. - omega6) * RM[M211] + omega6 * RMeq[M211];
        RMcoll[M121] = (1. - omega6) * RM[M121] + omega6 * RMeq[M121];
        RMcoll[M112] = (1. - omega6) * RM[M112] + omega6 * RMeq[M112];

        // Order 5
        RMcoll[M221] = (1. - omega7) * RM[M221] + omega7 * RMeq[M221];
        RMcoll[M212] = (1. - omega7) * RM[M212] + omega7 * RMeq[M212];
        RMcoll[M122] = (1. - omega7) * RM[M122] + omega7 * RMeq[M122];

        // Order 6
        RMcoll[M222] = (1. - omega8) * RM[M222] + omega8 * RMeq[M222];

        // Optimization based on symmetries between populations and their opposite counterpart
        cell[F000] = rho
                     * (1. - RMcoll[M200] - RMcoll[M020] - RMcoll[M002] + RMcoll[M220]
                        + RMcoll[M202] + RMcoll[M022] - RMcoll[M222]);

        cell[FP00] = 0.5 * rho
                     * (u[0] + RMcoll[M200] - RMcoll[M120] - RMcoll[M102] - RMcoll[M220]
                        - RMcoll[M202] + RMcoll[M122] + RMcoll[M222]);
        cell[FM00] = rho * (-u[0] + RMcoll[M120] + RMcoll[M102] - RMcoll[M122]) + cell[FP00];

        cell[F0P0] = 0.5 * rho
                     * (u[1] + RMcoll[M020] - RMcoll[M210] - RMcoll[M012] - RMcoll[M220]
                        - RMcoll[M022] + RMcoll[M212] + RMcoll[M222]);
        cell[F0M0] = rho * (-u[1] + RMcoll[M210] + RMcoll[M012] - RMcoll[M212]) + cell[F0P0];

        cell[F00P] = 0.5 * rho
                     * (u[2] + RMcoll[M002] - RMcoll[M201] - RMcoll[M021] - RMcoll[M202]
                        - RMcoll[M022] + RMcoll[M221] + RMcoll[M222]);
        cell[F00M] = rho * (-u[2] + RMcoll[M201] + RMcoll[M021] - RMcoll[M221]) + cell[F00P];

        cell[FPP0] = 0.25 * rho
                     * (RMcoll[M110] + RMcoll[M210] + RMcoll[M120] - RMcoll[M112] + RMcoll[M220]
                        - RMcoll[M212] - RMcoll[M122] - RMcoll[M222]);
        cell[FMP0] =
            0.5 * rho * (-RMcoll[M110] - RMcoll[M120] + RMcoll[M112] + RMcoll[M122]) + cell[FPP0];
        cell[FPM0] =
            0.5 * rho * (-RMcoll[M110] - RMcoll[M210] + RMcoll[M112] + RMcoll[M212]) + cell[FPP0];
        cell[FMM0] =
            0.5 * rho * (-RMcoll[M210] - RMcoll[M120] + RMcoll[M212] + RMcoll[M122]) + cell[FPP0];

        cell[FP0P] = 0.25 * rho
                     * (RMcoll[M101] + RMcoll[M201] + RMcoll[M102] - RMcoll[M121] + RMcoll[M202]
                        - RMcoll[M221] - RMcoll[M122] - RMcoll[M222]);
        cell[FM0P] =
            0.5 * rho * (-RMcoll[M101] - RMcoll[M102] + RMcoll[M121] + RMcoll[M122]) + cell[FP0P];
        cell[FP0M] =
            0.5 * rho * (-RMcoll[M101] - RMcoll[M201] + RMcoll[M121] + RMcoll[M221]) + cell[FP0P];
        cell[FM0M] =
            0.5 * rho * (-RMcoll[M201] - RMcoll[M102] + RMcoll[M221] + RMcoll[M122]) + cell[FP0P];

        cell[F0PP] = 0.25 * rho
                     * (RMcoll[M011] + RMcoll[M021] + RMcoll[M012] - RMcoll[M211] + RMcoll[M022]
                        - RMcoll[M221] - RMcoll[M212] - RMcoll[M222]);
        cell[F0MP] =
            0.5 * rho * (-RMcoll[M011] - RMcoll[M012] + RMcoll[M211] + RMcoll[M212]) + cell[F0PP];
        cell[F0PM] =
            0.5 * rho * (-RMcoll[M011] - RMcoll[M021] + RMcoll[M211] + RMcoll[M221]) + cell[F0PP];
        cell[F0MM] =
            0.5 * rho * (-RMcoll[M021] - RMcoll[M012] + RMcoll[M221] + RMcoll[M212]) + cell[F0PP];

        cell[FPPP] = 0.125 * rho
                     * (RMcoll[M111] + RMcoll[M211] + RMcoll[M121] + RMcoll[M112] + RMcoll[M221]
                        + RMcoll[M212] + RMcoll[M122] + RMcoll[M222]);
        cell[FMPP] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M121] - RMcoll[M112] - RMcoll[M122]) + cell[FPPP];
        cell[FPMP] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M211] - RMcoll[M112] - RMcoll[M212]) + cell[FPPP];
        cell[FPPM] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M211] - RMcoll[M121] - RMcoll[M221]) + cell[FPPP];
        cell[FMMP] =
            0.25 * rho * (-RMcoll[M211] - RMcoll[M121] - RMcoll[M212] - RMcoll[M122]) + cell[FPPP];
        cell[FMPM] =
            0.25 * rho * (-RMcoll[M211] - RMcoll[M112] - RMcoll[M221] - RMcoll[M122]) + cell[FPPP];
        cell[FPMM] =
            0.25 * rho * (-RMcoll[M121] - RMcoll[M112] - RMcoll[M221] - RMcoll[M212]) + cell[FPPP];
        cell[FMMM] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M221] - RMcoll[M212] - RMcoll[M122]) + cell[FPPP];

        for (int i = 0; i < 27; ++i) {
            cell[i] -= D::SkordosFactor() * D::t[i];
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////
    // Hermite Moments Formalism (Equilibrium is computed through raw moments formalism) //
    ///////////////////////////////////////////////////////////////////////////////////////

    /**
     * // General way to compute HMs
     * static void HMcomputeMoments(Array<T,D::q> const& cell, Array<T, D::q>& HM) {
     *
     *     Array<T, D::q> f;
     *     for (int i = 0; i<27; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *         HM[i] = 0.;
     *     }
     *
     *     for (int i = 0; i<27; ++i) {
     *
     *         Hxx = D::c[i][0] * D::c[i][0] - D::cs2;
     *         Hyy = D::c[i][1] * D::c[i][1] - D::cs2;
     *         Hzz = D::c[i][2] * D::c[i][2] - D::cs2;
     *
     *         // Order 0
     *         HM[M000] += f[i];
     *         // Order 1
     *         HM[M100] += D::c[i][0] * f[i];
     *         HM[M010] += D::c[i][1] * f[i];
     *         HM[M001] += D::c[i][2] * f[i];
     *         // Order 2
     *         HM[M200] += Hxx * f[i];
     *         HM[M020] += Hyy * f[i];
     *         HM[M002] += Hzz * f[i];
     *         HM[M110] += D::c[i][0] * D::c[i][1] * f[i];
     *         HM[M101] += D::c[i][0] * D::c[i][2] * f[i];
     *         HM[M011] += D::c[i][1] * D::c[i][2] * f[i];
     *         // Order 3
     *         HM[M210] += Hxx * D::c[i][1] * f[i];
     *         HM[M201] += Hxx * D::c[i][2] * f[i];
     *         HM[M021] += Hyy * D::c[i][2] * f[i];
     *         HM[M120] += D::c[i][0] * Hyy * f[i];
     *         HM[M102] += D::c[i][0] * Hzz * f[i];
     *         HM[M012] += D::c[i][1] * Hzz * f[i];
     *         HM[M111] += D::c[i][0] * D::c[i][1] * D::c[i][2] * f[i];
     *         // Order 4
     *         HM[M220] += Hxx * Hyy * f[i];
     *         HM[M202] += Hxx * Hzz * f[i];
     *         HM[M022] += Hyy * Hzz * f[i];
     *         HM[M211] += Hxx * D::c[i][1] * D::c[i][2] * f[i];
     *         HM[M121] += D::c[i][0] * Hyy * D::c[i][2] * f[i];
     *         HM[M112] += D::c[i][0] * D::c[i][1] * Hzz * f[i];
     *         // Order 5
     *         HM[M221] += Hxx * Hyy * D::c[i][2] * f[i];
     *         HM[M212] += Hxx * D::c[i][1] * Hzz * f[i];
     *         HM[M122] += D::c[i][0] * Hyy * Hzz * f[i];
     *         // Order 6
     *         HM[M222] += Hxx * Hyy * Hzz * f[i];
     *     }
     *
     *     T invRho = 1. / HM[M000];
     *     for (int i = 0; i<27; ++i) {
     *         HM[i] *= invRho;
     *     }
     * };
     *
     * // Optimized way to compute HMs based on Palabos ordering of discrete velocities
     * static void HMcomputeMoments(Array<T,D::q> const& cell, Array<T, D::q>& HM, T& rho) {
     *
     *     Array<T, D::q> f;
     *     for (int i = 0; i<27; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *     }
     *
     *     T a1 = 1./3. ;// D::cs2
     *     T a2 = 2./3. ;// 1.-D::cs2
     *
     *     T b1 = 1./9. ;// (D::cs2)*(D::cs2)
     *     T b2 = 2./9. ;// (1.-D::cs2)*(D::cs2)
     *     T b3 = 4./9. ;// (1.-D::cs2)*(1.-D::cs2)
     *
     *     T c1 = 1./27.;// (D::cs2)*(D::cs2)*(D::cs2)
     *     T c2 = 2./27.;// (1.-D::cs2)*(D::cs2)*(D::cs2)
     *     T c3 = 4./27.;// (1.-D::cs2)*(1.-D::cs2)*(D::cs2)
     *     T c4 = 8./27.;// (1.-D::cs2)*(1.-D::cs2)*(1.-D::cs2)
     *
     *     // Order 0
     *     rho =        f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[10]
     * + f[11] + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18] + f[19] + f[20] + f[21] +
     * f[22] + f[23] + f[24] + f[25] + f[26]; HM[M000] = 1.0; T invRho = 1./rho;
     *     // Order 1
     *     HM[M100] = invRho * (- f[1] - f[4] - f[5] - f[6] - f[7] - f[10] - f[11] - f[12] - f[13] +
     * f[14] + f[17] + f[18] + f[19] + f[20] + f[23] + f[24] + f[25] + f[26]); HM[M010] = invRho *
     * (- f[2] - f[4] + f[5] - f[8] - f[9] - f[10] - f[11] + f[12] + f[13] + f[15] + f[17] - f[18] +
     * f[21] + f[22] + f[23] + f[24] - f[25] - f[26]); HM[M001] = invRho * (- f[3] - f[6] + f[7] -
     * f[8] + f[9] - f[10] + f[11] - f[12] + f[13] + f[16] + f[19] - f[20] + f[21] - f[22] + f[23] -
     * f[24] + f[25] - f[26]);
     *     // Order 2
     *     HM[M200] = invRho * (- a1*f[0] + a2*f[1] - a1*f[2] - a1*f[3] + a2*f[4] + a2*f[5] +
     * a2*f[6] + a2*f[7] - a1*f[8] - a1*f[9] + a2*f[10] + a2*f[11] + a2*f[12] + a2*f[13] + a2*f[14]
     * - a1*f[15] - a1*f[16] + a2*f[17] + a2*f[18] + a2*f[19] + a2*f[20] - a1*f[21] - a1*f[22] +
     * a2*f[23] + a2*f[24] + a2*f[25] + a2*f[26]); HM[M020] = invRho * (- a1*f[0] - a1*f[1] +
     * a2*f[2] - a1*f[3] + a2*f[4] + a2*f[5] - a1*f[6] - a1*f[7] + a2*f[8] + a2*f[9] + a2*f[10] +
     * a2*f[11] + a2*f[12] + a2*f[13] - a1*f[14] + a2*f[15] - a1*f[16] + a2*f[17] + a2*f[18] -
     * a1*f[19] - a1*f[20] + a2*f[21] + a2*f[22] + a2*f[23] + a2*f[24] + a2*f[25] + a2*f[26]);
     *     HM[M002] = invRho * (- a1*f[0] - a1*f[1] - a1*f[2] + a2*f[3] - a1*f[4] - a1*f[5] +
     * a2*f[6] + a2*f[7] + a2*f[8] + a2*f[9] + a2*f[10] + a2*f[11] + a2*f[12] + a2*f[13] - a1*f[14]
     * - a1*f[15] + a2*f[16] - a1*f[17] - a1*f[18] + a2*f[19] + a2*f[20] + a2*f[21] + a2*f[22] +
     * a2*f[23] + a2*f[24] + a2*f[25] + a2*f[26]); HM[M110] = invRho * (f[4] - f[5] + f[10] + f[11]
     * - f[12] - f[13] + f[17] - f[18] + f[23] + f[24] - f[25] - f[26]); HM[M101] = invRho * (f[6] -
     * f[7] + f[10] - f[11] + f[12] - f[13] + f[19] - f[20] + f[23] - f[24] + f[25] - f[26]);
     *     HM[M011] = invRho * (f[8] - f[9] + f[10] - f[11] - f[12] + f[13] + f[21] - f[22] + f[23]
     * - f[24] - f[25] + f[26]);
     *     // Order 3
     *     HM[M210] = invRho * (a1*f[2] - a2*f[4] + a2*f[5] + a1*f[8] + a1*f[9] - a2*f[10] -
     * a2*f[11] + a2*f[12] + a2*f[13] - a1*f[15] + a2*f[17] - a2*f[18] - a1*f[21] - a1*f[22] +
     * a2*f[23] + a2*f[24] - a2*f[25] - a2*f[26]); HM[M201] = invRho * (a1*f[3] - a2*f[6] + a2*f[7]
     * + a1*f[8] - a1*f[9] - a2*f[10] + a2*f[11] - a2*f[12] + a2*f[13] - a1*f[16] + a2*f[19] -
     * a2*f[20] - a1*f[21] + a1*f[22] + a2*f[23] - a2*f[24] + a2*f[25] - a2*f[26]); HM[M021] =
     * invRho * (a1*f[3] + a1*f[6] - a1*f[7] - a2*f[8] + a2*f[9] - a2*f[10] + a2*f[11] - a2*f[12] +
     * a2*f[13] - a1*f[16] - a1*f[19] + a1*f[20] + a2*f[21] - a2*f[22] + a2*f[23] - a2*f[24] +
     * a2*f[25] - a2*f[26]); HM[M120] = invRho * (a1*f[1] - a2*f[4] - a2*f[5] + a1*f[6] + a1*f[7] -
     * a2*f[10] - a2*f[11] - a2*f[12] - a2*f[13] - a1*f[14] + a2*f[17] + a2*f[18] - a1*f[19] -
     * a1*f[20] + a2*f[23] + a2*f[24] + a2*f[25] + a2*f[26]); HM[M102] = invRho * (a1*f[1] + a1*f[4]
     * + a1*f[5] - a2*f[6] - a2*f[7] - a2*f[10] - a2*f[11] - a2*f[12] - a2*f[13] - a1*f[14] -
     * a1*f[17] - a1*f[18] + a2*f[19] + a2*f[20] + a2*f[23] + a2*f[24] + a2*f[25] + a2*f[26]);
     *     HM[M012] = invRho * (a1*f[2] + a1*f[4] - a1*f[5] - a2*f[8] - a2*f[9] - a2*f[10] -
     * a2*f[11] + a2*f[12] + a2*f[13] - a1*f[15] - a1*f[17] + a1*f[18] + a2*f[21] + a2*f[22] +
     * a2*f[23] + a2*f[24] - a2*f[25] - a2*f[26]); HM[M111] = invRho * (- f[10] + f[11] + f[12] -
     * f[13] + f[23] - f[24] - f[25] + f[26]);
     *     // Order 4
     *     HM[M220] = invRho * (b1*f[0] - b2*f[1] - b2*f[2] + b1*f[3] + b3*f[4] + b3*f[5] - b2*f[6]
     * - b2*f[7] - b2*f[8] - b2*f[9] + b3*f[10] + b3*f[11] + b3*f[12] + b3*f[13] - b2*f[14] -
     * b2*f[15] + b1*f[16] + b3*f[17] + b3*f[18] - b2*f[19] - b2*f[20] - b2*f[21] - b2*f[22] +
     * b3*f[23] + b3*f[24] + b3*f[25] + b3*f[26]); HM[M202] = invRho * (b1*f[0] - b2*f[1] + b1*f[2]
     * - b2*f[3] - b2*f[4] - b2*f[5] + b3*f[6] + b3*f[7] - b2*f[8] - b2*f[9] + b3*f[10] + b3*f[11] +
     * b3*f[12] + b3*f[13] - b2*f[14] + b1*f[15] - b2*f[16] - b2*f[17] - b2*f[18] + b3*f[19] +
     * b3*f[20] - b2*f[21] - b2*f[22] + b3*f[23] + b3*f[24] + b3*f[25] + b3*f[26]); HM[M022] =
     * invRho * (b1*f[0] + b1*f[1] - b2*f[2] - b2*f[3] - b2*f[4] - b2*f[5] - b2*f[6] - b2*f[7] +
     * b3*f[8] + b3*f[9] + b3*f[10] + b3*f[11] + b3*f[12] + b3*f[13] + b1*f[14] - b2*f[15] -
     * b2*f[16] - b2*f[17] - b2*f[18] - b2*f[19] - b2*f[20] + b3*f[21] + b3*f[22] + b3*f[23] +
     * b3*f[24] + b3*f[25] + b3*f[26]); HM[M211] = invRho * (- a1*f[8] + a1*f[9] + a2*f[10] -
     * a2*f[11] - a2*f[12] + a2*f[13] - a1*f[21] + a1*f[22] + a2*f[23] - a2*f[24] - a2*f[25] +
     * a2*f[26]); HM[M121] = invRho * (- a1*f[6] + a1*f[7] + a2*f[10] - a2*f[11] + a2*f[12] -
     * a2*f[13] - a1*f[19] + a1*f[20] + a2*f[23] - a2*f[24] + a2*f[25] - a2*f[26]); HM[M112] =
     * invRho * (- a1*f[4] + a1*f[5] + a2*f[10] + a2*f[11] - a2*f[12] - a2*f[13] - a1*f[17] +
     * a1*f[18] + a2*f[23] + a2*f[24] - a2*f[25] - a2*f[26]);
     *     // Order 5
     *     HM[M221] = invRho * (- b1*f[3] + b2*f[6] - b2*f[7] + b2*f[8] - b2*f[9] - b3*f[10] +
     * b3*f[11] - b3*f[12] + b3*f[13] + b1*f[16] - b2*f[19] + b2*f[20] - b2*f[21] + b2*f[22] +
     * b3*f[23] - b3*f[24] + b3*f[25] - b3*f[26]); HM[M212] = invRho * (- b1*f[2] + b2*f[4] -
     * b2*f[5] + b2*f[8] + b2*f[9] - b3*f[10] - b3*f[11] + b3*f[12] + b3*f[13] + b1*f[15] - b2*f[17]
     * + b2*f[18] - b2*f[21] - b2*f[22] + b3*f[23] + b3*f[24] - b3*f[25] - b3*f[26]); HM[M122] =
     * invRho * (- b1*f[1] + b2*f[4] + b2*f[5] + b2*f[6] + b2*f[7] - b3*f[10] - b3*f[11] - b3*f[12]
     * - b3*f[13] + b1*f[14] - b2*f[17] - b2*f[18] - b2*f[19] - b2*f[20] + b3*f[23] + b3*f[24] +
     * b3*f[25] + b3*f[26]);
     *     // Order 6
     *     HM[M222] = invRho * (- c1*f[0] + c2*f[1] + c2*f[2] + c2*f[3] - c3*f[4] - c3*f[5] -
     * c3*f[6] - c3*f[7] - c3*f[8] - c3*f[9] + c4*f[10] + c4*f[11] + c4*f[12] + c4*f[13] + c2*f[14]
     * + c2*f[15] + c2*f[16] - c3*f[17] - c3*f[18] - c3*f[19] - c3*f[20] - c3*f[21] - c3*f[22] +
     * c4*f[23] + c4*f[24] + c4*f[25] + c4*f[26]);
     * };
     */

    // Optimized way to compute HMs based on the general ordering of discrete velocities
    static void HMcomputeMoments(Array<T, D::q> const &cell, Array<T, D::q> &HM, T &rho)
    {
        Array<T, D::q> f;
        for (int i = 0; i < 27; ++i) {
            f[i] = cell[i] + D::SkordosFactor() * D::t[i];
        }

        T A1 = f[FMMM] + f[FMMP] + f[FMPM] + f[FMPP];
        T A2 = f[FPPP] + f[FPPM] + f[FPMP] + f[FPMM];
        T A3 = f[FMMM] + f[FMMP] - f[FMPM] - f[FMPP];
        T A4 = f[FPPP] + f[FPPM] - f[FPMP] - f[FPMM];
        T A5 = f[FMMM] - f[FMMP] + f[FMPM] - f[FMPP];
        T A6 = f[FPPP] - f[FPPM] + f[FPMP] - f[FPMM];
        T A7 = f[FMMM] - f[FMMP] - f[FMPM] + f[FMPP];
        T A8 = f[FPPP] - f[FPPM] - f[FPMP] + f[FPMM];

        T X_M1 = f[FM00] + f[FMM0] + f[FMP0] + f[FM0M] + f[FM0P] + A1;
        T X_P1 = f[FP00] + f[FPP0] + f[FPM0] + f[FP0P] + f[FP0M] + A2;
        T X_0 =
            f[F0M0] + f[F00M] + f[F0MM] + f[F0MP] + f[F000] + f[F0P0] + f[F00P] + f[F0PP] + f[F0PM];

        T Y_M1 =
            f[F0M0] + f[FMM0] + f[F0MM] + f[F0MP] + f[FMMM] + f[FMMP] + f[FPM0] + f[FPMP] + f[FPMM];
        T Y_P1 =
            f[F0P0] + f[FPP0] + f[F0PP] + f[F0PM] + f[FPPP] + f[FPPM] + f[FMP0] + f[FMPM] + f[FMPP];

        T Z_M1 =
            f[F00M] + f[FM0M] + f[F0MM] + f[FMMM] + f[FMPM] + f[FP0M] + f[F0PM] + f[FPPM] + f[FPMM];
        T Z_P1 =
            f[F00P] + f[FP0P] + f[F0PP] + f[FPPP] + f[FPMP] + f[FM0P] + f[F0MP] + f[FMMP] + f[FMPP];

        rho = X_M1 + X_P1 + X_0;
        T invRho = 1. / rho;

        // Order 6
        HM[M222] = invRho * (A1 + A2);
        // Order 5
        HM[M221] = invRho * (-A5 + A6);
        HM[M212] = invRho * (-A3 + A4);
        HM[M122] = invRho * (-A1 + A2);
        // Order 4
        HM[M220] = invRho * (f[FMM0] + f[FMP0] + f[FPP0] + f[FPM0]) + HM[M222];
        HM[M202] = invRho * (f[FM0M] + f[FM0P] + f[FP0P] + f[FP0M]) + HM[M222];
        HM[M022] = invRho * (f[F0MM] + f[F0MP] + f[F0PP] + f[F0PM]) + HM[M222];
        HM[M211] = invRho * (A7 + A8);
        HM[M121] = invRho * (A5 + A6);
        HM[M112] = invRho * (A3 + A4);
        // Order 3
        HM[M210] = invRho * (-f[FMM0] + f[FMP0] + f[FPP0] - f[FPM0]) + HM[M212];
        HM[M201] = invRho * (-f[FM0M] + f[FM0P] + f[FP0P] - f[FP0M]) + HM[M221];
        HM[M021] = invRho * (-f[F0MM] + f[F0MP] + f[F0PP] - f[F0PM]) + HM[M221];
        HM[M120] = invRho * (-f[FMM0] - f[FMP0] + f[FPP0] + f[FPM0]) + HM[M122];
        HM[M102] = invRho * (-f[FM0M] - f[FM0P] + f[FP0P] + f[FP0M]) + HM[M122];
        HM[M012] = invRho * (-f[F0MM] - f[F0MP] + f[F0PP] + f[F0PM]) + HM[M212];
        HM[M111] = invRho * (-A7 + A8);
        // Order 2
        HM[M200] = invRho * (X_P1 + X_M1);
        HM[M020] = invRho * (Y_P1 + Y_M1);
        HM[M002] = invRho * (Z_P1 + Z_M1);
        HM[M110] = invRho * (f[FMM0] - f[FMP0] + f[FPP0] - f[FPM0]) + HM[M112];
        HM[M101] = invRho * (f[FM0M] - f[FM0P] + f[FP0P] - f[FP0M]) + HM[M121];
        HM[M011] = invRho * (f[F0MM] - f[F0MP] + f[F0PP] - f[F0PM]) + HM[M211];
        // Order 1
        HM[M100] = invRho * (X_P1 - X_M1);
        HM[M010] = invRho * (Y_P1 - Y_M1);
        HM[M001] = invRho * (Z_P1 - Z_M1);
        // Order 0
        HM[M000] = 1.;

        // We come back to Hermite moments
        T cs4 = D::cs2 * D::cs2;
        HM[M200] -= D::cs2;
        HM[M020] -= D::cs2;
        HM[M002] -= D::cs2;

        HM[M210] -= D::cs2 * HM[M010];
        HM[M201] -= D::cs2 * HM[M001];
        HM[M021] -= D::cs2 * HM[M001];
        HM[M120] -= D::cs2 * HM[M100];
        HM[M102] -= D::cs2 * HM[M100];
        HM[M012] -= D::cs2 * HM[M010];

        HM[M220] -= (D::cs2 * (HM[M200] + HM[M020]) + cs4);
        HM[M202] -= (D::cs2 * (HM[M200] + HM[M002]) + cs4);
        HM[M022] -= (D::cs2 * (HM[M020] + HM[M002]) + cs4);
        HM[M211] -= (D::cs2 * HM[M011]);
        HM[M121] -= (D::cs2 * HM[M101]);
        HM[M112] -= (D::cs2 * HM[M110]);

        HM[M221] -= (D::cs2 * (HM[M201] + HM[M021]) + cs4 * HM[M001]);
        HM[M212] -= (D::cs2 * (HM[M210] + HM[M012]) + cs4 * HM[M010]);
        HM[M122] -= (D::cs2 * (HM[M120] + HM[M102]) + cs4 * HM[M100]);

        HM[M222] -=
            (D::cs2 * (HM[M220] + HM[M202] + HM[M022]) + cs4 * (HM[M200] + HM[M020] + HM[M002])
             + D::cs2 * cs4);
    };

    static void HMcomputeEquilibriumMoments(Array<T, D::d> const &u, Array<T, D::q> &HMeq)
    {
        // Order 0
        HMeq[M000] = 1.;
        // Order 1
        HMeq[M100] = u[0];
        HMeq[M010] = u[1];
        HMeq[M001] = u[2];
        // Order 2
        HMeq[M200] = u[0] * u[0];
        HMeq[M020] = u[1] * u[1];
        HMeq[M002] = u[2] * u[2];
        HMeq[M110] = u[0] * u[1];
        HMeq[M101] = u[0] * u[2];
        HMeq[M011] = u[1] * u[2];
        // Order 3
        HMeq[M210] = HMeq[M200] * u[1];
        HMeq[M201] = HMeq[M200] * u[2];
        HMeq[M021] = HMeq[M020] * u[2];
        HMeq[M120] = HMeq[M020] * u[0];
        HMeq[M102] = HMeq[M002] * u[0];
        HMeq[M012] = HMeq[M002] * u[1];
        HMeq[M111] = HMeq[M110] * u[2];
        // Order 4
        HMeq[M220] = HMeq[M200] * HMeq[M020];
        HMeq[M202] = HMeq[M200] * HMeq[M002];
        HMeq[M022] = HMeq[M020] * HMeq[M002];
        HMeq[M211] = HMeq[M200] * HMeq[M011];
        HMeq[M121] = HMeq[M020] * HMeq[M101];
        HMeq[M112] = HMeq[M002] * HMeq[M110];
        // Order 5
        HMeq[M221] = HMeq[M220] * HMeq[M001];
        HMeq[M212] = HMeq[M202] * HMeq[M010];
        HMeq[M122] = HMeq[M022] * HMeq[M100];
        // Order 6
        HMeq[M222] = HMeq[M220] * HMeq[M002];
    };

    // Equilibrium populations based on 27 moments can be computed using either RM, HM, CM, CHM or
    // Gauss-Hermite formalisms. All formulations are equivalent for the D3Q27 so we use the fastest
    // one (RMs)
    static void HMcomputeEquilibrium(T rho, Array<T, D::q> const &HMeq, Array<T, D::q> &eq)
    {
        Array<T, D::d> u(HMeq[1], HMeq[2], HMeq[3]);
        Array<T, D::q> RMeq;
        // Order 2
        RMeq[M200] = u[0] * u[0] + D::cs2;
        RMeq[M020] = u[1] * u[1] + D::cs2;
        RMeq[M002] = u[2] * u[2] + D::cs2;
        RMeq[M110] = u[0] * u[1];
        RMeq[M101] = u[0] * u[2];
        RMeq[M011] = u[1] * u[2];
        // Order 3
        RMeq[M210] = RMeq[M200] * u[1];
        RMeq[M201] = RMeq[M200] * u[2];
        RMeq[M021] = RMeq[M020] * u[2];
        RMeq[M120] = RMeq[M020] * u[0];
        RMeq[M102] = RMeq[M002] * u[0];
        RMeq[M012] = RMeq[M002] * u[1];
        RMeq[M111] = RMeq[M110] * u[2];
        // Order 4
        RMeq[M220] = RMeq[M200] * RMeq[M020];
        RMeq[M202] = RMeq[M200] * RMeq[M002];
        RMeq[M022] = RMeq[M020] * RMeq[M002];
        RMeq[M211] = RMeq[M200] * RMeq[M011];
        RMeq[M121] = RMeq[M020] * RMeq[M101];
        RMeq[M112] = RMeq[M002] * RMeq[M110];
        // Order 5
        RMeq[M221] = RMeq[M220] * u[2];
        RMeq[M212] = RMeq[M202] * u[1];
        RMeq[M122] = RMeq[M022] * u[0];
        // Order 6
        RMeq[M222] = RMeq[M220] * RMeq[M002];

        // Optimization based on symmetries between populations and their opposite counterpart
        eq[F000] = rho
                   * (1. - RMeq[M200] - RMeq[M020] - RMeq[M002] + RMeq[M220] + RMeq[M202]
                      + RMeq[M022] - RMeq[M222]);

        eq[FP00] = 0.5 * rho
                   * (u[0] + RMeq[M200] - RMeq[M120] - RMeq[M102] - RMeq[M220] - RMeq[M202]
                      + RMeq[M122] + RMeq[M222]);
        eq[FM00] = rho * (-u[0] + RMeq[M120] + RMeq[M102] - RMeq[M122]) + eq[FP00];

        eq[F0P0] = 0.5 * rho
                   * (u[1] + RMeq[M020] - RMeq[M210] - RMeq[M012] - RMeq[M220] - RMeq[M022]
                      + RMeq[M212] + RMeq[M222]);
        eq[F0M0] = rho * (-u[1] + RMeq[M210] + RMeq[M012] - RMeq[M212]) + eq[F0P0];

        eq[F00P] = 0.5 * rho
                   * (u[2] + RMeq[M002] - RMeq[M201] - RMeq[M021] - RMeq[M202] - RMeq[M022]
                      + RMeq[M221] + RMeq[M222]);
        eq[F00M] = rho * (-u[2] + RMeq[M201] + RMeq[M021] - RMeq[M221]) + eq[F00P];

        eq[FPP0] = 0.25 * rho
                   * (RMeq[M110] + RMeq[M210] + RMeq[M120] - RMeq[M112] + RMeq[M220] - RMeq[M212]
                      - RMeq[M122] - RMeq[M222]);
        eq[FMP0] = 0.5 * rho * (-RMeq[M110] - RMeq[M120] + RMeq[M112] + RMeq[M122]) + eq[FPP0];
        eq[FPM0] = 0.5 * rho * (-RMeq[M110] - RMeq[M210] + RMeq[M112] + RMeq[M212]) + eq[FPP0];
        eq[FMM0] = 0.5 * rho * (-RMeq[M210] - RMeq[M120] + RMeq[M212] + RMeq[M122]) + eq[FPP0];

        eq[FP0P] = 0.25 * rho
                   * (RMeq[M101] + RMeq[M201] + RMeq[M102] - RMeq[M121] + RMeq[M202] - RMeq[M221]
                      - RMeq[M122] - RMeq[M222]);
        eq[FM0P] = 0.5 * rho * (-RMeq[M101] - RMeq[M102] + RMeq[M121] + RMeq[M122]) + eq[FP0P];
        eq[FP0M] = 0.5 * rho * (-RMeq[M101] - RMeq[M201] + RMeq[M121] + RMeq[M221]) + eq[FP0P];
        eq[FM0M] = 0.5 * rho * (-RMeq[M201] - RMeq[M102] + RMeq[M221] + RMeq[M122]) + eq[FP0P];

        eq[F0PP] = 0.25 * rho
                   * (RMeq[M011] + RMeq[M021] + RMeq[M012] - RMeq[M211] + RMeq[M022] - RMeq[M221]
                      - RMeq[M212] - RMeq[M222]);
        eq[F0MP] = 0.5 * rho * (-RMeq[M011] - RMeq[M012] + RMeq[M211] + RMeq[M212]) + eq[F0PP];
        eq[F0PM] = 0.5 * rho * (-RMeq[M011] - RMeq[M021] + RMeq[M211] + RMeq[M221]) + eq[F0PP];
        eq[F0MM] = 0.5 * rho * (-RMeq[M021] - RMeq[M012] + RMeq[M221] + RMeq[M212]) + eq[F0PP];

        eq[FPPP] = 0.125 * rho
                   * (RMeq[M111] + RMeq[M211] + RMeq[M121] + RMeq[M112] + RMeq[M221] + RMeq[M212]
                      + RMeq[M122] + RMeq[M222]);
        eq[FMPP] = 0.25 * rho * (-RMeq[M111] - RMeq[M121] - RMeq[M112] - RMeq[M122]) + eq[FPPP];
        eq[FPMP] = 0.25 * rho * (-RMeq[M111] - RMeq[M211] - RMeq[M112] - RMeq[M212]) + eq[FPPP];
        eq[FPPM] = 0.25 * rho * (-RMeq[M111] - RMeq[M211] - RMeq[M121] - RMeq[M221]) + eq[FPPP];
        eq[FMMP] = 0.25 * rho * (-RMeq[M211] - RMeq[M121] - RMeq[M212] - RMeq[M122]) + eq[FPPP];
        eq[FMPM] = 0.25 * rho * (-RMeq[M211] - RMeq[M112] - RMeq[M221] - RMeq[M122]) + eq[FPPP];
        eq[FPMM] = 0.25 * rho * (-RMeq[M121] - RMeq[M112] - RMeq[M221] - RMeq[M212]) + eq[FPPP];
        eq[FMMM] = 0.25 * rho * (-RMeq[M111] - RMeq[M221] - RMeq[M212] - RMeq[M122]) + eq[FPPP];
    };

    static void HMcollide(
        Array<T, D::q> &cell, T rho, Array<T, D::d> const &u,
        Array<T, D::q> const &HM,    // Hermite moments
        Array<T, D::q> const &HMeq,  // Equilibrium moments (Hermite)
        Array<T, D::numRelaxationTimes> const &omega)
    {
        T omega1 = omega[0];
        T omega2 = omega[1];
        T omega3 = omega[2];
        T omega4 = omega[3];
        T omega5 = omega[4];
        T omega6 = omega[5];
        T omega7 = omega[6];
        T omega8 = omega[7];

        T omegaBulk = omega[8];
        T omegaPlus = (omegaBulk + 2. * omega1) / 3.;  // Notation used by Fei
        T omegaMinus = (omegaBulk - omega1) / 3.;      // Notation used by Fei

        T cs4 = D::cs2 * D::cs2;

        // Post-collision moments.
        Array<T, D::q> HMcoll;
        Array<T, D::q> RMcoll;

        // Collision in the Hermite moment space
        // Order 2 (non-diagonal collision so that we can easily modify the bulk viscosity)
        HMcoll[M200] = HM[M200] - omegaPlus * (HM[M200] - HMeq[M200])
                       - omegaMinus * (HM[M020] - HMeq[M020])
                       - omegaMinus * (HM[M002] - HMeq[M002]);
        HMcoll[M020] = HM[M020] - omegaMinus * (HM[M200] - HMeq[M200])
                       - omegaPlus * (HM[M020] - HMeq[M020]) - omegaMinus * (HM[M002] - HMeq[M002]);
        HMcoll[M002] = HM[M002] - omegaMinus * (HM[M200] - HMeq[M200])
                       - omegaMinus * (HM[M020] - HMeq[M020]) - omegaPlus * (HM[M002] - HMeq[M002]);
        HMcoll[M110] = (1. - omega2) * HM[M110] + omega2 * HMeq[M110];
        HMcoll[M101] = (1. - omega2) * HM[M101] + omega2 * HMeq[M101];
        HMcoll[M011] = (1. - omega2) * HM[M011] + omega2 * HMeq[M011];
        // Order 3
        HMcoll[M210] = (1. - omega3) * HM[M210] + omega3 * HMeq[M210];
        HMcoll[M201] = (1. - omega3) * HM[M201] + omega3 * HMeq[M201];
        HMcoll[M021] = (1. - omega3) * HM[M021] + omega3 * HMeq[M021];
        HMcoll[M120] = (1. - omega3) * HM[M120] + omega3 * HMeq[M120];
        HMcoll[M102] = (1. - omega3) * HM[M102] + omega3 * HMeq[M102];
        HMcoll[M012] = (1. - omega3) * HM[M012] + omega3 * HMeq[M012];
        HMcoll[M111] = (1. - omega4) * HM[M111] + omega4 * HMeq[M111];
        // Order 4
        HMcoll[M220] = (1. - omega5) * HM[M220] + omega5 * HMeq[M220];
        HMcoll[M202] = (1. - omega5) * HM[M202] + omega5 * HMeq[M202];
        HMcoll[M022] = (1. - omega5) * HM[M022] + omega5 * HMeq[M022];
        HMcoll[M211] = (1. - omega6) * HM[M211] + omega6 * HMeq[M211];
        HMcoll[M121] = (1. - omega6) * HM[M121] + omega6 * HMeq[M121];
        HMcoll[M112] = (1. - omega6) * HM[M112] + omega6 * HMeq[M112];
        // Order 5
        HMcoll[M221] = (1. - omega7) * HM[M221] + omega7 * HMeq[M221];
        HMcoll[M212] = (1. - omega7) * HM[M212] + omega7 * HMeq[M212];
        HMcoll[M122] = (1. - omega7) * HM[M122] + omega7 * HMeq[M122];
        // Order 6
        HMcoll[M222] = (1. - omega8) * HM[M222] + omega8 * HMeq[M222];

        // Come back to RMcoll using relationships between HMs and RMs
        RMcoll[M200] = HMcoll[M200] + D::cs2;
        RMcoll[M020] = HMcoll[M020] + D::cs2;
        RMcoll[M002] = HMcoll[M002] + D::cs2;

        RMcoll[M110] = HMcoll[M110];
        RMcoll[M101] = HMcoll[M101];
        RMcoll[M011] = HMcoll[M011];

        RMcoll[M210] = HMcoll[M210] + D::cs2 * u[1];
        RMcoll[M201] = HMcoll[M201] + D::cs2 * u[2];
        RMcoll[M021] = HMcoll[M021] + D::cs2 * u[2];
        RMcoll[M120] = HMcoll[M120] + D::cs2 * u[0];
        RMcoll[M102] = HMcoll[M102] + D::cs2 * u[0];
        RMcoll[M012] = HMcoll[M012] + D::cs2 * u[1];

        RMcoll[M111] = HMcoll[M111];

        RMcoll[M220] = HMcoll[M220] + D::cs2 * (HMcoll[M200] + HMcoll[M020]) + cs4;
        RMcoll[M202] = HMcoll[M202] + D::cs2 * (HMcoll[M200] + HMcoll[M002]) + cs4;
        RMcoll[M022] = HMcoll[M022] + D::cs2 * (HMcoll[M020] + HMcoll[M002]) + cs4;

        RMcoll[M211] = HMcoll[M211] + D::cs2 * HMcoll[M011];
        RMcoll[M121] = HMcoll[M121] + D::cs2 * HMcoll[M101];
        RMcoll[M112] = HMcoll[M112] + D::cs2 * HMcoll[M110];

        RMcoll[M221] = HMcoll[M221] + D::cs2 * (HMcoll[M201] + HMcoll[M021]) + cs4 * u[2];
        RMcoll[M212] = HMcoll[M212] + D::cs2 * (HMcoll[M210] + HMcoll[M012]) + cs4 * u[1];
        RMcoll[M122] = HMcoll[M122] + D::cs2 * (HMcoll[M120] + HMcoll[M102]) + cs4 * u[0];

        RMcoll[M222] = HMcoll[M222] + D::cs2 * (HMcoll[M220] + HMcoll[M202] + HMcoll[M022])
                       + cs4 * (HMcoll[M200] + HMcoll[M020] + HMcoll[M002]) + D::cs2 * cs4;

        // Compute post collision populations from RM
        // Optimization based on symmetries between populations and their opposite counterpart
        cell[F000] = rho
                     * (1. - RMcoll[M200] - RMcoll[M020] - RMcoll[M002] + RMcoll[M220]
                        + RMcoll[M202] + RMcoll[M022] - RMcoll[M222]);

        cell[FP00] = 0.5 * rho
                     * (u[0] + RMcoll[M200] - RMcoll[M120] - RMcoll[M102] - RMcoll[M220]
                        - RMcoll[M202] + RMcoll[M122] + RMcoll[M222]);
        cell[FM00] = rho * (-u[0] + RMcoll[M120] + RMcoll[M102] - RMcoll[M122]) + cell[FP00];

        cell[F0P0] = 0.5 * rho
                     * (u[1] + RMcoll[M020] - RMcoll[M210] - RMcoll[M012] - RMcoll[M220]
                        - RMcoll[M022] + RMcoll[M212] + RMcoll[M222]);
        cell[F0M0] = rho * (-u[1] + RMcoll[M210] + RMcoll[M012] - RMcoll[M212]) + cell[F0P0];

        cell[F00P] = 0.5 * rho
                     * (u[2] + RMcoll[M002] - RMcoll[M201] - RMcoll[M021] - RMcoll[M202]
                        - RMcoll[M022] + RMcoll[M221] + RMcoll[M222]);
        cell[F00M] = rho * (-u[2] + RMcoll[M201] + RMcoll[M021] - RMcoll[M221]) + cell[F00P];

        cell[FPP0] = 0.25 * rho
                     * (RMcoll[M110] + RMcoll[M210] + RMcoll[M120] - RMcoll[M112] + RMcoll[M220]
                        - RMcoll[M212] - RMcoll[M122] - RMcoll[M222]);
        cell[FMP0] =
            0.5 * rho * (-RMcoll[M110] - RMcoll[M120] + RMcoll[M112] + RMcoll[M122]) + cell[FPP0];
        cell[FPM0] =
            0.5 * rho * (-RMcoll[M110] - RMcoll[M210] + RMcoll[M112] + RMcoll[M212]) + cell[FPP0];
        cell[FMM0] =
            0.5 * rho * (-RMcoll[M210] - RMcoll[M120] + RMcoll[M212] + RMcoll[M122]) + cell[FPP0];

        cell[FP0P] = 0.25 * rho
                     * (RMcoll[M101] + RMcoll[M201] + RMcoll[M102] - RMcoll[M121] + RMcoll[M202]
                        - RMcoll[M221] - RMcoll[M122] - RMcoll[M222]);
        cell[FM0P] =
            0.5 * rho * (-RMcoll[M101] - RMcoll[M102] + RMcoll[M121] + RMcoll[M122]) + cell[FP0P];
        cell[FP0M] =
            0.5 * rho * (-RMcoll[M101] - RMcoll[M201] + RMcoll[M121] + RMcoll[M221]) + cell[FP0P];
        cell[FM0M] =
            0.5 * rho * (-RMcoll[M201] - RMcoll[M102] + RMcoll[M221] + RMcoll[M122]) + cell[FP0P];

        cell[F0PP] = 0.25 * rho
                     * (RMcoll[M011] + RMcoll[M021] + RMcoll[M012] - RMcoll[M211] + RMcoll[M022]
                        - RMcoll[M221] - RMcoll[M212] - RMcoll[M222]);
        cell[F0MP] =
            0.5 * rho * (-RMcoll[M011] - RMcoll[M012] + RMcoll[M211] + RMcoll[M212]) + cell[F0PP];
        cell[F0PM] =
            0.5 * rho * (-RMcoll[M011] - RMcoll[M021] + RMcoll[M211] + RMcoll[M221]) + cell[F0PP];
        cell[F0MM] =
            0.5 * rho * (-RMcoll[M021] - RMcoll[M012] + RMcoll[M221] + RMcoll[M212]) + cell[F0PP];

        cell[FPPP] = 0.125 * rho
                     * (RMcoll[M111] + RMcoll[M211] + RMcoll[M121] + RMcoll[M112] + RMcoll[M221]
                        + RMcoll[M212] + RMcoll[M122] + RMcoll[M222]);
        cell[FMPP] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M121] - RMcoll[M112] - RMcoll[M122]) + cell[FPPP];
        cell[FPMP] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M211] - RMcoll[M112] - RMcoll[M212]) + cell[FPPP];
        cell[FPPM] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M211] - RMcoll[M121] - RMcoll[M221]) + cell[FPPP];
        cell[FMMP] =
            0.25 * rho * (-RMcoll[M211] - RMcoll[M121] - RMcoll[M212] - RMcoll[M122]) + cell[FPPP];
        cell[FMPM] =
            0.25 * rho * (-RMcoll[M211] - RMcoll[M112] - RMcoll[M221] - RMcoll[M122]) + cell[FPPP];
        cell[FPMM] =
            0.25 * rho * (-RMcoll[M121] - RMcoll[M112] - RMcoll[M221] - RMcoll[M212]) + cell[FPPP];
        cell[FMMM] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M221] - RMcoll[M212] - RMcoll[M122]) + cell[FPPP];

        for (int i = 0; i < 27; ++i) {
            cell[i] -= D::SkordosFactor() * D::t[i];
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////
    // Central Moments Formalism (Equilibrium is computed through raw moments formalism) //
    ///////////////////////////////////////////////////////////////////////////////////////

    /**
     * // General way to compute CMs
     * static void CMcomputeMoments(Array<T,D::q> const& cell, Array<T, D::q>& CM, T& rho,
     * Array<T,D::d>& u) {
     *
     *     Array<T, D::q> f;
     *     for (int i = 0; i<27; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *         CM[i] = 0.;
     *     }
     *
     *     rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[10] + f[11]
     * + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18] + f[19] + f[20] + f[21] + f[22] +
     * f[23] + f[24] + f[25] + f[26]; CM[M000] = 1.0; T invRho = 1./rho; u[0] = invRho * (- f[1] -
     * f[4] - f[5] - f[6] - f[7] - f[10] - f[11] - f[12] - f[13] + f[14] + f[17] + f[18] + f[19] +
     * f[20] + f[23] + f[24] + f[25] + f[26]); u[1] = invRho * (- f[2] - f[4] + f[5] - f[8] - f[9] -
     * f[10] - f[11] + f[12] + f[13] + f[15] + f[17] - f[18] + f[21] + f[22] + f[23] + f[24] - f[25]
     * - f[26]); u[2] = invRho * (- f[3] - f[6] + f[7] - f[8] + f[9] - f[10] + f[11] - f[12] + f[13]
     * + f[16] + f[19] - f[20] + f[21] - f[22] + f[23] - f[24] + f[25] - f[26]); T cMux = 0.; T cMuy
     * = 0.; T cMuz = 0.;
     *
     *     for (int i = 0; i<27; ++i) {
     *         cMux = D::c[i][0]- u[0];
     *         cMuy = D::c[i][1]- u[1];
     *         cMuz = D::c[i][2]- u[2];
     *         // Order 1
     *         CM[M100] += cMux * f[i];
     *         CM[M010] += cMuy * f[i];
     *         CM[M001] += cMuz * f[i];
     *         // Order 2
     *         CM[M200] += cMux * cMux * f[i];
     *         CM[M020] += cMuy * cMuy * f[i];
     *         CM[M002] += cMuz * cMuz * f[i];
     *         CM[M110] += cMux * cMuy * f[i];
     *         CM[M101] += cMux * cMuz * f[i];
     *         CM[M011] += cMuy * cMuz * f[i];
     *         // Order 3
     *         CM[M210] += cMux * cMux * cMuy * f[i];
     *         CM[M201] += cMux * cMux * cMuz * f[i];
     *         CM[M021] += cMuy * cMuy * cMuz * f[i];
     *         CM[M120] += cMux * cMuy * cMuy * f[i];
     *         CM[M102] += cMux * cMuz * cMuz * f[i];
     *         CM[M012] += cMuy * cMuz * cMuz * f[i];
     *         CM[M111] += cMux * cMuy * cMuz * f[i];
     *         // Order 4
     *         CM[M220] += cMux * cMux * cMuy * cMuy * f[i];
     *         CM[M202] += cMux * cMux * cMuz * cMuz * f[i];
     *         CM[M022] += cMuy * cMuy * cMuz * cMuz * f[i];
     *         CM[M211] += cMux * cMux * cMuy * cMuz * f[i];
     *         CM[M121] += cMux * cMuy * cMuy * cMuz * f[i];
     *         CM[M112] += cMux * cMuy * cMuz * cMuz * f[i];
     *         // Order 5
     *         CM[M221] += cMux * cMux * cMuy * cMuy * cMuz * f[i];
     *         CM[M212] += cMux * cMux * cMuy * cMuz * cMuz * f[i];
     *         CM[M122] += cMux * cMuy * cMuy * cMuz * cMuz * f[i];
     *         // Order 6
     *         CM[M222] += cMux * cMux * cMuy * cMuy * cMuz * cMuz * f[i];
     *     }
     *
     *     for (int i = 1; i<27; ++i) {
     *         CM[i] *= invRho;
     *     }
     * };
     */

    // Optimized way to compute CMs based on the general ordering of discrete velocities
    static void CMcomputeMoments(
        Array<T, D::q> const &cell, Array<T, D::q> &CM, T &rho, Array<T, D::d> &u)
    {
        Array<T, D::q> f;
        for (int i = 0; i < 27; ++i) {
            f[i] = cell[i] + D::SkordosFactor() * D::t[i];
            CM[i] = 0.;
        }
        T A1 = f[FMMM] + f[FMMP] + f[FMPM] + f[FMPP];
        T A2 = f[FPPP] + f[FPPM] + f[FPMP] + f[FPMM];
        T A3 = f[FMMM] + f[FMMP] - f[FMPM] - f[FMPP];
        T A4 = f[FPPP] + f[FPPM] - f[FPMP] - f[FPMM];
        T A5 = f[FMMM] - f[FMMP] + f[FMPM] - f[FMPP];
        T A6 = f[FPPP] - f[FPPM] + f[FPMP] - f[FPMM];
        T A7 = f[FMMM] - f[FMMP] - f[FMPM] + f[FMPP];
        T A8 = f[FPPP] - f[FPPM] - f[FPMP] + f[FPMM];

        T X_M1 = f[FM00] + f[FMM0] + f[FMP0] + f[FM0M] + f[FM0P] + A1;
        T X_P1 = f[FP00] + f[FPP0] + f[FPM0] + f[FP0P] + f[FP0M] + A2;
        T X_0 =
            f[F0M0] + f[F00M] + f[F0MM] + f[F0MP] + f[F000] + f[F0P0] + f[F00P] + f[F0PP] + f[F0PM];

        T Y_M1 =
            f[F0M0] + f[FMM0] + f[F0MM] + f[F0MP] + f[FMMM] + f[FMMP] + f[FPM0] + f[FPMP] + f[FPMM];
        T Y_P1 =
            f[F0P0] + f[FPP0] + f[F0PP] + f[F0PM] + f[FPPP] + f[FPPM] + f[FMP0] + f[FMPM] + f[FMPP];

        T Z_M1 =
            f[F00M] + f[FM0M] + f[F0MM] + f[FMMM] + f[FMPM] + f[FP0M] + f[F0PM] + f[FPPM] + f[FPMM];
        T Z_P1 =
            f[F00P] + f[FP0P] + f[F0PP] + f[FPPP] + f[FPMP] + f[FM0P] + f[F0MP] + f[FMMP] + f[FMPP];

        rho = X_M1 + X_P1 + X_0;
        T invRho = 1. / rho;

        // Order 6
        CM[M222] = invRho * (A1 + A2);
        // Order 5
        CM[M221] = invRho * (-A5 + A6);
        CM[M212] = invRho * (-A3 + A4);
        CM[M122] = invRho * (-A1 + A2);
        // Order 4
        CM[M220] = invRho * (f[FMM0] + f[FMP0] + f[FPP0] + f[FPM0]) + CM[M222];
        CM[M202] = invRho * (f[FM0M] + f[FM0P] + f[FP0P] + f[FP0M]) + CM[M222];
        CM[M022] = invRho * (f[F0MM] + f[F0MP] + f[F0PP] + f[F0PM]) + CM[M222];
        CM[M211] = invRho * (A7 + A8);
        CM[M121] = invRho * (A5 + A6);
        CM[M112] = invRho * (A3 + A4);
        // Order 3
        CM[M210] = invRho * (-f[FMM0] + f[FMP0] + f[FPP0] - f[FPM0]) + CM[M212];
        CM[M201] = invRho * (-f[FM0M] + f[FM0P] + f[FP0P] - f[FP0M]) + CM[M221];
        CM[M021] = invRho * (-f[F0MM] + f[F0MP] + f[F0PP] - f[F0PM]) + CM[M221];
        CM[M120] = invRho * (-f[FMM0] - f[FMP0] + f[FPP0] + f[FPM0]) + CM[M122];
        CM[M102] = invRho * (-f[FM0M] - f[FM0P] + f[FP0P] + f[FP0M]) + CM[M122];
        CM[M012] = invRho * (-f[F0MM] - f[F0MP] + f[F0PP] + f[F0PM]) + CM[M212];
        CM[M111] = invRho * (-A7 + A8);
        // Order 2
        CM[M200] = invRho * (X_P1 + X_M1);
        CM[M020] = invRho * (Y_P1 + Y_M1);
        CM[M002] = invRho * (Z_P1 + Z_M1);
        CM[M110] = invRho * (f[FMM0] - f[FMP0] + f[FPP0] - f[FPM0]) + CM[M112];
        CM[M101] = invRho * (f[FM0M] - f[FM0P] + f[FP0P] - f[FP0M]) + CM[M121];
        CM[M011] = invRho * (f[F0MM] - f[F0MP] + f[F0PP] - f[F0PM]) + CM[M211];
        // Order 1
        CM[M100] = 0.;
        CM[M010] = 0.;
        CM[M001] = 0.;
        // Order 0
        CM[M000] = 1.;

        // Compute CMs from RMs using binomial formulas
        u[0] = invRho * (X_P1 - X_M1);
        u[1] = invRho * (Y_P1 - Y_M1);
        u[2] = invRho * (Z_P1 - Z_M1);
        T ux2 = u[0] * u[0];
        T uy2 = u[1] * u[1];
        T uz2 = u[2] * u[2];
        T uxyz = u[0] * u[1] * u[2];

        CM[M200] -= (ux2);
        CM[M020] -= (uy2);
        CM[M002] -= (uz2);

        CM[M110] -= (u[0] * u[1]);
        CM[M101] -= (u[0] * u[2]);
        CM[M011] -= (u[1] * u[2]);

        CM[M210] -= (u[1] * CM[M200] + 2. * u[0] * CM[M110] + ux2 * u[1]);
        CM[M201] -= (u[2] * CM[M200] + 2. * u[0] * CM[M101] + ux2 * u[2]);
        CM[M021] -= (u[2] * CM[M020] + 2. * u[1] * CM[M011] + uy2 * u[2]);
        CM[M120] -= (u[0] * CM[M020] + 2. * u[1] * CM[M110] + u[0] * uy2);
        CM[M102] -= (u[0] * CM[M002] + 2. * u[2] * CM[M101] + u[0] * uz2);
        CM[M012] -= (u[1] * CM[M002] + 2. * u[2] * CM[M011] + u[1] * uz2);

        CM[M111] -= (u[2] * CM[M110] + u[1] * CM[M101] + u[0] * CM[M011] + uxyz);

        CM[M220] -=
            (2. * u[1] * CM[M210] + 2. * u[0] * CM[M120] + uy2 * CM[M200] + ux2 * CM[M020]
             + 4. * u[0] * u[1] * CM[M110] + ux2 * uy2);
        CM[M202] -=
            (2. * u[2] * CM[M201] + 2. * u[0] * CM[M102] + uz2 * CM[M200] + ux2 * CM[M002]
             + 4. * u[0] * u[2] * CM[M101] + ux2 * uz2);
        CM[M022] -=
            (2. * u[2] * CM[M021] + 2. * u[1] * CM[M012] + uz2 * CM[M020] + uy2 * CM[M002]
             + 4. * u[1] * u[2] * CM[M011] + uy2 * uz2);

        CM[M211] -=
            (u[2] * CM[M210] + u[1] * CM[M201] + 2. * u[0] * CM[M111] + u[1] * u[2] * CM[M200]
             + 2. * u[0] * u[2] * CM[M110] + 2. * u[0] * u[1] * CM[M101] + ux2 * CM[M011]
             + ux2 * u[1] * u[2]);
        CM[M121] -=
            (u[2] * CM[M120] + u[0] * CM[M021] + 2. * u[1] * CM[M111] + u[0] * u[2] * CM[M020]
             + 2. * u[1] * u[2] * CM[M110] + 2. * u[0] * u[1] * CM[M011] + uy2 * CM[M101]
             + u[0] * uy2 * u[2]);
        CM[M112] -=
            (u[1] * CM[M102] + u[0] * CM[M012] + 2. * u[2] * CM[M111] + u[0] * u[1] * CM[M002]
             + 2. * u[1] * u[2] * CM[M101] + 2. * u[0] * u[2] * CM[M011] + uz2 * CM[M110]
             + u[0] * u[1] * uz2);

        CM[M221] -=
            (u[2] * CM[M220] + 2. * u[1] * CM[M211] + 2. * u[0] * CM[M121]
             + 2. * u[1] * u[2] * CM[M210] + uy2 * CM[M201] + ux2 * CM[M021]
             + 2. * u[0] * u[2] * CM[M120] + 4. * u[0] * u[1] * CM[M111] + uy2 * u[2] * CM[M200]
             + ux2 * u[2] * CM[M020] + 4. * uxyz * CM[M110] + 2. * u[0] * uy2 * CM[M101]
             + 2. * ux2 * u[1] * CM[M011] + ux2 * uy2 * u[2]);
        CM[M212] -=
            (u[1] * CM[M202] + 2. * u[2] * CM[M211] + 2. * u[0] * CM[M112]
             + 2. * u[1] * u[2] * CM[M201] + uz2 * CM[M210] + ux2 * CM[M012]
             + 2. * u[0] * u[1] * CM[M102] + 4. * u[0] * u[2] * CM[M111] + u[1] * uz2 * CM[M200]
             + ux2 * u[1] * CM[M002] + 4. * uxyz * CM[M101] + 2. * u[0] * uz2 * CM[M110]
             + 2. * ux2 * u[2] * CM[M011] + ux2 * u[1] * uz2);
        CM[M122] -=
            (u[0] * CM[M022] + 2. * u[2] * CM[M121] + 2. * u[1] * CM[M112]
             + 2. * u[0] * u[2] * CM[M021] + uz2 * CM[M120] + uy2 * CM[M102]
             + 2. * u[0] * u[1] * CM[M012] + 4. * u[1] * u[2] * CM[M111] + u[0] * uz2 * CM[M020]
             + u[0] * uy2 * CM[M002] + 4. * uxyz * CM[M011] + 2. * u[1] * uz2 * CM[M110]
             + 2. * uy2 * u[2] * CM[M101] + u[0] * uy2 * uz2);

        CM[M222] -=
            (2. * u[2] * CM[M221] + 2. * u[1] * CM[M212] + 2. * u[0] * CM[M122] + uz2 * CM[M220]
             + uy2 * CM[M202] + ux2 * CM[M022] + 4. * u[1] * u[2] * CM[M211]
             + 4. * u[0] * u[2] * CM[M121] + 4. * u[0] * u[1] * CM[M112]
             + 2. * u[1] * uz2 * CM[M210] + 2. * uy2 * u[2] * CM[M201] + 2. * ux2 * u[2] * CM[M021]
             + 2. * u[0] * uz2 * CM[M120] + 2. * u[0] * uy2 * CM[M102] + 2. * ux2 * u[1] * CM[M012]
             + 8. * uxyz * CM[M111] + uy2 * uz2 * CM[M200] + ux2 * uz2 * CM[M020]
             + ux2 * uy2 * CM[M002] + 4. * u[0] * u[1] * uz2 * CM[M110]
             + 4. * u[0] * uy2 * u[2] * CM[M101] + 4. * ux2 * u[1] * u[2] * CM[M011]
             + ux2 * uy2 * uz2);
    };

    static void CMcomputeEquilibriumMoments(Array<T, D::q> &CMeq)
    {
        // Order 0
        CMeq[M000] = 1.;
        // Order 1
        CMeq[M100] = 0.;
        CMeq[M010] = 0.;
        CMeq[M001] = 0.;
        // Order 2
        CMeq[M200] = D::cs2;
        CMeq[M020] = D::cs2;
        CMeq[M002] = D::cs2;
        CMeq[M110] = 0.;
        CMeq[M101] = 0.;
        CMeq[M011] = 0.;
        // Order 3
        CMeq[M210] = 0.;
        CMeq[M201] = 0.;
        CMeq[M021] = 0.;
        CMeq[M120] = 0.;
        CMeq[M102] = 0.;
        CMeq[M012] = 0.;
        CMeq[M111] = 0.;
        // Order 4
        CMeq[M220] = D::cs2 * D::cs2;
        CMeq[M202] = D::cs2 * D::cs2;
        CMeq[M022] = D::cs2 * D::cs2;
        CMeq[M211] = 0.;
        CMeq[M121] = 0.;
        CMeq[M112] = 0.;
        // Order 5
        CMeq[M221] = 0.;
        CMeq[M212] = 0.;
        CMeq[M122] = 0.;
        // Order 6
        CMeq[M222] = D::cs2 * D::cs2 * D::cs2;
    };

    // Equilibrium populations based on 27 moments can be computed using either RM, HM, CM, CHM or
    // Gauss-Hermite formalisms. All formulations are equivalent for the D3Q27 so we use the fastest
    // one (RMs)
    static void CMcomputeEquilibrium(
        T rho, Array<T, D::d> const &u, Array<T, D::q> const &CMeq, Array<T, D::q> &eq)
    {
        Array<T, D::q> RMeq;
        // Order 2
        RMeq[M200] = u[0] * u[0] + D::cs2;
        RMeq[M020] = u[1] * u[1] + D::cs2;
        RMeq[M002] = u[2] * u[2] + D::cs2;
        RMeq[M110] = u[0] * u[1];
        RMeq[M101] = u[0] * u[2];
        RMeq[M011] = u[1] * u[2];
        // Order 3
        RMeq[M210] = RMeq[M200] * u[1];
        RMeq[M201] = RMeq[M200] * u[2];
        RMeq[M021] = RMeq[M020] * u[2];
        RMeq[M120] = RMeq[M020] * u[0];
        RMeq[M102] = RMeq[M002] * u[0];
        RMeq[M012] = RMeq[M002] * u[1];
        RMeq[M111] = RMeq[M110] * u[2];
        // Order 4
        RMeq[M220] = RMeq[M200] * RMeq[M020];
        RMeq[M202] = RMeq[M200] * RMeq[M002];
        RMeq[M022] = RMeq[M020] * RMeq[M002];
        RMeq[M211] = RMeq[M200] * RMeq[M011];
        RMeq[M121] = RMeq[M020] * RMeq[M101];
        RMeq[M112] = RMeq[M002] * RMeq[M110];
        // Order 5
        RMeq[M221] = RMeq[M220] * u[2];
        RMeq[M212] = RMeq[M202] * u[1];
        RMeq[M122] = RMeq[M022] * u[0];
        // Order 6
        RMeq[M222] = RMeq[M220] * RMeq[M002];

        // Optimization based on symmetries between populations and their opposite counterpart
        eq[F000] = rho
                   * (1. - RMeq[M200] - RMeq[M020] - RMeq[M002] + RMeq[M220] + RMeq[M202]
                      + RMeq[M022] - RMeq[M222]);

        eq[FP00] = 0.5 * rho
                   * (u[0] + RMeq[M200] - RMeq[M120] - RMeq[M102] - RMeq[M220] - RMeq[M202]
                      + RMeq[M122] + RMeq[M222]);
        eq[FM00] = rho * (-u[0] + RMeq[M120] + RMeq[M102] - RMeq[M122]) + eq[FP00];

        eq[F0P0] = 0.5 * rho
                   * (u[1] + RMeq[M020] - RMeq[M210] - RMeq[M012] - RMeq[M220] - RMeq[M022]
                      + RMeq[M212] + RMeq[M222]);
        eq[F0M0] = rho * (-u[1] + RMeq[M210] + RMeq[M012] - RMeq[M212]) + eq[F0P0];

        eq[F00P] = 0.5 * rho
                   * (u[2] + RMeq[M002] - RMeq[M201] - RMeq[M021] - RMeq[M202] - RMeq[M022]
                      + RMeq[M221] + RMeq[M222]);
        eq[F00M] = rho * (-u[2] + RMeq[M201] + RMeq[M021] - RMeq[M221]) + eq[F00P];

        eq[FPP0] = 0.25 * rho
                   * (RMeq[M110] + RMeq[M210] + RMeq[M120] - RMeq[M112] + RMeq[M220] - RMeq[M212]
                      - RMeq[M122] - RMeq[M222]);
        eq[FMP0] = 0.5 * rho * (-RMeq[M110] - RMeq[M120] + RMeq[M112] + RMeq[M122]) + eq[FPP0];
        eq[FPM0] = 0.5 * rho * (-RMeq[M110] - RMeq[M210] + RMeq[M112] + RMeq[M212]) + eq[FPP0];
        eq[FMM0] = 0.5 * rho * (-RMeq[M210] - RMeq[M120] + RMeq[M212] + RMeq[M122]) + eq[FPP0];

        eq[FP0P] = 0.25 * rho
                   * (RMeq[M101] + RMeq[M201] + RMeq[M102] - RMeq[M121] + RMeq[M202] - RMeq[M221]
                      - RMeq[M122] - RMeq[M222]);
        eq[FM0P] = 0.5 * rho * (-RMeq[M101] - RMeq[M102] + RMeq[M121] + RMeq[M122]) + eq[FP0P];
        eq[FP0M] = 0.5 * rho * (-RMeq[M101] - RMeq[M201] + RMeq[M121] + RMeq[M221]) + eq[FP0P];
        eq[FM0M] = 0.5 * rho * (-RMeq[M201] - RMeq[M102] + RMeq[M221] + RMeq[M122]) + eq[FP0P];

        eq[F0PP] = 0.25 * rho
                   * (RMeq[M011] + RMeq[M021] + RMeq[M012] - RMeq[M211] + RMeq[M022] - RMeq[M221]
                      - RMeq[M212] - RMeq[M222]);
        eq[F0MP] = 0.5 * rho * (-RMeq[M011] - RMeq[M012] + RMeq[M211] + RMeq[M212]) + eq[F0PP];
        eq[F0PM] = 0.5 * rho * (-RMeq[M011] - RMeq[M021] + RMeq[M211] + RMeq[M221]) + eq[F0PP];
        eq[F0MM] = 0.5 * rho * (-RMeq[M021] - RMeq[M012] + RMeq[M221] + RMeq[M212]) + eq[F0PP];

        eq[FPPP] = 0.125 * rho
                   * (RMeq[M111] + RMeq[M211] + RMeq[M121] + RMeq[M112] + RMeq[M221] + RMeq[M212]
                      + RMeq[M122] + RMeq[M222]);
        eq[FMPP] = 0.25 * rho * (-RMeq[M111] - RMeq[M121] - RMeq[M112] - RMeq[M122]) + eq[FPPP];
        eq[FPMP] = 0.25 * rho * (-RMeq[M111] - RMeq[M211] - RMeq[M112] - RMeq[M212]) + eq[FPPP];
        eq[FPPM] = 0.25 * rho * (-RMeq[M111] - RMeq[M211] - RMeq[M121] - RMeq[M221]) + eq[FPPP];
        eq[FMMP] = 0.25 * rho * (-RMeq[M211] - RMeq[M121] - RMeq[M212] - RMeq[M122]) + eq[FPPP];
        eq[FMPM] = 0.25 * rho * (-RMeq[M211] - RMeq[M112] - RMeq[M221] - RMeq[M122]) + eq[FPPP];
        eq[FPMM] = 0.25 * rho * (-RMeq[M121] - RMeq[M112] - RMeq[M221] - RMeq[M212]) + eq[FPPP];
        eq[FMMM] = 0.25 * rho * (-RMeq[M111] - RMeq[M221] - RMeq[M212] - RMeq[M122]) + eq[FPPP];
    };

    static void CMcollide(
        Array<T, D::q> &cell, T rho, Array<T, D::d> const &u,
        Array<T, D::q> const &CM,    // Central moments
        Array<T, D::q> const &CMeq,  // Equilibrium moments (central)
        Array<T, D::numRelaxationTimes> const &omega)
    {
        T omega1 = omega[0];
        T omega2 = omega[1];
        T omega3 = omega[2];
        T omega4 = omega[3];
        T omega5 = omega[4];
        T omega6 = omega[5];
        T omega7 = omega[6];
        T omega8 = omega[7];

        T omegaBulk = omega[8];
        T omegaPlus = (omegaBulk + 2. * omega1) / 3.;  // Notation used by Fei
        T omegaMinus = (omegaBulk - omega1) / 3.;      // Notation used by Fei

        T ux2 = u[0] * u[0];
        T uy2 = u[1] * u[1];
        T uz2 = u[2] * u[2];
        T uxyz = u[0] * u[1] * u[2];

        // Post-collision moments.
        Array<T, D::q> CMcoll;
        Array<T, D::q> RMcoll;

        // Collision in the central moment space
        // Order 2 (non-diagonal collision so that we can easily modify the bulk viscosity)
        CMcoll[M200] = CM[M200] - omegaPlus * (CM[M200] - CMeq[M200])
                       - omegaMinus * (CM[M020] - CMeq[M020])
                       - omegaMinus * (CM[M002] - CMeq[M002]);
        CMcoll[M020] = CM[M020] - omegaMinus * (CM[M200] - CMeq[M200])
                       - omegaPlus * (CM[M020] - CMeq[M020]) - omegaMinus * (CM[M002] - CMeq[M002]);
        CMcoll[M002] = CM[M002] - omegaMinus * (CM[M200] - CMeq[M200])
                       - omegaMinus * (CM[M020] - CMeq[M020]) - omegaPlus * (CM[M002] - CMeq[M002]);

        CMcoll[M110] = (1. - omega2) * CM[M110] + omega2 * CMeq[M110];
        CMcoll[M101] = (1. - omega2) * CM[M101] + omega2 * CMeq[M101];
        CMcoll[M011] = (1. - omega2) * CM[M011] + omega2 * CMeq[M011];

        // Order 3
        CMcoll[M210] = (1. - omega3) * CM[M210] + omega3 * CMeq[M210];
        CMcoll[M201] = (1. - omega3) * CM[M201] + omega3 * CMeq[M201];
        CMcoll[M021] = (1. - omega3) * CM[M021] + omega3 * CMeq[M021];
        CMcoll[M120] = (1. - omega3) * CM[M120] + omega3 * CMeq[M120];
        CMcoll[M102] = (1. - omega3) * CM[M102] + omega3 * CMeq[M102];
        CMcoll[M012] = (1. - omega3) * CM[M012] + omega3 * CMeq[M012];

        CMcoll[M111] = (1. - omega4) * CM[M111] + omega4 * CMeq[M111];

        // Order 4
        CMcoll[M220] = (1. - omega5) * CM[M220] + omega5 * CMeq[M220];
        CMcoll[M202] = (1. - omega5) * CM[M202] + omega5 * CMeq[M202];
        CMcoll[M022] = (1. - omega5) * CM[M022] + omega5 * CMeq[M022];

        CMcoll[M211] = (1. - omega6) * CM[M211] + omega6 * CMeq[M211];
        CMcoll[M121] = (1. - omega6) * CM[M121] + omega6 * CMeq[M121];
        CMcoll[M112] = (1. - omega6) * CM[M112] + omega6 * CMeq[M112];

        // Order 5
        CMcoll[M221] = (1. - omega7) * CM[M221] + omega7 * CMeq[M221];
        CMcoll[M212] = (1. - omega7) * CM[M212] + omega7 * CMeq[M212];
        CMcoll[M122] = (1. - omega7) * CM[M122] + omega7 * CMeq[M122];

        // Order 6
        CMcoll[M222] = (1. - omega8) * CM[M222] + omega8 * CMeq[M222];

        // Come back to RMcoll using binomial formulas
        RMcoll[M200] = CMcoll[M200] + ux2;
        RMcoll[M020] = CMcoll[M020] + uy2;
        RMcoll[M002] = CMcoll[M002] + uz2;

        RMcoll[M110] = CMcoll[M110] + u[0] * u[1];
        RMcoll[M101] = CMcoll[M101] + u[0] * u[2];
        RMcoll[M011] = CMcoll[M011] + u[1] * u[2];

        RMcoll[M210] = CMcoll[M210] + u[1] * CMcoll[M200] + 2. * u[0] * CMcoll[M110] + ux2 * u[1];
        RMcoll[M201] = CMcoll[M201] + u[2] * CMcoll[M200] + 2. * u[0] * CMcoll[M101] + ux2 * u[2];
        RMcoll[M021] = CMcoll[M021] + u[2] * CMcoll[M020] + 2. * u[1] * CMcoll[M011] + uy2 * u[2];
        RMcoll[M120] = CMcoll[M120] + u[0] * CMcoll[M020] + 2. * u[1] * CMcoll[M110] + u[0] * uy2;
        RMcoll[M102] = CMcoll[M102] + u[0] * CMcoll[M002] + 2. * u[2] * CMcoll[M101] + u[0] * uz2;
        RMcoll[M012] = CMcoll[M012] + u[1] * CMcoll[M002] + 2. * u[2] * CMcoll[M011] + u[1] * uz2;

        RMcoll[M111] =
            CMcoll[M111] + u[2] * CMcoll[M110] + u[1] * CMcoll[M101] + u[0] * CMcoll[M011] + uxyz;

        RMcoll[M220] = CMcoll[M220] + 2. * u[1] * CMcoll[M210] + 2. * u[0] * CMcoll[M120]
                       + uy2 * CMcoll[M200] + ux2 * CMcoll[M020] + 4. * u[0] * u[1] * CMcoll[M110]
                       + ux2 * uy2;
        RMcoll[M202] = CMcoll[M202] + 2. * u[2] * CMcoll[M201] + 2. * u[0] * CMcoll[M102]
                       + uz2 * CMcoll[M200] + ux2 * CMcoll[M002] + 4. * u[0] * u[2] * CMcoll[M101]
                       + ux2 * uz2;
        RMcoll[M022] = CMcoll[M022] + 2. * u[2] * CMcoll[M021] + 2. * u[1] * CMcoll[M012]
                       + uz2 * CMcoll[M020] + uy2 * CMcoll[M002] + 4. * u[1] * u[2] * CMcoll[M011]
                       + uy2 * uz2;

        RMcoll[M211] = CMcoll[M211] + u[2] * CMcoll[M210] + u[1] * CMcoll[M201]
                       + 2. * u[0] * CMcoll[M111] + u[1] * u[2] * CMcoll[M200]
                       + 2. * u[0] * u[2] * CMcoll[M110] + 2. * u[0] * u[1] * CMcoll[M101]
                       + ux2 * CMcoll[M011] + ux2 * u[1] * u[2];
        RMcoll[M121] = CMcoll[M121] + u[2] * CMcoll[M120] + u[0] * CMcoll[M021]
                       + 2. * u[1] * CMcoll[M111] + u[0] * u[2] * CMcoll[M020]
                       + 2. * u[1] * u[2] * CMcoll[M110] + 2. * u[0] * u[1] * CMcoll[M011]
                       + uy2 * CMcoll[M101] + u[0] * uy2 * u[2];
        RMcoll[M112] = CMcoll[M112] + u[1] * CMcoll[M102] + u[0] * CMcoll[M012]
                       + 2. * u[2] * CMcoll[M111] + u[0] * u[1] * CMcoll[M002]
                       + 2. * u[1] * u[2] * CMcoll[M101] + 2. * u[0] * u[2] * CMcoll[M011]
                       + uz2 * CMcoll[M110] + u[0] * u[1] * uz2;

        RMcoll[M221] =
            CMcoll[M221] + u[2] * CMcoll[M220] + 2. * u[1] * CMcoll[M211] + 2. * u[0] * CMcoll[M121]
            + 2. * u[1] * u[2] * CMcoll[M210] + uy2 * CMcoll[M201] + ux2 * CMcoll[M021]
            + 2. * u[0] * u[2] * CMcoll[M120] + 4. * u[0] * u[1] * CMcoll[M111]
            + uy2 * u[2] * CMcoll[M200] + ux2 * u[2] * CMcoll[M020] + 4. * uxyz * CMcoll[M110]
            + 2. * u[0] * uy2 * CMcoll[M101] + 2. * ux2 * u[1] * CMcoll[M011] + ux2 * uy2 * u[2];
        RMcoll[M212] =
            CMcoll[M212] + u[1] * CMcoll[M202] + 2. * u[2] * CMcoll[M211] + 2. * u[0] * CMcoll[M112]
            + 2. * u[1] * u[2] * CMcoll[M201] + uz2 * CMcoll[M210] + ux2 * CMcoll[M012]
            + 2. * u[0] * u[1] * CMcoll[M102] + 4. * u[0] * u[2] * CMcoll[M111]
            + u[1] * uz2 * CMcoll[M200] + ux2 * u[1] * CMcoll[M002] + 4. * uxyz * CMcoll[M101]
            + 2. * u[0] * uz2 * CMcoll[M110] + 2. * ux2 * u[2] * CMcoll[M011] + ux2 * u[1] * uz2;
        RMcoll[M122] =
            CMcoll[M122] + u[0] * CMcoll[M022] + 2. * u[2] * CMcoll[M121] + 2. * u[1] * CMcoll[M112]
            + 2. * u[0] * u[2] * CMcoll[M021] + uz2 * CMcoll[M120] + uy2 * CMcoll[M102]
            + 2. * u[0] * u[1] * CMcoll[M012] + 4. * u[1] * u[2] * CMcoll[M111]
            + u[0] * uz2 * CMcoll[M020] + u[0] * uy2 * CMcoll[M002] + 4. * uxyz * CMcoll[M011]
            + 2. * u[1] * uz2 * CMcoll[M110] + 2. * uy2 * u[2] * CMcoll[M101] + u[0] * uy2 * uz2;

        RMcoll[M222] =
            CMcoll[M222] + 2. * u[2] * CMcoll[M221] + 2. * u[1] * CMcoll[M212]
            + 2. * u[0] * CMcoll[M122] + uz2 * CMcoll[M220] + uy2 * CMcoll[M202]
            + ux2 * CMcoll[M022] + 4. * u[1] * u[2] * CMcoll[M211] + 4. * u[0] * u[2] * CMcoll[M121]
            + 4. * u[0] * u[1] * CMcoll[M112] + 2. * u[1] * uz2 * CMcoll[M210]
            + 2. * uy2 * u[2] * CMcoll[M201] + 2. * ux2 * u[2] * CMcoll[M021]
            + 2. * u[0] * uz2 * CMcoll[M120] + 2. * u[0] * uy2 * CMcoll[M102]
            + 2. * ux2 * u[1] * CMcoll[M012] + 8. * uxyz * CMcoll[M111] + uy2 * uz2 * CMcoll[M200]
            + ux2 * uz2 * CMcoll[M020] + ux2 * uy2 * CMcoll[M002]
            + 4. * u[0] * u[1] * uz2 * CMcoll[M110] + 4. * u[0] * uy2 * u[2] * CMcoll[M101]
            + 4. * ux2 * u[1] * u[2] * CMcoll[M011] + ux2 * uy2 * uz2;

        // Compute post collision populations from RM
        // Optimization based on symmetries between populations and their opposite counterpart
        cell[F000] = rho
                     * (1. - RMcoll[M200] - RMcoll[M020] - RMcoll[M002] + RMcoll[M220]
                        + RMcoll[M202] + RMcoll[M022] - RMcoll[M222]);

        cell[FP00] = 0.5 * rho
                     * (u[0] + RMcoll[M200] - RMcoll[M120] - RMcoll[M102] - RMcoll[M220]
                        - RMcoll[M202] + RMcoll[M122] + RMcoll[M222]);
        cell[FM00] = rho * (-u[0] + RMcoll[M120] + RMcoll[M102] - RMcoll[M122]) + cell[FP00];

        cell[F0P0] = 0.5 * rho
                     * (u[1] + RMcoll[M020] - RMcoll[M210] - RMcoll[M012] - RMcoll[M220]
                        - RMcoll[M022] + RMcoll[M212] + RMcoll[M222]);
        cell[F0M0] = rho * (-u[1] + RMcoll[M210] + RMcoll[M012] - RMcoll[M212]) + cell[F0P0];

        cell[F00P] = 0.5 * rho
                     * (u[2] + RMcoll[M002] - RMcoll[M201] - RMcoll[M021] - RMcoll[M202]
                        - RMcoll[M022] + RMcoll[M221] + RMcoll[M222]);
        cell[F00M] = rho * (-u[2] + RMcoll[M201] + RMcoll[M021] - RMcoll[M221]) + cell[F00P];

        cell[FPP0] = 0.25 * rho
                     * (RMcoll[M110] + RMcoll[M210] + RMcoll[M120] - RMcoll[M112] + RMcoll[M220]
                        - RMcoll[M212] - RMcoll[M122] - RMcoll[M222]);
        cell[FMP0] =
            0.5 * rho * (-RMcoll[M110] - RMcoll[M120] + RMcoll[M112] + RMcoll[M122]) + cell[FPP0];
        cell[FPM0] =
            0.5 * rho * (-RMcoll[M110] - RMcoll[M210] + RMcoll[M112] + RMcoll[M212]) + cell[FPP0];
        cell[FMM0] =
            0.5 * rho * (-RMcoll[M210] - RMcoll[M120] + RMcoll[M212] + RMcoll[M122]) + cell[FPP0];

        cell[FP0P] = 0.25 * rho
                     * (RMcoll[M101] + RMcoll[M201] + RMcoll[M102] - RMcoll[M121] + RMcoll[M202]
                        - RMcoll[M221] - RMcoll[M122] - RMcoll[M222]);
        cell[FM0P] =
            0.5 * rho * (-RMcoll[M101] - RMcoll[M102] + RMcoll[M121] + RMcoll[M122]) + cell[FP0P];
        cell[FP0M] =
            0.5 * rho * (-RMcoll[M101] - RMcoll[M201] + RMcoll[M121] + RMcoll[M221]) + cell[FP0P];
        cell[FM0M] =
            0.5 * rho * (-RMcoll[M201] - RMcoll[M102] + RMcoll[M221] + RMcoll[M122]) + cell[FP0P];

        cell[F0PP] = 0.25 * rho
                     * (RMcoll[M011] + RMcoll[M021] + RMcoll[M012] - RMcoll[M211] + RMcoll[M022]
                        - RMcoll[M221] - RMcoll[M212] - RMcoll[M222]);
        cell[F0MP] =
            0.5 * rho * (-RMcoll[M011] - RMcoll[M012] + RMcoll[M211] + RMcoll[M212]) + cell[F0PP];
        cell[F0PM] =
            0.5 * rho * (-RMcoll[M011] - RMcoll[M021] + RMcoll[M211] + RMcoll[M221]) + cell[F0PP];
        cell[F0MM] =
            0.5 * rho * (-RMcoll[M021] - RMcoll[M012] + RMcoll[M221] + RMcoll[M212]) + cell[F0PP];

        cell[FPPP] = 0.125 * rho
                     * (RMcoll[M111] + RMcoll[M211] + RMcoll[M121] + RMcoll[M112] + RMcoll[M221]
                        + RMcoll[M212] + RMcoll[M122] + RMcoll[M222]);
        cell[FMPP] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M121] - RMcoll[M112] - RMcoll[M122]) + cell[FPPP];
        cell[FPMP] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M211] - RMcoll[M112] - RMcoll[M212]) + cell[FPPP];
        cell[FPPM] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M211] - RMcoll[M121] - RMcoll[M221]) + cell[FPPP];
        cell[FMMP] =
            0.25 * rho * (-RMcoll[M211] - RMcoll[M121] - RMcoll[M212] - RMcoll[M122]) + cell[FPPP];
        cell[FMPM] =
            0.25 * rho * (-RMcoll[M211] - RMcoll[M112] - RMcoll[M221] - RMcoll[M122]) + cell[FPPP];
        cell[FPMM] =
            0.25 * rho * (-RMcoll[M121] - RMcoll[M112] - RMcoll[M221] - RMcoll[M212]) + cell[FPPP];
        cell[FMMM] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M221] - RMcoll[M212] - RMcoll[M122]) + cell[FPPP];

        for (int i = 0; i < 27; ++i) {
            cell[i] -= D::SkordosFactor() * D::t[i];
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Central Hermite Moments Formalism (Equilibrium is computed through raw moments formalism) //
    ///////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * // General way to compute CHMs
     * static void CHMcomputeMoments(Array<T,D::q> const& cell, Array<T, D::q>& CHM, T& rho,
     * Array<T,D::d>& u) {
     *
     *     Array<T, D::q> f;
     *     for (int i = 0; i<27; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *         CHM[i] = 0.;
     *     }
     *
     *     rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[10] + f[11]
     * + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18] + f[19] + f[20] + f[21] + f[22] +
     * f[23] + f[24] + f[25] + f[26]; CHM[M000] = 1.0; T invRho = 1./rho; u[0] = invRho * (- f[1] -
     * f[4] - f[5] - f[6] - f[7] - f[10] - f[11] - f[12] - f[13] + f[14] + f[17] + f[18] + f[19] +
     * f[20] + f[23] + f[24] + f[25] + f[26]); u[1] = invRho * (- f[2] - f[4] + f[5] - f[8] - f[9] -
     * f[10] - f[11] + f[12] + f[13] + f[15] + f[17] - f[18] + f[21] + f[22] + f[23] + f[24] - f[25]
     * - f[26]); u[2] = invRho * (- f[3] - f[6] + f[7] - f[8] + f[9] - f[10] + f[11] - f[12] + f[13]
     * + f[16] + f[19] - f[20] + f[21] - f[22] + f[23] - f[24] + f[25] - f[26]); T cMux = 0.; T cMuy
     * = 0.; T cMuz = 0.; T Hxx = 0.; T Hyy = 0.; T Hzz = 0.;
     *
     *     for (int i = 0; i<27; ++i) {
     *
     *         cMux = D::c[i][0]- u[0];
     *         cMuy = D::c[i][1]- u[1];
     *         cMuz = D::c[i][2]- u[2];
     *
     *         Hxx = cMux * cMux - D::cs2;
     *         Hyy = cMuy * cMuy - D::cs2;
     *         Hzz = cMuz * cMuz - D::cs2;
     *
     *         // Order 1
     *         CHM[M100] += cMux * f[i];
     *         CHM[M010] += cMuy * f[i];
     *         CHM[M001] += cMuz * f[i];
     *         // Order 2
     *         CHM[M200] += Hxx * f[i];
     *         CHM[M020] += Hyy * f[i];
     *         CHM[M002] += Hzz * f[i];
     *         CHM[M110] += cMux * cMuy * f[i];
     *         CHM[M101] += cMux * cMuz * f[i];
     *         CHM[M011] += cMuy * cMuz * f[i];
     *         // Order 3
     *         CHM[M210] += Hxx * cMuy * f[i];
     *         CHM[M201] += Hxx * cMuz * f[i];
     *         CHM[M021] += Hyy * cMuz * f[i];
     *         CHM[M120] += cMux * Hyy * f[i];
     *         CHM[M102] += cMux * Hzz * f[i];
     *         CHM[M012] += cMuy * Hzz * f[i];
     *         CHM[M111] += cMux * cMuy * cMuz * f[i];
     *         // Order 4
     *         CHM[M220] += Hxx * Hyy * f[i];
     *         CHM[M202] += Hxx * Hzz * f[i];
     *         CHM[M022] += Hyy * Hzz * f[i];
     *         CHM[M211] += Hxx * cMuy * cMuz * f[i];
     *         CHM[M121] += cMux * Hyy * cMuz * f[i];
     *         CHM[M112] += cMux * cMuy * Hzz * f[i];
     *         // Order 5
     *         CHM[M221] += Hxx * Hyy * cMuz * f[i];
     *         CHM[M212] += Hxx * cMuy * Hzz * f[i];
     *         CHM[M122] += cMux * Hyy * Hzz * f[i];
     *         // Order 6
     *         CHM[M222] += Hxx * Hyy * Hzz * f[i];
     *     }
     *
     *     for (int i = 1; i<27; ++i) {
     *         CHM[i] *= invRho;
     *     }
     * };
     */

    // Optimized way to compute CHMs based on the general ordering of discrete velocities
    static void CHMcomputeMoments(
        Array<T, D::q> const &cell, Array<T, D::q> &CHM, T &rho, Array<T, D::d> &u)
    {
        Array<T, D::q> f;
        for (int i = 0; i < 27; ++i) {
            f[i] = cell[i] + D::SkordosFactor() * D::t[i];
            CHM[i] = 0.;
        }
        T A1 = f[FMMM] + f[FMMP] + f[FMPM] + f[FMPP];
        T A2 = f[FPPP] + f[FPPM] + f[FPMP] + f[FPMM];
        T A3 = f[FMMM] + f[FMMP] - f[FMPM] - f[FMPP];
        T A4 = f[FPPP] + f[FPPM] - f[FPMP] - f[FPMM];
        T A5 = f[FMMM] - f[FMMP] + f[FMPM] - f[FMPP];
        T A6 = f[FPPP] - f[FPPM] + f[FPMP] - f[FPMM];
        T A7 = f[FMMM] - f[FMMP] - f[FMPM] + f[FMPP];
        T A8 = f[FPPP] - f[FPPM] - f[FPMP] + f[FPMM];

        T X_M1 = f[FM00] + f[FMM0] + f[FMP0] + f[FM0M] + f[FM0P] + A1;
        T X_P1 = f[FP00] + f[FPP0] + f[FPM0] + f[FP0P] + f[FP0M] + A2;
        T X_0 =
            f[F0M0] + f[F00M] + f[F0MM] + f[F0MP] + f[F000] + f[F0P0] + f[F00P] + f[F0PP] + f[F0PM];

        T Y_M1 =
            f[F0M0] + f[FMM0] + f[F0MM] + f[F0MP] + f[FMMM] + f[FMMP] + f[FPM0] + f[FPMP] + f[FPMM];
        T Y_P1 =
            f[F0P0] + f[FPP0] + f[F0PP] + f[F0PM] + f[FPPP] + f[FPPM] + f[FMP0] + f[FMPM] + f[FMPP];

        T Z_M1 =
            f[F00M] + f[FM0M] + f[F0MM] + f[FMMM] + f[FMPM] + f[FP0M] + f[F0PM] + f[FPPM] + f[FPMM];
        T Z_P1 =
            f[F00P] + f[FP0P] + f[F0PP] + f[FPPP] + f[FPMP] + f[FM0P] + f[F0MP] + f[FMMP] + f[FMPP];

        rho = X_M1 + X_P1 + X_0;
        T invRho = 1. / rho;

        // Order 6
        CHM[M222] = invRho * (A1 + A2);
        // Order 5
        CHM[M221] = invRho * (-A5 + A6);
        CHM[M212] = invRho * (-A3 + A4);
        CHM[M122] = invRho * (-A1 + A2);
        // Order 4
        CHM[M220] = invRho * (f[FMM0] + f[FMP0] + f[FPP0] + f[FPM0]) + CHM[M222];
        CHM[M202] = invRho * (f[FM0M] + f[FM0P] + f[FP0P] + f[FP0M]) + CHM[M222];
        CHM[M022] = invRho * (f[F0MM] + f[F0MP] + f[F0PP] + f[F0PM]) + CHM[M222];
        CHM[M211] = invRho * (A7 + A8);
        CHM[M121] = invRho * (A5 + A6);
        CHM[M112] = invRho * (A3 + A4);
        // Order 3
        CHM[M210] = invRho * (-f[FMM0] + f[FMP0] + f[FPP0] - f[FPM0]) + CHM[M212];
        CHM[M201] = invRho * (-f[FM0M] + f[FM0P] + f[FP0P] - f[FP0M]) + CHM[M221];
        CHM[M021] = invRho * (-f[F0MM] + f[F0MP] + f[F0PP] - f[F0PM]) + CHM[M221];
        CHM[M120] = invRho * (-f[FMM0] - f[FMP0] + f[FPP0] + f[FPM0]) + CHM[M122];
        CHM[M102] = invRho * (-f[FM0M] - f[FM0P] + f[FP0P] + f[FP0M]) + CHM[M122];
        CHM[M012] = invRho * (-f[F0MM] - f[F0MP] + f[F0PP] + f[F0PM]) + CHM[M212];
        CHM[M111] = invRho * (-A7 + A8);
        // Order 2
        CHM[M200] = invRho * (X_P1 + X_M1);
        CHM[M020] = invRho * (Y_P1 + Y_M1);
        CHM[M002] = invRho * (Z_P1 + Z_M1);
        CHM[M110] = invRho * (f[FMM0] - f[FMP0] + f[FPP0] - f[FPM0]) + CHM[M112];
        CHM[M101] = invRho * (f[FM0M] - f[FM0P] + f[FP0P] - f[FP0M]) + CHM[M121];
        CHM[M011] = invRho * (f[F0MM] - f[F0MP] + f[F0PP] - f[F0PM]) + CHM[M211];
        // Order 1
        CHM[M100] = 0.;
        CHM[M010] = 0.;
        CHM[M001] = 0.;
        // Order 0
        CHM[M000] = 1.;

        // Compute CMs from RMs using binomial formulas
        u[0] = invRho * (X_P1 - X_M1);
        u[1] = invRho * (Y_P1 - Y_M1);
        u[2] = invRho * (Z_P1 - Z_M1);
        T ux2 = u[0] * u[0];
        T uy2 = u[1] * u[1];
        T uz2 = u[2] * u[2];
        T uxyz = u[0] * u[1] * u[2];

        CHM[M200] -= (ux2);
        CHM[M020] -= (uy2);
        CHM[M002] -= (uz2);

        CHM[M110] -= (u[0] * u[1]);
        CHM[M101] -= (u[0] * u[2]);
        CHM[M011] -= (u[1] * u[2]);

        CHM[M210] -= (u[1] * CHM[M200] + 2. * u[0] * CHM[M110] + ux2 * u[1]);
        CHM[M201] -= (u[2] * CHM[M200] + 2. * u[0] * CHM[M101] + ux2 * u[2]);
        CHM[M021] -= (u[2] * CHM[M020] + 2. * u[1] * CHM[M011] + uy2 * u[2]);
        CHM[M120] -= (u[0] * CHM[M020] + 2. * u[1] * CHM[M110] + u[0] * uy2);
        CHM[M102] -= (u[0] * CHM[M002] + 2. * u[2] * CHM[M101] + u[0] * uz2);
        CHM[M012] -= (u[1] * CHM[M002] + 2. * u[2] * CHM[M011] + u[1] * uz2);

        CHM[M111] -= (u[2] * CHM[M110] + u[1] * CHM[M101] + u[0] * CHM[M011] + uxyz);

        CHM[M220] -=
            (2. * u[1] * CHM[M210] + 2. * u[0] * CHM[M120] + uy2 * CHM[M200] + ux2 * CHM[M020]
             + 4. * u[0] * u[1] * CHM[M110] + ux2 * uy2);
        CHM[M202] -=
            (2. * u[2] * CHM[M201] + 2. * u[0] * CHM[M102] + uz2 * CHM[M200] + ux2 * CHM[M002]
             + 4. * u[0] * u[2] * CHM[M101] + ux2 * uz2);
        CHM[M022] -=
            (2. * u[2] * CHM[M021] + 2. * u[1] * CHM[M012] + uz2 * CHM[M020] + uy2 * CHM[M002]
             + 4. * u[1] * u[2] * CHM[M011] + uy2 * uz2);

        CHM[M211] -=
            (u[2] * CHM[M210] + u[1] * CHM[M201] + 2. * u[0] * CHM[M111] + u[1] * u[2] * CHM[M200]
             + 2. * u[0] * u[2] * CHM[M110] + 2. * u[0] * u[1] * CHM[M101] + ux2 * CHM[M011]
             + ux2 * u[1] * u[2]);
        CHM[M121] -=
            (u[2] * CHM[M120] + u[0] * CHM[M021] + 2. * u[1] * CHM[M111] + u[0] * u[2] * CHM[M020]
             + 2. * u[1] * u[2] * CHM[M110] + 2. * u[0] * u[1] * CHM[M011] + uy2 * CHM[M101]
             + u[0] * uy2 * u[2]);
        CHM[M112] -=
            (u[1] * CHM[M102] + u[0] * CHM[M012] + 2. * u[2] * CHM[M111] + u[0] * u[1] * CHM[M002]
             + 2. * u[1] * u[2] * CHM[M101] + 2. * u[0] * u[2] * CHM[M011] + uz2 * CHM[M110]
             + u[0] * u[1] * uz2);

        CHM[M221] -=
            (u[2] * CHM[M220] + 2. * u[1] * CHM[M211] + 2. * u[0] * CHM[M121]
             + 2. * u[1] * u[2] * CHM[M210] + uy2 * CHM[M201] + ux2 * CHM[M021]
             + 2. * u[0] * u[2] * CHM[M120] + 4. * u[0] * u[1] * CHM[M111] + uy2 * u[2] * CHM[M200]
             + ux2 * u[2] * CHM[M020] + 4. * uxyz * CHM[M110] + 2. * u[0] * uy2 * CHM[M101]
             + 2. * ux2 * u[1] * CHM[M011] + ux2 * uy2 * u[2]);
        CHM[M212] -=
            (u[1] * CHM[M202] + 2. * u[2] * CHM[M211] + 2. * u[0] * CHM[M112]
             + 2. * u[1] * u[2] * CHM[M201] + uz2 * CHM[M210] + ux2 * CHM[M012]
             + 2. * u[0] * u[1] * CHM[M102] + 4. * u[0] * u[2] * CHM[M111] + u[1] * uz2 * CHM[M200]
             + ux2 * u[1] * CHM[M002] + 4. * uxyz * CHM[M101] + 2. * u[0] * uz2 * CHM[M110]
             + 2. * ux2 * u[2] * CHM[M011] + ux2 * u[1] * uz2);
        CHM[M122] -=
            (u[0] * CHM[M022] + 2. * u[2] * CHM[M121] + 2. * u[1] * CHM[M112]
             + 2. * u[0] * u[2] * CHM[M021] + uz2 * CHM[M120] + uy2 * CHM[M102]
             + 2. * u[0] * u[1] * CHM[M012] + 4. * u[1] * u[2] * CHM[M111] + u[0] * uz2 * CHM[M020]
             + u[0] * uy2 * CHM[M002] + 4. * uxyz * CHM[M011] + 2. * u[1] * uz2 * CHM[M110]
             + 2. * uy2 * u[2] * CHM[M101] + u[0] * uy2 * uz2);

        CHM[M222] -=
            (2. * u[2] * CHM[M221] + 2. * u[1] * CHM[M212] + 2. * u[0] * CHM[M122] + uz2 * CHM[M220]
             + uy2 * CHM[M202] + ux2 * CHM[M022] + 4. * u[1] * u[2] * CHM[M211]
             + 4. * u[0] * u[2] * CHM[M121] + 4. * u[0] * u[1] * CHM[M112]
             + 2. * u[1] * uz2 * CHM[M210] + 2. * uy2 * u[2] * CHM[M201]
             + 2. * ux2 * u[2] * CHM[M021] + 2. * u[0] * uz2 * CHM[M120]
             + 2. * u[0] * uy2 * CHM[M102] + 2. * ux2 * u[1] * CHM[M012] + 8. * uxyz * CHM[M111]
             + uy2 * uz2 * CHM[M200] + ux2 * uz2 * CHM[M020] + ux2 * uy2 * CHM[M002]
             + 4. * u[0] * u[1] * uz2 * CHM[M110] + 4. * u[0] * uy2 * u[2] * CHM[M101]
             + 4. * ux2 * u[1] * u[2] * CHM[M011] + ux2 * uy2 * uz2);

        // Compute CHMs from CMs
        T cs4 = D::cs2 * D::cs2;

        CHM[M200] -= (D::cs2);
        CHM[M020] -= (D::cs2);
        CHM[M002] -= (D::cs2);

        CHM[M220] -= (D::cs2 * (CHM[M200] + CHM[M020]) + cs4);
        CHM[M202] -= (D::cs2 * (CHM[M200] + CHM[M002]) + cs4);
        CHM[M022] -= (D::cs2 * (CHM[M020] + CHM[M002]) + cs4);

        CHM[M211] -= (D::cs2 * CHM[M011]);
        CHM[M121] -= (D::cs2 * CHM[M101]);
        CHM[M112] -= (D::cs2 * CHM[M110]);

        CHM[M221] -= (D::cs2 * (CHM[M201] + CHM[M021]));
        CHM[M212] -= (D::cs2 * (CHM[M210] + CHM[M012]));
        CHM[M122] -= (D::cs2 * (CHM[M120] + CHM[M102]));

        CHM[M222] -=
            (D::cs2 * (CHM[M220] + CHM[M202] + CHM[M022])
             + cs4 * (CHM[M200] + CHM[M020] + CHM[M002]) + D::cs2 * cs4);
    };

    static void CHMcomputeEquilibriumMoments(Array<T, D::q> &CHMeq)
    {
        // Order 0
        CHMeq[M000] = 1.;
        // Order 1
        CHMeq[M100] = 0.;
        CHMeq[M010] = 0.;
        CHMeq[M001] = 0.;
        // Order 2
        CHMeq[M200] = 0.;
        CHMeq[M020] = 0.;
        CHMeq[M002] = 0.;
        CHMeq[M110] = 0.;
        CHMeq[M101] = 0.;
        CHMeq[M011] = 0.;
        // Order 3
        CHMeq[M210] = 0.;
        CHMeq[M201] = 0.;
        CHMeq[M021] = 0.;
        CHMeq[M120] = 0.;
        CHMeq[M102] = 0.;
        CHMeq[M012] = 0.;
        CHMeq[M111] = 0.;
        // Order 4
        CHMeq[M220] = 0.;
        CHMeq[M202] = 0.;
        CHMeq[M022] = 0.;
        CHMeq[M211] = 0.;
        CHMeq[M121] = 0.;
        CHMeq[M112] = 0.;
        // Order 5
        CHMeq[M221] = 0.;
        CHMeq[M212] = 0.;
        CHMeq[M122] = 0.;
        // Order 6
        CHMeq[M222] = 0.;
    };

    // Equilibrium populations based on 27 moments can be computed using either RM, HM, CM, CHM or
    // Gauss-Hermite formalisms. All formulations are equivalent for the D3Q27 so we use the fastest
    // one (RMs)
    static void CHMcomputeEquilibrium(
        T rho, Array<T, D::d> const &u, Array<T, D::q> const &CHMeq, Array<T, D::q> &eq)
    {
        Array<T, D::q> RMeq;
        // Order 2
        RMeq[M200] = u[0] * u[0] + D::cs2;
        RMeq[M020] = u[1] * u[1] + D::cs2;
        RMeq[M002] = u[2] * u[2] + D::cs2;
        RMeq[M110] = u[0] * u[1];
        RMeq[M101] = u[0] * u[2];
        RMeq[M011] = u[1] * u[2];
        // Order 3
        RMeq[M210] = RMeq[M200] * u[1];
        RMeq[M201] = RMeq[M200] * u[2];
        RMeq[M021] = RMeq[M020] * u[2];
        RMeq[M120] = RMeq[M020] * u[0];
        RMeq[M102] = RMeq[M002] * u[0];
        RMeq[M012] = RMeq[M002] * u[1];
        RMeq[M111] = RMeq[M110] * u[2];
        // Order 4
        RMeq[M220] = RMeq[M200] * RMeq[M020];
        RMeq[M202] = RMeq[M200] * RMeq[M002];
        RMeq[M022] = RMeq[M020] * RMeq[M002];
        RMeq[M211] = RMeq[M200] * RMeq[M011];
        RMeq[M121] = RMeq[M020] * RMeq[M101];
        RMeq[M112] = RMeq[M002] * RMeq[M110];
        // Order 5
        RMeq[M221] = RMeq[M220] * u[2];
        RMeq[M212] = RMeq[M202] * u[1];
        RMeq[M122] = RMeq[M022] * u[0];
        // Order 6
        RMeq[M222] = RMeq[M220] * RMeq[M002];

        // Optimization based on symmetries between populations and their opposite counterpart
        eq[F000] = rho
                   * (1. - RMeq[M200] - RMeq[M020] - RMeq[M002] + RMeq[M220] + RMeq[M202]
                      + RMeq[M022] - RMeq[M222]);

        eq[FP00] = 0.5 * rho
                   * (u[0] + RMeq[M200] - RMeq[M120] - RMeq[M102] - RMeq[M220] - RMeq[M202]
                      + RMeq[M122] + RMeq[M222]);
        eq[FM00] = rho * (-u[0] + RMeq[M120] + RMeq[M102] - RMeq[M122]) + eq[FP00];

        eq[F0P0] = 0.5 * rho
                   * (u[1] + RMeq[M020] - RMeq[M210] - RMeq[M012] - RMeq[M220] - RMeq[M022]
                      + RMeq[M212] + RMeq[M222]);
        eq[F0M0] = rho * (-u[1] + RMeq[M210] + RMeq[M012] - RMeq[M212]) + eq[F0P0];

        eq[F00P] = 0.5 * rho
                   * (u[2] + RMeq[M002] - RMeq[M201] - RMeq[M021] - RMeq[M202] - RMeq[M022]
                      + RMeq[M221] + RMeq[M222]);
        eq[F00M] = rho * (-u[2] + RMeq[M201] + RMeq[M021] - RMeq[M221]) + eq[F00P];

        eq[FPP0] = 0.25 * rho
                   * (RMeq[M110] + RMeq[M210] + RMeq[M120] - RMeq[M112] + RMeq[M220] - RMeq[M212]
                      - RMeq[M122] - RMeq[M222]);
        eq[FMP0] = 0.5 * rho * (-RMeq[M110] - RMeq[M120] + RMeq[M112] + RMeq[M122]) + eq[FPP0];
        eq[FPM0] = 0.5 * rho * (-RMeq[M110] - RMeq[M210] + RMeq[M112] + RMeq[M212]) + eq[FPP0];
        eq[FMM0] = 0.5 * rho * (-RMeq[M210] - RMeq[M120] + RMeq[M212] + RMeq[M122]) + eq[FPP0];

        eq[FP0P] = 0.25 * rho
                   * (RMeq[M101] + RMeq[M201] + RMeq[M102] - RMeq[M121] + RMeq[M202] - RMeq[M221]
                      - RMeq[M122] - RMeq[M222]);
        eq[FM0P] = 0.5 * rho * (-RMeq[M101] - RMeq[M102] + RMeq[M121] + RMeq[M122]) + eq[FP0P];
        eq[FP0M] = 0.5 * rho * (-RMeq[M101] - RMeq[M201] + RMeq[M121] + RMeq[M221]) + eq[FP0P];
        eq[FM0M] = 0.5 * rho * (-RMeq[M201] - RMeq[M102] + RMeq[M221] + RMeq[M122]) + eq[FP0P];

        eq[F0PP] = 0.25 * rho
                   * (RMeq[M011] + RMeq[M021] + RMeq[M012] - RMeq[M211] + RMeq[M022] - RMeq[M221]
                      - RMeq[M212] - RMeq[M222]);
        eq[F0MP] = 0.5 * rho * (-RMeq[M011] - RMeq[M012] + RMeq[M211] + RMeq[M212]) + eq[F0PP];
        eq[F0PM] = 0.5 * rho * (-RMeq[M011] - RMeq[M021] + RMeq[M211] + RMeq[M221]) + eq[F0PP];
        eq[F0MM] = 0.5 * rho * (-RMeq[M021] - RMeq[M012] + RMeq[M221] + RMeq[M212]) + eq[F0PP];

        eq[FPPP] = 0.125 * rho
                   * (RMeq[M111] + RMeq[M211] + RMeq[M121] + RMeq[M112] + RMeq[M221] + RMeq[M212]
                      + RMeq[M122] + RMeq[M222]);
        eq[FMPP] = 0.25 * rho * (-RMeq[M111] - RMeq[M121] - RMeq[M112] - RMeq[M122]) + eq[FPPP];
        eq[FPMP] = 0.25 * rho * (-RMeq[M111] - RMeq[M211] - RMeq[M112] - RMeq[M212]) + eq[FPPP];
        eq[FPPM] = 0.25 * rho * (-RMeq[M111] - RMeq[M211] - RMeq[M121] - RMeq[M221]) + eq[FPPP];
        eq[FMMP] = 0.25 * rho * (-RMeq[M211] - RMeq[M121] - RMeq[M212] - RMeq[M122]) + eq[FPPP];
        eq[FMPM] = 0.25 * rho * (-RMeq[M211] - RMeq[M112] - RMeq[M221] - RMeq[M122]) + eq[FPPP];
        eq[FPMM] = 0.25 * rho * (-RMeq[M121] - RMeq[M112] - RMeq[M221] - RMeq[M212]) + eq[FPPP];
        eq[FMMM] = 0.25 * rho * (-RMeq[M111] - RMeq[M221] - RMeq[M212] - RMeq[M122]) + eq[FPPP];
    };

    static void CHMcollide(
        Array<T, D::q> &cell, T rho, Array<T, D::d> const &u,
        Array<T, D::q> const &CHM,    // Central Hermite moments
        Array<T, D::q> const &CHMeq,  // Equilibrium moments (central Hermite)
        Array<T, D::numRelaxationTimes> const &omega)
    {
        T omega1 = omega[0];
        T omega2 = omega[1];
        T omega3 = omega[2];
        T omega4 = omega[3];
        T omega5 = omega[4];
        T omega6 = omega[5];
        T omega7 = omega[6];
        T omega8 = omega[7];

        T omegaBulk = omega[8];
        T omegaPlus = (omegaBulk + 2. * omega1) / 3.;  // Notation used by Fei
        T omegaMinus = (omegaBulk - omega1) / 3.;      // Notation used by Fei

        T ux2 = u[0] * u[0];
        T uy2 = u[1] * u[1];
        T uz2 = u[2] * u[2];
        T uxyz = u[0] * u[1] * u[2];
        T cs4 = D::cs2 * D::cs2;

        // Post-collision moments.
        Array<T, D::q> CHMcoll;
        Array<T, D::q> HMcoll;
        Array<T, D::q> RMcoll;

        // Collision in the central Hermite moment space
        // Order 2 (non-diagonal collision so that we can easily modify the bulk viscosity)
        CHMcoll[M200] = CHM[M200] - omegaPlus * (CHM[M200] - CHMeq[M200])
                        - omegaMinus * (CHM[M020] - CHMeq[M020])
                        - omegaMinus * (CHM[M002] - CHMeq[M002]);
        CHMcoll[M020] = CHM[M020] - omegaMinus * (CHM[M200] - CHMeq[M200])
                        - omegaPlus * (CHM[M020] - CHMeq[M020])
                        - omegaMinus * (CHM[M002] - CHMeq[M002]);
        CHMcoll[M002] = CHM[M002] - omegaMinus * (CHM[M200] - CHMeq[M200])
                        - omegaMinus * (CHM[M020] - CHMeq[M020])
                        - omegaPlus * (CHM[M002] - CHMeq[M002]);

        CHMcoll[M110] = (1. - omega2) * CHM[M110] + omega2 * CHMeq[M110];
        CHMcoll[M101] = (1. - omega2) * CHM[M101] + omega2 * CHMeq[M101];
        CHMcoll[M011] = (1. - omega2) * CHM[M011] + omega2 * CHMeq[M011];

        // Order 3
        CHMcoll[M210] = (1. - omega3) * CHM[M210] + omega3 * CHMeq[M210];
        CHMcoll[M201] = (1. - omega3) * CHM[M201] + omega3 * CHMeq[M201];
        CHMcoll[M021] = (1. - omega3) * CHM[M021] + omega3 * CHMeq[M021];
        CHMcoll[M120] = (1. - omega3) * CHM[M120] + omega3 * CHMeq[M120];
        CHMcoll[M102] = (1. - omega3) * CHM[M102] + omega3 * CHMeq[M102];
        CHMcoll[M012] = (1. - omega3) * CHM[M012] + omega3 * CHMeq[M012];

        CHMcoll[M111] = (1. - omega4) * CHM[M111] + omega4 * CHMeq[M111];

        // Order 4
        CHMcoll[M220] = (1. - omega5) * CHM[M220] + omega5 * CHMeq[M220];
        CHMcoll[M202] = (1. - omega5) * CHM[M202] + omega5 * CHMeq[M202];
        CHMcoll[M022] = (1. - omega5) * CHM[M022] + omega5 * CHMeq[M022];

        CHMcoll[M211] = (1. - omega6) * CHM[M211] + omega6 * CHMeq[M211];
        CHMcoll[M121] = (1. - omega6) * CHM[M121] + omega6 * CHMeq[M121];
        CHMcoll[M112] = (1. - omega6) * CHM[M112] + omega6 * CHMeq[M112];

        // Order 5
        CHMcoll[M221] = (1. - omega7) * CHM[M221] + omega7 * CHMeq[M221];
        CHMcoll[M212] = (1. - omega7) * CHM[M212] + omega7 * CHMeq[M212];
        CHMcoll[M122] = (1. - omega7) * CHM[M122] + omega7 * CHMeq[M122];

        // Order 6
        CHMcoll[M222] = (1. - omega8) * CHM[M222] + omega8 * CHMeq[M222];

        // Come back to HMcoll using relationships between CHMs and HMs
        HMcoll[M200] = CHMcoll[M200] + ux2;
        HMcoll[M020] = CHMcoll[M020] + uy2;
        HMcoll[M002] = CHMcoll[M002] + uz2;

        HMcoll[M110] = CHMcoll[M110] + u[0] * u[1];
        HMcoll[M101] = CHMcoll[M101] + u[0] * u[2];
        HMcoll[M011] = CHMcoll[M011] + u[1] * u[2];

        HMcoll[M210] =
            CHMcoll[M210] + u[1] * CHMcoll[M200] + 2. * u[0] * CHMcoll[M110] + ux2 * u[1];
        HMcoll[M201] =
            CHMcoll[M201] + u[2] * CHMcoll[M200] + 2. * u[0] * CHMcoll[M101] + ux2 * u[2];
        HMcoll[M021] =
            CHMcoll[M021] + u[2] * CHMcoll[M020] + 2. * u[1] * CHMcoll[M011] + uy2 * u[2];
        HMcoll[M120] =
            CHMcoll[M120] + u[0] * CHMcoll[M020] + 2. * u[1] * CHMcoll[M110] + u[0] * uy2;
        HMcoll[M102] =
            CHMcoll[M102] + u[0] * CHMcoll[M002] + 2. * u[2] * CHMcoll[M101] + u[0] * uz2;
        HMcoll[M012] =
            CHMcoll[M012] + u[1] * CHMcoll[M002] + 2. * u[2] * CHMcoll[M011] + u[1] * uz2;

        HMcoll[M111] = CHMcoll[M111] + u[2] * CHMcoll[M110] + u[1] * CHMcoll[M101]
                       + u[0] * CHMcoll[M011] + uxyz;

        HMcoll[M220] = CHMcoll[M220] + 2. * u[1] * CHMcoll[M210] + 2. * u[0] * CHMcoll[M120]
                       + uy2 * CHMcoll[M200] + ux2 * CHMcoll[M020]
                       + 4. * u[0] * u[1] * CHMcoll[M110] + ux2 * uy2;
        HMcoll[M202] = CHMcoll[M202] + 2. * u[2] * CHMcoll[M201] + 2. * u[0] * CHMcoll[M102]
                       + uz2 * CHMcoll[M200] + ux2 * CHMcoll[M002]
                       + 4. * u[0] * u[2] * CHMcoll[M101] + ux2 * uz2;
        HMcoll[M022] = CHMcoll[M022] + 2. * u[2] * CHMcoll[M021] + 2. * u[1] * CHMcoll[M012]
                       + uz2 * CHMcoll[M020] + uy2 * CHMcoll[M002]
                       + 4. * u[1] * u[2] * CHMcoll[M011] + uy2 * uz2;

        HMcoll[M211] = CHMcoll[M211] + u[2] * CHMcoll[M210] + u[1] * CHMcoll[M201]
                       + 2. * u[0] * CHMcoll[M111] + u[1] * u[2] * CHMcoll[M200]
                       + 2. * u[0] * u[2] * CHMcoll[M110] + 2. * u[0] * u[1] * CHMcoll[M101]
                       + ux2 * CHMcoll[M011] + ux2 * u[1] * u[2];
        HMcoll[M121] = CHMcoll[M121] + u[2] * CHMcoll[M120] + u[0] * CHMcoll[M021]
                       + 2. * u[1] * CHMcoll[M111] + u[0] * u[2] * CHMcoll[M020]
                       + 2. * u[1] * u[2] * CHMcoll[M110] + 2. * u[0] * u[1] * CHMcoll[M011]
                       + uy2 * CHMcoll[M101] + u[0] * uy2 * u[2];
        HMcoll[M112] = CHMcoll[M112] + u[1] * CHMcoll[M102] + u[0] * CHMcoll[M012]
                       + 2. * u[2] * CHMcoll[M111] + u[0] * u[1] * CHMcoll[M002]
                       + 2. * u[1] * u[2] * CHMcoll[M101] + 2. * u[0] * u[2] * CHMcoll[M011]
                       + uz2 * CHMcoll[M110] + u[0] * u[1] * uz2;

        HMcoll[M221] = CHMcoll[M221] + u[2] * CHMcoll[M220] + 2. * u[1] * CHMcoll[M211]
                       + 2. * u[0] * CHMcoll[M121] + 2. * u[1] * u[2] * CHMcoll[M210]
                       + uy2 * CHMcoll[M201] + ux2 * CHMcoll[M021]
                       + 2. * u[0] * u[2] * CHMcoll[M120] + 4. * u[0] * u[1] * CHMcoll[M111]
                       + uy2 * u[2] * CHMcoll[M200] + ux2 * u[2] * CHMcoll[M020]
                       + 4. * uxyz * CHMcoll[M110] + 2. * u[0] * uy2 * CHMcoll[M101]
                       + 2. * ux2 * u[1] * CHMcoll[M011] + ux2 * uy2 * u[2];
        HMcoll[M212] = CHMcoll[M212] + u[1] * CHMcoll[M202] + 2. * u[2] * CHMcoll[M211]
                       + 2. * u[0] * CHMcoll[M112] + 2. * u[1] * u[2] * CHMcoll[M201]
                       + uz2 * CHMcoll[M210] + ux2 * CHMcoll[M012]
                       + 2. * u[0] * u[1] * CHMcoll[M102] + 4. * u[0] * u[2] * CHMcoll[M111]
                       + u[1] * uz2 * CHMcoll[M200] + ux2 * u[1] * CHMcoll[M002]
                       + 4. * uxyz * CHMcoll[M101] + 2. * u[0] * uz2 * CHMcoll[M110]
                       + 2. * ux2 * u[2] * CHMcoll[M011] + ux2 * u[1] * uz2;
        HMcoll[M122] = CHMcoll[M122] + u[0] * CHMcoll[M022] + 2. * u[2] * CHMcoll[M121]
                       + 2. * u[1] * CHMcoll[M112] + 2. * u[0] * u[2] * CHMcoll[M021]
                       + uz2 * CHMcoll[M120] + uy2 * CHMcoll[M102]
                       + 2. * u[0] * u[1] * CHMcoll[M012] + 4. * u[1] * u[2] * CHMcoll[M111]
                       + u[0] * uz2 * CHMcoll[M020] + u[0] * uy2 * CHMcoll[M002]
                       + 4. * uxyz * CHMcoll[M011] + 2. * u[1] * uz2 * CHMcoll[M110]
                       + 2. * uy2 * u[2] * CHMcoll[M101] + u[0] * uy2 * uz2;

        HMcoll[M222] = CHMcoll[M222] + 2. * u[2] * CHMcoll[M221] + 2. * u[1] * CHMcoll[M212]
                       + 2. * u[0] * CHMcoll[M122] + uz2 * CHMcoll[M220] + uy2 * CHMcoll[M202]
                       + ux2 * CHMcoll[M022] + 4. * u[1] * u[2] * CHMcoll[M211]
                       + 4. * u[0] * u[2] * CHMcoll[M121] + 4. * u[0] * u[1] * CHMcoll[M112]
                       + 2. * u[1] * uz2 * CHMcoll[M210] + 2. * uy2 * u[2] * CHMcoll[M201]
                       + 2. * ux2 * u[2] * CHMcoll[M021] + 2. * u[0] * uz2 * CHMcoll[M120]
                       + 2. * u[0] * uy2 * CHMcoll[M102] + 2. * ux2 * u[1] * CHMcoll[M012]
                       + 8. * uxyz * CHMcoll[M111] + uy2 * uz2 * CHMcoll[M200]
                       + ux2 * uz2 * CHMcoll[M020] + ux2 * uy2 * CHMcoll[M002]
                       + 4. * u[0] * u[1] * uz2 * CHMcoll[M110]
                       + 4. * u[0] * uy2 * u[2] * CHMcoll[M101]
                       + 4. * ux2 * u[1] * u[2] * CHMcoll[M011] + ux2 * uy2 * uz2;

        // Come back to RMcoll using relationships between HMs and RMs
        RMcoll[M200] = HMcoll[M200] + D::cs2;
        RMcoll[M020] = HMcoll[M020] + D::cs2;
        RMcoll[M002] = HMcoll[M002] + D::cs2;

        RMcoll[M110] = HMcoll[M110];
        RMcoll[M101] = HMcoll[M101];
        RMcoll[M011] = HMcoll[M011];

        RMcoll[M210] = HMcoll[M210] + D::cs2 * u[1];
        RMcoll[M201] = HMcoll[M201] + D::cs2 * u[2];
        RMcoll[M021] = HMcoll[M021] + D::cs2 * u[2];
        RMcoll[M120] = HMcoll[M120] + D::cs2 * u[0];
        RMcoll[M102] = HMcoll[M102] + D::cs2 * u[0];
        RMcoll[M012] = HMcoll[M012] + D::cs2 * u[1];

        RMcoll[M111] = HMcoll[M111];

        RMcoll[M220] = HMcoll[M220] + D::cs2 * (HMcoll[M200] + HMcoll[M020]) + cs4;
        RMcoll[M202] = HMcoll[M202] + D::cs2 * (HMcoll[M200] + HMcoll[M002]) + cs4;
        RMcoll[M022] = HMcoll[M022] + D::cs2 * (HMcoll[M020] + HMcoll[M002]) + cs4;

        RMcoll[M211] = HMcoll[M211] + D::cs2 * HMcoll[M011];
        RMcoll[M121] = HMcoll[M121] + D::cs2 * HMcoll[M101];
        RMcoll[M112] = HMcoll[M112] + D::cs2 * HMcoll[M110];

        RMcoll[M221] = HMcoll[M221] + D::cs2 * (HMcoll[M201] + HMcoll[M021]) + cs4 * u[2];
        RMcoll[M212] = HMcoll[M212] + D::cs2 * (HMcoll[M210] + HMcoll[M012]) + cs4 * u[1];
        RMcoll[M122] = HMcoll[M122] + D::cs2 * (HMcoll[M120] + HMcoll[M102]) + cs4 * u[0];

        RMcoll[M222] = HMcoll[M222] + D::cs2 * (HMcoll[M220] + HMcoll[M202] + HMcoll[M022])
                       + cs4 * (HMcoll[M200] + HMcoll[M020] + HMcoll[M002]) + D::cs2 * cs4;

        // Compute post collision populations from RM
        // Optimization based on symmetries between populations and their opposite counterpart
        cell[F000] = rho
                     * (1. - RMcoll[M200] - RMcoll[M020] - RMcoll[M002] + RMcoll[M220]
                        + RMcoll[M202] + RMcoll[M022] - RMcoll[M222]);

        cell[FP00] = 0.5 * rho
                     * (u[0] + RMcoll[M200] - RMcoll[M120] - RMcoll[M102] - RMcoll[M220]
                        - RMcoll[M202] + RMcoll[M122] + RMcoll[M222]);
        cell[FM00] = rho * (-u[0] + RMcoll[M120] + RMcoll[M102] - RMcoll[M122]) + cell[FP00];

        cell[F0P0] = 0.5 * rho
                     * (u[1] + RMcoll[M020] - RMcoll[M210] - RMcoll[M012] - RMcoll[M220]
                        - RMcoll[M022] + RMcoll[M212] + RMcoll[M222]);
        cell[F0M0] = rho * (-u[1] + RMcoll[M210] + RMcoll[M012] - RMcoll[M212]) + cell[F0P0];

        cell[F00P] = 0.5 * rho
                     * (u[2] + RMcoll[M002] - RMcoll[M201] - RMcoll[M021] - RMcoll[M202]
                        - RMcoll[M022] + RMcoll[M221] + RMcoll[M222]);
        cell[F00M] = rho * (-u[2] + RMcoll[M201] + RMcoll[M021] - RMcoll[M221]) + cell[F00P];

        cell[FPP0] = 0.25 * rho
                     * (RMcoll[M110] + RMcoll[M210] + RMcoll[M120] - RMcoll[M112] + RMcoll[M220]
                        - RMcoll[M212] - RMcoll[M122] - RMcoll[M222]);
        cell[FMP0] =
            0.5 * rho * (-RMcoll[M110] - RMcoll[M120] + RMcoll[M112] + RMcoll[M122]) + cell[FPP0];
        cell[FPM0] =
            0.5 * rho * (-RMcoll[M110] - RMcoll[M210] + RMcoll[M112] + RMcoll[M212]) + cell[FPP0];
        cell[FMM0] =
            0.5 * rho * (-RMcoll[M210] - RMcoll[M120] + RMcoll[M212] + RMcoll[M122]) + cell[FPP0];

        cell[FP0P] = 0.25 * rho
                     * (RMcoll[M101] + RMcoll[M201] + RMcoll[M102] - RMcoll[M121] + RMcoll[M202]
                        - RMcoll[M221] - RMcoll[M122] - RMcoll[M222]);
        cell[FM0P] =
            0.5 * rho * (-RMcoll[M101] - RMcoll[M102] + RMcoll[M121] + RMcoll[M122]) + cell[FP0P];
        cell[FP0M] =
            0.5 * rho * (-RMcoll[M101] - RMcoll[M201] + RMcoll[M121] + RMcoll[M221]) + cell[FP0P];
        cell[FM0M] =
            0.5 * rho * (-RMcoll[M201] - RMcoll[M102] + RMcoll[M221] + RMcoll[M122]) + cell[FP0P];

        cell[F0PP] = 0.25 * rho
                     * (RMcoll[M011] + RMcoll[M021] + RMcoll[M012] - RMcoll[M211] + RMcoll[M022]
                        - RMcoll[M221] - RMcoll[M212] - RMcoll[M222]);
        cell[F0MP] =
            0.5 * rho * (-RMcoll[M011] - RMcoll[M012] + RMcoll[M211] + RMcoll[M212]) + cell[F0PP];
        cell[F0PM] =
            0.5 * rho * (-RMcoll[M011] - RMcoll[M021] + RMcoll[M211] + RMcoll[M221]) + cell[F0PP];
        cell[F0MM] =
            0.5 * rho * (-RMcoll[M021] - RMcoll[M012] + RMcoll[M221] + RMcoll[M212]) + cell[F0PP];

        cell[FPPP] = 0.125 * rho
                     * (RMcoll[M111] + RMcoll[M211] + RMcoll[M121] + RMcoll[M112] + RMcoll[M221]
                        + RMcoll[M212] + RMcoll[M122] + RMcoll[M222]);
        cell[FMPP] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M121] - RMcoll[M112] - RMcoll[M122]) + cell[FPPP];
        cell[FPMP] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M211] - RMcoll[M112] - RMcoll[M212]) + cell[FPPP];
        cell[FPPM] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M211] - RMcoll[M121] - RMcoll[M221]) + cell[FPPP];
        cell[FMMP] =
            0.25 * rho * (-RMcoll[M211] - RMcoll[M121] - RMcoll[M212] - RMcoll[M122]) + cell[FPPP];
        cell[FMPM] =
            0.25 * rho * (-RMcoll[M211] - RMcoll[M112] - RMcoll[M221] - RMcoll[M122]) + cell[FPPP];
        cell[FPMM] =
            0.25 * rho * (-RMcoll[M121] - RMcoll[M112] - RMcoll[M221] - RMcoll[M212]) + cell[FPPP];
        cell[FMMM] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M221] - RMcoll[M212] - RMcoll[M122]) + cell[FPPP];

        for (int i = 0; i < 27; ++i) {
            cell[i] -= D::SkordosFactor() * D::t[i];
        }
    };

    ////////////////////////////////////////////////////////////////////////////////
    // Cumulant Formalism (Equilibrium is computed through raw moments formalism) //
    ////////////////////////////////////////////////////////////////////////////////

    /**
     * // General way to compute Ks
     * static void KcomputeMoments(Array<T,D::q> const& cell, Array<T, D::q>& K, T& rho,
     * Array<T,D::d>& u) {
     *
     *     Array<T, D::q> f;
     *     Array<T,D::q> CM;
     *     for (int i = 0; i<27; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *         CM[i] = 0.;
     *     }
     *
     *     rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[10] + f[11]
     * + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18] + f[19] + f[20] + f[21] + f[22] +
     * f[23] + f[24] + f[25] + f[26]; CM[M000] = 1.0; T invRho = 1./rho; u[0] = invRho * (- f[1] -
     * f[4] - f[5] - f[6] - f[7] - f[10] - f[11] - f[12] - f[13] + f[14] + f[17] + f[18] + f[19] +
     * f[20] + f[23] + f[24] + f[25] + f[26]); u[1] = invRho * (- f[2] - f[4] + f[5] - f[8] - f[9] -
     * f[10] - f[11] + f[12] + f[13] + f[15] + f[17] - f[18] + f[21] + f[22] + f[23] + f[24] - f[25]
     * - f[26]); u[2] = invRho * (- f[3] - f[6] + f[7] - f[8] + f[9] - f[10] + f[11] - f[12] + f[13]
     * + f[16] + f[19] - f[20] + f[21] - f[22] + f[23] - f[24] + f[25] - f[26]); T cMux = 0.; T cMuy
     * = 0.; T cMuz = 0.;
     *
     *     // Computation of central moments in a first time
     *     for (int i = 0; i<27; ++i) {
     *         cMux = D::c[i][0]- u[0];
     *         cMuy = D::c[i][1]- u[1];
     *         cMuz = D::c[i][2]- u[2];
     *         // Order 1
     *         CM[M100] += cMux * f[i];
     *         CM[M010] += cMuy * f[i];
     *         CM[M001] += cMuz * f[i];
     *         // Order 2
     *         CM[M200] += cMux * cMux * f[i];
     *         CM[M020] += cMuy * cMuy * f[i];
     *         CM[M002] += cMuz * cMuz * f[i];
     *         CM[M110] += cMux * cMuy * f[i];
     *         CM[M101] += cMux * cMuz * f[i];
     *         CM[M011] += cMuy * cMuz * f[i];
     *         // Order 3
     *         CM[M210] += cMux * cMux * cMuy * f[i];
     *         CM[M201] += cMux * cMux * cMuz * f[i];
     *         CM[M021] += cMuy * cMuy * cMuz * f[i];
     *         CM[M120] += cMux * cMuy * cMuy * f[i];
     *         CM[M102] += cMux * cMuz * cMuz * f[i];
     *         CM[M012] += cMuy * cMuz * cMuz * f[i];
     *         CM[M111] += cMux * cMuy * cMuz * f[i];
     *         // Order 4
     *         CM[M220] += cMux * cMux * cMuy * cMuy * f[i];
     *         CM[M202] += cMux * cMux * cMuz * cMuz * f[i];
     *         CM[M022] += cMuy * cMuy * cMuz * cMuz * f[i];
     *         CM[M211] += cMux * cMux * cMuy * cMuz * f[i];
     *         CM[M121] += cMux * cMuy * cMuy * cMuz * f[i];
     *         CM[M112] += cMux * cMuy * cMuz * cMuz * f[i];
     *         // Order 5
     *         CM[M221] += cMux * cMux * cMuy * cMuy * cMuz * f[i];
     *         CM[M212] += cMux * cMux * cMuy * cMuz * cMuz * f[i];
     *         CM[M122] += cMux * cMuy * cMuy * cMuz * cMuz * f[i];
     *         // Order 6
     *         CM[M222] += cMux * cMux * cMuy * cMuy * cMuz * cMuz * f[i];
     *     }
     *
     *     // Normalize before the computation of cumulants !
     *     for (int i = 1; i<27; ++i) {
     *         CM[i] *= invRho;
     *     }
     *
     *     // Computation of cumulants through central moments
     *     K[M000] = CM[M000];
     *     K[M100] = CM[M100] + u[0];
     *     K[M010] = CM[M010] + u[1];
     *     K[M001] = CM[M001] + u[2];
     *     K[M200] = CM[M200];
     *     K[M020] = CM[M020];
     *     K[M002] = CM[M002];
     *     K[M110] = CM[M110];
     *     K[M101] = CM[M101];
     *     K[M011] = CM[M011];
     *     K[M210] = CM[M210];
     *     K[M201] = CM[M201];
     *     K[M021] = CM[M021];
     *     K[M120] = CM[M120];
     *     K[M102] = CM[M102];
     *     K[M012] = CM[M012];
     *     K[M111] = CM[M111];
     *
     *     K[M220] = CM[M220] - CM[M200]*CM[M020] - 2.*CM[M110]*CM[M110];
     *     K[M202] = CM[M202] - CM[M200]*CM[M002] - 2.*CM[M101]*CM[M101];
     *     K[M022] = CM[M022] - CM[M020]*CM[M002] - 2.*CM[M011]*CM[M011];
     *     K[M211] = CM[M211] - CM[M200]*CM[M011] - 2.*CM[M110]*CM[M101];
     *     K[M121] = CM[M121] - CM[M020]*CM[M101] - 2.*CM[M110]*CM[M011];
     *     K[M112] = CM[M112] - CM[M002]*CM[M110] - 2.*CM[M101]*CM[M011];
     *
     *     K[M221] = CM[M221] - CM[M201]*CM[M020] - CM[M021]*CM[M200] - 2.*CM[M210]*CM[M011]
     * - 2.*CM[M120]*CM[M101] - 4.*CM[M111]*CM[M110]; K[M212] = CM[M212] - CM[M210]*CM[M002] -
     * CM[M012]*CM[M200] - 2.*CM[M201]*CM[M011] - 2.*CM[M102]*CM[M110] - 4.*CM[M111]*CM[M101];
     *     K[M122] = CM[M122] - CM[M120]*CM[M002] - CM[M102]*CM[M020] - 2.*CM[M012]*CM[M110]
     * - 2.*CM[M021]*CM[M101] - 4.*CM[M111]*CM[M011];
     *
     *     K[M222] = CM[M222] - CM[M220]*CM[M002] - CM[M202]*CM[M020] - CM[M022]*CM[M200]
     * - 4.*(CM[M211]*CM[M011] + CM[M121]*CM[M101] + CM[M112]*CM[M110]) - 2.*(CM[M210]*CM[M012] +
     * CM[M201]*CM[M021] + CM[M120]*CM[M102]) - 4.*CM[M111]*CM[M111]
     * + 4.*(CM[M200]*CM[M011]*CM[M011] + CM[M020]*CM[M101]*CM[M101] + CM[M002]*CM[M110]*CM[M110])
     * + 16.*CM[M110]*CM[M101]*CM[M011] + 2.*CM[M200]*CM[M020]*CM[M002];
     * };
     */

    // Optimized way to compute Ks based on the general ordering of discrete velocities
    static void KcomputeMoments(
        Array<T, D::q> const &cell, Array<T, D::q> &K, T &rho, Array<T, D::d> &u)
    {
        Array<T, D::q> f;
        for (int i = 0; i < 27; ++i) {
            f[i] = cell[i] + D::SkordosFactor() * D::t[i];
            K[i] = 0.;
        }
        T A1 = f[FMMM] + f[FMMP] + f[FMPM] + f[FMPP];
        T A2 = f[FPPP] + f[FPPM] + f[FPMP] + f[FPMM];
        T A3 = f[FMMM] + f[FMMP] - f[FMPM] - f[FMPP];
        T A4 = f[FPPP] + f[FPPM] - f[FPMP] - f[FPMM];
        T A5 = f[FMMM] - f[FMMP] + f[FMPM] - f[FMPP];
        T A6 = f[FPPP] - f[FPPM] + f[FPMP] - f[FPMM];
        T A7 = f[FMMM] - f[FMMP] - f[FMPM] + f[FMPP];
        T A8 = f[FPPP] - f[FPPM] - f[FPMP] + f[FPMM];

        T X_M1 = f[FM00] + f[FMM0] + f[FMP0] + f[FM0M] + f[FM0P] + A1;
        T X_P1 = f[FP00] + f[FPP0] + f[FPM0] + f[FP0P] + f[FP0M] + A2;
        T X_0 =
            f[F0M0] + f[F00M] + f[F0MM] + f[F0MP] + f[F000] + f[F0P0] + f[F00P] + f[F0PP] + f[F0PM];

        T Y_M1 =
            f[F0M0] + f[FMM0] + f[F0MM] + f[F0MP] + f[FMMM] + f[FMMP] + f[FPM0] + f[FPMP] + f[FPMM];
        T Y_P1 =
            f[F0P0] + f[FPP0] + f[F0PP] + f[F0PM] + f[FPPP] + f[FPPM] + f[FMP0] + f[FMPM] + f[FMPP];

        T Z_M1 =
            f[F00M] + f[FM0M] + f[F0MM] + f[FMMM] + f[FMPM] + f[FP0M] + f[F0PM] + f[FPPM] + f[FPMM];
        T Z_P1 =
            f[F00P] + f[FP0P] + f[F0PP] + f[FPPP] + f[FPMP] + f[FM0P] + f[F0MP] + f[FMMP] + f[FMPP];

        rho = X_M1 + X_P1 + X_0;
        T invRho = 1. / rho;

        // Order 6
        K[M222] = invRho * (A1 + A2);
        // Order 5
        K[M221] = invRho * (-A5 + A6);
        K[M212] = invRho * (-A3 + A4);
        K[M122] = invRho * (-A1 + A2);
        // Order 4
        K[M220] = invRho * (f[FMM0] + f[FMP0] + f[FPP0] + f[FPM0]) + K[M222];
        K[M202] = invRho * (f[FM0M] + f[FM0P] + f[FP0P] + f[FP0M]) + K[M222];
        K[M022] = invRho * (f[F0MM] + f[F0MP] + f[F0PP] + f[F0PM]) + K[M222];
        K[M211] = invRho * (A7 + A8);
        K[M121] = invRho * (A5 + A6);
        K[M112] = invRho * (A3 + A4);
        // Order 3
        K[M210] = invRho * (-f[FMM0] + f[FMP0] + f[FPP0] - f[FPM0]) + K[M212];
        K[M201] = invRho * (-f[FM0M] + f[FM0P] + f[FP0P] - f[FP0M]) + K[M221];
        K[M021] = invRho * (-f[F0MM] + f[F0MP] + f[F0PP] - f[F0PM]) + K[M221];
        K[M120] = invRho * (-f[FMM0] - f[FMP0] + f[FPP0] + f[FPM0]) + K[M122];
        K[M102] = invRho * (-f[FM0M] - f[FM0P] + f[FP0P] + f[FP0M]) + K[M122];
        K[M012] = invRho * (-f[F0MM] - f[F0MP] + f[F0PP] + f[F0PM]) + K[M212];
        K[M111] = invRho * (-A7 + A8);
        // Order 2
        K[M200] = invRho * (X_P1 + X_M1);
        K[M020] = invRho * (Y_P1 + Y_M1);
        K[M002] = invRho * (Z_P1 + Z_M1);
        K[M110] = invRho * (f[FMM0] - f[FMP0] + f[FPP0] - f[FPM0]) + K[M112];
        K[M101] = invRho * (f[FM0M] - f[FM0P] + f[FP0P] - f[FP0M]) + K[M121];
        K[M011] = invRho * (f[F0MM] - f[F0MP] + f[F0PP] - f[F0PM]) + K[M211];
        // Order 1
        K[M100] = 0.;
        K[M010] = 0.;
        K[M001] = 0.;
        // Order 0
        K[M000] = 1.;

        // Compute CMs from RMs using binomial formulas
        u[0] = invRho * (X_P1 - X_M1);
        u[1] = invRho * (Y_P1 - Y_M1);
        u[2] = invRho * (Z_P1 - Z_M1);
        T ux2 = u[0] * u[0];
        T uy2 = u[1] * u[1];
        T uz2 = u[2] * u[2];
        T uxyz = u[0] * u[1] * u[2];

        K[M200] -= (ux2);
        K[M020] -= (uy2);
        K[M002] -= (uz2);

        K[M110] -= (u[0] * u[1]);
        K[M101] -= (u[0] * u[2]);
        K[M011] -= (u[1] * u[2]);

        K[M210] -= (u[1] * K[M200] + 2. * u[0] * K[M110] + ux2 * u[1]);
        K[M201] -= (u[2] * K[M200] + 2. * u[0] * K[M101] + ux2 * u[2]);
        K[M021] -= (u[2] * K[M020] + 2. * u[1] * K[M011] + uy2 * u[2]);
        K[M120] -= (u[0] * K[M020] + 2. * u[1] * K[M110] + u[0] * uy2);
        K[M102] -= (u[0] * K[M002] + 2. * u[2] * K[M101] + u[0] * uz2);
        K[M012] -= (u[1] * K[M002] + 2. * u[2] * K[M011] + u[1] * uz2);

        K[M111] -= (u[2] * K[M110] + u[1] * K[M101] + u[0] * K[M011] + uxyz);

        K[M220] -=
            (2. * u[1] * K[M210] + 2. * u[0] * K[M120] + uy2 * K[M200] + ux2 * K[M020]
             + 4. * u[0] * u[1] * K[M110] + ux2 * uy2);
        K[M202] -=
            (2. * u[2] * K[M201] + 2. * u[0] * K[M102] + uz2 * K[M200] + ux2 * K[M002]
             + 4. * u[0] * u[2] * K[M101] + ux2 * uz2);
        K[M022] -=
            (2. * u[2] * K[M021] + 2. * u[1] * K[M012] + uz2 * K[M020] + uy2 * K[M002]
             + 4. * u[1] * u[2] * K[M011] + uy2 * uz2);

        K[M211] -=
            (u[2] * K[M210] + u[1] * K[M201] + 2. * u[0] * K[M111] + u[1] * u[2] * K[M200]
             + 2. * u[0] * u[2] * K[M110] + 2. * u[0] * u[1] * K[M101] + ux2 * K[M011]
             + ux2 * u[1] * u[2]);
        K[M121] -=
            (u[2] * K[M120] + u[0] * K[M021] + 2. * u[1] * K[M111] + u[0] * u[2] * K[M020]
             + 2. * u[1] * u[2] * K[M110] + 2. * u[0] * u[1] * K[M011] + uy2 * K[M101]
             + u[0] * uy2 * u[2]);
        K[M112] -=
            (u[1] * K[M102] + u[0] * K[M012] + 2. * u[2] * K[M111] + u[0] * u[1] * K[M002]
             + 2. * u[1] * u[2] * K[M101] + 2. * u[0] * u[2] * K[M011] + uz2 * K[M110]
             + u[0] * u[1] * uz2);

        K[M221] -=
            (u[2] * K[M220] + 2. * u[1] * K[M211] + 2. * u[0] * K[M121] + 2. * u[1] * u[2] * K[M210]
             + uy2 * K[M201] + ux2 * K[M021] + 2. * u[0] * u[2] * K[M120]
             + 4. * u[0] * u[1] * K[M111] + uy2 * u[2] * K[M200] + ux2 * u[2] * K[M020]
             + 4. * uxyz * K[M110] + 2. * u[0] * uy2 * K[M101] + 2. * ux2 * u[1] * K[M011]
             + ux2 * uy2 * u[2]);
        K[M212] -=
            (u[1] * K[M202] + 2. * u[2] * K[M211] + 2. * u[0] * K[M112] + 2. * u[1] * u[2] * K[M201]
             + uz2 * K[M210] + ux2 * K[M012] + 2. * u[0] * u[1] * K[M102]
             + 4. * u[0] * u[2] * K[M111] + u[1] * uz2 * K[M200] + ux2 * u[1] * K[M002]
             + 4. * uxyz * K[M101] + 2. * u[0] * uz2 * K[M110] + 2. * ux2 * u[2] * K[M011]
             + ux2 * u[1] * uz2);
        K[M122] -=
            (u[0] * K[M022] + 2. * u[2] * K[M121] + 2. * u[1] * K[M112] + 2. * u[0] * u[2] * K[M021]
             + uz2 * K[M120] + uy2 * K[M102] + 2. * u[0] * u[1] * K[M012]
             + 4. * u[1] * u[2] * K[M111] + u[0] * uz2 * K[M020] + u[0] * uy2 * K[M002]
             + 4. * uxyz * K[M011] + 2. * u[1] * uz2 * K[M110] + 2. * uy2 * u[2] * K[M101]
             + u[0] * uy2 * uz2);

        K[M222] -=
            (2. * u[2] * K[M221] + 2. * u[1] * K[M212] + 2. * u[0] * K[M122] + uz2 * K[M220]
             + uy2 * K[M202] + ux2 * K[M022] + 4. * u[1] * u[2] * K[M211]
             + 4. * u[0] * u[2] * K[M121] + 4. * u[0] * u[1] * K[M112] + 2. * u[1] * uz2 * K[M210]
             + 2. * uy2 * u[2] * K[M201] + 2. * ux2 * u[2] * K[M021] + 2. * u[0] * uz2 * K[M120]
             + 2. * u[0] * uy2 * K[M102] + 2. * ux2 * u[1] * K[M012] + 8. * uxyz * K[M111]
             + uy2 * uz2 * K[M200] + ux2 * uz2 * K[M020] + ux2 * uy2 * K[M002]
             + 4. * u[0] * u[1] * uz2 * K[M110] + 4. * u[0] * uy2 * u[2] * K[M101]
             + 4. * ux2 * u[1] * u[2] * K[M011] + ux2 * uy2 * uz2);

        // Computation of Ks through non-linear transformations of CMs
        // Here I am only using one array for all moments so we need to
        // first compute higher-order Ks before lower-order ones!
        K[M222] -=
            (K[M220] * K[M002] + K[M202] * K[M020] + K[M022] * K[M200]
             + 4. * (K[M211] * K[M011] + K[M121] * K[M101] + K[M112] * K[M110])
             + 2. * (K[M210] * K[M012] + K[M201] * K[M021] + K[M120] * K[M102])
             + 4. * K[M111] * K[M111]
             - 4.
                   * (K[M200] * K[M011] * K[M011] + K[M020] * K[M101] * K[M101]
                      + K[M002] * K[M110] * K[M110])
             - 16. * K[M110] * K[M101] * K[M011] - 2. * K[M200] * K[M020] * K[M002]);

        K[M221] -=
            (K[M201] * K[M020] + K[M021] * K[M200] + 2. * K[M210] * K[M011] + 2. * K[M120] * K[M101]
             + 4. * K[M111] * K[M110]);
        K[M212] -=
            (K[M210] * K[M002] + K[M012] * K[M200] + 2. * K[M201] * K[M011] + 2. * K[M102] * K[M110]
             + 4. * K[M111] * K[M101]);
        K[M122] -=
            (K[M120] * K[M002] + K[M102] * K[M020] + 2. * K[M012] * K[M110] + 2. * K[M021] * K[M101]
             + 4. * K[M111] * K[M011]);

        K[M220] -= (K[M200] * K[M020] + 2. * K[M110] * K[M110]);
        K[M202] -= (K[M200] * K[M002] + 2. * K[M101] * K[M101]);
        K[M022] -= (K[M020] * K[M002] + 2. * K[M011] * K[M011]);
        K[M211] -= (K[M200] * K[M011] + 2. * K[M110] * K[M101]);
        K[M121] -= (K[M020] * K[M101] + 2. * K[M110] * K[M011]);
        K[M112] -= (K[M002] * K[M110] + 2. * K[M101] * K[M011]);
    };

    static void KcomputeEquilibriumMoments(Array<T, D::d> const &u, Array<T, D::q> &Keq)
    {
        // Order 0
        Keq[M000] = 1.;
        // Order 1
        Keq[M100] = u[0];
        Keq[M010] = u[1];
        Keq[M001] = u[2];
        // Order 2
        Keq[M200] = D::cs2;
        Keq[M020] = D::cs2;
        Keq[M002] = D::cs2;
        Keq[M110] = 0.;
        Keq[M101] = 0.;
        Keq[M011] = 0.;
        // Order 3
        Keq[M210] = 0.;
        Keq[M201] = 0.;
        Keq[M021] = 0.;
        Keq[M120] = 0.;
        Keq[M102] = 0.;
        Keq[M012] = 0.;
        Keq[M111] = 0.;
        // Order 4
        Keq[M220] = 0.;
        Keq[M202] = 0.;
        Keq[M022] = 0.;
        Keq[M211] = 0.;
        Keq[M121] = 0.;
        Keq[M112] = 0.;
        // Order 5
        Keq[M221] = 0.;
        Keq[M212] = 0.;
        Keq[M122] = 0.;
        // Order 6
        Keq[M222] = 0.;
    };

    // Equilibrium populations based on 27 moments can be computed using either RM, HM, CM, CHM or
    // Gauss-Hermite formalisms. All formulations are equivalent for the D3Q27 so we use the fastest
    // one (RMs)
    static void KcomputeEquilibrium(
        T rho, Array<T, D::d> const &u, Array<T, D::q> const &Keq, Array<T, D::q> &eq)
    {
        Array<T, D::q> RMeq;
        // Order 2
        RMeq[M200] = u[0] * u[0] + D::cs2;
        RMeq[M020] = u[1] * u[1] + D::cs2;
        RMeq[M002] = u[2] * u[2] + D::cs2;
        RMeq[M110] = u[0] * u[1];
        RMeq[M101] = u[0] * u[2];
        RMeq[M011] = u[1] * u[2];
        // Order 3
        RMeq[M210] = RMeq[M200] * u[1];
        RMeq[M201] = RMeq[M200] * u[2];
        RMeq[M021] = RMeq[M020] * u[2];
        RMeq[M120] = RMeq[M020] * u[0];
        RMeq[M102] = RMeq[M002] * u[0];
        RMeq[M012] = RMeq[M002] * u[1];
        RMeq[M111] = RMeq[M110] * u[2];
        // Order 4
        RMeq[M220] = RMeq[M200] * RMeq[M020];
        RMeq[M202] = RMeq[M200] * RMeq[M002];
        RMeq[M022] = RMeq[M020] * RMeq[M002];
        RMeq[M211] = RMeq[M200] * RMeq[M011];
        RMeq[M121] = RMeq[M020] * RMeq[M101];
        RMeq[M112] = RMeq[M002] * RMeq[M110];
        // Order 5
        RMeq[M221] = RMeq[M220] * u[2];
        RMeq[M212] = RMeq[M202] * u[1];
        RMeq[M122] = RMeq[M022] * u[0];
        // Order 6
        RMeq[M222] = RMeq[M220] * RMeq[M002];

        // Optimization based on symmetries between populations and their opposite counterpart
        eq[F000] = rho
                   * (1. - RMeq[M200] - RMeq[M020] - RMeq[M002] + RMeq[M220] + RMeq[M202]
                      + RMeq[M022] - RMeq[M222]);

        eq[FP00] = 0.5 * rho
                   * (u[0] + RMeq[M200] - RMeq[M120] - RMeq[M102] - RMeq[M220] - RMeq[M202]
                      + RMeq[M122] + RMeq[M222]);
        eq[FM00] = rho * (-u[0] + RMeq[M120] + RMeq[M102] - RMeq[M122]) + eq[FP00];

        eq[F0P0] = 0.5 * rho
                   * (u[1] + RMeq[M020] - RMeq[M210] - RMeq[M012] - RMeq[M220] - RMeq[M022]
                      + RMeq[M212] + RMeq[M222]);
        eq[F0M0] = rho * (-u[1] + RMeq[M210] + RMeq[M012] - RMeq[M212]) + eq[F0P0];

        eq[F00P] = 0.5 * rho
                   * (u[2] + RMeq[M002] - RMeq[M201] - RMeq[M021] - RMeq[M202] - RMeq[M022]
                      + RMeq[M221] + RMeq[M222]);
        eq[F00M] = rho * (-u[2] + RMeq[M201] + RMeq[M021] - RMeq[M221]) + eq[F00P];

        eq[FPP0] = 0.25 * rho
                   * (RMeq[M110] + RMeq[M210] + RMeq[M120] - RMeq[M112] + RMeq[M220] - RMeq[M212]
                      - RMeq[M122] - RMeq[M222]);
        eq[FMP0] = 0.5 * rho * (-RMeq[M110] - RMeq[M120] + RMeq[M112] + RMeq[M122]) + eq[FPP0];
        eq[FPM0] = 0.5 * rho * (-RMeq[M110] - RMeq[M210] + RMeq[M112] + RMeq[M212]) + eq[FPP0];
        eq[FMM0] = 0.5 * rho * (-RMeq[M210] - RMeq[M120] + RMeq[M212] + RMeq[M122]) + eq[FPP0];

        eq[FP0P] = 0.25 * rho
                   * (RMeq[M101] + RMeq[M201] + RMeq[M102] - RMeq[M121] + RMeq[M202] - RMeq[M221]
                      - RMeq[M122] - RMeq[M222]);
        eq[FM0P] = 0.5 * rho * (-RMeq[M101] - RMeq[M102] + RMeq[M121] + RMeq[M122]) + eq[FP0P];
        eq[FP0M] = 0.5 * rho * (-RMeq[M101] - RMeq[M201] + RMeq[M121] + RMeq[M221]) + eq[FP0P];
        eq[FM0M] = 0.5 * rho * (-RMeq[M201] - RMeq[M102] + RMeq[M221] + RMeq[M122]) + eq[FP0P];

        eq[F0PP] = 0.25 * rho
                   * (RMeq[M011] + RMeq[M021] + RMeq[M012] - RMeq[M211] + RMeq[M022] - RMeq[M221]
                      - RMeq[M212] - RMeq[M222]);
        eq[F0MP] = 0.5 * rho * (-RMeq[M011] - RMeq[M012] + RMeq[M211] + RMeq[M212]) + eq[F0PP];
        eq[F0PM] = 0.5 * rho * (-RMeq[M011] - RMeq[M021] + RMeq[M211] + RMeq[M221]) + eq[F0PP];
        eq[F0MM] = 0.5 * rho * (-RMeq[M021] - RMeq[M012] + RMeq[M221] + RMeq[M212]) + eq[F0PP];

        eq[FPPP] = 0.125 * rho
                   * (RMeq[M111] + RMeq[M211] + RMeq[M121] + RMeq[M112] + RMeq[M221] + RMeq[M212]
                      + RMeq[M122] + RMeq[M222]);
        eq[FMPP] = 0.25 * rho * (-RMeq[M111] - RMeq[M121] - RMeq[M112] - RMeq[M122]) + eq[FPPP];
        eq[FPMP] = 0.25 * rho * (-RMeq[M111] - RMeq[M211] - RMeq[M112] - RMeq[M212]) + eq[FPPP];
        eq[FPPM] = 0.25 * rho * (-RMeq[M111] - RMeq[M211] - RMeq[M121] - RMeq[M221]) + eq[FPPP];
        eq[FMMP] = 0.25 * rho * (-RMeq[M211] - RMeq[M121] - RMeq[M212] - RMeq[M122]) + eq[FPPP];
        eq[FMPM] = 0.25 * rho * (-RMeq[M211] - RMeq[M112] - RMeq[M221] - RMeq[M122]) + eq[FPPP];
        eq[FPMM] = 0.25 * rho * (-RMeq[M121] - RMeq[M112] - RMeq[M221] - RMeq[M212]) + eq[FPPP];
        eq[FMMM] = 0.25 * rho * (-RMeq[M111] - RMeq[M221] - RMeq[M212] - RMeq[M122]) + eq[FPPP];
    };

    static void Kcollide(
        Array<T, D::q> &cell, T rho, Array<T, D::d> const &u, Array<T, D::q> const &K,  // Cumulants
        Array<T, D::q> const &Keq,  // Equilibrium cumulants
        Array<T, D::numRelaxationTimes> const &omega)
    {
        T omega1 = omega[0];
        T omega2 = omega[1];
        T omega3 = omega[2];
        T omega4 = omega[3];
        T omega5 = omega[4];
        T omega6 = omega[5];
        T omega7 = omega[6];
        T omega8 = omega[7];

        T omegaBulk = omega[8];
        T omegaPlus = (omegaBulk + 2. * omega1) / 3.;  // Notation used by Fei
        T omegaMinus = (omegaBulk - omega1) / 3.;      // Notation used by Fei

        T ux2 = u[0] * u[0];
        T uy2 = u[1] * u[1];
        T uz2 = u[2] * u[2];
        T uxyz = u[0] * u[1] * u[2];

        // Post-collision moments.
        Array<T, D::q> Kcoll;
        Array<T, D::q> CMcoll;
        Array<T, D::q> RMcoll;

        // Collision in the cumulant space
        // Order 2 (non-diagonal collision so that we can easily modify the bulk viscosity)
        Kcoll[M200] = K[M200] - omegaPlus * (K[M200] - Keq[M200])
                      - omegaMinus * (K[M020] - Keq[M020]) - omegaMinus * (K[M002] - Keq[M002]);
        Kcoll[M020] = K[M020] - omegaMinus * (K[M200] - Keq[M200])
                      - omegaPlus * (K[M020] - Keq[M020]) - omegaMinus * (K[M002] - Keq[M002]);
        Kcoll[M002] = K[M002] - omegaMinus * (K[M200] - Keq[M200])
                      - omegaMinus * (K[M020] - Keq[M020]) - omegaPlus * (K[M002] - Keq[M002]);

        Kcoll[M110] = (1. - omega2) * K[M110] + omega2 * Keq[M110];
        Kcoll[M101] = (1. - omega2) * K[M101] + omega2 * Keq[M101];
        Kcoll[M011] = (1. - omega2) * K[M011] + omega2 * Keq[M011];

        // Order 3
        Kcoll[M210] = (1. - omega3) * K[M210] + omega3 * Keq[M210];
        Kcoll[M201] = (1. - omega3) * K[M201] + omega3 * Keq[M201];
        Kcoll[M021] = (1. - omega3) * K[M021] + omega3 * Keq[M021];
        Kcoll[M120] = (1. - omega3) * K[M120] + omega3 * Keq[M120];
        Kcoll[M102] = (1. - omega3) * K[M102] + omega3 * Keq[M102];
        Kcoll[M012] = (1. - omega3) * K[M012] + omega3 * Keq[M012];

        Kcoll[M111] = (1. - omega4) * K[M111] + omega4 * Keq[M111];

        // Order 4
        Kcoll[M220] = (1. - omega5) * K[M220] + omega5 * Keq[M220];
        Kcoll[M202] = (1. - omega5) * K[M202] + omega5 * Keq[M202];
        Kcoll[M022] = (1. - omega5) * K[M022] + omega5 * Keq[M022];

        Kcoll[M211] = (1. - omega6) * K[M211] + omega6 * Keq[M211];
        Kcoll[M121] = (1. - omega6) * K[M121] + omega6 * Keq[M121];
        Kcoll[M112] = (1. - omega6) * K[M112] + omega6 * Keq[M112];

        // Order 5
        Kcoll[M221] = (1. - omega7) * K[M221] + omega7 * Keq[M221];
        Kcoll[M212] = (1. - omega7) * K[M212] + omega7 * Keq[M212];
        Kcoll[M122] = (1. - omega7) * K[M122] + omega7 * Keq[M122];

        // Order 6
        Kcoll[M222] = (1. - omega8) * K[M222] + omega8 * Keq[M222];

        // Come back to CMcoll modifying fourth- and higher-order post-collision cumulants
        CMcoll[M200] = Kcoll[M200];
        CMcoll[M020] = Kcoll[M020];
        CMcoll[M002] = Kcoll[M002];
        CMcoll[M110] = Kcoll[M110];
        CMcoll[M101] = Kcoll[M101];
        CMcoll[M011] = Kcoll[M011];

        CMcoll[M210] = Kcoll[M210];
        CMcoll[M201] = Kcoll[M201];
        CMcoll[M021] = Kcoll[M021];
        CMcoll[M120] = Kcoll[M120];
        CMcoll[M102] = Kcoll[M102];
        CMcoll[M012] = Kcoll[M012];
        CMcoll[M111] = Kcoll[M111];

        CMcoll[M220] = Kcoll[M220] + Kcoll[M200] * Kcoll[M020] + 2. * Kcoll[M110] * Kcoll[M110];
        CMcoll[M202] = Kcoll[M202] + Kcoll[M200] * Kcoll[M002] + 2. * Kcoll[M101] * Kcoll[M101];
        CMcoll[M022] = Kcoll[M022] + Kcoll[M020] * Kcoll[M002] + 2. * Kcoll[M011] * Kcoll[M011];
        CMcoll[M211] = Kcoll[M211] + Kcoll[M200] * Kcoll[M011] + 2. * Kcoll[M110] * Kcoll[M101];
        CMcoll[M121] = Kcoll[M121] + Kcoll[M020] * Kcoll[M101] + 2. * Kcoll[M110] * Kcoll[M011];
        CMcoll[M112] = Kcoll[M112] + Kcoll[M002] * Kcoll[M110] + 2. * Kcoll[M101] * Kcoll[M011];

        CMcoll[M221] = Kcoll[M221] + Kcoll[M201] * Kcoll[M020] + Kcoll[M021] * Kcoll[M200]
                       + 2. * Kcoll[M210] * Kcoll[M011] + 2. * Kcoll[M120] * Kcoll[M101]
                       + 4. * Kcoll[M111] * Kcoll[M110];
        CMcoll[M212] = Kcoll[M212] + Kcoll[M210] * Kcoll[M002] + Kcoll[M012] * Kcoll[M200]
                       + 2. * Kcoll[M201] * Kcoll[M011] + 2. * Kcoll[M102] * Kcoll[M110]
                       + 4. * Kcoll[M111] * Kcoll[M101];
        CMcoll[M122] = Kcoll[M122] + Kcoll[M120] * Kcoll[M002] + Kcoll[M102] * Kcoll[M020]
                       + 2. * Kcoll[M012] * Kcoll[M110] + 2. * Kcoll[M021] * Kcoll[M101]
                       + 4. * Kcoll[M111] * Kcoll[M011];

        CMcoll[M222] = Kcoll[M222] + Kcoll[M220] * Kcoll[M002] + Kcoll[M202] * Kcoll[M020]
                       + Kcoll[M022] * Kcoll[M200]
                       + 4.
                             * (Kcoll[M211] * Kcoll[M011] + Kcoll[M121] * Kcoll[M101]
                                + Kcoll[M112] * Kcoll[M110])
                       + 2.
                             * (Kcoll[M210] * Kcoll[M012] + Kcoll[M201] * Kcoll[M021]
                                + Kcoll[M120] * Kcoll[M102])
                       + 4. * Kcoll[M111] * Kcoll[M111]
                       + 2.
                             * (Kcoll[M200] * Kcoll[M011] * Kcoll[M011]
                                + Kcoll[M020] * Kcoll[M101] * Kcoll[M101]
                                + Kcoll[M002] * Kcoll[M110] * Kcoll[M110])
                       + 8. * Kcoll[M110] * Kcoll[M101] * Kcoll[M011]
                       + Kcoll[M200] * Kcoll[M020] * Kcoll[M002];

        // Come back to RMcoll using binomial formulas
        RMcoll[M200] = CMcoll[M200] + ux2;
        RMcoll[M020] = CMcoll[M020] + uy2;
        RMcoll[M002] = CMcoll[M002] + uz2;

        RMcoll[M110] = CMcoll[M110] + u[0] * u[1];
        RMcoll[M101] = CMcoll[M101] + u[0] * u[2];
        RMcoll[M011] = CMcoll[M011] + u[1] * u[2];

        RMcoll[M210] = CMcoll[M210] + u[1] * CMcoll[M200] + 2. * u[0] * CMcoll[M110] + ux2 * u[1];
        RMcoll[M201] = CMcoll[M201] + u[2] * CMcoll[M200] + 2. * u[0] * CMcoll[M101] + ux2 * u[2];
        RMcoll[M021] = CMcoll[M021] + u[2] * CMcoll[M020] + 2. * u[1] * CMcoll[M011] + uy2 * u[2];
        RMcoll[M120] = CMcoll[M120] + u[0] * CMcoll[M020] + 2. * u[1] * CMcoll[M110] + u[0] * uy2;
        RMcoll[M102] = CMcoll[M102] + u[0] * CMcoll[M002] + 2. * u[2] * CMcoll[M101] + u[0] * uz2;
        RMcoll[M012] = CMcoll[M012] + u[1] * CMcoll[M002] + 2. * u[2] * CMcoll[M011] + u[1] * uz2;

        RMcoll[M111] =
            CMcoll[M111] + u[2] * CMcoll[M110] + u[1] * CMcoll[M101] + u[0] * CMcoll[M011] + uxyz;

        RMcoll[M220] = CMcoll[M220] + 2. * u[1] * CMcoll[M210] + 2. * u[0] * CMcoll[M120]
                       + uy2 * CMcoll[M200] + ux2 * CMcoll[M020] + 4. * u[0] * u[1] * CMcoll[M110]
                       + ux2 * uy2;
        RMcoll[M202] = CMcoll[M202] + 2. * u[2] * CMcoll[M201] + 2. * u[0] * CMcoll[M102]
                       + uz2 * CMcoll[M200] + ux2 * CMcoll[M002] + 4. * u[0] * u[2] * CMcoll[M101]
                       + ux2 * uz2;
        RMcoll[M022] = CMcoll[M022] + 2. * u[2] * CMcoll[M021] + 2. * u[1] * CMcoll[M012]
                       + uz2 * CMcoll[M020] + uy2 * CMcoll[M002] + 4. * u[1] * u[2] * CMcoll[M011]
                       + uy2 * uz2;

        RMcoll[M211] = CMcoll[M211] + u[2] * CMcoll[M210] + u[1] * CMcoll[M201]
                       + 2. * u[0] * CMcoll[M111] + u[1] * u[2] * CMcoll[M200]
                       + 2. * u[0] * u[2] * CMcoll[M110] + 2. * u[0] * u[1] * CMcoll[M101]
                       + ux2 * CMcoll[M011] + ux2 * u[1] * u[2];
        RMcoll[M121] = CMcoll[M121] + u[2] * CMcoll[M120] + u[0] * CMcoll[M021]
                       + 2. * u[1] * CMcoll[M111] + u[0] * u[2] * CMcoll[M020]
                       + 2. * u[1] * u[2] * CMcoll[M110] + 2. * u[0] * u[1] * CMcoll[M011]
                       + uy2 * CMcoll[M101] + u[0] * uy2 * u[2];
        RMcoll[M112] = CMcoll[M112] + u[1] * CMcoll[M102] + u[0] * CMcoll[M012]
                       + 2. * u[2] * CMcoll[M111] + u[0] * u[1] * CMcoll[M002]
                       + 2. * u[1] * u[2] * CMcoll[M101] + 2. * u[0] * u[2] * CMcoll[M011]
                       + uz2 * CMcoll[M110] + u[0] * u[1] * uz2;

        RMcoll[M221] =
            CMcoll[M221] + u[2] * CMcoll[M220] + 2. * u[1] * CMcoll[M211] + 2. * u[0] * CMcoll[M121]
            + 2. * u[1] * u[2] * CMcoll[M210] + uy2 * CMcoll[M201] + ux2 * CMcoll[M021]
            + 2. * u[0] * u[2] * CMcoll[M120] + 4. * u[0] * u[1] * CMcoll[M111]
            + uy2 * u[2] * CMcoll[M200] + ux2 * u[2] * CMcoll[M020] + 4. * uxyz * CMcoll[M110]
            + 2. * u[0] * uy2 * CMcoll[M101] + 2. * ux2 * u[1] * CMcoll[M011] + ux2 * uy2 * u[2];
        RMcoll[M212] =
            CMcoll[M212] + u[1] * CMcoll[M202] + 2. * u[2] * CMcoll[M211] + 2. * u[0] * CMcoll[M112]
            + 2. * u[1] * u[2] * CMcoll[M201] + uz2 * CMcoll[M210] + ux2 * CMcoll[M012]
            + 2. * u[0] * u[1] * CMcoll[M102] + 4. * u[0] * u[2] * CMcoll[M111]
            + u[1] * uz2 * CMcoll[M200] + ux2 * u[1] * CMcoll[M002] + 4. * uxyz * CMcoll[M101]
            + 2. * u[0] * uz2 * CMcoll[M110] + 2. * ux2 * u[2] * CMcoll[M011] + ux2 * u[1] * uz2;
        RMcoll[M122] =
            CMcoll[M122] + u[0] * CMcoll[M022] + 2. * u[2] * CMcoll[M121] + 2. * u[1] * CMcoll[M112]
            + 2. * u[0] * u[2] * CMcoll[M021] + uz2 * CMcoll[M120] + uy2 * CMcoll[M102]
            + 2. * u[0] * u[1] * CMcoll[M012] + 4. * u[1] * u[2] * CMcoll[M111]
            + u[0] * uz2 * CMcoll[M020] + u[0] * uy2 * CMcoll[M002] + 4. * uxyz * CMcoll[M011]
            + 2. * u[1] * uz2 * CMcoll[M110] + 2. * uy2 * u[2] * CMcoll[M101] + u[0] * uy2 * uz2;

        RMcoll[M222] =
            CMcoll[M222] + 2. * u[2] * CMcoll[M221] + 2. * u[1] * CMcoll[M212]
            + 2. * u[0] * CMcoll[M122] + uz2 * CMcoll[M220] + uy2 * CMcoll[M202]
            + ux2 * CMcoll[M022] + 4. * u[1] * u[2] * CMcoll[M211] + 4. * u[0] * u[2] * CMcoll[M121]
            + 4. * u[0] * u[1] * CMcoll[M112] + 2. * u[1] * uz2 * CMcoll[M210]
            + 2. * uy2 * u[2] * CMcoll[M201] + 2. * ux2 * u[2] * CMcoll[M021]
            + 2. * u[0] * uz2 * CMcoll[M120] + 2. * u[0] * uy2 * CMcoll[M102]
            + 2. * ux2 * u[1] * CMcoll[M012] + 8. * uxyz * CMcoll[M111] + uy2 * uz2 * CMcoll[M200]
            + ux2 * uz2 * CMcoll[M020] + ux2 * uy2 * CMcoll[M002]
            + 4. * u[0] * u[1] * uz2 * CMcoll[M110] + 4. * u[0] * uy2 * u[2] * CMcoll[M101]
            + 4. * ux2 * u[1] * u[2] * CMcoll[M011] + ux2 * uy2 * uz2;

        // Compute post collision populations from RM
        // Optimization based on symmetries between populations and their opposite counterpart
        cell[F000] = rho
                     * (1. - RMcoll[M200] - RMcoll[M020] - RMcoll[M002] + RMcoll[M220]
                        + RMcoll[M202] + RMcoll[M022] - RMcoll[M222]);

        cell[FP00] = 0.5 * rho
                     * (u[0] + RMcoll[M200] - RMcoll[M120] - RMcoll[M102] - RMcoll[M220]
                        - RMcoll[M202] + RMcoll[M122] + RMcoll[M222]);
        cell[FM00] = rho * (-u[0] + RMcoll[M120] + RMcoll[M102] - RMcoll[M122]) + cell[FP00];

        cell[F0P0] = 0.5 * rho
                     * (u[1] + RMcoll[M020] - RMcoll[M210] - RMcoll[M012] - RMcoll[M220]
                        - RMcoll[M022] + RMcoll[M212] + RMcoll[M222]);
        cell[F0M0] = rho * (-u[1] + RMcoll[M210] + RMcoll[M012] - RMcoll[M212]) + cell[F0P0];

        cell[F00P] = 0.5 * rho
                     * (u[2] + RMcoll[M002] - RMcoll[M201] - RMcoll[M021] - RMcoll[M202]
                        - RMcoll[M022] + RMcoll[M221] + RMcoll[M222]);
        cell[F00M] = rho * (-u[2] + RMcoll[M201] + RMcoll[M021] - RMcoll[M221]) + cell[F00P];

        cell[FPP0] = 0.25 * rho
                     * (RMcoll[M110] + RMcoll[M210] + RMcoll[M120] - RMcoll[M112] + RMcoll[M220]
                        - RMcoll[M212] - RMcoll[M122] - RMcoll[M222]);
        cell[FMP0] =
            0.5 * rho * (-RMcoll[M110] - RMcoll[M120] + RMcoll[M112] + RMcoll[M122]) + cell[FPP0];
        cell[FPM0] =
            0.5 * rho * (-RMcoll[M110] - RMcoll[M210] + RMcoll[M112] + RMcoll[M212]) + cell[FPP0];
        cell[FMM0] =
            0.5 * rho * (-RMcoll[M210] - RMcoll[M120] + RMcoll[M212] + RMcoll[M122]) + cell[FPP0];

        cell[FP0P] = 0.25 * rho
                     * (RMcoll[M101] + RMcoll[M201] + RMcoll[M102] - RMcoll[M121] + RMcoll[M202]
                        - RMcoll[M221] - RMcoll[M122] - RMcoll[M222]);
        cell[FM0P] =
            0.5 * rho * (-RMcoll[M101] - RMcoll[M102] + RMcoll[M121] + RMcoll[M122]) + cell[FP0P];
        cell[FP0M] =
            0.5 * rho * (-RMcoll[M101] - RMcoll[M201] + RMcoll[M121] + RMcoll[M221]) + cell[FP0P];
        cell[FM0M] =
            0.5 * rho * (-RMcoll[M201] - RMcoll[M102] + RMcoll[M221] + RMcoll[M122]) + cell[FP0P];

        cell[F0PP] = 0.25 * rho
                     * (RMcoll[M011] + RMcoll[M021] + RMcoll[M012] - RMcoll[M211] + RMcoll[M022]
                        - RMcoll[M221] - RMcoll[M212] - RMcoll[M222]);
        cell[F0MP] =
            0.5 * rho * (-RMcoll[M011] - RMcoll[M012] + RMcoll[M211] + RMcoll[M212]) + cell[F0PP];
        cell[F0PM] =
            0.5 * rho * (-RMcoll[M011] - RMcoll[M021] + RMcoll[M211] + RMcoll[M221]) + cell[F0PP];
        cell[F0MM] =
            0.5 * rho * (-RMcoll[M021] - RMcoll[M012] + RMcoll[M221] + RMcoll[M212]) + cell[F0PP];

        cell[FPPP] = 0.125 * rho
                     * (RMcoll[M111] + RMcoll[M211] + RMcoll[M121] + RMcoll[M112] + RMcoll[M221]
                        + RMcoll[M212] + RMcoll[M122] + RMcoll[M222]);
        cell[FMPP] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M121] - RMcoll[M112] - RMcoll[M122]) + cell[FPPP];
        cell[FPMP] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M211] - RMcoll[M112] - RMcoll[M212]) + cell[FPPP];
        cell[FPPM] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M211] - RMcoll[M121] - RMcoll[M221]) + cell[FPPP];
        cell[FMMP] =
            0.25 * rho * (-RMcoll[M211] - RMcoll[M121] - RMcoll[M212] - RMcoll[M122]) + cell[FPPP];
        cell[FMPM] =
            0.25 * rho * (-RMcoll[M211] - RMcoll[M112] - RMcoll[M221] - RMcoll[M122]) + cell[FPPP];
        cell[FPMM] =
            0.25 * rho * (-RMcoll[M121] - RMcoll[M112] - RMcoll[M221] - RMcoll[M212]) + cell[FPPP];
        cell[FMMM] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M221] - RMcoll[M212] - RMcoll[M122]) + cell[FPPP];

        for (int i = 0; i < 27; ++i) {
            cell[i] -= D::SkordosFactor() * D::t[i];
        }
    };

    /////////////////////////////////////////////////////////////////////////////////
    // Gauss-Hermite Formalism (Equilibrium is computed through the GH  formalism) //
    /////////////////////////////////////////////////////////////////////////////////

    /**
     * // General way to compute GHs
     * static void GHcomputeMoments(Array<T,D::q> const& cell, Array<T, D::q>& GH, T& rho) {
     *
     *     Array<T, D::q> f;
     *     for (int i = 0; i<27; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *         GH[i] = 0.;
     *     }
     *
     *     T Hxx = 0.;
     *     T Hyy = 0.;
     *     T Hzz = 0.;
     *
     *     for (int i = 0; i<27; ++i) {
     *
     *         Hxx = D::c[i][0] * D::c[i][0] - D::cs2;
     *         Hyy = D::c[i][1] * D::c[i][1] - D::cs2;
     *         Hzz = D::c[i][2] * D::c[i][2] - D::cs2;
     *
     *         // Order 0
     *         GH[M000] += f[i];
     *
     *         // Order 1
     *         GH[M100] += D::c[i][0] * f[i];
     *         GH[M010] += D::c[i][1] * f[i];
     *         GH[M001] += D::c[i][2] * f[i];
     *
     *         // Order 2
     *         GH[M200] += Hxx * f[i];
     *         GH[M020] += Hyy * f[i];
     *         GH[M002] += Hzz * f[i];
     *         GH[M110] += D::c[i][0] * D::c[i][1] * f[i];
     *         GH[M101] += D::c[i][0] * D::c[i][2] * f[i];
     *         GH[M011] += D::c[i][1] * D::c[i][2] * f[i];
     *
     *         // Order 3
     *         GH[M210] += Hxx * D::c[i][1] * f[i];
     *         GH[M201] += Hxx * D::c[i][2] * f[i];
     *         GH[M021] += Hyy * D::c[i][2] * f[i];
     *         GH[M120] += D::c[i][0] * Hyy * f[i];
     *         GH[M102] += D::c[i][0] * Hzz * f[i];
     *         GH[M012] += D::c[i][1] * Hzz * f[i];
     *         GH[M111] += D::c[i][0] * D::c[i][1] * D::c[i][2] * f[i];
     *
     *         // Order 4
     *         GH[M220] += Hxx * Hyy * f[i];
     *         GH[M202] += Hxx * Hzz * f[i];
     *         GH[M022] += Hyy * Hzz * f[i];
     *         GH[M211] += Hxx * D::c[i][1] * D::c[i][2] * f[i];
     *         GH[M121] += D::c[i][0] * Hyy * D::c[i][2] * f[i];
     *         GH[M112] += D::c[i][0] * D::c[i][1] * Hzz * f[i];
     *
     *         // Order 5
     *         GH[M221] += Hxx * Hyy * D::c[i][2] * f[i];
     *         GH[M212] += Hxx * D::c[i][1] * Hzz * f[i];
     *         GH[M122] += D::c[i][0] * Hyy * Hzz * f[i];
     *
     *         // Order 6
     *         GH[M222] += Hxx * Hyy * Hzz * f[i];
     *     }
     *
     *     rho = GH[M000];
     *     T invRho = 1. / rho;
     *     for (int i = 0; i<27; ++i) {
     *         GH[i] *= invRho;
     *     }
     * };
     */

    // Optimized way to compute GHs based on the general ordering of discrete velocities
    static void GHcomputeMoments(Array<T, D::q> const &cell, Array<T, D::q> &GH, T &rho)
    {
        Array<T, D::q> f;
        for (int i = 0; i < 27; ++i) {
            f[i] = cell[i] + D::SkordosFactor() * D::t[i];
        }

        T A1 = f[FMMM] + f[FMMP] + f[FMPM] + f[FMPP];
        T A2 = f[FPPP] + f[FPPM] + f[FPMP] + f[FPMM];
        T A3 = f[FMMM] + f[FMMP] - f[FMPM] - f[FMPP];
        T A4 = f[FPPP] + f[FPPM] - f[FPMP] - f[FPMM];
        T A5 = f[FMMM] - f[FMMP] + f[FMPM] - f[FMPP];
        T A6 = f[FPPP] - f[FPPM] + f[FPMP] - f[FPMM];
        T A7 = f[FMMM] - f[FMMP] - f[FMPM] + f[FMPP];
        T A8 = f[FPPP] - f[FPPM] - f[FPMP] + f[FPMM];

        T X_M1 = f[FM00] + f[FMM0] + f[FMP0] + f[FM0M] + f[FM0P] + A1;
        T X_P1 = f[FP00] + f[FPP0] + f[FPM0] + f[FP0P] + f[FP0M] + A2;
        T X_0 =
            f[F0M0] + f[F00M] + f[F0MM] + f[F0MP] + f[F000] + f[F0P0] + f[F00P] + f[F0PP] + f[F0PM];

        T Y_M1 =
            f[F0M0] + f[FMM0] + f[F0MM] + f[F0MP] + f[FMMM] + f[FMMP] + f[FPM0] + f[FPMP] + f[FPMM];
        T Y_P1 =
            f[F0P0] + f[FPP0] + f[F0PP] + f[F0PM] + f[FPPP] + f[FPPM] + f[FMP0] + f[FMPM] + f[FMPP];

        T Z_M1 =
            f[F00M] + f[FM0M] + f[F0MM] + f[FMMM] + f[FMPM] + f[FP0M] + f[F0PM] + f[FPPM] + f[FPMM];
        T Z_P1 =
            f[F00P] + f[FP0P] + f[F0PP] + f[FPPP] + f[FPMP] + f[FM0P] + f[F0MP] + f[FMMP] + f[FMPP];

        rho = X_M1 + X_P1 + X_0;
        T invRho = 1. / rho;

        // Order 6
        GH[M222] = invRho * (A1 + A2);
        // Order 5
        GH[M221] = invRho * (-A5 + A6);
        GH[M212] = invRho * (-A3 + A4);
        GH[M122] = invRho * (-A1 + A2);
        // Order 4
        GH[M220] = invRho * (f[FMM0] + f[FMP0] + f[FPP0] + f[FPM0]) + GH[M222];
        GH[M202] = invRho * (f[FM0M] + f[FM0P] + f[FP0P] + f[FP0M]) + GH[M222];
        GH[M022] = invRho * (f[F0MM] + f[F0MP] + f[F0PP] + f[F0PM]) + GH[M222];
        GH[M211] = invRho * (A7 + A8);
        GH[M121] = invRho * (A5 + A6);
        GH[M112] = invRho * (A3 + A4);
        // Order 3
        GH[M210] = invRho * (-f[FMM0] + f[FMP0] + f[FPP0] - f[FPM0]) + GH[M212];
        GH[M201] = invRho * (-f[FM0M] + f[FM0P] + f[FP0P] - f[FP0M]) + GH[M221];
        GH[M021] = invRho * (-f[F0MM] + f[F0MP] + f[F0PP] - f[F0PM]) + GH[M221];
        GH[M120] = invRho * (-f[FMM0] - f[FMP0] + f[FPP0] + f[FPM0]) + GH[M122];
        GH[M102] = invRho * (-f[FM0M] - f[FM0P] + f[FP0P] + f[FP0M]) + GH[M122];
        GH[M012] = invRho * (-f[F0MM] - f[F0MP] + f[F0PP] + f[F0PM]) + GH[M212];
        GH[M111] = invRho * (-A7 + A8);
        // Order 2
        GH[M200] = invRho * (X_P1 + X_M1);
        GH[M020] = invRho * (Y_P1 + Y_M1);
        GH[M002] = invRho * (Z_P1 + Z_M1);
        GH[M110] = invRho * (f[FMM0] - f[FMP0] + f[FPP0] - f[FPM0]) + GH[M112];
        GH[M101] = invRho * (f[FM0M] - f[FM0P] + f[FP0P] - f[FP0M]) + GH[M121];
        GH[M011] = invRho * (f[F0MM] - f[F0MP] + f[F0PP] - f[F0PM]) + GH[M211];
        // Order 1
        GH[M100] = invRho * (X_P1 - X_M1);
        GH[M010] = invRho * (Y_P1 - Y_M1);
        GH[M001] = invRho * (Z_P1 - Z_M1);
        // Order 0
        GH[M000] = 1.;

        // We come back to Hermite moments
        T cs4 = D::cs2 * D::cs2;
        GH[M200] -= D::cs2;
        GH[M020] -= D::cs2;
        GH[M002] -= D::cs2;

        GH[M210] -= D::cs2 * GH[M010];
        GH[M201] -= D::cs2 * GH[M001];
        GH[M021] -= D::cs2 * GH[M001];
        GH[M120] -= D::cs2 * GH[M100];
        GH[M102] -= D::cs2 * GH[M100];
        GH[M012] -= D::cs2 * GH[M010];

        GH[M220] -= (D::cs2 * (GH[M200] + GH[M020]) + cs4);
        GH[M202] -= (D::cs2 * (GH[M200] + GH[M002]) + cs4);
        GH[M022] -= (D::cs2 * (GH[M020] + GH[M002]) + cs4);
        GH[M211] -= (D::cs2 * GH[M011]);
        GH[M121] -= (D::cs2 * GH[M101]);
        GH[M112] -= (D::cs2 * GH[M110]);

        GH[M221] -= (D::cs2 * (GH[M201] + GH[M021]) + cs4 * GH[M001]);
        GH[M212] -= (D::cs2 * (GH[M210] + GH[M012]) + cs4 * GH[M010]);
        GH[M122] -= (D::cs2 * (GH[M120] + GH[M102]) + cs4 * GH[M100]);

        GH[M222] -=
            (D::cs2 * (GH[M220] + GH[M202] + GH[M022]) + cs4 * (GH[M200] + GH[M020] + GH[M002])
             + D::cs2 * cs4);
    };

    static void GHcomputeEquilibriumMoments(Array<T, D::d> const &u, Array<T, D::q> &GHeq)
    {
        // Order 0
        GHeq[M000] = 1.;
        // Order 1
        GHeq[M100] = u[0];
        GHeq[M010] = u[1];
        GHeq[M001] = u[2];
        // Order 2
        GHeq[M200] = u[0] * u[0];
        GHeq[M020] = u[1] * u[1];
        GHeq[M002] = u[2] * u[2];
        GHeq[M110] = u[0] * u[1];
        GHeq[M101] = u[0] * u[2];
        GHeq[M011] = u[1] * u[2];
        // Order 3
        GHeq[M210] = GHeq[M200] * u[1];
        GHeq[M201] = GHeq[M200] * u[2];
        GHeq[M021] = GHeq[M020] * u[2];
        GHeq[M120] = GHeq[M020] * u[0];
        GHeq[M102] = GHeq[M002] * u[0];
        GHeq[M012] = GHeq[M002] * u[1];
        GHeq[M111] = GHeq[M110] * u[2];
        // Order 4
        GHeq[M220] = GHeq[M200] * GHeq[M020];
        GHeq[M202] = GHeq[M200] * GHeq[M002];
        GHeq[M022] = GHeq[M020] * GHeq[M002];
        GHeq[M211] = GHeq[M200] * GHeq[M011];
        GHeq[M121] = GHeq[M020] * GHeq[M101];
        GHeq[M112] = GHeq[M002] * GHeq[M110];
        // Order 5
        GHeq[M221] = GHeq[M220] * GHeq[M001];
        GHeq[M212] = GHeq[M202] * GHeq[M010];
        GHeq[M122] = GHeq[M022] * GHeq[M100];
        // Order 6
        GHeq[M222] = GHeq[M220] * GHeq[M002];
    };

    // Equilibrium populations based on 27 moments can be computed using either RM, HM, CM, CHM or
    // Gauss-Hermite formalisms. All formulations are equivalent for the D3Q27 so we use the fastest
    // one (RMs)
    static void GHcomputeEquilibrium(T rho, Array<T, D::q> const &GHeq, Array<T, D::q> &eq)
    {
        Array<T, D::d> u(GHeq[1], GHeq[2], GHeq[3]);
        Array<T, D::q> RMeq;
        // Order 2
        RMeq[M200] = u[0] * u[0] + D::cs2;
        RMeq[M020] = u[1] * u[1] + D::cs2;
        RMeq[M002] = u[2] * u[2] + D::cs2;
        RMeq[M110] = u[0] * u[1];
        RMeq[M101] = u[0] * u[2];
        RMeq[M011] = u[1] * u[2];
        // Order 3
        RMeq[M210] = RMeq[M200] * u[1];
        RMeq[M201] = RMeq[M200] * u[2];
        RMeq[M021] = RMeq[M020] * u[2];
        RMeq[M120] = RMeq[M020] * u[0];
        RMeq[M102] = RMeq[M002] * u[0];
        RMeq[M012] = RMeq[M002] * u[1];
        RMeq[M111] = RMeq[M110] * u[2];
        // Order 4
        RMeq[M220] = RMeq[M200] * RMeq[M020];
        RMeq[M202] = RMeq[M200] * RMeq[M002];
        RMeq[M022] = RMeq[M020] * RMeq[M002];
        RMeq[M211] = RMeq[M200] * RMeq[M011];
        RMeq[M121] = RMeq[M020] * RMeq[M101];
        RMeq[M112] = RMeq[M002] * RMeq[M110];
        // Order 5
        RMeq[M221] = RMeq[M220] * u[2];
        RMeq[M212] = RMeq[M202] * u[1];
        RMeq[M122] = RMeq[M022] * u[0];
        // Order 6
        RMeq[M222] = RMeq[M220] * RMeq[M002];

        // Optimization based on symmetries between populations and their opposite counterpart
        eq[F000] = rho
                   * (1. - RMeq[M200] - RMeq[M020] - RMeq[M002] + RMeq[M220] + RMeq[M202]
                      + RMeq[M022] - RMeq[M222]);

        eq[FP00] = 0.5 * rho
                   * (u[0] + RMeq[M200] - RMeq[M120] - RMeq[M102] - RMeq[M220] - RMeq[M202]
                      + RMeq[M122] + RMeq[M222]);
        eq[FM00] = rho * (-u[0] + RMeq[M120] + RMeq[M102] - RMeq[M122]) + eq[FP00];

        eq[F0P0] = 0.5 * rho
                   * (u[1] + RMeq[M020] - RMeq[M210] - RMeq[M012] - RMeq[M220] - RMeq[M022]
                      + RMeq[M212] + RMeq[M222]);
        eq[F0M0] = rho * (-u[1] + RMeq[M210] + RMeq[M012] - RMeq[M212]) + eq[F0P0];

        eq[F00P] = 0.5 * rho
                   * (u[2] + RMeq[M002] - RMeq[M201] - RMeq[M021] - RMeq[M202] - RMeq[M022]
                      + RMeq[M221] + RMeq[M222]);
        eq[F00M] = rho * (-u[2] + RMeq[M201] + RMeq[M021] - RMeq[M221]) + eq[F00P];

        eq[FPP0] = 0.25 * rho
                   * (RMeq[M110] + RMeq[M210] + RMeq[M120] - RMeq[M112] + RMeq[M220] - RMeq[M212]
                      - RMeq[M122] - RMeq[M222]);
        eq[FMP0] = 0.5 * rho * (-RMeq[M110] - RMeq[M120] + RMeq[M112] + RMeq[M122]) + eq[FPP0];
        eq[FPM0] = 0.5 * rho * (-RMeq[M110] - RMeq[M210] + RMeq[M112] + RMeq[M212]) + eq[FPP0];
        eq[FMM0] = 0.5 * rho * (-RMeq[M210] - RMeq[M120] + RMeq[M212] + RMeq[M122]) + eq[FPP0];

        eq[FP0P] = 0.25 * rho
                   * (RMeq[M101] + RMeq[M201] + RMeq[M102] - RMeq[M121] + RMeq[M202] - RMeq[M221]
                      - RMeq[M122] - RMeq[M222]);
        eq[FM0P] = 0.5 * rho * (-RMeq[M101] - RMeq[M102] + RMeq[M121] + RMeq[M122]) + eq[FP0P];
        eq[FP0M] = 0.5 * rho * (-RMeq[M101] - RMeq[M201] + RMeq[M121] + RMeq[M221]) + eq[FP0P];
        eq[FM0M] = 0.5 * rho * (-RMeq[M201] - RMeq[M102] + RMeq[M221] + RMeq[M122]) + eq[FP0P];

        eq[F0PP] = 0.25 * rho
                   * (RMeq[M011] + RMeq[M021] + RMeq[M012] - RMeq[M211] + RMeq[M022] - RMeq[M221]
                      - RMeq[M212] - RMeq[M222]);
        eq[F0MP] = 0.5 * rho * (-RMeq[M011] - RMeq[M012] + RMeq[M211] + RMeq[M212]) + eq[F0PP];
        eq[F0PM] = 0.5 * rho * (-RMeq[M011] - RMeq[M021] + RMeq[M211] + RMeq[M221]) + eq[F0PP];
        eq[F0MM] = 0.5 * rho * (-RMeq[M021] - RMeq[M012] + RMeq[M221] + RMeq[M212]) + eq[F0PP];

        eq[FPPP] = 0.125 * rho
                   * (RMeq[M111] + RMeq[M211] + RMeq[M121] + RMeq[M112] + RMeq[M221] + RMeq[M212]
                      + RMeq[M122] + RMeq[M222]);
        eq[FMPP] = 0.25 * rho * (-RMeq[M111] - RMeq[M121] - RMeq[M112] - RMeq[M122]) + eq[FPPP];
        eq[FPMP] = 0.25 * rho * (-RMeq[M111] - RMeq[M211] - RMeq[M112] - RMeq[M212]) + eq[FPPP];
        eq[FPPM] = 0.25 * rho * (-RMeq[M111] - RMeq[M211] - RMeq[M121] - RMeq[M221]) + eq[FPPP];
        eq[FMMP] = 0.25 * rho * (-RMeq[M211] - RMeq[M121] - RMeq[M212] - RMeq[M122]) + eq[FPPP];
        eq[FMPM] = 0.25 * rho * (-RMeq[M211] - RMeq[M112] - RMeq[M221] - RMeq[M122]) + eq[FPPP];
        eq[FPMM] = 0.25 * rho * (-RMeq[M121] - RMeq[M112] - RMeq[M221] - RMeq[M212]) + eq[FPPP];
        eq[FMMM] = 0.25 * rho * (-RMeq[M111] - RMeq[M221] - RMeq[M212] - RMeq[M122]) + eq[FPPP];
    };

    static void GHcollide(
        Array<T, D::q> &cell, T rho, Array<T, D::d> const &u,
        Array<T, D::q> const &GH,    // Hermite moments
        Array<T, D::q> const &GHeq,  // Equilibrium moments (Hermite)
        Array<T, D::numRelaxationTimes> const &omega)
    {
        T omega1 = omega[0];
        T omega2 = omega[1];
        T omega3 = omega[2];
        T omega4 = omega[3];
        T omega5 = omega[4];
        T omega6 = omega[5];
        T omega7 = omega[6];
        T omega8 = omega[7];

        T omegaBulk = omega[8];
        T omegaPlus = (omegaBulk + 2. * omega1) / 3.;  // Notation used by Fei
        T omegaMinus = (omegaBulk - omega1) / 3.;      // Notation used by Fei

        // Post-collision moments.
        Array<T, D::q> GHcoll;

        // Collision in the Hermite moment space
        // Order 2 (non-diagonal collision so that we can easily modify the bulk viscosity)
        GHcoll[M200] = GH[M200] - omegaPlus * (GH[M200] - GHeq[M200])
                       - omegaMinus * (GH[M020] - GHeq[M020])
                       - omegaMinus * (GH[M002] - GHeq[M002]);
        GHcoll[M020] = GH[M020] - omegaMinus * (GH[M200] - GHeq[M200])
                       - omegaPlus * (GH[M020] - GHeq[M020]) - omegaMinus * (GH[M002] - GHeq[M002]);
        GHcoll[M002] = GH[M002] - omegaMinus * (GH[M200] - GHeq[M200])
                       - omegaMinus * (GH[M020] - GHeq[M020]) - omegaPlus * (GH[M002] - GHeq[M002]);

        GHcoll[M110] = (1. - omega2) * GH[M110] + omega2 * GHeq[M110];
        GHcoll[M101] = (1. - omega2) * GH[M101] + omega2 * GHeq[M101];
        GHcoll[M011] = (1. - omega2) * GH[M011] + omega2 * GHeq[M011];

        // Order 3
        GHcoll[M210] = (1. - omega3) * GH[M210] + omega3 * GHeq[M210];
        GHcoll[M201] = (1. - omega3) * GH[M201] + omega3 * GHeq[M201];
        GHcoll[M021] = (1. - omega3) * GH[M021] + omega3 * GHeq[M021];
        GHcoll[M120] = (1. - omega3) * GH[M120] + omega3 * GHeq[M120];
        GHcoll[M102] = (1. - omega3) * GH[M102] + omega3 * GHeq[M102];
        GHcoll[M012] = (1. - omega3) * GH[M012] + omega3 * GHeq[M012];

        GHcoll[M111] = (1. - omega4) * GH[M111] + omega4 * GHeq[M111];

        // Order 4
        GHcoll[M220] = (1. - omega5) * GH[M220] + omega5 * GHeq[M220];
        GHcoll[M202] = (1. - omega5) * GH[M202] + omega5 * GHeq[M202];
        GHcoll[M022] = (1. - omega5) * GH[M022] + omega5 * GHeq[M022];

        GHcoll[M211] = (1. - omega6) * GH[M211] + omega6 * GHeq[M211];
        GHcoll[M121] = (1. - omega6) * GH[M121] + omega6 * GHeq[M121];
        GHcoll[M112] = (1. - omega6) * GH[M112] + omega6 * GHeq[M112];

        // Order 5
        GHcoll[M221] = (1. - omega7) * GH[M221] + omega7 * GHeq[M221];
        GHcoll[M212] = (1. - omega7) * GH[M212] + omega7 * GHeq[M212];
        GHcoll[M122] = (1. - omega7) * GH[M122] + omega7 * GHeq[M122];

        // Order 6
        GHcoll[M222] = (1. - omega8) * GH[M222] + omega8 * GHeq[M222];

        // Compute post collision populations from Gauss-Hermite formalism
        cell[F000] = (rho * 8. / 27.)
                     * (1. - 3. / 2. * GHcoll[M200] - 3. / 2. * GHcoll[M020]
                        - 3. / 2. * GHcoll[M002] + 9. / 4. * GHcoll[M220] + 9. / 4. * GHcoll[M202]
                        + 9. / 4. * GHcoll[M022] - 27. / 8. * GHcoll[M222]);

        cell[FP00] = (rho * 2. / 27.)
                     * (1. + 3. * u[0] + 3. * GHcoll[M200] - 3. / 2. * GHcoll[M020]
                        - 3. / 2. * GHcoll[M002] - 9. / 2. * GHcoll[M120] - 9. / 2. * GHcoll[M102]
                        - 9. / 2. * GHcoll[M220] - 9. / 2. * GHcoll[M202] + 9. / 4. * GHcoll[M022]
                        + 27. / 4. * GHcoll[M122] + 27. / 4. * GHcoll[M222]);
        cell[FM00] = (rho * 2. / 27.)
                     * (1. - 3. * u[0] + 3. * GHcoll[M200] - 3. / 2. * GHcoll[M020]
                        - 3. / 2. * GHcoll[M002] + 9. / 2. * GHcoll[M120] + 9. / 2. * GHcoll[M102]
                        - 9. / 2. * GHcoll[M220] - 9. / 2. * GHcoll[M202] + 9. / 4. * GHcoll[M022]
                        - 27. / 4. * GHcoll[M122] + 27. / 4. * GHcoll[M222]);
        cell[F0P0] = (rho * 2. / 27.)
                     * (1. + 3. * u[1] - 3. / 2. * GHcoll[M200] + 3. * GHcoll[M020]
                        - 3. / 2. * GHcoll[M002] - 9. / 2. * GHcoll[M210] - 9. / 2. * GHcoll[M012]
                        - 9. / 2. * GHcoll[M220] + 9. / 4. * GHcoll[M202] - 9. / 2. * GHcoll[M022]
                        + 27. / 4. * GHcoll[M212] + 27. / 4. * GHcoll[M222]);
        cell[F0M0] = (rho * 2. / 27.)
                     * (1. - 3. * u[1] - 3. / 2. * GHcoll[M200] + 3. * GHcoll[M020]
                        - 3. / 2. * GHcoll[M002] + 9. / 2. * GHcoll[M210] + 9. / 2. * GHcoll[M012]
                        - 9. / 2. * GHcoll[M220] + 9. / 4. * GHcoll[M202] - 9. / 2. * GHcoll[M022]
                        - 27. / 4. * GHcoll[M212] + 27. / 4. * GHcoll[M222]);
        cell[F00P] = (rho * 2. / 27.)
                     * (1. + 3. * u[2] - 3. / 2. * GHcoll[M200] - 3. / 2. * GHcoll[M020]
                        + 3. * GHcoll[M002] - 9. / 2. * GHcoll[M201] - 9. / 2. * GHcoll[M021]
                        + 9. / 4. * GHcoll[M220] - 9. / 2. * GHcoll[M202] - 9. / 2. * GHcoll[M022]
                        + 27. / 4. * GHcoll[M221] + 27. / 4. * GHcoll[M222]);
        cell[F00M] = (rho * 2. / 27.)
                     * (1. - 3. * u[2] - 3. / 2. * GHcoll[M200] - 3. / 2. * GHcoll[M020]
                        + 3. * GHcoll[M002] + 9. / 2. * GHcoll[M201] + 9. / 2. * GHcoll[M021]
                        + 9. / 4. * GHcoll[M220] - 9. / 2. * GHcoll[M202] - 9. / 2. * GHcoll[M022]
                        - 27. / 4. * GHcoll[M221] + 27. / 4. * GHcoll[M222]);

        cell[FPP0] =
            (rho / 54.)
            * (1. + 3. * u[0] + 3. * u[1] + 3. * GHcoll[M200] + 3. * GHcoll[M020]
               - 3. / 2. * GHcoll[M002] + 9. * GHcoll[M110] + 9. * GHcoll[M210] + 9. * GHcoll[M120]
               - 9. / 2. * GHcoll[M102] - 9. / 2. * GHcoll[M012] + 9. * GHcoll[M220]
               - 9. / 2. * GHcoll[M202] - 9. / 2. * GHcoll[M022] - 27. / 2. * GHcoll[M112]
               - 27. / 2. * GHcoll[M212] - 27. / 2. * GHcoll[M122] - 27. / 2. * GHcoll[M222]);
        cell[FMP0] =
            (rho / 54.)
            * (1. - 3. * u[0] + 3. * u[1] + 3. * GHcoll[M200] + 3. * GHcoll[M020]
               - 3. / 2. * GHcoll[M002] - 9. * GHcoll[M110] + 9. * GHcoll[M210] - 9. * GHcoll[M120]
               + 9. / 2. * GHcoll[M102] - 9. / 2. * GHcoll[M012] + 9. * GHcoll[M220]
               - 9. / 2. * GHcoll[M202] - 9. / 2. * GHcoll[M022] + 27. / 2. * GHcoll[M112]
               - 27. / 2. * GHcoll[M212] + 27. / 2. * GHcoll[M122] - 27. / 2. * GHcoll[M222]);
        cell[FPM0] =
            (rho / 54.)
            * (1. + 3. * u[0] - 3. * u[1] + 3. * GHcoll[M200] + 3. * GHcoll[M020]
               - 3. / 2. * GHcoll[M002] - 9. * GHcoll[M110] - 9. * GHcoll[M210] + 9. * GHcoll[M120]
               - 9. / 2. * GHcoll[M102] + 9. / 2. * GHcoll[M012] + 9. * GHcoll[M220]
               - 9. / 2. * GHcoll[M202] - 9. / 2. * GHcoll[M022] + 27. / 2. * GHcoll[M112]
               + 27. / 2. * GHcoll[M212] - 27. / 2. * GHcoll[M122] - 27. / 2. * GHcoll[M222]);
        cell[FMM0] =
            (rho / 54.)
            * (1. - 3. * u[0] - 3. * u[1] + 3. * GHcoll[M200] + 3. * GHcoll[M020]
               - 3. / 2. * GHcoll[M002] + 9. * GHcoll[M110] - 9. * GHcoll[M210] - 9. * GHcoll[M120]
               + 9. / 2. * GHcoll[M102] + 9. / 2. * GHcoll[M012] + 9. * GHcoll[M220]
               - 9. / 2. * GHcoll[M202] - 9. / 2. * GHcoll[M022] - 27. / 2. * GHcoll[M112]
               + 27. / 2. * GHcoll[M212] + 27. / 2. * GHcoll[M122] - 27. / 2. * GHcoll[M222]);
        cell[FP0P] =
            (rho / 54.)
            * (1. + 3. * u[0] + 3. * u[2] + 3. * GHcoll[M200] - 3. / 2. * GHcoll[M020]
               + 3. * GHcoll[M002] + 9. * GHcoll[M101] + 9. * GHcoll[M201] - 9. / 2. * GHcoll[M021]
               - 9. / 2. * GHcoll[M120] + 9. * GHcoll[M102] - 9. / 2. * GHcoll[M220]
               + 9. * GHcoll[M202] - 9. / 2. * GHcoll[M022] - 27. / 2. * GHcoll[M121]
               - 27. / 2. * GHcoll[M221] - 27. / 2. * GHcoll[M122] - 27. / 2. * GHcoll[M222]);
        cell[FM0P] =
            (rho / 54.)
            * (1. - 3. * u[0] + 3. * u[2] + 3. * GHcoll[M200] - 3. / 2. * GHcoll[M020]
               + 3. * GHcoll[M002] - 9. * GHcoll[M101] + 9. * GHcoll[M201] - 9. / 2. * GHcoll[M021]
               + 9. / 2. * GHcoll[M120] - 9. * GHcoll[M102] - 9. / 2. * GHcoll[M220]
               + 9. * GHcoll[M202] - 9. / 2. * GHcoll[M022] + 27. / 2. * GHcoll[M121]
               - 27. / 2. * GHcoll[M221] + 27. / 2. * GHcoll[M122] - 27. / 2. * GHcoll[M222]);
        cell[FP0M] =
            (rho / 54.)
            * (1. + 3. * u[0] - 3. * u[2] + 3. * GHcoll[M200] - 3. / 2. * GHcoll[M020]
               + 3. * GHcoll[M002] - 9. * GHcoll[M101] - 9. * GHcoll[M201] + 9. / 2. * GHcoll[M021]
               - 9. / 2. * GHcoll[M120] + 9. * GHcoll[M102] - 9. / 2. * GHcoll[M220]
               + 9. * GHcoll[M202] - 9. / 2. * GHcoll[M022] + 27. / 2. * GHcoll[M121]
               + 27. / 2. * GHcoll[M221] - 27. / 2. * GHcoll[M122] - 27. / 2. * GHcoll[M222]);
        cell[FM0M] =
            (rho / 54.)
            * (1. - 3. * u[0] - 3. * u[2] + 3. * GHcoll[M200] - 3. / 2. * GHcoll[M020]
               + 3. * GHcoll[M002] + 9. * GHcoll[M101] - 9. * GHcoll[M201] + 9. / 2. * GHcoll[M021]
               + 9. / 2. * GHcoll[M120] - 9. * GHcoll[M102] - 9. / 2. * GHcoll[M220]
               + 9. * GHcoll[M202] - 9. / 2. * GHcoll[M022] - 27. / 2. * GHcoll[M121]
               + 27. / 2. * GHcoll[M221] + 27. / 2. * GHcoll[M122] - 27. / 2. * GHcoll[M222]);
        cell[F0PP] = (rho / 54.)
                     * (1. + 3. * u[1] + 3. * u[2] - 3. / 2. * GHcoll[M200] + 3. * GHcoll[M020]
                        + 3. * GHcoll[M002] + 9. * GHcoll[M011] - 9. / 2. * GHcoll[M210]
                        - 9. / 2. * GHcoll[M201] + 9. * GHcoll[M021] + 9. * GHcoll[M012]
                        - 9. / 2. * GHcoll[M220] - 9. / 2. * GHcoll[M202] + 9. * GHcoll[M022]
                        - 27. / 2. * GHcoll[M211] - 27. / 2. * GHcoll[M221]
                        - 27. / 2. * GHcoll[M212] - 27. / 2. * GHcoll[M222]);
        cell[F0MP] = (rho / 54.)
                     * (1. - 3. * u[1] + 3. * u[2] - 3. / 2. * GHcoll[M200] + 3. * GHcoll[M020]
                        + 3. * GHcoll[M002] - 9. * GHcoll[M011] + 9. / 2. * GHcoll[M210]
                        - 9. / 2. * GHcoll[M201] + 9. * GHcoll[M021] - 9. * GHcoll[M012]
                        - 9. / 2. * GHcoll[M220] - 9. / 2. * GHcoll[M202] + 9. * GHcoll[M022]
                        + 27. / 2. * GHcoll[M211] - 27. / 2. * GHcoll[M221]
                        + 27. / 2. * GHcoll[M212] - 27. / 2. * GHcoll[M222]);
        cell[F0PM] = (rho / 54.)
                     * (1. + 3. * u[1] - 3. * u[2] - 3. / 2. * GHcoll[M200] + 3. * GHcoll[M020]
                        + 3. * GHcoll[M002] - 9. * GHcoll[M011] - 9. / 2. * GHcoll[M210]
                        + 9. / 2. * GHcoll[M201] - 9. * GHcoll[M021] + 9. * GHcoll[M012]
                        - 9. / 2. * GHcoll[M220] - 9. / 2. * GHcoll[M202] + 9. * GHcoll[M022]
                        + 27. / 2. * GHcoll[M211] + 27. / 2. * GHcoll[M221]
                        - 27. / 2. * GHcoll[M212] - 27. / 2. * GHcoll[M222]);
        cell[F0MM] = (rho / 54.)
                     * (1. - 3. * u[1] - 3. * u[2] - 3. / 2. * GHcoll[M200] + 3. * GHcoll[M020]
                        + 3. * GHcoll[M002] + 9. * GHcoll[M011] + 9. / 2. * GHcoll[M210]
                        + 9. / 2. * GHcoll[M201] - 9. * GHcoll[M021] - 9. * GHcoll[M012]
                        - 9. / 2. * GHcoll[M220] - 9. / 2. * GHcoll[M202] + 9. * GHcoll[M022]
                        - 27. / 2. * GHcoll[M211] + 27. / 2. * GHcoll[M221]
                        + 27. / 2. * GHcoll[M212] - 27. / 2. * GHcoll[M222]);

        cell[FPPP] =
            (rho / 216.)
            * (1. + 3. * u[0] + 3. * u[1] + 3. * u[2] + 3. * GHcoll[M200] + 3. * GHcoll[M020]
               + 3. * GHcoll[M002] + 9. * GHcoll[M110] + 9. * GHcoll[M101] + 9. * GHcoll[M011]
               + 9. * GHcoll[M210] + 9. * GHcoll[M201] + 9. * GHcoll[M021] + 9. * GHcoll[M120]
               + 9. * GHcoll[M102] + 9. * GHcoll[M012] + 27. * GHcoll[M111] + 9. * GHcoll[M220]
               + 9. * GHcoll[M202] + 9. * GHcoll[M022] + 27. * GHcoll[M211] + 27. * GHcoll[M121]
               + 27. * GHcoll[M112] + 27. * GHcoll[M221] + 27. * GHcoll[M212] + 27. * GHcoll[M122]
               + 27. * GHcoll[M222]);
        cell[FMPP] =
            (rho / 216.)
            * (1. - 3. * u[0] + 3. * u[1] + 3. * u[2] + 3. * GHcoll[M200] + 3. * GHcoll[M020]
               + 3. * GHcoll[M002] - 9. * GHcoll[M110] - 9. * GHcoll[M101] + 9. * GHcoll[M011]
               + 9. * GHcoll[M210] + 9. * GHcoll[M201] + 9. * GHcoll[M021] - 9. * GHcoll[M120]
               - 9. * GHcoll[M102] + 9. * GHcoll[M012] - 27. * GHcoll[M111] + 9. * GHcoll[M220]
               + 9. * GHcoll[M202] + 9. * GHcoll[M022] + 27. * GHcoll[M211] - 27. * GHcoll[M121]
               - 27. * GHcoll[M112] + 27. * GHcoll[M221] + 27. * GHcoll[M212] - 27. * GHcoll[M122]
               + 27. * GHcoll[M222]);
        cell[FPMP] =
            (rho / 216.)
            * (1. + 3. * u[0] - 3. * u[1] + 3. * u[2] + 3. * GHcoll[M200] + 3. * GHcoll[M020]
               + 3. * GHcoll[M002] - 9. * GHcoll[M110] + 9. * GHcoll[M101] - 9. * GHcoll[M011]
               - 9. * GHcoll[M210] + 9. * GHcoll[M201] + 9. * GHcoll[M021] + 9. * GHcoll[M120]
               + 9. * GHcoll[M102] - 9. * GHcoll[M012] - 27. * GHcoll[M111] + 9. * GHcoll[M220]
               + 9. * GHcoll[M202] + 9. * GHcoll[M022] - 27. * GHcoll[M211] + 27. * GHcoll[M121]
               - 27. * GHcoll[M112] + 27. * GHcoll[M221] - 27. * GHcoll[M212] + 27. * GHcoll[M122]
               + 27. * GHcoll[M222]);
        cell[FPPM] =
            (rho / 216.)
            * (1. + 3. * u[0] + 3. * u[1] - 3. * u[2] + 3. * GHcoll[M200] + 3. * GHcoll[M020]
               + 3. * GHcoll[M002] + 9. * GHcoll[M110] - 9. * GHcoll[M101] - 9. * GHcoll[M011]
               + 9. * GHcoll[M210] - 9. * GHcoll[M201] - 9. * GHcoll[M021] + 9. * GHcoll[M120]
               + 9. * GHcoll[M102] + 9. * GHcoll[M012] - 27. * GHcoll[M111] + 9. * GHcoll[M220]
               + 9. * GHcoll[M202] + 9. * GHcoll[M022] - 27. * GHcoll[M211] - 27. * GHcoll[M121]
               + 27. * GHcoll[M112] - 27. * GHcoll[M221] + 27. * GHcoll[M212] + 27. * GHcoll[M122]
               + 27. * GHcoll[M222]);
        cell[FMMP] =
            (rho / 216.)
            * (1. - 3. * u[0] - 3. * u[1] + 3. * u[2] + 3. * GHcoll[M200] + 3. * GHcoll[M020]
               + 3. * GHcoll[M002] + 9. * GHcoll[M110] - 9. * GHcoll[M101] - 9. * GHcoll[M011]
               - 9. * GHcoll[M210] + 9. * GHcoll[M201] + 9. * GHcoll[M021] - 9. * GHcoll[M120]
               - 9. * GHcoll[M102] - 9. * GHcoll[M012] + 27. * GHcoll[M111] + 9. * GHcoll[M220]
               + 9. * GHcoll[M202] + 9. * GHcoll[M022] - 27. * GHcoll[M211] - 27. * GHcoll[M121]
               + 27. * GHcoll[M112] + 27. * GHcoll[M221] - 27. * GHcoll[M212] - 27. * GHcoll[M122]
               + 27. * GHcoll[M222]);
        cell[FMPM] =
            (rho / 216.)
            * (1. - 3. * u[0] + 3. * u[1] - 3. * u[2] + 3. * GHcoll[M200] + 3. * GHcoll[M020]
               + 3. * GHcoll[M002] - 9. * GHcoll[M110] + 9. * GHcoll[M101] - 9. * GHcoll[M011]
               + 9. * GHcoll[M210] - 9. * GHcoll[M201] - 9. * GHcoll[M021] - 9. * GHcoll[M120]
               - 9. * GHcoll[M102] + 9. * GHcoll[M012] + 27. * GHcoll[M111] + 9. * GHcoll[M220]
               + 9. * GHcoll[M202] + 9. * GHcoll[M022] - 27. * GHcoll[M211] + 27. * GHcoll[M121]
               - 27. * GHcoll[M112] - 27. * GHcoll[M221] + 27. * GHcoll[M212] - 27. * GHcoll[M122]
               + 27. * GHcoll[M222]);
        cell[FPMM] =
            (rho / 216.)
            * (1. + 3. * u[0] - 3. * u[1] - 3. * u[2] + 3. * GHcoll[M200] + 3. * GHcoll[M020]
               + 3. * GHcoll[M002] - 9. * GHcoll[M110] - 9. * GHcoll[M101] + 9. * GHcoll[M011]
               - 9. * GHcoll[M210] - 9. * GHcoll[M201] - 9. * GHcoll[M021] + 9. * GHcoll[M120]
               + 9. * GHcoll[M102] - 9. * GHcoll[M012] + 27. * GHcoll[M111] + 9. * GHcoll[M220]
               + 9. * GHcoll[M202] + 9. * GHcoll[M022] + 27. * GHcoll[M211] - 27. * GHcoll[M121]
               - 27. * GHcoll[M112] - 27. * GHcoll[M221] - 27. * GHcoll[M212] + 27. * GHcoll[M122]
               + 27. * GHcoll[M222]);
        cell[FMMM] =
            (rho / 216.)
            * (1. - 3. * u[0] - 3. * u[1] - 3. * u[2] + 3. * GHcoll[M200] + 3. * GHcoll[M020]
               + 3. * GHcoll[M002] + 9. * GHcoll[M110] + 9. * GHcoll[M101] + 9. * GHcoll[M011]
               - 9. * GHcoll[M210] - 9. * GHcoll[M201] - 9. * GHcoll[M021] - 9. * GHcoll[M120]
               - 9. * GHcoll[M102] - 9. * GHcoll[M012] - 27. * GHcoll[M111] + 9. * GHcoll[M220]
               + 9. * GHcoll[M202] + 9. * GHcoll[M022] + 27. * GHcoll[M211] + 27. * GHcoll[M121]
               + 27. * GHcoll[M112] - 27. * GHcoll[M221] - 27. * GHcoll[M212] - 27. * GHcoll[M122]
               + 27. * GHcoll[M222]);

        for (int i = 0; i < 27; ++i) {
            cell[i] -= D::SkordosFactor() * D::t[i];
        }
    };

    ///////////////////////////////////////////////////////////////////////////////
    // Recursive Regularization (RR) approach based on Gauss-Hermite formulation //
    ///////////////////////////////////////////////////////////////////////////////

    /**
     * // General way to compute RRs (we only need 2nd-order RRs)
     * static void RRcomputeMoments(Array<T,D::q> const& cell, Array<T, D::q>& RR, T& rho) {
     *
     *     Array<T, D::q> f;
     *     for (int i = 0; i<27; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *         RR[i] = 0.;
     *     }
     *
     *     T Hxx = 0.;
     *     T Hyy = 0.;
     *     T Hzz = 0.;
     *
     *     for (int i = 0; i<27; ++i) {
     *         Hxx = D::c[i][0] * D::c[i][0] - D::cs2;
     *         Hyy = D::c[i][1] * D::c[i][1] - D::cs2;
     *         Hzz = D::c[i][2] * D::c[i][2] - D::cs2;
     *
     *         // Order 0
     *         RR[M000] += f[i];
     *         // Order 1
     *         RR[M100] += D::c[i][0] * f[i];
     *         RR[M010] += D::c[i][1] * f[i];
     *         RR[M001] += D::c[i][2] * f[i];
     *         // Order 2
     *         RR[M200] += Hxx * f[i];
     *         RR[M020] += Hyy * f[i];
     *         RR[M002] += Hzz * f[i];
     *         RR[M110] += D::c[i][0] * D::c[i][1] * f[i];
     *         RR[M101] += D::c[i][0] * D::c[i][2] * f[i];
     *         RR[M011] += D::c[i][1] * D::c[i][2] * f[i];
     *     }
     *
     *     rho = RR[M000];
     *     T invRho = 1. / rho;
     *     for (int i = 0; i<27; ++i) {
     *         RR[i] *= invRho;
     *     }
     * };
     */

    // Optimized way to compute RMs based on the general ordering of discrete velocities
    // (we only need up to second order RRs since high-order ones are computed recursively)
    static void RRcomputeMoments(Array<T, D::q> const &cell, Array<T, D::q> &RR, T &rho)
    {
        Array<T, D::q> f;
        for (int i = 0; i < 27; ++i) {
            f[i] = cell[i] + D::SkordosFactor() * D::t[i];
            RR[i] = 0.;
        }

        T X_M1 =
            f[FM00] + f[FMM0] + f[FMP0] + f[FM0M] + f[FM0P] + f[FMMM] + f[FMMP] + f[FMPM] + f[FMPP];
        T X_P1 =
            f[FP00] + f[FPP0] + f[FPM0] + f[FP0P] + f[FP0M] + f[FPPP] + f[FPPM] + f[FPMP] + f[FPMM];
        T X_0 =
            f[F0M0] + f[F00M] + f[F0MM] + f[F0MP] + f[F000] + f[F0P0] + f[F00P] + f[F0PP] + f[F0PM];

        T Y_M1 =
            f[F0M0] + f[FMM0] + f[F0MM] + f[F0MP] + f[FMMM] + f[FMMP] + f[FPM0] + f[FPMP] + f[FPMM];
        T Y_P1 =
            f[F0P0] + f[FPP0] + f[F0PP] + f[F0PM] + f[FPPP] + f[FPPM] + f[FMP0] + f[FMPM] + f[FMPP];

        T Z_M1 =
            f[F00M] + f[FM0M] + f[F0MM] + f[FMMM] + f[FMPM] + f[FP0M] + f[F0PM] + f[FPPM] + f[FPMM];
        T Z_P1 =
            f[F00P] + f[FP0P] + f[F0PP] + f[FPPP] + f[FPMP] + f[FM0P] + f[F0MP] + f[FMMP] + f[FMPP];

        rho = X_M1 + X_P1 + X_0;
        T invRho = 1. / rho;

        // Order 0
        RR[M000] = 1.;
        // Order 1
        RR[M100] = invRho * (X_P1 - X_M1);
        RR[M010] = invRho * (Y_P1 - Y_M1);
        RR[M001] = invRho * (Z_P1 - Z_M1);
        // Order 2
        RR[M200] = invRho * (X_P1 + X_M1);
        RR[M020] = invRho * (Y_P1 + Y_M1);
        RR[M002] = invRho * (Z_P1 + Z_M1);
        RR[M110] = invRho
                   * (f[FMM0] - f[FMP0] + f[FMMM] + f[FMMP] - f[FMPM] - f[FMPP] + f[FPP0] - f[FPM0]
                      + f[FPPP] + f[FPPM] - f[FPMP] - f[FPMM]);
        RR[M101] = invRho
                   * (f[FM0M] - f[FM0P] + f[FMMM] - f[FMMP] + f[FMPM] - f[FMPP] + f[FP0P] - f[FP0M]
                      + f[FPPP] - f[FPPM] + f[FPMP] - f[FPMM]);
        RR[M011] = invRho
                   * (f[F0MM] - f[F0MP] + f[FMMM] - f[FMMP] - f[FMPM] + f[FMPP] + f[F0PP] - f[F0PM]
                      + f[FPPP] - f[FPPM] - f[FPMP] + f[FPMM]);

        // We come back to Hermite moments
        RR[M200] -= D::cs2;
        RR[M020] -= D::cs2;
        RR[M002] -= D::cs2;
    };

    static void RRcomputeEquilibriumMoments(Array<T, D::d> const &u, Array<T, D::q> &RReq)
    {
        // Order 0
        RReq[M000] = 1.;
        // Order 1
        RReq[M100] = u[0];
        RReq[M010] = u[1];
        RReq[M001] = u[2];
        // Order 2
        RReq[M200] = u[0] * u[0];
        RReq[M020] = u[1] * u[1];
        RReq[M002] = u[2] * u[2];
        RReq[M110] = u[0] * u[1];
        RReq[M101] = u[0] * u[2];
        RReq[M011] = u[1] * u[2];
        // Order 3
        RReq[M210] = RReq[M200] * u[1];
        RReq[M201] = RReq[M200] * u[2];
        RReq[M021] = RReq[M020] * u[2];
        RReq[M120] = RReq[M020] * u[0];
        RReq[M102] = RReq[M002] * u[0];
        RReq[M012] = RReq[M002] * u[1];
        RReq[M111] = RReq[M110] * u[2];
        // Order 4
        RReq[M220] = RReq[M200] * RReq[M020];
        RReq[M202] = RReq[M200] * RReq[M002];
        RReq[M022] = RReq[M020] * RReq[M002];
        RReq[M211] = RReq[M200] * RReq[M011];
        RReq[M121] = RReq[M020] * RReq[M101];
        RReq[M112] = RReq[M002] * RReq[M110];
        // Order 5
        RReq[M221] = RReq[M220] * RReq[M001];
        RReq[M212] = RReq[M202] * RReq[M010];
        RReq[M122] = RReq[M022] * RReq[M100];
        // Order 6
        RReq[M222] = RReq[M220] * RReq[M002];
    };

    // Equilibrium populations based on 27 moments can be computed using either RM, HM, CM, CHM or
    // Gauss-Hermite formalisms. All formulations are equivalent for the D3Q27 so we use the fastest
    // one (RMs)
    static void RRcomputeEquilibrium(T rho, Array<T, D::q> const &RReq, Array<T, D::q> &eq)
    {
        Array<T, D::d> u(RReq[1], RReq[2], RReq[3]);
        Array<T, D::q> RMeq;
        // Order 2
        RMeq[M200] = u[0] * u[0] + D::cs2;
        RMeq[M020] = u[1] * u[1] + D::cs2;
        RMeq[M002] = u[2] * u[2] + D::cs2;
        RMeq[M110] = u[0] * u[1];
        RMeq[M101] = u[0] * u[2];
        RMeq[M011] = u[1] * u[2];
        // Order 3
        RMeq[M210] = RMeq[M200] * u[1];
        RMeq[M201] = RMeq[M200] * u[2];
        RMeq[M021] = RMeq[M020] * u[2];
        RMeq[M120] = RMeq[M020] * u[0];
        RMeq[M102] = RMeq[M002] * u[0];
        RMeq[M012] = RMeq[M002] * u[1];
        RMeq[M111] = RMeq[M110] * u[2];
        // Order 4
        RMeq[M220] = RMeq[M200] * RMeq[M020];
        RMeq[M202] = RMeq[M200] * RMeq[M002];
        RMeq[M022] = RMeq[M020] * RMeq[M002];
        RMeq[M211] = RMeq[M200] * RMeq[M011];
        RMeq[M121] = RMeq[M020] * RMeq[M101];
        RMeq[M112] = RMeq[M002] * RMeq[M110];
        // Order 5
        RMeq[M221] = RMeq[M220] * u[2];
        RMeq[M212] = RMeq[M202] * u[1];
        RMeq[M122] = RMeq[M022] * u[0];
        // Order 6
        RMeq[M222] = RMeq[M220] * RMeq[M002];

        // Optimization based on symmetries between populations and their opposite counterpart
        eq[F000] = rho
                   * (1. - RMeq[M200] - RMeq[M020] - RMeq[M002] + RMeq[M220] + RMeq[M202]
                      + RMeq[M022] - RMeq[M222]);

        eq[FP00] = 0.5 * rho
                   * (u[0] + RMeq[M200] - RMeq[M120] - RMeq[M102] - RMeq[M220] - RMeq[M202]
                      + RMeq[M122] + RMeq[M222]);
        eq[FM00] = rho * (-u[0] + RMeq[M120] + RMeq[M102] - RMeq[M122]) + eq[FP00];

        eq[F0P0] = 0.5 * rho
                   * (u[1] + RMeq[M020] - RMeq[M210] - RMeq[M012] - RMeq[M220] - RMeq[M022]
                      + RMeq[M212] + RMeq[M222]);
        eq[F0M0] = rho * (-u[1] + RMeq[M210] + RMeq[M012] - RMeq[M212]) + eq[F0P0];

        eq[F00P] = 0.5 * rho
                   * (u[2] + RMeq[M002] - RMeq[M201] - RMeq[M021] - RMeq[M202] - RMeq[M022]
                      + RMeq[M221] + RMeq[M222]);
        eq[F00M] = rho * (-u[2] + RMeq[M201] + RMeq[M021] - RMeq[M221]) + eq[F00P];

        eq[FPP0] = 0.25 * rho
                   * (RMeq[M110] + RMeq[M210] + RMeq[M120] - RMeq[M112] + RMeq[M220] - RMeq[M212]
                      - RMeq[M122] - RMeq[M222]);
        eq[FMP0] = 0.5 * rho * (-RMeq[M110] - RMeq[M120] + RMeq[M112] + RMeq[M122]) + eq[FPP0];
        eq[FPM0] = 0.5 * rho * (-RMeq[M110] - RMeq[M210] + RMeq[M112] + RMeq[M212]) + eq[FPP0];
        eq[FMM0] = 0.5 * rho * (-RMeq[M210] - RMeq[M120] + RMeq[M212] + RMeq[M122]) + eq[FPP0];

        eq[FP0P] = 0.25 * rho
                   * (RMeq[M101] + RMeq[M201] + RMeq[M102] - RMeq[M121] + RMeq[M202] - RMeq[M221]
                      - RMeq[M122] - RMeq[M222]);
        eq[FM0P] = 0.5 * rho * (-RMeq[M101] - RMeq[M102] + RMeq[M121] + RMeq[M122]) + eq[FP0P];
        eq[FP0M] = 0.5 * rho * (-RMeq[M101] - RMeq[M201] + RMeq[M121] + RMeq[M221]) + eq[FP0P];
        eq[FM0M] = 0.5 * rho * (-RMeq[M201] - RMeq[M102] + RMeq[M221] + RMeq[M122]) + eq[FP0P];

        eq[F0PP] = 0.25 * rho
                   * (RMeq[M011] + RMeq[M021] + RMeq[M012] - RMeq[M211] + RMeq[M022] - RMeq[M221]
                      - RMeq[M212] - RMeq[M222]);
        eq[F0MP] = 0.5 * rho * (-RMeq[M011] - RMeq[M012] + RMeq[M211] + RMeq[M212]) + eq[F0PP];
        eq[F0PM] = 0.5 * rho * (-RMeq[M011] - RMeq[M021] + RMeq[M211] + RMeq[M221]) + eq[F0PP];
        eq[F0MM] = 0.5 * rho * (-RMeq[M021] - RMeq[M012] + RMeq[M221] + RMeq[M212]) + eq[F0PP];

        eq[FPPP] = 0.125 * rho
                   * (RMeq[M111] + RMeq[M211] + RMeq[M121] + RMeq[M112] + RMeq[M221] + RMeq[M212]
                      + RMeq[M122] + RMeq[M222]);
        eq[FMPP] = 0.25 * rho * (-RMeq[M111] - RMeq[M121] - RMeq[M112] - RMeq[M122]) + eq[FPPP];
        eq[FPMP] = 0.25 * rho * (-RMeq[M111] - RMeq[M211] - RMeq[M112] - RMeq[M212]) + eq[FPPP];
        eq[FPPM] = 0.25 * rho * (-RMeq[M111] - RMeq[M211] - RMeq[M121] - RMeq[M221]) + eq[FPPP];
        eq[FMMP] = 0.25 * rho * (-RMeq[M211] - RMeq[M121] - RMeq[M212] - RMeq[M122]) + eq[FPPP];
        eq[FMPM] = 0.25 * rho * (-RMeq[M211] - RMeq[M112] - RMeq[M221] - RMeq[M122]) + eq[FPPP];
        eq[FPMM] = 0.25 * rho * (-RMeq[M121] - RMeq[M112] - RMeq[M221] - RMeq[M212]) + eq[FPPP];
        eq[FMMM] = 0.25 * rho * (-RMeq[M111] - RMeq[M221] - RMeq[M212] - RMeq[M122]) + eq[FPPP];
    };

    static void RRcollide(
        Array<T, D::q> &cell, T rho, Array<T, D::d> const &u,
        Array<T, D::q> const &RR,    // Hermite moments
        Array<T, D::q> const &RReq,  // Equilibrium moments (Hermite)
        Array<T, D::numRelaxationTimes> const &omega)
    {
        T omega1 = omega[0];
        T omega2 = omega[1];
        T omega3 = omega[2];
        T omega4 = omega[3];
        T omega5 = omega[4];
        T omega6 = omega[5];
        T omega7 = omega[6];
        T omega8 = omega[7];

        T omegaBulk = omega[8];
        T omegaPlus = (omegaBulk + 2. * omega1) / 3.;  // Notation used by Fei
        T omegaMinus = (omegaBulk - omega1) / 3.;      // Notation used by Fei

        // Post-collision and Nonequilibrium moments.
        Array<T, D::q> RRneq;
        Array<T, D::q> RRcoll;
        Array<T, D::q> RMcoll;

        // Recursive computation of nonequilibrium Hermite moments
        // Order 2 (standard way to compute them)
        RRneq[M200] = RR[M200] - RReq[M200];
        RRneq[M020] = RR[M020] - RReq[M020];
        RRneq[M002] = RR[M002] - RReq[M002];
        RRneq[M110] = RR[M110] - RReq[M110];
        RRneq[M101] = RR[M101] - RReq[M101];
        RRneq[M011] = RR[M011] - RReq[M011];

        // Order 3 (reconstruction using Chapman-Enskog formulas)
        RRneq[M210] = u[1] * RRneq[M200] + 2. * u[0] * RRneq[M110];
        RRneq[M201] = u[2] * RRneq[M200] + 2. * u[0] * RRneq[M101];
        RRneq[M021] = u[2] * RRneq[M020] + 2. * u[1] * RRneq[M011];
        RRneq[M120] = u[0] * RRneq[M020] + 2. * u[1] * RRneq[M110];
        RRneq[M102] = u[0] * RRneq[M002] + 2. * u[2] * RRneq[M101];
        RRneq[M012] = u[1] * RRneq[M002] + 2. * u[2] * RRneq[M011];
        RRneq[M111] = u[2] * RRneq[M110] + u[1] * RRneq[M101] + u[0] * RRneq[M011];

        // Order 4 (reconstruction using Chapman-Enskog formulas)
        RRneq[M220] =
            u[1] * u[1] * RRneq[M200] + u[0] * u[0] * RRneq[M020] + 4. * u[0] * u[1] * RRneq[M110];
        RRneq[M202] =
            u[2] * u[2] * RRneq[M200] + u[0] * u[0] * RRneq[M002] + 4. * u[0] * u[2] * RRneq[M101];
        RRneq[M022] =
            u[2] * u[2] * RRneq[M020] + u[1] * u[1] * RRneq[M002] + 4. * u[1] * u[2] * RRneq[M011];
        RRneq[M211] = u[1] * u[2] * RRneq[M200] + 2. * u[0] * u[2] * RRneq[M110]
                      + 2. * u[0] * u[1] * RRneq[M101] + u[0] * u[0] * RRneq[M011];
        RRneq[M121] = u[0] * u[2] * RRneq[M020] + 2. * u[1] * u[2] * RRneq[M110]
                      + u[1] * u[1] * RRneq[M101] + 2. * u[0] * u[1] * RRneq[M011];
        RRneq[M112] = u[0] * u[1] * RRneq[M002] + u[2] * u[2] * RRneq[M110]
                      + 2. * u[1] * u[2] * RRneq[M101] + 2. * u[0] * u[2] * RRneq[M011];

        // Order 5 (reconstruction using Chapman-Enskog formulas)
        RRneq[M221] = u[1] * u[1] * u[2] * RRneq[M200] + u[0] * u[0] * u[2] * RRneq[M020]
                      + 4. * u[0] * u[1] * u[2] * RRneq[M110]
                      + 2. * u[0] * u[1] * u[1] * RRneq[M101]
                      + 2. * u[0] * u[0] * u[1] * RRneq[M011];
        RRneq[M212] = u[2] * u[2] * u[1] * RRneq[M200] + u[0] * u[0] * u[1] * RRneq[M002]
                      + 2. * u[0] * u[2] * u[2] * RRneq[M110]
                      + 4. * u[0] * u[1] * u[2] * RRneq[M101]
                      + 2. * u[0] * u[0] * u[2] * RRneq[M011];
        RRneq[M122] = u[2] * u[2] * u[0] * RRneq[M020] + u[1] * u[1] * u[0] * RRneq[M002]
                      + 2. * u[1] * u[2] * u[2] * RRneq[M110]
                      + 2. * u[1] * u[1] * u[2] * RRneq[M101]
                      + 4. * u[0] * u[1] * u[2] * RRneq[M011];

        // Order 6 (reconstruction using Chapman-Enskog formulas)
        RRneq[M222] =
            u[1] * u[1] * u[2] * u[2] * RRneq[M200] + u[0] * u[0] * u[2] * u[2] * RRneq[M020]
            + u[0] * u[0] * u[1] * u[1] * RRneq[M002] + 4. * u[0] * u[1] * u[2] * u[2] * RRneq[M110]
            + 4. * u[0] * u[1] * u[1] * u[2] * RRneq[M101]
            + 4. * u[0] * u[0] * u[1] * u[2] * RRneq[M011];

        // Collision in the Hermite moment space
        // Order 2 (non-diagonal collision so that we can easily modify the bulk viscosity)
        RRcoll[M200] = RR[M200] - omegaPlus * RRneq[M200] - omegaMinus * RRneq[M020]
                       - omegaMinus * RRneq[M002];
        RRcoll[M020] = RR[M020] - omegaMinus * RRneq[M200] - omegaPlus * RRneq[M020]
                       - omegaMinus * RRneq[M002];
        RRcoll[M002] = RR[M002] - omegaMinus * RRneq[M200] - omegaMinus * RRneq[M020]
                       - omegaPlus * RRneq[M002];

        RRcoll[M110] = (1. - omega2) * RRneq[M110] + RReq[M110];
        RRcoll[M101] = (1. - omega2) * RRneq[M101] + RReq[M101];
        RRcoll[M011] = (1. - omega2) * RRneq[M011] + RReq[M011];

        // Order 3
        RRcoll[M210] = (1. - omega3) * RRneq[M210] + RReq[M210];
        RRcoll[M201] = (1. - omega3) * RRneq[M201] + RReq[M201];
        RRcoll[M021] = (1. - omega3) * RRneq[M021] + RReq[M021];
        RRcoll[M120] = (1. - omega3) * RRneq[M120] + RReq[M120];
        RRcoll[M102] = (1. - omega3) * RRneq[M102] + RReq[M102];
        RRcoll[M012] = (1. - omega3) * RRneq[M012] + RReq[M012];

        RRcoll[M111] = (1. - omega4) * RRneq[M111] + RReq[M111];

        // Order 4
        RRcoll[M220] = (1. - omega5) * RRneq[M220] + RReq[M220];
        RRcoll[M202] = (1. - omega5) * RRneq[M202] + RReq[M202];
        RRcoll[M022] = (1. - omega5) * RRneq[M022] + RReq[M022];

        RRcoll[M211] = (1. - omega6) * RRneq[M211] + RReq[M211];
        RRcoll[M121] = (1. - omega6) * RRneq[M121] + RReq[M121];
        RRcoll[M112] = (1. - omega6) * RRneq[M112] + RReq[M112];

        // Order 5
        RRcoll[M221] = (1. - omega7) * RRneq[M221] + RReq[M221];
        RRcoll[M212] = (1. - omega7) * RRneq[M212] + RReq[M212];
        RRcoll[M122] = (1. - omega7) * RRneq[M122] + RReq[M122];

        // Order 6
        RRcoll[M222] = (1. - omega8) * RRneq[M222] + RReq[M222];

        // Come back to RMcoll using relationships between GHs and RMs
        T cs4 = D::cs2 * D::cs2;

        RMcoll[M200] = RRcoll[M200] + D::cs2;
        RMcoll[M020] = RRcoll[M020] + D::cs2;
        RMcoll[M002] = RRcoll[M002] + D::cs2;

        RMcoll[M110] = RRcoll[M110];
        RMcoll[M101] = RRcoll[M101];
        RMcoll[M011] = RRcoll[M011];

        RMcoll[M210] = RRcoll[M210] + D::cs2 * u[1];
        RMcoll[M201] = RRcoll[M201] + D::cs2 * u[2];
        RMcoll[M021] = RRcoll[M021] + D::cs2 * u[2];
        RMcoll[M120] = RRcoll[M120] + D::cs2 * u[0];
        RMcoll[M102] = RRcoll[M102] + D::cs2 * u[0];
        RMcoll[M012] = RRcoll[M012] + D::cs2 * u[1];

        RMcoll[M111] = RRcoll[M111];

        RMcoll[M220] = RRcoll[M220] + D::cs2 * (RRcoll[M200] + RRcoll[M020]) + cs4;
        RMcoll[M202] = RRcoll[M202] + D::cs2 * (RRcoll[M200] + RRcoll[M002]) + cs4;
        RMcoll[M022] = RRcoll[M022] + D::cs2 * (RRcoll[M020] + RRcoll[M002]) + cs4;

        RMcoll[M211] = RRcoll[M211] + D::cs2 * RRcoll[M011];
        RMcoll[M121] = RRcoll[M121] + D::cs2 * RRcoll[M101];
        RMcoll[M112] = RRcoll[M112] + D::cs2 * RRcoll[M110];

        RMcoll[M221] = RRcoll[M221] + D::cs2 * (RRcoll[M201] + RRcoll[M021]) + cs4 * u[2];
        RMcoll[M212] = RRcoll[M212] + D::cs2 * (RRcoll[M210] + RRcoll[M012]) + cs4 * u[1];
        RMcoll[M122] = RRcoll[M122] + D::cs2 * (RRcoll[M120] + RRcoll[M102]) + cs4 * u[0];

        RMcoll[M222] = RRcoll[M222] + D::cs2 * (RRcoll[M220] + RRcoll[M202] + RRcoll[M022])
                       + cs4 * (RRcoll[M200] + RRcoll[M020] + RRcoll[M002]) + D::cs2 * cs4;

        // Compute post collision populations from RM
        // Optimization based on symmetries between populations and their opposite counterpart
        cell[F000] = rho
                     * (1. - RMcoll[M200] - RMcoll[M020] - RMcoll[M002] + RMcoll[M220]
                        + RMcoll[M202] + RMcoll[M022] - RMcoll[M222]);

        cell[FP00] = 0.5 * rho
                     * (u[0] + RMcoll[M200] - RMcoll[M120] - RMcoll[M102] - RMcoll[M220]
                        - RMcoll[M202] + RMcoll[M122] + RMcoll[M222]);
        cell[FM00] = rho * (-u[0] + RMcoll[M120] + RMcoll[M102] - RMcoll[M122]) + cell[FP00];

        cell[F0P0] = 0.5 * rho
                     * (u[1] + RMcoll[M020] - RMcoll[M210] - RMcoll[M012] - RMcoll[M220]
                        - RMcoll[M022] + RMcoll[M212] + RMcoll[M222]);
        cell[F0M0] = rho * (-u[1] + RMcoll[M210] + RMcoll[M012] - RMcoll[M212]) + cell[F0P0];

        cell[F00P] = 0.5 * rho
                     * (u[2] + RMcoll[M002] - RMcoll[M201] - RMcoll[M021] - RMcoll[M202]
                        - RMcoll[M022] + RMcoll[M221] + RMcoll[M222]);
        cell[F00M] = rho * (-u[2] + RMcoll[M201] + RMcoll[M021] - RMcoll[M221]) + cell[F00P];

        cell[FPP0] = 0.25 * rho
                     * (RMcoll[M110] + RMcoll[M210] + RMcoll[M120] - RMcoll[M112] + RMcoll[M220]
                        - RMcoll[M212] - RMcoll[M122] - RMcoll[M222]);
        cell[FMP0] =
            0.5 * rho * (-RMcoll[M110] - RMcoll[M120] + RMcoll[M112] + RMcoll[M122]) + cell[FPP0];
        cell[FPM0] =
            0.5 * rho * (-RMcoll[M110] - RMcoll[M210] + RMcoll[M112] + RMcoll[M212]) + cell[FPP0];
        cell[FMM0] =
            0.5 * rho * (-RMcoll[M210] - RMcoll[M120] + RMcoll[M212] + RMcoll[M122]) + cell[FPP0];

        cell[FP0P] = 0.25 * rho
                     * (RMcoll[M101] + RMcoll[M201] + RMcoll[M102] - RMcoll[M121] + RMcoll[M202]
                        - RMcoll[M221] - RMcoll[M122] - RMcoll[M222]);
        cell[FM0P] =
            0.5 * rho * (-RMcoll[M101] - RMcoll[M102] + RMcoll[M121] + RMcoll[M122]) + cell[FP0P];
        cell[FP0M] =
            0.5 * rho * (-RMcoll[M101] - RMcoll[M201] + RMcoll[M121] + RMcoll[M221]) + cell[FP0P];
        cell[FM0M] =
            0.5 * rho * (-RMcoll[M201] - RMcoll[M102] + RMcoll[M221] + RMcoll[M122]) + cell[FP0P];

        cell[F0PP] = 0.25 * rho
                     * (RMcoll[M011] + RMcoll[M021] + RMcoll[M012] - RMcoll[M211] + RMcoll[M022]
                        - RMcoll[M221] - RMcoll[M212] - RMcoll[M222]);
        cell[F0MP] =
            0.5 * rho * (-RMcoll[M011] - RMcoll[M012] + RMcoll[M211] + RMcoll[M212]) + cell[F0PP];
        cell[F0PM] =
            0.5 * rho * (-RMcoll[M011] - RMcoll[M021] + RMcoll[M211] + RMcoll[M221]) + cell[F0PP];
        cell[F0MM] =
            0.5 * rho * (-RMcoll[M021] - RMcoll[M012] + RMcoll[M221] + RMcoll[M212]) + cell[F0PP];

        cell[FPPP] = 0.125 * rho
                     * (RMcoll[M111] + RMcoll[M211] + RMcoll[M121] + RMcoll[M112] + RMcoll[M221]
                        + RMcoll[M212] + RMcoll[M122] + RMcoll[M222]);
        cell[FMPP] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M121] - RMcoll[M112] - RMcoll[M122]) + cell[FPPP];
        cell[FPMP] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M211] - RMcoll[M112] - RMcoll[M212]) + cell[FPPP];
        cell[FPPM] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M211] - RMcoll[M121] - RMcoll[M221]) + cell[FPPP];
        cell[FMMP] =
            0.25 * rho * (-RMcoll[M211] - RMcoll[M121] - RMcoll[M212] - RMcoll[M122]) + cell[FPPP];
        cell[FMPM] =
            0.25 * rho * (-RMcoll[M211] - RMcoll[M112] - RMcoll[M221] - RMcoll[M122]) + cell[FPPP];
        cell[FPMM] =
            0.25 * rho * (-RMcoll[M121] - RMcoll[M112] - RMcoll[M221] - RMcoll[M212]) + cell[FPPP];
        cell[FMMM] =
            0.25 * rho * (-RMcoll[M111] - RMcoll[M221] - RMcoll[M212] - RMcoll[M122]) + cell[FPPP];

        for (int i = 0; i < 27; ++i) {
            cell[i] -= D::SkordosFactor() * D::t[i];
        }
    };

};  // struct comprehensiveDynamicsTemplatesImpl<T, descriptors::D3Q27DescriptorBase<T> >

// Efficient specialization for D3Q19 lattice
template <typename T>
struct comprehensiveDynamicsTemplatesImpl<T, descriptors::D3Q19DescriptorBase<T> > {
    typedef descriptors::D3Q19DescriptorBase<T> D;
    // Same order as in E1.
    enum {
        // Order 0
        M000 = 0,

        // Order 1
        M100 = 1,
        M010 = 2,
        M001 = 3,

        // Order 2
        M200 = 4,
        M020 = 5,
        M002 = 6,
        M110 = 7,
        M101 = 8,
        M011 = 9,

        // Order 3
        M210 = 10,
        M201 = 11,
        M021 = 12,
        M120 = 13,
        M102 = 14,
        M012 = 15,

        // Order 4
        M220 = 16,
        M202 = 17,
        M022 = 18,

    };

    enum {
        F000 = 0,
        FM00 = 1,
        F0M0 = 2,
        F00M = 3,
        FMM0 = 4,
        FMP0 = 5,
        FM0M = 6,
        FM0P = 7,
        F0MM = 8,
        F0MP = 9,

        FP00 = 10,
        F0P0 = 11,
        F00P = 12,
        FPP0 = 13,
        FPM0 = 14,
        FP0P = 15,
        FP0M = 16,
        F0PP = 17,
        F0PM = 18
    };

    /**
     * // General way to compute RMs
     * static void RMcomputeMoments(Array<T,D::q> const& cell, Array<T, D::q>& RM, T& rho) {
     *
     *     Array<T, D::q> f;
     *     for (int i = 0; i<19; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *         RM[i] = 0.;
     *     }
     *
     *     for (int i = 0; i<19; ++i) {
     *         // Order 0
     *         RM[M000] += f[i];
     *         // Order 1
     *         RM[M100] += D::c[i][0] * f[i];
     *         RM[M010] += D::c[i][1] * f[i];
     *         RM[M001] += D::c[i][2] * f[i];
     *         // Order 2
     *         RM[M200] += D::c[i][0] * D::c[i][0] * f[i];
     *         RM[M020] += D::c[i][1] * D::c[i][1] * f[i];
     *         RM[M002] += D::c[i][2] * D::c[i][2] * f[i];
     *         RM[M110] += D::c[i][0] * D::c[i][1] * f[i];
     *         RM[M101] += D::c[i][0] * D::c[i][2] * f[i];
     *         RM[M011] += D::c[i][1] * D::c[i][2] * f[i];
     *         // Order 3
     *         RM[M210] += D::c[i][0] * D::c[i][0] * D::c[i][1] * f[i];
     *         RM[M201] += D::c[i][0] * D::c[i][0] * D::c[i][2] * f[i];
     *         RM[M021] += D::c[i][1] * D::c[i][1] * D::c[i][2] * f[i];
     *         RM[M120] += D::c[i][0] * D::c[i][1] * D::c[i][1] * f[i];
     *         RM[M102] += D::c[i][0] * D::c[i][2] * D::c[i][2] * f[i];
     *         RM[M012] += D::c[i][1] * D::c[i][2] * D::c[i][2] * f[i];
     *         // Order 4
     *         RM[M220] += D::c[i][0] * D::c[i][0] * D::c[i][1] * D::c[i][1] * f[i];
     *         RM[M202] += D::c[i][0] * D::c[i][0] * D::c[i][2] * D::c[i][2] * f[i];
     *         RM[M022] += D::c[i][1] * D::c[i][1] * D::c[i][2] * D::c[i][2] * f[i];
     *     }
     *
     *     rho = RM[M000];
     *     T invRho = 1. / RM[M000];
     *     for (int i = 0; i<19; ++i) {
     *         RM[i] *= invRho;
     *     }
     * };
     *
     * // Optimized way to compute RMs based on Palabos ordering of discrete velocities
     * static void RMcomputeMoments(Array<T,D::q> const& cell, Array<T, D::q>& RM, T& rho) {
     *
     *     Array<T, D::q> f;
     *     for (int i = 0; i<19; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *     }
     *
     *     // Order 0
     *     rho =   f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[10] +
     * f[11] + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18]; RM[M000] = 1.0; T invRho
     * = 1./rho;
     *
     *     // Order 1
     *     RM[M100] = invRho * (- f[1] - f[4] - f[5] - f[6] - f[7] + f[10] + f[13] + f[14] + f[15] +
     * f[16]); RM[M010] = invRho * (- f[2] - f[4] + f[5] - f[8] - f[9] + f[11] + f[13] - f[14] +
     * f[17] + f[18]); RM[M001] = invRho * (- f[3] - f[6] + f[7] - f[8] + f[9] + f[12] + f[15] -
     * f[16] + f[17] - f[18]);
     *      // Order 2
     *     RM[M200] = invRho * (  f[1] + f[4] + f[5] + f[6] + f[7] + f[10] + f[13] + f[14] + f[15] +
     * f[16]); RM[M020] = invRho * (  f[2] + f[4] + f[5] + f[8] + f[9] + f[11] + f[13] + f[14] +
     * f[17] + f[18]); RM[M002] = invRho * (  f[3] + f[6] + f[7] + f[8] + f[9] + f[12] + f[15] +
     * f[16] + f[17] + f[18]);
     *
     *     RM[M110] = invRho * (  f[4] - f[5] + f[13] - f[14]);
     *     RM[M101] = invRho * (  f[6] - f[7] + f[15] - f[16]);
     *     RM[M011] = invRho * (  f[8] - f[9] + f[17] - f[18]);
     *      // Order 3
     *     RM[M210] = invRho * (- f[4] + f[5] + f[13] - f[14]);
     *     RM[M201] = invRho * (- f[6] + f[7] + f[15] - f[16]);
     *     RM[M021] = invRho * (- f[8] + f[9] + f[17] - f[18]);
     *     RM[M120] = invRho * (- f[4] - f[5] + f[13] + f[14]);
     *     RM[M102] = invRho * (- f[6] - f[7] + f[15] + f[16]);
     *     RM[M012] = invRho * (- f[8] - f[9] + f[17] + f[18]);
     *      // Order 4
     *     RM[M220] = invRho * (  f[4] + f[5] + f[13] + f[14]);
     *     RM[M202] = invRho * (  f[6] + f[7] + f[15] + f[16]);
     *     RM[M022] = invRho * (  f[8] + f[9] + f[17] + f[18]);
     * };
     */

    // Optimized way to compute RMs based on the general ordering of discrete velocities
    static void RMcomputeMoments(Array<T, D::q> const &cell, Array<T, D::q> &RM, T &rho)
    {
        Array<T, D::q> f;
        for (int i = 0; i < 19; ++i) {
            f[i] = cell[i] + D::SkordosFactor() * D::t[i];
        }
        T X_M1 = f[FM00] + f[FMM0] + f[FMP0] + f[FM0M] + f[FM0P];
        T X_P1 = f[FP00] + f[FPP0] + f[FPM0] + f[FP0P] + f[FP0M];
        T X_0 =
            f[F000] + f[F0M0] + f[F00M] + f[F0MM] + f[F0MP] + f[F0P0] + f[F00P] + f[F0PP] + f[F0PM];
        T Y_M1 = f[F0M0] + f[FMM0] + f[F0MM] + f[F0MP] + f[FPM0];
        T Y_P1 = f[FMP0] + f[F0P0] + f[FPP0] + f[F0PP] + f[F0PM];
        T Z_M1 = f[F00M] + f[FM0M] + f[F0MM] + f[FP0M] + f[F0PM];
        T Z_P1 = f[FM0P] + f[F0MP] + f[F00P] + f[FP0P] + f[F0PP];

        rho = X_M1 + X_P1 + X_0;
        T invRho = 1. / rho;
        T two_invRho = 2. * invRho;

        // Order 0
        RM[M000] = 1.;
        // Order 1
        RM[M100] = invRho * (X_P1 - X_M1);
        RM[M010] = invRho * (Y_P1 - Y_M1);
        RM[M001] = invRho * (Z_P1 - Z_M1);
        // Order 2
        RM[M200] = invRho * (X_M1 + X_P1);
        RM[M020] = invRho * (Y_M1 + Y_P1);
        RM[M002] = invRho * (Z_M1 + Z_P1);
        RM[M110] = invRho * (f[FMM0] - f[FMP0] + f[FPP0] - f[FPM0]);
        RM[M101] = invRho * (f[FM0M] - f[FM0P] + f[FP0P] - f[FP0M]);
        RM[M011] = invRho * (f[F0MM] - f[F0MP] + f[F0PP] - f[F0PM]);
        // Order 3
        RM[M210] = RM[M110] - two_invRho * (f[FMM0] - f[FMP0]);
        RM[M201] = RM[M101] - two_invRho * (f[FM0M] - f[FM0P]);
        RM[M021] = RM[M011] - two_invRho * (f[F0MM] - f[F0MP]);
        RM[M120] = RM[M110] - two_invRho * (f[FMM0] - f[FPM0]);
        RM[M102] = RM[M101] - two_invRho * (f[FM0M] - f[FP0M]);
        RM[M012] = RM[M011] - two_invRho * (f[F0MM] - f[F0PM]);
        // Order 4
        RM[M220] = RM[M110] + two_invRho * (f[FMP0] + f[FPM0]);
        RM[M202] = RM[M101] + two_invRho * (f[FM0P] + f[FP0M]);
        RM[M022] = RM[M011] + two_invRho * (f[F0MP] + f[F0PM]);
    };

    static void RMcomputeEquilibriumMoments(Array<T, D::d> const &u, Array<T, D::q> &RMeq)
    {
        // Order 0
        RMeq[M000] = 1.;
        // Order 1
        RMeq[M100] = u[0];
        RMeq[M010] = u[1];
        RMeq[M001] = u[2];
        // Order 2
        RMeq[M200] = u[0] * u[0] + D::cs2;
        RMeq[M020] = u[1] * u[1] + D::cs2;
        RMeq[M002] = u[2] * u[2] + D::cs2;
        RMeq[M110] = u[0] * u[1];
        RMeq[M101] = u[0] * u[2];
        RMeq[M011] = u[1] * u[2];
        // Order 3
        RMeq[M210] = RMeq[M200] * u[1];
        RMeq[M201] = RMeq[M200] * u[2];
        RMeq[M021] = RMeq[M020] * u[2];
        RMeq[M120] = RMeq[M020] * u[0];
        RMeq[M102] = RMeq[M002] * u[0];
        RMeq[M012] = RMeq[M002] * u[1];
        // Order 4
        RMeq[M220] = RMeq[M200] * RMeq[M020];
        RMeq[M202] = RMeq[M200] * RMeq[M002];
        RMeq[M022] = RMeq[M020] * RMeq[M002];
    };

    // Equilibrium populations based on 19 moments can be computed using either RM, HM, CM, CHM or
    // Gauss-Hermite formalisms. Here we use raw moments (RMs)
    static void RMcomputeEquilibrium(T rho, Array<T, D::q> const &RMeq, Array<T, D::q> &eq)
    {
        Array<T, D::d> u(RMeq[1], RMeq[2], RMeq[3]);

        eq[F000] =
            rho
            * (1. - RMeq[M200] - RMeq[M020] - RMeq[M002] + RMeq[M220] + RMeq[M202] + RMeq[M022]);

        eq[FP00] =
            0.5 * rho * (u[0] + RMeq[M200] - RMeq[M120] - RMeq[M102] - RMeq[M220] - RMeq[M202]);
        eq[FM00] = rho * (-u[0] + RMeq[M120] + RMeq[M102]) + eq[FP00];

        eq[F0P0] =
            0.5 * rho * (u[1] + RMeq[M020] - RMeq[M210] - RMeq[M012] - RMeq[M220] - RMeq[M022]);
        eq[F0M0] = rho * (-u[1] + RMeq[M210] + RMeq[M012]) + eq[F0P0];

        eq[F00P] =
            0.5 * rho * (u[2] + RMeq[M002] - RMeq[M201] - RMeq[M021] - RMeq[M202] - RMeq[M022]);
        eq[F00M] = rho * (-u[2] + RMeq[M201] + RMeq[M021]) + eq[F00P];

        eq[FPP0] = 0.25 * rho * (RMeq[M110] + RMeq[M210] + RMeq[M120] + RMeq[M220]);
        eq[FMP0] = 0.5 * rho * (-RMeq[M110] - RMeq[M120]) + eq[FPP0];
        eq[FPM0] = 0.5 * rho * (-RMeq[M110] - RMeq[M210]) + eq[FPP0];
        eq[FMM0] = 0.5 * rho * (-RMeq[M210] - RMeq[M120]) + eq[FPP0];

        eq[FP0P] = 0.25 * rho * (RMeq[M101] + RMeq[M201] + RMeq[M102] + RMeq[M202]);
        eq[FM0P] = 0.5 * rho * (-RMeq[M101] - RMeq[M102]) + eq[FP0P];
        eq[FP0M] = 0.5 * rho * (-RMeq[M101] - RMeq[M201]) + eq[FP0P];
        eq[FM0M] = 0.5 * rho * (-RMeq[M201] - RMeq[M102]) + eq[FP0P];

        eq[F0PP] = 0.25 * rho * (RMeq[M011] + RMeq[M021] + RMeq[M012] + RMeq[M022]);
        eq[F0MP] = 0.5 * rho * (-RMeq[M011] - RMeq[M012]) + eq[F0PP];
        eq[F0PM] = 0.5 * rho * (-RMeq[M011] - RMeq[M021]) + eq[F0PP];
        eq[F0MM] = 0.5 * rho * (-RMeq[M021] - RMeq[M012]) + eq[F0PP];
    };

    static void RMcollide(
        Array<T, D::q> &cell, T rho, Array<T, D::d> const &u,
        Array<T, D::q> const &RM,    // Raw moments
        Array<T, D::q> const &RMeq,  // Equilibrium moments (raw)
        Array<T, D::numRelaxationTimes> const &omega)
    {
        T omega1 = omega[0];
        T omega2 = omega[1];
        T omega3 = omega[2];
        T omega4 = omega[3];
        T omegaBulk = omega[4];
        T omegaPlus = (omegaBulk + 2. * omega1) / 3.;  // Notation used by Fei
        T omegaMinus = (omegaBulk - omega1) / 3.;      // Notation used by Fei

        // Post-collision moments.
        Array<T, D::q> RMcoll;

        // Order 2 (non-diagonal collision so that we can easily modify the bulk viscosity)
        RMcoll[M200] = RM[M200] - omegaPlus * (RM[M200] - RMeq[M200])
                       - omegaMinus * (RM[M020] - RMeq[M020])
                       - omegaMinus * (RM[M002] - RMeq[M002]);
        RMcoll[M020] = RM[M020] - omegaMinus * (RM[M200] - RMeq[M200])
                       - omegaPlus * (RM[M020] - RMeq[M020]) - omegaMinus * (RM[M002] - RMeq[M002]);
        RMcoll[M002] = RM[M002] - omegaMinus * (RM[M200] - RMeq[M200])
                       - omegaMinus * (RM[M020] - RMeq[M020]) - omegaPlus * (RM[M002] - RMeq[M002]);

        RMcoll[M110] = (1. - omega2) * RM[M110] + omega2 * RMeq[M110];
        RMcoll[M101] = (1. - omega2) * RM[M101] + omega2 * RMeq[M101];
        RMcoll[M011] = (1. - omega2) * RM[M011] + omega2 * RMeq[M011];

        // Order 3
        RMcoll[M210] = (1. - omega3) * RM[M210] + omega3 * RMeq[M210];
        RMcoll[M201] = (1. - omega3) * RM[M201] + omega3 * RMeq[M201];
        RMcoll[M021] = (1. - omega3) * RM[M021] + omega3 * RMeq[M021];
        RMcoll[M120] = (1. - omega3) * RM[M120] + omega3 * RMeq[M120];
        RMcoll[M102] = (1. - omega3) * RM[M102] + omega3 * RMeq[M102];
        RMcoll[M012] = (1. - omega3) * RM[M012] + omega3 * RMeq[M012];

        // Order 4
        RMcoll[M220] = (1. - omega4) * RM[M220] + omega4 * RMeq[M220];
        RMcoll[M202] = (1. - omega4) * RM[M202] + omega4 * RMeq[M202];
        RMcoll[M022] = (1. - omega4) * RM[M022] + omega4 * RMeq[M022];

        // Compute post collision populations from RM
        // Optimization based on symmetries between populations and their opposite counterpart
        cell[F000] = rho
                     * (1. - RMcoll[M200] - RMcoll[M020] - RMcoll[M002] + RMcoll[M220]
                        + RMcoll[M202] + RMcoll[M022]);

        cell[FP00] =
            0.5 * rho
            * (u[0] + RMcoll[M200] - RMcoll[M120] - RMcoll[M102] - RMcoll[M220] - RMcoll[M202]);
        cell[FM00] = rho * (-u[0] + RMcoll[M120] + RMcoll[M102]) + cell[FP00];

        cell[F0P0] =
            0.5 * rho
            * (u[1] + RMcoll[M020] - RMcoll[M210] - RMcoll[M012] - RMcoll[M220] - RMcoll[M022]);
        cell[F0M0] = rho * (-u[1] + RMcoll[M210] + RMcoll[M012]) + cell[F0P0];

        cell[F00P] =
            0.5 * rho
            * (u[2] + RMcoll[M002] - RMcoll[M201] - RMcoll[M021] - RMcoll[M202] - RMcoll[M022]);
        cell[F00M] = rho * (-u[2] + RMcoll[M201] + RMcoll[M021]) + cell[F00P];

        cell[FPP0] = 0.25 * rho * (RMcoll[M110] + RMcoll[M210] + RMcoll[M120] + RMcoll[M220]);
        cell[FMP0] = 0.5 * rho * (-RMcoll[M110] - RMcoll[M120]) + cell[FPP0];
        cell[FPM0] = 0.5 * rho * (-RMcoll[M110] - RMcoll[M210]) + cell[FPP0];
        cell[FMM0] = 0.5 * rho * (-RMcoll[M210] - RMcoll[M120]) + cell[FPP0];

        cell[FP0P] = 0.25 * rho * (RMcoll[M101] + RMcoll[M201] + RMcoll[M102] + RMcoll[M202]);
        cell[FM0P] = 0.5 * rho * (-RMcoll[M101] - RMcoll[M102]) + cell[FP0P];
        cell[FP0M] = 0.5 * rho * (-RMcoll[M101] - RMcoll[M201]) + cell[FP0P];
        cell[FM0M] = 0.5 * rho * (-RMcoll[M201] - RMcoll[M102]) + cell[FP0P];

        cell[F0PP] = 0.25 * rho * (RMcoll[M011] + RMcoll[M021] + RMcoll[M012] + RMcoll[M022]);
        cell[F0MP] = 0.5 * rho * (-RMcoll[M011] - RMcoll[M012]) + cell[F0PP];
        cell[F0PM] = 0.5 * rho * (-RMcoll[M011] - RMcoll[M021]) + cell[F0PP];
        cell[F0MM] = 0.5 * rho * (-RMcoll[M021] - RMcoll[M012]) + cell[F0PP];

        for (int i = 0; i < 19; ++i) {
            cell[i] -= D::SkordosFactor() * D::t[i];
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////
    // Hermite Moments Formalism (Equilibrium is computed through raw moments formalism) //
    ///////////////////////////////////////////////////////////////////////////////////////

    /**
     * // General way to compute HMs
     * static void HMcomputeMoments(Array<T,D::q> const& cell, Array<T, D::q>& HM, T& rho) {
     *
     *     Array<T, D::q> f;
     *     for (int i = 0; i<19; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *         HM[i] = 0.;
     *     }
     *
     *     for (int i = 0; i<19; ++i) {
     *         T Hxx = D::c[i][0] * D::c[i][0] - D::cs2;
     *         T Hyy = D::c[i][1] * D::c[i][1] - D::cs2;
     *         T Hzz = D::c[i][2] * D::c[i][2] - D::cs2;
     *
     *         // Order 0
     *         HM[M000] += f[i];
     *         // Order 1
     *         HM[M100] += D::c[i][0] * f[i];
     *         HM[M010] += D::c[i][1] * f[i];
     *         HM[M001] += D::c[i][2] * f[i];
     *         // Order 2
     *         HM[M200] += Hxx * f[i];
     *         HM[M020] += Hyy * f[i];
     *         HM[M002] += Hzz * f[i];
     *         HM[M110] += D::c[i][0] * D::c[i][1] * f[i];
     *         HM[M101] += D::c[i][0] * D::c[i][2] * f[i];
     *         HM[M011] += D::c[i][1] * D::c[i][2] * f[i];
     *         // Order 3
     *         HM[M210] += Hxx * D::c[i][1] * f[i];
     *         HM[M201] += Hxx * D::c[i][2] * f[i];
     *         HM[M021] += Hyy * D::c[i][2] * f[i];
     *         HM[M120] += D::c[i][0] * Hyy * f[i];
     *         HM[M102] += D::c[i][0] * Hzz * f[i];
     *         HM[M012] += D::c[i][1] * Hzz * f[i];
     *         // Order 4
     *         HM[M220] += Hxx * Hyy * f[i];
     *         HM[M202] += Hxx * Hzz * f[i];
     *         HM[M022] += Hyy * Hzz * f[i];
     *     }
     *     rho = HM[M000];
     *     T invRho = 1. / rho;
     *     for (int i = 0; i<19; ++i) {
     *         HM[i] *= invRho;
     *     }
     * }
     *
     * // Optimized way to compute HMs based on Palabos ordering of discrete velocities
     * static void HMcomputeMoments(Array<T,D::q> const& cell, Array<T, D::q>& HM, T& rho) {
     *
     *     Array<T, D::q> f;
     *     for (int i = 0; i<19; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *     }
     *
     *     T a1 = 1./3. ;
     *     T a2 = 2./3. ;
     *
     *     T b1 = 1./9. ;
     *     T b2 = 2./9. ;
     *     T b3 = 4./9. ;
     *
     *     // Order 0
     *     rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[10] + f[11]
     * + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18]; HM[M000] = 1.0; T invRho = 1./rho;
     *
     *     // Order 1
     *     HM[M100] = invRho * ( - f[1] - f[4] - f[5] - f[6] - f[7] + f[10] + f[13] + f[14] + f[15]
     * + f[16] ); HM[M010] = invRho * ( - f[2] - f[4] + f[5] - f[8] - f[9] + f[11] + f[13] - f[14] +
     * f[17] + f[18] ); HM[M001] = invRho * ( - f[3] - f[6] + f[7] - f[8] + f[9] + f[12] + f[15] -
     * f[16] + f[17] - f[18] );
     *     // Order 2
     *     HM[M200] = invRho * ( - a1*f[0] + a2*f[1] - a1*f[2] - a1*f[3] + a2*f[4] + a2*f[5] +
     * a2*f[6] + a2*f[7] - a1*f[8] - a1*f[9] + a2*f[10] - a1*f[11] - a1*f[12] + a2*f[13] + a2*f[14]
     * + a2*f[15] + a2*f[16] - a1*f[17] - a1*f[18] ); HM[M020] = invRho * ( - a1*f[0] - a1*f[1] +
     * a2*f[2] - a1*f[3] + a2*f[4] + a2*f[5] - a1*f[6] - a1*f[7] + a2*f[8] + a2*f[9] - a1*f[10] +
     * a2*f[11] - a1*f[12] + a2*f[13] + a2*f[14] - a1*f[15] - a1*f[16] + a2*f[17] + a2*f[18] );
     *     HM[M002] = invRho * ( - a1*f[0] - a1*f[1] - a1*f[2] + a2*f[3] - a1*f[4] - a1*f[5] +
     * a2*f[6] + a2*f[7] + a2*f[8] + a2*f[9] - a1*f[10] - a1*f[11] + a2*f[12] - a1*f[13] - a1*f[14]
     * + a2*f[15] + a2*f[16] + a2*f[17] + a2*f[18] ); HM[M110] = invRho * ( f[4] - f[5] + f[13] -
     * f[14] ); HM[M101] = invRho * ( f[6] - f[7] + f[15] - f[16] ); HM[M011] = invRho * ( f[8] -
     * f[9] + f[17] - f[18] );
     *     // Order 3
     *     HM[M210] = invRho * ( a1*f[2] - a2*f[4] + a2*f[5] + a1*f[8] + a1*f[9] - a1*f[11] +
     * a2*f[13] - a2*f[14] - a1*f[17] - a1*f[18] ); HM[M201] = invRho * ( a1*f[3] - a2*f[6] +
     * a2*f[7] + a1*f[8] - a1*f[9] - a1*f[12] + a2*f[15] - a2*f[16] - a1*f[17] + a1*f[18] );
     *     HM[M021] = invRho * ( a1*f[3] + a1*f[6] - a1*f[7] - a2*f[8] + a2*f[9] - a1*f[12] -
     * a1*f[15] + a1*f[16] + a2*f[17] - a2*f[18] ); HM[M120] = invRho * ( a1*f[1] - a2*f[4] -
     * a2*f[5] + a1*f[6] + a1*f[7] - a1*f[10] + a2*f[13] + a2*f[14] - a1*f[15] - a1*f[16] );
     *     HM[M102] = invRho * ( a1*f[1] + a1*f[4] + a1*f[5] - a2*f[6] - a2*f[7] - a1*f[10] -
     * a1*f[13] - a1*f[14] + a2*f[15] + a2*f[16] ); HM[M012] = invRho * ( a1*f[2] + a1*f[4] -
     * a1*f[5] - a2*f[8] - a2*f[9] - a1*f[11] - a1*f[13] + a1*f[14] + a2*f[17] + a2*f[18] );
     *     // Order 4
     *     HM[M220] = invRho * ( b1*f[0] - b2*f[1] - b2*f[2] + b1*f[3] + b3*f[4] + b3*f[5] - b2*f[6]
     * - b2*f[7] - b2*f[8] - b2*f[9] - b2*f[10] - b2*f[11] + b1*f[12] + b3*f[13] + b3*f[14] -
     * b2*f[15] - b2*f[16] - b2*f[17] - b2*f[18] ); HM[M202] = invRho * ( b1*f[0] - b2*f[1] +
     * b1*f[2] - b2*f[3] - b2*f[4] - b2*f[5] + b3*f[6] + b3*f[7] - b2*f[8] - b2*f[9] - b2*f[10] +
     * b1*f[11] - b2*f[12] - b2*f[13] - b2*f[14] + b3*f[15] + b3*f[16] - b2*f[17] - b2*f[18] );
     *     HM[M022] = invRho * ( b1*f[0] + b1*f[1] - b2*f[2] - b2*f[3] - b2*f[4] - b2*f[5] - b2*f[6]
     * - b2*f[7] + b3*f[8] + b3*f[9] + b1*f[10] - b2*f[11] - b2*f[12] - b2*f[13] - b2*f[14] -
     * b2*f[15] - b2*f[16] + b3*f[17] + b3*f[18] );
     * };
     */

    // Optimized way to compute HMs based on the general ordering of discrete velocities
    static void HMcomputeMoments(Array<T, D::q> const &cell, Array<T, D::q> &HM, T &rho)
    {
        Array<T, D::q> f;
        for (int i = 0; i < 19; ++i) {
            f[i] = cell[i] + D::SkordosFactor() * D::t[i];
        }
        T X_M1 = f[FM00] + f[FMM0] + f[FMP0] + f[FM0M] + f[FM0P];
        T X_P1 = f[FP00] + f[FPP0] + f[FPM0] + f[FP0P] + f[FP0M];
        T X_0 =
            f[F000] + f[F0M0] + f[F00M] + f[F0MM] + f[F0MP] + f[F0P0] + f[F00P] + f[F0PP] + f[F0PM];
        T Y_M1 = f[F0M0] + f[FMM0] + f[F0MM] + f[F0MP] + f[FPM0];
        T Y_P1 = f[FMP0] + f[F0P0] + f[FPP0] + f[F0PP] + f[F0PM];
        T Z_M1 = f[F00M] + f[FM0M] + f[F0MM] + f[FP0M] + f[F0PM];
        T Z_P1 = f[FM0P] + f[F0MP] + f[F00P] + f[FP0P] + f[F0PP];

        rho = X_M1 + X_P1 + X_0;
        T invRho = 1. / rho;
        T two_invRho = 2. * invRho;

        // Order 0
        HM[M000] = 1.;
        // Order 1
        HM[M100] = invRho * (X_P1 - X_M1);
        HM[M010] = invRho * (Y_P1 - Y_M1);
        HM[M001] = invRho * (Z_P1 - Z_M1);
        // Order 2
        HM[M200] = invRho * (X_M1 + X_P1);
        HM[M020] = invRho * (Y_M1 + Y_P1);
        HM[M002] = invRho * (Z_M1 + Z_P1);
        HM[M110] = invRho * (f[FMM0] - f[FMP0] + f[FPP0] - f[FPM0]);
        HM[M101] = invRho * (f[FM0M] - f[FM0P] + f[FP0P] - f[FP0M]);
        HM[M011] = invRho * (f[F0MM] - f[F0MP] + f[F0PP] - f[F0PM]);
        // Order 3
        HM[M210] = HM[M110] - two_invRho * (f[FMM0] - f[FMP0]);
        HM[M201] = HM[M101] - two_invRho * (f[FM0M] - f[FM0P]);
        HM[M021] = HM[M011] - two_invRho * (f[F0MM] - f[F0MP]);
        HM[M120] = HM[M110] - two_invRho * (f[FMM0] - f[FPM0]);
        HM[M102] = HM[M101] - two_invRho * (f[FM0M] - f[FP0M]);
        HM[M012] = HM[M011] - two_invRho * (f[F0MM] - f[F0PM]);
        // Order 4
        HM[M220] = HM[M110] + two_invRho * (f[FMP0] + f[FPM0]);
        HM[M202] = HM[M101] + two_invRho * (f[FM0P] + f[FP0M]);
        HM[M022] = HM[M011] + two_invRho * (f[F0MP] + f[F0PM]);

        // We come back to Hermite moments
        T cs4 = D::cs2 * D::cs2;
        HM[M200] -= D::cs2;
        HM[M020] -= D::cs2;
        HM[M002] -= D::cs2;

        HM[M210] -= D::cs2 * HM[M010];
        HM[M201] -= D::cs2 * HM[M001];
        HM[M021] -= D::cs2 * HM[M001];
        HM[M120] -= D::cs2 * HM[M100];
        HM[M102] -= D::cs2 * HM[M100];
        HM[M012] -= D::cs2 * HM[M010];

        HM[M220] -= (D::cs2 * (HM[M200] + HM[M020]) + cs4);
        HM[M202] -= (D::cs2 * (HM[M200] + HM[M002]) + cs4);
        HM[M022] -= (D::cs2 * (HM[M020] + HM[M002]) + cs4);
    };

    static void HMcomputeEquilibriumMoments(Array<T, D::d> const &u, Array<T, D::q> &HMeq)
    {
        // Order 0
        HMeq[M000] = 1.;
        // Order 1
        HMeq[M100] = u[0];
        HMeq[M010] = u[1];
        HMeq[M001] = u[2];
        // Order 2
        HMeq[M200] = u[0] * u[0];
        HMeq[M020] = u[1] * u[1];
        HMeq[M002] = u[2] * u[2];
        HMeq[M110] = u[0] * u[1];
        HMeq[M101] = u[0] * u[2];
        HMeq[M011] = u[1] * u[2];
        // Order 3
        HMeq[M210] = HMeq[M200] * u[1];
        HMeq[M201] = HMeq[M200] * u[2];
        HMeq[M021] = HMeq[M020] * u[2];
        HMeq[M120] = HMeq[M020] * u[0];
        HMeq[M102] = HMeq[M002] * u[0];
        HMeq[M012] = HMeq[M002] * u[1];
        // Order 4
        HMeq[M220] = HMeq[M200] * HMeq[M020];
        HMeq[M202] = HMeq[M200] * HMeq[M002];
        HMeq[M022] = HMeq[M020] * HMeq[M002];
    };

    // Equilibrium populations based on 19 moments can be computed using either RM, HM, CM, CHM or
    // Gauss-Hermite formalisms. All formulations are equivalent for the D3Q19 as long as we stick
    // to non-weighted formulas (i.e, no Gauss-Hermite) so we use the fastest one (RMs)
    static void HMcomputeEquilibrium(T rho, Array<T, D::q> const &HMeq, Array<T, D::q> &eq)
    {
        Array<T, D::d> u(HMeq[1], HMeq[2], HMeq[3]);
        Array<T, D::q> RMeq;
        // Order 2
        RMeq[M200] = u[0] * u[0] + D::cs2;
        RMeq[M020] = u[1] * u[1] + D::cs2;
        RMeq[M002] = u[2] * u[2] + D::cs2;
        RMeq[M110] = u[0] * u[1];
        RMeq[M101] = u[0] * u[2];
        RMeq[M011] = u[1] * u[2];
        // Order 3
        RMeq[M210] = RMeq[M200] * u[1];
        RMeq[M201] = RMeq[M200] * u[2];
        RMeq[M021] = RMeq[M020] * u[2];
        RMeq[M120] = RMeq[M020] * u[0];
        RMeq[M102] = RMeq[M002] * u[0];
        RMeq[M012] = RMeq[M002] * u[1];
        // Order 4
        RMeq[M220] = RMeq[M200] * RMeq[M020];
        RMeq[M202] = RMeq[M200] * RMeq[M002];
        RMeq[M022] = RMeq[M020] * RMeq[M002];

        // Optimization based on symmetries between populations and their opposite counterpart
        eq[F000] =
            rho
            * (1. - RMeq[M200] - RMeq[M020] - RMeq[M002] + RMeq[M220] + RMeq[M202] + RMeq[M022]);

        eq[FP00] =
            0.5 * rho * (u[0] + RMeq[M200] - RMeq[M120] - RMeq[M102] - RMeq[M220] - RMeq[M202]);
        eq[FM00] = rho * (-u[0] + RMeq[M120] + RMeq[M102]) + eq[FP00];

        eq[F0P0] =
            0.5 * rho * (u[1] + RMeq[M020] - RMeq[M210] - RMeq[M012] - RMeq[M220] - RMeq[M022]);
        eq[F0M0] = rho * (-u[1] + RMeq[M210] + RMeq[M012]) + eq[F0P0];

        eq[F00P] =
            0.5 * rho * (u[2] + RMeq[M002] - RMeq[M201] - RMeq[M021] - RMeq[M202] - RMeq[M022]);
        eq[F00M] = rho * (-u[2] + RMeq[M201] + RMeq[M021]) + eq[F00P];

        eq[FPP0] = 0.25 * rho * (RMeq[M110] + RMeq[M210] + RMeq[M120] + RMeq[M220]);
        eq[FMP0] = 0.5 * rho * (-RMeq[M110] - RMeq[M120]) + eq[FPP0];
        eq[FPM0] = 0.5 * rho * (-RMeq[M110] - RMeq[M210]) + eq[FPP0];
        eq[FMM0] = 0.5 * rho * (-RMeq[M210] - RMeq[M120]) + eq[FPP0];

        eq[FP0P] = 0.25 * rho * (RMeq[M101] + RMeq[M201] + RMeq[M102] + RMeq[M202]);
        eq[FM0P] = 0.5 * rho * (-RMeq[M101] - RMeq[M102]) + eq[FP0P];
        eq[FP0M] = 0.5 * rho * (-RMeq[M101] - RMeq[M201]) + eq[FP0P];
        eq[FM0M] = 0.5 * rho * (-RMeq[M201] - RMeq[M102]) + eq[FP0P];

        eq[F0PP] = 0.25 * rho * (RMeq[M011] + RMeq[M021] + RMeq[M012] + RMeq[M022]);
        eq[F0MP] = 0.5 * rho * (-RMeq[M011] - RMeq[M012]) + eq[F0PP];
        eq[F0PM] = 0.5 * rho * (-RMeq[M011] - RMeq[M021]) + eq[F0PP];
        eq[F0MM] = 0.5 * rho * (-RMeq[M021] - RMeq[M012]) + eq[F0PP];
    };

    static void HMcollide(
        Array<T, D::q> &cell, T rho, Array<T, D::d> const &u,
        Array<T, D::q> const &HM,    // Hermite moments
        Array<T, D::q> const &HMeq,  // Equilibrium moments (Hermite)
        Array<T, D::numRelaxationTimes> const &omega)
    {
        T omega1 = omega[0];
        T omega2 = omega[1];
        T omega3 = omega[2];
        T omega4 = omega[3];
        T omegaBulk = omega[4];
        T omegaPlus = (omegaBulk + 2. * omega1) / 3.;  // Notation used by Fei
        T omegaMinus = (omegaBulk - omega1) / 3.;      // Notation used by Fei

        T cs4 = D::cs2 * D::cs2;

        // Post-collision moments.
        Array<T, D::q> HMcoll;
        Array<T, D::q> RMcoll;

        // Collision in the Hermite moment space
        // Order 2 (non-diagonal collision so that we can easily modify the bulk viscosity)
        HMcoll[M200] = HM[M200] - omegaPlus * (HM[M200] - HMeq[M200])
                       - omegaMinus * (HM[M020] - HMeq[M020])
                       - omegaMinus * (HM[M002] - HMeq[M002]);
        HMcoll[M020] = HM[M020] - omegaMinus * (HM[M200] - HMeq[M200])
                       - omegaPlus * (HM[M020] - HMeq[M020]) - omegaMinus * (HM[M002] - HMeq[M002]);
        HMcoll[M002] = HM[M002] - omegaMinus * (HM[M200] - HMeq[M200])
                       - omegaMinus * (HM[M020] - HMeq[M020]) - omegaPlus * (HM[M002] - HMeq[M002]);

        HMcoll[M110] = (1. - omega2) * HM[M110] + omega2 * HMeq[M110];
        HMcoll[M101] = (1. - omega2) * HM[M101] + omega2 * HMeq[M101];
        HMcoll[M011] = (1. - omega2) * HM[M011] + omega2 * HMeq[M011];

        // Order 3
        HMcoll[M210] = (1. - omega3) * HM[M210] + omega3 * HMeq[M210];
        HMcoll[M201] = (1. - omega3) * HM[M201] + omega3 * HMeq[M201];
        HMcoll[M021] = (1. - omega3) * HM[M021] + omega3 * HMeq[M021];
        HMcoll[M120] = (1. - omega3) * HM[M120] + omega3 * HMeq[M120];
        HMcoll[M102] = (1. - omega3) * HM[M102] + omega3 * HMeq[M102];
        HMcoll[M012] = (1. - omega3) * HM[M012] + omega3 * HMeq[M012];

        // Order 4
        HMcoll[M220] = (1. - omega4) * HM[M220] + omega4 * HMeq[M220];
        HMcoll[M202] = (1. - omega4) * HM[M202] + omega4 * HMeq[M202];
        HMcoll[M022] = (1. - omega4) * HM[M022] + omega4 * HMeq[M022];

        // Come back to RMcoll using relationships between HMs and RMs
        RMcoll[M200] = HMcoll[M200] + D::cs2;
        RMcoll[M020] = HMcoll[M020] + D::cs2;
        RMcoll[M002] = HMcoll[M002] + D::cs2;

        RMcoll[M110] = HMcoll[M110];
        RMcoll[M101] = HMcoll[M101];
        RMcoll[M011] = HMcoll[M011];

        RMcoll[M210] = HMcoll[M210] + D::cs2 * u[1];
        RMcoll[M201] = HMcoll[M201] + D::cs2 * u[2];
        RMcoll[M021] = HMcoll[M021] + D::cs2 * u[2];
        RMcoll[M120] = HMcoll[M120] + D::cs2 * u[0];
        RMcoll[M102] = HMcoll[M102] + D::cs2 * u[0];
        RMcoll[M012] = HMcoll[M012] + D::cs2 * u[1];

        RMcoll[M220] = HMcoll[M220] + D::cs2 * (HMcoll[M200] + HMcoll[M020]) + cs4;
        RMcoll[M202] = HMcoll[M202] + D::cs2 * (HMcoll[M200] + HMcoll[M002]) + cs4;
        RMcoll[M022] = HMcoll[M022] + D::cs2 * (HMcoll[M020] + HMcoll[M002]) + cs4;

        // Compute post collision populations from RM
        // Optimization based on symmetries between populations and their opposite counterpart
        cell[F000] = rho
                     * (1. - RMcoll[M200] - RMcoll[M020] - RMcoll[M002] + RMcoll[M220]
                        + RMcoll[M202] + RMcoll[M022]);

        cell[FP00] =
            0.5 * rho
            * (u[0] + RMcoll[M200] - RMcoll[M120] - RMcoll[M102] - RMcoll[M220] - RMcoll[M202]);
        cell[FM00] = rho * (-u[0] + RMcoll[M120] + RMcoll[M102]) + cell[FP00];

        cell[F0P0] =
            0.5 * rho
            * (u[1] + RMcoll[M020] - RMcoll[M210] - RMcoll[M012] - RMcoll[M220] - RMcoll[M022]);
        cell[F0M0] = rho * (-u[1] + RMcoll[M210] + RMcoll[M012]) + cell[F0P0];

        cell[F00P] =
            0.5 * rho
            * (u[2] + RMcoll[M002] - RMcoll[M201] - RMcoll[M021] - RMcoll[M202] - RMcoll[M022]);
        cell[F00M] = rho * (-u[2] + RMcoll[M201] + RMcoll[M021]) + cell[F00P];

        cell[FPP0] = 0.25 * rho * (RMcoll[M110] + RMcoll[M210] + RMcoll[M120] + RMcoll[M220]);
        cell[FMP0] = 0.5 * rho * (-RMcoll[M110] - RMcoll[M120]) + cell[FPP0];
        cell[FPM0] = 0.5 * rho * (-RMcoll[M110] - RMcoll[M210]) + cell[FPP0];
        cell[FMM0] = 0.5 * rho * (-RMcoll[M210] - RMcoll[M120]) + cell[FPP0];

        cell[FP0P] = 0.25 * rho * (RMcoll[M101] + RMcoll[M201] + RMcoll[M102] + RMcoll[M202]);
        cell[FM0P] = 0.5 * rho * (-RMcoll[M101] - RMcoll[M102]) + cell[FP0P];
        cell[FP0M] = 0.5 * rho * (-RMcoll[M101] - RMcoll[M201]) + cell[FP0P];
        cell[FM0M] = 0.5 * rho * (-RMcoll[M201] - RMcoll[M102]) + cell[FP0P];

        cell[F0PP] = 0.25 * rho * (RMcoll[M011] + RMcoll[M021] + RMcoll[M012] + RMcoll[M022]);
        cell[F0MP] = 0.5 * rho * (-RMcoll[M011] - RMcoll[M012]) + cell[F0PP];
        cell[F0PM] = 0.5 * rho * (-RMcoll[M011] - RMcoll[M021]) + cell[F0PP];
        cell[F0MM] = 0.5 * rho * (-RMcoll[M021] - RMcoll[M012]) + cell[F0PP];

        for (int i = 0; i < 19; ++i) {
            cell[i] -= D::SkordosFactor() * D::t[i];
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////
    // Central Moments Formalism (Equilibrium is computed through raw moments formalism) //
    ///////////////////////////////////////////////////////////////////////////////////////

    /**
     * General way to compute CMs
     * static void CMcomputeMoments(Array<T,D::q> const& cell, Array<T, D::q>& CM, T& rho,
     * Array<T,D::d>& u) {
     *
     *     Array<T, D::q> f;
     *     for (int i = 0; i<19; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *         CM[i] = 0.;
     *     }
     *
     *     rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[10] + f[11]
     * + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18]; CM[M000] = 1.0; T invRho = 1./rho;
     *
     *     u[0] = invRho * ( - f[1] - f[4] - f[5] - f[6] - f[7] + f[10] + f[13] + f[14] + f[15] +
     * f[16] ); u[1] = invRho * ( - f[2] - f[4] + f[5] - f[8] - f[9] + f[11] + f[13] - f[14] + f[17]
     * + f[18] ); u[2] = invRho * ( - f[3] - f[6] + f[7] - f[8] + f[9] + f[12] + f[15] - f[16] +
     * f[17] - f[18] );
     *
     *     T cMux = 0.;
     *     T cMuy = 0.;
     *     T cMuz = 0.;
     *
     *     for (int i = 0; i<19; ++i) {
     *
     *         cMux = D::c[i][0]- u[0];
     *         cMuy = D::c[i][1]- u[1];
     *         cMuz = D::c[i][2]- u[2];
     *
     *         // // Order 0
     *         // CM[M000] += f[i];
     *
     *         // Order 1
     *         CM[M100] += cMux * f[i];
     *         CM[M010] += cMuy * f[i];
     *         CM[M001] += cMuz * f[i];
     *
     *         // Order 2
     *         CM[M200] += cMux * cMux * f[i];
     *         CM[M020] += cMuy * cMuy * f[i];
     *         CM[M002] += cMuz * cMuz * f[i];
     *         CM[M110] += cMux * cMuy * f[i];
     *         CM[M101] += cMux * cMuz * f[i];
     *         CM[M011] += cMuy * cMuz * f[i];
     *
     *         // Order 3
     *         CM[M210] += cMux * cMux * cMuy * f[i];
     *         CM[M201] += cMux * cMux * cMuz * f[i];
     *         CM[M021] += cMuy * cMuy * cMuz * f[i];
     *         CM[M120] += cMux * cMuy * cMuy * f[i];
     *         CM[M102] += cMux * cMuz * cMuz * f[i];
     *         CM[M012] += cMuy * cMuz * cMuz * f[i];
     *
     *         // Order 4
     *         CM[M220] += cMux * cMux * cMuy * cMuy * f[i];
     *         CM[M202] += cMux * cMux * cMuz * cMuz * f[i];
     *         CM[M022] += cMuy * cMuy * cMuz * cMuz * f[i];
     *     }
     *
     *     for (int i = 1; i<19; ++i) {
     *         CM[i] *= invRho;
     *     }
     * };
     */

    // Optimized way to compute CMs based on the general ordering of discrete velocities
    static void CMcomputeMoments(
        Array<T, D::q> const &cell, Array<T, D::q> &CM, T &rho, Array<T, D::d> &u)
    {
        Array<T, D::q> f;
        for (int i = 0; i < 19; ++i) {
            f[i] = cell[i] + D::SkordosFactor() * D::t[i];
            CM[i] = 0.;
        }
        T X_M1 = f[FM00] + f[FMM0] + f[FMP0] + f[FM0M] + f[FM0P];
        T X_P1 = f[FP00] + f[FPP0] + f[FPM0] + f[FP0P] + f[FP0M];
        T X_0 =
            f[F000] + f[F0M0] + f[F00M] + f[F0MM] + f[F0MP] + f[F0P0] + f[F00P] + f[F0PP] + f[F0PM];
        T Y_M1 = f[F0M0] + f[FMM0] + f[F0MM] + f[F0MP] + f[FPM0];
        T Y_P1 = f[FMP0] + f[F0P0] + f[FPP0] + f[F0PP] + f[F0PM];
        T Z_M1 = f[F00M] + f[FM0M] + f[F0MM] + f[FP0M] + f[F0PM];
        T Z_P1 = f[FM0P] + f[F0MP] + f[F00P] + f[FP0P] + f[F0PP];

        rho = X_M1 + X_P1 + X_0;
        T invRho = 1. / rho;
        T two_invRho = 2. * invRho;

        // Order 0
        CM[M000] = 1.;
        // Order 1
        CM[M100] = 0.;
        CM[M010] = 0.;
        CM[M001] = 0.;
        // Order 2
        CM[M200] = invRho * (X_M1 + X_P1);
        CM[M020] = invRho * (Y_M1 + Y_P1);
        CM[M002] = invRho * (Z_M1 + Z_P1);
        CM[M110] = invRho * (f[FMM0] - f[FMP0] + f[FPP0] - f[FPM0]);
        CM[M101] = invRho * (f[FM0M] - f[FM0P] + f[FP0P] - f[FP0M]);
        CM[M011] = invRho * (f[F0MM] - f[F0MP] + f[F0PP] - f[F0PM]);
        // Order 3
        CM[M210] = CM[M110] - two_invRho * (f[FMM0] - f[FMP0]);
        CM[M201] = CM[M101] - two_invRho * (f[FM0M] - f[FM0P]);
        CM[M021] = CM[M011] - two_invRho * (f[F0MM] - f[F0MP]);
        CM[M120] = CM[M110] - two_invRho * (f[FMM0] - f[FPM0]);
        CM[M102] = CM[M101] - two_invRho * (f[FM0M] - f[FP0M]);
        CM[M012] = CM[M011] - two_invRho * (f[F0MM] - f[F0PM]);
        // Order 4
        CM[M220] = CM[M110] + two_invRho * (f[FMP0] + f[FPM0]);
        CM[M202] = CM[M101] + two_invRho * (f[FM0P] + f[FP0M]);
        CM[M022] = CM[M011] + two_invRho * (f[F0MP] + f[F0PM]);

        // Compute CMs from RMs using binomial formulas
        u[0] = invRho * (X_P1 - X_M1);
        u[1] = invRho * (Y_P1 - Y_M1);
        u[2] = invRho * (Z_P1 - Z_M1);
        T ux2 = u[0] * u[0];
        T uy2 = u[1] * u[1];
        T uz2 = u[2] * u[2];

        CM[M200] -= (ux2);
        CM[M020] -= (uy2);
        CM[M002] -= (uz2);

        CM[M110] -= (u[0] * u[1]);
        CM[M101] -= (u[0] * u[2]);
        CM[M011] -= (u[1] * u[2]);

        CM[M210] -= (u[1] * CM[M200] + 2. * u[0] * CM[M110] + ux2 * u[1]);
        CM[M201] -= (u[2] * CM[M200] + 2. * u[0] * CM[M101] + ux2 * u[2]);
        CM[M021] -= (u[2] * CM[M020] + 2. * u[1] * CM[M011] + uy2 * u[2]);
        CM[M120] -= (u[0] * CM[M020] + 2. * u[1] * CM[M110] + u[0] * uy2);
        CM[M102] -= (u[0] * CM[M002] + 2. * u[2] * CM[M101] + u[0] * uz2);
        CM[M012] -= (u[1] * CM[M002] + 2. * u[2] * CM[M011] + u[1] * uz2);

        CM[M220] -=
            (2. * u[1] * CM[M210] + 2. * u[0] * CM[M120] + uy2 * CM[M200] + ux2 * CM[M020]
             + 4. * u[0] * u[1] * CM[M110] + ux2 * uy2);
        CM[M202] -=
            (2. * u[2] * CM[M201] + 2. * u[0] * CM[M102] + uz2 * CM[M200] + ux2 * CM[M002]
             + 4. * u[0] * u[2] * CM[M101] + ux2 * uz2);
        CM[M022] -=
            (2. * u[2] * CM[M021] + 2. * u[1] * CM[M012] + uz2 * CM[M020] + uy2 * CM[M002]
             + 4. * u[1] * u[2] * CM[M011] + uy2 * uz2);
    };

    static void CMcomputeEquilibriumMoments(Array<T, D::q> &CMeq)
    {
        // Order 0
        CMeq[M000] = 1.;
        // Order 1
        CMeq[M100] = 0.;
        CMeq[M010] = 0.;
        CMeq[M001] = 0.;
        // Order 2
        CMeq[M200] = D::cs2;
        CMeq[M020] = D::cs2;
        CMeq[M002] = D::cs2;
        CMeq[M110] = 0.;
        CMeq[M101] = 0.;
        CMeq[M011] = 0.;
        // Order 3
        CMeq[M210] = 0.;
        CMeq[M201] = 0.;
        CMeq[M021] = 0.;
        CMeq[M120] = 0.;
        CMeq[M102] = 0.;
        CMeq[M012] = 0.;
        // Order 4
        CMeq[M220] = D::cs2 * D::cs2;
        CMeq[M202] = D::cs2 * D::cs2;
        CMeq[M022] = D::cs2 * D::cs2;
    };

    // Equilibrium populations based on 19 moments can be computed using either RM, HM, CM, CHM or
    // Gauss-Hermite formalisms. All formulations are equivalent for the D3Q19 as long as we stick
    // to non-weighted formulas (i.e, no Gauss-Hermite) so we use the fastest one (RMs)
    static void CMcomputeEquilibrium(
        T rho, Array<T, D::d> const &u, Array<T, D::q> const &CMeq, Array<T, D::q> &eq)
    {
        Array<T, D::q> RMeq;
        // Order 2
        RMeq[M200] = u[0] * u[0] + D::cs2;
        RMeq[M020] = u[1] * u[1] + D::cs2;
        RMeq[M002] = u[2] * u[2] + D::cs2;
        RMeq[M110] = u[0] * u[1];
        RMeq[M101] = u[0] * u[2];
        RMeq[M011] = u[1] * u[2];
        // Order 3
        RMeq[M210] = RMeq[M200] * u[1];
        RMeq[M201] = RMeq[M200] * u[2];
        RMeq[M021] = RMeq[M020] * u[2];
        RMeq[M120] = RMeq[M020] * u[0];
        RMeq[M102] = RMeq[M002] * u[0];
        RMeq[M012] = RMeq[M002] * u[1];
        // Order 4
        RMeq[M220] = RMeq[M200] * RMeq[M020];
        RMeq[M202] = RMeq[M200] * RMeq[M002];
        RMeq[M022] = RMeq[M020] * RMeq[M002];

        // Optimization based on symmetries between populations and their opposite counterpart
        eq[F000] =
            rho
            * (1. - RMeq[M200] - RMeq[M020] - RMeq[M002] + RMeq[M220] + RMeq[M202] + RMeq[M022]);

        eq[FP00] =
            0.5 * rho * (u[0] + RMeq[M200] - RMeq[M120] - RMeq[M102] - RMeq[M220] - RMeq[M202]);
        eq[FM00] = rho * (-u[0] + RMeq[M120] + RMeq[M102]) + eq[FP00];

        eq[F0P0] =
            0.5 * rho * (u[1] + RMeq[M020] - RMeq[M210] - RMeq[M012] - RMeq[M220] - RMeq[M022]);
        eq[F0M0] = rho * (-u[1] + RMeq[M210] + RMeq[M012]) + eq[F0P0];

        eq[F00P] =
            0.5 * rho * (u[2] + RMeq[M002] - RMeq[M201] - RMeq[M021] - RMeq[M202] - RMeq[M022]);
        eq[F00M] = rho * (-u[2] + RMeq[M201] + RMeq[M021]) + eq[F00P];

        eq[FPP0] = 0.25 * rho * (RMeq[M110] + RMeq[M210] + RMeq[M120] + RMeq[M220]);
        eq[FMP0] = 0.5 * rho * (-RMeq[M110] - RMeq[M120]) + eq[FPP0];
        eq[FPM0] = 0.5 * rho * (-RMeq[M110] - RMeq[M210]) + eq[FPP0];
        eq[FMM0] = 0.5 * rho * (-RMeq[M210] - RMeq[M120]) + eq[FPP0];

        eq[FP0P] = 0.25 * rho * (RMeq[M101] + RMeq[M201] + RMeq[M102] + RMeq[M202]);
        eq[FM0P] = 0.5 * rho * (-RMeq[M101] - RMeq[M102]) + eq[FP0P];
        eq[FP0M] = 0.5 * rho * (-RMeq[M101] - RMeq[M201]) + eq[FP0P];
        eq[FM0M] = 0.5 * rho * (-RMeq[M201] - RMeq[M102]) + eq[FP0P];

        eq[F0PP] = 0.25 * rho * (RMeq[M011] + RMeq[M021] + RMeq[M012] + RMeq[M022]);
        eq[F0MP] = 0.5 * rho * (-RMeq[M011] - RMeq[M012]) + eq[F0PP];
        eq[F0PM] = 0.5 * rho * (-RMeq[M011] - RMeq[M021]) + eq[F0PP];
        eq[F0MM] = 0.5 * rho * (-RMeq[M021] - RMeq[M012]) + eq[F0PP];
    };

    static void CMcollide(
        Array<T, D::q> &cell, T rho, Array<T, D::d> const &u,
        Array<T, D::q> const &CM,    // Central moments
        Array<T, D::q> const &CMeq,  // Equilibrium moments (central)
        Array<T, D::numRelaxationTimes> const &omega)
    {
        T omega1 = omega[0];
        T omega2 = omega[1];
        T omega3 = omega[2];
        T omega4 = omega[3];
        T omegaBulk = omega[4];
        T omegaPlus = (omegaBulk + 2. * omega1) / 3.;  // Notation used by Fei
        T omegaMinus = (omegaBulk - omega1) / 3.;      // Notation used by Fei

        T ux2 = u[0] * u[0];
        T uy2 = u[1] * u[1];
        T uz2 = u[2] * u[2];

        // Post-collision moments.
        Array<T, D::q> CMcoll;
        Array<T, D::q> RMcoll;

        // Collision in the central moment space
        // Order 2 (non-diagonal collision so that we can easily modify the bulk viscosity)
        CMcoll[M200] = CM[M200] - omegaPlus * (CM[M200] - CMeq[M200])
                       - omegaMinus * (CM[M020] - CMeq[M020])
                       - omegaMinus * (CM[M002] - CMeq[M002]);
        CMcoll[M020] = CM[M020] - omegaMinus * (CM[M200] - CMeq[M200])
                       - omegaPlus * (CM[M020] - CMeq[M020]) - omegaMinus * (CM[M002] - CMeq[M002]);
        CMcoll[M002] = CM[M002] - omegaMinus * (CM[M200] - CMeq[M200])
                       - omegaMinus * (CM[M020] - CMeq[M020]) - omegaPlus * (CM[M002] - CMeq[M002]);

        CMcoll[M110] = (1. - omega2) * CM[M110] + omega2 * CMeq[M110];
        CMcoll[M101] = (1. - omega2) * CM[M101] + omega2 * CMeq[M101];
        CMcoll[M011] = (1. - omega2) * CM[M011] + omega2 * CMeq[M011];

        // Order 3
        CMcoll[M210] = (1. - omega3) * CM[M210] + omega3 * CMeq[M210];
        CMcoll[M201] = (1. - omega3) * CM[M201] + omega3 * CMeq[M201];
        CMcoll[M021] = (1. - omega3) * CM[M021] + omega3 * CMeq[M021];
        CMcoll[M120] = (1. - omega3) * CM[M120] + omega3 * CMeq[M120];
        CMcoll[M102] = (1. - omega3) * CM[M102] + omega3 * CMeq[M102];
        CMcoll[M012] = (1. - omega3) * CM[M012] + omega3 * CMeq[M012];

        // Order 4
        CMcoll[M220] = (1. - omega4) * CM[M220] + omega4 * CMeq[M220];
        CMcoll[M202] = (1. - omega4) * CM[M202] + omega4 * CMeq[M202];
        CMcoll[M022] = (1. - omega4) * CM[M022] + omega4 * CMeq[M022];

        // Come back to RMcoll using binomial formulas
        RMcoll[M200] = CMcoll[M200] + ux2;
        RMcoll[M020] = CMcoll[M020] + uy2;
        RMcoll[M002] = CMcoll[M002] + uz2;

        RMcoll[M110] = CMcoll[M110] + u[0] * u[1];
        RMcoll[M101] = CMcoll[M101] + u[0] * u[2];
        RMcoll[M011] = CMcoll[M011] + u[1] * u[2];

        RMcoll[M210] = CMcoll[M210] + u[1] * CMcoll[M200] + 2. * u[0] * CMcoll[M110] + ux2 * u[1];
        RMcoll[M201] = CMcoll[M201] + u[2] * CMcoll[M200] + 2. * u[0] * CMcoll[M101] + ux2 * u[2];
        RMcoll[M021] = CMcoll[M021] + u[2] * CMcoll[M020] + 2. * u[1] * CMcoll[M011] + uy2 * u[2];
        RMcoll[M120] = CMcoll[M120] + u[0] * CMcoll[M020] + 2. * u[1] * CMcoll[M110] + u[0] * uy2;
        RMcoll[M102] = CMcoll[M102] + u[0] * CMcoll[M002] + 2. * u[2] * CMcoll[M101] + u[0] * uz2;
        RMcoll[M012] = CMcoll[M012] + u[1] * CMcoll[M002] + 2. * u[2] * CMcoll[M011] + u[1] * uz2;

        RMcoll[M220] = CMcoll[M220] + 2. * u[1] * CMcoll[M210] + 2. * u[0] * CMcoll[M120]
                       + uy2 * CMcoll[M200] + ux2 * CMcoll[M020] + 4. * u[0] * u[1] * CMcoll[M110]
                       + ux2 * uy2;
        RMcoll[M202] = CMcoll[M202] + 2. * u[2] * CMcoll[M201] + 2. * u[0] * CMcoll[M102]
                       + uz2 * CMcoll[M200] + ux2 * CMcoll[M002] + 4. * u[0] * u[2] * CMcoll[M101]
                       + ux2 * uz2;
        RMcoll[M022] = CMcoll[M022] + 2. * u[2] * CMcoll[M021] + 2. * u[1] * CMcoll[M012]
                       + uz2 * CMcoll[M020] + uy2 * CMcoll[M002] + 4. * u[1] * u[2] * CMcoll[M011]
                       + uy2 * uz2;

        // Compute post collision populations from RM
        // Optimization based on symmetries between populations and their opposite counterpart
        cell[F000] = rho
                     * (1. - RMcoll[M200] - RMcoll[M020] - RMcoll[M002] + RMcoll[M220]
                        + RMcoll[M202] + RMcoll[M022]);

        cell[FP00] =
            0.5 * rho
            * (u[0] + RMcoll[M200] - RMcoll[M120] - RMcoll[M102] - RMcoll[M220] - RMcoll[M202]);
        cell[FM00] = rho * (-u[0] + RMcoll[M120] + RMcoll[M102]) + cell[FP00];

        cell[F0P0] =
            0.5 * rho
            * (u[1] + RMcoll[M020] - RMcoll[M210] - RMcoll[M012] - RMcoll[M220] - RMcoll[M022]);
        cell[F0M0] = rho * (-u[1] + RMcoll[M210] + RMcoll[M012]) + cell[F0P0];

        cell[F00P] =
            0.5 * rho
            * (u[2] + RMcoll[M002] - RMcoll[M201] - RMcoll[M021] - RMcoll[M202] - RMcoll[M022]);
        cell[F00M] = rho * (-u[2] + RMcoll[M201] + RMcoll[M021]) + cell[F00P];

        cell[FPP0] = 0.25 * rho * (RMcoll[M110] + RMcoll[M210] + RMcoll[M120] + RMcoll[M220]);
        cell[FMP0] = 0.5 * rho * (-RMcoll[M110] - RMcoll[M120]) + cell[FPP0];
        cell[FPM0] = 0.5 * rho * (-RMcoll[M110] - RMcoll[M210]) + cell[FPP0];
        cell[FMM0] = 0.5 * rho * (-RMcoll[M210] - RMcoll[M120]) + cell[FPP0];

        cell[FP0P] = 0.25 * rho * (RMcoll[M101] + RMcoll[M201] + RMcoll[M102] + RMcoll[M202]);
        cell[FM0P] = 0.5 * rho * (-RMcoll[M101] - RMcoll[M102]) + cell[FP0P];
        cell[FP0M] = 0.5 * rho * (-RMcoll[M101] - RMcoll[M201]) + cell[FP0P];
        cell[FM0M] = 0.5 * rho * (-RMcoll[M201] - RMcoll[M102]) + cell[FP0P];

        cell[F0PP] = 0.25 * rho * (RMcoll[M011] + RMcoll[M021] + RMcoll[M012] + RMcoll[M022]);
        cell[F0MP] = 0.5 * rho * (-RMcoll[M011] - RMcoll[M012]) + cell[F0PP];
        cell[F0PM] = 0.5 * rho * (-RMcoll[M011] - RMcoll[M021]) + cell[F0PP];
        cell[F0MM] = 0.5 * rho * (-RMcoll[M021] - RMcoll[M012]) + cell[F0PP];

        for (int i = 0; i < 19; ++i) {
            cell[i] -= D::SkordosFactor() * D::t[i];
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Central Hermite Moments Formalism (Equilibrium is computed through raw moments formalism) //
    ///////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * // General way to compute CHMs
     * static void CHMcomputeMoments(Array<T,D::q> const& cell, Array<T, D::q>& CHM, T& rho,
     * Array<T,D::d>& u) {
     *
     *     Array<T, D::q> f;
     *     for (int i = 0; i<19; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *         CHM[i] = 0.;
     *     }
     *
     *     rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[10] + f[11]
     * + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18]; CHM[M000] = 1.0; T invRho = 1./rho;
     *     u[0] = invRho * ( - f[1] - f[4] - f[5] - f[6] - f[7] + f[10] + f[13] + f[14] + f[15] +
     * f[16] ); u[1] = invRho * ( - f[2] - f[4] + f[5] - f[8] - f[9] + f[11] + f[13] - f[14] + f[17]
     * + f[18] ); u[2] = invRho * ( - f[3] - f[6] + f[7] - f[8] + f[9] + f[12] + f[15] - f[16] +
     * f[17] - f[18] ); T cMux = 0.; T cMuy = 0.; T cMuz = 0.; T Hxx = 0.; T Hyy = 0.; T Hzz = 0.;
     *
     *     for (int i = 0; i<19; ++i) {
     *
     *         cMux = D::c[i][0]- u[0];
     *         cMuy = D::c[i][1]- u[1];
     *         cMuz = D::c[i][2]- u[2];
     *
     *         Hxx = cMux * cMux - D::cs2;
     *         Hyy = cMuy * cMuy - D::cs2;
     *         Hzz = cMuz * cMuz - D::cs2;
     *
     *         // // Order 0
     *         // CHM[M000] += f[i];
     *
     *         // Order 1
     *         CHM[M100] += cMux * f[i];
     *         CHM[M010] += cMuy * f[i];
     *         CHM[M001] += cMuz * f[i];
     *
     *         // Order 2
     *         CHM[M200] += Hxx * f[i];
     *         CHM[M020] += Hyy * f[i];
     *         CHM[M002] += Hzz * f[i];
     *         CHM[M110] += cMux * cMuy * f[i];
     *         CHM[M101] += cMux * cMuz * f[i];
     *         CHM[M011] += cMuy * cMuz * f[i];
     *
     *         // Order 3
     *         CHM[M210] += Hxx * cMuy * f[i];
     *         CHM[M201] += Hxx * cMuz * f[i];
     *         CHM[M021] += Hyy * cMuz * f[i];
     *         CHM[M120] += cMux * Hyy * f[i];
     *         CHM[M102] += cMux * Hzz * f[i];
     *         CHM[M012] += cMuy * Hzz * f[i];
     *
     *         // Order 4
     *         CHM[M220] += Hxx * Hyy * f[i];
     *         CHM[M202] += Hxx * Hzz * f[i];
     *         CHM[M022] += Hyy * Hzz * f[i];
     *     }
     *
     *     for (int i = 1; i<19; ++i) {
     *         CHM[i] *= invRho;
     *     }
     * };
     */

    // Optimized way to compute CHMs based on the general ordering of discrete velocities
    static void CHMcomputeMoments(
        Array<T, D::q> const &cell, Array<T, D::q> &CHM, T &rho, Array<T, D::d> &u)
    {
        Array<T, D::q> f;
        for (int i = 0; i < 19; ++i) {
            f[i] = cell[i] + D::SkordosFactor() * D::t[i];
            CHM[i] = 0.;
        }
        T X_M1 = f[FM00] + f[FMM0] + f[FMP0] + f[FM0M] + f[FM0P];
        T X_P1 = f[FP00] + f[FPP0] + f[FPM0] + f[FP0P] + f[FP0M];
        T X_0 =
            f[F000] + f[F0M0] + f[F00M] + f[F0MM] + f[F0MP] + f[F0P0] + f[F00P] + f[F0PP] + f[F0PM];
        T Y_M1 = f[F0M0] + f[FMM0] + f[F0MM] + f[F0MP] + f[FPM0];
        T Y_P1 = f[FMP0] + f[F0P0] + f[FPP0] + f[F0PP] + f[F0PM];
        T Z_M1 = f[F00M] + f[FM0M] + f[F0MM] + f[FP0M] + f[F0PM];
        T Z_P1 = f[FM0P] + f[F0MP] + f[F00P] + f[FP0P] + f[F0PP];

        rho = X_M1 + X_P1 + X_0;
        T invRho = 1. / rho;
        T two_invRho = 2. * invRho;

        // Order 0
        CHM[M000] = 1.;
        // Order 1
        CHM[M100] = 0.;
        CHM[M010] = 0.;
        CHM[M001] = 0.;
        // Order 2
        CHM[M200] = invRho * (X_M1 + X_P1);
        CHM[M020] = invRho * (Y_M1 + Y_P1);
        CHM[M002] = invRho * (Z_M1 + Z_P1);
        CHM[M110] = invRho * (f[FMM0] - f[FMP0] + f[FPP0] - f[FPM0]);
        CHM[M101] = invRho * (f[FM0M] - f[FM0P] + f[FP0P] - f[FP0M]);
        CHM[M011] = invRho * (f[F0MM] - f[F0MP] + f[F0PP] - f[F0PM]);
        // Order 3
        CHM[M210] = CHM[M110] - two_invRho * (f[FMM0] - f[FMP0]);
        CHM[M201] = CHM[M101] - two_invRho * (f[FM0M] - f[FM0P]);
        CHM[M021] = CHM[M011] - two_invRho * (f[F0MM] - f[F0MP]);
        CHM[M120] = CHM[M110] - two_invRho * (f[FMM0] - f[FPM0]);
        CHM[M102] = CHM[M101] - two_invRho * (f[FM0M] - f[FP0M]);
        CHM[M012] = CHM[M011] - two_invRho * (f[F0MM] - f[F0PM]);
        // Order 4
        CHM[M220] = CHM[M110] + two_invRho * (f[FMP0] + f[FPM0]);
        CHM[M202] = CHM[M101] + two_invRho * (f[FM0P] + f[FP0M]);
        CHM[M022] = CHM[M011] + two_invRho * (f[F0MP] + f[F0PM]);

        // Compute CMs from RMs using binomial formulas
        u[0] = invRho * (X_P1 - X_M1);
        u[1] = invRho * (Y_P1 - Y_M1);
        u[2] = invRho * (Z_P1 - Z_M1);
        T ux2 = u[0] * u[0];
        T uy2 = u[1] * u[1];
        T uz2 = u[2] * u[2];

        CHM[M200] -= (ux2);
        CHM[M020] -= (uy2);
        CHM[M002] -= (uz2);

        CHM[M110] -= (u[0] * u[1]);
        CHM[M101] -= (u[0] * u[2]);
        CHM[M011] -= (u[1] * u[2]);

        CHM[M210] -= (u[1] * CHM[M200] + 2. * u[0] * CHM[M110] + ux2 * u[1]);
        CHM[M201] -= (u[2] * CHM[M200] + 2. * u[0] * CHM[M101] + ux2 * u[2]);
        CHM[M021] -= (u[2] * CHM[M020] + 2. * u[1] * CHM[M011] + uy2 * u[2]);
        CHM[M120] -= (u[0] * CHM[M020] + 2. * u[1] * CHM[M110] + u[0] * uy2);
        CHM[M102] -= (u[0] * CHM[M002] + 2. * u[2] * CHM[M101] + u[0] * uz2);
        CHM[M012] -= (u[1] * CHM[M002] + 2. * u[2] * CHM[M011] + u[1] * uz2);

        CHM[M220] -=
            (2. * u[1] * CHM[M210] + 2. * u[0] * CHM[M120] + uy2 * CHM[M200] + ux2 * CHM[M020]
             + 4. * u[0] * u[1] * CHM[M110] + ux2 * uy2);
        CHM[M202] -=
            (2. * u[2] * CHM[M201] + 2. * u[0] * CHM[M102] + uz2 * CHM[M200] + ux2 * CHM[M002]
             + 4. * u[0] * u[2] * CHM[M101] + ux2 * uz2);
        CHM[M022] -=
            (2. * u[2] * CHM[M021] + 2. * u[1] * CHM[M012] + uz2 * CHM[M020] + uy2 * CHM[M002]
             + 4. * u[1] * u[2] * CHM[M011] + uy2 * uz2);

        // Compute CHMs from CMs
        T cs4 = D::cs2 * D::cs2;

        CHM[M200] -= D::cs2;
        CHM[M020] -= D::cs2;
        CHM[M002] -= D::cs2;

        CHM[M220] -= (D::cs2 * (CHM[M200] + CHM[M020]) + cs4);
        CHM[M202] -= (D::cs2 * (CHM[M200] + CHM[M002]) + cs4);
        CHM[M022] -= (D::cs2 * (CHM[M020] + CHM[M002]) + cs4);
    };

    static void CHMcomputeEquilibriumMoments(Array<T, D::q> &CHMeq)
    {
        // Order 0
        CHMeq[M000] = 1.;
        // Order 1
        CHMeq[M100] = 0.;
        CHMeq[M010] = 0.;
        CHMeq[M001] = 0.;
        // Order 2
        CHMeq[M200] = 0.;
        CHMeq[M020] = 0.;
        CHMeq[M002] = 0.;
        CHMeq[M110] = 0.;
        CHMeq[M101] = 0.;
        CHMeq[M011] = 0.;
        // Order 3
        CHMeq[M210] = 0.;
        CHMeq[M201] = 0.;
        CHMeq[M021] = 0.;
        CHMeq[M120] = 0.;
        CHMeq[M102] = 0.;
        CHMeq[M012] = 0.;
        // Order 4
        CHMeq[M220] = 0.;
        CHMeq[M202] = 0.;
        CHMeq[M022] = 0.;
    };

    // Equilibrium populations based on 19 moments can be computed using either RM, HM, CM, CHM or
    // Gauss-Hermite formalisms. All formulations are equivalent for the D3Q19 as long as we stick
    // to non-weighted formulas (i.e, no Gauss-Hermite) so we use the fastest one (RMs)
    static void CHMcomputeEquilibrium(
        T rho, Array<T, D::d> const &u, Array<T, D::q> const &CHMeq, Array<T, D::q> &eq)
    {
        Array<T, D::q> RMeq;
        // Order 2
        RMeq[M200] = u[0] * u[0] + D::cs2;
        RMeq[M020] = u[1] * u[1] + D::cs2;
        RMeq[M002] = u[2] * u[2] + D::cs2;
        RMeq[M110] = u[0] * u[1];
        RMeq[M101] = u[0] * u[2];
        RMeq[M011] = u[1] * u[2];
        // Order 3
        RMeq[M210] = RMeq[M200] * u[1];
        RMeq[M201] = RMeq[M200] * u[2];
        RMeq[M021] = RMeq[M020] * u[2];
        RMeq[M120] = RMeq[M020] * u[0];
        RMeq[M102] = RMeq[M002] * u[0];
        RMeq[M012] = RMeq[M002] * u[1];
        // Order 4
        RMeq[M220] = RMeq[M200] * RMeq[M020];
        RMeq[M202] = RMeq[M200] * RMeq[M002];
        RMeq[M022] = RMeq[M020] * RMeq[M002];

        // Optimization based on symmetries between populations and their opposite counterpart
        eq[F000] =
            rho
            * (1. - RMeq[M200] - RMeq[M020] - RMeq[M002] + RMeq[M220] + RMeq[M202] + RMeq[M022]);

        eq[FP00] =
            0.5 * rho * (u[0] + RMeq[M200] - RMeq[M120] - RMeq[M102] - RMeq[M220] - RMeq[M202]);
        eq[FM00] = rho * (-u[0] + RMeq[M120] + RMeq[M102]) + eq[FP00];

        eq[F0P0] =
            0.5 * rho * (u[1] + RMeq[M020] - RMeq[M210] - RMeq[M012] - RMeq[M220] - RMeq[M022]);
        eq[F0M0] = rho * (-u[1] + RMeq[M210] + RMeq[M012]) + eq[F0P0];

        eq[F00P] =
            0.5 * rho * (u[2] + RMeq[M002] - RMeq[M201] - RMeq[M021] - RMeq[M202] - RMeq[M022]);
        eq[F00M] = rho * (-u[2] + RMeq[M201] + RMeq[M021]) + eq[F00P];

        eq[FPP0] = 0.25 * rho * (RMeq[M110] + RMeq[M210] + RMeq[M120] + RMeq[M220]);
        eq[FMP0] = 0.5 * rho * (-RMeq[M110] - RMeq[M120]) + eq[FPP0];
        eq[FPM0] = 0.5 * rho * (-RMeq[M110] - RMeq[M210]) + eq[FPP0];
        eq[FMM0] = 0.5 * rho * (-RMeq[M210] - RMeq[M120]) + eq[FPP0];

        eq[FP0P] = 0.25 * rho * (RMeq[M101] + RMeq[M201] + RMeq[M102] + RMeq[M202]);
        eq[FM0P] = 0.5 * rho * (-RMeq[M101] - RMeq[M102]) + eq[FP0P];
        eq[FP0M] = 0.5 * rho * (-RMeq[M101] - RMeq[M201]) + eq[FP0P];
        eq[FM0M] = 0.5 * rho * (-RMeq[M201] - RMeq[M102]) + eq[FP0P];

        eq[F0PP] = 0.25 * rho * (RMeq[M011] + RMeq[M021] + RMeq[M012] + RMeq[M022]);
        eq[F0MP] = 0.5 * rho * (-RMeq[M011] - RMeq[M012]) + eq[F0PP];
        eq[F0PM] = 0.5 * rho * (-RMeq[M011] - RMeq[M021]) + eq[F0PP];
        eq[F0MM] = 0.5 * rho * (-RMeq[M021] - RMeq[M012]) + eq[F0PP];
    };

    static void CHMcollide(
        Array<T, D::q> &cell, T rho, Array<T, D::d> const &u,
        Array<T, D::q> const &CHM,    // Central Hermite moments
        Array<T, D::q> const &CHMeq,  // Equilibrium moments (central Hermite)
        Array<T, D::numRelaxationTimes> const &omega)
    {
        T omega1 = omega[0];
        T omega2 = omega[1];
        T omega3 = omega[2];
        T omega4 = omega[3];
        T omegaBulk = omega[4];
        T omegaPlus = (omegaBulk + 2. * omega1) / 3.;  // Notation used by Fei
        T omegaMinus = (omegaBulk - omega1) / 3.;      // Notation used by Fei

        T ux2 = u[0] * u[0];
        T uy2 = u[1] * u[1];
        T uz2 = u[2] * u[2];
        T cs4 = D::cs2 * D::cs2;

        // Post-collision moments.
        Array<T, D::q> CHMcoll;
        Array<T, D::q> HMcoll;
        Array<T, D::q> RMcoll;

        // Collision in the central Hermite moment space
        // Order 2 (non-diagonal collision so that we can easily modify the bulk viscosity)
        CHMcoll[M200] = CHM[M200] - omegaPlus * (CHM[M200] - CHMeq[M200])
                        - omegaMinus * (CHM[M020] - CHMeq[M020])
                        - omegaMinus * (CHM[M002] - CHMeq[M002]);
        CHMcoll[M020] = CHM[M020] - omegaMinus * (CHM[M200] - CHMeq[M200])
                        - omegaPlus * (CHM[M020] - CHMeq[M020])
                        - omegaMinus * (CHM[M002] - CHMeq[M002]);
        CHMcoll[M002] = CHM[M002] - omegaMinus * (CHM[M200] - CHMeq[M200])
                        - omegaMinus * (CHM[M020] - CHMeq[M020])
                        - omegaPlus * (CHM[M002] - CHMeq[M002]);

        CHMcoll[M110] = (1. - omega2) * CHM[M110] + omega2 * CHMeq[M110];
        CHMcoll[M101] = (1. - omega2) * CHM[M101] + omega2 * CHMeq[M101];
        CHMcoll[M011] = (1. - omega2) * CHM[M011] + omega2 * CHMeq[M011];

        // Order 3
        CHMcoll[M210] = (1. - omega3) * CHM[M210] + omega3 * CHMeq[M210];
        CHMcoll[M201] = (1. - omega3) * CHM[M201] + omega3 * CHMeq[M201];
        CHMcoll[M021] = (1. - omega3) * CHM[M021] + omega3 * CHMeq[M021];
        CHMcoll[M120] = (1. - omega3) * CHM[M120] + omega3 * CHMeq[M120];
        CHMcoll[M102] = (1. - omega3) * CHM[M102] + omega3 * CHMeq[M102];
        CHMcoll[M012] = (1. - omega3) * CHM[M012] + omega3 * CHMeq[M012];

        // Order 4
        CHMcoll[M220] = (1. - omega4) * CHM[M220] + omega4 * CHMeq[M220];
        CHMcoll[M202] = (1. - omega4) * CHM[M202] + omega4 * CHMeq[M202];
        CHMcoll[M022] = (1. - omega4) * CHM[M022] + omega4 * CHMeq[M022];

        // Come back to HMcoll using relationships between CHMs and HMs
        HMcoll[M200] = CHMcoll[M200] + ux2;
        HMcoll[M020] = CHMcoll[M020] + uy2;
        HMcoll[M002] = CHMcoll[M002] + uz2;

        HMcoll[M110] = CHMcoll[M110] + u[0] * u[1];
        HMcoll[M101] = CHMcoll[M101] + u[0] * u[2];
        HMcoll[M011] = CHMcoll[M011] + u[1] * u[2];

        HMcoll[M210] =
            CHMcoll[M210] + u[1] * CHMcoll[M200] + 2. * u[0] * CHMcoll[M110] + ux2 * u[1];
        HMcoll[M201] =
            CHMcoll[M201] + u[2] * CHMcoll[M200] + 2. * u[0] * CHMcoll[M101] + ux2 * u[2];
        HMcoll[M021] =
            CHMcoll[M021] + u[2] * CHMcoll[M020] + 2. * u[1] * CHMcoll[M011] + uy2 * u[2];
        HMcoll[M120] =
            CHMcoll[M120] + u[0] * CHMcoll[M020] + 2. * u[1] * CHMcoll[M110] + u[0] * uy2;
        HMcoll[M102] =
            CHMcoll[M102] + u[0] * CHMcoll[M002] + 2. * u[2] * CHMcoll[M101] + u[0] * uz2;
        HMcoll[M012] =
            CHMcoll[M012] + u[1] * CHMcoll[M002] + 2. * u[2] * CHMcoll[M011] + u[1] * uz2;

        HMcoll[M220] = CHMcoll[M220] + 2. * u[1] * CHMcoll[M210] + 2. * u[0] * CHMcoll[M120]
                       + uy2 * CHMcoll[M200] + ux2 * CHMcoll[M020]
                       + 4. * u[0] * u[1] * CHMcoll[M110] + ux2 * uy2;
        HMcoll[M202] = CHMcoll[M202] + 2. * u[2] * CHMcoll[M201] + 2. * u[0] * CHMcoll[M102]
                       + uz2 * CHMcoll[M200] + ux2 * CHMcoll[M002]
                       + 4. * u[0] * u[2] * CHMcoll[M101] + ux2 * uz2;
        HMcoll[M022] = CHMcoll[M022] + 2. * u[2] * CHMcoll[M021] + 2. * u[1] * CHMcoll[M012]
                       + uz2 * CHMcoll[M020] + uy2 * CHMcoll[M002]
                       + 4. * u[1] * u[2] * CHMcoll[M011] + uy2 * uz2;

        // Come back to RMcoll using relationships between HMs and RMs
        RMcoll[M200] = HMcoll[M200] + D::cs2;
        RMcoll[M020] = HMcoll[M020] + D::cs2;
        RMcoll[M002] = HMcoll[M002] + D::cs2;

        RMcoll[M110] = HMcoll[M110];
        RMcoll[M101] = HMcoll[M101];
        RMcoll[M011] = HMcoll[M011];

        RMcoll[M210] = HMcoll[M210] + D::cs2 * u[1];
        RMcoll[M201] = HMcoll[M201] + D::cs2 * u[2];
        RMcoll[M021] = HMcoll[M021] + D::cs2 * u[2];
        RMcoll[M120] = HMcoll[M120] + D::cs2 * u[0];
        RMcoll[M102] = HMcoll[M102] + D::cs2 * u[0];
        RMcoll[M012] = HMcoll[M012] + D::cs2 * u[1];

        RMcoll[M220] = HMcoll[M220] + D::cs2 * (HMcoll[M200] + HMcoll[M020]) + cs4;
        RMcoll[M202] = HMcoll[M202] + D::cs2 * (HMcoll[M200] + HMcoll[M002]) + cs4;
        RMcoll[M022] = HMcoll[M022] + D::cs2 * (HMcoll[M020] + HMcoll[M002]) + cs4;

        // Compute post collision populations from RM
        // Optimization based on symmetries between populations and their opposite counterpart
        cell[F000] = rho
                     * (1. - RMcoll[M200] - RMcoll[M020] - RMcoll[M002] + RMcoll[M220]
                        + RMcoll[M202] + RMcoll[M022]);

        cell[FP00] =
            0.5 * rho
            * (u[0] + RMcoll[M200] - RMcoll[M120] - RMcoll[M102] - RMcoll[M220] - RMcoll[M202]);
        cell[FM00] = rho * (-u[0] + RMcoll[M120] + RMcoll[M102]) + cell[FP00];

        cell[F0P0] =
            0.5 * rho
            * (u[1] + RMcoll[M020] - RMcoll[M210] - RMcoll[M012] - RMcoll[M220] - RMcoll[M022]);
        cell[F0M0] = rho * (-u[1] + RMcoll[M210] + RMcoll[M012]) + cell[F0P0];

        cell[F00P] =
            0.5 * rho
            * (u[2] + RMcoll[M002] - RMcoll[M201] - RMcoll[M021] - RMcoll[M202] - RMcoll[M022]);
        cell[F00M] = rho * (-u[2] + RMcoll[M201] + RMcoll[M021]) + cell[F00P];

        cell[FPP0] = 0.25 * rho * (RMcoll[M110] + RMcoll[M210] + RMcoll[M120] + RMcoll[M220]);
        cell[FMP0] = 0.5 * rho * (-RMcoll[M110] - RMcoll[M120]) + cell[FPP0];
        cell[FPM0] = 0.5 * rho * (-RMcoll[M110] - RMcoll[M210]) + cell[FPP0];
        cell[FMM0] = 0.5 * rho * (-RMcoll[M210] - RMcoll[M120]) + cell[FPP0];

        cell[FP0P] = 0.25 * rho * (RMcoll[M101] + RMcoll[M201] + RMcoll[M102] + RMcoll[M202]);
        cell[FM0P] = 0.5 * rho * (-RMcoll[M101] - RMcoll[M102]) + cell[FP0P];
        cell[FP0M] = 0.5 * rho * (-RMcoll[M101] - RMcoll[M201]) + cell[FP0P];
        cell[FM0M] = 0.5 * rho * (-RMcoll[M201] - RMcoll[M102]) + cell[FP0P];

        cell[F0PP] = 0.25 * rho * (RMcoll[M011] + RMcoll[M021] + RMcoll[M012] + RMcoll[M022]);
        cell[F0MP] = 0.5 * rho * (-RMcoll[M011] - RMcoll[M012]) + cell[F0PP];
        cell[F0PM] = 0.5 * rho * (-RMcoll[M011] - RMcoll[M021]) + cell[F0PP];
        cell[F0MM] = 0.5 * rho * (-RMcoll[M021] - RMcoll[M012]) + cell[F0PP];

        for (int i = 0; i < 19; ++i) {
            cell[i] -= D::SkordosFactor() * D::t[i];
        }
    };

    ////////////////////////////////////////////////////////////////////////////////
    // Cumulant Formalism (Equilibrium is computed through raw moments formalism) //
    ////////////////////////////////////////////////////////////////////////////////

    /**
     * // General way to compute Ks
     * static void KcomputeMoments(Array<T,D::q> const& cell, Array<T, D::q>& K, T& rho,
     * Array<T,D::d>& u) {
     *
     *     Array<T, D::q> f;
     *     Array<T,D::q> CM;
     *     for (int i = 0; i<19; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *         CM[i] = 0.;
     *     }
     *     rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[10] + f[11]
     * + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18]; CM[M000] = 1.0; T invRho = 1./rho;
     *     u[0] = invRho * ( - f[1] - f[4] - f[5] - f[6] - f[7] + f[10] + f[13] + f[14] + f[15] +
     * f[16] ); u[1] = invRho * ( - f[2] - f[4] + f[5] - f[8] - f[9] + f[11] + f[13] - f[14] + f[17]
     * + f[18] ); u[2] = invRho * ( - f[3] - f[6] + f[7] - f[8] + f[9] + f[12] + f[15] - f[16] +
     * f[17] - f[18] ); T cMux = 0.; T cMuy = 0.; T cMuz = 0.;
     *
     *     // Computation of central moments in a first time
     *     for (int i = 0; i<19; ++i) {
     *
     *         cMux = D::c[i][0]- u[0];
     *         cMuy = D::c[i][1]- u[1];
     *         cMuz = D::c[i][2]- u[2];
     *
     *         // // Order 0
     *         // CM[M000] += f[i];
     *
     *         // Order 1
     *         CM[M100] += cMux * f[i];
     *         CM[M010] += cMuy * f[i];
     *         CM[M001] += cMuz * f[i];
     *
     *         // Order 2
     *         CM[M200] += cMux * cMux * f[i];
     *         CM[M020] += cMuy * cMuy * f[i];
     *         CM[M002] += cMuz * cMuz * f[i];
     *         CM[M110] += cMux * cMuy * f[i];
     *         CM[M101] += cMux * cMuz * f[i];
     *         CM[M011] += cMuy * cMuz * f[i];
     *
     *         // Order 3
     *         CM[M210] += cMux * cMux * cMuy * f[i];
     *         CM[M201] += cMux * cMux * cMuz * f[i];
     *         CM[M021] += cMuy * cMuy * cMuz * f[i];
     *         CM[M120] += cMux * cMuy * cMuy * f[i];
     *         CM[M102] += cMux * cMuz * cMuz * f[i];
     *         CM[M012] += cMuy * cMuz * cMuz * f[i];
     *
     *         // Order 4
     *         CM[M220] += cMux * cMux * cMuy * cMuy * f[i];
     *         CM[M202] += cMux * cMux * cMuz * cMuz * f[i];
     *         CM[M022] += cMuy * cMuy * cMuz * cMuz * f[i];
     *     }
     *
     *     // Normalize before the computation of cumulants !
     *     for (int i = 1; i<19; ++i) {
     *         CM[i] *= invRho;
     *     }
     *
     *     // Computation of cumulants through central moments
     *     K[M000] = CM[M000];
     *     K[M100] = CM[M100] + u[0];
     *     K[M010] = CM[M010] + u[1];
     *     K[M001] = CM[M001] + u[2];
     *     K[M200] = CM[M200];
     *     K[M020] = CM[M020];
     *     K[M002] = CM[M002];
     *     K[M110] = CM[M110];
     *     K[M101] = CM[M101];
     *     K[M011] = CM[M011];
     *     K[M210] = CM[M210];
     *     K[M201] = CM[M201];
     *     K[M021] = CM[M021];
     *     K[M120] = CM[M120];
     *     K[M102] = CM[M102];
     *     K[M012] = CM[M012];
     *     K[M220] = CM[M220] - CM[M200]*CM[M020] - 2.*CM[M110]*CM[M110];
     *     K[M202] = CM[M202] - CM[M200]*CM[M002] - 2.*CM[M101]*CM[M101];
     *     K[M022] = CM[M022] - CM[M020]*CM[M002] - 2.*CM[M011]*CM[M011];
     * };
     */

    // Optimized way to compute Ks based on the general ordering of discrete velocities
    static void KcomputeMoments(
        Array<T, D::q> const &cell, Array<T, D::q> &K, T &rho, Array<T, D::d> &u)
    {
        Array<T, D::q> f;
        for (int i = 0; i < 19; ++i) {
            f[i] = cell[i] + D::SkordosFactor() * D::t[i];
            K[i] = 0.;
        }
        T X_M1 = f[FM00] + f[FMM0] + f[FMP0] + f[FM0M] + f[FM0P];
        T X_P1 = f[FP00] + f[FPP0] + f[FPM0] + f[FP0P] + f[FP0M];
        T X_0 =
            f[F000] + f[F0M0] + f[F00M] + f[F0MM] + f[F0MP] + f[F0P0] + f[F00P] + f[F0PP] + f[F0PM];
        T Y_M1 = f[F0M0] + f[FMM0] + f[F0MM] + f[F0MP] + f[FPM0];
        T Y_P1 = f[FMP0] + f[F0P0] + f[FPP0] + f[F0PP] + f[F0PM];
        T Z_M1 = f[F00M] + f[FM0M] + f[F0MM] + f[FP0M] + f[F0PM];
        T Z_P1 = f[FM0P] + f[F0MP] + f[F00P] + f[FP0P] + f[F0PP];

        rho = X_M1 + X_P1 + X_0;
        T invRho = 1. / rho;
        T two_invRho = 2. * invRho;

        // Order 0
        K[M000] = 1.;
        // Order 1
        K[M100] = invRho * (X_P1 - X_M1);
        K[M010] = invRho * (Y_P1 - Y_M1);
        K[M001] = invRho * (Z_P1 - Z_M1);
        // Order 2
        K[M200] = invRho * (X_M1 + X_P1);
        K[M020] = invRho * (Y_M1 + Y_P1);
        K[M002] = invRho * (Z_M1 + Z_P1);
        K[M110] = invRho * (f[FMM0] - f[FMP0] + f[FPP0] - f[FPM0]);
        K[M101] = invRho * (f[FM0M] - f[FM0P] + f[FP0P] - f[FP0M]);
        K[M011] = invRho * (f[F0MM] - f[F0MP] + f[F0PP] - f[F0PM]);
        // Order 3
        K[M210] = K[M110] - two_invRho * (f[FMM0] - f[FMP0]);
        K[M201] = K[M101] - two_invRho * (f[FM0M] - f[FM0P]);
        K[M021] = K[M011] - two_invRho * (f[F0MM] - f[F0MP]);
        K[M120] = K[M110] - two_invRho * (f[FMM0] - f[FPM0]);
        K[M102] = K[M101] - two_invRho * (f[FM0M] - f[FP0M]);
        K[M012] = K[M011] - two_invRho * (f[F0MM] - f[F0PM]);
        // Order 4
        K[M220] = K[M110] + two_invRho * (f[FMP0] + f[FPM0]);
        K[M202] = K[M101] + two_invRho * (f[FM0P] + f[FP0M]);
        K[M022] = K[M011] + two_invRho * (f[F0MP] + f[F0PM]);

        // Compute CMs from RMs using binomial formulas
        u[0] = K[M100];
        u[1] = K[M010];
        u[2] = K[M001];
        T ux2 = u[0] * u[0];
        T uy2 = u[1] * u[1];
        T uz2 = u[2] * u[2];

        K[M200] -= (ux2);
        K[M020] -= (uy2);
        K[M002] -= (uz2);

        K[M110] -= (u[0] * u[1]);
        K[M101] -= (u[0] * u[2]);
        K[M011] -= (u[1] * u[2]);

        K[M210] -= (u[1] * K[M200] + 2. * u[0] * K[M110] + ux2 * u[1]);
        K[M201] -= (u[2] * K[M200] + 2. * u[0] * K[M101] + ux2 * u[2]);
        K[M021] -= (u[2] * K[M020] + 2. * u[1] * K[M011] + uy2 * u[2]);
        K[M120] -= (u[0] * K[M020] + 2. * u[1] * K[M110] + u[0] * uy2);
        K[M102] -= (u[0] * K[M002] + 2. * u[2] * K[M101] + u[0] * uz2);
        K[M012] -= (u[1] * K[M002] + 2. * u[2] * K[M011] + u[1] * uz2);

        K[M220] -=
            (2. * u[1] * K[M210] + 2. * u[0] * K[M120] + uy2 * K[M200] + ux2 * K[M020]
             + 4. * u[0] * u[1] * K[M110] + ux2 * uy2);
        K[M202] -=
            (2. * u[2] * K[M201] + 2. * u[0] * K[M102] + uz2 * K[M200] + ux2 * K[M002]
             + 4. * u[0] * u[2] * K[M101] + ux2 * uz2);
        K[M022] -=
            (2. * u[2] * K[M021] + 2. * u[1] * K[M012] + uz2 * K[M020] + uy2 * K[M002]
             + 4. * u[1] * u[2] * K[M011] + uy2 * uz2);

        // Computation of Ks through CMs
        K[M220] -= (K[M200] * K[M020] + 2. * K[M110] * K[M110]);
        K[M202] -= (K[M200] * K[M002] + 2. * K[M101] * K[M101]);
        K[M022] -= (K[M020] * K[M002] + 2. * K[M011] * K[M011]);
    };

    static void KcomputeEquilibriumMoments(Array<T, D::d> const &u, Array<T, D::q> &Keq)
    {
        // Order 0
        Keq[M000] = 1.;
        // Order 1
        Keq[M100] = u[0];
        Keq[M010] = u[1];
        Keq[M001] = u[2];
        // Order 2
        Keq[M200] = D::cs2;
        Keq[M020] = D::cs2;
        Keq[M002] = D::cs2;
        Keq[M110] = 0.;
        Keq[M101] = 0.;
        Keq[M011] = 0.;
        // Order 3
        Keq[M210] = 0.;
        Keq[M201] = 0.;
        Keq[M021] = 0.;
        Keq[M120] = 0.;
        Keq[M102] = 0.;
        Keq[M012] = 0.;
        // Order 4
        Keq[M220] = 0.;
        Keq[M202] = 0.;
        Keq[M022] = 0.;
    };

    // Equilibrium populations based on 19 moments can be computed using either RM, HM, CM, CHM or
    // Gauss-Hermite formalisms. All formulations are equivalent for the D3Q19 as long as we stick
    // to non-weighted formulas (i.e, no Gauss-Hermite) so we use the fastest one (RMs)
    static void KcomputeEquilibrium(
        T rho, Array<T, D::d> const &u, Array<T, D::q> const &Keq, Array<T, D::q> &eq)
    {
        Array<T, D::q> RMeq;
        // Order 2
        RMeq[M200] = u[0] * u[0] + D::cs2;
        RMeq[M020] = u[1] * u[1] + D::cs2;
        RMeq[M002] = u[2] * u[2] + D::cs2;
        RMeq[M110] = u[0] * u[1];
        RMeq[M101] = u[0] * u[2];
        RMeq[M011] = u[1] * u[2];
        // Order 3
        RMeq[M210] = RMeq[M200] * u[1];
        RMeq[M201] = RMeq[M200] * u[2];
        RMeq[M021] = RMeq[M020] * u[2];
        RMeq[M120] = RMeq[M020] * u[0];
        RMeq[M102] = RMeq[M002] * u[0];
        RMeq[M012] = RMeq[M002] * u[1];
        // Order 4
        RMeq[M220] = RMeq[M200] * RMeq[M020];
        RMeq[M202] = RMeq[M200] * RMeq[M002];
        RMeq[M022] = RMeq[M020] * RMeq[M002];

        // Optimization based on symmetries between populations and their opposite counterpart
        eq[F000] =
            rho
            * (1. - RMeq[M200] - RMeq[M020] - RMeq[M002] + RMeq[M220] + RMeq[M202] + RMeq[M022]);

        eq[FP00] =
            0.5 * rho * (u[0] + RMeq[M200] - RMeq[M120] - RMeq[M102] - RMeq[M220] - RMeq[M202]);
        eq[FM00] = rho * (-u[0] + RMeq[M120] + RMeq[M102]) + eq[FP00];

        eq[F0P0] =
            0.5 * rho * (u[1] + RMeq[M020] - RMeq[M210] - RMeq[M012] - RMeq[M220] - RMeq[M022]);
        eq[F0M0] = rho * (-u[1] + RMeq[M210] + RMeq[M012]) + eq[F0P0];

        eq[F00P] =
            0.5 * rho * (u[2] + RMeq[M002] - RMeq[M201] - RMeq[M021] - RMeq[M202] - RMeq[M022]);
        eq[F00M] = rho * (-u[2] + RMeq[M201] + RMeq[M021]) + eq[F00P];

        eq[FPP0] = 0.25 * rho * (RMeq[M110] + RMeq[M210] + RMeq[M120] + RMeq[M220]);
        eq[FMP0] = 0.5 * rho * (-RMeq[M110] - RMeq[M120]) + eq[FPP0];
        eq[FPM0] = 0.5 * rho * (-RMeq[M110] - RMeq[M210]) + eq[FPP0];
        eq[FMM0] = 0.5 * rho * (-RMeq[M210] - RMeq[M120]) + eq[FPP0];

        eq[FP0P] = 0.25 * rho * (RMeq[M101] + RMeq[M201] + RMeq[M102] + RMeq[M202]);
        eq[FM0P] = 0.5 * rho * (-RMeq[M101] - RMeq[M102]) + eq[FP0P];
        eq[FP0M] = 0.5 * rho * (-RMeq[M101] - RMeq[M201]) + eq[FP0P];
        eq[FM0M] = 0.5 * rho * (-RMeq[M201] - RMeq[M102]) + eq[FP0P];

        eq[F0PP] = 0.25 * rho * (RMeq[M011] + RMeq[M021] + RMeq[M012] + RMeq[M022]);
        eq[F0MP] = 0.5 * rho * (-RMeq[M011] - RMeq[M012]) + eq[F0PP];
        eq[F0PM] = 0.5 * rho * (-RMeq[M011] - RMeq[M021]) + eq[F0PP];
        eq[F0MM] = 0.5 * rho * (-RMeq[M021] - RMeq[M012]) + eq[F0PP];
    };

    static void Kcollide(
        Array<T, D::q> &cell, T rho, Array<T, D::d> const &u, Array<T, D::q> const &K,  // Cumulants
        Array<T, D::q> const &Keq,  // Equilibrium cumulants
        Array<T, D::numRelaxationTimes> const &omega)
    {
        T omega1 = omega[0];
        T omega2 = omega[1];
        T omega3 = omega[2];
        T omega4 = omega[3];
        T omegaBulk = omega[4];
        T omegaPlus = (omegaBulk + 2. * omega1) / 3.;  // Notation used by Fei
        T omegaMinus = (omegaBulk - omega1) / 3.;      // Notation used by Fei

        T ux2 = u[0] * u[0];
        T uy2 = u[1] * u[1];
        T uz2 = u[2] * u[2];

        // Post-collision moments.
        Array<T, D::q> Kcoll;
        Array<T, D::q> CMcoll;
        Array<T, D::q> RMcoll;

        // Collision in the cumulant space
        // Order 2 (non-diagonal collision so that we can easily modify the bulk viscosity)
        Kcoll[M200] = K[M200] - omegaPlus * (K[M200] - Keq[M200])
                      - omegaMinus * (K[M020] - Keq[M020]) - omegaMinus * (K[M002] - Keq[M002]);
        Kcoll[M020] = K[M020] - omegaMinus * (K[M200] - Keq[M200])
                      - omegaPlus * (K[M020] - Keq[M020]) - omegaMinus * (K[M002] - Keq[M002]);
        Kcoll[M002] = K[M002] - omegaMinus * (K[M200] - Keq[M200])
                      - omegaMinus * (K[M020] - Keq[M020]) - omegaPlus * (K[M002] - Keq[M002]);

        Kcoll[M110] = (1. - omega2) * K[M110] + omega2 * Keq[M110];
        Kcoll[M101] = (1. - omega2) * K[M101] + omega2 * Keq[M101];
        Kcoll[M011] = (1. - omega2) * K[M011] + omega2 * Keq[M011];

        // Order 3
        Kcoll[M210] = (1. - omega3) * K[M210] + omega3 * Keq[M210];
        Kcoll[M201] = (1. - omega3) * K[M201] + omega3 * Keq[M201];
        Kcoll[M021] = (1. - omega3) * K[M021] + omega3 * Keq[M021];
        Kcoll[M120] = (1. - omega3) * K[M120] + omega3 * Keq[M120];
        Kcoll[M102] = (1. - omega3) * K[M102] + omega3 * Keq[M102];
        Kcoll[M012] = (1. - omega3) * K[M012] + omega3 * Keq[M012];

        // Order 4
        Kcoll[M220] = (1. - omega4) * K[M220] + omega4 * Keq[M220];
        Kcoll[M202] = (1. - omega4) * K[M202] + omega4 * Keq[M202];
        Kcoll[M022] = (1. - omega4) * K[M022] + omega4 * Keq[M022];

        // Come back to CMcoll using modifying fourth- and higher-order post-collision cumulants
        CMcoll[M200] = Kcoll[M200];
        CMcoll[M020] = Kcoll[M020];
        CMcoll[M002] = Kcoll[M002];
        CMcoll[M110] = Kcoll[M110];
        CMcoll[M101] = Kcoll[M101];
        CMcoll[M011] = Kcoll[M011];

        CMcoll[M210] = Kcoll[M210];
        CMcoll[M201] = Kcoll[M201];
        CMcoll[M021] = Kcoll[M021];
        CMcoll[M120] = Kcoll[M120];
        CMcoll[M102] = Kcoll[M102];
        CMcoll[M012] = Kcoll[M012];

        CMcoll[M220] = Kcoll[M220] + Kcoll[M200] * Kcoll[M020] + 2. * Kcoll[M110] * Kcoll[M110];
        CMcoll[M202] = Kcoll[M202] + Kcoll[M200] * Kcoll[M002] + 2. * Kcoll[M101] * Kcoll[M101];
        CMcoll[M022] = Kcoll[M022] + Kcoll[M020] * Kcoll[M002] + 2. * Kcoll[M011] * Kcoll[M011];

        // Come back to RMcoll using binomial formulas
        RMcoll[M200] = CMcoll[M200] + ux2;
        RMcoll[M020] = CMcoll[M020] + uy2;
        RMcoll[M002] = CMcoll[M002] + uz2;

        RMcoll[M110] = CMcoll[M110] + u[0] * u[1];
        RMcoll[M101] = CMcoll[M101] + u[0] * u[2];
        RMcoll[M011] = CMcoll[M011] + u[1] * u[2];

        RMcoll[M210] = CMcoll[M210] + u[1] * CMcoll[M200] + 2. * u[0] * CMcoll[M110] + ux2 * u[1];
        RMcoll[M201] = CMcoll[M201] + u[2] * CMcoll[M200] + 2. * u[0] * CMcoll[M101] + ux2 * u[2];
        RMcoll[M021] = CMcoll[M021] + u[2] * CMcoll[M020] + 2. * u[1] * CMcoll[M011] + uy2 * u[2];
        RMcoll[M120] = CMcoll[M120] + u[0] * CMcoll[M020] + 2. * u[1] * CMcoll[M110] + u[0] * uy2;
        RMcoll[M102] = CMcoll[M102] + u[0] * CMcoll[M002] + 2. * u[2] * CMcoll[M101] + u[0] * uz2;
        RMcoll[M012] = CMcoll[M012] + u[1] * CMcoll[M002] + 2. * u[2] * CMcoll[M011] + u[1] * uz2;

        RMcoll[M220] = CMcoll[M220] + 2. * u[1] * CMcoll[M210] + 2. * u[0] * CMcoll[M120]
                       + uy2 * CMcoll[M200] + ux2 * CMcoll[M020] + 4. * u[0] * u[1] * CMcoll[M110]
                       + ux2 * uy2;
        RMcoll[M202] = CMcoll[M202] + 2. * u[2] * CMcoll[M201] + 2. * u[0] * CMcoll[M102]
                       + uz2 * CMcoll[M200] + ux2 * CMcoll[M002] + 4. * u[0] * u[2] * CMcoll[M101]
                       + ux2 * uz2;
        RMcoll[M022] = CMcoll[M022] + 2. * u[2] * CMcoll[M021] + 2. * u[1] * CMcoll[M012]
                       + uz2 * CMcoll[M020] + uy2 * CMcoll[M002] + 4. * u[1] * u[2] * CMcoll[M011]
                       + uy2 * uz2;

        // Compute post collision populations from RM
        // Optimization based on symmetries between populations and their opposite counterpart
        cell[F000] = rho
                     * (1. - RMcoll[M200] - RMcoll[M020] - RMcoll[M002] + RMcoll[M220]
                        + RMcoll[M202] + RMcoll[M022]);

        cell[FP00] =
            0.5 * rho
            * (u[0] + RMcoll[M200] - RMcoll[M120] - RMcoll[M102] - RMcoll[M220] - RMcoll[M202]);
        cell[FM00] = rho * (-u[0] + RMcoll[M120] + RMcoll[M102]) + cell[FP00];

        cell[F0P0] =
            0.5 * rho
            * (u[1] + RMcoll[M020] - RMcoll[M210] - RMcoll[M012] - RMcoll[M220] - RMcoll[M022]);
        cell[F0M0] = rho * (-u[1] + RMcoll[M210] + RMcoll[M012]) + cell[F0P0];

        cell[F00P] =
            0.5 * rho
            * (u[2] + RMcoll[M002] - RMcoll[M201] - RMcoll[M021] - RMcoll[M202] - RMcoll[M022]);
        cell[F00M] = rho * (-u[2] + RMcoll[M201] + RMcoll[M021]) + cell[F00P];

        cell[FPP0] = 0.25 * rho * (RMcoll[M110] + RMcoll[M210] + RMcoll[M120] + RMcoll[M220]);
        cell[FMP0] = 0.5 * rho * (-RMcoll[M110] - RMcoll[M120]) + cell[FPP0];
        cell[FPM0] = 0.5 * rho * (-RMcoll[M110] - RMcoll[M210]) + cell[FPP0];
        cell[FMM0] = 0.5 * rho * (-RMcoll[M210] - RMcoll[M120]) + cell[FPP0];

        cell[FP0P] = 0.25 * rho * (RMcoll[M101] + RMcoll[M201] + RMcoll[M102] + RMcoll[M202]);
        cell[FM0P] = 0.5 * rho * (-RMcoll[M101] - RMcoll[M102]) + cell[FP0P];
        cell[FP0M] = 0.5 * rho * (-RMcoll[M101] - RMcoll[M201]) + cell[FP0P];
        cell[FM0M] = 0.5 * rho * (-RMcoll[M201] - RMcoll[M102]) + cell[FP0P];

        cell[F0PP] = 0.25 * rho * (RMcoll[M011] + RMcoll[M021] + RMcoll[M012] + RMcoll[M022]);
        cell[F0MP] = 0.5 * rho * (-RMcoll[M011] - RMcoll[M012]) + cell[F0PP];
        cell[F0PM] = 0.5 * rho * (-RMcoll[M011] - RMcoll[M021]) + cell[F0PP];
        cell[F0MM] = 0.5 * rho * (-RMcoll[M021] - RMcoll[M012]) + cell[F0PP];

        for (int i = 0; i < 19; ++i) {
            cell[i] -= D::SkordosFactor() * D::t[i];
        }
    };

    /////////////////////////////////////////////////////////////////////////////////
    // Gauss-Hermite Formalism (Equilibrium is computed through the GH  formalism) //
    /////////////////////////////////////////////////////////////////////////////////

    /**
     * // General way to compute GHs
     * static void GHcomputeMoments(Array<T,D::q> const& cell, Array<T, D::q>& GH, T& rho) {
     *
     *     Array<T, D::q> f;
     *     for (int i = 0; i<19; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *         GH[i] = 0.;
     *     }
     *
     *     for (int i = 0; i<19; ++i) {
     *         T Hxx = D::c[i][0] * D::c[i][0] - D::cs2;
     *         T Hyy = D::c[i][1] * D::c[i][1] - D::cs2;
     *         T Hzz = D::c[i][2] * D::c[i][2] - D::cs2;
     *
     *         // Order 0
     *         GH[M000] += f[i];
     *         // Order 1
     *         GH[M100] += D::c[i][0] * f[i];
     *         GH[M010] += D::c[i][1] * f[i];
     *         GH[M001] += D::c[i][2] * f[i];
     *         // Order 2
     *         GH[M200] += Hxx * f[i];
     *         GH[M020] += Hyy * f[i];
     *         GH[M002] += Hzz * f[i];
     *         GH[M110] += D::c[i][0] * D::c[i][1] * f[i];
     *         GH[M101] += D::c[i][0] * D::c[i][2] * f[i];
     *         GH[M011] += D::c[i][1] * D::c[i][2] * f[i];
     *         // Order 3
     *         GH[M210] += Hxx * D::c[i][1] * f[i];
     *         GH[M201] += Hxx * D::c[i][2] * f[i];
     *         GH[M021] += Hyy * D::c[i][2] * f[i];
     *         GH[M120] += D::c[i][0] * Hyy * f[i];
     *         GH[M102] += D::c[i][0] * Hzz * f[i];
     *         GH[M012] += D::c[i][1] * Hzz * f[i];
     *         // Order 4
     *         GH[M220] += Hxx * Hyy * f[i];
     *         GH[M202] += Hxx * Hzz * f[i];
     *         GH[M022] += Hyy * Hzz * f[i];
     *     }
     *     rho = GH[M000];
     *     T invRho = 1. / rho;
     *     for (int i = 0; i<19; ++i) {
     *         GH[i] *= invRho;
     *     }
     * };
     *
     * // Optimized way to compute HMs based on Palabos ordering of discrete velocities
     * static void GHcomputeMoments(Array<T,D::q> const& cell, Array<T, D::q>& GH, T& rho) {
     *
     *     Array<T, D::q> f;
     *     for (int i = 0; i<19; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *     }
     *
     *     T a1 = 1./3. ;
     *     T a2 = 2./3. ;
     *
     *     T b1 = 1./9. ;
     *     T b2 = 2./9. ;
     *     T b3 = 4./9. ;
     *
     *     // Order 0
     *     rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[10] + f[11]
     * + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18]; GH[M000] = 1.0; T invRho = 1./rho;
     *
     *     // Order 1
     *     GH[M100] = invRho * ( - f[1] - f[4] - f[5] - f[6] - f[7] + f[10] + f[13] + f[14] + f[15]
     * + f[16] ); GH[M010] = invRho * ( - f[2] - f[4] + f[5] - f[8] - f[9] + f[11] + f[13] - f[14] +
     * f[17] + f[18] ); GH[M001] = invRho * ( - f[3] - f[6] + f[7] - f[8] + f[9] + f[12] + f[15] -
     * f[16] + f[17] - f[18] );
     *     // Order 2
     *     GH[M200] = invRho * ( - a1*f[0] + a2*f[1] - a1*f[2] - a1*f[3] + a2*f[4] + a2*f[5] +
     * a2*f[6] + a2*f[7] - a1*f[8] - a1*f[9] + a2*f[10] - a1*f[11] - a1*f[12] + a2*f[13] + a2*f[14]
     * + a2*f[15] + a2*f[16] - a1*f[17] - a1*f[18] ); GH[M020] = invRho * ( - a1*f[0] - a1*f[1] +
     * a2*f[2] - a1*f[3] + a2*f[4] + a2*f[5] - a1*f[6] - a1*f[7] + a2*f[8] + a2*f[9] - a1*f[10] +
     * a2*f[11] - a1*f[12] + a2*f[13] + a2*f[14] - a1*f[15] - a1*f[16] + a2*f[17] + a2*f[18] );
     *     GH[M002] = invRho * ( - a1*f[0] - a1*f[1] - a1*f[2] + a2*f[3] - a1*f[4] - a1*f[5] +
     * a2*f[6] + a2*f[7] + a2*f[8] + a2*f[9] - a1*f[10] - a1*f[11] + a2*f[12] - a1*f[13] - a1*f[14]
     * + a2*f[15] + a2*f[16] + a2*f[17] + a2*f[18] ); GH[M110] = invRho * ( f[4] - f[5] + f[13] -
     * f[14] ); GH[M101] = invRho * ( f[6] - f[7] + f[15] - f[16] ); GH[M011] = invRho * ( f[8] -
     * f[9] + f[17] - f[18] );
     *     // Order 3
     *     GH[M210] = invRho * ( a1*f[2] - a2*f[4] + a2*f[5] + a1*f[8] + a1*f[9] - a1*f[11] +
     * a2*f[13] - a2*f[14] - a1*f[17] - a1*f[18] ); GH[M201] = invRho * ( a1*f[3] - a2*f[6] +
     * a2*f[7] + a1*f[8] - a1*f[9] - a1*f[12] + a2*f[15] - a2*f[16] - a1*f[17] + a1*f[18] );
     *     GH[M021] = invRho * ( a1*f[3] + a1*f[6] - a1*f[7] - a2*f[8] + a2*f[9] - a1*f[12] -
     * a1*f[15] + a1*f[16] + a2*f[17] - a2*f[18] ); GH[M120] = invRho * ( a1*f[1] - a2*f[4] -
     * a2*f[5] + a1*f[6] + a1*f[7] - a1*f[10] + a2*f[13] + a2*f[14] - a1*f[15] - a1*f[16] );
     *     GH[M102] = invRho * ( a1*f[1] + a1*f[4] + a1*f[5] - a2*f[6] - a2*f[7] - a1*f[10] -
     * a1*f[13] - a1*f[14] + a2*f[15] + a2*f[16] ); GH[M012] = invRho * ( a1*f[2] + a1*f[4] -
     * a1*f[5] - a2*f[8] - a2*f[9] - a1*f[11] - a1*f[13] + a1*f[14] + a2*f[17] + a2*f[18] );
     *     // Order 4
     *     GH[M220] = invRho * ( b1*f[0] - b2*f[1] - b2*f[2] + b1*f[3] + b3*f[4] + b3*f[5] - b2*f[6]
     * - b2*f[7] - b2*f[8] - b2*f[9] - b2*f[10] - b2*f[11] + b1*f[12] + b3*f[13] + b3*f[14] -
     * b2*f[15] - b2*f[16] - b2*f[17] - b2*f[18] ); GH[M202] = invRho * ( b1*f[0] - b2*f[1] +
     * b1*f[2] - b2*f[3] - b2*f[4] - b2*f[5] + b3*f[6] + b3*f[7] - b2*f[8] - b2*f[9] - b2*f[10] +
     * b1*f[11] - b2*f[12] - b2*f[13] - b2*f[14] + b3*f[15] + b3*f[16] - b2*f[17] - b2*f[18] );
     *     GH[M022] = invRho * ( b1*f[0] + b1*f[1] - b2*f[2] - b2*f[3] - b2*f[4] - b2*f[5] - b2*f[6]
     * - b2*f[7] + b3*f[8] + b3*f[9] + b1*f[10] - b2*f[11] - b2*f[12] - b2*f[13] - b2*f[14] -
     * b2*f[15] - b2*f[16] + b3*f[17] + b3*f[18] );
     * };
     */

    // Optimized way to compute GHs based on the general ordering of discrete velocities
    static void GHcomputeMoments(Array<T, D::q> const &cell, Array<T, D::q> &GH, T &rho)
    {
        Array<T, D::q> f;
        for (int i = 0; i < 19; ++i) {
            f[i] = cell[i] + D::SkordosFactor() * D::t[i];
        }
        T X_M1 = f[FM00] + f[FMM0] + f[FMP0] + f[FM0M] + f[FM0P];
        T X_P1 = f[FP00] + f[FPP0] + f[FPM0] + f[FP0P] + f[FP0M];
        T X_0 =
            f[F000] + f[F0M0] + f[F00M] + f[F0MM] + f[F0MP] + f[F0P0] + f[F00P] + f[F0PP] + f[F0PM];
        T Y_M1 = f[F0M0] + f[FMM0] + f[F0MM] + f[F0MP] + f[FPM0];
        T Y_P1 = f[FMP0] + f[F0P0] + f[FPP0] + f[F0PP] + f[F0PM];
        T Z_M1 = f[F00M] + f[FM0M] + f[F0MM] + f[FP0M] + f[F0PM];
        T Z_P1 = f[FM0P] + f[F0MP] + f[F00P] + f[FP0P] + f[F0PP];

        rho = X_M1 + X_P1 + X_0;
        T invRho = 1. / rho;
        T two_invRho = 2. * invRho;

        // Order 0
        GH[M000] = 1.;
        // Order 1
        GH[M100] = invRho * (X_P1 - X_M1);
        GH[M010] = invRho * (Y_P1 - Y_M1);
        GH[M001] = invRho * (Z_P1 - Z_M1);
        // Order 2
        GH[M200] = invRho * (X_M1 + X_P1);
        GH[M020] = invRho * (Y_M1 + Y_P1);
        GH[M002] = invRho * (Z_M1 + Z_P1);
        GH[M110] = invRho * (f[FMM0] - f[FMP0] + f[FPP0] - f[FPM0]);
        GH[M101] = invRho * (f[FM0M] - f[FM0P] + f[FP0P] - f[FP0M]);
        GH[M011] = invRho * (f[F0MM] - f[F0MP] + f[F0PP] - f[F0PM]);
        // Order 3
        GH[M210] = GH[M110] - two_invRho * (f[FMM0] - f[FMP0]);
        GH[M201] = GH[M101] - two_invRho * (f[FM0M] - f[FM0P]);
        GH[M021] = GH[M011] - two_invRho * (f[F0MM] - f[F0MP]);
        GH[M120] = GH[M110] - two_invRho * (f[FMM0] - f[FPM0]);
        GH[M102] = GH[M101] - two_invRho * (f[FM0M] - f[FP0M]);
        GH[M012] = GH[M011] - two_invRho * (f[F0MM] - f[F0PM]);
        // Order 4
        GH[M220] = GH[M110] + two_invRho * (f[FMP0] + f[FPM0]);
        GH[M202] = GH[M101] + two_invRho * (f[FM0P] + f[FP0M]);
        GH[M022] = GH[M011] + two_invRho * (f[F0MP] + f[F0PM]);

        // We come back to Hermite moments
        T cs4 = D::cs2 * D::cs2;
        GH[M200] -= D::cs2;
        GH[M020] -= D::cs2;
        GH[M002] -= D::cs2;

        GH[M210] -= D::cs2 * GH[M010];
        GH[M201] -= D::cs2 * GH[M001];
        GH[M021] -= D::cs2 * GH[M001];
        GH[M120] -= D::cs2 * GH[M100];
        GH[M102] -= D::cs2 * GH[M100];
        GH[M012] -= D::cs2 * GH[M010];

        GH[M220] -= (D::cs2 * (GH[M200] + GH[M020]) + cs4);
        GH[M202] -= (D::cs2 * (GH[M200] + GH[M002]) + cs4);
        GH[M022] -= (D::cs2 * (GH[M020] + GH[M002]) + cs4);
    };

    static void GHcomputeEquilibriumMoments(Array<T, D::d> const &u, Array<T, D::q> &GHeq)
    {
        // Order 0
        GHeq[M000] = 1.;
        // Order 1
        GHeq[M100] = u[0];
        GHeq[M010] = u[1];
        GHeq[M001] = u[2];
        // Order 2
        GHeq[M200] = u[0] * u[0];
        GHeq[M020] = u[1] * u[1];
        GHeq[M002] = u[2] * u[2];
        GHeq[M110] = u[0] * u[1];
        GHeq[M101] = u[0] * u[2];
        GHeq[M011] = u[1] * u[2];
        // Order 3
        GHeq[M210] = GHeq[M200] * u[1];
        GHeq[M201] = GHeq[M200] * u[2];
        GHeq[M021] = GHeq[M020] * u[2];
        GHeq[M120] = GHeq[M020] * u[0];
        GHeq[M102] = GHeq[M002] * u[0];
        GHeq[M012] = GHeq[M002] * u[1];
        // Order 4
        GHeq[M220] = GHeq[M200] * GHeq[M020];
        GHeq[M202] = GHeq[M200] * GHeq[M002];
        GHeq[M022] = GHeq[M020] * GHeq[M002];
    };

    // Equilibrium populations based on 19 moments can be computed using either RM, HM, CM, CHM or
    // Gauss-Hermite formalisms. Here we use the Gauss-Hermite formalism that requires the
    // orthogonalisation of moments due to spurious couplings between third-order moments, and
    // second/fourth order moments.
    static void GHcomputeEquilibrium(T rho, Array<T, D::q> const &GHeq, Array<T, D::q> &eq)
    {
        Array<T, D::d> u(GHeq[1], GHeq[2], GHeq[3]);

        // Orthogonalization
        // Order 3
        T GHeq_M210 = GHeq[M210] + GHeq[M012];
        T GHeq_M201 = GHeq[M201] + GHeq[M021];
        T GHeq_M120 = GHeq[M120] + GHeq[M102];
        T GHeq_M012 = GHeq[M012] - GHeq[M210];
        T GHeq_M021 = GHeq[M021] - GHeq[M201];
        T GHeq_M102 = GHeq[M102] - GHeq[M120];
        // 4 ortho with 2
        T GHeq_M220 = GHeq[M220] + (1. / 6.) * GHeq[M002];
        T GHeq_M202 = GHeq[M202] + (1. / 6.) * GHeq[M020];
        T GHeq_M022 = GHeq[M022] + (1. / 6.) * GHeq[M200];
        // 4 ortho with 4
        GHeq_M202 = GHeq[M202] + (2. / 7.) * GHeq_M220;
        GHeq_M022 = GHeq[M022] + (2. / 7.) * GHeq_M220 + (2. / 5.) * GHeq_M202;

        eq[F000] = (rho * D::t[F000])
                   * (1. - 1.5 * (GHeq[M200] + GHeq[M020] + GHeq[M002]) + (9. / 7.) * GHeq_M220
                      + (9. / 5.) * GHeq_M202 + 3. * GHeq_M022);

        eq[FP00] = (rho * D::t[FP00])
                   * (1. + 3. * u[0] + 1.5 * (2. * GHeq[M200] - GHeq[M020] - GHeq[M002])
                      - 9. * GHeq_M120 - (45. / 7.) * GHeq_M220 - 9. * GHeq_M202);
        eq[FM00] = (rho * D::t[FM00])
                   * (1. - 3. * u[0] + 1.5 * (2. * GHeq[M200] - GHeq[M020] - GHeq[M002])
                      + 9. * GHeq_M120 - (45. / 7.) * GHeq_M220 - 9. * GHeq_M202);

        eq[F0P0] =
            (rho * D::t[F0P0])
            * (1. + 3. * u[1] + 1.5 * (2. * GHeq[M020] - GHeq[M200] - GHeq[M002]) - 9. * GHeq_M210
               - (45. / 7.) * GHeq_M220 + (18. / 5.) * GHeq_M202 - 9. * GHeq_M022);
        eq[F0M0] =
            (rho * D::t[F0M0])
            * (1. - 3. * u[1] + 1.5 * (2. * GHeq[M020] - GHeq[M200] - GHeq[M002]) + 9. * GHeq_M210
               - (45. / 7.) * GHeq_M220 + (18. / 5.) * GHeq_M202 - 9. * GHeq_M022);

        eq[F00P] =
            (rho * D::t[F00P])
            * (1. + 3. * u[2] + 1.5 * (2. * GHeq[M002] - GHeq[M200] - GHeq[M020]) - 9. * GHeq_M201
               + (36. / 7.) * GHeq_M220 - (27. / 5.) * GHeq_M202 - 9. * GHeq_M022);
        eq[F00M] =
            (rho * D::t[F00M])
            * (1. - 3. * u[2] + 1.5 * (2. * GHeq[M002] - GHeq[M200] - GHeq[M020]) + 9. * GHeq_M201
               + (36. / 7.) * GHeq_M220 - (27. / 5.) * GHeq_M202 - 9. * GHeq_M022);

        eq[FPP0] = (rho * D::t[FPP0])
                   * (1. + 3. * (+u[0] + u[1])
                      + 1.5 * (2. * GHeq[M200] + 2. * GHeq[M020] - GHeq[M002]) + 9. * GHeq[M110]
                      + 4.5 * (GHeq_M210 + GHeq_M120 - GHeq_M102 - GHeq_M012) + 9. * GHeq_M220);
        eq[FMP0] = (rho * D::t[FMP0])
                   * (1. + 3. * (-u[0] + u[1])
                      + 1.5 * (2. * GHeq[M200] + 2. * GHeq[M020] - GHeq[M002]) - 9. * GHeq[M110]
                      + 4.5 * (GHeq_M210 - GHeq_M120 + GHeq_M102 - GHeq_M012) + 9. * GHeq_M220);
        eq[FPM0] = (rho * D::t[FPM0])
                   * (1. + 3. * (+u[0] - u[1])
                      + 1.5 * (2. * GHeq[M200] + 2. * GHeq[M020] - GHeq[M002]) - 9. * GHeq[M110]
                      + 4.5 * (-GHeq_M210 + GHeq_M120 - GHeq_M102 + GHeq_M012) + 9. * GHeq_M220);
        eq[FMM0] = (rho * D::t[FMM0])
                   * (1. + 3. * (-u[0] - u[1])
                      + 1.5 * (2. * GHeq[M200] + 2. * GHeq[M020] - GHeq[M002]) + 9. * GHeq[M110]
                      + 4.5 * (-GHeq_M210 - GHeq_M120 + GHeq_M102 + GHeq_M012) + 9. * GHeq_M220);

        eq[FP0P] =
            (rho * D::t[FP0P])
            * (1. + 3. * (+u[0] + u[2]) + 1.5 * (2. * GHeq[M200] + 2. * GHeq[M002] - GHeq[M020])
               + 9. * GHeq[M101] + 4.5 * (GHeq_M201 + GHeq_M102 + GHeq_M120 - GHeq_M021)
               - (18. / 7.) * GHeq_M220 + 9. * GHeq_M202);
        eq[FM0P] =
            (rho * D::t[FM0P])
            * (1. + 3. * (-u[0] + u[2]) + 1.5 * (2. * GHeq[M200] + 2. * GHeq[M002] - GHeq[M020])
               - 9. * GHeq[M101] + 4.5 * (GHeq_M201 - GHeq_M102 - GHeq_M120 - GHeq_M021)
               - (18. / 7.) * GHeq_M220 + 9. * GHeq_M202);
        eq[FP0M] =
            (rho * D::t[FP0M])
            * (1. + 3. * (+u[0] - u[2]) + 1.5 * (2. * GHeq[M200] + 2. * GHeq[M002] - GHeq[M020])
               - 9. * GHeq[M101] + 4.5 * (-GHeq_M201 + GHeq_M102 + GHeq_M120 + GHeq_M021)
               - (18. / 7.) * GHeq_M220 + 9. * GHeq_M202);
        eq[FM0M] =
            (rho * D::t[FM0M])
            * (1. + 3. * (-u[0] - u[2]) + 1.5 * (2. * GHeq[M200] + 2. * GHeq[M002] - GHeq[M020])
               + 9. * GHeq[M101] + 4.5 * (-GHeq_M201 - GHeq_M102 - GHeq_M120 + GHeq_M021)
               - (18. / 7.) * GHeq_M220 + 9. * GHeq_M202);

        eq[F0PP] =
            (rho * D::t[F0PP])
            * (1. + 3. * (+u[1] + u[2]) + 1.5 * (2. * GHeq[M020] + 2. * GHeq[M002] - GHeq[M200])
               + 9. * GHeq[M011] + 4.5 * (GHeq_M021 + GHeq_M012 + GHeq_M210 + GHeq_M201)
               - (18. / 7.) * GHeq_M220 - (18. / 5.) * GHeq_M202 + 9. * GHeq_M022);
        eq[F0MP] =
            (rho * D::t[F0MP])
            * (1. + 3. * (-u[1] + u[2]) + 1.5 * (2. * GHeq[M020] + 2. * GHeq[M002] - GHeq[M200])
               - 9. * GHeq[M011] + 4.5 * (GHeq_M021 - GHeq_M012 - GHeq_M210 + GHeq_M201)
               - (18. / 7.) * GHeq_M220 - (18. / 5.) * GHeq_M202 + 9. * GHeq_M022);
        eq[F0PM] =
            (rho * D::t[F0PM])
            * (1. + 3. * (+u[1] - u[2]) + 1.5 * (2. * GHeq[M020] + 2. * GHeq[M002] - GHeq[M200])
               - 9. * GHeq[M011] + 4.5 * (-GHeq_M021 + GHeq_M012 + GHeq_M210 - GHeq_M201)
               - (18. / 7.) * GHeq_M220 - (18. / 5.) * GHeq_M202 + 9. * GHeq_M022);
        eq[F0MM] =
            (rho * D::t[F0MM])
            * (1. + 3. * (-u[1] - u[2]) + 1.5 * (2. * GHeq[M020] + 2. * GHeq[M002] - GHeq[M200])
               + 9. * GHeq[M011] + 4.5 * (-GHeq_M021 - GHeq_M012 - GHeq_M210 - GHeq_M201)
               - (18. / 7.) * GHeq_M220 - (18. / 5.) * GHeq_M202 + 9. * GHeq_M022);
    };

    /////////// Full Ortho ///////////
    static void GHcollide(
        Array<T, D::q> &cell, T rho, Array<T, D::d> const &u,
        Array<T, D::q> const &GH,    // Hermite moments
        Array<T, D::q> const &GHeq,  // Equilibrium moments (Hermite)
        Array<T, D::numRelaxationTimes> const &omega)
    {
        T omega1 = omega[0];
        T omega2 = omega[1];
        T omega3 = omega[2];
        T omega4 = omega[3];
        T omegaBulk = omega[4];
        T omegaPlus = (omegaBulk + 2. * omega1) / 3.;  // Notation used by Fei
        T omegaMinus = (omegaBulk - omega1) / 3.;      // Notation used by Fei

        // Post-collision moments.
        Array<T, D::q> GHcoll;

        // Collision in the Hermite moment space
        // Order 2 (non-diagonal collision so that we can easily modify the bulk viscosity)
        GHcoll[M200] = GH[M200] - omegaPlus * (GH[M200] - GHeq[M200])
                       - omegaMinus * (GH[M020] - GHeq[M020])
                       - omegaMinus * (GH[M002] - GHeq[M002]);
        GHcoll[M020] = GH[M020] - omegaMinus * (GH[M200] - GHeq[M200])
                       - omegaPlus * (GH[M020] - GHeq[M020]) - omegaMinus * (GH[M002] - GHeq[M002]);
        GHcoll[M002] = GH[M002] - omegaMinus * (GH[M200] - GHeq[M200])
                       - omegaMinus * (GH[M020] - GHeq[M020]) - omegaPlus * (GH[M002] - GHeq[M002]);

        GHcoll[M110] = (1. - omega2) * GH[M110] + omega2 * GHeq[M110];
        GHcoll[M101] = (1. - omega2) * GH[M101] + omega2 * GHeq[M101];
        GHcoll[M011] = (1. - omega2) * GH[M011] + omega2 * GHeq[M011];

        // Orthogonalization
        // Order 3
        T GHeq_M210 = GHeq[M210] + GHeq[M012];
        T GHeq_M201 = GHeq[M201] + GHeq[M021];
        T GHeq_M120 = GHeq[M120] + GHeq[M102];
        T GHeq_M012 = GHeq[M012] - GHeq[M210];
        T GHeq_M021 = GHeq[M021] - GHeq[M201];
        T GHeq_M102 = GHeq[M102] - GHeq[M120];

        T GH_M210 = GH[M210] + GH[M012];
        T GH_M201 = GH[M201] + GH[M021];
        T GH_M120 = GH[M120] + GH[M102];
        T GH_M012 = GH[M012] - GH[M210];
        T GH_M021 = GH[M021] - GH[M201];
        T GH_M102 = GH[M102] - GH[M120];

        // 4 ortho with 2
        T GHeq_M220 = GHeq[M220] + (1. / 6.) * GHeq[M002];
        T GHeq_M202 = GHeq[M202] + (1. / 6.) * GHeq[M020];
        T GHeq_M022 = GHeq[M022] + (1. / 6.) * GHeq[M200];

        T GH_M220 = GH[M220] + (1. / 6.) * GH[M002];
        T GH_M202 = GH[M202] + (1. / 6.) * GH[M020];
        T GH_M022 = GH[M022] + (1. / 6.) * GH[M200];

        // 4 ortho with 4
        GHeq_M202 = GHeq[M202] + (2. / 7.) * GHeq_M220;
        GHeq_M022 = GHeq[M022] + (2. / 7.) * GHeq_M220 + (2. / 5.) * GHeq_M202;

        GH_M202 = GH[M202] + (2. / 7.) * GH_M220;
        GH_M022 = GH[M022] + (2. / 7.) * GH_M220 + (2. / 5.) * GH_M202;

        // Order 3
        GHcoll[M210] = (1. - omega3) * GH_M210 + omega3 * GHeq_M210;
        GHcoll[M201] = (1. - omega3) * GH_M201 + omega3 * GHeq_M201;
        GHcoll[M021] = (1. - omega3) * GH_M021 + omega3 * GHeq_M021;
        GHcoll[M120] = (1. - omega3) * GH_M120 + omega3 * GHeq_M120;
        GHcoll[M102] = (1. - omega3) * GH_M102 + omega3 * GHeq_M102;
        GHcoll[M012] = (1. - omega3) * GH_M012 + omega3 * GHeq_M012;

        // Order 4
        GHcoll[M220] = (1. - omega4) * GH_M220 + omega4 * GHeq_M220;
        GHcoll[M202] = (1. - omega4) * GH_M202 + omega4 * GHeq_M202;
        GHcoll[M022] = (1. - omega4) * GH_M022 + omega4 * GHeq_M022;

        cell[F000] = (rho * D::t[F000])
                     * (1. - 1.5 * (GHcoll[M200] + GHcoll[M020] + GHcoll[M002])
                        + (9. / 7.) * GHcoll[M220] + (9. / 5.) * GHcoll[M202] + 3. * GHcoll[M022]);

        cell[FP00] = (rho * D::t[FP00])
                     * (1. + 3. * u[0] + 1.5 * (2. * GHcoll[M200] - GHcoll[M020] - GHcoll[M002])
                        - 9. * GHcoll[M120] - (45. / 7.) * GHcoll[M220] - 9. * GHcoll[M202]);
        cell[FM00] = (rho * D::t[FM00])
                     * (1. - 3. * u[0] + 1.5 * (2. * GHcoll[M200] - GHcoll[M020] - GHcoll[M002])
                        + 9. * GHcoll[M120] - (45. / 7.) * GHcoll[M220] - 9. * GHcoll[M202]);

        cell[F0P0] = (rho * D::t[F0P0])
                     * (1. + 3. * u[1] + 1.5 * (2. * GHcoll[M020] - GHcoll[M200] - GHcoll[M002])
                        - 9. * GHcoll[M210] - (45. / 7.) * GHcoll[M220] + (18. / 5.) * GHcoll[M202]
                        - 9. * GHcoll[M022]);
        cell[F0M0] = (rho * D::t[F0M0])
                     * (1. - 3. * u[1] + 1.5 * (2. * GHcoll[M020] - GHcoll[M200] - GHcoll[M002])
                        + 9. * GHcoll[M210] - (45. / 7.) * GHcoll[M220] + (18. / 5.) * GHcoll[M202]
                        - 9. * GHcoll[M022]);

        cell[F00P] = (rho * D::t[F00P])
                     * (1. + 3. * u[2] + 1.5 * (2. * GHcoll[M002] - GHcoll[M200] - GHcoll[M020])
                        - 9. * GHcoll[M201] + (36. / 7.) * GHcoll[M220] - (27. / 5.) * GHcoll[M202]
                        - 9. * GHcoll[M022]);
        cell[F00M] = (rho * D::t[F00M])
                     * (1. - 3. * u[2] + 1.5 * (2. * GHcoll[M002] - GHcoll[M200] - GHcoll[M020])
                        + 9. * GHcoll[M201] + (36. / 7.) * GHcoll[M220] - (27. / 5.) * GHcoll[M202]
                        - 9. * GHcoll[M022]);

        cell[FPP0] =
            (rho * D::t[FPP0])
            * (1. + 3. * (+u[0] + u[1])
               + 1.5 * (2. * GHcoll[M200] + 2. * GHcoll[M020] - GHcoll[M002]) + 9. * GHcoll[M110]
               + 4.5 * (GHcoll[M210] + GHcoll[M120] - GHcoll[M102] - GHcoll[M012])
               + 9. * GHcoll[M220]);
        cell[FMP0] =
            (rho * D::t[FMP0])
            * (1. + 3. * (-u[0] + u[1])
               + 1.5 * (2. * GHcoll[M200] + 2. * GHcoll[M020] - GHcoll[M002]) - 9. * GHcoll[M110]
               + 4.5 * (GHcoll[M210] - GHcoll[M120] + GHcoll[M102] - GHcoll[M012])
               + 9. * GHcoll[M220]);
        cell[FPM0] =
            (rho * D::t[FPM0])
            * (1. + 3. * (+u[0] - u[1])
               + 1.5 * (2. * GHcoll[M200] + 2. * GHcoll[M020] - GHcoll[M002]) - 9. * GHcoll[M110]
               + 4.5 * (-GHcoll[M210] + GHcoll[M120] - GHcoll[M102] + GHcoll[M012])
               + 9. * GHcoll[M220]);
        cell[FMM0] =
            (rho * D::t[FMM0])
            * (1. + 3. * (-u[0] - u[1])
               + 1.5 * (2. * GHcoll[M200] + 2. * GHcoll[M020] - GHcoll[M002]) + 9. * GHcoll[M110]
               + 4.5 * (-GHcoll[M210] - GHcoll[M120] + GHcoll[M102] + GHcoll[M012])
               + 9. * GHcoll[M220]);

        cell[FP0P] =
            (rho * D::t[FP0P])
            * (1. + 3. * (+u[0] + u[2])
               + 1.5 * (2. * GHcoll[M200] + 2. * GHcoll[M002] - GHcoll[M020]) + 9. * GHcoll[M101]
               + 4.5 * (GHcoll[M201] + GHcoll[M102] + GHcoll[M120] - GHcoll[M021])
               - (18. / 7.) * GHcoll[M220] + 9. * GHcoll[M202]);
        cell[FM0P] =
            (rho * D::t[FM0P])
            * (1. + 3. * (-u[0] + u[2])
               + 1.5 * (2. * GHcoll[M200] + 2. * GHcoll[M002] - GHcoll[M020]) - 9. * GHcoll[M101]
               + 4.5 * (GHcoll[M201] - GHcoll[M102] - GHcoll[M120] - GHcoll[M021])
               - (18. / 7.) * GHcoll[M220] + 9. * GHcoll[M202]);
        cell[FP0M] =
            (rho * D::t[FP0M])
            * (1. + 3. * (+u[0] - u[2])
               + 1.5 * (2. * GHcoll[M200] + 2. * GHcoll[M002] - GHcoll[M020]) - 9. * GHcoll[M101]
               + 4.5 * (-GHcoll[M201] + GHcoll[M102] + GHcoll[M120] + GHcoll[M021])
               - (18. / 7.) * GHcoll[M220] + 9. * GHcoll[M202]);
        cell[FM0M] =
            (rho * D::t[FM0M])
            * (1. + 3. * (-u[0] - u[2])
               + 1.5 * (2. * GHcoll[M200] + 2. * GHcoll[M002] - GHcoll[M020]) + 9. * GHcoll[M101]
               + 4.5 * (-GHcoll[M201] - GHcoll[M102] - GHcoll[M120] + GHcoll[M021])
               - (18. / 7.) * GHcoll[M220] + 9. * GHcoll[M202]);

        cell[F0PP] =
            (rho * D::t[F0PP])
            * (1. + 3. * (+u[1] + u[2])
               + 1.5 * (2. * GHcoll[M020] + 2. * GHcoll[M002] - GHcoll[M200]) + 9. * GHcoll[M011]
               + 4.5 * (GHcoll[M021] + GHcoll[M012] + GHcoll[M210] + GHcoll[M201])
               - (18. / 7.) * GHcoll[M220] - (18. / 5.) * GHcoll[M202] + 9. * GHcoll[M022]);
        cell[F0MP] =
            (rho * D::t[F0MP])
            * (1. + 3. * (-u[1] + u[2])
               + 1.5 * (2. * GHcoll[M020] + 2. * GHcoll[M002] - GHcoll[M200]) - 9. * GHcoll[M011]
               + 4.5 * (GHcoll[M021] - GHcoll[M012] - GHcoll[M210] + GHcoll[M201])
               - (18. / 7.) * GHcoll[M220] - (18. / 5.) * GHcoll[M202] + 9. * GHcoll[M022]);
        cell[F0PM] =
            (rho * D::t[F0PM])
            * (1. + 3. * (+u[1] - u[2])
               + 1.5 * (2. * GHcoll[M020] + 2. * GHcoll[M002] - GHcoll[M200]) - 9. * GHcoll[M011]
               + 4.5 * (-GHcoll[M021] + GHcoll[M012] + GHcoll[M210] - GHcoll[M201])
               - (18. / 7.) * GHcoll[M220] - (18. / 5.) * GHcoll[M202] + 9. * GHcoll[M022]);
        cell[F0MM] =
            (rho * D::t[F0MM])
            * (1. + 3. * (-u[1] - u[2])
               + 1.5 * (2. * GHcoll[M020] + 2. * GHcoll[M002] - GHcoll[M200]) + 9. * GHcoll[M011]
               + 4.5 * (-GHcoll[M021] - GHcoll[M012] - GHcoll[M210] - GHcoll[M201])
               - (18. / 7.) * GHcoll[M220] - (18. / 5.) * GHcoll[M202] + 9. * GHcoll[M022]);

        for (int i = 0; i < 19; ++i) {
            cell[i] -= D::SkordosFactor() * D::t[i];
        }
    };

    //////////////////////////////////////////////////////////////////////////////////
    // Recursive Regularization (RR) approach based on the NON-WEIGHTED formulation //
    //////////////////////////////////////////////////////////////////////////////////

    /**
     * // General way to compute RRs (we only need 2nd-order RRs)
     * static void RRcomputeMoments(Array<T,D::q> const& cell, Array<T, D::q>& RR, T& rho) {
     *
     *     Array<T, D::q> f;
     *     for (int i = 0; i<19; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *         RR[i] = 0.;
     *     }
     *
     *     T Hxx = 0.;
     *     T Hyy = 0.;
     *     T Hzz = 0.;
     *
     *     for (int i = 0; i<19; ++i) {
     *
     *         Hxx = D::c[i][0] * D::c[i][0] - D::cs2;
     *         Hyy = D::c[i][1] * D::c[i][1] - D::cs2;
     *         Hzz = D::c[i][2] * D::c[i][2] - D::cs2;
     *
     *         // Order 0
     *         RR[M000] += f[i];
     *
     *         // Order 1
     *         RR[M100] += D::c[i][0] * f[i];
     *         RR[M010] += D::c[i][1] * f[i];
     *         RR[M001] += D::c[i][2] * f[i];
     *
     *         // Order 2
     *         RR[M200] += Hxx * f[i];
     *         RR[M020] += Hyy * f[i];
     *         RR[M002] += Hzz * f[i];
     *         RR[M110] += D::c[i][0] * D::c[i][1] * f[i];
     *         RR[M101] += D::c[i][0] * D::c[i][2] * f[i];
     *         RR[M011] += D::c[i][1] * D::c[i][2] * f[i];
     *     }
     *
     *     rho = RR[M000];
     *     T invRho = 1. / rho;
     *     for (int i = 0; i<19; ++i) {
     *         RR[i] *= invRho;
     *     }
     * }
     */

    // Optimized way to compute RRs based on the general ordering of discrete velocities
    static void RRcomputeMoments(Array<T, D::q> const &cell, Array<T, D::q> &RR, T &rho)
    {
        Array<T, D::q> f;
        for (int i = 0; i < 19; ++i) {
            f[i] = cell[i] + D::SkordosFactor() * D::t[i];
        }
        T X_M1 = f[FM00] + f[FMM0] + f[FMP0] + f[FM0M] + f[FM0P];
        T X_P1 = f[FP00] + f[FPP0] + f[FPM0] + f[FP0P] + f[FP0M];
        T X_0 =
            f[F000] + f[F0M0] + f[F00M] + f[F0MM] + f[F0MP] + f[F0P0] + f[F00P] + f[F0PP] + f[F0PM];
        T Y_M1 = f[F0M0] + f[FMM0] + f[F0MM] + f[F0MP] + f[FPM0];
        T Y_P1 = f[FMP0] + f[F0P0] + f[FPP0] + f[F0PP] + f[F0PM];
        T Z_M1 = f[F00M] + f[FM0M] + f[F0MM] + f[FP0M] + f[F0PM];
        T Z_P1 = f[FM0P] + f[F0MP] + f[F00P] + f[FP0P] + f[F0PP];

        rho = X_M1 + X_P1 + X_0;
        T invRho = 1. / rho;

        // Order 0
        RR[M000] = 1.;
        // Order 1
        RR[M100] = invRho * (X_P1 - X_M1);
        RR[M010] = invRho * (Y_P1 - Y_M1);
        RR[M001] = invRho * (Z_P1 - Z_M1);
        // Order 2
        RR[M200] = invRho * (X_M1 + X_P1);
        RR[M020] = invRho * (Y_M1 + Y_P1);
        RR[M002] = invRho * (Z_M1 + Z_P1);
        RR[M110] = invRho * (f[FMM0] - f[FMP0] + f[FPP0] - f[FPM0]);
        RR[M101] = invRho * (f[FM0M] - f[FM0P] + f[FP0P] - f[FP0M]);
        RR[M011] = invRho * (f[F0MM] - f[F0MP] + f[F0PP] - f[F0PM]);

        // We come back to Hermite moments
        RR[M200] -= D::cs2;
        RR[M020] -= D::cs2;
        RR[M002] -= D::cs2;
    }

    static void RRcomputeEquilibriumMoments(Array<T, D::d> const &u, Array<T, D::q> &RReq)
    {
        // Order 0
        RReq[M000] = 1.;
        // Order 1
        RReq[M100] = u[0];
        RReq[M010] = u[1];
        RReq[M001] = u[2];
        // Order 2
        RReq[M200] = u[0] * u[0];
        RReq[M020] = u[1] * u[1];
        RReq[M002] = u[2] * u[2];
        RReq[M110] = u[0] * u[1];
        RReq[M101] = u[0] * u[2];
        RReq[M011] = u[1] * u[2];
        // Order 3
        RReq[M210] = RReq[M200] * u[1];
        RReq[M201] = RReq[M200] * u[2];
        RReq[M021] = RReq[M020] * u[2];
        RReq[M120] = RReq[M020] * u[0];
        RReq[M102] = RReq[M002] * u[0];
        RReq[M012] = RReq[M002] * u[1];
        // Order 4
        RReq[M220] = RReq[M200] * RReq[M020];
        RReq[M202] = RReq[M200] * RReq[M002];
        RReq[M022] = RReq[M020] * RReq[M002];
    };

    // Equilibrium populations based on 19 moments can be computed using either RM, HM, CM, CHM or
    // Gauss-Hermite formalisms. Here we use Hermite moments (RRs)
    static void RRcomputeEquilibrium(T rho, Array<T, D::q> const &RReq, Array<T, D::q> &eq)
    {
        Array<T, D::d> u(RReq[1], RReq[2], RReq[3]);
        Array<T, D::q> RMeq;
        // Order 2
        RMeq[M200] = u[0] * u[0] + D::cs2;
        RMeq[M020] = u[1] * u[1] + D::cs2;
        RMeq[M002] = u[2] * u[2] + D::cs2;
        RMeq[M110] = u[0] * u[1];
        RMeq[M101] = u[0] * u[2];
        RMeq[M011] = u[1] * u[2];
        // Order 3
        RMeq[M210] = RMeq[M200] * u[1];
        RMeq[M201] = RMeq[M200] * u[2];
        RMeq[M021] = RMeq[M020] * u[2];
        RMeq[M120] = RMeq[M020] * u[0];
        RMeq[M102] = RMeq[M002] * u[0];
        RMeq[M012] = RMeq[M002] * u[1];
        // Order 4
        RMeq[M220] = RMeq[M200] * RMeq[M020];
        RMeq[M202] = RMeq[M200] * RMeq[M002];
        RMeq[M022] = RMeq[M020] * RMeq[M002];

        eq[F000] =
            rho
            * (1. - RMeq[M200] - RMeq[M020] - RMeq[M002] + RMeq[M220] + RMeq[M202] + RMeq[M022]);

        eq[FP00] =
            0.5 * rho * (u[0] + RMeq[M200] - RMeq[M120] - RMeq[M102] - RMeq[M220] - RMeq[M202]);
        eq[FM00] =
            0.5 * rho * (-u[0] + RMeq[M200] + RMeq[M120] + RMeq[M102] - RMeq[M220] - RMeq[M202]);

        eq[F0P0] =
            0.5 * rho * (u[1] + RMeq[M020] - RMeq[M210] - RMeq[M012] - RMeq[M220] - RMeq[M022]);
        eq[F0M0] =
            0.5 * rho * (-u[1] + RMeq[M020] + RMeq[M210] + RMeq[M012] - RMeq[M220] - RMeq[M022]);

        eq[F00P] =
            0.5 * rho * (u[2] + RMeq[M002] - RMeq[M201] - RMeq[M021] - RMeq[M202] - RMeq[M022]);
        eq[F00M] =
            0.5 * rho * (-u[2] + RMeq[M002] + RMeq[M201] + RMeq[M021] - RMeq[M202] - RMeq[M022]);

        eq[FPP0] = 0.25 * rho * (RMeq[M110] + RMeq[M210] + RMeq[M120] + RMeq[M220]);
        eq[FMP0] = 0.25 * rho * (-RMeq[M110] + RMeq[M210] - RMeq[M120] + RMeq[M220]);
        eq[FPM0] = 0.25 * rho * (-RMeq[M110] - RMeq[M210] + RMeq[M120] + RMeq[M220]);
        eq[FMM0] = 0.25 * rho * (RMeq[M110] - RMeq[M210] - RMeq[M120] + RMeq[M220]);

        eq[FP0P] = 0.25 * rho * (RMeq[M101] + RMeq[M201] + RMeq[M102] + RMeq[M202]);
        eq[FM0P] = 0.25 * rho * (-RMeq[M101] + RMeq[M201] - RMeq[M102] + RMeq[M202]);
        eq[FP0M] = 0.25 * rho * (-RMeq[M101] - RMeq[M201] + RMeq[M102] + RMeq[M202]);
        eq[FM0M] = 0.25 * rho * (RMeq[M101] - RMeq[M201] - RMeq[M102] + RMeq[M202]);

        eq[F0PP] = 0.25 * rho * (RMeq[M011] + RMeq[M021] + RMeq[M012] + RMeq[M022]);
        eq[F0MP] = 0.25 * rho * (-RMeq[M011] + RMeq[M021] - RMeq[M012] + RMeq[M022]);
        eq[F0PM] = 0.25 * rho * (-RMeq[M011] - RMeq[M021] + RMeq[M012] + RMeq[M022]);
        eq[F0MM] = 0.25 * rho * (RMeq[M011] - RMeq[M021] - RMeq[M012] + RMeq[M022]);
    };

    static void RRcollide(
        Array<T, D::q> &cell, T rho, Array<T, D::d> const &u,
        Array<T, D::q> const &RR,    // Hermite moments
        Array<T, D::q> const &RReq,  // Equilibrium moments (Hermite)
        Array<T, D::numRelaxationTimes> const &omega)
    {
        T omega1 = omega[0];
        T omega2 = omega[1];
        T omega3 = omega[2];
        T omega4 = omega[3];
        T omegaBulk = omega[4];
        T omegaPlus = (omegaBulk + 2. * omega1) / 3.;  // Notation used by Fei
        T omegaMinus = (omegaBulk - omega1) / 3.;      // Notation used by Fei

        T cs4 = D::cs2 * D::cs2;

        // Post-collision and Nonequilibrium moments.
        Array<T, D::q> RRneq;
        Array<T, D::q> RRcoll;
        Array<T, D::q> RMcoll;

        // Recursive computation of nonequilibrium Hermite moments
        // Order 2 (standard way to compute them)
        RRneq[M200] = RR[M200] - RReq[M200];
        RRneq[M020] = RR[M020] - RReq[M020];
        RRneq[M002] = RR[M002] - RReq[M002];
        RRneq[M110] = RR[M110] - RReq[M110];
        RRneq[M101] = RR[M101] - RReq[M101];
        RRneq[M011] = RR[M011] - RReq[M011];

        // Order 3 (reconstruction using Chapman-Enskog formulas)
        RRneq[M210] = u[1] * RRneq[M200] + 2. * u[0] * RRneq[M110];
        RRneq[M201] = u[2] * RRneq[M200] + 2. * u[0] * RRneq[M101];
        RRneq[M021] = u[2] * RRneq[M020] + 2. * u[1] * RRneq[M011];
        RRneq[M120] = u[0] * RRneq[M020] + 2. * u[1] * RRneq[M110];
        RRneq[M102] = u[0] * RRneq[M002] + 2. * u[2] * RRneq[M101];
        RRneq[M012] = u[1] * RRneq[M002] + 2. * u[2] * RRneq[M011];

        // Order 4 (reconstruction using Chapman-Enskog formulas)
        RRneq[M220] =
            u[1] * u[1] * RRneq[M200] + u[0] * u[0] * RRneq[M020] + 4. * u[0] * u[1] * RRneq[M110];
        RRneq[M202] =
            u[2] * u[2] * RRneq[M200] + u[0] * u[0] * RRneq[M002] + 4. * u[0] * u[2] * RRneq[M101];
        RRneq[M022] =
            u[2] * u[2] * RRneq[M020] + u[1] * u[1] * RRneq[M002] + 4. * u[1] * u[2] * RRneq[M011];

        // Collision in the Hermite moment space
        // Order 2 (non-diagonal collision so that we can easily modify the bulk viscosity)
        RRcoll[M200] = RR[M200] - omegaPlus * RRneq[M200] - omegaMinus * RRneq[M020]
                       - omegaMinus * RRneq[M002];
        RRcoll[M020] = RR[M020] - omegaMinus * RRneq[M200] - omegaPlus * RRneq[M020]
                       - omegaMinus * RRneq[M002];
        RRcoll[M002] = RR[M002] - omegaMinus * RRneq[M200] - omegaMinus * RRneq[M020]
                       - omegaPlus * RRneq[M002];

        RRcoll[M110] = (1. - omega2) * RRneq[M110] + RReq[M110];
        RRcoll[M101] = (1. - omega2) * RRneq[M101] + RReq[M101];
        RRcoll[M011] = (1. - omega2) * RRneq[M011] + RReq[M011];

        // Order 3
        RRcoll[M210] = (1. - omega3) * RRneq[M210] + RReq[M210];
        RRcoll[M201] = (1. - omega3) * RRneq[M201] + RReq[M201];
        RRcoll[M021] = (1. - omega3) * RRneq[M021] + RReq[M021];
        RRcoll[M120] = (1. - omega3) * RRneq[M120] + RReq[M120];
        RRcoll[M102] = (1. - omega3) * RRneq[M102] + RReq[M102];
        RRcoll[M012] = (1. - omega3) * RRneq[M012] + RReq[M012];

        // Order 4
        RRcoll[M220] = (1. - omega4) * RRneq[M220] + RReq[M220];
        RRcoll[M202] = (1. - omega4) * RRneq[M202] + RReq[M202];
        RRcoll[M022] = (1. - omega4) * RRneq[M022] + RReq[M022];

        // Come back to RMcoll using relationships between RRs and RMs
        RMcoll[M200] = RRcoll[M200] + D::cs2;
        RMcoll[M020] = RRcoll[M020] + D::cs2;
        RMcoll[M002] = RRcoll[M002] + D::cs2;

        RMcoll[M110] = RRcoll[M110];
        RMcoll[M101] = RRcoll[M101];
        RMcoll[M011] = RRcoll[M011];

        RMcoll[M210] = RRcoll[M210] + D::cs2 * u[1];
        RMcoll[M201] = RRcoll[M201] + D::cs2 * u[2];
        RMcoll[M021] = RRcoll[M021] + D::cs2 * u[2];
        RMcoll[M120] = RRcoll[M120] + D::cs2 * u[0];
        RMcoll[M102] = RRcoll[M102] + D::cs2 * u[0];
        RMcoll[M012] = RRcoll[M012] + D::cs2 * u[1];

        RMcoll[M220] = RRcoll[M220] + D::cs2 * (RRcoll[M200] + RRcoll[M020]) + cs4;
        RMcoll[M202] = RRcoll[M202] + D::cs2 * (RRcoll[M200] + RRcoll[M002]) + cs4;
        RMcoll[M022] = RRcoll[M022] + D::cs2 * (RRcoll[M020] + RRcoll[M002]) + cs4;

        // Compute post collision populations from RM
        // Optimization based on symmetries between populations and their opposite counterpart
        cell[F000] = rho
                     * (1. - RMcoll[M200] - RMcoll[M020] - RMcoll[M002] + RMcoll[M220]
                        + RMcoll[M202] + RMcoll[M022]);

        cell[FP00] =
            0.5 * rho
            * (u[0] + RMcoll[M200] - RMcoll[M120] - RMcoll[M102] - RMcoll[M220] - RMcoll[M202]);
        cell[FM00] = rho * (-u[0] + RMcoll[M120] + RMcoll[M102]) + cell[FP00];

        cell[F0P0] =
            0.5 * rho
            * (u[1] + RMcoll[M020] - RMcoll[M210] - RMcoll[M012] - RMcoll[M220] - RMcoll[M022]);
        cell[F0M0] = rho * (-u[1] + RMcoll[M210] + RMcoll[M012]) + cell[F0P0];

        cell[F00P] =
            0.5 * rho
            * (u[2] + RMcoll[M002] - RMcoll[M201] - RMcoll[M021] - RMcoll[M202] - RMcoll[M022]);
        cell[F00M] = rho * (-u[2] + RMcoll[M201] + RMcoll[M021]) + cell[F00P];

        cell[FPP0] = 0.25 * rho * (RMcoll[M110] + RMcoll[M210] + RMcoll[M120] + RMcoll[M220]);
        cell[FMP0] = 0.5 * rho * (-RMcoll[M110] - RMcoll[M120]) + cell[FPP0];
        cell[FPM0] = 0.5 * rho * (-RMcoll[M110] - RMcoll[M210]) + cell[FPP0];
        cell[FMM0] = 0.5 * rho * (-RMcoll[M210] - RMcoll[M120]) + cell[FPP0];

        cell[FP0P] = 0.25 * rho * (RMcoll[M101] + RMcoll[M201] + RMcoll[M102] + RMcoll[M202]);
        cell[FM0P] = 0.5 * rho * (-RMcoll[M101] - RMcoll[M102]) + cell[FP0P];
        cell[FP0M] = 0.5 * rho * (-RMcoll[M101] - RMcoll[M201]) + cell[FP0P];
        cell[FM0M] = 0.5 * rho * (-RMcoll[M201] - RMcoll[M102]) + cell[FP0P];

        cell[F0PP] = 0.25 * rho * (RMcoll[M011] + RMcoll[M021] + RMcoll[M012] + RMcoll[M022]);
        cell[F0MP] = 0.5 * rho * (-RMcoll[M011] - RMcoll[M012]) + cell[F0PP];
        cell[F0PM] = 0.5 * rho * (-RMcoll[M011] - RMcoll[M021]) + cell[F0PP];
        cell[F0MM] = 0.5 * rho * (-RMcoll[M021] - RMcoll[M012]) + cell[F0PP];

        for (int i = 0; i < 19; ++i) {
            cell[i] -= D::SkordosFactor() * D::t[i];
        }
    };

};  // struct comprehensiveDynamicsTemplatesImpl<T, descriptors::D3Q19DescriptorBase<T> >

}  // namespace plb

#endif  // COMPREHENSIVE_MODELS_TEMPLATES_3D_H
