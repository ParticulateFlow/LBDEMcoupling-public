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
 * 2D specialization of dynamicsTemplates functions.
 *
 * Theoretical background about these collision models can be found in
 * Coreixas et al. 'Comprehensive comparison of collision models in the
 * lattice Boltzmann framework: Theoretical investigations', PRE, 2019.
 *
 * Stability domains and preliminary accuracy comparisons are available in
 * Coreixas et al. 'Impact of collision models on the physical properties
 * and the stability of lattice Boltzmann methods', PTRSA, 2020.
 */

#ifndef COMPREHENSIVE_MODELS_TEMPLATES_2D_H
#define COMPREHENSIVE_MODELS_TEMPLATES_2D_H

#include "core/globalDefs.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"

namespace plb {
// Efficient specialization for D2Q9 lattice
template <typename T>
struct comprehensiveDynamicsTemplatesImpl<T, descriptors::D2Q9DescriptorBase<T> > {
    typedef descriptors::D2Q9DescriptorBase<T> D;
    // Same order as in E1.
    enum {
        // Order 0
        M00 = 0,

        // Order 1
        M10 = 1,
        M01 = 2,

        // Order 2
        M20 = 3,
        M02 = 4,
        M11 = 5,

        // Order 3
        M21 = 6,
        M12 = 7,

        // Order 4
        M22 = 8,
    };

    /**
     * // General way to compute RMs
     * static void RMcomputeMoments(Array<T,D::q> const& cell, Array<T, D::q>& RM, T& rho) {
     *
     *     Array<T, D::q> f;
     *     for (int i = 0; i<D::q; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *         RM[i] = 0.;
     *     }
     *
     *     for (int i = 0; i<D::q; ++i) {
     *         // Order 0
     *         RM[M00] += f[i];
     *
     *         // Order 1
     *         RM[M10] += D::c[i][0] * f[i];
     *         RM[M01] += D::c[i][1] * f[i];
     *
     *         // Order 2
     *         RM[M20] += D::c[i][0] * D::c[i][0] * f[i];
     *         RM[M02] += D::c[i][1] * D::c[i][1] * f[i];
     *         RM[M11] += D::c[i][0] * D::c[i][1] * f[i];
     *
     *         // Order 3
     *         RM[M21] += D::c[i][0] * D::c[i][0] * D::c[i][1] * f[i];
     *         RM[M12] += D::c[i][0] * D::c[i][1] * D::c[i][1] * f[i];
     *
     *         // Order 4
     *         RM[M22] += D::c[i][0] * D::c[i][0] * D::c[i][1] * D::c[i][1] * f[i];
     *     }
     *
     *     rho = RM[M00];
     *     T invRho = 1. / rho;
     *     for (int i = 0; i<D::q; ++i) {
     *         RM[i] *= invRho;
     *     }
     * };
     */

    // Optimized way to compute RMs based on Palabos ordering of discrete velocities
    static void RMcomputeMoments(Array<T, D::q> const &cell, Array<T, D::q> &RM, T &rho)
    {
        Array<T, D::q> f;
        for (int i = 0; i < D::q; ++i) {
            f[i] = cell[i] + D::SkordosFactor() * D::t[i];
            RM[i] = 0.;
        }

        double X_M1 = f[1] + f[2] + f[3];
        double X_P1 = f[5] + f[6] + f[7];
        double X_0 = f[0] + f[4] + f[8];

        double Y_M1 = f[3] + f[4] + f[5];
        double Y_P1 = f[1] + f[7] + f[8];

        // Order 0
        RM[M00] = X_M1 + X_P1 + X_0;

        // Order 1
        RM[M10] = X_P1 - X_M1;
        RM[M01] = Y_P1 - Y_M1;

        // Order 2
        RM[M20] = X_P1 + X_M1;
        RM[M02] = Y_P1 + Y_M1;
        RM[M11] = -f[1] + f[3] - f[5] + f[7];

        // Order 3
        RM[M21] = f[1] - f[3] - f[5] + f[7];
        RM[M12] = -f[1] - f[3] + f[5] + f[7];

        // Order 4
        RM[M22] = f[1] + f[3] + f[5] + f[7];

        rho = RM[M00];
        T invRho = 1. / rho;
        for (int i = 0; i < D::q; ++i) {
            RM[i] *= invRho;
        }
    };

    static void RMcomputeEquilibriumMoments(Array<T, D::d> const &u, Array<T, D::q> &RMeq)
    {
        // Order 0
        RMeq[M00] = 1.;

        // Order 1
        RMeq[M10] = u[0];
        RMeq[M01] = u[1];

        // Order 2
        RMeq[M20] = u[0] * u[0] + D::cs2;
        RMeq[M02] = u[1] * u[1] + D::cs2;
        RMeq[M11] = u[0] * u[1];

        // Order 3
        RMeq[M21] = RMeq[M20] * u[1];
        RMeq[M12] = RMeq[M02] * u[0];

        // Order 4
        RMeq[M22] = RMeq[M20] * RMeq[M02];
    };

    enum { F00 = 0, FMP = 1, FM0 = 2, FMM = 3, F0M = 4, FPM = 5, FP0 = 6, FPP = 7, F0P = 8 };

    // Equilibrium populations based on 9 moments can be computed using either RM, HM, CM, CHM or
    // Gauss-Hermite formalisms. Here we use raw moments (RMs)
    static void RMcomputeEquilibrium(T rho, Array<T, D::q> const &RMeq, Array<T, D::q> &eq)
    {
        Array<T, D::d> u(RMeq[1], RMeq[2]);
        eq[F00] = rho * (1. - RMeq[M20] - RMeq[M02] + RMeq[M22]);

        eq[FP0] = 0.5 * rho * (u[0] + RMeq[M20] - RMeq[M12] - RMeq[M22]);
        eq[FM0] = 0.5 * rho * (-u[0] + RMeq[M20] + RMeq[M12] - RMeq[M22]);

        eq[F0P] = 0.5 * rho * (u[1] + RMeq[M02] - RMeq[M21] - RMeq[M22]);
        eq[F0M] = 0.5 * rho * (-u[1] + RMeq[M02] + RMeq[M21] - RMeq[M22]);

        eq[FPP] = 0.25 * rho * (RMeq[M11] + RMeq[M21] + RMeq[M12] + RMeq[M22]);
        eq[FMP] = 0.25 * rho * (-RMeq[M11] + RMeq[M21] - RMeq[M12] + RMeq[M22]);
        eq[FPM] = 0.25 * rho * (-RMeq[M11] - RMeq[M21] + RMeq[M12] + RMeq[M22]);
        eq[FMM] = 0.25 * rho * (RMeq[M11] - RMeq[M21] - RMeq[M12] + RMeq[M22]);
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
        T omegaPlus = (omegaBulk + omega1) / 2.;   // Notation used by Fei
        T omegaMinus = (omegaBulk - omega1) / 2.;  // Notation used by Fei

        // Post-collision moments.
        Array<T, D::q> RMcoll;

        // Order 2 (non-diagonal collision so that we can easily modify the bulk viscosity)
        RMcoll[M20] =
            RM[M20] - omegaPlus * (RM[M20] - RMeq[M20]) - omegaMinus * (RM[M02] - RMeq[M02]);
        RMcoll[M02] =
            RM[M02] - omegaMinus * (RM[M20] - RMeq[M20]) - omegaPlus * (RM[M02] - RMeq[M02]);

        RMcoll[M11] = (1. - omega2) * RM[M11] + omega2 * RMeq[M11];

        // Order 3
        RMcoll[M21] = (1. - omega3) * RM[M21] + omega3 * RMeq[M21];
        RMcoll[M12] = (1. - omega3) * RM[M12] + omega3 * RMeq[M12];

        // Order 4
        RMcoll[M22] = (1. - omega4) * RM[M22] + omega4 * RMeq[M22];

        cell[F00] = rho * (1. - RMcoll[M20] - RMcoll[M02] + RMcoll[M22]);

        cell[FP0] = 0.5 * rho * (u[0] + RMcoll[M20] - RMcoll[M12] - RMcoll[M22]);
        cell[FM0] = 0.5 * rho * (-u[0] + RMcoll[M20] + RMcoll[M12] - RMcoll[M22]);

        cell[F0P] = 0.5 * rho * (u[1] + RMcoll[M02] - RMcoll[M21] - RMcoll[M22]);
        cell[F0M] = 0.5 * rho * (-u[1] + RMcoll[M02] + RMcoll[M21] - RMcoll[M22]);

        cell[FPP] = 0.25 * rho * (RMcoll[M11] + RMcoll[M21] + RMcoll[M12] + RMcoll[M22]);
        cell[FMP] = 0.25 * rho * (-RMcoll[M11] + RMcoll[M21] - RMcoll[M12] + RMcoll[M22]);
        cell[FPM] = 0.25 * rho * (-RMcoll[M11] - RMcoll[M21] + RMcoll[M12] + RMcoll[M22]);
        cell[FMM] = 0.25 * rho * (RMcoll[M11] - RMcoll[M21] - RMcoll[M12] + RMcoll[M22]);
        for (int i = 0; i < D::q; ++i) {
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
     *     for (int i = 0; i<D::q; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *         HM[i] = 0.;
     *     }
     *
     *     T Hxx = 0.;
     *     T Hyy = 0.;
     *
     *     for (int i = 0; i<D::q; ++i) {
     *
     *         Hxx = D::c[i][0] * D::c[i][0] - D::cs2;
     *         Hyy = D::c[i][1] * D::c[i][1] - D::cs2;
     *
     *         // Order 0
     *         HM[M00] += f[i];
     *
     *         // Order 1
     *         HM[M10] += D::c[i][0] * f[i];
     *         HM[M01] += D::c[i][1] * f[i];
     *
     *         // Order 2
     *         HM[M20] += Hxx * f[i];
     *         HM[M02] += Hyy * f[i];
     *         HM[M11] += D::c[i][0] * D::c[i][1] * f[i];
     *
     *         // Order 3
     *         HM[M21] += Hxx * D::c[i][1] * f[i];
     *         HM[M12] += D::c[i][0] * Hyy * f[i];
     *
     *         // Order 4
     *         HM[M22] += Hxx * Hyy * f[i];
     *     }
     *
     *     rho = HM[M00];
     *     T invRho = 1. / rho;
     *     for (int i = 0; i<D::q; ++i) {
     *         HM[i] *= invRho;
     *     }
     * };
     */

    // Optimized way to compute HMs based on Palabos ordering of discrete velocities
    static void HMcomputeMoments(Array<T, D::q> const &cell, Array<T, D::q> &HM, T &rho)
    {
        Array<T, D::q> f;
        for (int i = 0; i < D::q; ++i) {
            f[i] = cell[i] + D::SkordosFactor() * D::t[i];
            HM[i] = 0.;
        }

        double X_M1 = f[1] + f[2] + f[3];
        double X_P1 = f[5] + f[6] + f[7];
        double X_0 = f[0] + f[4] + f[8];

        double Y_M1 = f[3] + f[4] + f[5];
        double Y_P1 = f[1] + f[7] + f[8];

        // Order 0
        HM[M00] = X_M1 + X_P1 + X_0;

        // Order 1
        HM[M10] = X_P1 - X_M1;
        HM[M01] = Y_P1 - Y_M1;

        // Order 2
        HM[M20] = X_P1 + X_M1;
        HM[M02] = Y_P1 + Y_M1;
        HM[M11] = -f[1] + f[3] - f[5] + f[7];

        // Order 3
        HM[M21] = f[1] - f[3] - f[5] + f[7];
        HM[M12] = -f[1] - f[3] + f[5] + f[7];

        // Order 4
        HM[M22] = f[1] + f[3] + f[5] + f[7];

        rho = HM[M00];
        T invRho = 1. / rho;
        for (int i = 0; i < D::q; ++i) {
            HM[i] *= invRho;
        }

        // We come back to Hermite moments
        T cs4 = D::cs2 * D::cs2;
        HM[M20] -= D::cs2;
        HM[M02] -= D::cs2;

        HM[M21] -= D::cs2 * HM[M01];
        HM[M12] -= D::cs2 * HM[M10];

        HM[M22] -= (D::cs2 * (HM[M20] + HM[M02]) + cs4);
    };

    static void HMcomputeEquilibriumMoments(Array<T, D::d> const &u, Array<T, D::q> &HMeq)
    {
        // Order 0
        HMeq[M00] = 1.;

        // Order 1
        HMeq[M10] = u[0];
        HMeq[M01] = u[1];

        // Order 2
        HMeq[M20] = u[0] * u[0];
        HMeq[M02] = u[1] * u[1];
        HMeq[M11] = u[0] * u[1];

        // Order 3
        HMeq[M21] = HMeq[M20] * u[1];
        HMeq[M12] = HMeq[M02] * u[0];

        // Order 4
        HMeq[M22] = HMeq[M20] * HMeq[M02];
    };

    // Equilibrium populations based on 9 moments can be computed using either RM, HM, CM, CHM or
    // Gauss-Hermite formalisms. Here we use hermite moments (HMs)
    static void HMcomputeEquilibrium(T rho, Array<T, D::q> const &HMeq, Array<T, D::q> &eq)
    {
        Array<T, D::d> u(HMeq[1], HMeq[2]);
        Array<T, D::q> RMeq;
        T cs4 = D::cs2 * D::cs2;

        RMeq[M20] = HMeq[M20] + D::cs2;
        RMeq[M02] = HMeq[M02] + D::cs2;

        RMeq[M11] = HMeq[M11];

        RMeq[M21] = HMeq[M21] + D::cs2 * u[1];
        RMeq[M12] = HMeq[M12] + D::cs2 * u[0];

        RMeq[M22] = HMeq[M22] + D::cs2 * (HMeq[M20] + HMeq[M02]) + cs4;

        eq[F00] = rho * (1. - RMeq[M20] - RMeq[M02] + RMeq[M22]);

        eq[FP0] = 0.5 * rho * (u[0] + RMeq[M20] - RMeq[M12] - RMeq[M22]);
        eq[FM0] = 0.5 * rho * (-u[0] + RMeq[M20] + RMeq[M12] - RMeq[M22]);

        eq[F0P] = 0.5 * rho * (u[1] + RMeq[M02] - RMeq[M21] - RMeq[M22]);
        eq[F0M] = 0.5 * rho * (-u[1] + RMeq[M02] + RMeq[M21] - RMeq[M22]);

        eq[FPP] = 0.25 * rho * (RMeq[M11] + RMeq[M21] + RMeq[M12] + RMeq[M22]);
        eq[FMP] = 0.25 * rho * (-RMeq[M11] + RMeq[M21] - RMeq[M12] + RMeq[M22]);
        eq[FPM] = 0.25 * rho * (-RMeq[M11] - RMeq[M21] + RMeq[M12] + RMeq[M22]);
        eq[FMM] = 0.25 * rho * (RMeq[M11] - RMeq[M21] - RMeq[M12] + RMeq[M22]);
    };

    static void HMcollide(
        Array<T, D::q> &cell, T rho, Array<T, D::d> const &u,
        Array<T, D::q> const &HM,    // Hermite moments
        Array<T, D::q> const &HMeq,  // Equilibrium moments (hermite)
        Array<T, D::numRelaxationTimes> const &omega)
    {
        T omega1 = omega[0];
        T omega2 = omega[1];
        T omega3 = omega[2];
        T omega4 = omega[3];

        T omegaBulk = omega[4];
        T omegaPlus = (omegaBulk + omega1) / 2.;   // Notation used by Fei
        T omegaMinus = (omegaBulk - omega1) / 2.;  // Notation used by Fei

        T cs4 = D::cs2 * D::cs2;

        // Post-collision moments.
        Array<T, D::q> HMcoll;
        Array<T, D::q> RMcoll;

        // Collision in the Hermite moment space
        // Order 2 (non-diagonal collision so that we can easily modify the bulk viscosity)
        HMcoll[M20] =
            HM[M20] - omegaPlus * (HM[M20] - HMeq[M20]) - omegaMinus * (HM[M02] - HMeq[M02]);
        HMcoll[M02] =
            HM[M02] - omegaMinus * (HM[M20] - HMeq[M20]) - omegaPlus * (HM[M02] - HMeq[M02]);

        HMcoll[M11] = (1. - omega2) * HM[M11] + omega2 * HMeq[M11];

        // Order 3
        HMcoll[M21] = (1. - omega3) * HM[M21] + omega3 * HMeq[M21];
        HMcoll[M12] = (1. - omega3) * HM[M12] + omega3 * HMeq[M12];

        // Order 4
        HMcoll[M22] = (1. - omega4) * HM[M22] + omega4 * HMeq[M22];

        // Come back to RMcoll using relationships between HMs and RMs
        RMcoll[M20] = HMcoll[M20] + D::cs2;
        RMcoll[M02] = HMcoll[M02] + D::cs2;

        RMcoll[M11] = HMcoll[M11];

        RMcoll[M21] = HMcoll[M21] + D::cs2 * u[1];
        RMcoll[M12] = HMcoll[M12] + D::cs2 * u[0];

        RMcoll[M22] = HMcoll[M22] + D::cs2 * (HMcoll[M20] + HMcoll[M02]) + cs4;

        // Compute post collision populations from RM
        cell[F00] = rho * (1. - RMcoll[M20] - RMcoll[M02] + RMcoll[M22]);

        cell[FP0] = 0.5 * rho * (u[0] + RMcoll[M20] - RMcoll[M12] - RMcoll[M22]);
        cell[FM0] = 0.5 * rho * (-u[0] + RMcoll[M20] + RMcoll[M12] - RMcoll[M22]);

        cell[F0P] = 0.5 * rho * (u[1] + RMcoll[M02] - RMcoll[M21] - RMcoll[M22]);
        cell[F0M] = 0.5 * rho * (-u[1] + RMcoll[M02] + RMcoll[M21] - RMcoll[M22]);

        cell[FPP] = 0.25 * rho * (RMcoll[M11] + RMcoll[M21] + RMcoll[M12] + RMcoll[M22]);
        cell[FMP] = 0.25 * rho * (-RMcoll[M11] + RMcoll[M21] - RMcoll[M12] + RMcoll[M22]);
        cell[FPM] = 0.25 * rho * (-RMcoll[M11] - RMcoll[M21] + RMcoll[M12] + RMcoll[M22]);
        cell[FMM] = 0.25 * rho * (RMcoll[M11] - RMcoll[M21] - RMcoll[M12] + RMcoll[M22]);

        for (int i = 0; i < D::q; ++i) {
            cell[i] -= D::SkordosFactor() * D::t[i];
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////
    // Central Moments Formalism (Equilibrium is computed through raw moments formalism) //
    ///////////////////////////////////////////////////////////////////////////////////////

    /**
     * // General way to compute CMs (computations of velocity components are based on Palabos
     * ordering of discrete velocities !) static void CMcomputeMoments(Array<T,D::q> const& cell,
     * Array<T, D::q>& CM, T& rho, Array<T,D::d>& u) {
     *
     *     Array<T, D::q> f;
     *     for (int i = 0; i<D::q; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *         CM[i] = 0.;
     *     }
     *
     *     rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
     *     CM[M00] = 1.0;
     *     T invRho = 1./rho;
     *
     *     u[0] = invRho * ( - f[1] - f[2] - f[3] + f[5] + f[6] + f[7]);
     *     u[1] = invRho * (   f[1] - f[3] - f[4] - f[5] + f[7] + f[8]);
     *
     *     T cMux = 0.;
     *     T cMuy = 0.;
     *
     *     for (int i = 0; i<D::q; ++i) {
     *
     *         cMux = D::c[i][0]- u[0];
     *         cMuy = D::c[i][1]- u[1];
     *
     *         // Order 1
     *         CM[M10] += cMux * f[i];
     *         CM[M01] += cMuy * f[i];
     *
     *         // Order 2
     *         CM[M20] += cMux * cMux * f[i];
     *         CM[M02] += cMuy * cMuy * f[i];
     *         CM[M11] += cMux * cMuy * f[i];
     *
     *
     *         // Order 3
     *         CM[M21] += cMux * cMux * cMuy * f[i];
     *         CM[M12] += cMux * cMuy * cMuy * f[i];
     *
     *         // Order 4
     *         CM[M22] += cMux * cMux * cMuy * cMuy * f[i];
     *     }
     *
     *     for (int i = 1; i<D::q; ++i) {
     *         CM[i] *= invRho;
     *     }
     * };
     */

    // Optimized way to compute CMs based on Palabos ordering of discrete velocities
    static void CMcomputeMoments(
        Array<T, D::q> const &cell, Array<T, D::q> &CM, T &rho, Array<T, D::d> &u)
    {
        Array<T, D::q> f;
        for (int i = 0; i < D::q; ++i) {
            f[i] = cell[i] + D::SkordosFactor() * D::t[i];
            CM[i] = 0.;
        }

        double X_M1 = f[1] + f[2] + f[3];
        double X_P1 = f[5] + f[6] + f[7];
        double X_0 = f[0] + f[4] + f[8];

        double Y_M1 = f[3] + f[4] + f[5];
        double Y_P1 = f[1] + f[7] + f[8];

        // Order 0
        CM[M00] = X_M1 + X_P1 + X_0;

        // Order 1
        CM[M10] = 0.;
        CM[M01] = 0.;

        // Order 2
        CM[M20] = X_P1 + X_M1;
        CM[M02] = Y_P1 + Y_M1;
        CM[M11] = -f[1] + f[3] - f[5] + f[7];

        // Order 3
        CM[M21] = f[1] - f[3] - f[5] + f[7];
        CM[M12] = -f[1] - f[3] + f[5] + f[7];

        // Order 4
        CM[M22] = f[1] + f[3] + f[5] + f[7];

        rho = CM[M00];
        T invRho = 1. / rho;
        for (int i = 0; i < D::q; ++i) {
            CM[i] *= invRho;
        }

        // Compute CMs from RMs using binomial formulas
        u[0] = invRho * (X_P1 - X_M1);
        u[1] = invRho * (Y_P1 - Y_M1);
        T ux2 = u[0] * u[0];
        T uy2 = u[1] * u[1];

        CM[M20] -= (ux2);
        CM[M02] -= (uy2);
        CM[M11] -= (u[0] * u[1]);

        CM[M21] -= (u[1] * CM[M20] + 2. * u[0] * CM[M11] + ux2 * u[1]);
        CM[M12] -= (u[0] * CM[M02] + 2. * u[1] * CM[M11] + u[0] * uy2);

        CM[M22] -=
            (2. * u[1] * CM[M21] + 2. * u[0] * CM[M12] + uy2 * CM[M20] + ux2 * CM[M02]
             + 4. * u[0] * u[1] * CM[M11] + ux2 * uy2);
    };

    static void CMcomputeEquilibriumMoments(Array<T, D::q> &CMeq)
    {
        // Order 0
        CMeq[M00] = 1.;

        // Order 1
        CMeq[M10] = 0.;
        CMeq[M01] = 0.;

        // Order 2
        CMeq[M20] = D::cs2;
        CMeq[M02] = D::cs2;
        CMeq[M11] = 0.;

        // Order 3
        CMeq[M21] = 0.;
        CMeq[M12] = 0.;

        // Order 4
        CMeq[M22] = D::cs2 * D::cs2;
    };

    // Equilibrium populations based on 9 moments can be computed using either RM, HM, CM, CHM or
    // Gauss-Hermite formalisms. Here we use hermite moments (CMs)
    static void CMcomputeEquilibrium(
        T rho, Array<T, D::d> const &u, Array<T, D::q> const &CMeq, Array<T, D::q> &eq)
    {
        Array<T, D::q> RMeq;
        T ux2 = u[0] * u[0];
        T uy2 = u[1] * u[1];

        RMeq[M20] = CMeq[M20] + ux2;
        RMeq[M02] = CMeq[M02] + uy2;

        RMeq[M11] = CMeq[M11] + u[0] * u[1];

        RMeq[M21] = CMeq[M21] + u[1] * CMeq[M20] + 2. * u[0] * CMeq[M11] + ux2 * u[1];
        RMeq[M12] = CMeq[M12] + u[0] * CMeq[M02] + 2. * u[1] * CMeq[M11] + u[0] * uy2;

        RMeq[M22] = CMeq[M22] + 2. * u[1] * CMeq[M21] + 2. * u[0] * CMeq[M12] + uy2 * CMeq[M20]
                    + ux2 * CMeq[M02] + 4. * u[0] * u[1] * CMeq[M11] + ux2 * uy2;

        eq[F00] = rho * (1. - RMeq[M20] - RMeq[M02] + RMeq[M22]);

        eq[FP0] = 0.5 * rho * (u[0] + RMeq[M20] - RMeq[M12] - RMeq[M22]);
        eq[FM0] = 0.5 * rho * (-u[0] + RMeq[M20] + RMeq[M12] - RMeq[M22]);

        eq[F0P] = 0.5 * rho * (u[1] + RMeq[M02] - RMeq[M21] - RMeq[M22]);
        eq[F0M] = 0.5 * rho * (-u[1] + RMeq[M02] + RMeq[M21] - RMeq[M22]);

        eq[FPP] = 0.25 * rho * (RMeq[M11] + RMeq[M21] + RMeq[M12] + RMeq[M22]);
        eq[FMP] = 0.25 * rho * (-RMeq[M11] + RMeq[M21] - RMeq[M12] + RMeq[M22]);
        eq[FPM] = 0.25 * rho * (-RMeq[M11] - RMeq[M21] + RMeq[M12] + RMeq[M22]);
        eq[FMM] = 0.25 * rho * (RMeq[M11] - RMeq[M21] - RMeq[M12] + RMeq[M22]);
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
        T omegaPlus = (omegaBulk + omega1) / 2.;   // Notation used by Fei
        T omegaMinus = (omegaBulk - omega1) / 2.;  // Notation used by Fei

        T ux2 = u[0] * u[0];
        T uy2 = u[1] * u[1];

        // Post-collision moments.
        Array<T, D::q> CMcoll;
        Array<T, D::q> RMcoll;

        // Collision in the central moment space
        // Order 2 (non-diagonal collision so that we can easily modify the bulk viscosity)
        CMcoll[M20] =
            CM[M20] - omegaPlus * (CM[M20] - CMeq[M20]) - omegaMinus * (CM[M02] - CMeq[M02]);
        CMcoll[M02] =
            CM[M02] - omegaMinus * (CM[M20] - CMeq[M20]) - omegaPlus * (CM[M02] - CMeq[M02]);

        CMcoll[M11] = (1. - omega2) * CM[M11] + omega2 * CMeq[M11];

        // Order 3
        CMcoll[M21] = (1. - omega3) * CM[M21] + omega3 * CMeq[M21];
        CMcoll[M12] = (1. - omega3) * CM[M12] + omega3 * CMeq[M12];

        // Order 4
        CMcoll[M22] = (1. - omega4) * CM[M22] + omega4 * CMeq[M22];

        // Come back to RMcoll using binomial formulas
        RMcoll[M20] = CMcoll[M20] + ux2;
        RMcoll[M02] = CMcoll[M02] + uy2;

        RMcoll[M11] = CMcoll[M11] + u[0] * u[1];

        RMcoll[M21] = CMcoll[M21] + u[1] * CMcoll[M20] + 2. * u[0] * CMcoll[M11] + ux2 * u[1];
        RMcoll[M12] = CMcoll[M12] + u[0] * CMcoll[M02] + 2. * u[1] * CMcoll[M11] + u[0] * uy2;

        RMcoll[M22] = CMcoll[M22] + 2. * u[1] * CMcoll[M21] + 2. * u[0] * CMcoll[M12]
                      + uy2 * CMcoll[M20] + ux2 * CMcoll[M02] + 4. * u[0] * u[1] * CMcoll[M11]
                      + ux2 * uy2;

        // Compute post collision populations from RM
        cell[F00] = rho * (1. - RMcoll[M20] - RMcoll[M02] + RMcoll[M22]);

        cell[FP0] = 0.5 * rho * (u[0] + RMcoll[M20] - RMcoll[M12] - RMcoll[M22]);
        cell[FM0] = 0.5 * rho * (-u[0] + RMcoll[M20] + RMcoll[M12] - RMcoll[M22]);

        cell[F0P] = 0.5 * rho * (u[1] + RMcoll[M02] - RMcoll[M21] - RMcoll[M22]);
        cell[F0M] = 0.5 * rho * (-u[1] + RMcoll[M02] + RMcoll[M21] - RMcoll[M22]);

        cell[FPP] = 0.25 * rho * (RMcoll[M11] + RMcoll[M21] + RMcoll[M12] + RMcoll[M22]);
        cell[FMP] = 0.25 * rho * (-RMcoll[M11] + RMcoll[M21] - RMcoll[M12] + RMcoll[M22]);
        cell[FPM] = 0.25 * rho * (-RMcoll[M11] - RMcoll[M21] + RMcoll[M12] + RMcoll[M22]);
        cell[FMM] = 0.25 * rho * (RMcoll[M11] - RMcoll[M21] - RMcoll[M12] + RMcoll[M22]);

        for (int i = 0; i < D::q; ++i) {
            cell[i] -= D::SkordosFactor() * D::t[i];
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Central Hermite Moments Formalism (Equilibrium is computed through raw moments formalism) //
    ///////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * // General way to compute CHMs (computations of velocity components are based on Palabos
     * ordering of discrete velocities !) static void CHMcomputeMoments(Array<T,D::q> const& cell,
     * Array<T, D::q>& CHM, T& rho, Array<T,D::d>& u) {
     *
     *     Array<T, D::q> f;
     *     for (int i = 0; i<D::q; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *         CHM[i] = 0.;
     *     }
     *
     *     rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
     *     CHM[M00] = 1.0;
     *     T invRho = 1./rho;
     *     u[0] = invRho * ( - f[1] - f[2] - f[3] + f[5] + f[6] + f[7]);
     *     u[1] = invRho * (   f[1] - f[3] - f[4] - f[5] + f[7] + f[8]);
     *
     *     T cMux = 0.;
     *     T cMuy = 0.;
     *     T Hxx = 0.;
     *     T Hyy = 0.;
     *
     *     for (int i = 0; i<D::q; ++i) {
     *
     *         cMux = D::c[i][0] - u[0];
     *         cMuy = D::c[i][1] - u[1];
     *
     *         Hxx = cMux * cMux - D::cs2;
     *         Hyy = cMuy * cMuy - D::cs2;
     *
     *         // // Order 0
     *         // CHM[M000] += f[i];
     *
     *         // Order 1
     *         CHM[M10] += cMux * f[i];
     *         CHM[M01] += cMuy * f[i];
     *
     *         // Order 2
     *         CHM[M20] += Hxx * f[i];
     *         CHM[M02] += Hyy * f[i];
     *         CHM[M11] += cMux * cMuy * f[i];
     *
     *         // Order 3
     *         CHM[M21] += Hxx * cMuy * f[i];
     *         CHM[M12] += cMux * Hyy * f[i];
     *
     *         // Order 4
     *         CHM[M22] += Hxx * Hyy * f[i];
     *     }
     *
     *     for (int i = 1; i<D::q; ++i) {
     *         CHM[i] *= invRho;
     *     }
     * };
     */

    // Optimized way to compute CHs based on Palabos ordering of discrete velocities
    static void CHMcomputeMoments(
        Array<T, D::q> const &cell, Array<T, D::q> &CHM, T &rho, Array<T, D::d> &u)
    {
        Array<T, D::q> f;
        for (int i = 0; i < D::q; ++i) {
            f[i] = cell[i] + D::SkordosFactor() * D::t[i];
            CHM[i] = 0.;
        }

        double X_M1 = f[1] + f[2] + f[3];
        double X_P1 = f[5] + f[6] + f[7];
        double X_0 = f[0] + f[4] + f[8];

        double Y_M1 = f[3] + f[4] + f[5];
        double Y_P1 = f[1] + f[7] + f[8];

        // Order 0
        CHM[M00] = X_M1 + X_P1 + X_0;

        // Order 1
        CHM[M10] = 0.;
        CHM[M01] = 0.;

        // Order 2
        CHM[M20] = X_P1 + X_M1;
        CHM[M02] = Y_P1 + Y_M1;
        CHM[M11] = -f[1] + f[3] - f[5] + f[7];

        // Order 3
        CHM[M21] = f[1] - f[3] - f[5] + f[7];
        CHM[M12] = -f[1] - f[3] + f[5] + f[7];

        // Order 4
        CHM[M22] = f[1] + f[3] + f[5] + f[7];

        rho = CHM[M00];
        T invRho = 1. / rho;
        for (int i = 0; i < D::q; ++i) {
            CHM[i] *= invRho;
        }

        // Compute CMs from RMs using binomial formulas
        u[0] = invRho * (X_P1 - X_M1);
        u[1] = invRho * (Y_P1 - Y_M1);
        T ux2 = u[0] * u[0];
        T uy2 = u[1] * u[1];

        CHM[M20] -= (ux2);
        CHM[M02] -= (uy2);
        CHM[M11] -= (u[0] * u[1]);

        CHM[M21] -= (u[1] * CHM[M20] + 2. * u[0] * CHM[M11] + ux2 * u[1]);
        CHM[M12] -= (u[0] * CHM[M02] + 2. * u[1] * CHM[M11] + u[0] * uy2);

        CHM[M22] -=
            (2. * u[1] * CHM[M21] + 2. * u[0] * CHM[M12] + uy2 * CHM[M20] + ux2 * CHM[M02]
             + 4. * u[0] * u[1] * CHM[M11] + ux2 * uy2);

        // Compute CHMs from CMs
        T cs4 = D::cs2 * D::cs2;

        CHM[M20] -= (D::cs2);
        CHM[M02] -= (D::cs2);

        CHM[M22] -= (D::cs2 * (CHM[M20] + CHM[M02]) + cs4);
    };

    static void CHMcomputeEquilibriumMoments(Array<T, D::q> &CHMeq)
    {
        // Order 0
        CHMeq[M00] = 1.;

        // Order 1
        CHMeq[M10] = 0.;
        CHMeq[M01] = 0.;

        // Order 2
        CHMeq[M20] = 0.;
        CHMeq[M02] = 0.;
        CHMeq[M11] = 0.;

        // Order 3
        CHMeq[M21] = 0.;
        CHMeq[M12] = 0.;

        // Order 4
        CHMeq[M22] = 0.;
    };

    // Equilibrium populations based on 9 moments can be computed using either RM, HM, CM, CHM or
    // Gauss-Hermite formalisms. Here we use central hermite moments (CHMs)
    static void CHMcomputeEquilibrium(
        T rho, Array<T, D::d> const &u, Array<T, D::q> const &CHMeq, Array<T, D::q> &eq)
    {
        Array<T, D::q> RMeq;
        Array<T, D::q> HMeq;
        T ux2 = u[0] * u[0];
        T uy2 = u[1] * u[1];
        T cs4 = D::cs2 * D::cs2;

        // Come back to HMeq using relationships between HMs and CHMs
        HMeq[M20] = CHMeq[M20] + ux2;
        HMeq[M02] = CHMeq[M02] + uy2;

        HMeq[M11] = CHMeq[M11] + u[0] * u[1];

        HMeq[M21] = CHMeq[M21] + u[1] * CHMeq[M20] + 2. * u[0] * CHMeq[M11] + ux2 * u[1];
        HMeq[M12] = CHMeq[M12] + u[0] * CHMeq[M02] + 2. * u[1] * CHMeq[M11] + u[0] * uy2;

        HMeq[M22] = CHMeq[M22] + 2. * u[1] * CHMeq[M21] + 2. * u[0] * CHMeq[M12] + uy2 * CHMeq[M20]
                    + ux2 * CHMeq[M02] + 4. * u[0] * u[1] * CHMeq[M11] + ux2 * uy2;

        // Come back to RMeq using relationships between HMs and RMs
        RMeq[M20] = HMeq[M20] + D::cs2;
        RMeq[M02] = HMeq[M02] + D::cs2;

        RMeq[M11] = HMeq[M11];

        RMeq[M21] = HMeq[M21] + D::cs2 * u[1];
        RMeq[M12] = HMeq[M12] + D::cs2 * u[0];

        RMeq[M22] = HMeq[M22] + D::cs2 * (HMeq[M20] + HMeq[M02]) + cs4;

        eq[F00] = rho * (1. - RMeq[M20] - RMeq[M02] + RMeq[M22]);

        eq[FP0] = 0.5 * rho * (u[0] + RMeq[M20] - RMeq[M12] - RMeq[M22]);
        eq[FM0] = 0.5 * rho * (-u[0] + RMeq[M20] + RMeq[M12] - RMeq[M22]);

        eq[F0P] = 0.5 * rho * (u[1] + RMeq[M02] - RMeq[M21] - RMeq[M22]);
        eq[F0M] = 0.5 * rho * (-u[1] + RMeq[M02] + RMeq[M21] - RMeq[M22]);

        eq[FPP] = 0.25 * rho * (RMeq[M11] + RMeq[M21] + RMeq[M12] + RMeq[M22]);
        eq[FMP] = 0.25 * rho * (-RMeq[M11] + RMeq[M21] - RMeq[M12] + RMeq[M22]);
        eq[FPM] = 0.25 * rho * (-RMeq[M11] - RMeq[M21] + RMeq[M12] + RMeq[M22]);
        eq[FMM] = 0.25 * rho * (RMeq[M11] - RMeq[M21] - RMeq[M12] + RMeq[M22]);
    };

    static void CHMcollide(
        Array<T, D::q> &cell, T rho, Array<T, D::d> const &u,
        Array<T, D::q> const &CHM,    // Central hermite moments
        Array<T, D::q> const &CHMeq,  // Equilibrium moments (central hermite)
        Array<T, D::numRelaxationTimes> const &omega)
    {
        T omega1 = omega[0];
        T omega2 = omega[1];
        T omega3 = omega[2];
        T omega4 = omega[3];

        T omegaBulk = omega[4];
        T omegaPlus = (omegaBulk + omega1) / 2.;   // Notation used by Fei
        T omegaMinus = (omegaBulk - omega1) / 2.;  // Notation used by Fei

        T ux2 = u[0] * u[0];
        T uy2 = u[1] * u[1];
        T cs4 = D::cs2 * D::cs2;

        // Post-collision moments.
        Array<T, D::q> CHMcoll;
        Array<T, D::q> HMcoll;
        Array<T, D::q> RMcoll;

        // Collision in the central Hermite moment space
        // Order 2 (non-diagonal collision so that we can easily modify the bulk viscosity)
        CHMcoll[M20] =
            CHM[M20] - omegaPlus * (CHM[M20] - CHMeq[M20]) - omegaMinus * (CHM[M02] - CHMeq[M02]);
        CHMcoll[M02] =
            CHM[M02] - omegaMinus * (CHM[M20] - CHMeq[M20]) - omegaPlus * (CHM[M02] - CHMeq[M02]);

        CHMcoll[M11] = (1. - omega2) * CHM[M11] + omega2 * CHMeq[M11];

        // Order 3
        CHMcoll[M21] = (1. - omega3) * CHM[M21] + omega3 * CHMeq[M21];
        CHMcoll[M12] = (1. - omega3) * CHM[M12] + omega3 * CHMeq[M12];

        // Order 4
        CHMcoll[M22] = (1. - omega4) * CHM[M22] + omega4 * CHMeq[M22];

        // Come back to HMcoll using relationships between CHMs and HMs
        HMcoll[M20] = CHMcoll[M20] + ux2;
        HMcoll[M02] = CHMcoll[M02] + uy2;

        HMcoll[M11] = CHMcoll[M11] + u[0] * u[1];

        HMcoll[M21] = CHMcoll[M21] + u[1] * CHMcoll[M20] + 2. * u[0] * CHMcoll[M11] + ux2 * u[1];
        HMcoll[M12] = CHMcoll[M12] + u[0] * CHMcoll[M02] + 2. * u[1] * CHMcoll[M11] + u[0] * uy2;

        HMcoll[M22] = CHMcoll[M22] + 2. * u[1] * CHMcoll[M21] + 2. * u[0] * CHMcoll[M12]
                      + uy2 * CHMcoll[M20] + ux2 * CHMcoll[M02] + 4. * u[0] * u[1] * CHMcoll[M11]
                      + ux2 * uy2;

        // Come back to RMcoll using relationships between HMs and RMs
        RMcoll[M20] = HMcoll[M20] + D::cs2;
        RMcoll[M02] = HMcoll[M02] + D::cs2;

        RMcoll[M11] = HMcoll[M11];

        RMcoll[M21] = HMcoll[M21] + D::cs2 * u[1];
        RMcoll[M12] = HMcoll[M12] + D::cs2 * u[0];

        RMcoll[M22] = HMcoll[M22] + D::cs2 * (HMcoll[M20] + HMcoll[M02]) + cs4;

        // Compute post collision populations from RM
        cell[F00] = rho * (1. - RMcoll[M20] - RMcoll[M02] + RMcoll[M22]);

        cell[FP0] = 0.5 * rho * (u[0] + RMcoll[M20] - RMcoll[M12] - RMcoll[M22]);
        cell[FM0] = 0.5 * rho * (-u[0] + RMcoll[M20] + RMcoll[M12] - RMcoll[M22]);

        cell[F0P] = 0.5 * rho * (u[1] + RMcoll[M02] - RMcoll[M21] - RMcoll[M22]);
        cell[F0M] = 0.5 * rho * (-u[1] + RMcoll[M02] + RMcoll[M21] - RMcoll[M22]);

        cell[FPP] = 0.25 * rho * (RMcoll[M11] + RMcoll[M21] + RMcoll[M12] + RMcoll[M22]);
        cell[FMP] = 0.25 * rho * (-RMcoll[M11] + RMcoll[M21] - RMcoll[M12] + RMcoll[M22]);
        cell[FPM] = 0.25 * rho * (-RMcoll[M11] - RMcoll[M21] + RMcoll[M12] + RMcoll[M22]);
        cell[FMM] = 0.25 * rho * (RMcoll[M11] - RMcoll[M21] - RMcoll[M12] + RMcoll[M22]);

        for (int i = 0; i < D::q; ++i) {
            cell[i] -= D::SkordosFactor() * D::t[i];
        }
    };

    ////////////////////////////////////////////////////////////////////////////////
    // Cumulant Formalism (Equilibrium is computed through raw moments formalism) //
    ////////////////////////////////////////////////////////////////////////////////

    /**
     * // General way to compute Ks (computations of velocity components are based on Palabos
     * ordering of discrete velocities !) static void KcomputeMoments(Array<T,D::q> const& cell,
     * Array<T, D::q>& K, T& rho, Array<T,D::d>& u) {
     *
     *     Array<T, D::q> f;
     *     Array<T,D::q> CM;
     *     for (int i = 0; i<D::q; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *         CM[i] = 0.;
     *     }
     *
     *     rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
     *     CM[M00] = 1.0;
     *     T invRho = 1./rho;
     *     u[0] = invRho * ( - f[1] - f[2] - f[3] + f[5] + f[6] + f[7]);
     *     u[1] = invRho * (   f[1] - f[3] - f[4] - f[5] + f[7] + f[8]);
     *
     *     T cMux = 0.;
     *     T cMuy = 0.;
     *
     *     // Computation of central moments in a first time
     *     for (int i = 0; i<D::q; ++i) {
     *
     *         cMux = D::c[i][0]- u[0];
     *         cMuy = D::c[i][1]- u[1];
     *
     *         // // Order 0
     *         // CM[M000] += f[i];
     *
     *         // Order 1
     *         CM[M10] += cMux * f[i];
     *         CM[M01] += cMuy * f[i];
     *
     *         // Order 2
     *         CM[M20] += cMux * cMux * f[i];
     *         CM[M02] += cMuy * cMuy * f[i];
     *         CM[M11] += cMux * cMuy * f[i];
     *
     *         // Order 3
     *         CM[M21] += cMux * cMux * cMuy * f[i];
     *         CM[M12] += cMux * cMuy * cMuy * f[i];
     *
     *         // Order 4
     *         CM[M22] += cMux * cMux * cMuy * cMuy * f[i];
     *     }
     *
     *     // Normalize before the computation of cumulants !
     *     for (int i = 1; i<D::q; ++i) {
     *         CM[i] *= invRho;
     *     }
     *
     *     // Computation of cumulants through central moments
     *     K[M00] = CM[M00];
     *     K[M10] = CM[M10] + u[0];
     *     K[M01] = CM[M01] + u[1];
     *     K[M20] = CM[M20];
     *     K[M02] = CM[M02];
     *     K[M11] = CM[M11];
     *     K[M21] = CM[M21];
     *     K[M12] = CM[M12];
     *     K[M22] = CM[M22] - CM[M20]*CM[M02] - 2.*CM[M11]*CM[M11];
     * };
     */

    // Optimized way to compute Kcs based on Palabos ordering of discrete velocities
    static void KcomputeMoments(
        Array<T, D::q> const &cell, Array<T, D::q> &K, T &rho, Array<T, D::d> &u)
    {
        Array<T, D::q> f;
        for (int i = 0; i < D::q; ++i) {
            f[i] = cell[i] + D::SkordosFactor() * D::t[i];
            K[i] = 0.;
        }

        double X_M1 = f[1] + f[2] + f[3];
        double X_P1 = f[5] + f[6] + f[7];
        double X_0 = f[0] + f[4] + f[8];

        double Y_M1 = f[3] + f[4] + f[5];
        double Y_P1 = f[1] + f[7] + f[8];

        // Order 0
        K[M00] = X_M1 + X_P1 + X_0;

        // Order 1
        K[M10] = X_P1 - X_M1;
        K[M01] = Y_P1 - Y_M1;

        // Order 2
        K[M20] = X_P1 + X_M1;
        K[M02] = Y_P1 + Y_M1;
        K[M11] = -f[1] + f[3] - f[5] + f[7];

        // Order 3
        K[M21] = f[1] - f[3] - f[5] + f[7];
        K[M12] = -f[1] - f[3] + f[5] + f[7];

        // Order 4
        K[M22] = f[1] + f[3] + f[5] + f[7];

        rho = K[M00];
        T invRho = 1. / rho;
        for (int i = 0; i < D::q; ++i) {
            K[i] *= invRho;
        }

        // Compute CMs from RMs using binomial formulas
        u[0] = K[M10];
        u[1] = K[M01];
        T ux2 = u[0] * u[0];
        T uy2 = u[1] * u[1];

        K[M20] -= (ux2);
        K[M02] -= (uy2);
        K[M11] -= (u[0] * u[1]);

        K[M21] -= (u[1] * K[M20] + 2. * u[0] * K[M11] + ux2 * u[1]);
        K[M12] -= (u[0] * K[M02] + 2. * u[1] * K[M11] + u[0] * uy2);

        K[M22] -=
            (2. * u[1] * K[M21] + 2. * u[0] * K[M12] + uy2 * K[M20] + ux2 * K[M02]
             + 4. * u[0] * u[1] * K[M11] + ux2 * uy2);

        // Computation of Ks through non-linear transformations of CMs
        K[M22] -= (K[M20] * K[M02] + 2. * K[M11] * K[M11]);
    };

    static void KcomputeEquilibriumMoments(Array<T, D::d> const &u, Array<T, D::q> &Keq)
    {
        // Order 0
        Keq[M00] = 1.;

        // Order 1
        Keq[M10] = u[0];
        Keq[M01] = u[1];

        // Order 2
        Keq[M20] = D::cs2;
        Keq[M02] = D::cs2;
        Keq[M11] = 0.;

        // Order 3
        Keq[M21] = 0.;
        Keq[M12] = 0.;

        // Order 4
        Keq[M22] = 0.;
    };

    // Equilibrium populations based on 9 moments can be computed using either RM, HM, CM, CHM or
    // Gauss-Hermite formalisms. Here we use cumulants (Ks)
    static void KcomputeEquilibrium(
        T rho, Array<T, D::d> const &u, Array<T, D::q> const &Keq, Array<T, D::q> &eq)
    {
        Array<T, D::q> CMeq;
        Array<T, D::q> RMeq;
        T ux2 = u[0] * u[0];
        T uy2 = u[1] * u[1];

        CMeq[M20] = Keq[M20];
        CMeq[M02] = Keq[M02];
        CMeq[M11] = Keq[M11];

        CMeq[M21] = Keq[M21];
        CMeq[M12] = Keq[M12];

        CMeq[M22] = Keq[M22] + Keq[M20] * Keq[M02] + 2. * Keq[M11] * Keq[M11];

        RMeq[M20] = CMeq[M20] + ux2;
        RMeq[M02] = CMeq[M02] + uy2;

        RMeq[M11] = CMeq[M11] + u[0] * u[1];

        RMeq[M21] = CMeq[M21] + u[1] * CMeq[M20] + 2. * u[0] * CMeq[M11] + ux2 * u[1];
        RMeq[M12] = CMeq[M12] + u[0] * CMeq[M02] + 2. * u[1] * CMeq[M11] + u[0] * uy2;

        RMeq[M22] = CMeq[M22] + 2. * u[1] * CMeq[M21] + 2. * u[0] * CMeq[M12] + uy2 * CMeq[M20]
                    + ux2 * CMeq[M02] + 4. * u[0] * u[1] * CMeq[M11] + ux2 * uy2;

        eq[F00] = rho * (1. - RMeq[M20] - RMeq[M02] + RMeq[M22]);

        eq[FP0] = 0.5 * rho * (u[0] + RMeq[M20] - RMeq[M12] - RMeq[M22]);
        eq[FM0] = 0.5 * rho * (-u[0] + RMeq[M20] + RMeq[M12] - RMeq[M22]);

        eq[F0P] = 0.5 * rho * (u[1] + RMeq[M02] - RMeq[M21] - RMeq[M22]);
        eq[F0M] = 0.5 * rho * (-u[1] + RMeq[M02] + RMeq[M21] - RMeq[M22]);

        eq[FPP] = 0.25 * rho * (RMeq[M11] + RMeq[M21] + RMeq[M12] + RMeq[M22]);
        eq[FMP] = 0.25 * rho * (-RMeq[M11] + RMeq[M21] - RMeq[M12] + RMeq[M22]);
        eq[FPM] = 0.25 * rho * (-RMeq[M11] - RMeq[M21] + RMeq[M12] + RMeq[M22]);
        eq[FMM] = 0.25 * rho * (RMeq[M11] - RMeq[M21] - RMeq[M12] + RMeq[M22]);
    };

    static void Kcollide(
        Array<T, D::q> &cell, T rho, Array<T, D::d> const &u,
        Array<T, D::q> const &K,    // Central moments
        Array<T, D::q> const &Keq,  // Equilibrium moments (central)
        Array<T, D::numRelaxationTimes> const &omega)
    {
        T omega1 = omega[0];
        T omega2 = omega[1];
        T omega3 = omega[2];
        T omega4 = omega[3];

        T omegaBulk = omega[4];
        T omegaPlus = (omegaBulk + omega1) / 2.;   // Notation used by Fei
        T omegaMinus = (omegaBulk - omega1) / 2.;  // Notation used by Fei

        T ux2 = u[0] * u[0];
        T uy2 = u[1] * u[1];

        // Post-collision moments.
        Array<T, D::q> Kcoll;
        Array<T, D::q> CMcoll;
        Array<T, D::q> RMcoll;

        // Collision in the cumulant space
        // Order 2 (non-diagonal collision so that we can easily modify the bulk viscosity)
        Kcoll[M20] = K[M20] - omegaPlus * (K[M20] - Keq[M20]) - omegaMinus * (K[M02] - Keq[M02]);
        Kcoll[M02] = K[M02] - omegaMinus * (K[M20] - Keq[M20]) - omegaPlus * (K[M02] - Keq[M02]);

        Kcoll[M11] = (1. - omega2) * K[M11] + omega2 * Keq[M11];

        // Order 3
        Kcoll[M21] = (1. - omega3) * K[M21] + omega3 * Keq[M21];
        Kcoll[M12] = (1. - omega3) * K[M12] + omega3 * Keq[M12];

        // Order 4
        Kcoll[M22] = (1. - omega4) * K[M22] + omega4 * Keq[M22];

        // Come back to CMcoll using modifying fourth- and higher-order post-collision cumulants
        CMcoll[M20] = Kcoll[M20];
        CMcoll[M02] = Kcoll[M02];
        CMcoll[M11] = Kcoll[M11];

        CMcoll[M21] = Kcoll[M21];
        CMcoll[M12] = Kcoll[M12];

        CMcoll[M22] = Kcoll[M22] + Kcoll[M20] * Kcoll[M02] + 2. * Kcoll[M11] * Kcoll[M11];

        // Come back to RMcoll using binomial formulas
        RMcoll[M20] = CMcoll[M20] + ux2;
        RMcoll[M02] = CMcoll[M02] + uy2;

        RMcoll[M11] = CMcoll[M11] + u[0] * u[1];

        RMcoll[M21] = CMcoll[M21] + u[1] * CMcoll[M20] + 2. * u[0] * CMcoll[M11] + ux2 * u[1];
        RMcoll[M12] = CMcoll[M12] + u[0] * CMcoll[M02] + 2. * u[1] * CMcoll[M11] + u[0] * uy2;

        RMcoll[M22] = CMcoll[M22] + 2. * u[1] * CMcoll[M21] + 2. * u[0] * CMcoll[M12]
                      + uy2 * CMcoll[M20] + ux2 * CMcoll[M02] + 4. * u[0] * u[1] * CMcoll[M11]
                      + ux2 * uy2;

        // Compute post collision populations from RM
        cell[F00] = rho * (1. - RMcoll[M20] - RMcoll[M02] + RMcoll[M22]);

        cell[FP0] = 0.5 * rho * (u[0] + RMcoll[M20] - RMcoll[M12] - RMcoll[M22]);
        cell[FM0] = 0.5 * rho * (-u[0] + RMcoll[M20] + RMcoll[M12] - RMcoll[M22]);

        cell[F0P] = 0.5 * rho * (u[1] + RMcoll[M02] - RMcoll[M21] - RMcoll[M22]);
        cell[F0M] = 0.5 * rho * (-u[1] + RMcoll[M02] + RMcoll[M21] - RMcoll[M22]);

        cell[FPP] = 0.25 * rho * (RMcoll[M11] + RMcoll[M21] + RMcoll[M12] + RMcoll[M22]);
        cell[FMP] = 0.25 * rho * (-RMcoll[M11] + RMcoll[M21] - RMcoll[M12] + RMcoll[M22]);
        cell[FPM] = 0.25 * rho * (-RMcoll[M11] - RMcoll[M21] + RMcoll[M12] + RMcoll[M22]);
        cell[FMM] = 0.25 * rho * (RMcoll[M11] - RMcoll[M21] - RMcoll[M12] + RMcoll[M22]);

        for (int i = 0; i < D::q; ++i) {
            cell[i] -= D::SkordosFactor() * D::t[i];
        }
    };

    ////////////////////////////////////////////////////////////////////////////////
    // Gauss-Hermite Formalism (Equilibrium is computed through the GH formalism) //
    ////////////////////////////////////////////////////////////////////////////////

    /**
     * // General way to compute GHs
     * static void GHcomputeMoments(Array<T,D::q> const& cell, Array<T, D::q>& GH, T& rho) {
     *
     *     Array<T, D::q> f;
     *     for (int i = 0; i<D::q; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *         GH[i] = 0.;
     *     }
     *
     *     T Hxx = 0.;
     *     T Hyy = 0.;
     *
     *     for (int i = 0; i<D::q; ++i) {
     *
     *         Hxx = D::c[i][0] * D::c[i][0] - D::cs2;
     *         Hyy = D::c[i][1] * D::c[i][1] - D::cs2;
     *
     *         // Order 0
     *         GH[M00] += f[i];
     *
     *         // Order 1
     *         GH[M10] += D::c[i][0] * f[i];
     *         GH[M01] += D::c[i][1] * f[i];
     *
     *         // Order 2
     *         GH[M20] += Hxx * f[i];
     *         GH[M02] += Hyy * f[i];
     *         GH[M11] += D::c[i][0] * D::c[i][1] * f[i];
     *
     *         // Order 3
     *         GH[M21] += Hxx * D::c[i][1] * f[i];
     *         GH[M12] += D::c[i][0] * Hyy * f[i];
     *
     *         // Order 4
     *         GH[M22] += Hxx * Hyy * f[i];
     *     }
     *
     *     rho = GH[M00];
     *     T invRho = 1. / rho;
     *     for (int i = 0; i<D::q; ++i) {
     *         GH[i] *= invRho;
     *     }
     * };
     */

    // Optimized way to compute GHs based on Palabos ordering of discrete velocities
    static void GHcomputeMoments(Array<T, D::q> const &cell, Array<T, D::q> &GH, T &rho)
    {
        Array<T, D::q> f;
        for (int i = 0; i < D::q; ++i) {
            f[i] = cell[i] + D::SkordosFactor() * D::t[i];
            GH[i] = 0.;
        }

        double X_M1 = f[1] + f[2] + f[3];
        double X_P1 = f[5] + f[6] + f[7];
        double X_0 = f[0] + f[4] + f[8];

        double Y_M1 = f[3] + f[4] + f[5];
        double Y_P1 = f[1] + f[7] + f[8];

        // Order 0
        GH[M00] = X_M1 + X_P1 + X_0;

        // Order 1
        GH[M10] = X_P1 - X_M1;
        GH[M01] = Y_P1 - Y_M1;

        // Order 2
        GH[M20] = X_P1 + X_M1;
        GH[M02] = Y_P1 + Y_M1;
        GH[M11] = -f[1] + f[3] - f[5] + f[7];

        // Order 3
        GH[M21] = f[1] - f[3] - f[5] + f[7];
        GH[M12] = -f[1] - f[3] + f[5] + f[7];

        // Order 4
        GH[M22] = f[1] + f[3] + f[5] + f[7];

        rho = GH[M00];
        T invRho = 1. / rho;
        for (int i = 0; i < D::q; ++i) {
            GH[i] *= invRho;
        }

        // We come back to Hermite moments
        T cs4 = D::cs2 * D::cs2;
        GH[M20] -= D::cs2;
        GH[M02] -= D::cs2;

        GH[M21] -= D::cs2 * GH[M01];
        GH[M12] -= D::cs2 * GH[M10];

        GH[M22] -= (D::cs2 * (GH[M20] + GH[M02]) + cs4);
    };

    static void GHcomputeEquilibriumMoments(Array<T, D::d> const &u, Array<T, D::q> &GHeq)
    {
        // Order 0
        GHeq[M00] = 1.;

        // Order 1
        GHeq[M10] = u[0];
        GHeq[M01] = u[1];

        // Order 2
        GHeq[M20] = u[0] * u[0];
        GHeq[M02] = u[1] * u[1];
        GHeq[M11] = u[0] * u[1];

        // Order 3
        GHeq[M21] = GHeq[M20] * u[1];
        GHeq[M12] = GHeq[M02] * u[0];

        // Order 4
        GHeq[M22] = GHeq[M20] * GHeq[M02];
    };

    // Equilibrium populations based on 9 moments can be computed using either RM, HM, CM, CHM or
    // Gauss-Hermite formalisms. Here we use hermite moments (GHs)
    static void GHcomputeEquilibrium(T rho, Array<T, D::q> const &GHeq, Array<T, D::q> &eq)
    {
        Array<T, D::d> u(GHeq[1], GHeq[2]);

        eq[F00] =
            (rho * 4. / 9.) * (1. - (3. / 2.) * (GHeq[M20] + GHeq[M02]) + (9. / 4.) * (GHeq[M22]));

        eq[FP0] = (rho / 9.)
                  * (1. + 3. * u[0] + (3. / 2.) * (2. * GHeq[M20] - GHeq[M02])
                     - (9. / 2.) * (GHeq[M12]) - (9. / 4.) * (2. * GHeq[M22]));
        eq[FM0] = (rho / 9.)
                  * (1. - 3. * u[0] + (3. / 2.) * (2. * GHeq[M20] - GHeq[M02])
                     + (9. / 2.) * (GHeq[M12]) - (9. / 4.) * (2. * GHeq[M22]));
        eq[F0P] = (rho / 9.)
                  * (1. + 3. * u[1] + (3. / 2.) * (2. * GHeq[M02] - GHeq[M20])
                     - (9. / 2.) * (GHeq[M21]) - (9. / 4.) * (2. * GHeq[M22]));
        eq[F0M] = (rho / 9.)
                  * (1. - 3. * u[1] + (3. / 2.) * (2. * GHeq[M02] - GHeq[M20])
                     + (9. / 2.) * (GHeq[M21]) - (9. / 4.) * (2. * GHeq[M22]));

        eq[FPP] = (rho / 36.)
                  * (1. + 3. * (+u[0] + u[1]) + (3. / 2.) * (2. * GHeq[M20] + 2. * GHeq[M02])
                     + 9. * GHeq[M11] + (9. / 2.) * (2. * GHeq[M21] + 2. * GHeq[M12])
                     + (9. / 2.) * (2. * GHeq[M22]));
        eq[FMP] = (rho / 36.)
                  * (1. + 3. * (-u[0] + u[1]) + (3. / 2.) * (2. * GHeq[M20] + 2. * GHeq[M02])
                     - 9. * GHeq[M11] + (9. / 2.) * (2. * GHeq[M21] - 2. * GHeq[M12])
                     + (9. / 2.) * (2. * GHeq[M22]));
        eq[FPM] = (rho / 36.)
                  * (1. + 3. * (+u[0] - u[1]) + (3. / 2.) * (2. * GHeq[M20] + 2. * GHeq[M02])
                     - 9. * GHeq[M11] + (9. / 2.) * (-2. * GHeq[M21] + 2. * GHeq[M12])
                     + (9. / 2.) * (2. * GHeq[M22]));
        eq[FMM] = (rho / 36.)
                  * (1. + 3. * (-u[0] - u[1]) + (3. / 2.) * (2. * GHeq[M20] + 2. * GHeq[M02])
                     + 9. * GHeq[M11] + (9. / 2.) * (-2. * GHeq[M21] - 2. * GHeq[M12])
                     + (9. / 2.) * (2. * GHeq[M22]));
    };

    static void GHcollide(
        Array<T, D::q> &cell, T rho, Array<T, D::d> const &u,
        Array<T, D::q> const &GH,    // Hermite moments
        Array<T, D::q> const &GHeq,  // Equilibrium moments (hermite)
        Array<T, D::numRelaxationTimes> const &omega)
    {
        T omega1 = omega[0];
        T omega2 = omega[1];
        T omega3 = omega[2];
        T omega4 = omega[3];

        T omegaBulk = omega[4];
        T omegaPlus = (omegaBulk + omega1) / 2.;   // Notation used by Fei
        T omegaMinus = (omegaBulk - omega1) / 2.;  // Notation used by Fei

        // Post-collision moments.
        Array<T, D::q> GHcoll;

        // Collision in the Hermite moment space
        // Order 2 (non-diagonal collision so that we can easily modify the bulk viscosity)
        GHcoll[M20] =
            GH[M20] - omegaPlus * (GH[M20] - GHeq[M20]) - omegaMinus * (GH[M02] - GHeq[M02]);
        GHcoll[M02] =
            GH[M02] - omegaMinus * (GH[M20] - GHeq[M20]) - omegaPlus * (GH[M02] - GHeq[M02]);

        GHcoll[M11] = (1. - omega2) * GH[M11] + omega2 * GHeq[M11];

        // Order 3
        GHcoll[M21] = (1. - omega3) * GH[M21] + omega3 * GHeq[M21];
        GHcoll[M12] = (1. - omega3) * GH[M12] + omega3 * GHeq[M12];

        // Order 4
        GHcoll[M22] = (1. - omega4) * GH[M22] + omega4 * GHeq[M22];

        // Compute post collision populations from RM
        cell[F00] = (rho * 4. / 9.)
                    * (1. - (3. / 2.) * (GHcoll[M20] + GHcoll[M02]) + (9. / 4.) * (GHcoll[M22]));

        cell[FP0] = (rho / 9.)
                    * (1. + 3. * u[0] + (3. / 2.) * (2. * GHcoll[M20] - GHcoll[M02])
                       - (9. / 2.) * (GHcoll[M12]) - (9. / 4.) * (2. * GHcoll[M22]));
        cell[FM0] = (rho / 9.)
                    * (1. - 3. * u[0] + (3. / 2.) * (2. * GHcoll[M20] - GHcoll[M02])
                       + (9. / 2.) * (GHcoll[M12]) - (9. / 4.) * (2. * GHcoll[M22]));
        cell[F0P] = (rho / 9.)
                    * (1. + 3. * u[1] + (3. / 2.) * (2. * GHcoll[M02] - GHcoll[M20])
                       - (9. / 2.) * (GHcoll[M21]) - (9. / 4.) * (2. * GHcoll[M22]));
        cell[F0M] = (rho / 9.)
                    * (1. - 3. * u[1] + (3. / 2.) * (2. * GHcoll[M02] - GHcoll[M20])
                       + (9. / 2.) * (GHcoll[M21]) - (9. / 4.) * (2. * GHcoll[M22]));

        cell[FPP] = (rho / 36.)
                    * (1. + 3. * (+u[0] + u[1]) + (3. / 2.) * (2. * GHcoll[M20] + 2. * GHcoll[M02])
                       + 9. * GHcoll[M11] + (9. / 2.) * (2. * GHcoll[M21] + 2. * GHcoll[M12])
                       + (9. / 2.) * (2. * GHcoll[M22]));
        cell[FMP] = (rho / 36.)
                    * (1. + 3. * (-u[0] + u[1]) + (3. / 2.) * (2. * GHcoll[M20] + 2. * GHcoll[M02])
                       - 9. * GHcoll[M11] + (9. / 2.) * (2. * GHcoll[M21] - 2. * GHcoll[M12])
                       + (9. / 2.) * (2. * GHcoll[M22]));
        cell[FPM] = (rho / 36.)
                    * (1. + 3. * (+u[0] - u[1]) + (3. / 2.) * (2. * GHcoll[M20] + 2. * GHcoll[M02])
                       - 9. * GHcoll[M11] + (9. / 2.) * (-2. * GHcoll[M21] + 2. * GHcoll[M12])
                       + (9. / 2.) * (2. * GHcoll[M22]));
        cell[FMM] = (rho / 36.)
                    * (1. + 3. * (-u[0] - u[1]) + (3. / 2.) * (2. * GHcoll[M20] + 2. * GHcoll[M02])
                       + 9. * GHcoll[M11] + (9. / 2.) * (-2. * GHcoll[M21] - 2. * GHcoll[M12])
                       + (9. / 2.) * (2. * GHcoll[M22]));

        for (int i = 0; i < D::q; ++i) {
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
     *     for (int i = 0; i<D::q; ++i) {
     *         f[i] = cell[i] + D::SkordosFactor() * D::t[i];
     *         RR[i] = 0.;
     *     }
     *
     *     T Hxx = 0.;
     *     T Hyy = 0.;
     *
     *     for (int i = 0; i<D::q; ++i) {
     *
     *         Hxx = D::c[i][0] * D::c[i][0] - D::cs2;
     *         Hyy = D::c[i][1] * D::c[i][1] - D::cs2;
     *
     *         // Order 0
     *         RR[M00] += f[i];
     *
     *         // Order 1
     *         RR[M10] += D::c[i][0] * f[i];
     *         RR[M01] += D::c[i][1] * f[i];
     *
     *         // Order 2
     *         RR[M20] += Hxx * f[i];
     *         RR[M02] += Hyy * f[i];
     *         RR[M11] += D::c[i][0] * D::c[i][1] * f[i];
     *     }
     *
     *     rho = RR[M00];
     *     T invRho = 1. / rho;
     *     for (int i = 0; i<D::q; ++i) {
     *         RR[i] *= invRho;
     *     }
     * };
     */

    // Optimized way to compute RRs based on Palabos ordering of discrete velocities
    // (we only need up to second order RRs since high-order ones are computed recursively)
    static void RRcomputeMoments(Array<T, D::q> const &cell, Array<T, D::q> &RR, T &rho)
    {
        Array<T, D::q> f;
        for (int i = 0; i < D::q; ++i) {
            f[i] = cell[i] + D::SkordosFactor() * D::t[i];
            RR[i] = 0.;
        }

        double X_M1 = f[1] + f[2] + f[3];
        double X_P1 = f[5] + f[6] + f[7];
        double X_0 = f[0] + f[4] + f[8];

        double Y_M1 = f[3] + f[4] + f[5];
        double Y_P1 = f[1] + f[7] + f[8];

        // Order 0
        RR[M00] = X_M1 + X_P1 + X_0;

        // Order 1
        RR[M10] = X_P1 - X_M1;
        RR[M01] = Y_P1 - Y_M1;

        // Order 2
        RR[M20] = X_P1 + X_M1;
        RR[M02] = Y_P1 + Y_M1;
        RR[M11] = -f[1] + f[3] - f[5] + f[7];

        rho = RR[M00];
        T invRho = 1. / rho;
        for (int i = 0; i < D::q; ++i) {
            RR[i] *= invRho;
        }

        // We come back to Hermite moments
        RR[M20] -= D::cs2;
        RR[M02] -= D::cs2;
    };

    static void RRcomputeEquilibriumMoments(Array<T, D::d> const &u, Array<T, D::q> &RReq)
    {
        // Order 0
        RReq[M00] = 1.;

        // Order 1
        RReq[M10] = u[0];
        RReq[M01] = u[1];

        // Order 2
        RReq[M20] = u[0] * u[0];
        RReq[M02] = u[1] * u[1];
        RReq[M11] = u[0] * u[1];

        // Order 3
        RReq[M21] = RReq[M20] * u[1];
        RReq[M12] = RReq[M02] * u[0];

        // Order 4
        RReq[M22] = RReq[M20] * RReq[M02];
    };

    // Equilibrium populations based on 9 moments can be computed using either RM, HM, CM, CHM or
    // Gauss-Hermite formalisms. Here we use hermite moments (RRs)
    static void RRcomputeEquilibrium(T rho, Array<T, D::q> const &RReq, Array<T, D::q> &eq)
    {
        Array<T, D::d> u(RReq[1], RReq[2]);

        eq[F00] =
            (rho * 4. / 9.) * (1. - (3. / 2.) * (RReq[M20] + RReq[M02]) + (9. / 4.) * (RReq[M22]));

        eq[FP0] = (rho / 9.)
                  * (1. + 3. * u[0] + (3. / 2.) * (2. * RReq[M20] - RReq[M02])
                     - (9. / 2.) * (RReq[M12]) - (9. / 4.) * (2. * RReq[M22]));
        eq[FM0] = (rho / 9.)
                  * (1. - 3. * u[0] + (3. / 2.) * (2. * RReq[M20] - RReq[M02])
                     + (9. / 2.) * (RReq[M12]) - (9. / 4.) * (2. * RReq[M22]));
        eq[F0P] = (rho / 9.)
                  * (1. + 3. * u[1] + (3. / 2.) * (2. * RReq[M02] - RReq[M20])
                     - (9. / 2.) * (RReq[M21]) - (9. / 4.) * (2. * RReq[M22]));
        eq[F0M] = (rho / 9.)
                  * (1. - 3. * u[1] + (3. / 2.) * (2. * RReq[M02] - RReq[M20])
                     + (9. / 2.) * (RReq[M21]) - (9. / 4.) * (2. * RReq[M22]));

        eq[FPP] = (rho / 36.)
                  * (1. + 3. * (+u[0] + u[1]) + (3. / 2.) * (2. * RReq[M20] + 2. * RReq[M02])
                     + 9. * RReq[M11] + (9. / 2.) * (2. * RReq[M21] + 2. * RReq[M12])
                     + (9. / 2.) * (2. * RReq[M22]));
        eq[FMP] = (rho / 36.)
                  * (1. + 3. * (-u[0] + u[1]) + (3. / 2.) * (2. * RReq[M20] + 2. * RReq[M02])
                     - 9. * RReq[M11] + (9. / 2.) * (2. * RReq[M21] - 2. * RReq[M12])
                     + (9. / 2.) * (2. * RReq[M22]));
        eq[FPM] = (rho / 36.)
                  * (1. + 3. * (+u[0] - u[1]) + (3. / 2.) * (2. * RReq[M20] + 2. * RReq[M02])
                     - 9. * RReq[M11] + (9. / 2.) * (-2. * RReq[M21] + 2. * RReq[M12])
                     + (9. / 2.) * (2. * RReq[M22]));
        eq[FMM] = (rho / 36.)
                  * (1. + 3. * (-u[0] - u[1]) + (3. / 2.) * (2. * RReq[M20] + 2. * RReq[M02])
                     + 9. * RReq[M11] + (9. / 2.) * (-2. * RReq[M21] - 2. * RReq[M12])
                     + (9. / 2.) * (2. * RReq[M22]));
    };

    // Compute post collision populations from RM instead of RR
    // (this is 5-10% more efficient than with Gauss-Hermite formulas)
    static void RRcollide(
        Array<T, D::q> &cell, T rho, Array<T, D::d> const &u,
        Array<T, D::q> const &RR,    // Hermite moments
        Array<T, D::q> const &RReq,  // Equilibrium moments (hermite)
        Array<T, D::numRelaxationTimes> const &omega)
    {
        T omega1 = omega[0];
        T omega2 = omega[1];
        T omega3 = omega[2];
        T omega4 = omega[3];

        T omegaBulk = omega[4];
        T omegaPlus = (omegaBulk + omega1) / 2.;   // Notation used by Fei
        T omegaMinus = (omegaBulk - omega1) / 2.;  // Notation used by Fei

        // Post-collision and Nonequilibrium moments.
        Array<T, D::q> RRneq;
        Array<T, D::q> RRcoll;

        // Recursive computation of nonequilibrium Hermite moments
        // Order 2 (standard way to compute them)
        RRneq[M20] = RR[M20] - RReq[M20];
        RRneq[M02] = RR[M02] - RReq[M02];
        RRneq[M11] = RR[M11] - RReq[M11];
        // Order 3 (reconstruction using Chapman-Enskog formulas)
        RRneq[M21] = u[1] * RRneq[M20] + 2. * u[0] * RRneq[M11];
        RRneq[M12] = u[0] * RRneq[M02] + 2. * u[1] * RRneq[M11];
        // Order 4 (reconstruction using Chapman-Enskog formulas)
        RRneq[M22] =
            u[1] * u[1] * RRneq[M20] + u[0] * u[0] * RRneq[M02] + 4. * u[0] * u[1] * RRneq[M11];

        // Collision in the Hermite moment space
        // Order 2 (non-diagonal collision so that we can easily modify the bulk viscosity)
        RRcoll[M20] = RR[M20] - omegaPlus * RRneq[M20] - omegaMinus * RRneq[M02];
        RRcoll[M02] = RR[M02] - omegaMinus * RRneq[M20] - omegaPlus * RRneq[M02];

        RRcoll[M11] = (1. - omega2) * RRneq[M11] + RReq[M11];

        // Order 3
        RRcoll[M21] = (1. - omega3) * RRneq[M21] + RReq[M21];
        RRcoll[M12] = (1. - omega3) * RRneq[M12] + RReq[M12];

        // Order 4
        RRcoll[M22] = (1. - omega4) * RRneq[M22] + RReq[M22];

        // Come back to RMcoll using relationships between RRs and RMs
        T cs4 = D::cs2 * D::cs2;
        Array<T, D::q> RMcoll;

        RMcoll[M20] = RRcoll[M20] + D::cs2;
        RMcoll[M02] = RRcoll[M02] + D::cs2;

        RMcoll[M11] = RRcoll[M11];

        RMcoll[M21] = RRcoll[M21] + D::cs2 * u[1];
        RMcoll[M12] = RRcoll[M12] + D::cs2 * u[0];

        RMcoll[M22] = RRcoll[M22] + D::cs2 * (RRcoll[M20] + RRcoll[M02]) + cs4;

        // Compute post collision populations from RM
        cell[F00] = rho * (1. - RMcoll[M20] - RMcoll[M02] + RMcoll[M22]);

        cell[FP0] = 0.5 * rho * (u[0] + RMcoll[M20] - RMcoll[M12] - RMcoll[M22]);
        cell[FM0] = 0.5 * rho * (-u[0] + RMcoll[M20] + RMcoll[M12] - RMcoll[M22]);

        cell[F0P] = 0.5 * rho * (u[1] + RMcoll[M02] - RMcoll[M21] - RMcoll[M22]);
        cell[F0M] = 0.5 * rho * (-u[1] + RMcoll[M02] + RMcoll[M21] - RMcoll[M22]);

        cell[FPP] = 0.25 * rho * (RMcoll[M11] + RMcoll[M21] + RMcoll[M12] + RMcoll[M22]);
        cell[FMP] = 0.25 * rho * (-RMcoll[M11] + RMcoll[M21] - RMcoll[M12] + RMcoll[M22]);
        cell[FPM] = 0.25 * rho * (-RMcoll[M11] - RMcoll[M21] + RMcoll[M12] + RMcoll[M22]);
        cell[FMM] = 0.25 * rho * (RMcoll[M11] - RMcoll[M21] - RMcoll[M12] + RMcoll[M22]);
        for (int i = 0; i < D::q; ++i) {
            cell[i] -= D::SkordosFactor() * D::t[i];
        }
    };

};  // struct comprehensiveDynamicsTemplatesImpl<T, descriptors::D2Q9DescriptorBase<T> >

}  // namespace plb

#endif  // COMPREHENSIVE_MODELS_TEMPLATES_2D_H
