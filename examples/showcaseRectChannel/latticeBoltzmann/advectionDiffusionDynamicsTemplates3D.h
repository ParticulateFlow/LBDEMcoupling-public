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
 * Helper functions for the implementation of LB dynamics. This file is all
 * about efficiency. The generic template code is specialized for commonly
 * used Lattices, so that a maximum performance can be taken out of each
 * case.
 */
#ifndef ADVECTION_DIFFUSION_DYNAMICS_TEMPLATES_3D_H
#define ADVECTION_DIFFUSION_DYNAMICS_TEMPLATES_3D_H

namespace plb {

/// All helper functions are inside this structure
template <typename T>
struct advectionDiffusionDynamicsTemplatesImpl<T, descriptors::D3Q7DescriptorBase<T> > {
    typedef descriptors::D3Q7DescriptorBase<T> D;

    static T bgk_ma1_equilibrium(plint iPop, T rhoBar, Array<T, D::d> const &jEq)
    {
        return D::t[iPop]
               * (rhoBar
                  + D::invCs2
                        * (D::c[iPop][0] * jEq[0] + D::c[iPop][1] * jEq[1]
                           + D::c[iPop][2] * jEq[2]));
    }

    /// Regularization
    static void regularize(
        Array<T, D::q> &f, T rhoBar, Array<T, D::d> const &jAdvDiff, Array<T, D::d> const &jEq)
    {
        f[0] = D::t[0] * rhoBar;

        f[1] = D::t[1] * (rhoBar - D::invCs2 * jAdvDiff[0]);
        f[2] = D::t[2] * (rhoBar - D::invCs2 * jAdvDiff[1]);
        f[3] = D::t[3] * (rhoBar - D::invCs2 * jAdvDiff[2]);
        f[4] = D::t[4] * (rhoBar + D::invCs2 * jAdvDiff[0]);
        f[5] = D::t[5] * (rhoBar + D::invCs2 * jAdvDiff[1]);
        f[6] = D::t[6] * (rhoBar + D::invCs2 * jAdvDiff[2]);
    }

    static T no_corr_bgk_collision(Array<T, D::q> &f, T rhoBar, Array<T, D::d> const &jEq, T omega)
    {
        T invRho = D::invRho(rhoBar);
        const T jSqr = jEq[0] * jEq[0] + jEq[1] * jEq[1] + jEq[2] * jEq[2];

        const T oneMinusOmega = (T)1 - omega;
        const T halfOmega = (T)0.5 * omega;
        const T cs2RhoBar = D::cs2 * rhoBar;

        f[0] = oneMinusOmega * f[0] + omega * ((T)1 - (T)3 * D::cs2) * rhoBar;

        f[1] = oneMinusOmega * f[1] + halfOmega * (cs2RhoBar - jEq[0]);
        f[2] = oneMinusOmega * f[2] + halfOmega * (cs2RhoBar - jEq[1]);
        f[3] = oneMinusOmega * f[3] + halfOmega * (cs2RhoBar - jEq[2]);

        f[4] = oneMinusOmega * f[4] + halfOmega * (cs2RhoBar + jEq[0]);
        f[5] = oneMinusOmega * f[5] + halfOmega * (cs2RhoBar + jEq[1]);
        f[6] = oneMinusOmega * f[6] + halfOmega * (cs2RhoBar + jEq[2]);

        return jSqr * invRho * invRho;
    }

    static T no_corr_bgk_collision(
        Array<T, D::q> &f, T rhoBar, Array<T, D::d> const &jEq, T omega, T source)
    {
        T invRho = D::invRho(rhoBar);
        const T jSqr = jEq[0] * jEq[0] + jEq[1] * jEq[1] + jEq[2] * jEq[2];

        const T oneMinusOmega = (T)1 - omega;
        const T halfOmega = (T)0.5 * omega;
        const T cs2RhoBar = D::cs2 * rhoBar;
        const T halfSourceCs2 = (T)0.5 * source * D::cs2;

        f[0] = oneMinusOmega * f[0] + ((T)1 - (T)3 * D::cs2) * (omega * rhoBar + source);

        f[1] = oneMinusOmega * f[1] + halfOmega * (cs2RhoBar - jEq[0]) + halfSourceCs2;
        f[2] = oneMinusOmega * f[2] + halfOmega * (cs2RhoBar - jEq[1]) + halfSourceCs2;
        f[3] = oneMinusOmega * f[3] + halfOmega * (cs2RhoBar - jEq[2]) + halfSourceCs2;

        f[4] = oneMinusOmega * f[4] + halfOmega * (cs2RhoBar + jEq[0]) + halfSourceCs2;
        f[5] = oneMinusOmega * f[5] + halfOmega * (cs2RhoBar + jEq[1]) + halfSourceCs2;
        f[6] = oneMinusOmega * f[6] + halfOmega * (cs2RhoBar + jEq[2]) + halfSourceCs2;

        return jSqr * invRho * invRho;
    }

    static T no_corr_rlb_collision(
        Array<T, D::q> &f, T rhoBar, Array<T, D::d> const &jEq, Array<T, D::d> const &jNeq, T omega)
    {
        T invRho = D::invRho(rhoBar);
        const T jSqr = jEq[0] * jEq[0] + jEq[1] * jEq[1] + jEq[2] * jEq[2];

        const T oneHalfMinusHalfOmega = (T)0.5 - (T)0.5 * omega;
        const T cs2RhoBar = D::cs2 * rhoBar;

        const T jNeqTerm_0 = oneHalfMinusHalfOmega * jNeq[0];
        const T jNeqTerm_1 = oneHalfMinusHalfOmega * jNeq[1];
        const T jNeqTerm_2 = oneHalfMinusHalfOmega * jNeq[2];

        f[0] = ((T)1 - (T)3 * D::cs2) * rhoBar;

        f[1] = -jNeqTerm_0 + (T)0.5 * (cs2RhoBar - jEq[0]);
        f[2] = -jNeqTerm_1 + (T)0.5 * (cs2RhoBar - jEq[1]);
        f[3] = -jNeqTerm_2 + (T)0.5 * (cs2RhoBar - jEq[2]);

        f[4] = -f[1] + cs2RhoBar;
        f[5] = -f[2] + cs2RhoBar;
        f[6] = -f[3] + cs2RhoBar;

        return jSqr * invRho * invRho;
    }

    static T no_corr_rlb_collision(
        Array<T, D::q> &f, T rhoBar, Array<T, D::d> const &jEq, Array<T, D::d> const &jNeq, T omega,
        T source)
    {
        T invRho = D::invRho(rhoBar);
        const T jSqr = jEq[0] * jEq[0] + jEq[1] * jEq[1] + jEq[2] * jEq[2];

        const T oneHalfMinusHalfOmega = (T)0.5 - (T)0.5 * omega;
        const T cs2RhoBar = D::cs2 * rhoBar;
        const T halfSourceCs2 = (T)0.5 * source * D::cs2;

        const T jNeqTerm_0 = oneHalfMinusHalfOmega * jNeq[0];
        const T jNeqTerm_1 = oneHalfMinusHalfOmega * jNeq[1];
        const T jNeqTerm_2 = oneHalfMinusHalfOmega * jNeq[2];

        f[0] = ((T)1 - (T)3 * D::cs2) * (rhoBar + source);

        f[1] = -jNeqTerm_0 + (T)0.5 * (cs2RhoBar - jEq[0]);
        f[2] = -jNeqTerm_1 + (T)0.5 * (cs2RhoBar - jEq[1]);
        f[3] = -jNeqTerm_2 + (T)0.5 * (cs2RhoBar - jEq[2]);

        f[4] = -f[1] + cs2RhoBar;
        f[5] = -f[2] + cs2RhoBar;
        f[6] = -f[3] + cs2RhoBar;

        f[1] += halfSourceCs2;
        f[2] += halfSourceCs2;
        f[3] += halfSourceCs2;
        f[4] += halfSourceCs2;
        f[5] += halfSourceCs2;
        f[6] += halfSourceCs2;

        return jSqr * invRho * invRho;
    }

    static void complete_bgk_ma2_regularized_collision(
        Array<T, D::q> &f, T rhoPhiBar, T rhoBar, Array<T, D::d> const &jEq,
        Array<T, D::d> const &jNeq, const Array<T, SymmetricTensorImpl<T, D::d>::n> &piNeq, T omega,
        T omegaNonPhys, T omegaFluid, T omegaFluidNonPhys)
    {
        PLB_ASSERT(false && "method not implemented for d3q7");
    }

};  // struct advectionDiffusionDynamicsTemplatesImpl

template <typename T>
struct advectionDiffusionDynamicsTemplatesImpl<T, descriptors::D3Q15DescriptorBase<T> > {
    typedef descriptors::D3Q15DescriptorBase<T> D;

    static T computePsiComplete(T omega)
    {
        PLB_ASSERT(false && "method not implemented for d3q15");
    }

    static T bgk_ma1_equilibrium(plint iPop, T rhoBar, Array<T, D::d> const &jEq)
    {
        return dynamicsTemplatesImpl<T, D>::bgk_ma2_equilibrium(
            iPop, rhoBar, D::invRho(rhoBar), jEq, VectorTemplateImpl<T, D::d>::normSqr(jEq));
    }

    static void regularize(
        Array<T, D::q> &f, T rhoBar, Array<T, D::d> const &jAdvDiff, Array<T, D::d> const &jEq)
    {
        T invRho = D::invRho(rhoBar);
        T jSqr = jEq[0] * jEq[0] + jEq[1] * jEq[1] + jEq[2] * jEq[2];
        dynamicsTemplatesImpl<T, D>::bgk_ma2_equilibria(rhoBar, invRho, jEq, jSqr, f);
        Array<T, D::d> jNeq = jAdvDiff - jEq;
        for (plint iPop = 0; iPop < D::q; ++iPop) {
            f[iPop] +=
                D::t[iPop] * D::invCs2
                * (D::c[iPop][0] * jNeq[0] + D::c[iPop][1] * jNeq[1] + D::c[iPop][2] * jNeq[2]);
        }
    }

    static T no_corr_bgk_collision(Array<T, D::q> &f, T rhoBar, Array<T, D::d> const &jEq, T omega)
    {
        T invRho = D::invRho(rhoBar);
        T jSqr = dynamicsTemplatesImpl<T, D>::bgk_ma2_collision(f, rhoBar, jEq, omega);

        return jSqr * invRho * invRho;
    }

    static T no_corr_bgk_collision(
        Array<T, D::q> &f, T rhoBar, Array<T, D::d> const &jEq, T omega, T source)
    {
        PLB_ASSERT(
            false && "no correction bgk collision with source not implemented yet for d3q15.");
    }

    static T no_corr_rlb_collision(
        Array<T, D::q> &f, T rhoBar, Array<T, D::d> const &jEq, Array<T, D::d> const &jNeq, T omega)
    {
        T invRho = D::invRho(rhoBar);
        T jSqr = jEq[0] * jEq[0] + jEq[1] * jEq[1] + jEq[2] * jEq[2];
        dynamicsTemplatesImpl<T, D>::bgk_ma2_equilibria(rhoBar, invRho, jEq, jSqr, f);
        for (plint iPop = 0; iPop < D::q; ++iPop) {
            f[iPop] +=
                ((T)1 - omega)
                * (D::t[iPop] * D::invCs2
                   * (D::c[iPop][0] * jNeq[0] + D::c[iPop][1] * jNeq[1] + D::c[iPop][2] * jNeq[2]));
        }

        return jSqr * invRho * invRho;
    }

    static T no_corr_rlb_collision(
        Array<T, D::q> &f, T rhoBar, Array<T, D::d> const &jEq, Array<T, D::d> const &jNeq, T omega,
        T source)
    {
        T invRho = D::invRho(rhoBar);
        T jSqr = jEq[0] * jEq[0] + jEq[1] * jEq[1] + jEq[2] * jEq[2];
        dynamicsTemplatesImpl<T, D>::bgk_ma2_equilibria(rhoBar, invRho, jEq, jSqr, f);
        for (plint iPop = 0; iPop < D::q; ++iPop) {
            f[iPop] += ((T)1 - omega)
                           * (D::t[iPop] * D::invCs2
                              * (D::c[iPop][0] * jNeq[0] + D::c[iPop][1] * jNeq[1]
                                 + D::c[iPop][2] * jNeq[2]))
                       + D::t[iPop] * source;
        }

        return jSqr * invRho * invRho;
    }

    static void bgk_ma2_off_equilibra(
        T phi, Array<T, D::d> const &u, Array<T, D::d> const &jNeq,
        const Array<T, SymmetricTensorImpl<T, D::d>::n> &piNeq, T omega, T omegaFluid,
        Array<T, D::q> &fNeq)
    {
        typedef SymmetricTensorImpl<T, D::d> S;

        T omegaRatio = omegaFluid / omega;

        Array<T, SymmetricTensorImpl<T, D::d>::n> Psi;
        Psi[S::xx] = phi * piNeq[S::xx] * omegaRatio + 2 * jNeq[0] * u[0];
        Psi[S::xy] = phi * piNeq[S::xy] * omegaRatio + u[0] * jNeq[1] + u[1] * jNeq[0];
        Psi[S::xz] = phi * piNeq[S::xz] * omegaRatio + jNeq[2] * u[0] + jNeq[0] * u[2];
        Psi[S::yy] = phi * piNeq[S::yy] * omegaRatio + 2 * jNeq[1] * u[1];
        Psi[S::yz] = phi * piNeq[S::yz] * omegaRatio + jNeq[2] * u[1] + jNeq[1] * u[2];
        Psi[S::zz] = phi * piNeq[S::zz] * omegaRatio + 2 * jNeq[2] * u[2];

        fNeq[0] = D::t[0] * (-1.5 * Psi[S::xx] - 1.5 * Psi[S::yy] - 1.5 * Psi[S::zz]);
        fNeq[1] = D::t[1] * (-3 * jNeq[0] + 3 * Psi[S::xx] - 1.5 * Psi[S::yy] - 1.5 * Psi[S::zz]);
        fNeq[2] = D::t[2] * (-3 * jNeq[1] - 1.5 * Psi[S::xx] + 3 * Psi[S::yy] - 1.5 * Psi[S::zz]);
        fNeq[3] = D::t[3] * (-3 * jNeq[2] - 1.5 * Psi[S::xx] - 1.5 * Psi[S::yy] + 3 * Psi[S::zz]);
        fNeq[4] = D::t[4]
                  * (-3 * jNeq[0] - 3 * jNeq[1] - 3 * jNeq[2] + 3 * Psi[S::xx] + 9 * Psi[S::xy]
                     + 9 * Psi[S::xz] + 3 * Psi[S::yy] + 9 * Psi[S::yz] + 3 * Psi[S::zz]);
        fNeq[5] = D::t[5]
                  * (-3 * jNeq[0] - 3 * jNeq[1] + 3 * jNeq[2] + 3 * Psi[S::xx] + 9 * Psi[S::xy]
                     - 9 * Psi[S::xz] + 3 * Psi[S::yy] - 9 * Psi[S::yz] + 3 * Psi[S::zz]);
        fNeq[6] = D::t[6]
                  * (-3 * jNeq[0] + 3 * jNeq[1] - 3 * jNeq[2] + 3 * Psi[S::xx] - 9 * Psi[S::xy]
                     + 9 * Psi[S::xz] + 3 * Psi[S::yy] - 9 * Psi[S::yz] + 3 * Psi[S::zz]);
        fNeq[7] = D::t[7]
                  * (-3 * jNeq[0] + 3 * jNeq[1] + 3 * jNeq[2] + 3 * Psi[S::xx] - 9 * Psi[S::xy]
                     - 9 * Psi[S::xz] + 3 * Psi[S::yy] + 9 * Psi[S::yz] + 3 * Psi[S::zz]);
        fNeq[8] = D::t[8] * (3 * jNeq[0] + 3 * Psi[S::xx] - 1.5 * Psi[S::yy] - 1.5 * Psi[S::zz]);
        fNeq[9] = D::t[9] * (3 * jNeq[1] - 1.5 * Psi[S::xx] + 3 * Psi[S::yy] - 1.5 * Psi[S::zz]);
        fNeq[10] = D::t[10] * (3 * jNeq[2] - 1.5 * Psi[S::xx] - 1.5 * Psi[S::yy] + 3 * Psi[S::zz]);
        fNeq[11] = D::t[11]
                   * (3 * jNeq[0] + 3 * jNeq[1] + 3 * jNeq[2] + 3 * Psi[S::xx] + 9 * Psi[S::xy]
                      + 9 * Psi[S::xz] + 3 * Psi[S::yy] + 9 * Psi[S::yz] + 3 * Psi[S::zz]);
        fNeq[12] = D::t[12]
                   * (3 * jNeq[0] + 3 * jNeq[1] - 3 * jNeq[2] + 3 * Psi[S::xx] + 9 * Psi[S::xy]
                      - 9 * Psi[S::xz] + 3 * Psi[S::yy] - 9 * Psi[S::yz] + 3 * Psi[S::zz]);
        fNeq[13] = D::t[13]
                   * (3 * jNeq[0] - 3 * jNeq[1] + 3 * jNeq[2] + 3 * Psi[S::xx] - 9 * Psi[S::xy]
                      + 9 * Psi[S::xz] + 3 * Psi[S::yy] - 9 * Psi[S::yz] + 3 * Psi[S::zz]);
        fNeq[14] = D::t[14]
                   * (3 * jNeq[0] - 3 * jNeq[1] - 3 * jNeq[2] + 3 * Psi[S::xx] - 9 * Psi[S::xy]
                      - 9 * Psi[S::xz] + 3 * Psi[S::yy] + 9 * Psi[S::yz] + 3 * Psi[S::zz]);
    }

    static void complete_bgk_ma2_regularize(
        Array<T, D::q> &f, T rhoPhiBar, T rhoBar, Array<T, D::d> const &jEq,
        Array<T, D::d> const &jNeq, const Array<T, SymmetricTensorImpl<T, D::d>::n> &piNeq, T omega,
        T omegaNonPhys, T omegaFluid, T omegaFluidNonPhys)
    {  // TODO: should check omegaNonPhys, and omegaNonPhysFluid that are unused. Looks suspicious.
        T invRho = D::invRho(rhoBar);
        T phi = invRho * D::fullRho(rhoPhiBar);
        T invRhoPhi = D::invRho(rhoPhiBar);
        T jSqr = jEq[0] * jEq[0] + jEq[1] * jEq[1] + jEq[2] * jEq[2];
        dynamicsTemplatesImpl<T, D>::bgk_ma2_equilibria(rhoPhiBar, invRhoPhi, jEq, jSqr, f);
        Array<T, D::q> fNeq;
        bgk_ma2_off_equilibra(phi, jEq * invRhoPhi, jNeq, piNeq, omega, omegaFluid, fNeq);

        for (plint iPop = 0; iPop < D::q; ++iPop) {
            f[iPop] += fNeq[iPop];
        }
    }

    static T complete_bgk_ma2_regularized_collision(
        Array<T, D::q> &f, T rhoPhiBar, T rhoBar, Array<T, D::d> const &jEq,
        Array<T, D::d> const &jNeq, const Array<T, SymmetricTensorImpl<T, D::d>::n> &piNeq, T omega,
        T omegaNonPhys, T omegaFluid, T omegaFluidNonPhys)
    {
        T invRho = D::invRho(rhoBar);
        T phi = invRho * D::fullRho(rhoPhiBar);
        T invRhoPhi = D::invRho(rhoPhiBar);
        T jSqr = jEq[0] * jEq[0] + jEq[1] * jEq[1] + jEq[2] * jEq[2];
        dynamicsTemplatesImpl<T, D>::bgk_ma2_equilibria(rhoPhiBar, invRhoPhi, jEq, jSqr, f);
        Array<T, D::q> fNeq;
        bgk_ma2_off_equilibra(
            phi, jEq * invRhoPhi, jNeq, piNeq, omega, omegaNonPhys, omegaFluid, omegaFluidNonPhys,
            fNeq);

        for (plint iPop = 0; iPop < D::q; ++iPop) {
            f[iPop] += ((T)1 - omega) * fNeq[iPop];
        }
        return jSqr * invRhoPhi * invRhoPhi;
    }

};  // struct advectionDiffusionDynamicsTemplatesImpl

// TODO: It seems there is the error that the advection-diffusion templates
//       for D3Q15 and D3Q19 have been specialized only partially.
//       The with-source version is missing in the specialization of the RLB
//       case at least. For codes with D3Q19 to compile, we temporarily
//       comment-out this specialization for now. At some point we must add
//       back the RLB-with-source case into the D3Q19 specialization (and in
//       any other it might be missing from).

template <typename T>
struct advectionDiffusionDynamicsTemplatesImpl<T, descriptors::D3Q19DescriptorBase<T> > {
    typedef descriptors::D3Q19DescriptorBase<T> D;

    static T computePsiComplete(T omega)
    {
        PLB_ASSERT(false && "method not implemented for d3q19");
    }

    static void bgk_ma2_off_equilibra(
        T phi, Array<T, D::d> const &u, Array<T, D::d> const &jNeq,
        const Array<T, SymmetricTensorImpl<T, D::d>::n> &piNeq, T omega, T omegaFluid,
        Array<T, D::q> &fNeq)
    {
        typedef SymmetricTensorImpl<T, D::d> S;

        T omegaRatio = omegaFluid / omega;

        Array<T, SymmetricTensorImpl<T, D::d>::n> Psi;
        Psi[S::xx] = phi * piNeq[S::xx] * omegaRatio + 2 * jNeq[0] * u[0];
        Psi[S::xy] = phi * piNeq[S::xy] * omegaRatio + u[0] * jNeq[1] + u[1] * jNeq[0];
        Psi[S::xz] = phi * piNeq[S::xz] * omegaRatio + jNeq[2] * u[0] + jNeq[0] * u[2];
        Psi[S::yy] = phi * piNeq[S::yy] * omegaRatio + 2 * jNeq[1] * u[1];
        Psi[S::yz] = phi * piNeq[S::yz] * omegaRatio + jNeq[2] * u[1] + jNeq[1] * u[2];
        Psi[S::zz] = phi * piNeq[S::zz] * omegaRatio + 2 * jNeq[2] * u[2];

        fNeq[0] = D::t[0] * (-(T)1.5 * Psi[S::xx] - (T)1.5 * Psi[S::yy] - (T)1.5 * Psi[S::zz]);
        fNeq[1] =
            D::t[1] * (-3 * jNeq[0] + 3 * Psi[S::xx] - (T)1.5 * Psi[S::yy] - (T)1.5 * Psi[S::zz]);
        fNeq[2] =
            D::t[2] * (-3 * jNeq[1] - (T)1.5 * Psi[S::xx] + 3 * Psi[S::yy] - (T)1.5 * Psi[S::zz]);
        fNeq[3] =
            D::t[3] * (-3 * jNeq[2] - (T)1.5 * Psi[S::xx] - (T)1.5 * Psi[S::yy] + 3 * Psi[S::zz]);
        fNeq[4] = D::t[4]
                  * (-3 * jNeq[0] - 3 * jNeq[1] + 3 * Psi[S::xx] + 9 * Psi[S::xy] + 3 * Psi[S::yy]
                     - (T)1.5 * Psi[S::zz]);
        fNeq[5] = D::t[5]
                  * (-3 * jNeq[0] + 3 * jNeq[1] + 3 * Psi[S::xx] - 9 * Psi[S::xy] + 3 * Psi[S::yy]
                     - (T)1.5 * Psi[S::zz]);
        fNeq[6] = D::t[6]
                  * (-3 * jNeq[0] - 3 * jNeq[2] + 3 * Psi[S::xx] + 9 * Psi[S::xz]
                     - (T)1.5 * Psi[S::yy] + 3 * Psi[S::zz]);
        fNeq[7] = D::t[7]
                  * (-3 * jNeq[0] + 3 * jNeq[2] + 3 * Psi[S::xx] - 9 * Psi[S::xz]
                     - (T)1.5 * Psi[S::yy] + 3 * Psi[S::zz]);
        fNeq[8] = D::t[8]
                  * (-3 * jNeq[1] - 3 * jNeq[2] - (T)1.5 * Psi[S::xx] + 3 * Psi[S::yy]
                     + 9 * Psi[S::yz] + 3 * Psi[S::zz]);
        fNeq[9] = D::t[9]
                  * (-3 * jNeq[1] + 3 * jNeq[2] - (T)1.5 * Psi[S::xx] + 3 * Psi[S::yy]
                     - 9 * Psi[S::yz] + 3 * Psi[S::zz]);
        fNeq[10] =
            D::t[10] * (3 * jNeq[0] + 3 * Psi[S::xx] - (T)1.5 * Psi[S::yy] - (T)1.5 * Psi[S::zz]);
        fNeq[11] =
            D::t[11] * (3 * jNeq[1] - (T)1.5 * Psi[S::xx] + 3 * Psi[S::yy] - (T)1.5 * Psi[S::zz]);
        fNeq[12] =
            D::t[12] * (3 * jNeq[2] - (T)1.5 * Psi[S::xx] - (T)1.5 * Psi[S::yy] + 3 * Psi[S::zz]);
        fNeq[13] = D::t[13]
                   * (3 * jNeq[0] + 3 * jNeq[1] + 3 * Psi[S::xx] + 9 * Psi[S::xy] + 3 * Psi[S::yy]
                      - (T)1.5 * Psi[S::zz]);
        fNeq[14] = D::t[14]
                   * (3 * jNeq[0] - 3 * jNeq[1] + 3 * Psi[S::xx] - 9 * Psi[S::xy] + 3 * Psi[S::yy]
                      - (T)1.5 * Psi[S::zz]);
        fNeq[15] = D::t[15]
                   * (3 * jNeq[0] + 3 * jNeq[2] + 3 * Psi[S::xx] + 9 * Psi[S::xz]
                      - (T)1.5 * Psi[S::yy] + 3 * Psi[S::zz]);
        fNeq[16] = D::t[16]
                   * (3 * jNeq[0] - 3 * jNeq[2] + 3 * Psi[S::xx] - 9 * Psi[S::xz]
                      - (T)1.5 * Psi[S::yy] + 3 * Psi[S::zz]);
        fNeq[17] = D::t[17]
                   * (3 * jNeq[1] + 3 * jNeq[2] - (T)1.5 * Psi[S::xx] + 3 * Psi[S::yy]
                      + 9 * Psi[S::yz] + 3 * Psi[S::zz]);
        fNeq[18] = D::t[18]
                   * (3 * jNeq[1] - 3 * jNeq[2] - (T)1.5 * Psi[S::xx] + 3 * Psi[S::yy]
                      - 9 * Psi[S::yz] + 3 * Psi[S::zz]);
    }

    static T bgk_ma1_equilibrium(plint iPop, T rhoBar, Array<T, D::d> const &jEq)
    {
        return D::t[iPop]
               * (rhoBar
                  + D::invCs2
                        * (D::c[iPop][0] * jEq[0] + D::c[iPop][1] * jEq[1]
                           + D::c[iPop][2] * jEq[2]));
    }

    static void regularize(
        Array<T, D::q> &f, T rhoBar, Array<T, D::d> const &jAdvDiff, Array<T, D::d> const &jEq)
    {
        T invRho = D::invRho(rhoBar);
        T jSqr = jEq[0] * jEq[0] + jEq[1] * jEq[1] + jEq[2] * jEq[2];
        dynamicsTemplatesImpl<T, D>::bgk_ma2_equilibria(rhoBar, invRho, jEq, jSqr, f);
        Array<T, D::d> jNeq = jAdvDiff - jEq;
        for (plint iPop = 0; iPop < D::q; ++iPop) {
            f[iPop] +=
                D::t[iPop] * D::invCs2
                * (D::c[iPop][0] * jNeq[0] + D::c[iPop][1] * jNeq[1] + D::c[iPop][2] * jNeq[2]);
        }
    }

    static T no_corr_bgk_collision(Array<T, D::q> &f, T rhoBar, Array<T, D::d> const &jEq, T omega)
    {
        T invRho = D::invRho(rhoBar);
        T jSqr = dynamicsTemplatesImpl<T, D>::bgk_ma2_collision(f, rhoBar, jEq, omega);

        return jSqr * invRho * invRho;
    }

    static T no_corr_bgk_collision(
        Array<T, D::q> &f, T rhoBar, Array<T, D::d> const &jEq, T omega, T source)
    {
        PLB_ASSERT(
            false && "no correction bgk collision with source not implemented yet for d3q19.");
    }

    static T no_corr_rlb_collision(
        Array<T, D::q> &f, T rhoBar, Array<T, D::d> const &jEq, Array<T, D::d> const &jNeq, T omega)
    {
        T invRho = D::invRho(rhoBar);
        T jSqr = jEq[0] * jEq[0] + jEq[1] * jEq[1] + jEq[2] * jEq[2];
        dynamicsTemplatesImpl<T, D>::bgk_ma2_equilibria(rhoBar, invRho, jEq, jSqr, f);
        for (plint iPop = 0; iPop < D::q; ++iPop) {
            f[iPop] +=
                ((T)1 - omega)
                * (D::t[iPop] * D::invCs2
                   * (D::c[iPop][0] * jNeq[0] + D::c[iPop][1] * jNeq[1] + D::c[iPop][2] * jNeq[2]));
        }

        return jSqr * invRho * invRho;
    }

    static T no_corr_rlb_collision(
        Array<T, D::q> &f, T rhoBar, Array<T, D::d> const &jEq, Array<T, D::d> const &jNeq, T omega,
        T source)
    {
        T invRho = D::invRho(rhoBar);
        T jSqr = jEq[0] * jEq[0] + jEq[1] * jEq[1] + jEq[2] * jEq[2];
        dynamicsTemplatesImpl<T, D>::bgk_ma2_equilibria(rhoBar, invRho, jEq, jSqr, f);
        for (plint iPop = 0; iPop < D::q; ++iPop) {
            f[iPop] += ((T)1 - omega)
                           * (D::t[iPop] * D::invCs2
                              * (D::c[iPop][0] * jNeq[0] + D::c[iPop][1] * jNeq[1]
                                 + D::c[iPop][2] * jNeq[2]))
                       + D::t[iPop] * source;
        }

        return jSqr * invRho * invRho;
    }

    static void complete_bgk_ma2_regularize(
        Array<T, D::q> &f, T rhoPhiBar, T rhoBar, Array<T, D::d> const &jEq,
        Array<T, D::d> const &jNeq, const Array<T, SymmetricTensorImpl<T, D::d>::n> &piNeq, T omega,
        T omegaNonPhys, T omegaFluid, T omegaFluidNonPhys)
    {  // TODO: check both unused seem suspicious (and misleading)
        T invRho = D::invRho(rhoBar);
        T phi = invRho * D::fullRho(rhoPhiBar);
        T invRhoPhi = D::invRho(rhoPhiBar);
        T jSqr = jEq[0] * jEq[0] + jEq[1] * jEq[1] + jEq[2] * jEq[2];
        dynamicsTemplatesImpl<T, D>::bgk_ma2_equilibria(rhoPhiBar, invRhoPhi, jEq, jSqr, f);
        Array<T, D::q> fNeq;
        bgk_ma2_off_equilibra(phi, jEq * invRhoPhi, jNeq, piNeq, omega, omegaFluid, fNeq);

        for (plint iPop = 0; iPop < D::q; ++iPop) {
            f[iPop] += fNeq[iPop];
        }
    }

    static T complete_bgk_ma2_regularized_collision(
        Array<T, D::q> &f, T rhoPhiBar, T rhoBar, Array<T, D::d> const &jEq,
        Array<T, D::d> const &jNeq, const Array<T, SymmetricTensorImpl<T, D::d>::n> &piNeq, T omega,
        T omegaNonPhys, T omegaFluid, T omegaFluidNonPhys)
    {  // TODO: check both unused omegas. Misleading.
        T invRho = D::invRho(rhoBar);
        T phi = invRho * D::fullRho(rhoPhiBar);
        T invRhoPhi = D::invRho(rhoPhiBar);
        T jSqr = jEq[0] * jEq[0] + jEq[1] * jEq[1] + jEq[2] * jEq[2];
        dynamicsTemplatesImpl<T, D>::bgk_ma2_equilibria(rhoPhiBar, invRhoPhi, jEq, jSqr, f);
        Array<T, D::q> fNeq;
        bgk_ma2_off_equilibra(phi, jEq * invRhoPhi, jNeq, piNeq, omega, omegaFluid, fNeq);

        for (plint iPop = 0; iPop < D::q; ++iPop) {
            f[iPop] += ((T)1 - omega) * fNeq[iPop];
        }
        return jSqr * invRhoPhi * invRhoPhi;
    }

};  // struct advectionDiffusionDynamicsTemplatesImpl

}  // namespace plb

#endif
