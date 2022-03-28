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
#ifndef ADVECTION_DIFFUSION_DYNAMICS_TEMPLATES_2D_H
#define ADVECTION_DIFFUSION_DYNAMICS_TEMPLATES_2D_H

namespace plb {

/// All helper functions are inside this structure
template <typename T>
struct advectionDiffusionDynamicsTemplatesImpl<T, descriptors::D2Q5DescriptorBase<T> > {
    typedef descriptors::D2Q5DescriptorBase<T> Descriptor;

    static T bgk_ma1_equilibrium(plint iPop, T rhoBar, Array<T, Descriptor::d> const &jEq)
    {
        return Descriptor::t[iPop]
               * (rhoBar
                  + Descriptor::invCs2
                        * (Descriptor::c[iPop][0] * jEq[0] + Descriptor::c[iPop][1] * jEq[1]));
    }

    /// Regularization
    static void regularize(
        Array<T, Descriptor::q> &f, T rhoBar, Array<T, Descriptor::d> const &jAdvDiff,
        Array<T, Descriptor::d> const &jEq)
    {
        f[0] = Descriptor::t[0] * rhoBar;

        f[1] = Descriptor::t[1] * (rhoBar - Descriptor::invCs2 * jAdvDiff[0]);
        f[2] = Descriptor::t[2] * (rhoBar - Descriptor::invCs2 * jAdvDiff[1]);
        f[3] = Descriptor::t[3] * (rhoBar + Descriptor::invCs2 * jAdvDiff[0]);
        f[4] = Descriptor::t[4] * (rhoBar + Descriptor::invCs2 * jAdvDiff[1]);
    }

    static T no_corr_bgk_collision(
        Array<T, Descriptor::q> &f, T rhoBar, Array<T, Descriptor::d> const &jEq, T omega)
    {
        T invRho = Descriptor::invRho(rhoBar);
        const T jSqr = jEq[0] * jEq[0] + jEq[1] * jEq[1];

        const T oneMinusOmega = (T)1 - omega;
        const T halfOmega = (T)0.5 * omega;
        const T cs2RhoBar = Descriptor::cs2 * rhoBar;

        f[0] = oneMinusOmega * f[0] + omega * ((T)1 - 2 * Descriptor::cs2) * rhoBar;

        f[1] = oneMinusOmega * f[1] + halfOmega * (cs2RhoBar - jEq[0]);
        f[2] = oneMinusOmega * f[2] + halfOmega * (cs2RhoBar - jEq[1]);
        f[3] = oneMinusOmega * f[3] + halfOmega * (cs2RhoBar + jEq[0]);
        f[4] = oneMinusOmega * f[4] + halfOmega * (cs2RhoBar + jEq[1]);

        return jSqr * invRho * invRho;
    }

    static T no_corr_bgk_collision(
        Array<T, Descriptor::q> &f, T rhoBar, Array<T, Descriptor::d> const &jEq, T omega, T source)
    {
        T invRho = Descriptor::invRho(rhoBar);
        const T jSqr = jEq[0] * jEq[0] + jEq[1] * jEq[1];

        const T oneMinusOmega = (T)1 - omega;
        const T halfOmega = (T)0.5 * omega;
        const T cs2RhoBar = Descriptor::cs2 * rhoBar;
        const T halfSourceCs2 = (T)0.5 * source * Descriptor::cs2;

        f[0] = oneMinusOmega * f[0] + ((T)1 - 2 * Descriptor::cs2) * (source + omega * rhoBar);

        f[1] = oneMinusOmega * f[1] + halfOmega * (cs2RhoBar - jEq[0]) + halfSourceCs2;
        f[2] = oneMinusOmega * f[2] + halfOmega * (cs2RhoBar - jEq[1]) + halfSourceCs2;
        f[3] = oneMinusOmega * f[3] + halfOmega * (cs2RhoBar + jEq[0]) + halfSourceCs2;
        f[4] = oneMinusOmega * f[4] + halfOmega * (cs2RhoBar + jEq[1]) + halfSourceCs2;

        return jSqr * invRho * invRho;
    }

    static T no_corr_rlb_collision(
        Array<T, Descriptor::q> &f, T rhoBar, Array<T, Descriptor::d> const &jEq,
        Array<T, Descriptor::d> const &jNeq, T omega)
    {
        T invRho = Descriptor::invRho(rhoBar);
        const T jSqr = jEq[0] * jEq[0] + jEq[1] * jEq[1];

        const T oneHalfMinusHalfOmega = (T)0.5 - (T)0.5 * omega;
        const T cs2RhoBar = Descriptor::cs2 * rhoBar;

        const T jNeqTerm_0 = oneHalfMinusHalfOmega * jNeq[0];
        const T jNeqTerm_1 = oneHalfMinusHalfOmega * jNeq[1];

        f[0] = ((T)1 - (T)2 * Descriptor::cs2) * rhoBar;

        f[1] = -jNeqTerm_0 + (T)0.5 * (cs2RhoBar - jEq[0]);
        f[2] = -jNeqTerm_1 + (T)0.5 * (cs2RhoBar - jEq[1]);

        f[3] = -f[1] + cs2RhoBar;
        f[4] = -f[2] + cs2RhoBar;

        return jSqr * invRho * invRho;
    }

    static T no_corr_rlb_collision(
        Array<T, Descriptor::q> &f, T rhoBar, Array<T, Descriptor::d> const &jEq,
        Array<T, Descriptor::d> const &jNeq, T omega, T source)
    {
        T invRho = Descriptor::invRho(rhoBar);
        const T jSqr = jEq[0] * jEq[0] + jEq[1] * jEq[1];

        const T oneHalfMinusHalfOmega = (T)0.5 - (T)0.5 * omega;
        const T cs2RhoBar = Descriptor::cs2 * rhoBar;
        const T halfSourceCs2 = (T)0.5 * source * Descriptor::cs2;

        const T jNeqTerm_0 = oneHalfMinusHalfOmega * jNeq[0];
        const T jNeqTerm_1 = oneHalfMinusHalfOmega * jNeq[1];

        f[0] = ((T)1 - (T)2 * Descriptor::cs2) * (rhoBar + source);

        f[1] = -jNeqTerm_0 + (T)0.5 * (cs2RhoBar - jEq[0]);
        f[2] = -jNeqTerm_1 + (T)0.5 * (cs2RhoBar - jEq[1]);

        f[3] = -f[1] + cs2RhoBar;
        f[4] = -f[2] + cs2RhoBar;

        f[1] += halfSourceCs2;
        f[2] += halfSourceCs2;
        f[3] += halfSourceCs2;
        f[4] += halfSourceCs2;

        return jSqr * invRho * invRho;
    }

};  // struct advectionDiffusionDynamicsTemplatesImpl

template <typename T>
struct advectionDiffusionDynamicsTemplatesImpl<T, descriptors::D2Q9DescriptorBase<T> > {
    typedef descriptors::D2Q9DescriptorBase<T> D;

    static T computePsiComplete(T omega)
    {
        T psi = T();
        if (omega > (T)1.75) {
            psi = -629.1520418797601 - 562.5634533117104 * omega * omega
                  + 102.9178439049484 * omega * omega * omega + 1029.031566819295 * omega;
        } else {
            psi = -2.400317064381656 * omega * omega + .6470123480620908 * omega * omega * omega
                  + 2.435047526713970 * omega;
        }
        if (psi > (T)2)
            psi = (T)2;
        return psi;
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
        T jSqr = jEq[0] * jEq[0] + jEq[1] * jEq[1];
        dynamicsTemplatesImpl<T, D>::bgk_ma2_equilibria(rhoBar, invRho, jEq, jSqr, f);
        Array<T, D::d> jNeq = jAdvDiff - jEq;
        for (plint iPop = 0; iPop < D::q; ++iPop) {
            f[iPop] += D::t[iPop] * D::invCs2 * (D::c[iPop][0] * jNeq[0] + D::c[iPop][1] * jNeq[1]);
        }
    }

    static T no_corr_bgk_collision(Array<T, D::q> &f, T rhoBar, Array<T, D::d> const &jEq, T omega)
    {
        T invRho = D::invRho(rhoBar);
        T jSqr = dynamicsTemplatesImpl<T, D>::bgk_ma2_collision(f, rhoBar, jEq, omega);

        return jSqr * invRho * invRho;
    }

    static T no_corr_rlb_collision(
        Array<T, D::q> &f, T rhoBar, Array<T, D::d> const &jEq, Array<T, D::d> const &jNeq, T omega)
    {
        T invRho = D::invRho(rhoBar);
        T jSqr = jEq[0] * jEq[0] + jEq[1] * jEq[1];
        dynamicsTemplatesImpl<T, D>::bgk_ma2_equilibria(rhoBar, invRho, jEq, jSqr, f);
        for (plint iPop = 0; iPop < D::q; ++iPop) {
            f[iPop] +=
                ((T)1 - omega)
                * (D::t[iPop] * D::invCs2 * (D::c[iPop][0] * jNeq[0] + D::c[iPop][1] * jNeq[1]));
        }

        return jSqr * invRho * invRho;
    }

    // TODO: this omegaFluidNonPhys looks suspicious. Should check it.
    static void bgk_ma2_off_equilibra(
        T phi, Array<T, D::d> const &u, Array<T, D::d> const &jNeq,
        const Array<T, SymmetricTensorImpl<T, D::d>::n> &piNeq, T omega, T omegaNonPhys,
        T omegaFluid, T omegaFluidNonPhys, Array<T, D::q> &fNeq)
    {
        T ux2 = u[0] * u[0];
        T uy2 = u[1] * u[1];
        Array<T, SymmetricTensorImpl<T, D::d>::n> psiNeq;

        psiNeq[0] =
            2 * omega * u[0] * jNeq[0] / omegaNonPhys + phi * omegaFluid * piNeq[0] / omegaNonPhys;
        psiNeq[1] = omega * u[0] * jNeq[1] / omegaNonPhys + omega * u[1] * jNeq[0] / omegaNonPhys
                    + phi * omegaFluid * piNeq[1] / omegaNonPhys;
        psiNeq[2] =
            2 * omega * u[1] * jNeq[1] / omegaNonPhys + phi * omegaFluid * piNeq[2] / omegaNonPhys;

        Array<T, D::d> q;
        // q[0] = psiNeq[1]*u[1] + psiNeq[2]*u[0] - u[0]*u[1]*jNeq[1];
        q[0] = -2 * u[0] * u[1] * jNeq[1] - uy2 * jNeq[0] + u[0] * psiNeq[2] + 2 * u[1] * psiNeq[1];

        // q[1] = psiNeq[1]*u[0] + psiNeq[0]*u[1] - u[0]*u[1]*jNeq[0];
        q[1] = -ux2 * jNeq[1] - 2 * u[0] * u[1] * jNeq[0] + 2 * u[0] * psiNeq[1] + u[1] * psiNeq[0];
        T chi = q[1] * u[0] + q[0] * u[1] - ux2 * u[1] * jNeq[1] - u[0] * uy2 * jNeq[0];

        fNeq[0] = -(T)2 * D::cs2 * psiNeq[0] - (T)2 * D::cs2 * psiNeq[2] + chi;

        fNeq[2] = -D::cs2 * jNeq[0] + D::cs2 * psiNeq[0] - (T)0.5 * D::cs2 * psiNeq[2]
                  + (T)0.5 * q[0] - (T)0.5 * chi;
        fNeq[6] = D::cs2 * jNeq[0] + D::cs2 * psiNeq[0] - (T)0.5 * D::cs2 * psiNeq[2]
                  - (T)0.5 * q[0] - (T)0.5 * chi;

        fNeq[4] = -D::cs2 * jNeq[1] - (T)0.5 * D::cs2 * psiNeq[0] + D::cs2 * psiNeq[2]
                  + (T)0.5 * q[1] - (T)0.5 * chi;
        fNeq[8] = D::cs2 * jNeq[1] - (T)0.5 * D::cs2 * psiNeq[0] + D::cs2 * psiNeq[2]
                  - (T)0.5 * q[1] - (T)0.5 * chi;

        fNeq[1] = -(T)0.25 * D::cs2 * jNeq[0] + (T)0.25 * D::cs2 * jNeq[1]
                  + (T)0.25 * D::cs2 * psiNeq[0] - (T)0.25 * psiNeq[1]
                  + (T)0.25 * D::cs2 * psiNeq[2] - (T)0.25 * q[0] + (T)0.25 * q[1] + (T)0.25 * chi;
        fNeq[5] = (T)0.25 * D::cs2 * jNeq[0] - (T)0.25 * D::cs2 * jNeq[1]
                  + (T)0.25 * D::cs2 * psiNeq[0] - (T)0.25 * psiNeq[1]
                  + (T)0.25 * D::cs2 * psiNeq[2] + (T)0.25 * q[0] - (T)0.25 * q[1] + (T)0.25 * chi;

        fNeq[3] = -(T)0.25 * D::cs2 * jNeq[0] - (T)0.25 * D::cs2 * jNeq[1]
                  + (T)0.25 * D::cs2 * psiNeq[0] + (T)0.25 * psiNeq[1]
                  + (T)0.25 * D::cs2 * psiNeq[2] - (T)0.25 * q[0] - (T)0.25 * q[1] + (T)0.25 * chi;
        fNeq[7] = (T)0.25 * D::cs2 * jNeq[0] + (T)0.25 * D::cs2 * jNeq[1]
                  + (T)0.25 * D::cs2 * psiNeq[0] + (T)0.25 * psiNeq[1]
                  + (T)0.25 * D::cs2 * psiNeq[2] + (T)0.25 * q[0] + (T)0.25 * q[1] + (T)0.25 * chi;
    }

    static void complete_bgk_ma2_regularize(
        Array<T, D::q> &f, T rhoPhiBar, T rhoBar, Array<T, D::d> const &jEq,
        Array<T, D::d> const &jNeq, const Array<T, SymmetricTensorImpl<T, D::d>::n> &piNeq, T omega,
        T omegaNonPhys, T omegaFluid, T omegaFluidNonPhys)
    {
        T invRhoPhi = D::invRho(rhoPhiBar);
        T jEqSqr = jEq[0] * jEq[0] + jEq[1] * jEq[1];
        dynamicsTemplatesImpl<T, D>::complete_bgk_ma2_equilibria(
            rhoPhiBar, invRhoPhi, jEq, jEqSqr, f);

        T phi = D::fullRho(rhoPhiBar) * D::invRho(rhoBar);

        Array<T, D::d> u = jEq * invRhoPhi;
        Array<T, D::q> fNeq;
        bgk_ma2_off_equilibra(
            phi, u, jNeq, piNeq, omega, omegaNonPhys, omegaFluid, omegaFluidNonPhys, fNeq);

        f += fNeq;
    }

    static T complete_bgk_ma2_regularized_collision(
        Array<T, D::q> &f, T rhoPhiBar, T rhoBar, Array<T, D::d> const &jEq,
        Array<T, D::d> const &jNeq, const Array<T, SymmetricTensorImpl<T, D::d>::n> &piNeq, T omega,
        T omegaNonPhys, T omegaFluid, T omegaFluidNonPhys)
    {
        T invRho = D::invRho(rhoBar);
        T phi = invRho * D::fullRho(rhoPhiBar);
        T invRhoPhi = D::invRho(rhoPhiBar);
        T jSqr = jEq[0] * jEq[0] + jEq[1] * jEq[1];
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

}  // namespace plb

#endif
