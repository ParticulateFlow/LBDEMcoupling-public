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
 * functions of the header file mrtTemplates.h, for the D2Q9 grid.
 */

#ifndef MRT_TEMPLATES_2D_H
#define MRT_TEMPLATES_2D_H

#include "core/globalDefs.h"

namespace plb {

// Efficient specialization for D2Q9 lattice
template <typename T>
struct mrtTemplatesImpl<T, descriptors::MRTD2Q9DescriptorBase<T> > {
    typedef descriptors::D2Q9DescriptorBase<T> Descriptor;
    typedef descriptors::MRTD2Q9DescriptorBase<T> MRTDescriptor;

    /// Computation of all equilibrium distribution (in moments space)
    static void computeEquilibriumMoments(
        Array<T, Descriptor::q> &momentsEq, T rhoBar, Array<T, 2> const &j, T jSqr)
    {
        T invRho = Descriptor::invRho(rhoBar);
        momentsEq[0] = rhoBar;
        momentsEq[1] = (T)3 * jSqr * invRho - 2 * rhoBar;
        momentsEq[2] = -(T)3 * jSqr * invRho + rhoBar;
        momentsEq[3] = j[0];
        momentsEq[4] = -j[0];
        momentsEq[5] = j[1];
        momentsEq[6] = -j[1];
        momentsEq[7] = (j[0] * j[0] - j[1] * j[1]) * invRho;
        momentsEq[8] = j[1] * j[0] * invRho;
    }

    /// Computation of all moments (specialized for d2q9)
    static void computeMoments(Array<T, Descriptor::q> &moments, const Array<T, Descriptor::q> &f)
    {
        T f1_f3 = f[1] - f[3];
        T f1pf3 = f[1] + f[3];

        T f2pf4 = f[2] + f[4];
        T f5_f7 = f[5] - f[7];
        T f5pf7 = f[5] + f[7];
        T f6pf8 = f[6] + f[8];

        moments[0] = f[0] + f1pf3 + f2pf4 + f5pf7 + f6pf8;
        moments[1] = -(T)4 * f[0] + (T)2 * (f1pf3 + f5pf7) - f2pf4 - f6pf8;
        moments[2] = (T)4 * f[0] + f1pf3 - (T)2 * (f2pf4 + f6pf8) + f5pf7;
        moments[3] = -f1pf3 - f[2] + f5pf7 + f[6];
        moments[4] = -f1pf3 + (T)2 * (f[2] - f[6]) + f5pf7;
        moments[5] = f1_f3 - f[4] - f5_f7 + f[8];
        moments[6] = f1_f3 + (T)2 * (f[4] - f[8]) - f5_f7;
        moments[7] = f[2] - f[4] + f[6] - f[8];
        moments[8] = -f1_f3 - f5_f7;
    }

    static void computeMneqInPlace(Array<T, 9> &moments, const Array<T, 9> &momentsEq)
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
    }

    static void computef_InvM_Smoments(Array<T, 9> &f, const Array<T, 9> &moments, const T &omega)
    {
        T mom0 = moments[0] * MRTDescriptor::S[0];
        T mom1 = moments[1] * MRTDescriptor::S[1];
        T mom2 = moments[2] * MRTDescriptor::S[2];
        T mom3 = moments[3] * MRTDescriptor::S[3];
        T mom4 = moments[4] * MRTDescriptor::S[4];
        T mom5 = moments[5] * MRTDescriptor::S[5];
        T mom6 = moments[6] * MRTDescriptor::S[6];
        T mom7 = moments[7] * omega;
        T mom8 = moments[8] * omega;

        const T oneOverSix = (T)1 / (T)6;
        mom3 *= oneOverSix;
        mom5 *= oneOverSix;
        const T oneOverNine = (T)1 / (T)9;
        mom0 *= oneOverNine;
        const T oneOverTwelve = oneOverSix * (T)0.5;
        const T oneOverEighteen = oneOverNine * (T)0.5;
        const T oneOverThirtysix = oneOverEighteen * (T)0.5;
        mom7 *= (T)0.25;
        mom8 *= (T)0.25;

        f[0] -= mom0 - oneOverNine * mom1 + oneOverNine * mom2;
        T mom1tmp = oneOverThirtysix * mom1;
        T mom2tmp = oneOverEighteen * mom2;
        T mom4tmp = oneOverSix * mom4;
        T mom6tmp = oneOverSix * mom6;
        f[2] -= mom0 - mom1tmp - mom2tmp - mom3 + mom4tmp + mom7;
        f[4] -= mom0 - mom1tmp - mom2tmp - mom5 + mom6tmp - mom7;
        f[6] -= mom0 - mom1tmp - mom2tmp + mom3 - mom4tmp + mom7;
        f[8] -= mom0 - mom1tmp - mom2tmp + mom5 - mom6tmp - mom7;

        mom1tmp = oneOverEighteen * mom1;
        mom2tmp = oneOverThirtysix * mom2;
        mom4tmp = oneOverTwelve * mom4;
        mom6tmp = oneOverTwelve * mom6;
        f[1] -= mom0 + mom1tmp + mom2tmp - mom3 - mom4tmp + mom5 + mom6tmp - mom8;
        f[3] -= mom0 + mom1tmp + mom2tmp - mom3 - mom4tmp - mom5 - mom6tmp + mom8;
        f[5] -= mom0 + mom1tmp + mom2tmp + mom3 + mom4tmp - mom5 - mom6tmp - mom8;
        f[7] -= mom0 + mom1tmp + mom2tmp + mom3 + mom4tmp + mom5 + mom6tmp + mom8;
    }

    /// MRT collision step imposed rhoBar and j
    static T mrtCollision(Array<T, Descriptor::q> &f, const T &omega)
    {
        Array<T, 9> moments, momentsEq;

        computeMoments(moments, f);
        T rhoBar = moments[0];
        Array<T, 2> j(
            moments[MRTDescriptor::momentumIndexes[0]], moments[MRTDescriptor::momentumIndexes[1]]);
        T jSqr = VectorTemplateImpl<T, 2>::normSqr(j);

        computeEquilibriumMoments(momentsEq, rhoBar, j, jSqr);
        computeMneqInPlace(moments, momentsEq);
        computef_InvM_Smoments(f, moments, omega);

        return jSqr;
    }

    /// MRT collision step imposed rhoBar and j
    static T mrtCollision(
        Array<T, Descriptor::q> &f, const T &rhoBar, const Array<T, 2> &j, const T &omega)
    {
        Array<T, 9> moments, momentsEq;

        computeMoments(moments, f);
        T jSqr = VectorTemplateImpl<T, 2>::normSqr(j);

        computeEquilibriumMoments(momentsEq, rhoBar, j, jSqr);
        computeMneqInPlace(moments, momentsEq);
        computef_InvM_Smoments(f, moments, omega);

        return jSqr;
    }

    static void addGuoForce(
        Array<T, Descriptor::q> &f, const Array<T, Descriptor::d> &g,
        Array<T, Descriptor::d> const &u, const T &omega, T amplitude)
    {
        Array<T, Descriptor::q> forceMom, force;

        T gx = amplitude * g[0];
        T gy = amplitude * g[1];
        T g_u = gx * u[0] + gy * u[1];

        forceMom[0] = 0;
        forceMom[1] = (T)6 * g_u;
        forceMom[2] = -(T)6 * g_u;
        forceMom[3] = gx;
        forceMom[4] = -gx;
        forceMom[5] = gy;
        forceMom[6] = -gy;
        forceMom[7] = (T)2 * (gx * u[0] - gy * u[1]);
        forceMom[8] = gx * u[1] + gy * u[0];

        forceMom[0] *= (T)0.5;
        forceMom[1] *= (T)0.5;
        forceMom[2] *= (T)0.5;
        forceMom[3] *= (T)0.5;
        forceMom[4] *= (T)0.5;
        forceMom[5] *= (T)0.5;
        forceMom[6] *= (T)0.5;
        forceMom[7] *= (T)0.5;
        forceMom[8] *= (T)0.5;

        computef_InvM_Smoments(f, forceMom, omega);

        T gxuy_pgyux = (gx * u[1] + gy * u[0]);
        force[0] = -(T)4 * Descriptor::cs2 * g_u;

        T prefactor = (T)0.25 * Descriptor::cs2;
        force[1] = prefactor * ((T)2 * g_u - (T)3 * gxuy_pgyux - gx + gy);
        force[3] = prefactor * ((T)2 * g_u + (T)3 * gxuy_pgyux - gx - gy);
        force[5] = prefactor * ((T)2 * g_u - (T)3 * gxuy_pgyux + gx - gy);
        force[7] = prefactor * ((T)2 * g_u + (T)3 * gxuy_pgyux + gx + gy);

        force[2] = Descriptor::cs2 * ((T)2 * gx * u[0] - gy * u[1] - gx);
        force[4] = -Descriptor::cs2 * (gx * u[0] - 2 * gy * u[1] + gy);
        force[6] = Descriptor::cs2 * ((T)2 * gx * u[0] - gy * u[1] + gx);
        force[8] = -Descriptor::cs2 * (gx * u[0] - 2 * gy * u[1] - gy);

        f[0] += force[0];
        f[1] += force[1];
        f[2] += force[2];
        f[3] += force[3];
        f[4] += force[4];
        f[5] += force[5];
        f[6] += force[6];
        f[7] += force[7];
        f[8] += force[8];
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

        T c_u_0 = 0;  //(0,0)

        T c_u_1 = -uLB[0] + uLB[1];  //{-1, 1}
        T c_u_2 = -uLB[0];           // {-1, 0}
        T c_u_3 = -uLB[0] - uLB[1];  // {-1,-1}
        T c_u_4 = -uLB[1];           // { 0,-1}

        T c_u_5 = uLB[0] - uLB[1];  // { 1,-1}
        T c_u_6 = uLB[0];           //{ 1, 0}

        T c_u_7 = uLB[0] + uLB[1];  // { 1, 1}
        T c_u_8 = uLB[1];           // { 0, 1}

        // common terms fo ug_i calculations
        T uLB0_force0 = uLB[0] * force[0];
        T uLB1_force1 = uLB[1] * force[1];

        T mines_force_0 = (T(-1) - uLB[0]) * force[0];
        T mines_force_1 = (T(-1) - uLB[1]) * force[1];

        T plus_force_0 = (T(1) - uLB[0]) * force[0];
        T plus_force_1 = (T(1) - uLB[1]) * force[1];

        T ug_0 = -uLB0_force0 - uLB1_force1;     //(0,0)
        T ug_1 = mines_force_0 + plus_force_1;   //{-1, 1}
        T ug_2 = mines_force_0 - uLB1_force1;    //{-1, 0},
        T ug_3 = mines_force_0 + mines_force_1;  //{-1,-1}
        T ug_4 = -uLB0_force0 + mines_force_1;   //{ 0,-1}

        T ug_5 = plus_force_0 + mines_force_1;  //{ 1,-1}
        T ug_6 = plus_force_0 - uLB1_force1;    //{ 1, 0},

        T ug_7 = plus_force_0 + plus_force_1;  //{ 1, 1}
        T ug_8 = -uLB0_force0 + plus_force_1;  // { 0, 1}

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

        Array<T, 9> f_full, forceMoments, forcing;

        forcing[0] = invRho * Descriptor::invCs2 * ug_0 * eqContribution_0;
        forcing[1] = invRho * Descriptor::invCs2 * ug_1 * eqContribution_1;
        forcing[2] = invRho * Descriptor::invCs2 * ug_2 * eqContribution_2;
        forcing[3] = invRho * Descriptor::invCs2 * ug_3 * eqContribution_3;
        forcing[4] = invRho * Descriptor::invCs2 * ug_4 * eqContribution_4;
        forcing[5] = invRho * Descriptor::invCs2 * ug_5 * eqContribution_5;
        forcing[6] = invRho * Descriptor::invCs2 * ug_6 * eqContribution_6;
        forcing[7] = invRho * Descriptor::invCs2 * ug_7 * eqContribution_7;
        forcing[8] = invRho * Descriptor::invCs2 * ug_8 * eqContribution_8;

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

        f[0] = f_full[0] - Descriptor::SkordosFactor() * Descriptor::t[0];
        f[1] = f_full[1] - Descriptor::SkordosFactor() * Descriptor::t[1];
        f[2] = f_full[2] - Descriptor::SkordosFactor() * Descriptor::t[2];
        f[3] = f_full[3] - Descriptor::SkordosFactor() * Descriptor::t[3];
        f[4] = f_full[4] - Descriptor::SkordosFactor() * Descriptor::t[4];
        f[5] = f_full[5] - Descriptor::SkordosFactor() * Descriptor::t[5];
        f[6] = f_full[6] - Descriptor::SkordosFactor() * Descriptor::t[6];
        f[7] = f_full[7] - Descriptor::SkordosFactor() * Descriptor::t[7];
        f[8] = f_full[8] - Descriptor::SkordosFactor() * Descriptor::t[8];
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

    /// MRT collision step
    static T incMrtCollisionWithForce(
        Array<T, Descriptor::q> &f, const T &rhoBar, const Array<T, Descriptor::d> &u,
        T invM_S[Descriptor::q][Descriptor::q], const Array<T, Descriptor::d> &force, T amplitude)
    {
        Array<T, Descriptor::d> j = Descriptor::fullRho(rhoBar) * u;
        T jSqr = incMrtCollision(f, rhoBar, j, invM_S);
        addGuoForce(f, force, u, invM_S, amplitude);

        return jSqr;
    }

    /// Computation of all equilibrium distribution (in moments space)
    static void computeIncEquilibriumMoments(
        Array<T, Descriptor::q> &momentsEq, T rhoBar, Array<T, 2> const &j, T jSqr)
    {
        momentsEq[0] = rhoBar;
        momentsEq[1] = (T)3 * jSqr - 2 * rhoBar;
        momentsEq[2] = -(T)3 * jSqr + rhoBar;
        momentsEq[3] = j[0];
        momentsEq[4] = -j[0];
        momentsEq[5] = j[1];
        momentsEq[6] = -j[1];
        momentsEq[7] = (j[0] * j[0] - j[1] * j[1]);
        momentsEq[8] = j[1] * j[0];
    }

    /// MRT collision step
    static T incMrtCollision(
        Array<T, Descriptor::q> &f, const T &rhoBar, const Array<T, 2> &j, const T &omega)
    {
        Array<T, 9> moments, momentsEq;

        computeMoments(moments, f);
        T jSqr = VectorTemplateImpl<T, 2>::normSqr(j);

        computeIncEquilibriumMoments(momentsEq, rhoBar, j, jSqr);
        computeMneqInPlace(moments, momentsEq);
        computef_InvM_Smoments(f, moments, omega);

        return jSqr;
    }

    /// MRT collision step
    static T incMrtCollision(Array<T, Descriptor::q> &f, const T &omega)
    {
        Array<T, 9> moments, momentsEq;

        computeMoments(moments, f);
        T rhoBar = moments[0];
        Array<T, 2> j(
            moments[MRTDescriptor::momentumIndexes[0]], moments[MRTDescriptor::momentumIndexes[1]]);
        T jSqr = VectorTemplateImpl<T, 2>::normSqr(j);

        computeIncEquilibriumMoments(momentsEq, rhoBar, j, jSqr);
        computeMneqInPlace(moments, momentsEq);
        computef_InvM_Smoments(f, moments, omega);

        return jSqr;
    }
};

}  // namespace plb

#endif
