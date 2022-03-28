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
 */

#ifndef DYNAMICS_TEMPLATES_2D_H
#define DYNAMICS_TEMPLATES_2D_H

#include "core/globalDefs.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"

namespace plb {

/// Compute Pi tensor efficiently on D2Q9 lattice
template <typename T>
struct neqPiD2Q9 {
    typedef SymmetricTensorImpl<T, 2> S;

    static T fromPiToFneq0(Array<T, 3> const &pi)
    {
        return (T)2 * (-(T)1 / (T)3 * pi[S::xx] - (T)1 / (T)3 * pi[S::yy]);
    }

    static T fromPiToFneq1(Array<T, 3> const &pi)
    {
        return (T)1 / (T)4 * ((T)1 / (T)3 * pi[S::xx] + (T)1 / (T)3 * pi[S::yy] - pi[S::xy]);
    }

    static T fromPiToFneq2(Array<T, 3> const &pi)
    {
        return (T)1 / (T)2 * ((T)2 / (T)3 * pi[S::xx] - (T)1 / (T)3 * pi[S::yy]);
    }

    static T fromPiToFneq3(Array<T, 3> const &pi)
    {
        return (T)1 / (T)4 * ((T)1 / (T)3 * pi[S::xx] + (T)1 / (T)3 * pi[S::yy] + pi[S::xy]);
    }

    static T fromPiToFneq4(Array<T, 3> const &pi)
    {
        return (T)1 / (T)2 * (-(T)1 / (T)3 * pi[S::xx] + (T)2 / (T)3 * pi[S::yy]);
    }
};  // struct neqPiD2Q9

// Efficient specialization for D2Q9 base lattice
template <typename T>
struct dynamicsTemplatesImpl<T, descriptors::D2Q9DescriptorBase<T> > {
    typedef descriptors::D2Q9DescriptorBase<T> D;

    static T computePsiComplete(T omega)
    {
        T psi = T();
        if (omega > (T)1.6) {
            psi = 121.9507329271408 * omega * omega * omega * omega
                  - 896.7378885747224 * omega * omega * omega + 2463.451616583476 * omega * omega
                  - 2992.597971569829 * omega + 1356.080857677411;
        } else {
            psi = .2026712402141036 * omega * omega * omega * omega
                  + 0.8142762886362326e-1 * omega * omega * omega
                  - 1.816012019668248 * omega * omega + 2.203381432056469 * omega;
        }
        if (psi > (T)2)
            psi = (T)2;
        return psi;
    }

    static T computePsiTruncated(T omega)
    {
        T psi = T();
        if (omega > (T)1.6) {
            psi = 72.39513224423209 * omega * omega * omega * omega
                  - 531.2772094802919 * omega * omega * omega + 1457.098269731003 * omega * omega
                  - 1766.422775230620 * omega + 798.3480312798796;
        } else {
            psi = .2328939082036888 * omega * omega * omega * omega
                  + 0.4056579590482348e-1 * omega * omega * omega
                  - 1.820935751061464 * omega * omega + 2.218496607915061 * omega;
        }
        if (psi > (T)2)
            psi = (T)2;
        return psi;
    }

    static T bgk_ma2_equilibrium(plint iPop, T rhoBar, T invRho, Array<T, 2> const &j, T jSqr)
    {
        typedef descriptors::D2Q9DescriptorBase<T> L;
        T c_j = L::c[iPop][0] * j[0] + L::c[iPop][1] * j[1];
        return L::t[iPop] * (rhoBar + (T)3 * c_j + invRho * (4.5 * c_j * c_j - 1.5 * jSqr));
    }

    static T complete_bgk_ma2_equilibrium(
        plint iPop, T rhoBar, T invRho, Array<T, 2> const &j, T jSqr)
    {
        typedef descriptors::D2Q9DescriptorBase<T> L;
        T c_j = L::c[iPop][0] * j[0] + L::c[iPop][1] * j[1];
        T ux = j[0] * invRho;
        T uy = j[1] * invRho;
        T ux2 = ux * ux;
        T uy2 = uy * uy;

        T v6 = (T)L::c[iPop][0] * (util::sqr((T)L::c[iPop][1]) - D::cs2);
        T v7 = (T)L::c[iPop][1] * (util::sqr((T)L::c[iPop][0]) - D::cs2);
        T v8 = (util::sqr((T)L::c[iPop][0]) - D::cs2) * (util::sqr((T)L::c[iPop][1]) - D::cs2);

        return L::t[iPop]
               * (rhoBar + (T)3 * c_j + invRho * (4.5 * c_j * c_j - 1.5 * jSqr)
                  + (T)27 / (T)2 * (v6 * j[0] * uy2 + v7 * j[1] * ux2)
                  + (T)81 / (T)4 * v8 * j[0] * ux * uy2);
    }

    static void bgk_ma2_equilibria(
        T rhoBar, T invRho, Array<T, D::d> const &j, T jSqr, Array<T, D::q> &eqPop)
    {
        T t0 = D::t[0];
        T t1 = D::t[1];
        T t2 = D::t[2];

        T kx = (T)3 * j[0];
        T ky = (T)3 * j[1];
        T kxSqr_ = invRho / (T)2 * kx * kx;
        T kySqr_ = invRho / (T)2 * ky * ky;
        T kxky_ = invRho * kx * ky;

        T C1 = rhoBar + invRho * (T)3 * jSqr;
        T C2, C3;

        // i=0
        C3 = -kxSqr_ - kySqr_;
        eqPop[0] = t0 * (C1 + C3);

        // i=1 and i=5
        C2 = -kx + ky;
        C3 = -kxky_;
        eqPop[1] = t1 * (C1 + C2 + C3);
        eqPop[5] = t1 * (C1 - C2 + C3);

        // i=2 and i=6
        C2 = -kx;
        C3 = -kySqr_;
        eqPop[2] = t2 * (C1 + C2 + C3);
        eqPop[6] = t2 * (C1 - C2 + C3);

        // i=3 and i=7
        C2 = -kx - ky;
        C3 = kxky_;
        eqPop[3] = t1 * (C1 + C2 + C3);
        eqPop[7] = t1 * (C1 - C2 + C3);

        // i=4 and i=8
        C2 = -ky;
        C3 = -kxSqr_;
        eqPop[4] = t2 * (C1 + C2 + C3);
        eqPop[8] = t2 * (C1 - C2 + C3);
    }

    static void complete_bgk_ma2_equilibria(
        T rhoBar, T invRho, Array<T, D::d> const &j, T jSqr, Array<T, D::q> &eqPop)
    {
        // T t0 = D::t[0];
        // T t1 = D::t[1];
        // T t2 = D::t[2];

        // T kx     = (T)3 * j[0];
        // T ky     = (T)3 * j[1];
        // T kxSqr_ = invRho / (T)2 * kx*kx;
        // T kySqr_ = invRho / (T)2 * ky*ky;
        // T kxky_  = invRho * kx*ky;

        // T C1 = rhoBar + invRho*(T)3*jSqr;
        // T C2, C3;

        // T ux = j[0]*invRho; T uy = j[1]*invRho;
        // T ux2 = ux*ux; T uy2 = uy*uy;

        // // i=0
        // C3 = -kxSqr_ - kySqr_;
        // eqPop[0] = t0 * (C1+C3) + j[0]*ux*uy2;

        // // i=1 and i=5
        // C2 = -kx + ky;
        // C3 = -kxky_;
        // eqPop[1] = t1 * (C1+C2+C3) + (T)0.25*j[0]*uy*(ux*uy+ux-uy);
        // eqPop[5] = t1 * (C1-C2+C3) + (T)0.25*j[0]*uy*(ux*uy-ux+uy);

        // // i=2 and i=6
        // C2 = -kx;
        // C3 = -kySqr_;
        // eqPop[2] = t2 * (C1+C2+C3) - (T)0.5*j[0]*uy2*(ux-(T)1);
        // eqPop[6] = t2 * (C1-C2+C3) - (T)0.5*j[0]*uy2*(ux+(T)1);

        // // i=3 and i=7
        // C2 = -kx - ky;
        // C3 = kxky_;
        // eqPop[3] = t1 * (C1+C2+C3) + (T)0.25*j[0]*uy*(ux*uy-ux-uy);
        // eqPop[7] = t1 * (C1-C2+C3) + (T)0.25*j[0]*uy*(ux*uy+ux+uy);

        // // i=4 and i=8
        // C2 = -ky;
        // C3 = -kxSqr_;
        // eqPop[4] = t2 * (C1+C2+C3) - (T)0.5*j[1]*ux2*(uy-(T)1);
        // eqPop[8] = t2 * (C1-C2+C3) - (T)0.5*j[1]*ux2*(uy+(T)1);

        T t0 = (T)0.25 * D::t[0];
        T t1 = D::t[1];
        T t2 = (T)0.5 * D::t[2];

        T rho = D::fullRho(rhoBar);
        T ux = j[0] * invRho;
        T uy = j[1] * invRho;
        T ux2 = ux * ux;
        T uy2 = uy * uy;

        T Cx = (T)3 * ux2 - 2;
        T Cy = 3 * uy2 - 2;
        eqPop[0] = t0 * (rho * Cy * Cx) - D::t[0];

        T k1x = (T)3 * (ux2 - ux) + (T)1;
        T k2x = k1x + 6 * ux;
        T k1y = (T)3 * (uy2 - uy) + (T)1;
        T k2y = k1y + 6 * uy;

        eqPop[2] = -t2 * (rho * Cy * k1x) - D::t[2];
        eqPop[4] = -t2 * (rho * k1y * Cx) - D::t[2];
        eqPop[6] = -t2 * (rho * Cy * k2x) - D::t[2];
        eqPop[8] = -t2 * (rho * k2y * Cx) - D::t[2];

        eqPop[1] = t1 * (rho * k2y * k1x) - t1;
        eqPop[5] = t1 * (rho * k1y * k2x) - t1;
        eqPop[3] = t1 * (rho * k1y * k1x) - t1;
        eqPop[7] = t1 * (rho * k2y * k2x) - t1;
    }

    static void complete_ma2_equilibrium_moments(
        T rhoBar, T invRho, Array<T, D::d> const &j, Array<T, D::q> &momEq, plint order)
    {
        momEq[0] = rhoBar;
        momEq[1] = j[0];
        momEq[2] = j[1];
        if (order >= 2) {
            momEq[3] = j[0] * j[0] * invRho;
            momEq[4] = j[0] * j[1] * invRho;
            momEq[5] = j[1] * j[1] * invRho;
            if (order >= 3) {
                momEq[6] = j[0] * momEq[5] * invRho;
                momEq[7] = momEq[3] * j[1] * invRho;
                if (order >= 4) {
                    momEq[8] = momEq[3] * momEq[5] * invRho;
                }
            }
        } else {
            PLB_ASSERT(order >= 2 && "Order must be greater than 2.");
        }
    }

    static void truncated_ma2_equilibrium_moments(
        T rhoBar, T invRho, Array<T, D::d> const &j, Array<T, D::q> &momEq)
    {
        momEq[0] = rhoBar;
        momEq[1] = j[0];
        momEq[2] = j[1];
        momEq[3] = j[0] * j[0] * invRho;
        momEq[4] = j[0] * j[1] * invRho;
        momEq[5] = j[1] * j[1] * invRho;
        momEq[6] = T();
        momEq[7] = T();
        momEq[8] = T();
    }

    static void complete_ma2_populations(Array<T, D::q> &a, Array<T, D::q> &f)
    {
        f[0] = a[0] - (T)1.5 * a[3] - (T)1.5 * a[5] + (T)2.25 * a[8];
        f[1] = a[0] - 3 * a[1] + 3 * a[2] + 3 * a[3] - 9 * a[4] + 3 * a[5] - 9 * a[6] + 9 * a[7]
               + 9 * a[8];
        f[2] = a[0] - 3 * a[1] + 3 * a[3] - (T)1.5 * a[5] + (T)4.5 * a[6] - (T)4.5 * a[8];

        f[3] = a[0] - 3 * a[1] - 3 * a[2] + 3 * a[3] + 9 * a[4] + 3 * a[5] - 9 * a[6] - 9 * a[7]
               + 9 * a[8];
        f[4] = a[0] - 3 * a[2] - (T)1.5 * a[3] + 3 * a[5] + (T)4.5 * a[7] - (T)4.5 * a[8];
        f[5] = a[0] + 3 * a[1] - 3 * a[2] + 3 * a[3] - 9 * a[4] + 3 * a[5] + 9 * a[6] - 9 * a[7]
               + 9 * a[8];
        f[6] = a[0] + 3 * a[1] + 3 * a[3] - (T)1.5 * a[5] - (T)4.5 * a[6] - (T)4.5 * a[8];
        f[7] = a[0] + 3 * a[1] + 3 * a[2] + 3 * a[3] + 9 * a[4] + 3 * a[5] + 9 * a[6] + 9 * a[7]
               + 9 * a[8];
        f[8] = a[0] + 3 * a[2] - (T)1.5 * a[3] + 3 * a[5] - (T)4.5 * a[7] - (T)4.5 * a[8];

        for (plint iPop = 0; iPop < D::q; ++iPop)
            f[iPop] *= D::t[iPop];
    }

    static void complete_ma2_moments(Array<T, D::q> &f, Array<T, D::q> &mom)
    {
        T f15 = f[1] + f[5];
        T f13 = f[1] + f[3];
        T f1_3 = f[1] - f[3];
        T f37 = f[3] + f[7];
        T f57 = f[5] + f[7];
        T f5_7 = f[5] - f[7];
        T f26 = f[2] + f[6];
        T f2_6 = f[2] - f[6];
        T f48 = f[4] + f[8];
        T f4_8 = f[4] - f[8];

        mom[0] = f[0] + f15 + f26 + f37 + f48;
        mom[1] = -f13 - f2_6 + f57;
        mom[2] = f1_3 - f4_8 - f5_7;

        mom[3] = D::cs2 * ((T)2 * (f15 + f37 + f26) - f[0] - f48);

        mom[4] = -f15 + f37;

        mom[5] = D::cs2 * ((T)2 * (f15 + f37 + f48) - f[0] - f26);

        mom[6] = D::cs2 * ((T)2 * (-f13 + f57) + f2_6);

        mom[7] = D::cs2 * ((T)2 * (f1_3 - f5_7) + f4_8);

        mom[8] = D::cs2 * D::cs2 * (f[0] - (T)2 * (f26 + f48) + (T)4 * (f13 + f57));
    }

    static void computeInvMmNeq(const Array<T, D::q> &mNeq, Array<T, D::q> &f)
    {
        f[0] = D::t[0] * (mNeq[0] - (T)1.5 * (mNeq[3] + mNeq[5]) + (T)2.25 * mNeq[8]);

        T sym = mNeq[0] + (T)3 * mNeq[3] - (T)9 * mNeq[4] + (T)3 * mNeq[5] + (T)9 * mNeq[8];
        T asym = -(T)3 * mNeq[1] + (T)3 * mNeq[2] - (T)9 * mNeq[6] + (T)9 * mNeq[7];
        f[1] = D::t[1] * (sym + asym);
        f[5] = D::t[1] * (sym - asym);
        sym = mNeq[0] + (T)3 * mNeq[3] + (T)9 * mNeq[4] + (T)3 * mNeq[5] + (T)9 * mNeq[8];
        asym = -(T)3 * mNeq[1] - (T)3 * mNeq[2] - (T)9 * mNeq[6] - (T)9 * mNeq[7];
        f[3] = D::t[1] * (sym + asym);
        f[7] = D::t[1] * (sym - asym);

        sym = mNeq[0] + (T)3 * mNeq[3] - (T)1.5 * mNeq[5] - (T)4.5 * mNeq[8];
        asym = -(T)3 * mNeq[1] + (T)4.5 * mNeq[6];
        f[2] = D::t[2] * (sym + asym);
        f[6] = D::t[2] * (sym - asym);

        sym = mNeq[0] - (T)1.5 * mNeq[3] + (T)3 * mNeq[5] - (T)4.5 * mNeq[8];
        asym = -(T)3 * mNeq[2] + (T)4.5 * mNeq[7];
        f[4] = D::t[2] * (sym + asym);
        f[8] = D::t[2] * (sym - asym);
    }

    static void complete_neq_ma2_moments_from_phys_moments(
        Array<T, D::q> &mNeq, T rhoBar, T invRho, const Array<T, D::d> &j,
        const Array<T, SymmetricTensorImpl<T, D::d>::n> &piNeq, plint order, T omega,
        T omegaNonPhys)
    {
        T omegaRatio = omega / omegaNonPhys;
        mNeq[0] = T();

        mNeq[1] = T();
        mNeq[2] = T();

        if (order >= 2) {
            mNeq[3] = piNeq[0];
            mNeq[4] = piNeq[1];
            mNeq[5] = piNeq[2];
            if (order >= 3) {
                mNeq[6] = invRho * omegaRatio * (piNeq[2] * j[0] + (T)2 * piNeq[1] * j[1]);
                mNeq[7] = invRho * omegaRatio * (piNeq[0] * j[1] + (T)2 * piNeq[1] * j[0]);
                if (order >= 4) {
                    mNeq[8] = invRho * invRho * omegaRatio
                              * (j[0] * j[0] * piNeq[2] + (T)4 * j[0] * j[1] * piNeq[1]
                                 + j[1] * j[1] * piNeq[0]);
                }
            }
        } else {
            PLB_ASSERT(order >= 2 && "Order must be greater than 2.");
        }
    }

    static T complete_mrt_ma2_collision_base(
        Array<T, D::q> &f, plint order, T omega, T omegaNonPhys, plint iPhys)
    {
        Array<T, D::q> mNeq, mEq;
        complete_ma2_moments(f, mNeq);
        T invRho = D::invRho(mNeq[0]);
        complete_ma2_equilibrium_moments(
            mNeq[0], invRho, Array<T, D::d>(mNeq[1], mNeq[2]), mEq, order);

        T jSqr = mNeq[1] * mNeq[1] + mNeq[2] * mNeq[2];

        for (plint iPop = 1 + D::d; iPop <= iPhys; ++iPop) {
            mNeq[iPop] -= mEq[iPop];
            mNeq[iPop] *= -omega;
        }

        for (plint iPop = iPhys + 1; iPop < D::q; ++iPop) {
            mNeq[iPop] -= mEq[iPop];
            mNeq[iPop] *= -omegaNonPhys;
        }

        T m35 = mNeq[3] + mNeq[5];
        T m46 = mNeq[4] + mNeq[6];
        T m4_6 = mNeq[4] - mNeq[6];
        T m68 = mNeq[6] + mNeq[8];
        T m6_8 = mNeq[6] - mNeq[8];
        T m78 = mNeq[7] + mNeq[8];
        T m7_8 = mNeq[7] - mNeq[8];

        f[0] += -(T)2 * D::cs2 * m35 + mNeq[8];

        f[1] += (T)0.25 * (D::cs2 * m35 - m46 + m78);
        f[3] += (T)0.25 * (D::cs2 * m35 + m4_6 - m7_8);
        f[5] += (T)0.25 * (D::cs2 * m35 - m4_6 - m7_8);
        f[7] += (T)0.25 * (D::cs2 * m35 + m46 + m78);

        f[2] += D::cs2 * (mNeq[3] - (T)0.5 * mNeq[5]) + (T)0.5 * m6_8;
        f[6] += D::cs2 * (mNeq[3] - (T)0.5 * mNeq[5]) - (T)0.5 * m68;

        f[4] += -D::cs2 * ((T)0.5 * mNeq[3] - mNeq[5]) + (T)0.5 * m7_8;
        f[8] += -D::cs2 * ((T)0.5 * mNeq[3] - mNeq[5]) - (T)0.5 * m78;

        return invRho * invRho * jSqr;
    }

    static T complete_regularized_mrt_ma2_collision_base(
        Array<T, D::q> &f, plint order, T rhoBar, const Array<T, D::d> &j,
        const Array<T, SymmetricTensorImpl<T, D::d>::n> &piNeq, T omega, T omegaNonPhys,
        plint iPhys)
    {
        Array<T, D::q> mNeq, mEq;
        // T rhoBar;
        // Array<T,D::d> j;
        // Array<T,SymmetricTensorImpl<T,D::d>::n> piNeq;
        // momentTemplatesImpl<T,D>::compute_rhoBar_j_PiNeq(f,rhoBar,j,piNeq);
        T jSqr = j[0] * j[0] + j[1] * j[1];
        T invRho = D::invRho(rhoBar);

        complete_neq_ma2_moments_from_phys_moments(
            mNeq, rhoBar, invRho, j, piNeq, order, omega, omegaNonPhys);

        complete_ma2_equilibrium_moments(rhoBar, invRho, j, mEq, order + 1);

        for (plint iPop = 0; iPop < 1 + D::d; ++iPop) {
            mNeq[iPop] = mEq[iPop];
        }

        T one_m_omega = ((T)1 - omega);
        for (plint iPop = 1 + D::d; iPop <= iPhys; ++iPop) {
            mNeq[iPop] *= one_m_omega;
            mNeq[iPop] += mEq[iPop];
        }

        one_m_omega = ((T)1 - omegaNonPhys);
        for (plint iPop = iPhys + 1; iPop < D::q; ++iPop) {
            mNeq[iPop] *= one_m_omega;
            mNeq[iPop] += mEq[iPop];
        }

        computeInvMmNeq(mNeq, f);

        return invRho * invRho * jSqr;
    }

    static T truncated_mrt_ma2_collision_base(
        Array<T, D::q> &f, T omega, T omegaNonPhys, plint iPhys)
    {
        Array<T, D::q> mNeq, mEq;
        complete_ma2_moments(f, mNeq);
        T invRho = D::invRho(mNeq[0]);
        truncated_ma2_equilibrium_moments(mNeq[0], invRho, Array<T, D::d>(mNeq[1], mNeq[2]), mEq);

        T jSqr = mNeq[1] * mNeq[1] + mNeq[2] * mNeq[2];

        for (plint iPop = 3; iPop <= iPhys; ++iPop) {
            mNeq[iPop] -= mEq[iPop];
            mNeq[iPop] *= -omega;
        }

        for (plint iPop = iPhys + 1; iPop < D::q; ++iPop) {
            mNeq[iPop] -= mEq[iPop];
            mNeq[iPop] *= -omegaNonPhys;
        }

        T m35 = mNeq[3] + mNeq[5];
        T m46 = mNeq[4] + mNeq[6];
        T m4_6 = mNeq[4] - mNeq[6];
        T m68 = mNeq[6] + mNeq[8];
        T m6_8 = mNeq[6] - mNeq[8];
        T m78 = mNeq[7] + mNeq[8];
        T m7_8 = mNeq[7] - mNeq[8];

        f[0] += -(T)2 * D::cs2 * m35 + mNeq[8];

        f[1] += (T)0.25 * (D::cs2 * m35 - m46 + m78);
        f[3] += (T)0.25 * (D::cs2 * m35 + m4_6 - m7_8);
        f[5] += (T)0.25 * (D::cs2 * m35 - m4_6 - m7_8);
        f[7] += (T)0.25 * (D::cs2 * m35 + m46 + m78);

        f[2] += D::cs2 * (mNeq[3] - (T)0.5 * mNeq[5]) + (T)0.5 * m6_8;
        f[6] += D::cs2 * (mNeq[3] - (T)0.5 * mNeq[5]) - (T)0.5 * m68;

        f[4] += -D::cs2 * ((T)0.5 * mNeq[3] - mNeq[5]) + (T)0.5 * m7_8;
        f[8] += -D::cs2 * ((T)0.5 * mNeq[3] - mNeq[5]) - (T)0.5 * m78;

        return invRho * invRho * jSqr;
    }

    static T complete_mrt_ma2_ext_rhoBar_j_collision_base(
        Array<T, D::q> &f, T rhoBar, Array<T, D::d> const &j, plint order, T omega, T omegaNonPhys,
        plint iPhys)
    {
        Array<T, D::q> mNeq, mEq;
        complete_ma2_moments(f, mNeq);
        T invRho = D::invRho(rhoBar);
        complete_ma2_equilibrium_moments(rhoBar, invRho, j, mEq, order + 1);

        T jSqr = j[0] * j[0] + j[1] * j[1];

        //     // max moment to relax with omega
        for (plint iPop = 0; iPop <= iPhys; ++iPop) {
            mNeq[iPop] -= mEq[iPop];
            mNeq[iPop] *= -omega;
        }

        for (plint iPop = iPhys + 1; iPop < D::q; ++iPop) {
            mNeq[iPop] -= mEq[iPop];
            mNeq[iPop] *= -omegaNonPhys;
        }

        T m35 = mNeq[3] + mNeq[5];
        T m46 = mNeq[4] + mNeq[6];
        T m4_6 = mNeq[4] - mNeq[6];
        T m68 = mNeq[6] + mNeq[8];
        T m6_8 = mNeq[6] - mNeq[8];
        T m78 = mNeq[7] + mNeq[8];
        T m7_8 = mNeq[7] - mNeq[8];

        f[0] += -(T)2 * D::cs2 * m35 + mNeq[8];

        f[1] += (T)0.25 * (D::cs2 * m35 - m46 + m78);
        f[3] += (T)0.25 * (D::cs2 * m35 + m4_6 - m7_8);
        f[5] += (T)0.25 * (D::cs2 * m35 - m4_6 - m7_8);
        f[7] += (T)0.25 * (D::cs2 * m35 + m46 + m78);

        f[2] += D::cs2 * (mNeq[3] - (T)0.5 * mNeq[5]) + (T)0.5 * m6_8;
        f[6] += D::cs2 * (mNeq[3] - (T)0.5 * mNeq[5]) - (T)0.5 * m68;

        f[4] += -D::cs2 * ((T)0.5 * mNeq[3] - mNeq[5]) + (T)0.5 * m7_8;
        f[8] += -D::cs2 * ((T)0.5 * mNeq[3] - mNeq[5]) - (T)0.5 * m78;

        mNeq[0] *= D::cs2 * D::cs2;
        f[0] += (T)4 * mNeq[0];

        mNeq[1] *= D::cs2;
        mNeq[2] *= D::cs2;

        f[2] += mNeq[0] - mNeq[1];
        f[6] += mNeq[0] + mNeq[1];
        f[4] += mNeq[0] - mNeq[2];
        f[8] += mNeq[0] + mNeq[2];

        mNeq[0] *= (T)0.25;
        mNeq[1] *= (T)0.25;
        mNeq[2] *= (T)0.25;

        f[1] += mNeq[0] - mNeq[1] + mNeq[2];
        f[5] += mNeq[0] + mNeq[1] - mNeq[2];
        f[3] += mNeq[0] - mNeq[1] - mNeq[2];
        f[7] += mNeq[0] + mNeq[1] + mNeq[2];

        return invRho * invRho * jSqr;
    }

    static T truncated_mrt_ma2_ext_rhoBar_j_collision_base(
        Array<T, D::q> &f, T rhoBar, Array<T, D::d> const &j, T omega, T omegaNonPhys, plint iPhys)
    {
        Array<T, D::q> mNeq, mEq;
        complete_ma2_moments(f, mNeq);
        T invRho = D::invRho(rhoBar);
        truncated_ma2_equilibrium_moments(rhoBar, invRho, j, mEq);

        T jSqr = j[0] * j[0] + j[1] * j[1];

        //     // max moment to relax with omega
        for (plint iPop = 0; iPop <= iPhys; ++iPop) {
            mNeq[iPop] -= mEq[iPop];
            mNeq[iPop] *= -omega;
        }

        for (plint iPop = iPhys + 1; iPop < D::q; ++iPop) {
            mNeq[iPop] -= mEq[iPop];
            mNeq[iPop] *= -omegaNonPhys;
        }

        T m35 = mNeq[3] + mNeq[5];
        T m46 = mNeq[4] + mNeq[6];
        T m4_6 = mNeq[4] - mNeq[6];
        T m68 = mNeq[6] + mNeq[8];
        T m6_8 = mNeq[6] - mNeq[8];
        T m78 = mNeq[7] + mNeq[8];
        T m7_8 = mNeq[7] - mNeq[8];

        f[0] += -(T)2 * D::cs2 * m35 + mNeq[8];

        f[1] += (T)0.25 * (D::cs2 * m35 - m46 + m78);
        f[3] += (T)0.25 * (D::cs2 * m35 + m4_6 - m7_8);
        f[5] += (T)0.25 * (D::cs2 * m35 - m4_6 - m7_8);
        f[7] += (T)0.25 * (D::cs2 * m35 + m46 + m78);

        f[2] += D::cs2 * (mNeq[3] - (T)0.5 * mNeq[5]) + (T)0.5 * m6_8;
        f[6] += D::cs2 * (mNeq[3] - (T)0.5 * mNeq[5]) - (T)0.5 * m68;

        f[4] += -D::cs2 * ((T)0.5 * mNeq[3] - mNeq[5]) + (T)0.5 * m7_8;
        f[8] += -D::cs2 * ((T)0.5 * mNeq[3] - mNeq[5]) - (T)0.5 * m78;

        mNeq[0] *= D::cs2 * D::cs2;
        f[0] += (T)4 * mNeq[0];

        mNeq[1] *= D::cs2;
        mNeq[2] *= D::cs2;

        f[2] += mNeq[0] - mNeq[1];
        f[6] += mNeq[0] + mNeq[1];
        f[4] += mNeq[0] - mNeq[2];
        f[8] += mNeq[0] + mNeq[2];

        mNeq[0] *= (T)0.25;
        mNeq[1] *= (T)0.25;
        mNeq[2] *= (T)0.25;

        f[1] += mNeq[0] - mNeq[1] + mNeq[2];
        f[5] += mNeq[0] + mNeq[1] - mNeq[2];
        f[3] += mNeq[0] - mNeq[1] - mNeq[2];
        f[7] += mNeq[0] + mNeq[1] + mNeq[2];

        return invRho * invRho * jSqr;
    }

    static void complete_bgk_ma2_regularize(
        Array<T, D::q> &f, T rhoBar, T invRho, Array<T, D::d> const &j, T jSqr,
        Array<T, SymmetricTensorImpl<T, D::d>::n> const &piNeq, T omega, T omegaNonPhys,
        plint iPhys)
    {
        complete_bgk_ma2_equilibria(rhoBar, invRho, j, jSqr, f);
        T omegaRatio = omega / omegaNonPhys;

        Array<T, D::q> aOne;

        aOne.resetToZero();

        aOne[6] = invRho * (piNeq[2] * j[0] + (T)2 * piNeq[1] * j[1]);
        aOne[7] = invRho * (piNeq[0] * j[1] + (T)2 * piNeq[1] * j[0]);
        // aOne[8] = invRho * invRho *
        // (j[0]*j[0]*piNeq[2]+(T)4*j[0]*j[1]*piNeq[1]+j[1]*j[1]*piNeq[0]);
        aOne[8] = invRho * (aOne[6] * j[0] + aOne[7] * j[1]);

        for (plint iPop = iPhys + 1; iPop < D::q; ++iPop)
            aOne[iPop] *= omegaRatio;

        f[0] += -((T)2 / (T)3) * (piNeq[0] + piNeq[2]) + aOne[8];

        T sym = (T)0.25 * (D::cs2 * (piNeq[0] + piNeq[2]) - piNeq[1] + aOne[8]);
        T asym = (T)0.25 * (aOne[6] - aOne[7]);
        f[1] += sym - asym;
        f[5] += sym + asym;

        sym = D::cs2 * (piNeq[0] - (T)0.5 * piNeq[2]) - (T)0.5 * aOne[8];
        asym = (T)0.5 * aOne[6];
        f[2] += sym + asym;
        f[6] += sym - asym;

        sym = (T)0.25 * (D::cs2 * (piNeq[0] + piNeq[2]) + piNeq[1] + aOne[8]);
        asym = (T)0.25 * (aOne[6] + aOne[7]);

        f[3] += sym - asym;
        f[7] += sym + asym;

        sym = D::cs2 * (-(T)0.5 * piNeq[0] + piNeq[2]) - (T)0.5 * aOne[8];
        asym = (T)0.5 * aOne[7];
        f[4] += sym + asym;
        f[8] += sym - asym;
    }

    static T bgk_ma2_collision_base(
        Array<T, D::q> &f, T rhoBar, Array<T, 2> const &j, T omega, T invRho)
    {
        T one_m_omega = (T)1 - omega;
        T t0_omega = D::t[0] * omega;
        T t1_omega = D::t[1] * omega;
        T t2_omega = D::t[2] * omega;

        T jSqr = j[0] * j[0] + j[1] * j[1];
        T kx = (T)3 * j[0];
        T ky = (T)3 * j[1];
        T kxSqr_ = invRho / (T)2 * kx * kx;
        T kySqr_ = invRho / (T)2 * ky * ky;
        T kxky_ = invRho * kx * ky;

        T C1 = rhoBar + invRho * (T)3 * jSqr;
        T C2, C3;

        // i=0
        C3 = -kxSqr_ - kySqr_;
        f[0] *= one_m_omega;
        f[0] += t0_omega * (C1 + C3);

        // i=1 and i=5
        C2 = -kx + ky;
        C3 = -kxky_;
        f[1] *= one_m_omega;
        f[1] += t1_omega * (C1 + C2 + C3);
        f[5] *= one_m_omega;
        f[5] += t1_omega * (C1 - C2 + C3);

        // i=2 and i=6
        C2 = -kx;
        C3 = -kySqr_;
        f[2] *= one_m_omega;
        f[2] += t2_omega * (C1 + C2 + C3);
        f[6] *= one_m_omega;
        f[6] += t2_omega * (C1 - C2 + C3);

        // i=3 and i=7
        C2 = -kx - ky;
        C3 = kxky_;
        f[3] *= one_m_omega;
        f[3] += t1_omega * (C1 + C2 + C3);
        f[7] *= one_m_omega;
        f[7] += t1_omega * (C1 - C2 + C3);

        // i=4 and i=8
        C2 = -ky;
        C3 = -kxSqr_;
        f[4] *= one_m_omega;
        f[4] += t2_omega * (C1 + C2 + C3);
        f[8] *= one_m_omega;
        f[8] += t2_omega * (C1 - C2 + C3);

        return invRho * invRho * jSqr;
    }

    static T complete_bgk_ma2_collision_base(
        Array<T, D::q> &f, T rhoBar, T invRho, Array<T, 2> const &j, T omega)
    {
        T one_m_omega = (T)1 - omega;
        T t0_omega = D::t[0] * omega;
        T t1_omega = D::t[1] * omega;
        T t2_omega = D::t[2] * omega;

        T jSqr = j[0] * j[0] + j[1] * j[1];
        T kx = (T)3 * j[0];
        T ky = (T)3 * j[1];
        T kxSqr_ = invRho / (T)2 * kx * kx;
        T kySqr_ = invRho / (T)2 * ky * ky;
        T kxky_ = invRho * kx * ky;

        T C1 = rhoBar + invRho * (T)3 * jSqr;
        T C2, C3;

        T ux = j[0] * invRho;
        T uy = j[1] * invRho;

        T ux2 = ux * ux;
        T uy2 = uy * uy;

        // i=0
        C3 = -kxSqr_ - kySqr_;
        f[0] *= one_m_omega;
        f[0] += t0_omega * (C1 + C3) + omega * j[0] * ux * uy2;
        //     f[0] *= one_m_omega; f[0] += t0_omega * (C1+C3+j[0]*ux*uy2);

        // i=1 and i=5
        C2 = -kx + ky;
        C3 = -kxky_;
        f[1] *= one_m_omega;
        f[1] += t1_omega * (C1 + C2 + C3) + omega * (T)0.25 * j[0] * uy * (ux * uy + ux - uy);
        f[5] *= one_m_omega;
        f[5] += t1_omega * (C1 - C2 + C3) + omega * (T)0.25 * j[0] * uy * (ux * uy - ux + uy);
        //     f[1] *= one_m_omega; f[1] += t1_omega * (C1+C2+C3+(T)0.25*j[0]*uy*(ux*uy+ux-uy));
        //     f[5] *= one_m_omega; f[5] += t1_omega * (C1-C2+C3+(T)0.25*j[0]*uy*(ux*uy-ux+uy));

        // i=2 and i=6
        C2 = -kx;
        C3 = -kySqr_;
        f[2] *= one_m_omega;
        f[2] += t2_omega * (C1 + C2 + C3) - omega * (T)0.5 * j[0] * uy2 * (ux - (T)1);
        f[6] *= one_m_omega;
        f[6] += t2_omega * (C1 - C2 + C3) - omega * (T)0.5 * j[0] * uy2 * (ux + (T)1);
        //     f[2] *= one_m_omega; f[2] += t2_omega * (C1+C2+C3 - (T)0.5*j[0]*uy2*(ux-(T)1));
        //     f[6] *= one_m_omega; f[6] += t2_omega * (C1-C2+C3 - (T)0.5*j[0]*uy2*(ux+(T)1));

        // i=3 and i=7
        C2 = -kx - ky;
        C3 = kxky_;
        f[3] *= one_m_omega;
        f[3] += t1_omega * (C1 + C2 + C3) + omega * (T)0.25 * j[0] * uy * (ux * uy - ux - uy);
        f[7] *= one_m_omega;
        f[7] += t1_omega * (C1 - C2 + C3) + omega * (T)0.25 * j[0] * uy * (ux * uy + ux + uy);
        //     f[3] *= one_m_omega; f[3] += t1_omega * (C1+C2+C3 + (T)0.25*j[0]*uy*(ux*uy-ux-uy));
        //     f[7] *= one_m_omega; f[7] += t1_omega * (C1-C2+C3 + (T)0.25*j[0]*uy*(ux*uy+ux+uy));

        // i=4 and i=8
        C2 = -ky;
        C3 = -kxSqr_;
        f[4] *= one_m_omega;
        f[4] += t2_omega * (C1 + C2 + C3) - omega * (T)0.5 * j[1] * ux2 * (uy - (T)1);
        f[8] *= one_m_omega;
        f[8] += t2_omega * (C1 - C2 + C3) - omega * (T)0.5 * j[1] * ux2 * (uy + (T)1);
        //     f[4] *= one_m_omega; f[4] += t2_omega * (C1+C2+C3 - (T)0.5*j[1]*ux2*(uy-(T)1));
        //     f[8] *= one_m_omega; f[8] += t2_omega * (C1-C2+C3 - (T)0.5*j[1]*ux2*(uy+(T)1));

        return invRho * invRho * jSqr;
    }

    static T complete_regularized_bgk_ma2_collision_base(
        Array<T, D::q> &f, T rhoBar, T invRho, Array<T, 2> const &j, Array<T, 3> const &piNeq,
        T omega)
    {
        T rho = D::fullRho(rhoBar);
        T t0 = (T)0.25 * D::t[0] * rho;
        T t1 = D::t[1] * rho;
        T t2 = (T)0.5 * D::t[2] * rho;
        T ux = j[0] * invRho;
        T uy = j[1] * invRho;
        T ux2 = ux * ux;
        T uy2 = uy * uy;
        T Cx = (T)3 * ux2 - 2;
        T Cy = 3 * uy2 - 2;
        T k1x = (T)3 * (ux2 - ux) + (T)1;
        T k2x = k1x + 6 * ux;
        T k1y = (T)3 * (uy2 - uy) + (T)1;
        T k2y = k1y + 6 * uy;

        T one_omega = (T)1 - omega;

        f[0] = t0 * Cy * Cx - D::t[0];
        // f[0] += one_omega*(D::cs2*(Cy*piNeq[0]+Cx*piNeq[2])+(T)4*ux*uy*piNeq[1]);

        f[2] = -t2 * Cy * k1x - D::t[2];
        // f[2] += one_omega*(-2*uy*(ux-0.5)*piNeq[1]-0.5*D::cs2*(Cy*piNeq[0]+k1x*piNeq[2]));
        f[6] = -t2 * Cy * k2x - D::t[2];
        // f[6] += one_omega*(-2*uy*(ux+0.5)*piNeq[1]-0.5*D::cs2*(Cy*piNeq[0]+k2x*piNeq[2]));
        f[4] = -t2 * k1y * Cx - D::t[2];
        // f[4] += one_omega*(-2*ux*(uy-0.5)*piNeq[1]-0.5*D::cs2*(k1y*piNeq[0]+Cx*piNeq[2]));
        f[8] = -t2 * k2y * Cx - D::t[2];
        // f[8] += one_omega*(-2*ux*(uy+0.5)*piNeq[1]-0.5*D::cs2*(k2y*piNeq[0]+Cx*piNeq[2]));

        f[1] = t1 * k2y * k1x - D::t[1];
        // f[1] += one_omega*((ux-0.5)*(uy+0.5)*piNeq[1]+0.25*D::cs2*(k2y*piNeq[0]+k1x*piNeq[2]));
        f[5] = t1 * k1y * k2x - D::t[1];
        // f[5] += one_omega*((ux+0.5)*(uy-0.5)*piNeq[1]+0.25*D::cs2*(k1y*piNeq[0]+k2x*piNeq[2]));
        f[3] = t1 * k1y * k1x - D::t[1];
        // f[3] += one_omega*((ux-0.5)*(uy-0.5)*piNeq[1]+0.25*D::cs2*(k1y*piNeq[0]+k1x*piNeq[2]));
        f[7] = t1 * k2y * k2x - D::t[1];
        // f[7] += one_omega*((ux+0.5)*(uy+0.5)*piNeq[1]+0.25*D::cs2*(k2y*piNeq[0]+k2x*piNeq[2]));

        T aOne6 = piNeq[2] * ux + (T)2 * piNeq[1] * uy;
        T aOne7 = piNeq[0] * uy + (T)2 * piNeq[1] * ux;
        T aOne8 = aOne6 * ux + aOne7 * uy;
        // T aOne8 = T();
        f[0] += one_omega * (-((T)2 * D::cs2) * (piNeq[0] + piNeq[2]) + aOne8);

        T sym = (T)0.25 * (D::cs2 * (piNeq[0] + piNeq[2]) - piNeq[1] + aOne8);
        T asym = (T)0.25 * (aOne6 - aOne7);
        f[1] += one_omega * (sym - asym);
        f[5] += one_omega * (sym + asym);

        sym = D::cs2 * (piNeq[0] - (T)0.5 * piNeq[2]) - (T)0.5 * aOne8;
        asym = (T)0.5 * aOne6;
        f[2] += one_omega * (sym + asym);
        f[6] += one_omega * (sym - asym);
        sym = (T)0.25 * (D::cs2 * (piNeq[0] + piNeq[2]) + piNeq[1] + aOne8);
        asym = (T)0.25 * (aOne6 + aOne7);
        f[3] += one_omega * (sym - asym);
        f[7] += one_omega * (sym + asym);

        sym = D::cs2 * (-(T)0.5 * piNeq[0] + piNeq[2]) - (T)0.5 * aOne8;
        asym = (T)0.5 * aOne7;
        f[4] += one_omega * (sym + asym);
        f[8] += one_omega * (sym - asym);

        return ux2 + uy2;
    }

    static T complete_bgk_ma2_collision(
        Array<T, D::q> &f, T rhoBar, T invRho, Array<T, 2> const &j, T omega)
    {
        return complete_bgk_ma2_collision_base(f, rhoBar, invRho, j, omega);
    }

    static T complete_regularized_bgk_ma2_collision(
        Array<T, D::q> &f, T rhoBar, Array<T, 2> const &j, Array<T, 3> const &piNeq, T omega)
    {
        return complete_regularized_bgk_ma2_collision_base(
            f, rhoBar, D::invRho(rhoBar), j, piNeq, omega);
    }

    static T complete_mrt_ma2_collision(
        Array<T, D::q> &f, plint order, T omega, T omegaNonPhys, plint iPhys)
    {
        return complete_mrt_ma2_collision_base(f, order, omega, omegaNonPhys, iPhys);
    }

    static T complete_regularized_mrt_ma2_collision(
        Array<T, D::q> &f, T rhoBar, const Array<T, D::d> &j,
        const Array<T, SymmetricTensorImpl<T, D::d>::n> &piNeq, plint order, T omega,
        T omegaNonPhys, plint iPhys)
    {
        return complete_regularized_mrt_ma2_collision_base(
            f, order, rhoBar, j, piNeq, omega, omegaNonPhys, iPhys);
    }

    static T truncated_mrt_ma2_collision(Array<T, D::q> &f, T omega, T omegaNonPhys, plint iPhys)
    {
        return truncated_mrt_ma2_collision_base(f, omega, omegaNonPhys, iPhys);
    }

    static T complete_mrt_ma2_ext_rhoBar_j_collision(
        Array<T, D::q> &f, T rhoBar, Array<T, D::d> const &j, plint order, T omega, T omegaNonPhys,
        plint iPhys)
    {
        return complete_mrt_ma2_ext_rhoBar_j_collision_base(
            f, rhoBar, j, order, omega, omegaNonPhys, iPhys);
    }

    static T truncated_mrt_ma2_ext_rhoBar_j_collision(
        Array<T, D::q> &f, T rhoBar, Array<T, D::d> const &j, T omega, T omegaNonPhys, plint iPhys)
    {
        return truncated_mrt_ma2_ext_rhoBar_j_collision_base(
            f, rhoBar, j, omega, omegaNonPhys, iPhys);
    }

    static T bgk_ma2_collision(Array<T, D::q> &f, T rhoBar, Array<T, 2> const &j, T omega)
    {
        return bgk_ma2_collision_base(f, rhoBar, j, omega, D::invRho(rhoBar));
    }

    static T bgk_inc_collision(Array<T, D::q> &f, T rhoBar, Array<T, 2> const &j, T omega)
    {
        // Incompressible: rho0=1
        static const T invRho0 = (T)1;
        return bgk_ma2_collision_base(f, rhoBar, j, omega, invRho0);
    }

    static T rlb_collision(
        Array<T, D::q> &f, T rhoBar, T invRho, Array<T, 2> const &j, Array<T, 3> const &PiNeq,
        T omega)
    {
        typedef dynamicsTemplatesImpl<T, descriptors::D2Q9DescriptorBase<T> > DH;
        const T jSqr = j[0] * j[0] + j[1] * j[1];

        T piNeq0 = neqPiD2Q9<T>::fromPiToFneq0(PiNeq);
        T piNeq1 = neqPiD2Q9<T>::fromPiToFneq1(PiNeq);
        T piNeq2 = neqPiD2Q9<T>::fromPiToFneq2(PiNeq);
        T piNeq3 = neqPiD2Q9<T>::fromPiToFneq3(PiNeq);
        T piNeq4 = neqPiD2Q9<T>::fromPiToFneq4(PiNeq);

        f[0] = DH::bgk_ma2_equilibrium(0, rhoBar, invRho, j, jSqr) + ((T)1 - omega) * piNeq0;
        f[1] = DH::bgk_ma2_equilibrium(1, rhoBar, invRho, j, jSqr) + ((T)1 - omega) * piNeq1;
        f[2] = DH::bgk_ma2_equilibrium(2, rhoBar, invRho, j, jSqr) + ((T)1 - omega) * piNeq2;
        f[3] = DH::bgk_ma2_equilibrium(3, rhoBar, invRho, j, jSqr) + ((T)1 - omega) * piNeq3;
        f[4] = DH::bgk_ma2_equilibrium(4, rhoBar, invRho, j, jSqr) + ((T)1 - omega) * piNeq4;
        f[5] = DH::bgk_ma2_equilibrium(5, rhoBar, invRho, j, jSqr) + ((T)1 - omega) * piNeq1;
        f[6] = DH::bgk_ma2_equilibrium(6, rhoBar, invRho, j, jSqr) + ((T)1 - omega) * piNeq2;
        f[7] = DH::bgk_ma2_equilibrium(7, rhoBar, invRho, j, jSqr) + ((T)1 - omega) * piNeq3;
        f[8] = DH::bgk_ma2_equilibrium(8, rhoBar, invRho, j, jSqr) + ((T)1 - omega) * piNeq4;
        return jSqr * invRho * invRho;
    }

    static T bgk_ma2_constRho_collision(
        Array<T, D::q> &f, T rhoBar, Array<T, 2> const &j, T ratioRho, T omega)
    {
        typedef descriptors::D2Q9DescriptorBase<T> L;
        T invRho = L::invRho(rhoBar);
        const T jSqr = j[0] * j[0] + j[1] * j[1];
        for (plint iPop = 0; iPop < L::q; ++iPop) {
            T feq = bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
            f[iPop] =
                ratioRho * feq + L::t[iPop] * (ratioRho - (T)1) + ((T)1 - omega) * (f[iPop] - feq);
        }
        return jSqr * invRho * invRho;
    }

    static T precond_bgk_ma2_equilibrium(
        plint iPop, T rhoBar, T invRho, Array<T, 2> const &j, T jSqr, T invGamma)
    {
        typedef descriptors::D2Q9DescriptorBase<T> L;
        T c_j = L::c[iPop][0] * j[0] + L::c[iPop][1] * j[1];
        return L::t[iPop]
               * (rhoBar + (T)3 * c_j + invGamma * invRho * (4.5 * c_j * c_j - 1.5 * jSqr));
    }

    static T precond_bgk_ma2_collision_base(
        Array<T, D::q> &f, T rhoBar, Array<T, 2> const &j, T omega, T invGamma, bool incompressible)
    {
        T invRho = incompressible ? (T)1 : D::invRho(rhoBar);
        T one_m_omega = (T)1 - omega;
        T t0_omega = D::t[0] * omega;
        T t1_omega = D::t[1] * omega;
        T t2_omega = D::t[2] * omega;

        T jSqr = j[0] * j[0] + j[1] * j[1];
        T kx = (T)3 * j[0];
        T ky = (T)3 * j[1];
        T kxSqr_ = invGamma * invRho / (T)2 * kx * kx;
        T kySqr_ = invGamma * invRho / (T)2 * ky * ky;
        T kxky_ = invGamma * invRho * kx * ky;

        T C1 = rhoBar + invGamma * invRho * (T)3 * jSqr;
        T C2, C3;

        // i=0
        C3 = -kxSqr_ - kySqr_;
        f[0] *= one_m_omega;
        f[0] += t0_omega * (C1 + C3);

        // i=1 and i=5
        C2 = -kx + ky;
        C3 = -kxky_;
        f[1] *= one_m_omega;
        f[1] += t1_omega * (C1 + C2 + C3);
        f[5] *= one_m_omega;
        f[5] += t1_omega * (C1 - C2 + C3);

        // i=2 and i=6
        C2 = -kx;
        C3 = -kySqr_;
        f[2] *= one_m_omega;
        f[2] += t2_omega * (C1 + C2 + C3);
        f[6] *= one_m_omega;
        f[6] += t2_omega * (C1 - C2 + C3);

        // i=3 and i=7
        C2 = -kx - ky;
        C3 = kxky_;
        f[3] *= one_m_omega;
        f[3] += t1_omega * (C1 + C2 + C3);
        f[7] *= one_m_omega;
        f[7] += t1_omega * (C1 - C2 + C3);

        // i=4 and i=8
        C2 = -ky;
        C3 = -kxSqr_;
        f[4] *= one_m_omega;
        f[4] += t2_omega * (C1 + C2 + C3);
        f[8] *= one_m_omega;
        f[8] += t2_omega * (C1 - C2 + C3);

        return invRho * invRho * jSqr;
    }

    static T precond_bgk_ma2_collision(
        Array<T, D::q> &f, T rhoBar, Array<T, 2> const &j, T omega, T invGamma)
    {
        return precond_bgk_ma2_collision_base(f, rhoBar, j, omega, invGamma, false);
    }

};  // struct dynamicsTemplatesImpl<D2Q9DescriptorBase>

}  // namespace plb

#endif  // DYNAMICS_TEMPLATES_2D_H
