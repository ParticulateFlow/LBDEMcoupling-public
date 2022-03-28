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
 * Helper functions for the implementation of Asinari's LW-AWC model.
 * This file is all about efficiency. The generic template code is
 * specialized for commonly used Lattices, so that a maximum performance
 * can be taken out of each case.
 */
#ifndef ASINARI_TEMPLATES_H
#define ASINARI_TEMPLATES_H

#include "core/cell.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/offEquilibriumTemplates.h"

namespace plb {

template <typename T, class Descriptor>
struct asinariTemplatesImpl;

/// This structure forwards the calls to the appropriate helper class
template <typename T, template <typename U> class Descriptor>
struct asinariTemplates {
    static T bgk_ma2_nonsymmetric_equilibrium(plint iPop, Array<T, Descriptor<T>::d> const &j)
    {
        return asinariTemplatesImpl<
            T, typename Descriptor<T>::BaseDescriptor>::bgk_ma2_nonsymmetric_equilibrium(iPop, j);
    }

    static T bgk_collision_stage1(
        Cell<T, Descriptor> &cell, T rhoBar, T invRho, Array<T, Descriptor<T>::d> const &j, T jSqr,
        T prefactor)
    {
        return asinariTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::
            bgk_collision_stage1(cell.getRawPopulations(), rhoBar, invRho, j, jSqr, prefactor);
    }

    static void bgk_collision_stage3(
        Cell<T, Descriptor> &cell, Array<T, Descriptor<T>::d> const &j, T prefactor)
    {
        asinariTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::bgk_collision_stage3(
            cell.getRawPopulations(), j, prefactor);
    }

};  // struct dynamicsTemplates

/// All helper functions are inside this structure
template <typename T, class Descriptor>
struct asinariTemplatesImpl {
    static T bgk_ma2_nonsymmetric_equilibrium(plint iPop, Array<T, Descriptor::d> const &j)
    {
        T c_j = Descriptor::c[iPop][0] * j[0];
        for (int iD = 1; iD < Descriptor::d; ++iD) {
            c_j += Descriptor::c[iPop][iD] * j[iD];
        }
        return Descriptor::t[iPop] * Descriptor::invCs2 * c_j;
    }

    static T bgk_collision_stage1(
        Array<T, Descriptor::q> &f, T rhoBar, T invRho, Array<T, Descriptor::d> const &j, T jSqr,
        T prefactor)
    {
        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            f[iPop] = dynamicsTemplatesImpl<T, Descriptor>::bgk_ma2_equilibrium(
                          iPop, rhoBar, invRho, j, jSqr)
                      - prefactor * bgk_ma2_nonsymmetric_equilibrium(iPop, j);
        }
        return jSqr * util::sqr(invRho);
    }

    static void bgk_collision_stage3(
        Array<T, Descriptor::q> &f, Array<T, Descriptor::d> const &j, T prefactor)
    {
        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            f[iPop] += prefactor * bgk_ma2_nonsymmetric_equilibrium(iPop, j);
        }
    }

};  // struct asinariTemplatesImpl

// Efficient specialization for D2Q9 base lattice
template <typename T>
struct asinariTemplatesImpl<T, descriptors::D2Q9DescriptorBase<T> > {
    typedef descriptors::D2Q9DescriptorBase<T> Descriptor;

    static T bgk_ma2_nonsymmetric_equilibrium(plint iPop, Array<T, 2> const &j)
    {
        T c_j = Descriptor::c[iPop][0] * j[0] + Descriptor::c[iPop][1] * j[1];
        return Descriptor::t[iPop] * 3. * c_j;
    }

    static T bgk_collision_stage1(
        Array<T, Descriptor::q> &f, T rhoBar, T invRho, Array<T, Descriptor::d> const &j, T jSqr,
        T prefactor)
    {
        T p = (T)1 - prefactor;

        T t0 = Descriptor::t[0];
        T t1 = Descriptor::t[1];
        T t2 = Descriptor::t[2];

        T kx = (T)3 * j[0];
        T ky = (T)3 * j[1];
        T kxSqr_ = invRho / (T)2 * kx * kx;
        T kySqr_ = invRho / (T)2 * ky * ky;
        T kxky_ = invRho * kx * ky;

        T C1 = rhoBar + invRho * (T)3 * jSqr;
        T C2, C3;

        // i=0
        C3 = -kxSqr_ - kySqr_;
        f[0] = t0 * (C1 + C3);

        // i=1 and i=5
        C2 = p * (-kx + ky);
        C3 = -kxky_;
        f[1] = t1 * (C1 + C2 + C3);
        f[5] = t1 * (C1 - C2 + C3);

        // i=2 and i=6
        C2 = p * (-kx);
        C3 = -kySqr_;
        f[2] = t2 * (C1 + C2 + C3);
        f[6] = t2 * (C1 - C2 + C3);

        // i=3 and i=7
        C2 = p * (-kx - ky);
        C3 = kxky_;
        f[3] = t1 * (C1 + C2 + C3);
        f[7] = t1 * (C1 - C2 + C3);

        // i=4 and i=8
        C2 = p * (-ky);
        C3 = -kxSqr_;
        f[4] = t2 * (C1 + C2 + C3);
        f[8] = t2 * (C1 - C2 + C3);

        return jSqr * util::sqr(invRho);
    }

    static void bgk_collision_stage3(
        Array<T, Descriptor::q> &f, Array<T, Descriptor::d> const &j, T prefactor)
    {
        T t1 = Descriptor::t[1];
        T t2 = Descriptor::t[2];

        T kx = (T)3 * j[0];
        T ky = (T)3 * j[1];

        T C2;

        // i=1 and i=5
        C2 = prefactor * (-kx + ky);
        f[1] += t1 * (C2);
        f[5] += t1 * (-C2);

        // i=2 and i=6
        C2 = prefactor * (-kx);
        f[2] += t2 * (C2);
        f[6] += t2 * (-C2);

        // i=3 and i=7
        C2 = prefactor * (-kx - ky);
        f[3] += t1 * (C2);
        f[7] += t1 * (-C2);

        // i=4 and i=8
        C2 = prefactor * (-ky);
        f[4] += t2 * (C2);
        f[8] += t2 * (-C2);
    }
};

// Efficient specialization for D3Q27 lattice
template <typename T>
struct asinariTemplatesImpl<T, descriptors::D3Q27DescriptorBase<T> > {
    typedef descriptors::D3Q27DescriptorBase<T> Descriptor;

    static T bgk_ma2_nonsymmetric_equilibrium(plint iPop, Array<T, 3> const &j)
    {
        T c_j = Descriptor::c[iPop][0] * j[0] + Descriptor::c[iPop][1] * j[1]
                + Descriptor::c[iPop][2] * j[2];
        return Descriptor::t[iPop] * 3. * c_j;
    }

    static T bgk_collision_stage1(
        Array<T, Descriptor::q> &f, T rhoBar, T invRho, Array<T, Descriptor::d> const &j, T jSqr,
        T prefactor)
    {
        T p = (T)1 - prefactor;

        T t0 = Descriptor::t[0];
        T t1 = Descriptor::t[1];
        T t4 = Descriptor::t[4];
        T t10 = Descriptor::t[10];

        T kx = (T)3 * j[0];
        T ky = (T)3 * j[1];
        T kz = (T)3 * j[2];
        T kxSqr_ = invRho / (T)2 * kx * kx;
        T kySqr_ = invRho / (T)2 * ky * ky;
        T kzSqr_ = invRho / (T)2 * kz * kz;
        T kxky_ = invRho * kx * ky;
        T kxkz_ = invRho * kx * kz;
        T kykz_ = invRho * ky * kz;

        T C1 = rhoBar + invRho * (T)3 * jSqr;
        T C2, C3;

        // i=0
        C3 = -kxSqr_ - kySqr_ - kzSqr_;
        f[0] = t0 * (C1 + C3);

        // i=1 and i=14
        C2 = p * (-kx);
        C3 = -kySqr_ - kzSqr_;
        f[1] = t1 * (C1 + C2 + C3);
        f[14] = t1 * (C1 - C2 + C3);

        // i=2 and i=15
        C2 = p * (-ky);
        C3 = -kxSqr_ - kzSqr_;
        f[2] = t1 * (C1 + C2 + C3);
        f[15] = t1 * (C1 - C2 + C3);

        // i=3 and i=16
        C2 = p * (-kz);
        C3 = -kxSqr_ - kySqr_;
        f[3] = t1 * (C1 + C2 + C3);
        f[16] = t1 * (C1 - C2 + C3);

        // i=4 and i=17
        C2 = p * (-kx - ky);
        C3 = kxky_ - kzSqr_;
        f[4] = t4 * (C1 + C2 + C3);
        f[17] = t4 * (C1 - C2 + C3);

        // i=5 and i=18
        C2 = p * (-kx + ky);
        C3 = -kxky_ - kzSqr_;
        f[5] = t4 * (C1 + C2 + C3);
        f[18] = t4 * (C1 - C2 + C3);

        // i=6 and i=19
        C2 = p * (-kx - kz);
        C3 = kxkz_ - kySqr_;
        f[6] = t4 * (C1 + C2 + C3);
        f[19] = t4 * (C1 - C2 + C3);

        // i=7 and i=20
        C2 = p * (-kx + kz);
        C3 = -kxkz_ - kySqr_;
        f[7] = t4 * (C1 + C2 + C3);
        f[20] = t4 * (C1 - C2 + C3);

        // i=8 and i=21
        C2 = p * (-ky - kz);
        C3 = kykz_ - kxSqr_;
        f[8] = t4 * (C1 + C2 + C3);
        f[21] = t4 * (C1 - C2 + C3);

        // i=9 and i=22
        C2 = p * (-ky + kz);
        C3 = -kykz_ - kxSqr_;
        f[9] = t4 * (C1 + C2 + C3);
        f[22] = t4 * (C1 - C2 + C3);

        // i=10 and i=23
        C2 = p * (-kx - ky - kz);
        C3 = kxky_ + kxkz_ + kykz_;
        f[10] = t10 * (C1 + C2 + C3);
        f[23] = t10 * (C1 - C2 + C3);

        // i=11 and i=24
        C2 = p * (-kx - ky + kz);
        C3 = kxky_ - kxkz_ - kykz_;
        f[11] = t10 * (C1 + C2 + C3);
        f[24] = t10 * (C1 - C2 + C3);

        // i=12 and i=25
        C2 = p * (-kx + ky - kz);
        C3 = -kxky_ + kxkz_ - kykz_;
        f[12] = t10 * (C1 + C2 + C3);
        f[25] = t10 * (C1 - C2 + C3);

        // i=13 and i=26
        C2 = p * (-kx + ky + kz);
        C3 = -kxky_ - kxkz_ + kykz_;
        f[13] = t10 * (C1 + C2 + C3);
        f[26] = t10 * (C1 - C2 + C3);

        return jSqr * util::sqr(invRho);
    }

    static void bgk_collision_stage3(
        Array<T, Descriptor::q> &f, Array<T, Descriptor::d> const &j, T prefactor)
    {
        T t1 = Descriptor::t[1];
        T t4 = Descriptor::t[4];
        T t10 = Descriptor::t[10];

        T kx = (T)3 * j[0];
        T ky = (T)3 * j[1];
        T kz = (T)3 * j[2];

        T C2;

        // i=1 and i=14
        C2 = prefactor * (-kx);
        f[1] += t1 * (C2);
        f[14] += t1 * (-C2);

        // i=2 and i=15
        C2 = prefactor * (-ky);
        f[2] += t1 * (C2);
        f[15] += t1 * (-C2);

        // i=3 and i=16
        C2 = prefactor * (-kz);
        f[3] += t1 * (C2);
        f[16] += t1 * (-C2);

        // i=4 and i=17
        C2 = prefactor * (-kx - ky);
        f[4] += t4 * (C2);
        f[17] += t4 * (-C2);

        // i=5 and i=18
        C2 = prefactor * (-kx + ky);
        f[5] += t4 * (C2);
        f[18] += t4 * (-C2);

        // i=6 and i=19
        C2 = prefactor * (-kx - kz);
        f[6] += t4 * (C2);
        f[19] += t4 * (-C2);

        // i=7 and i=20
        C2 = prefactor * (-kx + kz);
        f[7] += t4 * (C2);
        f[20] += t4 * (-C2);

        // i=8 and i=21
        C2 = prefactor * (-ky - kz);
        f[8] += t4 * (C2);
        f[21] += t4 * (-C2);

        // i=9 and i=22
        C2 = prefactor * (-ky + kz);
        f[9] += t4 * (C2);
        f[22] += t4 * (-C2);

        // i=10 and i=23
        C2 = prefactor * (-kx - ky - kz);
        f[10] += t10 * (C2);
        f[23] += t10 * (-C2);

        // i=11 and i=24
        C2 = prefactor * (-kx - ky + kz);
        f[11] += t10 * (C2);
        f[24] += t10 * (-C2);

        // i=12 and i=25
        C2 = prefactor * (-kx + ky - kz);
        f[12] += t10 * (C2);
        f[25] += t10 * (-C2);

        // i=13 and i=26
        C2 = prefactor * (-kx + ky + kz);
        f[13] += t10 * (C2);
        f[26] += t10 * (-C2);
    }
};

// Efficient specialization for D3Q19 lattice
template <typename T>
struct asinariTemplatesImpl<T, descriptors::D3Q19DescriptorBase<T> > {
    typedef descriptors::D3Q19DescriptorBase<T> Descriptor;

    static T bgk_ma2_nonsymmetric_equilibrium(plint iPop, Array<T, 3> const &j)
    {
        T c_j = Descriptor::c[iPop][0] * j[0] + Descriptor::c[iPop][1] * j[1]
                + Descriptor::c[iPop][2] * j[2];
        return Descriptor::t[iPop] * 3. * c_j;
    }

    static T bgk_collision_stage1(
        Array<T, Descriptor::q> &f, T rhoBar, T invRho, Array<T, Descriptor::d> const &j, T jSqr,
        T prefactor)
    {
        T p = (T)1 - prefactor;

        T t0 = Descriptor::t[0];
        T t1 = Descriptor::t[1];
        T t4 = Descriptor::t[4];

        T kx = (T)3 * j[0];
        T ky = (T)3 * j[1];
        T kz = (T)3 * j[2];
        T kxSqr_ = invRho / (T)2 * kx * kx;
        T kySqr_ = invRho / (T)2 * ky * ky;
        T kzSqr_ = invRho / (T)2 * kz * kz;
        T kxky_ = invRho * kx * ky;
        T kxkz_ = invRho * kx * kz;
        T kykz_ = invRho * ky * kz;

        T C1 = rhoBar + invRho * (T)3 * jSqr;
        T C2, C3;

        // i=0
        C3 = -kxSqr_ - kySqr_ - kzSqr_;
        f[0] = t0 * (C1 + C3);

        // i=1 and i=10
        C2 = p * (-kx);
        C3 = -kySqr_ - kzSqr_;
        f[1] = t1 * (C1 + C2 + C3);
        f[10] = t1 * (C1 - C2 + C3);

        // i=2 and i=11
        C2 = p * (-ky);
        C3 = -kxSqr_ - kzSqr_;
        f[2] = t1 * (C1 + C2 + C3);
        f[11] = t1 * (C1 - C2 + C3);

        // i=3 and i=12
        C2 = p * (-kz);
        C3 = -kxSqr_ - kySqr_;
        f[3] = t1 * (C1 + C2 + C3);
        f[12] = t1 * (C1 - C2 + C3);

        // i=4 and i=13
        C2 = p * (-kx - ky);
        C3 = kxky_ - kzSqr_;
        f[4] = t4 * (C1 + C2 + C3);
        f[13] = t4 * (C1 - C2 + C3);

        // i=5 and i=14
        C2 = p * (-kx + ky);
        C3 = -kxky_ - kzSqr_;
        f[5] = t4 * (C1 + C2 + C3);
        f[14] = t4 * (C1 - C2 + C3);

        // i=6 and i=15
        C2 = p * (-kx - kz);
        C3 = kxkz_ - kySqr_;
        f[6] = t4 * (C1 + C2 + C3);
        f[15] = t4 * (C1 - C2 + C3);

        // i=7 and i=16
        C2 = p * (-kx + kz);
        C3 = -kxkz_ - kySqr_;
        f[7] = t4 * (C1 + C2 + C3);
        f[16] = t4 * (C1 - C2 + C3);

        // i=8 and i=17
        C2 = p * (-ky - kz);
        C3 = kykz_ - kxSqr_;
        f[8] = t4 * (C1 + C2 + C3);
        f[17] = t4 * (C1 - C2 + C3);

        // i=9 and i=18
        C2 = p * (-ky + kz);
        C3 = -kykz_ - kxSqr_;
        f[9] = t4 * (C1 + C2 + C3);
        f[18] = t4 * (C1 - C2 + C3);

        return jSqr * util::sqr(invRho);
    }

    static void bgk_collision_stage3(
        Array<T, Descriptor::q> &f, Array<T, Descriptor::d> const &j, T prefactor)
    {
        T t1 = Descriptor::t[1];
        T t4 = Descriptor::t[4];

        T kx = (T)3 * j[0];
        T ky = (T)3 * j[1];
        T kz = (T)3 * j[2];

        T C2;

        // i=1 and i=10
        C2 = prefactor * (-kx);
        f[1] += t1 * (C2);
        f[10] += t1 * (-C2);

        // i=2 and i=11
        C2 = prefactor * (-ky);
        f[2] += t1 * (C2);
        f[11] += t1 * (-C2);

        // i=3 and i=12
        C2 = prefactor * (-kz);
        f[3] += t1 * (C2);
        f[12] += t1 * (-C2);

        // i=4 and i=13
        C2 = prefactor * (-kx - ky);
        f[4] += t4 * (C2);
        f[13] += t4 * (-C2);

        // i=5 and i=14
        C2 = prefactor * (-kx + ky);
        f[5] += t4 * (C2);
        f[14] += t4 * (-C2);

        // i=6 and i=15
        C2 = prefactor * (-kx - kz);
        f[6] += t4 * (C2);
        f[15] += t4 * (-C2);

        // i=7 and i=16
        C2 = prefactor * (-kx + kz);
        f[7] += t4 * (C2);
        f[16] += t4 * (-C2);

        // i=8 and i=17
        C2 = prefactor * (-ky - kz);
        f[8] += t4 * (C2);
        f[17] += t4 * (-C2);

        // i=9 and i=18
        C2 = prefactor * (-ky + kz);
        f[9] += t4 * (C2);
        f[18] += t4 * (-C2);
    }
};

// Efficient specialization for D3Q15 lattice
template <typename T>
struct asinariTemplatesImpl<T, descriptors::D3Q15DescriptorBase<T> > {
    typedef descriptors::D3Q15DescriptorBase<T> Descriptor;

    static T bgk_ma2_nonsymmetric_equilibrium(plint iPop, Array<T, 3> const &j)
    {
        T c_j = Descriptor::c[iPop][0] * j[0] + Descriptor::c[iPop][1] * j[1]
                + Descriptor::c[iPop][2] * j[2];
        return Descriptor::t[iPop] * 3. * c_j;
    }

    static T bgk_collision_stage1(
        Array<T, Descriptor::q> &f, T rhoBar, T invRho, Array<T, Descriptor::d> const &j, T jSqr,
        T prefactor)
    {
        T p = (T)1 - prefactor;

        T t0 = Descriptor::t[0];
        T t1 = Descriptor::t[1];
        T t4 = Descriptor::t[4];

        T kx = (T)3 * j[0];
        T ky = (T)3 * j[1];
        T kz = (T)3 * j[2];
        T kxSqr_ = invRho / (T)2 * kx * kx;
        T kySqr_ = invRho / (T)2 * ky * ky;
        T kzSqr_ = invRho / (T)2 * kz * kz;
        T kxky_ = invRho * kx * ky;
        T kxkz_ = invRho * kx * kz;
        T kykz_ = invRho * ky * kz;

        T C1 = rhoBar + invRho * (T)3 * jSqr;
        T C2, C3;

        // i=0
        C3 = -kxSqr_ - kySqr_ - kzSqr_;
        f[0] = t0 * (C1 + C3);

        // i=1 and i=8
        C2 = p * (-kx);
        C3 = -kySqr_ - kzSqr_;
        f[1] = t1 * (C1 + C2 + C3);
        f[8] = t1 * (C1 - C2 + C3);

        // i=2 and i=9
        C2 = p * (-ky);
        C3 = -kxSqr_ - kzSqr_;
        f[2] = t1 * (C1 + C2 + C3);
        f[9] = t1 * (C1 - C2 + C3);

        // i=3 and i=10
        C2 = p * (-kz);
        C3 = -kxSqr_ - kySqr_;
        f[3] = t1 * (C1 + C2 + C3);
        f[10] = t1 * (C1 - C2 + C3);

        // i=4 and i=11
        C2 = p * (-kx - ky - kz);
        C3 = kxky_ + kxkz_ + kykz_;
        f[4] = t4 * (C1 + C2 + C3);
        f[11] = t4 * (C1 - C2 + C3);

        // i=5 and i=12
        C2 = p * (-kx - ky + kz);
        C3 = kxky_ - kxkz_ - kykz_;
        f[5] = t4 * (C1 + C2 + C3);
        f[12] = t4 * (C1 - C2 + C3);

        // i=6 and i=13
        C2 = p * (-kx + ky - kz);
        C3 = -kxky_ + kxkz_ - kykz_;
        f[6] = t4 * (C1 + C2 + C3);
        f[13] = t4 * (C1 - C2 + C3);

        // i=7 and i=14
        C2 = p * (-kx + ky + kz);
        C3 = -kxky_ - kxkz_ + kykz_;
        f[7] = t4 * (C1 + C2 + C3);
        f[14] = t4 * (C1 - C2 + C3);

        return jSqr * util::sqr(invRho);
    }

    static void bgk_collision_stage3(
        Array<T, Descriptor::q> &f, Array<T, Descriptor::d> const &j, T prefactor)
    {
        T t1 = Descriptor::t[1];
        T t4 = Descriptor::t[4];

        T kx = (T)3 * j[0];
        T ky = (T)3 * j[1];
        T kz = (T)3 * j[2];

        T C2;

        // i=1 and i=8
        C2 = prefactor * (-kx);
        f[1] += t1 * (C2);
        f[8] += t1 * (-C2);

        // i=2 and i=9
        C2 = prefactor * (-ky);
        f[2] += t1 * (C2);
        f[9] += t1 * (-C2);

        // i=3 and i=10
        C2 = prefactor * (-kz);
        f[3] += t1 * (C2);
        f[10] += t1 * (-C2);

        // i=4 and i=11
        C2 = prefactor * (-kx - ky - kz);
        f[4] += t4 * (C2);
        f[11] += t4 * (-C2);

        // i=5 and i=12
        C2 = prefactor * (-kx - ky + kz);
        f[5] += t4 * (C2);
        f[12] += t4 * (-C2);

        // i=6 and i=13
        C2 = prefactor * (-kx + ky - kz);
        f[6] += t4 * (C2);
        f[13] += t4 * (-C2);

        // i=7 and i=14
        C2 = prefactor * (-kx + ky + kz);
        f[7] += t4 * (C2);
        f[14] += t4 * (-C2);
    }
};

}  // namespace plb

#endif  // ASINARI_TEMPLATES_H
