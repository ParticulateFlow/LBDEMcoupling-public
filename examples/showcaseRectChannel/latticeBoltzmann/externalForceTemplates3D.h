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
 * 3D specialization of externalForceTemplates functions
 */
#ifndef EXTERNAL_FORCE_TEMPLATES_3D_H
#define EXTERNAL_FORCE_TEMPLATES_3D_H

namespace plb {

template <typename T>
struct externalForceTemplatesImpl<T, descriptors::ForcedD3Q19Descriptor<T> > {
    static void addNaiveForce(
        Array<T, descriptors::ForcedD3Q19Descriptor<T>::q> &f, T *externalScalars)
    {
        static const int forceBeginsAt =
            descriptors::ForcedD3Q19Descriptor<T>::ExternalField::forceBeginsAt;
        T *force = externalScalars + forceBeginsAt;
        T &fx = force[0];
        T &fy = force[1];
        T &fz = force[2];

        f[1] += (-fx) * (T)1 / (T)6;
        f[2] += (-fy) * (T)1 / (T)6;
        f[3] += (-fz) * (T)1 / (T)6;
        f[4] += (-fx - fy) * (T)1 / (T)12;
        f[5] += (-fx + fy) * (T)1 / (T)12;
        f[6] += (-fx - fz) * (T)1 / (T)12;
        f[7] += (-fx + fz) * (T)1 / (T)12;
        f[8] += (-fy - fz) * (T)1 / (T)12;
        f[9] += (-fy + fz) * (T)1 / (T)12;

        f[10] += (fx) * (T)1 / (T)6;
        f[11] += (fy) * (T)1 / (T)6;
        f[12] += (fz) * (T)1 / (T)6;
        f[13] += (fx + fy) * (T)1 / (T)12;
        f[14] += (fx - fy) * (T)1 / (T)12;
        f[15] += (fx + fz) * (T)1 / (T)12;
        f[16] += (fx - fz) * (T)1 / (T)12;
        f[17] += (fy + fz) * (T)1 / (T)12;
        f[18] += (fy - fz) * (T)1 / (T)12;
    }

    static void addGuoForce(
        Array<T, descriptors::ForcedD3Q19Descriptor<T>::q> &f, T *externalScalars,
        Array<T, descriptors::ForcedD3Q19Descriptor<T>::d> const &u, T omega, T amplitude)
    {
        static const int forceBeginsAt =
            descriptors::ForcedD3Q19Descriptor<T>::ExternalField::forceBeginsAt;
        T *force = externalScalars + forceBeginsAt;
        T mu = amplitude * ((T)1 - omega / (T)2);

        static const T oneOver6 = (T)1 / (T)6;
        static const T oneOver12 = (T)1 / (T)12;

        f[0] += -mu * (force[0] * u[0] + force[1] * u[1] + force[2] * u[2]);

        f[1] += oneOver6 * mu * (force[0] * (-(T)1 + 2 * u[0]) - force[1] * u[1] - force[2] * u[2]);

        f[2] += -oneOver6 * mu * (force[0] * u[0] + force[1] * ((T)1 - 2 * u[1]) + force[2] * u[2]);

        f[3] += -oneOver6 * mu * (force[0] * u[0] + force[1] * u[1] + force[2] * ((T)1 - 2 * u[2]));

        f[4] += oneOver12 * mu
                * (force[0] * (-(T)1 + 2 * u[0] + 3 * u[1])
                   + force[1] * (-(T)1 + 2 * u[1] + 3 * u[0]) - force[2] * u[2]);

        f[5] += oneOver12 * mu
                * (force[0] * (-(T)1 + 2 * u[0] - 3 * u[1])
                   + force[1] * ((T)1 + 2 * u[1] - 3 * u[0]) - force[2] * u[2]);

        f[6] += oneOver12 * mu
                * (force[0] * (-(T)1 + 2 * u[0] + 3 * u[2]) - force[1] * u[1]
                   + force[2] * (-(T)1 + 2 * u[2] + 3 * u[0]));

        f[7] += oneOver12 * mu
                * (force[0] * (-(T)1 + 2 * u[0] - 3 * u[2]) - force[1] * u[1]
                   + force[2] * ((T)1 + 2 * u[2] - 3 * u[0]));

        f[8] += -oneOver12 * mu
                * (force[0] * u[0] + force[1] * ((T)1 - 2 * u[1] - 3 * u[2])
                   + force[2] * ((T)1 - 2 * u[2] - 3 * u[1]));

        f[9] += -oneOver12 * mu
                * (force[0] * u[0] + force[1] * ((T)1 - 2 * u[1] + 3 * u[2])
                   + force[2] * (-(T)1 - 2 * u[2] + 3 * u[1]));

        f[10] += oneOver6 * mu * (force[0] * ((T)1 + 2 * u[0]) - force[1] * u[1] - force[2] * u[2]);

        f[11] +=
            -oneOver6 * mu * (force[0] * u[0] + force[1] * (-(T)1 - 2 * u[1]) + force[2] * u[2]);

        f[12] +=
            -oneOver6 * mu * (force[0] * u[0] + force[1] * u[1] + force[2] * (-(T)1 - 2 * u[2]));

        f[13] += oneOver12 * mu
                 * (force[0] * ((T)1 + 2 * u[0] + 3 * u[1])
                    + force[1] * ((T)1 + 2 * u[1] + 3 * u[0]) - force[2] * u[2]);

        f[14] += oneOver12 * mu
                 * (force[0] * ((T)1 + 2 * u[0] - 3 * u[1])
                    + force[1] * (-(T)1 + 2 * u[1] - 3 * u[0]) - force[2] * u[2]);

        f[15] += oneOver12 * mu
                 * (force[0] * ((T)1 + 2 * u[0] + 3 * u[2]) - force[1] * u[1]
                    + force[2] * ((T)1 + 2 * u[2] + 3 * u[0]));

        f[16] += oneOver12 * mu
                 * (force[0] * ((T)1 + 2 * u[0] - 3 * u[2]) - force[1] * u[1]
                    + force[2] * (-(T)1 + 2 * u[2] - 3 * u[0]));

        f[17] += -oneOver12 * mu
                 * (force[0] * u[0] + force[1] * (-(T)1 - 2 * u[1] - 3 * u[2])
                    + force[2] * (-(T)1 - 2 * u[2] - 3 * u[1]));

        f[18] += -oneOver12 * mu
                 * (force[0] * u[0] + force[1] * (-(T)1 - 2 * u[1] + 3 * u[2])
                    + force[2] * ((T)1 - 2 * u[2] + 3 * u[1]));
    }

    static T heForcedBGKCollision(
        Array<T, descriptors::ForcedD3Q19Descriptor<T>::q> &f, T *externalScalars, T rhoBar,
        Array<T, descriptors::ForcedD3Q19Descriptor<T>::d> const &j, T omega)
    {
        static const int forceBeginsAt =
            descriptors::ForcedD3Q19Descriptor<T>::ExternalField::forceBeginsAt;
        T *force = externalScalars + forceBeginsAt;
        T invRho = descriptors::ForcedD3Q19Descriptor<T>::invRho(rhoBar);

        Array<T, descriptors::ForcedD3Q19Descriptor<T>::d> uLB(
            j[0] * invRho, j[1] * invRho, j[2] * invRho);
        const T uSqrLB =
            VectorTemplateImpl<T, descriptors::ForcedD3Q19Descriptor<T>::d>::normSqr(uLB);

        dynamicsTemplatesImpl<T, descriptors::D3Q19DescriptorBase<T> >::bgk_ma2_collision(
            f, rhoBar, j, omega);

        for (plint iPop = 0; iPop < descriptors::ForcedD3Q19Descriptor<T>::q; ++iPop) {
            T ug = (descriptors::ForcedD3Q19Descriptor<T>::c[iPop][0] - uLB[0]) * force[0]
                   + (descriptors::ForcedD3Q19Descriptor<T>::c[iPop][1] - uLB[1]) * force[1]
                   + (descriptors::ForcedD3Q19Descriptor<T>::c[iPop][2] - uLB[2]) * force[2];

            f[iPop] += (1 - omega / (T)2) * descriptors::ForcedD3Q19Descriptor<T>::invCs2 * ug
                       * dynamicsTemplatesImpl<T, descriptors::ForcedD3Q19Descriptor<T> >::
                           bgk_ma2_equilibrium(iPop, (T)1, (T)1, uLB, uSqrLB);
        }
        T uSqr = util::sqr(uLB[0]) + util::sqr(uLB[1]) + util::sqr(uLB[2]);
        return uSqr;
    }
};

}  // namespace plb

#endif
