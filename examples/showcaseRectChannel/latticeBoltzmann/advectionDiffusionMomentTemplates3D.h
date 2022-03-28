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
 * Helper functions for the computation of velocity moments for the f's.
 * This file is all about efficiency. The generic template code is specialized
 * for commonly used Lattices, so that a maximum performance can be taken out
 * of each case.
 */
#ifndef ADVECTION_DIFFUSION_MOMENTS_TEMPLATES_3D_H
#define ADVECTION_DIFFUSION_MOMENTS_TEMPLATES_3D_H

namespace plb {

// This structure forwards the calls to the appropriate helper class
template <typename T>
struct advectionDiffusionMomentTemplatesImpl<T, descriptors::D3Q7DescriptorBase<T> > {
    typedef descriptors::D3Q7DescriptorBase<T> Descriptor;

    static void get_rhoBar_jEq(
        Array<T, Descriptor::q> const &f, T &rhoBar, Array<T, Descriptor::d> &jEq,
        const T u[Descriptor::d])
    {
        rhoBar = momentTemplatesImpl<T, Descriptor>::get_rhoBar(f);
        T rho = Descriptor::fullRho(rhoBar);
        jEq[0] = rho * u[0];
        jEq[1] = rho * u[1];
        jEq[2] = rho * u[2];
    }

    static void get_jEq(const T &rhoBar, Array<T, Descriptor::d> &jEq, const T *u)
    {
        T rho = Descriptor::fullRho(rhoBar);
        jEq[0] = rho * u[0];
        jEq[1] = rho * u[1];
        jEq[2] = rho * u[2];
    }

    static void get_rhoBar_jEq_jNeq_linear(
        Array<T, Descriptor::q> const &f, T &rhoBar, Array<T, Descriptor::d> &jEq,
        Array<T, Descriptor::d> &jNeq, const T *u)
    {
        rhoBar = momentTemplatesImpl<T, Descriptor>::get_rhoBar(f);
        jEq[0] = u[0];
        jEq[1] = u[1];
        jEq[2] = u[2];

        jNeq[0] = -f[1] + f[4] - jEq[0];
        jNeq[1] = -f[2] + f[5] - jEq[1];
        jNeq[2] = -f[3] + f[6] - jEq[2];
    }

    static void get_rhoBar_jEq_jNeq(
        Array<T, Descriptor::q> const &f, T &rhoBar, Array<T, Descriptor::d> &jEq,
        Array<T, Descriptor::d> &jNeq, const T *u)
    {
        rhoBar = momentTemplatesImpl<T, Descriptor>::get_rhoBar(f);
        T rho = Descriptor::fullRho(rhoBar);
        jEq[0] = rho * u[0];
        jEq[1] = rho * u[1];
        jEq[2] = rho * u[2];

        jNeq[0] = -f[1] + f[4] - jEq[0];
        jNeq[1] = -f[2] + f[5] - jEq[1];
        jNeq[2] = -f[3] + f[6] - jEq[2];
    }

};  // struct advectionDiffusionMomentTemplatesImpl

}  // namespace plb

#endif
