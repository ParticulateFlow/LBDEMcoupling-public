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
#ifndef ADVECTION_DIFFUSION_MOMENTS_TEMPLATES_H
#define ADVECTION_DIFFUSION_MOMENTS_TEMPLATES_H

#include "core/cell.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/roundOffPolicy.h"

namespace plb {

template <typename T, class Descriptor>
struct advectionDiffusionMomentTemplatesImpl;

// This structure forwards the calls to the appropriate helper class
template <typename T, template <typename U> class Descriptor>
struct advectionDiffusionMomentTemplates {
    static void get_rhoBar_jEq(
        Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &jEq)
    {
        advectionDiffusionMomentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::
            get_rhoBar_jEq(
                cell.getRawPopulations(), rhoBar, jEq,
                cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));
    }

    static void get_jEq(
        Cell<T, Descriptor> const &cell, const T &rhoBar, Array<T, Descriptor<T>::d> &jEq)
    {
        advectionDiffusionMomentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::get_jEq(
            rhoBar, jEq, cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));
    }

    static void get_rhoBar_jEq_jNeq(
        Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &jEq,
        Array<T, Descriptor<T>::d> &jNeq)
    {
        advectionDiffusionMomentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::
            get_rhoBar_jEq_jNeq(
                cell.getRawPopulations(), rhoBar, jEq, jNeq,
                cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));
    }

    static void get_rhoBar_jEq_jNeq_linear(
        Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &jEq,
        Array<T, Descriptor<T>::d> &jNeq)
    {
        advectionDiffusionMomentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::
            get_rhoBar_jEq_jNeq_linear(
                cell.getRawPopulations(), rhoBar, jEq, jNeq,
                cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));
    }

};  // struct advectionDiffusionMomentTemplates

// This structure forwards the calls to the appropriate helper class
template <typename T, class Descriptor>
struct advectionDiffusionMomentTemplatesImpl {
    static void get_rhoBar_jEq(
        Array<T, Descriptor::q> const &f, T &rhoBar, Array<T, Descriptor::d> &jEq,
        const T u[Descriptor::d])
    {
        rhoBar = momentTemplatesImpl<T, Descriptor>::get_rhoBar(f);
        T rho = Descriptor::fullRho(rhoBar);
        for (plint iD = 0; iD < Descriptor::d; ++iD) {
            jEq[iD] = rho * u[iD];
        }
    }

    static void get_jEq(const T &rhoBar, Array<T, Descriptor::d> &jEq, const T u[Descriptor::d])
    {
        T rho = Descriptor::fullRho(rhoBar);
        for (plint iD = 0; iD < Descriptor::d; ++iD) {
            jEq[iD] = rho * u[iD];
        }
    }

    static void get_rhoBar_jEq_jNeq(
        Array<T, Descriptor::q> const &f, T &rhoBar, Array<T, Descriptor::d> &jEq,
        Array<T, Descriptor::d> &jNeq, const T u[Descriptor::d])
    {
        rhoBar = momentTemplatesImpl<T, Descriptor>::get_rhoBar(f);
        T rho = Descriptor::fullRho(rhoBar);
        Array<T, Descriptor::d> jReal;  // sum f_i*c_i
        for (plint iD = 0; iD < Descriptor::d; ++iD) {
            jReal[iD] = (T)Descriptor::c[0][iD] * f[0];
            for (plint iPop = 1; iPop < Descriptor::q; ++iPop) {
                jReal[iD] += (T)Descriptor::c[iPop][iD] * f[iPop];
            }
            jEq[iD] = rho * u[iD];
            jNeq[iD] = jReal[iD] - jEq[iD];
        }
    }

    static void get_rhoBar_jEq_jNeq_linear(
        Array<T, Descriptor::q> const &f, T &rhoBar, Array<T, Descriptor::d> &jEq,
        Array<T, Descriptor::d> &jNeq, const T u[Descriptor::d])
    {
        rhoBar = momentTemplatesImpl<T, Descriptor>::get_rhoBar(f);
        T rho = Descriptor::fullRho(rhoBar);
        Array<T, Descriptor::d> jReal;  // sum f_i*c_i
        for (plint iD = 0; iD < Descriptor::d; ++iD) {
            jReal[iD] = (T)Descriptor::c[0][iD] * f[0];
            for (plint iPop = 1; iPop < Descriptor::q; ++iPop) {
                jReal[iD] += (T)Descriptor::c[iPop][iD] * f[iPop];
            }
            jEq[iD] = u[iD];
            jNeq[iD] = jReal[iD] - jEq[iD];
        }
    }

};  // struct advectionDiffusionMomentTemplatesImpl

}  // namespace plb

#include "latticeBoltzmann/advectionDiffusionMomentTemplates2D.h"
#include "latticeBoltzmann/advectionDiffusionMomentTemplates3D.h"

#endif
