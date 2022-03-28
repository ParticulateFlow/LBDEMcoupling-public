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
 * Helper functions for the implementation of force terms in LB dynamics.
 * This file is all about efficiency. The generic template code is specialized
 * for commonly used Lattices, so that a maximum performance can be taken out of
 * each case.
 */
#ifndef EXTERNAL_FORCE_TEMPLATES_H
#define EXTERNAL_FORCE_TEMPLATES_H

#include "core/cell.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

template <typename T, class Descriptor>
struct externalForceTemplatesImpl;

template <typename T, template <typename U> class Descriptor>
struct externalForceTemplates {
    /// Add a force term after BGK collision, according to the "old", linear algorithm.
    static void addNaiveForce(Cell<T, Descriptor> &cell)
    {
        externalForceTemplatesImpl<T, Descriptor<T> >::addNaiveForce(
            cell.getRawPopulations(), cell.getExternal(0));
    }

    /// Add a force term after BGK collision, according to the Guo algorithm
    static void addGuoForce(
        Cell<T, Descriptor> &cell, Array<T, Descriptor<T>::d> const &u, T omega, T amplitude)
    {
        externalForceTemplatesImpl<T, Descriptor<T> >::addGuoForce(
            cell.getRawPopulations(), cell.getExternal(0), u, omega, amplitude);
    }

    /// Add a force term after BGK collision, according to the Shan algorithm (with third order
    /// terms in hermite polynomial)
    static void addGuoCompressibleForce(
        Cell<T, Descriptor> &cell, T rho, Array<T, Descriptor<T>::d> const &u, T uSqr, T thetaBar,
        T omega, T amplitude)
    {
        externalForceTemplatesImpl<T, Descriptor<T> >::addGuoCompressibleForce(
            cell.getRawPopulations(), cell.getExternal(0), rho, u, uSqr, thetaBar, omega,
            amplitude);
    }

    /// BGK collision with Shan/Chen external force.
    static T shanChenForcedBGKCollision(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T omega)
    {
        return externalForceTemplatesImpl<T, Descriptor<T> >::shanChenForcedBGKCollision(
            cell.getRawPopulations(), cell.getExternal(0), rhoBar, j, omega);
    }

    /// BGK collision with He et al. external force.
    static T heForcedBGKCollision(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T omega)
    {
        return externalForceTemplatesImpl<T, Descriptor<T> >::heForcedBGKCollision(
            cell.getRawPopulations(), cell.getExternal(0), rhoBar, j, omega);
    }
};

template <typename T, class Descriptor>
struct externalForceTemplatesImpl {
    static void addNaiveForce(Array<T, Descriptor::q> &f, T *externalScalars)
    {
        static const int forceBeginsAt = Descriptor::ExternalField::forceBeginsAt;
        T *force = externalScalars + forceBeginsAt;
        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            T c_f = T();
            for (int iD = 0; iD < Descriptor::d; ++iD) {
                c_f += Descriptor::c[iPop][iD] * force[iD];
            }
            f[iPop] += c_f * Descriptor::t[iPop] * Descriptor::invCs2;
        }
    }

    static void addGuoForce(
        Array<T, Descriptor::q> &f, T *externalScalars, Array<T, Descriptor::d> const &u, T omega,
        T amplitude)
    {
        static const int forceBeginsAt = Descriptor::ExternalField::forceBeginsAt;
        T *force = externalScalars + forceBeginsAt;
        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            T c_u = T();
            for (int iD = 0; iD < Descriptor::d; ++iD) {
                c_u += Descriptor::c[iPop][iD] * u[iD];
            }
            c_u *= Descriptor::invCs2 * Descriptor::invCs2;
            T forceTerm = T();
            for (int iD = 0; iD < Descriptor::d; ++iD) {
                forceTerm += (((T)Descriptor::c[iPop][iD] - u[iD]) * Descriptor::invCs2
                              + c_u * (T)Descriptor::c[iPop][iD])
                             * force[iD];
            }
            forceTerm *= Descriptor::t[iPop];
            forceTerm *= 1 - omega / (T)2;
            forceTerm *= amplitude;
            f[iPop] += forceTerm;
        }
    }

    static void addGuoCompressibleForce(
        Array<T, Descriptor::q> &f, T *externalScalars, T rho, Array<T, Descriptor::d> const &u,
        T uSqr, T thetaBar, T omega, T amplitude)
    {
        static const int forceBeginsAt = Descriptor::ExternalField::forceBeginsAt;
        T *force = externalScalars + forceBeginsAt;
        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            T c_u = T();
            T f_u = T();
            T c_f = T();
            for (int iD = 0; iD < Descriptor::d; ++iD) {
                c_u += Descriptor::c[iPop][iD] * u[iD];
                f_u += force[iD] * u[iD];
                f_u += Descriptor::c[iPop][iD] * force[iD];
            }

            T forceTerm = c_f * Descriptor::invCs2;  // Hermite order 1
            forceTerm += Descriptor::invCs2 * Descriptor::invCs2 / (T)2
                         * (c_u * c_f - Descriptor::cs2 * f_u);  // Hermite order 2
            forceTerm +=
                Descriptor::invCs2 * Descriptor::invCs2 * Descriptor::invCs2 / (T)6
                * (c_f * (c_u * c_u + Descriptor::cs2 * (Descriptor::cNormSqr[iPop] * thetaBar))
                   - Descriptor::cs2
                         * (c_u * ((T)2 * f_u + uSqr + Descriptor::cs2 * Descriptor::d * thetaBar)
                            + (T)2 * Descriptor::cs2 * thetaBar * c_f));  // Hermite order 3

            forceTerm *= Descriptor::t[iPop];
            forceTerm *= 1 - omega / (T)2;
            forceTerm *= rho;
            forceTerm *= amplitude;

            f[iPop] += forceTerm;
        }
    }

    static T heForcedBGKCollision(
        Array<T, Descriptor::q> &f, T *externalScalars, T rhoBar, Array<T, Descriptor::d> const &j,
        T omega)
    {
        static const int forceBeginsAt = Descriptor::ExternalField::forceBeginsAt;
        T *force = externalScalars + forceBeginsAt;
        T invRho = Descriptor::invRho(rhoBar);
        const T jSqr = VectorTemplateImpl<T, Descriptor::d>::normSqr(j);

        Array<T, Descriptor::d> uLB;
        for (plint iD = 0; iD < Descriptor::d; ++iD) {
            uLB[iD] = invRho * j[iD];
        }
        const T uSqrLB = VectorTemplateImpl<T, Descriptor::d>::normSqr(uLB);

        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            T ug = 0.;
            for (int iD = 0; iD < Descriptor::d; ++iD) {
                ug += (Descriptor::c[iPop][iD] - uLB[iD]) * force[iD];
            }
            f[iPop] *= (T)1 - omega;
            f[iPop] += omega
                       * dynamicsTemplatesImpl<T, Descriptor>::bgk_ma2_equilibrium(
                           iPop, rhoBar, invRho, j, jSqr);

            f[iPop] += (1 - omega / (T)2) * Descriptor::invCs2 * ug
                       * dynamicsTemplatesImpl<T, Descriptor>::bgk_ma2_equilibrium(
                           iPop, (T)1, (T)1, uLB, uSqrLB);
        }

        T uSqr = T();
        for (int iD = 0; iD < Descriptor::d; ++iD) {
            uSqr += util::sqr(uLB[iD]);
        }

        return uSqr;
    }

};  // struct externalForceTemplates

}  // namespace plb

#include "latticeBoltzmann/externalForceTemplates2D.h"
#include "latticeBoltzmann/externalForceTemplates3D.h"

#endif
