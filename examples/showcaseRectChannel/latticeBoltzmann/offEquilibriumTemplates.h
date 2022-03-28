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
 * Specialized helper functions for advanced techniques around LB
 * implementations. They implement the physics of the first-order terms
 * of the Chapman-Enskog expansion and are useful whenever a transition
 * from hydrodynamical variables (rho, u) to kinetic variables (f) si to
 * be implemented. Additionally, they are used for the implementation of
 * the stable RLB dynamics.
 *
 * This file is all about efficiency. The generic
 * template code is specialized for commonly used Lattices, so that a
 * maximum performance can be taken out of each case.
 */
#ifndef OFF_EQUILIBRIUM_TEMPLATES_H
#define OFF_EQUILIBRIUM_TEMPLATES_H

#include "core/cell.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "geometricOperationTemplates.h"
#include "hermitePolynomialsTemplates.h"

namespace plb {

template <typename T, class Descriptor>
struct offEquilibriumTemplatesImpl;

/// General first-order functions
template <typename T, template <typename U> class Descriptor>
struct offEquilibriumTemplates {
    /// Compute off-equilibrium part of the f's from the stress tensor Pi.
    /** Implements the following formula (with Einstein index contraction):
     * /f[ f_i^{neq} = t_i / (2 c_s^4) *
     *                 (c_{ia} c_{ib} - c_s^2 \delta_{ab}) \Pi_{ab} /f]
     * By Pi we mean the tensor computed from the off-equilibrium functions:
     * /f[ \Pi = \sum c_i c_i f_i^{neq}
     *         = \sum c_i c_i f_i - \rho u u - c_s^2 \rho\ Id /f]
     */
    static T fromPiToFneq(plint iPop, Array<T, SymmetricTensor<T, Descriptor>::n> const &pi)
    {
        return offEquilibriumTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::fromPiToFneq(
            iPop, pi);
    }

    static T fromPiAndQtoFneq(
        plint iPop, Array<T, SymmetricTensor<T, Descriptor>::n> const &pi,
        Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> const &q)
    {
        return offEquilibriumTemplatesImpl<
            T, typename Descriptor<T>::BaseDescriptor>::fromPiAndQtoFneq(iPop, pi, q);
    }

    /// Compute off-equilibrium part of the f's from the strain rate tensor S.
    /** Implements the following formula:
     * /f[ f_i^{neq} = - t_i / (c_s^2\omega) *
     *                 (c_{ia} c_{ib} - c_s^2 \delta_{ab}) S_{ab} /f]
     * By S we mean the tensor computed from the velocity gradients:
     * /f[ S_{\alpha\beta} = 1/2 (
     *     \partial_\alpha(\rho u_\beta) + \partial_\beta(\rho u_\alpha) ) /f]
     */
    static T fromStrainToFneq(
        plint iPop, Array<T, SymmetricTensor<T, Descriptor>::n> const &S, T density, T omega)
    {
        return offEquilibriumTemplatesImpl<
            T, typename Descriptor<T>::BaseDescriptor>::fromStrainToFneq(iPop, S, density, omega);
    }

};  // struct offEquilibriumTemplates

template <typename T, class Descriptor>
struct offEquilibriumTemplatesImpl {
    static T fromPiToFneq(
        plint iPop, Array<T, SymmetricTensorImpl<T, Descriptor::d>::n> const &PiNeq)
    {
        typedef Descriptor L;
        T fNeq = T();
        plint iPi = 0;
        // Iterate only over superior triangle + diagonal, and add
        // the elements under the diagonal by symmetry
        for (int iAlpha = 0; iAlpha < L::d; ++iAlpha) {
            // Treat diagonal term first
            fNeq += PiNeq[iPi] * (L::c[iPop][iAlpha] * L::c[iPop][iAlpha] - L::cs2);
            ++iPi;
            // Then, treat off-diagonal terms
            for (int iBeta = iAlpha + 1; iBeta < L::d; ++iBeta) {
                // Multiply off-diagonal elements by 2 because
                // the Q tensor is symmetric
                fNeq += PiNeq[iPi] * (T)2 * L::c[iPop][iAlpha] * L::c[iPop][iBeta];
                ++iPi;
            }
        }
        fNeq *= L::t[iPop] * L::invCs2 * L::invCs2 / (T)2;
        return fNeq;
    }

    static T fromPiAndQtoFneq(
        plint iPop, Array<T, SymmetricTensorImpl<T, Descriptor::d>::n> const &PiNeq,
        Array<T, SymmetricRankThreeTensorImpl<T, Descriptor::d>::n> const &Q)
    {
        typedef Descriptor L;

        Array<T, SymmetricTensorImpl<T, Descriptor::d>::n> H2 =
            HermiteTemplateImpl<T, Descriptor, Descriptor::d>::order2(iPop);
        Array<T, SymmetricRankThreeTensorImpl<T, Descriptor::d>::n> H3 =
            HermiteTemplateImpl<T, Descriptor, Descriptor::d>::order3(iPop);

        T factor = L::t[iPop] * L::invCs2 * L::invCs2 / (T)2;

        T fNeqPi = factor * SymmetricTensorImpl<T, Descriptor::d>::contractIndexes(H2, PiNeq);
        T fNeqQ = factor * L::invCs2 / (T)3
                  * SymmetricRankThreeTensorImpl<T, Descriptor::d>::contractIndexes(H3, Q);

        return fNeqPi + fNeqQ;
    }

    /// Compute off-equilibrium part of the f's from the strain rate tensor S.
    /** Implements the following formula:
     * /f[ f_i^{neq} = - t_i / (c_s^2\omega) *
     *                 (c_{ia} c_{ib} - c_s^2 \delta_{ab}) S_{ab} /f]
     * By S we mean the tensor computed from the velocity gradients:
     * /f[ S_{\alpha\beta} = 1/2 (
     *     \partial_\alpha(\rho u_\beta) + \partial_\beta(\rho u_\alpha) ) /f]
     */
    static T fromStrainToFneq(
        plint iPop, Array<T, SymmetricTensorImpl<T, Descriptor::d>::n> const &S, T density, T omega)
    {
        typedef Descriptor L;
        T fNeq = fromPiToFneq(iPop, S) * (-(T)2 * density * L::cs2 / omega);
        return fNeq;
    }

};  // struct offEquilibriumTemplates

}  // namespace plb

#include "latticeBoltzmann/offEquilibriumTemplates2D.h"
#include "latticeBoltzmann/offEquilibriumTemplates3D.h"

#endif
