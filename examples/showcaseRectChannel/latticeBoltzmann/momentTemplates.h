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
#ifndef MOMENT_TEMPLATES_H
#define MOMENT_TEMPLATES_H

#include "core/cell.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/roundOffPolicy.h"

namespace plb {

template <typename T, class Descriptor>
struct momentTemplatesImpl;

// This structure forwards the calls to the appropriate helper class
template <typename T, template <typename U> class Descriptor>
struct momentTemplates {
    static T get_rhoBar(Cell<T, Descriptor> const &cell)
    {
        return momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::get_rhoBar(
            cell.getRawPopulations());
    }

    static void get_j(Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &j)
    {
        momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::get_j(
            cell.getRawPopulations(), j);
    }

    static T get_eBar(Cell<T, Descriptor> const &cell)
    {
        return momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::get_eBar(
            cell.getRawPopulations());
    }

    static void get_rhoBar_j(
        Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j)
    {
        momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::get_rhoBar_j(
            cell.getRawPopulations(), rhoBar, j);
    }

    static void get_rhoBar_j_thetaBar(
        Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j, T &thetaBar)
    {
        momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::get_rhoBar_j_thetaBar(
            cell.getRawPopulations(), rhoBar, j, thetaBar);
    }

    static T compute_rho(Cell<T, Descriptor> const &cell)
    {
        return momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::compute_rho(
            cell.getRawPopulations());
    }

    /// Get order-1 moment of f's, divided by rho ("lattice-boltzmann-velocity", or "uLb")
    /** In many cases, such as the plain BGK model, this is equal to to physical
     * velocity, but in other cases not. In presence of a body force g for example,
     * the velocity is uLb + g/2.
     **/
    static void compute_uLb(Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &uLb)
    {
        momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::compute_uLb(
            cell.getRawPopulations(), uLb);
    }

    static void compute_rho_uLb(
        Cell<T, Descriptor> const &cell, T &rho, Array<T, Descriptor<T>::d> &uLb)
    {
        momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::compute_rho_uLb(
            cell.getRawPopulations(), rho, uLb);
    }

    static T compute_e(Cell<T, Descriptor> const &cell)
    {
        return momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::compute_e(
            cell.getRawPopulations());
    }

    static T compute_rhoThetaBar(Cell<T, Descriptor> const &cell, T rhoBar, T jSqr)
    {
        return momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::compute_rhoThetaBar(
            cell.getRawPopulations(), rhoBar, jSqr);
    }

    static void compute_rho_rhoThetaBar(Cell<T, Descriptor> const &cell, T &rho, T &rhoThetaBar)
    {
        momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::compute_rho_rhoThetaBar(
            cell.getRawPopulations(), rho, rhoThetaBar);
    }

    static T compute_theta(Cell<T, Descriptor> const &cell, T rhoBar, T jSqr)
    {
        return momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::compute_theta(
            cell.getRawPopulations(), rhoBar, jSqr);
    }

    static T compute_rhoEpsilon(Cell<T, Descriptor> const &cell, T rhoBar, T jSqr)
    {
        return momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::compute_rhoEpsilon(
            cell.getRawPopulations(), rhoBar, jSqr);
    }

    static void compute_PiNeq(
        Cell<T, Descriptor> const &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq)
    {
        momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::compute_PiNeq(
            cell.getRawPopulations(), rhoBar, j, PiNeq, Descriptor<T>::invRho(rhoBar));
    }

    static void compute_PiNeq(
        Cell<T, Descriptor> const &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq, T invRho)
    {
        momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::compute_PiNeq(
            cell.getRawPopulations(), rhoBar, j, PiNeq, invRho);
    }

    static void compute_thermal_PiNeq(
        Cell<T, Descriptor> const &cell, T rhoBar, T thetaBar, Array<T, Descriptor<T>::d> const &j,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq)
    {
        momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::compute_thermal_PiNeq(
            cell.getRawPopulations(), rhoBar, thetaBar, j, PiNeq);
    }

    static void compute_rhoBar_j_PiNeq(
        Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq)
    {
        momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::compute_rhoBar_j_PiNeq(
            cell.getRawPopulations(), rhoBar, j, PiNeq);
    }

    static void compute_rhoBar_j_PiNeq(
        Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq, T invRho)
    {
        momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::compute_rhoBar_j_PiNeq(
            cell.getRawPopulations(), rhoBar, j, PiNeq, invRho);
    }

    static void compute_rhoBar_thetaBar_j_PiNeq(
        Cell<T, Descriptor> const &cell, T &rhoBar, T &thetaBar, Array<T, Descriptor<T>::d> &j,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq)
    {
        momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::
            compute_rhoBar_thetaBar_j_PiNeq(cell.getRawPopulations(), rhoBar, thetaBar, j, PiNeq);
    }

    static void compute_rhoBar_thetaBar_j_PiNeq_qNeq(
        Cell<T, Descriptor> const &cell, T &rhoBar, T &thetaBar, Array<T, Descriptor<T>::d> &j,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq,
        Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> &qNeq)
    {
        momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::
            compute_rhoBar_thetaBar_j_PiNeq_qNeq(
                cell.getRawPopulations(), rhoBar, thetaBar, j, PiNeq, qNeq);
    }

    /// Get local, order-2 moment: sum_i (c_i-uLb)(c_i-uLb) f_i = -rho uLb uLb + sum_i c_i c_i f_i
    /** The full stress tensor Pi is equal to P + rho u u. The deviatoric stress tensor sigma
     *  is equal to P - c_s^2 rho I (\sa compute_Pi_neq)
     **/
    static void compute_P(
        Cell<T, Descriptor> const &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j,
        Array<T, SymmetricTensor<T, Descriptor>::n> &P)
    {
        momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::compute_P(
            cell.getRawPopulations(), rhoBar, j, P);
    }

    static void modifyJ(T &cell, Array<T, Descriptor<T>::d> const &newJ)
    {
        momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::modifyVelocity(
            cell.getRawPopulations(), newJ);
    }

    static void compute_Qneq(
        Cell<T, Descriptor> const &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> &qNeq)
    {
        momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::compute_Qneq(
            cell.getRawPopulations(), rhoBar, j, thetaBar, qNeq);
    }

    static void compute_heat_flux(
        Cell<T, Descriptor> const &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        Array<T, Descriptor<T>::d> &q)
    {
        momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::compute_heat_flux(
            cell.getRawPopulations(), rhoBar, j, thetaBar, q);
    }

};  // struct momentTemplates

// This structure forwards the calls to the appropriate helper class
template <typename T, class Descriptor>
struct momentTemplatesImpl {
    static T get_rhoBar(Array<T, Descriptor::q> const &f)
    {
        T rhoBar = f[0];
        for (plint iPop = 1; iPop < Descriptor::q; ++iPop) {
            rhoBar += f[iPop];
        }
        return rhoBar;
    }

    static void get_j(Array<T, Descriptor::q> const &f, Array<T, Descriptor::d> &j)
    {
        for (int iD = 0; iD < Descriptor::d; ++iD) {
            j[iD] = f[0] * Descriptor::c[0][iD];
        }
        for (plint iPop = 1; iPop < Descriptor::q; ++iPop) {
            for (int iD = 0; iD < Descriptor::d; ++iD) {
                j[iD] += f[iPop] * Descriptor::c[iPop][iD];
            }
        }
    }

    static T get_eBar(Array<T, Descriptor::q> const &f)
    {
        T eBar = f[0] * Descriptor::cNormSqr[0];
        for (plint iPop = 1; iPop < Descriptor::q; ++iPop) {
            eBar += f[iPop] * Descriptor::cNormSqr[iPop];
        }
        return eBar;
    }

    static void get_rhoBar_j(
        Array<T, Descriptor::q> const &f, T &rhoBar, Array<T, Descriptor::d> &j)
    {
        rhoBar = get_rhoBar(f);
        get_j(f, j);
    }

    static void get_rhoBar_j_thetaBar(
        Array<T, Descriptor::q> const &f, T &rhoBar, Array<T, Descriptor::d> &j, T &thetaBar)
    {
        get_rhoBar_j(f, rhoBar, j);
        T jSqr = VectorTemplateImpl<T, Descriptor::d>::normSqr(j);
        T invRho = Descriptor::invRho(rhoBar);
        thetaBar = compute_rhoThetaBar(f, rhoBar, jSqr) * invRho;
    }

    static T compute_rho(Array<T, Descriptor::q> const &f)
    {
        return Descriptor::fullRho(get_rhoBar(f));
    }

    static void compute_uLb(Array<T, Descriptor::q> const &f, Array<T, Descriptor::d> &uLb)
    {
        get_j(f, uLb);
        T invRho = Descriptor::invRho(get_rhoBar(f));
        for (int iD = 0; iD < Descriptor::d; ++iD) {
            uLb[iD] *= invRho;
        }
    }

    static void compute_rho_uLb(
        Array<T, Descriptor::q> const &f, T &rho, Array<T, Descriptor::d> &uLb)
    {
        get_j(f, uLb);
        T rhoBar = get_rhoBar(f);
        T invRho = Descriptor::invRho(rhoBar);
        rho = Descriptor::fullRho(rhoBar);
        for (int iD = 0; iD < Descriptor::d; ++iD) {
            uLb[iD] *= invRho;
        }
    }

    static T compute_e(Array<T, Descriptor::q> const &f)
    {
        return get_eBar(f) + Descriptor::SkordosFactor() * Descriptor::d * Descriptor::cs2;
    }

    static T compute_rhoThetaBar(Array<T, Descriptor::q> const &f, T rhoBar, T jSqr)
    {
        T invRho = Descriptor::invRho(rhoBar);
        return Descriptor::invCs2 * Descriptor::invD * (get_eBar(f) - invRho * jSqr) - rhoBar;
    }

    static void compute_rho_rhoThetaBar(Array<T, Descriptor::q> const &f, T &rho, T &rhoThetaBar)
    {
        T rhoBar, j[Descriptor::d];
        get_rhoBar_j(f, rhoBar, j);
        T jSqr = VectorTemplateImpl<T, Descriptor::d>::normSqr(j);
        rho = Descriptor::fullRho(rhoBar);
        rhoThetaBar = compute_rhoThetaBar(f, rhoBar, jSqr);
    }

    static T compute_theta(Array<T, Descriptor::q> const &f, T rhoBar, T jSqr)
    {
        T invRho = Descriptor::invRho(rhoBar);
        T e = compute_e(f);
        return invRho * Descriptor::invD * Descriptor::invCs2 * (e - invRho * jSqr);
    }

    static T compute_rhoEpsilon(Array<T, Descriptor::q> const &f, T rhoBar, T jSqr)
    {
        T invRho = Descriptor::invRho(rhoBar);
        T e = compute_e(f);
        return (e - invRho * jSqr) / (T)2;
    }

    static void compute_PiNeq(
        Array<T, Descriptor::q> const &f, T rhoBar, Array<T, Descriptor::d> const &j,
        Array<T, SymmetricTensorImpl<T, Descriptor::d>::n> &PiNeq)
    {
        T invRho = Descriptor::invRho(rhoBar);
        compute_PiNeq(f, rhoBar, j, PiNeq, invRho);
    }

    static void compute_PiNeq(
        Array<T, Descriptor::q> const &f, T rhoBar, Array<T, Descriptor::d> const &j,
        Array<T, SymmetricTensorImpl<T, Descriptor::d>::n> &PiNeq, T invRho)
    {
        int iPi = 0;
        for (int iAlpha = 0; iAlpha < Descriptor::d; ++iAlpha) {
            int iDiagonal = iPi;
            for (int iBeta = iAlpha; iBeta < Descriptor::d; ++iBeta) {
                PiNeq[iPi] = Descriptor::c[0][iAlpha] * Descriptor::c[0][iBeta] * f[0];
                for (plint iPop = 1; iPop < Descriptor::q; ++iPop) {
                    PiNeq[iPi] +=
                        Descriptor::c[iPop][iAlpha] * Descriptor::c[iPop][iBeta] * f[iPop];
                }
                // Stripe off relative velocity
                PiNeq[iPi] -= invRho * j[iAlpha] * j[iBeta];
                ++iPi;
            }
            // Stripe off diagonal term
            PiNeq[iDiagonal] -= Descriptor::cs2 * rhoBar;
        }
    }

    static void compute_PiNeq(
        Array<T, Descriptor::q> const &fNeq,
        Array<T, SymmetricTensorImpl<T, Descriptor::d>::n> &PiNeq)
    {
        int iPi = 0;
        for (int iAlpha = 0; iAlpha < Descriptor::d; ++iAlpha) {
            for (int iBeta = iAlpha; iBeta < Descriptor::d; ++iBeta) {
                PiNeq[iPi] = Descriptor::c[0][iAlpha] * Descriptor::c[0][iBeta] * fNeq[0];
                for (plint iPop = 1; iPop < Descriptor::q; ++iPop) {
                    PiNeq[iPi] +=
                        Descriptor::c[iPop][iAlpha] * Descriptor::c[iPop][iBeta] * fNeq[iPop];
                }
                ++iPi;
            }
        }
    }

    static void compute_thermal_PiNeq(
        Array<T, Descriptor::q> const &f, T rhoBar, T thetaBar, Array<T, Descriptor::d> const &j,
        Array<T, SymmetricTensorImpl<T, Descriptor::d>::n> &PiNeq)
    {
        // rhoTheta_bar == rho*theta - 1
        T rhoTheta_bar = rhoBar * thetaBar + rhoBar + Descriptor::SkordosFactor() * thetaBar;

        T invRho = Descriptor::invRho(rhoBar);
        int iPi = 0;
        for (int iAlpha = 0; iAlpha < Descriptor::d; ++iAlpha) {
            int iDiagonal = iPi;
            for (int iBeta = iAlpha; iBeta < Descriptor::d; ++iBeta) {
                PiNeq[iPi] = Descriptor::c[0][iAlpha] * Descriptor::c[0][iBeta] * f[0];
                for (plint iPop = 1; iPop < Descriptor::q; ++iPop) {
                    PiNeq[iPi] +=
                        Descriptor::c[iPop][iAlpha] * Descriptor::c[iPop][iBeta] * f[iPop];
                }
                // Stripe off relative velocity
                PiNeq[iPi] -= invRho * j[iAlpha] * j[iBeta];
                ++iPi;
            }
            // Stripe off diagonal term
            PiNeq[iDiagonal] -= Descriptor::cs2 * rhoTheta_bar;
        }
    }

    static void compute_rhoBar_j_PiNeq(
        Array<T, Descriptor::q> const &f, T &rhoBar, Array<T, Descriptor::d> &j,
        Array<T, SymmetricTensorImpl<T, Descriptor::d>::n> &PiNeq)
    {
        get_rhoBar_j(f, rhoBar, j);
        compute_PiNeq(f, rhoBar, j, PiNeq);
    }

    static void compute_rhoBar_j_PiNeq(
        Array<T, Descriptor::q> const &f, T &rhoBar, Array<T, Descriptor::d> &j,
        Array<T, SymmetricTensorImpl<T, Descriptor::d>::n> &PiNeq, T invRho)
    {
        get_rhoBar_j(f, rhoBar, j);
        compute_PiNeq(f, rhoBar, j, PiNeq, invRho);
    }

    static void compute_rhoBar_thetaBar_j_PiNeq(
        Array<T, Descriptor::q> const &f, T &rhoBar, T &thetaBar, Array<T, Descriptor::d> &j,
        Array<T, SymmetricTensorImpl<T, Descriptor::d>::n> &PiNeq)
    {
        get_rhoBar_j_thetaBar(f, rhoBar, j, thetaBar);
        compute_thermal_PiNeq(f, rhoBar, thetaBar, j, PiNeq);
    }

    static void compute_rhoBar_thetaBar_j_PiNeq_qNeq(
        Array<T, Descriptor::q> const &f, T &rhoBar, T &thetaBar, Array<T, Descriptor::d> &j,
        Array<T, SymmetricTensorImpl<T, Descriptor::d>::n> &PiNeq,
        Array<T, SymmetricRankThreeTensorImpl<T, Descriptor::d>::n> &qNeq)
    {
        compute_rhoBar_thetaBar_j_PiNeq(f, rhoBar, thetaBar, j, PiNeq);
        compute_Qneq(f, rhoBar, j, thetaBar, qNeq);
    }

    static void compute_P(
        Array<T, Descriptor::q> const &f, T rhoBar, Array<T, Descriptor::d> const &j,
        Array<T, SymmetricTensorImpl<T, Descriptor::d>::n> &P)
    {
        T invRho = Descriptor::invRho(rhoBar);
        plint iP = 0;
        for (int iAlpha = 0; iAlpha < Descriptor::d; ++iAlpha) {
            for (int iBeta = iAlpha; iBeta < Descriptor::d; ++iBeta) {
                P[iP] = Descriptor::c[0][iAlpha] * Descriptor::c[0][iBeta] * f[0];
                for (plint iPop = 1; iPop < Descriptor::q; ++iPop) {
                    P[iP] += Descriptor::c[iPop][iAlpha] * Descriptor::c[iPop][iBeta] * f[iPop];
                }
                // Stripe off relative velocity
                P[iP] -= invRho * j[iAlpha] * j[iBeta];
                ++iP;
            }
        }
    }

    static void modifyJ(T *f, Array<T, Descriptor::d> const &newJ)
    {
        T rhoBar;
        Array<T, Descriptor::d> oldJ;
        get_rhoBar_j(f, rhoBar, oldJ);
        T invRho = Descriptor::invRho(rhoBar);
        const T oldJSqr = VectorTemplateImpl<T, Descriptor::d>::normSqr(oldJ);
        const T newJSqr = VectorTemplateImpl<T, Descriptor::d>::normSqr(newJ);
        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            f[iPop] = f[iPop] - bgk_ma2_equilibrium(iPop, rhoBar, invRho, oldJ, oldJSqr)
                      + bgk_ma2_equilibrium(iPop, rhoBar, invRho, newJ, newJSqr);
        }
    }

    static void compute_Qneq(
        Array<T, Descriptor::q> const &f, T rhoBar, Array<T, Descriptor::d> const &j, T thetaBar,
        Array<T, SymmetricRankThreeTensorImpl<T, Descriptor::d>::n> &qNeq)
    {
        typedef Descriptor L;
        T invRho = L::invRho(rhoBar);
        T factor = L::cs2 + thetaBar;
        plint iQ = 0;
        for (plint iA = 0; iA < L::d; ++iA) {
            for (plint iB = iA; iB < L::d; ++iB) {
                for (plint iC = iB; iC < L::d; ++iC) {
                    qNeq[iQ] = -j[iA] * j[iB] * j[iC] * invRho * invRho;
                    for (plint iPop = 0; iPop < L::q; ++iPop) {
                        qNeq[iQ] = L::c[iPop][iA] * L::c[iPop][iB] * L::c[iPop][iC] * f[iPop];
                    }
                    if (iA == iB && iB == iC) {
                        qNeq[iQ] -= (T)3 * factor * j[iA];
                    } else if (iA == iB && iB != iC) {
                        qNeq[iQ] -= factor * j[iC];
                    } else if (iA == iC && iC != iB) {
                        qNeq[iQ] -= factor * j[iB];
                    } else if (iB == iC && iC != iA) {
                        qNeq[iQ] -= factor * j[iA];
                    }

                    ++iQ;
                }
            }
        }
    }

    static void compute_heat_flux(
        Array<T, Descriptor::q> const &f, T rhoBar, Array<T, Descriptor::d> const &j, T thetaBar,
        Array<T, Descriptor::d> &q)
    {
        typedef Descriptor L;
        Array<T, SymmetricRankThreeTensorImpl<T, Descriptor::d>::n> qNeq;
        compute_Qneq(f, rhoBar, j, thetaBar, qNeq);
        SymmetricRankThreeTensorImpl<T, L::d>::contractLastTwoIndexes(qNeq, q);
    }

};  // struct momentTemplatesImpl

}  // namespace plb

#include "latticeBoltzmann/momentTemplates2D.h"
#include "latticeBoltzmann/momentTemplates3D.h"

#endif  // MOMENT_TEMPLATES_H
