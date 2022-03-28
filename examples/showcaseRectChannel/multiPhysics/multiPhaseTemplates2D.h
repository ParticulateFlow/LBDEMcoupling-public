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
 * Helper functions for the implementation of multi-phase models. This file is all
 * about efficiency. The generic template code is specialized for commonly
 * used Lattices, so that a maximum performance can be taken out of each
 * case.
 */
#ifndef MULTI_PHASE_TEMPLATES_2D_H
#define MULTI_PHASE_TEMPLATES_2D_H

#include "core/cell.h"
#include "core/globalDefs.h"
#include "multiPhysics/shanChenLattices2D.h"

namespace plb {

/// Helper functions with full-lattice access
template <typename T, template <typename U> class Descriptor>
struct multiPhaseTemplates2D {
    static void shanChenInteraction(
        BlockLattice2D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> &rhoContribution,
        plint iX, plint iY)
    {
        enum { densityOffset = Descriptor<T>::ExternalField::densityBeginsAt };

        rhoContribution.resetToZero();
        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            plint nextX = iX + Descriptor<T>::c[iPop][0];
            plint nextY = iY + Descriptor<T>::c[iPop][1];
            Cell<T, Descriptor> const &cell = lattice.get(nextX, nextY);
            T rho = *cell.getExternal(densityOffset);
            for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
                rhoContribution[iD] += Descriptor<T>::t[iPop] * rho * Descriptor<T>::c[iPop][iD];
            }
        }
    }
};

template <typename T, template <typename U> class Descriptor>
struct multiPhaseHelpers2D {
    typedef Descriptor<T> D;

    static void computeRhoAndJ(BlockLattice2D<T, Descriptor> &lattice, plint iX, plint iY)
    {
        enum {
            densityOffset = D::ExternalField::densityBeginsAt,
            momentumOffset = D::ExternalField::momentumBeginsAt
        };

        Cell<T, Descriptor> &cell = lattice.get(iX, iY);
        Array<T, Descriptor<T>::d> j;
        T rhoBar;
        cell.getDynamics().computeRhoBarJ(cell, rhoBar, j);
        momentTemplates<T, Descriptor>::get_j(cell, j);
        *cell.getExternal(densityOffset) = Descriptor<T>::fullRho(rhoBar);
        j.to_cArray(cell.getExternal(momentumOffset));
    }

    static void computeUstar(
        std::vector<BlockLattice2D<T, Descriptor> *> lattices, plint iX, plint iY,
        Array<T, Descriptor<T>::d> &uStar)
    {
        enum {
            densityOffset = D::ExternalField::densityBeginsAt,
            momentumOffset = D::ExternalField::momentumBeginsAt
        };

        plint numSpecies = (plint)lattices.size();
        std::vector<T> omega(numSpecies);
        T rhoStar = T();
        for (plint iSpecies = 0; iSpecies < numSpecies; ++iSpecies) {
            Cell<T, Descriptor> const &cell = lattices[iSpecies]->get(iX, iY);
            omega[iSpecies] = cell.getDynamics().getOmega();
            rhoStar += omega[iSpecies] * (*cell.getExternal(densityOffset));
        }
        // Computation of the common velocity, shared among all populations.
        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            uStar[iD] = T();
            for (plint iSpecies = 0; iSpecies < numSpecies; ++iSpecies) {
                T *momentum = lattices[iSpecies]->get(iX, iY).getExternal(momentumOffset);
                uStar[iD] += momentum[iD] * omega[iSpecies];
            }
            uStar[iD] /= rhoStar;
        }
    }

    static void computeUphys(
        std::vector<BlockLattice2D<T, Descriptor> *> lattices, plint iX, plint iY,
        Array<T, Descriptor<T>::d> &uPhys)
    {
        enum {
            densityOffset = D::ExternalField::densityBeginsAt,
            momentumOffset = D::ExternalField::momentumBeginsAt,
            gOffset = D::ExternalField::GbeginsAt,
            interfOffset = D::ExternalField::interactionFbeginsAt
        };
        plint numSpecies = (plint)lattices.size();
        T rhoPhys = T();
        for (plint iSpecies = 0; iSpecies < numSpecies; ++iSpecies) {
            Cell<T, Descriptor> const &cell = lattices[iSpecies]->get(iX, iY);
            rhoPhys += *cell.getExternal(densityOffset);
        }
        // Computation of the common velocity, shared among all populations.
        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            uPhys[iD] = T();
            for (plint iSpecies = 0; iSpecies < numSpecies; ++iSpecies) {
                Cell<T, Descriptor> const &cell = lattices[iSpecies]->get(iX, iY);
                T const *momentum = cell.getExternal(momentumOffset);
                T const *interf = cell.getExternal(interfOffset);
                T rho = *cell.getExternal(densityOffset);
                uPhys[iD] += momentum[iD] + 0.5 * rho * interf[iD];
            }
            uPhys[iD] /= rhoPhys;
            T *externalf = lattices[0]->get(iX, iY).getExternal(gOffset);
            uPhys[iD] += 0.5 * externalf[iD];
        }
    }
};

template <typename T>
struct multiPhaseTemplates2D<T, descriptors::ForcedShanChenD2Q9Descriptor> {
    typedef descriptors::ForcedShanChenD2Q9Descriptor<T> D;

    static void shanChenInteraction(
        BlockLattice2D<T, descriptors::ForcedShanChenD2Q9Descriptor> &lattice,
        Array<T, D::d> &rhoContribution, plint iX, plint iY)
    {
        enum { densityOffset = D::ExternalField::densityBeginsAt };

        T rho;
        rho = *lattice.get(iX - 1, iY + 1).getExternal(densityOffset);
        rhoContribution[0] = -D::t[1] * rho;
        rhoContribution[1] = D::t[1] * rho;
        rho = *lattice.get(iX - 1, iY).getExternal(densityOffset);
        rhoContribution[0] -= D::t[2] * rho;
        rho = *lattice.get(iX - 1, iY - 1).getExternal(densityOffset);
        rhoContribution[0] -= D::t[3] * rho;
        rhoContribution[1] -= D::t[3] * rho;
        rho = *lattice.get(iX, iY - 1).getExternal(densityOffset);
        rhoContribution[1] -= D::t[4] * rho;
        rho = *lattice.get(iX + 1, iY - 1).getExternal(densityOffset);
        rhoContribution[0] += D::t[5] * rho;
        rhoContribution[1] -= D::t[5] * rho;
        rho = *lattice.get(iX + 1, iY).getExternal(densityOffset);
        rhoContribution[0] += D::t[6] * rho;
        rho = *lattice.get(iX + 1, iY + 1).getExternal(densityOffset);
        rhoContribution[0] += D::t[7] * rho;
        rhoContribution[1] += D::t[7] * rho;
        rho = *lattice.get(iX, iY + 1).getExternal(densityOffset);
        rhoContribution[1] += D::t[8] * rho;
    }
};

}  // namespace plb

#endif  // MULTI_PHASE_TEMPLATES_2D_H
