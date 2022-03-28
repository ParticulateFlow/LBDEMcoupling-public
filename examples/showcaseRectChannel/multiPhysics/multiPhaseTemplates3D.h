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
#ifndef MULTI_PHASE_TEMPLATES_3D_H
#define MULTI_PHASE_TEMPLATES_3D_H

#include "core/cell.h"
#include "core/globalDefs.h"
#include "multiPhysics/shanChenLattices3D.h"

namespace plb {

/// Helper functions with full-lattice access
template <typename T, template <typename U> class Descriptor>
struct multiPhaseTemplates3D {
    static void shanChenInteraction(
        BlockLattice3D<T, Descriptor> &lattice, Array<T, Descriptor<T>::d> &rhoContribution,
        plint iX, plint iY, plint iZ)
    {
        enum { densityOffset = Descriptor<T>::ExternalField::densityBeginsAt };

        rhoContribution.resetToZero();
        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            plint nextX = iX + Descriptor<T>::c[iPop][0];
            plint nextY = iY + Descriptor<T>::c[iPop][1];
            plint nextZ = iZ + Descriptor<T>::c[iPop][2];
            Cell<T, Descriptor> const &cell = lattice.get(nextX, nextY, nextZ);
            T rho = *cell.getExternal(densityOffset);
            for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
                rhoContribution[iD] += Descriptor<T>::t[iPop] * rho * Descriptor<T>::c[iPop][iD];
            }
        }
    }

    static void shanChenInteraction(
        BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &rhoBar,
        Array<T, Descriptor<T>::d> &rhoContribution, plint iX, plint iY, plint iZ)
    {
        Dot3D ofs = computeRelativeDisplacement(lattice, rhoBar);
        rhoContribution.resetToZero();
        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            plint nextX = iX + Descriptor<T>::c[iPop][0];
            plint nextY = iY + Descriptor<T>::c[iPop][1];
            plint nextZ = iZ + Descriptor<T>::c[iPop][2];
            T rho = Descriptor<T>::fullRho(rhoBar.get(nextX + ofs.x, nextY + ofs.y, nextZ + ofs.z));
            for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
                rhoContribution[iD] += Descriptor<T>::t[iPop] * rho * Descriptor<T>::c[iPop][iD];
            }
        }
    }
};

template <typename T>
struct multiPhaseTemplates3D<T, descriptors::ForcedShanChenD3Q19Descriptor> {
    typedef descriptors::ForcedShanChenD3Q19Descriptor<T> D;

    static void shanChenInteraction(
        BlockLattice3D<T, descriptors::ForcedShanChenD3Q19Descriptor> &lattice,
        Array<T, D::d> &rhoContribution, plint iX, plint iY, plint iZ)
    {
        enum { densityOffset = D::ExternalField::densityBeginsAt };

        T rho;
        rho = *lattice.get(iX - 1, iY, iZ).getExternal(densityOffset);
        rhoContribution[0] = -D::t[1] * rho;
        rho = *lattice.get(iX, iY - 1, iZ).getExternal(densityOffset);
        rhoContribution[1] = -D::t[2] * rho;
        rho = *lattice.get(iX, iY, iZ - 1).getExternal(densityOffset);
        rhoContribution[2] = -D::t[3] * rho;
        rho = *lattice.get(iX - 1, iY - 1, iZ).getExternal(densityOffset);
        rhoContribution[0] -= D::t[4] * rho;
        rhoContribution[1] -= D::t[4] * rho;
        rho = *lattice.get(iX - 1, iY + 1, iZ).getExternal(densityOffset);
        rhoContribution[0] -= D::t[5] * rho;
        rhoContribution[1] += D::t[5] * rho;
        rho = *lattice.get(iX - 1, iY, iZ - 1).getExternal(densityOffset);
        rhoContribution[0] -= D::t[6] * rho;
        rhoContribution[2] -= D::t[6] * rho;
        rho = *lattice.get(iX - 1, iY, iZ + 1).getExternal(densityOffset);
        rhoContribution[0] -= D::t[7] * rho;
        rhoContribution[2] += D::t[7] * rho;
        rho = *lattice.get(iX, iY - 1, iZ - 1).getExternal(densityOffset);
        rhoContribution[1] -= D::t[8] * rho;
        rhoContribution[2] -= D::t[8] * rho;
        rho = *lattice.get(iX, iY - 1, iZ + 1).getExternal(densityOffset);
        rhoContribution[1] -= D::t[9] * rho;
        rhoContribution[2] += D::t[9] * rho;

        rho = *lattice.get(iX + 1, iY, iZ).getExternal(densityOffset);
        rhoContribution[0] += D::t[10] * rho;
        rho = *lattice.get(iX, iY + 1, iZ).getExternal(densityOffset);
        rhoContribution[1] += D::t[11] * rho;
        rho = *lattice.get(iX, iY, iZ + 1).getExternal(densityOffset);
        rhoContribution[2] += D::t[12] * rho;
        rho = *lattice.get(iX + 1, iY + 1, iZ).getExternal(densityOffset);
        rhoContribution[0] += D::t[13] * rho;
        rhoContribution[1] += D::t[13] * rho;
        rho = *lattice.get(iX + 1, iY - 1, iZ).getExternal(densityOffset);
        rhoContribution[0] += D::t[14] * rho;
        rhoContribution[1] -= D::t[14] * rho;
        rho = *lattice.get(iX + 1, iY, iZ + 1).getExternal(densityOffset);
        rhoContribution[0] += D::t[15] * rho;
        rhoContribution[2] += D::t[15] * rho;
        rho = *lattice.get(iX + 1, iY, iZ - 1).getExternal(densityOffset);
        rhoContribution[0] += D::t[16] * rho;
        rhoContribution[2] -= D::t[16] * rho;
        rho = *lattice.get(iX, iY + 1, iZ + 1).getExternal(densityOffset);
        rhoContribution[1] += D::t[17] * rho;
        rhoContribution[2] += D::t[17] * rho;
        rho = *lattice.get(iX, iY + 1, iZ - 1).getExternal(densityOffset);
        rhoContribution[1] += D::t[18] * rho;
        rhoContribution[2] -= D::t[18] * rho;
    }
};

}  // namespace plb

#endif  // MULTI_PHASE_TEMPLATES_3D_H
