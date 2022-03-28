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

#ifndef ADVECTION_DIFFUSION_PROCESSOR_3D_HH
#define ADVECTION_DIFFUSION_PROCESSOR_3D_HH

#include "complexDynamics/advectionDiffusionProcessor3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
void CompleteAdvectionDiffusionEdgeBoundaryFunctional3D<T, Descriptor, plane, normal1, normal2>::
    process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    PLB_ASSERT(
        (plane == 2 && domain.x0 == domain.x1 && domain.y0 == domain.y1)
        || (plane == 1 && domain.x0 == domain.x1 && domain.z0 == domain.z1)
        || (plane == 0 && domain.y0 == domain.y1 && domain.z0 == domain.z1));

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                processCell(iX, iY, iZ, lattice);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
void CompleteAdvectionDiffusionEdgeBoundaryFunctional3D<T, Descriptor, plane, normal1, normal2>::
    processCell(plint iX, plint iY, plint iZ, BlockLattice3D<T, Descriptor> &lattice)
{
    Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
    T phi = cell.computeDensity();
    Dynamics<T, Descriptor> &dyn = cell.getDynamics();

    Array<T, 3> d_phi;
    interpolateGradients<plane, 0>(lattice, d_phi[plane], iX, iY, iZ);
    interpolateGradients<direction1, normal1>(lattice, d_phi[direction1], iX, iY, iZ);
    interpolateGradients<direction2, normal2>(lattice, d_phi[direction2], iX, iY, iZ);

    T rho = Descriptor<T>::fullRho(*cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt));
    Array<T, Descriptor<T>::d> jEq;
    jEq.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));
    jEq *= phi;

    Array<T, Descriptor<T>::d> jOne = -Descriptor<T>::cs2 * d_phi * rho / dyn.getOmega();

    Array<T, SymmetricTensor<T, Descriptor>::n> dummyPiNeq;
    dyn.regularize(cell, Descriptor<T>::rhoBar(rho * phi), jEq + jOne, T(), dummyPiNeq);
}

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
template <int deriveDirection, int orientation>
void CompleteAdvectionDiffusionEdgeBoundaryFunctional3D<T, Descriptor, plane, normal1, normal2>::
    interpolateGradients(
        BlockLattice3D<T, Descriptor> const &lattice, T &phiDeriv, plint iX, plint iY, plint iZ)
{
    fd::DirectedGradients3D<
        T, Descriptor, deriveDirection, orientation, deriveDirection,
        deriveDirection != plane>::o1_densityDerivative(phiDeriv, lattice, iX, iY, iZ);
}

// This data processor uses asymmetric finite differences to compute a gradient.
template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
void CompleteAdvectionDiffusionCornerBoundaryFunctional3D<
    T, Descriptor, xNormal, yNormal,
    zNormal>::process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
                Dynamics<T, Descriptor> &dyn = cell.getDynamics();
                Cell<T, Descriptor> &cell100 =
                    lattice.get(iX - 1 * xNormal, iY - 0 * yNormal, iZ - 0 * zNormal);
                Cell<T, Descriptor> &cell010 =
                    lattice.get(iX - 0 * xNormal, iY - 1 * yNormal, iZ - 0 * zNormal);
                Cell<T, Descriptor> &cell001 =
                    lattice.get(iX - 0 * xNormal, iY - 0 * yNormal, iZ - 1 * zNormal);

                T phi100 = cell100.computeDensity();
                T phi010 = cell010.computeDensity();
                T phi001 = cell001.computeDensity();
                T phi = cell.computeDensity();

                T dPhiDx = -xNormal * (phi100 - phi);
                T dPhiDy = -yNormal * (phi010 - phi);
                T dPhiDz = -zNormal * (phi001 - phi);

                T rho = Descriptor<T>::fullRho(
                    *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt));
                Array<T, Descriptor<T>::d> jEq;
                jEq.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));
                jEq *= phi;

                Array<T, Descriptor<T>::d> jOne(
                    -Descriptor<T>::cs2 * dPhiDx * rho / dyn.getOmega(),
                    -Descriptor<T>::cs2 * dPhiDy * rho / dyn.getOmega(),
                    -Descriptor<T>::cs2 * dPhiDz * rho / dyn.getOmega());

                Array<T, SymmetricTensor<T, Descriptor>::n> dummyPiNeq;
                dyn.regularize(cell, Descriptor<T>::rhoBar(rho * phi), jEq + jOne, T(), dummyPiNeq);
            }
        }
    }
}

// This data processor uses symmetric finite differences to compute a gradient.
// It cannot be applied on any part of the boundary of the global simulation
// domain, except if this boundary is periodic.
template <typename T, template <typename U> class Descriptor>
void SetEffectiveDiffusivity3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &rhoBarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, rhoBarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint iX2 = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint iY2 = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint iZ2 = iZ + offset.z;
                T gradX =
                    0.5 * (rhoBarField.get(iX2 + 1, iY2, iZ2) - rhoBarField.get(iX2 - 1, iY2, iZ2));
                T gradY =
                    0.5 * (rhoBarField.get(iX2, iY2 + 1, iZ2) - rhoBarField.get(iX2, iY2 - 1, iZ2));
                T gradZ =
                    0.5 * (rhoBarField.get(iX2, iY2, iZ2 + 1) - rhoBarField.get(iX2, iY2, iZ2 - 1));
                T normGradT = std::sqrt(gradX * gradX + gradY * gradY + gradZ * gradZ);

                // Model: d = d0 ( 1 + C^2 h |gradT|/T0 )
                T omega = omega0 / (1. + (1. - 0.5 * omega0) * cSmagoSqr * invT0 * normGradT);

                lattice.get(iX, iY, iZ).getDynamics().setOmega(omega);
            }
        }
    }
}

}  // namespace plb

#endif  // ADVECTION_DIFFUSION_PROCESSOR_3D_HH
