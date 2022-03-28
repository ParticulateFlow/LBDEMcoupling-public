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

#ifndef ADVECTION_DIFFUSION_PROCESSOR_2D_HH
#define ADVECTION_DIFFUSION_PROCESSOR_2D_HH

#include "complexDynamics/advectionDiffusionProcessor2D.h"

namespace plb {

// This data processor uses asymmetric finite differences to compute a gradient.
template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal>
void AdvectionDiffusionCornerBoundaryFunctional2D<T, Descriptor, xNormal, yNormal>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Cell<T, Descriptor> &cell = lattice.get(iX, iY);
            Dynamics<T, Descriptor> &dyn = cell.getDynamics();

            T phi10 = lattice.get(iX - 1 * xNormal, iY - 0 * yNormal).computeDensity();
            T phi01 = lattice.get(iX - 0 * xNormal, iY - 1 * yNormal).computeDensity();
            T phi = cell.computeDensity();

            T dPhiDx = -xNormal * (phi10 - phi);
            T dPhiDy = -yNormal * (phi01 - phi);

            Array<T, Descriptor<T>::d> jEq;
            jEq.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));
            jEq *= phi;

            Array<T, Descriptor<T>::d> jOne(
                -Descriptor<T>::cs2 * dPhiDx / dyn.getOmega(),
                -Descriptor<T>::cs2 * dPhiDy / dyn.getOmega());

            Array<T, SymmetricTensor<T, Descriptor>::n> dummyPiNeq;
            dyn.regularize(cell, Descriptor<T>::rhoBar(phi), jEq + jOne, T(), dummyPiNeq);
        }
    }
}

// This data processor uses asymmetric finite differences to compute a gradient.
template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal>
void CompleteAdvectionDiffusionCornerBoundaryFunctional2D<T, Descriptor, xNormal, yNormal>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Cell<T, Descriptor> &cell = lattice.get(iX, iY);
            Dynamics<T, Descriptor> &dyn = cell.getDynamics();

            T phi10 = lattice.get(iX - 1 * xNormal, iY - 0 * yNormal).computeDensity();
            T phi01 = lattice.get(iX - 0 * xNormal, iY - 1 * yNormal).computeDensity();
            T phi = cell.computeDensity();

            T dPhiDx = -xNormal * (phi10 - phi);
            T dPhiDy = -yNormal * (phi01 - phi);

            T rho = Descriptor<T>::fullRho(
                *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt));
            Array<T, Descriptor<T>::d> jEq;
            jEq.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));
            jEq *= phi;

            Array<T, Descriptor<T>::d> jOne(
                -Descriptor<T>::cs2 * dPhiDx * rho / dyn.getOmega(),
                -Descriptor<T>::cs2 * dPhiDy * rho / dyn.getOmega());

            Array<T, SymmetricTensor<T, Descriptor>::n> dummyPiNeq;
            dyn.regularize(cell, Descriptor<T>::rhoBar(rho * phi), jEq + jOne, T(), dummyPiNeq);
        }
    }
}

}  // namespace plb

#endif  // ADVECTION_DIFFUSION_PROCESSOR_2D_HH
