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
 * Neumann boundary conditions for temperature -- generic implementation.
 */
#ifndef ADIABATIC_BOUNDARY_PROCESSOR_2D_HH
#define ADIABATIC_BOUNDARY_PROCESSOR_2D_HH

#include "complexDynamics/adiabaticBoundaryProcessor2D.h"
#include "finiteDifference/fdStencils1D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void FlatAdiabaticBoundaryFunctional2D<T, Descriptor, direction, orientation>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint iX_prev = iX + ((direction == 0) ? (-orientation) : 0);
            plint iY_prev = iY + ((direction == 1) ? (-orientation) : 0);

            T temperature_1 = lattice.get(iX_prev, iY_prev).computeDensity();

            lattice.get(iX, iY).defineDensity(temperature_1);
        }
    }
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
FlatAdiabaticBoundaryFunctional2D<T, Descriptor, direction, orientation>
    *FlatAdiabaticBoundaryFunctional2D<T, Descriptor, direction, orientation>::clone() const
{
    return new FlatAdiabaticBoundaryFunctional2D<T, Descriptor, direction, orientation>(*this);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void FlatAdiabaticBoundaryFunctional2D<T, Descriptor, direction, orientation>::
    getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
BlockDomain::DomainT
    FlatAdiabaticBoundaryFunctional2D<T, Descriptor, direction, orientation>::appliesTo() const
{
    return BlockDomain::bulk;
}

}  // namespace plb

#endif  // NEUMANN_CONDITION_2D_HH
