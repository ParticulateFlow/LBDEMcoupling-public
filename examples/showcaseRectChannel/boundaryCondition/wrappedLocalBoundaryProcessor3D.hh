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

#ifndef WRAPPED_LOCAL_BOUNDARY_PROCESSOR_3D_HH
#define WRAPPED_LOCAL_BOUNDARY_PROCESSOR_3D_HH

#include <typeinfo>

#include "atomicBlock/blockLattice3D.h"
#include "boundaryCondition/wrappedLocalBoundaryProcessor3D.h"
#include "core/processorIdentifiers3D.h"

namespace plb {

///////////  WrappedLocalBoundaryFunctional3D ///////////////////////////////////

template <typename T, template <typename U> class Descriptor>
const int WrappedLocalBoundaryFunctional3D<T, Descriptor>::staticId =
    meta::registerProcessor3D<WrappedLocalBoundaryFunctional3D<T, Descriptor>, T, Descriptor>(
        std::string("WrappedLocalBoundary3D"));

template <typename T, template <typename U> class Descriptor>
void WrappedLocalBoundaryFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    PLB_ASSERT(domain.x0 == domain.x1 || domain.y0 == domain.y1 || domain.z0 == domain.z1);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
                Dynamics<T, Descriptor> &dynamics = cell.getDynamics();
                // Make a safety check to avoid errors in case the user replaces
                //   the original composite dynamics by something else; for
                //   example, when bounce-back is used as an "eraser tool" at
                //   inlets/outlets.
                if (dynamics.isComposite()) {
                    dynamic_cast<CompositeDynamics<T, Descriptor> &>(dynamics).prepareCollision(
                        cell);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
WrappedLocalBoundaryFunctional3D<T, Descriptor>
    *WrappedLocalBoundaryFunctional3D<T, Descriptor>::clone() const
{
    return new WrappedLocalBoundaryFunctional3D<T, Descriptor>(*this);
}

}  // namespace plb

#endif  // WRAPPED_LOCAL_BOUNDARY_PROCESSOR_3D_HH
