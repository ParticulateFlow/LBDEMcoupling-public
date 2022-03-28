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

#ifndef CREATE_BUBBLES_3D_H
#define CREATE_BUBBLES_3D_H

#include <limits>

#include "atomicBlock/atomicContainerBlock3D.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/reductiveDataProcessingFunctional3D.h"
#include "core/globalDefs.h"
#include "multiPhysics/freeSurfaceModel3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
void punchSphere(
    FreeSurfaceFields3D<T, Descriptor> &fields, Array<T, 3> const &center, T radius, T rhoEmpty,
    T rho0, Dynamics<T, Descriptor> &dynamics);

template <typename T, template <typename U> class Descriptor>
void analyticalPunchSphere(
    FreeSurfaceFields3D<T, Descriptor> &fields, Array<T, 3> const &center, T radius, T rhoEmpty,
    T rho0, plint subDivision, Dynamics<T, Descriptor> &dynamics);

template <typename T, template <typename U> class Descriptor>
T computeAverageSphereDensity(
    FreeSurfaceFields3D<T, Descriptor> &fields, Array<T, 3> const &center, T radius);

template <typename T, template <typename U> class Descriptor>
void punchSphere(
    FreeSurfaceSetup<T, Descriptor> &setup, Array<T, 3> const &center, T radius, T rhoEmpty, T rho0,
    Dynamics<T, Descriptor> &dynamics);

template <typename T, template <typename U> class Descriptor>
void analyticalPunchSphere(
    FreeSurfaceSetup<T, Descriptor> &setup, Array<T, 3> const &center, T radius, T rhoEmpty, T rho0,
    plint subDivision, Dynamics<T, Descriptor> &dynamics);

template <typename T, template <typename U> class Descriptor>
T computeAverageSphereDensity(
    FreeSurfaceSetup<T, Descriptor> &setup, Array<T, 3> const &center, T radius);

template <typename T, template <typename U> class Descriptor>
void punchSphere(
    FreeSurfaceWrapper<T, Descriptor> &wrapper, Array<T, 3> const &center, T radius, T rhoEmpty,
    T rho0, Dynamics<T, Descriptor> &dynamics);

template <typename T, template <typename U> class Descriptor>
void analyticalPunchSphere(
    FreeSurfaceWrapper<T, Descriptor> &wrapper, Array<T, 3> const &center, T radius, T rhoEmpty,
    T rho0, plint subDivision, Dynamics<T, Descriptor> &dynamics);

template <typename T, template <typename U> class Descriptor>
T computeAverageSphereDensity(
    FreeSurfaceWrapper<T, Descriptor> &wrapper, Array<T, 3> const &center, T radius);

}  // namespace plb

#endif  // CREATE_BUBBLES_3D_H
