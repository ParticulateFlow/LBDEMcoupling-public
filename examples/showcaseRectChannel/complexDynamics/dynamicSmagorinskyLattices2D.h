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
 * Descriptor for nearest-neighbor 2D lattices. In principle, thanks
 * to static genericity of the code, definition of the lattice and
 * implementation of the dynamics are independent. Of course, there
 * remain a few dependencies. For example, a dynamics which uses higher-
 * order moments such as heat flux needs a lattice with an extended
 * neighborhood. Another example is the D3Q13 lattice, which works only
 * with its very own dynamics.
 *
 * Lattice descriptors for dynamic Smagorinsky models -- header file
 */
#ifndef DYNAMIC_SMAGORINSKY_LATTICES_2D_H
#define DYNAMIC_SMAGORINSKY_LATTICES_2D_H

#include "core/globalDefs.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"

namespace plb {

namespace descriptors {

struct ExternalOmegaDescriptor2d {
    static const int numScalars = 1;
    static const int numSpecies = 1;
    static const int omegaBeginsAt = 0;
    static const int sizeOfOmega = 1;
};

struct ExternalOmegaAndForceDescriptor2d {
    static const int numScalars = 3;
    static const int numSpecies = 2;
    static const int omegaBeginsAt = 0;
    static const int sizeOfOmega = 1;
    static const int forceBeginsAt = 1;
    static const int sizeOfForce = 2;
};

template <typename T>
struct ExternalOmegaD2Q9Descriptor :
    public D2Q9DescriptorBase<T>,
    public ExternalOmegaDescriptor2d {
    static const char name[];
};

template <typename T>
struct ForcedExternalOmegaD2Q9Descriptor :
    public D2Q9DescriptorBase<T>,
    public ExternalOmegaAndForceDescriptor2d {
    static const char name[];
};

template <typename T>
const char ExternalOmegaD2Q9Descriptor<T>::name[] = "ExternalOmegaD2Q9";

template <typename T>
const char ExternalOmegaAndForceD2Q9Descriptor<T>::name[] = "ExternalOmegaAndForceD2Q9";

}  // namespace descriptors

}  // namespace plb

#endif  // DYNAMIC_SMAGORINSKY_LATTICES_2D_H
