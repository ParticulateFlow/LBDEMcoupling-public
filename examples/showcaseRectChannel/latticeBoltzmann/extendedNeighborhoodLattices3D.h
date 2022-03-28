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
 * Descriptor for nearest-neighbor 3D lattices. In principle, thanks
 * to static genericity of the code, definition of the lattice and
 * implementation of the dynamics are independent. Of course, there
 * remain a few dependencies. For example, a dynamics which uses higher-
 * order moments such as heat flux needs a lattice with an extended
 * neighborhood. Another example is the D3Q13 lattice, which works only
 * with its very own dynamics.
 *
 *  -- header file
 */
#ifndef EXTENDED_NEIGHBORHOOD_LATTICES_3D_H
#define EXTENDED_NEIGHBORHOOD_LATTICES_3D_H

#include <vector>

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "latticeBoltzmann/externalFields.h"
#include "latticeBoltzmann/roundOffPolicy.h"

namespace plb {

/// Descriptors for the 2D and 3D lattices.
/** \warning Attention: The lattice directions must always be ordered in
 * such a way that c[i] = -c[i+(q-1)/2] for i=1..(q-1)/2, and c[0] = 0 must
 * be the rest velocity. Furthermore, the velocities c[i] for i=1..(q-1)/2
 * must verify
 *  - in 2D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *  - in 3D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *                       || (c[i][0]==0 && c[i][1]==0 && c[i][2]<0)
 * Otherwise some of the code will work erroneously, because the
 * aformentioned relations are taken as given to enable a few
 * optimizations.
 *
 * Givan that this constraint is difficult to verify by hand, we copy-paste
 * a Python code at the end of extendedNeighborhoodLattices3D.hh which does the work.
 */
namespace descriptors {

/// d3q121 lattice constants
template <typename T>
struct D3Q121Constants {
    enum { d = 3, q = 121 };       ///< number of dimensions/distr. functions
    static const T invD;           ///< 1 / (number of dimensions)
    static const int vicinity;     ///< size of neighborhood
    static const int c[q][d];      ///< lattice directions
    static const int cNormSqr[q];  ///< norm-square of the vector c
    static const T t[q];           ///< lattice weights
    static const T cs2;     ///< lattice constant cs2 (in BGK, this is the square-speed-of-sound)
    static const T invCs2;  ///< 1 / cs2
};

template <typename T>
struct D3Q121DescriptorBase : public D3Q121Constants<T>, public DefaultRoundOffPolicy<T> {
    typedef D3Q121DescriptorBase<T> BaseDescriptor;
    enum { numPop = D3Q121Constants<T>::q };
};

template <typename T>
struct D3Q121Descriptor : public D3Q121DescriptorBase<T>, public NoExternalFieldBase {
    static const char name[];
};

template <typename T>
struct ForcedD3Q121Descriptor : public D3Q121DescriptorBase<T>, public Force3dDescriptorBase {
    static const char name[];
};

/// d3q39 lattice constants
template <typename T>
struct D3Q39Constants {
    enum { d = 3, q = 39 };        ///< number of dimensions/distr. functions
    static const T invD;           ///< 1 / (number of dimensions)
    static const int vicinity;     ///< size of neighborhood
    static const int c[q][d];      ///< lattice directions
    static const int cNormSqr[q];  ///< norm-square of the vector c
    static const T t[q];           ///< lattice weights
    static const T cs2;     ///< lattice constant cs2 (in BGK, this is the square-speed-of-sound)
    static const T invCs2;  ///< 1 / cs2
};

template <typename T>
struct D3Q39DescriptorBase : public D3Q39Constants<T>, public DefaultRoundOffPolicy<T> {
    typedef D3Q39DescriptorBase<T> BaseDescriptor;
    enum { numPop = D3Q39Constants<T>::q };
};

template <typename T>
struct D3Q39Descriptor : public D3Q39DescriptorBase<T>, public NoExternalFieldBase {
    static const char name[];
};

template <typename T>
struct ForcedD3Q39Descriptor : public D3Q39DescriptorBase<T>, public Force3dDescriptorBase {
    static const char name[];
};

}  // namespace descriptors

}  // namespace plb

#endif
