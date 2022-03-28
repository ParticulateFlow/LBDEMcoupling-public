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

/* The original version of this file was written by Orestis Malaspinas
 * and Andrea Parmigiani.
 */

/** \file
 *  Lattice descriptors for 3D Shan/Chen multi-component flows -- header file
 */
#ifndef SHAN_CHEN_LATTICES_3D_H
#define SHAN_CHEN_LATTICES_3D_H

#include "core/globalDefs.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"

namespace plb {

namespace descriptors {

/// Density and Momentum as external scalars
struct DensityMomentumNoForceExternals3D {
    static const int numScalars = 4;
    static const int numSpecies = 2;

    static const int densityBeginsAt = 0;
    static const int sizeOfDensity = 1;

    static const int momentumBeginsAt = 1;
    static const int sizeOfMomentum = 3;

    static const int forceBeginsAt = 0;
    static const int sizeOfForce = 0;
};

/// Density, Momentum and Force as external scalars
struct DensityMomentumForceExternals3D {
    static const int numScalars = 7;
    static const int numSpecies = 3;

    static const int densityBeginsAt = 0;
    static const int sizeOfDensity = 1;

    static const int momentumBeginsAt = 1;
    static const int sizeOfMomentum = 3;

    static const int forceBeginsAt = 4;
    static const int sizeOfForce = 3;
};

struct ShanChenExternalBase3D {
    typedef DensityMomentumNoForceExternals3D ExternalField;
};

struct ForcedShanChenExternalBase3D {
    typedef DensityMomentumForceExternals3D ExternalField;
};

/// D3Q19 lattice for Shan-Chen model
template <typename T>
struct ShanChenD3Q19Descriptor : public D3Q19DescriptorBase<T>, public ShanChenExternalBase3D {
    static const char name[];
};

template <typename T>
const char ShanChenD3Q19Descriptor<T>::name[] = "ShanChenD3Q19";

/// D3Q19 lattice for Shan-Chen model with force
template <typename T>
struct ForcedShanChenD3Q19Descriptor :
    public D3Q19DescriptorBase<T>,
    public ForcedShanChenExternalBase3D {
    static const char name[];
};

template <typename T>
const char ForcedShanChenD3Q19Descriptor<T>::name[] = "ForcedShanChenD3Q19";

}  // namespace descriptors

}  // namespace plb

#endif  // SHAN_CHEN_LATTICES_3D_H
