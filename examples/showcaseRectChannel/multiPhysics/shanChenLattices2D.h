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
 *  Lattice descriptors for 2D Shan/Chen multi-component flows -- header file
 */
#ifndef SHAN_CHEN_LATTICES_2D_H
#define SHAN_CHEN_LATTICES_2D_H

#include "core/globalDefs.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"

namespace plb {

namespace descriptors {

/// Density and Momentum as external scalars
struct DensityMomentumNoForceExternals2D {
    static const int numScalars = 3;
    static const int numSpecies = 2;

    static const int densityBeginsAt = 0;
    static const int sizeOfDensity = 1;

    static const int momentumBeginsAt = 1;
    static const int sizeOfMomentum = 2;

    static const int forceBeginsAt = 0;
    static const int sizeOfForce = 0;
};

/// Density, Momentum and Force as external scalars
struct DensityMomentumForceExternals2D {
    static const int numScalars = 5;
    static const int numSpecies = 3;

    static const int densityBeginsAt = 0;
    static const int sizeOfDensity = 1;

    static const int momentumBeginsAt = 1;
    static const int sizeOfMomentum = 2;

    static const int forceBeginsAt = 3;
    static const int sizeOfForce = 2;
};

struct ShanChenExternalBase2D {
    typedef DensityMomentumNoForceExternals2D ExternalField;
};

struct ForcedShanChenExternalBase2D {
    typedef DensityMomentumForceExternals2D ExternalField;
};

/// D2Q9 lattice for Shan-Chen model
template <typename T>
struct ShanChenD2Q9Descriptor : public D2Q9DescriptorBase<T>, public ShanChenExternalBase2D {
    static const char name[];
};

template <typename T>
const char ShanChenD2Q9Descriptor<T>::name[] = "ShanChenD2Q9";

/// D2Q9 lattice for Shan-Chen model with force
template <typename T>
struct ForcedShanChenD2Q9Descriptor :
    public D2Q9DescriptorBase<T>,
    public ForcedShanChenExternalBase2D {
    static const char name[];
};

template <typename T>
const char ForcedShanChenD2Q9Descriptor<T>::name[] = "ForcedShanChenD2Q9";

}  // namespace descriptors

}  // namespace plb

#endif  // SHAN_CHEN_LATTICES_2D_H
