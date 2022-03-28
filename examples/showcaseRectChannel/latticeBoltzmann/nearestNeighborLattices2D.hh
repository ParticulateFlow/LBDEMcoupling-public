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
 * Descriptor for nearest-neighbor 2D lattices -- generic code.
 **/
#ifndef NEAREST_NEIGHBOR_LATTICES_2D_HH
#define NEAREST_NEIGHBOR_LATTICES_2D_HH

#include "latticeBoltzmann/nearestNeighborLattices2D.h"

// clang-format off

namespace plb {

namespace descriptors {

    // D2Q9 ////////////////////////////////////////////////////////////

    template<typename T>
    const T D2Q9Constants<T>::invD = (T)1 / (T) d;

    template<typename T>
    const int D2Q9Constants<T>::vicinity = 1;

    template<typename T>
    const int D2Q9Constants<T>::c
        [D2Q9Constants<T>::q][D2Q9Constants<T>::d] =
        {
            { 0, 0},
            {-1, 1}, {-1, 0}, {-1,-1}, { 0,-1},
            { 1,-1}, { 1, 0}, { 1, 1}, { 0, 1}
        };

    template<typename T>
    const int D2Q9Constants<T>::cNormSqr[D2Q9Constants<T>::q] =
        { 0, 2, 1, 2, 1, 2, 1, 2, 1 };

    template<typename T>
    const T D2Q9Constants<T>::t[D2Q9Constants<T>::q] =
        {
            (T)4/(T)9, (T)1/(T)36, (T)1/(T)9, (T)1/(T)36, (T)1/(T)9,
                       (T)1/(T)36, (T)1/(T)9, (T)1/(T)36, (T)1/(T)9
        };

    template<typename T>
    const T D2Q9Constants<T>::cs2 = (T)1 / (T)3;

    template<typename T>
    const T D2Q9Constants<T>::invCs2 = (T)3;

    template<typename T>
    const char D2Q9Descriptor<T>::name[] = "D2Q9";

    template<typename T>
    const char ForcedD2Q9Descriptor<T>::name[] = "ForcedD2Q9";

    template<typename T>
    const char RhoBarJD2Q9Descriptor<T>::name[] = "RhoBarJD2Q9";
    
    template<typename T>
    const char RhoBarVelocityPiNeqOmegaD2Q9Descriptor<T>::name[] = "RhoBarVelocityPiNeqOmegaD2Q9";
    
    template<typename T>
    const char VelocityD2Q9Descriptor<T>::name[] = "VelocityD2Q9";

    template<typename T>
    const char Tau1_D2Q9Descriptor<T>::name[] = "Tau1_D2Q9";

    template<typename T>
    const char AbsorbingWaveD2Q9Descriptor<T>::name[] = "AbsorbingWave_D2Q9";

    template<typename T>
    const char ExtSymTensorD2Q9Descriptor<T>::name[] = "ExtSymTensor_D2Q9";

}  // namespace descriptors

}  // namespace plb

// clang-format on

#endif
