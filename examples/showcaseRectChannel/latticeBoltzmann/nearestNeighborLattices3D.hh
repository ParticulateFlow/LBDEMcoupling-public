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
 * Descriptor for nearest-neighbor 3D lattices -- generic code.
 **/
#ifndef NEAREST_NEIGHBOR_LATTICES_3D_HH
#define NEAREST_NEIGHBOR_LATTICES_3D_HH

#include "latticeBoltzmann/nearestNeighborLattices3D.h"

// clang-format off

namespace plb {

namespace descriptors {

    // D3Q13 ///////////////////////////////////////////////////////////

    template<typename T>
    const T D3Q13Constants<T>::invD = (T)1 / (T) d;

    template<typename T>
    const int D3Q13Constants<T>::vicinity = 1;

    template<typename T>
    const int D3Q13Constants<T>::c
        [D3Q13Constants<T>::q][D3Q13Constants<T>::d] =
        {
            { 0, 0, 0},

            {-1,-1, 0}, {-1, 1, 0}, {-1, 0,-1},
            {-1, 0, 1}, { 0,-1,-1}, { 0,-1, 1},

            { 1, 1, 0}, { 1,-1, 0}, { 1, 0, 1},
            { 1, 0,-1}, { 0, 1, 1}, { 0, 1,-1}
        }; 

    template<typename T>
    const int D3Q13Constants<T>::cNormSqr[D3Q13Constants<T>::q] =
    { 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 };

    template<typename T>
    const T D3Q13Constants<T>::t[D3Q13Constants<T>::q] =
        {
            (T)1/(T)2,

            (T)1/(T)24, (T)1/(T)24, (T)1/(T)24,
            (T)1/(T)24, (T)1/(T)24, (T)1/(T)24,

            (T)1/(T)24, (T)1/(T)24, (T)1/(T)24,
            (T)1/(T)24, (T)1/(T)24, (T)1/(T)24
        };

    /** This parameter is chosen to enhance numerical stability */
    template<typename T>
    const T D3Q13Constants<T>::cs2 = (T)1 / (T)3;

    /** This parameter is chosen to enhance numerical stability */
    template<typename T>
    const T D3Q13Constants<T>::invCs2 = (T)3;

    /** This parameter is chosen to enhance numerical stability */
    template<typename T>
    const T D3Q13Constants<T>::lambda_e = (T)1.5;

    /** This parameter is chosen to enhance numerical stability */
    template<typename T>
    const T D3Q13Constants<T>::lambda_h = (T)1.8;

    template<typename T>
    const char D3Q13Descriptor<T>::name[] = "D3Q13";

    template<typename T>
    const char ForcedD3Q13Descriptor<T>::name[] = "ForcedD3Q13";


    // D3Q15 ///////////////////////////////////////////////////////////

    template<typename T>
    const T D3Q15Constants<T>::invD = (T)1 / (T) d;

    template<typename T>
    const int D3Q15Constants<T>::vicinity = 1;

    template<typename T>
    const int D3Q15Constants<T>::c
        [D3Q15Constants<T>::q][D3Q15Constants<T>::d] =
        {
            { 0, 0, 0},

            {-1, 0, 0}, { 0,-1, 0}, { 0, 0,-1},
            {-1,-1,-1}, {-1,-1, 1}, {-1, 1,-1}, {-1, 1, 1},

            { 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1},
            { 1, 1, 1}, { 1, 1,-1}, { 1,-1, 1}, { 1,-1,-1},

        }; 

    template<typename T>
    const int D3Q15Constants<T>::cNormSqr
        [D3Q15Constants<T>::q] =
        { 0,  1, 1, 1,  3, 3, 3, 3, 
              1, 1, 1,  3, 3, 3, 3 };
 
    template<typename T>
    const T D3Q15Constants<T>::t[D3Q15Constants<T>::q] =
        {
            (T)2/(T)9,

            (T)1/(T)9, (T)1/(T)9, (T)1/(T)9, 
            (T)1/(T)72, (T)1/(T)72, (T)1/(T)72, (T)1/(T)72,

            (T)1/(T)9, (T)1/(T)9, (T)1/(T)9, 
            (T)1/(T)72, (T)1/(T)72, (T)1/(T)72, (T)1/(T)72
        };

    template<typename T>
    const T D3Q15Constants<T>::cs2 = (T)1 / (T)3;

    template<typename T>
    const T D3Q15Constants<T>::invCs2 = (T)3;

    template<typename T>
    const char D3Q15Descriptor<T>::name[] = "D3Q15";

    template<typename T>
    const char ForcedD3Q15Descriptor<T>::name[] = "ForcedD3Q15";

    template<typename T>
    const char RhoBarJD3Q15Descriptor<T>::name[] = "RhoBarJD3Q15";

    template<typename T>
    const char RhoBarVelocityPiNeqOmegaD3Q15Descriptor<T>::name[] = "RhoBarVelocityPiNeqOmegaD3Q15";

    template<typename T>
    const char ForcedRhoBarJD3Q15Descriptor<T>::name[] = "ForcedRhoBarJD3Q15";

    template<typename T>
    const char AbsorbingWaveD3Q15Descriptor<T>::name[] = "AbsorbingWave_D3Q15";


    // D3Q19 ///////////////////////////////////////////////////////////

    template<typename T>
    const T D3Q19Constants<T>::invD = (T)1 / (T) d;

    template<typename T>
    const int D3Q19Constants<T>::vicinity = 1;

    template<typename T>
    const int D3Q19Constants<T>::c
        [D3Q19Constants<T>::q][D3Q19Constants<T>::d] =
        {
            { 0, 0, 0},

            {-1, 0, 0}, { 0,-1, 0}, { 0, 0,-1},
            {-1,-1, 0}, {-1, 1, 0}, {-1, 0,-1},
            {-1, 0, 1}, { 0,-1,-1}, { 0,-1, 1},

            { 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1},
            { 1, 1, 0}, { 1,-1, 0}, { 1, 0, 1},
            { 1, 0,-1}, { 0, 1, 1}, { 0, 1,-1}
        }; 

    template<typename T>
    const int D3Q19Constants<T>::cNormSqr
        [D3Q19Constants<T>::q] =
        { 0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 
             1, 1, 1, 2, 2, 2, 2, 2, 2 };
 
    template<typename T>
    const T D3Q19Constants<T>::t[D3Q19Constants<T>::q] =
        {
            (T)1/(T)3,

            (T)1/(T)18, (T)1/(T)18, (T)1/(T)18, 
            (T)1/(T)36, (T)1/(T)36, (T)1/(T)36,
            (T)1/(T)36, (T)1/(T)36, (T)1/(T)36,

            (T)1/(T)18, (T)1/(T)18, (T)1/(T)18, 
            (T)1/(T)36, (T)1/(T)36, (T)1/(T)36,
            (T)1/(T)36, (T)1/(T)36, (T)1/(T)36
        };

    template<typename T>
    const T D3Q19Constants<T>::cs2 = (T)1 / (T)3;

    template<typename T>
    const T D3Q19Constants<T>::invCs2 = (T)3;

    template<typename T>
    const char D3Q19Descriptor<T>::name[] = "D3Q19";

    template<typename T>
    const char ForcedD3Q19Descriptor<T>::name[] = "ForcedD3Q19";

    template<typename T>
    const char RhoBarJD3Q19Descriptor<T>::name[] = "RhoBarJD3Q19";

    template<typename T>
    const char RhoBarVelocityPiNeqOmegaD3Q19Descriptor<T>::name[] = "RhoBarVelocityPiNeqOmegaD3Q19";

    template<typename T>
    const char RhoBarVelocityPiNeqOmegaWithSourceD3Q19Descriptor<T>::name[] = "RhoBarVelocityPiNeqOmegaWithSourceD3Q19";

    template<typename T>
    const char ForcedRhoBarJD3Q19Descriptor<T>::name[] = "ForcedRhoBarJD3Q19";
    
    template<typename T>
    const char VelocityD3Q19Descriptor<T>::name[] = "VelocityD3Q19";
    
    template<typename T>
    const char SymTensorD3Q19Descriptor<T>::name[] = "SymTensorD3Q19";

    template<typename T>
    const char AbsorbingWaveD3Q19Descriptor<T>::name[] = "AbsorbingWave_D3Q19";
    
    // D3Q27 ///////////////////////////////////////////////////////////

    template<typename T>
    const T D3Q27Constants<T>::invD = (T)1 / (T) d;

    template<typename T>
    const int D3Q27Constants<T>::vicinity = 1;

    template<typename T>
    const int D3Q27Constants<T>::c
        [D3Q27Constants<T>::q][D3Q27Constants<T>::d] =
        {
            { 0, 0, 0},

            {-1, 0, 0}, { 0,-1, 0}, { 0, 0,-1},
            {-1,-1, 0}, {-1, 1, 0}, {-1, 0,-1},
            {-1, 0, 1}, { 0,-1,-1}, { 0,-1, 1},
            {-1,-1,-1}, {-1,-1, 1}, {-1, 1,-1}, {-1, 1, 1},

            { 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1},
            { 1, 1, 0}, { 1,-1, 0}, { 1, 0, 1},
            { 1, 0,-1}, { 0, 1, 1}, { 0, 1,-1},
            { 1, 1, 1}, { 1, 1,-1}, { 1,-1, 1}, { 1,-1,-1}
        }; 
 

    template<typename T>
    const int D3Q27Constants<T>::cNormSqr[D3Q27Constants<T>::q] =
    {
        0,
        1, 1, 1,  2, 2, 2,  2, 2, 2,  3, 3, 3, 3,
        1, 1, 1,  2, 2, 2,  2, 2, 2,  3, 3, 3, 3
    };

    template<typename T>
    const T D3Q27Constants<T>::t[D3Q27Constants<T>::q] =
        {
            (T)8/(T)27,

            (T)2/(T)27, (T)2/(T)27, (T)2/(T)27, 
            (T)1/(T)54, (T)1/(T)54, (T)1/(T)54,
            (T)1/(T)54, (T)1/(T)54, (T)1/(T)54,
            (T)1/(T)216, (T)1/(T)216, (T)1/(T)216, (T)1/(T)216,

            (T)2/(T)27, (T)2/(T)27, (T)2/(T)27, 
            (T)1/(T)54, (T)1/(T)54, (T)1/(T)54,
            (T)1/(T)54, (T)1/(T)54, (T)1/(T)54,
            (T)1/(T)216, (T)1/(T)216, (T)1/(T)216, (T)1/(T)216
        };

    template<typename T>
    const T D3Q27Constants<T>::cs2 = (T)1 / (T)3;

    template<typename T>
    const T D3Q27Constants<T>::invCs2 = (T)3;

    template<typename T>
    const char D3Q27Descriptor<T>::name[] = "D3Q27";

    template<typename T>
    const char ForcedD3Q27Descriptor<T>::name[] = "ForcedD3Q27";

    template<typename T>
    const char RhoBarJD3Q27Descriptor<T>::name[] = "RhoBarJD3Q27";

    template<typename T>
    const char ForcedRhoBarJD3Q27Descriptor<T>::name[] = "ForcedRhoBarJD3Q27";
    
    template<typename T>
    const char AbsorbingWaveD3Q27Descriptor<T>::name[] = "AbsorbingWave_D3Q27";

    

}  // namespace descriptors

}  // namespace plb

// clang-format on

#endif
