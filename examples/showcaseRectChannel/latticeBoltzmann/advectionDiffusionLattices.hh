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
 * Descriptor for all types of 2D and 3D lattices. In principle, thanks
 * to the fact that the Palabos code is generic, it is sufficient to
 * write a new descriptor when a new type of lattice is to be used.
 *  -- generic code
 */
#ifndef ADVECTION_DIFFUSION_LATTICES_HH
#define ADVECTION_DIFFUSION_LATTICES_HH

#include "latticeBoltzmann/advectionDiffusionLattices.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"

namespace plb {

namespace descriptors {

// AdvectionDiffusion D2Q5 //////////////////////////////////////////////

template <typename T>
const T D2Q5Constants<T>::invD = (T)1 / (T)d;

template <typename T>
const int D2Q5Constants<T>::vicinity = 1;

template <typename T>
const int D2Q5Constants<T>::c[D2Q5Constants<T>::q][D2Q5Constants<T>::d] = {
    {0, 0}, {-1, 0}, {0, -1}, {1, 0}, {0, 1}};

template <typename T>
const int D2Q5Constants<T>::cNormSqr[D2Q5Constants<T>::q] = {0, 1, 1, 1, 1};

template <typename T>
const T D2Q5Constants<T>::t[D2Q5Constants<T>::q] = {
    (T)1 - (T)2 / invCs2, (T)1 / (invCs2 * (T)2), (T)1 / (invCs2 * (T)2), (T)1 / (invCs2 * (T)2),
    (T)1 / (invCs2 * (T)2)};

template <typename T>
const T D2Q5Constants<T>::cs2 = (T)1 / (T)3;

template <typename T>
const T D2Q5Constants<T>::invCs2 = (T)3;

template <typename T>
const char AdvectionDiffusionD2Q5Descriptor<T>::name[] = "AdvectionDiffusionD2Q5";

template <typename T>
const char AdvectionDiffusionWithSourceD2Q5Descriptor<T>::name[] =
    "AdvectionDiffusionWithSourceD2Q5";

// AdvectionDiffusion D3Q7 ////////////////////////////////////////////////////

template <typename T>
const T D3Q7Constants<T>::invD = (T)1 / (T)d;

template <typename T>
const int D3Q7Constants<T>::vicinity = 1;

template <typename T>
const int D3Q7Constants<T>::c[D3Q7Constants<T>::q][D3Q7Constants<T>::d] = {
    {0, 0, 0},

    {-1, 0, 0}, {0, -1, 0}, {0, 0, -1}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
};

template <typename T>
const int D3Q7Constants<T>::cNormSqr[D3Q7Constants<T>::q] = {0, 1, 1, 1, 1, 1, 1};

template <typename T>
const T D3Q7Constants<T>::t[D3Q7Constants<T>::q] = {(T)1 - (T)3 / invCs2,

                                                    (T)1 / (invCs2 * (T)2), (T)1 / (invCs2 * (T)2),
                                                    (T)1 / (invCs2 * (T)2), (T)1 / (invCs2 * (T)2),
                                                    (T)1 / (invCs2 * (T)2), (T)1 / (invCs2 * (T)2)};

template <typename T>
const T D3Q7Constants<T>::cs2 = (T)1 / (T)4;

template <typename T>
const T D3Q7Constants<T>::invCs2 = (T)4;

template <typename T>
const char AdvectionDiffusionD3Q7Descriptor<T>::name[] = "AdvectionDiffusionD3Q7";

template <typename T>
const char AdvectionDiffusionWithSourceD3Q7Descriptor<T>::name[] =
    "AdvectionDiffusionWithSourceD3Q7";

// AdvectionDiffusion D3Q19 ////////////////////////////////////////////////////

template <typename T>
const char AdvectionDiffusionD3Q19Descriptor<T>::name[] = "AdvectionDiffusionD3Q19";

template <typename T>
const char AdvectionDiffusionWithSourceD3Q19Descriptor<T>::name[] =
    "AdvectionDiffusionWithSourceD3Q19";

// AdvectionDiffusion D3Q27 ////////////////////////////////////////////////////

template <typename T>
const char AdvectionDiffusionD3Q27Descriptor<T>::name[] = "AdvectionDiffusionD3Q27";

template <typename T>
const char AdvectionDiffusionWithSourceD3Q27Descriptor<T>::name[] =
    "AdvectionDiffusionWithSourceD3Q27";

}  // namespace descriptors

}  // namespace plb

#endif
