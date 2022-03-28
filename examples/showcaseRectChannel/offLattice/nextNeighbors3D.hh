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

#ifndef NEXT_NEIGHBORS_3D_HH
#define NEXT_NEIGHBORS_3D_HH

#include <cmath>

#include "latticeBoltzmann/externalFieldAccess.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "offLattice/nextNeighbors3D.h"

namespace plb {

template <typename T>
const int NextNeighbor<T>::numNeighbors = 26;

template <typename T>
const int NextNeighbor<T>::c[26][3] = {
    {0, 0, 1},   {0, 0, -1},  {0, 1, 0},   {0, 1, 1},   {0, 1, -1},
    {0, -1, 0},  {0, -1, 1},  {0, -1, -1},

    {1, 0, 0},   {1, 0, 1},   {1, 0, -1},  {1, 1, 0},   {1, 1, 1},
    {1, 1, -1},  {1, -1, 0},  {1, -1, 1},  {1, -1, -1},

    {-1, 0, 0},  {-1, 0, 1},  {-1, 0, -1}, {-1, 1, 0},  {-1, 1, 1},
    {-1, 1, -1}, {-1, -1, 0}, {-1, -1, 1}, {-1, -1, -1}};

template <typename T>
const T NextNeighbor<T>::d1 = (T)1;
template <typename T>
const T NextNeighbor<T>::d2 = std::sqrt((T)2);
template <typename T>
const T NextNeighbor<T>::d3 = std::sqrt((T)3);

template <typename T>
const T NextNeighbor<T>::d[26] = {
    (T)1,
    (T)1,
    (T)1,
    std::sqrt((T)2),
    std::sqrt((T)2),
    (T)1,
    std::sqrt((T)2),
    std::sqrt((T)2),
    (T)1,
    std::sqrt((T)2),
    std::sqrt((T)2),
    std::sqrt((T)2),
    std::sqrt((T)3),
    std::sqrt((T)3),
    std::sqrt((T)2),
    std::sqrt((T)3),
    std::sqrt((T)3),
    (T)1,
    std::sqrt((T)2),
    std::sqrt((T)2),
    std::sqrt((T)2),
    std::sqrt((T)3),
    std::sqrt((T)3),
    std::sqrt((T)2),
    std::sqrt((T)3),
    std::sqrt((T)3)};

template <typename T>
const T NextNeighbor<T>::id1 = (T)1;
template <typename T>
const T NextNeighbor<T>::id2 = (T)1 / std::sqrt((T)2);
template <typename T>
const T NextNeighbor<T>::id3 = (T)1 / std::sqrt((T)3);

template <typename T>
const T NextNeighbor<T>::invD[26] = {
    (T)1,
    (T)1,
    (T)1,
    (T)1 / std::sqrt((T)2),
    (T)1 / std::sqrt((T)2),
    (T)1,
    (T)1 / std::sqrt((T)2),
    (T)1 / std::sqrt((T)2),
    (T)1,
    (T)1 / std::sqrt((T)2),
    (T)1 / std::sqrt((T)2),
    (T)1 / std::sqrt((T)2),
    (T)1 / std::sqrt((T)3),
    (T)1 / std::sqrt((T)3),
    (T)1 / std::sqrt((T)2),
    (T)1 / std::sqrt((T)3),
    (T)1 / std::sqrt((T)3),
    (T)1,
    (T)1 / std::sqrt((T)2),
    (T)1 / std::sqrt((T)2),
    (T)1 / std::sqrt((T)2),
    (T)1 / std::sqrt((T)3),
    (T)1 / std::sqrt((T)3),
    (T)1 / std::sqrt((T)2),
    (T)1 / std::sqrt((T)3),
    (T)1 / std::sqrt((T)3)};

template <typename T, template <typename U> class Descriptor>
NextNeighborPop<T, Descriptor>::NextNeighborPop()
{
    typedef Descriptor<T> D;
    for (plint iNeighbor = 0; iNeighbor < NextNeighbor<T>::numNeighbors; ++iNeighbor) {
        int const *c = NextNeighbor<T>::c[iNeighbor];
        ids[iNeighbor] = -1;
        for (plint iPop = 0; iPop < D::q; ++iPop) {
            if (D::c[iPop][0] == c[0] && D::c[iPop][1] == c[1] && D::c[iPop][2] == c[2]) {
                ids[iNeighbor] = iPop;
                break;
            }
        }
    }
}

}  // namespace plb

#endif  // NEXT_NEIGHBORS_3D_HH
