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

#ifndef FD_STENCILS_2D_H
#define FD_STENCILS_2D_H

#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField2D.h"
#include "atomicBlock/dataField3D.h"
#include "core/globalDefs.h"

namespace plb {

namespace fd {

/// Fourth order central 2d-interpolation (a_p1p1=a(x+1/2,y+1/2), a_m2m2=a(x-3/2,y-3/2)), returns
/// a(x)
template <typename T>
inline T centralSixteenPointsInterp(
    T a_m2m2, T a_m2m1, T a_m2p1, T a_m2p2, T a_m1m2, T a_m1m1, T a_m1p1, T a_m1p2, T a_p1m2,
    T a_p1m1, T a_p1p1, T a_p1p2, T a_p2m2, T a_p2m1, T a_p2p1, T a_p2p2)
{
    //   o --- o --- o --- o // a_m2p2  a_m1p2    a_p1p2  a_p2p2

    //   o --- o --- o --- o // a_m2p1  a_m1p1    a_p1p1  a_p2p1
    //   ---------x--------- //
    //   o --- o --- o --- o // a_m2m1  a_m1m1    a_p1m1  a_p2m1

    //   o --- o --- o --- o // a_m2m2  a_m1m2    a_p1m2  a_p2m2

    T distTwo =
        -(T)9 / (T)256 * (a_m2m1 + a_m2p1 + a_m1m2 + a_m1p2 + a_p1m2 + a_p1p2 + a_p2m1 + a_p2p1);
    T distThree = (T)1 / (T)256 * (a_m2m2 + a_m2p2 + a_p2p2 + a_p2m2);
    T distOne = (T)81 / (T)256 * (a_p1p1 + a_m1p1 + a_m1m1 + a_p1m1);

    return distTwo + distThree + distOne;
}

/// Third order central 2d-interpolation (a_p1p1=a(x+1/2,y+1/2), a_m2p1=a(x-3/2,y+1/2)), returns
/// a(x)
template <typename T>
inline T asymTwelvePointsInterp(
    T a_m2m1, T a_m2p1, T a_m2p2, T a_m1m1, T a_m1p1, T a_m1p2, T a_p1m1, T a_p1p1, T a_p1p2,
    T a_p2m1, T a_p2p1, T a_p2p2)
{
    //   o --- o --- o --- o // a_m2p2  a_m1p2    a_p1p2  a_p2p2

    //   o --- o --- o --- o // a_m2p1  a_m1p1    a_p1p1  a_p2p1
    //   ---------x--------- //
    //   o --- o --- o --- o // a_m2m1  a_m1m1    a_p1m1  a_p2m1

    return (T)1 / (T)128 * (a_m2p2 + a_p2p2) - (T)3 / (T)128 * (a_m2m1 + a_p2m1)
           - (T)3 / (T)64 * (a_m2p1 + a_p2p1) - (T)9 / (T)128 * (a_m1p2 + a_p1p2)
           + (T)27 / (T)128 * (a_m1m1 + a_p1m1) + (T)27 / (T)64 * (a_m1p1 + a_p1p1);
}

/// Third order central 2d-interpolation (a_p1p1=a(x+1/2,y+1/2), a_m2p1=a(x-3/2,y+1/2)), returns
/// a(x)
template <typename T>
inline T asymNinePointsInterp(
    T a_m1m1, T a_m1p1, T a_m1p2, T a_p1m1, T a_p1p1, T a_p1p2, T a_p2m1, T a_p2p1, T a_p2p2)
{
    //   o --- o --- o // a_m1p2    a_p1p2  a_p2p2

    //   o --- o --- o // a_m1p1    a_p1p1  a_p2p1
    //   ---x--------- //
    //   o --- o --- o // a_m1m1    a_p1m1  a_p2m1

    return (T)9 / (T)64 * a_m1m1 + (T)9 / (T)32 * (a_m1p1 + a_p1m1)
           - (T)3 / (T)64 * (a_m1p2 + a_p2m1) + (T)9 / (T)16 * a_p1p1
           - (T)3 / (T)32 * (a_p1p2 + a_p2p1) + (T)1 / (T)64 * a_p2p2;
}

}  // namespace fd

}  // namespace plb

#endif  // FD_STENCILS_2D_H
