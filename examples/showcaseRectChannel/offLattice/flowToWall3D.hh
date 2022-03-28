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

#ifndef FLOW_TO_WALL_3D_HH
#define FLOW_TO_WALL_3D_HH

#include <algorithm>

#include "core/globalDefs.h"
#include "core/util.h"
#include "offLattice/flowToWall3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor, class ValueT>
ValueT optimalInterpolate(
    Lattice<T, Descriptor> &lattice, T x, T y, T z,
    ExtrapolateValue3D<T, Descriptor, ValueT> const &extrapValue, bool &success, float &quality)
{
    Array<ValueT, 4> values;
    Array<int, 4> ok;
    values[0] = extrapValue(lattice, -2, x, y, z, ok[0]);
    values[1] = extrapValue(lattice, -1, x, y, z, ok[1]);
    values[2] = extrapValue(lattice, +1, x, y, z, ok[2]);
    values[3] = extrapValue(lattice, +2, x, y, z, ok[3]);

    ValueT value = ValueT();
    quality = 0;              // By default, success failed.
    if (ok[1]) {              //  ? x ? ?
        if (ok[2]) {          //  ? x x ?
            if (ok[3]) {      //  ? x x x
                if (ok[0]) {  //  x x x x ####
                    // O(4) interpolate: excellent quality.
                    quality = 4;
                } else {  //  . x x x ####
                    // O(3) interpolate: good quality.
                    quality = 3;
                }
            } else {          //  ? x x .
                if (ok[0]) {  //  x x x . ####
                    // O(3) interpolate: good quality.
                    quality = 3;
                } else {  //  . x x . ####
                    // O(2) interpolate: fair quality.
                    quality = 2;
                }
            }
        } else {          //  ? x . ?
            if (ok[0]) {  //  x x . ? ####
                // Extrapolate: bad quality.
                quality = 1;
            }
        }
    } else {                   //  ? . ? ?
        if (ok[2] && ok[3]) {  //  ? . x x ####
            // Extrapolate: bad quality.
            quality = 1;
        }
    }
    return value;
}

template <typename T, template <typename U> class Descriptor, class ValueT>
ExtrapolateValue3D<T, Descriptor, ValueT>::ExtrapolateValue3D(
    ComputeValue3D<T, Descriptor, ValueT> *computeValue_) :
    computeValue(computeValue_)
{ }

template <typename T, template <typename U> class Descriptor, class ValueT>
ExtrapolateValue3D<T, Descriptor, ValueT>::ExtrapolateValue3D(
    ExtrapolateValue3D<T, Descriptor, ValueT> const &rhs) computeValue(rhs.computeValue->clone())
{ }

template <typename T, template <typename U> class Descriptor, class ValueT>
ExtrapolateValue3D<T, Descriptor, ValueT> &ExtrapolateValue3D<T, Descriptor, ValueT>::operator=(
    ExtrapolateValue3D<T, Descriptor, ValueT> const &rhs);
{
    ExtrapolateValue3D<T, Descriptor, ValueT>(rhs).swap(*this);
    return *this;
}

template <typename T, template <typename U> class Descriptor, class ValueT>
void ExtrapolateValue3D<T, Descriptor, ValueT>::swap(ExtrapolateValue3D<T, Descriptor, ValueT> &rhs)
{
    std::swap(computeValue, rhs.computeValue);
}

template <typename T, template <typename U> class Descriptor, class ValueT>
ExtrapolateValue3D<T, Descriptor, ValueT>::~ExtrapolateValue3D()
{
    delete computeValue;
}

template <typename T, template <typename U> class Descriptor, class ValueT>
TerminalExtrapolateValue3D::TerminalExtrapolateValue3D(
    ComputeValue3D<T, Descriptor, ValueT> *computeValue_) :
    ExtrapolateValue3D<T, Descriptor, ValueT>(computeValue_)
{ }

template <typename T, template <typename U> class Descriptor, class ValueT>
std::vector<T> TerminalExtrapolateValue3D::operator()(
    BlockLattice3D<T, Descriptor> const &lattice, plint increment, T x, T y, T z) const
{
    PLB_ASSERT(increment == 0);
    plint iX = util::roundToInt(x);
    plint iY = util::roundToInt(y);
    plint iZ = util::roundToInt(z);
    this->compute(lattice, iX, iY, iZ);
}

template <typename T, template <typename U> class Descriptor, class ValueT>
TerminalExtrapolateValue3D<T, Descriptor, ValueT>
    *TerminalExtrapolateValue3D<T, Descriptor, ValueT>::clone() const
{
    return new TerminalExtrapolateValue3D<T, Descriptor, ValueT>(*this);
}

template <typename T, template <typename U> class Descriptor, class ValueT>
X_ExtrapolateValue3D<T, Descriptor, ValueT>::X_ExtrapolateValue3D(
    ComputeValue3D<T, Descriptor, ValueT> *computeValue_,
    X_ExtrapolateValue3D<T, Descriptor, ValueT> *recurseValue_) :
    ExtrapolateValue3D<T, Descriptor, ValueT>(computeValue_), recurseValue(recurseValue_)
{ }

template <typename T, template <typename U> class Descriptor, class ValueT>
X_ExtrapolateValue3D<T, Descriptor, ValueT>::X_ExtrapolateValue3D(
    X_ExtrapolateValue3D<T, Descriptor, ValueT> const &rhs) :
    ExtrapolateValue3D<T, Descriptor, ValueT>(rhs), recurseValue(rhs.recurseValue->clone())
{ }

template <typename T, template <typename U> class Descriptor, class ValueT>
X_ExtrapolateValue3D<T, Descriptor, ValueT> &X_ExtrapolateValue3D<T, Descriptor, ValueT>::operator=(
    X_ExtrapolateValue3D<T, Descriptor, ValueT> const &rhs)
{
    X_ExtrapolateValue3D<T, Descriptor, ValueT>(rhs).swap(*this);
    return *this;
}

template <typename T, template <typename U> class Descriptor, class ValueT>
void X_ExtrapolateValue3D<T, Descriptor, ValueT>::swap(
    X_ExtrapolateValue3D<T, Descriptor, ValueT> &rhs)
{
    ExtrapolateValue3D<T, Descriptor, ValueT>::swap(rhs);
    std::swap(recurseValue, rhs.recurseValue);
}

template <typename T, template <typename U> class Descriptor, class ValueT>
X_ExtrapolateValue3D<T, Descriptor, ValueT>::~X_ExtrapolateValue3D()
{
    delete recurseValue;
}

template <typename T, template <typename U> class Descriptor, class ValueT>
std::vector<T> X_ExtrapolateValue3D<T, Descriptor, ValueT>::operator()(
    BlockLattice3D<T, Descriptor> const &lattice, plint increment, T x, T y, T z) const
{
    plint iX = util::roundDown(x);
    switch (increment) {
    case -1:
        break;
    case -2:
        iX -= 1;
        break;
    case +1:
        iX += 1;
        break;
    case +2:
        iX += 2;
        break;
defaul:
        PLB_ASSERT(false);
    }
    BLA
}

template <typename T, template <typename U> class Descriptor, class ValueT>
X_ExtrapolateValue3D<T, Descriptor, ValueT> *X_ExtrapolateValue3D<T, Descriptor, ValueT>::clone()
    const
{
    return new X_ExtrapolateValue3D<T, Descriptor, ValueT>(*this);
}

}  // namespace plb

#endif  // FLOW_TO_WALL_3D_HH
