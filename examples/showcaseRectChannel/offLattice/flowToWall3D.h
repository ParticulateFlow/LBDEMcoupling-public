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

#ifndef FLOW_TO_WALL_3D_H
#define FLOW_TO_WALL_3D_H

#include "atomicBlock/blockLattice3D.h"
#include "core/globalDefs.h"

namespace plb {

template <typename T, template <typename U> class Descriptor, class ValueT>
class ComputeValue3D {
public:
    virtual ~ComputeValue3D() { }
    virtual ValueT operator()(
        BlockLattice3D<T, Descriptor> const &lattice, T x, T y, T z) const = 0;
    virtual ComputeValue3D<T, Descriptor, ValueT> *clone() const = 0;
};

template <typename T, template <typename U> class Descriptor, class ValueT>
class ExtrapolateValue3D {
public:
    ExtrapolateValue3D(ComputeValue3D<T, Descriptor, ValueT> *computeValue_);
    ExtrapolateValue3D(ExtrapolateValue3D<T, Descriptor, ValueT> const &rhs);
    ExtrapolateValue3D<T, Descriptor, ValueT> &operator=(
        ExtrapolateValue3D<T, Descriptor, ValueT> const &rhs);
    void swap(ExtrapolateValue3D<T, Descriptor, ValueT> &rhs);
    virtual ~ExtrapolateValue3D();
    virtual ValueT operator()(
        BlockLattice3D<T, Descriptor> const &lattice, plint increment, T x, T y, T z) const = 0;
    virtual ExtrapolateValue3D<T, Descriptor, ValueT> *clone() const = 0;
    void compute(BlockLattice3D<T, Descriptor> const &lattice, plint iX, plint iY, plint iZ) const
    {
        computeValue->compute(lattice, iX, iY, iZ);
    }

private:
    ComputeValue3D<T, Descriptor, ValueT> *computeValue;
};

template <typename T, template <typename U> class Descriptor, class ValueT>
ValueT optimalInterpolate(
    Lattice<T, Descriptor> &lattice, T x, T y, T z,
    ExtrapolateValue3D<T, Descriptor, ValueT> const &extrapValue);

template <typename T, template <typename U> class Descriptor, class ValueT>
class TerminalExtrapolateValue3D {
public:
    TerminalExtrapolateValue3D(ComputeValue3D<T, Descriptor, ValueT> *computeValue_);
    virtual ValueT operator()(
        BlockLattice3D<T, Descriptor> const &lattice, plint increment, T x, T y, T z) const;
    virtual TerminalExtrapolateValue3D<T, Descriptor, ValueT> *clone() const;
};

template <typename T, template <typename U> class Descriptor, class ValueT>
class X_ExtrapolateValue3D {
public:
    X_ExtrapolateValue3D(
        ComputeValue3D<T, Descriptor, ValueT> *computeValue_,
        X_ExtrapolateValue3D<T, Descriptor, ValueT> *recurseValue_ =
            new TerminalExtrapolateValue3D<T, Descriptor, ValueT>);
    X_ExtrapolateValue3D(X_ExtrapolateValue3D<T, Descriptor, ValueT> const &rhs);
    X_ExtrapolateValue3D<T, Descriptor, ValueT> &operator=(
        X_ExtrapolateValue3D<T, Descriptor, ValueT> const &rhs);
    void swap(X_ExtrapolateValue3D<T, Descriptor, ValueT> &rhs);
    virtual ~X_ExtrapolateValue3D();
    virtual ValueT operator()(
        BlockLattice3D<T, Descriptor> const &lattice, plint increment, T x, T y, T z) const;
    virtual X_ExtrapolateValue3D<T, Descriptor, ValueT> *clone() const;

private:
    ExtrapolateValue3D<T, Descriptor, ValueT> *recurseValue;
};

}  // namespace plb

#endif  // FLOW_TO_WALL_3D_H
