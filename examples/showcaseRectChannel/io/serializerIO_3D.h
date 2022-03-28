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

#ifndef SERIALIZER_IO_3D_H
#define SERIALIZER_IO_3D_H

#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "core/block3D.h"
#include "core/globalDefs.h"
#include "io/parallelIO.h"
#include "io/serializerIO.h"
#include "multiBlock/multiBlock3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"

namespace plb {

/// Save the content of a Block3D into a Base64 encoded binary file.
/** The content includes external scalars in the case of a BlockLattice3D. Only raw
 *  data is written, and structural information such as length and width of the
 *  block is lost.
 *
 *  Index-ordering for this operation can be chosen through a call to
 *  global::IOpolicy().setIndexOrderingForStreams(IndexOrdering::OrderingT).
 */
void saveBinaryBlock(Block3D const &block, std::string fName, bool enforceUint = false);

/// Load the content of a Block3D from a Base64 encoded binary file.
/** The content includes external scalars in the case of a BlockLattice3D. Only raw
 *  data is written, and structural information such as length and width of the
 *  block is lost.
 *
 *  Index-ordering for this operation can be chosen through a call to
 *  global::IOpolicy().setIndexOrderingForStreams(IndexOrdering::OrderingT).
 */
void loadBinaryBlock(Block3D &block, std::string fName, bool enforceUint = false);

/// Flush the content of a Block3D into a generic C++ stream with space-separated ASCII words.
/** The content includes external scalars in the case of a BlockLattice3D. Only raw
 *  data is written, and structural information such as length and width of the
 *  block is lost.
 *
 *  Index-ordering for this operation can be chosen through a call to
 *  global::IOpolicy().setIndexOrderingForStreams(IndexOrdering::OrderingT).
 *
 *  This file format is not exact and should be used for data post-processing only,
 *  and not for checkpointing.
 */
template <typename T>
std::ostream &block2ostream(std::ostream &ostr, Block3D const &block);

template <typename T>
std::ostream &operator<<(std::ostream &ostr, MultiScalarField3D<T> const &block)
{
    return block2ostream<T>(ostr, block);
}

template <typename T>
std::ostream &operator<<(std::ostream &ostr, MultiNTensorField3D<T> const &block)
{
    return block2ostream<T>(ostr, block);
}

template <typename T, int nDim>
std::ostream &operator<<(std::ostream &ostr, MultiTensorField3D<T, nDim> const &block)
{
    return block2ostream<T>(ostr, block);
}

template <typename T, template <typename U> class Descriptor>
std::ostream &operator<<(std::ostream &ostr, MultiBlockLattice3D<T, Descriptor> const &block)
{
    return block2ostream<T>(ostr, block);
}

template <typename T>
std::ostream &operator<<(std::ostream &ostr, ScalarField3D<T> const &block)
{
    return block2ostream<T>(ostr, block);
}

template <typename T>
std::ostream &operator<<(std::ostream &ostr, NTensorField3D<T> const &block)
{
    return block2ostream<T>(ostr, block);
}

template <typename T, int nDim>
std::ostream &operator<<(std::ostream &ostr, TensorField3D<T, nDim> const &block)
{
    return block2ostream<T>(ostr, block);
}

template <typename T, template <typename U> class Descriptor>
std::ostream &operator<<(std::ostream &ostr, BlockLattice3D<T, Descriptor> const &block)
{
    return block2ostream<T>(ostr, block);
}

/// Flush the content of a generic C++ stream with ASCII content into a Block3D.
/** The content includes external scalars in the case of a BlockLattice3D.
 *
 *  Index-ordering for this operation can be chosen through a call to
 *  global::IOpolicy().setIndexOrderingForStreams(IndexOrdering::OrderingT).
 *
 *  This file format is not exact and should be used for data post-processing only,
 *  and not for checkpointing.
 */
template <typename T>
std::istream &istream2block(std::istream &istr, Block3D &block);

template <typename T>
std::istream &operator>>(std::istream &istr, MultiScalarField3D<T> &block)
{
    return istream2block<T>(istr, block);
}

template <typename T>
std::istream &operator>>(std::istream &istr, MultiNTensorField3D<T> &block)
{
    return istream2block<T>(istr, block);
}

template <typename T, int nDim>
std::istream &operator>>(std::istream &istr, MultiTensorField3D<T, nDim> &block)
{
    return istream2block<T>(istr, block);
}

template <typename T, template <typename U> class Descriptor>
std::istream &operator>>(std::istream &istr, MultiBlockLattice3D<T, Descriptor> &block)
{
    return istream2block<T>(istr, block);
}

template <typename T>
std::istream &operator>>(std::istream &istr, ScalarField3D<T> &block)
{
    return istream2block<T>(istr, block);
}

template <typename T>
std::istream &operator>>(std::istream &istr, NTensorField3D<T> &block)
{
    return istream2block<T>(istr, block);
}

template <typename T, int nDim>
std::istream &operator>>(std::istream &istr, TensorField3D<T, nDim> &block)
{
    return istream2block<T>(istr, block);
}

template <typename T, template <typename U> class Descriptor>
std::istream &operator>>(std::istream &istr, BlockLattice3D<T, Descriptor> &block)
{
    return istream2block<T>(istr, block);
}

}  // namespace plb

#endif  // SERIALIZER_IO_3D_H
