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
 * Serializer and UnSerializer for multi blocks -- header file.
 */
#ifndef MULTI_BLOCK_SERIALIZER_3D_H
#define MULTI_BLOCK_SERIALIZER_3D_H

#include "core/globalDefs.h"
#include "core/serializer.h"
#include "core/util.h"
#include "multiBlock/multiBlock3D.h"

namespace plb {

class MultiBlockSerializer3D : public DataSerializer {
public:
    MultiBlockSerializer3D(MultiBlock3D const &multiBlock_, IndexOrdering::OrderingT ordering_);
    MultiBlockSerializer3D(
        MultiBlock3D const &multiBlock_, Box3D domain_, IndexOrdering::OrderingT ordering_);
    virtual MultiBlockSerializer3D *clone() const;
    virtual pluint getSize() const;
    virtual const char *getNextDataBuffer(pluint &bufferSize) const;
    virtual bool isEmpty() const;

private:
    SparseBlockStructure3D const &getSparseBlockStructure() const;
    bool isLocal(plint blockId) const;
    void computeBufferAlongX(plint nextBlockId, plint nextChunkSize) const;
    void computeBufferAlongY(plint nextBlockId, plint nextChunkSize) const;
    void computeBufferAlongZ(plint nextBlockId, plint nextChunkSize) const;
    void communicateBuffer(plint bufferSize, plint fromBlockId, bool isAllocated) const;
    void fillBufferWithZeros(plint nextChunkSize) const;

private:
    MultiBlock3D const &multiBlock;
    IndexOrdering::OrderingT ordering;
    Box3D domain;
    mutable plint iX, iY, iZ;
    mutable std::vector<char> buffer;
};

class MultiBlockUnSerializer3D : public DataUnSerializer {
public:
    MultiBlockUnSerializer3D(MultiBlock3D &multiBlock_, IndexOrdering::OrderingT ordering_);
    MultiBlockUnSerializer3D(
        MultiBlock3D &multiBlock_, Box3D domain_, IndexOrdering::OrderingT ordering_);
    virtual MultiBlockUnSerializer3D *clone() const;
    virtual pluint getSize() const;
    virtual char *getNextDataBuffer(pluint &bufferSize);
    virtual void commitData();
    virtual bool isFull() const;

private:
    SparseBlockStructure3D const &getSparseBlockStructure() const;
    bool isLocal(plint blockId) const;
    void fillBufferAlongX(plint nextBlockId, plint nextChunkSize);
    void fillBufferAlongY(plint nextBlockId, plint nextChunkSize);
    void fillBufferAlongZ(plint nextBlockId, plint nextChunkSize);
    void communicateBuffer(plint bufferSize, plint toBlockId, bool isAllocated) const;

private:
    MultiBlock3D &multiBlock;
    IndexOrdering::OrderingT ordering;
    Box3D domain;
    mutable plint iX, iY, iZ;
    mutable std::vector<char> buffer;
};

class MultiBlockFastSerializer3D : public DataSerializer {
public:
    MultiBlockFastSerializer3D(MultiBlock3D const &multiBlock_, IndexOrdering::OrderingT ordering_);
    MultiBlockFastSerializer3D(
        MultiBlock3D const &multiBlock_, Box3D domain_, IndexOrdering::OrderingT ordering_);
    virtual MultiBlockFastSerializer3D *clone() const;
    virtual pluint getSize() const;
    virtual const char *getNextDataBuffer(pluint &bufferSize) const;
    virtual bool isEmpty() const;

private:
    pluint computeSlice() const;

private:
    MultiBlock3D const &multiBlock;
    IndexOrdering::OrderingT ordering;
    Box3D domain;
    mutable plint pos;
    mutable std::vector<char> buffer;
};

class MultiBlockFastUnSerializer3D : public DataUnSerializer {
public:
    MultiBlockFastUnSerializer3D(MultiBlock3D &multiBlock_, IndexOrdering::OrderingT ordering_);
    MultiBlockFastUnSerializer3D(
        MultiBlock3D &multiBlock_, Box3D domain_, IndexOrdering::OrderingT ordering_);
    virtual MultiBlockFastUnSerializer3D *clone() const;
    virtual pluint getSize() const;
    virtual char *getNextDataBuffer(pluint &bufferSize);
    virtual void commitData();
    virtual bool isFull() const;

private:
    MultiBlock3D &multiBlock;
    IndexOrdering::OrderingT ordering;
    Box3D domain;
    mutable plint pos;
    mutable std::vector<char> buffer;
};

}  //  namespace plb

#endif  // MULTI_BLOCK_SERIALIZER_3D_H
