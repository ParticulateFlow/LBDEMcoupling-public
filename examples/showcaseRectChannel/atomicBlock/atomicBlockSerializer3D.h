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
 * Serializer and UnSerializer for atomic blocks -- header file.
 */
#ifndef ATOMIC_BLOCK_SERIALIZER_3D_H
#define ATOMIC_BLOCK_SERIALIZER_3D_H

#include "atomicBlock/atomicBlock3D.h"
#include "core/globalDefs.h"
#include "core/serializer.h"
#include "core/util.h"

namespace plb {

class AtomicBlockSerializer3D : public DataSerializer {
public:
    AtomicBlockSerializer3D(AtomicBlock3D const &block_, IndexOrdering::OrderingT ordering_);
    AtomicBlockSerializer3D(
        AtomicBlock3D const &block_, Box3D domain_, IndexOrdering::OrderingT ordering_);
    virtual AtomicBlockSerializer3D *clone() const;
    virtual pluint getSize() const;
    virtual const char *getNextDataBuffer(pluint &bufferSize) const;
    virtual bool isEmpty() const;

private:
    AtomicBlock3D const &block;
    IndexOrdering::OrderingT ordering;
    Box3D domain;
    mutable std::vector<char> buffer;
    mutable plint iX, iY, iZ;
};

class AtomicBlockUnSerializer3D : public DataUnSerializer {
public:
    AtomicBlockUnSerializer3D(AtomicBlock3D &block_, IndexOrdering::OrderingT ordering_);
    AtomicBlockUnSerializer3D(
        AtomicBlock3D &block_, Box3D domain_, IndexOrdering::OrderingT ordering_);
    virtual AtomicBlockUnSerializer3D *clone() const;
    virtual pluint getSize() const;
    virtual char *getNextDataBuffer(pluint &bufferSize);
    virtual void commitData();
    virtual bool isFull() const;

private:
    AtomicBlock3D &block;
    IndexOrdering::OrderingT ordering;
    Box3D domain;
    mutable std::vector<char> buffer;
    mutable plint iX, iY, iZ;
};

}  //  namespace plb

#endif  // ATOMIC_BLOCK_SERIALIZER_3D_H
