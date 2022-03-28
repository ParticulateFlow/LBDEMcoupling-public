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
 * Serializer and UnSerializer for atomic blocks -- generic implementation.
 */

#include "atomicBlock/atomicBlockSerializer2D.h"

#include "core/plbDebug.h"

namespace plb {

////////// class AtomicBlockSerializer2D ////////////////////////////

AtomicBlockSerializer2D::AtomicBlockSerializer2D(
    AtomicBlock2D const &block_, IndexOrdering::OrderingT ordering_) :
    block(block_), ordering(ordering_), domain(block.getBoundingBox()), iX(domain.x0), iY(domain.y0)
{ }

AtomicBlockSerializer2D::AtomicBlockSerializer2D(
    AtomicBlock2D const &block_, Box2D domain_, IndexOrdering::OrderingT ordering_) :
    block(block_), ordering(ordering_), domain(domain_), iX(domain.x0), iY(domain.y0)
{ }

AtomicBlockSerializer2D *AtomicBlockSerializer2D::clone() const
{
    return new AtomicBlockSerializer2D(*this);
}

pluint AtomicBlockSerializer2D::getSize() const
{
    return domain.nCells() * block.getDataTransfer().staticCellSize();
}

const char *AtomicBlockSerializer2D::getNextDataBuffer(pluint &bufferSize) const
{
    PLB_PRECONDITION(!isEmpty());
    if (ordering == IndexOrdering::forward || ordering == IndexOrdering::memorySaving) {
        bufferSize = domain.getNy() * block.getDataTransfer().staticCellSize();
        buffer.resize(bufferSize);
        block.getDataTransfer().send(
            Box2D(iX, iX, domain.y0, domain.y1), buffer, modif::staticVariables);
        ++iX;
    } else {
        bufferSize = domain.getNx() * block.getDataTransfer().staticCellSize();
        buffer.resize(bufferSize);
        block.getDataTransfer().send(
            Box2D(domain.x0, domain.x1, iY, iY), buffer, modif::staticVariables);
        ++iY;
    }
    return &buffer[0];
}

bool AtomicBlockSerializer2D::isEmpty() const
{
    if (ordering == IndexOrdering::forward || ordering == IndexOrdering::memorySaving) {
        return iX > domain.x1;
    } else {
        return iY > domain.y1;
    }
}

////////// class AtomicBlockUnSerializer2D ////////////////////////////

AtomicBlockUnSerializer2D::AtomicBlockUnSerializer2D(
    AtomicBlock2D &block_, IndexOrdering::OrderingT ordering_) :
    block(block_), ordering(ordering_), domain(block.getBoundingBox()), iX(domain.x0), iY(domain.y0)
{ }

AtomicBlockUnSerializer2D::AtomicBlockUnSerializer2D(
    AtomicBlock2D &block_, Box2D domain_, IndexOrdering::OrderingT ordering_) :
    block(block_), ordering(ordering_), domain(domain_), iX(domain.x0), iY(domain.y0)
{ }

AtomicBlockUnSerializer2D *AtomicBlockUnSerializer2D::clone() const
{
    return new AtomicBlockUnSerializer2D(*this);
}

pluint AtomicBlockUnSerializer2D::getSize() const
{
    return domain.nCells() * block.getDataTransfer().staticCellSize();
}

char *AtomicBlockUnSerializer2D::getNextDataBuffer(pluint &bufferSize)
{
    PLB_PRECONDITION(!isFull());
    if (ordering == IndexOrdering::forward || ordering == IndexOrdering::memorySaving) {
        bufferSize = domain.getNy() * block.getDataTransfer().staticCellSize();
    } else {
        bufferSize = domain.getNx() * block.getDataTransfer().staticCellSize();
    }
    buffer.resize(bufferSize);
    return &buffer[0];
}

void AtomicBlockUnSerializer2D::commitData()
{
    PLB_PRECONDITION(!isFull());
    if (ordering == IndexOrdering::forward || ordering == IndexOrdering::memorySaving) {
        block.getDataTransfer().receive(
            Box2D(iX, iX, domain.y0, domain.y1), buffer, modif::staticVariables);
        ++iX;
    } else {
        block.getDataTransfer().receive(
            Box2D(domain.x0, domain.x1, iY, iY), buffer, modif::staticVariables);
        ++iY;
    }
}

bool AtomicBlockUnSerializer2D::isFull() const
{
    if (ordering == IndexOrdering::forward || ordering == IndexOrdering::memorySaving) {
        return iX > domain.x1;
    } else {
        return iY > domain.y1;
    }
}

}  //  namespace plb
