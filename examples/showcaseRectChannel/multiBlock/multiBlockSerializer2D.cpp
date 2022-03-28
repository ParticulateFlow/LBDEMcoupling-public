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
#include "multiBlock/multiBlockSerializer2D.h"

#include "atomicBlock/atomicBlock2D.h"
#include "core/plbDebug.h"
#include "io/parallelIO.h"
#include "parallelism/mpiManager.h"

namespace plb {

////////// class MultiBlockSerializer2D ////////////////////////////

MultiBlockSerializer2D::MultiBlockSerializer2D(
    MultiBlock2D const &multiBlock_, IndexOrdering::OrderingT ordering_) :
    multiBlock(multiBlock_),
    ordering(ordering_),
    domain(multiBlock.getBoundingBox()),
    iX(domain.x0),
    iY(domain.y0),
    buffer(1)  // this avoids buffer of size 0 which one cannot point to
{ }

MultiBlockSerializer2D::MultiBlockSerializer2D(
    MultiBlock2D const &multiBlock_, Box2D domain_, IndexOrdering::OrderingT ordering_) :
    multiBlock(multiBlock_),
    ordering(ordering_),
    domain(domain_),
    iX(domain.x0),
    iY(domain.y0),
    buffer(1)  // this avoids buffer of size 0 which one cannot point to
{ }

MultiBlockSerializer2D *MultiBlockSerializer2D::clone() const
{
    return new MultiBlockSerializer2D(*this);
}

pluint MultiBlockSerializer2D::getSize() const
{
    if (ordering == IndexOrdering::memorySaving) {
        return getSparseBlockStructure().getNumBulkCells() * multiBlock.sizeOfCell();
    } else {
        return domain.nCells() * multiBlock.sizeOfCell();
    }
}

const char *MultiBlockSerializer2D::getNextDataBuffer(pluint &bufferSize) const
{
    PLB_PRECONDITION(!isEmpty());
    EuclideanIterator2D iterator(getSparseBlockStructure());
    if (ordering == IndexOrdering::forward || ordering == IndexOrdering::memorySaving) {
        plint nextBlockId, nextChunkSize;
        bool hasData = iterator.getNextChunkY(iX, iY, nextBlockId, nextChunkSize);
        if (iY + nextChunkSize > domain.y1 + 1) {
            nextChunkSize = domain.y1 - iY + 1;
        }
        if (hasData) {
            computeBufferAlongY(nextBlockId, nextChunkSize);
            bufferSize = nextChunkSize * multiBlock.sizeOfCell();
        } else {
            if (ordering == IndexOrdering::forward) {
                fillBufferWithZeros(nextChunkSize);
                bufferSize = nextChunkSize * multiBlock.sizeOfCell();
            } else {
                bufferSize = 0;
            }
        }
        iY += nextChunkSize;
        if (iY > domain.y1) {
            PLB_ASSERT(iY == domain.y1 + 1);
            iY = domain.y0;
            ++iX;
        }
    } else {
        plint nextBlockId, nextChunkSize;
        bool hasData = iterator.getNextChunkX(iX, iY, nextBlockId, nextChunkSize);
        if (iX + nextChunkSize > domain.x1 + 1) {
            nextChunkSize = domain.x1 - iX + 1;
        }
        if (hasData) {
            computeBufferAlongX(nextBlockId, nextChunkSize);
            bufferSize = nextChunkSize * multiBlock.sizeOfCell();
        } else {
            fillBufferWithZeros(nextChunkSize);
            bufferSize = nextChunkSize * multiBlock.sizeOfCell();
        }
        iX += nextChunkSize;
        if (iX > domain.x1) {
            PLB_ASSERT(iX == domain.x1 + 1);
            iX = domain.x0;
            ++iY;
        }
    }
    if (global::mpi().isMainProcessor()) {
        return &buffer[0];
    } else {
        return 0;
    }
}

bool MultiBlockSerializer2D::isEmpty() const
{
    if (ordering == IndexOrdering::forward || ordering == IndexOrdering::memorySaving) {
        return iX > domain.x1;
    } else {
        return iY > domain.y1;
    }
}

SparseBlockStructure2D const &MultiBlockSerializer2D::getSparseBlockStructure() const
{
    return multiBlock.getMultiBlockManagement().getSparseBlockStructure();
}

bool MultiBlockSerializer2D::isLocal(plint blockId) const
{
    return multiBlock.getMultiBlockManagement().getThreadAttribution().isLocal(blockId);
}

void MultiBlockSerializer2D::computeBufferAlongX(plint nextBlockId, plint nextChunkSize) const
{
    plint bufferSize = nextChunkSize * multiBlock.sizeOfCell();
    bool blockIsLocal = isLocal(nextBlockId);
    if (blockIsLocal) {
        SmartBulk2D bulk(multiBlock.getMultiBlockManagement(), nextBlockId);
        // Avoid pointing to a buffer of size 0, as this leads to undefined behavior.
        PLB_ASSERT(bufferSize > 0);
        buffer.resize(bufferSize);
        plint localX = bulk.toLocalX(iX);
        plint localY = bulk.toLocalY(iY);
        AtomicBlock2D const &nextBlock = multiBlock.getComponent(nextBlockId);
        nextBlock.getDataTransfer().send(
            Box2D(localX, localX + nextChunkSize - 1, localY, localY), buffer,
            modif::staticVariables);
    }
    communicateBuffer(bufferSize, nextBlockId, blockIsLocal);
}

void MultiBlockSerializer2D::computeBufferAlongY(plint nextBlockId, plint nextChunkSize) const
{
    plint bufferSize = nextChunkSize * multiBlock.sizeOfCell();
    bool blockIsLocal = isLocal(nextBlockId);
    if (blockIsLocal) {
        SmartBulk2D bulk(multiBlock.getMultiBlockManagement(), nextBlockId);
        // Avoid pointing to a buffer of size 0, as this leads to undefined behavior.
        PLB_ASSERT(bufferSize > 0);
        buffer.resize(bufferSize);
        plint localX = bulk.toLocalX(iX);
        plint localY = bulk.toLocalY(iY);
        AtomicBlock2D const &nextBlock = multiBlock.getComponent(nextBlockId);
        nextBlock.getDataTransfer().send(
            Box2D(localX, localX, localY, localY + nextChunkSize - 1), buffer,
            modif::staticVariables);
    }
    communicateBuffer(bufferSize, nextBlockId, blockIsLocal);
}

void MultiBlockSerializer2D::communicateBuffer(
    plint bufferSize, plint fromBlockId, bool isAllocated) const
{
#ifdef PLB_MPI_PARALLEL
    // Exchanging a dummy message ahead of the main data has the effect
    //   of synchronizing sender and receiver and therefore avoiding network
    //   jams: given that only processor 0 is receiving and treating the data,
    //   it does a poor job handling all the requests from the other processors
    //   at the same time.
    plint dummyMessage;
    if (isAllocated && !global::mpi().isMainProcessor()) {
        global::mpi().receive(&dummyMessage, 1, 0);
        global::mpi().rSend(&buffer[0], bufferSize, 0);
    }
    if (!isAllocated && global::mpi().isMainProcessor()) {
        int fromProc =
            multiBlock.getMultiBlockManagement().getThreadAttribution().getMpiProcess(fromBlockId);
        // Avoid pointing to a buffer of size 0, as this leads to undefined behavior.
        PLB_ASSERT(bufferSize > 0);
        buffer.resize(bufferSize);
        MPI_Request request;
        MPI_Status status;
        global::mpi().iRecv(&buffer[0], bufferSize, fromProc, &request);
        global::mpi().send(&dummyMessage, 1, fromProc);
        global::mpi().wait(&request, &status);
    }
#endif  // PLB_MPI_PARALLEL
}

void MultiBlockSerializer2D::fillBufferWithZeros(plint nextChunkSize) const
{
    plint bufferSize = nextChunkSize * multiBlock.sizeOfCell();
    if (global::mpi().isMainProcessor()) {
        // Avoid pointing to a buffer of size 0, as this leads to undefined behavior.
        PLB_ASSERT(bufferSize > 0);
        buffer.resize(bufferSize);
        for (plint iBuffer = 0; iBuffer < bufferSize; ++iBuffer) {
            buffer[iBuffer] = 0;
        }
    }
}

////////// class MultiBlockUnSerializer2D ////////////////////////////

MultiBlockUnSerializer2D::MultiBlockUnSerializer2D(
    MultiBlock2D &multiBlock_, IndexOrdering::OrderingT ordering_) :
    multiBlock(multiBlock_),
    ordering(ordering_),
    domain(multiBlock.getBoundingBox()),
    iX(domain.x0),
    iY(domain.y0),
    buffer(1)  // this avoids buffer of size 0 which one cannot point to
{ }

MultiBlockUnSerializer2D::MultiBlockUnSerializer2D(
    MultiBlock2D &multiBlock_, Box2D domain_, IndexOrdering::OrderingT ordering_) :
    multiBlock(multiBlock_),
    ordering(ordering_),
    domain(domain_),
    iX(domain.x0),
    iY(domain.y0),
    buffer(1)  // this avoids buffer of size 0 which one cannot point to
{ }

MultiBlockUnSerializer2D *MultiBlockUnSerializer2D::clone() const
{
    return new MultiBlockUnSerializer2D(*this);
}

pluint MultiBlockUnSerializer2D::getSize() const
{
    if (ordering == IndexOrdering::memorySaving) {
        return getSparseBlockStructure().getNumBulkCells() * multiBlock.sizeOfCell();
    } else {
        return domain.nCells() * multiBlock.sizeOfCell();
    }
}

char *MultiBlockUnSerializer2D::getNextDataBuffer(pluint &bufferSize)
{
    PLB_PRECONDITION(!isFull());
    EuclideanIterator2D iterator(getSparseBlockStructure());
    plint nextBlockId, nextChunkSize;
    if (ordering == IndexOrdering::forward || ordering == IndexOrdering::memorySaving) {
        bool hasData = iterator.getNextChunkY(iX, iY, nextBlockId, nextChunkSize);
        if (iY + nextChunkSize > domain.y1 + 1) {
            nextChunkSize = domain.y1 - iY + 1;
        }
        bufferSize = nextChunkSize * multiBlock.sizeOfCell();
        if (ordering == IndexOrdering::memorySaving && !hasData) {
            bufferSize = 0;
        }
    } else {
        iterator.getNextChunkX(iX, iY, nextBlockId, nextChunkSize);
        if (iX + nextChunkSize > domain.x1 + 1) {
            nextChunkSize = domain.x1 - iX + 1;
        }
        bufferSize = nextChunkSize * multiBlock.sizeOfCell();
    }
    if (global::mpi().isMainProcessor()) {
        // Avoid pointing to a buffer of size 0, as this leads to undefined behavior.
        PLB_ASSERT(bufferSize > 0);
        buffer.resize(bufferSize);
    }
    if (global::mpi().isMainProcessor()) {
        return &buffer[0];
    } else {
        return 0;
    }
}

void MultiBlockUnSerializer2D::commitData()
{
    PLB_PRECONDITION(!isFull());
    EuclideanIterator2D iterator(getSparseBlockStructure());
    if (ordering == IndexOrdering::forward || ordering == IndexOrdering::memorySaving) {
        plint nextBlockId, nextChunkSize;
        bool hasData = iterator.getNextChunkY(iX, iY, nextBlockId, nextChunkSize);
        if (iY + nextChunkSize > domain.y1 + 1) {
            nextChunkSize = domain.y1 - iY + 1;
        }
        if (hasData) {
            fillBufferAlongY(nextBlockId, nextChunkSize);
        }
        iY += nextChunkSize;
        if (iY > domain.y1) {
            iY = domain.y0;
            ++iX;
        }
    } else {
        plint nextBlockId, nextChunkSize;
        bool hasData = iterator.getNextChunkX(iX, iY, nextBlockId, nextChunkSize);
        if (iX + nextChunkSize > domain.x1 + 1) {
            nextChunkSize = domain.x1 - iX + 1;
        }
        if (hasData) {
            fillBufferAlongX(nextBlockId, nextChunkSize);
        }
        iX += nextChunkSize;
        if (iX > domain.x1) {
            iX = domain.x0;
            ++iY;
        }
    }
    // At the end of unserialization, duplicate overlaps.
    if (isFull()) {
        multiBlock.duplicateOverlaps(modif::staticVariables);
    }
}

bool MultiBlockUnSerializer2D::isFull() const
{
    if (ordering == IndexOrdering::forward || ordering == IndexOrdering::memorySaving) {
        return iX > domain.x1;
    } else {
        return iY > domain.y1;
    }
}

SparseBlockStructure2D const &MultiBlockUnSerializer2D::getSparseBlockStructure() const
{
    return multiBlock.getMultiBlockManagement().getSparseBlockStructure();
}

bool MultiBlockUnSerializer2D::isLocal(plint blockId) const
{
    return multiBlock.getMultiBlockManagement().getThreadAttribution().isLocal(blockId);
}

void MultiBlockUnSerializer2D::fillBufferAlongX(plint nextBlockId, plint nextChunkSize)
{
    plint bufferSize = nextChunkSize * multiBlock.sizeOfCell();
    bool blockIsLocal = isLocal(nextBlockId);
    communicateBuffer(bufferSize, nextBlockId, blockIsLocal);
    if (blockIsLocal) {
        SmartBulk2D bulk(multiBlock.getMultiBlockManagement(), nextBlockId);
        plint localX = bulk.toLocalX(iX);
        plint localY = bulk.toLocalY(iY);
        AtomicBlock2D &nextBlock = multiBlock.getComponent(nextBlockId);
        nextBlock.getDataTransfer().receive(
            Box2D(localX, localX + nextChunkSize - 1, localY, localY), buffer,
            modif::staticVariables);
    }
}

void MultiBlockUnSerializer2D::fillBufferAlongY(plint nextBlockId, plint nextChunkSize)
{
    plint bufferSize = nextChunkSize * multiBlock.sizeOfCell();
    bool blockIsLocal = isLocal(nextBlockId);
    communicateBuffer(bufferSize, nextBlockId, blockIsLocal);
    if (blockIsLocal) {
        SmartBulk2D bulk(multiBlock.getMultiBlockManagement(), nextBlockId);
        plint localX = bulk.toLocalX(iX);
        plint localY = bulk.toLocalY(iY);
        AtomicBlock2D &nextBlock = multiBlock.getComponent(nextBlockId);
        nextBlock.getDataTransfer().receive(
            Box2D(localX, localX, localY, localY + nextChunkSize - 1), buffer,
            modif::staticVariables);
    }
}

void MultiBlockUnSerializer2D::communicateBuffer(
    plint bufferSize, plint toBlockId, bool isAllocated) const
{
#ifdef PLB_MPI_PARALLEL
    if (isAllocated && !global::mpi().isMainProcessor()) {
        // Avoid pointing to a buffer of size 0, as this leads to undefined behavior.
        PLB_ASSERT(bufferSize > 0);
        buffer.resize(bufferSize);
        global::mpi().receive(&buffer[0], bufferSize, 0);
    }
    if (!isAllocated && global::mpi().isMainProcessor()) {
        int toProc =
            multiBlock.getMultiBlockManagement().getThreadAttribution().getMpiProcess(toBlockId);
        global::mpi().send(&buffer[0], bufferSize, toProc);
    }
#endif
}

MultiBlockFastSerializer2D::MultiBlockFastSerializer2D(
    MultiBlock2D const &multiBlock_, IndexOrdering::OrderingT ordering_) :
    multiBlock(multiBlock_), ordering(ordering_), domain(multiBlock.getBoundingBox())
{
    if (ordering == IndexOrdering::forward) {
        pos = domain.x0;
    } else if (ordering == IndexOrdering::backward) {
        pos = domain.y0;
    } else {
        // Sparse ordering not implemented.
        PLB_ASSERT(false);
    }
}

MultiBlockFastSerializer2D::MultiBlockFastSerializer2D(
    MultiBlock2D const &multiBlock_, Box2D domain_, IndexOrdering::OrderingT ordering_) :
    multiBlock(multiBlock_), ordering(ordering_), domain(domain_)
{
    if (ordering == IndexOrdering::forward) {
        pos = domain.x0;
    } else if (ordering == IndexOrdering::backward) {
        pos = domain.y0;
    } else {
        // Sparse ordering not implemented.
        PLB_ASSERT(false);
    }
}

MultiBlockFastSerializer2D *MultiBlockFastSerializer2D::clone() const
{
    return new MultiBlockFastSerializer2D(*this);
}

pluint MultiBlockFastSerializer2D::getSize() const
{
    return domain.nCells() * multiBlock.sizeOfCell();
}

const char *MultiBlockFastSerializer2D::getNextDataBuffer(pluint &bufferSize) const
{
    bufferSize = computeSlice();
    if (buffer.empty())
        buffer.resize(1);
    return &buffer[0];
}

bool MultiBlockFastSerializer2D::isEmpty() const
{
    if (ordering == IndexOrdering::forward) {
        return pos > domain.x1;
    } else {
        return pos > domain.y1;
    }
}

pluint MultiBlockFastSerializer2D::computeSlice() const
{
    Box2D slice(domain);
    if (ordering == IndexOrdering::forward) {
        slice.x0 = slice.x1 = pos;
    } else if (ordering == IndexOrdering::backward) {
        slice.y0 = slice.y1 = pos;
    } else {
        // Sparse ordering not implemented.
        PLB_ASSERT(false);
    }
    SparseBlockStructure2D blockStructure(multiBlock.getBoundingBox());
    blockStructure.addBlock(slice, 0);
    plint envelopeWidth = 1;
    MultiBlockManagement2D serialMultiBlockManagement(
        blockStructure, new OneToOneThreadAttribution, envelopeWidth);
    MultiBlock2D *multiSerialBlock = multiBlock.clone(serialMultiBlockManagement);

    pluint bufferSize = slice.nCells() * multiBlock.sizeOfCell();
    if (global::mpi().isMainProcessor()) {
        AtomicBlock2D &atomicSerialBlock = multiSerialBlock->getComponent(0);
        SmartBulk2D oneBlockSlice(serialMultiBlockManagement, 0);
        Box2D localSlice(oneBlockSlice.toLocal(slice));
        atomicSerialBlock.getDataTransfer().send(localSlice, buffer, modif::staticVariables);
        PLB_ASSERT(bufferSize == buffer.size());
    }
    delete multiSerialBlock;
    ++pos;
    return bufferSize;
}

MultiBlockFastUnSerializer2D::MultiBlockFastUnSerializer2D(
    MultiBlock2D &multiBlock_, IndexOrdering::OrderingT ordering_) :
    multiBlock(multiBlock_), ordering(ordering_), domain(multiBlock.getBoundingBox())
{
    if (ordering == IndexOrdering::forward) {
        pos = domain.x0;
    } else if (ordering == IndexOrdering::backward) {
        pos = domain.y0;
    } else {
        // Sparse ordering not implemented.
        PLB_ASSERT(false);
    }
}

MultiBlockFastUnSerializer2D::MultiBlockFastUnSerializer2D(
    MultiBlock2D &multiBlock_, Box2D domain_, IndexOrdering::OrderingT ordering_) :
    multiBlock(multiBlock_), ordering(ordering_), domain(domain_)
{
    if (ordering == IndexOrdering::forward) {
        pos = domain.x0;
    } else if (ordering == IndexOrdering::backward) {
        pos = domain.y0;
    } else {
        // Sparse ordering not implemented.
        PLB_ASSERT(false);
    }
}

MultiBlockFastUnSerializer2D *MultiBlockFastUnSerializer2D::clone() const
{
    return new MultiBlockFastUnSerializer2D(*this);
}

pluint MultiBlockFastUnSerializer2D::getSize() const
{
    return domain.nCells() * multiBlock.sizeOfCell();
}

char *MultiBlockFastUnSerializer2D::getNextDataBuffer(pluint &bufferSize)
{
    if (ordering == IndexOrdering::forward) {
        bufferSize = domain.getNy() * multiBlock.sizeOfCell();
    } else {
        bufferSize = domain.getNx() * multiBlock.sizeOfCell();
    }

    if (global::mpi().isMainProcessor()) {
        buffer.resize(bufferSize);
    } else {
        buffer.clear();
    }
    if (buffer.empty()) {
        buffer.resize(1);
    }
    return &buffer[0];
}

void MultiBlockFastUnSerializer2D::commitData()
{
    Box2D slice(domain);
    if (ordering == IndexOrdering::forward) {
        slice.x0 = slice.x1 = pos;
    } else if (ordering == IndexOrdering::backward) {
        slice.y0 = slice.y1 = pos;
    } else {
        // Sparse ordering not implemented.
        PLB_ASSERT(false);
    }
    SparseBlockStructure2D blockStructure(multiBlock.getBoundingBox());
    blockStructure.addBlock(slice, 0);
    plint envelopeWidth = 1;
    MultiBlockManagement2D serialMultiBlockManagement(
        blockStructure, new OneToOneThreadAttribution, envelopeWidth);
    MultiBlock2D *multiSerialBlock = multiBlock.clone(serialMultiBlockManagement);
    if (global::mpi().isMainProcessor()) {
        AtomicBlock2D &atomicSerialBlock = multiSerialBlock->getComponent(0);
        SmartBulk2D oneBlockBulk(serialMultiBlockManagement, 0);
        Box2D localSlice(oneBlockBulk.toLocal(slice));
        atomicSerialBlock.getDataTransfer().receive(localSlice, buffer, modif::staticVariables);
    }
    multiBlock.copyReceive(*multiSerialBlock, slice, slice);
    delete multiSerialBlock;
    ++pos;
}

bool MultiBlockFastUnSerializer2D::isFull() const
{
    if (ordering == IndexOrdering::forward) {
        return pos > domain.x1;
    } else {
        return pos > domain.y1;
    }
}

}  //  namespace plb
