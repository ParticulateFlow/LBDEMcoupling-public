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
#include "multiBlock/multiBlockSerializer3D.h"

#include "atomicBlock/atomicBlock3D.h"
#include "core/plbDebug.h"
#include "parallelism/mpiManager.h"

namespace plb {

////////// class MultiBlockSerializer3D ////////////////////////////

MultiBlockSerializer3D::MultiBlockSerializer3D(
    MultiBlock3D const &multiBlock_, IndexOrdering::OrderingT ordering_) :
    multiBlock(multiBlock_),
    ordering(ordering_),
    domain(multiBlock.getBoundingBox()),
    iX(domain.x0),
    iY(domain.y0),
    iZ(domain.z0),
    buffer(1)  // this avoids buffer of size 0 which one cannot point to
{ }

MultiBlockSerializer3D::MultiBlockSerializer3D(
    MultiBlock3D const &multiBlock_, Box3D domain_, IndexOrdering::OrderingT ordering_) :
    multiBlock(multiBlock_),
    ordering(ordering_),
    domain(domain_),
    iX(domain.x0),
    iY(domain.y0),
    iZ(domain.z0),
    buffer(1)  // this avoids buffer of size 0 which one cannot point to
{ }

MultiBlockSerializer3D *MultiBlockSerializer3D::clone() const
{
    return new MultiBlockSerializer3D(*this);
}

pluint MultiBlockSerializer3D::getSize() const
{
    if (ordering == IndexOrdering::memorySaving) {
        return getSparseBlockStructure().getNumBulkCells() * multiBlock.sizeOfCell();
    } else {
        return domain.nCells() * multiBlock.sizeOfCell();
    }
}

const char *MultiBlockSerializer3D::getNextDataBuffer(pluint &bufferSize) const
{
    PLB_PRECONDITION(!isEmpty());
    EuclideanIterator3D iterator(getSparseBlockStructure());
    if (ordering == IndexOrdering::forward || ordering == IndexOrdering::memorySaving) {
        plint nextBlockId, nextChunkSize;
        bool hasData = iterator.getNextChunkZ(iX, iY, iZ, nextBlockId, nextChunkSize);
        if (iZ + nextChunkSize > domain.z1 + 1) {
            nextChunkSize = domain.z1 - iZ + 1;
        }
        if (hasData) {
            computeBufferAlongZ(nextBlockId, nextChunkSize);
            bufferSize = nextChunkSize * multiBlock.sizeOfCell();
        } else {
            if (ordering == IndexOrdering::forward) {
                fillBufferWithZeros(nextChunkSize);
                bufferSize = nextChunkSize * multiBlock.sizeOfCell();
            } else {
                bufferSize = 0;
            }
        }
        iZ += nextChunkSize;
        if (iZ > domain.z1) {
            PLB_ASSERT(iZ == domain.z1 + 1);
            iZ = domain.z0;
            ++iY;
            if (iY > domain.y1) {
                PLB_ASSERT(iY == domain.y1 + 1);
                iY = domain.y0;
                ++iX;
            }
        }
    } else {
        plint nextBlockId, nextChunkSize;
        bool hasData = iterator.getNextChunkX(iX, iY, iZ, nextBlockId, nextChunkSize);
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
            if (iY > domain.y1) {
                PLB_ASSERT(iY == domain.y1 + 1);
                iY = domain.y0;
                ++iZ;
            }
        }
    }
    if (global::mpi().isMainProcessor()) {
        return &buffer[0];
    } else {
        return 0;
    }
}

bool MultiBlockSerializer3D::isEmpty() const
{
    if (ordering == IndexOrdering::forward || ordering == IndexOrdering::memorySaving) {
        return iX > domain.x1;
    } else {
        return iZ > domain.z1;
    }
}

SparseBlockStructure3D const &MultiBlockSerializer3D::getSparseBlockStructure() const
{
    return multiBlock.getMultiBlockManagement().getSparseBlockStructure();
}

bool MultiBlockSerializer3D::isLocal(plint blockId) const
{
    return multiBlock.getMultiBlockManagement().getThreadAttribution().isLocal(blockId);
}

void MultiBlockSerializer3D::computeBufferAlongX(plint nextBlockId, plint nextChunkSize) const
{
    plint bufferSize = nextChunkSize * multiBlock.sizeOfCell();
    bool blockIsLocal = isLocal(nextBlockId);
    if (blockIsLocal) {
        SmartBulk3D bulk(multiBlock.getMultiBlockManagement(), nextBlockId);
        // Avoid pointing to a buffer of size 0, as this leads to undefined behavior.
        PLB_ASSERT(bufferSize > 0);
        buffer.resize(bufferSize);
        plint localX = bulk.toLocalX(iX);
        plint localY = bulk.toLocalY(iY);
        plint localZ = bulk.toLocalZ(iZ);
        AtomicBlock3D const &nextBlock = multiBlock.getComponent(nextBlockId);
        nextBlock.getDataTransfer().send(
            Box3D(localX, localX + nextChunkSize - 1, localY, localY, localZ, localZ), buffer,
            modif::staticVariables);
    }
    communicateBuffer(bufferSize, nextBlockId, blockIsLocal);
}

void MultiBlockSerializer3D::computeBufferAlongY(plint nextBlockId, plint nextChunkSize) const
{
    plint bufferSize = nextChunkSize * multiBlock.sizeOfCell();
    bool blockIsLocal = isLocal(nextBlockId);
    if (blockIsLocal) {
        SmartBulk3D bulk(multiBlock.getMultiBlockManagement(), nextBlockId);
        // Avoid pointing to a buffer of size 0, as this leads to undefined behavior.
        PLB_ASSERT(bufferSize > 0);
        buffer.resize(bufferSize);
        plint localX = bulk.toLocalX(iX);
        plint localY = bulk.toLocalY(iY);
        plint localZ = bulk.toLocalZ(iZ);
        AtomicBlock3D const &nextBlock = multiBlock.getComponent(nextBlockId);
        nextBlock.getDataTransfer().send(
            Box3D(localX, localX, localY, localY + nextChunkSize - 1, localZ, localZ), buffer,
            modif::staticVariables);
    }
    communicateBuffer(bufferSize, nextBlockId, blockIsLocal);
}

void MultiBlockSerializer3D::computeBufferAlongZ(plint nextBlockId, plint nextChunkSize) const
{
    plint bufferSize = nextChunkSize * multiBlock.sizeOfCell();
    bool blockIsLocal = isLocal(nextBlockId);
    if (blockIsLocal) {
        SmartBulk3D bulk(multiBlock.getMultiBlockManagement(), nextBlockId);
        // Avoid pointing to a buffer of size 0, as this leads to undefined behavior.
        PLB_ASSERT(bufferSize > 0);
        buffer.resize(bufferSize);
        plint localX = bulk.toLocalX(iX);
        plint localY = bulk.toLocalY(iY);
        plint localZ = bulk.toLocalZ(iZ);
        AtomicBlock3D const &nextBlock = multiBlock.getComponent(nextBlockId);
        nextBlock.getDataTransfer().send(
            Box3D(localX, localX, localY, localY, localZ, localZ + nextChunkSize - 1), buffer,
            modif::staticVariables);
    }
    communicateBuffer(bufferSize, nextBlockId, blockIsLocal);
}

void MultiBlockSerializer3D::communicateBuffer(
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

void MultiBlockSerializer3D::fillBufferWithZeros(plint nextChunkSize) const
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

////////// class MultiBlockUnSerializer3D ////////////////////////////

MultiBlockUnSerializer3D::MultiBlockUnSerializer3D(
    MultiBlock3D &multiBlock_, IndexOrdering::OrderingT ordering_) :
    multiBlock(multiBlock_),
    ordering(ordering_),
    domain(multiBlock.getBoundingBox()),
    iX(domain.x0),
    iY(domain.y0),
    iZ(domain.z0),
    buffer(1)  // this avoids buffer of size 0 which one cannot point to
{ }

MultiBlockUnSerializer3D::MultiBlockUnSerializer3D(
    MultiBlock3D &multiBlock_, Box3D domain_, IndexOrdering::OrderingT ordering_) :
    multiBlock(multiBlock_),
    ordering(ordering_),
    domain(domain_),
    iX(domain.x0),
    iY(domain.y0),
    iZ(domain.z0),
    buffer(1)  // this avoids buffer of size 0 which one cannot point to
{ }

MultiBlockUnSerializer3D *MultiBlockUnSerializer3D::clone() const
{
    return new MultiBlockUnSerializer3D(*this);
}

pluint MultiBlockUnSerializer3D::getSize() const
{
    if (ordering == IndexOrdering::memorySaving) {
        return getSparseBlockStructure().getNumBulkCells() * multiBlock.sizeOfCell();
    } else {
        return domain.nCells() * multiBlock.sizeOfCell();
    }
}

char *MultiBlockUnSerializer3D::getNextDataBuffer(pluint &bufferSize)
{
    PLB_PRECONDITION(!isFull());
    EuclideanIterator3D iterator(getSparseBlockStructure());
    plint nextBlockId, nextChunkSize;
    if (ordering == IndexOrdering::forward || ordering == IndexOrdering::memorySaving) {
        bool hasData = iterator.getNextChunkZ(iX, iY, iZ, nextBlockId, nextChunkSize);
        if (iZ + nextChunkSize > domain.z1 + 1) {
            nextChunkSize = domain.z1 - iZ + 1;
        }
        bufferSize = nextChunkSize * multiBlock.sizeOfCell();
        if (ordering == IndexOrdering::memorySaving && !hasData) {
            bufferSize = 0;
        }
    } else {
        iterator.getNextChunkX(iX, iY, iZ, nextBlockId, nextChunkSize);
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

void MultiBlockUnSerializer3D::commitData()
{
    PLB_PRECONDITION(!isFull());
    EuclideanIterator3D iterator(getSparseBlockStructure());
    if (ordering == IndexOrdering::forward || ordering == IndexOrdering::memorySaving) {
        plint nextBlockId, nextChunkSize;
        bool hasData = iterator.getNextChunkZ(iX, iY, iZ, nextBlockId, nextChunkSize);
        if (iZ + nextChunkSize > domain.z1 + 1) {
            nextChunkSize = domain.z1 - iZ + 1;
        }
        if (hasData) {
            fillBufferAlongZ(nextBlockId, nextChunkSize);
        }
        iZ += nextChunkSize;
        if (iZ > domain.z1) {
            iZ = domain.z0;
            ++iY;
            if (iY > domain.y1) {
                iY = domain.y0;
                ++iX;
            }
        }
    } else {
        plint nextBlockId, nextChunkSize;
        bool hasData = iterator.getNextChunkX(iX, iY, iZ, nextBlockId, nextChunkSize);
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
            if (iY > domain.y1) {
                iY = domain.y0;
                ++iZ;
            }
        }
    }
    // At the end of unserialization, duplicate overlaps.
    if (isFull()) {
        multiBlock.duplicateOverlaps(modif::staticVariables);
    }
}

bool MultiBlockUnSerializer3D::isFull() const
{
    if (ordering == IndexOrdering::forward || ordering == IndexOrdering::memorySaving) {
        return iX > domain.x1;
    } else {
        return iZ > domain.z1;
    }
}

SparseBlockStructure3D const &MultiBlockUnSerializer3D::getSparseBlockStructure() const
{
    return multiBlock.getMultiBlockManagement().getSparseBlockStructure();
}

bool MultiBlockUnSerializer3D::isLocal(plint blockId) const
{
    return multiBlock.getMultiBlockManagement().getThreadAttribution().isLocal(blockId);
}

void MultiBlockUnSerializer3D::fillBufferAlongX(plint nextBlockId, plint nextChunkSize)
{
    plint bufferSize = nextChunkSize * multiBlock.sizeOfCell();
    bool blockIsLocal = isLocal(nextBlockId);
    communicateBuffer(bufferSize, nextBlockId, blockIsLocal);
    if (blockIsLocal) {
        SmartBulk3D bulk(multiBlock.getMultiBlockManagement(), nextBlockId);
        plint localX = bulk.toLocalX(iX);
        plint localY = bulk.toLocalY(iY);
        plint localZ = bulk.toLocalZ(iZ);
        AtomicBlock3D &nextBlock = multiBlock.getComponent(nextBlockId);
        nextBlock.getDataTransfer().receive(
            Box3D(localX, localX + nextChunkSize - 1, localY, localY, localZ, localZ), buffer,
            modif::staticVariables);
    }
}

void MultiBlockUnSerializer3D::fillBufferAlongY(plint nextBlockId, plint nextChunkSize)
{
    plint bufferSize = nextChunkSize * multiBlock.sizeOfCell();
    bool blockIsLocal = isLocal(nextBlockId);
    communicateBuffer(bufferSize, nextBlockId, blockIsLocal);
    if (blockIsLocal) {
        SmartBulk3D bulk(multiBlock.getMultiBlockManagement(), nextBlockId);
        plint localX = bulk.toLocalX(iX);
        plint localY = bulk.toLocalY(iY);
        plint localZ = bulk.toLocalZ(iZ);
        AtomicBlock3D &nextBlock = multiBlock.getComponent(nextBlockId);
        nextBlock.getDataTransfer().receive(
            Box3D(localX, localX, localY, localY + nextChunkSize - 1, localZ, localZ), buffer,
            modif::staticVariables);
    }
}

void MultiBlockUnSerializer3D::fillBufferAlongZ(plint nextBlockId, plint nextChunkSize)
{
    plint bufferSize = nextChunkSize * multiBlock.sizeOfCell();
    bool blockIsLocal = isLocal(nextBlockId);
    communicateBuffer(bufferSize, nextBlockId, blockIsLocal);
    if (blockIsLocal) {
        SmartBulk3D bulk(multiBlock.getMultiBlockManagement(), nextBlockId);
        plint localX = bulk.toLocalX(iX);
        plint localY = bulk.toLocalY(iY);
        plint localZ = bulk.toLocalZ(iZ);
        AtomicBlock3D &nextBlock = multiBlock.getComponent(nextBlockId);
        nextBlock.getDataTransfer().receive(
            Box3D(localX, localX, localY, localY, localZ, localZ + nextChunkSize - 1), buffer,
            modif::staticVariables);
    }
}

void MultiBlockUnSerializer3D::communicateBuffer(
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

MultiBlockFastSerializer3D::MultiBlockFastSerializer3D(
    MultiBlock3D const &multiBlock_, IndexOrdering::OrderingT ordering_) :
    multiBlock(multiBlock_), ordering(ordering_), domain(multiBlock.getBoundingBox())
{
    if (ordering == IndexOrdering::forward) {
        pos = domain.x0;
    } else if (ordering == IndexOrdering::backward) {
        pos = domain.z0;
    } else {
        // Sparse ordering not implemented.
        PLB_ASSERT(false);
    }
}

MultiBlockFastSerializer3D::MultiBlockFastSerializer3D(
    MultiBlock3D const &multiBlock_, Box3D domain_, IndexOrdering::OrderingT ordering_) :
    multiBlock(multiBlock_), ordering(ordering_), domain(domain_)
{
    if (ordering == IndexOrdering::forward) {
        pos = domain.x0;
    } else if (ordering == IndexOrdering::backward) {
        pos = domain.z0;
    } else {
        // Sparse ordering not implemented.
        PLB_ASSERT(false);
    }
}

MultiBlockFastSerializer3D *MultiBlockFastSerializer3D::clone() const
{
    return new MultiBlockFastSerializer3D(*this);
}

pluint MultiBlockFastSerializer3D::getSize() const
{
    return domain.nCells() * multiBlock.sizeOfCell();
}

const char *MultiBlockFastSerializer3D::getNextDataBuffer(pluint &bufferSize) const
{
    bufferSize = computeSlice();
    if (buffer.empty())
        buffer.resize(1);
    return &buffer[0];
}

bool MultiBlockFastSerializer3D::isEmpty() const
{
    if (ordering == IndexOrdering::forward) {
        return pos > domain.x1;
    } else {
        return pos > domain.z1;
    }
}

pluint MultiBlockFastSerializer3D::computeSlice() const
{
    Box3D slice(domain);
    if (ordering == IndexOrdering::forward) {
        slice.x0 = slice.x1 = pos;
    } else if (ordering == IndexOrdering::backward) {
        slice.z0 = slice.z1 = pos;
    } else {
        // Sparse ordering not implemented.
        PLB_ASSERT(false);
    }
    SparseBlockStructure3D blockStructure(multiBlock.getBoundingBox());
    blockStructure.addBlock(slice, 0);
    plint envelopeWidth = 1;
    MultiBlockManagement3D serialMultiBlockManagement(
        blockStructure, new OneToOneThreadAttribution, envelopeWidth);
    MultiBlock3D *multiSerialBlock = multiBlock.clone(serialMultiBlockManagement);

    pluint bufferSize = slice.nCells() * multiBlock.sizeOfCell();
    if (global::mpi().isMainProcessor()) {
        AtomicBlock3D &atomicSerialBlock = multiSerialBlock->getComponent(0);
        if (ordering == IndexOrdering::forward) {
            SmartBulk3D oneBlockSlice(serialMultiBlockManagement, 0);
            Box3D localSlice(oneBlockSlice.toLocal(slice));
            atomicSerialBlock.getDataTransfer().send(localSlice, buffer, modif::staticVariables);
        } else {
            buffer.clear();
            std::vector<char> tmpBuf;
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                Box3D subSlice(domain.x0, domain.x1, iY, iY, pos, pos);
                SmartBulk3D oneBlockSlice(serialMultiBlockManagement, 0);
                Box3D localSlice(oneBlockSlice.toLocal(subSlice));
                atomicSerialBlock.getDataTransfer().send(
                    localSlice, tmpBuf, modif::staticVariables);
                buffer.insert(buffer.end(), tmpBuf.begin(), tmpBuf.end());
            }
        }
        PLB_ASSERT(bufferSize == buffer.size());
    }
    delete multiSerialBlock;
    ++pos;
    return bufferSize;
}

MultiBlockFastUnSerializer3D::MultiBlockFastUnSerializer3D(
    MultiBlock3D &multiBlock_, IndexOrdering::OrderingT ordering_) :
    multiBlock(multiBlock_), ordering(ordering_), domain(multiBlock.getBoundingBox())
{
    if (ordering == IndexOrdering::forward) {
        pos = domain.x0;
    } else if (ordering == IndexOrdering::backward) {
        pos = domain.z0;
    } else {
        // Sparse ordering not implemented.
        PLB_ASSERT(false);
    }
}

MultiBlockFastUnSerializer3D::MultiBlockFastUnSerializer3D(
    MultiBlock3D &multiBlock_, Box3D domain_, IndexOrdering::OrderingT ordering_) :
    multiBlock(multiBlock_), ordering(ordering_), domain(domain_)
{
    if (ordering == IndexOrdering::forward) {
        pos = domain.x0;
    } else if (ordering == IndexOrdering::backward) {
        pos = domain.z0;
    } else {
        // Sparse ordering not implemented.
        PLB_ASSERT(false);
    }
}

MultiBlockFastUnSerializer3D *MultiBlockFastUnSerializer3D::clone() const
{
    return new MultiBlockFastUnSerializer3D(*this);
}

pluint MultiBlockFastUnSerializer3D::getSize() const
{
    return domain.nCells() * multiBlock.sizeOfCell();
}

char *MultiBlockFastUnSerializer3D::getNextDataBuffer(pluint &bufferSize)
{
    if (ordering == IndexOrdering::forward) {
        bufferSize = domain.getNy() * domain.getNz() * multiBlock.sizeOfCell();
    } else {
        bufferSize = domain.getNy() * domain.getNx() * multiBlock.sizeOfCell();
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

void MultiBlockFastUnSerializer3D::commitData()
{
    Box3D slice(domain);
    if (ordering == IndexOrdering::forward) {
        slice.x0 = slice.x1 = pos;
    } else if (ordering == IndexOrdering::backward) {
        slice.z0 = slice.z1 = pos;
    } else {
        // Sparse ordering not implemented.
        PLB_ASSERT(false);
    }
    SparseBlockStructure3D blockStructure(multiBlock.getBoundingBox());
    blockStructure.addBlock(slice, 0);
    plint envelopeWidth = 1;
    MultiBlockManagement3D serialMultiBlockManagement(
        blockStructure, new OneToOneThreadAttribution, envelopeWidth);
    MultiBlock3D *multiSerialBlock = multiBlock.clone(serialMultiBlockManagement);
    if (global::mpi().isMainProcessor()) {
        AtomicBlock3D &atomicSerialBlock = multiSerialBlock->getComponent(0);
        if (ordering == IndexOrdering::forward) {
            SmartBulk3D oneBlockBulk(serialMultiBlockManagement, 0);
            Box3D localSlice(oneBlockBulk.toLocal(slice));
            atomicSerialBlock.getDataTransfer().receive(localSlice, buffer, modif::staticVariables);
        } else {
            plint sizeOfSubSlice = domain.getNx() * multiBlock.sizeOfCell();
            std::vector<char> subBuffer(sizeOfSubSlice);
            plint subPos = 0;
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                std::copy(
                    buffer.begin() + subPos, buffer.begin() + subPos + sizeOfSubSlice,
                    subBuffer.begin());
                subPos += sizeOfSubSlice;
                Box3D subSlice(domain.x0, domain.x1, iY, iY, pos, pos);
                SmartBulk3D oneBlockBulk(serialMultiBlockManagement, 0);
                Box3D localSubSlice(oneBlockBulk.toLocal(subSlice));
                atomicSerialBlock.getDataTransfer().receive(
                    localSubSlice, subBuffer, modif::staticVariables);
            }
        }
    }
    multiBlock.copyReceive(*multiSerialBlock, slice, slice);
    delete multiSerialBlock;
    ++pos;
}

bool MultiBlockFastUnSerializer3D::isFull() const
{
    if (ordering == IndexOrdering::forward) {
        return pos > domain.x1;
    } else {
        return pos > domain.z1;
    }
}

}  //  namespace plb
