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
 * Geometry specifications for 3D multiblocks -- header file.
 */

#ifdef PLB_MPI_PARALLEL

#include "multiBlock/nonLocalTransfer2D_3D.h"

#include "core/globalDefs.h"
#include "core/plbProfiler.h"

namespace plb {

Box3D to3d(Box2D const &box)
{
    return Box3D(box.x0, box.x1, box.y0, box.y1, 0, 0);
}

Box2D to2d(Box3D const &box)
{
    return Box2D(box.x0, box.x1, box.y0, box.y1);
}

Dot2D to2d(Dot3D const &dot3d)
{
    return Dot2D(dot3d.x, dot3d.y);
}

std::vector<Overlap3D> copyDomainDataTransfer(
    SparseBlockStructure2D const &block1, Box2D block1Domain, SparseBlockStructure3D const &block2,
    Box3D block2Domain)
{
    PLB_PRECONDITION(block1Domain.getNx() == block2Domain.getNx());
    PLB_PRECONDITION(block1Domain.getNy() == block2Domain.getNy());
    PLB_PRECONDITION(block2Domain.getNz() == 1);
    plint shiftX = block1Domain.x0 - block2Domain.x0;
    plint shiftY = block1Domain.y0 - block2Domain.y0;
    plint shiftZ = 0 - block2Domain.z0;
    std::vector<plint> block1Ids, block2Ids;
    std::vector<Box2D> block1Components;
    std::vector<Box3D> block2Inters;
    std::vector<Overlap3D> dataTransfer;  // The return value.
    block1.intersect(block1Domain, block1Ids, block1Components);
    PLB_ASSERT(block1Ids.size() == block1Components.size());
    for (pluint iComp1 = 0; iComp1 < block1Ids.size(); ++iComp1) {
        block2Ids.clear();
        block2Inters.clear();
        block2.intersect(
            to3d(block1Components[iComp1]).shift(-shiftX, -shiftY, -shiftZ), block2Ids,
            block2Inters);
        for (pluint iInters = 0; iInters < block2Inters.size(); ++iInters) {
            dataTransfer.push_back(Overlap3D(
                block1Ids[iComp1], block2Ids[iInters],
                block2Inters[iInters].shift(shiftX, shiftY, shiftZ), shiftX, shiftY, shiftZ));
        }
    }

    return dataTransfer;
}

std::vector<Overlap3D> copyDomainDataTransfer(
    SparseBlockStructure3D const &block1, Box3D block1Domain, SparseBlockStructure2D const &block2,
    Box2D block2Domain)
{
    PLB_PRECONDITION(block1Domain.getNx() == block2Domain.getNx());
    PLB_PRECONDITION(block1Domain.getNy() == block2Domain.getNy());
    PLB_PRECONDITION(block1Domain.getNz() == 1);
    plint shiftX = block1Domain.x0 - block2Domain.x0;
    plint shiftY = block1Domain.y0 - block2Domain.y0;
    plint shiftZ = block1Domain.z0;
    std::vector<plint> block1Ids, block2Ids;
    std::vector<Box3D> block1Components;
    std::vector<Box2D> block2Inters;
    std::vector<Overlap3D> dataTransfer;  // The return value.
    block1.intersect(block1Domain, block1Ids, block1Components);
    PLB_ASSERT(block1Ids.size() == block1Components.size());
    for (pluint iComp1 = 0; iComp1 < block1Ids.size(); ++iComp1) {
        block2Ids.clear();
        block2Inters.clear();
        block2.intersect(
            to2d(block1Components[iComp1]).shift(-shiftX, -shiftY), block2Ids, block2Inters);
        for (pluint iInters = 0; iInters < block2Inters.size(); ++iInters) {
            dataTransfer.push_back(Overlap3D(
                block1Ids[iComp1], block2Ids[iInters],
                to3d(block2Inters[iInters]).shift(shiftX, shiftY, shiftZ), shiftX, shiftY, shiftZ));
        }
    }

    return dataTransfer;
}

CommunicationStructure2D_3D::CommunicationStructure2D_3D(
    std::vector<Overlap3D> const &overlaps, MultiBlockManagement2D const &originManagement,
    MultiBlockManagement3D const &destinationManagement, plint sizeOfCell)
{
    plint fromEnvelopeWidth = originManagement.getEnvelopeWidth();
    plint toEnvelopeWidth = destinationManagement.getEnvelopeWidth();
    SparseBlockStructure2D const &fromSparseBlock = originManagement.getSparseBlockStructure();
    SparseBlockStructure3D const &toSparseBlock = destinationManagement.getSparseBlockStructure();

    SendRecvPool sendPool, recvPool;
    for (pluint iOverlap = 0; iOverlap < overlaps.size(); ++iOverlap) {
        Overlap3D const &overlap = overlaps[iOverlap];
        CommunicationInfo3D info;

        info.fromBlockId = overlap.getOriginalId();
        info.toBlockId = overlap.getOverlapId();

        SmartBulk2D originalBulk(fromSparseBlock, fromEnvelopeWidth, info.fromBlockId);
        SmartBulk3D overlapBulk(toSparseBlock, toEnvelopeWidth, info.toBlockId);

        Box3D originalCoordinates(overlap.getOriginalCoordinates());
        Box3D overlapCoordinates(overlap.getOverlapCoordinates());
        info.fromDomain = to3d(originalBulk.toLocal(to2d(originalCoordinates)));
        info.toDomain = overlapBulk.toLocal(overlapCoordinates);
        info.absoluteOffset = Dot3D(
            overlapCoordinates.x0 - originalCoordinates.x0,
            overlapCoordinates.y0 - originalCoordinates.y0,
            overlapCoordinates.z0 - originalCoordinates.z0);

        plint lx = info.fromDomain.x1 - info.fromDomain.x0 + 1;
        plint ly = info.fromDomain.y1 - info.fromDomain.y0 + 1;
        plint lz = info.fromDomain.z1 - info.fromDomain.z0 + 1;
        PLB_PRECONDITION(lx == info.toDomain.x1 - info.toDomain.x0 + 1);
        PLB_PRECONDITION(ly == info.toDomain.y1 - info.toDomain.y0 + 1);
        PLB_PRECONDITION(lz == info.toDomain.z1 - info.toDomain.z0 + 1);

        plint numberOfCells = lx * ly * lz;

        ThreadAttribution const &fromAttribution = originManagement.getThreadAttribution();
        ThreadAttribution const &toAttribution = destinationManagement.getThreadAttribution();
        info.fromProcessId = fromAttribution.getMpiProcess(info.fromBlockId);
        info.toProcessId = toAttribution.getMpiProcess(info.toBlockId);

        if (fromAttribution.isLocal(info.fromBlockId) && toAttribution.isLocal(info.toBlockId)) {
            sendRecvPackage.push_back(info);
        } else if (fromAttribution.isLocal(info.fromBlockId)) {
            sendPackage.push_back(info);
            sendPool.subscribeMessage(info.toProcessId, numberOfCells * sizeOfCell);
        } else if (toAttribution.isLocal(info.toBlockId)) {
            recvPackage.push_back(info);
            recvPool.subscribeMessage(info.fromProcessId, numberOfCells * sizeOfCell);
        }
    }

    sendComm = SendPoolCommunicator(sendPool);
    recvComm = RecvPoolCommunicator(recvPool);
}

CommunicationStructure3D_2D::CommunicationStructure3D_2D(
    std::vector<Overlap3D> const &overlaps, MultiBlockManagement3D const &originManagement,
    MultiBlockManagement2D const &destinationManagement, plint sizeOfCell)
{
    plint fromEnvelopeWidth = originManagement.getEnvelopeWidth();
    plint toEnvelopeWidth = destinationManagement.getEnvelopeWidth();
    SparseBlockStructure3D const &fromSparseBlock = originManagement.getSparseBlockStructure();
    SparseBlockStructure2D const &toSparseBlock = destinationManagement.getSparseBlockStructure();

    SendRecvPool sendPool, recvPool;
    for (pluint iOverlap = 0; iOverlap < overlaps.size(); ++iOverlap) {
        Overlap3D const &overlap = overlaps[iOverlap];
        CommunicationInfo3D info;

        info.fromBlockId = overlap.getOriginalId();
        info.toBlockId = overlap.getOverlapId();

        SmartBulk3D originalBulk(fromSparseBlock, fromEnvelopeWidth, info.fromBlockId);
        SmartBulk2D overlapBulk(toSparseBlock, toEnvelopeWidth, info.toBlockId);

        Box3D originalCoordinates(overlap.getOriginalCoordinates());
        Box3D overlapCoordinates(overlap.getOverlapCoordinates());
        info.fromDomain = originalBulk.toLocal(originalCoordinates);
        info.toDomain = to3d(overlapBulk.toLocal(to2d(overlapCoordinates)));
        info.absoluteOffset = Dot3D(
            overlapCoordinates.x0 - originalCoordinates.x0,
            overlapCoordinates.y0 - originalCoordinates.y0,
            overlapCoordinates.z0 - originalCoordinates.z0);

        plint lx = info.fromDomain.x1 - info.fromDomain.x0 + 1;
        plint ly = info.fromDomain.y1 - info.fromDomain.y0 + 1;
        plint lz = info.fromDomain.z1 - info.fromDomain.z0 + 1;
        PLB_PRECONDITION(lx == info.toDomain.x1 - info.toDomain.x0 + 1);
        PLB_PRECONDITION(ly == info.toDomain.y1 - info.toDomain.y0 + 1);
        PLB_PRECONDITION(lz == info.toDomain.z1 - info.toDomain.z0 + 1);

        plint numberOfCells = lx * ly * lz;

        ThreadAttribution const &fromAttribution = originManagement.getThreadAttribution();
        ThreadAttribution const &toAttribution = destinationManagement.getThreadAttribution();
        info.fromProcessId = fromAttribution.getMpiProcess(info.fromBlockId);
        info.toProcessId = toAttribution.getMpiProcess(info.toBlockId);

        if (fromAttribution.isLocal(info.fromBlockId) && toAttribution.isLocal(info.toBlockId)) {
            sendRecvPackage.push_back(info);
        } else if (fromAttribution.isLocal(info.fromBlockId)) {
            sendPackage.push_back(info);
            sendPool.subscribeMessage(info.toProcessId, numberOfCells * sizeOfCell);
        } else if (toAttribution.isLocal(info.toBlockId)) {
            recvPackage.push_back(info);
            recvPool.subscribeMessage(info.fromProcessId, numberOfCells * sizeOfCell);
        }
    }

    sendComm = SendPoolCommunicator(sendPool);
    recvComm = RecvPoolCommunicator(recvPool);
}

void communicate(
    CommunicationStructure2D_3D &communication, MultiBlock2D const &originMultiBlock,
    MultiBlock3D &destinationMultiBlock, modif::ModifT whichData)
{
    global::profiler().start("mpiCommunication");
    bool staticMessage = whichData == modif::staticVariables;
    // 1. Non-blocking receives.
    communication.recvComm.startBeingReceptive(staticMessage);

    // 2. Non-blocking sends.
    for (unsigned iSend = 0; iSend < communication.sendPackage.size(); ++iSend) {
        CommunicationInfo3D const &info = communication.sendPackage[iSend];
        AtomicBlock2D const &fromBlock = originMultiBlock.getComponent(info.fromBlockId);
        fromBlock.getDataTransfer().send(
            to2d(info.fromDomain), communication.sendComm.getSendBuffer(info.toProcessId),
            whichData);
        communication.sendComm.acceptMessage(info.toProcessId, staticMessage);
    }

    // 3. Local copies which require no communication.
    for (unsigned iSendRecv = 0; iSendRecv < communication.sendRecvPackage.size(); ++iSendRecv) {
        CommunicationInfo3D const &info = communication.sendRecvPackage[iSendRecv];
        AtomicBlock2D const &fromBlock = originMultiBlock.getComponent(info.fromBlockId);
        AtomicBlock3D &toBlock = destinationMultiBlock.getComponent(info.toBlockId);
        std::vector<char> tmpBuf;
        fromBlock.getDataTransfer().send(to2d(info.fromDomain), tmpBuf, whichData);
        toBlock.getDataTransfer().receive(info.toDomain, tmpBuf, whichData);
    }

    // 4. Finalize the receives.
    for (unsigned iRecv = 0; iRecv < communication.recvPackage.size(); ++iRecv) {
        CommunicationInfo3D const &info = communication.recvPackage[iRecv];
        AtomicBlock3D &toBlock = destinationMultiBlock.getComponent(info.toBlockId);
        toBlock.getDataTransfer().receive(
            info.toDomain, communication.recvComm.receiveMessage(info.fromProcessId, staticMessage),
            whichData, info.absoluteOffset);
    }

    // 5. Finalize the sends.
    communication.sendComm.finalize(staticMessage);
    global::profiler().stop("mpiCommunication");
}

void communicate(
    std::vector<Overlap3D> const &overlaps, MultiBlock2D const &originMultiBlock,
    MultiBlock3D &destinationMultiBlock, modif::ModifT whichData)
{
    PLB_PRECONDITION(originMultiBlock.sizeOfCell() == destinationMultiBlock.sizeOfCell());

    CommunicationStructure2D_3D communication(
        overlaps, originMultiBlock.getMultiBlockManagement(),
        destinationMultiBlock.getMultiBlockManagement(), originMultiBlock.sizeOfCell());
    communicate(communication, originMultiBlock, destinationMultiBlock, whichData);
}

void copy_generic(
    MultiBlock2D const &from, Box2D const &fromDomain, MultiBlock3D &to, Box3D const &toDomain,
    modif::ModifT typeOfModif)
{
    Box3D fromDomain_(to3d(fromDomain));
    Box3D toDomain_(toDomain);
    adjustEqualSize(fromDomain_, toDomain_);
    std::vector<Overlap3D> dataTransfer = copyDomainDataTransfer(
        from.getMultiBlockManagement().getSparseBlockStructure(), to2d(fromDomain_),
        to.getMultiBlockManagement().getSparseBlockStructure(), toDomain_);
    communicate(dataTransfer, from, to, typeOfModif);
    to.getBlockCommunicator().duplicateOverlaps(to, typeOfModif);
}

void communicate(
    CommunicationStructure3D_2D &communication, MultiBlock3D const &originMultiBlock,
    MultiBlock2D &destinationMultiBlock, modif::ModifT whichData)
{
    global::profiler().start("mpiCommunication");
    bool staticMessage = whichData == modif::staticVariables;
    // 1. Non-blocking receives.
    communication.recvComm.startBeingReceptive(staticMessage);

    // 2. Non-blocking sends.
    for (unsigned iSend = 0; iSend < communication.sendPackage.size(); ++iSend) {
        CommunicationInfo3D const &info = communication.sendPackage[iSend];
        AtomicBlock3D const &fromBlock = originMultiBlock.getComponent(info.fromBlockId);
        fromBlock.getDataTransfer().send(
            info.fromDomain, communication.sendComm.getSendBuffer(info.toProcessId), whichData);
        communication.sendComm.acceptMessage(info.toProcessId, staticMessage);
    }

    // 3. Local copies which require no communication.
    for (unsigned iSendRecv = 0; iSendRecv < communication.sendRecvPackage.size(); ++iSendRecv) {
        CommunicationInfo3D const &info = communication.sendRecvPackage[iSendRecv];
        AtomicBlock3D const &fromBlock = originMultiBlock.getComponent(info.fromBlockId);
        AtomicBlock2D &toBlock = destinationMultiBlock.getComponent(info.toBlockId);
        std::vector<char> tmpBuf;
        fromBlock.getDataTransfer().send(info.fromDomain, tmpBuf, whichData);
        toBlock.getDataTransfer().receive(to2d(info.toDomain), tmpBuf, whichData);
    }

    // 4. Finalize the receives.
    for (unsigned iRecv = 0; iRecv < communication.recvPackage.size(); ++iRecv) {
        CommunicationInfo3D const &info = communication.recvPackage[iRecv];
        AtomicBlock2D &toBlock = destinationMultiBlock.getComponent(info.toBlockId);
        toBlock.getDataTransfer().receive(
            to2d(info.toDomain),
            communication.recvComm.receiveMessage(info.fromProcessId, staticMessage), whichData,
            to2d(info.absoluteOffset));
    }

    // 5. Finalize the sends.
    communication.sendComm.finalize(staticMessage);
    global::profiler().stop("mpiCommunication");
}

void communicate(
    std::vector<Overlap3D> const &overlaps, MultiBlock3D const &originMultiBlock,
    MultiBlock2D &destinationMultiBlock, modif::ModifT whichData)
{
    PLB_PRECONDITION(originMultiBlock.sizeOfCell() == destinationMultiBlock.sizeOfCell());

    CommunicationStructure3D_2D communication(
        overlaps, originMultiBlock.getMultiBlockManagement(),
        destinationMultiBlock.getMultiBlockManagement(), originMultiBlock.sizeOfCell());
    communicate(communication, originMultiBlock, destinationMultiBlock, whichData);
}

void copy_generic(
    MultiBlock3D const &from, Box3D const &fromDomain, MultiBlock2D &to, Box2D const &toDomain,
    modif::ModifT typeOfModif)
{
    Box3D fromDomain_(fromDomain);
    Box3D toDomain_(to3d(toDomain));
    adjustEqualSize(fromDomain_, toDomain_);
    std::vector<Overlap3D> dataTransfer = copyDomainDataTransfer(
        from.getMultiBlockManagement().getSparseBlockStructure(), fromDomain_,
        to.getMultiBlockManagement().getSparseBlockStructure(), to2d(toDomain_));
    communicate(dataTransfer, from, to, typeOfModif);
    to.getBlockCommunicator().duplicateOverlaps(to, typeOfModif);
}

}  // namespace plb

#endif  // PLB_MPI_PARALLEL
