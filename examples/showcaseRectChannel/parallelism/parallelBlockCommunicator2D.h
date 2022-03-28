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
 * Helper classes for parallel 2D multiblock lattice -- header file.
 */

#ifndef PARALLEL_BLOCK_COMMUNICATOR_2D_H
#define PARALLEL_BLOCK_COMMUNICATOR_2D_H

#include <vector>

#include "core/globalDefs.h"
#include "multiBlock/blockCommunicator2D.h"
#include "multiBlock/multiBlock2D.h"
#include "multiBlock/multiBlockManagement2D.h"
#include "parallelism/communicationPackage2D.h"
#include "parallelism/sendRecvPool.h"

namespace plb {

#ifdef PLB_MPI_PARALLEL

struct CommunicationStructure2D {
    CommunicationStructure2D(
        std::vector<Overlap2D> const &overlaps, MultiBlockManagement2D const &originManagement,
        MultiBlockManagement2D const &destinationManagement, plint sizeOfCell);
    CommunicationPackage2D sendPackage;
    CommunicationPackage2D recvPackage;
    CommunicationPackage2D sendRecvPackage;
    SendPoolCommunicator sendComm;
    RecvPoolCommunicator recvComm;
};

class ParallelBlockCommunicator2D : public BlockCommunicator2D {
public:
    ParallelBlockCommunicator2D();
    ~ParallelBlockCommunicator2D();
    ParallelBlockCommunicator2D(ParallelBlockCommunicator2D const &rhs);
    ParallelBlockCommunicator2D &operator=(ParallelBlockCommunicator2D const &rhs);
    void swap(ParallelBlockCommunicator2D &rhs);
    virtual ParallelBlockCommunicator2D *clone() const;
    virtual void duplicateOverlaps(MultiBlock2D &multiBlock, modif::ModifT whichData) const;
    virtual void communicate(
        std::vector<Overlap2D> const &overlaps, MultiBlock2D const &originMultiBlock,
        MultiBlock2D &destinationMultiBlock, modif::ModifT whichData) const;
    virtual void signalPeriodicity() const;

private:
    void communicate(
        CommunicationStructure2D &communication, MultiBlock2D const &originMultiBlock,
        MultiBlock2D &destinationMultiBlock, modif::ModifT whichData) const;
    void subscribeOverlap(
        Overlap2D const &overlap, MultiBlockManagement2D const &multiBlockManagement,
        SendRecvPool &sendPool, SendRecvPool &recvPool, plint sizeOfCell) const;

private:
    mutable bool overlapsModified;
    mutable CommunicationStructure2D *communication;
};

#endif  // PLB_MPI_PARALLEL

}  // namespace plb

#endif  // PARALLEL_BLOCK_COMMUNICATOR_2D_H
