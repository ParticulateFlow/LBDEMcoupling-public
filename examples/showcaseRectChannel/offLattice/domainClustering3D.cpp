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

#include "offLattice/domainClustering3D.h"

#include <algorithm>
#include <iterator>
#include <numeric>

#include "core/globalDefs.h"
#include "io/parallelIO.h"
#include "parallelism/mpiManager.h"

namespace plb {

#ifdef PLB_MPI_PARALLEL

class BlockWithContact {
public:
    BlockWithContact() : blockId(0), numContacts(0) { }
    BlockWithContact(plint blockId_, plint numContacts_) :
        blockId(blockId_), numContacts(numContacts_)
    { }
    bool operator<(BlockWithContact const &rhs) const
    {
        return numContacts < rhs.numContacts;
    }
    plint getBlockId() const
    {
        return blockId;
    }

private:
    plint blockId, numContacts;
};

ExplicitThreadAttribution *optimalThreadAttribution(
    SparseBlockStructure3D const &sparseBlock, ThreadAttribution const &oldThreadAttribution)
{
    std::map<plint, Box3D> const &domains = sparseBlock.getBulks();
    std::map<plint, Box3D>::const_iterator it = domains.begin();
    std::vector<plint> localBlocks;
    for (; it != domains.end(); ++it) {
        plint id = it->first;
        if (oldThreadAttribution.isLocal(id)) {
            localBlocks.push_back(id);
        }
    }

    std::vector<plint> optimalLocalBlocks = optimalRepartition(sparseBlock, localBlocks);
    std::vector<plint> numBlocks(global::mpi().getSize(), 0);
    numBlocks[global::mpi().getRank()] = (plint)optimalLocalBlocks.size();
    global::mpi().allReduceVect(numBlocks, MPI_SUM);

    std::vector<plint> sumNumBlocks(numBlocks.size() + 1);
    sumNumBlocks[0] = 0;
    std::partial_sum(numBlocks.begin(), numBlocks.end(), sumNumBlocks.begin() + 1);
    plint totalNumBlocks = sumNumBlocks.back();

    std::vector<plint> everyOnesBlocks(totalNumBlocks, 0);
    plint startBlock = sumNumBlocks[global::mpi().getRank()];
    for (pluint iBlock = 0; iBlock < optimalLocalBlocks.size(); ++iBlock) {
        everyOnesBlocks[startBlock + iBlock] = optimalLocalBlocks[iBlock];
    }
    global::mpi().allReduceVect(everyOnesBlocks, MPI_SUM);

    ExplicitThreadAttribution *newThreadAttribution = new ExplicitThreadAttribution;
    for (plint iThread = 0; iThread < global::mpi().getSize(); ++iThread) {
        plint startBlock = sumNumBlocks[iThread];
        plint endBlock = sumNumBlocks[iThread + 1];
        for (plint iBlock = startBlock; iBlock < endBlock; ++iBlock) {
            newThreadAttribution->addBlock(everyOnesBlocks[iBlock], iThread);
        }
    }
    return newThreadAttribution;
}

std::vector<plint> optimalRepartition(
    SparseBlockStructure3D const &sparseBlock, std::vector<plint> localBlocks)
{
    plint numIter = 6;
    std::vector<plint> selectedOptimalRepartition;
    plint bestFitness = std::numeric_limits<plint>::min();
    for (plint iter = 0; iter < numIter; ++iter) {
        std::vector<plint> nextOptimalRepartition = tryOptimalRepartition(sparseBlock, localBlocks);

        plint numLocalNeighbors = getNumLocalNeighbors(
            sparseBlock,
            std::set<plint>(nextOptimalRepartition.begin(), nextOptimalRepartition.end()));
        global::mpi().reduceAndBcast(numLocalNeighbors, MPI_SUM);
        pcout << "num local neighbors=" << numLocalNeighbors << std::endl;
        if (numLocalNeighbors > bestFitness) {
            bestFitness = numLocalNeighbors;
            selectedOptimalRepartition = nextOptimalRepartition;
        }
    }
    return selectedOptimalRepartition;
}

std::vector<plint> tryOptimalRepartition(
    SparseBlockStructure3D const &sparseBlock, std::vector<plint> localBlocks)
{
    // Switch to a set representation to find out quickly if a given neighbor is local or not.
    std::set<plint> localBlockSet(localBlocks.begin(), localBlocks.end());
    if (global::mpi().isMainProcessor()) {
        return mainProcessorRepartition(sparseBlock, localBlockSet);
    }
    const plint maxIter = 2 * sparseBlock.getNumBlocks();
    for (plint iter = 0; iter < maxIter; ++iter) {
        plint myPartner;
        global::mpi().receive(&myPartner, 1, 0);
        bool exchangeAnyway = iter < maxIter / 3;
        repartitionIter(sparseBlock, localBlockSet, myPartner, exchangeAnyway);
    }
    return std::vector<plint>(localBlockSet.begin(), localBlockSet.end());
}

std::vector<plint> mainProcessorRepartition(
    SparseBlockStructure3D const &sparseBlock, std::set<plint> localBlocks)
{
    const plint maxIter = 2 * sparseBlock.getNumBlocks();
    std::vector<int> processIds(global::mpi().getSize());
    for (int i = 0; i < global::mpi().getSize(); ++i) {
        processIds[i] = i;
    }
    for (plint iter = 0; iter < maxIter; ++iter) {
        std::vector<int> shuffledIds(processIds);
        std::random_shuffle(shuffledIds.begin(), shuffledIds.end());
        int myPartner = 0;
        for (int i = 0; i < global::mpi().getSize() - 1; i += 2) {
            int partner1 = shuffledIds[i];
            int partner2 = shuffledIds[i + 1];
            if (partner1 == global::mpi().bossId()) {
                myPartner = partner2;
            } else {
                global::mpi().send(&partner2, 1, partner1);
            }
            if (partner2 == global::mpi().bossId()) {
                myPartner = partner1;
            } else {
                global::mpi().send(&partner1, 1, partner2);
            }
        }
        bool exchangeAnyway = iter < maxIter / 3;
        repartitionIter(sparseBlock, localBlocks, myPartner, exchangeAnyway);
    }
    return std::vector<plint>(localBlocks.begin(), localBlocks.end());
}

void repartitionIter(
    SparseBlockStructure3D const &sparseBlock, std::set<plint> &localBlocks, int myPartner,
    bool exchangeAnyway)
{
    double epsilon = 0.1;
    if (exchangeAnyway) {
        epsilon = 1.0;
    }
    std::vector<BlockWithContact> contacts;
    plint totalContactsBefore;
    evaluateNeighbors(sparseBlock, localBlocks, contacts, totalContactsBefore);
    std::sort(contacts.begin(), contacts.end());

    plint numConsideredBlocks = (plint)(localBlocks.size() * epsilon);
    if (numConsideredBlocks == 0)
        numConsideredBlocks = 1;
    plint choice = (plint)((double)rand() / ((double)RAND_MAX + 1.) * (double)numConsideredBlocks);
    plint sacrificedBlock = contacts[choice].getBlockId();

    plint receivedBlock;
    global::mpi().sendRecv(&sacrificedBlock, &receivedBlock, 1, myPartner, myPartner);
    plint partnerTotalContactsBefore;
    global::mpi().sendRecv(
        &totalContactsBefore, &partnerTotalContactsBefore, 1, myPartner, myPartner);

    plint totalContactsAfter, partnerTotalContactsAfter;
    localBlocks.erase(sacrificedBlock);
    localBlocks.insert(receivedBlock);
    evaluateNeighbors(sparseBlock, localBlocks, contacts, totalContactsAfter);
    global::mpi().sendRecv(
        &totalContactsAfter, &partnerTotalContactsAfter, 1, myPartner, myPartner);
    // The exchange didn't offer any advantage; let's switch back.
    if (totalContactsAfter < totalContactsBefore
        && partnerTotalContactsAfter < partnerTotalContactsBefore && !exchangeAnyway)
    {
        localBlocks.erase(receivedBlock);
        localBlocks.insert(sacrificedBlock);
    }
}

void evaluateNeighbors(
    SparseBlockStructure3D const &sparseBlock, std::set<plint> const &localBlocks,
    std::vector<BlockWithContact> &contacts, plint &totalContacts)
{
    contacts.resize(localBlocks.size());
    totalContacts = 0;

    std::set<plint>::const_iterator it = localBlocks.begin();
    for (plint iBlock = 0; it != localBlocks.end(); ++it, ++iBlock) {
        plint blockId = *it;
        plint numContacts = getNumLocalNeighbors(sparseBlock, localBlocks, blockId);
        totalContacts += numContacts;
        contacts[iBlock] = BlockWithContact(blockId, numContacts);
    }
}

plint getNumLocalNeighbors(
    SparseBlockStructure3D const &sparseBlock, std::set<plint> const &localBlocks)
{
    std::set<plint>::const_iterator it = localBlocks.begin();
    plint totalContacts = 0;
    for (; it != localBlocks.end(); ++it) {
        plint blockId = *it;
        plint numContacts = getNumLocalNeighbors(sparseBlock, localBlocks, blockId);
        totalContacts += numContacts;
    }
    return totalContacts;
}

plint getNumLocalNeighbors(
    SparseBlockStructure3D const &sparseBlock, std::set<plint> const &localBlocks, plint blockId)
{
    Box3D myBulk;
    sparseBlock.getBulk(blockId, myBulk);
    myBulk = myBulk.enlarge(1);
    std::vector<plint> neighbors;
    sparseBlock.findNeighbors(blockId, 1, neighbors);
    plint numContacts = 0;
    for (pluint i = 0; i < neighbors.size(); ++i) {
        plint neighborBlock = neighbors[i];
        std::set<plint>::const_iterator neighbIt = localBlocks.find(neighborBlock);
        if (neighbIt != localBlocks.end()) {
            Box3D neighborBulk, intersection;
            sparseBlock.getBulk(neighborBlock, neighborBulk);
            if (intersect(myBulk, neighborBulk, intersection)) {
                numContacts += intersection.nCells();
            }
        }
    }
    return numContacts;
}

#endif  // PLB_MPI_PARALLEL

}  // namespace plb
