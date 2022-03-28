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
 * Wrapper functions that simplify the use of MPI
 */

#include "parallelism/sendRecvPool.h"

#include <numeric>

#include "core/plbDebug.h"
#include "core/plbProfiler.h"
#include "parallelism/mpiManager.h"

namespace plb {

#ifdef PLB_MPI_PARALLEL

SendPoolCommunicator::SendPoolCommunicator(SendRecvPool const &pool) :
    subscriptions(pool.begin(), pool.end())
{
    // PLB_PRECONDITION(!pool.empty());
}

std::vector<char> &SendPoolCommunicator::getSendBuffer(int toProc)
{
    std::map<int, CommunicatorEntry>::iterator entryPtr = subscriptions.find(toProc);
    PLB_ASSERT(entryPtr != subscriptions.end());
    CommunicatorEntry &entry = entryPtr->second;
    PLB_ASSERT(entry.currentMessage < (int)entry.messages.size());
    return entry.messages[entry.currentMessage];
}

void SendPoolCommunicator::acceptMessage(int toProc, bool staticMessage)
{
    std::map<int, CommunicatorEntry>::iterator entryPtr = subscriptions.find(toProc);
    PLB_ASSERT(entryPtr != subscriptions.end());
    CommunicatorEntry &entry = entryPtr->second;
    PLB_ASSERT(entry.currentMessage < (int)entry.messages.size());
    // If communication is static, make sure that the message has
    //   the right size.
    PLB_ASSERT(
        !staticMessage
        || ((int)entry.messages[entry.currentMessage].size()
            == entry.lengths[entry.currentMessage]));
    entry.currentMessage++;

    if (entry.currentMessage == (int)entry.lengths.size()) {
        startCommunication(toProc, staticMessage);
        entry.reset();
    }
}

void SendPoolCommunicator::finalize(bool staticMessage)
{
    // PLB_ASSERT( !subscriptions.empty() );
    std::map<int, CommunicatorEntry>::iterator iter = subscriptions.begin();
    for (; iter != subscriptions.end(); ++iter) {
        CommunicatorEntry &entry = iter->second;
        if (!staticMessage) {
            global::mpi().wait(&entry.sizeRequest, &entry.sizeStatus);
        }
        // Empty messages are neither sent nor received.
        if (!entry.data.empty()) {
            global::mpi().wait(&entry.messageRequest, &entry.messageStatus);
        }
    }
}

void SendPoolCommunicator::startCommunication(int toProc, bool staticMessage)
{
    std::map<int, CommunicatorEntry>::iterator entryPtr = subscriptions.find(toProc);
    PLB_ASSERT(entryPtr != subscriptions.end());
    CommunicatorEntry &entry = entryPtr->second;
    if (staticMessage) {
        entry.data.resize(entry.cumDataLength);
    } else {
        // If the communicated data is non-static, the overall size of transmitted
        //   data must be computed.
        int dynamicDataLength = 0;
        entry.dynamicDataSizes.resize(entry.messages.size());
        for (pluint iMessage = 0; iMessage < entry.messages.size(); ++iMessage) {
            dynamicDataLength += entry.messages[iMessage].size();
            entry.dynamicDataSizes[iMessage] = entry.messages[iMessage].size();
        }
        entry.data.resize(dynamicDataLength);
    }
    // Merge the individual messages into a single vector.
    int pos = 0;
    for (pluint iMessage = 0; iMessage < entry.messages.size(); ++iMessage) {
        PLB_ASSERT(
            !staticMessage || ((int)entry.messages[iMessage].size() == entry.lengths[iMessage]));
        PLB_ASSERT(pos + entry.messages[iMessage].size() <= entry.data.size());
        if (!entry.messages[iMessage].empty() && !entry.data.empty()) {
            std::copy(
                entry.messages[iMessage].begin(), entry.messages[iMessage].end(),
                entry.data.begin() + pos);
        }
        pos += entry.messages[iMessage].size();
    }
    if (!staticMessage) {
        PLB_ASSERT(entry.dynamicDataSizes.size() > 0);
        global::profiler().increment("mpiSendChar", (plint)entry.dynamicDataSizes.size());
        global::mpi().iSend(
            &entry.dynamicDataSizes[0], entry.dynamicDataSizes.size(), toProc, &entry.sizeRequest);
    }
    // Empty messages are neither sent nor received.
    if (!entry.data.empty()) {
        global::profiler().increment("mpiSendChar", (plint)entry.data.size());
        global::mpi().iSend(&entry.data[0], entry.data.size(), toProc, &entry.messageRequest);
    }
}

RecvPoolCommunicator::RecvPoolCommunicator(SendRecvPool const &pool) :
    subscriptions(pool.begin(), pool.end())
{ }

void RecvPoolCommunicator::startBeingReceptive(bool staticMessage)
{
    // If the message has dynamic content, the receives cannot be intantiated
    //   at this point, because the message size is unknown. The message size
    //   is being transmitted by MPI communication.
    if (!staticMessage) {
        return;
    }
    std::map<int, CommunicatorEntry>::iterator iter = subscriptions.begin();
    for (; iter != subscriptions.end(); ++iter) {
        int fromProc = iter->first;
        CommunicatorEntry &entry = iter->second;
        entry.data.resize(entry.cumDataLength);
        // Empty messages are neither sent nor received.
        if (!entry.data.empty()) {
            global::profiler().increment("mpiReceiveChar", (plint)entry.data.size());
            global::mpi().iRecv(&entry.data[0], entry.data.size(), fromProc, &entry.messageRequest);
        }
    }
}

std::vector<char> const &RecvPoolCommunicator::receiveMessage(int fromProc, bool staticMessage)
{
    std::map<int, CommunicatorEntry>::iterator entryPtr = subscriptions.find(fromProc);
    PLB_ASSERT(entryPtr != subscriptions.end());
    CommunicatorEntry &entry = entryPtr->second;
    PLB_ASSERT(entry.currentMessage < (int)entry.messages.size());
    if (entry.currentMessage == 0) {
        if (staticMessage) {
            finalizeStatic(fromProc);
        } else {
            receiveDynamic(fromProc);
        }
    }
    std::vector<char> const &message = entry.messages[entry.currentMessage];
    entry.currentMessage++;
    if (entry.currentMessage == (int)entry.lengths.size()) {
        entry.reset();
    }
    return message;
}

void RecvPoolCommunicator::receiveDynamic(int fromProc)
{
    std::map<int, CommunicatorEntry>::iterator entryPtr = subscriptions.find(fromProc);
    PLB_ASSERT(entryPtr != subscriptions.end());
    CommunicatorEntry &entry = entryPtr->second;

    // 1. In a first MPI communication the individual message sizes
    //    are obtained.
    pluint numMessages = entry.messages.size();
    std::vector<int> messageSizes(numMessages);
    PLB_ASSERT(numMessages > 0);
    global::mpi().receive(&messageSizes[0], numMessages, fromProc);

    // 2. All messages are received in a single MPI communication.
    int totalSize = std::accumulate(messageSizes.begin(), messageSizes.end(), 0);
    entry.data.resize(totalSize);
    // Empty messages are neither sent nor received.
    if (!entry.data.empty()) {
        global::profiler().increment("mpiReceiveChar", (plint)totalSize);
        global::mpi().receive(&entry.data[0], totalSize, fromProc);
    }

    // 3. The message package is split into individual messages.
    int pos = 0;
    for (pluint iMessage = 0; iMessage < entry.messages.size(); ++iMessage) {
        int length = messageSizes[iMessage];
        entry.messages[iMessage].resize(length);
        PLB_ASSERT(pos + length <= (int)entry.data.size());
        if (!entry.data.empty() && !entry.messages[iMessage].empty()) {
            std::copy(
                entry.data.begin() + pos, entry.data.begin() + pos + length,
                entry.messages[iMessage].begin());
        }
        pos += length;
    }
}

void RecvPoolCommunicator::finalizeStatic(int fromProc)
{
    std::map<int, CommunicatorEntry>::iterator entryPtr = subscriptions.find(fromProc);
    PLB_ASSERT(entryPtr != subscriptions.end());
    CommunicatorEntry &entry = entryPtr->second;

    // Empty messages are neither sent nor received.
    if (!entry.data.empty()) {
        // 1. Make sure the package of messages has been received.
        global::mpi().wait(&entry.messageRequest, &entry.messageStatus);

        // 2. The message package is split into individual messages.
        int pos = 0;
        PLB_ASSERT(entry.messages.size() == entry.lengths.size());
        for (pluint iMessage = 0; iMessage < entry.messages.size(); ++iMessage) {
            int length = entry.lengths[iMessage];
            entry.messages[iMessage].resize(length);
            PLB_ASSERT(pos + length <= (int)entry.data.size());
            if (!entry.messages[iMessage].empty()) {
                std::copy(
                    entry.data.begin() + pos, entry.data.begin() + pos + length,
                    entry.messages[iMessage].begin());
            }
            pos += length;
        }
    }
}

#endif  // PLB_MPI_PARALLEL

}  // namespace plb
