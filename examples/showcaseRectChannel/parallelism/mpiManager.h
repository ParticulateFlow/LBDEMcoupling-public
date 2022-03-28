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

#ifndef MPI_MANAGER_H
#define MPI_MANAGER_H

#include "core/globalDefs.h"

#ifdef PLB_MPI_PARALLEL
// DISABLE_WARNING_PUSH
// DISABLE_WARNING_CAST_FUNCTION_TYPE
#include <mpi.h>
// DISABLE_WARNING_POP
#include <string>
#include <vector>
#endif

namespace plb {

namespace global {

#ifdef PLB_MPI_PARALLEL

/// Wrapper functions that simplify the use of MPI
class MpiManager {
public:
    /// Initializes the MPI manager and the MPI machine.
    void init(int *argc, char ***argv, bool verbous = false);
    /// Initializes the MPI manager, but assumes that the MPI
    ///   machine is handled by another instance.
    void init(MPI_Comm globalCommunicator_);
    /// Initializes the MPI manager, but assumes that the MPI
    ///   machine is handled by another instance.
    void init();
    /// Returns the number of processes
    int getSize() const;
    /// Returns the process ID
    int getRank() const;
    /// Returns process ID of main processor
    int bossId() const;
    /// Tells whether current processor is main processor
    bool isMainProcessor() const;
    /// Returns universal MPI-time in seconds
    double getTime() const;
    /// Returns the global communicator for this program or library instance.
    MPI_Comm getGlobalCommunicator() const;

    /// Synchronizes the processes
    void barrier();

    /// Sends data at *buf, blocking
    template <typename T>
    void send(T *buf, int count, int dest, int tag = 0);

    /// Sends data at *buf, non blocking
    template <typename T>
    void iSend(T *buf, int count, int dest, MPI_Request *request, int tag = 0);

    /// Sends data at *buf, assuming that receiver is ready.
    template <typename T>
    void rSend(T *buf, int count, int dest, int tag = 0);

    /// Sends data at *buf, non blocking and request free
    template <typename T>
    void iSendRequestFree(T *buf, int count, int dest, int tag = 0);

    /// Receives data at *buf, blocking
    template <typename T>
    void receive(T *buf, int count, int source, int tag = 0);

    /// Receives data at *buf, non blocking
    template <typename T>
    void iRecv(T *buf, int count, int source, MPI_Request *request, int tag = 0);

    /// Send and receive data between two partners
    template <typename T>
    void sendRecv(T *sendBuf, T *recvBuf, int count, int dest, int source, int tag = 0);

    /// Sends data to master processor
    template <typename T>
    void sendToMaster(T *sendBuf, int sendCount, bool iAmRoot);

    /// Special case for sending strings. Memory handling is automatic.
    void sendToMaster(std::string &message, bool iAmRoot);

    /// Scatter data from one processor over multiple processors
    template <typename T>
    void scatterV(T *sendBuf, T *recvBuf, int *sendCounts, int root = 0);

    /// Gather data from multiple processors to one processor
    template <typename T>
    void gatherV(T *sendBuf, T *recvBuf, int *recvCounts, int root = 0);

    /// Broadcast data from one processor to multiple processors
    template <typename T>
    void bCast(T *sendBuf, int sendCount, int root = 0);

    /// Special case for broadcasting strings. Memory handling is automatic.
    void bCast(std::string &message, int root = 0);

    /// Broadcast data when root is unknown to other processors
    template <typename T>
    void bCastThroughMaster(T *sendBuf, int sendCount, bool iAmRoot);

    /// Reduction operation toward one processor
    template <typename T>
    void reduce(T sendVal, T &recvVal, MPI_Op op, int root = 0);

    /// Element-per-element reduction of a vector of data
    template <typename T>
    void reduceVect(std::vector<T> &sendVal, std::vector<T> &recvVal, MPI_Op op, int root = 0);

    /// Inplace element-per-element reduction of a vector of data; result
    ///   available on all MPI threads.
    template <typename T>
    void allReduceVect(std::vector<T> &sendRecvVal, MPI_Op op);

    /// Reduction operation, followed by a broadcast
    template <typename T>
    void reduceAndBcast(T &reductVal, MPI_Op op, int root = 0);

    /// Scan operation (partial reductions) of data on a collection of processes
    template <typename T>
    void scan(T sendVal, T &recvVal, MPI_Op op);

    /// Complete a non-blocking MPI operation
    void wait(MPI_Request *request, MPI_Status *status);

private:
    /// Implementation code for Scatter
    template <typename T>
    void scatterv_impl(
        T *sendBuf, int *sendCounts, int *displs, T *recvBuf, int recvCount, int root);

    /// Implementation code for Gather
    template <typename T>
    void gatherv_impl(
        T *sendBuf, int sendCount, T *recvBuf, int *recvCounts, int *displs, int root);

private:
    MpiManager();
    ~MpiManager();

private:
    int numTasks, taskId;
    bool ok;
    bool responsibleForMpiMachine;
    MPI_Comm globalCommunicator;

    friend MpiManager &mpi();
};

#else  // #ifdef PLB_MPI_PARALLEL

class MpiManager {
public:
    /// Initializes the mpi manager
    void init(int *argc, char ***argv, bool verbous = false) { }
    /// Initializes the MPI manager, but assumes that the MPI
    ///   machine is handled by another instance.
    void init() { }
    /// Returns the number of processes
    int getSize() const
    {
        return 1;
    }
    /// Returns the process ID
    int getRank() const
    {
        return 0;
    }
    /// Returns process ID of main processor
    int bossId() const
    {
        return 0;
    }
    /// Tells whether current processor is main processor
    bool isMainProcessor() const
    {
        return true;
    }
    /// Broadcast data from one processor to multiple processors
    template <typename T>
    void bCast(T *sendBuf, int sendCount, int root = 0)
    { }
    /// Special case for broadcasting strings. Memory handling is automatic.
    void bCast(std::string &message, int root = 0) { }
    /// Special case for sending strings. Memory handling is automatic.
    void sendToMaster(std::string &message, bool iAmRoot) { }
    /// Synchronizes the processes
    void barrier() { }

    friend MpiManager &mpi();
};

#endif  // PLB_MPI_PARALLEL

inline MpiManager &mpi()
{
    static MpiManager instance;
    return instance;
}

}  // namespace global

}  // namespace plb

#endif  // MPI_MANAGER_H
