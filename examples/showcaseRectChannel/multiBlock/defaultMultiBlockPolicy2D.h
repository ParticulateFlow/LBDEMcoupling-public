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
 * Default policy for instantiating serial/parallel multi-blocks
 */
#ifndef DEFAULT_MULTI_BLOCK_POLICY_2D_H
#define DEFAULT_MULTI_BLOCK_POLICY_2D_H

#include <cmath>

#include "core/globalDefs.h"
#include "multiBlock/combinedStatistics.h"
#include "multiBlock/multiBlockManagement2D.h"
#include "multiBlock/serialBlockCommunicator2D.h"
#include "multiBlock/serialMultiBlockLattice2D.h"
#include "multiBlock/serialMultiDataField2D.h"
#include "multiBlock/serialMultiDataField2D.hh"
#include "multiBlock/staticRepartitions2D.h"
#include "multiBlock/threadAttribution.h"
#include "parallelism/mpiManager.h"
#include "parallelism/parallelBlockCommunicator2D.h"
#include "parallelism/parallelMultiBlockLattice2D.h"
#include "parallelism/parallelMultiDataField2D.h"
#include "parallelism/parallelStatistics.h"

namespace plb {

class DefaultMultiBlockPolicy2D {
public:
    void toggleBlockingCommunication(bool useBlockingCommunication_)
    {
        useBlockingCommunication = useBlockingCommunication_;
    }

    BlockCommunicator2D *getBlockCommunicator()
    {
#ifdef PLB_MPI_PARALLEL
        if (useBlockingCommunication) {
            return new ParallelBlockCommunicator2D();
        } else {
            return new ParallelBlockCommunicator2D();
        }
#else
        return new SerialBlockCommunicator2D();
#endif
    }

    CombinedStatistics *getCombinedStatistics()
    {
#ifdef PLB_MPI_PARALLEL
        return new ParallelCombinedStatistics();
#else
        return new SerialCombinedStatistics();
#endif
    }

    template <typename T, template <typename U> class Descriptor>
    MultiCellAccess2D<T, Descriptor> *getMultiCellAccess()
    {
#ifdef PLB_MPI_PARALLEL
        return new ParallelCellAccess2D<T, Descriptor>();
#else
        return new SerialCellAccess2D<T, Descriptor>();
#endif
    }

    template <typename T>
    MultiScalarAccess2D<T> *getMultiScalarAccess()
    {
#ifdef PLB_MPI_PARALLEL
        return new ParallelScalarAccess2D<T>();
#else
        return new SerialScalarAccess2D<T>();
#endif
    }

    template <typename T, int nDim>
    MultiTensorAccess2D<T, nDim> *getMultiTensorAccess()
    {
#ifdef PLB_MPI_PARALLEL
        return new ParallelTensorAccess2D<T, nDim>();
#else
        return new SerialTensorAccess2D<T, nDim>();
#endif
    }

    template <typename T>
    MultiNTensorAccess2D<T> *getMultiNTensorAccess()
    {
#ifdef PLB_MPI_PARALLEL
        return new ParallelNTensorAccess2D<T>();
#else
        return new SerialNTensorAccess2D<T>();
#endif
    }

    ThreadAttribution *getThreadAttribution()
    {
#ifdef PLB_MPI_PARALLEL
        return new OneToOneThreadAttribution();
#else
        return new SerialThreadAttribution();
#endif
    }

    MultiBlockManagement2D getMultiBlockManagement(Box2D const &domain, plint envelopeWidth)
    {
        return MultiBlockManagement2D(
            createRegularDistribution2D(domain, numProcesses), getThreadAttribution(),
            envelopeWidth);
    }

    MultiBlockManagement2D getMultiBlockManagement(plint nx, plint ny, plint envelopeWidth)
    {
        return MultiBlockManagement2D(
            createRegularDistribution2D(nx, ny, numProcesses), getThreadAttribution(),
            envelopeWidth);
    }

    void setNumGridPoints(plint numGridPoints_)
    {
        numGridPoints = numGridPoints_;
        numGridPointsSpecified = true;
    }

    plint getNumGridPoints() const
    {
        return numGridPoints;
    }

    void setNumProcesses(int numProcesses_)
    {
        numProcesses = numProcesses_;
        if (!numGridPointsSpecified) {
            numGridPoints = numProcesses;
        }
    }

    int getNumProcesses() const
    {
        return numProcesses;
    }

private:
    DefaultMultiBlockPolicy2D() :
        numProcesses(global::mpi().getSize()),
        numGridPointsSpecified(false),
        useBlockingCommunication(false)
    {
        numGridPoints = numProcesses;
    }
    friend DefaultMultiBlockPolicy2D &defaultMultiBlockPolicy2D();

private:
    int numProcesses;
    plint numGridPoints;
    bool numGridPointsSpecified;
    bool useBlockingCommunication;
};

inline DefaultMultiBlockPolicy2D &defaultMultiBlockPolicy2D()
{
    static DefaultMultiBlockPolicy2D singleton;
    return singleton;
}

}  // namespace plb

#endif  // DEFAULT_MULTI_BLOCK_POLICY_2D_H
