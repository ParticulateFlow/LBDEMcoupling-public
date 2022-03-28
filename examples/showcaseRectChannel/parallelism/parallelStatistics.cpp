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
 * The CombinedStatistics class -- generic implementation.
 */
#include "parallelism/parallelStatistics.h"

#include <cmath>

#include "parallelism/mpiManager.h"

namespace plb {

#ifdef PLB_MPI_PARALLEL

ParallelCombinedStatistics *ParallelCombinedStatistics::clone() const
{
    return new ParallelCombinedStatistics(*this);
}

void ParallelCombinedStatistics::reduceStatistics(
    std::vector<double> &averageObservables, std::vector<double> &sumWeights,
    std::vector<double> &sumObservables, std::vector<double> &maxObservables,
    std::vector<plint> &intSumObservables) const
{
    // Averages
    for (pluint iAverage = 0; iAverage < averageObservables.size(); ++iAverage) {
        double globalAverage, globalWeight;
        global::mpi().reduce(
            averageObservables[iAverage] * sumWeights[iAverage], globalAverage, MPI_SUM);
        global::mpi().reduce(sumWeights[iAverage], globalWeight, MPI_SUM);
        if (global::mpi().isMainProcessor() && std::fabs(globalWeight) > 0.5) {
            globalAverage /= globalWeight;
        }
        global::mpi().bCast(&globalAverage, 1);
        averageObservables[iAverage] = globalAverage;
    }

    // Sum
    for (pluint iSum = 0; iSum < sumObservables.size(); ++iSum) {
        double globalSum;
        global::mpi().reduce(sumObservables[iSum], globalSum, MPI_SUM);
        global::mpi().bCast(&globalSum, 1);
        sumObservables[iSum] = globalSum;
    }

    // Max
    for (pluint iMax = 0; iMax < maxObservables.size(); ++iMax) {
        double globalMax;
        global::mpi().reduce(maxObservables[iMax], globalMax, MPI_MAX);
        global::mpi().bCast(&globalMax, 1);
        maxObservables[iMax] = globalMax;
    }

    // Integer sum
    for (pluint iSum = 0; iSum < intSumObservables.size(); ++iSum) {
        plint globalSum;
        global::mpi().reduce(intSumObservables[iSum], globalSum, MPI_SUM);
        global::mpi().bCast(&globalSum, 1);
        intSumObservables[iSum] = globalSum;
    }
}

#endif  // PLB_MPI_PARALLEL

}  // namespace plb
