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
 * The CombinedStatistics class -- implementation.
 */
#include "multiBlock/combinedStatistics.h"

#include <cmath>
#include <limits>
#include <numeric>

namespace plb {

CombinedStatistics::~CombinedStatistics() { }

void CombinedStatistics::computeLocalAverage(
    std::vector<BlockStatistics const *> const &individualStatistics,
    std::vector<double> &averageObservables, std::vector<double> &sumWeights) const
{
    // For each "average observable"
    for (pluint iAverage = 0; iAverage < averageObservables.size(); ++iAverage) {
        averageObservables[iAverage] = 0.;
        sumWeights[iAverage] = 0.;
        // Compute local average
        for (pluint iStat = 0; iStat < individualStatistics.size(); ++iStat) {
            double newElement = individualStatistics[iStat]->getAverage(iAverage);
            double newWeight = individualStatistics[iStat]->getNumCells();
            averageObservables[iAverage] += newWeight * newElement;
            sumWeights[iAverage] += newWeight;
        }
        // Avoid division by zero
        if (std::fabs(sumWeights[iAverage]) > 0.5) {
            averageObservables[iAverage] /= sumWeights[iAverage];
        }
    }
}

void CombinedStatistics::computeLocalSum(
    std::vector<BlockStatistics const *> const &individualStatistics,
    std::vector<double> &sumObservables) const
{
    // For each "sum observable"
    for (pluint iSum = 0; iSum < sumObservables.size(); ++iSum) {
        sumObservables[iSum] = 0.;
        // Compute local sum
        for (pluint iStat = 0; iStat < individualStatistics.size(); ++iStat) {
            sumObservables[iSum] += individualStatistics[iStat]->getSum(iSum);
        }
    }
}

void CombinedStatistics::computeLocalMax(
    std::vector<BlockStatistics const *> const &individualStatistics,
    std::vector<double> &maxObservables) const
{
    // For each "max observable"
    for (pluint iMax = 0; iMax < maxObservables.size(); ++iMax) {
        // Use -max() instead of min(), because min<float> yields a positive value close to zero.
        maxObservables[iMax] = -std::numeric_limits<double>::max();
        // Compute local max
        for (pluint iStat = 0; iStat < individualStatistics.size(); ++iStat) {
            double newMax = individualStatistics[iStat]->getMax(iMax);
            if (newMax > maxObservables[iMax]) {
                maxObservables[iMax] = newMax;
            }
        }
    }
}

void CombinedStatistics::computeLocalIntSum(
    std::vector<BlockStatistics const *> const &individualStatistics,
    std::vector<plint> &intSumObservables) const
{
    // For each integer "sum observable"
    for (pluint iSum = 0; iSum < intSumObservables.size(); ++iSum) {
        intSumObservables[iSum] = 0;
        // Compute local sum
        for (pluint iStat = 0; iStat < individualStatistics.size(); ++iStat) {
            intSumObservables[iSum] += individualStatistics[iStat]->getIntSum(iSum);
        }
    }
}

void CombinedStatistics::combine(
    std::vector<BlockStatistics const *> &individualStatistics, BlockStatistics &result) const
{
    // Local averages
    std::vector<double> averageObservables(result.getAverageVect().size());
    std::vector<double> sumWeights(result.getAverageVect().size());
    computeLocalAverage(individualStatistics, averageObservables, sumWeights);

    // Local sums
    std::vector<double> sumObservables(result.getSumVect().size());
    computeLocalSum(individualStatistics, sumObservables);

    // Local maxima
    std::vector<double> maxObservables(result.getMaxVect().size());
    computeLocalMax(individualStatistics, maxObservables);

    // Local integer sums
    std::vector<plint> intSumObservables(result.getIntSumVect().size());
    computeLocalIntSum(individualStatistics, intSumObservables);

    // Compute global, cross-core statistics
    this->reduceStatistics(
        averageObservables, sumWeights, sumObservables, maxObservables, intSumObservables);

    // Update public statistics in resulting block
    result.evaluate(averageObservables, sumObservables, maxObservables, intSumObservables, 0);
}

SerialCombinedStatistics *SerialCombinedStatistics::clone() const
{
    return new SerialCombinedStatistics(*this);
}

void SerialCombinedStatistics::reduceStatistics(
    std::vector<double> &averageObservables, std::vector<double> &sumWeights,
    std::vector<double> &sumObservables, std::vector<double> &maxObservables,
    std::vector<plint> &intSumObservables) const
{
    // Do nothing in serial case
}

}  // namespace plb
