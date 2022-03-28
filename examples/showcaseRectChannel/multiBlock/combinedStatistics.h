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
 * The CombinedStatistics class -- header file.
 */
#ifndef COMBINED_STATISTICS_H
#define COMBINED_STATISTICS_H

#include <vector>

#include "core/blockStatistics.h"
#include "core/globalDefs.h"

namespace plb {

class CombinedStatistics {
public:
    virtual ~CombinedStatistics();
    virtual CombinedStatistics *clone() const = 0;
    void combine(
        std::vector<BlockStatistics const *> &individualStatistics, BlockStatistics &result) const;

protected:
    virtual void reduceStatistics(
        std::vector<double> &averageObservables, std::vector<double> &sumWeights,
        std::vector<double> &sumObservables, std::vector<double> &maxObservables,
        std::vector<plint> &intSumObservables) const = 0;

private:
    void computeLocalAverage(
        std::vector<BlockStatistics const *> const &individualStatistics,
        std::vector<double> &averageObservables, std::vector<double> &sumWeights) const;
    void computeLocalSum(
        std::vector<BlockStatistics const *> const &individualStatistics,
        std::vector<double> &sumObservables) const;
    void computeLocalMax(
        std::vector<BlockStatistics const *> const &individualStatistics,
        std::vector<double> &maxObservables) const;
    void computeLocalIntSum(
        std::vector<BlockStatistics const *> const &individualStatistics,
        std::vector<plint> &intSumObservables) const;
};

class SerialCombinedStatistics : public CombinedStatistics {
public:
    virtual SerialCombinedStatistics *clone() const;

protected:
    virtual void reduceStatistics(
        std::vector<double> &averageObservables, std::vector<double> &sumWeights,
        std::vector<double> &sumObservables, std::vector<double> &maxObservables,
        std::vector<plint> &intSumObservables) const;
};

}  // namespace plb

#endif
