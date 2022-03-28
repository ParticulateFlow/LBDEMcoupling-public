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

#include "algorithm/statistics.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "core/util.h"

namespace plb {

namespace util {

Stats::Stats(std::vector<double> const &data)
{
    if (data.empty()) {
        min = max = mean = stddev = 0.;
        return;
    }

    min = std::numeric_limits<double>::max();
    max = -std::numeric_limits<double>::max();
    mean = 0.;
    stddev = 0.;
    for (pluint i = 0; i < data.size(); ++i) {
        mean += data[i];
        if (data[i] < min)
            min = data[i];
        if (data[i] > max)
            max = data[i];
    }
    mean /= (double)data.size();
    if (data.size() == 1) {
        stddev = 0.;
    } else {
        for (pluint i = 0; i < data.size(); ++i) {
            stddev += sqr(data[i] - mean);
        }
        stddev = std::sqrt(stddev / ((double)(data.size() - 1)));
    }
}

double Stats::getMean() const
{
    return mean;
}

double Stats::getStddev() const
{
    return stddev;
}

double Stats::getMin() const
{
    return min;
}

double Stats::getMax() const
{
    return max;
}

}  // namespace util

}  // namespace plb
