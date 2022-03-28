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

#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <vector>

#include "core/globalDefs.h"

namespace plb {

template <typename T>
class GaussLegendreQuadrature {
public:
    typedef T (*IntegralKernel)(T x, void *data);

public:
    GaussLegendreQuadrature(plint n_, plint maxNumOfIterations_ = 64);
    GaussLegendreQuadrature<T> *clone() const;
    T evaluateIntegral(IntegralKernel integralKernel, void *data, T x0, T x1) const;
    plint getN() const
    {
        return n;
    }
    std::vector<T> const &getNodes() const
    {
        return nodes;
    };
    std::vector<T> const &getWeights() const
    {
        return weights;
    };

private:
    void evaluateLegendrePolynomialAndDerivative(plint k, T x, T &L_k, T &dL_k) const;
    void evaluateNodesAndWeights();

private:
    plint n;
    std::vector<T> nodes;
    std::vector<T> weights;
    plint maxNumOfIterations;
};

}  // namespace plb

#endif  // QUADRATURE_H
