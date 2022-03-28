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

#ifndef QUADRATURE_HH
#define QUADRATURE_HH

#include <cmath>
#include <limits>
#include <vector>

#include "algorithm/quadrature.h"
#include "core/globalDefs.h"
#include "core/plbDebug.h"

namespace plb {

/* ***************** class GaussLegendreQuadrature ***************************************** */

template <typename T>
GaussLegendreQuadrature<T>::GaussLegendreQuadrature(plint n_, plint maxNumOfIterations_) :
    n(n_), maxNumOfIterations(maxNumOfIterations_)
{
    PLB_ASSERT(n >= 0);
    PLB_ASSERT(maxNumOfIterations > 0);
    nodes.resize(n);
    weights.resize(n);
    evaluateNodesAndWeights();
}

template <typename T>
GaussLegendreQuadrature<T> *GaussLegendreQuadrature<T>::clone() const
{
    return new GaussLegendreQuadrature<T>(*this);
}

template <typename T>
T GaussLegendreQuadrature<T>::evaluateIntegral(
    IntegralKernel integralKernel, void *data, T x0, T x1) const
{
    T a = 0.5 * (x1 - x0);
    T b = 0.5 * (x0 + x1);

    T integral = 0.0;
    for (plint i = 0; i < n; i++) {
        T x_i = a * nodes[i] + b;
        integral += weights[i] * integralKernel(x_i, data);
    }
    integral *= a;

    return integral;
}

template <typename T>
void GaussLegendreQuadrature<T>::evaluateLegendrePolynomialAndDerivative(
    plint k, T x, T &L_k, T &dL_k) const
{
    if (k == 0) {
        L_k = 1.0;
        dL_k = 0.0;
    } else if (k == 1) {
        L_k = x;
        dL_k = 1.0;
    } else {
        T L_m1 = 1.0;
        T L_m = x;
        T dL_m1 = 0.0;
        T dL_m = 1.0;
        T L = 0.0;
        T dL = 0.0;
        for (plint m = 1; m < k; m++) {
            L = (2.0 * (T)m + 1.0) * x * L_m - (T)m * L_m1;
            L /= (T)m + 1.0;
            dL = (2.0 * (T)m + 1.0) * (L_m + x * dL_m) - (T)m * dL_m1;
            dL /= (T)m + 1.0;
            L_m1 = L_m;
            L_m = L;
            dL_m1 = dL_m;
            dL_m = dL;
        }
        L_k = L;
        dL_k = dL;
    }
}

template <typename T>
void GaussLegendreQuadrature<T>::evaluateNodesAndWeights()
{
    static T tol = getEpsilon<T>((T)4);
    static T pi = std::acos((T)-1);

    if (n == 0) {
        nodes[0] = 0.0;
        weights[0] = 2.0;
    } else if (n == 1) {
        nodes[0] = -std::sqrt((T)1.0 / (T)3.0);
        weights[0] = 1.0;
        nodes[1] = -nodes[0];
        weights[1] = weights[0];
    } else {
        for (plint j = 0; j < n / 2; j++) {
            T x_j = -std::cos(((2.0 * (T)j + 1.0) / (2.0 * (T)n)) * pi);
            if (j > 0) {
                x_j = 0.5 * (x_j + nodes[j - 1]);
            }
            T L = 0.0;
            T dL = 0.0;
            for (plint k = 0; k < maxNumOfIterations; k++) {
                T s = 0.0;
                for (plint l = 0; l < j; l++) {
                    s += 1.0 / (x_j - nodes[l]);
                }
                evaluateLegendrePolynomialAndDerivative(n, x_j, L, dL);
                T delta = -L / (dL - L * s);
                x_j += delta;
                if (std::fabs(delta) <= tol * std::fabs(x_j)) {
                    break;
                }
            }
            evaluateLegendrePolynomialAndDerivative(n, x_j, L, dL);
            nodes[j] = x_j;
            nodes[n - 1 - j] = -x_j;
            weights[j] = 2.0 / ((1.0 - x_j * x_j) * dL * dL);
            weights[n - 1 - j] = weights[j];
        }
    }

    if ((n + 1) % 2 == 0) {
        T L = 0.0;
        T dL = 0.0;
        evaluateLegendrePolynomialAndDerivative(n, 0.0, L, dL);
        nodes[n / 2] = 0.0;
        weights[n / 2] = 2.0 / (dL * dL);
    }
}

}  // namespace plb

#endif  // QUADRATURE_HH
