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

#ifndef FD_WEIGHTS_HH
#define FD_WEIGHTS_HH

#include <algorithm>
#include <cstdlib>

#include "core/util.h"
#include "finiteDifference/fdWeights.h"

namespace plb {

/* ***************** class FdWeights ***************************************** */

template <typename T, int order, int maxWidth>
FdWeights<T, order, maxWidth>::FdWeights()
{
    int n = maxWidth * (maxWidth + 1) / 2;
    weights = (T **)calloc(n, sizeof(T *));
    PLB_ASSERT(weights != 0);
}

template <typename T, int order, int maxWidth>
FdWeights<T, order, maxWidth>::FdWeights(FdWeights<T, order, maxWidth> const &rhs)
{
    int n = maxWidth * (maxWidth + 1) / 2;
    weights = (T **)calloc(n, sizeof(T *));
    PLB_ASSERT(weights != 0);

    for (int width = 1; width <= maxWidth; width++) {
        for (int position = 0; position < width; position++) {
            int index = (width - 1) * width / 2 + position;
            if (rhs.weights[index] != 0) {
                weights[index] = (T *)malloc(width * sizeof(T));
                PLB_ASSERT(weights[index] != 0);
                for (int i = 0; i < width; i++) {
                    weights[index][i] = rhs.weights[index][i];
                }
            }
        }
    }
}

template <typename T, int order, int maxWidth>
void FdWeights<T, order, maxWidth>::swap(FdWeights<T, order, maxWidth> &rhs)
{
    std::swap(weights, rhs.weights);
}

template <typename T, int order, int maxWidth>
FdWeights<T, order, maxWidth> &FdWeights<T, order, maxWidth>::operator=(
    FdWeights<T, order, maxWidth> const &rhs)
{
    FdWeights<T, order, maxWidth> tmp(rhs);
    swap(tmp);
    return *this;
}

template <typename T, int order, int maxWidth>
FdWeights<T, order, maxWidth> *FdWeights<T, order, maxWidth>::clone() const
{
    return new FdWeights<T, order, maxWidth>(*this);
}

template <typename T, int order, int maxWidth>
FdWeights<T, order, maxWidth>::~FdWeights()
{
    int n = maxWidth * (maxWidth + 1) / 2;
    for (int i = 0; i < n; i++) {
        free(weights[i]);
    }
    free(weights);
}

template <typename T, int order, int maxWidth>
T const *FdWeights<T, order, maxWidth>::getWeights(int width, int position) const
{
    PLB_ASSERT(width > 0 && width <= maxWidth);
    PLB_ASSERT(position >= 0 && position < width);

    int index = (width - 1) * width / 2 + position;
    if (weights[index] == 0) {
        weights[index] = (T *)malloc(width * sizeof(T));
        PLB_ASSERT(weights[index] != 0);
        computeWeights(width, position, weights[index]);
    }

    return (weights[index]);
}

template <typename T, int order, int maxWidth>
void FdWeights<T, order, maxWidth>::getWeights(int width, T z, T *w) const
{
    PLB_ASSERT(width > 0);
    PLB_ASSERT(util::greaterEqual_abs<T>(z, (T)0) && util::lessEqual_abs<T>(z, (T)(width - 1)));
    PLB_ASSERT(w != 0);

    computeWeights(width, z, w);
}

/* This algorithm for computing the finite difference weights is based on
 * the article:
 *
 * Bengt Fornberg, "Calculation of weights in finite difference formulas",
 * SIAM Rev., Vol. 40, No. 3, pp. 685-691, September 1998.
 */
template <typename T, int order, int maxWidth>
void FdWeights<T, order, maxWidth>::computeWeights(int width, T z, T *w) const
{
    int n = width - 1;
    int m = order;

    T *c = (T *)calloc((n + 1) * (m + 1), sizeof(T));
    PLB_ASSERT(c != 0);

    int mp1 = m + 1;

#define PLB_IND(i, j) ((i)*mp1 + (j))

    T c1 = 1.0;
    T c4 = (T)0 - z;

    c[0] = 1.0;
    for (int i = 1; i <= n; i++) {
        int mn = std::min(i, m);
        T c2 = 1.0;
        T c5 = c4;
        c4 = (T)i - z;
        for (int j = 0; j <= i - 1; j++) {
            T c3 = i - j;
            c2 *= c3;
            if (j == i - 1) {
                for (int k = mn; k >= 1; k--) {
                    c[PLB_IND(i, k)] =
                        c1 * ((T)k * c[PLB_IND(i - 1, k - 1)] - c5 * c[PLB_IND(i - 1, k)]) / c2;
                }
                c[PLB_IND(i, 0)] = -c1 * c5 * c[PLB_IND(i - 1, 0)] / c2;
            }
            for (int k = mn; k >= 1; k--) {
                c[PLB_IND(j, k)] = (c4 * c[PLB_IND(j, k)] - (T)k * c[PLB_IND(j, k - 1)]) / c3;
            }
            c[PLB_IND(j, 0)] *= c4 / c3;
        }
        c1 = c2;
    }

    for (int i = 0; i < width; i++) {
        w[i] = c[PLB_IND(i, order)];
    }

#undef PLB_IND

    free(c);
}

}  // namespace plb

#endif  // FD_WEIGHTS_HH
