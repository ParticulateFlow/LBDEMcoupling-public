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

#ifndef LINEAR_ALGEBRA_HH
#define LINEAR_ALGEBRA_HH

#include <algorithm>
#include <cmath>

#include "algorithm/linearAlgebra.h"

namespace plb {

template <typename T>
void gramSchmidt(Array<T, 3> const &v1Unit, Array<T, 3> &v2Unit, Array<T, 3> &v3Unit)
{
    Array<T, 3> e1((T)1, T(), T()), e2(T(), (T)1, T()), e3(T(), T(), (T)1);
    T v11 = dot(v1Unit, e1), av11 = std::fabs(v11);
    T v12 = dot(v1Unit, e2), av12 = std::fabs(v12);
    T v13 = dot(v1Unit, e3), av13 = std::fabs(v13);

    if (av11 > av12) {
        if (av11 > av13) {  // e1 best aligned with v1Unit.
            v2Unit = e2 - v12 * v1Unit;
            v3Unit = e3 - v13 * v1Unit;
        } else {  // e3 best aligned with v1Unit.
            v2Unit = e1 - v11 * v1Unit;
            v3Unit = e2 - v12 * v1Unit;
        }
    } else {                // av12 >= av11
        if (av12 > av13) {  // e2 best aligned with v1Unit.
            v2Unit = e1 - v11 * v1Unit;
            v3Unit = e3 - v13 * v1Unit;
        } else {  // e3 best aligned with v1Unit.
            v2Unit = e1 - v11 * v1Unit;
            v3Unit = e2 - v12 * v1Unit;
        }
    }

    v2Unit /= norm(v2Unit);
    v3Unit -= dot(v2Unit, v3Unit) * v2Unit;
    v3Unit /= norm(v3Unit);
}

template <typename T>
T hypot2(T x, T y)
{
    return std::sqrt(x * x + y * y);
}

// Symmetric Householder reduction to tridiagonal form.
template <typename T>
void tred2(Array<Array<T, 3>, 3> &V, Array<T, 3> &d, Array<T, 3> &e)
{
    //  This is derived from the Algol procedures tred2 by
    //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutine in EISPACK.

    static const int n = 3;
    for (int j = 0; j < n; j++) {
        d[j] = V[n - 1][j];
    }

    // Householder reduction to tridiagonal form.

    for (int i = n - 1; i > 0; i--) {
        // Scale to avoid under/overflow.

        T scale = 0.0;
        T h = 0.0;
        for (int k = 0; k < i; k++) {
            scale = scale + std::fabs(d[k]);
        }
        if (scale == 0.0) {
            e[i] = d[i - 1];
            for (int j = 0; j < i; j++) {
                d[j] = V[i - 1][j];
                V[i][j] = 0.0;
                V[j][i] = 0.0;
            }
        } else {
            // Generate Householder vector.

            for (int k = 0; k < i; k++) {
                d[k] /= scale;
                h += d[k] * d[k];
            }
            T f = d[i - 1];
            T g = std::sqrt(h);
            if (f > 0) {
                g = -g;
            }
            e[i] = scale * g;
            h = h - f * g;
            d[i - 1] = f - g;
            for (int j = 0; j < i; j++) {
                e[j] = 0.0;
            }

            // Apply similarity transformation to remaining columns.

            for (int j = 0; j < i; j++) {
                f = d[j];
                V[j][i] = f;
                g = e[j] + V[j][j] * f;
                for (int k = j + 1; k <= i - 1; k++) {
                    g += V[k][j] * d[k];
                    e[k] += V[k][j] * f;
                }
                e[j] = g;
            }
            f = 0.0;
            for (int j = 0; j < i; j++) {
                e[j] /= h;
                f += e[j] * d[j];
            }
            T hh = f / (h + h);
            for (int j = 0; j < i; j++) {
                e[j] -= hh * d[j];
            }
            for (int j = 0; j < i; j++) {
                f = d[j];
                g = e[j];
                for (int k = j; k <= i - 1; k++) {
                    V[k][j] -= (f * e[k] + g * d[k]);
                }
                d[j] = V[i - 1][j];
                V[i][j] = 0.0;
            }
        }
        d[i] = h;
    }

    // Accumulate transformations.

    for (int i = 0; i < n - 1; i++) {
        V[n - 1][i] = V[i][i];
        V[i][i] = 1.0;
        T h = d[i + 1];
        if (h != 0.0) {
            for (int k = 0; k <= i; k++) {
                d[k] = V[k][i + 1] / h;
            }
            for (int j = 0; j <= i; j++) {
                T g = 0.0;
                for (int k = 0; k <= i; k++) {
                    g += V[k][i + 1] * V[k][j];
                }
                for (int k = 0; k <= i; k++) {
                    V[k][j] -= g * d[k];
                }
            }
        }
        for (int k = 0; k <= i; k++) {
            V[k][i + 1] = 0.0;
        }
    }
    for (int j = 0; j < n; j++) {
        d[j] = V[n - 1][j];
        V[n - 1][j] = 0.0;
    }
    V[n - 1][n - 1] = 1.0;
    e[0] = 0.0;
}

// Symmetric tridiagonal QL algorithm.

template <typename T>
void tql2(Array<Array<T, 3>, 3> &V, Array<T, 3> &d, Array<T, 3> &e)
{
    //  This is derived from the Algol procedures tql2, by
    //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutine in EISPACK.

    static const int n = 3;
    for (int i = 1; i < n; i++) {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;

    T f = 0.0;
    T tst1 = 0.0;
    T eps = std::pow((T)2.0, (T)-52.0);
    for (int l = 0; l < n; l++) {
        // Find small subdiagonal element

        tst1 = std::max(tst1, std::fabs(d[l]) + std::fabs(e[l]));
        int m = l;
        while (m < n) {
            if (std::fabs(e[m]) <= eps * tst1) {
                break;
            }
            m++;
        }
        PLB_ASSERT(m < n);  // e[n-1] is always zero so there is never any exit through the bottom
                            // of the while-loop above.

        // If m == l, d[l] is an eigenvalue,
        // otherwise, iterate.

        if (m > l) {
            int iter = 0;
            do {
                iter = iter + 1;  // (Could check iteration count here.)

                // Compute implicit shift

                T g = d[l];
                T p = (d[l + 1] - g) / (2.0 * e[l]);
                T r = hypot2(p, (T)1.0);
                if (p < 0) {
                    r = -r;
                }
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                T dl1 = d[l + 1];
                T h = g - d[l];
                for (int i = l + 2; i < n; i++) {
                    d[i] -= h;
                }
                f = f + h;

                // Implicit QL transformation.

                p = d[m];
                T c = 1.0;
                T c2 = c;
                T c3 = c;
                T el1 = e[l + 1];
                T s = 0.0;
                T s2 = 0.0;
                for (int i = m - 1; i >= l; i--) {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = hypot2(p, e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);

                    // Accumulate transformation.

                    for (int k = 0; k < n; k++) {
                        h = V[k][i + 1];
                        V[k][i + 1] = s * V[k][i] + c * h;
                        V[k][i] = c * V[k][i] - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;

                // Check for convergence.

            } while (std::fabs(e[l]) > eps * tst1);
        }
        d[l] = d[l] + f;
        e[l] = 0.0;
    }

    // Sort eigenvalues and corresponding vectors.

    for (int i = 0; i < n - 1; i++) {
        int k = i;
        T p = d[i];
        for (int j = i + 1; j < n; j++) {
            if (d[j] < p) {
                k = j;
                p = d[j];
            }
        }
        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (int j = 0; j < n; j++) {
                p = V[j][i];
                V[j][i] = V[j][k];
                V[j][k] = p;
            }
        }
    }
}

template <typename T>
void eigenDecomposition(Array<Array<T, 3>, 3> const &A, Array<Array<T, 3>, 3> &V, Array<T, 3> &d)
{
    static const int n = 3;
    Array<T, n> e;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            V[i][j] = A[i][j];
        }
    }
    tred2(V, d, e);
    tql2(V, d, e);
}

}  // namespace plb

#endif  // LINEAR_ALGEBRA_HH
