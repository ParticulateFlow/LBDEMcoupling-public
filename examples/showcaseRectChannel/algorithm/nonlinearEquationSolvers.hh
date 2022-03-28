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

#ifndef NONLINEAR_EQUATION_SOLVERS_HH
#define NONLINEAR_EQUATION_SOLVERS_HH

#include <cmath>

#include "algorithm/nonlinearEquationSolvers.h"
#include "core/globalDefs.h"
#include "core/plbDebug.h"

namespace plb {

/* ***************** class NewtonRaphsonMethod ***************************************** */

template <typename T>
NewtonRaphsonMethod<T>::NewtonRaphsonMethod(T tolerance_, plint maxNumOfIterations_) :
    tolerance(tolerance_), maxNumOfIterations(maxNumOfIterations_)
{
    PLB_ASSERT(tolerance > 0.0);
    PLB_ASSERT(maxNumOfIterations > 0);
}

template <typename T>
NewtonRaphsonMethod<T> *NewtonRaphsonMethod<T>::clone() const
{
    return new NewtonRaphsonMethod<T>(*this);
}

template <typename T>
bool NewtonRaphsonMethod<T>::solve(
    Function f, Function df, void *data, T initialGuess, T &solution, T &absoluteError) const
{
    bool converged = false;
    T xOld = initialGuess;

    T xNew = xOld;
    T error = -1.0;
    for (plint i = 0; i < maxNumOfIterations; i++) {
        xNew = xOld - f(xOld, data) / df(xOld, data);
        error = std::fabs(xNew - xOld);
        if (error <= tolerance) {
            converged = true;
            break;
        }
        xOld = xNew;
    }

    solution = xNew;
    absoluteError = error;
    return (converged);
}

/* ***************** Bisection Method ***************************************** */

template <typename T, class Function>
bool bisect(Function const &function, T x0, T x1, T xacc, plint maxIter, T &result)
{
    PLB_ASSERT(maxIter > 0);
    PLB_ASSERT(xacc > T());

    T dx, f, fmid, xmid;

    f = function(x0);
    fmid = function(x1);
    if (f * fmid >= T()) {
        result = T();
        return false;
    }

    // Orient the search so that f>0 lies at x+dx.
    if (f < T()) {
        dx = x1 - x0;
        result = x0;
    } else {
        dx = x0 - x1;
        result = x1;
    }

    for (int j = 0; j < maxIter; ++j) {  // Bisection loop.
        dx *= 0.5;
        xmid = result + dx;
        fmid = function(xmid);

        if (fmid <= T()) {
            result = xmid;
        }

        if (std::fabs(dx) < xacc || fmid == T()) {
            return true;
        }
    }

    result = T();
    return false;
}

template <typename T, class Function>
bool brentSolve(Function const &func, T x0, T x1, T xacc, plint maxIter, T &result)
{
    result = x0;
    plint iter;
    T a = x0, b = x1, c = x1, d = 0., e = 0., min1, min2;
    T fa = func(a), fb = func(b), fc, p, q, r, s, tol1, xm;
    static const T eps = getEpsilon<T>();
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
        return false;
    fc = fb;
    for (iter = 0; iter < maxIter; iter++) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c = a;
            fc = fa;
            e = d = b - a;
        }
        if (std::fabs(fc) < std::fabs(fb)) {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        tol1 = 2.0 * eps * std::fabs(b) + 0.5 * xacc;
        xm = 0.5 * (c - b);
        if (std::fabs(xm) <= tol1 || fb == 0.0) {
            result = b;
            return true;
        }
        if (std::fabs(e) >= tol1 && std::fabs(fa) > std::fabs(fb)) {
            s = fb / fa;
            if (a == c) {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0)
                q = -q;
            p = std::fabs(p);
            min1 = 3.0 * xm * q - std::fabs(tol1 * q);
            min2 = std::fabs(e * q);
            if (2.0 * p < (min1 < min2 ? min1 : min2)) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }
        a = b;
        fa = fb;
        if (std::fabs(d) > tol1)
            b += d;
        else
            b += xm >= 0.0 ? std::fabs(tol1) : -std::fabs(tol1);
        fb = func(b);
    }
    return false;
}

}  // namespace plb

#endif  // NONLINEAR_EQUATION_SOLVERS_HH
