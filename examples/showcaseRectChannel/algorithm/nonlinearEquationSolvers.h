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

#ifndef NONLINEAR_EQUATION_SOLVERS_H
#define NONLINEAR_EQUATION_SOLVERS_H

#include "core/globalDefs.h"

namespace plb {

template <typename T>
class NewtonRaphsonMethod {
public:
    typedef T (*Function)(T x, void *data);

public:
    NewtonRaphsonMethod(T tolerance_, plint maxNumOfIterations_);
    NewtonRaphsonMethod<T> *clone() const;
    // "f" is the real function that defines the equation "f = 0" to be solved.
    // "df" is the derivative of f with respect to the uknown variable.
    // "data" is a pointer to the data of the f and df functions.
    // "initialGuess" is obviously the initial guess to start the iterations.
    // "solution" is a reference to the solution. When this function returns, the value of
    // "solution" will be the last approximation of the solution by the Newton-Raphson algorithm.
    // This function returns "true" if the absolute error level "tolerance" was achieved
    // before "maxNumOfIterations" iterations were performed. The actual absolute error
    // for the "solution" value which is "returned", is "absoluteError".
    bool solve(
        Function f, Function df, void *data, T initialGuess, T &solution, T &absoluteError) const;
    T getTolerance() const
    {
        return tolerance;
    }
    T getMaxNumOfIterations() const
    {
        return maxNumOfIterations;
    }

private:
    T tolerance;
    plint maxNumOfIterations;
};

template <typename T, class Function>
bool bisect(Function const &function, T x0, T x1, T xacc, plint maxIter, T &result);

template <typename T, class Function>
bool brentSolve(Function const &func, T x0, T x1, T xacc, plint maxIter, T &result);

}  // namespace plb

#endif  // NONLINEAR_EQUATION_SOLVERS_H
