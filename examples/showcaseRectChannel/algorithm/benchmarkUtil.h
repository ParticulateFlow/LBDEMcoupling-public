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

#ifndef BENCHMARK_UTIL_H
#define BENCHMARK_UTIL_H

#include <deque>

#include "core/globalDefs.h"

namespace plb {

namespace util {

/// Check time-convergence of a scalar.
/** This class is useful, for example to check convergence of
 * the velocity field for the simulation of a stationary flow.
 * Convergence is claimed when the standard deviation of the
 * monitored value is smaller than epsilon times the average.
 * The statistics are taken over a macroscopic time scale of the
 * system.
 */
template <typename T>
class ValueTracer {
public:
    /// The only constructor.
    /** \param u The characteristic velocity of the system, for
     *          computation of the characteristic time scale.
     * \param L The characteristic length of the system, for
     *          computation of the characteristic time scale.
     * \param _epsilon Precision of the convergence.
     */
    ValueTracer(T u, T L, T epsilon);
    /// Change values of u and L to update characteristic scales of the system.
    void resetScale(T u, T L);
    /// reinitializes the values
    void resetValues();
    /// Get characteristic time scale.
    plint getDeltaT() const;
    /// Feed the object with a new measured scalar.
    void takeValue(T val, bool doPrint = false);
    /// Test for convergence, with respect to stdDev.
    bool hasConverged() const;
    /// Test for convergence, with respect to difference between min and max value;
    bool hasConvergedMinMax() const;
    T computeAverage() const;
    T computeStdDev(T average) const;
    void setEpsilon(T epsilon_);

private:
    plint deltaT;
    T epsilon;
    plint t;
    bool converged;
    std::deque<T> values;
};

/// Propose successive test values of a scalar (e.g. Re) to check stability of a system.
/** At first, the stability limit is explored by constant
 * increments/decrements of the scalar, and then, by successive
 * bisection.
 */
template <typename T>
class BisectStepper {
public:
    /// The only constructor.
    /** \param _iniVal Initial guess for the stability limit.
     * \param _step   Step size at which the value is initially
     *                incremented/decremented.
     */
    BisectStepper(T _iniVal, T _step = 0.);
    /// Get new value, and indicate if the previous value yielded a stable system or not.
    T getVal(bool stable, bool doPrint = false);
    /// Test for convergence.
    bool hasConverged(T epsilon) const;

private:
    T iniVal, currentVal, lowerVal, upperVal;
    T step;
    enum { first, up, down, bisect } state;
};

}  // namespace util

}  // namespace plb

#endif
