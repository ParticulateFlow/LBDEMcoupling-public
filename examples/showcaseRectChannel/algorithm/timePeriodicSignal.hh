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

#ifndef TIME_PERIODIC_SIGNAL_HH
#define TIME_PERIODIC_SIGNAL_HH

#include <limits>

#include "algorithm/spline.h"
#include "algorithm/timePeriodicSignal.h"
#include "core/globalDefs.h"
#include "core/runTimeDiagnostics.h"
#include "core/util.h"

namespace plb {

template <typename T>
TimePeriodicSignal<T>::TimePeriodicSignal(std::string fname) : signal(fname)
{
    std::vector<T> const &t = signal.getAbscissae();
    std::vector<T> const &x = signal.getOrdinates();
    plint n = (plint)t.size();

    PLB_ASSERT(util::isZero(t[0]));
    PLB_ASSERT(util::fpequal(x[0], x[n - 1]));

    period = t[n - 1] - t[0];

    PLB_ASSERT(util::greaterThan(period, (T)0));
}

template <typename T>
TimePeriodicSignal<T>::TimePeriodicSignal(std::vector<T> const &t, std::vector<T> const &x) :
    signal(t, x)
{
    plint n = (plint)t.size();

    PLB_ASSERT(util::isZero(t[0]));
    PLB_ASSERT(util::fpequal(x[0], x[n - 1]));

    period = t[n - 1] - t[0];

    PLB_ASSERT(util::greaterThan(period, (T)0));
}

template <typename T>
TimePeriodicSignal<T> *TimePeriodicSignal<T>::clone() const
{
    return new TimePeriodicSignal<T>(*this);
}

template <typename T>
T TimePeriodicSignal<T>::getSignalValue(T t) const
{
    PLB_ASSERT(util::greaterEqual(t, (T)0));
    T trel = t - (plint)(t / period) * period;
    return signal.getFunctionValue(trel);
}

template <typename T>
T TimePeriodicSignal<T>::getDerivativeValue(T t) const
{
    PLB_ASSERT(util::greaterEqual(t, (T)0));
    T trel = t - (plint)(t / period) * period;
    return signal.getDerivativeValue(trel);
}

template <typename T>
T TimePeriodicSignal<T>::getSecondDerivativeValue(T t) const
{
    PLB_ASSERT(util::greaterEqual(t, (T)0));
    T trel = t - (plint)(t / period) * period;
    return signal.getSecondDerivativeValue(trel);
}

template <typename T>
T TimePeriodicSignal<T>::getThirdDerivativeValue(T t) const
{
    PLB_ASSERT(util::greaterEqual(t, (T)0));
    T trel = t - (plint)(t / period) * period;
    return signal.getThirdDerivativeValue(trel);
}

template <typename T>
T TimePeriodicSignal<T>::getIntegralValue() const
{
    return signal.getIntegralValue();
}

template <typename T>
T TimePeriodicSignal<T>::getIntegralValue(T tmin, T tmax) const
{
    PLB_ASSERT(util::lessEqual(tmin, tmax) && util::greaterEqual(tmin, (T)0));

    if (util::fpequal(tmin, tmax)) {
        return ((T)0);
    }

    T dt = tmax - tmin;
    plint k = (plint)(dt / period);

    T integral = (k == 0) ? (T)0 : k * signal.getIntegralValue();

    T tmin_rel = tmin - (plint)(tmin / period) * period;
    T tmax_rel = tmax - (plint)(tmax / period) * period;

    integral += signal.getIntegralValue(tmin_rel, tmax_rel);

    return integral;
}

}  // namespace plb

#endif  // TIME_PERIODIC_SIGNAL_HH
