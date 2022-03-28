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

#ifndef SPLINE_HH
#define SPLINE_HH

#include <cstdio>
#include <limits>

#include "algorithm/spline.h"
#include "core/globalDefs.h"
#include "core/util.h"

namespace plb {

/* ***************** class Spline ***************************************** */

template <typename T>
Spline<T>::Spline(std::string fname)
{
    FILE *fp = fopen(fname.c_str(), "r");
    PLB_ASSERT(fp != 0);

    long double tmp0, tmp1;
    while (!feof(fp))
        if (fscanf(fp, "%Lf%Lf", &tmp0, &tmp1) != EOF) {
            x.push_back((T)tmp0);
            y.push_back((T)tmp1);
        }

    fclose(fp);

    PLB_ASSERT(x.size() >= 2);
    PLB_ASSERT(x.size() == y.size());
#ifdef PLB_DEBUG
    for (plint i = 1; i < (plint)x.size(); i++) {
        PLB_ASSERT(util::greaterThan(x[i], x[i - 1]));
    }
#endif
}

template <typename T>
Spline<T>::Spline(std::vector<T> const &x_, std::vector<T> const &y_) : x(x_), y(y_)
{
    PLB_ASSERT(x.size() >= 2);
    PLB_ASSERT(x.size() == y.size());
#ifdef PLB_DEBUG
    for (plint i = 1; i < (plint)x.size(); i++) {
        PLB_ASSERT(util::greaterThan(x[i], x[i - 1]));
    }
#endif
}

/* ***************** class NaturalCubicSpline ***************************************** */

template <typename T>
NaturalCubicSpline<T>::NaturalCubicSpline(std::string fname) : CubicSpline<T>(fname)
{
    icache = 0;
    constructSpline();
}

template <typename T>
NaturalCubicSpline<T>::NaturalCubicSpline(std::vector<T> const &x_, std::vector<T> const &y_) :
    CubicSpline<T>(x_, y_)
{
    icache = 0;
    constructSpline();
}

template <typename T>
NaturalCubicSpline<T> *NaturalCubicSpline<T>::clone() const
{
    return new NaturalCubicSpline<T>(*this);
}

template <typename T>
void NaturalCubicSpline<T>::constructSpline()
{
    std::vector<T> const &x = this->getAbscissae();
    std::vector<T> const &y = this->getOrdinates();
    plint n = (plint)x.size();
    std::vector<T> a1(n - 1), a2(n - 1), a3(n), a4(n), a5(n);

    y1.resize(n - 1);
    y2.resize(n);
    y3.resize(n - 1);

    // Calculate the spline coefficients for polynomials of the form:
    // S_j(x) = y_j + y1_j * (x - x_j) + y2_j * (x - x_j)^2 + y3_j * (x - x_j)^3

    for (plint i = 0; i < n - 1; i++)
        a1[i] = x[i + 1] - x[i];

    for (plint i = 1; i < n - 1; i++)
        a2[i] = 3.0 * (y[i + 1] - y[i]) / a1[i] - 3.0 * (y[i] - y[i - 1]) / a1[i - 1];

    a3[0] = 1.0;
    a4[0] = 0.0;
    a5[0] = 0.0;

    for (plint i = 1; i < n - 1; i++) {
        a3[i] = 2.0 * (x[i + 1] - x[i - 1]) - a1[i - 1] * a4[i - 1];
        a4[i] = a1[i] / a3[i];
        a5[i] = (a2[i] - a1[i - 1] * a5[i - 1]) / a3[i];
    }

    a3[n - 1] = 1.0;
    a5[n - 1] = 0.0;
    y2[n - 1] = 0.0;

    for (plint i = n - 2; i > -1; i--) {
        y2[i] = a5[i] - a4[i] * y2[i + 1];
        y1[i] = (y[i + 1] - y[i]) / a1[i] - a1[i] * (y2[i + 1] + 2.0 * y2[i]) / 3.0;
        y3[i] = (y2[i + 1] - y2[i]) / (3.0 * a1[i]);
    }

    y2.resize(n - 1);
}

// Return an index i in the range from il to ih for which x[i] <= t < x[i + 1]
template <typename T>
plint NaturalCubicSpline<T>::bsrch(T t, plint il, plint ih) const
{
    std::vector<T> const &x = this->getAbscissae();
    plint ilo = il;
    plint ihi = ih;

    plint i;
    while (ihi > (ilo + 1)) {
        i = (ihi + ilo) / 2;

        if (x[i] > t)
            ihi = i;
        else
            ilo = i;
    }

    return ilo;
}

template <typename T>
void NaturalCubicSpline<T>::locate(T &t) const
{
    std::vector<T> const &x = this->getAbscissae();
    plint n = (plint)x.size();
    if (util::lessEqual(t, x[0])) {
        t = x[0];
        icache = 0;
    } else if (util::greaterEqual(t, x[n - 1])) {
        t = x[n - 1];
        icache = n - 2;
    } else {
        if (t < x[icache]) {
            icache = bsrch(t, 0, icache);
        } else if (t >= x[icache + 1]) {
            icache = bsrch(t, icache, n - 1);
        }
    }
}

template <typename T>
T NaturalCubicSpline<T>::getFunctionValue(T t) const
{
    std::vector<T> const &x = this->getAbscissae();
    std::vector<T> const &y = this->getOrdinates();

    locate(t);

    T xtmp = t - x[icache];
    T ytmp = y[icache] + xtmp * (y1[icache] + xtmp * (y2[icache] + xtmp * (y3[icache])));
    return ytmp;
}

template <typename T>
T NaturalCubicSpline<T>::getDerivativeValue(T t) const
{
    std::vector<T> const &x = this->getAbscissae();

    locate(t);

    T xtmp = t - x[icache];
    T ydtmp = y1[icache] + xtmp * (2.0 * y2[icache] + xtmp * (3.0 * y3[icache]));
    return ydtmp;
}

template <typename T>
T NaturalCubicSpline<T>::getSecondDerivativeValue(T t) const
{
    std::vector<T> const &x = this->getAbscissae();

    locate(t);

    T xtmp = t - x[icache];
    T yd2tmp = 2.0 * y2[icache] + 6.0 * y3[icache] * xtmp;
    return yd2tmp;
}

template <typename T>
T NaturalCubicSpline<T>::getThirdDerivativeValue(T t) const
{
    locate(t);

    T yd3tmp = 6.0 * y3[icache];
    return yd3tmp;
}

template <typename T>
T NaturalCubicSpline<T>::getIntegralValue() const
{
    std::vector<T> const &x = this->getAbscissae();
    std::vector<T> const &y = this->getOrdinates();
    plint n = (plint)x.size();

    T integral = (T)0;
    for (plint i = 0; i < n - 1; i++) {
        T dx = x[i + 1] - x[i];
        integral += dx * (y[i] + dx * (y1[i] / 2.0 + dx * (y2[i] / 3.0 + dx * (y3[i] / 4.0))));
    }

    return integral;
}

template <typename T>
T NaturalCubicSpline<T>::getIntegralValue(T t1, T t2) const
{
    if (util::fpequal(t1, t2)) {
        return ((T)0);
    }

    std::vector<T> const &x = this->getAbscissae();
    std::vector<T> const &y = this->getOrdinates();
    plint n = (plint)x.size();

    T sign = (T)1;
    if (util::lessThan(t2, t1)) {
        std::swap(t1, t2);
        sign = (T)-1;
    }

    T leftIntegral = (T)0;
    if (util::lessEqual(t1, x[0])) {
        if (util::lessEqual(t2, x[0])) {
            return sign * (t2 - t1) * y[0];
        } else {
            leftIntegral = (x[0] - t1) * y[0];
        }
    }

    T rightIntegral = (T)0;
    if (util::greaterEqual(t2, x[n - 1])) {
        if (util::greaterEqual(t1, x[n - 1])) {
            return sign * (t2 - t1) * y[n - 1];
        } else {
            rightIntegral = (t2 - x[n - 1]) * y[n - 1];
        }
    }

    locate(t1);
    plint imin = icache;

    locate(t2);
    plint imax = icache;

    T middleIntegral = (T)0;

    if (imin == imax) {
        plint i = imin;
        T dxmin = t1 - x[i];
        T dxmax = t2 - x[i];
        middleIntegral =
            y[i] * (t2 - t1) + (y1[i] / 2.0) * (dxmax * dxmax - dxmin * dxmin)
            + (y2[i] / 3.0) * (dxmax * dxmax * dxmax - dxmin * dxmin * dxmin)
            + (y3[i] / 4.0) * (dxmax * dxmax * dxmax * dxmax - dxmin * dxmin * dxmin * dxmin);
    } else {
        plint i = imin;
        T dxmin = t1 - x[i];
        T dxmax = x[i + 1] - x[i];
        middleIntegral =
            y[i] * (x[i + 1] - t1) + (y1[i] / 2.0) * (dxmax * dxmax - dxmin * dxmin)
            + (y2[i] / 3.0) * (dxmax * dxmax * dxmax - dxmin * dxmin * dxmin)
            + (y3[i] / 4.0) * (dxmax * dxmax * dxmax * dxmax - dxmin * dxmin * dxmin * dxmin);

        T dx;
        for (i = imin + 1; i < imax; i++) {
            dx = x[i + 1] - x[i];
            middleIntegral +=
                dx * (y[i] + dx * (y1[i] / 2.0 + dx * (y2[i] / 3.0 + dx * (y3[i] / 4.0))));
        }

        i = imax;
        dx = t2 - x[i];
        middleIntegral +=
            dx * (y[i] + dx * (y1[i] / 2.0 + dx * (y2[i] / 3.0 + dx * (y3[i] / 4.0))));
    }

    return sign * (leftIntegral + middleIntegral + rightIntegral);
}

}  // namespace plb

#endif  // SPLINE_HH
