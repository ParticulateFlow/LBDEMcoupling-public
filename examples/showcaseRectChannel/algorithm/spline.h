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

#ifndef SPLINE_H
#define SPLINE_H

#include <string>
#include <vector>

namespace plb {

template <typename T>
class Spline {
public:
    Spline() { }
    Spline(std::string fname);
    Spline(std::vector<T> const &x_, std::vector<T> const &y_);
    virtual ~Spline() { }
    virtual Spline<T> *clone() const = 0;
    std::vector<T> const &getAbscissae() const
    {
        return x;
    }
    std::vector<T> &getAbscissae()
    {
        return x;
    }
    std::vector<T> const &getOrdinates() const
    {
        return y;
    }
    std::vector<T> &getOrdinates()
    {
        return y;
    }
    virtual T getFunctionValue(T t) const = 0;
    virtual T getIntegralValue() const = 0;
    virtual T getIntegralValue(T t1, T t2) const = 0;

private:
    std::vector<T> x, y;
};

template <typename T>
class CubicSpline : public Spline<T> {
public:
    CubicSpline() : Spline<T>() { }
    CubicSpline(std::string fname) : Spline<T>(fname) { }
    CubicSpline(std::vector<T> const &x_, std::vector<T> const &y_) : Spline<T>(x_, y_) { }
    virtual ~CubicSpline() { }
    virtual CubicSpline<T> *clone() const = 0;
    virtual T getFunctionValue(T t) const = 0;
    virtual T getDerivativeValue(T t) const = 0;
    virtual T getSecondDerivativeValue(T t) const = 0;
    virtual T getThirdDerivativeValue(T t) const = 0;
    virtual T getIntegralValue() const = 0;
    virtual T getIntegralValue(T t1, T t2) const = 0;
};

template <typename T>
class NaturalCubicSpline : public CubicSpline<T> {
public:
    NaturalCubicSpline() : CubicSpline<T>(), icache(0) { }
    NaturalCubicSpline(std::string fname);
    NaturalCubicSpline(std::vector<T> const &x_, std::vector<T> const &y_);
    virtual ~NaturalCubicSpline() { }
    virtual NaturalCubicSpline<T> *clone() const;
    virtual T getFunctionValue(T t) const;
    virtual T getDerivativeValue(T t) const;
    virtual T getSecondDerivativeValue(T t) const;
    virtual T getThirdDerivativeValue(T t) const;
    virtual T getIntegralValue() const;
    virtual T getIntegralValue(T t1, T t2) const;

private:
    void constructSpline();
    plint bsrch(T t, plint il, plint ih) const;
    void locate(T &t) const;

private:
    mutable plint icache;
    std::vector<T> y1, y2, y3;
};

}  // namespace plb

#endif  // SPLINE_H
