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

#ifndef INTERPOLATION_2D_H
#define INTERPOLATION_2D_H

#include <string>
#include <vector>

namespace plb {

template <typename T>
class Interpolation2D {
public:
    Interpolation2D() { }
    Interpolation2D(std::string fname);
    Interpolation2D(
        std::vector<T> const &x_, std::vector<T> const &y_, std::vector<std::vector<T> > const &f_);
    virtual ~Interpolation2D() { }
    virtual Interpolation2D<T> *clone() const = 0;
    std::vector<T> const &getCoordinatesX() const
    {
        return x;
    }
    std::vector<T> &getCoordinatesX()
    {
        return x;
    }
    std::vector<T> const &getCoordinatesY() const
    {
        return y;
    }
    std::vector<T> &getCoordinatesY()
    {
        return y;
    }
    std::vector<std::vector<T> > const &getFunctionValues() const
    {
        return f;
    }
    std::vector<std::vector<T> > &getFunctionValues()
    {
        return f;
    }
    void transformCoordinatesX(T scale, T offset);
    void transformCoordinatesY(T scale, T offset);
    void transformFunctionValues(T scale, T offset);
    void exportInXYZFormat(std::string fname, int numDecimalDigits = 10) const;
    virtual T getFunctionValue(T x0, T y0) const = 0;

private:
    std::vector<T> x, y;
    std::vector<std::vector<T> > f;
};

template <typename T>
class LinearInterpolation2D : public Interpolation2D<T> {
public:
    LinearInterpolation2D() : Interpolation2D<T>(), icache(0), jcache(0) { }
    LinearInterpolation2D(std::string fname) : Interpolation2D<T>(fname), icache(0), jcache(0) { }
    LinearInterpolation2D(
        std::vector<T> const &x_, std::vector<T> const &y_,
        std::vector<std::vector<T> > const &f_) :
        Interpolation2D<T>(x_, y_, f_), icache(0), jcache(0)
    { }
    virtual ~LinearInterpolation2D() { }
    virtual LinearInterpolation2D<T> *clone() const;
    virtual T getFunctionValue(T x0, T y0) const;

private:
    plint bsrch(std::vector<T> const &v, T t, plint il, plint ih) const;
    void locate(std::vector<T> const &v, T &t, plint &cache) const;

private:
    mutable plint icache, jcache;
};

}  // namespace plb

#endif  // INTERPOLATION_2D_H
