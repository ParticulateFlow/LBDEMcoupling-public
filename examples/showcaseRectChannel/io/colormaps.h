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

#ifndef COLORMAPS_H
#define COLORMAPS_H

#include <string>
#include <vector>

#include "core/globalDefs.h"

namespace plb {

struct ScalarFunction {
    virtual ~ScalarFunction() { }
    virtual double operator()(double x) const = 0;
    virtual ScalarFunction *clone() const = 0;
};

class LinearFunction : public ScalarFunction {
public:
    LinearFunction(double x1_, double x2_, double y1_, double y2_);
    virtual double operator()(double x) const;
    virtual LinearFunction *clone() const;

private:
    double x1, x2, y1, y2;
};

class PowerLawFunction : public ScalarFunction {
public:
    PowerLawFunction(double x1_, double x2_, double y1_, double y2_, double b_);
    virtual double operator()(double x) const;
    virtual PowerLawFunction *clone() const;

private:
    double x1, x2, y1, y2;
    double b;
};

struct Piece {
    Piece(double closedBegin_, double openEnd_) : closedBegin(closedBegin_), openEnd(openEnd_) { }
    double closedBegin, openEnd;
};

class PiecewiseFunction : public ScalarFunction {
public:
    PiecewiseFunction() { }
    ~PiecewiseFunction();
    PiecewiseFunction(PiecewiseFunction const &rhs);
    PiecewiseFunction &operator=(PiecewiseFunction const &rhs);
    void swap(PiecewiseFunction &rhs);
    void addPiece(Piece piece, ScalarFunction *f);
    virtual double operator()(double x) const;
    virtual PiecewiseFunction *clone() const;

private:
    std::vector<Piece> pieces;
    std::vector<ScalarFunction *> functions;
};

struct rgb {
    rgb(double r_, double g_, double b_) : r(r_), g(g_), b(b_) { }
    double r, g, b;
};

class ColorMap {
public:
    ColorMap(
        PiecewiseFunction const &red_, PiecewiseFunction const &green_,
        PiecewiseFunction const &blue_);
    rgb get(double x) const;

private:
    PiecewiseFunction red, green, blue;
};

namespace mapGenerators {

PiecewiseFunction generateEarthRed();
PiecewiseFunction generateEarthGreen();
PiecewiseFunction generateEarthBlue();

PiecewiseFunction generateWaterRed();
PiecewiseFunction generateWaterGreen();
PiecewiseFunction generateWaterBlue();

PiecewiseFunction generateAirRed();
PiecewiseFunction generateAirGreen();
PiecewiseFunction generateAirBlue();

PiecewiseFunction generateFireRed();
PiecewiseFunction generateFireGreen();
PiecewiseFunction generateFireBlue();

PiecewiseFunction generateLeeLooRed();
PiecewiseFunction generateLeeLooGreen();
PiecewiseFunction generateLeeLooBlue();

PiecewiseFunction generateRedBlueRed();
PiecewiseFunction generateRedBlueGreen();
PiecewiseFunction generateRedBlueBlue();

ColorMap generateMap(std::string mapName);

}  // namespace mapGenerators

}  // namespace plb

#endif  // COLORMAPS_H
