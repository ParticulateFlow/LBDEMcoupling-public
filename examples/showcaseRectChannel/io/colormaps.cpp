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

#include "io/colormaps.h"

#include <cmath>

namespace plb {

LinearFunction::LinearFunction(double x1_, double x2_, double y1_, double y2_) :
    x1(x1_), x2(x2_), y1(y1_), y2(y2_)
{ }

double LinearFunction::operator()(double x) const
{
    return ((y2 - y1) * x + x2 * y1 - x1 * y2) / (x2 - x1);
}

LinearFunction *LinearFunction::clone() const
{
    return new LinearFunction(*this);
}

PowerLawFunction::PowerLawFunction(double x1_, double x2_, double y1_, double y2_, double b_) :
    x1(std::pow(x1_, b_)), x2(std::pow(x2_, b_)), y1(y1_), y2(y2_), b(b_)
{ }

double PowerLawFunction::operator()(double x) const
{
    return ((y2 - y1) * std::pow(x, b) + x2 * y1 - x1 * y2) / (x2 - x1);
}

PowerLawFunction *PowerLawFunction::clone() const
{
    return new PowerLawFunction(*this);
}

PiecewiseFunction::PiecewiseFunction(PiecewiseFunction const &rhs) :
    pieces(rhs.pieces), functions(rhs.functions.size())
{
    for (pluint iF = 0; iF < functions.size(); ++iF) {
        functions[iF] = rhs.functions[iF]->clone();
    }
}

PiecewiseFunction &PiecewiseFunction::operator=(PiecewiseFunction const &rhs)
{
    PiecewiseFunction(rhs).swap(*this);
    return *this;
}

PiecewiseFunction::~PiecewiseFunction()
{
    for (pluint iF = 0; iF < functions.size(); ++iF) {
        delete functions[iF];
    }
}

void PiecewiseFunction::swap(PiecewiseFunction &rhs)
{
    pieces.swap(rhs.pieces);
    functions.swap(rhs.functions);
}

void PiecewiseFunction::addPiece(Piece piece, ScalarFunction *f)
{
    std::vector<Piece>::iterator pieceIt = pieces.begin();
    std::vector<ScalarFunction *>::iterator fIt = functions.begin();
    while (pieceIt != pieces.end() && piece.closedBegin >= pieceIt->closedBegin) {
        ++pieceIt;
        ++fIt;
    }
    pieces.insert(pieceIt, piece);
    functions.insert(fIt, f);
}

double PiecewiseFunction::operator()(double x) const
{
    if (pieces.empty() || x < pieces[0].closedBegin) {
        return double();
    }
    pluint iPiece = 0;
    while (iPiece != pieces.size() && x >= pieces[iPiece].openEnd) {
        ++iPiece;
    }
    if (iPiece == pieces.size() || x < pieces[iPiece].closedBegin) {
        return double();
    }
    return (*functions[iPiece])(x);
}

PiecewiseFunction *PiecewiseFunction::clone() const
{
    return new PiecewiseFunction(*this);
}

ColorMap::ColorMap(
    PiecewiseFunction const &red_, PiecewiseFunction const &green_,
    PiecewiseFunction const &blue_) :
    red(red_), green(green_), blue(blue_)
{ }

rgb ColorMap::get(double x) const
{
    return rgb(red(x), green(x), blue(x));
}

namespace mapGenerators {

PiecewiseFunction generateEarthRed()
{
    double p0 = 0.;
    double p1 = 3. / 8.;
    double p2 = 6. / 8.;
    double p3 = 1.;

    PiecewiseFunction earthRed;
    earthRed.addPiece(Piece(p0, p1), new PowerLawFunction(p0, p1, 0., 0.8, 0.6));
    earthRed.addPiece(Piece(p1, p2), new PowerLawFunction(p1, p2, 0.8, 0.9, 0.9));
    earthRed.addPiece(Piece(p2, p3), new PowerLawFunction(p2, p3, 0.9, 1.0, 0.2));

    return earthRed;
}

PiecewiseFunction generateEarthGreen()
{
    double p0 = 0.;
    double p1 = 3. / 8.;
    double p2 = 6. / 8.;
    double p3 = 1.;

    PiecewiseFunction earthGreen;
    earthGreen.addPiece(Piece(p0, p1), new PowerLawFunction(p0, p1, 0., 0.5, 0.6));
    earthGreen.addPiece(Piece(p1, p2), new PowerLawFunction(p1, p2, 0.5, 0.9, 0.2));
    earthGreen.addPiece(Piece(p2, p3), new PowerLawFunction(p2, p3, 0.9, 1.0, 0.2));

    return earthGreen;
}

PiecewiseFunction generateEarthBlue()
{
    double p0 = 0.;
    double p1 = 3. / 8.;
    double p2 = 6. / 8.;
    double p3 = 1.;

    PiecewiseFunction earthBlue;
    earthBlue.addPiece(Piece(p0, p1), new PowerLawFunction(p0, p1, 0., 0.5, 0.6));
    earthBlue.addPiece(Piece(p1, p2), new PowerLawFunction(p1, p2, 0.5, 0.7, 0.2));
    earthBlue.addPiece(Piece(p2, p3), new PowerLawFunction(p2, p3, 0.7, 1.0, 0.2));

    return earthBlue;
}

PiecewiseFunction generateWaterRed()
{
    double p0 = 0.;
    double p1 = 3. / 8.;
    double p2 = 6. / 8.;
    double p3 = 1.;

    PiecewiseFunction waterRed;
    waterRed.addPiece(Piece(p0, p1), new PowerLawFunction(p0, p1, 0., 0.5, 0.6));
    waterRed.addPiece(Piece(p1, p2), new PowerLawFunction(p1, p2, 0.5, 0.7, 0.2));
    waterRed.addPiece(Piece(p2, p3), new PowerLawFunction(p2, p3, 0.7, 1.0, 0.2));

    return waterRed;
}

PiecewiseFunction generateWaterGreen()
{
    double p0 = 0.;
    double p1 = 3. / 8.;
    double p2 = 6. / 8.;
    double p3 = 1.;

    PiecewiseFunction waterGreen;
    waterGreen.addPiece(Piece(p0, p1), new PowerLawFunction(p0, p1, 0., 0.5, 0.6));
    waterGreen.addPiece(Piece(p1, p2), new PowerLawFunction(p1, p2, 0.5, 0.9, 0.2));
    waterGreen.addPiece(Piece(p2, p3), new PowerLawFunction(p2, p3, 0.9, 1.0, 0.2));

    return waterGreen;
}

PiecewiseFunction generateWaterBlue()
{
    double p0 = 0.;
    double p1 = 3. / 8.;
    double p2 = 6. / 8.;
    double p3 = 1.;

    PiecewiseFunction waterBlue;
    waterBlue.addPiece(Piece(p0, p1), new PowerLawFunction(p0, p1, 0., 0.8, 0.6));
    waterBlue.addPiece(Piece(p1, p2), new PowerLawFunction(p1, p2, 0.8, 0.9, 0.9));
    waterBlue.addPiece(Piece(p2, p3), new PowerLawFunction(p2, p3, 0.9, 1.0, 0.2));

    return waterBlue;
}

PiecewiseFunction generateAirRed()
{
    double p0 = 0.;
    double p1 = 1.;

    PiecewiseFunction airRed;
    airRed.addPiece(Piece(p0, p1), new LinearFunction(p0, p1, 0., 1.));

    return airRed;
}

PiecewiseFunction generateAirGreen()
{
    double p0 = 0.;
    double p1 = 1.;

    PiecewiseFunction airGreen;
    airGreen.addPiece(Piece(p0, p1), new LinearFunction(p0, p1, 1., 0.));

    return airGreen;
}

PiecewiseFunction generateAirBlue()
{
    double p0 = 0.;
    double p1 = 1.;

    PiecewiseFunction airBlue;
    airBlue.addPiece(Piece(p0, p1), new LinearFunction(p0, p1, 1., 1.));

    return airBlue;
}

PiecewiseFunction generateFireRed()
{
    double p0 = 0.;
    double p1 = 0.36;
    double p3 = 1.;

    PiecewiseFunction fireRed;
    fireRed.addPiece(Piece(p0, p1), new LinearFunction(p0, p1, 0., 1.));
    fireRed.addPiece(Piece(p1, p3), new LinearFunction(p1, p3, 1., 1.));

    return fireRed;
}

PiecewiseFunction generateFireGreen()
{
    double p0 = 0.;
    double p1 = 0.36;
    double p2 = 0.75;
    double p3 = 1.;

    PiecewiseFunction fireGreen;
    fireGreen.addPiece(Piece(p0, p1), new LinearFunction(p0, p1, 0., 0.));
    fireGreen.addPiece(Piece(p1, p2), new LinearFunction(p1, p2, 0., 1.));
    fireGreen.addPiece(Piece(p2, p3), new LinearFunction(p2, p3, 1., 1.));

    return fireGreen;
}

PiecewiseFunction generateFireBlue()
{
    double p0 = 0.;
    double p2 = 0.75;
    double p3 = 1.;

    PiecewiseFunction fireBlue;
    fireBlue.addPiece(Piece(p0, p2), new LinearFunction(p0, p2, 0., 0.));
    fireBlue.addPiece(Piece(p2, p3), new LinearFunction(p2, p3, 0., 1.));

    return fireBlue;
}

PiecewiseFunction generateLeeLooRed()
{
    double p0 = 0.;
    double p2 = 3. / 8.;
    double p3 = 5. / 8.;
    double p4 = 7. / 8.;
    double p5 = 1.;
    double p6 = 9. / 8.;

    PiecewiseFunction leeLooRed;
    leeLooRed.addPiece(Piece(p0, p2), new LinearFunction(p0, p2, 0., 0.));
    leeLooRed.addPiece(Piece(p2, p3), new LinearFunction(p2, p3, 0., 1.));
    leeLooRed.addPiece(Piece(p3, p4), new LinearFunction(p3, p4, 1., 1.));
    leeLooRed.addPiece(Piece(p4, p5), new LinearFunction(p4, p6, 1., 0.));

    return leeLooRed;
}

PiecewiseFunction generateLeeLooGreen()
{
    double p0 = 0.;
    double p1 = 1. / 8.;
    double p2 = 3. / 8.;
    double p3 = 5. / 8.;
    double p4 = 7. / 8.;
    double p5 = 1.;
    double p6 = 9 / 8;

    PiecewiseFunction leeLooGreen;
    leeLooGreen.addPiece(Piece(p0, p1), new LinearFunction(p0, p1, 0., 0.));
    leeLooGreen.addPiece(Piece(p1, p2), new LinearFunction(p1, p2, 0., 1.));
    leeLooGreen.addPiece(Piece(p2, p3), new LinearFunction(p2, p3, 1., 1.));
    leeLooGreen.addPiece(Piece(p3, p4), new LinearFunction(p3, p4, 1., 0.));
    leeLooGreen.addPiece(Piece(p4, p5), new LinearFunction(p4, p6, 0., 0.));

    return leeLooGreen;
}

PiecewiseFunction generateLeeLooBlue()
{
    double pm1 = -1. / 8.;
    double p0 = 0.;
    double p1 = 1. / 8.;
    double p2 = 3. / 8.;
    double p3 = 5. / 8.;
    double p5 = 1.;

    PiecewiseFunction leeLooBlue;
    leeLooBlue.addPiece(Piece(p0, p1), new LinearFunction(pm1, p1, 0., 1.));
    leeLooBlue.addPiece(Piece(p1, p2), new LinearFunction(p1, p2, 1., 1.));
    leeLooBlue.addPiece(Piece(p2, p3), new LinearFunction(p2, p3, 1., 0.));
    leeLooBlue.addPiece(Piece(p3, p5), new LinearFunction(p3, p5, 0., 0.));

    return leeLooBlue;
}

PiecewiseFunction generateRedBlueRed()
{
    double p0 = 0.;
    double p1 = 0.45;
    double p2 = 0.5;
    double p3 = 0.55;
    double p4 = 1.;
    double damping = 0.4;

    PiecewiseFunction rbRed;
    rbRed.addPiece(Piece(p0, p1), new LinearFunction(p0, p1, 0., 0.45));
    rbRed.addPiece(Piece(p1, p2), new PowerLawFunction(p1, p2, 0.45, damping * 0.5, 50.0));
    rbRed.addPiece(Piece(p2, p3), new PowerLawFunction(p2, p3, damping * 0.5, 0.55, -25.0));
    rbRed.addPiece(Piece(p3, p4), new LinearFunction(p3, p4, 0.55, 1.));

    return rbRed;
}

PiecewiseFunction generateRedBlueGreen()
{
    double p0 = 0.;
    double p1 = 1.;

    PiecewiseFunction rbGreen;
    rbGreen.addPiece(Piece(p0, p1), new LinearFunction(p0, p1, 0., 0.));

    return rbGreen;
}

PiecewiseFunction generateRedBlueBlue()
{
    double p0 = 0.;
    double p1 = 0.45;
    double p2 = 0.5;
    double p3 = 0.55;
    double p4 = 1.;
    double damping = 0.4;

    PiecewiseFunction rbBlue;
    rbBlue.addPiece(Piece(p0, p1), new LinearFunction(p0, p1, 1., 0.55));
    rbBlue.addPiece(Piece(p1, p2), new PowerLawFunction(p1, p2, 0.55, damping * 0.5, 25.0));
    rbBlue.addPiece(Piece(p2, p3), new PowerLawFunction(p2, p3, damping * 0.5, 0.45, -50.0));
    rbBlue.addPiece(Piece(p3, p4), new LinearFunction(p3, p4, 0.45, 0.));

    return rbBlue;
}

ColorMap generateMap(std::string mapName)
{
    if (mapName == "earth") {
        return ColorMap(generateEarthRed(), generateEarthGreen(), generateEarthBlue());
    } else if (mapName == "water") {
        return ColorMap(generateWaterRed(), generateWaterGreen(), generateWaterBlue());
    } else if (mapName == "air") {
        return ColorMap(generateAirRed(), generateAirGreen(), generateAirBlue());
    } else if (mapName == "fire") {
        return ColorMap(generateFireRed(), generateFireGreen(), generateFireBlue());
    } else if (mapName == "leeloo") {
        return ColorMap(generateLeeLooRed(), generateLeeLooGreen(), generateLeeLooBlue());
    } else if (mapName == "redblue") {
        return ColorMap(generateRedBlueRed(), generateRedBlueGreen(), generateRedBlueBlue());
    }
    return ColorMap(generateLeeLooRed(), generateLeeLooGreen(), generateLeeLooBlue());
}

}  // namespace mapGenerators

}  // namespace plb
