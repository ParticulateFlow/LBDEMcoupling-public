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

#include "algorithm/empiricalData.h"

#include <cmath>

#include "algorithm/nonlinearEquationSolvers.h"
#include "algorithm/nonlinearEquationSolvers.hh"
#include "io/parallelIO.h"

namespace plb {

double empirical_sphere_drag(double Re)
{
    double w = std::log10(Re);
    double CD = 1.0;
    if (Re <= 0.01) {
        CD = 9. / 2. + 24. / Re;
    } else if (Re <= 20.) {
        CD = 24. / Re * (1.0 + 0.1315 * std::pow(Re, 0.82 - 0.05 * w));
    } else if (Re <= 260) {
        CD = 24. / Re * (1.0 + 0.1935 * std::pow(Re, 0.6305));
    } else if (Re <= 1.5e3) {
        CD = std::pow(10.0, 1.6435 - 1.1242 * w + 0.1558 * w * w);
    } else if (Re <= 1.2e4) {
        CD = std::pow(10.0, -2.4571 + 2.5558 * w - 0.9295 * w * w + 0.1049 * w * w * w);
    } else if (Re <= 4.4e4) {
        CD = std::pow(10.0, -1.9181 + 0.6370 * w - 0.0636 * w * w);
    } else if (Re <= 3.38e5) {
        CD = std::pow(10.0, -4.3390 + 1.5809 * w - 0.1546 * w * w);
    } else if (Re <= 4e5) {
        CD = 29.78 - 5.3 * w;
    } else if (Re <= 1e6) {
        CD = 0.1 * w - 0.49;
    } else if (Re > 1e6) {
        CD = 0.19 - 8e4 / Re;
    }
    return CD;
}

class TerminalVelocity {
public:
    TerminalVelocity(
        double densityRatio_, double volume_, double r_, double kinematicViscosity_,
        double gravity_) :
        densityRatio(densityRatio_),
        volume(volume_),
        r(r_),
        kinematicViscosity(kinematicViscosity_),
        gravity(gravity_)
    { }
    double operator()(double vel) const
    {
        static const double pi = std::acos(-1.0);
        // Reynolds number is defined with respect to the sphere diameter.
        double Re = 2. * r * vel / kinematicViscosity;
        return vel
               - sqrt(
                   densityRatio * 2. / empirical_sphere_drag(Re) * volume * gravity / (r * r * pi));
    }

private:
    double densityRatio;
    double volume;
    double r;
    double kinematicViscosity;
    double gravity;
};

double computeTerminalVelocity(
    double densityRatio, double vSphere, double r, double kinematicViscosity, double gravity,
    bool doOutput)
{
    double u0 = 1.e-4;
    double u1 = 1.e4;
    double u_acc = 1.e-10;
    plint maxIter = 10000;
    double result;
    brentSolve(
        TerminalVelocity(densityRatio, vSphere, r, kinematicViscosity, gravity), u0, u1, u_acc,
        maxIter, result);
    if (doOutput) {
        pcout << "Brent solver for terminal velocity: error is "
              << std::fabs(TerminalVelocity(
                     densityRatio, vSphere, r, kinematicViscosity, gravity)(result))
              << std::endl;
    }
    return result;
}

class OptimalRadius {
public:
    OptimalRadius(
        double targetRe_, double densityRatio_, double kinematicViscosity_, double gravity_) :
        targetRe(targetRe_),
        densityRatio(densityRatio_),
        kinematicViscosity(kinematicViscosity_),
        gravity(gravity_)
    { }
    double operator()(double r) const
    {
        static const double pi = acos(-1.);
        double volume = 4. / 3. * pi * std::pow(r, 3.0);
        bool doOutput = false;
        return targetRe
               - 2. * r
                     * computeTerminalVelocity(
                         densityRatio, volume, r, kinematicViscosity, gravity, doOutput)
                     / kinematicViscosity;
    }

private:
    double targetRe;
    double densityRatio;
    double kinematicViscosity;
    double gravity;
};

double computeOptimalRadius(
    double targetRe, double densityRatio, double kinematicViscosity, double gravity, bool doOutput)
{
    double r0 = 1.e-8;
    double r1 = 1.0;
    double r_acc = 1.e-10;
    plint maxIter = 10000;
    double result;
    brentSolve(
        OptimalRadius(targetRe, densityRatio, kinematicViscosity, gravity), r0, r1, r_acc, maxIter,
        result);
    if (doOutput) {
        pcout << "Brent solver for terminal radius: error is "
              << std::fabs(
                     OptimalRadius(targetRe, densityRatio, kinematicViscosity, gravity)(result))
              << std::endl;
    }
    return result;
}

}  // namespace plb
