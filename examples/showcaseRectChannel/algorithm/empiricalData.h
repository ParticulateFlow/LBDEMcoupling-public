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

#ifndef EMPIRICAL_DATA_H
#define EMPIRICAL_DATA_H

#include "core/globalDefs.h"

namespace plb {

/// Returns the drag coefficient.
/// CD = 2 Fd/(rho u^2 r^2 pi)
/// Input: Re = 2 r u / nu
/// The paper to cite for this is:
/// @book{clift_bubbles_1978,
///	title = {Bubbles, {Drops}, and {Particles}},
///	isbn = {978-0-12-176950-5},
///	url = {https://books.google.ch/books?id=n8gRAQAAIAAJ},
///	publisher = {Academic Press},
///	author = {Clift, R. and Grace, J.R. and Weber, M.E. and Weber, M.F.},
///	year = {1978},
///	lccn = {77006592},
///	keywords = {bench-data}
///}
/// But I used the easier to read document
/// Drag on spherical particles and steady settling velocities
/// R. Shankar Subramanian
/// Department of Chemical and Biomolecular Engineering
/// Clarkson University

double empirical_sphere_drag(double Re);

double computeTerminalVelocity(
    double densityRatio, double vSphere, double r, double kinematicViscosity, double gravity = 9.8,
    bool doOutput = true);

double computeOptimalRadius(
    double targetRe, double densityRatio, double kinematicViscosity, double gravity = 9.8,
    bool doOutput = true);

}  // namespace plb

#endif  // EMPIRICAL_DATA_H
