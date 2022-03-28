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

/** \file Constants and definitions for boundary conditions -- header file.  */

#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

namespace plb {

namespace boundary {
/** dirichlet:     impose a velocity (* obsolete ...or a density *).
 *  neumann:       zero-gradient for all velocity components or for the density.
 *  freeslip:      zero-gradient for tangential velocity components, and zero for normal ones.
 *  density:       fixed density, Neumann for tangential velocity components.
 *  (*obsolete*) outflow:       zero-gradient for all velocity components.
 *  (*obsolete*) normalOutflow: zero-gradient for normal velocity components, and zero for
 *tangential ones.
 **/
typedef enum { dirichlet, neumann, freeslip, density, outflow, normalOutflow } BcType;
}  // namespace boundary

}  // namespace plb

#endif
