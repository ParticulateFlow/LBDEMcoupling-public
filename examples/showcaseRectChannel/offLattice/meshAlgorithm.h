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

#ifndef MESH_ALGORITHM_H
#define MESH_ALGORITHM_H

/*
 * Algorithms, applied to generic meshes.
 */

#include <string>
#include <vector>

#include "core/globalDefs.h"
#include "core/units.h"
#include "offLattice/triangleMesh.h"
#include "offLattice/triangleSet.h"

namespace plb {

/// Convert to lattice units, do not recenter.
template <typename T>
void toLatticeUnits(
    TriangleMesh<T> &mesh, T dx, Array<T, 3> const &offset = Array<T, 3>(T(), T(), T()));

/// Convert to lattice units, according to the physical and LBM domains of the Units3D object.
template <typename T>
void toLatticeUnits(TriangleMesh<T> &mesh, Units3D const &units);

/// Convert to lattice units, according to the physical and LBM domains of the Units3D object.
template <typename T>
void toLatticeUnits(TriangleSet<T> &triangleSet, Units3D const &units);

/// Convert to lattice units and move the center-of-mass to the specified position.
template <typename T>
void toLatticeUnits(
    TriangleMesh<T> &mesh, Array<T, 3> const &center, T dx,
    Array<T, 3> const &offset = Array<T, 3>(T(), T(), T()));

template <typename T>
void toPhysicalUnits(TriangleMesh<T> &mesh, T dx, Array<T, 3> offset = Array<T, 3>(T(), T(), T()));

/// Convert from lattice to physical units, according to the physical and LBM domains of the Units3D
/// object.
template <typename T>
void toPhysicalUnits(TriangleMesh<T> &mesh, Units3D const &units);

/// Convert from lattice to physical units, according to the physical and LBM domains of the Units3D
/// object.
template <typename T>
void toPhysicalUnits(TriangleSet<T> &triangleSet, Units3D const &units);

template <typename T>
void smooth(
    RawConnectedTriangleMesh<T> &triangleMesh, plint maxiter = 1, T relax = (T)1.0,
    bool isMeasureWeighted = false);

template <typename T>
Array<T, 3> centerOfMass(TriangleMesh<T> &mesh);

/// Executes a rotation to the obstacle which, if applied to the
///   x- and y- axes, would apply them to the new x- and y- axes
///   as provided.
template <typename T>
void rotateAxesTo(TriangleMesh<T> &mesh, Array<T, 3> const &ex, Array<T, 3> const &ey);

/// Executes a rotation to the obstacle which, if applied to the
///   provided new x- and y- axes, would rotate them back to the
///   x- and y-axes of the origin.
template <typename T>
void rotateAxesFrom(TriangleMesh<T> &mesh, Array<T, 3> const &ex, Array<T, 3> const &ey);

}  // namespace plb

#endif  // MESH_ALGORITHM_H
