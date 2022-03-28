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

#ifndef TRIANGLE_SET_GENERATOR_H
#define TRIANGLE_SET_GENERATOR_H

#include <vector>

#include "core/array.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "multiBlock/multiBlock3D.h"
#include "offLattice/triangleSet.h"

namespace plb {

/// Create and return an ellipsoid as a set of triangles. The center and the half-axes (a, b, c)
///   of the ellipsoid must be given as arguments. A minimum number of triangles for the
///   triangulation must be provided as well. This number is suggestive for the resolution. The
///   actual number of triangles can be greater than the one provided.
template <typename T>
TriangleSet<T> constructEllipsoid(
    Array<T, 3> const &center, T a, T b, T c, plint minNumOfTriangles, T eps = getEpsilon<T>());

/// Create and return a sphere as a set of triangles. The center and radius of the sphere
///   must be given as arguments. A minimum number of triangles for the sphere triangulation
///   must be provided as well. This number is suggestive for the resolution. The
///   actual number of triangles can be greater than the one provided.
template <typename T>
TriangleSet<T> constructSphere(
    Array<T, 3> const &center, T radius, plint minNumOfTriangles, T eps = getEpsilon<T>());

/// Create and return a half-sphere as a set of triangles. The center and radius of the half-sphere
///   must be given as arguments. A minimum number of triangles for the half-sphere triangulation
///   must be provided as well. This number is suggestive for the resolution. The
///   actual number of triangles can be greater than the one provided.
template <typename T>
TriangleSet<T> constructHalfSphere(
    Array<T, 3> const &center, T radius, plint minNumOfTriangles, T eps = getEpsilon<T>());

/// Create and return a 2D disk as a set of triangles. The center and radius of the disk
///   must be given as arguments. A minimum number of triangles for the disk triangulation
///   must be provided as well. This number is suggestive for the resolution. The
///   actual number of triangles can be greater than the one provided.
template <typename T>
TriangleSet<T> constructDisk(
    Array<T, 3> const &center, T radius, plint minNumOfTriangles, T eps = getEpsilon<T>());

/// Create and return a circular tapered cylinder. The axis of the cylinder is parallel
///   to the x-axis. The center of the inlet disk must be given as argument, as well as
///   the inlet and outer radii and the length of the object. "nAxial" is the
///   number of points in the axial direction before triangulation, and "nCirc"
///   is the number of points in the circumferential direction before triangulation.
///   The final number of points will be 2*nAxial-1 and 2*nCirc in the axial
///   and circumferential directions, respectively.
template <typename T>
TriangleSet<T> constructCylinder(
    Array<T, 3> const &inletCenter, T inletRadius, T outletRadius, T length, plint nAxial,
    plint nCirc, T eps = getEpsilon<T>());

template <typename T>
TriangleSet<T> constructCylinder(
    Array<T, 3> const &inletCenter, T inletRadius, T outletRadius, T length, plint nAxial,
    plint nCirc, std::vector<Array<T, 3> > &inletPoints, T eps = getEpsilon<T>());

template <typename T>
TriangleSet<T> constructCylinder(
    Array<T, 3> const &inletCenter, Array<T, 3> const &axis, T inletRadius, T outletRadius,
    T length, plint nAxial, plint nCirc, std::vector<Array<T, 3> > &inletPoints,
    T eps = getEpsilon<T>());

template <typename T>
TriangleSet<T> constructCuboidWithInletOutlet(
    Array<T, 3> const &lowerCorner, Array<T, 3> const &upperCorner,
    Array<plint, 3> const &nSegments, plint inletOutletDirection, T eps = getEpsilon<T>());

template <typename T>
TriangleSet<T> constructCuboidWithInletOutlet(
    Cuboid<T> const &cuboid, Array<plint, 3> const &nSegments, plint inletOutletDirection,
    T eps = getEpsilon<T>());

template <typename T>
TriangleSet<T> constructCuboid(
    Array<T, 3> const &lowerCorner, Array<T, 3> const &upperCorner,
    Array<plint, 3> const &nSegments, T eps = getEpsilon<T>());

template <typename T>
TriangleSet<T> constructCuboid(
    Cuboid<T> const &cuboid, Array<plint, 3> const &nSegments, T eps = getEpsilon<T>());

template <typename T>
TriangleSet<T> patchTubes(
    TriangleSet<T> const &geometryWithOpenings, plint sortDirection,
    std::vector<T> const &patchLengths);

/// Create and return a rectangle. The rectangle is on the x-y plane, and its lower left
///   corner is at the origin of the axes. Its sides have length "lx" and "ly", while
///   the number of points for the triangulation are "nx" and "ny" on the x and y axis,
///   respectively. This means that the total number of triangles is 2*(nx-1)*(ny-1).
template <typename T>
TriangleSet<T> constructRectangle(T lx, T ly, plint nx, plint ny, T eps = getEpsilon<T>());

/// Create and return a generically placed rectangle. The center of the rectangle is placed
///   at "center", and its normal is "normal". Its sides have length "lx" and "ly", while
///   the number of points for the triangulation are "nx" and "ny" on the x and y axis,
///   respectively. This means that the total number of triangles is 2*(nx-1)*(ny-1).
template <typename T>
TriangleSet<T> constructGenericRectangle(
    T lx, T ly, plint nx, plint ny, Array<T, 3> const &center, Array<T, 3> const &normal,
    T eps = getEpsilon<T>());

/// Create a strip of triangles. There are given two sets of points "from" and "to". These
///   must contain the same number of points. This function creates a set of triangles
///   by connecting points from "from" with points from "to".
template <typename T>
TriangleSet<T> constructStrip(
    std::vector<Array<T, 3> > const &from, std::vector<Array<T, 3> > const &to,
    T eps = getEpsilon<T>());

/// Create an annulus (closed strip) of triangles. The normal of the annulus is parallel
///   to the positive x-axis. The center must be given as argument, as well as the inner
///   and outer radii of the object. "nCirc" is the number of points in the circumferential
///   direction. Remember that if the annulus is constructed to fit a cylinder, then the "nCirc"
///   parameter given to this function must be twice the corresponding one given to the
///   constructCylinder function. If innerRadius = 0, then a disk is constructed.
template <typename T>
TriangleSet<T> constructAnnulus(
    Array<T, 3> const &center, T innerRadius, T outerRadius, plint nCirc, T eps = getEpsilon<T>());

/// Create a TriangleSet from the SparseBlockStructure of a given block. If highResolution is false,
///   then each bulk will be represented by a cuboid with one segment per coordinate direction,
///   otherwise the number of segments will reflect the number of cells in each bulk.
template <typename T>
TriangleSet<T> sparseBlockStructureToTriangleSet(
    MultiBlock3D const &block, bool highResolution = false, T eps = getEpsilon<T>());

}  // namespace plb

#endif  // TRIANGLE_SET_GENERATOR_H
