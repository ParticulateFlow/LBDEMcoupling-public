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

#ifndef TRIANGLE_UTIL_H
#define TRIANGLE_UTIL_H

#include "offLattice/triangleMesh.h"

namespace plb {

/// Function to compute the intersection between a triangle and a line segment
///   between points "point1" and "point2" (if flag = 0), or between a triangle
///   and a half-line starting at "point1" and containing "point2" (if flag = 1),
///   or between a triangle and a whole line containing both "point1" and
///   "point2" (if flag = 2). "intersection", "normal" and "distance" are
///   objects whose states are changed by this function. These states are undefined,
///   and cannot be used by the caller function, when the return value of this
///   function is not 1.
///   Returns 1 if an intersection is found and the intersection is inside
///   the triangle, 0 if there is no intersection or the intersection is
///   outside the triangle, and -1 is the line is parallel to the triangle
///   and intersects with an infinity of points.
template <typename T>
int pointOnTriangle(
    Array<T, 3> const &point1, Array<T, 3> const &point2, int flag,
    typename TriangleMesh<T>::Triangle const &triangle, Array<T, 3> &intersection,
    Array<T, 3> &normal, T &distance);

/// The following function is a specialized version of the previous one. It simply checks
///   if a line segment between points "point1" and "point2" intersects a triangle. This
///   function is written for optimization purposes, and returns "true" if an intersection
///   is found and the intersection is inside the triangle, or "false" if there is no
///   intersection or the intersection is outside the triangle or if the line segment
///   belongs to the triangle and intersects with an infinity of points.
template <typename T>
bool segmentIntersectsTriangle(
    Array<T, 3> const &point1, Array<T, 3> const &point2,
    typename TriangleMesh<T>::Triangle const &triangle);

/// Compute the distance between a point and a line which goes through a
///   given edge. Returns also a flag indicating whether the location on
///   the line which is closest to the point is inside the edge or not.
template <typename T>
void distanceToEdgeLine(
    Array<T, 3> const &point, typename TriangleMesh<T>::Triangle const &triangle, plint whichEdge,
    T &distance, bool &intersectionIsInside);

/// Compute the distance between a point and a plane which goes through a
///   given triangle. Returns two flags. The first indicates whether the
///   location on the plane which is closest to the point is inside the
///   triangle or not. The second indicates whether the point is "behind"
///   the plane, i.e. on the side which is opposed to the triangle normal.
template <typename T>
void distanceToTrianglePlane(
    Array<T, 3> const &point, typename TriangleMesh<T>::Triangle const &triangle, T &distance,
    bool &intersectionIsInside, bool &pointIsBehind);

/// Compute the distance between the point and a triangle. The location
///   on the triangle the distance to which is being computed is chosen
///   as follows. The point is first projected on the plane defined by
///   the triangle. If the obtained location is inside the triangle it is
///   being selected. Otherwise, it is projected on the three lines defined
///   by the edges, and the nearest one is picked out. If it is inside the
///   edge it is being selected. Otherwise, the answer consists of the
///   shortest distance between the point and the three vertices. Note
///   that none of the three cases (plane, edge, vertex) is degenerate:
///   there exists an infinity of points in space for which the nearest
///   distance to the triangle is on a vertex.
template <typename T>
void distanceToTriangle(
    Array<T, 3> const &point, typename TriangleMesh<T>::Triangle const &triangle, T &distance,
    bool &pointIsBehind);

}  // namespace plb

#endif  // TRIANGLE_UTIL_H
