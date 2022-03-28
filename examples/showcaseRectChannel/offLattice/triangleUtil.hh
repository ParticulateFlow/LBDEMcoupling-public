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

#ifndef TRIANGLE_UTIL_HH
#define TRIANGLE_UTIL_HH

#include <limits>

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "core/util.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "offLattice/triangleUtil.h"

namespace plb {

template <typename T>
struct precision {
    static const T eps1;
};

template <typename T>
const T precision<T>::eps1 = (sizeof(T) == sizeof(float))
                                 ? std::numeric_limits<float>::epsilon()
                                 : (T)100.0 * std::numeric_limits<T>::epsilon();

template <typename T>
int pointOnTriangle(
    Array<T, 3> const &point1, Array<T, 3> const &point2, int flag,
    typename TriangleMesh<T>::Triangle const &triangle, Array<T, 3> &intersection,
    Array<T, 3> &normal, T &distance)
{
    PLB_ASSERT((flag == 0 || flag == 1 || flag == 2));

    Array<T, 3> v0 = triangle[0];
    Array<T, 3> v1 = triangle[1];
    Array<T, 3> v2 = triangle[2];

    Array<T, 3> e0 = v1 - v0;  // Triangle edge vectors starting at v0
    Array<T, 3> e1 = v2 - v0;

    crossProduct(e0, e1, normal);  // Triangle unit normal
    T normN = norm(normal);
    if (util::isZero(normN)) {
        return 0;  // The triangle has zero area.
    }
    normal /= normN;

    Array<T, 3> direction = point2 - point1;

    T num = dot(v0, normal) - dot(point1, normal);
    T denom = dot(direction, normal);

    T t = num / denom;

    // The function pointOnTriangle is essentially to verify crossings
    //   through the surface from inside to outside or vice versa. When
    //   one of the end-points of a segment is right on top of the surface
    //   this creates an awkward ambiguity. To remove this ambiguity, the
    //   triangle is slightly shifted in direction of its normal vector.
    if (((flag == 0 || flag == 1) && util::isZero(t, precision<T>::eps1))
        || ((flag == 0) && util::isOne(t, precision<T>::eps1)))
    {
        // Shift the vertex, and recompute everything that needs to be
        //   recomputed. v1 and v2 are not used any more, so there is no
        //   need to recompute them.
        v0 += (T)2 * precision<T>::eps1 * normal;
        num = dot(v0, normal) - dot(point1, normal);
        t = num / denom;
    }

    if (util::isZero(denom, precision<T>::eps1)) {
        if (util::isZero(num, precision<T>::eps1)) {
            return -1;  // Line belongs to the plane
        } else {
            return 0;  // Line does not intersect the plane
        }
    }

    if (flag == 0) {  // For intersection with a line segment
        if (util::lessThan(t, (T)0, precision<T>::eps1)
            || util::greaterThan(t, (T)1, precision<T>::eps1)) {
            return 0;
        }
    } else if (flag == 1) {  // For intersection with a half-line
        if (util::lessThan(t, (T)0, precision<T>::eps1)) {
            return 0;
        }
    }

    intersection = point1 + direction * t;  // Intersection point with the plane

    T a[2][2];
    a[0][0] = dot(e0, e0);
    a[0][1] = dot(e0, e1);
    a[1][0] = a[0][1];
    a[1][1] = dot(e1, e1);

    Array<T, 3> tmp((T)0.0, (T)0.0, (T)0.0);
    tmp[0] = intersection[0] - v0[0];
    tmp[1] = intersection[1] - v0[1];
    tmp[2] = intersection[2] - v0[2];

    T b[2];
    b[0] = dot(tmp, e0);
    b[1] = dot(tmp, e1);

    T det = a[0][0] * a[1][1] - a[0][1] * a[1][0];

    T u = (a[1][1] * b[0] - a[0][1] * b[1]) / det;
    T v = (a[0][0] * b[1] - a[1][0] * b[0]) / det;

    T upv = u + v;

    bool ueq0 = util::isZero(u, precision<T>::eps1);
    bool ueq1 = util::isOne(u, precision<T>::eps1);
    bool veq0 = util::isZero(v, precision<T>::eps1);
    bool veq1 = util::isOne(v, precision<T>::eps1);
    bool upveq1 = util::isOne(upv, precision<T>::eps1);

    int is_vertex = -1;
    int is_in_edge = -1;

    if ((u < (T)0.0 && !ueq0) || (u > (T)1.0 && !ueq1) || (v < (T)0.0 && !veq0)
        || (v > (T)1.0 && !veq1) || (upv > (T)1.0 && !upveq1))
        return 0;  // The point does not belong to the triangle
    else {         // The point belongs to the triangle (boundary or interior)
        distance = std::fabs(t) * norm(direction);
        // Distance between intersection and point1

        if (!ueq0 && !veq0 && !upveq1)
            return 1;  // Point is in the triangle interior
        else if (ueq0 && veq0)
            is_vertex = 0;
        else if (ueq1 && veq0)
            is_vertex = 1;
        else if (ueq0 && veq1)
            is_vertex = 2;
        else if (veq0 && !ueq0 && !ueq1)
            is_in_edge = 0;
        else if (ueq0 && !veq0 && !veq1)
            is_in_edge = 2;
        else if (upveq1 && !ueq1 && !veq1)
            is_in_edge = 1;
    }

    if (is_in_edge != -1) {  // The point belongs to an edge
        normal = triangle.edgeNormal(is_in_edge);
    } else if (is_vertex != -1) {  // The point is a vertex
        normal = triangle.vertex(is_vertex)->normal();
    }

    return 1;
}

template <typename T>
bool segmentIntersectsTriangle(
    Array<T, 3> const &p1, Array<T, 3> const &p2,
    typename TriangleMesh<T>::Triangle const &triangle)
{
    Array<T, 3> v0 = triangle[0];
    Array<T, 3> v1 = triangle[1];
    Array<T, 3> v2 = triangle[2];

    // Triangle edge vectors starting at v0
    Array<T, 3> e0(v1 - v0);
    Array<T, 3> e1(v2 - v0);

    Array<T, 3> normal(crossProduct(e0, e1));
    Array<T, 3> direction(p2 - p1);
    T denom = dot(direction, normal);

    if (util::isZero(denom, precision<T>::eps1))
        return false;  // The segment belongs to the plane or it is parallel to it
                       //   and does not intersect it.

    T num = dot((v0 - p1), normal);
    T t = num / denom;

    // The function segmentIntersectsTriangle is essentially to verify crossings
    //   through the surface from inside to outside or vice versa. When
    //   one of the end-points of a segment is right on top of the surface
    //   this creates an awkward ambiguity. To remove this ambiguity, the
    //   triangle is slightly shifted in direction of its normal vector.
    if ((util::isZero(t, precision<T>::eps1)) || (util::isOne(t, precision<T>::eps1))) {
        normal /= norm(normal);
        // Shift the vertex, and recompute everything that needs to be
        //   recomputed. v1 and v2 are not used any more, so there is no
        //   need to recompute them.
        T tmp = (T)2 * precision<T>::eps1;
        v0 += tmp * normal;
        num = dot((v0 - p1), normal);
        t = num / denom;
    }

    if (util::lessThan(t, (T)0, precision<T>::eps1)
        || util::greaterThan(t, (T)1, precision<T>::eps1))
        return false;

    Array<T, 3> intersection(p1 + direction * t);

    T a[2][2];
    a[0][0] = e0[0] * e0[0] + e0[1] * e0[1] + e0[2] * e0[2];
    a[0][1] = e0[0] * e1[0] + e0[1] * e1[1] + e0[2] * e1[2];
    a[1][0] = a[0][1];
    a[1][1] = e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2];

    Array<T, 3> tmp(intersection - v0);

    T b[2];
    b[0] = tmp[0] * e0[0] + tmp[1] * e0[1] + tmp[2] * e0[2];
    b[1] = tmp[0] * e1[0] + tmp[1] * e1[1] + tmp[2] * e1[2];

    T det = a[0][0] * a[1][1] - a[0][1] * a[1][0];

    T u = (a[1][1] * b[0] - a[0][1] * b[1]) / det;
    T v = (a[0][0] * b[1] - a[1][0] * b[0]) / det;

    T upv = u + v;

    if (util::lessThan(u, (T)0, precision<T>::eps1)
        || util::greaterThan(u, (T)1, precision<T>::eps1)
        || util::lessThan(v, (T)0, precision<T>::eps1)
        || util::greaterThan(v, (T)1, precision<T>::eps1)
        || util::greaterThan(upv, (T)1, precision<T>::eps1))
        return false;  // The point does not belong to the triangle

    return true;
}

template <typename T>
void distanceToEdgeLine(
    Array<T, 3> const &point, typename TriangleMesh<T>::Triangle const &triangle, plint whichEdge,
    T &distance, bool &intersectionIsInside)
{
    PLB_ASSERT((whichEdge == 0 || whichEdge == 1 || whichEdge == 2));

    plint iVertex1 = whichEdge;
    plint iVertex2 = (whichEdge + 1) % 3;
    Array<T, 3> vertex1(triangle[iVertex1]);
    Array<T, 3> vertex2(triangle[iVertex2]);
    Array<T, 3> pv(point - vertex1);
    Array<T, 3> e(vertex2 - vertex1);

    T u = dot(pv, e) / dot(e, e);
    Array<T, 3> x = vertex1 + u * e;
    distance = norm(point - x);

    intersectionIsInside = util::greaterEqual(u, (T)0, precision<T>::eps1)
                           && util::lessEqual(u, (T)1, precision<T>::eps1);
}

template <typename T>
void distanceToTrianglePlane(
    Array<T, 3> const &point, typename TriangleMesh<T>::Triangle const &triangle, T &distance,
    bool &intersectionIsInside, bool &pointIsBehind)
{
    Array<T, 3> normal = triangle.normal();
    Array<T, 3> intersection((T)0, (T)0, (T)0);
    int flag = 2;  // Line.
    intersectionIsInside =
        pointOnTriangle(point, point + normal, flag, triangle, intersection, normal, distance) == 1;
    T projection = dot(normal, point - intersection);

    pointIsBehind = util::lessEqual(projection, (T)0, precision<T>::eps1);
}

template <typename T>
void distanceToTriangle(
    Array<T, 3> const &point, typename TriangleMesh<T>::Triangle const &triangle, T &distance,
    bool &pointIsBehind)
{
    bool intersectionIsInside;
    // First possibility: The projection of the point is inside the triangle.
    distanceToTrianglePlane(point, triangle, distance, intersectionIsInside, pointIsBehind);
    if (intersectionIsInside) {
        return;
    }
    // Second possibility: The projection of the point is inside an edge.
    bool intersectsWithEdge = false;
    for (plint iEdge = 0; iEdge <= 2; ++iEdge) {
        T newDistance;
        distanceToEdgeLine(point, triangle, iEdge, newDistance, intersectionIsInside);
        if (iEdge == 0 || newDistance < distance) {
            distance = newDistance;
            // The edge-line to which the point is closest should be selected if
            //   and only if the intersection is inside the edge. Otherwise
            //   the closest distance to the triangle is on one of the vertices
            //   which is selected further down.
            intersectsWithEdge = intersectionIsInside;
        }
    }
    if (intersectsWithEdge) {
        return;
    }
    // Default: Compute the closest distance to one of the vertices.
    T d0Sqr = normSqr(point - triangle[0]);
    T d1Sqr = normSqr(point - triangle[1]);
    T d2Sqr = normSqr(point - triangle[2]);
    T minDistSqr = std::min(d0Sqr, std::min(d1Sqr, d2Sqr));
    distance = std::sqrt(minDistSqr);
}

}  // namespace plb

#endif  // TRIANGLE_UTIL_HH
