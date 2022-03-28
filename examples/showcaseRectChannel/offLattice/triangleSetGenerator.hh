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

#ifndef TRIANGLE_SET_GENERATOR_HH
#define TRIANGLE_SET_GENERATOR_HH

#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "offLattice/triangleBoundary3D.h"
#include "offLattice/triangleSetGenerator.h"

namespace plb {

template <typename T>
TriangleSet<T> constructEllipsoid(
    Array<T, 3> const &center, T a, T b, T c, plint minNumOfTriangles, T eps)
{
    TriangleSet<T> ellipsoid(constructSphere<T>(center, (T)1, minNumOfTriangles, eps));
    ellipsoid.scale(a, b, c);
    return ellipsoid;
}

template <typename T>
TriangleSet<T> constructSphere(Array<T, 3> const &center, T radius, plint minNumOfTriangles, T eps)
{
    PLB_ASSERT(util::greaterThan(radius, (T)0) && minNumOfTriangles >= 8);

    typedef typename TriangleSet<T>::Triangle Triangle;

    std::vector<Triangle> triangles;

    // Create a triangularized unit sphere

    // Initial 6 vertices

    Array<T, 3> va;
    va[0] = (T)1;
    va[1] = (T)0;
    va[2] = (T)0;

    Array<T, 3> vb;
    vb[0] = (T)0;
    vb[1] = (T)1;
    vb[2] = (T)0;

    Array<T, 3> vc;
    vc[0] = (T)-1;
    vc[1] = (T)0;
    vc[2] = (T)0;

    Array<T, 3> vd;
    vd[0] = (T)0;
    vd[1] = (T)-1;
    vd[2] = (T)0;

    Array<T, 3> ve;
    ve[0] = (T)0;
    ve[1] = (T)0;
    ve[2] = (T)1;

    Array<T, 3> vf;
    vf[0] = (T)0;
    vf[1] = (T)0;
    vf[2] = (T)-1;

    // Initial 8 triangles

    Triangle triangle;

    triangle[0] = ve;
    triangle[1] = va;
    triangle[2] = vb;
    triangles.push_back(triangle);

    triangle[0] = ve;
    triangle[1] = vb;
    triangle[2] = vc;
    triangles.push_back(triangle);

    triangle[0] = ve;
    triangle[1] = vc;
    triangle[2] = vd;
    triangles.push_back(triangle);

    triangle[0] = ve;
    triangle[1] = vd;
    triangle[2] = va;
    triangles.push_back(triangle);

    triangle[0] = vf;
    triangle[1] = vb;
    triangle[2] = va;
    triangles.push_back(triangle);

    triangle[0] = vf;
    triangle[1] = vc;
    triangle[2] = vb;
    triangles.push_back(triangle);

    triangle[0] = vf;
    triangle[1] = vd;
    triangle[2] = vc;
    triangles.push_back(triangle);

    triangle[0] = vf;
    triangle[1] = va;
    triangle[2] = vd;
    triangles.push_back(triangle);

    // Perform refinement iterations

    plint size = 0;
    while ((size = triangles.size()) < minNumOfTriangles) {
        for (plint i = 0; i < size; i++) {
            va = triangles[i][0];
            vb = triangles[i][1];
            vc = triangles[i][2];

            vd = (T)0.5 * (va + vb);
            ve = (T)0.5 * (vb + vc);
            vf = (T)0.5 * (vc + va);

            vd /= norm(vd);
            ve /= norm(ve);
            vf /= norm(vf);

            triangles[i][0] = vd;
            triangles[i][1] = ve;
            triangles[i][2] = vf;

            triangle[0] = va;
            triangle[1] = vd;
            triangle[2] = vf;
            triangles.push_back(triangle);

            triangle[0] = vd;
            triangle[1] = vb;
            triangle[2] = ve;
            triangles.push_back(triangle);

            triangle[0] = vf;
            triangle[1] = ve;
            triangle[2] = vc;
            triangles.push_back(triangle);
        }
    }

    // Scale and translate the mesh

    TriangleSet<T> triangleSet(triangles, eps);

    triangleSet.scale(radius);
    triangleSet.translate(center);

    return triangleSet;
}

template <typename T>
TriangleSet<T> constructHalfSphere(
    Array<T, 3> const &center, T radius, plint minNumOfTriangles, T eps)
{
    PLB_ASSERT(util::greaterThan(radius, (T)0) && minNumOfTriangles >= 4);

    typedef typename TriangleSet<T>::Triangle Triangle;

    std::vector<Triangle> triangles;

    // Create a triangularized unit half-sphere

    // Initial 5 vertices

    Array<T, 3> va;
    va[0] = (T)1;
    va[1] = (T)0;
    va[2] = (T)0;

    Array<T, 3> vb;
    vb[0] = (T)0;
    vb[1] = (T)1;
    vb[2] = (T)0;

    Array<T, 3> vc;
    vc[0] = (T)-1;
    vc[1] = (T)0;
    vc[2] = (T)0;

    Array<T, 3> vd;
    vd[0] = (T)0;
    vd[1] = (T)-1;
    vd[2] = (T)0;

    Array<T, 3> ve;
    ve[0] = (T)0;
    ve[1] = (T)0;
    ve[2] = (T)1;

    // Initial 4 triangles

    Triangle triangle;

    triangle[0] = ve;
    triangle[1] = va;
    triangle[2] = vb;
    triangles.push_back(triangle);

    triangle[0] = ve;
    triangle[1] = vb;
    triangle[2] = vc;
    triangles.push_back(triangle);

    triangle[0] = ve;
    triangle[1] = vc;
    triangle[2] = vd;
    triangles.push_back(triangle);

    triangle[0] = ve;
    triangle[1] = vd;
    triangle[2] = va;
    triangles.push_back(triangle);

    // Perform refinement iterations

    Array<T, 3> vf;
    plint size = 0;
    while ((size = triangles.size()) < minNumOfTriangles) {
        for (plint i = 0; i < size; i++) {
            va = triangles[i][0];
            vb = triangles[i][1];
            vc = triangles[i][2];

            vd = (T)0.5 * (va + vb);
            ve = (T)0.5 * (vb + vc);
            vf = (T)0.5 * (vc + va);

            vd /= norm(vd);
            ve /= norm(ve);
            vf /= norm(vf);

            triangles[i][0] = vd;
            triangles[i][1] = ve;
            triangles[i][2] = vf;

            triangle[0] = va;
            triangle[1] = vd;
            triangle[2] = vf;
            triangles.push_back(triangle);

            triangle[0] = vd;
            triangle[1] = vb;
            triangle[2] = ve;
            triangles.push_back(triangle);

            triangle[0] = vf;
            triangle[1] = ve;
            triangle[2] = vc;
            triangles.push_back(triangle);
        }
    }

    // Scale and translate the mesh

    TriangleSet<T> triangleSet(triangles, eps);

    triangleSet.scale(radius);
    triangleSet.translate(center);

    return triangleSet;
}

template <typename T>
static void sphereIsomorphToCenter(Array<T, 3> &p, Array<T, 3> const &center, T radius)
{
    static T pi = std::acos((T)-1);
    Array<T, 3> v(p - center);
    T a = norm(v);
    if (a > getEpsilon<T>((T)100)) {
        v /= a;
        p = center + v * radius * ((T)1 - (T)2 / pi * std::acos(a / radius));
    }
}

template <typename T>
TriangleSet<T> constructDisk(Array<T, 3> const &center, T radius, plint minNumOfTriangles, T eps)
{
    typedef typename TriangleSet<T>::Triangle Triangle;
    TriangleSet<T> halfSphere(constructHalfSphere(center, radius, minNumOfTriangles, eps));
    std::vector<Triangle> const &triangles = halfSphere.getTriangles();
    std::vector<Triangle> diskTriangles;
    for (pluint i = 0; i < triangles.size(); ++i) {
        Triangle triangle = triangles[i];
        triangle[0][2] = center[2];
        triangle[1][2] = center[2];
        triangle[2][2] = center[2];
        sphereIsomorphToCenter(triangle[0], center, radius);
        sphereIsomorphToCenter(triangle[1], center, radius);
        sphereIsomorphToCenter(triangle[2], center, radius);
        diskTriangles.push_back(triangle);
    }
    return TriangleSet<T>(diskTriangles, eps);
}

template <typename T>
TriangleSet<T> constructCylinder(
    Array<T, 3> const &inletCenter, T inletRadius, T outletRadius, T length, plint nAxial,
    plint nCirc, T eps)
{
    PLB_ASSERT(
        util::greaterThan(inletRadius, (T)0) && util::greaterThan(outletRadius, (T)0)
        && util::greaterThan(length, (T)0) && nAxial >= 2 && nCirc >= 2);

    // This is needed to be backwards compatible.

    nAxial = 2 * nAxial - 1;
    nCirc *= 2;

    // Construction of the cylindrical grid and triangulation

    typedef typename TriangleSet<T>::Triangle Triangle;

    std::vector<Triangle> triangles;

    static T pi = std::acos((T)-1);

    T dtheta = (T)2 * pi / nCirc;
    T dx = length / (T)(nAxial - 1);

    T x0 = inletCenter[0];
    T y0 = inletCenter[1];
    T z0 = inletCenter[2];

    std::vector<Array<T, 3> > circle0(nCirc);
    for (plint i = 0; i < nCirc; i++) {
        T r = inletRadius;
        T theta = (T)i * dtheta;
        circle0[i] = Array<T, 3>(x0, y0 + r * std::cos(theta), z0 + r * std::sin(theta));
    }

    std::vector<Array<T, 3> > circle1(nCirc);
    for (plint i = 0; i < nAxial - 1; i++) {
        T x = x0 + (T)(i + 1) * dx;
        T r = (outletRadius - inletRadius) * (x - x0) / length + inletRadius;
        for (plint i = 0; i < nCirc; i++) {
            T theta = (T)i * dtheta;
            circle1[i] = Array<T, 3>(x, y0 + r * std::cos(theta), z0 + r * std::sin(theta));
        }

        for (plint i = 0; i < nCirc; i++) {
            Array<T, 3> va = circle0[i];
            Array<T, 3> vb = circle1[i];
            Array<T, 3> vc = i != nCirc - 1 ? circle1[i + 1] : circle1[0];
            Array<T, 3> vd = i != nCirc - 1 ? circle0[i + 1] : circle0[0];

            Triangle triangle;

            triangle[0] = va;
            triangle[1] = vd;
            triangle[2] = vb;
            triangles.push_back(triangle);

            triangle[0] = vb;
            triangle[1] = vd;
            triangle[2] = vc;
            triangles.push_back(triangle);
        }

        std::swap(circle0, circle1);
    }

    return TriangleSet<T>(triangles, eps);
}

template <typename T>
TriangleSet<T> constructCylinder(
    Array<T, 3> const &inletCenter, T inletRadius, T outletRadius, T length, plint nAxial,
    plint nCirc, std::vector<Array<T, 3> > &inletPoints, T eps)
{
    PLB_ASSERT(
        util::greaterThan(inletRadius, (T)0) && util::greaterThan(outletRadius, (T)0)
        && util::greaterThan(length, (T)0) && nAxial >= 2 && nCirc >= 2);

    // This is needed to be backwards compatible.

    nAxial = 2 * nAxial - 1;
    nCirc *= 2;

    // Construction of the cylindrical grid and triangulation

    typedef typename TriangleSet<T>::Triangle Triangle;

    std::vector<Triangle> triangles;
    inletPoints.clear();

    static T pi = std::acos((T)-1);

    T dtheta = (T)2 * pi / nCirc;
    T dx = length / (T)(nAxial - 1);

    T x0 = inletCenter[0];
    T y0 = inletCenter[1];
    T z0 = inletCenter[2];

    std::vector<Array<T, 3> > circle0(nCirc);
    for (plint i = 0; i < nCirc; i++) {
        T r = inletRadius;
        T theta = (T)i * dtheta;
        circle0[i] = Array<T, 3>(x0, y0 + r * std::cos(theta), z0 + r * std::sin(theta));
    }
    inletPoints = circle0;

    std::vector<Array<T, 3> > circle1(nCirc);
    for (plint i = 0; i < nAxial - 1; i++) {
        T x = x0 + (T)(i + 1) * dx;
        T r = (outletRadius - inletRadius) * (x - x0) / length + inletRadius;
        for (plint i = 0; i < nCirc; i++) {
            T theta = (T)i * dtheta;
            circle1[i] = Array<T, 3>(x, y0 + r * std::cos(theta), z0 + r * std::sin(theta));
        }

        for (plint i = 0; i < nCirc; i++) {
            Array<T, 3> va = circle0[i];
            Array<T, 3> vb = circle1[i];
            Array<T, 3> vc = i != nCirc - 1 ? circle1[i + 1] : circle1[0];
            Array<T, 3> vd = i != nCirc - 1 ? circle0[i + 1] : circle0[0];

            Triangle triangle;

            triangle[0] = va;
            triangle[1] = vd;
            triangle[2] = vb;
            triangles.push_back(triangle);

            triangle[0] = vb;
            triangle[1] = vd;
            triangle[2] = vc;
            triangles.push_back(triangle);
        }

        std::swap(circle0, circle1);
    }

    return TriangleSet<T>(triangles, eps);
}

template <typename T>
TriangleSet<T> constructCylinder(
    Array<T, 3> const &inletCenter, Array<T, 3> const &axis, T inletRadius, T outletRadius,
    T length, plint nAxial, plint nCirc, std::vector<Array<T, 3> > &inletPoints, T eps)
{
    TriangleSet<T> cylinder(constructCylinder(
        Array<T, 3>(T(), T(), T()), inletRadius, outletRadius, length, nAxial, nCirc, inletPoints,
        eps));
    Array<T, 3> xAxis((T)1, T(), T());
    Array<T, 3> rotAxis;
    T alignment = dot(xAxis, axis);
    if (!util::isOne(alignment)) {
        crossProduct<T>(xAxis, axis, rotAxis);
        rotAxis /= norm(rotAxis);
        T angle = angleBetweenVectors(xAxis, axis);
        cylinder.rotateAtOrigin(rotAxis, angle);
        for (pluint i = 0; i < inletPoints.size(); ++i) {
            inletPoints[i] = rotateAtOrigin(inletPoints[i], rotAxis, angle);
        }
    }
    cylinder.translate(inletCenter);
    for (pluint i = 0; i < inletPoints.size(); ++i) {
        inletPoints[i] += inletCenter;
    }
    return cylinder;
}

template <typename T>
static void addSurface(
    Array<T, 3> const &lowerCorner, Array<T, 3> const &delta1, plint n1, Array<T, 3> const &delta2,
    plint n2, std::vector<typename TriangleSet<T>::Triangle> &triangles)
{
    Array<T, 3> pos1(lowerCorner);
    for (plint i1 = 0; i1 < n1; ++i1, pos1 += delta1) {
        Array<T, 3> pos2(pos1);
        for (plint i2 = 0; i2 < n2; ++i2, pos2 += delta2) {
            typename TriangleSet<T>::Triangle triangle;
            triangle[0] = pos2;
            triangle[1] = pos2 + delta1;
            triangle[2] = pos2 + delta2;
            triangles.push_back(triangle);
            triangle[0] += delta1 + delta2;
            std::swap(triangle[1], triangle[2]);
            triangles.push_back(triangle);
        }
    }
}

template <typename T>
TriangleSet<T> constructCuboidWithInletOutlet(
    Array<T, 3> const &lowerCorner, Array<T, 3> const &upperCorner,
    Array<plint, 3> const &nSegments, plint inletOutletDirection, T eps)
{
    PLB_ASSERT(inletOutletDirection == 0 || inletOutletDirection == 1 || inletOutletDirection == 2);
    std::vector<typename TriangleSet<T>::Triangle> triangles;
    T lx = upperCorner[0] - lowerCorner[0];
    T ly = upperCorner[1] - lowerCorner[1];
    T lz = upperCorner[2] - lowerCorner[2];
    Array<T, 3> deltaX(lx / (T)nSegments[0], T(), T());
    Array<T, 3> deltaY(T(), ly / (T)nSegments[1], T());
    Array<T, 3> deltaZ(T(), T(), lz / (T)nSegments[2]);

    if (inletOutletDirection != 0) {
        addSurface(lowerCorner, deltaZ, nSegments[2], deltaY, nSegments[1], triangles);
        addSurface(
            lowerCorner + Array<T, 3>(lx, T(), T()), deltaY, nSegments[1], deltaZ, nSegments[2],
            triangles);
    }

    if (inletOutletDirection != 1) {
        addSurface(lowerCorner, deltaX, nSegments[0], deltaZ, nSegments[2], triangles);
        addSurface(
            lowerCorner + Array<T, 3>(T(), ly, T()), deltaZ, nSegments[2], deltaX, nSegments[0],
            triangles);
    }

    if (inletOutletDirection != 2) {
        addSurface(lowerCorner, deltaY, nSegments[1], deltaX, nSegments[0], triangles);
        addSurface(
            lowerCorner + Array<T, 3>(T(), T(), lz), deltaX, nSegments[0], deltaY, nSegments[1],
            triangles);
    }

    return TriangleSet<T>(triangles, eps);
}

template <typename T>
TriangleSet<T> constructCuboidWithInletOutlet(
    Cuboid<T> const &cuboid, Array<plint, 3> const &nSegments, plint inletOutletDirection, T eps)
{
    return constructCuboidWithInletOutlet<T>(
        cuboid.lowerLeftCorner, cuboid.upperRightCorner, nSegments, inletOutletDirection, eps);
}

template <typename T>
TriangleSet<T> constructCuboid(
    Array<T, 3> const &lowerCorner, Array<T, 3> const &upperCorner,
    Array<plint, 3> const &nSegments, T eps)
{
    std::vector<typename TriangleSet<T>::Triangle> triangles;
    T lx = upperCorner[0] - lowerCorner[0];
    T ly = upperCorner[1] - lowerCorner[1];
    T lz = upperCorner[2] - lowerCorner[2];
    Array<T, 3> deltaX(lx / (T)nSegments[0], T(), T());
    Array<T, 3> deltaY(T(), ly / (T)nSegments[1], T());
    Array<T, 3> deltaZ(T(), T(), lz / (T)nSegments[2]);

    addSurface(lowerCorner, deltaZ, nSegments[2], deltaY, nSegments[1], triangles);
    addSurface(
        lowerCorner + Array<T, 3>(lx, T(), T()), deltaY, nSegments[1], deltaZ, nSegments[2],
        triangles);

    addSurface(lowerCorner, deltaX, nSegments[0], deltaZ, nSegments[2], triangles);
    addSurface(
        lowerCorner + Array<T, 3>(T(), ly, T()), deltaZ, nSegments[2], deltaX, nSegments[0],
        triangles);

    addSurface(lowerCorner, deltaY, nSegments[1], deltaX, nSegments[0], triangles);
    addSurface(
        lowerCorner + Array<T, 3>(T(), T(), lz), deltaX, nSegments[0], deltaY, nSegments[1],
        triangles);

    return TriangleSet<T>(triangles, eps);
}

template <typename T>
TriangleSet<T> constructCuboid(Cuboid<T> const &cuboid, Array<plint, 3> const &nSegments, T eps)
{
    return constructCuboid<T>(cuboid.lowerLeftCorner, cuboid.upperRightCorner, nSegments, eps);
}

template <typename T>
TriangleSet<T> patchTubes(
    TriangleSet<T> const &geometryWithOpenings, plint sortDirection,
    std::vector<T> const &patchLengths)
{
    typedef typename TriangleSet<T>::Triangle Triangle;

    std::vector<Triangle> fullGeometry(geometryWithOpenings.getTriangles());

    DEFscaledMesh<T> *defMesh = new DEFscaledMesh<T>(geometryWithOpenings);
    TriangularSurfaceMesh<T> &mesh = defMesh->getMesh();

    std::vector<Lid> holes = mesh.closeHoles();
    std::sort(holes.begin(), holes.end(), LidLessThan<T>(sortDirection, mesh));

    PLB_ASSERT(holes.size() == patchLengths.size());

    for (pluint iHole = 0; iHole < holes.size(); ++iHole) {
        Array<T, 3> baryCenter = computeGeometricCenter(mesh, holes[iHole]);
        plint numHoleVertices = (plint)holes[iHole].boundaryVertices.size();

        Array<T, 3> normal = computeNormal(mesh, holes[iHole]);
        T aveRadius = computeGeometricRadius(mesh, holes[iHole]);

        Array<T, 3> nextCenter = baryCenter + normal * aveRadius * (T)10 / (T)numHoleVertices;

        plint numInletPoints = numHoleVertices;
        bool oddNumber = numHoleVertices % 2 == 1;
        if (oddNumber)
            numInletPoints--;  // Must be even for cylinder construction algorithm.
        std::vector<Array<T, 3> > inletPoints;
        plint numPointsOnLength = numInletPoints * patchLengths[iHole] / aveRadius / 8;
        if (numPointsOnLength < 3)
            numPointsOnLength = 3;
        TriangleSet<T> piece = constructCylinder(
            nextCenter, normal, aveRadius, aveRadius, patchLengths[iHole], numPointsOnLength,
            numInletPoints / 2, inletPoints, geometryWithOpenings.getEpsilon());
        std::vector<Triangle> pieceTriangles = piece.getTriangles();

        plint newId = 0;
        T minDistance = std::numeric_limits<T>::max();
        plint minDistanceId = -1;
        for (plint i = 0; i < numHoleVertices; ++i) {
            plint iVertex = holes[iHole].boundaryVertices[i];
            Array<T, 3> p1 = mesh.getVertex(iVertex);
            T nextDistance = norm(inletPoints[newId] - p1);
            if (nextDistance < minDistance) {
                minDistance = nextDistance;
                minDistanceId = i;
            }
        }
        plint oldId = minDistanceId;
        plint newId_p1 = 0;
        for (plint i = 0; i < numInletPoints; ++i) {
            plint newId_p1 = (newId + 1) % numInletPoints;
            plint oldId_p1 = oldId - 1;
            if (oldId_p1 < 0)
                oldId_p1 = numHoleVertices - 1;

            plint oldVertex1 = holes[iHole].boundaryVertices[oldId];
            plint oldVertex2 = holes[iHole].boundaryVertices[oldId_p1];
            Array<T, 3> p1 = mesh.getVertex(oldVertex1);
            Array<T, 3> p2 = inletPoints[newId];
            Array<T, 3> p3 = mesh.getVertex(oldVertex2);
            Array<T, 3> p4 = inletPoints[newId_p1];

            pieceTriangles.push_back(Triangle(p1, p3, p2));
            pieceTriangles.push_back(Triangle(p2, p3, p4));

            std::swap(newId, newId_p1);
            std::swap(oldId, oldId_p1);
        }

        if (oddNumber) {
            plint id_a = holes[iHole].boundaryVertices[oldId];
            plint oldId_p2 = oldId - 1;
            if (oldId_p2 < 0)
                oldId_p2 = numHoleVertices - 1;
            plint id_b = holes[iHole].boundaryVertices[oldId_p2];
            plint id_c = newId_p1;
            Array<T, 3> a = mesh.getVertex(id_a);
            Array<T, 3> b = mesh.getVertex(id_b);
            Array<T, 3> c = inletPoints[id_c];
            pieceTriangles.push_back(Triangle(a, b, c));
        }

        fullGeometry.insert(fullGeometry.end(), pieceTriangles.begin(), pieceTriangles.end());
    }

    return TriangleSet<T>(fullGeometry, geometryWithOpenings.getEpsilon());
}

template <typename T>
TriangleSet<T> constructRectangle(T lx, T ly, plint nx, plint ny, T eps)
{
    PLB_ASSERT(util::greaterThan(lx, (T)0) && util::greaterThan(ly, (T)0) && nx >= 2 && ny >= 2);

    T dx = lx / (T)(nx - 1);
    T dy = ly / (T)(ny - 1);

    typedef typename TriangleSet<T>::Triangle Triangle;

    std::vector<Triangle> triangles;
    T z = (T)0;
    for (plint i = 0; i < nx - 1; i++) {
        T x0 = i * dx;
        T x1 = x0 + dx;
        for (plint j = 0; j < ny - 1; j++) {
            T y0 = j * dy;
            T y1 = y0 + dy;

            Array<T, 3> v0, v1, v2;
            v0 = Array<T, 3>(x0, y0, z);
            v1 = Array<T, 3>(x1, y0, z);
            v2 = Array<T, 3>(x1, y1, z);

            Triangle triangle;
            triangle[0] = v0;
            triangle[1] = v1;
            triangle[2] = v2;

            triangles.push_back(triangle);

            v1 = Array<T, 3>(x0, y1, z);

            triangle[1] = v2;
            triangle[2] = v1;

            triangles.push_back(triangle);
        }
    }

    return TriangleSet<T>(triangles, eps);
}

template <typename T>
TriangleSet<T> constructGenericRectangle(
    T lx, T ly, plint nx, plint ny, Array<T, 3> const &center, Array<T, 3> const &normal, T eps)
{
    PLB_ASSERT(util::greaterThan(lx, (T)0) && util::greaterThan(ly, (T)0) && nx >= 2 && ny >= 2);

    static T pi = std::acos((T)-1);

    TriangleSet<T> rectangle = constructRectangle<T>(lx, ly, nx, ny, eps);
    Array<T, 3> oldCenter((T)0.5 * lx, (T)0.5 * ly, (T)0.0);
    rectangle.translate(-oldCenter);

    T r = norm(normal);
    PLB_ASSERT(!util::isZero(r));
    Array<T, 3> unitNormal = normal / r;
    T theta = std::acos(unitNormal[2]);
    Array<T, 3> zAxis((T)0.0, (T)0.0, (T)1.0);
    Array<T, 3> normedAxis = crossProduct<T>(zAxis, unitNormal);
    T l = norm(normedAxis);
    if (util::isZero(l)) {
        normedAxis = Array<T, 3>((T)1.0, (T)0.0, (T)0.0);
        if (dot<T>(zAxis, unitNormal) > 0.0) {
            theta = 0.0;
        } else {
            theta = pi;
        }
    } else {
        normedAxis /= l;
    }
    rectangle.rotateAtOrigin(normedAxis, theta);
    rectangle.translate(center);

    return rectangle;
}

template <typename T>
TriangleSet<T> constructStrip(
    std::vector<Array<T, 3> > const &from, std::vector<Array<T, 3> > const &to, T eps)
{
    PLB_ASSERT(from.size() >= 2);
    PLB_ASSERT(from.size() == to.size());
    plint size = from.size();

    typedef typename TriangleSet<T>::Triangle Triangle;

    std::vector<Triangle> triangles;

    plint count = 0;
    for (plint i = 0; i < size - 1; i++) {
        Array<T, 3> from0 = from[i];
        Array<T, 3> from1 = from[i + 1];

        Array<T, 3> to0 = to[i];
        Array<T, 3> to1 = to[i + 1];

        Triangle triangle;
        if (count % 2 == 0) {
            triangle[0] = from0;
            triangle[1] = from1;
            triangle[2] = to0;
            triangles.push_back(triangle);

            triangle[0] = from1;
            triangle[1] = to1;
            triangle[2] = to0;
            triangles.push_back(triangle);
        } else {
            triangle[0] = from0;
            triangle[1] = from1;
            triangle[2] = to1;
            triangles.push_back(triangle);

            triangle[0] = from0;
            triangle[1] = to1;
            triangle[2] = to0;
            triangles.push_back(triangle);
        }
        count++;
    }

    return TriangleSet<T>(triangles, eps);
}

template <typename T>
TriangleSet<T> constructAnnulus(
    Array<T, 3> const &center, T innerRadius, T outerRadius, plint nCirc, T eps)
{
    PLB_ASSERT(
        util::greaterEqual(innerRadius, (T)0) && util::greaterThan(outerRadius, (T)0)
        && !util::fpequal(innerRadius, outerRadius) && nCirc >= 3);

    static T pi = std::acos((T)-1);

    T dtheta = (T)2 * pi / nCirc;

    T x0 = center[0];
    T y0 = center[1];
    T z0 = center[2];

    if (!util::isZero(innerRadius)) {
        std::vector<Array<T, 3> > from(nCirc + 1), to(nCirc + 1);
        for (plint i = 0; i < nCirc; i++) {
            T theta = -i * dtheta;  // The minus is to have +x orientation.

            T cosTheta = std::cos(theta);
            T sinTheta = std::sin(theta);

            from[i][0] = x0;
            from[i][1] = y0 + innerRadius * cosTheta;
            from[i][2] = z0 + innerRadius * sinTheta;

            to[i][0] = x0;
            to[i][1] = y0 + outerRadius * cosTheta;
            to[i][2] = z0 + outerRadius * sinTheta;
        }
        from[nCirc] = from[0];
        to[nCirc] = to[0];

        return constructStrip<T>(from, to, eps);
    } else {
        std::vector<Array<T, 3> > to(nCirc + 1);
        for (plint i = 0; i < nCirc; i++) {
            T theta = i * dtheta;

            T cosTheta = std::cos(theta);
            T sinTheta = std::sin(theta);

            to[i][0] = x0;
            to[i][1] = y0 + outerRadius * cosTheta;
            to[i][2] = z0 + outerRadius * sinTheta;
        }
        to[nCirc] = to[0];

        typedef typename TriangleSet<T>::Triangle Triangle;
        std::vector<Triangle> triangles;
        for (plint i = 0; i < nCirc; i++) {
            Triangle triangle;
            triangle[0] = center;
            triangle[1] = to[i];
            triangle[2] = to[i + 1];
            triangles.push_back(triangle);
        }

        return TriangleSet<T>(triangles, eps);
    }
}

template <typename T>
TriangleSet<T> sparseBlockStructureToTriangleSet(
    MultiBlock3D const &block, bool highResolution, T eps)
{
    std::map<plint, Box3D> const &bulks = block.getSparseBlockStructure().getBulks();
    Array<plint, 3> nSegments(1, 1, 1);

    TriangleSet<T> bulkSet(eps);
    std::map<plint, Box3D>::const_iterator it;
    for (it = bulks.begin(); it != bulks.end(); ++it) {
        Box3D const &bulk = it->second;
        Array<T, 3> lowerCorner((T)bulk.x0, (T)bulk.y0, (T)bulk.z0);
        Array<T, 3> upperCorner((T)bulk.x1, (T)bulk.y1, (T)bulk.z1);
        if (highResolution) {
            nSegments[0] = bulk.x1 - bulk.x0;
            nSegments[1] = bulk.y1 - bulk.y0;
            nSegments[2] = bulk.z1 - bulk.z0;
        }
        TriangleSet<T> set = constructCuboid(lowerCorner, upperCorner, nSegments, eps);
        bulkSet.append(set);
    }

    return bulkSet;
}

}  // namespace plb

#endif  // TRIANGLE_SET_GENERATOR_HH
