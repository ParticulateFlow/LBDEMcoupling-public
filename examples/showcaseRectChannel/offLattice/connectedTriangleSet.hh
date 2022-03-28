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

#ifndef CONNECTED_TRIANGLE_SET_HH
#define CONNECTED_TRIANGLE_SET_HH

#include <cmath>
#include <cstdio>
#include <limits>

#include "core/array.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "offLattice/connectedTriangleSet.h"
#include "offLattice/triangleSet.h"

namespace plb {

template <typename T>
ConnectedTriangleSet<T>::ConnectedTriangleSet(TriangleSet<T> const &triangleSet)
{
    std::vector<Array<Array<T, 3>, 3> > oldTriangles = triangleSet.getTriangles();

    numTriangles = (plint)oldTriangles.size();
    if (numTriangles == 0) {
        return;
    }

    triangles.resize(numTriangles);

    T epsilon = triangleSet.getEpsilon();
    PositionLessThan3D<T, VertexSetNode> lessThan(epsilon);
    VertexSet vertexSet(lessThan);

    // Unify duplicated vertices.
    numVertices = 0;
    plint globalVertex = 0;
    for (plint iTriangle = 0; iTriangle < numTriangles; iTriangle++) {
        for (plint localVertex = 0; localVertex < 3; localVertex++) {
            Array<T, 3> *vertex = &oldTriangles[iTriangle][localVertex];
            VertexSetIterator it = vertexSet.find(VertexSetNode(-1, vertex));
            if (it == vertexSet.end()) {
                globalVertex = numVertices;
                numVertices++;
                VertexSetNode newNode(globalVertex, vertex);
                vertexSet.insert(newNode);
            } else {
                globalVertex = it->i;
            }
            triangles[iTriangle][localVertex] = globalVertex;
        }
    }
    PLB_ASSERT(numVertices >= 3);

    // Order and copy vertex coordinates.
    vertices.resize(numVertices);
    VertexSetConstIterator it = vertexSet.begin();
    for (; it != vertexSet.end(); ++it) {
        vertices[it->i] = *(it->vertex);
    }

    // Create the connectivity information on the triangle global indices
    // that meet on a unified vertex.
    trianglesOnVertex.resize(numVertices);
    for (plint iTriangle = 0; iTriangle < numTriangles; iTriangle++) {
        for (plint localVertex = 0; localVertex < 3; localVertex++) {
            plint globalVertex = triangles[iTriangle][localVertex];
            trianglesOnVertex[globalVertex].push_back(iTriangle);
        }
    }
}

template <typename T>
void ConnectedTriangleSet<T>::swapGeometry(std::vector<Array<T, 3> > &newVertices)
{
    PLB_ASSERT(newVertices.size() == vertices.size());
    vertices.swap(newVertices);
}

template <typename T>
void ConnectedTriangleSet<T>::computeVertexAreaAndUnitNormal(
    plint iVertex, T &area, Array<T, 3> &unitNormal, std::vector<Array<T, 3> > *newVertices,
    plint indexOffset) const
{
    PLB_ASSERT(iVertex >= 0 && iVertex < numVertices);
    area = (T)0;
    unitNormal.resetToZero();
    for (plint iTriangle = 0; iTriangle < (plint)trianglesOnVertex[iVertex].size(); iTriangle++) {
        plint globalTriangle = trianglesOnVertex[iVertex][iTriangle];
        T triangleArea = (T)0;
        Array<T, 3> triangleNormal((T)0, (T)0, (T)0);
        computeTriangleAreaAndUnitNormal(
            globalTriangle, triangleArea, triangleNormal, newVertices, indexOffset);
        area += triangleArea;
        unitNormal += triangleNormal;
    }
    area /= (T)3;

    T normNormal = norm(unitNormal);
    if (!util::isZero(normNormal)) {
        unitNormal /= normNormal;
    } else {
        unitNormal.resetToZero();
    }
}

template <typename T>
void ConnectedTriangleSet<T>::computeTriangleAreaAndUnitNormal(
    plint iTriangle, T &area, Array<T, 3> &unitNormal, std::vector<Array<T, 3> > *newVertices,
    plint indexOffset) const
{
    PLB_ASSERT(iTriangle >= 0 && iTriangle < numTriangles);

    std::vector<Array<T, 3> > const *workingVertices = newVertices;
    if (workingVertices == 0) {
        workingVertices = &vertices;
        indexOffset = 0;
    }

    plint v0 = triangles[iTriangle][0] + indexOffset;
    plint v1 = triangles[iTriangle][1] + indexOffset;
    plint v2 = triangles[iTriangle][2] + indexOffset;

    Array<T, 3> const &a0 = (*workingVertices)[v0];
    Array<T, 3> const &a1 = (*workingVertices)[v1];
    Array<T, 3> const &a2 = (*workingVertices)[v2];

    plb::computeTriangleAreaAndUnitNormal(a0, a1, a2, area, unitNormal);
}

template <typename T>
void ConnectedTriangleSet<T>::writeOFF(
    std::string fname, std::vector<Array<T, 3> > *newVertices, plint indexOffset,
    int numDecimalDigits) const
{
    if (global::mpi().isMainProcessor()) {
        if (numTriangles == 0) {
            return;
        }
        FILE *fp = fopen(fname.c_str(), "w");  // Only ASCII OFF files are written for the moment.
        PLB_ASSERT(fp != 0);

        std::vector<Array<T, 3> > const *workingVertices = newVertices;
        if (workingVertices == 0) {
            workingVertices = &vertices;
            indexOffset = 0;
        }

        fprintf(fp, "OFF\n");
        long numVerticesL = (long)numVertices;
        long numTrianglesL = (long)numTriangles;
        long zero = (long)0;
        fprintf(fp, "%ld %ld %ld\n", numVerticesL, numTrianglesL, zero);

        for (plint iVertex = 0; iVertex < numVertices; iVertex++) {
            double v0 = (double)(*workingVertices)[iVertex + indexOffset][0];
            double v1 = (double)(*workingVertices)[iVertex + indexOffset][1];
            double v2 = (double)(*workingVertices)[iVertex + indexOffset][2];
            fprintf(
                fp, "% .*e % .*e % .*e\n", numDecimalDigits, v0, numDecimalDigits, v1,
                numDecimalDigits, v2);
        }

        for (plint iTriangle = 0; iTriangle < numTriangles; iTriangle++) {
            long three = (long)3;
            long t0 = (long)triangles[iTriangle][0];
            long t1 = (long)triangles[iTriangle][1];
            long t2 = (long)triangles[iTriangle][2];
            fprintf(fp, "%ld %ld %ld %ld\n", three, t0, t1, t2);
        }

        fclose(fp);
    }
}

template <typename T>
TriangleSet<T> *ConnectedTriangleSet<T>::toTriangleSet(
    Precision precision, std::vector<Array<T, 3> > const *newVertices, plint indexOffset) const
{
    if (numTriangles == 0) {
        return 0;
    }

    std::vector<Array<T, 3> > const *workingVertices = newVertices;
    if (workingVertices == 0) {
        workingVertices = &vertices;
        indexOffset = 0;
    }

    typedef typename TriangleSet<T>::Triangle Triangle;
    std::vector<Triangle> vectorOfTriangles;

    for (plint iTriangle = 0; iTriangle < numTriangles; iTriangle++) {
        plint iVertex0 = triangles[iTriangle][0] + indexOffset;
        plint iVertex1 = triangles[iTriangle][1] + indexOffset;
        plint iVertex2 = triangles[iTriangle][2] + indexOffset;
        Array<T, 3> vertex0 = (*workingVertices)[iVertex0];
        Array<T, 3> vertex1 = (*workingVertices)[iVertex1];
        Array<T, 3> vertex2 = (*workingVertices)[iVertex2];
        Triangle triangle(vertex0, vertex1, vertex2);
        vectorOfTriangles.push_back(triangle);
    }

    return new TriangleSet<T>(vectorOfTriangles, precision);
}

template <typename T>
TriangleSet<T> *ConnectedTriangleSet<T>::toTriangleSet(
    T eps, std::vector<Array<T, 3> > const *newVertices, plint indexOffset) const
{
    if (numTriangles == 0) {
        return 0;
    }

    std::vector<Array<T, 3> > const *workingVertices = newVertices;
    if (workingVertices == 0) {
        workingVertices = &vertices;
        indexOffset = 0;
    }

    typedef typename TriangleSet<T>::Triangle Triangle;
    std::vector<Triangle> vectorOfTriangles;

    for (plint iTriangle = 0; iTriangle < numTriangles; iTriangle++) {
        plint iVertex0 = triangles[iTriangle][0] + indexOffset;
        plint iVertex1 = triangles[iTriangle][1] + indexOffset;
        plint iVertex2 = triangles[iTriangle][2] + indexOffset;
        Array<T, 3> vertex0 = (*workingVertices)[iVertex0];
        Array<T, 3> vertex1 = (*workingVertices)[iVertex1];
        Array<T, 3> vertex2 = (*workingVertices)[iVertex2];
        Triangle triangle(vertex0, vertex1, vertex2);
        vectorOfTriangles.push_back(triangle);
    }

    return new TriangleSet<T>(vectorOfTriangles, eps);
}

}  // namespace plb

#endif  // CONNECTED_TRIANGLE_SET_HH
