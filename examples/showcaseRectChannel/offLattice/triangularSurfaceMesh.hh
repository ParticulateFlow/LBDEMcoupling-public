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

#ifndef TRIANGULAR_SURFACE_MESH_HH
#define TRIANGULAR_SURFACE_MESH_HH

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "core/util.h"
#include "io/parallelIO.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "offLattice/triangleSet.h"
#include "offLattice/triangularSurfaceMesh.h"
#include "parallelism/mpiManager.h"

namespace plb {

template <typename T>
const T TriangularSurfaceMesh<T>::eps1 = (sizeof(T) == sizeof(float))
                                             ? std::numeric_limits<float>::epsilon()
                                             : (T)100.0 * std::numeric_limits<T>::epsilon();

template <typename T>
TriangularSurfaceMesh<T>::TriangularSurfaceMesh(
    std::vector<Array<T, 3> > &vertexList_, std::vector<plint> &emanatingEdgeList_,
    std::vector<Edge> &edgeList_, plint numVertices_) :
    vertexList(&vertexList_),
    emanatingEdgeList(&emanatingEdgeList_),
    edgeList(&edgeList_),
    numTriangles((plint)edges().size() / (plint)3),
    numVertices(numVertices_ >= 0 ? numVertices_ : ((plint)vertices().size()))
{
    avoidIntegerPositions();
}

template <typename T>
void TriangularSurfaceMesh<T>::replaceVertex(plint iVertex, Array<T, 3> const &newPosition)
{
    PLB_PRECONDITION(iVertex < getNumVertices());
    (*vertexList)[iVertex] = newPosition;
    // avoidIntegerPosition(iVertex);
}

template <typename T>
void TriangularSurfaceMesh<T>::resetVertices(Array<T, 3> const &defaultVertex)
{
    for (plint iVertex = 0; iVertex < getNumVertices(); ++iVertex) {
        vertices()[iVertex] = defaultVertex;
    }
}

template <typename T>
inline void TriangularSurfaceMesh<T>::assertVertex(Array<T, 3> const &vertex) const
{
#ifdef PLB_DEBUG
    if (global::IOpolicy().stlFilesHaveLowerBound()) {
        double bound = global::IOpolicy().getLowerBoundForStlFiles();
        PLB_ASSERT(!(vertex[0] < bound && vertex[1] < bound && vertex[2] < bound));
    }
#endif  // PLB_DEBUG
}

template <typename T>
inline Array<T, 3> const &TriangularSurfaceMesh<T>::getVertex(
    plint iTriangle, int localVertex) const
{
    PLB_ASSERT(
        iTriangle >= 0 && iTriangle < numTriangles
        && (localVertex == 0 || localVertex == 1 || localVertex == 2));
    Array<T, 3> const &vertex(
        vertices()[edges()[3 * iTriangle + ((localVertex == 0) ? 2 : localVertex - 1)].pv]);
    assertVertex(vertex);
    return vertex;
}

template <typename T>
inline Array<T, 3> &TriangularSurfaceMesh<T>::getVertex(plint iTriangle, int localVertex)
{
    PLB_ASSERT(
        iTriangle >= 0 && iTriangle < numTriangles
        && (localVertex == 0 || localVertex == 1 || localVertex == 2));
    Array<T, 3> &vertex(
        vertices()[edges()[3 * iTriangle + ((localVertex == 0) ? 2 : localVertex - 1)].pv]);
    assertVertex(vertex);
    return vertex;
}

template <typename T>
bool TriangularSurfaceMesh<T>::isValidVertex(plint iTriangle, int localVertex) const
{
    PLB_ASSERT(
        iTriangle >= 0 && iTriangle < numTriangles
        && (localVertex == 0 || localVertex == 1 || localVertex == 2));
    if (global::IOpolicy().stlFilesHaveLowerBound()) {
        double bound = global::IOpolicy().getLowerBoundForStlFiles();
        Array<T, 3> const &vertex(
            vertices()[edges()[3 * iTriangle + ((localVertex == 0) ? 2 : localVertex - 1)].pv]);
        return !(vertex[0] < bound && vertex[1] < bound && vertex[2] < bound);
    } else {
        return true;
    }
}

template <typename T>
inline Array<T, 3> const &TriangularSurfaceMesh<T>::getVertex(plint iVertex) const
{
    PLB_ASSERT(iVertex >= 0 && iVertex < numVertices);
    Array<T, 3> const &vertex(vertices()[iVertex]);
    assertVertex(vertex);
    return vertex;
}

template <typename T>
inline Array<T, 3> &TriangularSurfaceMesh<T>::getVertex(plint iVertex)
{
    PLB_ASSERT(iVertex >= 0 && iVertex < numVertices);
    Array<T, 3> &vertex(vertices()[iVertex]);
    assertVertex(vertex);
    return vertex;
}

template <typename T>
inline bool TriangularSurfaceMesh<T>::isValidVertex(plint iVertex) const
{
    PLB_ASSERT(iVertex >= 0 && iVertex < numVertices);
    if (global::IOpolicy().stlFilesHaveLowerBound()) {
        Array<T, 3> const &vertex(vertices()[iVertex]);
        double bound = global::IOpolicy().getLowerBoundForStlFiles();
        return !(vertex[0] < bound && vertex[1] < bound && vertex[2] < bound);
    } else {
        return true;
    }
}

template <typename T>
void TriangularSurfaceMesh<T>::computeBoundingBox(
    Array<T, 2> &xRange, Array<T, 2> &yRange, Array<T, 2> &zRange) const
{
    T maxVal = std::numeric_limits<T>::max();
    xRange = Array<T, 2>(maxVal, -maxVal);
    yRange = Array<T, 2>(maxVal, -maxVal);
    zRange = Array<T, 2>(maxVal, -maxVal);
    for (plint iVertex = 0; iVertex < getNumVertices(); ++iVertex) {
        Array<T, 3> const &vertex = getVertex(iVertex);
        xRange[0] = std::min(xRange[0], vertex[0]);
        xRange[1] = std::max(xRange[1], vertex[0]);
        yRange[0] = std::min(yRange[0], vertex[1]);
        yRange[1] = std::max(yRange[1], vertex[1]);
        zRange[0] = std::min(zRange[0], vertex[2]);
        zRange[1] = std::max(zRange[1], vertex[2]);
    }
}

template <typename T>
void TriangularSurfaceMesh<T>::translate(Array<T, 3> const &vector)
{
    if (util::isZero(norm(vector))) {
        return;
    }

    for (plint i = 0; i < numVertices; i++) {
        getVertex(i) += vector;
    }
}

template <typename T>
void TriangularSurfaceMesh<T>::scale(T alpha)
{
    PLB_ASSERT(!util::isZero(alpha));
    if (util::isOne(alpha)) {
        return;
    }

    for (plint i = 0; i < numVertices; i++) {
        getVertex(i) *= alpha;
    }
}

template <typename T>
void TriangularSurfaceMesh<T>::rotate(T phi, T theta, T psi)
{
#ifdef PLB_DEBUG
    static const T pi = std::acos((T)-1);
#endif
    PLB_ASSERT(util::greaterEqual(theta, (T)0) && util::lessEqual(theta, pi));

    T cosPhi = std::cos(phi);
    T sinPhi = std::sin(phi);
    T cosTheta = std::cos(theta);
    T sinTheta = std::sin(theta);
    T cosPsi = std::cos(psi);
    T sinPsi = std::sin(psi);

    T a[3][3];
    a[0][0] = (T)1.0;
    a[0][1] = (T)0.0;
    a[0][2] = (T)0.0;
    a[1][0] = (T)0.0;
    a[1][1] = cosTheta;
    a[1][2] = -sinTheta;
    a[2][0] = (T)0.0;
    a[2][1] = sinTheta;
    a[2][2] = cosTheta;

    T b[3][3];
    b[0][0] = cosPhi;
    b[0][1] = -sinPhi;
    b[0][2] = (T)0.0;
    b[1][0] = sinPhi;
    b[1][1] = cosPhi;
    b[1][2] = (T)0.0;
    b[2][0] = (T)0.0;
    b[2][1] = (T)0.0;
    b[2][2] = (T)1.0;

    T c[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            c[i][j] = (T)0.0;
            for (int k = 0; k < 3; k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }

    b[0][0] = cosPsi;
    b[0][1] = -sinPsi;
    b[0][2] = (T)0.0;
    b[1][0] = sinPsi;
    b[1][1] = cosPsi;
    b[1][2] = (T)0.0;
    b[2][0] = (T)0.0;
    b[2][1] = (T)0.0;
    b[2][2] = (T)1.0;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            a[i][j] = (T)0.0;
            for (int k = 0; k < 3; k++) {
                a[i][j] += b[i][k] * c[k][j];
            }
        }
    }

    for (plint iVertex = 0; iVertex < numVertices; iVertex++) {
        Array<T, 3> x = getVertex(iVertex);
        for (int i = 0; i < 3; i++) {
            getVertex(iVertex)[i] = (T)0.0;
            for (int j = 0; j < 3; j++) {
                getVertex(iVertex)[i] += a[i][j] * x[j];
            }
        }
    }
}

template <typename T>
void TriangularSurfaceMesh<T>::smooth(plint maxiter, T relax, bool isMeasureWeighted)
{
    PLB_ASSERT(maxiter >= 0);
    PLB_ASSERT(util::greaterEqual(relax, (T)0) && util::lessEqual(relax, (T)1));

    if (maxiter <= 0) {
        return;
    }

    std::vector<Array<T, 3> > *bp0;
    bp0 = vertexList;

    std::vector<Array<T, 3> > buf;
    buf.resize(numVertices);

    std::vector<Array<T, 3> > *bp1;
    bp1 = &buf;

    if (isMeasureWeighted) {
        for (plint iter = 0; iter < maxiter; iter++) {
            for (plint iVertex = 0; iVertex < numVertices; iVertex++) {
                Array<T, 3> iVertexPos = (*bp0)[iVertex];
                std::vector<plint> neighborVertexIds = getNeighborVertexIds(iVertex);

                if (isInteriorVertex(iVertex)) {
                    pluint sizeIds = neighborVertexIds.size();
                    T area = (T)0.0;
                    Array<T, 3> tmp((T)0.0, (T)0.0, (T)0.0);
                    for (pluint i = 0; i < sizeIds; i++) {
                        plint jVertex = neighborVertexIds[i];
                        plint kVertex =
                            (i + 1 < sizeIds) ? neighborVertexIds[i + 1] : neighborVertexIds[0];

                        Array<T, 3> jVertexPos = (*bp0)[jVertex];
                        Array<T, 3> kVertexPos = (*bp0)[kVertex];

                        Array<T, 3> middleVertexPos =
                            (iVertexPos + jVertexPos + kVertexPos) / (T)3.0;

                        Array<T, 3> eij = jVertexPos - iVertexPos;
                        Array<T, 3> eik = kVertexPos - iVertexPos;
                        Array<T, 3> n;
                        crossProduct(eij, eik, n);
                        T locArea = (T)0.5 * norm(n);

                        area += locArea;
                        tmp += middleVertexPos * locArea;
                    }
                    (*bp1)[iVertex] = ((T)1.0 - relax) * (*bp0)[iVertex] + relax * tmp / area;
                } else {
                    plint jVertex = -1, kVertex = -1;
                    int counter = 0;
                    for (pluint i = 0; i < neighborVertexIds.size(); i++) {
                        plint vertexId = neighborVertexIds[i];
                        if (isBoundaryVertex(vertexId)) {
                            if (isBoundaryEdge(iVertex, vertexId)) {
                                counter++;
                                if (counter == 1) {
                                    jVertex = vertexId;
                                } else if (counter == 2) {
                                    kVertex = vertexId;
                                } else {
                                    PLB_ASSERT(
                                        false);  // Problem with the boundary of the surface mesh.
                                }
                            }
                        }
                    }

                    PLB_ASSERT(counter == 2);  // Problem with the boundary of the surface mesh.

                    Array<T, 3> jVertexPos = (*bp0)[jVertex];
                    Array<T, 3> kVertexPos = (*bp0)[kVertex];

                    Array<T, 3> middleVertexPos0 = (iVertexPos + jVertexPos) / (T)2.0;
                    Array<T, 3> middleVertexPos1 = (iVertexPos + kVertexPos) / (T)2.0;

                    T length0 = norm(jVertexPos - iVertexPos);
                    T length1 = norm(kVertexPos - iVertexPos);

                    (*bp1)[iVertex] =
                        ((T)1.0 - relax) * (*bp0)[iVertex]
                        + relax * (middleVertexPos0 * length0 + middleVertexPos1 * length1)
                              / (length0 + length1);
                }
            }
            // Swap buffers.
            std::vector<Array<T, 3> > *tmp = bp0;
            bp0 = bp1;
            bp1 = tmp;
        }
    } else {
        for (plint iter = 0; iter < maxiter; iter++) {
            for (plint iVertex = 0; iVertex < numVertices; iVertex++) {
                std::vector<plint> neighborVertexIds = getNeighborVertexIds(iVertex);

                if (isInteriorVertex(iVertex)) {
                    pluint sizeIds = neighborVertexIds.size();
                    Array<T, 3> tmp((T)0.0, (T)0.0, (T)0.0);
                    for (pluint i = 0; i < sizeIds; i++) {
                        tmp += (*bp0)[neighborVertexIds[i]];
                    }
                    (*bp1)[iVertex] = ((T)1.0 - relax) * (*bp0)[iVertex] + relax * tmp / (T)sizeIds;
                } else {
                    plint jVertex = -1, kVertex = -1;
                    int counter = 0;
                    for (pluint i = 0; i < neighborVertexIds.size(); i++) {
                        plint vertexId = neighborVertexIds[i];
                        if (isBoundaryVertex(vertexId)) {
                            if (isBoundaryEdge(iVertex, vertexId)) {
                                counter++;
                                if (counter == 1) {
                                    jVertex = vertexId;
                                } else if (counter == 2) {
                                    kVertex = vertexId;
                                } else {
                                    PLB_ASSERT(
                                        false);  // Problem with the boundary of the surface mesh.
                                }
                            }
                        }
                    }

                    PLB_ASSERT(counter == 2);  // Problem with the boundary of the surface mesh.

                    (*bp1)[iVertex] = ((T)1.0 - relax) * (*bp0)[iVertex]
                                      + relax * ((*bp0)[jVertex] + (*bp0)[kVertex]) / (T)2.0;
                }
            }
            // Swap buffers.
            std::vector<Array<T, 3> > *tmp = bp0;
            bp0 = bp1;
            bp1 = tmp;
        }
    }

    // Final copy of vertex positions.
    if (vertexList != bp0)
        for (plint iVertex = 0; iVertex < numVertices; iVertex++)
            (*vertexList)[iVertex] = (*bp0)[iVertex];
}

template <typename T>
inline plint TriangularSurfaceMesh<T>::getVertexId(plint iTriangle, plint localVertex) const
{
    PLB_ASSERT(
        iTriangle >= 0 && iTriangle < numTriangles
        && (localVertex == 0 || localVertex == 1 || localVertex == 2));
    return edges()[3 * iTriangle + ((localVertex == 0) ? 2 : localVertex - 1)].pv;
}

template <typename T>
std::vector<plint> TriangularSurfaceMesh<T>::getNeighborVertexIds(plint iVertex) const
{
    PLB_ASSERT(iVertex >= 0 && iVertex < numVertices);

    std::vector<plint> neighborVertexIds;

    plint ee = emanatingEdges()[iVertex];
    neighborVertexIds.push_back(edges()[ee].pv);

    plint pe = prev(ee);
    plint e = pe;
    if ((e = edges()[e].ne) < 0) {
        neighborVertexIds.push_back(edges()[prev(pe)].pv);
        e = ee;
    }

    while (e != ee) {
        neighborVertexIds.push_back(edges()[e].pv);

        e = pe = prev(e);
        if ((e = edges()[e].ne) < 0) {
            neighborVertexIds.push_back(edges()[prev(pe)].pv);
            e = ee;
        }
    }

    PLB_ASSERT(neighborVertexIds.size() > 1);  // Problem with the topology of the surface mesh.

    return neighborVertexIds;
}

template <typename T>
std::vector<plint> TriangularSurfaceMesh<T>::getNeighborVertexIds(
    plint iVertex, plint jVertex) const
{
    std::vector<plint> neighborVertexIds;
    std::vector<plint> adjacentTriangleIds = getAdjacentTriangleIds(iVertex, jVertex);
    std::vector<plint>::iterator tit = adjacentTriangleIds.begin();
    for (; tit != adjacentTriangleIds.end(); ++tit) {
        plint iTriangle = *tit;
        plint id0 = getVertexId(iTriangle, 0);
        plint id1 = getVertexId(iTriangle, 1);
        plint id2 = getVertexId(iTriangle, 2);
        plint kVertex = (iVertex != id0 && jVertex != id0)
                            ? id0
                            : ((iVertex != id1 && jVertex != id1) ? id1 : id2);
        neighborVertexIds.push_back(kVertex);
    }

#ifdef PLB_DEBUG
    plint size = neighborVertexIds.size();
#endif
    PLB_ASSERT(size == 1 || size == 2);  // Problem with the topology of the surface mesh.

    return neighborVertexIds;
}

template <typename T>
std::vector<plint> TriangularSurfaceMesh<T>::getNeighborTriangleIds(plint iVertex) const
{
    PLB_ASSERT(iVertex >= 0 && iVertex < numVertices);

    std::vector<plint> neighborTriangleIds;

    plint ee = emanatingEdges()[iVertex];
    neighborTriangleIds.push_back(ee / 3);

    plint e = prev(ee);
    if ((e = edges()[e].ne) < 0)
        e = ee;

    while (e != ee) {
        neighborTriangleIds.push_back(e / 3);

        e = prev(e);
        if ((e = edges()[e].ne) < 0)
            e = ee;
    }

    PLB_ASSERT(neighborTriangleIds.size() != 0);  // Problem with the topology of the surface mesh.

    return neighborTriangleIds;
}

template <typename T>
std::vector<plint> TriangularSurfaceMesh<T>::getAdjacentTriangleIds(plint iTriangle) const
{
    PLB_ASSERT(iTriangle >= 0 && iTriangle < numTriangles);

    std::vector<plint> adjacentTriangleIds;
    for (int localEdge = 0; localEdge < 3; localEdge++) {
        plint e = edges()[3 * iTriangle + localEdge].ne;
        if (e >= 0)
            adjacentTriangleIds.push_back(e / 3);
    }
    return adjacentTriangleIds;
}

template <typename T>
std::vector<plint> TriangularSurfaceMesh<T>::getAdjacentTriangleIds(
    plint iVertex, plint jVertex) const
{
    PLB_ASSERT(
        iVertex >= 0 && iVertex < numVertices && jVertex >= 0 && jVertex < numVertices
        && iVertex != jVertex);

#ifdef PLB_DEBUG
    std::vector<plint> neighborVertexIds = getNeighborVertexIds(iVertex);
    std::vector<plint>::iterator vit;
    vit = find(neighborVertexIds.begin(), neighborVertexIds.end(), jVertex);
    PLB_ASSERT(vit != neighborVertexIds.end());  // Vertices do not belong to the same edge.
#endif                                           // PLB_DEBUG

    std::vector<plint> adjacentTriangleIds;
    std::vector<plint> neighborTriangleIds = getNeighborTriangleIds(iVertex);
    std::vector<plint>::iterator tit = neighborTriangleIds.begin();
    for (; tit != neighborTriangleIds.end(); ++tit) {
        plint iTriangle = *tit;
        plint id0 = getVertexId(iTriangle, 0);
        plint id1 = getVertexId(iTriangle, 1);
        plint id2 = getVertexId(iTriangle, 2);
        if ((iVertex == id0 || iVertex == id1 || iVertex == id2)
            && (jVertex == id0 || jVertex == id1 || jVertex == id2))
            adjacentTriangleIds.push_back(iTriangle);
    }

    PLB_ASSERT(
        adjacentTriangleIds.size() == 1
        || adjacentTriangleIds.size() == 2);  // Problem with the topology of the surface mesh.

    return adjacentTriangleIds;
}

template <typename T>
Array<T, 3> TriangularSurfaceMesh<T>::computeTriangleNormal(
    plint iTriangle, bool isAreaWeighted) const
{
    PLB_ASSERT(iTriangle >= 0 && iTriangle < numTriangles);

    Array<T, 3> const &v0 = getVertex(iTriangle, 0);
    Array<T, 3> const &v1 = getVertex(iTriangle, 1);
    Array<T, 3> const &v2 = getVertex(iTriangle, 2);

    return plb::computeTriangleNormal(v0, v1, v2, isAreaWeighted);
}

template <typename T>
Array<T, 3> TriangularSurfaceMesh<T>::computeTriangleNormal(
    plint iVertex, plint jVertex, plint kVertex, bool isAreaWeighted) const
{
    PLB_ASSERT(
        iVertex >= 0 && iVertex < numVertices && jVertex >= 0 && jVertex < numVertices
        && kVertex >= 0 && kVertex < numVertices && iVertex != jVertex && iVertex != kVertex
        && jVertex != kVertex);

    plint id0 = iVertex;
    plint id1, id2;

    std::vector<plint> neighborVertexIds = getNeighborVertexIds(iVertex);
    std::vector<plint>::iterator it;
    it = find(neighborVertexIds.begin(), neighborVertexIds.end(), jVertex);
    PLB_ASSERT(it != neighborVertexIds.end());  // Vertices do not belong to the same triangle.

    plint prevVertex, nextVertex;
    if (jVertex == neighborVertexIds[0]) {
        prevVertex = neighborVertexIds[neighborVertexIds.size() - 1];
        nextVertex = neighborVertexIds[1];
    } else if (jVertex == neighborVertexIds[neighborVertexIds.size() - 1]) {
        prevVertex = neighborVertexIds[neighborVertexIds.size() - 2];
        nextVertex = neighborVertexIds[0];
    } else {
        prevVertex = *(it - 1);
        nextVertex = *(it + 1);
    }

    if (kVertex == prevVertex) {
        id1 = kVertex;
        id2 = jVertex;
    } else if (kVertex == nextVertex) {
        id1 = jVertex;
        id2 = kVertex;
    } else {
        PLB_ASSERT(false);  // Vertices do not belong to the same triangle.
    }

    Array<T, 3> const &v0 = getVertex(id0);
    Array<T, 3> const &v1 = getVertex(id1);
    Array<T, 3> const &v2 = getVertex(id2);

    return plb::computeTriangleNormal(v0, v1, v2, isAreaWeighted);
}

template <typename T>
Array<T, 3> TriangularSurfaceMesh<T>::computeEdgeNormal(
    plint iVertex, plint jVertex, bool isAreaWeighted) const
{
    Array<T, 3> n;
    std::vector<plint> adjacentTriangleIds = getAdjacentTriangleIds(iVertex, jVertex);
    std::vector<plint>::iterator tit = adjacentTriangleIds.begin();
    for (n.resetToZero(); tit != adjacentTriangleIds.end(); ++tit) {
        n += computeTriangleNormal(*tit, isAreaWeighted);
    }

    T normN = norm(n);
    if (!util::isZero(normN)) {
        n /= normN;
    } else {
        n.resetToZero();
    }
    return n;
}

template <typename T>
Array<T, 3> TriangularSurfaceMesh<T>::computeVertexNormal(plint iVertex, bool isAreaWeighted) const
{
    std::vector<plint> neighborTriangleIds = getNeighborTriangleIds(iVertex);
    Array<T, 3> n;
    std::vector<plint>::iterator it = neighborTriangleIds.begin();
    for (n.resetToZero(); it != neighborTriangleIds.end(); ++it) {
        n += computeTriangleNormal(*it, isAreaWeighted);
    }

    T normN = norm(n);
    if (!util::isZero(normN)) {
        n /= normN;
    } else {
        n.resetToZero();
    }
    return n;
}

template <typename T>
Array<T, 3> TriangularSurfaceMesh<T>::computeContinuousNormal(
    Array<T, 3> const &p, plint iTriangle, bool isAreaWeighted) const
{
    T area = computeTriangleArea(iTriangle);
    if (util::isZero(area)) {
        return Array<T, 3>((T)0, (T)0, (T)0);
    }

    plint id0 = getVertexId(iTriangle, 0);
    plint id1 = getVertexId(iTriangle, 1);
    plint id2 = getVertexId(iTriangle, 2);

    Array<T, 3> v0 = getVertex(id0);
    Array<T, 3> v1 = getVertex(id1);
    Array<T, 3> v2 = getVertex(id2);

    Array<T, 3> ep0 = v0 - p;
    Array<T, 3> ep1 = v1 - p;
    Array<T, 3> ep2 = v2 - p;

    Array<T, 3> n;
    crossProduct(ep1, ep2, n);
    T area0 = (T)0.5 * norm(n);
    crossProduct(ep2, ep0, n);
    T area1 = (T)0.5 * norm(n);

    T u = area0 / area;
    T v = area1 / area;

    Array<T, 3> n0 = computeVertexNormal(id0, isAreaWeighted);
    Array<T, 3> n1 = computeVertexNormal(id1, isAreaWeighted);
    Array<T, 3> n2 = computeVertexNormal(id2, isAreaWeighted);

    n = u * n0 + v * n1 + ((T)1. - u - v) * n2;
    T normN = norm(n);
    PLB_ASSERT(!util::isZero(normN));
    n /= normN;
    return n;
}

template <typename T>
T TriangularSurfaceMesh<T>::computeTriangleArea(plint iTriangle) const
{
    Array<T, 3> const &v0 = getVertex(iTriangle, 0);
    Array<T, 3> const &v1 = getVertex(iTriangle, 1);
    Array<T, 3> const &v2 = getVertex(iTriangle, 2);

    return plb::computeTriangleArea(v0, v1, v2);
}

template <typename T>
T TriangularSurfaceMesh<T>::computeTriangleArea(plint iVertex, plint jVertex, plint kVertex) const
{
    PLB_ASSERT(
        iVertex >= 0 && iVertex < numVertices && jVertex >= 0 && jVertex < numVertices
        && kVertex >= 0 && kVertex < numVertices && iVertex != jVertex && iVertex != kVertex
        && jVertex != kVertex);

#ifdef PLB_DEBUG
    std::vector<plint> neighborVertexIds = getNeighborVertexIds(iVertex);
    std::vector<plint>::iterator it;
    it = find(neighborVertexIds.begin(), neighborVertexIds.end(), jVertex);
    PLB_ASSERT(it != neighborVertexIds.end());  // Vertices do not belong to the same triangle.

    plint prevVertex, nextVertex;
    if (jVertex == neighborVertexIds[0]) {
        prevVertex = neighborVertexIds[neighborVertexIds.size() - 1];
        nextVertex = neighborVertexIds[1];
    } else if (jVertex == neighborVertexIds[neighborVertexIds.size() - 1]) {
        prevVertex = neighborVertexIds[neighborVertexIds.size() - 2];
        nextVertex = neighborVertexIds[0];
    } else {
        prevVertex = *(it - 1);
        nextVertex = *(it + 1);
    }

    PLB_ASSERT(
        kVertex == prevVertex
        || kVertex == nextVertex);  // Vertices do not belong to the same triangle.
#endif                              // PLB_DEBUG

    Array<T, 3> const &v0 = getVertex(iVertex);
    Array<T, 3> const &v1 = getVertex(jVertex);
    Array<T, 3> const &v2 = getVertex(kVertex);

    return plb::computeTriangleArea(v0, v1, v2);
}

template <typename T>
T TriangularSurfaceMesh<T>::computeEdgeArea(plint iVertex, plint jVertex) const
{
    T area = T();
    std::vector<plint> adjacentTriangleIds = getAdjacentTriangleIds(iVertex, jVertex);
    std::vector<plint>::iterator tit = adjacentTriangleIds.begin();
    for (; tit != adjacentTriangleIds.end(); ++tit)
        area += computeTriangleArea(*tit);

    return area / 3.0;
}

template <typename T>
T TriangularSurfaceMesh<T>::computeVertexArea(plint iVertex) const
{
    std::vector<plint> neighborTriangleIds = getNeighborTriangleIds(iVertex);
    T area = T();
    std::vector<plint>::iterator it = neighborTriangleIds.begin();
    for (; it != neighborTriangleIds.end(); ++it) {
        area += computeTriangleArea(*it);
    }

    return area / 3.0;
}

template <typename T>
T TriangularSurfaceMesh<T>::computeEdgeLength(plint iVertex, plint jVertex) const
{
#ifdef PLB_DEBUG
    std::vector<plint> neighborVertexIds = getNeighborVertexIds(iVertex);
    std::vector<plint>::iterator it;
    it = find(neighborVertexIds.begin(), neighborVertexIds.end(), jVertex);
    PLB_ASSERT(it != neighborVertexIds.end());  // Vertices do not belong to the same edge.
#endif                                          // PLB_DEBUG

    return norm(getVertex(iVertex) - getVertex(jVertex));
}

template <typename T>
T TriangularSurfaceMesh<T>::computeDihedralAngle(plint iVertex, plint jVertex) const
{
    std::vector<plint> adjacentTriangleIds = getAdjacentTriangleIds(iVertex, jVertex);
    if (adjacentTriangleIds.size() == 1)
        return T();

    return angleBetweenVectors(
        computeTriangleNormal(adjacentTriangleIds[0]),
        computeTriangleNormal(adjacentTriangleIds[1]));
}

template <typename T>
T TriangularSurfaceMesh<T>::computeEdgeTileSpan(plint iVertex, plint jVertex) const
{
    std::vector<plint> neighborVertexIds = getNeighborVertexIds(iVertex, jVertex);

    Array<T, 3> v0 = getVertex(neighborVertexIds[0]);
    Array<T, 3> v1 = getVertex(iVertex);
    Array<T, 3> v2 = getVertex(jVertex);

    Array<T, 3> v01 = v1 - v0;
    Array<T, 3> v21 = v1 - v2;
    T angle_012 = angleBetweenVectors(v21, v01);

    T span = std::fabs(std::sin(angle_012)) * norm(v01);

    if (neighborVertexIds.size() == 1)
        return span / 6.0;

    Array<T, 3> v3 = getVertex(neighborVertexIds[1]);

    Array<T, 3> v32 = v2 - v3;
    Array<T, 3> v12 = -v21;
    T angle_321 = angleBetweenVectors(v12, v32);

    span += std::fabs(std::sin(angle_321)) * norm(v32);

    return span / 6.0;
}

template <typename T>
TriangleSet<T> TriangularSurfaceMesh<T>::toTriangleSet(Precision precision) const
{
    typedef typename TriangleSet<T>::Triangle Triangle;
    std::vector<Triangle> triangles(numTriangles);
    for (plint i = 0; i < numTriangles; i++) {
        Triangle triangle;
        triangle[0] = getVertex(i, 0);
        triangle[1] = getVertex(i, 1);
        triangle[2] = getVertex(i, 2);
        triangles[i] = triangle;
    }
    return TriangleSet<T>(triangles, precision);
}

template <typename T>
TriangleSet<T> TriangularSurfaceMesh<T>::toTriangleSet(T eps) const
{
    typedef typename TriangleSet<T>::Triangle Triangle;
    std::vector<Triangle> triangles(numTriangles);
    for (plint i = 0; i < numTriangles; i++) {
        Triangle triangle;
        triangle[0] = getVertex(i, 0);
        triangle[1] = getVertex(i, 1);
        triangle[2] = getVertex(i, 2);
        triangles[i] = triangle;
    }
    return TriangleSet<T>(triangles, eps);
}

template <typename T>
void TriangularSurfaceMesh<T>::writeAsciiSTL(std::string fname, T dx, int numDecimalDigits) const
{
    Array<T, 3> location((T)0, (T)0, (T)0);
    writeAsciiSTL(fname, dx, location, numDecimalDigits);
}

template <typename T>
void TriangularSurfaceMesh<T>::writeAsciiSTL(
    std::string fname, T dx, Array<T, 3> location, int numDecimalDigits) const
{
    // Output only from one MPI process.
    if (!global::mpi().isMainProcessor()) {
        return;
    }
    if (numTriangles == 0) {
        return;
    }
    FILE *fp = fopen(fname.c_str(), "w");
    PLB_ASSERT(fp != 0);

    fprintf(fp, "solid plb\n");
    bool isAreaWeighted = false;
    int d = numDecimalDigits;
    for (plint i = 0; i < numTriangles; i++) {
        Array<double, 3> n = computeTriangleNormal(i, isAreaWeighted);
        Array<double, 3> v;
        fprintf(fp, "  facet normal % .*e % .*e % .*e\n", d, n[0], d, n[1], d, n[2]);
        fprintf(fp, "    outer loop\n");
        v = getVertex(i, 0) * dx + location;
        fprintf(fp, "      vertex % .*e % .*e % .*e\n", d, v[0], d, v[1], d, v[2]);
        v = getVertex(i, 1) * dx + location;
        fprintf(fp, "      vertex % .*e % .*e % .*e\n", d, v[0], d, v[1], d, v[2]);
        v = getVertex(i, 2) * dx + location;
        fprintf(fp, "      vertex % .*e % .*e % .*e\n", d, v[0], d, v[1], d, v[2]);
        fprintf(fp, "    endloop\n");
        fprintf(fp, "  endfacet\n");
    }
    fprintf(fp, "endsolid plb\n");

    fclose(fp);
}

template <typename T>
void TriangularSurfaceMesh<T>::writeBinarySTL(std::string fname, T dx) const
{
    Array<T, 3> location((T)0, (T)0, (T)0);
    writeBinarySTL(fname, dx, location);
}

template <typename T>
void TriangularSurfaceMesh<T>::writeBinarySTL(std::string fname, T dx, Array<T, 3> location) const
{
    // Output only from one MPI process.
    if (!global::mpi().isMainProcessor()) {
        return;
    }
    if (numTriangles == 0) {
        return;
    }
    FILE *fp = fopen(fname.c_str(), "wb");
    PLB_ASSERT(fp != 0);

    unsigned int nt = (unsigned int)numTriangles;
    unsigned short abc = 0;
    char buf[80];

    for (int i = 0; i < 80; i++)
        buf[i] = '\0';

    fwrite(buf, sizeof(char), 80, fp);
    fwrite(&nt, sizeof(unsigned int), 1, fp);
    bool isAreaWeighted = false;
    for (plint i = 0; i < numTriangles; i++) {
        Array<T, 3> vertex;
        Array<T, 3> normal = computeTriangleNormal(i, isAreaWeighted);
        float n[3];
        n[0] = normal[0];
        n[1] = normal[1];
        n[2] = normal[2];
        fwrite((void *)n, sizeof(float), 3, fp);
        vertex = getVertex(i, 0) * dx + location;
        float v[3];
        v[0] = vertex[0];
        v[1] = vertex[1];
        v[2] = vertex[2];
        fwrite((void *)v, sizeof(float), 3, fp);
        vertex = getVertex(i, 1) * dx + location;
        v[0] = vertex[0];
        v[1] = vertex[1];
        v[2] = vertex[2];
        fwrite((void *)v, sizeof(float), 3, fp);
        vertex = getVertex(i, 2) * dx + location;
        v[0] = vertex[0];
        v[1] = vertex[1];
        v[2] = vertex[2];
        fwrite((void *)v, sizeof(float), 3, fp);
        fwrite(&abc, sizeof(unsigned short), 1, fp);
    }

    fclose(fp);
}

template <typename T>
inline bool TriangularSurfaceMesh<T>::isBoundaryVertex(plint iVertex) const
{
    PLB_ASSERT(iVertex >= 0 && iVertex < numVertices);
    if (edges()[emanatingEdges()[iVertex]].ne < 0) {
        return true;
    }
    return false;
}

template <typename T>
inline bool TriangularSurfaceMesh<T>::isInteriorVertex(plint iVertex) const
{
    PLB_ASSERT(iVertex >= 0 && iVertex < numVertices);
    if (edges()[emanatingEdges()[iVertex]].ne >= 0)
        return true;
    return false;
}

template <typename T>
inline bool TriangularSurfaceMesh<T>::isBoundaryEdge(plint iVertex, plint jVertex) const
{
    if (getAdjacentTriangleIds(iVertex, jVertex).size() == 2)
        return false;

    return true;
}

template <typename T>
inline bool TriangularSurfaceMesh<T>::isInteriorEdge(plint iVertex, plint jVertex) const
{
    return !isBoundaryEdge(iVertex, jVertex);
}

template <typename T>
int TriangularSurfaceMesh<T>::pointOnTriangle(
    Array<T, 3> const &point1, Array<T, 3> const &point2, int flag, plint iTriangle,
    Array<T, 3> &intersection, Array<T, 3> &normal, T &distance) const
{
    if (!((flag == 0 || flag == 1 || flag == 2) && iTriangle >= 0 && iTriangle < numTriangles)) {
        pcout << "iTriangle=" << iTriangle << std::endl;
        pcout << "numTriangles=" << numTriangles << std::endl;
        pcout << "flag=" << flag << std::endl;
    }
    PLB_ASSERT((flag == 0 || flag == 1 || flag == 2) && iTriangle >= 0 && iTriangle < numTriangles);

    Array<T, 3> v0 = getVertex(iTriangle, 0);  // Triangle vertex coordinates
    Array<T, 3> v1 = getVertex(iTriangle, 1);
    Array<T, 3> v2 = getVertex(iTriangle, 2);

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
    if (((flag == 0 || flag == 1) && util::isZero(t, eps1))
        || ((flag == 0) && util::isOne(t, eps1))) {
        // Shift the vertex, and recompute everything that needs to be
        //   recomputed. v1 and v2 are not used any more, so there is no
        //   need to recompute them.
        v0 += (T)2 * eps1 * normal;
        num = dot(v0, normal) - dot(point1, normal);
        t = num / denom;
    }

    if (util::isZero(denom, eps1)) {
        if (util::isZero(num, eps1)) {
            return -1;  // Line belongs to the plane
        } else {
            return 0;  // Line does not intersect the plane
        }
    }

    if (flag == 0) {  // For intersection with a line segment
        if (util::lessThan(t, (T)0, eps1) || util::greaterThan(t, (T)1, eps1)) {
            return 0;
        }
    } else if (flag == 1) {  // For intersection with a half-line
        if (util::lessThan(t, (T)0, eps1)) {
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

    bool ueq0 = util::isZero(u, eps1);
    bool ueq1 = util::isOne(u, eps1);
    bool veq0 = util::isZero(v, eps1);
    bool veq1 = util::isOne(v, eps1);
    bool upveq1 = util::isOne(upv, eps1);

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
        plint iVertex = getVertexId(iTriangle, is_in_edge);
        plint jVertex =
            (is_in_edge == 2) ? getVertexId(iTriangle, 0) : getVertexId(iTriangle, is_in_edge + 1);

        normal = computeEdgeNormal(iVertex, jVertex);
    } else if (is_vertex != -1) {  // The point is a vertex
        normal = computeVertexNormal(getVertexId(iTriangle, is_vertex));
    }

    return 1;
}

template <typename T>
bool TriangularSurfaceMesh<T>::segmentIntersectsTriangle(
    Array<T, 3> const &point1, Array<T, 3> const &point2, plint iTriangle) const
{
    PLB_ASSERT(iTriangle >= 0 && iTriangle < numTriangles);

    T v0[3], v1[3], v2[3];  // Triangle vertex coordinates

    Array<T, 3> tv0 = getVertex(iTriangle, 0);
    tv0.to_cArray(v0);
    Array<T, 3> tv1 = getVertex(iTriangle, 1);
    tv1.to_cArray(v1);
    Array<T, 3> tv2 = getVertex(iTriangle, 2);
    tv2.to_cArray(v2);

    T e0[3], e1[3];  // Triangle edge vectors starting at v0

    e0[0] = v1[0] - v0[0];
    e0[1] = v1[1] - v0[1];
    e0[2] = v1[2] - v0[2];

    e1[0] = v2[0] - v0[0];
    e1[1] = v2[1] - v0[1];
    e1[2] = v2[2] - v0[2];

    T normal[3];  // Triangle normal

    normal[0] = e0[1] * e1[2] - e0[2] * e1[1];
    normal[1] = e0[2] * e1[0] - e0[0] * e1[2];
    normal[2] = e0[0] * e1[1] - e0[1] * e1[0];

    T direction[3], p1[3], p2[3];

    point1.to_cArray(p1);
    point2.to_cArray(p2);

    direction[0] = p2[0] - p1[0];
    direction[1] = p2[1] - p1[1];
    direction[2] = p2[2] - p1[2];

    T denom = direction[0] * normal[0] + direction[1] * normal[1] + direction[2] * normal[2];

    if (util::isZero(denom, eps1))
        return false;  // The segment belongs to the plane or it is parallel to it
                       //   and does not intersect it.

    T num = (v0[0] - p1[0]) * normal[0] + (v0[1] - p1[1]) * normal[1] + (v0[2] - p1[2]) * normal[2];

    T t = num / denom;

    // The function segmentIntersectsTriangle is essentially to verify crossings
    //   through the surface from inside to outside or vice versa. When
    //   one of the end-points of a segment is right on top of the surface
    //   this creates an awkward ambiguity. To remove this ambiguity, the
    //   triangle is slightly shifted in direction of its normal vector.
    if ((util::isZero(t, eps1)) || (util::isOne(t, eps1))) {
        T norm_normal =
            std::sqrt(util::sqr(normal[0]) + util::sqr(normal[1]) + util::sqr(normal[2]));
        normal[0] /= norm_normal;
        normal[1] /= norm_normal;
        normal[2] /= norm_normal;
        // Shift the vertex, and recompute everything that needs to be
        //   recomputed. v1 and v2 are not used any more, so there is no
        //   need to recompute them.
        T tmp = (T)2 * eps1;
        v0[0] += tmp * normal[0];
        v0[1] += tmp * normal[1];
        v0[2] += tmp * normal[2];
        num =
            (v0[0] - p1[0]) * normal[0] + (v0[1] - p1[1]) * normal[1] + (v0[2] - p1[2]) * normal[2];
        t = num / denom;
    }

    if (util::lessThan(t, (T)0, eps1) || util::greaterThan(t, (T)1, eps1))
        return false;

    T intersection[3];

    intersection[0] = p1[0] + direction[0] * t;  // Intersection point with the plane
    intersection[1] = p1[1] + direction[1] * t;
    intersection[2] = p1[2] + direction[2] * t;

    T a[2][2];
    a[0][0] = e0[0] * e0[0] + e0[1] * e0[1] + e0[2] * e0[2];
    a[0][1] = e0[0] * e1[0] + e0[1] * e1[1] + e0[2] * e1[2];
    a[1][0] = a[0][1];
    a[1][1] = e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2];

    T tmp[3];
    tmp[0] = intersection[0] - v0[0];
    tmp[1] = intersection[1] - v0[1];
    tmp[2] = intersection[2] - v0[2];

    T b[2];
    b[0] = tmp[0] * e0[0] + tmp[1] * e0[1] + tmp[2] * e0[2];
    b[1] = tmp[0] * e1[0] + tmp[1] * e1[1] + tmp[2] * e1[2];

    T det = a[0][0] * a[1][1] - a[0][1] * a[1][0];

    T u = (a[1][1] * b[0] - a[0][1] * b[1]) / det;
    T v = (a[0][0] * b[1] - a[1][0] * b[0]) / det;

    T upv = u + v;

    if (util::lessThan(u, (T)0, eps1) || util::greaterThan(u, (T)1, eps1)
        || util::lessThan(v, (T)0, eps1) || util::greaterThan(v, (T)1, eps1)
        || util::greaterThan(upv, (T)1, eps1))
        return false;  // The point does not belong to the triangle

    return true;
}

template <typename T>
void TriangularSurfaceMesh<T>::distanceToEdgeLine(
    Array<T, 3> const &point, plint iTriangle, plint whichEdge, T &distance,
    bool &intersectionIsInside) const
{
    PLB_ASSERT(
        iTriangle >= 0 && iTriangle < numTriangles
        && (whichEdge == 0 || whichEdge == 1 || whichEdge == 2));

    plint iVertex1 = edges()[3 * iTriangle + ((whichEdge == 0) ? 2 : whichEdge - 1)].pv;
    plint iVertex2 = edges()[3 * iTriangle + whichEdge].pv;
    Array<T, 3> vertex1 = getVertex(iVertex1);
    Array<T, 3> vertex2 = getVertex(iVertex2);
    Array<T, 3> pv = point - vertex1;
    Array<T, 3> e = vertex2 - vertex1;

    T u = dot(pv, e) / dot(e, e);
    Array<T, 3> x = vertex1 + u * e;
    distance = norm(point - x);

    intersectionIsInside = util::greaterEqual(u, (T)0, eps1) && util::lessEqual(u, (T)1, eps1);
}

template <typename T>
void TriangularSurfaceMesh<T>::distanceToTrianglePlane(
    Array<T, 3> const &point, plint iTriangle, T &distance, bool &intersectionIsInside,
    bool &pointIsBehind) const
{
    Array<T, 3> normal = computeTriangleNormal(iTriangle);
    Array<T, 3> intersection((T)0, (T)0, (T)0);
    int flag = 2;  // Line.
    intersectionIsInside =
        pointOnTriangle(point, point + normal, flag, iTriangle, intersection, normal, distance)
        == 1;
    T projection = dot(normal, point - intersection);

    pointIsBehind = util::lessEqual(projection, (T)0, eps1);
}

template <typename T>
void TriangularSurfaceMesh<T>::distanceToTriangle(
    Array<T, 3> const &point, plint iTriangle, T &distance, bool &pointIsBehind) const
{
    bool intersectionIsInside;
    // First possibility: The projection of the point is inside the triangle.
    distanceToTrianglePlane(point, iTriangle, distance, intersectionIsInside, pointIsBehind);
    if (intersectionIsInside) {
        return;
    }
    // Second possibility: The projection of the point is inside an edge.
    bool intersectsWithEdge = false;
    for (plint iEdge = 0; iEdge <= 2; ++iEdge) {
        T newDistance;
        distanceToEdgeLine(point, iTriangle, iEdge, newDistance, intersectionIsInside);
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
    T d0Sqr = normSqr(point - getVertex(iTriangle, 0));
    T d1Sqr = normSqr(point - getVertex(iTriangle, 1));
    T d2Sqr = normSqr(point - getVertex(iTriangle, 2));
    T minDistSqr = std::min(d0Sqr, std::min(d1Sqr, d2Sqr));
    distance = std::sqrt(minDistSqr);
}

template <typename T>
void TriangularSurfaceMesh<T>::reverseOrientation()
{
    // Rearrange the emanating edge list.
    for (plint iVertex = 0; iVertex < numVertices; iVertex++)
        if (isInteriorVertex(iVertex)) {
            emanatingEdges()[iVertex] = changeEdgeId(edges()[emanatingEdges()[iVertex]].ne);
        } else {
            plint ee = emanatingEdges()[iVertex];
            plint pe = prev(ee);
            plint e = pe;
            if ((e = edges()[e].ne) < 0) {
                emanatingEdges()[iVertex] = changeEdgeId(pe);
                e = ee;
            }
            while (e != ee) {
                e = pe = prev(e);
                if ((e = edges()[e].ne) < 0) {
                    emanatingEdges()[iVertex] = changeEdgeId(pe);
                    e = ee;
                }
            }
        }

    // Rearrange the neighboring edge list.
    std::vector<Edge> tmp = edges();
    for (plint iEdge = 0; iEdge < 3 * numTriangles; iEdge++)
        if (edges()[iEdge].ne < 0) {
            tmp[changeEdgeId(iEdge)].ne = -1;
        } else {
            tmp[changeEdgeId(iEdge)].ne = changeEdgeId(edges()[iEdge].ne);
        }

    edges().swap(tmp);

    // Rearrange the pointing vertex list.
    for (plint iTriangle = 0; iTriangle < numTriangles; iTriangle++) {
        plint id0 = getVertexId(iTriangle, 0);
        plint id1 = getVertexId(iTriangle, 1);
        plint id2 = getVertexId(iTriangle, 2);

        edges()[3 * iTriangle].pv = id2;
        edges()[3 * iTriangle + 1].pv = id1;
        edges()[3 * iTriangle + 2].pv = id0;
    }
}

template <typename T>
std::vector<Lid> TriangularSurfaceMesh<T>::closeHoles()
{
    std::vector<std::vector<plint> > holes = detectHoles();
    std::vector<Lid> lids(holes.size());
    for (pluint iHole = 0; iHole < holes.size(); ++iHole) {
        lids[iHole] = closeHole(holes[iHole]);
    }
    return lids;
}

template <typename T>
void TriangularSurfaceMesh<T>::avoidIntegerPositions()
{
    for (plint iVertex = 0; iVertex < getNumVertices(); ++iVertex) {
        avoidIntegerPosition(iVertex);
    }
}

template <typename T>
void TriangularSurfaceMesh<T>::avoidIntegerPosition(plint iVertex)
{
    Array<T, 3> &vertex = getVertex(iVertex);
    Array<T, 3> normal = computeVertexNormal(iVertex);
    if (vertex[0] - (plint)vertex[0] < 1.e-12 || vertex[1] - (plint)vertex[1] < 1.e-12
        || vertex[2] - (plint)vertex[2] < 1.e-12)
    {
        vertex += (T)1.e-12 * normal;
    }
}

template <typename T>
void TriangularSurfaceMesh<T>::inflate(T amount)
{
    for (plint iVertex = 0; iVertex < getNumVertices(); ++iVertex) {
        Array<T, 3> &vertex = getVertex(iVertex);
        Array<T, 3> normal = computeVertexNormal(iVertex);
        vertex += amount * normal;
    }
}

template <typename T>
inline plint TriangularSurfaceMesh<T>::prev(plint iEdge) const
{
    PLB_ASSERT(iEdge >= 0 && iEdge < 3 * numTriangles);
    return (iEdge % 3 == 0) ? iEdge + 2 : iEdge - 1;
}

template <typename T>
inline plint TriangularSurfaceMesh<T>::next(plint iEdge) const
{
    PLB_ASSERT(iEdge >= 0 && iEdge < 3 * numTriangles);
    return (iEdge % 3 == 2) ? iEdge - 2 : iEdge + 1;
}

template <typename T>
inline plint TriangularSurfaceMesh<T>::changeEdgeId(plint iEdge) const
{
    plint iTriangle = iEdge / 3;
    plint mod = iEdge % 3;
    return (mod == 0) ? 3 * iTriangle + 2 : ((mod == 1) ? iEdge : 3 * iTriangle);
}

template <typename T>
std::vector<std::vector<plint> > TriangularSurfaceMesh<T>::detectHoles()
{
    std::vector<std::vector<plint> > holes;
    std::vector<bool> check(getNumVertices());
    // Whenever a vertex is identified as part of a given hole boundary, it is
    //   opted out in the "check" array to avoid that it is mistakingly
    //   identified later on as the starting vertex for a new hole.
    std::fill(check.begin(), check.end(), false);
    // Check all vertices ...
    for (plint iVertex = 0; iVertex < numVertices; ++iVertex) {
        //                               ... unless they've been previously opted out.
        if (isBoundaryVertex(iVertex) && !check[iVertex]) {
            check[iVertex] = true;
            std::vector<plint> newHole;
            newHole.push_back(iVertex);
            Edge nextEdge = edges()[emanatingEdges()[iVertex]];
            // Loop around the hole, until we're back to the first vertex.
            while (nextEdge.pv != iVertex) {
                plint nextVertex = nextEdge.pv;
                check[nextVertex] = true;
                newHole.push_back(nextVertex);
                nextEdge = edges()[emanatingEdges()[nextVertex]];
            }
            holes.push_back(newHole);
        }
    }
    return holes;
}

template <typename T>
Lid TriangularSurfaceMesh<T>::closeHole(std::vector<plint> const &hole)
{
    plint numHoleVertices = (plint)hole.size();

    // 1. Create a new vertex, corresponding to the barycenter.
    Array<T, 3> baryCenter;
    baryCenter.resetToZero();
    for (plint iHole = 0; iHole < numHoleVertices; ++iHole) {
        plint iVertex = hole[iHole];
        baryCenter += getVertex(iVertex);
    }
    baryCenter /= (T)numHoleVertices;
    vertices().push_back(baryCenter);
    plint iCenter = numVertices;
    ++numVertices;  // Don't forget to update this class variable.

    // 2. Create the new triangles.
    plint edgeIndex = edges().size() - 1;
    plint firstAddedEdge = edgeIndex + 1;
    // There is exactly one new triangle for each vertex along the hole.
    for (plint iTriangle = 0; iTriangle < numHoleVertices; ++iTriangle) {
        plint iVertex = hole[iTriangle];
        plint nextVertex = hole[iTriangle == numHoleVertices - 1 ? 0 : iTriangle + 1];

        // a) Add the new edge which is horizontal to the hole boundary.
        Edge counterEdge;
        ++edgeIndex;
        counterEdge.pv = iVertex;
        counterEdge.ne = emanatingEdges()[iVertex];
        edges().push_back(counterEdge);
        // Convert this former boundary edge to a bulk edge.
        edges()[counterEdge.ne].ne = edgeIndex;

        // b) Add the edge that goes from a boundary vertex to the barycenter.
        Edge firstEdge;
        ++edgeIndex;
        firstEdge.pv = iCenter;
        if (iTriangle > 0) {
            firstEdge.ne = edgeIndex - 2;
            // Convert this former boundary edge to a bulk edge.
            edges()[firstEdge.ne].ne = edgeIndex;
        }
        edges().push_back(firstEdge);

        // c) Add the edge that goes from the barycenter to a boundary vertex.
        Edge secondEdge;
        ++edgeIndex;
        secondEdge.pv = nextVertex;
        if (iTriangle == numHoleVertices - 1) {
            secondEdge.ne = firstAddedEdge + 1;
            // Convert this former boundary edge to a bulk edge.
            edges()[secondEdge.ne].ne = edgeIndex;
        }
        edges().push_back(secondEdge);

        ++numTriangles;  // Don't forget to update this class variable.
    }
    // For the newly created barycenter, define the emanating edge to
    //   be the one corresponding to the first added triangle.
    emanatingEdges().push_back(firstAddedEdge + 2);

    // 3. Make sure this will not become a degenerate case.
    avoidIntegerPosition(iCenter);

    // 4. Keep track of the newly created triangles with help of a "Lid" structure.
    Lid lid;
    lid.firstTriangle = firstAddedEdge / 3;
    lid.numTriangles = numHoleVertices;
    lid.boundaryVertices = hole;
    lid.centerVertex = iCenter;
    return lid;
}

template <typename T>
void toLatticeUnits(
    TriangularSurfaceMesh<T> &mesh, plint resolution, plint referenceDirection,
    Array<T, 3> &location, T &dx)
{
    PLB_ASSERT(referenceDirection >= 0 && referenceDirection <= 2);
    Array<T, 2> xRange, yRange, zRange;
    mesh.computeBoundingBox(xRange, yRange, zRange);
    T deltaX = T();
    switch (referenceDirection) {
    case 0:
        deltaX = xRange[1] - xRange[0];
        break;
    case 1:
        deltaX = yRange[1] - yRange[0];
        break;
    case 2:
        deltaX = zRange[1] - zRange[0];
        break;
    }

    T scalingFactor = (T)(resolution) / deltaX;

    // Transform mesh into the frame of coordinate of the multiScalarField.
    location = Array<T, 3>(xRange[0], yRange[0], zRange[0]);
    mesh.translate(Array<T, 3>(-location));
    mesh.scale(scalingFactor);
    dx = 1. / scalingFactor;
}

/* ******* Lid operations ************************************************** */

template <typename T>
Array<T, 3> computeBaryCenter(TriangularSurfaceMesh<T> const &mesh, Lid const &lid)
{
    Array<T, 3> baryCenter;
    baryCenter.resetToZero();
    for (plint iBoundary = 0; iBoundary < (plint)lid.boundaryVertices.size(); ++iBoundary) {
        plint iVertex = lid.boundaryVertices[iBoundary];
        baryCenter += mesh.getVertex(iVertex);
    }
    baryCenter /= (T)lid.boundaryVertices.size();
    return baryCenter;
}

template <typename T>
Array<T, 3> computeGeometricCenter(TriangularSurfaceMesh<T> const &mesh, Lid const &lid)
{
    typedef typename TriangularSurfaceMesh<T>::Edge Edge;
    Array<T, 3> center;
    center.resetToZero();
    T circumference = T();
    for (plint iBoundary = 0; iBoundary < (plint)lid.boundaryVertices.size(); ++iBoundary) {
        plint iVertex = lid.boundaryVertices[iBoundary];
        Edge edge = mesh.edges()[mesh.emanatingEdges()[iVertex]];
        plint iNextVertex = edge.pv;
        Array<T, 3> v1 = mesh.getVertex(iVertex);
        Array<T, 3> v2 = mesh.getVertex(iNextVertex);
        T l = norm(v2 - v1);
        center += (v1 + v2) * l;
        circumference += l;
    }
    center /= (T)2. * circumference;
    return center;
}

template <typename T>
T computeGeometricRadius(TriangularSurfaceMesh<T> const &mesh, Lid const &lid)
{
    typedef typename TriangularSurfaceMesh<T>::Edge Edge;
    Array<T, 3> center = computeGeometricCenter(mesh, lid);
    T radius = T();
    T circumference = T();
    for (plint iBoundary = 0; iBoundary < (plint)lid.boundaryVertices.size(); ++iBoundary) {
        plint iVertex = lid.boundaryVertices[iBoundary];
        Edge edge = mesh.edges()[mesh.emanatingEdges()[iVertex]];
        plint iNextVertex = edge.pv;
        Array<T, 3> v1 = mesh.getVertex(iVertex);
        Array<T, 3> v2 = mesh.getVertex(iNextVertex);
        T l = norm(v2 - v1);
        radius += norm(v1 - center) * l;
        circumference += l;
    }
    radius /= circumference;
    return radius;
}

template <typename T>
void computeBoundingBox(
    TriangularSurfaceMesh<T> const &mesh, Lid const &lid, Array<T, 2> &xLim, Array<T, 2> &yLim,
    Array<T, 2> &zLim)
{
    xLim[0] = yLim[0] = zLim[0] = std::numeric_limits<T>::max();
    xLim[1] = yLim[1] = zLim[1] = -std::numeric_limits<T>::max();
    for (plint iBoundary = 0; iBoundary < (plint)lid.boundaryVertices.size(); ++iBoundary) {
        plint iVertex = lid.boundaryVertices[iBoundary];
        Array<T, 3> const &vertex = mesh.getVertex(iVertex);
        xLim[0] = std::min(xLim[0], vertex[0]);
        xLim[1] = std::max(xLim[1], vertex[0]);
        yLim[0] = std::min(yLim[0], vertex[1]);
        yLim[1] = std::max(yLim[1], vertex[1]);
        zLim[0] = std::min(zLim[0], vertex[2]);
        zLim[1] = std::max(zLim[1], vertex[2]);
    }
}

template <typename T>
T computeInnerRadius(TriangularSurfaceMesh<T> const &mesh, Lid const &lid)
{
    if (lid.boundaryVertices.empty()) {
        return T();
    } else {
        T innerRadius = mesh.computeEdgeLength(lid.centerVertex, lid.boundaryVertices[0]);
        for (plint iBoundary = 1; iBoundary < (plint)lid.boundaryVertices.size(); ++iBoundary) {
            T nextLength =
                mesh.computeEdgeLength(lid.centerVertex, lid.boundaryVertices[iBoundary]);
            innerRadius = std::min(innerRadius, nextLength);
        }
        return innerRadius;
    }
}

template <typename T>
T computeOuterRadius(TriangularSurfaceMesh<T> const &mesh, Lid const &lid)
{
    if (lid.boundaryVertices.empty()) {
        return T();
    } else {
        T innerRadius = mesh.computeEdgeLength(lid.centerVertex, lid.boundaryVertices[0]);
        for (plint iBoundary = 1; iBoundary < (plint)lid.boundaryVertices.size(); ++iBoundary) {
            T nextLength =
                mesh.computeEdgeLength(lid.centerVertex, lid.boundaryVertices[iBoundary]);
            innerRadius = std::max(innerRadius, nextLength);
        }
        return innerRadius;
    }
}

template <typename T>
T computeArea(TriangularSurfaceMesh<T> const &mesh, Lid const &lid)
{
    T area = T();
    plint numVertices = (plint)lid.boundaryVertices.size();
    plint centerVertex = lid.centerVertex;
    for (plint iBoundary = 0; iBoundary < numVertices; ++iBoundary) {
        plint vertex1 = lid.boundaryVertices[iBoundary];
        plint vertex2 = lid.boundaryVertices[(iBoundary + 1) % numVertices];
        area += mesh.computeTriangleArea(centerVertex, vertex1, vertex2);
    }
    return area;
}

template <typename T>
Array<T, 3> computeNormal(TriangularSurfaceMesh<T> const &mesh, Lid const &lid)
{
    return mesh.computeVertexNormal(lid.centerVertex);
}

template <typename T>
void reCenter(TriangularSurfaceMesh<T> &mesh, Lid const &lid)
{
    Array<T, 3> newCenter = computeBaryCenter(mesh, lid);
    mesh.replaceVertex(lid.centerVertex, newCenter);
}

template <typename T>
void TriangularSurfaceMesh<T>::writeHTML(std::string fname)
{
    writeHTML(fname, "Palabos html output", T(), Array<T, 3>((T)0., (T)0., (T)0.));
}

// TODO: title unused could we remove it?
template <typename T>
void TriangularSurfaceMesh<T>::writeHTML(
    std::string fname, std::string title, T phys_dx, Array<T, 3> phys_location)
{
    if (!global::mpi().isMainProcessor()) {
        return;
    }
    if (numTriangles == 0) {
        return;
    }
    std::vector<Array<T, 3> > posVect(getNumVertices());
    T maxLimit = std::numeric_limits<T>::max();
    T minLimit = -std::numeric_limits<T>::max();
    Array<T, 3> minPos(maxLimit, maxLimit, maxLimit), maxPos(minLimit, minLimit, minLimit);
    for (plint i = 0; i < getNumVertices(); ++i) {
        posVect[i] = (*vertexList)[i];
        if (posVect[i][0] < minPos[0])
            minPos[0] = posVect[i][0];
        if (posVect[i][1] < minPos[1])
            minPos[1] = posVect[i][1];
        if (posVect[i][2] < minPos[2])
            minPos[2] = posVect[i][2];
        if (posVect[i][0] > maxPos[0])
            maxPos[0] = posVect[i][0];
        if (posVect[i][1] > maxPos[1])
            maxPos[1] = posVect[i][1];
        if (posVect[i][2] > maxPos[2])
            maxPos[2] = posVect[i][2];
    }

    T diffx = maxPos[0] - minPos[0];
    T diffy = maxPos[1] - minPos[1];
    T diffz = maxPos[2] - minPos[2];
    T diff = std::max(diffx, std::max(diffy, diffz));
    T dx_nonDim = (T)2. / diff;
    Array<T, 3> offset_nonDim(
        -((T)1.0 + minPos[0] / diff), -((T)1.0 + minPos[1] / diff), -((T)1.0 + minPos[2] / diff));
    for (plint i = 0; i < getNumVertices(); ++i) {
        posVect[i] *= dx_nonDim;
        posVect[i] += offset_nonDim;
    }
    T dx_ndToPhys = phys_dx / dx_nonDim;
    Array<T, 3> location_ndToPhys = phys_location - dx_ndToPhys * offset_nonDim;

    std::ofstream ofile(fname.c_str());
    ofile << std::fixed << std::setprecision(5);
    ofile << "<!DOCTYPE HTML>\n";
    ofile << "<html>\n";
    ofile << "    <head>\n";
    ofile << "        <title>"
          << "Palabos Geometry"
          << "</title>\n";
    ofile << "        <link rel=\"stylesheet\" type=\"text/css\" "
             "href=\"http://www.x3dom.org/download/x3dom.css\">\n";
    ofile << "        <script type=\"text/javascript\" "
             "src=\"http://www.x3dom.org/download/1.4/x3dom-full.js\">\n";
    ofile << "        </script>\n";

    ofile << "        <script type=\"text/javascript\">\n";
    ofile << "            var scale = " << std::setprecision(5) << dx_ndToPhys << ";\n";
    ofile << "            var offset = new Array(" << std::setprecision(5) << location_ndToPhys[0]
          << "," << location_ndToPhys[1] << "," << location_ndToPhys[2] << ");\n";
    ofile << "        </script>\n";
    ofile << "    </head>\n";
    ofile << "    <x3d style=\"width: 100%; height :100%\">\n";
    ofile << "        <scene>\n";
    ofile << "            <shape>\n";
    ofile << "            <indexedtriangleset index=\"";

    for (plint iTriangle = 0; iTriangle < getNumTriangles(); ++iTriangle) {
        plint i0 = getVertexId(iTriangle, 0);
        plint i1 = getVertexId(iTriangle, 1);
        plint i2 = getVertexId(iTriangle, 2);
        ofile << i0 << " " << i1 << " " << i2 << " ";
    }
    ofile << "\" solid=\"false\">\n";
    ofile << "            <coordinate point=\"";
    for (plint iVertex = 0; iVertex < getNumVertices(); ++iVertex) {
        ofile << posVect[iVertex][0] << " " << posVect[iVertex][1] << " " << posVect[iVertex][2];
        if (iVertex < getNumVertices() - 1) {
            ofile << " ";
        }
    }
    ofile << "\"> </coordinate>\n";
    ofile << "            <normal vector=\"";
    for (plint iVertex = 0; iVertex < getNumVertices(); ++iVertex) {
        Array<T, 3> normal = computeVertexNormal(iVertex);
        ofile << normal[0] << " " << normal[1] << " " << normal[2];
        if (iVertex < getNumVertices() - 1) {
            ofile << " ";
        }
    }
    ofile << "\"> </normal>\n";
    ofile << "            </indexedtriangleset>\n";
    ofile << "            </shape>\n";
    ofile << "            <viewpoint position=\'0.0 0.0 4.0\' orientation=\'0.0 0.0 0.0 "
             "0.0\'></viewpoint>\n";
    ofile << "        </scene>\n";
    ofile << "    </x3d>\n";
    ofile << "</html>\n";
}

template <typename T>
void TriangularSurfaceMesh<T>::writeX3D(std::string fname)
{
    writeX3D(fname, "Palabos X3D output", T(), Array<T, 3>((T)0., (T)0., (T)0.));
}

// TODO: title, phys_dx and phys_location unused could we remove them?
template <typename T>
void TriangularSurfaceMesh<T>::writeX3D(
    std::string fname, std::string title, T phys_dx, Array<T, 3> phys_location)
{
    if (!global::mpi().isMainProcessor()) {
        return;
    }
    if (numTriangles == 0) {
        return;
    }
    std::vector<Array<T, 3> > posVect(getNumVertices());
    T maxLimit = std::numeric_limits<T>::max();
    T minLimit = -std::numeric_limits<T>::max();
    Array<T, 3> minPos(maxLimit, maxLimit, maxLimit), maxPos(minLimit, minLimit, minLimit);
    for (plint i = 0; i < getNumVertices(); ++i) {
        posVect[i] = (*vertexList)[i];
        if (posVect[i][0] < minPos[0])
            minPos[0] = posVect[i][0];
        if (posVect[i][1] < minPos[1])
            minPos[1] = posVect[i][1];
        if (posVect[i][2] < minPos[2])
            minPos[2] = posVect[i][2];
        if (posVect[i][0] > maxPos[0])
            maxPos[0] = posVect[i][0];
        if (posVect[i][1] > maxPos[1])
            maxPos[1] = posVect[i][1];
        if (posVect[i][2] > maxPos[2])
            maxPos[2] = posVect[i][2];
    }

    T diffx = maxPos[0] - minPos[0];
    T diffy = maxPos[1] - minPos[1];
    T diffz = maxPos[2] - minPos[2];
    T diff = std::max(diffx, std::max(diffy, diffz));
    T dx_nonDim = (T)2. / diff;
    Array<T, 3> offset_nonDim(
        -((T)1.0 + minPos[0] / diff), -((T)1.0 + minPos[1] / diff), -((T)1.0 + minPos[2] / diff));
    for (plint i = 0; i < getNumVertices(); ++i) {
        posVect[i] *= dx_nonDim;
        posVect[i] += offset_nonDim;
    }

    std::ofstream ofile(fname.c_str());
    ofile << std::fixed << std::setprecision(5);
    ofile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    ofile << "<!DOCTYPE X3D PUBLIC \"ISO//Web3D//DTD X3D 3.0//EN\" "
             "\"http://www.web3d.org/specifications/x3d-3.0.dtd\">\n";
    ofile << "<X3D version=\"3.0\" profile=\"Immersive\" "
             "xmlns:xsd=\"http://www.w3.org/2001/XMLSchema-instance\" "
             "xsd:noNamespaceSchemaLocation=\"http://www.web3d.org/specifications/x3d-3.0.xsd\">\n";
    ofile << "        <Scene>\n";
    ofile << "            <Shape>\n";
    ofile << "            <IndexedFaceSet coordIndex=\"";

    for (plint iTriangle = 0; iTriangle < getNumTriangles(); ++iTriangle) {
        plint i0 = getVertexId(iTriangle, 0);
        plint i1 = getVertexId(iTriangle, 1);
        plint i2 = getVertexId(iTriangle, 2);
        ofile << i0 << " " << i1 << " " << i2 << " -1 ";
    }
    ofile << "\" solid=\"false\">\n";
    ofile << "            <Coordinate point=\"";
    for (plint iVertex = 0; iVertex < getNumVertices(); ++iVertex) {
        ofile << posVect[iVertex][0] << " " << posVect[iVertex][1] << " " << posVect[iVertex][2];
        if (iVertex < getNumVertices() - 1) {
            ofile << " ";
        }
    }
    ofile << "\"/>\n";
    ofile << "            <Normal vector=\"";
    for (plint iVertex = 0; iVertex < getNumVertices(); ++iVertex) {
        Array<T, 3> normal = computeVertexNormal(iVertex);
        ofile << normal[0] << " " << normal[1] << " " << normal[2];
        if (iVertex < getNumVertices() - 1) {
            ofile << " ";
        }
    }
    ofile << "\"/>\n";
    ofile << "            </IndexedFaceSet>\n";
    ofile << "            </Shape>\n";
    ofile << "        </Scene>\n";
    ofile << "    </X3D>\n";
}

}  // namespace plb

#endif  // TRIANGULAR_SURFACE_MESH_HH
