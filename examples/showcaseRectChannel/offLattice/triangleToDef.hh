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

#ifndef TRIANGLE_TO_DEF_HH
#define TRIANGLE_TO_DEF_HH

#include <cstdlib>
#include <limits>
#include <queue>
#include <utility>

#include "core/geometry3D.h"
#include "core/util.h"
#include "offLattice/triangleToDef.h"

namespace plb {

template <typename T>
void TriangleToDef<T>::vsAdd(Array<T, 3> const &coord, plint &index, plint &count)
{
    VsNodeIt it = vertexSet.find(VertexSetNode(-1, &coord));
    if (it == vertexSet.end()) {
        index = count;
        ++count;
        VertexSetNode newNode(index, &coord);
        vertexSet.insert(newNode);
    } else {
        index = it->i;
    }
}

template <typename T>
void TriangleToDef<T>::vsOrder()
{
    VsNodeConstIt it = vertexSet.begin();
    for (; it != vertexSet.end(); ++it) {
        vertexList[it->i] = *(it->vertex);
    }
}

template <typename T>
plint TriangleToDef<T>::searchEdgeList(
    std::vector<typename TriangleToDef<T>::EdgeListNode> const &edgeList, plint maxv) const
{
    for (pluint iList = 0; iList < edgeList.size(); ++iList) {
        if (edgeList[iList].maxv == maxv) {
            return iList;
        }
    }
    return -1;
}

template <typename T>
typename TriangleToDef<T>::BvmNodeIt TriangleToDef<T>::bvmAdd(plint id)
{
    BvmNodeIt it = boundaryVertexMap.find(id);
    if (it == boundaryVertexMap.end()) {
        it = boundaryVertexMap.insert(std::make_pair(id, BoundaryVertexMapNode())).first;
    }
    return it;
}

template <typename T>
bool TriangleToDef<T>::bvmCheck() const
{
    BvmNodeConstIt it = boundaryVertexMap.begin();
    for (; it != boundaryVertexMap.end(); ++it) {
        BoundaryVertexMapNode const &node = it->second;
        if (node.v1 == -1 || node.v2 == -1 || node.t1 == -1 || node.t2 == -1 || node.counter != 2) {
#ifdef PLB_DEBUG
            Array<T, 3> v = vertexList[it->first];
            pcout << "There is a problem with the boundary of the triangular surface mesh."
                  << std::endl;
            pcout << "The problematic vertex is: [" << v[0] << ", " << v[1] << ", " << v[2] << "]"
                  << std::endl;
#endif
            return true;
        }
    }
    return false;
}

template <typename T>
void TriangleToDef<T>::bvmLabel()
{
    plint vertex2, triangle2, localEdge2;
    BvmNodeConstIt it = boundaryVertexMap.begin();
    for (; it != boundaryVertexMap.end(); ++it) {
        plint vertex = it->first;  // Index of the boundary vertex
        BoundaryVertexMapNode const &node = it->second;
        vertex2 = node.v2;    // Index of the boundary vertex to which
                              // the emanating boundary edge from the
                              // boundary vertex `vertex' points to
        triangle2 = node.t2;  // Index of the adjacent triangle to
                              // the boundary edge emanating from
                              // the boundary vertex `vertex'

        localEdge2 = -1;  // Local index (in the triangle t2) of the
                          // boundary edge that points from `vertex' to
                          // `vertex2'
        if (vertex == edgeList[3 * triangle2 + 2].pv && vertex2 == edgeList[3 * triangle2].pv) {
            localEdge2 = 0;
        } else if (
            vertex == edgeList[3 * triangle2].pv && vertex2 == edgeList[3 * triangle2 + 1].pv) {
            localEdge2 = 1;
        } else if (
            vertex == edgeList[3 * triangle2 + 1].pv && vertex2 == edgeList[3 * triangle2 + 2].pv) {
            localEdge2 = 2;
        } else {
            PLB_ASSERT(false);  // Problem with the boundary of the surface mesh.
        }

        if (emanatingEdgeList[vertex] == -1)
            emanatingEdgeList[vertex] = 3 * triangle2 + localEdge2;
        else {
            PLB_ASSERT(false);  // Problem with the boundary of the surface mesh.
        }
    }
}

template <typename T>
plint &TriangleToDef<T>::globalVertex(plint triangle, plint localVertex)
{
    return edgeList[3 * triangle + ((localVertex == 0) ? 2 : localVertex - 1)].pv;
}

template <typename T>
TriangleToDef<T>::TriangleToDef(std::vector<Triangle> const &triangles, T epsilon) :
    vertexSet(PositionLessThan3D<T, VertexSetNode>(epsilon))
{
    numTriangles = triangles.size();

    numVertices = uniqueVertices(triangles);
    vertexList.resize(numVertices);
    emanatingEdgeList.resize(numVertices);

    vsOrder();

    computePointingVertex();

    plint nbe = createEdgeTable();

    if (fixOrientation()) {
#ifdef PLB_DEBUG
        plint nbe_new = createEdgeTable();
#else
        (void)createEdgeTable();
#endif
        PLB_ASSERT(nbe == nbe_new);  // Problem with the surface mesh topology.
    }

    if (nbe != 0) {
        findBoundaryVertices();
    }

    computeNeighboringEdges();

    computeEmanatingEdges();
}

// Fix the orientation of all triangles, so that they are consistently oriented,
//   meaning that all triangle normals point either inwards or outwards.
template <typename T>
bool TriangleToDef<T>::fixOrientation()
{
    bool fixedSomething = false;
    std::queue<plint> trianglesToFixNeighbors;
    char *visitedTriangles = (char *)calloc(numTriangles, sizeof(char));
    for (plint iTriangle = 0; iTriangle < numTriangles; iTriangle++) {
        if (visitedTriangles[iTriangle] == 0) {
            trianglesToFixNeighbors.push(iTriangle);
            while (!trianglesToFixNeighbors.empty()) {
                plint triangle = trianglesToFixNeighbors.front();
                trianglesToFixNeighbors.pop();
                fixOrientationOfNeighbors(
                    triangle, trianglesToFixNeighbors, visitedTriangles, fixedSomething);
            }
        }
    }

    free(visitedTriangles);
    return fixedSomething;
}

// Fix the orientation of all neighboring triangles to a triangle of index iTriangle,
//   so that ther orientations are consistent to the one of the triangle iTriangle.
template <typename T>
void TriangleToDef<T>::fixOrientationOfNeighbors(
    plint iTriangle, std::queue<plint> &trianglesToFixNeighbors, char *visitedTriangles, bool &flag)
{
    visitedTriangles[iTriangle] = 1;

    plint ind[3] = {
        globalVertex(iTriangle, 0), globalVertex(iTriangle, 1), globalVertex(iTriangle, 2)};

    // Loop over local edges

    for (plint localEdge = 0; localEdge < 3; localEdge++) {
        plint ia = ind[localEdge];
        plint ib = ind[(localEdge == 2) ? 0 : (localEdge + 1)];

        plint min_i = std::min(ia, ib);
        plint max_i = std::max(ia, ib);
        plint jTriangle = -2;
        for (pluint iEdge = 0; iEdge < edgeTable[min_i].size(); iEdge++)
            if (edgeTable[min_i][iEdge].maxv == max_i) {
                plint t1 = edgeTable[min_i][iEdge].t1;
                plint t2 = edgeTable[min_i][iEdge].t2;

                if (iTriangle == t1) {
                    jTriangle = t2;
                } else if (iTriangle == t2) {
                    jTriangle = t1;
                } else {
                    PLB_ASSERT(false);  // Problem with the surface mesh edges.
                }

                break;
            }

        PLB_ASSERT(jTriangle != -2);  // Problem with the surface mesh.

        if (jTriangle != -1) {
            if (visitedTriangles[jTriangle] == 0) {
                plint j0 = globalVertex(jTriangle, 0);
                plint j1 = globalVertex(jTriangle, 1);
                plint j2 = globalVertex(jTriangle, 2);

                if (ia == j0 && ib == j1) {
                    std::swap(triangleIndices[jTriangle][0], triangleIndices[jTriangle][1]);
                    globalVertex(jTriangle, 0) = triangleIndices[jTriangle][0];
                    globalVertex(jTriangle, 1) = triangleIndices[jTriangle][1];
                    flag = true;
                } else if (ia == j1 && ib == j2) {
                    std::swap(triangleIndices[jTriangle][1], triangleIndices[jTriangle][2]);
                    globalVertex(jTriangle, 1) = triangleIndices[jTriangle][1];
                    globalVertex(jTriangle, 2) = triangleIndices[jTriangle][2];
                    flag = true;
                } else if (ia == j2 && ib == j0) {
                    std::swap(triangleIndices[jTriangle][2], triangleIndices[jTriangle][0]);
                    globalVertex(jTriangle, 2) = triangleIndices[jTriangle][2];
                    globalVertex(jTriangle, 0) = triangleIndices[jTriangle][0];
                    flag = true;
                } else if (
                    (ia == j1 && ib == j0) || (ia == j2 && ib == j1) || (ia == j0 && ib == j2)) {
                    flag = false;
                } else {
                    PLB_ASSERT(false);  // Problem with the surface mesh.
                }

                visitedTriangles[jTriangle] = 1;
                trianglesToFixNeighbors.push(jTriangle);
            }
        }
    }
}

// Make vertices unique through a unique labelling scheme.
//   Return the number of unique vertices.
template <typename T>
plint TriangleToDef<T>::uniqueVertices(std::vector<Triangle> const &triangles)
{
    edgeList.resize(3 * numTriangles);

    triangleIndices.resize(numTriangles);
    plint index = 0;
    plint count = 0;
    for (plint iTriangle = 0; iTriangle < numTriangles; ++iTriangle) {
        vsAdd(triangles[iTriangle][0], index, count);
        triangleIndices[iTriangle][0] = index;

        vsAdd(triangles[iTriangle][1], index, count);
        triangleIndices[iTriangle][1] = index;

        vsAdd(triangles[iTriangle][2], index, count);
        triangleIndices[iTriangle][2] = index;
    }
    return count;
}

template <typename T>
void TriangleToDef<T>::computePointingVertex()
{
    for (plint iTriangle = 0; iTriangle < numTriangles; ++iTriangle) {
        globalVertex(iTriangle, 0) = triangleIndices[iTriangle][0];
        globalVertex(iTriangle, 1) = triangleIndices[iTriangle][1];
        globalVertex(iTriangle, 2) = triangleIndices[iTriangle][2];
    }
}

/// Create hash-table of edges, distinguish boundary edges
/// and count boundary edges.
/// The purpose of this table is to identify the edges and
/// to separate between boundary and interior edges.
template <typename T>
plint TriangleToDef<T>::createEdgeTable()
{
    /* Create the edge hash table */
    edgeTable.resize(0);
    edgeTable.resize(vertexList.size());
    /* Local indices of edge endpoints */
    int epoint[3][2] = {{0, 1}, {1, 2}, {2, 0}};

    plint nbe = 0; /* Number of boundary edges */

    for (plint iTriangle = 0; iTriangle < numTriangles; ++iTriangle) {
        for (plint localEdge = 0; localEdge < 3; ++localEdge) {
            plint i0 = globalVertex(iTriangle, epoint[localEdge][0]);
            plint i1 = globalVertex(iTriangle, epoint[localEdge][1]);
            plint min_i = std::min(i0, i1);
            plint max_i = std::max(i0, i1);

            plint iList = searchEdgeList(edgeTable[min_i], max_i);
            if (iList == -1) {
                EdgeListNode newNode;
                newNode.maxv = max_i;
                newNode.t1 = iTriangle;
                newNode.t2 = -1;
                edgeTable[min_i].push_back(newNode);
                ++nbe;
            } else {
                EdgeListNode &node = edgeTable[min_i][iList];
                if (node.t2 == -1) {
                    node.t2 = iTriangle;
                    --nbe;
                } else {
                    PLB_ASSERT(false);  // The surface mesh contains non-manifold edges.
                }
            }
        }
    }
    return nbe;
}

// Create the map of boundary vertices with additional information on the
//   adjacent boundary edges and their orientation.
template <typename T>
void TriangleToDef<T>::findBoundaryVertices()
{
    BoundaryVertexMap boundaryVertexMap;
    BvmNodeIt node_i = boundaryVertexMap.end();
    BvmNodeIt node_j = boundaryVertexMap.end();
    plint jVertex, triangle;
    plint iVertex;
    int lock;
    for (iVertex = 0, lock = 1; iVertex < (plint)edgeTable.size(); ++iVertex, lock = 1) {
        for (pluint iEdge = 0; iEdge < edgeTable[iVertex].size(); ++iEdge) {
            if (edgeTable[iVertex][iEdge].t2 == -1) {
                if (lock) {
                    node_i = bvmAdd(iVertex);
                    lock = 0;
                }
                jVertex = edgeTable[iVertex][iEdge].maxv;
                triangle = edgeTable[iVertex][iEdge].t1;
                node_j = bvmAdd(jVertex);
                node_i->second.counter++;
                node_j->second.counter++;

                plint localVertex_i = -1;
                if (iVertex == globalVertex(triangle, 0)) {
                    localVertex_i = 0;
                } else if (iVertex == globalVertex(triangle, 1)) {
                    localVertex_i = 1;
                } else if (iVertex == globalVertex(triangle, 2)) {
                    localVertex_i = 2;
                } else {
                    PLB_ASSERT(false);  // Problem with the boundary of the surface mesh.
                }

                plint localVertex_j = -1;
                if (jVertex == globalVertex(triangle, 0)) {
                    localVertex_j = 0;
                } else if (jVertex == globalVertex(triangle, 1)) {
                    localVertex_j = 1;
                } else if (jVertex == globalVertex(triangle, 2)) {
                    localVertex_j = 2;
                } else {
                    PLB_ASSERT(false);  // Problem with the boundary of the surface mesh.
                }

                if ((localVertex_i == 0 && localVertex_j == 1)
                    || (localVertex_i == 1 && localVertex_j == 2)
                    || (localVertex_i == 2 && localVertex_j == 0))
                {
                    node_i->second.v2 = jVertex;
                    node_i->second.t2 = triangle;
                    node_j->second.v1 = iVertex;
                    node_j->second.t1 = triangle;
                } else if (
                    (localVertex_i == 1 && localVertex_j == 0)
                    || (localVertex_i == 2 && localVertex_j == 1)
                    || (localVertex_i == 0 && localVertex_j == 2))
                {
                    node_i->second.v1 = jVertex;
                    node_i->second.t1 = triangle;
                    node_j->second.v2 = iVertex;
                    node_j->second.t2 = triangle;
                }
            }
        }
    }
#ifdef PLB_DEBUG
    bool failure = bvmCheck();
#else
    (void)bvmCheck();
#endif
    PLB_ASSERT(!failure);  // Problem with the boundary of the surface mesh.
}

/// Use all previous information to fill in the missing fields in
///   the directed edge data.
template <typename T>
void TriangleToDef<T>::computeNeighboringEdges()
{
    for (plint iVertex = 0; iVertex < (plint)edgeTable.size(); ++iVertex) {
        for (plint iEdge = 0; iEdge < (plint)edgeTable[iVertex].size(); ++iEdge) {
            EdgeListNode &node = edgeTable[iVertex][iEdge];
            plint jVertex = node.maxv;
            plint triangle1 = node.t1;
            plint triangle2 = node.t2;
            plint localEdge1 = -1;

            plint id0 = globalVertex(triangle1, 0);
            plint id1 = globalVertex(triangle1, 1);
            plint id2 = globalVertex(triangle1, 2);

            if ((iVertex == id0 && jVertex == id1) || (iVertex == id1 && jVertex == id0)) {
                localEdge1 = 0;
            } else if ((iVertex == id1 && jVertex == id2) || (iVertex == id2 && jVertex == id1)) {
                localEdge1 = 1;
            } else if ((iVertex == id2 && jVertex == id0) || (iVertex == id0 && jVertex == id2)) {
                localEdge1 = 2;
            } else {
                PLB_ASSERT(false);  // Problem with the surface mesh connectivity.
            }

            if (triangle2 != -1) { /* Internal edge */
                plint va1 = globalVertex(triangle1, localEdge1);
                plint vb1 = globalVertex(triangle1, ((localEdge1 != 2) ? localEdge1 + 1 : 0));

                computeInternalNeighboringEdges(
                    iVertex, jVertex, triangle1, triangle2, localEdge1, va1, vb1);
            } else { /* Boundary edge */
                // By definition, on boundary edges, the neighboring
                //   edge has negative id equal to -1.
                edgeList[3 * triangle1 + localEdge1].ne = -1;
            }
        }
    }
}

template <typename T>
void TriangleToDef<T>::computeInternalNeighboringEdges(
    plint iVertex, plint jVertex, plint triangle1, plint triangle2, plint localEdge1, plint va1,
    plint vb1)
{
    plint localEdge2 = -1;

    plint id0 = globalVertex(triangle2, 0);
    plint id1 = globalVertex(triangle2, 1);
    plint id2 = globalVertex(triangle2, 2);

    if ((iVertex == id0 && jVertex == id1) || (iVertex == id1 && jVertex == id0)) {
        localEdge2 = 0;
    } else if ((iVertex == id1 && jVertex == id2) || (iVertex == id2 && jVertex == id1)) {
        localEdge2 = 1;
    } else if ((iVertex == id2 && jVertex == id0) || (iVertex == id0 && jVertex == id2)) {
        localEdge2 = 2;
    } else {
        PLB_ASSERT(false);  // Problem with the surface mesh connectivity.
    }

    plint va2 = globalVertex(triangle2, localEdge2);
    plint vb2 = globalVertex(triangle2, ((localEdge2 != 2) ? localEdge2 + 1 : 0));

    if (va1 != vb2 || vb1 != va2) {
        PLB_ASSERT(false);  // Problem with the surface mesh orientation.
    }

    edgeList[3 * triangle1 + localEdge1].ne = 3 * triangle2 + localEdge2;
    edgeList[3 * triangle2 + localEdge2].ne = 3 * triangle1 + localEdge1;
}

template <typename T>
void TriangleToDef<T>::computeEmanatingEdges()
{
    // Handle emanating edges.
    for (pluint iVertex = 0; iVertex < vertexList.size(); ++iVertex) {
        emanatingEdgeList[iVertex] = -1;
    }
    // For every boundary vertex choose as an emanating edge
    // the one that belongs to the boundary to facilitate
    // boundary operations.
    bvmLabel();

    for (plint iTriangle = 0; iTriangle < numTriangles; ++iTriangle) {
        for (plint localEdge = 0; localEdge < 3; ++localEdge) {
            plint &emanatingEdge = emanatingEdgeList[globalVertex(iTriangle, localEdge)];
            if (emanatingEdge == -1) {
                emanatingEdge = 3 * iTriangle + localEdge;
            }
        }
    }
}

template <typename T>
void TriangleToDef<T>::generateOnce(
    std::vector<Array<T, 3> > &vertexList_, std::vector<plint> &emanatingEdgeList_,
    std::vector<Edge> &edgeList_)
{
    vertexList.swap(vertexList_);
    emanatingEdgeList.swap(emanatingEdgeList_);
    edgeList.swap(edgeList_);
}

}  // namespace plb

#endif  // TRIANGLE_TO_DEF_HH
