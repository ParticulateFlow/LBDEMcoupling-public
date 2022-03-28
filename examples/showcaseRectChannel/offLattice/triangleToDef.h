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

#ifndef TRIANGLE_TO_DEF_H
#define TRIANGLE_TO_DEF_H

#include <map>
#include <queue>
#include <set>
#include <vector>

#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "offLattice/triangleSet.h"
#include "offLattice/triangularSurfaceMesh.h"

namespace plb {

template <typename T>
class TriangleToDef {
public:
    typedef typename TriangleSet<T>::Triangle Triangle;

private:  // This class should only be used by the function constructSurfaceMesh().
    TriangleToDef(std::vector<Triangle> const &triangles, T epsilon = getEpsilon<T>());
    void generateOnce(
        std::vector<Array<T, 3> > &vertexList_, std::vector<plint> &emanatingEdgeList_,
        std::vector<Edge> &edgeList_);

private:
    struct EdgeListNode {
        plint maxv;   /* Integer key value representing the maximum of the global
                         indices of the two vertices that constitute a mesh edge */
        plint t1, t2; /* Indices of the adjacent triangles to the specific
                         edge. The value -1 in t2 indicates that the edge has no
                         triangle neighbor and therefore is a boundary edge */
    };
    struct VertexSetNode {
        VertexSetNode(plint i_, Array<T, 3> const *vertex_) : i(i_), vertex(vertex_) { }
        Array<T, 3> const &getPosition() const
        {
            return *vertex;
        }
        plint i;                    // Global index of the vertex
        Array<T, 3> const *vertex;  // Pointer to vertex coordinates
    };
    struct BoundaryVertexMapNode {
        BoundaryVertexMapNode() : v1(-1), t1(-1), v2(-1), t2(-1), counter(0) { }
        plint v1; /* Global index of the vertex that the boundary edge
                     that points to the boundary vertex of index equal to
                     the map key starts from */
        plint t1; /* Index of the adjacent triangle of the boundary edge
                     that points to the boundary vertex of index equal to
                     the map key */

        plint v2; /* Global index of the vertex that the boundary edge
                     that starts from the boundary vertex of index equal
                     to the map key points to */
        plint t2; /* Index of the adjacent triangle of the boundary edge
                     that starts from the boundary vertex of index equal
                     to the map key */

        plint counter; /* Counter of boundary edges attached to the boundary
                          vertex of index equal to the map key. If `counter'
                          is not equal to two then there is a problem with
                          the mesh */
    };

    typedef std::set<VertexSetNode, PositionLessThan3D<T, VertexSetNode> > VertexSet;
    typedef typename VertexSet::iterator VsNodeIt;
    typedef typename VertexSet::const_iterator VsNodeConstIt;

    typedef std::map<plint, BoundaryVertexMapNode> BoundaryVertexMap;
    typedef typename BoundaryVertexMap::iterator BvmNodeIt;
    typedef typename BoundaryVertexMap::const_iterator BvmNodeConstIt;

private:
    void vsAdd(Array<T, 3> const &coord, plint &index, plint &count);
    void vsOrder();
    plint searchEdgeList(std::vector<EdgeListNode> const &edgeList, plint maxv) const;
    BvmNodeIt bvmAdd(plint id);
    bool bvmCheck() const;
    void bvmLabel();
    plint &globalVertex(plint triangle, plint localVertex);
    plint uniqueVertices(std::vector<Triangle> const &triangles);
    void computePointingVertex();
    plint createEdgeTable();
    void findBoundaryVertices();
    void computeNeighboringEdges();
    void computeInternalNeighboringEdges(
        plint iVertex, plint jVertex, plint triangle1, plint triangle2, plint localEdge1, plint va1,
        plint vb1);
    void computeEmanatingEdges();
    bool fixOrientation();
    void fixOrientationOfNeighbors(
        plint iTriangle, std::queue<plint> &trianglesToFixNeighbors, char *visitedTriangles,
        bool &flag);

private:
    std::vector<Array<plint, 3> > triangleIndices;
    VertexSet vertexSet;
    std::vector<std::vector<EdgeListNode> > edgeTable;
    BoundaryVertexMap boundaryVertexMap;

private:
    std::vector<Array<T, 3> > vertexList;
    std::vector<plint> emanatingEdgeList;
    std::vector<Edge> edgeList;
    plint numTriangles, numVertices;
    template <typename U>
    friend void constructSurfaceMesh(
        std::vector<typename TriangleToDef<U>::Triangle> const &triangles,
        std::vector<Array<U, 3> > &vertexList, std::vector<plint> &emanatingEdgeList,
        std::vector<Edge> &edgeList, U epsilon);
};

template <typename T>
void constructSurfaceMesh(
    std::vector<typename TriangleToDef<T>::Triangle> const &triangles,
    std::vector<Array<T, 3> > &vertexList, std::vector<plint> &emanatingEdgeList,
    std::vector<Edge> &edgeList, T epsilon = getEpsilon<T>())
{
    TriangleToDef<T>(triangles, epsilon).generateOnce(vertexList, emanatingEdgeList, edgeList);
}

}  // namespace plb

#endif  // TRIANGLE_TO_DEF_H
