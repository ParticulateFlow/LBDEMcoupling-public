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

#ifndef CONNECTED_TRIANGLE_UTIL_H
#define CONNECTED_TRIANGLE_UTIL_H

#include "core/geometry3D.h"
#include "offLattice/connectedTriangleMesh.h"

namespace plb {

/// Computes the area of all triangles, and stores them as a new triangle property.
/// If the triangle property already exists, it gets overwritten.
template <typename T>
void createTriangleAreaCache(ConnectedTriangleMesh<T> &mesh, std::string triangleProperty);

/// Computes the area of all vertices, and stores them as a new vertex property.
/// If the vertex property already exists, it gets overwritten.
template <typename T>
void createVertexAreaCache(ConnectedTriangleMesh<T> &mesh, std::string vertexProperty);

/// Extracts a part from a connected mesh (of any kind, could be the DEF one),
/// and create a RawConnectedTriangleMesh from it.
template <typename T>
RawConnectedTriangleMesh<T> extractConnectedPart(ConnectedTriangleMesh<T> &mesh, plint whichPart);

/// Returns a vector with the positions of all the unique vertices of a ConnectedTriangleMesh.
/// All parts of the mesh are included without distinction.
template <typename T>
std::vector<Array<T, 3> > getAllUniqueVertices(ConnectedTriangleMesh<T> &mesh);

template <typename T>
class MeshConnector {
public:
    MeshConnector(TriangleMesh<T> &triangleMesh_, T eps_ = getEpsilon<T>(DBL));
    MeshConnector(RawTriangleMesh<T> &rawTriangleMesh_);
    RawConnectedTriangleMesh<T> generateConnectedMesh() const;

private:
    void uniqueVertices();

private:
    // To unify duplicated vertices they need to be inserted in a std::set container.
    struct VertexSetNode {
        VertexSetNode(plint i_, Array<T, 3> const *vertex_) : i(i_), vertex(vertex_) { }
        Array<T, 3> const &getPosition() const
        {
            return *vertex;
        }
        plint i;                    // Global index of the vertex
        Array<T, 3> const *vertex;  // Pointer to vertex coordinates
    };

    typedef std::set<VertexSetNode, PositionLessThan3D<T, VertexSetNode> > VertexSet;
    typedef typename VertexSet::iterator VertexSetIterator;
    typedef typename VertexSet::const_iterator VertexSetConstIterator;

private:
    TriangleMesh<T> &triangleMesh;
    T eps;
    std::vector<Array<T, 3> > vertices;       // Positions of vertices. The vector
                                              // index is the global vertex index.
    std::vector<Array<plint, 3> > triangles;  // Global indices of the vertices that
                                              // constitute the triangle. The vector
                                              // index is the global triangle index.
    std::vector<plint> partTagging;
};

template <typename T>
void inflate(ConnectedTriangleMesh<T> &mesh, T distance);

template <typename T>
void transformVertexProperty(
    ConnectedTriangleMesh<T> &mesh, plint whichProperty, T scale, T offset = 0.0);

template <typename T>
void transformVertexProperty(
    ConnectedTriangleMesh<T> &mesh, std::string propertyName, T scale, T offset = 0.0);

}  // namespace plb

#endif  // CONNECTED_TRIANGLE_UTIL_H
