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

#ifndef OBSOLETE_FORMAT_WRAPPER_HH
#define OBSOLETE_FORMAT_WRAPPER_HH

#include "offLattice/obsoleteFormatWrapper.h"

namespace plb {

template <typename T>
RawTriangleMesh<T> triangleSetToRawTriangleMesh(TriangleSet<T> const &triangleSet, T eps)
{
    return RawTriangleMesh<T>(triangleSet.getTriangles(), eps);
}

template <typename T>
RawConnectedTriangleMesh<T> triangleSetToConnectedTriangleMesh(
    TriangleSet<T> const &triangleSet, T eps)
{
    RawTriangleMesh<T> rawMesh(triangleSet.getTriangles());
    return MeshConnector<T>(rawMesh, eps).generateConnectedMesh();
}

template <typename T>
TriangleSet<T> rawTriangleMeshToTriangleSet(
    RawTriangleMesh<T> const &triangleMesh, Precision precision)
{
    typedef typename RawTriangleMesh<T>::RawTriangle RawTriangle;
    std::vector<std::vector<RawTriangle> > const &triangles = triangleMesh.getTriangles();
    std::vector<RawTriangle> allTriangles;
    for (pluint i = 0; i < triangles.size(); ++i) {
        allTriangles.insert(allTriangles.end(), triangles[i].begin(), triangles[i].end());
    }
    return TriangleSet<T>(allTriangles, precision);
}

template <typename T>
TriangleSet<T> rawTriangleMeshToTriangleSet(RawTriangleMesh<T> const &triangleMesh, T eps)
{
    typedef typename RawTriangleMesh<T>::RawTriangle RawTriangle;
    std::vector<std::vector<RawTriangle> > const &triangles = triangleMesh.getTriangles();
    std::vector<RawTriangle> allTriangles;
    for (pluint i = 0; i < triangles.size(); ++i) {
        allTriangles.insert(allTriangles.end(), triangles[i].begin(), triangles[i].end());
    }
    return TriangleSet<T>(allTriangles, eps);
}

template <typename T>
RawConnectedTriangleMesh<T> def_to_ConnectedMesh(
    TriangleBoundary3D<T> &boundary, std::string partNaming)
{
    TriangularSurfaceMesh<T> &mesh = boundary.getMesh();
    plint numVertices = mesh.getNumVertices();
    plint numTriangles = mesh.getNumTriangles();

    std::vector<Array<T, 3> > vertices(numVertices);
    std::vector<Array<plint, 3> > triangles(numTriangles);
    std::vector<plint> partTagging(numTriangles);
    std::vector<std::string> nameOfParts;

    vertices.resize(numVertices);
    for (pluint i = 0; i < vertices.size(); ++i) {
        vertices[i] = mesh.getVertex(i);
    }

    plint numParts = 0;
    for (pluint i = 0; i < triangles.size(); ++i) {
        triangles[i][0] = mesh.getVertexId(i, 0);
        triangles[i][1] = mesh.getVertexId(i, 1);
        triangles[i][2] = mesh.getVertexId(i, 2);
        partTagging[i] = boundary.getTag(i);
        numParts = std::max(numParts, partTagging[i] + 1);
    }

    nameOfParts.push_back("Body");
    for (plint i = 1; i < numParts; ++i) {
        nameOfParts.push_back(partNaming + util::val2str(i - 1));
    }

    return RawConnectedTriangleMesh<T>(vertices, triangles, partTagging, nameOfParts);
}

}  // namespace plb

#endif  // OBSOLETE_FORMAT_WRAPPER_HH
