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

#ifndef CONNECTED_TRIANGLE_UTIL_HH
#define CONNECTED_TRIANGLE_UTIL_HH

#include "core/geometry3D.h"
#include "core/util.h"
#include "offLattice/connectedTriangleUtil.h"

namespace plb {

/* ******** Caching triangle and vertex properties.  ******************** */

template <typename T>
void createTriangleAreaCache(ConnectedTriangleMesh<T> &mesh, std::string triangleProperty)
{
    plint property = mesh.registerTriangleProperty(triangleProperty);
    RawConnectedTriangleMesh<double>::PTriangleIterator it = mesh.triangleIterator();
    while (!it->end()) {
        RawConnectedTriangleMesh<double>::PTriangle triangle = it->next();
        T area = triangle->area();
        triangle->setProperty(property, area);
    }
}

template <typename T>
void createVertexAreaCache(ConnectedTriangleMesh<T> &mesh, std::string vertexProperty)
{
    plint property = mesh.registerVertexProperty(vertexProperty);
    RawConnectedTriangleMesh<double>::PVertexIterator it = mesh.vertexIterator();
    while (!it->end()) {
        RawConnectedTriangleMesh<double>::PVertex vertex = it->next();
        T area = vertex->area();
        vertex->setProperty(property, area);
    }
}

/* ******** General Utilities ******************************************* */

template <typename T>
RawConnectedTriangleMesh<T> extractConnectedPart(ConnectedTriangleMesh<T> &mesh, plint whichPart)
{
    typedef typename ConnectedTriangleMesh<T>::PTriangleIterator PTriangleIterator;
    typedef typename ConnectedTriangleMesh<T>::PVertexIterator PVertexIterator;
    typedef typename ConnectedTriangleMesh<T>::PTriangle PTriangle;
    typedef typename ConnectedTriangleMesh<T>::PVertex PVertex;

    plint numVertexProperties = mesh.numVertexProperties();
    plint numTriangleProperties = mesh.numTriangleProperties();

    std::vector<std::vector<plint> > vertexTags, triangleTags;
    std::vector<std::vector<T> > vertexProperties(numVertexProperties),
        triangleProperties(numTriangleProperties);

    plint vertexIDtagging = mesh.getVertexTag("UniqueID");
    std::vector<Array<T, 3> > vertices;
    std::map<plint, plint> globalToLocalVertex;
    PTriangleIterator triangleIt = mesh.triangleIterator(whichPart);
    plint nextLocalVertexId = 0;
    while (!triangleIt->end()) {
        PTriangle triangle = triangleIt->next();
        // Store tags for each local triangle.
        for (pluint i = 0; i < triangleTags.size(); ++i) {
            // We happen to also store the tags for the unique-id and the part-tagging.
            // They are not needed and will be automatically discarded, but this way the code
            // is kept simpler.
            triangleTags[i].push_back(triangle->tag(i));
        }
        // Store properties for each local triangle.
        for (pluint i = 0; i < triangleProperties.size(); ++i) {
            triangleProperties[i].push_back(triangle->property(i));
        }
        for (plint iVert = 0; iVert < 3; ++iVert) {
            plint globalVertexId = triangle->vertex(iVert)->tag(vertexIDtagging);
            std::map<plint, plint>::const_iterator it = globalToLocalVertex.find(globalVertexId);
            if (it == globalToLocalVertex.end()) {  // A new local vertex has been spotted.
                // Store tags for each local vertex.
                for (pluint i = 0; i < vertexTags.size(); ++i) {
                    // We happen to also store the tags for the unique-id. They are not needed and
                    // will be automatically discarded, but this way the code is kept simpler.
                    vertexTags[i].push_back(triangle->vertex(iVert)->tag(i));
                }
                // Store properties for each local vertex.
                for (pluint i = 0; i < vertexProperties.size(); ++i) {
                    vertexProperties[i].push_back(triangle->vertex(iVert)->property(i));
                }
                // Store a mapping from global to local vertex.
                globalToLocalVertex.insert(std::make_pair(globalVertexId, nextLocalVertexId++));
                vertices.push_back(triangle->vertex(iVert)->get());
            }
        }
    }

    // Loop a second time to create the triangles based on local vertex coordinates.
    std::vector<Array<plint, 3> > triangles;
    triangleIt = mesh.triangleIterator(whichPart);
    while (!triangleIt->end()) {
        PTriangle triangle = triangleIt->next();
        Array<plint, 3> newTriangle;
        for (plint iVert = 0; iVert < 3; ++iVert) {
            plint globalVertexId = triangle->vertex(iVert)->tag(vertexIDtagging);
            newTriangle[iVert] = globalToLocalVertex[globalVertexId];
        }
        triangles.push_back(newTriangle);
    }

    RawConnectedTriangleMesh<T> partMesh(vertices, triangles);
    for (plint i = 0; i < (plint)triangleTags.size(); ++i) {
        partMesh.registerTriangleTag(mesh.getTriangleTagName(i));
    }
    for (plint i = 0; i < (plint)triangleProperties.size(); ++i) {
        partMesh.registerTriangleProperty(mesh.getTrianglePropertyName(i));
    }
    for (plint i = 0; i < (plint)vertexTags.size(); ++i) {
        partMesh.registerVertexTag(mesh.getVertexTagName(i));
    }
    for (plint i = 0; i < (plint)vertexProperties.size(); ++i) {
        partMesh.registerVertexProperty(mesh.getVertexPropertyName(i));
    }

    PTriangleIterator partTriangleIt = partMesh.triangleIterator();
    plint iTriangle = 0;
    while (!partTriangleIt->end()) {
        PTriangle triangle = partTriangleIt->next();
        for (plint iTag = 0; iTag < (plint)triangleTags.size(); ++iTag) {
            triangle->setTag(iTag, triangleTags[iTag][iTriangle]);
        }
        for (plint iProp = 0; iProp < (plint)triangleProperties.size(); ++iProp) {
            triangle->setProperty(iProp, triangleProperties[iProp][iTriangle]);
        }
        ++iTriangle;
    }

    PVertexIterator partVertexIt = partMesh.vertexIterator();
    plint iVertex = 0;
    while (!partVertexIt->end()) {
        PVertex vertex = partVertexIt->next();
        for (plint iTag = 0; iTag < (plint)vertexTags.size(); ++iTag) {
            vertex->setTag(iTag, vertexTags[iTag][iVertex]);
        }
        for (plint iProp = 0; iProp < (plint)vertexProperties.size(); ++iProp) {
            vertex->setProperty(iProp, vertexProperties[iProp][iVertex]);
        }
        ++iVertex;
    }

    return partMesh;
}

template <typename T>
std::vector<Array<T, 3> > getAllUniqueVertices(ConnectedTriangleMesh<T> &mesh)
{
    typedef typename ConnectedTriangleMesh<T>::PVertexIterator PVertexIterator;
    typedef typename ConnectedTriangleMesh<T>::PVertex PVertex;

    plint numVertices = mesh.getNumVertices();
    std::vector<Array<T, 3> > vertices(numVertices);

    if (numVertices == 0) {
        return vertices;
    }

    plint vertexIDtagging = mesh.getVertexTag("UniqueID");
    PVertexIterator vertexIterator(mesh.vertexIterator());
    while (!vertexIterator->end()) {
        PVertex vertex(vertexIterator->next());
        vertices[vertex->tag(vertexIDtagging)] = vertex->get();
    }

    return vertices;
}

/* ******** Class MeshConnector ******************************************* */

template <typename T>
MeshConnector<T>::MeshConnector(TriangleMesh<T> &triangleMesh_, T eps_) :
    triangleMesh(triangleMesh_), eps(eps_)
{
    uniqueVertices();
}

template <typename T>
MeshConnector<T>::MeshConnector(RawTriangleMesh<T> &rawTriangleMesh_) :
    triangleMesh(rawTriangleMesh_), eps(rawTriangleMesh_.getEps())
{
    uniqueVertices();
}

template <typename T>
void MeshConnector<T>::uniqueVertices()
{
    typedef typename TriangleMesh<T>::PTriangleIterator PTriangleIterator;
    typedef typename TriangleMesh<T>::PTriangle PTriangle;
    typedef typename TriangleMesh<T>::PVertex PVertex;

    PTriangleIterator triangleIt = triangleMesh.triangleIterator();

    PositionLessThan3D<T, VertexSetNode> lessThan(eps);
    VertexSet vertexSet(lessThan);

    plint numVertices = 0;
    plint globalVertex = 0;
    triangles.clear();
    partTagging.clear();

    plint triangleIDtagging = triangleMesh.getTriangleTag("Part");
    std::vector<Array<T, 3> *> tmp;
    while (!triangleIt->end()) {
        triangles.push_back(Array<plint, 3>(0, 0, 0));
        PTriangle triangle = triangleIt->next();
        partTagging.push_back(triangle->tag(triangleIDtagging));
        for (plint i = 0; i < 3; ++i) {
            PVertex vertex = triangle->vertex(i);
            VertexSetIterator vertexIt = vertexSet.find(VertexSetNode(-1, &(vertex->get())));
            if (vertexIt == vertexSet.end()) {  // Vertex does not exist yet.
                globalVertex = numVertices++;
                Array<T, 3> &data = vertex->get();
                tmp.push_back(&data);
                VertexSetNode newNode(globalVertex, &data);
                vertexSet.insert(newNode);
            } else {  // Vertex has already been registered.
                globalVertex = vertexIt->i;
            }
            triangles.back()[i] = globalVertex;
        }
    }

    // Order and copy vertex coordinates.
    vertices.resize(numVertices);
    VertexSetConstIterator it = vertexSet.begin();
    for (; it != vertexSet.end(); ++it) {
        vertices[it->i] = *(it->vertex);
    }
}

template <typename T>
RawConnectedTriangleMesh<T> MeshConnector<T>::generateConnectedMesh() const
{
    std::vector<std::string> nameOfParts;
    for (plint i = 0; i < triangleMesh.numParts(); ++i) {
        nameOfParts.push_back(triangleMesh.partName(i));
    }
    return RawConnectedTriangleMesh<T>(vertices, triangles, partTagging, nameOfParts);
}

template <typename T>
void inflate(ConnectedTriangleMesh<T> &mesh, T distance)
{
    typedef typename TriangleMesh<T>::PVertexIterator PVertexIterator;
    typedef typename TriangleMesh<T>::PVertex PVertex;
    PVertexIterator vertexIt = mesh.vertexIterator();
    while (!vertexIt->end()) {
        PVertex vertex = vertexIt->next();
        Array<T, 3> normal = vertex->normal();
        Array<T, 3> &coordinates = vertex->get();
        coordinates += distance * normal;
    }
}

template <typename T>
void transformVertexProperty(ConnectedTriangleMesh<T> &mesh, plint whichProperty, T scale, T offset)
{
    if (whichProperty < 0 || whichProperty >= mesh.numVertexProperties()) {
        return;
    }
    typedef typename ConnectedTriangleMesh<T>::PVertexIterator PVertexIterator;
    typedef typename ConnectedTriangleMesh<T>::PVertex PVertex;
    PVertexIterator vertexIt = mesh.vertexIterator();
    while (!vertexIt->end()) {
        PVertex vertex = vertexIt->next();
        vertex->setProperty(whichProperty, scale * vertex->property(whichProperty) + offset);
    }
}

template <typename T>
void transformVertexProperty(
    ConnectedTriangleMesh<T> &mesh, std::string propertyName, T scale, T offset)
{
    plint whichProperty = mesh.getVertexProperty(propertyName);
    if (whichProperty >= 0 && whichProperty < mesh.numVertexProperties()) {
        transformVertexProperty(mesh, whichProperty, scale, offset);
    }
}

}  // namespace plb

#endif  // CONNECTED_TRIANGLE_UTIL_HH
