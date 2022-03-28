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

#ifndef RAW_CONNECTED_TRIANGLE_MESH_HH
#define RAW_CONNECTED_TRIANGLE_MESH_HH

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <istream>
#include <limits>
#include <ostream>
#include <string>
#include <vector>

#include "core/array.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "io/plbFiles.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "offLattice/connectedTriangleMesh.h"
#include "offLattice/connectedTriangleUtil.h"
#include "offLattice/offFileIO.h"
#include "offLattice/rawTriangleMesh.h"
#include "offLattice/stlFileIO.h"

namespace plb {

/* *************** Free Functions ******************************************* */

template <typename T>
RawConnectedTriangleMesh<T> stlToConnectedTriangleMesh(STLreader<T> const &reader)
{
    RawTriangleMesh<T> rawMesh(stlToRawTriangleMesh<T>(reader));
    return generateConnectedTriangleMesh<T>(rawMesh, reader.getEps());
}

template <typename T>
RawConnectedTriangleMesh<T> stlToConnectedTriangleMesh(FileName stlFileName, T eps)
{
    RawTriangleMesh<T> rawMesh(stlToRawTriangleMesh<T>(stlFileName, eps));
    return generateConnectedTriangleMesh<T>(rawMesh, eps);
}

template <typename T>
RawConnectedTriangleMesh<T> offToConnectedTriangleMesh(OFFreader<T> const &reader)
{
    std::vector<Array<T, 3> > const &vertices = reader.getVertices();
    std::vector<std::vector<plint> > const &facets = reader.getFacets();
    std::vector<Array<plint, 3> > triangles(facets.size());

    for (pluint iTriangle = 0; iTriangle < triangles.size(); iTriangle++) {
        PLB_ASSERT(
            facets[iTriangle].size() == 3);  // The surface in the OFF file must be triangulated.
        Array<plint, 3> globalVertices;
        for (plint localVertex = 0; localVertex < 3; localVertex++) {
            globalVertices[localVertex] = facets[iTriangle][localVertex];
        }
        triangles[iTriangle] = globalVertices;
    }

    return RawConnectedTriangleMesh<T>(vertices, triangles);
}

template <typename T>
RawConnectedTriangleMesh<T> offToConnectedTriangleMesh(FileName offFileName)
{
    OFFreader<T> reader(offFileName.get());
    return offToConnectedTriangleMesh<T>(reader);
}

template <typename T>
RawConnectedTriangleMesh<T> generateConnectedTriangleMesh(RawTriangleMesh<T> &rawMesh, T eps)
{
    return MeshConnector<T>(rawMesh, eps).generateConnectedMesh();
}

/* *************** Class ManualConnectedTriangleMesh ************************ */

template <typename T>
ManualConnectedTriangleMesh<T>::ManualConnectedTriangleMesh(
    plint numTriangles_, plint numVertices_) :
    numTriangles(numTriangles_), numVertices(numVertices_)
{
    // Add dummy triangle tags for the unique ID and for the part tagging.
    // In this case, there's only one part, and the tagging is useless.
    // And, the unique ID does not need to be stored anyway.
    // But triangleTags is a vector in which the index
    // is the ID for the tagging, so we're forced to take the slot anyway.
    triangleTags.push_back(std::vector<plint>());
    triangleTags.push_back(std::vector<plint>());
    triangleTagNames.push_back("Part");
    triangleTagNames.push_back("UniqueID");
    nameOfParts.push_back("Body");
    // Both triangleTags[0] and trianglesPerPart  remain empty, as there's only a
    // single part.

    // Add dummy vertex tags for the default unique numbering of vertices.
    vertexTags.push_back(std::vector<plint>());
    vertexTagNames.push_back("UniqueID");
}

template <typename T>
ManualConnectedTriangleMesh<T>::ManualConnectedTriangleMesh(
    plint numTriangles_, plint numVertices_, std::vector<plint> partTagging,
    std::vector<std::string> nameOfParts_) :
    numTriangles(numTriangles_), numVertices(numVertices_), nameOfParts(nameOfParts_)
{
    triangleTagNames.push_back("Part");
    triangleTagNames.push_back("UniqueID");
    PLB_ASSERT((plint)partTagging.size() == numTriangles);
    if (nameOfParts.size() == 1) {
        // If there's only one part, triangleTags[0] and trianglesPerPart
        // can remain empty.
        triangleTags.push_back(std::vector<plint>());
    } else {
        triangleTags.push_back(partTagging);
        // From the constructor, we got a triangle-to-part mapping, which we now
        // convert to a part-to-triangle map, so the user can iterate over parts.
        trianglesPerPart.resize(nameOfParts.size());
        for (plint triangle = 0; triangle < (plint)partTagging.size(); ++triangle) {
            plint part = partTagging[triangle];
            PLB_ASSERT(part < (plint)trianglesPerPart.size());
            trianglesPerPart[part].push_back(triangle);
        }
    }
    // Add dummy triangle tags for the unique ID. The unique ID does not
    // actually need to be stored.
    triangleTags.push_back(std::vector<plint>());

    // Add dummy vertex tags for the default unique numbering of vertices.
    vertexTags.push_back(std::vector<plint>());
    vertexTagNames.push_back("UniqueID");
}

template <typename T>
plint ManualConnectedTriangleMesh<T>::numParts() const
{
    return (plint)nameOfParts.size();
}

template <typename T>
std::string ManualConnectedTriangleMesh<T>::partName(plint iPart) const
{
    PLB_ASSERT(iPart < (plint)nameOfParts.size());
    return nameOfParts[iPart];
}

template <typename T>
plint ManualConnectedTriangleMesh<T>::partId(std::string partName) const
{
    std::vector<std::string>::const_iterator it =
        find(nameOfParts.begin(), nameOfParts.end(), partName);
    if (it == nameOfParts.end()) {
        return -1;
    } else {
        return it - nameOfParts.begin();
    }
}

template <typename T>
plint ManualConnectedTriangleMesh<T>::registerTriangleTag(std::string tagName)
{
    // Make sure this tag has not yet been registered.
    if (getTriangleTag(tagName) != -1) {
        return getTriangleTag(tagName);
    }
    plint tagID = (plint)triangleTags.size();
    triangleTags.resize(tagID + 1);
    triangleTags[tagID].resize(numTriangles);
    std::fill(triangleTags[tagID].begin(), triangleTags[tagID].end(), 0);
    triangleTagNames.push_back(tagName);
    PLB_ASSERT((plint)triangleTagNames.size() == tagID + 1);
    return tagID;
}

template <typename T>
plint ManualConnectedTriangleMesh<T>::numTriangleTags() const
{
    return (plint)triangleTags.size();
}

template <typename T>
plint ManualConnectedTriangleMesh<T>::registerTriangleProperty(std::string propertyName)
{
    // Make sure this property has not yet been registered.
    if (getTriangleProperty(propertyName) != -1) {
        return getTriangleProperty(propertyName);
    }
    plint propertyID = triangleProperties.size();
    triangleProperties.resize(propertyID + 1);
    triangleProperties[propertyID].resize(numTriangles);
    std::fill(triangleProperties[propertyID].begin(), triangleProperties[propertyID].end(), 0);
    trianglePropertyNames.push_back(propertyName);
    PLB_ASSERT((plint)trianglePropertyNames.size() == propertyID + 1);
    return propertyID;
}

template <typename T>
plint ManualConnectedTriangleMesh<T>::numTriangleProperties() const
{
    return (plint)triangleProperties.size();
}

template <typename T>
plint ManualConnectedTriangleMesh<T>::registerVertexTag(std::string tagName)
{
    // Make sure this tag has not yet been registered.
    if (getVertexTag(tagName) != -1) {
        return getVertexTag(tagName);
    }
    plint tagID = vertexTags.size();
    vertexTags.resize(tagID + 1);
    vertexTags[tagID].resize(numVertices);
    vertexTagNames.push_back(tagName);
    PLB_ASSERT((plint)vertexTagNames.size() == tagID + 1);
    return tagID;
}

template <typename T>
plint ManualConnectedTriangleMesh<T>::numVertexTags() const
{
    return (plint)vertexTags.size();
}

template <typename T>
plint ManualConnectedTriangleMesh<T>::registerVertexProperty(std::string propertyName)
{
    // Make sure this property has not yet been registered.
    if (getVertexProperty(propertyName) != -1) {
        return getVertexProperty(propertyName);
    }
    plint propertyID = vertexProperties.size();
    vertexProperties.resize(propertyID + 1);
    vertexProperties[propertyID].resize(numVertices);
    vertexPropertyNames.push_back(propertyName);
    PLB_ASSERT((plint)vertexPropertyNames.size() == propertyID + 1);
    return propertyID;
}

template <typename T>
plint ManualConnectedTriangleMesh<T>::numVertexProperties() const
{
    return (plint)vertexProperties.size();
}

template <typename T>
plint ManualConnectedTriangleMesh<T>::getTriangleTag(std::string tagName) const
{
    for (plint i = 0; i < (plint)triangleTagNames.size(); ++i) {
        if (triangleTagNames[i] == tagName) {
            return i;
        }
    }
    return -1;
}

template <typename T>
std::string ManualConnectedTriangleMesh<T>::getTriangleTagName(plint tag) const
{
    PLB_ASSERT(tag < (plint)triangleTagNames.size());
    return triangleTagNames[tag];
}

template <typename T>
plint ManualConnectedTriangleMesh<T>::getTriangleProperty(std::string propertyName) const
{
    for (plint i = 0; i < (plint)trianglePropertyNames.size(); ++i) {
        if (trianglePropertyNames[i] == propertyName) {
            return i;
        }
    }
    return -1;
}

template <typename T>
std::string ManualConnectedTriangleMesh<T>::getTrianglePropertyName(plint property) const
{
    PLB_ASSERT(property < (plint)trianglePropertyNames.size());
    return trianglePropertyNames[property];
}

template <typename T>
plint ManualConnectedTriangleMesh<T>::getVertexTag(std::string tagName) const
{
    for (plint i = 0; i < (plint)vertexTagNames.size(); ++i) {
        if (vertexTagNames[i] == tagName) {
            return i;
        }
    }
    return -1;
}

template <typename T>
std::string ManualConnectedTriangleMesh<T>::getVertexTagName(plint tag) const
{
    PLB_ASSERT(tag < (plint)vertexTagNames.size());
    return vertexTagNames[tag];
}

template <typename T>
plint ManualConnectedTriangleMesh<T>::getVertexProperty(std::string propertyName) const
{
    for (plint i = 0; i < (plint)vertexPropertyNames.size(); ++i) {
        if (vertexPropertyNames[i] == propertyName) {
            return i;
        }
    }
    return -1;
}

template <typename T>
std::string ManualConnectedTriangleMesh<T>::getVertexPropertyName(plint property) const
{
    PLB_ASSERT(property < (plint)vertexPropertyNames.size());
    return vertexPropertyNames[property];
}

template <typename T>
void ManualConnectedTriangleMesh<T>::clearTriangleTags()
{
    PLB_ASSERT(!triangleTags.empty());
    std::vector<plint> backupPartTagging = triangleTags[0];
    triangleTags.clear();
    triangleTags.push_back(backupPartTagging);
}

template <typename T>
void ManualConnectedTriangleMesh<T>::clearTriangleProperties()
{
    triangleProperties.clear();
}

template <typename T>
void ManualConnectedTriangleMesh<T>::clearVertexTags()
{
    vertexTags.clear();
}

template <typename T>
void ManualConnectedTriangleMesh<T>::clearVertexProperties()
{
    vertexProperties.clear();
}

template <typename T>
plint ManualConnectedTriangleMesh<T>::getNumTriangles() const
{
    return numTriangles;
}

template <typename T>
plint ManualConnectedTriangleMesh<T>::getNumVertices() const
{
    return numVertices;
}

/* *************** Class RawConnectedTriangleMesh::Vertex **************************** */

template <typename T>
RawConnectedTriangleMesh<T>::Vertex::Vertex(RawConnectedTriangleMesh<T> *mesh_, plint vertex_) :
    mesh(mesh_), vertex(vertex_)
{ }

template <typename T>
typename RawConnectedTriangleMesh<T>::Vertex *RawConnectedTriangleMesh<T>::Vertex::clone() const
{
    return new typename RawConnectedTriangleMesh<T>::Vertex(*this);
}

template <typename T>
Array<T, 3> const &RawConnectedTriangleMesh<T>::Vertex::get() const
{
    PLB_ASSERT(vertex < (plint)mesh->vertices.size());
    return mesh->vertices[vertex];
}

template <typename T>
Array<T, 3> &RawConnectedTriangleMesh<T>::Vertex::get()
{
    PLB_ASSERT(vertex < (plint)mesh->vertices.size());
    return mesh->vertices[vertex];
}

template <typename T>
T const &RawConnectedTriangleMesh<T>::Vertex::operator[](plint i) const
{
    PLB_ASSERT(vertex < (plint)mesh->vertices.size());
    return mesh->vertices[vertex][i];
}

template <typename T>
T &RawConnectedTriangleMesh<T>::Vertex::operator[](plint i)
{
    PLB_ASSERT(vertex < (plint)mesh->vertices.size());
    return mesh->vertices[vertex][i];
}

template <typename T>
Array<T, 3> RawConnectedTriangleMesh<T>::Vertex::normal(bool areaWeighted) const
{
    PLB_ASSERT(vertex < (plint)mesh->trianglesOnVertex.size());
    std::vector<plint> const &nbTriangles = mesh->trianglesOnVertex[vertex];
    Array<T, 3> n((T)0, (T)0, (T)0);
    for (pluint i = 0; i < nbTriangles.size(); ++i) {
        n += areaWeighted ? Triangle(mesh, nbTriangles[i]).normalTimesArea()
                          : Triangle(mesh, nbTriangles[i]).normal();
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
T RawConnectedTriangleMesh<T>::Vertex::area() const
{
    PLB_ASSERT(vertex < (plint)mesh->trianglesOnVertex.size());
    std::vector<plint> const &nbTriangles = mesh->trianglesOnVertex[vertex];
    T area = 0.0;
    for (pluint i = 0; i < nbTriangles.size(); ++i) {
        area += Triangle(mesh, nbTriangles[i]).area();
    }
    return area / (T)3;
}

template <typename T>
plint RawConnectedTriangleMesh<T>::Vertex::numAdjacentTriangles() const
{
    PLB_ASSERT(vertex < (plint)mesh->trianglesOnVertex.size());
    return mesh->trianglesOnVertex[vertex].size();
}

template <typename T>
typename RawConnectedTriangleMesh<T>::CPTriangle
    RawConnectedTriangleMesh<T>::Vertex::adjacentTriangle(plint iTriangle) const
{
    PLB_ASSERT(vertex < (plint)mesh->trianglesOnVertex.size());
    PLB_ASSERT(iTriangle < (plint)mesh->trianglesOnVertex[vertex].size());
    return CPTriangle(new Triangle(mesh, mesh->trianglesOnVertex[vertex][iTriangle]));
}

template <typename T>
typename RawConnectedTriangleMesh<T>::PTriangle
    RawConnectedTriangleMesh<T>::Vertex::adjacentTriangle(plint iTriangle)
{
    PLB_ASSERT(vertex < (plint)mesh->trianglesOnVertex.size());
    PLB_ASSERT(iTriangle < (plint)mesh->trianglesOnVertex[vertex].size());
    return PTriangle(new Triangle(mesh, mesh->trianglesOnVertex[vertex][iTriangle]));
}

template <typename T>
T RawConnectedTriangleMesh<T>::Vertex::property(plint whichProperty) const
{
    PLB_ASSERT(whichProperty >= 0);
    PLB_ASSERT(whichProperty < (plint)mesh->vertexProperties.size());
    PLB_ASSERT(vertex < (plint)mesh->vertexProperties[whichProperty].size());
    return mesh->vertexProperties[whichProperty][vertex];
}

template <typename T>
void RawConnectedTriangleMesh<T>::Vertex::setProperty(plint whichProperty, T value)
{
    PLB_ASSERT(whichProperty >= 0);
    PLB_ASSERT(whichProperty < (plint)mesh->vertexProperties.size());
    PLB_ASSERT(vertex < (plint)mesh->vertexProperties[whichProperty].size());
    mesh->vertexProperties[whichProperty][vertex] = value;
}

template <typename T>
plint RawConnectedTriangleMesh<T>::Vertex::tag(plint whichTag) const
{
    PLB_ASSERT(whichTag < (plint)mesh->vertexTags.size());
    // The 0-tag is reserved to return a unique number for each vertex.
    if (whichTag == 0) {
        return vertex;
    } else {
        PLB_ASSERT(vertex < (plint)mesh->vertexTags[whichTag].size());
        return mesh->vertexTags[whichTag][vertex];
    }
}

template <typename T>
void RawConnectedTriangleMesh<T>::Vertex::setTag(plint whichTag, plint value)
{
    PLB_ASSERT(whichTag >= 0);
    PLB_ASSERT(whichTag < (plint)mesh->vertexTags.size());
    PLB_ASSERT(vertex < (plint)mesh->vertexTags[whichTag].size());
    // You can't assign to the "tag-0" tag, which is reserved for the unique vertex numbering.
    if (whichTag > 0) {
        mesh->vertexTags[whichTag][vertex] = value;
    }
}

template <typename T>
bool RawConnectedTriangleMesh<T>::Vertex::isInterior() const
{
    std::vector<plint> neighbors = adjacentVertices();
    std::set<T> neighborSet(neighbors.begin(), neighbors.end());
    for (pluint i = 0; i < neighbors.size(); ++i) {
        std::vector<plint> neighborNeighbors = mesh->vertex(neighbors[i])->adjacentVertices();
        pluint nearbyNeighbors = 0;
        for (pluint j = 0; j < neighborNeighbors.size(); ++j) {
            typename std::set<T>::const_iterator it = neighborSet.find(neighborNeighbors[j]);
            if (it != neighborSet.end()) {
                ++nearbyNeighbors;
            }
        }
        if (nearbyNeighbors <= 1) {
            return false;
        }
    }
    return true;
}

template <typename T>
std::vector<plint> RawConnectedTriangleMesh<T>::Vertex::adjacentVertices() const
{
    plint uniqueIDTag = mesh->getVertexTag("UniqueID");
    plint numTriangles = numAdjacentTriangles();
    std::set<plint> neighborSet;
    for (plint iTriangle = 0; iTriangle < numTriangles; ++iTriangle) {
        CPTriangle triangle = adjacentTriangle(iTriangle);
        neighborSet.insert(triangle->vertex(0)->tag(uniqueIDTag));
        neighborSet.insert(triangle->vertex(1)->tag(uniqueIDTag));
        neighborSet.insert(triangle->vertex(2)->tag(uniqueIDTag));
    }
    std::vector<plint> vertices;
    plint myTag = this->tag(uniqueIDTag);
    for (std::set<plint>::iterator it = neighborSet.begin(); it != neighborSet.end(); ++it) {
        if (*it != myTag) {
            vertices.push_back(*it);
        }
    }
    return vertices;
    // return std::vector<plint>(neighborSet.begin(), neighborSet.end());
}

/* *************** Class RawConnectedTriangleMesh::Triangle **************************** */

template <typename T>
RawConnectedTriangleMesh<T>::Triangle::Triangle(
    RawConnectedTriangleMesh<T> *mesh_, plint triangle_) :
    mesh(mesh_), triangle(triangle_)
{ }

template <typename T>
typename RawConnectedTriangleMesh<T>::Triangle *RawConnectedTriangleMesh<T>::Triangle::clone() const
{
    return new typename RawConnectedTriangleMesh<T>::Triangle(*this);
}

template <typename T>
typename RawConnectedTriangleMesh<T>::CPVertex RawConnectedTriangleMesh<T>::Triangle::vertex(
    plint iVertex) const
{
    PLB_ASSERT(iVertex < 3);
    PLB_ASSERT(triangle < (plint)mesh->triangles.size());
    return CPVertex(new Vertex(mesh, mesh->triangles[triangle][iVertex]));
}

template <typename T>
typename RawConnectedTriangleMesh<T>::PVertex RawConnectedTriangleMesh<T>::Triangle::vertex(
    plint iVertex)
{
    PLB_ASSERT(iVertex < 3);
    PLB_ASSERT(triangle < (plint)mesh->triangles.size());
    return PVertex(new Vertex(mesh, mesh->triangles[triangle][iVertex]));
}

template <typename T>
Array<T, 3> const &RawConnectedTriangleMesh<T>::Triangle::operator[](plint iVertex) const
{
    PLB_ASSERT(iVertex < 3);
    PLB_ASSERT(triangle < (plint)mesh->triangles.size());
    return mesh->vertices[mesh->triangles[triangle][iVertex]];
}

template <typename T>
Array<T, 3> &RawConnectedTriangleMesh<T>::Triangle::operator[](plint iVertex)
{
    PLB_ASSERT(iVertex < 3);
    PLB_ASSERT(triangle < (plint)mesh->triangles.size());
    return mesh->vertices[mesh->triangles[triangle][iVertex]];
}

template <typename T>
T RawConnectedTriangleMesh<T>::Triangle::area() const
{
    return computeTriangleArea(
        mesh->vertices[mesh->triangles[triangle][0]], mesh->vertices[mesh->triangles[triangle][1]],
        mesh->vertices[mesh->triangles[triangle][2]]);
}

template <typename T>
Array<T, 3> RawConnectedTriangleMesh<T>::Triangle::normal() const
{
    Array<T, 3> n(normalTimesArea());
    T normN = norm(n);
    if (!util::isZero(normN)) {
        n /= normN;
    } else {
        n.resetToZero();
    }
    return n;
}

template <typename T>
Array<T, 3> RawConnectedTriangleMesh<T>::Triangle::normalTimesArea() const
{
    Array<T, 3> const &v0 = mesh->vertices[mesh->triangles[triangle][0]];
    Array<T, 3> const &v1 = mesh->vertices[mesh->triangles[triangle][1]];
    Array<T, 3> const &v2 = mesh->vertices[mesh->triangles[triangle][2]];
    bool isAreaWeighted = true;
    return plb::computeTriangleNormal(v0, v1, v2, isAreaWeighted);
}

template <typename T>
Array<T, 3> RawConnectedTriangleMesh<T>::Triangle::edgeNormal(plint iEdge, bool areaWeighted) const
{
    PLB_ASSERT(iEdge >= 0 && iEdge <= 2);
    Array<T, 3> result(areaWeighted ? normalTimesArea() : normal());
    CPTriangle triangle2(this->edgeNeighbor(iEdge));
    if (triangle2.get()) {
        result += areaWeighted ? triangle2->normalTimesArea() : triangle2->normal();
    }

    T normN = norm(result);
    if (!util::isZero(normN)) {
        result /= normN;
    } else {
        result.resetToZero();
    }

    return result;
}
template <typename T>
Array<T, 3> RawConnectedTriangleMesh<T>::Triangle::continuousNormal(
    Array<T, 3> const &p, bool areaWeighted) const
{
    T area = this->area();
    if (util::isZero(area)) {
        return Array<T, 3>((T)0, (T)0, (T)0);
    }

    Array<T, 3> v0((*this)[0]);
    Array<T, 3> v1((*this)[1]);
    Array<T, 3> v2((*this)[2]);

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

    Array<T, 3> n0 = vertex(0)->normal(areaWeighted);
    Array<T, 3> n1 = vertex(1)->normal(areaWeighted);
    Array<T, 3> n2 = vertex(2)->normal(areaWeighted);

    n = u * n0 + v * n1 + ((T)1. - u - v) * n2;
    T normN = norm(n);
    PLB_ASSERT(!util::isZero(normN));
    n /= normN;
    return n;
}

template <typename T>
typename RawConnectedTriangleMesh<T>::CPTriangle
    RawConnectedTriangleMesh<T>::Triangle::edgeNeighbor(plint iEdge) const
{
    PLB_ASSERT(iEdge >= 0 && iEdge <= 2);
    plint iVertex = iEdge;
    plint jVertex = (iEdge + 1) % 3;
    std::vector<plint> const &iVertexNeighbors =
        mesh->trianglesOnVertex[mesh->triangles[triangle][iVertex]];
    std::vector<plint> const &jVertexNeighbors =
        mesh->trianglesOnVertex[mesh->triangles[triangle][jVertex]];
    CPTriangle neighborTriangle(nullptr);
    plint count = 0;
    for (pluint iN = 0; iN < iVertexNeighbors.size(); iN++) {
        plint neighbor = iVertexNeighbors[iN];
        if (neighbor != triangle) {
            for (pluint jN = 0; jN < jVertexNeighbors.size(); jN++) {
                if (neighbor == jVertexNeighbors[jN]) {
                    ++count;
                    neighborTriangle = CPTriangle(new Triangle(mesh, neighbor));
                }
            }
        }
    }
    // If there is more than 1 neighbor, the edge is non-manifold. In this case
    // we return nothing, as the user should use the method "edgeNeighbors" instead.
    if (count > 1) {
        neighborTriangle = CPTriangle(nullptr);
    }
    return neighborTriangle;
}

template <typename T>
std::vector<plint> RawConnectedTriangleMesh<T>::Triangle::edgeNeighbors(plint iEdge) const
{
    PLB_ASSERT(iEdge >= 0 && iEdge <= 2);
    plint iVertex = iEdge;
    plint jVertex = (iEdge + 1) % 3;
    std::vector<plint> const &iVertexNeighbors =
        mesh->trianglesOnVertex[mesh->triangles[triangle][iVertex]];
    std::vector<plint> const &jVertexNeighbors =
        mesh->trianglesOnVertex[mesh->triangles[triangle][jVertex]];
    std::vector<plint> neighborTriangles;
    for (pluint iN = 0; iN < iVertexNeighbors.size(); iN++) {
        plint neighbor = iVertexNeighbors[iN];
        if (neighbor != triangle) {
            for (pluint jN = 0; jN < jVertexNeighbors.size(); jN++) {
                if (neighbor == jVertexNeighbors[jN]) {
                    neighborTriangles.push_back(neighbor);
                }
            }
        }
    }
    return neighborTriangles;
}

template <typename T>
plint RawConnectedTriangleMesh<T>::Triangle::numVertexNeighbors(plint iVertex) const
{
    PLB_ASSERT(iVertex >= 0 && iVertex <= 2);
    plint vertex = mesh->triangles[triangle][iVertex];
    return (plint)mesh->trianglesOnVertex[vertex].size();
}

template <typename T>
typename RawConnectedTriangleMesh<T>::CPTriangle
    RawConnectedTriangleMesh<T>::Triangle::vertexNeighbor(plint iVertex, plint iNeighbor) const
{
    plint numVertexNeighbors = this->numVertexNeighbors(iVertex);
    if (iNeighbor < 0 || iNeighbor >= numVertexNeighbors) {
        return CPTriangle(nullptr);
    }
    plint vertex = mesh->triangles[triangle][iVertex];
    return CPTriangle(new Triangle(mesh, mesh->trianglesOnVertex[vertex][iNeighbor]));
}

template <typename T>
T RawConnectedTriangleMesh<T>::Triangle::property(plint whichProperty) const
{
    PLB_ASSERT(whichProperty >= 0);
    PLB_ASSERT(whichProperty < (plint)mesh->triangleProperties.size());
    PLB_ASSERT(triangle < (plint)mesh->triangleProperties[whichProperty].size());
    return mesh->triangleProperties[whichProperty][triangle];
}

template <typename T>
void RawConnectedTriangleMesh<T>::Triangle::setProperty(plint whichProperty, T value)
{
    PLB_ASSERT(whichProperty >= 0);
    PLB_ASSERT(whichProperty < (plint)mesh->triangleProperties.size());
    PLB_ASSERT(triangle < (plint)mesh->triangleProperties[whichProperty].size());
    mesh->triangleProperties[whichProperty][triangle] = value;
}

template <typename T>
plint RawConnectedTriangleMesh<T>::Triangle::tag(plint whichTag) const
{
    PLB_ASSERT(whichTag >= 0);
    PLB_ASSERT(whichTag < (plint)mesh->triangleTags.size());
    if (whichTag == 0) {                             // Tag 0 is for geometry parts.
        if (mesh->triangleTags[whichTag].empty()) {  // There's only one part.
            return 0;
        } else {
            PLB_ASSERT(triangle < (plint)mesh->triangleTags[whichTag].size());
            return mesh->triangleTags[whichTag][triangle];
        }
    } else if (whichTag == 1) {  // Tag 1 is for unique triangle ID.
        return triangle;
    } else {
        PLB_ASSERT(triangle < (plint)mesh->triangleTags[whichTag].size());
        return mesh->triangleTags[whichTag][triangle];
    }
}

template <typename T>
void RawConnectedTriangleMesh<T>::Triangle::setTag(plint whichTag, plint value)
{
    PLB_ASSERT(whichTag < (plint)mesh->triangleTags.size());
    // Tag 0 is for the unique ID, and tag 1 is for geometry parts.
    if (whichTag >= 2) {
        PLB_ASSERT(triangle < (plint)mesh->triangleTags[whichTag].size());
        mesh->triangleTags[whichTag][triangle] = value;
    }
}

/* *************** Class RawConnectedTriangleMesh::AllTriangleIterator ****************************
 */

template <typename T>
RawConnectedTriangleMesh<T>::AllTriangleIterator::AllTriangleIterator(
    RawConnectedTriangleMesh<T> *mesh_) :
    mesh(mesh_), currentTriangle(0)
{ }

template <typename T>
typename RawConnectedTriangleMesh<T>::PTriangle
    RawConnectedTriangleMesh<T>::AllTriangleIterator::next()
{
    PLB_ASSERT(currentTriangle < (plint)mesh->triangles.size());

    Triangle *triangle = new Triangle(mesh, currentTriangle);
    ++currentTriangle;
    return PTriangle(triangle);
}

template <typename T>
bool RawConnectedTriangleMesh<T>::AllTriangleIterator::end() const
{
    return currentTriangle >= (plint)mesh->triangles.size();
}

template <typename T>
typename RawConnectedTriangleMesh<T>::PTriangleIterator
    RawConnectedTriangleMesh<T>::AllTriangleIterator::clone() const
{
    return PTriangleIterator(new AllTriangleIterator(*this));
}

/* *************** Class RawConnectedTriangleMesh::PartTriangleIterator ****************************
 */

template <typename T>
RawConnectedTriangleMesh<T>::PartTriangleIterator::PartTriangleIterator(
    RawConnectedTriangleMesh<T> *mesh_, plint partId_) :
    mesh(mesh_), partId(partId_), currentTriangle(0)
{
    PLB_ASSERT(partId < (plint)mesh->trianglesPerPart.size());
}

template <typename T>
typename RawConnectedTriangleMesh<T>::PTriangle
    RawConnectedTriangleMesh<T>::PartTriangleIterator::next()
{
    PLB_ASSERT(partId < (plint)mesh->trianglesPerPart.size());
    PLB_ASSERT(currentTriangle < (plint)mesh->trianglesPerPart[partId].size());

    Triangle *triangle = new Triangle(mesh, mesh->trianglesPerPart[partId][currentTriangle]);
    ++currentTriangle;
    return PTriangle(triangle);
}

template <typename T>
bool RawConnectedTriangleMesh<T>::PartTriangleIterator::end() const
{
    PLB_ASSERT(partId < (plint)mesh->trianglesPerPart.size());
    return currentTriangle >= (plint)mesh->trianglesPerPart[partId].size();
}

template <typename T>
typename RawConnectedTriangleMesh<T>::PTriangleIterator
    RawConnectedTriangleMesh<T>::PartTriangleIterator::clone() const
{
    return PTriangleIterator(new PartTriangleIterator(*this));
}

/* *************** Class RawConnectedTriangleMesh::VertexIterator **************************** */

template <typename T>
RawConnectedTriangleMesh<T>::VertexIterator::VertexIterator(RawConnectedTriangleMesh<T> *mesh_) :
    mesh(mesh_), currentVertex(0)
{ }

template <typename T>
typename RawConnectedTriangleMesh<T>::PVertex RawConnectedTriangleMesh<T>::VertexIterator::next()
{
    PLB_ASSERT(currentVertex < (plint)mesh->vertices.size());

    Vertex *vertex = new Vertex(mesh, currentVertex);
    ++currentVertex;
    return PVertex(vertex);
}

template <typename T>
bool RawConnectedTriangleMesh<T>::VertexIterator::end() const
{
    return currentVertex >= (plint)mesh->vertices.size();
}

template <typename T>
typename RawConnectedTriangleMesh<T>::PVertexIterator
    RawConnectedTriangleMesh<T>::VertexIterator::clone() const
{
    return PVertexIterator(new VertexIterator(*this));
}

/* *************** Class RawConnectedTriangleMesh  *************************************************
 */

template <typename T>
RawConnectedTriangleMesh<T>::RawConnectedTriangleMesh(
    std::vector<Array<T, 3> > const &vertices_, std::vector<Array<plint, 3> > const &triangles_) :
    ManualConnectedTriangleMesh<T>((plint)triangles_.size(), (plint)vertices_.size()),
    vertices(vertices_),
    triangles(triangles_)
{
    computeTrianglesOnVertex();
    signalVertexUpdate();
}

template <typename T>
RawConnectedTriangleMesh<T>::RawConnectedTriangleMesh(
    std::vector<Array<T, 3> > const &vertices_, std::vector<Array<plint, 3> > const &triangles_,
    std::vector<plint> const &partTagging, std::vector<std::string> const &nameOfParts) :
    ManualConnectedTriangleMesh<T>(
        (plint)triangles_.size(), (plint)vertices_.size(), partTagging, nameOfParts),
    vertices(vertices_),
    triangles(triangles_)
{
    computeTrianglesOnVertex();
    signalVertexUpdate();
}

template <typename T>
RawConnectedTriangleMesh<T> RawConnectedTriangleMesh<T>::merge(
    RawConnectedTriangleMesh<T> &rhs, T eps)
{
    RawTriangleMesh<T> thisRawTriangleMesh((TriangleMesh<T> &)*this, eps);
    RawTriangleMesh<T> rhsRawTriangleMesh(rhs, eps);
    RawTriangleMesh<T> mergedRawTriangleMesh = thisRawTriangleMesh.merge(rhsRawTriangleMesh);
    RawConnectedTriangleMesh<T> newMesh =
        MeshConnector<T>(mergedRawTriangleMesh).generateConnectedMesh();
    return newMesh;
}

template <typename T>
RawConnectedTriangleMesh<T> RawConnectedTriangleMesh<T>::refine(T eps) const
{
    RawTriangleMesh<T> rawTriangleMesh((TriangleMesh<T> &)*this, eps);
    RawTriangleMesh<T> refinedRawTriangleMesh = rawTriangleMesh.refine();
    RawConnectedTriangleMesh<T> newMesh =
        MeshConnector<T>(refinedRawTriangleMesh).generateConnectedMesh();
    return newMesh;
}

template <typename T>
RawConnectedTriangleMesh<T> RawConnectedTriangleMesh<T>::refineRecursively(
    T targetMaxEdgeLength, plint maxNumIterations, bool &success, T eps) const
{
    RawTriangleMesh<T> rawTriangleMesh((TriangleMesh<T> &)*this, eps);
    RawTriangleMesh<T> refinedRawTriangleMesh =
        rawTriangleMesh.refineRecursively(targetMaxEdgeLength, maxNumIterations, success);
    RawConnectedTriangleMesh<T> newMesh =
        MeshConnector<T>(refinedRawTriangleMesh).generateConnectedMesh();
    return newMesh;
}

template <typename T>
RawConnectedTriangleMesh<T> RawConnectedTriangleMesh<T>::select(
    TriangleSelector<T> const &selector, T eps) const
{
    RawTriangleMesh<T> rawTriangleMesh((TriangleMesh<T> &)*this, eps);
    RawTriangleMesh<T> subRawTriangleMesh = rawTriangleMesh.select(selector);
    RawConnectedTriangleMesh<T> newMesh =
        MeshConnector<T>(subRawTriangleMesh).generateConnectedMesh();
    return newMesh;
}

template <typename T>
RawConnectedTriangleMesh<T> RawConnectedTriangleMesh<T>::cutWithPlane(
    Array<T, 3> planePos, Array<T, 3> normal, T eps) const
{
    RawTriangleMesh<T> rawTriangleMesh((TriangleMesh<T> &)*this, eps);
    RawTriangleMesh<T> cutTriangleMesh = rawTriangleMesh.cutWithPlane(planePos, normal);
    RawConnectedTriangleMesh<T> newMesh = MeshConnector<T>(cutTriangleMesh).generateConnectedMesh();
    return newMesh;
}

template <typename T>
typename RawConnectedTriangleMesh<T>::PTriangleIterator
    RawConnectedTriangleMesh<T>::triangleIterator(plint partId)
{
    PLB_ASSERT(partId < this->numParts());
    if (partId == -1 || this->numParts() == 1) {
        return PTriangleIterator(new AllTriangleIterator(this));
    } else {
        return PTriangleIterator(new PartTriangleIterator(this, partId));
    }
}

template <typename T>
typename RawConnectedTriangleMesh<T>::PTriangle RawConnectedTriangleMesh<T>::triangle(
    plint uniqueID)
{
    if (uniqueID < (plint)this->triangles.size()) {
        return PTriangle(new Triangle(this, uniqueID));
    } else {
        return PTriangle(nullptr);
    }
}

template <typename T>
typename RawConnectedTriangleMesh<T>::PVertexIterator RawConnectedTriangleMesh<T>::vertexIterator()
{
    return PVertexIterator(new VertexIterator(this));
}

template <typename T>
typename RawConnectedTriangleMesh<T>::PVertex RawConnectedTriangleMesh<T>::vertex(plint uniqueID)
{
    if (uniqueID < (plint)this->vertices.size()) {
        return PVertex(new Vertex(this, uniqueID));
    } else {
        return PVertex(nullptr);
    }
}

template <typename T>
void RawConnectedTriangleMesh<T>::signalVertexUpdate()
{
    minMaxEdgeLengthValid = false;
    minMaxTriangleAreaValid = false;
    boundingCuboidValid = false;
}

template <typename T>
void RawConnectedTriangleMesh<T>::signalIsometricVertexUpdate()
{
    boundingCuboidValid = false;
}

template <typename T>
void RawConnectedTriangleMesh<T>::translate(Array<T, 3> const &vector)
{
    for (pluint i = 0; i < vertices.size(); ++i) {
        vertices[i] += vector;
    }
    boundingCuboid.lowerLeftCorner += vector;
    boundingCuboid.upperRightCorner += vector;
}

template <typename T>
void RawConnectedTriangleMesh<T>::scale(T alpha)
{
    for (pluint i = 0; i < vertices.size(); ++i) {
        vertices[i] *= alpha;
    }
    minMaxEdgeLengthValid = false;
    minMaxTriangleAreaValid = false;
    boundingCuboidValid = false;
}

template <typename T>
void RawConnectedTriangleMesh<T>::rotateAtOrigin(Array<T, 3> const &normedAxis, T theta)
{
    for (pluint i = 0; i < vertices.size(); ++i) {
        vertices[i] = plb::rotateAtOrigin(vertices[i], normedAxis, theta);
    }
    boundingCuboidValid = false;
}

template <typename T>
void RawConnectedTriangleMesh<T>::rotate(T phi, T theta, T psi)
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

    for (pluint iVertex = 0; iVertex < vertices.size(); ++iVertex) {
        Array<T, 3> x = vertices[iVertex];
        for (int i = 0; i < 3; i++) {
            vertices[iVertex][i] = (T)0.0;
            for (int j = 0; j < 3; j++) {
                vertices[iVertex][i] += a[i][j] * x[j];
            }
        }
    }

    boundingCuboidValid = false;
}

template <typename T>
void RawConnectedTriangleMesh<T>::reverseOrientation()
{
    for (pluint i = 0; i < triangles.size(); ++i) {
        std::swap(triangles[i][1], triangles[i][2]);
    }
}

template <typename T>
T RawConnectedTriangleMesh<T>::getMinEdgeLength() const
{
    computeMinMaxEdges();  // Uses lazy evaluation.
    return minEdgeLength;
}

template <typename T>
T RawConnectedTriangleMesh<T>::getMaxEdgeLength() const
{
    computeMinMaxEdges();  // Uses lazy evaluation.
    return maxEdgeLength;
}

template <typename T>
T RawConnectedTriangleMesh<T>::getMinTriangleArea() const
{
    computeMinMaxAreas();  // Uses lazy evaluation.
    return minTriangleArea;
}

template <typename T>
T RawConnectedTriangleMesh<T>::getMaxTriangleArea() const
{
    computeMinMaxAreas();  // Uses lazy evaluation.
    return maxTriangleArea;
}

template <typename T>
Cuboid<T> RawConnectedTriangleMesh<T>::getBoundingCuboid() const
{
    computeBoundingCuboid();  // Uses lazy evaluation.
    return boundingCuboid;
}

template <typename T>
void RawConnectedTriangleMesh<T>::verticesToFile(FileName fileName) const
{
    if (global::mpi().isMainProcessor()) {
        std::ofstream ofile(fileName.get().c_str(), std::ios_base::out | std::ios_base::binary);
        ofile.write((const char *)(&vertices[0][0]), vertices.size() * 3 * sizeof(T));
    }
}

template <typename T>
void RawConnectedTriangleMesh<T>::verticesFromFile(FileName fileName)
{
    std::ifstream ifile(fileName.get().c_str(), std::ios_base::in | std::ios_base::binary);
    int numchar = vertices.size() * 3 * sizeof(T);
    char *data = new char[numchar];
    ifile.read(data, numchar);
    memcpy(&vertices[0][0], data, numchar);
    signalVertexUpdate();
}

template <typename T>
void RawConnectedTriangleMesh<T>::computeTrianglesOnVertex()
{
    trianglesOnVertex.resize(vertices.size());
    for (pluint iTriangle = 0; iTriangle < triangles.size(); ++iTriangle) {
        for (plint i = 0; i < 3; ++i) {
            plint iVertex = triangles[iTriangle][i];
            trianglesOnVertex[iVertex].push_back(iTriangle);
        }
    }
}

template <typename T>
void RawConnectedTriangleMesh<T>::computeMinMaxEdges() const
{
    if (minMaxEdgeLengthValid)
        return;  // Lazy evaluation.
    minEdgeLength = std::numeric_limits<T>::max();
    maxEdgeLength = std::numeric_limits<T>::min();
    for (pluint i = 0; i < triangles.size(); ++i) {
        Array<T, 3> const &v0 = vertices[triangles[i][0]];
        Array<T, 3> const &v1 = vertices[triangles[i][1]];
        Array<T, 3> const &v2 = vertices[triangles[i][2]];
        T edge1 = norm(v1 - v0);
        T edge2 = norm(v2 - v1);
        T edge3 = norm(v0 - v2);
        T minEdge = std::min(edge1, std::min(edge2, edge3));
        T maxEdge = std::max(edge1, std::max(edge2, edge3));
        minEdgeLength = std::min(minEdgeLength, minEdge);
        maxEdgeLength = std::max(maxEdgeLength, maxEdge);
    }
    minMaxEdgeLengthValid = true;
}

template <typename T>
void RawConnectedTriangleMesh<T>::computeMinMaxAreas() const
{
    if (minMaxTriangleAreaValid)
        return;  // Lazy evaluation.
    minTriangleArea = std::numeric_limits<T>::max();
    maxTriangleArea = std::numeric_limits<T>::min();
    for (pluint i = 0; i < triangles.size(); ++i) {
        Array<T, 3> const &v0 = vertices[triangles[i][0]];
        Array<T, 3> const &v1 = vertices[triangles[i][1]];
        Array<T, 3> const &v2 = vertices[triangles[i][2]];
        T area = computeTriangleArea(v0, v1, v2);
        minTriangleArea = std::min(minTriangleArea, area);
        maxTriangleArea = std::max(maxTriangleArea, area);
    }
    minMaxTriangleAreaValid = true;
}

template <typename T>
void RawConnectedTriangleMesh<T>::computeBoundingCuboid() const
{
    if (boundingCuboidValid)
        return;  // Lazy evaluation.
    T xMin, yMin, zMin;
    T xMax, yMax, zMax;

    xMin = yMin = zMin = std::numeric_limits<T>::max();
    xMax = yMax = zMax = -std::numeric_limits<T>::max();
    for (pluint i = 0; i < triangles.size(); ++i) {
        Array<T, 3> const &v0 = vertices[triangles[i][0]];
        Array<T, 3> const &v1 = vertices[triangles[i][1]];
        Array<T, 3> const &v2 = vertices[triangles[i][2]];

        xMin = std::min(xMin, std::min(v0[0], std::min(v1[0], v2[0])));
        yMin = std::min(yMin, std::min(v0[1], std::min(v1[1], v2[1])));
        zMin = std::min(zMin, std::min(v0[2], std::min(v1[2], v2[2])));

        xMax = std::max(xMax, std::max(v0[0], std::max(v1[0], v2[0])));
        yMax = std::max(yMax, std::max(v0[1], std::max(v1[1], v2[1])));
        zMax = std::max(zMax, std::max(v0[2], std::max(v1[2], v2[2])));
    }
    boundingCuboid.lowerLeftCorner = Array<T, 3>(xMin, yMin, zMin);
    boundingCuboid.upperRightCorner = Array<T, 3>(xMax, yMax, zMax);
    boundingCuboidValid = true;
}

template <typename T>
RawConnectedTriangleMesh<T> disconnectedMerge(
    std::vector<RawConnectedTriangleMesh<T> const *> parts, std::vector<std::string> names,
    std::string vertexTagForParts)
{
    PLB_ASSERT(parts.size() == names.size());
    // First, join all vertices and triangles of the provided parts, and create a
    // RawConnectedTriangleMesh.
    std::vector<Array<T, 3> > allVertices;
    std::vector<Array<plint, 3> > allTriangles;
    std::vector<plint> triangleTags;
    std::vector<plint> vertexTags;
    plint vertexOffset = 0;
    for (plint i = 0; i < (plint)parts.size(); ++i) {
        std::vector<Array<T, 3> > const &vertices = parts[i]->getVertices();
        std::vector<Array<plint, 3> > triangles = parts[i]->getTriangles();
        for (pluint j = 0; j < triangles.size(); ++j) {
            triangles[j][0] += vertexOffset;
            triangles[j][1] += vertexOffset;
            triangles[j][2] += vertexOffset;
        }
        allVertices.insert(allVertices.end(), vertices.begin(), vertices.end());
        allTriangles.insert(allTriangles.end(), triangles.begin(), triangles.end());
        for (plint j = 0; j < (plint)triangles.size(); ++j) {
            triangleTags.push_back(i);
        }
        for (plint j = 0; j < (plint)vertices.size(); ++j) {
            vertexTags.push_back(i);
        }
        vertexOffset += vertices.size();
    }
    RawConnectedTriangleMesh<T> mesh(allVertices, allTriangles, triangleTags, names);

    // Now, RawConnectedTriangleMesh have part-tagging for triangles. But they don't
    // have part-tagging for vertices, which would not be possible in general,
    // as a vertex can belong to triangles in multiple parts.
    // In the present case, we know that each vertex belongs to a unique part,
    // because we performed a disjoint merge.
    // We can therefore create a part-tagging for the vertices manually.
    plint vertexPartTagID = mesh.registerVertexTag(vertexTagForParts);
    plint uniqueNumberingID = mesh.getVertexTag("UniqueID");
    RawConnectedTriangleMesh<double>::PVertexIterator it = mesh.vertexIterator();
    while (!it->end()) {
        RawConnectedTriangleMesh<double>::PVertex vertex = it->next();
        plint id = vertex->tag(uniqueNumberingID);
        vertex->setTag(vertexPartTagID, vertexTags[id]);
    }
    return mesh;
}

}  // namespace plb

#endif  // RAW_CONNECTED_TRIANGLE_MESH_HH
