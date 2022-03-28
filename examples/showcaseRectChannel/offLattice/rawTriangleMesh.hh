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

#ifndef RAW_TRIANGLE_MESH_HH
#define RAW_TRIANGLE_MESH_HH

#include <algorithm>
#include <string>
#include <vector>

#include "core/array.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "io/plbFiles.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "offLattice/offFileIO.h"
#include "offLattice/rawTriangleMesh.h"
#include "offLattice/stlFileIO.h"

namespace plb {

/* *************** Free Functions ******************************************* */

template <typename T>
RawTriangleMesh<T> stlToRawTriangleMesh(STLreader<T> const &reader)
{
    return RawTriangleMesh<T>(reader.getParts(), reader.getNames(), reader.getEps());
}

template <typename T>
RawTriangleMesh<T> stlToRawTriangleMesh(FileName stlFileName, T eps)
{
    return stlToRawTriangleMesh(STLreader<T>(stlFileName.get(), eps));
}

template <typename T>
RawTriangleMesh<T> offToRawTriangleMesh(OFFreader<T> const &reader)
{
    std::vector<Array<T, 3> > const &vertices = reader.getVertices();
    std::vector<std::vector<plint> > const &facets = reader.getFacets();
    typedef Array<Array<T, 3>, 3> RawTriangle;
    std::vector<RawTriangle> triangles(facets.size());

    for (pluint iTriangle = 0; iTriangle < triangles.size(); iTriangle++) {
        PLB_ASSERT(
            facets[iTriangle].size() == 3);  // The surface in the OFF file must be triangulated.
        RawTriangle triangle;
        triangle[0] = vertices[facets[iTriangle][0]];
        triangle[1] = vertices[facets[iTriangle][1]];
        triangle[2] = vertices[facets[iTriangle][2]];
        triangles[iTriangle] = triangle;
    }

    return RawTriangleMesh<T>(triangles);
}

template <typename T>
RawTriangleMesh<T> offToRawTriangleMesh(FileName offFileName)
{
    return offToRawTriangleMesh(OFFreader<T>(offFileName.get()));
}

/* *************** Class RawTriangleMesh::Vertex **************************** */

template <typename T>
RawTriangleMesh<T>::Vertex::Vertex(
    RawTriangleMesh<T> *mesh_, plint part_, plint triangle_, plint vertex_) :
    mesh(mesh_), part(part_), triangle(triangle_), vertex(vertex_)
{ }

template <typename T>
typename RawTriangleMesh<T>::Vertex *RawTriangleMesh<T>::Vertex::clone() const
{
    return new typename RawTriangleMesh<T>::Vertex(*this);
}

template <typename T>
Array<T, 3> const &RawTriangleMesh<T>::Vertex::get() const
{
    return mesh->triangles[part][triangle][vertex];
}

template <typename T>
Array<T, 3> &RawTriangleMesh<T>::Vertex::get()
{
    return mesh->triangles[part][triangle][vertex];
}

template <typename T>
T const &RawTriangleMesh<T>::Vertex::operator[](plint i) const
{
    PLB_ASSERT(i < 3);
    return mesh->triangles[part][triangle][vertex][i];
}

template <typename T>
T &RawTriangleMesh<T>::Vertex::operator[](plint i)
{
    PLB_ASSERT(i < 3);
    return mesh->triangles[part][triangle][vertex][i];
}

template <typename T>
Array<T, 3> RawTriangleMesh<T>::Vertex::normal(bool areaWeighted) const
{
    // Not implemented at the level of raw triangle meshes.
    PLB_ASSERT(false);
    return Array<T, 3>((T)0, (T)0, (T)0);
}

template <typename T>
T RawTriangleMesh<T>::Vertex::area() const
{
    // Not implemented at the level of raw triangle meshes.
    PLB_ASSERT(false);
    return ((T)0);
}

template <typename T>
plint RawTriangleMesh<T>::Vertex::numAdjacentTriangles() const
{
    // Not implemented at the level of raw triangle meshes.
    PLB_ASSERT(false);
    return -1;
}

template <typename T>
typename RawTriangleMesh<T>::CPTriangle RawTriangleMesh<T>::Vertex::adjacentTriangle(
    plint iTriangle) const
{
    // Not implemented at the level of raw triangle meshes.
    PLB_ASSERT(false);
    return typename RawTriangleMesh<T>::CPTriangle(nullptr);
}

template <typename T>
typename RawTriangleMesh<T>::PTriangle RawTriangleMesh<T>::Vertex::adjacentTriangle(plint iTriangle)
{
    // Not implemented at the level of raw triangle meshes.
    PLB_ASSERT(false);
    return typename RawTriangleMesh<T>::PTriangle(nullptr);
}

template <typename T>
bool RawTriangleMesh<T>::Vertex::isInterior() const
{
    // Not implemented at the level of raw triangle meshes.
    PLB_ASSERT(false);
    return false;
}

template <typename T>
std::vector<plint> RawTriangleMesh<T>::Vertex::adjacentVertices() const
{
    // Not implemented at the level of raw triangle meshes.
    PLB_ASSERT(false);
    return std::vector<plint>();
}

template <typename T>
T RawTriangleMesh<T>::Vertex::property(plint whichProperty) const
{
    // Not implemented at the level of raw triangle meshes.
    PLB_ASSERT(false);
    return ((T)0);
}

template <typename T>
void RawTriangleMesh<T>::Vertex::setProperty(plint whichProperty, T value)
{
    // Not implemented at the level of raw triangle meshes.
    PLB_ASSERT(false);
}

template <typename T>
plint RawTriangleMesh<T>::Vertex::tag(plint whichTag) const
{
    // Not implemented at the level of raw triangle meshes.
    PLB_ASSERT(false);
    return -1;
}

template <typename T>
void RawTriangleMesh<T>::Vertex::setTag(plint whichTag, plint value)
{
    // Not implemented at the level of raw triangle meshes.
    PLB_ASSERT(false);
}

/* *************** Class RawTriangleMesh::Triangle **************************** */

template <typename T>
RawTriangleMesh<T>::Triangle::Triangle(RawTriangleMesh<T> *mesh_, plint part_, plint triangle_) :
    mesh(mesh_), part(part_), triangle(triangle_)
{ }

template <typename T>
typename RawTriangleMesh<T>::Triangle *RawTriangleMesh<T>::Triangle::clone() const
{
    return new typename RawTriangleMesh<T>::Triangle(*this);
}

template <typename T>
typename RawTriangleMesh<T>::CPVertex RawTriangleMesh<T>::Triangle::vertex(plint iVertex) const
{
    PLB_ASSERT(iVertex < 3);
    return CPVertex(new Vertex(mesh, part, triangle, iVertex));
}

template <typename T>
typename RawTriangleMesh<T>::PVertex RawTriangleMesh<T>::Triangle::vertex(plint iVertex)
{
    PLB_ASSERT(iVertex < 3);
    return PVertex(new Vertex(mesh, part, triangle, iVertex));
}

template <typename T>
Array<T, 3> const &RawTriangleMesh<T>::Triangle::operator[](plint iVertex) const
{
    RawTriangle const &data = mesh->triangles[part][triangle];
    return data[iVertex];
}

template <typename T>
Array<T, 3> &RawTriangleMesh<T>::Triangle::operator[](plint iVertex)
{
    RawTriangle &data = mesh->triangles[part][triangle];
    return data[iVertex];
}

template <typename T>
T RawTriangleMesh<T>::Triangle::area() const
{
    RawTriangle const &data = mesh->triangles[part][triangle];
    return computeTriangleArea(data[0], data[1], data[2]);
}

template <typename T>
Array<T, 3> RawTriangleMesh<T>::Triangle::normal() const
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
Array<T, 3> RawTriangleMesh<T>::Triangle::normalTimesArea() const
{
    RawTriangle const &data = mesh->triangles[part][triangle];
    bool isAreaWeighted = true;
    return plb::computeTriangleNormal(data[0], data[1], data[2], isAreaWeighted);
}

template <typename T>
Array<T, 3> RawTriangleMesh<T>::Triangle::edgeNormal(plint iEdge, bool areaWeighted) const
{
    // Not implemented at the level of raw triangle meshes.
    PLB_ASSERT(false);
    return Array<T, 3>(T(), T(), T());
}
template <typename T>
Array<T, 3> RawTriangleMesh<T>::Triangle::continuousNormal(
    Array<T, 3> const &p, bool areaWeighted) const
{
    // Not implemented at the level of raw triangle meshes.
    PLB_ASSERT(false);
    return Array<T, 3>(T(), T(), T());
}

template <typename T>
typename RawTriangleMesh<T>::CPTriangle RawTriangleMesh<T>::Triangle::edgeNeighbor(
    plint iEdge) const
{
    // Not implemented at the level of raw triangle meshes.
    PLB_ASSERT(false);
    return CPTriangle(nullptr);
}

template <typename T>
std::vector<plint> RawTriangleMesh<T>::Triangle::edgeNeighbors(plint iEdge) const
{
    // Not implemented at the level of raw triangle meshes.
    PLB_ASSERT(false);
    return std::vector<plint>();
}

template <typename T>
plint RawTriangleMesh<T>::Triangle::numVertexNeighbors(plint iVertex) const
{
    // Not implemented at the level of raw triangle meshes.
    PLB_ASSERT(false);
    return 0;
}

template <typename T>
typename RawTriangleMesh<T>::CPTriangle RawTriangleMesh<T>::Triangle::vertexNeighbor(
    plint iVertex, plint iNeighbor) const
{
    // Not implemented at the level of raw triangle meshes.
    PLB_ASSERT(false);
    return CPTriangle(nullptr);
}

template <typename T>
T RawTriangleMesh<T>::Triangle::property(plint whichProperty) const
{
    // Not implemented at the level of raw triangle meshes.
    PLB_ASSERT(false);
    return ((T)0.0);
}

template <typename T>
void RawTriangleMesh<T>::Triangle::setProperty(plint whichProperty, T value)
{
    // Not implemented at the level of raw triangle meshes.
    PLB_ASSERT(false);
}

template <typename T>
plint RawTriangleMesh<T>::Triangle::tag(plint whichTag) const
{
    if (whichTag == 0) {
        return part;
    } else {
        // Not implemented at the level of raw triangle meshes.
        PLB_ASSERT(false);
        return -1;
    }
}

template <typename T>
void RawTriangleMesh<T>::Triangle::setTag(plint whichTag, plint value)
{
    if (whichTag == 0) {
        // Can't reassign part numbering.
        PLB_ASSERT(false);
    } else {
        // Not implemented at the level of raw triangle meshes.
        PLB_ASSERT(false);
    }
}

/* *************** Class RawTriangleMesh::AllTriangleIterator **************************** */

template <typename T>
RawTriangleMesh<T>::AllTriangleIterator::AllTriangleIterator(RawTriangleMesh<T> *mesh_) :
    mesh(mesh_), currentPart(0), currentTriangle(0)
{ }

template <typename T>
typename RawTriangleMesh<T>::PTriangle RawTriangleMesh<T>::AllTriangleIterator::next()
{
    PLB_ASSERT(currentPart < (plint)mesh->triangles.size());
    PLB_ASSERT(currentTriangle < (plint)mesh->triangles[currentPart].size());

    Triangle *triangle = new Triangle(mesh, currentPart, currentTriangle);
    ++currentTriangle;
    if (currentTriangle >= (plint)mesh->triangles[currentPart].size()) {
        currentTriangle = 0;
        ++currentPart;
    }
    return PTriangle(triangle);
}

template <typename T>
bool RawTriangleMesh<T>::AllTriangleIterator::end() const
{
    return currentPart >= (plint)mesh->triangles.size();
}

template <typename T>
typename RawTriangleMesh<T>::PTriangleIterator RawTriangleMesh<T>::AllTriangleIterator::clone()
    const
{
    return PTriangleIterator(new AllTriangleIterator(*this));
}

/* *************** Class RawTriangleMesh::PartTriangleIterator **************************** */

template <typename T>
RawTriangleMesh<T>::PartTriangleIterator::PartTriangleIterator(
    RawTriangleMesh<T> *mesh_, plint partId_) :
    mesh(mesh_), partId(partId_), currentTriangle(0)
{ }

template <typename T>
typename RawTriangleMesh<T>::PTriangle RawTriangleMesh<T>::PartTriangleIterator::next()
{
    PLB_ASSERT(partId < (plint)mesh->triangles.size());
    PLB_ASSERT(currentTriangle < (plint)mesh->triangles[partId].size());

    Triangle *triangle = new Triangle(mesh, partId, currentTriangle);
    ++currentTriangle;
    return PTriangle(triangle);
}

template <typename T>
bool RawTriangleMesh<T>::PartTriangleIterator::end() const
{
    return currentTriangle >= (plint)mesh->triangles[partId].size();
}

template <typename T>
typename RawTriangleMesh<T>::PTriangleIterator RawTriangleMesh<T>::PartTriangleIterator::clone()
    const
{
    return PTriangleIterator(new PartTriangleIterator(*this));
}

/* *************** Class RawTriangleMesh  ************************************************* */

template <typename T>
RawTriangleMesh<T>::RawTriangleMesh(
    std::vector<typename RawTriangleMesh<T>::RawTriangle> const &triangles_, T eps_) :
    eps(eps_)
{
    // We don't add empty parts, because this would make things more difficult for iterators.
    if (!triangles_.empty()) {
        triangles.push_back(triangles_);
    }
    partNames.push_back("Body");
    signalVertexUpdate();
}

template <typename T>
RawTriangleMesh<T>::RawTriangleMesh(
    std::vector<std::vector<typename RawTriangleMesh<T>::RawTriangle> > const &triangles_,
    std::vector<std::string> const &partNames_, T eps_) :
    eps(eps_)
{
    PLB_ASSERT(triangles_.size() == partNames_.size());
    for (pluint i = 0; i < triangles_.size(); ++i) {
        // We don't add empty parts, because this would make things more difficult for iterators.
        if (!triangles_[i].empty()) {
            triangles.push_back(triangles_[i]);
            partNames.push_back(partNames_[i]);
        }
    }
    signalVertexUpdate();
}

template <typename T>
RawTriangleMesh<T>::RawTriangleMesh(RawTriangleMesh<T> const &rhs) :
    TriangleMesh<T>(rhs),
    triangles(rhs.triangles),
    partNames(rhs.partNames),
    minEdgeLength(rhs.minEdgeLength),
    maxEdgeLength(rhs.maxEdgeLength),
    minTriangleArea(rhs.minTriangleArea),
    maxTriangleArea(rhs.maxTriangleArea),
    boundingCuboid(rhs.boundingCuboid),
    eps(rhs.eps)
{ }

template <typename T>
RawTriangleMesh<T>::RawTriangleMesh(TriangleMesh<T> &rhs, T eps_)
{
    eps = eps_;
    for (plint part = 0; part < rhs.numParts(); ++part) {
        std::vector<RawTriangle> nextTriangles;
        typename TriangleMesh<T>::PTriangleIterator it = rhs.triangleIterator(part);
        while (!it->end()) {
            typename TriangleMesh<T>::PTriangle triangle = it->next();
            nextTriangles.push_back(RawTriangle((*triangle)[0], (*triangle)[1], (*triangle)[2]));
        }
        triangles.push_back(nextTriangles);
        partNames.push_back(rhs.partName(part));
    }
    signalVertexUpdate();
}

template <typename T>
RawTriangleMesh<T> RawTriangleMesh<T>::merge(RawTriangleMesh<T> const &rhs) const
{
    std::vector<std::vector<RawTriangle> > newTriangles(triangles);
    newTriangles.insert(newTriangles.end(), rhs.triangles.begin(), rhs.triangles.end());
    std::vector<std::string> newPartNames;
    for (pluint i = 0; i < partNames.size(); ++i) {
        newPartNames.push_back("mesh1." + partNames[i]);
    }
    for (pluint i = 0; i < rhs.partNames.size(); ++i) {
        newPartNames.push_back("mesh2." + rhs.partNames[i]);
    }
    return RawTriangleMesh(newTriangles, newPartNames);
}

template <typename T>
RawTriangleMesh<T> RawTriangleMesh<T>::refine() const
{
    std::vector<std::vector<RawTriangle> > newTriangles;

    for (pluint part = 0; part < triangles.size(); ++part) {
        std::vector<RawTriangle> nextPart;
        std::vector<RawTriangle> newZeroAreaTriangles;
        for (pluint i = 0; i < triangles[part].size(); ++i) {
            RawTriangle const &triangle = triangles[part][i];

            Array<T, 3> v00 = triangle[0];
            Array<T, 3> v11 = triangle[1];
            Array<T, 3> v22 = triangle[2];

            Array<T, 3> v01 = (T)0.5 * (v00 + v11);
            Array<T, 3> v12 = (T)0.5 * (v11 + v22);
            Array<T, 3> v20 = (T)0.5 * (v22 + v00);

            RawTriangle newTriangle;

            newTriangle = RawTriangle(v01, v12, v20);
            if (triangleHasZeroArea(triangle)) {
                if (!triangleHasZeroLengthEdges(newTriangle)
                    && !containedAndErase(newTriangle, newZeroAreaTriangles))
                {
                    newZeroAreaTriangles.push_back(newTriangle);
                }
            } else {
                nextPart.push_back(newTriangle);
            }

            newTriangle = RawTriangle(v00, v01, v20);
            if (triangleHasZeroArea(triangle)) {
                if (!triangleHasZeroLengthEdges(newTriangle)
                    && !containedAndErase(newTriangle, newZeroAreaTriangles))
                {
                    newZeroAreaTriangles.push_back(newTriangle);
                }
            } else {
                nextPart.push_back(newTriangle);
            }

            newTriangle = RawTriangle(v01, v11, v12);
            if (triangleHasZeroArea(triangle)) {
                if (!triangleHasZeroLengthEdges(newTriangle)
                    && !containedAndErase(newTriangle, newZeroAreaTriangles))
                {
                    newZeroAreaTriangles.push_back(newTriangle);
                }
            } else {
                nextPart.push_back(newTriangle);
            }

            newTriangle = RawTriangle(v20, v12, v22);
            if (triangleHasZeroArea(triangle)) {
                if (!triangleHasZeroLengthEdges(newTriangle)
                    && !containedAndErase(newTriangle, newZeroAreaTriangles))
                {
                    newZeroAreaTriangles.push_back(newTriangle);
                }
            } else {
                nextPart.push_back(newTriangle);
            }
        }
        nextPart.insert(nextPart.end(), newZeroAreaTriangles.begin(), newZeroAreaTriangles.end());
        newTriangles.push_back(nextPart);
    }

    return RawTriangleMesh<T>(newTriangles, partNames);
}

template <typename T>
RawTriangleMesh<T> RawTriangleMesh<T>::refine(T edgeLengthThreshold) const
{
    PLB_ASSERT(util::greaterThan(edgeLengthThreshold, (T)0, eps));

    T thr2 = util::sqr(edgeLengthThreshold);

    std::vector<std::vector<RawTriangle> > allNewTriangles;

    for (pluint part = 0; part < triangles.size(); ++part) {
        std::vector<RawTriangle> newTriangles;
        std::vector<RawTriangle> newZeroAreaTriangles;

        RawTriangle newTriangle;

        for (pluint i = 0; i < triangles[part].size(); ++i) {
            RawTriangle const &triangle = triangles[part][i];

            Array<T, 3> v0 = triangle[0];
            Array<T, 3> v1 = triangle[1];
            Array<T, 3> v2 = triangle[2];

            T e0 = normSqr<T, 3>(v1 - v0);
            T e1 = normSqr<T, 3>(v2 - v1);
            T e2 = normSqr<T, 3>(v0 - v2);

            if (e0 >= thr2) {
                Array<T, 3> vA = (T)0.5 * (v0 + v1);
                if (e1 >= thr2) {
                    Array<T, 3> vB = (T)0.5 * (v1 + v2);
                    if (e2 >= thr2) {
                        Array<T, 3> vC = (T)0.5 * (v2 + v0);

                        newTriangle[0] = v0;
                        newTriangle[1] = vA;
                        newTriangle[2] = vC;
                        if (triangleHasZeroArea(triangle)) {
                            if (!triangleHasZeroLengthEdges(newTriangle)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }

                        newTriangle[0] = vA;
                        newTriangle[1] = v1;
                        newTriangle[2] = vB;
                        if (triangleHasZeroArea(triangle)) {
                            if (!triangleHasZeroLengthEdges(newTriangle)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }

                        newTriangle[0] = vC;
                        newTriangle[1] = vB;
                        newTriangle[2] = v2;
                        if (triangleHasZeroArea(triangle)) {
                            if (!triangleHasZeroLengthEdges(newTriangle)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }

                        newTriangle[0] = vA;
                        newTriangle[1] = vB;
                        newTriangle[2] = vC;
                        if (triangleHasZeroArea(triangle)) {
                            if (!triangleHasZeroLengthEdges(newTriangle)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }
                    } else {
                        newTriangle[0] = vA;
                        newTriangle[1] = v1;
                        newTriangle[2] = vB;
                        if (triangleHasZeroArea(triangle)) {
                            if (!triangleHasZeroLengthEdges(newTriangle)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }

                        T eA2 = normSqr<T, 3>(vA - v2);
                        T eB0 = normSqr<T, 3>(vB - v0);
                        if (eA2 < eB0) {
                            newTriangle[0] = v2;
                            newTriangle[1] = v0;
                            newTriangle[2] = vA;
                            if (triangleHasZeroArea(triangle)) {
                                if (!triangleHasZeroLengthEdges(newTriangle)
                                    && !containedAndErase(newTriangle, newZeroAreaTriangles))
                                {
                                    newZeroAreaTriangles.push_back(newTriangle);
                                }
                            } else {
                                newTriangles.push_back(newTriangle);
                            }

                            newTriangle[0] = v2;
                            newTriangle[1] = vA;
                            newTriangle[2] = vB;
                            if (triangleHasZeroArea(triangle)) {
                                if (!triangleHasZeroLengthEdges(newTriangle)
                                    && !containedAndErase(newTriangle, newZeroAreaTriangles))
                                {
                                    newZeroAreaTriangles.push_back(newTriangle);
                                }
                            } else {
                                newTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangle[0] = v0;
                            newTriangle[1] = vA;
                            newTriangle[2] = vB;
                            if (triangleHasZeroArea(triangle)) {
                                if (!triangleHasZeroLengthEdges(newTriangle)
                                    && !containedAndErase(newTriangle, newZeroAreaTriangles))
                                {
                                    newZeroAreaTriangles.push_back(newTriangle);
                                }
                            } else {
                                newTriangles.push_back(newTriangle);
                            }

                            newTriangle[0] = v0;
                            newTriangle[1] = vB;
                            newTriangle[2] = v2;
                            if (triangleHasZeroArea(triangle)) {
                                if (!triangleHasZeroLengthEdges(newTriangle)
                                    && !containedAndErase(newTriangle, newZeroAreaTriangles))
                                {
                                    newZeroAreaTriangles.push_back(newTriangle);
                                }
                            } else {
                                newTriangles.push_back(newTriangle);
                            }
                        }
                    }
                } else {
                    if (e2 >= thr2) {
                        Array<T, 3> vC = (T)0.5 * (v2 + v0);

                        newTriangle[0] = v0;
                        newTriangle[1] = vA;
                        newTriangle[2] = vC;
                        if (triangleHasZeroArea(triangle)) {
                            if (!triangleHasZeroLengthEdges(newTriangle)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }

                        T eA2 = normSqr<T, 3>(vA - v2);
                        T eC1 = normSqr<T, 3>(vC - v1);
                        if (eA2 < eC1) {
                            newTriangle[0] = v2;
                            newTriangle[1] = vC;
                            newTriangle[2] = vA;
                            if (triangleHasZeroArea(triangle)) {
                                if (!triangleHasZeroLengthEdges(newTriangle)
                                    && !containedAndErase(newTriangle, newZeroAreaTriangles))
                                {
                                    newZeroAreaTriangles.push_back(newTriangle);
                                }
                            } else {
                                newTriangles.push_back(newTriangle);
                            }

                            newTriangle[0] = v2;
                            newTriangle[1] = vA;
                            newTriangle[2] = v1;
                            if (triangleHasZeroArea(triangle)) {
                                if (!triangleHasZeroLengthEdges(newTriangle)
                                    && !containedAndErase(newTriangle, newZeroAreaTriangles))
                                {
                                    newZeroAreaTriangles.push_back(newTriangle);
                                }
                            } else {
                                newTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangle[0] = v1;
                            newTriangle[1] = v2;
                            newTriangle[2] = vC;
                            if (triangleHasZeroArea(triangle)) {
                                if (!triangleHasZeroLengthEdges(newTriangle)
                                    && !containedAndErase(newTriangle, newZeroAreaTriangles))
                                {
                                    newZeroAreaTriangles.push_back(newTriangle);
                                }
                            } else {
                                newTriangles.push_back(newTriangle);
                            }

                            newTriangle[0] = v1;
                            newTriangle[1] = vC;
                            newTriangle[2] = vA;
                            if (triangleHasZeroArea(triangle)) {
                                if (!triangleHasZeroLengthEdges(newTriangle)
                                    && !containedAndErase(newTriangle, newZeroAreaTriangles))
                                {
                                    newZeroAreaTriangles.push_back(newTriangle);
                                }
                            } else {
                                newTriangles.push_back(newTriangle);
                            }
                        }
                    } else {
                        newTriangle[0] = v2;
                        newTriangle[1] = v0;
                        newTriangle[2] = vA;
                        if (triangleHasZeroArea(triangle)) {
                            if (!triangleHasZeroLengthEdges(newTriangle)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }

                        newTriangle[0] = v2;
                        newTriangle[1] = vA;
                        newTriangle[2] = v1;
                        if (triangleHasZeroArea(triangle)) {
                            if (!triangleHasZeroLengthEdges(newTriangle)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }
                    }
                }
            } else {
                if (e1 >= thr2) {
                    Array<T, 3> vB = (T)0.5 * (v1 + v2);
                    if (e2 >= thr2) {
                        Array<T, 3> vC = (T)0.5 * (v2 + v0);

                        newTriangle[0] = v2;
                        newTriangle[1] = vC;
                        newTriangle[2] = vB;
                        if (triangleHasZeroArea(triangle)) {
                            if (!triangleHasZeroLengthEdges(newTriangle)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }

                        T eB0 = normSqr<T, 3>(vB - v0);
                        T eC1 = normSqr<T, 3>(vC - v1);
                        if (eB0 < eC1) {
                            newTriangle[0] = v0;
                            newTriangle[1] = v1;
                            newTriangle[2] = vB;
                            if (triangleHasZeroArea(triangle)) {
                                if (!triangleHasZeroLengthEdges(newTriangle)
                                    && !containedAndErase(newTriangle, newZeroAreaTriangles))
                                {
                                    newZeroAreaTriangles.push_back(newTriangle);
                                }
                            } else {
                                newTriangles.push_back(newTriangle);
                            }

                            newTriangle[0] = v0;
                            newTriangle[1] = vB;
                            newTriangle[2] = vC;
                            if (triangleHasZeroArea(triangle)) {
                                if (!triangleHasZeroLengthEdges(newTriangle)
                                    && !containedAndErase(newTriangle, newZeroAreaTriangles))
                                {
                                    newZeroAreaTriangles.push_back(newTriangle);
                                }
                            } else {
                                newTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangle[0] = v1;
                            newTriangle[1] = vC;
                            newTriangle[2] = v0;
                            if (triangleHasZeroArea(triangle)) {
                                if (!triangleHasZeroLengthEdges(newTriangle)
                                    && !containedAndErase(newTriangle, newZeroAreaTriangles))
                                {
                                    newZeroAreaTriangles.push_back(newTriangle);
                                }
                            } else {
                                newTriangles.push_back(newTriangle);
                            }

                            newTriangle[0] = v1;
                            newTriangle[1] = vB;
                            newTriangle[2] = vC;
                            if (triangleHasZeroArea(triangle)) {
                                if (!triangleHasZeroLengthEdges(newTriangle)
                                    && !containedAndErase(newTriangle, newZeroAreaTriangles))
                                {
                                    newZeroAreaTriangles.push_back(newTriangle);
                                }
                            } else {
                                newTriangles.push_back(newTriangle);
                            }
                        }
                    } else {
                        newTriangle[0] = v0;
                        newTriangle[1] = v1;
                        newTriangle[2] = vB;
                        if (triangleHasZeroArea(triangle)) {
                            if (!triangleHasZeroLengthEdges(newTriangle)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }

                        newTriangle[0] = v0;
                        newTriangle[1] = vB;
                        newTriangle[2] = v2;
                        if (triangleHasZeroArea(triangle)) {
                            if (!triangleHasZeroLengthEdges(newTriangle)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }
                    }
                } else {
                    if (e2 >= thr2) {
                        Array<T, 3> vC = (T)0.5 * (v2 + v0);

                        newTriangle[0] = v1;
                        newTriangle[1] = vC;
                        newTriangle[2] = v0;
                        if (triangleHasZeroArea(triangle)) {
                            if (!triangleHasZeroLengthEdges(newTriangle)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }

                        newTriangle[0] = v1;
                        newTriangle[1] = v2;
                        newTriangle[2] = vC;
                        if (triangleHasZeroArea(triangle)) {
                            if (!triangleHasZeroLengthEdges(newTriangle)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }
                    } else {
                        newTriangle[0] = v0;
                        newTriangle[1] = v1;
                        newTriangle[2] = v2;
                        if (triangleHasZeroArea(triangle)) {
                            if (!triangleHasZeroLengthEdges(newTriangle)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }
                    }
                }
            }
        }

        newTriangles.insert(
            newTriangles.end(), newZeroAreaTriangles.begin(), newZeroAreaTriangles.end());

        allNewTriangles.push_back(newTriangles);
    }

    return RawTriangleMesh<T>(allNewTriangles, partNames);
}

template <typename T>
RawTriangleMesh<T> RawTriangleMesh<T>::refineRecursively(
    T targetMaxEdgeLength, plint maxNumIterations, bool &success) const
{
    PLB_ASSERT(util::greaterThan(targetMaxEdgeLength, (T)0, eps));
    PLB_ASSERT(maxNumIterations > 0);

    RawTriangleMesh<T> result(*this);
    for (plint iter = 0; iter < maxNumIterations; ++iter) {
        result = result.refine(targetMaxEdgeLength);
        if (result.getMaxEdgeLength() < targetMaxEdgeLength)
            break;
    }

    success = result.getMaxEdgeLength() < targetMaxEdgeLength;
    return result;
}

template <typename T>
RawTriangleMesh<T> RawTriangleMesh<T>::refineByArea(T triangleAreaThreshold) const
{
    PLB_ASSERT(util::greaterThan(triangleAreaThreshold, (T)0, eps));

    std::vector<std::vector<RawTriangle> > newTriangles;

    for (pluint part = 0; part < triangles.size(); ++part) {
        std::vector<RawTriangle> nextPart;
        std::vector<RawTriangle> newZeroAreaTriangles;
        for (pluint i = 0; i < triangles[part].size(); ++i) {
            RawTriangle const &triangle = triangles[part][i];
            T area = computeTriangleArea(triangle[0], triangle[1], triangle[2]);
            if (area >= triangleAreaThreshold) {
                Array<T, 3> v00 = triangle[0];
                Array<T, 3> v11 = triangle[1];
                Array<T, 3> v22 = triangle[2];

                Array<T, 3> v01 = (T)0.5 * (v00 + v11);
                Array<T, 3> v12 = (T)0.5 * (v11 + v22);
                Array<T, 3> v20 = (T)0.5 * (v22 + v00);

                RawTriangle newTriangle;

                newTriangle = RawTriangle(v01, v12, v20);
                if (triangleHasZeroArea(triangle)) {
                    if (!triangleHasZeroLengthEdges(newTriangle)
                        && !containedAndErase(newTriangle, newZeroAreaTriangles))
                    {
                        newZeroAreaTriangles.push_back(newTriangle);
                    }
                } else {
                    nextPart.push_back(newTriangle);
                }

                newTriangle = RawTriangle(v00, v01, v20);
                if (triangleHasZeroArea(triangle)) {
                    if (!triangleHasZeroLengthEdges(newTriangle)
                        && !containedAndErase(newTriangle, newZeroAreaTriangles))
                    {
                        newZeroAreaTriangles.push_back(newTriangle);
                    }
                } else {
                    nextPart.push_back(newTriangle);
                }

                newTriangle = RawTriangle(v01, v11, v12);
                if (triangleHasZeroArea(triangle)) {
                    if (!triangleHasZeroLengthEdges(newTriangle)
                        && !containedAndErase(newTriangle, newZeroAreaTriangles))
                    {
                        newZeroAreaTriangles.push_back(newTriangle);
                    }
                } else {
                    nextPart.push_back(newTriangle);
                }

                newTriangle = RawTriangle(v20, v12, v22);
                if (triangleHasZeroArea(triangle)) {
                    if (!triangleHasZeroLengthEdges(newTriangle)
                        && !containedAndErase(newTriangle, newZeroAreaTriangles))
                    {
                        newZeroAreaTriangles.push_back(newTriangle);
                    }
                } else {
                    nextPart.push_back(newTriangle);
                }
            } else {
                nextPart.push_back(triangle);
            }
        }
        nextPart.insert(nextPart.end(), newZeroAreaTriangles.begin(), newZeroAreaTriangles.end());
        newTriangles.push_back(nextPart);
    }

    return RawTriangleMesh<T>(newTriangles, partNames);
}

template <typename T>
RawTriangleMesh<T> RawTriangleMesh<T>::refineByAreaRecursively(
    T targetMaxTriangleArea, plint maxNumIterations, bool &success) const
{
    PLB_ASSERT(util::greaterThan(targetMaxTriangleArea, (T)0, eps));
    PLB_ASSERT(maxNumIterations > 0);

    RawTriangleMesh<T> result(*this);
    for (plint iter = 0; iter < maxNumIterations; ++iter) {
        result = result.refineByArea(targetMaxTriangleArea);
        if (result.getMaxTriangleArea() < targetMaxTriangleArea)
            break;
    }

    success = result.getMaxTriangleArea() < targetMaxTriangleArea;
    return result;
}

template <typename T>
RawTriangleMesh<T> RawTriangleMesh<T>::select(TriangleSelector<T> const &selector) const
{
    std::vector<std::vector<RawTriangle> > newTriangles;

    for (pluint part = 0; part < triangles.size(); ++part) {
        std::vector<RawTriangle> nextPart;
        for (pluint i = 0; i < triangles[part].size(); ++i) {
            RawTriangle const &triangle = triangles[part][i];

            if (selector(triangle, part)) {
                nextPart.push_back(triangle);
            }
        }
        newTriangles.push_back(nextPart);
    }

    return RawTriangleMesh<T>(newTriangles, partNames);
}

template <typename T>
RawTriangleMesh<T> RawTriangleMesh<T>::cutWithPlane(Array<T, 3> planePos, Array<T, 3> normal) const
{
    std::vector<std::vector<RawTriangle> > newTriangles;
    std::vector<std::string> newPartNames;
    T norm_normal = norm(normal);
    PLB_ASSERT(
        !util::isZero(norm_normal));  // The cut plane normal vector cannot have zero magnitude.
    normal /= norm_normal;
    for (pluint iPart = 0; iPart < triangles.size(); ++iPart) {
        std::vector<RawTriangle> tmpTriangles;
        std::vector<RawTriangle> newZeroAreaTriangles;
        for (pluint iTriangle = 0; iTriangle < triangles[iPart].size(); iTriangle++) {
            cutTriangleWithPlane(
                triangles[iPart][iTriangle], planePos, normal, tmpTriangles, newZeroAreaTriangles);
        }
        tmpTriangles.insert(
            tmpTriangles.end(), newZeroAreaTriangles.begin(), newZeroAreaTriangles.end());
        if (!tmpTriangles.empty()) {
            newTriangles.push_back(tmpTriangles);
            newPartNames.push_back(partNames[iPart]);
        }
    }
    return RawTriangleMesh<T>(newTriangles, newPartNames);
}

template <typename T>
bool RawTriangleMesh<T>::triangleHasZeroArea(RawTriangle const &triangle) const
{
    T area = computeTriangleArea(triangle[0], triangle[1], triangle[2]);
    if (util::isZero(area, eps)) {
        return true;
    }
    return false;
}

template <typename T>
bool RawTriangleMesh<T>::triangleHasZeroLengthEdges(RawTriangle const &triangle) const
{
    if (util::isZero(norm(triangle[1] - triangle[0]), eps)
        || util::isZero(norm(triangle[2] - triangle[0]), eps)
        || util::isZero(norm(triangle[2] - triangle[1]), eps))
    {
        return true;
    }
    return false;
}

template <typename T>
bool RawTriangleMesh<T>::areTheSameTriangle(RawTriangle const &t1, RawTriangle const &t2) const
{
    Array<T, 3> const &v10 = t1[0];
    Array<T, 3> const &v11 = t1[1];
    Array<T, 3> const &v12 = t1[2];

    Array<T, 3> const &v20 = t2[0];
    Array<T, 3> const &v21 = t2[1];
    Array<T, 3> const &v22 = t2[2];

    T n1020 = norm(v10 - v20);
    T n1021 = norm(v10 - v21);
    T n1022 = norm(v10 - v22);
    if (!util::isZero(n1020, eps) && !util::isZero(n1021, eps) && !util::isZero(n1022, eps)) {
        return false;
    }

    T n1120 = norm(v11 - v20);
    T n1121 = norm(v11 - v21);
    T n1122 = norm(v11 - v22);
    if (!util::isZero(n1120, eps) && !util::isZero(n1121, eps) && !util::isZero(n1122, eps)) {
        return false;
    }

    T n1220 = norm(v12 - v20);
    T n1221 = norm(v12 - v21);
    T n1222 = norm(v12 - v22);
    if (!util::isZero(n1220, eps) && !util::isZero(n1221, eps) && !util::isZero(n1222, eps)) {
        return false;
    }

    return true;
}

template <typename T>
bool RawTriangleMesh<T>::containedAndErase(
    RawTriangle const &triangle, std::vector<RawTriangle> &triangles) const
{
    typename std::vector<RawTriangle>::iterator it;
    for (it = triangles.begin(); it != triangles.end(); ++it) {
        if (areTheSameTriangle(triangle, *it)) {
            triangles.erase(it);
            return true;
        }
    }
    return false;
}

template <typename T>
bool RawTriangleMesh<T>::cutTriangleWithPlane(
    RawTriangle const &triangle, Array<T, 3> const &planePos, Array<T, 3> const &normal,
    std::vector<RawTriangle> &newTriangles, std::vector<RawTriangle> &newZeroAreaTriangles) const
{
    int vertexTags[3];

    // Tag the triangle vertices.
    for (int iVertex = 0; iVertex < 3; iVertex++) {
        Array<T, 3> tmp = triangle[iVertex] - planePos;
        T norm_tmp = norm(tmp);
        if (!util::isZero(norm_tmp)) {
            tmp /= norm_tmp;
        } else {
            tmp[0] = tmp[1] = tmp[2] = (T)0.0;
        }
        T dotp = dot(tmp, normal);
        if (util::isZero(dotp, eps)) {
            vertexTags[iVertex] = 0;
        } else if (util::greaterThan(dotp, (T)0, eps)) {
            vertexTags[iVertex] = -1;
        } else if (util::lessThan(dotp, (T)0, eps)) {
            vertexTags[iVertex] = 1;
        } else {
            return false;
        }
    }

    // All three vertices belong to one side of the cut plane.
    if (vertexTags[0] == 1 && vertexTags[1] == 1 && vertexTags[2] == 1) {
        if (triangleHasZeroArea(triangle)) {
            if (!triangleHasZeroLengthEdges(triangle)
                && !containedAndErase(triangle, newZeroAreaTriangles)) {
                newZeroAreaTriangles.push_back(triangle);
            }
        } else {
            newTriangles.push_back(triangle);
        }
        return true;
    } else if (vertexTags[0] == -1 && vertexTags[1] == -1 && vertexTags[2] == -1) {
        return false;
    }

    // One vertex belongs to one side of the cut plane and the other two vertices
    //   belong to the other side.
    for (int i = 0; i < 3; i++) {
        int j = (i + 1) != 3 ? (i + 1) : 0;
        int k = (j + 1) != 3 ? (j + 1) : 0;

        if (vertexTags[i] == 1 && vertexTags[j] == -1 && vertexTags[k] == -1) {
            Array<T, 3> intersection_ij((T)0.0, (T)0.0, (T)0.0),
                intersection_ik((T)0.0, (T)0.0, (T)0.0);
            int rv = 0;
            rv = lineIntersectionWithPlane<T>(
                Plane<T>(planePos, normal), triangle[i], triangle[j], eps, intersection_ij);
            if (rv != 1) {
                return false;
            }
            rv = lineIntersectionWithPlane<T>(
                Plane<T>(planePos, normal), triangle[i], triangle[k], eps, intersection_ik);
            if (rv != 1) {
                return false;
            }
            RawTriangle newTriangle;
            newTriangle[0] = triangle[i];
            newTriangle[1] = intersection_ij;
            newTriangle[2] = intersection_ik;
            if (triangleHasZeroArea(triangle)) {
                if (!triangleHasZeroLengthEdges(newTriangle)
                    && !containedAndErase(newTriangle, newZeroAreaTriangles))
                {
                    newZeroAreaTriangles.push_back(newTriangle);
                }
            } else {
                newTriangles.push_back(newTriangle);
            }
            return true;
        } else if (vertexTags[i] == -1 && vertexTags[j] == 1 && vertexTags[k] == 1) {
            Array<T, 3> intersection_ij((T)0.0, (T)0.0, (T)0.0),
                intersection_ik((T)0.0, (T)0.0, (T)0.0);
            int rv = 0;
            rv = lineIntersectionWithPlane<T>(
                Plane<T>(planePos, normal), triangle[i], triangle[j], eps, intersection_ij);
            if (rv != 1) {
                return false;
            }
            rv = lineIntersectionWithPlane<T>(
                Plane<T>(planePos, normal), triangle[i], triangle[k], eps, intersection_ik);
            if (rv != 1) {
                return false;
            }
            RawTriangle newTriangle_0(triangle[k], intersection_ij, triangle[j]);
            RawTriangle newTriangle_1(triangle[k], intersection_ik, intersection_ij);
            if (triangleHasZeroArea(triangle)) {
                if (!triangleHasZeroLengthEdges(newTriangle_0)
                    && !containedAndErase(newTriangle_0, newZeroAreaTriangles))
                {
                    newZeroAreaTriangles.push_back(newTriangle_0);
                }
                if (!triangleHasZeroLengthEdges(newTriangle_1)
                    && !containedAndErase(newTriangle_1, newZeroAreaTriangles))
                {
                    newZeroAreaTriangles.push_back(newTriangle_1);
                }
            } else {
                newTriangles.push_back(newTriangle_0);
                newTriangles.push_back(newTriangle_1);
            }
            return true;
        }
    }

    // Only one vertex belongs to the cut plane.
    for (int i = 0; i < 3; i++) {
        int j = (i + 1) != 3 ? (i + 1) : 0;
        int k = (j + 1) != 3 ? (j + 1) : 0;

        if (vertexTags[i] == 0) {
            if (vertexTags[j] == 1 && vertexTags[k] == 1) {
                if (triangleHasZeroArea(triangle)) {
                    if (!triangleHasZeroLengthEdges(triangle)
                        && !containedAndErase(triangle, newZeroAreaTriangles)) {
                        newZeroAreaTriangles.push_back(triangle);
                    }
                } else {
                    newTriangles.push_back(triangle);
                }
                return true;
            } else if (vertexTags[j] == -1 && vertexTags[k] == -1) {
                return false;
            } else if (vertexTags[j] == 1 && vertexTags[k] == -1) {
                Array<T, 3> intersection((T)0.0, (T)0.0, (T)0.0);
                int rv = 0;
                rv = lineIntersectionWithPlane<T>(
                    Plane<T>(planePos, normal), triangle[j], triangle[k], eps, intersection);
                if (rv != 1) {
                    return false;
                }
                RawTriangle newTriangle(triangle[i], triangle[j], intersection);
                if (triangleHasZeroArea(triangle)) {
                    if (!triangleHasZeroLengthEdges(newTriangle)
                        && !containedAndErase(newTriangle, newZeroAreaTriangles))
                    {
                        newZeroAreaTriangles.push_back(newTriangle);
                    }
                } else {
                    newTriangles.push_back(newTriangle);
                }
                return true;
            } else if (vertexTags[j] == -1 && vertexTags[k] == 1) {
                Array<T, 3> intersection((T)0.0, (T)0.0, (T)0.0);
                int rv = 0;
                rv = lineIntersectionWithPlane<T>(
                    Plane<T>(planePos, normal), triangle[j], triangle[k], eps, intersection);
                if (rv != 1) {
                    return false;
                }
                RawTriangle newTriangle(triangle[i], intersection, triangle[k]);
                if (triangleHasZeroArea(triangle)) {
                    if (!triangleHasZeroLengthEdges(newTriangle)
                        && !containedAndErase(newTriangle, newZeroAreaTriangles))
                    {
                        newZeroAreaTriangles.push_back(newTriangle);
                    }
                } else {
                    newTriangles.push_back(newTriangle);
                }
                return true;
            }
        }
    }

    // Only two of the three vertices belong to the cut plane.
    for (int i = 0; i < 3; i++) {
        int j = (i + 1) != 3 ? (i + 1) : 0;
        int k = (j + 1) != 3 ? (j + 1) : 0;

        if (vertexTags[i] == 0 && vertexTags[j] == 0) {
            if (vertexTags[k] == 1) {
                if (triangleHasZeroArea(triangle)) {
                    if (!triangleHasZeroLengthEdges(triangle)
                        && !containedAndErase(triangle, newZeroAreaTriangles)) {
                        newZeroAreaTriangles.push_back(triangle);
                    }
                } else {
                    newTriangles.push_back(triangle);
                }
                return true;
            } else if (vertexTags[k] == -1) {
                return false;
            }
        }
    }

    // All 3 vertices belong to the cut plane.
    if (vertexTags[0] == 0 && vertexTags[1] == 0 && vertexTags[2] == 0) {
        if (triangleHasZeroArea(triangle)) {
            if (!triangleHasZeroLengthEdges(triangle)
                && !containedAndErase(triangle, newZeroAreaTriangles)) {
                newZeroAreaTriangles.push_back(triangle);
            }
        } else {
            newTriangles.push_back(triangle);
        }
        return true;
    }

    return false;
}

template <typename T>
typename RawTriangleMesh<T>::PTriangleIterator RawTriangleMesh<T>::triangleIterator(plint partId)
{
    if (partId == -1) {
        return PTriangleIterator(new AllTriangleIterator(this));
    } else {
        PLB_ASSERT(partId < numParts());
        return PTriangleIterator(new PartTriangleIterator(this, partId));
    }
}

template <typename T>
typename RawTriangleMesh<T>::PVertexIterator RawTriangleMesh<T>::vertexIterator()
{
    // Not implemented for vertices.
    PLB_ASSERT(false);
    return typename RawTriangleMesh<T>::PVertexIterator(nullptr);
}

template <typename T>
plint RawTriangleMesh<T>::numParts() const
{
    return (plint)triangles.size();
}

template <typename T>
std::string RawTriangleMesh<T>::partName(plint iPart) const
{
    PLB_ASSERT(iPart < (plint)partNames.size());
    return partNames[iPart];
}

template <typename T>
plint RawTriangleMesh<T>::partId(std::string partName) const
{
    std::vector<std::string>::const_iterator it =
        find(partNames.begin(), partNames.end(), partName);
    if (it == partNames.end()) {
        return -1;
    } else {
        return it - partNames.begin();
    }
}

template <typename T>
void RawTriangleMesh<T>::signalVertexUpdate()
{
    computeMinMaxEdges();
    computeMinMaxAreas();
    computeBoundingCuboid();
}

template <typename T>
void RawTriangleMesh<T>::signalIsometricVertexUpdate()
{
    computeBoundingCuboid();
}

template <typename T>
void RawTriangleMesh<T>::translate(Array<T, 3> const &vector)
{
    for (pluint part = 0; part < triangles.size(); ++part) {
        for (pluint i = 0; i < triangles[part].size(); i++) {
            for (plint j = 0; j < 3; j++) {
                triangles[part][i][j] += vector;
            }
        }
    }
    boundingCuboid.lowerLeftCorner += vector;
    boundingCuboid.upperRightCorner += vector;
}

template <typename T>
void RawTriangleMesh<T>::scale(T alpha)
{
    for (pluint part = 0; part < triangles.size(); ++part) {
        for (pluint i = 0; i < triangles[part].size(); i++) {
            for (plint j = 0; j < 3; j++) {
                triangles[part][i][j] *= alpha;
            }
        }
    }
    computeMinMaxEdges();
    computeMinMaxAreas();
    computeBoundingCuboid();
}

template <typename T>
void RawTriangleMesh<T>::rotateAtOrigin(Array<T, 3> const &normedAxis, T theta)
{
    for (pluint part = 0; part < triangles.size(); ++part) {
        for (pluint i = 0; i < triangles[part].size(); i++) {
            for (plint j = 0; j < 3; j++) {
                triangles[part][i][j] =
                    plb::rotateAtOrigin(triangles[part][i][j], normedAxis, theta);
            }
        }
    }

    computeBoundingCuboid();
}

template <typename T>
void RawTriangleMesh<T>::rotate(T phi, T theta, T psi)
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

    for (pluint part = 0; part < triangles.size(); ++part) {
        for (pluint iTriangle = 0; iTriangle < triangles[part].size(); iTriangle++) {
            for (pluint iVertex = 0; iVertex < 3; iVertex++) {
                Array<T, 3> x = triangles[part][iTriangle][iVertex];
                for (int i = 0; i < 3; i++) {
                    triangles[part][iTriangle][iVertex][i] = (T)0.0;
                    for (int j = 0; j < 3; j++) {
                        triangles[part][iTriangle][iVertex][i] += a[i][j] * x[j];
                    }
                }
            }
        }
    }

    computeBoundingCuboid();
}

template <typename T>
void RawTriangleMesh<T>::reverseOrientation()
{
    for (pluint part = 0; part < triangles.size(); ++part) {
        for (pluint i = 0; i < triangles[part].size(); i++) {
            std::swap(triangles[part][i][1], triangles[part][i][2]);
        }
    }
}

template <typename T>
T RawTriangleMesh<T>::getMinEdgeLength() const
{
    return minEdgeLength;
}

template <typename T>
T RawTriangleMesh<T>::getMaxEdgeLength() const
{
    return maxEdgeLength;
}

template <typename T>
T RawTriangleMesh<T>::getMinTriangleArea() const
{
    return minTriangleArea;
}

template <typename T>
T RawTriangleMesh<T>::getMaxTriangleArea() const
{
    return maxTriangleArea;
}

template <typename T>
Cuboid<T> RawTriangleMesh<T>::getBoundingCuboid() const
{
    return boundingCuboid;
}

template <typename T>
std::vector<std::vector<typename RawTriangleMesh<T>::RawTriangle> > const &
    RawTriangleMesh<T>::getTriangles() const
{
    return triangles;
}

template <typename T>
plint RawTriangleMesh<T>::getTriangleTag(std::string tagName) const
{
    if (tagName == "Part") {
        return 0;
    }
    return -1;
}
template <typename T>
std::string RawTriangleMesh<T>::getTriangleTagName(plint tag) const
{
    PLB_ASSERT(tag == 0);
    return "Part";
}

template <typename T>
plint RawTriangleMesh<T>::numTriangleTags() const
{
    return 1;
}

template <typename T>
void RawTriangleMesh<T>::computeMinMaxEdges()
{
    minEdgeLength = std::numeric_limits<T>::max();
    maxEdgeLength = std::numeric_limits<T>::min();
    for (pluint part = 0; part < triangles.size(); ++part) {
        for (pluint i = 0; i < triangles[part].size(); i++) {
            RawTriangle const &triangle = triangles[part][i];
            T edge1 = norm(triangle[1] - triangle[0]);
            T edge2 = norm(triangle[2] - triangle[1]);
            T edge3 = norm(triangle[0] - triangle[2]);
            T minEdge = std::min(edge1, std::min(edge2, edge3));
            T maxEdge = std::max(edge1, std::max(edge2, edge3));
            minEdgeLength = std::min(minEdgeLength, minEdge);
            maxEdgeLength = std::max(maxEdgeLength, maxEdge);
        }
    }
}

template <typename T>
void RawTriangleMesh<T>::computeMinMaxAreas()
{
    minTriangleArea = std::numeric_limits<T>::max();
    maxTriangleArea = std::numeric_limits<T>::min();
    for (pluint part = 0; part < triangles.size(); ++part) {
        for (pluint i = 0; i < triangles[part].size(); i++) {
            RawTriangle const &triangle = triangles[part][i];
            T area = computeTriangleArea(triangle[0], triangle[1], triangle[2]);
            minTriangleArea = std::min(minTriangleArea, area);
            maxTriangleArea = std::max(maxTriangleArea, area);
        }
    }
}

template <typename T>
void RawTriangleMesh<T>::computeBoundingCuboid()
{
    T xMin, yMin, zMin;
    T xMax, yMax, zMax;

    xMin = yMin = zMin = std::numeric_limits<T>::max();
    xMax = yMax = zMax = -std::numeric_limits<T>::max();
    for (pluint part = 0; part < triangles.size(); ++part) {
        for (pluint i = 0; i < triangles[part].size(); i++) {
            RawTriangle const &triangle = triangles[part][i];

            xMin =
                std::min(xMin, std::min(triangle[0][0], std::min(triangle[1][0], triangle[2][0])));
            yMin =
                std::min(yMin, std::min(triangle[0][1], std::min(triangle[1][1], triangle[2][1])));
            zMin =
                std::min(zMin, std::min(triangle[0][2], std::min(triangle[1][2], triangle[2][2])));

            xMax =
                std::max(xMax, std::max(triangle[0][0], std::max(triangle[1][0], triangle[2][0])));
            yMax =
                std::max(yMax, std::max(triangle[0][1], std::max(triangle[1][1], triangle[2][1])));
            zMax =
                std::max(zMax, std::max(triangle[0][2], std::max(triangle[1][2], triangle[2][2])));
        }
    }
    boundingCuboid.lowerLeftCorner = Array<T, 3>(xMin, yMin, zMin);
    boundingCuboid.upperRightCorner = Array<T, 3>(xMax, yMax, zMax);
}

template <typename T>
plint RawTriangleMesh<T>::getNumTriangles() const
{
    plint numTriangles = 0;
    for (pluint part = 0; part < triangles.size(); ++part) {
        numTriangles += triangles[part].size();
    }
    return numTriangles;
}

}  // namespace plb

#endif  // RAW_TRIANGLE_MESH_HH
