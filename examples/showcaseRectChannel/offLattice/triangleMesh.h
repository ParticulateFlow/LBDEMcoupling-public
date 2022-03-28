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

#ifndef TRIANGLE_MESH_H
#define TRIANGLE_MESH_H

#include <memory>
#include <vector>

#include "core/array.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"

/* Virtual bases for triangular mesh classes. The basic hierarchy is:
 *
 * TriangleMesh  --> Set of triangles without explicit connectivity
 *       |           information (like an STL file).
 *       |
 * ConnectedTriangleMesh  --> These meshes explicitly identify vertices
 *       |                    that are shared among triangles, but impose
 *       |                    nothing on orientation or anything else
 *       |                    (like an X3D or OFF file).
 *       |
 * ManifoldTriangleMesh  --> In this case, the surface must be orientable
 *                           and manifold. It can for example be represented
 *                           by a Directed-Edge format (DEF).
 *
 * These classes have the following fundamentally useful properties:
 * - They offer a generic way of using meshes of different degree of
 *   sophistication. For example, if you write an algorithm that works for
 *   a naive set of triangles, it will also work for a DEF mesh.
 * - They have a built-in mechanism for distringuishing parts of a geometry.
 * - Additionally to the parts mechanism, they have an arbitrary number
 *   of tags and properties.
 *
 * Iterators
 * =========
 * All meshes have simple iterators for looping over all triangles or
 * over all vertices. It is not possible to loop over arbitrary sub-sets,
 * but you can loop over the geometry parts individually, by requesting
 * a corresponding iterator.
 *
 * The TriangleMesh can loop only over triangles, and not over vertices.
 * Apart from that, all meshes have fully functional iterators. The
 * triangles and vertices returned by the iterators provide additional
 * functionalities (depending on the type of mesh), such as triangle
 * area or adjacent triangles.
 *
 * Tags and properties
 * ===================
 * All triangles and vertices (except for TriangleMesh which only offers
 * this ability for triangles) can have an unlimited number of associated
 * tags and properties (and all of them have names defined by a string).
 *
 * Properties are real-valued numbers which could be a color or a factor
 * of rigidity.
 *
 * Tags are integer-valued, and you can use them to assign a more complex
 * property (e.g. a boundary condition) to triangles and vertices. The
 * mesh classes provide no further support for doing so: you must by
 * your own means create a tag-to-property map, outside the mesh class.
 *
 * All vertices have a predefined tag "UniqueID". It is read-only:
 * redefining or overwriting it has no effect.
 * All triangles have two predefined tags: "Part" and "UniqueID". The are
 * read-only: redefining or overwriting them has no effect.
 *
 * Here's an example on how to get the tag "UniqueID" of a vertex:
 * plint vertexIDtagging = part.getVertexTag("UniqueID"); // Do this once for all.
 * PVertex vertex = ...
 * plint tag = vertex->tag(vertexIDtagging);
 *
 * Unique numering of connected meshes
 * ===================================
 * In a connected mesh, all vertices and all triangles have a predefined tag
 * "UniqueID", which provides a unique, numbering of the vertices/triangles.
 * This provides a convenient way of referencing vertices and triangles without
 * the heavy artillery of iterators.
 * You can get a triangle or vertex from their unique ID through the methods
 * - PTriangle triangle(plint uniqueID)
 * - Vertex vertex(plint uniqueID)
 * Non-connected, raw meshes don't offer a unique numbering. That's because
 * the RawTriangleMesh is multi-part, and each part has its own unique numbering,
 * which however is not globally unique: that would be too much of a hassle.
 * Attention: The IDs are unique but not necessarily continuous. This is not a
 * valid way to iterate through the triangles or vertices. Use iterators
 * instead.
 *
 * Parts
 * =====
 * One set of triangle-tags is always created by default. It is called
 * "Part" and labels different parts of the geometry, as they are
 * for example obtained from a CAD file. They can be used to assign
 * boundary conditions. Some useful information
 * - When requiring a triangleIterator, you can either iterate over the full
 *   mesh ( mesh.triangleIterator() ), or over a given part
 *   ( mesh.triangleIterator(partId) ).
 * - The RawConnectedTriangleMesh has a constructor which explicitly accepts
 *   multiple parts, defined through the vector "partTagging".
 * - You can extract a part and create a new mesh from it, using the function
 *   "extractConnectedPart".
 * - Attention: It is not guaranteed that the mesh is the union of the parts,
 *   or that the parts don't overlap. To iterate over the full mesh, you must
 *   request the corresponding iterator, and not iterate over all parts.
 *
 * How to convert meshes
 * =====================
 * - TriangleSet is the raw triangle mesh format from the old Palabos. It is equivalent
 *   to RawTriangleMesh, but it is mutable. This makes a crucial speed difference
 *   when you add up pieces into a mesh (e.g. create a streamline out of multiple
 *   spheres). Therefore, we keep TriangleSet in our framework for efficiency reasons,
 *   but the default thing to use should be RawTriangleMesh.
 * - TriangleBoundary is the DEF format from the old Palabos (obsolete).
 * - RawTriangleMesh is the new raw format.
 * - RawConnectedTriangleMesh is the new connected format.
 *
 * - TriangleBoundary to RawConnectedTriangleMesh.
 *   Use def_to_ConnectedMesh from obsoleteFormatWrapper.h:
 *   RawConnectedTriangleMesh<T> connectedMesh = def_to_ConnectedMesh(boundary);
 *   This generates a multi-part connected mesh, interpreting the "lids" in the DEF
 *   mesh as parts.
 *
 * - TriangleSet to RawTriangleMesh.
 *   Use triangleSetToRawTriangleMesh from obsoleteFormatWrapper.h:
 *   RawTriangleMesh<T> rawMesh = triangleSetToRawTriangleMesh(triangleSet);
 *
 * - TriangleSet to RawConnectedTriangleMesh.
 *   Use triangleSetToConnectedTriangleMesh from obsoleteFormatWrapper.h:
 *
 * - RawTriangleMesh to RawConnectedTriangleMesh:
 *   Use generateConnectedTriangleMesh from connectedTriangleMesh.h
 *
 * - RawConnectedTriangleMesh to RawTriangleMesh:
 *   Use the any-mesh to RawMesh constructor:
 *   RawTriangleMesh<T> rawMesh(connectedMesh);
 *
 * - RawTriangleMesh to TriangleSet:
 *   Use rawTriangleMeshToTriangleSet from obsoleteFormatWrapper:
 *   TriangleSet<T> triangleSet(rawTriangleMeshToTriangleSet(triangleMesh));
 *
 * How to read meshes
 * ==================
 *
 * - RawTriangleMesh:
 *      RawTriangleMesh<T> stlToRawTriangleMesh(STLreader<T> const& reader);
 *      RawTriangleMesh<T> stlToRawTriangleMesh(FileName stlFileName, T eps=getEpsilon<T>(DBL));
 *      RawTriangleMesh<T> offToRawTriangleMesh(OFFreader<T> const& reader);
 *      RawTriangleMesh<T> offToRawTriangleMesh(FileName otlFileName);
 *
 * - RawConnectedTriangleMesh:
 *      RawConnectedTriangleMesh<T> stlToConnectedTriangleMesh(STLreader<T> const& reader);
 *      RawConnectedTriangleMesh<T> stlToConnectedTriangleMesh(FileName stlFileName, T
 * eps=getEpsilon<T>(DBL)); RawConnectedTriangleMesh<T> offToConnectedTriangleMesh(OFFreader<T>
 * const& reader); RawConnectedTriangleMesh<T> offToConnectedTriangleMesh(FileName offFileName);
 */

namespace plb {

/// Set of triangles without explicit connectivity information (like an STL file).
template <typename T>
struct TriangleMesh {
public:
    struct Triangle;
    struct Vertex;
    typedef std::unique_ptr<Triangle> PTriangle;
    typedef std::unique_ptr<const Triangle> CPTriangle;
    typedef std::unique_ptr<Vertex> PVertex;
    typedef std::unique_ptr<const Vertex> CPVertex;
    /// A smart vertex with connectivity information.
    struct Vertex {
        virtual ~Vertex() { }
        virtual Vertex *clone() const = 0;
        virtual Array<T, 3> const &get() const = 0;      /// Raw coordinates.
        virtual Array<T, 3> &get() = 0;                  /// Raw coordinates.
        virtual T const &operator[](plint i) const = 0;  /// Access elements.
        virtual T &operator[](plint i) = 0;              /// Access elements.
        /// Compute the normal vector at a given vertex. If "areaWeighted" is false, then
        /// the vertex normal is calculated as the simple normalized sum of the individual
        /// unit normals of the triangles which share the specific vertex. If
        /// "areaWeighted" is true, then the vertex normal is calculated as the
        /// normalized sum of the individual area weighted normals of the triangles which
        /// share the vertex under consideration. An area weighted normal of a given triangle
        /// is the normal that has length equal to twice the area of the triangle.
        virtual Array<T, 3> normal(bool areaWeighted = false) const = 0;
        virtual T area()
            const = 0;  /// One third of the sum of the areas of the neighboring triangles.
        virtual plint numAdjacentTriangles() const = 0;  /// How many triangles share this vertex?
        virtual CPTriangle adjacentTriangle(
            plint iTriangle) const = 0;  /// Get any of the adjacent triangles.
        virtual PTriangle adjacentTriangle(
            plint iTriangle) = 0;  /// Get any of the adjacent triangles.
        virtual bool isInterior()
            const = 0;  /// Check if the vertex is an interior, non-boundary vertex.
        virtual std::vector<plint> adjacentVertices()
            const = 0;  /// Return a list of adjacent vertices.
        virtual T property(plint whichProperty)
            const = 0;  /// Get a vertex property (to know its name, ask the mesh).
        virtual void setProperty(plint whichProperty, T value) = 0;
        virtual plint tag(
            plint whichTag) const = 0;  /// Get a vertex tag (to know its name, ask the mesh).
                                        /// To get the tag for the unique vertex numbering,
                                        /// one must use getVertexTag("UniqueID")).
        virtual void setTag(plint whichTag, plint value) = 0;
    };
    /// A smart triangle with info like area and normal,
    /// and connectivity information to neighbouring triangles.
    struct Triangle {
        virtual ~Triangle() { }
        virtual Triangle *clone() const = 0;
        virtual CPVertex vertex(plint iVertex) const = 0;
        virtual PVertex vertex(plint iVertex) = 0;
        virtual Array<T, 3> const &operator[](plint iVertex) const = 0;
        virtual Array<T, 3> &operator[](plint iVertex) = 0;
        virtual T area() const = 0;
        virtual Array<T, 3> normal() const = 0;
        /// An area weighted normal of a given triangle is the normal that has length equal
        /// to twice the area of the triangle.
        virtual Array<T, 3> normalTimesArea() const = 0;
        /// If an edge has two neighboring triangles, this function returns the average of
        /// the normals of these two. If it has only one, ore more than two triangle, the
        /// function returns the normal of the current triangle.
        virtual Array<T, 3> edgeNormal(plint iEdge, bool areaWeighted = false) const = 0;
        /// Compute a "continuous" normal at a point "p". The computed normal vector is
        ///   normalized (has length one) and is continuous as its base point moves across
        ///   the triangles. The point "p" must belong to the interior or to the boundary of
        ///   the triangle. The currently implemented algorithm uses a barycentric
        ///   coordinate representation of the normal with respect to the three vertex
        ///   normals of the given triangle. If "isAreaWeighted" is true, then the vertex
        ///   normals used in the barycentric representation are area weighted.
        ///   If the triangle has zero area, then the normal is the zero vector.
        virtual Array<T, 3> continuousNormal(
            Array<T, 3> const &p, bool areaWeighted = false) const = 0;
        /// Returns the neighboring triangle that meets on the edge with local index "iEdge".
        /// "iNeighbor" must be greater equal to 0, and less than the size of edgeNeighbors(iEdge).
        /// Returns a Null pointer if the number of neighboring triangles is zero (the edge
        /// is a boundary edge) or if the number of neighboring triangles is larger than 1
        /// (the edge is non-manifold). In the latter case, you must use the method
        /// "egdeNeighbors" to get all triangles.
        virtual CPTriangle edgeNeighbor(plint iEdge) const = 0;
        /// For a manifold mesh with a boundary, there exist 0 or 1 edge-neighbor triangles
        /// for a given edge. For a non-manifold mesh, there can exist more. This function
        /// returns the global ID of all triangles that meet on a given edge. You can convert
        /// the global ID into a triangle with the mesh method "triangle".
        virtual std::vector<plint> edgeNeighbors(plint iEdge) const = 0;
        /// Returns the number of triangles that meet on the vertex with the specified
        /// local index.
        virtual plint numVertexNeighbors(plint iVertex) const = 0;
        /// Returns one of the triangles that meet on the vertex with local index "iVertex".
        /// "iNeighbor" must be greater equal to 0, and less than numVertexNeighbors(iVertex).
        /// If not, a NULL pointer is returned.
        virtual CPTriangle vertexNeighbor(plint iVertex, plint iNeighbor) const = 0;
        virtual T property(plint whichProperty) const = 0;  /// Any property (color, etc.).
        virtual void setProperty(plint whichProperty, T value) = 0;
        virtual plint tag(plint whichTag) const = 0;  /// Any tag (to get the tag for the part,
                                                      /// one must use getTriangleTag("Part")).
        virtual void setTag(plint whichTag, plint value) = 0;
    };
    struct TriangleIterator;
    typedef std::unique_ptr<TriangleIterator> PTriangleIterator;
    /// Iterate over all triangles, forward only.
    struct TriangleIterator {
        virtual ~TriangleIterator() { }
        virtual PTriangle next() = 0;
        virtual bool end() const = 0;
        virtual PTriangleIterator clone() const = 0;
    };
    struct VertexIterator;
    typedef std::unique_ptr<VertexIterator> PVertexIterator;
    /// Iterate over all vertices, forward only.
    struct VertexIterator {
        virtual ~VertexIterator() { }
        virtual PVertex next() = 0;
        virtual bool end() const = 0;
        virtual PVertexIterator clone() const = 0;
    };

public:
    virtual ~TriangleMesh() { }
    /// Get a forward iterator for triangles of the full mesh (if partId==-1)
    /// or of a part of the mesh (if partId>=0).
    virtual PTriangleIterator triangleIterator(plint partId = -1) = 0;
    /// Get a forward iterator for vertices of the full mesh. You can't
    /// get it for single parts, as you can for triangles.
    virtual PVertexIterator vertexIterator() = 0;
    /// Get number of parts in the geometry (they are always numbered
    /// continuously from 0 to numParts()-1.
    virtual plint numParts() const = 0;
    /// Get the name of any of the parts.
    virtual std::string partName(plint iPart) const = 0;
    /// Get the ID of any of the parts.
    virtual plint partId(std::string partName) const = 0;
    /// If an external function updates the vertex positions, it needs to
    /// call signalVertexUpdate() so the mesh can update its internal state.
    virtual void signalVertexUpdate() = 0;
    /// If an external function updates the vertex positions, but preserves the
    /// lengths between vertices, signalIsometricVertexUpdate() is more efficient than
    /// (but has equivalent result to) signalVertexUpdate().
    virtual void signalIsometricVertexUpdate() = 0;
    /// Translate all vertices of all parts by constant offset.
    virtual void translate(Array<T, 3> const &vector) = 0;
    /// Scale all vertices of all parts by constant factor.
    virtual void scale(T alpha) = 0;
    /// Rotate all parts of the geometry by a given angle, around an
    /// axis that passes through the origin. If you want to play with
    /// an arbitrary axis, do a translate-rotate-translate cycle.
    virtual void rotateAtOrigin(Array<T, 3> const &normedAxis, T theta) = 0;
    /// Rotate the triangular surface mesh.
    ///   The arguments of this function are the Euler angles in radians.
    ///   The so-called "x-convention" is used, in which the rotation is
    ///   given by the three angles (phi, theta, psi), where:
    ///   1.] The first rotation is by an angle phi about the z-axis,
    ///   2.] The second rotation is by an angle theta in [0, pi] about
    ///       the new x-axis,
    ///   3.] The third rotation is by an angle psi about the new z-axis.
    virtual void rotate(T phi, T theta, T psi) = 0;
    /// Reverse the orientation of all triangles in all parts.
    virtual void reverseOrientation() = 0;
    /// Return the length of the shortest edge (including all parts).
    virtual T getMinEdgeLength() const = 0;
    /// Return the length of the longest edge (including all parts).
    virtual T getMaxEdgeLength() const = 0;
    /// Return the minimum triangle area (including all parts).
    virtual T getMinTriangleArea() const = 0;
    /// Return the maximum triangle area (including all parts).
    virtual T getMaxTriangleArea() const = 0;
    /// Return the bounding box of the geometry (including all parts).
    virtual Cuboid<T> getBoundingCuboid() const = 0;

    virtual plint getTriangleTag(std::string tagName) const = 0;
    virtual std::string getTriangleTagName(plint tag) const = 0;
    virtual plint numTriangleTags() const = 0;

    virtual plint getNumTriangles() const = 0;
};

/// These meshes explicitly identifies vertices that are shared among
/// triangles, but imposes nothing on orientation or anything else
/// (just like an X3D or OFF file).
/// At this level, we start implementing tagging and property naming
/// (this didn't seem worth it at the TriangleMesh level, which had
/// no other ability than naming of parts).
template <typename T>
class ConnectedTriangleMesh : public TriangleMesh<T> {
public:
    typedef typename TriangleMesh<T>::PTriangle PTriangle;
    typedef typename TriangleMesh<T>::PVertex PVertex;

public:
    /// Get a triangle, as defined through its unique ID (corresponding
    /// to the tag "UniqueID" of the triangle).
    virtual PTriangle triangle(plint uniqueID) = 0;
    /// Get a vertex, as defined through its unique ID (corresponding
    /// to the tag "UniqueID" of the vertex).
    virtual PVertex vertex(plint uniqueID) = 0;

    virtual plint registerTriangleTag(std::string tagName) = 0;
    virtual plint getTriangleTag(std::string tagName) const = 0;
    virtual std::string getTriangleTagName(plint tag) const = 0;
    virtual plint numTriangleTags() const = 0;
    virtual void clearTriangleTags() = 0;

    virtual plint registerTriangleProperty(std::string propertyName) = 0;
    virtual plint getTriangleProperty(std::string propertyName) const = 0;
    virtual std::string getTrianglePropertyName(plint property) const = 0;
    virtual plint numTriangleProperties() const = 0;
    virtual void clearTriangleProperties() = 0;

    virtual plint registerVertexTag(std::string tagName) = 0;
    virtual plint getVertexTag(std::string tagName) const = 0;
    virtual std::string getVertexTagName(plint tag) const = 0;
    virtual plint numVertexTags() const = 0;
    virtual void clearVertexTags() = 0;

    virtual plint registerVertexProperty(std::string propertyName) = 0;
    virtual plint getVertexProperty(std::string propertyName) const = 0;
    virtual std::string getVertexPropertyName(plint property) const = 0;
    virtual plint numVertexProperties() const = 0;
    virtual void clearVertexProperties() = 0;

    virtual plint getNumTriangles() const = 0;
    virtual plint getNumVertices() const = 0;
};

/// In these meshes case, the surface must be orientable and manifold.
/// It can for example be represented by a Directed-Edge format.
template <typename T>
struct ManifoldTriangleMesh : public ConnectedTriangleMesh<T> {
};

/*
template<typename T>
class DEFtriangleMesh : public ManifoldTriangleMesh<T> {
};
*/

}  // namespace plb

#endif  // TRIANGLE_MESH_H
