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

#ifndef RAW_CONNECTED_TRIANGLE_MESH_H
#define RAW_CONNECTED_TRIANGLE_MESH_H

#include <string>
#include <vector>

#include "core/array.h"
#include "core/globalDefs.h"
#include "io/plbFiles.h"
#include "offLattice/offFileIO.h"
#include "offLattice/rawTriangleMesh.h"
#include "offLattice/stlFileIO.h"
#include "offLattice/triangleSelector.h"

namespace plb {

template <typename T>
class RawConnectedTriangleMesh;

template <typename T>
RawConnectedTriangleMesh<T> stlToConnectedTriangleMesh(STLreader<T> const &reader);

template <typename T>
RawConnectedTriangleMesh<T> stlToConnectedTriangleMesh(
    FileName stlFileName, T eps = getEpsilon<T>(DBL));

template <typename T>
RawConnectedTriangleMesh<T> offToConnectedTriangleMesh(OFFreader<T> const &reader);

template <typename T>
RawConnectedTriangleMesh<T> offToConnectedTriangleMesh(FileName offFileName);

template <typename T>
RawConnectedTriangleMesh<T> generateConnectedTriangleMesh(
    RawTriangleMesh<T> &rawMesh, T eps = getEpsilon<T>(DBL));

template <typename T>
class ManualConnectedTriangleMesh : public ConnectedTriangleMesh<T> {
public:
    /// In this constructor, it is assumed that there is only one part in
    /// the geometry, and that the tagging is therefore trivial.
    ManualConnectedTriangleMesh(plint numTriangles_, plint numVertices_);
    /// This constructor explicitly specifies the tagging for the geometry parts.
    ManualConnectedTriangleMesh(
        plint numTriangles_, plint numVertices_, std::vector<plint> partTagging,
        std::vector<std::string> nameOfParts_);

    virtual plint registerTriangleTag(std::string tagName);
    virtual plint getTriangleTag(std::string tagName) const;
    virtual std::string getTriangleTagName(plint tag) const;
    virtual plint numTriangleTags() const;
    virtual void clearTriangleTags();

    virtual plint registerTriangleProperty(std::string propertyName);
    virtual plint getTriangleProperty(std::string propertyName) const;
    virtual std::string getTrianglePropertyName(plint property) const;
    virtual plint numTriangleProperties() const;
    virtual void clearTriangleProperties();

    virtual plint registerVertexTag(std::string tagName);
    virtual plint getVertexTag(std::string tagName) const;
    virtual std::string getVertexTagName(plint tag) const;
    virtual plint numVertexTags() const;
    virtual void clearVertexTags();

    virtual plint registerVertexProperty(std::string propertyName);
    virtual plint getVertexProperty(std::string propertyName) const;
    virtual std::string getVertexPropertyName(plint property) const;
    virtual plint numVertexProperties() const;
    virtual void clearVertexProperties();

    virtual plint numParts() const;
    virtual std::string partName(plint iPart) const;
    virtual plint partId(std::string partName) const;

    virtual plint getNumTriangles() const;
    virtual plint getNumVertices() const;

protected:
    plint numTriangles, numVertices;
    std::vector<std::vector<plint> >
        trianglesPerPart;  /// Aggregates the index of all triangles
                           /// defined on a given part. We do this only for
                           /// the parts. None of the other triangle
                           /// properties have the luxury of having a
                           /// property-to-triangle map.
                           /// Note: trianglesPerPart.size()==nameOfParts.size(), except if
                           /// there's only one part. In that case trianglesPerPart is empty.
    std::vector<std::vector<plint> >
        triangleTags;  /// Note: triangleTags.size()==triangleTagNames.size() and
                       /// triangleTags[i].size()==triangles.size().
    std::vector<std::vector<T> >
        triangleProperties;  /// Note: triangleProperties.size()==trianglePropertyNames.size() and
                             /// triangleProperties[i].size()==triangles.size().
    std::vector<std::vector<plint> > vertexTags;  /// Note: vertexTags.size()==vertexTagNames.size()
                                                  /// and vertexTags[i].size()==vertices.size().
    std::vector<std::vector<T> >
        vertexProperties;  /// Note: vertexProperties.size()==vertexPropertyNames.size() and
                           /// vertexProperties[i].size()==vertices.size().
    std::vector<std::string> triangleTagNames, trianglePropertyNames;
    std::vector<std::string> vertexTagNames, vertexPropertyNames;
    std::vector<std::string> nameOfParts;
};

template <typename T>
class RawConnectedTriangleMesh : public ManualConnectedTriangleMesh<T> {
public:
    class Triangle;
    typedef Array<Array<T, 3>, 3> RawTriangle;
    typedef typename TriangleMesh<T>::PTriangle PTriangle;
    typedef typename TriangleMesh<T>::CPTriangle CPTriangle;
    typedef typename TriangleMesh<T>::PVertex PVertex;
    typedef typename TriangleMesh<T>::CPVertex CPVertex;
    typedef typename TriangleMesh<T>::PTriangleIterator PTriangleIterator;
    typedef typename TriangleMesh<T>::PVertexIterator PVertexIterator;

    class Vertex : public TriangleMesh<T>::Vertex {
    public:
        Vertex(RawConnectedTriangleMesh<T> *mesh_, plint vertex_);
        virtual Vertex *clone() const;
        virtual Array<T, 3> const &get() const;      /// Raw coordinates.
        virtual Array<T, 3> &get();                  /// Raw coordinates.
        virtual T const &operator[](plint i) const;  /// Access elements.
        virtual T &operator[](plint i);              /// Access elements.
        /// Compute the normal vector at a given vertex. If "areaWeighted" is false, then
        /// the vertex normal is calculated as the simple normalized sum of the individual
        /// unit normals of the triangles which share the specific vertex. If
        /// "areaWeighted" is true, then the vertex normal is calculated as the
        /// normalized sum of the individual area weighted normals of the triangles which
        /// share the vertex under consideration. An area weighted normal of a given triangle
        /// is the normal that has length equal to twice the area of the triangle.
        virtual Array<T, 3> normal(bool areaWeighted = false) const;
        virtual T area() const;  /// One third of the sum of the areas of the neighboring triangles.
        virtual plint numAdjacentTriangles() const;  /// How many triangles share this vertex?
        virtual CPTriangle adjacentTriangle(
            plint iTriangle) const;                           /// Get any of the adjacent triangles.
        virtual PTriangle adjacentTriangle(plint iTriangle);  /// Get any of the adjacent triangles.
        virtual T property(
            plint whichProperty) const;  /// Get a vertex property (to know its name, ask the mesh).
        virtual void setProperty(plint whichProperty, T value);
        virtual plint tag(
            plint whichTag) const;  /// Get a vertex tag (to know its name, ask the mesh).
                                    /// To get the tag for the unique vertex numbering,
                                    /// one must use getVertexTag("UniqueID")).
        virtual void setTag(plint whichTag, plint value);
        virtual bool isInterior()
            const;  /// Check if the vertex is an interior, non-boundary vertex.
        virtual std::vector<plint> adjacentVertices()
            const;  /// Return a list of adjacent vertices.
    private:
        mutable RawConnectedTriangleMesh<T> *mesh;
        plint vertex;
    };

    class Triangle : public TriangleMesh<T>::Triangle {
    public:
        Triangle(RawConnectedTriangleMesh<T> *mesh_, plint triangle_);
        virtual Triangle *clone() const;
        virtual CPVertex vertex(plint iVertex) const;
        virtual PVertex vertex(plint iVertex);
        virtual Array<T, 3> const &operator[](plint iVertex) const;
        virtual Array<T, 3> &operator[](plint iVertex);
        virtual T area() const;
        virtual Array<T, 3> normal() const;
        /// An area weighted normal of a given triangle is the normal that has length equal
        /// to twice the area of the triangle.
        virtual Array<T, 3> normalTimesArea() const;
        /// If an edge has two neighboring triangles, this function returns the average of
        /// the normals of these two. If it has only one, ore more than two triangle, the
        /// function returns the normal of the current triangle.
        virtual Array<T, 3> edgeNormal(plint iEdge, bool areaWeighted = false) const;
        /// Compute a "continuous" normal at a point "p". The computed normal vector is
        ///   normalized (has length one) and is continuous as its base point moves across
        ///   the triangles. The point "p" must belong to the interior or to the boundary of
        ///   the triangle. The currently implemented algorithm uses a barycentric
        ///   coordinate representation of the normal with respect to the three vertex
        ///   normals of the given triangle. If "isAreaWeighted" is true, then the vertex
        ///   normals used in the barycentric representation are area weighted.
        ///   If the triangle has zero area, then the normal is the zero vector.
        virtual Array<T, 3> continuousNormal(Array<T, 3> const &p, bool areaWeighted = false) const;
        /// Returns the neighboring triangle that meets on the edge with local index "iEdge".
        /// "iNeighbor" must be greater equal to 0, and less than the size of edgeNeighbors(iEdge).
        /// Returns a Null pointer if the number of neighboring triangles is zero (the edge
        /// is a boundary edge) or if the number of neighboring triangles is larger than 1
        /// (the edge is non-manifold). In the latter case, you must use the method
        /// "egdeNeighbors" to get all triangles.
        virtual CPTriangle edgeNeighbor(plint iEdge) const;
        /// For a manifold mesh with a boundary, there exist 0 or 1 edge-neighbor triangles
        /// for a given edge. For a non-manifold mesh, there can exist more. This function
        /// returns the global ID of all triangles that meet on a given edge. You can convert
        /// the global ID into a triangle with the mesh method "triangle".
        virtual std::vector<plint> edgeNeighbors(plint iEdge) const;
        /// Returns the number of triangles that meet on the vertex with the specified
        /// local index.
        virtual plint numVertexNeighbors(plint iVertex) const;
        /// Returns one of the triangles that meet on the vertex with local index "iVertex".
        /// "iNeighbor" must be greater equal to 0, and less than numVertexNeighbors(iVertex).
        /// If not, a NULL pointer is returned.
        virtual CPTriangle vertexNeighbor(plint iVertex, plint iNeighbor) const;
        virtual T property(plint whichProperty) const;  /// Any property (color, etc.).
        virtual void setProperty(plint whichProperty, T value);
        virtual plint tag(plint whichTag) const;  /// Get a triangle tag. Predefined tags are:
                                                  /// To get the tag for the part,
                                                  /// one must use getTriangleTag("Part"),
                                                  /// and to get the tag for unique triangle
                                                  /// numbering, one must use
                                                  /// getTriangleTag("UniqueID").
        virtual void setTag(plint whichTag, plint value);

    private:
        RawConnectedTriangleMesh<T> *mesh;
        plint triangle;
    };
    /// Iterate over all triangles, forward only.
    class AllTriangleIterator : public TriangleMesh<T>::TriangleIterator {
    public:
        AllTriangleIterator(RawConnectedTriangleMesh<T> *mesh_);
        virtual PTriangle next();
        virtual bool end() const;
        virtual PTriangleIterator clone() const;

    private:
        RawConnectedTriangleMesh<T> *mesh;
        plint currentPart;
        plint currentTriangle;
    };
    /// Iterate over all triangles, forward only.
    class PartTriangleIterator : public TriangleMesh<T>::TriangleIterator {
    public:
        PartTriangleIterator(RawConnectedTriangleMesh<T> *mesh, plint partId_);
        virtual PTriangle next();
        virtual bool end() const;
        virtual PTriangleIterator clone() const;

    private:
        RawConnectedTriangleMesh<T> *mesh;
        plint partId;
        plint currentTriangle;
    };
    /// Iterate over all vertices, forward only.
    class VertexIterator : public TriangleMesh<T>::VertexIterator {
    public:
        VertexIterator(RawConnectedTriangleMesh<T> *mesh_);
        virtual PVertex next();
        virtual bool end() const;
        virtual PVertexIterator clone() const;

    private:
        RawConnectedTriangleMesh<T> *mesh;
        plint currentVertex;
    };

public:
    /// Constructor for single part, with raw triangles.
    RawConnectedTriangleMesh(
        std::vector<Array<T, 3> > const &vertices_,
        std::vector<Array<plint, 3> > const &triangles_);
    RawConnectedTriangleMesh(
        std::vector<Array<T, 3> > const &vertices_, std::vector<Array<plint, 3> > const &triangles_,
        std::vector<plint> const &partTagging, std::vector<std::string> const &nameOfParts);

    RawConnectedTriangleMesh<T> merge(
        RawConnectedTriangleMesh<T> &mesh, T eps = getEpsilon<T>(DBL));
    RawConnectedTriangleMesh<T> refine(T eps = getEpsilon<T>(DBL)) const;
    RawConnectedTriangleMesh<T> refineRecursively(
        T targetMaxEdgeLength, plint maxNumIterations, bool &success,
        T eps = getEpsilon<T>(DBL)) const;
    RawConnectedTriangleMesh<T> select(
        TriangleSelector<T> const &selector, T eps = getEpsilon<T>(DBL)) const;
    RawConnectedTriangleMesh<T> cutWithPlane(
        Array<T, 3> planePos, Array<T, 3> normal, T eps = getEpsilon<T>(DBL)) const;

    virtual PTriangleIterator triangleIterator(plint partId = -1);
    virtual PTriangle triangle(plint uniqueID);
    virtual PVertexIterator vertexIterator();
    virtual PVertex vertex(plint uniqueID);

    virtual void signalVertexUpdate();
    virtual void signalIsometricVertexUpdate();
    virtual void translate(Array<T, 3> const &vector);
    virtual void scale(T alpha);
    virtual void rotateAtOrigin(Array<T, 3> const &normedAxis, T theta);
    virtual void rotate(T phi, T theta, T psi);
    virtual void reverseOrientation();
    virtual T getMinEdgeLength() const;
    virtual T getMaxEdgeLength() const;
    virtual T getMinTriangleArea() const;
    virtual T getMaxTriangleArea() const;
    virtual Cuboid<T> getBoundingCuboid() const;
    std::vector<Array<plint, 3> > const &getTriangles() const
    {
        return triangles;
    }
    std::vector<Array<T, 3> > const &getVertices() const
    {
        return vertices;
    }
    std::vector<std::vector<plint> > const &getTrianglesOnVertex() const
    {
        return trianglesOnVertex;
    }
    void verticesToFile(FileName fileName) const;
    void verticesFromFile(FileName fileName);

private:
    void computeTrianglesOnVertex();
    void computeMinMaxEdges() const;     // Has lazy evaluation.
    void computeMinMaxAreas() const;     // Has lazy evaluation.
    void computeBoundingCuboid() const;  // Has lazy evaluation.
private:
    /**** This part is for the actual implementation of the mesh *************/
    std::vector<Array<T, 3> > vertices;  /// Positions of vertices. The vector
                                         /// index is the global vertex index.
    std::vector<Array<plint, 3> >
        triangles;  /// Global indices of the vertices that constitute the
                    /// triangle. The vector index is the global triangle index.
    std::vector<std::vector<plint> > trianglesOnVertex;  /// Global indices of the triangles that
                                                         /// share a given vertex.
    /// Note: vertices.size()==trianglesOnVertex.size()
    mutable bool minMaxEdgeLengthValid;  /// Lazy evaluation for min-max edge-length.
    mutable T minEdgeLength, maxEdgeLength;
    mutable bool minMaxTriangleAreaValid;  /// Lazy evaluation for min-max triangle-area.
    mutable T minTriangleArea, maxTriangleArea;
    mutable bool boundingCuboidValid;  /// Lazy evaluation for bounding-cuboid.
    mutable Cuboid<T> boundingCuboid;
};

/// Merges multiple meshes into a single one.
//  It is assumed that all input meshes are disconnected. Coinciding vertices will not be merged.
//  If you need to merge vertices, use RawConnectedTriangleMesh::merge() instead.
//
//  The partwise subdivision of every input mesh is ignored. Every input mesh is treated as
//  a single part.
template <typename T>
RawConnectedTriangleMesh<T> disconnectedMerge(
    std::vector<RawConnectedTriangleMesh<T> const *> parts, std::vector<std::string> names,
    std::string vertexTagForParts = "Part");

}  // namespace plb

#endif  // RAW_CONNECTED_TRIANGLE_MESH_H
