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

#ifndef RAW_TRIANGLE_MESH_H
#define RAW_TRIANGLE_MESH_H

#include <string>
#include <vector>

#include "core/array.h"
#include "core/globalDefs.h"
#include "io/plbFiles.h"
#include "offLattice/offFileIO.h"
#include "offLattice/stlFileIO.h"
#include "offLattice/triangleMesh.h"
#include "offLattice/triangleSelector.h"

namespace plb {

template <typename T>
class RawTriangleMesh;
template <typename T>
class STLreader;

template <typename T>
RawTriangleMesh<T> stlToRawTriangleMesh(STLreader<T> const &reader);

template <typename T>
RawTriangleMesh<T> stlToRawTriangleMesh(FileName stlFileName, T eps = getEpsilon<T>(DBL));

template <typename T>
RawTriangleMesh<T> offToRawTriangleMesh(OFFreader<T> const &reader);

template <typename T>
RawTriangleMesh<T> offToRawTriangleMesh(FileName offFileName);

/// The raw triangle mesh can be subdivided into parts (in its internal
/// representation, the triangles are subdivided into sections that
/// reflect the geometry parts). Apart from that, there's no possibility
/// to assign properties and tags, though. For more advanced needs, you
/// must use the ConnectedTriangleMesh implementations.
template <typename T>
class RawTriangleMesh : public TriangleMesh<T> {
public:
    typedef Array<Array<T, 3>, 3> RawTriangle;
    typedef typename TriangleMesh<T>::PVertex PVertex;
    typedef typename TriangleMesh<T>::CPVertex CPVertex;
    typedef typename TriangleMesh<T>::PTriangle PTriangle;
    typedef typename TriangleMesh<T>::CPTriangle CPTriangle;
    typedef typename TriangleMesh<T>::PTriangleIterator PTriangleIterator;
    typedef typename TriangleMesh<T>::PVertexIterator PVertexIterator;
    class Triangle;

    class Vertex : public TriangleMesh<T>::Vertex {
    public:
        Vertex(RawTriangleMesh<T> *mesh_, plint part_, plint triangle_, plint vertex_);
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
        virtual bool isInterior()
            const;  /// Check if the vertex is an interior, non-boundary vertex.
        virtual std::vector<plint> adjacentVertices()
            const;  /// Return a list of adjacent vertices.
        virtual T property(
            plint whichProperty) const;  /// Get a vertex property (to know its name, ask the mesh).
        virtual void setProperty(plint whichProperty, T value);
        virtual plint tag(plint whichTag)
            const;  /// Get a vertex tag (to know its name, ask the mesh).
                    /// Tags with whichTag==0 are reserved for the unique vertex numbering.
        virtual void setTag(plint whichTag, plint value);

    private:
        RawTriangleMesh<T> *mesh;
        plint part, triangle, vertex;
    };

    class Triangle : public TriangleMesh<T>::Triangle {
    public:
        Triangle(RawTriangleMesh<T> *mesh_, plint part_, plint triangle_);
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
        /// Returns the number of triangles that meet on the vertex with the specified
        /// local index.
        virtual plint numVertexNeighbors(plint iVertex) const;
        /// Returns one of the triangles that meet on the vertex with local index "iVertex".
        /// "iNeighbor" must be greater equal to 0, and less than numVertexNeighbors(iVertex).
        /// If not, a NULL pointer is returned.
        virtual CPTriangle vertexNeighbor(plint iVertex, plint iNeighbor) const;
        virtual T property(plint id) const;  /// Any property (color, etc.).
        virtual void setProperty(plint whichProperty, T value);
        virtual plint tag(plint id) const;  /// Any tag (to get the tag for the part,
                                            /// one must use getTriangleTag("Part")).
        virtual void setTag(plint whichTag, plint value);

    private:
        RawTriangleMesh<T> *mesh;
        plint part, triangle;
    };
    /// Iterate over all triangles, forward only.
    class AllTriangleIterator : public TriangleMesh<T>::TriangleIterator {
    public:
        AllTriangleIterator(RawTriangleMesh<T> *mesh_);
        virtual PTriangle next();
        virtual bool end() const;
        virtual PTriangleIterator clone() const;

    private:
        RawTriangleMesh<T> *mesh;
        plint currentPart;
        plint currentTriangle;
    };
    /// Iterate over all triangles, forward only.
    class PartTriangleIterator : public TriangleMesh<T>::TriangleIterator {
    public:
        PartTriangleIterator(RawTriangleMesh<T> *mesh, plint partId_);
        virtual PTriangle next();
        virtual bool end() const;
        virtual PTriangleIterator clone() const;

    private:
        RawTriangleMesh<T> *mesh;
        plint partId;
        plint currentTriangle;
    };

public:
    /// Constructor for single part, with raw triangles.
    RawTriangleMesh(std::vector<RawTriangle> const &triangles_, T eps_ = getEpsilon<T>(DBL));
    /// Constructor for multiple parts, with raw triangles.
    RawTriangleMesh(
        std::vector<std::vector<RawTriangle> > const &triangles_,
        std::vector<std::string> const &partNames_, T eps_ = getEpsilon<T>(DBL));
    /// Copy constructor.
    RawTriangleMesh(RawTriangleMesh<T> const &rhs);
    /// Construct from any other type of triangle mesh.
    RawTriangleMesh(TriangleMesh<T> &rhs, T eps_ = getEpsilon<T>(DBL));
    RawTriangleMesh<T> merge(RawTriangleMesh<T> const &mesh) const;
    RawTriangleMesh<T> refine() const;
    RawTriangleMesh<T> refine(T edgeLengthThreshold) const;
    RawTriangleMesh<T> refineRecursively(
        T targetMaxEdgeLength, plint maxNumIterations, bool &success) const;
    RawTriangleMesh<T> refineByArea(T triangleAreaThreshold) const;
    RawTriangleMesh<T> refineByAreaRecursively(
        T targetMaxTriangleArea, plint maxNumIterations, bool &success) const;
    RawTriangleMesh<T> select(TriangleSelector<T> const &selector) const;
    RawTriangleMesh<T> cutWithPlane(Array<T, 3> planePos, Array<T, 3> orientation) const;
    T getEps() const
    {
        return eps;
    }

    virtual PTriangleIterator triangleIterator(plint partId = -1);
    virtual PVertexIterator vertexIterator();
    virtual plint numParts() const;
    virtual std::string partName(plint iPart) const;
    virtual plint partId(std::string partName) const;

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
    std::vector<std::vector<RawTriangle> > const &getTriangles() const;

    virtual plint getTriangleTag(std::string tagName) const;
    virtual std::string getTriangleTagName(plint tag) const;
    virtual plint numTriangleTags() const;

    virtual plint getNumTriangles() const;

private:
    bool triangleHasZeroArea(RawTriangle const &triangle) const;
    bool triangleHasZeroLengthEdges(RawTriangle const &triangle) const;
    bool areTheSameTriangle(RawTriangle const &t1, RawTriangle const &t2) const;
    bool containedAndErase(RawTriangle const &triangle, std::vector<RawTriangle> &triangles) const;
    bool cutTriangleWithPlane(
        RawTriangle const &triangle, Array<T, 3> const &planePos, Array<T, 3> const &normal,
        std::vector<RawTriangle> &newTriangles,
        std::vector<RawTriangle> &newZeroAreaTriangles) const;
    void computeMinMaxEdges();
    void computeMinMaxAreas();
    void computeBoundingCuboid();

private:
    std::vector<std::vector<RawTriangle> > triangles;
    std::vector<std::string> partNames;
    T minEdgeLength, maxEdgeLength;
    T minTriangleArea, maxTriangleArea;
    Cuboid<T> boundingCuboid;
    T eps;
};

}  // namespace plb

#endif  // RAW_TRIANGLE_MESH_H
