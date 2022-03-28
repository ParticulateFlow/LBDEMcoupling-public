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

#ifndef TRIANGLE_SET_H
#define TRIANGLE_SET_H

#include <cstdio>
#include <string>
#include <vector>

#include "core/array.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "offLattice/triangleSelector.h"

namespace plb {

template <typename T>
class TriangleSet {
public:
    typedef Array<Array<T, 3>, 3> Triangle;

public:
    TriangleSet(Precision precision_ = DBL);
    TriangleSet(T eps_);
    TriangleSet(
        std::vector<Triangle> const &triangles_, Precision precision_ = DBL,
        TriangleSelector<T> *selector = 0);
    TriangleSet(std::vector<Triangle> const &triangles_, T eps_, TriangleSelector<T> *selector = 0);
    // Currently only STL and OFF files are supported by this class.
    TriangleSet(
        std::string fname, Precision precision_ = DBL, SurfaceGeometryFileFormat fformat = STL,
        TriangleSelector<T> *selector = 0);
    TriangleSet(
        std::string fname, T eps_, SurfaceGeometryFileFormat fformat = STL,
        TriangleSelector<T> *selector = 0);

    TriangleSet<T> *clone() const
    {
        return new TriangleSet<T>(*this);
    }

    std::vector<Triangle> const &getTriangles() const;
    T getEpsilon() const
    {
        return eps;
    }
    void setPrecision(Precision precision_);
    void setEpsilon(T eps_);

    /// Translate the triangle set surface mesh.
    void translate(Array<T, 3> const &vector);
    /// Scale the triangle set surface mesh uniformly in all directions.
    void scale(T alpha);
    /// Scale the triangle set surface mesh differently in each direction.
    void scale(T alpha, T beta, T gamma);
    /// Rotate the triangle set surface mesh.
    ///   The arguments of this function are the Euler angles in radians.
    ///   The so-called "x-convention" is used, in which the rotation is
    ///   given by the three angles (phi, theta, psi), where:
    ///   1.] The first rotation is by an angle phi about the z-axis,
    ///   2.] The second rotation is by an angle theta in [0, pi] about
    ///       the new x-axis,
    ///   3.] The third rotation is by an angle psi about the new z-axis.
    void rotate(T phi, T theta, T psi);
    /// Rotate the triangle set by an angle theta around an axis that
    ///   traverses the origin and points in the direction normedAxis
    ///   (which must be normalized).
    void rotateAtOrigin(Array<T, 3> const &normedAxis, T theta);
    /// Executes a rotation to the obstacle which, if applied to the
    ///   x- and y- axes, would apply them to the new x- and y- axes
    ///   as provided.
    void rotateAxesTo(Array<T, 3> const &ex, Array<T, 3> const &ey);
    /// Executes a rotation to the obstacle which, if applied to the
    ///   provided new x- and y- axes, would rotate them back to the
    ///   x- and y-axes of the origin.
    void rotateAxesFrom(Array<T, 3> const &ex, Array<T, 3> const &ey);

    /// Clear the triangle mesh (free all the associated memory as well).
    void clear();

    /// Erase the current triangle set surface mesh, and merge into it the new meshes.
    ///   This function currently does not check for duplicate
    ///   triangles in the new meshes, and does not handle
    ///   cases when the resulting merged surface mesh is
    ///   non-manifold. It is in the user's jurisdiction to
    ///   make sure that the resulting mesh is well defined.
    ///   Practically, this can be achieved if the new triangle set
    ///   meshes given as arguments are mutually disjoint.
    void merge(std::vector<TriangleSet<T> *> meshes);

    /// Append to the current triangle mesh, the mesh that is passed as an argument.
    ///   This function currently does not check for duplicate
    ///   triangles in the two meshes, and does not handle
    ///   cases when the resulting merged surface mesh is
    ///   non-manifold. It is in the user's jurisdiction to
    ///   make sure that the resulting mesh is well defined.
    ///   Practically, this can be achieved if the two triangle set
    ///   meshes are mutually disjoint.
    void append(TriangleSet<T> const &mesh);

    /// Refine the current triangle set surface mesh by splitting each
    ///   original triangle into four smaller triangles constructed by
    ///   joining the middle points of the original edges. The old mesh
    ///   is deleted.
    void refine();
    /// Refine the current triangle set surface mesh by considering
    ///   edge lengths. If the length of an edge is greater than or
    ///   equal to the edgeLengthThreshold, then this edge is split
    ///   in two. The old mesh is deleted.
    void refine(T edgeLengthThreshold);
    /// Refine the current triangle set surface mesh by considering
    ///   edge lengths. If the length of an edge is greater than or
    ///   equal to the targetMaxEdgeLength, then this edge is split
    ///   in two. The procedure is repeated until the maximum edge length
    ///   is less than the one provided, or until the maxNumIterations
    ///   have been reached. The old mesh is deleted. The function returns
    ///   true upon success or false otherwise.
    bool refineRecursively(T targetMaxEdgeLength, plint maxNumIterations);
    /// Refine the current triangle set surface mesh by considering
    ///   triangle areas. If the area of a triangle is greater than or
    ///   equal to the triangleAreaThreshold, then this triangle is
    ///   split into four smaller triangles constructed by
    ///   joining the middle points of the original edges.
    ///   The old mesh is deleted.
    ///
    ///   CAUTION: The resulting triangle mesh usually is non-conforming.
    void refineByArea(T triangleAreaThreshold);
    /// Refine the current triangle set surface mesh by considering
    ///   triangle areas. If the area of a triangle is greater than or
    ///   equal to the targetMaxTriangleArea, then this triangle is
    ///   split into four smaller triangles constructed by
    ///   joining the middle points of the original edges. The procedure
    ///   is repeated until the maximum triangle area is less than the
    ///   one provided, or until the maxNumIterations have been reached.
    ///   The old mesh is deleted. The function returns true upon success
    ///   or false otherwise.
    ///
    ///   CAUTION: The resulting triangle mesh usually is non-conforming.
    bool refineByAreaRecursively(T targetMaxTriangleArea, plint maxNumIterations);

    /// Select a subset of all triangles by using a TriangleSelector.
    ///   The TriangleSet contains only the selected triangles after this operation.
    void select(TriangleSelector<T> const &selector);

    /// Remove all triangles with a normal n that verifies n*normal >= tolerance.
    void removeTrianglesWithOrientation(Array<T, 3> const &normal, T tolerance);

    /// A very simple orientation reversing function.
    void reverseOrientation();

    /// Cast all coordinates to float and then back to type T.
    void toFloatAndBack();

    /// Check if the triangle set has any triangles with zero area.
    bool hasZeroAreaTriangles() const;

    /// Find the number of triangles with zero area (if any).
    plint numZeroAreaTriangles() const;

    /// Check if the the triangle set mesh has a dependence on the floating point
    /// precision.
    bool hasFloatingPointPrecisionDependence() const;

    /// Export the mesh as an ASCII STL file.
    /// If mainProcOnly = false, then all processes will write, a fact which the
    /// caller must take under consideration when constructing the fname.
    void writeAsciiSTL(
        std::string fname, int numDecimalDigits = 10, bool mainProcOnly = true,
        TriangleSelector<T> *selector = 0) const;
    void writeAsciiSTL(
        std::string fname, T scale, Array<T, 3> const &offset, int numDecimalDigits = 10,
        bool mainProcOnly = true, TriangleSelector<T> *selector = 0) const;
    void writeAsciiContentSTL(
        FILE *fp, std::string solidName, T scale, Array<T, 3> const &offset, int numDecimalDigits,
        TriangleSelector<T> *selector) const;

    /// Export the mesh as an binary STL file.
    /// If mainProcOnly = false, then all processes will write, a fact which the
    /// caller must take under consideration when constructing the fname.
    void writeBinarySTL(
        std::string fname, bool mainProcOnly = true, TriangleSelector<T> *selector = 0) const;
    void writeBinarySTL(
        std::string fname, T scale, Array<T, 3> const &offset, bool mainProcOnly = true,
        TriangleSelector<T> *selector = 0) const;
    void writeBinaryContentSTL(
        FILE *fp, std::string solidName, T scale, Array<T, 3> const &offset,
        TriangleSelector<T> *selector) const;

    /// Cut the current triangle set mesh by a plane "plane" which is
    ///   defined by a point and a normal. This cutting operation will
    ///   produce a new triangle set "newTriangleSet" which is
    ///   the one that belongs to the half-space that the normal of the
    ///   cutting plane points outside of. The function returns 1 if
    ///   the cut was successful, 0 if there was no intersection between
    ///   the original triangle set and the plane provided, and -1 if an
    ///   error occurred.
    int cutWithPlane(Plane<T> const &plane, TriangleSet<T> &newTriangleSet) const;

    /// Cut the current triangle set mesh by a plane "plane" which is
    ///   defined by a point and a normal and a cuboid "cuboid" which is
    ///   defined by a lower left corner and a right upper corner. This
    ///   cutting operation will cut off all triangles, or parts of triangles
    ///   that are fully contained in the cuboid and are positioned in
    ///   the half-space that the normal of the cutting plane points into.
    ///   Caution is needed, because this function does not use the cuboid
    ///   to cut, but only to select parts of the original triangle set,
    ///   to be then cut by the plane. Obviously, at least a part of the
    ///   cutting plane must be contained in the cuboid, for the cutting
    ///   to make any sense. If not used wisely, this function can lead to
    ///   broken STL files. The function returns 1 if the cut was successful,
    ///   0 if there was no intersection between the original triangle set
    ///   the plane and the cuboid provided, and -1 if an error occurred.
    int cutWithPlane(
        Plane<T> const &plane, Cuboid<T> const &cuboid, TriangleSet<T> &newTriangleSet) const;

    T getMinEdgeLength() const
    {
        return minEdgeLength;
    }
    T getMaxEdgeLength() const
    {
        return maxEdgeLength;
    }

    T getMinTriangleArea() const
    {
        return minTriangleArea;
    }
    T getMaxTriangleArea() const
    {
        return maxTriangleArea;
    }

    Cuboid<T> getBoundingCuboid() const
    {
        return boundingCuboid;
    }
    Array<T, 3> getCentroid() const;
    Array<T, 3> getCenterOfMass() const;

private:
    void readSTL(std::string fname, TriangleSelector<T> *selector);
    bool isAsciiSTL(FILE *fp);
    void readAsciiSTL(FILE *fp, TriangleSelector<T> *selector);
    void readBinarySTL(FILE *fp, TriangleSelector<T> *selector);
    void readOFF(std::string fname, TriangleSelector<T> *selector);
    void readAsciiOFF(FILE *fp, TriangleSelector<T> *selector);
    bool triangleHasZeroArea(Triangle const &triangle, T epsilon) const;
    bool triangleHasZeroLengthEdges(Triangle const &triangle, T epsilon) const;
    void fixOrientation(Triangle &triangle, Array<T, 3> const &n) const;
    bool areTheSameTriangle(Triangle const &t1, Triangle const &t2, T epsilon) const;
    bool containedAndErase(
        Triangle const &triangle, std::vector<Triangle> &triangles, T epsilon) const;
    void computeMinMaxEdge(pluint iTriangle, T &minEdge, T &maxEdge) const;
    void computeMinMaxEdges();
    T computeArea(pluint iTriangle) const;
    void computeMinMaxAreas();
    void computeBoundingCuboid();
    /// Cut the current triangle by a plane "plane" which is
    ///   defined by a point and a normal. This cutting operation will
    ///   add or not one or more triangles to the vectors of triangles.
    ///   The triangles added belong to the half-space that the normal of
    ///   the plane given points outside of. The function returns 1 if the
    ///   cut was successful, 0 if there was no intersection between the
    ///   original triangle set and the plane and -1 if an error occurred.
    int cutTriangleWithPlane(
        Plane<T> const &plane, Triangle const &triangle, std::vector<Triangle> &newTriangles,
        std::vector<Triangle> &newZeroAreaTriangles) const;
    /// Skip nLines number of lines in the file.
    void skipLines(plint nLines, FILE *fp) const;
    /// Read the file character-by-character. This function returns 0 if a non-white-space
    /// character is found, or EOF if the end-of-file is found. It "consumes" the
    /// rest of the line if the "commentCharacter" is found. Obviously, it works only
    /// with text files.
    int readAhead(FILE *fp, char commentCharacter) const;
    /// Check if the buffer which holds an input line of text is full.
    bool checkForBufferOverflow(char *buf) const;

private:
    std::vector<Triangle> triangles;
    T minEdgeLength, maxEdgeLength;
    T minTriangleArea, maxTriangleArea;
    Cuboid<T> boundingCuboid;
    T eps;  // Tolerance for floating point comparisons.
};

}  // namespace plb

#endif  // TRIANGLE_SET_H
