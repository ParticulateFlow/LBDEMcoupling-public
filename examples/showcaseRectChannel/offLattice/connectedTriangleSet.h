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

#ifndef CONNECTED_TRIANGLE_SET_H
#define CONNECTED_TRIANGLE_SET_H

#include <set>
#include <string>
#include <vector>

#include "core/array.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "offLattice/triangleSet.h"

namespace plb {

/* ConnectedTriangleSet: This class is similar to the TriangleSet class with
 * the addition that also contains connectivity information. In this data
 * structure, a vertex is a point in space, a triangle is a set of 3 integers
 * that correspond to the global indices of vertices. There is also a list
 * of triangle indices that meet on a vertex. The construction of the data
 * of this class performs very few checks so that objects of this class can
 * be instantiated only with an object of type TriangleSet. What this means
 * is that the ConnectedTriangleSet represents a triangular mesh which does
 * not need to be well formed or even manifold.
 */
template <typename T>
class ConnectedTriangleSet {
public:
    ConnectedTriangleSet(TriangleSet<T> const &triangleSet);

    plint getNumVertices() const
    {
        return numVertices;
    }
    plint getNumTriangles() const
    {
        return numTriangles;
    }
    Array<T, 3> getVertex(plint iVertex) const
    {
        return vertices[iVertex];
    }
    Array<plint, 3> getTriangle(plint iTriangle) const
    {
        return triangles[iTriangle];
    }
    std::vector<plint> getTrianglesOnVertex(plint iVertex) const
    {
        return trianglesOnVertex[iVertex];
    }
    /*
     * This class contains geometrical and topological information in its data.
     * The only piece of geometrical information is the vector "vertices" which holds
     * the vertex positions. All the other vectors contain topological information.
     * It is possible in some occasions, that the geometry of a surface changes (from
     * a rigid body movement for example) but the topology remains the same. The "swapGeometry"
     * function gives the possibility to the user to change the geometry of the surface,
     * but retain its topology. It is a responsibility of the user to make sure that the
     * vector "newVertices" provided as an argument contains the vertices of the surface
     * in the same order as the "vertices" data member of the object (in other words: the surface
     * topology must be retained). This function will swap the vertex vectors, as its name
     * suggests.
     */
    void swapGeometry(std::vector<Array<T, 3> > &newVertices);
    /*
     * Relative to the comment above, quantities like the area and the unit normal can
     * be computed for any set of vertices that have the same topology. The next functions
     * (along with the output OFF function) take as an argument a pointer to a vector of
     * vertex positions. If this pointer is 0, then the default set of vertices (the one
     * contained in the object of the class) is used. If it is not 0, then this user-provided
     * vertex set is used for the computations. The optional vertex position vector can be
     * a larger vector than the one corresponding to the current ConnectedTriangleSet object.
     * The optional parameter "indexOffset" (meaningful only when "newVertices" is not 0),
     * is the index in the new vector which maps to the vertex with index 0 in the current
     * ConnectedTriangleSet object. To make this clearer, sometimes the user at his application
     * has the vertices of many surfaces put in a single vector. In such a case,
     * he has a list of index offsets to know where each surface starts in this large
     * vector of vertex positions. This "starting index" is the "indexOffset" in these utility
     * functions. If at the current ConnectedTriangleSet object, a vertex has a local index "v0",
     * then its corresponding index in the "newVertices" vector will be "v0 + indexOffset".
     * Needless to say that "iVertex" and "iTriangle" are the local ids of the current
     * ConnectedTriangleSet object.
     */
    void computeVertexAreaAndUnitNormal(
        plint iVertex, T &area, Array<T, 3> &unitNormal, std::vector<Array<T, 3> > *newVertices = 0,
        plint indexOffset = 0) const;
    void computeTriangleAreaAndUnitNormal(
        plint iTriangle, T &area, Array<T, 3> &unitNormal,
        std::vector<Array<T, 3> > *newVertices = 0, plint indexOffset = 0) const;
    void writeOFF(
        std::string fname, std::vector<Array<T, 3> > *newVertices = 0, plint indexOffset = 0,
        int numDecimalDigits = 10) const;
    TriangleSet<T> *toTriangleSet(
        Precision precision, std::vector<Array<T, 3> > const *newVertices = 0,
        plint indexOffset = 0) const;
    TriangleSet<T> *toTriangleSet(
        T eps, std::vector<Array<T, 3> > const *newVertices = 0, plint indexOffset = 0) const;

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
    plint numVertices, numTriangles;                     // Total number of unique vertices,
                                                         // and total number of triangles.
    std::vector<Array<T, 3> > vertices;                  // Positions of vertices. The vector
                                                         // index is the global vertex index.
    std::vector<Array<plint, 3> > triangles;             // Global indices of the vertices that
                                                         // constitute the triangle. The vector
                                                         // index is the global triangle index.
    std::vector<std::vector<plint> > trianglesOnVertex;  // List of the triangle global indices
                                                         // that meet on a vertex.
};

}  // namespace plb

#endif  // CONNECTED_TRIANGLE_SET_H
