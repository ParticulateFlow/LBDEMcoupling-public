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

#ifndef MESH_ALGORITHM_HH
#define MESH_ALGORITHM_HH

/*
 * Algorithms, applied to generic meshes.
 */

#include <algorithm>
#include <limits>
#include <set>

#include "core/array.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "offLattice/meshAlgorithm.h"

namespace plb {

template <typename T>
void toLatticeUnits(TriangleMesh<T> &mesh, T dx, Array<T, 3> const &offset)
{
    if (!util::isZero(norm(offset))) {
        mesh.translate(-offset);
    }
    if (!util::isOne(dx)) {
        mesh.scale(1.0 / dx);
    }
}

template <typename T>
void toLatticeUnits(TriangleMesh<T> &mesh, Units3D const &units)
{
    toLatticeUnits(mesh, units.dx(), units.physOffset());
}

template <typename T>
void toLatticeUnits(TriangleSet<T> &triangleSet, Units3D const &units)
{
    if (!util::isZero(norm(units.physOffset()))) {
        triangleSet.translate(-units.physOffset());
    }
    if (!util::isOne(units.dx())) {
        triangleSet.scale(1.0 / units.dx());
    }
}

template <typename T>
void toLatticeUnits(
    TriangleMesh<T> &mesh, Array<T, 3> const &center, T dx, Array<T, 3> const &offset)
{
    Cuboid<T> bCuboid = mesh.getBoundingCuboid();
    Array<T, 3> meshCenter = (T)0.5 * (bCuboid.lowerLeftCorner + bCuboid.upperRightCorner);
    if (!util::isZero(norm(center - meshCenter))) {
        mesh.translate(center - meshCenter);
    }
    toLatticeUnits(mesh, dx, offset);
}

template <typename T>
void toPhysicalUnits(TriangleMesh<T> &mesh, Units3D const &units)
{
    toPhysicalUnits(mesh, units.dx(), units.physOffset());
}

template <typename T>
void toPhysicalUnits(TriangleSet<T> &triangleSet, Units3D const &units)
{
    triangleSet.scale(units.dx());
    triangleSet.translate(units.physOffset());
}

template <typename T>
void toPhysicalUnits(TriangleMesh<T> &mesh, T dx, Array<T, 3> offset)
{
    mesh.scale(dx);
    mesh.translate(offset);
}

template <typename T>
void smooth(
    RawConnectedTriangleMesh<T> &triangleMesh, plint maxiter, T relax, bool isMeasureWeighted)
{
    typedef typename RawConnectedTriangleMesh<T>::PVertexIterator PVertexIterator;
    typedef typename RawConnectedTriangleMesh<T>::PVertex PVertex;
    plint uniqueIDTag = triangleMesh.getVertexTag("UniqueID");

    PLB_ASSERT(maxiter >= 0);
    PLB_ASSERT(util::greaterEqual(relax, (T)0) && util::lessEqual(relax, (T)1));

    if (maxiter <= 0)
        return;

    RawConnectedTriangleMesh<T> *bp0;
    bp0 = &triangleMesh;

    RawConnectedTriangleMesh<T> buf(triangleMesh);

    RawConnectedTriangleMesh<T> *bp1;
    bp1 = &buf;

    if (isMeasureWeighted) {
        for (plint iter = 0; iter < maxiter; iter++) {
            PVertexIterator vertexIterator(triangleMesh.vertexIterator());
            while (!vertexIterator->end()) {
                PVertex vertex(vertexIterator->next());
                plint iVertex = vertex->tag(uniqueIDTag);
                Array<T, 3> iVertexPos = (*bp0).vertex(iVertex)->get();
                std::vector<plint> neighborVertexIds = vertex->adjacentVertices();
                if (vertex->isInterior()) {
                    pluint sizeIds = neighborVertexIds.size();
                    T area = (T)0.0;
                    Array<T, 3> tmp((T)0.0, (T)0.0, (T)0.0);
                    for (pluint i = 0; i < sizeIds; i++) {
                        plint jVertex = neighborVertexIds[i];
                        plint kVertex =
                            (i + 1 < sizeIds) ? neighborVertexIds[i + 1] : neighborVertexIds[0];

                        Array<T, 3> jVertexPos = (*bp0).vertex(jVertex)->get();
                        Array<T, 3> kVertexPos = (*bp0).vertex(kVertex)->get();

                        Array<T, 3> middleVertexPos =
                            (iVertexPos + jVertexPos + kVertexPos) / (T)3.0;

                        Array<T, 3> eij = jVertexPos - iVertexPos;
                        Array<T, 3> eik = kVertexPos - iVertexPos;
                        Array<T, 3> n;
                        crossProduct(eij, eik, n);
                        T locArea = (T)0.5 * norm(n);

                        area += locArea;
                        tmp += middleVertexPos * locArea;
                    }
                    (*bp1).vertex(iVertex)->get() =
                        ((T)1.0 - relax) * (*bp0).vertex(iVertex)->get() + relax * tmp / area;
                } else {  // vertex->isInterior()
                    plint jVertex = -1, kVertex = -1;
                    int counter = 0;
                    for (pluint i = 0; i < neighborVertexIds.size(); i++) {
                        plint vertexId = neighborVertexIds[i];
                        PVertex neighborVertex(triangleMesh.vertex(vertexId));
                        if (!neighborVertex->isInterior()) {
                            counter++;
                            if (counter == 1) {
                                jVertex = vertexId;
                            } else if (counter == 2) {
                                kVertex = vertexId;
                            } else {
                                PLB_ASSERT(
                                    false);  // Problem with the boundary of the surface mesh.
                            }
                        }
                    }

                    PLB_ASSERT(counter == 2);  // Problem with the boundary of the surface mesh.

                    Array<T, 3> jVertexPos = (*bp0).vertex(jVertex)->get();
                    Array<T, 3> kVertexPos = (*bp0).vertex(kVertex)->get();

                    Array<T, 3> middleVertexPos0 = (iVertexPos + jVertexPos) / (T)2.0;
                    Array<T, 3> middleVertexPos1 = (iVertexPos + kVertexPos) / (T)2.0;

                    T length0 = norm(jVertexPos - iVertexPos);
                    T length1 = norm(kVertexPos - iVertexPos);

                    (*bp1).vertex(iVertex)->get() =
                        ((T)1.0 - relax) * (*bp0).vertex(iVertex)->get()
                        + relax * (middleVertexPos0 * length0 + middleVertexPos1 * length1)
                              / (length0 + length1);
                }
            }
            // Swap buffers.
            std::swap(bp0, bp1);
        }
    } else {
        for (plint iter = 0; iter < maxiter; iter++) {
            PVertexIterator vertexIterator(triangleMesh.vertexIterator());
            while (!vertexIterator->end()) {
                PVertex vertex(vertexIterator->next());
                plint iVertex = vertex->tag(uniqueIDTag);
                std::vector<plint> neighborVertexIds = vertex->adjacentVertices();

                if (vertex->isInterior()) {
                    pluint sizeIds = neighborVertexIds.size();
                    Array<T, 3> tmp((T)0.0, (T)0.0, (T)0.0);
                    for (pluint i = 0; i < sizeIds; i++) {
                        tmp += (*bp0).vertex(neighborVertexIds[i])->get();
                    }
                    (*bp1).vertex(iVertex)->get() =
                        ((T)1.0 - relax) * (*bp0).vertex(iVertex)->get() + relax * tmp / (T)sizeIds;
                } else  // vertex->isInterior()
                {
                    plint jVertex = -1, kVertex = -1;
                    int counter = 0;
                    for (pluint i = 0; i < neighborVertexIds.size(); i++) {
                        plint vertexId = neighborVertexIds[i];
                        PVertex neighborVertex(triangleMesh.vertex(vertexId));
                        if (!neighborVertex->isInterior()) {
                            counter++;
                            if (counter == 1) {
                                jVertex = vertexId;
                            } else if (counter == 2) {
                                kVertex = vertexId;
                            } else {
                                PLB_ASSERT(
                                    false);  // Problem with the boundary of the surface mesh.
                            }
                        }
                    }

                    PLB_ASSERT(counter == 2);  // Problem with the boundary of the surface mesh.

                    (*bp1).vertex(iVertex)->get() =
                        ((T)1.0 - relax) * (*bp0).vertex(iVertex)->get()
                        + relax * ((*bp0).vertex(jVertex)->get() + (*bp0).vertex(kVertex)->get())
                              / (T)2.0;
                }
            }
            // Swap buffers.
            std::swap(bp0, bp1);
        }
    }

    // Final copy of vertex positions.
    if (&triangleMesh != bp0) {
        triangleMesh = *bp0;
    }
    triangleMesh.signalVertexUpdate();
}

template <typename T>
Array<T, 3> centerOfMass(TriangleMesh<T> &mesh)
{
    typedef typename RawConnectedTriangleMesh<T>::PTriangleIterator PTriangleIterator;
    typedef typename RawConnectedTriangleMesh<T>::PTriangle PTriangle;
    PTriangleIterator iterator(mesh.triangleIterator());
    T totalVolume = T();
    Array<T, 3> center(0., 0., 0.);
    while (!iterator->end()) {
        PTriangle triangle(iterator->next());
        T currentVolume = ((*triangle)[0][0] * (*triangle)[1][1] * (*triangle)[2][2]
                           - (*triangle)[0][0] * (*triangle)[2][1] * (*triangle)[1][2]
                           - (*triangle)[1][0] * (*triangle)[0][1] * (*triangle)[2][2]
                           + (*triangle)[1][0] * (*triangle)[2][1] * (*triangle)[0][2]
                           + (*triangle)[2][0] * (*triangle)[0][1] * (*triangle)[1][2]
                           - (*triangle)[2][0] * (*triangle)[1][1] * (*triangle)[0][2])
                          / (T)6;
        totalVolume += currentVolume;
        center[0] +=
            (((*triangle)[0][0] + (*triangle)[1][0] + (*triangle)[2][0]) / (T)4) * currentVolume;
        center[1] +=
            (((*triangle)[0][1] + (*triangle)[1][1] + (*triangle)[2][1]) / (T)4) * currentVolume;
        center[2] +=
            (((*triangle)[0][2] + (*triangle)[1][2] + (*triangle)[2][2]) / (T)4) * currentVolume;
    }
    return center / totalVolume;
}

template <typename T>
void rotateAxesTo(TriangleMesh<T> &mesh, Array<T, 3> const &ex, Array<T, 3> const &ey)
{
    typedef typename RawConnectedTriangleMesh<T>::PTriangleIterator PTriangleIterator;
    typedef typename RawConnectedTriangleMesh<T>::PTriangle PTriangle;
    Array<T, 3> ez(crossProduct(ex, ey));
    PTriangleIterator iterator(mesh.triangleIterator());
    while (!iterator->end()) {
        PTriangle triangle(iterator->next());
        for (int iVertex = 0; iVertex < 3; iVertex++) {
            Array<T, 3> &v = (*triangle)[iVertex];
            v = Array<T, 3>(
                ex[0] * v[0] + ey[0] * v[1] + ez[0] * v[2],
                ex[1] * v[0] + ey[1] * v[1] + ez[1] * v[2],
                ex[2] * v[0] + ey[2] * v[1] + ez[2] * v[2]);
        }
    }

    mesh.signalVertexUpdate();
}

template <typename T>
void rotateAxesFrom(TriangleMesh<T> &mesh, Array<T, 3> const &ex, Array<T, 3> const &ey)
{
    typedef typename RawConnectedTriangleMesh<T>::PTriangleIterator PTriangleIterator;
    typedef typename RawConnectedTriangleMesh<T>::PTriangle PTriangle;
    Array<T, 3> ez(crossProduct(ex, ey));
    PTriangleIterator iterator(mesh.triangleIterator());
    while (!iterator->end()) {
        PTriangle triangle(iterator->next());
        for (int iVertex = 0; iVertex < 3; iVertex++) {
            Array<T, 3> &v = (*triangle)[iVertex];
            v = Array<T, 3>(
                ex[0] * v[0] + ex[1] * v[1] + ex[2] * v[2],
                ey[0] * v[0] + ey[1] * v[1] + ey[2] * v[2],
                ez[0] * v[0] + ez[1] * v[1] + ez[2] * v[2]);
        }
    }

    mesh.signalVertexUpdate();
}

}  // namespace plb

#endif  // MESH_ALGORITHM_HH
