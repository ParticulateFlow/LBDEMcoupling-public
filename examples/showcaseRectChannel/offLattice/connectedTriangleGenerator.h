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

#ifndef CONNECTED_TRIANGLE_GENERATOR_H
#define CONNECTED_TRIANGLE_GENERATOR_H

#include "offLattice/connectedTriangleMesh.h"

namespace plb {

// Avoid conflict with legacy functions of same name.
namespace connected {

template <typename T>
RawConnectedTriangleMesh<T> constructSphere(
    Array<T, 3> const &center, T radius, plint minNumOfTriangles, T eps = getEpsilon<T>());

template <typename T>
RawConnectedTriangleMesh<T> constructCylinder(
    Array<T, 3> const &inletCenter, T inletRadius, T outletRadius, T length, plint nAxial,
    plint nCirc, T eps = getEpsilon<T>());

template <typename T>
RawConnectedTriangleMesh<T> constructCylinder(
    Array<T, 3> const &inletCenter, T inletRadius, T outletRadius, T length, plint nAxial,
    plint nCirc, std::vector<Array<T, 3> > &inletPoints, T eps = getEpsilon<T>());

template <typename T>
RawConnectedTriangleMesh<T> constructCylinder(
    Array<T, 3> const &inletCenter, Array<T, 3> const &axis, T inletRadius, T outletRadius,
    T length, plint nAxial, plint nCirc, std::vector<Array<T, 3> > &inletPoints,
    T eps = getEpsilon<T>());

template <typename T>
RawConnectedTriangleMesh<T> constructCuboid(
    Array<T, 3> const &lowerCorner, Array<T, 3> const &upperCorner,
    Array<plint, 3> const &nSegments, T eps = getEpsilon<T>());

template <typename T>
RawConnectedTriangleMesh<T> patchTubes(
    RawTriangleMesh<T> &geometryWithOpenings, plint sortDirection, std::vector<T> patchLengths);

template <typename T>
RawConnectedTriangleMesh<T> constructRectangle(
    T lx, T ly, plint nx, plint ny, T eps = getEpsilon<T>());

template <typename T>
RawConnectedTriangleMesh<T> constructStrip(
    std::vector<Array<T, 3> > from, std::vector<Array<T, 3> > to, T eps = getEpsilon<T>());

template <typename T>
RawConnectedTriangleMesh<T> constructAnnulus(
    Array<T, 3> const &center, T innerRadius, T outerRadius, plint nCirc, T eps = getEpsilon<T>());

}  // namespace connected

}  // namespace plb

#endif  // CONNECTED_TRIANGLE_GENERATOR_HH
