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

#ifndef BOUNDARY_SHAPES_3D_HH
#define BOUNDARY_SHAPES_3D_HH

#include <limits>

#include "offLattice/boundaryShapes3D.h"

namespace plb {

template <typename T, class SurfaceData>
BoundaryShapeIntersection3D<T, SurfaceData>::BoundaryShapeIntersection3D()
{ }

template <typename T, class SurfaceData>
BoundaryShapeIntersection3D<T, SurfaceData>::BoundaryShapeIntersection3D(
    BoundaryShapeIntersection3D<T, SurfaceData> const &rhs)
{
    components.resize(rhs.components.size());
    for (pluint iComponent = 0; iComponent < components.size(); ++iComponent) {
        components[iComponent] = rhs.components[iComponent]->clone();
    }
}

template <typename T, class SurfaceData>
BoundaryShapeIntersection3D<T, SurfaceData> &BoundaryShapeIntersection3D<T, SurfaceData>::operator=(
    BoundaryShapeIntersection3D<T, SurfaceData> const &rhs)
{
    BoundaryShapeIntersection3D<T, SurfaceData>(rhs).swap(*this);
    return *this;
}

template <typename T, class SurfaceData>
BoundaryShapeIntersection3D<T, SurfaceData>::~BoundaryShapeIntersection3D()
{
    for (pluint iComponent = 0; iComponent < components.size(); ++iComponent) {
        delete components[iComponent];
    }
}

template <typename T, class SurfaceData>
void BoundaryShapeIntersection3D<T, SurfaceData>::swap(
    BoundaryShapeIntersection3D<T, SurfaceData> &rhs)
{
    components.swap(rhs.components);
}

template <typename T, class SurfaceData>
void BoundaryShapeIntersection3D<T, SurfaceData>::addComponent(
    BoundaryShape3D<T, SurfaceData> *component)
{
    components.push_back(component);
}

template <typename T, class SurfaceData>
bool BoundaryShapeIntersection3D<T, SurfaceData>::isInside(Array<T, 3> const &location) const
{
    bool isInside = true;
    for (pluint iComponent = 0; iComponent < components.size(); ++iComponent) {
        if (!components[iComponent]->isInside(location)) {
            isInside = false;
        }
    }
    return isInside;
}

template <typename T, class SurfaceData>
bool BoundaryShapeIntersection3D<T, SurfaceData>::isInside(
    Array<T, 3> const &location, T epsilon) const
{
    bool isInside = true;
    for (pluint iComponent = 0; iComponent < components.size(); ++iComponent) {
        if (!components[iComponent]->isInside(location, epsilon)) {
            isInside = false;
        }
    }
    return isInside;
}

template <typename T, class SurfaceData>
bool BoundaryShapeIntersection3D<T, SurfaceData>::pointOnSurface(
    Array<T, 3> const &fromPoint, Array<T, 3> const &direction, Array<T, 3> &locatedPoint,
    T &distance, Array<T, 3> &wallNormal, SurfaceData &surfaceData, OffBoundary::Type &bdType,
    plint &id) const
{
    bool doesIntersect = false;
    for (pluint iComponent = 0; iComponent < components.size(); ++iComponent) {
        Array<T, 3> positionOnWall;
        T newDistance;
        Array<T, 3> newWallNormal;
        SurfaceData newSurfaceData;
        OffBoundary::Type newBdType;
        if (components[iComponent]->pointOnSurface(
                fromPoint, direction, positionOnWall, newDistance, newWallNormal, newSurfaceData,
                newBdType, id))
        {
            if (!doesIntersect || (doesIntersect && newDistance < distance)) {
                locatedPoint = positionOnWall;
                distance = newDistance;
                wallNormal = newWallNormal;
                surfaceData = newSurfaceData;
                bdType = newBdType;
            }
            doesIntersect = true;
        }
    }
    return doesIntersect;
}

template <typename T, class SurfaceData>
bool BoundaryShapeIntersection3D<T, SurfaceData>::intersectsSurface(
    Array<T, 3> const &p1, Array<T, 3> const &p2, plint &id) const
{
    bool doesIntersect = false;
    for (pluint iComponent = 0; iComponent < components.size(); ++iComponent) {
        Array<T, 3> positionOnWall;
        T newDistance;
        Array<T, 3> newWallNormal;
        SurfaceData newSurfaceData;
        OffBoundary::Type newBdType;
        if (components[iComponent]->pointOnSurface(
                p1, p2 - p1, positionOnWall, newDistance, newWallNormal, newSurfaceData, newBdType,
                id))
        {
            doesIntersect = true;
        }
    }
    return doesIntersect;
}

template <typename T, class SurfaceData>
bool BoundaryShapeIntersection3D<T, SurfaceData>::distanceToSurface(
    Array<T, 3> const &point, T &distance, bool &isBehind) const
{
    bool surfaceFound = false;
    for (pluint iComponent = 0; iComponent < components.size(); ++iComponent) {
        T newDistance;
        bool newIsBehind;
        if (components[iComponent]->distanceToSurface(point, newDistance, newIsBehind)) {
            if (!surfaceFound || (surfaceFound && newDistance < distance)) {
                distance = newDistance;
                isBehind = newIsBehind;
            }
            surfaceFound = true;
        }
    }
    return surfaceFound;
}

template <typename T, class SurfaceData>
plint BoundaryShapeIntersection3D<T, SurfaceData>::getTag(plint id) const
{
    for (pluint iComponent = 0; iComponent < components.size(); ++iComponent) {
        plint tag = components[iComponent]->getTag(id);
        if (tag >= 0) {
            return tag;
        }
    }
    return -1;
}

template <typename T, class SurfaceData>
BoundaryShapeIntersection3D<T, SurfaceData> *BoundaryShapeIntersection3D<T, SurfaceData>::clone()
    const
{
    return new BoundaryShapeIntersection3D<T, SurfaceData>(*this);
}

}  // namespace plb

#endif  // BOUNDARY_SHAPES_3D_H
