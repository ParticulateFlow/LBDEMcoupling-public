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

#ifndef BOUNDARY_SHAPES_3D_H
#define BOUNDARY_SHAPES_3D_H

#include "atomicBlock/atomicBlock3D.h"
#include "core/array.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"

namespace plb {

struct OffBoundary {
    enum Type { dirichlet, neumann, freeSlip, constRhoInlet, densityNeumann, flux, isolation };
};

/// Description of a shape as the boundary of a given volume.
template <typename T, class SurfaceData>
struct BoundaryShape3D {
    virtual ~BoundaryShape3D() { }
    /// Decide whether a given discrete lattice node is inside the solid shape.
    /** The reason why there is an isInside function instead of a isOutside
     *  one is that while there exists only one voxelFlag::inside flag,
     *  there are several flags which count as outside: undetermined,
     *  outside, and borderline. This is particularly important for the
     *  undetermined flag, because the outer envelopes in a multi-block
     *  structure have no special meaning and are default-initialized to
     *  undetermined.
     **/
    virtual bool isInside(Dot3D const &location) const = 0;
    /// Decide whether a given discrete lattice node is outside the solid shape.
    /** The reason why there is an isOutside function now, is that
     *  there exists voxelFlag::outside, voxelFlag::outerBorder flags,
     *  there are several flags which count as outside: undetermined,
     *  outside, and borderline. This is particularly important for the
     *  undetermined flag, because the outer envelopes in a multi-block
     *  structure have no special meaning and are default-initialized to
     *  undetermined. In any case isInside or isOutside will be false in this case.
     **/
    virtual bool isOutside(Dot3D const &location) const = 0;
    /// Get the distance to the wall, and the velocity value on the wall,
    ///   from a given real-valued position (in lattice units), and along
    ///   a given direction. Returns true if there exists an intersection
    ///   along the indicated direction.
    ///   ATTENTION: The id is an in-and-out value. If you know the right
    ///   id of the surface intersection (e.g. the triangle ID in case of
    ///   a triangular mesh) you can provide it to speed up the calculation.
    ///   However, if you don't know it, you MUST provide the value -1,
    ///   because otherwise the result might be wrong.
    virtual bool pointOnSurface(
        Array<T, 3> const &fromPoint, Array<T, 3> const &direction, Array<T, 3> &locatedPoint,
        T &distance, Array<T, 3> &wallNormal, SurfaceData &surfaceData, OffBoundary::Type &bdType,
        plint &id) const = 0;
    /// Get the distance to the wall, and the data on the wall,
    ///   from a given discrete lattice node, and along a given direction.
    ///   Returns true if there exists an intersection along the indicated
    ///   direction.
    ///   ATTENTION: The id is an in-and-out value. If you know the right
    ///   id of the surface intersection (e.g. the triangle ID in case of
    ///   a triangular mesh) you can provide it to speed up the calculation.
    ///   However, if you don't know it, you MUST provide the value -1,
    ///   because otherwise the result might be wrong.
    virtual bool gridPointOnSurface(
        Dot3D const &fromPoint, Dot3D const &direction, Array<T, 3> &locatedPoint, T &distance,
        Array<T, 3> &wallNormal, SurfaceData &surfaceData, OffBoundary::Type &bdType,
        plint &id) const
    {
        return pointOnSurface(
            Array<T, 3>((T)fromPoint.x, (T)fromPoint.y, (T)fromPoint.z),
            Array<T, 3>((T)direction.x, (T)direction.y, (T)direction.z), locatedPoint, distance,
            wallNormal, surfaceData, bdType, id);
    }
    /// Given a point p on the surface of the shape, determine its "continuous normal".
    ///   If the shape is for example piecewise linear, the normal is adjusted to vary
    ///   continuously over the surface.
    virtual Array<T, 3> computeContinuousNormal(
        Array<T, 3> const &p, plint id, bool isAreaWeighted = false) const = 0;
    /// Say if a given segment intersects the surface.
    virtual bool intersectsSurface(
        Array<T, 3> const &p1, Array<T, 3> const &p2, plint &id) const = 0;
    /// Say if a given segment with integer-valued endpoints intersects the surface.
    virtual bool gridIntersectsSurface(Dot3D const &p1, Dot3D const &p2, plint &id) const
    {
        return intersectsSurface(
            Array<T, 3>((T)p1.x, (T)p1.y, (T)p1.z), Array<T, 3>((T)p2.x, (T)p2.y, (T)p2.z), id);
    }
    /// Get the shortest distance to the wall. Returns true in case of success.
    ///   The flag isBehind indicates if the point is "behind" the wall, i.e.
    ///   in the direction opposite to the wall normal.
    virtual bool distanceToSurface(Array<T, 3> const &point, T &distance, bool &isBehind) const = 0;
    /// Get the shortest distance to the wall. Returns true in case of success.
    ///   The flag isBehind indicates if the point is "behind" the wall, i.e.
    ///   in the direction opposite to the wall normal.
    virtual bool gridDistanceToSurface(Dot3D const &point, T &distance, bool &isBehind) const
    {
        return distanceToSurface(
            Array<T, 3>((T)point.x, (T)point.y, (T)point.z), distance, isBehind);
    }
    /// Return the tag (id of boundary-portion with specific boundary,
    ///   condition); the id is the one returned by pointOnSurface.
    virtual plint getTag(plint id) const = 0;
    /// Plain clone function.
    virtual BoundaryShape3D<T, SurfaceData> *clone() const = 0;
    /// In case the shape class needs additional meshed data, this clone function
    ///   offers the possibility to provide the data.
    virtual BoundaryShape3D<T, SurfaceData> *clone(std::vector<AtomicBlock3D *> args) const
    {
        return clone();
    }
};

template <typename T, class SurfaceData>
class BoundaryShapeIntersection3D : public BoundaryShape3D<T, SurfaceData> {
public:
    BoundaryShapeIntersection3D();
    virtual ~BoundaryShapeIntersection3D();
    BoundaryShapeIntersection3D(BoundaryShapeIntersection3D<T, SurfaceData> const &rhs);
    BoundaryShapeIntersection3D<T, SurfaceData> &operator=(
        BoundaryShapeIntersection3D<T, SurfaceData> const &rhs);
    void swap(BoundaryShapeIntersection3D<T, SurfaceData> &rhs);
    void addComponent(BoundaryShape3D<T, SurfaceData> *component);

    virtual bool isInside(Array<T, 3> const &location) const;
    virtual bool isInside(Array<T, 3> const &location, T epsilon) const;
    virtual bool pointOnSurface(
        Array<T, 3> const &fromPoint, Array<T, 3> const &direction, Array<T, 3> &locatedPoint,
        T &distance, Array<T, 3> &wallNormal, SurfaceData &surfaceData, OffBoundary::Type &bdType,
        plint &id) const;
    virtual bool intersectsSurface(Array<T, 3> const &p1, Array<T, 3> const &p2, plint &id) const;
    virtual bool distanceToSurface(Array<T, 3> const &point, T &distance, bool &isBehind) const;
    virtual plint getTag(plint id) const;
    virtual BoundaryShapeIntersection3D<T, SurfaceData> *clone() const;

private:
    std::vector<BoundaryShape3D<T, SurfaceData> *> components;
};

}  // namespace plb

#endif  // BOUNDARY_SHAPES_3D_H
