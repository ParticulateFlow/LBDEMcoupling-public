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

#ifndef TRIANGLE_SELECTOR_H
#define TRIANGLE_SELECTOR_H

#include <vector>

#include "core/array.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "multiBlock/multiBlockManagement3D.h"

namespace plb {

template <typename T>
struct TriangleSelector {
    virtual ~TriangleSelector() { }
    virtual bool operator()(Array<Array<T, 3>, 3> const &triangle, plint partId) const = 0;
    virtual TriangleSelector<T> *clone() const = 0;
};

// This selector chooses triangles that intersect with a given cuboid.
// It can be called from inside a data processor (with the cuboid
// being the domain provided to the data processor for example).
// This selector is "non-unique" in the sense that a given triangle
// may intersect with two neighbor cuboids, so it will be selected by
// both of them.
// Care must be taken so that the cuboid and the triangles are at
// the same system of units.
template <typename T>
class CuboidTriangleSelector : public TriangleSelector<T> {
public:
    CuboidTriangleSelector(Cuboid<T> const &cuboid_) : cuboid(cuboid_) { }

    virtual bool operator()(Array<Array<T, 3>, 3> const &triangle, plint partId) const
    {
        return doesIntersect(cuboid, triangle);
    }

    virtual CuboidTriangleSelector<T> *clone() const
    {
        return new CuboidTriangleSelector<T>(*this);
    }

private:
    Cuboid<T> cuboid;
};

// This selector chooses triangles with centroids that belong (uniquely) to a given cuboid.
// It can be called from inside a data processor (with the cuboid being the domain provided
// to the data processor enlarged by 0.5).
// Care must be taken so that the cuboid and the triangles are at the same system of units.
template <typename T>
class CuboidUniqueTriangleSelector : public TriangleSelector<T> {
public:
    CuboidUniqueTriangleSelector(Cuboid<T> const &cuboid_) :
        cuboid(cuboid_), oneOverThree((T)1 / (T)3)
    { }

    virtual bool operator()(Array<Array<T, 3>, 3> const &triangle, plint partId) const
    {
        Array<T, 3> centroid = oneOverThree * (triangle[0] + triangle[1] + triangle[2]);
        return isContained(centroid, cuboid);
    }

    virtual CuboidUniqueTriangleSelector<T> *clone() const
    {
        return new CuboidUniqueTriangleSelector<T>(*this);
    }

private:
    bool isContained(Array<T, 3> const &p, Cuboid<T> const &c) const
    {
        return (p[0] > c.x0()) && (p[0] <= c.x1()) && (p[1] > c.y0()) && (p[1] <= c.y1())
               && (p[2] > c.z0()) && (p[2] <= c.z1());
    }

private:
    Cuboid<T> cuboid;
    T oneOverThree;
};

// This selector chooses triangles that belong to at least one of the local blocks.
// The triangles are chosen based on their intersection with these blocks.
// This means that the same triangle can be selected by more than one threads.
// The blocks are enlarged by "width" before the intersection test with the triangles
// takes place (it is a kind of envelopeWidth for the selection of triangles).
// In general, "width" should be greater than or equal to 1.0.
template <typename T>
class LocalTriangleSelector : public TriangleSelector<T> {
public:
    LocalTriangleSelector(MultiBlockManagement3D const &management_, T width_, T deltaX_ = (T)1) :
        management(management_), width(width_), scale((T)1 / deltaX_), offset(Array<T, 3>::zero())
    {
        if (util::isOne(scale)) {
            scale = (T)1;
        }
        if (util::isZero(norm(offset))) {
            offset = Array<T, 3>::zero();
        }
    }

    LocalTriangleSelector(
        MultiBlockManagement3D const &management_, T width_, T deltaX_,
        Array<T, 3> const &offset_) :
        management(management_), width(width_), scale((T)1 / deltaX_), offset(offset_)
    {
        if (util::isOne(scale)) {
            scale = (T)1;
        }
        if (util::isZero(norm(offset))) {
            offset = Array<T, 3>::zero();
        }
    }

    virtual bool operator()(Array<Array<T, 3>, 3> const &triangle, plint partId) const
    {
        SparseBlockStructure3D const &sparseBlock = management.getSparseBlockStructure();
        ThreadAttribution const &attribution = management.getThreadAttribution();

        Array<Array<T, 3>, 3> triangleLU = toLatticeUnits(triangle);

        Box3D triangleBoundingBox(boxCover(computeBoundingCuboid(triangleLU)));
        std::vector<plint> neighbors;
        plint neighborhoodWidth = width + 1;  // To be on the safe side.
        sparseBlock.findNeighbors(triangleBoundingBox, neighborhoodWidth, neighbors);
        for (pluint iNeighbor = 0; iNeighbor < neighbors.size(); iNeighbor++) {
            plint blockId = neighbors[iNeighbor];
            if (attribution.isLocal(blockId)) {
                Box3D bulk;
                sparseBlock.getBulk(blockId, bulk);
                Cuboid<T> extendedBulk(Cuboid<T>(bulk).enlarge(width));
                if (doesIntersect(extendedBulk, triangleLU)) {
                    return true;
                }
            }
        }

        return false;
    }

    virtual LocalTriangleSelector<T> *clone() const
    {
        return new LocalTriangleSelector<T>(*this);
    }

private:
    Array<Array<T, 3>, 3> toLatticeUnits(Array<Array<T, 3>, 3> const &triangle) const
    {
        Array<Array<T, 3>, 3> newTriangle;
        newTriangle[0] = scale * (triangle[0] - offset);
        newTriangle[1] = scale * (triangle[1] - offset);
        newTriangle[2] = scale * (triangle[2] - offset);
        return newTriangle;
    }

private:
    MultiBlockManagement3D const &management;
    T width;
    T scale;
    Array<T, 3> offset;
};

// This selector chooses triangles with centroids that belong (uniquely) to one of the local blocks.
template <typename T>
class LocalUniqueTriangleSelector : public TriangleSelector<T> {
public:
    LocalUniqueTriangleSelector(MultiBlockManagement3D const &management_, T deltaX_ = (T)1) :
        management(management_),
        scale((T)1 / deltaX_),
        offset(Array<T, 3>::zero()),
        oneOverThree((T)1 / (T)3)
    {
        if (util::isOne(scale)) {
            scale = (T)1;
        }
        if (util::isZero(norm(offset))) {
            offset = Array<T, 3>::zero();
        }
    }

    LocalUniqueTriangleSelector(
        MultiBlockManagement3D const &management_, T deltaX_, Array<T, 3> const &offset_) :
        management(management_), scale((T)1 / deltaX_), offset(offset_), oneOverThree((T)1 / (T)3)
    {
        if (util::isOne(scale)) {
            scale = (T)1;
        }
        if (util::isZero(norm(offset))) {
            offset = Array<T, 3>::zero();
        }
    }

    virtual bool operator()(Array<Array<T, 3>, 3> const &triangle, plint partId) const
    {
        SparseBlockStructure3D const &sparseBlock = management.getSparseBlockStructure();
        ThreadAttribution const &attribution = management.getThreadAttribution();

        Array<Array<T, 3>, 3> triangleLU = toLatticeUnits(triangle);
        Array<T, 3> centroidLU = oneOverThree * (triangleLU[0] + triangleLU[1] + triangleLU[2]);

        Box3D triangleBoundingBox(boxCover(computeBoundingCuboid(triangleLU)));
        std::vector<plint> neighbors;
        plint neighborhoodWidth = 1;  // To be on the safe side.
        sparseBlock.findNeighbors(triangleBoundingBox, neighborhoodWidth, neighbors);
        for (pluint iNeighbor = 0; iNeighbor < neighbors.size(); iNeighbor++) {
            plint blockId = neighbors[iNeighbor];
            if (attribution.isLocal(blockId)) {
                Box3D bulk;
                sparseBlock.getBulk(blockId, bulk);
                if (isContained(centroidLU, bulk)) {
                    return true;
                }
            }
        }

        return false;
    }

    virtual LocalUniqueTriangleSelector<T> *clone() const
    {
        return new LocalUniqueTriangleSelector<T>(*this);
    }

private:
    Array<Array<T, 3>, 3> toLatticeUnits(Array<Array<T, 3>, 3> const &triangle) const
    {
        Array<Array<T, 3>, 3> newTriangle;
        newTriangle[0] = scale * (triangle[0] - offset);
        newTriangle[1] = scale * (triangle[1] - offset);
        newTriangle[2] = scale * (triangle[2] - offset);
        return newTriangle;
    }

    bool isContained(Array<T, 3> const &p, Box3D const &b) const
    {
        return (p[0] > (T)b.x0 - (T)0.5) && (p[0] <= (T)b.x1 + (T)0.5) && (p[1] > (T)b.y0 - (T)0.5)
               && (p[1] <= (T)b.y1 + (T)0.5) && (p[2] > (T)b.z0 - (T)0.5)
               && (p[2] <= (T)b.z1 + (T)0.5);
    }

private:
    MultiBlockManagement3D const &management;
    T scale;
    Array<T, 3> offset;
    T oneOverThree;
};

// This selector chooses triangles that belong to a given block. The triangles are chosen
// based on their intersection with the block. This means that the same triangle can be
// selected by more than one blocks. The blocks are enlarged by "width" before the
// intersection test with the triangles takes place (it is a kind of envelopeWidth for the
// selection of triangles). In general, "width" should be greater than or equal to 1.0.
template <typename T>
class BlockTriangleSelector : public TriangleSelector<T> {
public:
    BlockTriangleSelector(
        SparseBlockStructure3D const &sparseBlockStructure_, plint blockId_, T width_,
        T deltaX_ = (T)1) :
        sparseBlockStructure(sparseBlockStructure_),
        blockId(blockId_),
        width(width_),
        scale((T)1 / deltaX_),
        offset(Array<T, 3>::zero())
    {
        if (util::isOne(scale)) {
            scale = (T)1;
        }
        if (util::isZero(norm(offset))) {
            offset = Array<T, 3>::zero();
        }
    }

    BlockTriangleSelector(
        SparseBlockStructure3D const &sparseBlockStructure_, plint blockId_, T width_, T deltaX_,
        Array<T, 3> const &offset_) :
        sparseBlockStructure(sparseBlockStructure_),
        blockId(blockId_),
        width(width_),
        scale((T)1 / deltaX_),
        offset(offset_)
    {
        if (util::isOne(scale)) {
            scale = (T)1;
        }
        if (util::isZero(norm(offset))) {
            offset = Array<T, 3>::zero();
        }
    }

    virtual bool operator()(Array<Array<T, 3>, 3> const &triangle, plint partId) const
    {
        Array<Array<T, 3>, 3> triangleLU = toLatticeUnits(triangle);

        Box3D bulk;
        sparseBlockStructure.getBulk(blockId, bulk);
        Cuboid<T> extendedBulk(Cuboid<T>(bulk).enlarge(width));
        if (doesIntersect(extendedBulk, triangleLU)) {
            return true;
        }

        return false;
    }

    virtual BlockTriangleSelector<T> *clone() const
    {
        return new BlockTriangleSelector<T>(*this);
    }

private:
    Array<Array<T, 3>, 3> toLatticeUnits(Array<Array<T, 3>, 3> const &triangle) const
    {
        Array<Array<T, 3>, 3> newTriangle;
        newTriangle[0] = scale * (triangle[0] - offset);
        newTriangle[1] = scale * (triangle[1] - offset);
        newTriangle[2] = scale * (triangle[2] - offset);
        return newTriangle;
    }

private:
    SparseBlockStructure3D const &sparseBlockStructure;
    plint blockId;
    T width;
    T scale;
    Array<T, 3> offset;
};

// This selector chooses triangles with centroids that belong (uniquely) to the given block.
template <typename T>
class BlockUniqueTriangleSelector : public TriangleSelector<T> {
public:
    BlockUniqueTriangleSelector(
        SparseBlockStructure3D const &sparseBlockStructure_, plint blockId_, T deltaX_ = (T)1) :
        sparseBlockStructure(sparseBlockStructure_),
        blockId(blockId_),
        scale((T)1 / deltaX_),
        offset(Array<T, 3>::zero()),
        oneOverThree((T)1 / (T)3)
    {
        if (util::isOne(scale)) {
            scale = (T)1;
        }
        if (util::isZero(norm(offset))) {
            offset = Array<T, 3>::zero();
        }
    }

    BlockUniqueTriangleSelector(
        SparseBlockStructure3D const &sparseBlockStructure_, plint blockId_, T deltaX_,
        Array<T, 3> const &offset_) :
        sparseBlockStructure(sparseBlockStructure_),
        blockId(blockId_),
        scale((T)1 / deltaX_),
        offset(offset_),
        oneOverThree((T)1 / (T)3)
    {
        if (util::isOne(scale)) {
            scale = (T)1;
        }
        if (util::isZero(norm(offset))) {
            offset = Array<T, 3>::zero();
        }
    }

    virtual bool operator()(Array<Array<T, 3>, 3> const &triangle, plint partId) const
    {
        Array<Array<T, 3>, 3> triangleLU = toLatticeUnits(triangle);
        Array<T, 3> centroidLU = oneOverThree * (triangleLU[0] + triangleLU[1] + triangleLU[2]);

        Box3D bulk;
        sparseBlockStructure.getBulk(blockId, bulk);
        if (isContained(centroidLU, bulk)) {
            return true;
        }

        return false;
    }

    virtual BlockUniqueTriangleSelector<T> *clone() const
    {
        return new BlockUniqueTriangleSelector<T>(*this);
    }

private:
    Array<Array<T, 3>, 3> toLatticeUnits(Array<Array<T, 3>, 3> const &triangle) const
    {
        Array<Array<T, 3>, 3> newTriangle;
        newTriangle[0] = scale * (triangle[0] - offset);
        newTriangle[1] = scale * (triangle[1] - offset);
        newTriangle[2] = scale * (triangle[2] - offset);
        return newTriangle;
    }

    bool isContained(Array<T, 3> const &p, Box3D const &b) const
    {
        return (p[0] > (T)b.x0 - (T)0.5) && (p[0] <= (T)b.x1 + (T)0.5) && (p[1] > (T)b.y0 - (T)0.5)
               && (p[1] <= (T)b.y1 + (T)0.5) && (p[2] > (T)b.z0 - (T)0.5)
               && (p[2] <= (T)b.z1 + (T)0.5);
    }

private:
    SparseBlockStructure3D const &sparseBlockStructure;
    plint blockId;
    T scale;
    Array<T, 3> offset;
    T oneOverThree;
};

}  // namespace plb

#endif  // TRIANGLE_SELECTOR_H
