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

#ifndef BOX_LOGIC_3D_H
#define BOX_LOGIC_3D_H

#include <string>
#include <vector>

#include "core/array.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "gridRefinement/octree.h"
#include "gridRefinement/octreeGridGenerator.h"

namespace plb {

namespace boxLogic {

// A Plane is defined by a Box3D which has
// two coordinates which are equal (either x0=x1, y0=y1, or z0=z1)
struct Plane {
    Plane(Box3D bb_) : bb(bb_)
    {
        PLB_ASSERT(
            ((bb.x0 == bb.x1) || (bb.y0 == bb.y1) || (bb.z0 == bb.z1))
            && "The box must define a plane therefore two coordinates must be the same.");
    }

    std::string toStr() const
    {
        return "Plane = (" + util::val2str(bb.x0) + ", " + util::val2str(bb.x1) + ", "
               + util::val2str(bb.y0) + ", " + util::val2str(bb.y1) + ", " + util::val2str(bb.z0)
               + ", " + util::val2str(bb.z1) + ")";
    }

    plint nCells() const
    {
        return bb.nCells();
    }

    bool operator==(Plane const &rhs) const
    {
        return bb == rhs.bb;
    }

    Box3D bb;
};

// A directed plane is defined by a Plane and its normal vector
// (the nomral is along
// one of thed principal axes: x=0, y=1, z=2).
struct DirectedPlane : public Plane {
    DirectedPlane(Box3D bb_, int direction_, int orientation_) :
        Plane(bb_), direction(direction_), orientation(orientation_)
    {
        PLB_ASSERT(
            (direction == 0 || direction == 1 || direction == 2)
            && "Direction must be 0, 1, or 2 (x, y, z).");

        PLB_ASSERT((orientation == -1 || orientation == +1) && "Orientation must be -1 or +1.");
    }

    std::string toStr() const
    {
        return Plane::toStr() + " , direction = " + util::val2str(direction)
               + " , orientation = " + util::val2str(orientation) + ".";
    }

    plint nCells() const
    {
        return bb.nCells();
    }

    bool operator==(DirectedPlane const &rhs) const
    {
        return bb == rhs.bb && direction == rhs.direction && orientation == rhs.orientation;
    }

    int direction, orientation;
};

// An Edge is defined by a Box3D which has
// four coordinates which are equal (either x0=x1 && y0=y1,
// or x0=x1 && z0=z1, or y0=y1 && z0=z1)
struct Edge {
    Edge(Box3D bb_) : bb(bb_)
    {
        PLB_ASSERT(
            ((bb.x0 == bb.x1 && bb.y0 == bb.y1) || (bb.x0 == bb.x1 && bb.z0 == bb.z1)
             || (bb.y0 == bb.y1 && bb.z0 == bb.z1))
            && "The box must define an edge therefore 4 coordinates must be the same.");
    }

    std::string toStr() const
    {
        return "Edge = (" + util::val2str(bb.x0) + ", " + util::val2str(bb.x1) + ", "
               + util::val2str(bb.y0) + ", " + util::val2str(bb.y1) + ", " + util::val2str(bb.z0)
               + ", " + util::val2str(bb.z1) + ")";
    }

    plint nCells() const
    {
        return bb.nCells();
    }

    bool operator==(Edge const &rhs) const
    {
        return bb == rhs.bb;
    }

    Box3D bb;
};

// A directed Edge is an edge where one defines the plane it is adjacent to (planeDir).
// Furthermore information is given about its "direction" which is the dir which
// is null. Finally the non null dir is the orientation pointing on the "outside"
// of the edge.
struct DirectedEdge : public Edge {
    DirectedEdge(Box3D bb_, int planeDir_, int dir1_, int dir2_) :
        Edge(bb_), planeDir(planeDir_), dir1(dir1_), dir2(dir2_)
    {
        PLB_ASSERT(
            (planeDir == 0 || planeDir == 1 || planeDir == 2)
            && "The planeDir must be 0, 1, or 2.");
        PLB_ASSERT((dir1 == 0 || dir2 == 0) && "The dir1 or dir2 must be 0.");
        PLB_ASSERT(
            ((dir1 == -1 || dir1 == +1) || (dir2 == -1 || dir2 == +1))
            && "The dir1 or dir2 must be -1 or +1.");
    }

    plint getDirectionAlongEdge()
    {
        if (planeDir == 0) {
            if (dir1 == 0) {
                return 1;
            } else {
                return 2;
            }
        } else if (planeDir == 1) {
            if (dir1 == 0) {
                return 2;
            } else {
                return 0;
            }

        } else if (planeDir == 2) {
            if (dir1 == 0) {
                return 0;
            } else {
                return 1;
            }
        } else {
            PLB_ASSERT(false && "Dir must enter the if statement otherwise edge misformed.")
            return -1;
        }
    }

    std::string toStr() const
    {
        return Edge::toStr() + " , planeDir = " + util::val2str(planeDir)
               + " , dir1 = " + util::val2str(dir1) + " , dir2 = " + util::val2str(dir2);
    }

    plint nCells() const
    {
        return bb.nCells();
    }

    bool operator==(DirectedEdge const &rhs) const
    {
        return bb == rhs.bb && planeDir == rhs.planeDir && dir1 == rhs.dir1 && dir2 == rhs.dir2;
    }

    int planeDir, dir1, dir2;
};

// A Corner is defined by a Box3D which has
// six coordinates which are equal (or x0=x1 && y0=y1 && z0=z1)
struct Corner {
    Corner(Box3D bb_) : bb(bb_)
    {
        PLB_ASSERT(
            (bb.x0 == bb.x1 && bb.y0 == bb.y1 && bb.z0 == bb.z1)
            && "The box must define a corner therefore all coordinates must be the same.");
    }

    std::string toStr() const
    {
        return "Corner = (" + util::val2str(bb.x0) + ", " + util::val2str(bb.x1) + ", "
               + util::val2str(bb.y0) + ", " + util::val2str(bb.y1) + ", " + util::val2str(bb.z0)
               + ", " + util::val2str(bb.z1) + ")";
    }

    plint nCells() const
    {
        return bb.nCells();
    }

    bool operator==(Corner const &rhs) const
    {
        return bb == rhs.bb;
    }

    Box3D bb;
};

// A DirectedCorner is a Corner where one defines the plane it is adjacent to (planeDir).
// The directions are the orientations pointing on the "outside" of the corner.
struct DirectedCorner : public Corner {
    DirectedCorner(Box3D bb_, int planeDir_, int dir1_, int dir2_) :
        Corner(bb_), planeDir(planeDir_), dir1(dir1_), dir2(dir2_)
    {
        PLB_ASSERT(
            (dir1 == -1 || dir1 == +1) && (dir2 == -1 || dir2 == +1)
            && "The dir1 and dir2 must be -1 or +1.");
    }

    std::string toStr() const
    {
        return Corner::toStr() + " , planeDir = " + util::val2str(planeDir)
               + " , dir1 = " + util::val2str(dir1) + " , dir2 = " + util::val2str(dir2);
    }

    plint nCells() const
    {
        return bb.nCells();
    }

    bool operator==(DirectedCorner const &rhs) const
    {
        return bb == rhs.bb && planeDir == rhs.planeDir && dir1 == rhs.dir1 && dir2 == rhs.dir2;
    }

    int planeDir, dir1, dir2;
};

// does a vector of boxes have two identical boxes?
// if the size of the box is of 1, then an exception can be made
template <class Boundary>
bool containsDuplicates(const std::vector<Boundary> &boxes, bool exceptSizeOne);
template <class Box>
bool containsNoDuplicates(const std::vector<Box> &boxes, bool exceptSizeOne);
template <class Surf>
std::vector<Surf> merge(const std::vector<Surf> &origBoxes);
/// If the two planes can be merged into one, do it, and write the result
///   into the first plane. Return value states if merging was successful.
bool merge(Plane &p1, Plane const &p2);
/// If the two edges can be merged into one, do it, and write the result
///   into the first edge. Return value states if merging was successful.
bool merge(Edge &e1, Edge const &e2);
/// If the two corners can be merged into one, do it, and write the result
///   into the first corner. Return value states if merging was successful.
bool merge(Corner &c1, Corner const &c2);
/// If the two directed planes can be merged into one, do it, and write the result
///   into the first directed plane. Return value states if merging was successful.
bool merge(DirectedPlane &p1, DirectedPlane const &p2);
/// If the two directed edges can be merged into one, do it, and write the result
///   into the first directed edge. Return value states if merging was successful.
bool merge(DirectedEdge &e1, DirectedEdge const &e2);
/// If the two directed corners can be merged into one, do it, and write the result
///   into the first directed corner. Return value states if merging was successful.
bool merge(DirectedCorner &c1, DirectedCorner const &c2);

/// The global ordering on boxes has no specific interpretation, but
///   is useful to hold boxes in ordered data structures like sets.
bool operator<(std::pair<plint, Box3D> const &box1, std::pair<plint, Box3D> const &box2);

// returns the corners of a box (excluding the corner points)
std::vector<Corner> getCorners(const Box3D &box);
// returns the edges of a box (excluding the corner points)
std::vector<Edge> getEdges(const Box3D &box);
bool thicknessOne(std::vector<Box3D> &b);
std::vector<Box3D> removeOverlapsOfBoxes(const std::vector<Box3D> &origBoxes);
// returns the planes of a box (excluding the edges and corners)
// the neighbor array says if in one of the 3 axes and 2 directions there
// is a neighbor with the same grid level or no grid level (in this case neighbor is true)
std::vector<DirectedPlane> getDirectedPlanes(
    const Box3D &box, const Array<bool, 26> &neighbor, const Array<bool, 26> &bcNeighbor,
    const Array<bool, 26> &allocNeighbor);
// returns the planes of a box (excluding the edges and corners)
// the neighbor array says if in one of the 3 axes and 2 directions there
// is a neighbor with the same grid level or no grid level (in this case neighbor is true)
std::vector<DirectedPlane> getCompleteDirectedPlanes(
    const Box3D &box, const Array<bool, 26> &neighbor, const Array<bool, 26> &bcNeighbor,
    const Array<bool, 26> &allocNeighbor);
// returns the planes of a box (excluding the edges and corners)
std::vector<DirectedEdge> getDirectedEdges(
    const Box3D &box, const Array<bool, 26> &neighbor, const Array<bool, 26> &bcNeighbor,
    const Array<bool, 26> &allocNeighbor);
// returns the planes of a box (excluding the edges and corners)
std::vector<Edge> getBoundaryEdges(
    const Box3D &box, const Array<bool, 26> &neighbor, const Array<bool, 26> &bcNeighbor,
    const Array<bool, 26> &allocNeighbor);
// returns the planes of a box (excluding the edges and corners)
std::vector<DirectedEdge> getBoundaryDirectedEdges(
    const Box3D &box, const Array<bool, 26> &neighbor, const Array<bool, 26> &bcNeighbor,
    const Array<bool, 26> &allocNeighbor);
// returns the corners of a box (excluding the corner points)
std::vector<DirectedCorner> getDirectedCorners(
    const Box3D &box, const Array<bool, 26> &neighbor, const Array<bool, 26> &bcNeighbor,
    const Array<bool, 26> &allocNeighbor);
// returns the corners of a box (excluding the corner points)
std::vector<DirectedCorner> getBoundaryDirectedCorners(
    const Box3D &box, const Array<bool, 26> &neighbor, const Array<bool, 26> &bcNeighbor,
    const Array<bool, 26> &allocNeighbor);

bool hasFinerNeighbor(OctreeGridStructure const &ogs, plint blockId);

bool hasCoarserNeighbor(OctreeGridStructure const &ogs, plint blockId);

Array<bool, 26> getNeighbors(OctreeGridStructure const &ogs, plint blockId);

Array<bool, 26> getBcNeighbors(OctreeGridStructure const &ogs, plint blockId);

Array<bool, 26> getAllocatedNeighbors(OctreeGridStructure const &ogs, plint blockId);

std::pair<plint, plint> getDirections(plint dir);

std::pair<Array<plint, 2>, Array<plint, 2> > getEdgesIndices(
    plint dir, plint ori, const Array<bool, 26> &neighbors, const Array<bool, 26> &bcNeighbors,
    plint l1 = 1, plint l2 = 1);

void printBox(Box3D const &b, std::string name);

void printBoxes(
    std::vector<Box3D> const &b, std::string name, plint level, double dx,
    Array<double, 3> const &pos);

void makeIntersectionsConsistent(
    const Box3D &box, const std::vector<Box3D> &boxes, std::vector<Box3D> &intBoxes);

void getInterfaces(
    const Box3D &box, const Array<bool, 26> &neighbor, const Array<bool, 26> &bcNeighbor,
    const Array<bool, 26> &allocNeighbor, plint overlapWidthCoarseUnits,
    std::vector<boxLogic::Plane>
        &coarseToFinePlanesFineUnits,  // used for recompose of   fine interface
    std::vector<boxLogic::DirectedPlane>
        &coarseToFinePlanesCoarseUnits,  // used for interpolation Processing Functionals (PF)
    std::vector<boxLogic::DirectedEdge> &coarseToFineEdgesCoarseUnits,  // used for interpolation PF
    std::vector<boxLogic::DirectedCorner>
        &coarseToFineCornersCoarseUnits,  // used for interpolation PF
    std::vector<boxLogic::DirectedPlane>
        &coarseToFinePlanesCoarseUnitsExtended,       // used for decompose of coarse interface
    std::vector<Box3D> &fineToCoarsePlanesFineUnits,  // used for decompose of fine interfare
    std::vector<boxLogic::Plane>
        &fineToCoarsePlanesCoarseUnits,  // used for filtering and recompose of coarse interface

    std::vector<boxLogic::DirectedEdge>
        &coarseToFineBoundaryEdgesCoarseUnits,  // special inepolation for BC edge
    std::vector<boxLogic::DirectedCorner>
        &coarseToFineBoundaryCornersCoarseUnits,  // special inepolation for BC corner
    std::vector<boxLogic::Edge> &fineToCoarseBoundaryEdgesCoarseUnits  // used for interpolation PF)
);

}  // namespace boxLogic

}  // namespace plb

#endif  // BOX_LOGIC_3D_H
