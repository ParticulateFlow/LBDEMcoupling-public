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

#include "gridRefinement/boxLogic3D.h"

#include <algorithm>
#include <string>
#include <vector>

#include "core/array.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "gridRefinement/boxLogic3D.hh"
#include "offLattice/triangleSet.h"
#include "offLattice/triangleSet.hh"
#include "offLattice/triangleSetGenerator.h"
#include "offLattice/triangleSetGenerator.hh"

namespace plb {

namespace boxLogic {

typedef OctreeTables OT;

bool merge(Plane &p1, Plane const &p2)
{
    return merge(p1.bb, p2.bb);
}

bool merge(Edge &e1, Edge const &e2)
{
    return merge(e1.bb, e2.bb);
}

bool merge(Corner &c1, Corner const &c2)
{
    return merge(c1.bb, c2.bb);
}

bool merge(DirectedPlane &p1, DirectedPlane const &p2)
{
    if ((p1.direction == p2.direction) && (p1.orientation == p2.orientation)) {
        return merge(p1.bb, p2.bb);
    }
    return false;
}

bool merge(DirectedEdge &e1, DirectedEdge const &e2)
{
    if ((e1.planeDir == e2.planeDir) && (e1.dir1 == e2.dir1) && (e1.dir2 == e2.dir2)) {
        return merge(e1.bb, e2.bb);
    }
    return false;
}

bool merge(DirectedCorner &c1, DirectedCorner const &c2)
{
    if ((c1.planeDir == c2.planeDir) && (c1.dir1 == c2.dir1) && (c1.dir2 == c2.dir2)) {
        return merge(c1.bb, c2.bb);
    }
    return false;
}

bool operator<(std::pair<plint, Box3D> const &box1, std::pair<plint, Box3D> const &box2)
{
    return (box1.second < box2.second);
}

// returns the corners of a box (excluding the corner points)
std::vector<Corner> getCorners(const Box3D &box)
{
    std::vector<Corner> corners;

    corners.push_back(Corner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0)));
    corners.push_back(Corner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1)));
    corners.push_back(Corner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0)));
    corners.push_back(Corner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1)));
    corners.push_back(Corner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0)));
    corners.push_back(Corner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1)));
    corners.push_back(Corner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0)));
    corners.push_back(Corner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1)));

    PLB_ASSERT(containsNoDuplicates(corners, false) && "Corners must be unique.");
    PLB_ASSERT(corners.size() == 8 && "There are 8 corners in a box.");

    return corners;
}

// returns the edges of a box (excluding the corner points)
std::vector<Edge> getEdges(const Box3D &box)
{
    std::vector<Edge> edges;

    // z edges (must be 4 of them)
    edges.push_back(Edge(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0 + 1, box.z1 - 1)));
    edges.push_back(Edge(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0 + 1, box.z1 - 1)));
    edges.push_back(Edge(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0 + 1, box.z1 - 1)));
    edges.push_back(Edge(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0 + 1, box.z1 - 1)));

    // y edges (must be 4 of them)
    edges.push_back(Edge(Box3D(box.x0, box.x0, box.y0 + 1, box.y1 - 1, box.z0, box.z0)));
    edges.push_back(Edge(Box3D(box.x0, box.x0, box.y0 + 1, box.y1 - 1, box.z1, box.z1)));
    edges.push_back(Edge(Box3D(box.x1, box.x1, box.y0 + 1, box.y1 - 1, box.z0, box.z0)));
    edges.push_back(Edge(Box3D(box.x1, box.x1, box.y0 + 1, box.y1 - 1, box.z1, box.z1)));

    // x edges (must be 4 of them)
    edges.push_back(Edge(Box3D(box.x0 + 1, box.x1 - 1, box.y0, box.y0, box.z0, box.z0)));
    edges.push_back(Edge(Box3D(box.x0 + 1, box.x1 - 1, box.y0, box.y0, box.z1, box.z1)));
    edges.push_back(Edge(Box3D(box.x0 + 1, box.x1 - 1, box.y1, box.y1, box.z0, box.z0)));
    edges.push_back(Edge(Box3D(box.x0 + 1, box.x1 - 1, box.y1, box.y1, box.z1, box.z1)));

    PLB_ASSERT(containsNoDuplicates(edges, false) && "Edges must be unique.");
    PLB_ASSERT(edges.size() == 12 && "There are 12 edges in a box.");

    return edges;
}

bool thicknessOne(std::vector<Box3D> &b)
{
    for (std::vector<Box3D>::iterator it = b.begin(); it != b.end(); ++it) {
        bool one = ((it->x0 == it->x1) || (it->y0 == it->y1) || (it->z0 == it->z1));
        if (one) {
            return true;
        }
    }
    return false;
}

std::vector<Box3D> removeOverlapsOfBoxes(const std::vector<Box3D> &origBoxes)
{
    std::vector<Box3D> result;

    result = origBoxes;
    plint iA = 0;
    bool allIntRemoved = false;

    while (!allIntRemoved && origBoxes.size() > 0) {
        // pcout << "before = " << iA << ", " << result.size() << std::endl;

        std::vector<Box3D> newBoxes;
        newBoxes.insert(newBoxes.begin(), result.begin(), result.begin() + iA + 1);
        for (plint iB = iA + 1; iB < (plint)result.size(); ++iB) {
            std::vector<Box3D> resultTmp;
            except(result[iB], result[iA], resultTmp);
            // if (thicknessOne(resultTmp)) {
            // resultTmp.clear();
            // except(result[iA], result[iB], resultTmp);
            // }
            newBoxes.insert(newBoxes.end(), resultTmp.begin(), resultTmp.end());
        }
        result = newBoxes;
        ++iA;
        // pcout << "after = " << iA << ", " << result.size() << std::endl;
        if ((plint)result.size() == iA) {
            allIntRemoved = true;
        }
    }
    return result;
}

std::vector<DirectedPlane> getDirectedPlanes(
    const Box3D &box, const Array<bool, 26> &neighbor, const Array<bool, 26> &bcNeighbor,
    const Array<bool, 26> &allocNeighbor)
{
    std::vector<DirectedPlane> planes;

    // x planes (must be 2 of them)
    if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
        && !bcNeighbor[OT::surface0N()]) {
        planes.push_back(DirectedPlane(
            Box3D(box.x0, box.x0, box.y0 + 1, box.y1 - 1, box.z0 + 1, box.z1 - 1), 0, -1));
    }
    if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
        && !bcNeighbor[OT::surface0P()]) {
        planes.push_back(DirectedPlane(
            Box3D(box.x1, box.x1, box.y0 + 1, box.y1 - 1, box.z0 + 1, box.z1 - 1), 0, +1));
    }

    // y planes (must be 2 of them)
    if ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
        && !bcNeighbor[OT::surface1N()]) {
        planes.push_back(DirectedPlane(
            Box3D(box.x0 + 1, box.x1 - 1, box.y0, box.y0, box.z0 + 1, box.z1 - 1), 1, -1));
    }
    if ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
        && !bcNeighbor[OT::surface1P()]) {
        planes.push_back(DirectedPlane(
            Box3D(box.x0 + 1, box.x1 - 1, box.y1, box.y1, box.z0 + 1, box.z1 - 1), 1, +1));
    }

    // z planes (must be 2 of them)
    if ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
        && !bcNeighbor[OT::surface2N()]) {
        planes.push_back(DirectedPlane(
            Box3D(box.x0 + 1, box.x1 - 1, box.y0 + 1, box.y1 - 1, box.z0, box.z0), 2, -1));
    }
    if ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
        && !bcNeighbor[OT::surface2P()]) {
        planes.push_back(DirectedPlane(
            Box3D(box.x0 + 1, box.x1 - 1, box.y0 + 1, box.y1 - 1, box.z1, box.z1), 2, +1));
    }

    return planes;
}

std::vector<DirectedPlane> getCompleteDirectedPlanes(
    const Box3D &box, const Array<bool, 26> &neighbor, const Array<bool, 26> &bcNeighbor,
    const Array<bool, 26> &allocNeighbor)
{
    std::vector<DirectedPlane> planes;

    // x planes (must be 2 of them)
    if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
        && !bcNeighbor[OT::surface0N()]) {
        std::vector<Box3D> result;
        result.push_back(Box3D(box.x0, box.x0, box.y0, box.y1, box.z0, box.z1));

        if ((!neighbor[OT::edge0NN()] && allocNeighbor[OT::edge0NN()]) && neighbor[OT::surface1N()]
            && neighbor[OT::surface2N()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x0, box.x0, box.y0, box.y0 + 1, box.z0, box.z0 + 1),
                    resultTmp);
            }
            result = resultTmp;
        }

        if ((!neighbor[OT::edge0NP()] && allocNeighbor[OT::edge0NP()]) && neighbor[OT::surface1N()]
            && neighbor[OT::surface2P()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x0, box.x0, box.y0, box.y0 + 1, box.z1 - 1, box.z1),
                    resultTmp);
            }
            result = resultTmp;
        }

        if ((!neighbor[OT::edge0PN()] && allocNeighbor[OT::edge0PN()]) && neighbor[OT::surface1P()]
            && neighbor[OT::surface2N()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x0, box.x0, box.y1 - 1, box.y1, box.z0, box.z0 + 1),
                    resultTmp);
            }
            result = resultTmp;
        }

        if ((!neighbor[OT::edge0PP()] && allocNeighbor[OT::edge0PP()]) && neighbor[OT::surface1P()]
            && neighbor[OT::surface2P()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x0, box.x0, box.y1 - 1, box.y1, box.z1 - 1, box.z1),
                    resultTmp);
            }
            result = resultTmp;
        }

        for (plint iA = 0; iA < (plint)result.size(); ++iA) {
            planes.push_back(DirectedPlane(result[iA], 0, -1));
        }
    }
    if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
        && !bcNeighbor[OT::surface0P()]) {
        std::vector<Box3D> result;
        result.push_back(Box3D(box.x1, box.x1, box.y0, box.y1, box.z0, box.z1));

        if ((!neighbor[OT::edge0NN()] && allocNeighbor[OT::edge0NN()]) && neighbor[OT::surface1N()]
            && neighbor[OT::surface2N()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x1, box.x1, box.y0, box.y0 + 1, box.z0, box.z0 + 1),
                    resultTmp);
            }
            result = resultTmp;
        }

        if ((!neighbor[OT::edge0NP()] && allocNeighbor[OT::edge0NP()]) && neighbor[OT::surface1N()]
            && neighbor[OT::surface2P()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x1, box.x1, box.y0, box.y0 + 1, box.z1 - 1, box.z1),
                    resultTmp);
            }
            result = resultTmp;
        }

        if ((!neighbor[OT::edge0PN()] && allocNeighbor[OT::edge0PN()]) && neighbor[OT::surface1P()]
            && neighbor[OT::surface2N()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x1, box.x1, box.y1 - 1, box.y1, box.z0, box.z0 + 1),
                    resultTmp);
            }
            result = resultTmp;
        }

        if ((!neighbor[OT::edge0PP()] && allocNeighbor[OT::edge0PP()]) && neighbor[OT::surface1P()]
            && neighbor[OT::surface2P()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x1, box.x1, box.y1 - 1, box.y1, box.z1 - 1, box.z1),
                    resultTmp);
            }
            result = resultTmp;
        }

        for (plint iA = 0; iA < (plint)result.size(); ++iA) {
            planes.push_back(DirectedPlane(result[iA], 0, +1));
        }
    }

    // y planes (must be 2 of them)
    if ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
        && !bcNeighbor[OT::surface1N()]) {
        std::vector<Box3D> result;
        result.push_back(Box3D(box.x0, box.x1, box.y0, box.y0, box.z0, box.z1));

        if ((!neighbor[OT::edge1NN()] && allocNeighbor[OT::edge1NN()]) && neighbor[OT::surface0N()]
            && neighbor[OT::surface2N()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x0, box.x0 + 1, box.y0, box.y0, box.z0, box.z0 + 1),
                    resultTmp);
            }
            result = resultTmp;
        }

        if ((!neighbor[OT::edge1NP()] && allocNeighbor[OT::edge1NP()]) && neighbor[OT::surface0P()]
            && neighbor[OT::surface2N()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x1 - 1, box.x1, box.y0, box.y0, box.z0, box.z0 + 1),
                    resultTmp);
            }
            result = resultTmp;
        }

        if ((!neighbor[OT::edge1PN()] && allocNeighbor[OT::edge1PN()]) && neighbor[OT::surface0N()]
            && neighbor[OT::surface2P()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x0, box.x0, box.y0, box.y0, box.z1 - 1, box.z1),
                    resultTmp);
            }
            result = resultTmp;
        }

        if ((!neighbor[OT::edge1PP()] && allocNeighbor[OT::edge1PP()]) && neighbor[OT::surface0P()]
            && neighbor[OT::surface2P()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x1 - 1, box.x1, box.y0, box.y0, box.z1 - 1, box.z1),
                    resultTmp);
            }
            result = resultTmp;
        }

        for (plint iA = 0; iA < (plint)result.size(); ++iA) {
            planes.push_back(DirectedPlane(result[iA], 1, -1));
        }
    }
    if ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
        && !bcNeighbor[OT::surface1P()]) {
        std::vector<Box3D> result;
        result.push_back(Box3D(box.x0, box.x1, box.y1, box.y1, box.z0, box.z1));

        if ((!neighbor[OT::edge1NN()] && allocNeighbor[OT::edge1NN()]) && neighbor[OT::surface0N()]
            && neighbor[OT::surface2N()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x0, box.x0 + 1, box.y1, box.y1, box.z0, box.z0 + 1),
                    resultTmp);
            }
            result = resultTmp;
        }

        if ((!neighbor[OT::edge1NP()] && allocNeighbor[OT::edge1NP()]) && neighbor[OT::surface0P()]
            && neighbor[OT::surface2N()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x1 - 1, box.x1, box.y1, box.y1, box.z0, box.z0 + 1),
                    resultTmp);
            }
            result = resultTmp;
        }

        if ((!neighbor[OT::edge1PN()] && allocNeighbor[OT::edge1PN()]) && neighbor[OT::surface0N()]
            && neighbor[OT::surface2P()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x0, box.x0, box.y1, box.y1, box.z1 - 1, box.z1),
                    resultTmp);
            }
            result = resultTmp;
        }

        if ((!neighbor[OT::edge1PP()] && allocNeighbor[OT::edge1PP()]) && neighbor[OT::surface0P()]
            && neighbor[OT::surface2P()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x1 - 1, box.x1, box.y1, box.y1, box.z1 - 1, box.z1),
                    resultTmp);
            }
            result = resultTmp;
        }

        for (plint iA = 0; iA < (plint)result.size(); ++iA) {
            planes.push_back(DirectedPlane(result[iA], 1, +1));
        }
    }

    // z planes (must be 2 of them)
    if ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
        && !bcNeighbor[OT::surface2N()]) {
        std::vector<Box3D> result;
        result.push_back(Box3D(box.x0, box.x1, box.y0, box.y1, box.z0, box.z0));

        if ((!neighbor[OT::edge2NN()] && allocNeighbor[OT::edge2NN()]) && neighbor[OT::surface0N()]
            && neighbor[OT::surface1N()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x0, box.x0 + 1, box.y0, box.y0 + 1, box.z0, box.z0),
                    resultTmp);
            }
            result = resultTmp;
        }
        if ((!neighbor[OT::edge2NP()] && allocNeighbor[OT::edge2NP()]) && neighbor[OT::surface0N()]
            && neighbor[OT::surface1P()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x0, box.x0 + 1, box.y1 - 1, box.y1, box.z0, box.z0),
                    resultTmp);
            }
            result = resultTmp;
        }

        if ((!neighbor[OT::edge2PN()] && allocNeighbor[OT::edge2PN()]) && neighbor[OT::surface0P()]
            && neighbor[OT::surface1N()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x1 - 1, box.x1, box.y0, box.y0 + 1, box.z0, box.z0),
                    resultTmp);
            }
            result = resultTmp;
        }

        if ((!neighbor[OT::edge2PP()] && allocNeighbor[OT::edge2PP()]) && neighbor[OT::surface0P()]
            && neighbor[OT::surface1P()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x1 - 1, box.x1, box.y1 - 1, box.y1, box.z0, box.z0),
                    resultTmp);
            }
            result = resultTmp;
        }

        for (plint iA = 0; iA < (plint)result.size(); ++iA) {
            planes.push_back(DirectedPlane(result[iA], 2, -1));
        }
    }
    if ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
        && !bcNeighbor[OT::surface2P()]) {
        std::vector<Box3D> result;
        result.push_back(Box3D(box.x0, box.x1, box.y0, box.y1, box.z1, box.z1));

        if ((!neighbor[OT::edge2NN()] && allocNeighbor[OT::edge2NN()]) && neighbor[OT::surface0N()]
            && neighbor[OT::surface1N()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x0, box.x0 + 1, box.y0, box.y0 + 1, box.z1, box.z1),
                    resultTmp);
            }
            result = resultTmp;
        }
        if ((!neighbor[OT::edge2NP()] && allocNeighbor[OT::edge2NP()]) && neighbor[OT::surface0N()]
            && neighbor[OT::surface1P()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x0, box.x0 + 1, box.y1 - 1, box.y1, box.z1, box.z1),
                    resultTmp);
            }
            result = resultTmp;
        }

        if ((!neighbor[OT::edge2PN()] && allocNeighbor[OT::edge2PN()]) && neighbor[OT::surface0P()]
            && neighbor[OT::surface1N()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x1 - 1, box.x1, box.y0, box.y0 + 1, box.z1, box.z1),
                    resultTmp);
            }
            result = resultTmp;
        }

        if ((!neighbor[OT::edge2PP()] && allocNeighbor[OT::edge2PP()]) && neighbor[OT::surface0P()]
            && neighbor[OT::surface1P()])
        {
            std::vector<Box3D> resultTmp;
            for (plint iA = 0; iA < (plint)result.size(); ++iA) {
                except(
                    result[iA], Box3D(box.x1 - 1, box.x1, box.y1 - 1, box.y1, box.z1, box.z1),
                    resultTmp);
            }
            result = resultTmp;
        }

        for (plint iA = 0; iA < (plint)result.size(); ++iA) {
            planes.push_back(DirectedPlane(result[iA], 2, +1));
        }
    }

    return planes;
}

std::vector<DirectedEdge> getDirectedEdges(
    const Box3D &box, const Array<bool, 26> &neighbor, const Array<bool, 26> &bcNeighbor,
    const Array<bool, 26> &allocNeighbor)
{
    std::vector<DirectedEdge> edges;

    // ====================== External edges ====================== //
    // z edges (must be 8 of them)
    if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
        && (!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
        && !bcNeighbor[OT::surface0N()] && !bcNeighbor[OT::surface1N()])
    {
        edges.push_back(
            DirectedEdge(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0 + 1, box.z1 - 1), 0, -1, 0));
        edges.push_back(
            DirectedEdge(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0 + 1, box.z1 - 1), 1, 0, -1));
    }
    if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
        && (!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
        && !bcNeighbor[OT::surface0N()] && !bcNeighbor[OT::surface1P()])
    {
        edges.push_back(
            DirectedEdge(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0 + 1, box.z1 - 1), 0, +1, 0));
        edges.push_back(
            DirectedEdge(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0 + 1, box.z1 - 1), 1, 0, -1));
    }

    if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
        && (!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
        && !bcNeighbor[OT::surface0P()] && !bcNeighbor[OT::surface1N()])
    {
        edges.push_back(
            DirectedEdge(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0 + 1, box.z1 - 1), 0, -1, 0));
        edges.push_back(
            DirectedEdge(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0 + 1, box.z1 - 1), 1, 0, +1));
    }

    if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
        && (!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
        && !bcNeighbor[OT::surface0P()] && !bcNeighbor[OT::surface1P()])
    {
        edges.push_back(
            DirectedEdge(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0 + 1, box.z1 - 1), 0, +1, 0));
        edges.push_back(
            DirectedEdge(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0 + 1, box.z1 - 1), 1, 0, +1));
    }

    // y edges (must be 8 of them)
    if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
        && (!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
        && !bcNeighbor[OT::surface0N()] && !bcNeighbor[OT::surface2N()])
    {
        edges.push_back(
            DirectedEdge(Box3D(box.x0, box.x0, box.y0 + 1, box.y1 - 1, box.z0, box.z0), 0, 0, -1));
        edges.push_back(
            DirectedEdge(Box3D(box.x0, box.x0, box.y0 + 1, box.y1 - 1, box.z0, box.z0), 2, -1, 0));
    }

    if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
        && (!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
        && !bcNeighbor[OT::surface0N()] && !bcNeighbor[OT::surface2P()])
    {
        edges.push_back(
            DirectedEdge(Box3D(box.x0, box.x0, box.y0 + 1, box.y1 - 1, box.z1, box.z1), 0, 0, +1));
        edges.push_back(
            DirectedEdge(Box3D(box.x0, box.x0, box.y0 + 1, box.y1 - 1, box.z1, box.z1), 2, -1, 0));
    }

    if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
        && (!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
        && !bcNeighbor[OT::surface0P()] && !bcNeighbor[OT::surface2N()])
    {
        edges.push_back(
            DirectedEdge(Box3D(box.x1, box.x1, box.y0 + 1, box.y1 - 1, box.z0, box.z0), 0, 0, -1));
        edges.push_back(
            DirectedEdge(Box3D(box.x1, box.x1, box.y0 + 1, box.y1 - 1, box.z0, box.z0), 2, +1, 0));
    }

    if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
        && (!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
        && !bcNeighbor[OT::surface0P()] && !bcNeighbor[OT::surface2P()])
    {
        edges.push_back(
            DirectedEdge(Box3D(box.x1, box.x1, box.y0 + 1, box.y1 - 1, box.z1, box.z1), 0, 0, +1));
        edges.push_back(
            DirectedEdge(Box3D(box.x1, box.x1, box.y0 + 1, box.y1 - 1, box.z1, box.z1), 2, +1, 0));
    }

    // x edges (must be 8 of them)
    if ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
        && (!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
        && !bcNeighbor[OT::surface1N()] && !bcNeighbor[OT::surface2N()])
    {
        edges.push_back(
            DirectedEdge(Box3D(box.x0 + 1, box.x1 - 1, box.y0, box.y0, box.z0, box.z0), 1, -1, 0));
        edges.push_back(
            DirectedEdge(Box3D(box.x0 + 1, box.x1 - 1, box.y0, box.y0, box.z0, box.z0), 2, 0, -1));
    }

    if ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
        && (!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
        && !bcNeighbor[OT::surface1N()] && !bcNeighbor[OT::surface2P()])
    {
        edges.push_back(
            DirectedEdge(Box3D(box.x0 + 1, box.x1 - 1, box.y0, box.y0, box.z1, box.z1), 1, +1, 0));
        edges.push_back(
            DirectedEdge(Box3D(box.x0 + 1, box.x1 - 1, box.y0, box.y0, box.z1, box.z1), 2, 0, -1));
    }

    if ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
        && (!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
        && !bcNeighbor[OT::surface1P()] && !bcNeighbor[OT::surface2N()])
    {
        edges.push_back(
            DirectedEdge(Box3D(box.x0 + 1, box.x1 - 1, box.y1, box.y1, box.z0, box.z0), 1, -1, 0));
        edges.push_back(
            DirectedEdge(Box3D(box.x0 + 1, box.x1 - 1, box.y1, box.y1, box.z0, box.z0), 2, 0, +1));
    }

    if ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
        && (!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
        && !bcNeighbor[OT::surface1P()] && !bcNeighbor[OT::surface2P()])
    {
        edges.push_back(
            DirectedEdge(Box3D(box.x0 + 1, box.x1 - 1, box.y1, box.y1, box.z1, box.z1), 1, +1, 0));
        edges.push_back(
            DirectedEdge(Box3D(box.x0 + 1, box.x1 - 1, box.y1, box.y1, box.z1, box.z1), 2, 0, +1));
    }

    // ====================== Internal edges ====================== //
    // z edges along internal "corners" of the refinement (must be at most 8 of them)
    if (neighbor[OT::surface0N()] && neighbor[OT::surface1N()]
        && (!neighbor[OT::edge2NN()] && allocNeighbor[OT::edge2NN()]))
    {
        plint dz0 = neighbor[OT::surface2N()]
                    && ((!neighbor[OT::edge0NN()] && allocNeighbor[OT::edge0NN()])
                        || (!neighbor[OT::edge1NN()]
                            && allocNeighbor[OT::edge1NN()])) /* && !((!neighbor[OT::edge0NN()] &&
                                                                 allocNeighbor[OT::edge0NN()]) &&
                                                                 (!neighbor[OT::edge1NN()] &&
                                                                 allocNeighbor[OT::edge1NN()]))*/
            ;
        plint dz1 = neighbor[OT::surface2P()]
                    && ((!neighbor[OT::edge0NP()] && allocNeighbor[OT::edge0NP()])
                        || (!neighbor[OT::edge1PN()]
                            && allocNeighbor[OT::edge1PN()])) /* && !((!neighbor[OT::edge0NP()] &&
                                                                 allocNeighbor[OT::edge0NP()]) &&
                                                                 (!neighbor[OT::edge1PN()] &&
                                                                 allocNeighbor[OT::edge1PN()]))*/
            ;

        edges.push_back(DirectedEdge(
            Box3D(box.x0, box.x0, box.y0, box.y0, box.z0 + 1 + dz0, box.z1 - 1 - dz1), 0, +1, 0));
        edges.push_back(DirectedEdge(
            Box3D(box.x0, box.x0, box.y0, box.y0, box.z0 + 1 + dz0, box.z1 - 1 - dz1), 1, 0, +1));
    }
    if (neighbor[OT::surface0N()] && neighbor[OT::surface1P()]
        && (!neighbor[OT::edge2NP()] && allocNeighbor[OT::edge2NP()]))
    {
        plint dz0 = neighbor[OT::surface2N()]
                    && ((!neighbor[OT::edge0PN()] && allocNeighbor[OT::edge0PN()])
                        || (!neighbor[OT::edge1NN()]
                            && allocNeighbor[OT::edge1NN()])) /* && !((!neighbor[OT::edge0PN()] &&
                                                                 allocNeighbor[OT::edge0PN()]) &&
                                                                 (!neighbor[OT::edge1NN()] &&
                                                                 allocNeighbor[OT::edge1NN()]))*/
            ;
        plint dz1 = neighbor[OT::surface2P()]
                    && ((!neighbor[OT::edge0PP()] && allocNeighbor[OT::edge0PP()])
                        || (!neighbor[OT::edge1PN()]
                            && allocNeighbor[OT::edge1PN()])) /* && !((!neighbor[OT::edge0PP()] &&
                                                                 allocNeighbor[OT::edge0PP()]) &&
                                                                 (!neighbor[OT::edge1PN()] &&
                                                                 allocNeighbor[OT::edge1PN()]))*/
            ;

        edges.push_back(DirectedEdge(
            Box3D(box.x0, box.x0, box.y1, box.y1, box.z0 + 1 + dz0, box.z1 - 1 - dz1), 0, -1, 0));
        edges.push_back(DirectedEdge(
            Box3D(box.x0, box.x0, box.y1, box.y1, box.z0 + 1 + dz0, box.z1 - 1 - dz1), 1, 0, +1));
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1N()]
        && (!neighbor[OT::edge2PN()] && allocNeighbor[OT::edge2PN()]))
    {
        plint dz0 = neighbor[OT::surface2N()]
                    && ((!neighbor[OT::edge0NN()] && allocNeighbor[OT::edge0NN()])
                        || (!neighbor[OT::edge1NP()]
                            && allocNeighbor[OT::edge1NP()])) /* && !((!neighbor[OT::edge0NN()] &&
                                                                 allocNeighbor[OT::edge0NN()]) &&
                                                                 (!neighbor[OT::edge1NP()] &&
                                                                 allocNeighbor[OT::edge1NP()]))*/
            ;
        plint dz1 = neighbor[OT::surface2P()]
                    && ((!neighbor[OT::edge0NP()] && allocNeighbor[OT::edge0NP()])
                        || (!neighbor[OT::edge1PP()]
                            && allocNeighbor[OT::edge1PP()])) /* && !((!neighbor[OT::edge0NP()] &&
                                                                 allocNeighbor[OT::edge0NP()]) &&
                                                                 (!neighbor[OT::edge1PP()] &&
                                                                 allocNeighbor[OT::edge1PP()]))*/
            ;

        edges.push_back(DirectedEdge(
            Box3D(box.x1, box.x1, box.y0, box.y0, box.z0 + 1 + dz0, box.z1 - 1 - dz1), 0, +1, 0));
        edges.push_back(DirectedEdge(
            Box3D(box.x1, box.x1, box.y0, box.y0, box.z0 + 1 + dz0, box.z1 - 1 - dz1), 1, 0, -1));
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1P()]
        && (!neighbor[OT::edge2PP()] && allocNeighbor[OT::edge2PP()]))
    {
        plint dz0 = neighbor[OT::surface2N()]
                    && ((!neighbor[OT::edge0PN()] && allocNeighbor[OT::edge0PN()])
                        || (!neighbor[OT::edge1NP()]
                            && allocNeighbor[OT::edge1NP()])) /* && !((!neighbor[OT::edge0PN()] &&
                                                                 allocNeighbor[OT::edge0PN()]) &&
                                                                 (!neighbor[OT::edge1NP()] &&
                                                                 allocNeighbor[OT::edge1NP()]))*/
            ;
        plint dz1 = neighbor[OT::surface2P()]
                    && ((!neighbor[OT::edge0PP()] && allocNeighbor[OT::edge0PP()])
                        || (!neighbor[OT::edge1PP()]
                            && allocNeighbor[OT::edge1PP()])) /* && !((!neighbor[OT::edge0PP()] &&
                                                                 allocNeighbor[OT::edge0PP()]) &&
                                                                 (!neighbor[OT::edge1PP()] &&
                                                                 allocNeighbor[OT::edge1PP()]))*/
            ;

        edges.push_back(DirectedEdge(
            Box3D(box.x1, box.x1, box.y1, box.y1, box.z0 + 1 + dz0, box.z1 - 1 - dz1), 0, -1, 0));
        edges.push_back(DirectedEdge(
            Box3D(box.x1, box.x1, box.y1, box.y1, box.z0 + 1 + dz0, box.z1 - 1 - dz1), 1, 0, -1));
    }

    // y edges (must be 8 of them)
    if (neighbor[OT::surface0N()] && neighbor[OT::surface2N()]
        && (!neighbor[OT::edge1NN()] && allocNeighbor[OT::edge1NN()]))
    {
        plint dy0 = neighbor[OT::surface1N()]
                    && ((!neighbor[OT::edge0NN()] && allocNeighbor[OT::edge0NN()])
                        || (!neighbor[OT::edge2NN()]
                            && allocNeighbor[OT::edge2NN()])) /* && !((!neighbor[OT::edge0NN()] &&
                                                                 allocNeighbor[OT::edge0NN()]) &&
                                                                 (!neighbor[OT::edge2NN()] &&
                                                                 allocNeighbor[OT::edge2NN()]))*/
            ;
        plint dy1 = neighbor[OT::surface1P()]
                    && ((!neighbor[OT::edge0PN()] && allocNeighbor[OT::edge0PN()])
                        || (!neighbor[OT::edge2NP()]
                            && allocNeighbor[OT::edge2NP()])) /* && !((!neighbor[OT::edge0PN()] &&
                                                                 allocNeighbor[OT::edge0PN()]) &&
                                                                 (!neighbor[OT::edge2NP()] &&
                                                                 allocNeighbor[OT::edge2NP()]))*/
            ;

        edges.push_back(DirectedEdge(
            Box3D(box.x0, box.x0, box.y0 + 1 + dy0, box.y1 - 1 - dy1, box.z0, box.z0), 0, 0, +1));
        edges.push_back(DirectedEdge(
            Box3D(box.x0, box.x0, box.y0 + 1 + dy0, box.y1 - 1 - dy1, box.z0, box.z0), 2, +1, 0));
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface2P()]
        && (!neighbor[OT::edge1PN()] && allocNeighbor[OT::edge1PN()]))
    {
        plint dy0 = neighbor[OT::surface1N()]
                    && ((!neighbor[OT::edge0NP()] && allocNeighbor[OT::edge0NP()])
                        || (!neighbor[OT::edge2NN()]
                            && allocNeighbor[OT::edge2NN()])) /* && !((!neighbor[OT::edge0NP()] &&
                                                                 allocNeighbor[OT::edge0NP()]) &&
                                                                 (!neighbor[OT::edge2NN()] &&
                                                                 allocNeighbor[OT::edge2NN()]))*/
            ;
        plint dy1 = neighbor[OT::surface1P()]
                    && ((!neighbor[OT::edge0PP()] && allocNeighbor[OT::edge0PP()])
                        || (!neighbor[OT::edge2NP()]
                            && allocNeighbor[OT::edge2NP()])) /* && !((!neighbor[OT::edge0PP()] &&
                                                                 allocNeighbor[OT::edge0PP()]) &&
                                                                 (!neighbor[OT::edge2NP()] &&
                                                                 allocNeighbor[OT::edge2NP()]))*/
            ;

        edges.push_back(DirectedEdge(
            Box3D(box.x0, box.x0, box.y0 + 1 + dy0, box.y1 - 1 - dy1, box.z1, box.z1), 0, 0, -1));
        edges.push_back(DirectedEdge(
            Box3D(box.x0, box.x0, box.y0 + 1 + dy0, box.y1 - 1 - dy1, box.z1, box.z1), 2, +1, 0));
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface2N()]
        && (!neighbor[OT::edge1NP()] && allocNeighbor[OT::edge1NP()]))
    {
        plint dy0 = neighbor[OT::surface1N()]
                    && ((!neighbor[OT::edge0NN()] && allocNeighbor[OT::edge0NN()])
                        || (!neighbor[OT::edge2PN()]
                            && allocNeighbor[OT::edge2PN()])) /* && !((!neighbor[OT::edge0NN()] &&
                                                                 allocNeighbor[OT::edge0NN()]) &&
                                                                 (!neighbor[OT::edge2PN()] &&
                                                                 allocNeighbor[OT::edge2PN()]))*/
            ;
        plint dy1 = neighbor[OT::surface1P()]
                    && ((!neighbor[OT::edge0PN()] && allocNeighbor[OT::edge0PN()])
                        || (!neighbor[OT::edge2PP()]
                            && allocNeighbor[OT::edge2PP()])) /* && !((!neighbor[OT::edge0PN()] &&
                                                                 allocNeighbor[OT::edge0PN()]) &&
                                                                 (!neighbor[OT::edge2PP()] &&
                                                                 allocNeighbor[OT::edge2PP()]))*/
            ;

        edges.push_back(DirectedEdge(
            Box3D(box.x1, box.x1, box.y0 + 1 + dy0, box.y1 - 1 - dy1, box.z0, box.z0), 0, 0, +1));
        edges.push_back(DirectedEdge(
            Box3D(box.x1, box.x1, box.y0 + 1 + dy0, box.y1 - 1 - dy1, box.z0, box.z0), 2, -1, 0));
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface2P()]
        && (!neighbor[OT::edge1PP()] && allocNeighbor[OT::edge1PP()]))
    {
        plint dy0 = neighbor[OT::surface1N()]
                    && ((!neighbor[OT::edge0NP()] && allocNeighbor[OT::edge0NP()])
                        || (!neighbor[OT::edge2PN()]
                            && allocNeighbor[OT::edge2PN()])) /* && !((!neighbor[OT::edge0NP()] &&
                                                                 allocNeighbor[OT::edge0NP()]) &&
                                                                 (!neighbor[OT::edge2PN()] &&
                                                                 allocNeighbor[OT::edge2PN()]))*/
            ;
        plint dy1 = neighbor[OT::surface1P()]
                    && ((!neighbor[OT::edge0PP()] && allocNeighbor[OT::edge0PP()])
                        || (!neighbor[OT::edge2PP()]
                            && allocNeighbor[OT::edge2PP()])) /* && !((!neighbor[OT::edge0PP()] &&
                                                                 allocNeighbor[OT::edge0PP()]) &&
                                                                 (!neighbor[OT::edge2PP()] &&
                                                                 allocNeighbor[OT::edge2PP()]))*/
            ;

        edges.push_back(DirectedEdge(
            Box3D(box.x1, box.x1, box.y0 + 1 + dy0, box.y1 - 1 - dy1, box.z1, box.z1), 0, 0, -1));
        edges.push_back(DirectedEdge(
            Box3D(box.x1, box.x1, box.y0 + 1 + dy0, box.y1 - 1 - dy1, box.z1, box.z1), 2, -1, 0));
    }

    // x edges (must be 8 of them)
    if (neighbor[OT::surface1N()] && neighbor[OT::surface2N()]
        && (!neighbor[OT::edge0NN()] && allocNeighbor[OT::edge0NN()]))
    {
        plint dx0 = neighbor[OT::surface0N()]
                    && ((!neighbor[OT::edge1NN()] && allocNeighbor[OT::edge1NN()])
                        || (!neighbor[OT::edge2NN()]
                            && allocNeighbor[OT::edge2NN()])) /* && !((!neighbor[OT::edge1NN()] &&
                                                                 allocNeighbor[OT::edge1NN()]) &&
                                                                 (!neighbor[OT::edge2NN()] &&
                                                                 allocNeighbor[OT::edge2NN()]))*/
            ;
        plint dx1 = neighbor[OT::surface0P()]
                    && ((!neighbor[OT::edge1NP()] && allocNeighbor[OT::edge1NP()])
                        || (!neighbor[OT::edge2PN()]
                            && allocNeighbor[OT::edge2PN()])) /* && !((!neighbor[OT::edge1NP()] &&
                                                                 allocNeighbor[OT::edge1NP()]) &&
                                                                 (!neighbor[OT::edge2PN()] &&
                                                                 allocNeighbor[OT::edge2PN()]))*/
            ;

        edges.push_back(DirectedEdge(
            Box3D(box.x0 + 1 + dx0, box.x1 - 1 - dx1, box.y0, box.y0, box.z0, box.z0), 1, +1, 0));
        edges.push_back(DirectedEdge(
            Box3D(box.x0 + 1 + dx0, box.x1 - 1 - dx1, box.y0, box.y0, box.z0, box.z0), 2, 0, +1));
    }

    if (neighbor[OT::surface1N()] && neighbor[OT::surface2P()]
        && (!neighbor[OT::edge0NP()] && allocNeighbor[OT::edge0NP()]))
    {
        plint dx0 = neighbor[OT::surface0N()]
                    && ((!neighbor[OT::edge1PN()] && allocNeighbor[OT::edge1PN()])
                        || (!neighbor[OT::edge2NN()]
                            && allocNeighbor[OT::edge2NN()])) /* && !((!neighbor[OT::edge1PN()] &&
                                                                 allocNeighbor[OT::edge1PN()]) &&
                                                                 (!neighbor[OT::edge2NN()] &&
                                                                 allocNeighbor[OT::edge2NN()]))*/
            ;
        plint dx1 = neighbor[OT::surface0P()]
                    && ((!neighbor[OT::edge1PP()] && allocNeighbor[OT::edge1PP()])
                        || (!neighbor[OT::edge2PN()]
                            && allocNeighbor[OT::edge2PN()])) /* && !((!neighbor[OT::edge1PP()] &&
                                                                 allocNeighbor[OT::edge1PP()]) &&
                                                                 (!neighbor[OT::edge2PN()] &&
                                                                 allocNeighbor[OT::edge2PN()]))*/
            ;

        edges.push_back(DirectedEdge(
            Box3D(box.x0 + 1 + dx0, box.x1 - 1 - dx1, box.y0, box.y0, box.z1, box.z1), 1, -1, 0));
        edges.push_back(DirectedEdge(
            Box3D(box.x0 + 1 + dx0, box.x1 - 1 - dx1, box.y0, box.y0, box.z1, box.z1), 2, 0, +1));
    }

    if (neighbor[OT::surface1P()] && neighbor[OT::surface2N()]
        && (!neighbor[OT::edge0PN()] && allocNeighbor[OT::edge0PN()]))
    {
        plint dx0 = neighbor[OT::surface0N()]
                    && ((!neighbor[OT::edge1NN()] && allocNeighbor[OT::edge1NN()])
                        || (!neighbor[OT::edge2NP()]
                            && allocNeighbor[OT::edge2NP()])) /* && !((!neighbor[OT::edge1NN()] &&
                                                                 allocNeighbor[OT::edge1NN()]) &&
                                                                 (!neighbor[OT::edge2NP()] &&
                                                                 allocNeighbor[OT::edge2NP()]))*/
            ;
        plint dx1 = neighbor[OT::surface0P()]
                    && ((!neighbor[OT::edge1NP()] && allocNeighbor[OT::edge1NP()])
                        || (!neighbor[OT::edge2PP()]
                            && allocNeighbor[OT::edge2PP()])) /* && !((!neighbor[OT::edge1NP()] &&
                                                                 allocNeighbor[OT::edge1NP()]) &&
                                                                 (!neighbor[OT::edge2PP()] &&
                                                                 allocNeighbor[OT::edge2PP()]))*/
            ;

        edges.push_back(DirectedEdge(
            Box3D(box.x0 + 1 + dx0, box.x1 - 1 - dx1, box.y1, box.y1, box.z0, box.z0), 1, +1, 0));
        edges.push_back(DirectedEdge(
            Box3D(box.x0 + 1 + dx0, box.x1 - 1 - dx1, box.y1, box.y1, box.z0, box.z0), 2, 0, -1));
    }

    if (neighbor[OT::surface1P()] && neighbor[OT::surface2P()]
        && (!neighbor[OT::edge0PP()] && allocNeighbor[OT::edge0PP()]))
    {
        plint dx0 = neighbor[OT::surface0N()]
                    && ((!neighbor[OT::edge1PN()] && allocNeighbor[OT::edge1PN()])
                        || (!neighbor[OT::edge2NP()]
                            && allocNeighbor[OT::edge2NP()])) /* && !((!neighbor[OT::edge1PN()] &&
                                                                 allocNeighbor[OT::edge1PN()]) &&
                                                                 (!neighbor[OT::edge2NP()] &&
                                                                 allocNeighbor[OT::edge2NP()]))*/
            ;
        plint dx1 = neighbor[OT::surface0P()]
                    && ((!neighbor[OT::edge1PP()] && allocNeighbor[OT::edge1PP()])
                        || (!neighbor[OT::edge2PP()]
                            && allocNeighbor[OT::edge2PP()])) /* && !((!neighbor[OT::edge1PP()] &&
                                                                 allocNeighbor[OT::edge1PP()]) &&
                                                                 (!neighbor[OT::edge2PP()] &&
                                                                 allocNeighbor[OT::edge2PP()]))*/
            ;

        edges.push_back(DirectedEdge(
            Box3D(box.x0 + 1 + dx0, box.x1 - 1 - dx1, box.y1, box.y1, box.z1, box.z1), 1, -1, 0));
        edges.push_back(DirectedEdge(
            Box3D(box.x0 + 1 + dx0, box.x1 - 1 - dx1, box.y1, box.y1, box.z1, box.z1), 2, 0, -1));
    }

    PLB_ASSERT(containsNoDuplicates(edges, false) && "Dirdected edges must be unique.");
    // PLB_ASSERT(edges.size() == 24 && "There are 24 directed edges in a box.");

    return edges;
}

std::vector<Edge> getBoundaryEdges(
    const Box3D &box, const Array<bool, 26> &neighbor, const Array<bool, 26> &bcNeighbor,
    const Array<bool, 26> &allocNeighbor)
{
    std::vector<Edge> edges;

    if (bcNeighbor[OT::surface0N()]) {
        {
            plint dz0 = 0;
            plint dz1 = 0;
            if ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
                && !bcNeighbor[OT::surface1N()]) {
                dz1 = 2 * neighbor[OT::edge0NP()];
                dz0 = 2 * neighbor[OT::edge0NN()];
                edges.push_back(
                    Edge(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0 - dz0, box.z1 + dz1)));
            }
            if ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
                && !bcNeighbor[OT::surface1P()]) {
                dz1 = 2 * neighbor[OT::edge0PP()];
                dz0 = 2 * neighbor[OT::edge0PN()];
                edges.push_back(
                    Edge(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0 - dz0, box.z1 + dz1)));
            }
        }

        {
            plint dy0 = 0;
            plint dy1 = 0;
            if ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
                && !bcNeighbor[OT::surface2N()]) {
                dy1 = 2 * neighbor[OT::edge0PN()];
                dy0 = 2 * neighbor[OT::edge0NN()];
                edges.push_back(
                    Edge(Box3D(box.x0, box.x0, box.y0 - dy0, box.y1 + dy1, box.z0, box.z0)));
            }
            if ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
                && !bcNeighbor[OT::surface2P()]) {
                dy1 = 2 * neighbor[OT::edge0PP()];
                dy0 = 2 * neighbor[OT::edge0NP()];
                edges.push_back(
                    Edge(Box3D(box.x0, box.x0, box.y0 - dy0, box.y1 + dy1, box.z1, box.z1)));
            }
        }
    }

    if (bcNeighbor[OT::surface0P()]) {
        {
            plint dz0 = 0;
            plint dz1 = 0;
            if ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
                && !bcNeighbor[OT::surface1N()]) {
                dz1 = 2 * neighbor[OT::edge0NP()];
                dz0 = 2 * neighbor[OT::edge0NN()];
                edges.push_back(
                    Edge(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0 - dz0, box.z1 + dz1)));
            }

            if ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
                && !bcNeighbor[OT::surface1P()]) {
                dz1 = 2 * neighbor[OT::edge0PP()];
                dz0 = 2 * neighbor[OT::edge0PN()];
                edges.push_back(
                    Edge(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0 - dz0, box.z1 + dz1)));
            }
        }

        {
            plint dy0 = 0;
            plint dy1 = 0;
            if ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
                && !bcNeighbor[OT::surface2N()]) {
                dy1 = 2 * neighbor[OT::edge0PN()];
                dy0 = 2 * neighbor[OT::edge0NN()];

                edges.push_back(
                    Edge(Box3D(box.x1, box.x1, box.y0 - dy0, box.y1 + dy1, box.z0, box.z0)));
            }
            if ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
                && !bcNeighbor[OT::surface2P()]) {
                dy1 = 2 * neighbor[OT::edge0PP()];
                dy0 = 2 * neighbor[OT::edge0NP()];
                edges.push_back(
                    Edge(Box3D(box.x1, box.x1, box.y0 - dy0, box.y1 + dy1, box.z1, box.z1)));
            }
        }
    }

    if (bcNeighbor[OT::surface1N()]) {
        {
            plint dz0 = 0;
            plint dz1 = 0;
            if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
                && !bcNeighbor[OT::surface0N()]) {
                dz1 = 2 * neighbor[OT::edge1PN()];
                dz0 = 2 * neighbor[OT::edge1NN()];
                edges.push_back(
                    Edge(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0 - dz0, box.z1 + dz1)));
            }
            if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
                && !bcNeighbor[OT::surface0P()]) {
                dz1 = 2 * neighbor[OT::edge1PP()];
                dz0 = 2 * neighbor[OT::edge1NP()];
                edges.push_back(
                    Edge(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0 - dz0, box.z1 + dz1)));
            }
        }

        {
            plint dx0 = 0;
            plint dx1 = 0;
            if ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
                && !bcNeighbor[OT::surface2N()]) {
                dx1 = 2 * neighbor[OT::edge1NP()];
                dx0 = 2 * neighbor[OT::edge1NN()];
                edges.push_back(
                    Edge(Box3D(box.x0 - dx0, box.x1 + dx1, box.y0, box.y0, box.z0, box.z0)));
            }
            if ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
                && !bcNeighbor[OT::surface2P()]) {
                dx1 = 2 * neighbor[OT::edge1PP()];
                dx0 = 2 * neighbor[OT::edge1PN()];
                edges.push_back(
                    Edge(Box3D(box.x0 - dx0, box.x1 + dx1, box.y0, box.y0, box.z1, box.z1)));
            }
        }
    }

    if (bcNeighbor[OT::surface1P()]) {
        {
            plint dz0 = 0;
            plint dz1 = 0;
            if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
                && !bcNeighbor[OT::surface0N()]) {
                dz1 = 2 * neighbor[OT::edge1PN()];
                dz0 = 2 * neighbor[OT::edge1NN()];
                edges.push_back(
                    Edge(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0 - dz0, box.z1 + dz1)));
            }
            if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
                && !bcNeighbor[OT::surface0P()]) {
                dz1 = 2 * neighbor[OT::edge1PP()];
                dz0 = 2 * neighbor[OT::edge1NP()];
                edges.push_back(
                    Edge(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0 - dz0, box.z1 + dz1)));
            }
        }

        {
            plint dx0 = 0;
            plint dx1 = 0;
            if ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
                && !bcNeighbor[OT::surface2N()]) {
                dx1 = 2 * neighbor[OT::edge1NP()];
                dx0 = 2 * neighbor[OT::edge1NN()];
                edges.push_back(
                    Edge(Box3D(box.x0 - dx0, box.x1 + dx1, box.y1, box.y1, box.z0, box.z0)));
            }
            if ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
                && !bcNeighbor[OT::surface2P()]) {
                dx1 = 2 * neighbor[OT::edge1PP()];
                dx0 = 2 * neighbor[OT::edge1PN()];
                edges.push_back(
                    Edge(Box3D(box.x0 - dx0, box.x1 + dx1, box.y1, box.y1, box.z1, box.z1)));
            }
        }
    }

    if (bcNeighbor[OT::surface2N()]) {
        {
            plint dy0 = 0;
            plint dy1 = 0;
            if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
                && !bcNeighbor[OT::surface0N()]) {
                dy1 = 2 * neighbor[OT::edge2NP()];
                dy0 = 2 * neighbor[OT::edge2NN()];
                edges.push_back(
                    Edge(Box3D(box.x0, box.x0, box.y0 - dy0, box.y1 + dy1, box.z0, box.z0)));
            }
            if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
                && !bcNeighbor[OT::surface0P()]) {
                dy1 = 2 * neighbor[OT::edge2PP()];
                dy0 = 2 * neighbor[OT::edge2PN()];
                edges.push_back(
                    Edge(Box3D(box.x1, box.x1, box.y0 - dy0, box.y1 + dy1, box.z0, box.z0)));
            }
        }
        {
            plint dx0 = 0;
            plint dx1 = 0;
            if ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
                && !bcNeighbor[OT::surface1N()]) {
                dx1 = 2 * neighbor[OT::edge2PN()];
                dx0 = 2 * neighbor[OT::edge2NN()];
                edges.push_back(
                    Edge(Box3D(box.x0 - dx0, box.x1 + dx1, box.y0, box.y0, box.z0, box.z0)));
            }
            if ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
                && !bcNeighbor[OT::surface1P()]) {
                dx1 = 2 * neighbor[OT::edge2PP()];
                dx0 = 2 * neighbor[OT::edge2NP()];
                edges.push_back(
                    Edge(Box3D(box.x0 - dx0, box.x1 + dx1, box.y1, box.y1, box.z0, box.z0)));
            }
        }
    }

    if (bcNeighbor[OT::surface2P()]) {
        {
            plint dy0 = 0;
            plint dy1 = 0;
            if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
                && !bcNeighbor[OT::surface0N()]) {
                dy1 = 2 * neighbor[OT::edge2NP()];
                dy0 = 2 * neighbor[OT::edge2NN()];
                edges.push_back(
                    Edge(Box3D(box.x0, box.x0, box.y0 - dy0, box.y1 + dy1, box.z1, box.z1)));
            }
            if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
                && !bcNeighbor[OT::surface0P()]) {
                dy1 = 2 * neighbor[OT::edge2PP()];
                dy0 = 2 * neighbor[OT::edge2PN()];
                edges.push_back(
                    Edge(Box3D(box.x1, box.x1, box.y0 - dy0, box.y1 + dy1, box.z1, box.z1)));
            }
        }
        {
            plint dx0 = 0;
            plint dx1 = 0;
            if ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
                && !bcNeighbor[OT::surface1N()]) {
                dx1 = 2 * neighbor[OT::edge2PN()];
                dx0 = 2 * neighbor[OT::edge2NN()];
                edges.push_back(
                    Edge(Box3D(box.x0 - dx0, box.x1 + dx1, box.y0, box.y0, box.z1, box.z1)));
            }
            if ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
                && !bcNeighbor[OT::surface1P()]) {
                dx1 = 2 * neighbor[OT::edge2PP()];
                dx0 = 2 * neighbor[OT::edge2NP()];
                edges.push_back(
                    Edge(Box3D(box.x0 - dx0, box.x1 + dx1, box.y1, box.y1, box.z1, box.z1)));
            }
        }
    }

    PLB_ASSERT(containsNoDuplicates(edges, true) && "Edges must be unique.");

    return edges;
}

std::vector<DirectedEdge> getBoundaryDirectedEdges(
    const Box3D &box, const Array<bool, 26> &neighbor, const Array<bool, 26> &bcNeighbor,
    const Array<bool, 26> &allocNeighbor)
{
    std::vector<DirectedEdge> edges;

    if (bcNeighbor[OT::surface0N()]) {
        if ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
            && !bcNeighbor[OT::surface1N()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x0, box.x0, box.y0, box.y0,
                    box.z0 + ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])),
                    box.z1 - ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()]))),
                1, 0, -1));
        }
        if ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
            && !bcNeighbor[OT::surface1P()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x0, box.x0, box.y1, box.y1,
                    box.z0 + ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])),
                    box.z1 - ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()]))),
                1, 0, -1));
        }
        if ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
            && !bcNeighbor[OT::surface2N()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x0, box.x0,
                    box.y0 + ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])),
                    box.y1 - ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])),
                    box.z0, box.z0),
                2, -1, 0));
        }
        if ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
            && !bcNeighbor[OT::surface2P()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x0, box.x0,
                    box.y0 + ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])),
                    box.y1 - ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])),
                    box.z1, box.z1),
                2, -1, 0));
        }
    }

    if (bcNeighbor[OT::surface0P()]) {
        if ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
            && !bcNeighbor[OT::surface1N()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x1, box.x1, box.y0, box.y0,
                    box.z0 + ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])),
                    box.z1 - ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()]))),
                1, 0, +1));
        }
        if ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
            && !bcNeighbor[OT::surface1P()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x1, box.x1, box.y1, box.y1,
                    box.z0 + ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])),
                    box.z1 - ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()]))),
                1, 0, +1));
        }
        if ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
            && !bcNeighbor[OT::surface2N()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x1, box.x1,
                    box.y0 + ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])),
                    box.y1 - ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])),
                    box.z0, box.z0),
                2, +1, 0));
        }
        if ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
            && !bcNeighbor[OT::surface2P()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x1, box.x1,
                    box.y0 + ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])),
                    box.y1 - ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])),
                    box.z1, box.z1),
                2, +1, 0));
        }
    }

    if (bcNeighbor[OT::surface1N()]) {
        if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
            && !bcNeighbor[OT::surface0N()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x0, box.x0, box.y0, box.y0,
                    box.z0 + ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])),
                    box.z1 - ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()]))),
                0, -1, 0));
        }
        if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
            && !bcNeighbor[OT::surface0P()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x1, box.x1, box.y0, box.y0,
                    box.z0 + ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])),
                    box.z1 - ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()]))),
                0, -1, 0));
        }
        if ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
            && !bcNeighbor[OT::surface2N()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x0 + ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])),
                    box.x1 - ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])),
                    box.y0, box.y0, box.z0, box.z0),
                2, 0, -1));
        }
        if ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
            && !bcNeighbor[OT::surface2P()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x0 + ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])),
                    box.x1 - ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])),
                    box.y0, box.y0, box.z1, box.z1),
                2, 0, -1));
        }
    }

    if (bcNeighbor[OT::surface1P()]) {
        if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
            && !bcNeighbor[OT::surface0N()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x0, box.x0, box.y1, box.y1,
                    box.z0 + ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])),
                    box.z1 - ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()]))),
                0, +1, 0));
        }
        if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
            && !bcNeighbor[OT::surface0P()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x1, box.x1, box.y1, box.y1,
                    box.z0 + ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])),
                    box.z1 - ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()]))),
                0, +1, 0));
        }
        if ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
            && !bcNeighbor[OT::surface2N()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x0 + ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])),
                    box.x1 - ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])),
                    box.y1, box.y1, box.z0, box.z0),
                2, 0, +1));
        }
        if ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
            && !bcNeighbor[OT::surface2P()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x0 + ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])),
                    box.x1 - ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])),
                    box.y1, box.y1, box.z1, box.z1),
                2, 0, +1));
        }
    }

    if (bcNeighbor[OT::surface2N()]) {
        if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
            && !bcNeighbor[OT::surface0N()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x0, box.x0,
                    box.y0 + ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])),
                    box.y1 - ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])),
                    box.z0, box.z0),
                0, 0, -1));
        }
        if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
            && !bcNeighbor[OT::surface0P()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x1, box.x1,
                    box.y0 + ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])),
                    box.y1 - ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])),
                    box.z0, box.z0),
                0, 0, -1));
        }
        if ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
            && !bcNeighbor[OT::surface1N()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x0 + ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])),
                    box.x1 - ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])),
                    box.y0, box.y0, box.z0, box.z0),
                1, -1, 0));
        }
        if ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
            && !bcNeighbor[OT::surface1P()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x0 + ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])),
                    box.x1 - ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])),
                    box.y1, box.y1, box.z0, box.z0),
                1, -1, 0));
        }
    }

    if (bcNeighbor[OT::surface2P()]) {
        if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
            && !bcNeighbor[OT::surface0N()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x0, box.x0,
                    box.y0 + ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])),
                    box.y1 - ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])),
                    box.z1, box.z1),
                0, 0, +1));
        }
        if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
            && !bcNeighbor[OT::surface0P()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x1, box.x1,
                    box.y0 + ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])),
                    box.y1 - ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])),
                    box.z1, box.z1),
                0, 0, +1));
        }
        if ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
            && !bcNeighbor[OT::surface1N()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x0 + ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])),
                    box.x1 - ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])),
                    box.y0, box.y0, box.z1, box.z1),
                1, +1, 0));
        }
        if ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
            && !bcNeighbor[OT::surface1P()]) {
            edges.push_back(DirectedEdge(
                Box3D(
                    box.x0 + ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])),
                    box.x1 - ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])),
                    box.y1, box.y1, box.z1, box.z1),
                1, +1, 0));
        }
    }

    PLB_ASSERT(containsNoDuplicates(edges, false) && "Dirdected edges must be unique.");

    return edges;
}

std::vector<DirectedCorner> getDirectedCorners(
    const Box3D &box, const Array<bool, 26> &neighbor, const Array<bool, 26> &bcNeighbor,
    const Array<bool, 26> &allocNeighbor)
{
    std::vector<DirectedCorner> corners;

    if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
        && (!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
        && (!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
        && !bcNeighbor[OT::surface0N()] && !bcNeighbor[OT::surface1N()]
        && !bcNeighbor[OT::surface2N()])
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 0, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 1, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 2, -1, -1));
    }

    if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
        && (!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
        && (!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
        && !bcNeighbor[OT::surface0N()] && !bcNeighbor[OT::surface1N()]
        && !bcNeighbor[OT::surface2P()])
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 0, -1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 1, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 2, -1, -1));
    }

    if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
        && (!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
        && (!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
        && !bcNeighbor[OT::surface0N()] && !bcNeighbor[OT::surface1P()]
        && !bcNeighbor[OT::surface2N()])
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 0, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 1, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 2, -1, +1));
    }

    if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
        && (!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
        && (!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
        && !bcNeighbor[OT::surface0N()] && !bcNeighbor[OT::surface1P()]
        && !bcNeighbor[OT::surface2P()])
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 0, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 1, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 2, -1, +1));
    }

    if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
        && (!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
        && (!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
        && !bcNeighbor[OT::surface0P()] && !bcNeighbor[OT::surface1N()]
        && !bcNeighbor[OT::surface2N()])
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 0, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 1, -1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 2, +1, -1));
    }

    if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
        && (!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
        && (!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
        && !bcNeighbor[OT::surface0P()] && !bcNeighbor[OT::surface1N()]
        && !bcNeighbor[OT::surface2P()])
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 0, -1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 1, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 2, +1, -1));
    }

    if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
        && (!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
        && (!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
        && !bcNeighbor[OT::surface0P()] && !bcNeighbor[OT::surface1P()]
        && !bcNeighbor[OT::surface2N()])
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 0, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 1, -1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 2, +1, +1));
    }

    if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
        && (!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
        && (!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
        && !bcNeighbor[OT::surface0P()] && !bcNeighbor[OT::surface1P()]
        && !bcNeighbor[OT::surface2P()])
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 0, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 1, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 2, +1, +1));
    }

    // Corners must be added on the end of internal edges
    // z edges along internal "corners" of the refinement (must be at most 8 of them)
    if (neighbor[OT::surface0N()] && neighbor[OT::surface1N()]
        && (!neighbor[OT::edge2NN()] && allocNeighbor[OT::edge2NN()])
        && !bcNeighbor[OT::surface0N()] && !bcNeighbor[OT::surface1N()])
    {
        if ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
            && !bcNeighbor[OT::surface2P()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 0, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 1, +1, +1));

            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 2, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 2, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 2, +1, -1));
        }
        if ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
            && !bcNeighbor[OT::surface2N()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 0, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 1, -1, +1));

            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 2, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 2, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 2, +1, -1));
        }
    }
    if (neighbor[OT::surface0N()] && neighbor[OT::surface1P()]
        && (!neighbor[OT::edge2NP()] && allocNeighbor[OT::edge2NP()])
        && !bcNeighbor[OT::surface0N()] && !bcNeighbor[OT::surface1P()])
    {
        if ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
            && !bcNeighbor[OT::surface2P()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 0, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 1, +1, -1));

            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 2, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 2, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 2, -1, -1));
        }
        if ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
            && !bcNeighbor[OT::surface2N()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 0, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 1, -1, +1));

            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 2, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 2, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 2, +1, +1));
        }
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1N()]
        && (!neighbor[OT::edge2PN()] && allocNeighbor[OT::edge2PN()])
        && !bcNeighbor[OT::surface0P()] && !bcNeighbor[OT::surface1N()])
    {
        if ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
            && !bcNeighbor[OT::surface2P()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 0, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 1, +1, -1));

            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 2, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 2, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 2, +1, +1));
        }
        if ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
            && !bcNeighbor[OT::surface2N()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 0, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 1, -1, -1));

            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 2, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 2, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 2, +1, +1));
        }
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1P()]
        && (!neighbor[OT::edge2PP()] && allocNeighbor[OT::edge2PP()])
        && !bcNeighbor[OT::surface0P()] && !bcNeighbor[OT::surface1P()])
    {
        if ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
            && !bcNeighbor[OT::surface2P()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 0, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 1, +1, -1));

            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 2, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 2, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 2, +1, -1));
        }
        if ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
            && !bcNeighbor[OT::surface2N()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 0, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 1, -1, -1));

            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 2, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 2, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 2, +1, -1));
        }
    }

    // y edges (must be 8 of them)
    if (neighbor[OT::surface0N()] && neighbor[OT::surface2N()]
        && (!neighbor[OT::edge1NN()] && allocNeighbor[OT::edge1NN()])
        && !bcNeighbor[OT::surface0N()] && !bcNeighbor[OT::surface2N()])
    {
        if ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
            && !bcNeighbor[OT::surface1P()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 0, +1, +1));

            corners.push_back(DirectedCorner(
                Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 1, -1,
                -1));  // CHANGED RECENTLY (16.2.17) NOT 100% SURE IT'S NEEDED (TODO)
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 1, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 1, -1, +1));

            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 2, +1, +1));
        }
        if ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
            && !bcNeighbor[OT::surface1N()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 0, -1, +1));

            corners.push_back(DirectedCorner(
                Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 1, -1,
                -1));  // CHANGED RECENTLY (16.2.17) NOT 100% SURE IT'S NEEDED (TODO)
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 1, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 1, -1, +1));

            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 2, +1, -1));
        }
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface2P()]
        && (!neighbor[OT::edge1PN()] && allocNeighbor[OT::edge1PN()])
        && !bcNeighbor[OT::surface0N()] && !bcNeighbor[OT::surface2P()])
    {
        if ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
            && !bcNeighbor[OT::surface1P()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 0, +1, -1));

            corners.push_back(DirectedCorner(
                Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 1, +1,
                -1));  // CHANGED RECENTLY (16.2.17) NOT 100% SURE IT'S NEEDED (TODO)
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 1, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 1, -1, -1));

            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 2, +1, +1));
        }
        if ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
            && !bcNeighbor[OT::surface1N()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 0, -1, -1));

            corners.push_back(DirectedCorner(
                Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 1, +1,
                -1));  // CHANGED RECENTLY (16.2.17) NOT 100% SURE IT'S NEEDED (TODO)
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 1, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 1, -1, -1));

            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 2, +1, -1));
        }
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface2N()]
        && (!neighbor[OT::edge1NP()] && allocNeighbor[OT::edge1NP()])
        && !bcNeighbor[OT::surface0P()] && !bcNeighbor[OT::surface2N()])
    {
        if ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
            && !bcNeighbor[OT::surface1P()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 0, +1, +1));

            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 1, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 1, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 1, -1, -1));

            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 2, -1, +1));
        }
        if ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
            && !bcNeighbor[OT::surface1N()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 0, -1, +1));

            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 1, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 1, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 1, -1, -1));

            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 2, -1, -1));
        }
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface2P()]
        && (!neighbor[OT::edge1PP()] && allocNeighbor[OT::edge1PP()])
        && !bcNeighbor[OT::surface0P()] && !bcNeighbor[OT::surface2P()])
    {
        if ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
            && !bcNeighbor[OT::surface1P()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 0, +1, -1));

            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 1, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 1, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 1, -1, +1));

            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 2, -1, +1));
        }
        if ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
            && !bcNeighbor[OT::surface1N()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 0, -1, -1));

            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 1, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 1, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 1, -1, +1));

            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 2, -1, -1));
        }
    }

    // x edges (must be 8 of them)
    if (neighbor[OT::surface1N()] && neighbor[OT::surface2N()]
        && (!neighbor[OT::edge0NN()] && allocNeighbor[OT::edge0NN()])
        && !bcNeighbor[OT::surface1N()] && !bcNeighbor[OT::surface2N()])
    {
        if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
            && !bcNeighbor[OT::surface0P()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 0, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 0, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 0, -1, +1));

            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 1, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 2, +1, +1));
        }
        if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
            && !bcNeighbor[OT::surface0N()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 0, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 0, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 0, -1, +1));

            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 1, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 2, -1, +1));
        }
    }

    if (neighbor[OT::surface1N()] && neighbor[OT::surface2P()]
        && (!neighbor[OT::edge0NP()] && allocNeighbor[OT::edge0NP()])
        && !bcNeighbor[OT::surface1N()] && !bcNeighbor[OT::surface2P()])
    {
        if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
            && !bcNeighbor[OT::surface0P()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 0, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 0, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 0, -1, +1));

            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 1, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 2, +1, +1));
        }
        if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
            && !bcNeighbor[OT::surface0N()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 0, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 0, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 0, -1, +1));

            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 1, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 2, -1, +1));
        }
    }

    if (neighbor[OT::surface1P()] && neighbor[OT::surface2N()]
        && (!neighbor[OT::edge0PN()] && allocNeighbor[OT::edge0PN()])
        && !bcNeighbor[OT::surface1P()] && !bcNeighbor[OT::surface2N()])
    {
        if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
            && !bcNeighbor[OT::surface0P()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 0, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 0, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 0, -1, -1));

            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 1, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 2, +1, -1));
        }
        if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
            && !bcNeighbor[OT::surface0N()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 0, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 0, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 0, -1, -1));

            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 1, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 2, -1, -1));
        }
    }

    if (neighbor[OT::surface1P()] && neighbor[OT::surface2P()]
        && (!neighbor[OT::edge0PP()] && allocNeighbor[OT::edge0PP()])
        && !bcNeighbor[OT::surface1P()] && !bcNeighbor[OT::surface2P()])
    {
        if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
            && !bcNeighbor[OT::surface0P()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 0, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 0, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 0, -1, +1));

            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 1, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 2, +1, -1));
        }
        if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
            && !bcNeighbor[OT::surface0N()]) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 0, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 0, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 0, -1, +1));

            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 1, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 2, -1, -1));
        }
    }

    // Corners must be added on the end of internal corners
    if (neighbor[OT::surface0N()] && neighbor[OT::surface1N()] && neighbor[OT::surface2N()]
        && (!neighbor[OT::cornerNNN()] && allocNeighbor[OT::cornerNNN()]))
    {
        Box3D bb = Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0);

        if (neighbor[OT::edge0NN()] && (!neighbor[OT::edge1NN()] && allocNeighbor[OT::edge1NN()])
            && (!neighbor[OT::edge2NN()] && allocNeighbor[OT::edge2NN()]))
        {
            corners.push_back(DirectedCorner(bb, 0, +1, +1));
            corners.push_back(DirectedCorner(bb, 0, -1, +1));
            corners.push_back(DirectedCorner(bb, 0, +1, -1));

            corners.push_back(DirectedCorner(bb, 1, -1, +1));
            corners.push_back(DirectedCorner(bb, 2, +1, -1));
        }

        if ((!neighbor[OT::edge0NN()] && allocNeighbor[OT::edge0NN()]) && neighbor[OT::edge1NN()]
            && (!neighbor[OT::edge2NN()] && allocNeighbor[OT::edge2NN()]))
        {
            corners.push_back(DirectedCorner(bb, 1, +1, +1));
            corners.push_back(DirectedCorner(bb, 1, -1, +1));
            corners.push_back(DirectedCorner(bb, 1, +1, -1));

            corners.push_back(DirectedCorner(bb, 0, +1, -1));
            corners.push_back(DirectedCorner(bb, 2, -1, +1));
        }

        if ((!neighbor[OT::edge0NN()] && allocNeighbor[OT::edge0NN()])
            && (!neighbor[OT::edge1NN()] && allocNeighbor[OT::edge1NN()])
            && neighbor[OT::edge2NN()])
        {
            corners.push_back(DirectedCorner(bb, 2, -1, +1));
            corners.push_back(DirectedCorner(bb, 2, +1, +1));
            corners.push_back(DirectedCorner(bb, 2, +1, -1));

            corners.push_back(DirectedCorner(bb, 0, -1, +1));
            corners.push_back(DirectedCorner(bb, 1, +1, -1));
        }
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1N()] && neighbor[OT::surface2P()]
        && (!neighbor[OT::cornerNNP()] && allocNeighbor[OT::cornerNNP()]))
    {
        Box3D bb = Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1);

        if (neighbor[OT::edge0NP()] && (!neighbor[OT::edge1PN()] && allocNeighbor[OT::edge1PN()])
            && (!neighbor[OT::edge2NN()] && allocNeighbor[OT::edge2NN()]))
        {
            corners.push_back(DirectedCorner(bb, 0, +1, -1));
            corners.push_back(DirectedCorner(bb, 0, -1, -1));
            corners.push_back(DirectedCorner(bb, 0, +1, +1));

            corners.push_back(DirectedCorner(bb, 1, +1, +1));
            corners.push_back(DirectedCorner(bb, 2, +1, -1));
        }

        if ((!neighbor[OT::edge0NP()] && allocNeighbor[OT::edge0NP()]) && neighbor[OT::edge1PN()]
            && (!neighbor[OT::edge2NN()] && allocNeighbor[OT::edge2NN()]))
        {
            corners.push_back(DirectedCorner(bb, 1, +1, +1));
            corners.push_back(DirectedCorner(bb, 1, -1, +1));
            corners.push_back(DirectedCorner(bb, 1, -1, -1));

            corners.push_back(DirectedCorner(bb, 0, +1, +1));
            corners.push_back(DirectedCorner(bb, 2, -1, +1));
        }

        if ((!neighbor[OT::edge0NP()] && allocNeighbor[OT::edge0NP()])
            && (!neighbor[OT::edge1PN()] && allocNeighbor[OT::edge1PN()])
            && neighbor[OT::edge2NN()])
        {
            corners.push_back(DirectedCorner(bb, 2, -1, +1));
            corners.push_back(DirectedCorner(bb, 2, +1, +1));
            corners.push_back(DirectedCorner(bb, 2, +1, -1));

            corners.push_back(DirectedCorner(bb, 0, -1, -1));
            corners.push_back(DirectedCorner(bb, 1, -1, -1));
        }
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1P()] && neighbor[OT::surface2N()]
        && (!neighbor[OT::cornerNPN()] && allocNeighbor[OT::cornerNPN()]))
    {
        Box3D bb = Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0);

        if (neighbor[OT::edge0PN()] && (!neighbor[OT::edge1NN()] && allocNeighbor[OT::edge1NN()])
            && (!neighbor[OT::edge2NP()] && allocNeighbor[OT::edge2NP()]))
        {
            corners.push_back(DirectedCorner(bb, 0, -1, +1));
            corners.push_back(DirectedCorner(bb, 0, +1, +1));
            corners.push_back(DirectedCorner(bb, 0, -1, -1));

            corners.push_back(DirectedCorner(bb, 1, -1, +1));
            corners.push_back(DirectedCorner(bb, 2, +1, +1));
        }

        if ((!neighbor[OT::edge0PN()] && allocNeighbor[OT::edge0PN()]) && neighbor[OT::edge1NN()]
            && (!neighbor[OT::edge2NP()] && allocNeighbor[OT::edge2NP()]))
        {
            corners.push_back(DirectedCorner(bb, 1, +1, +1));
            corners.push_back(DirectedCorner(bb, 1, -1, +1));
            corners.push_back(DirectedCorner(bb, 1, +1, -1));

            corners.push_back(DirectedCorner(bb, 0, -1, -1));
            corners.push_back(DirectedCorner(bb, 2, -1, -1));
        }

        if ((!neighbor[OT::edge0PN()] && allocNeighbor[OT::edge0PN()])
            && (!neighbor[OT::edge1NN()] && allocNeighbor[OT::edge1NN()])
            && neighbor[OT::edge2NP()])
        {
            corners.push_back(DirectedCorner(bb, 2, +1, -1));
            corners.push_back(DirectedCorner(bb, 2, +1, +1));
            corners.push_back(DirectedCorner(bb, 2, -1, -1));

            corners.push_back(DirectedCorner(bb, 0, +1, +1));
            corners.push_back(DirectedCorner(bb, 1, +1, -1));
        }
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1P()] && neighbor[OT::surface2P()]
        && (!neighbor[OT::cornerNPP()] && allocNeighbor[OT::cornerNPP()]))
    {
        Box3D bb = Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1);

        if (neighbor[OT::edge0PP()] && (!neighbor[OT::edge1PN()] && allocNeighbor[OT::edge1PN()])
            && (!neighbor[OT::edge2NP()] && allocNeighbor[OT::edge2NP()]))
        {
            corners.push_back(DirectedCorner(bb, 0, -1, -1));
            corners.push_back(DirectedCorner(bb, 0, +1, -1));
            corners.push_back(DirectedCorner(bb, 0, -1, +1));

            corners.push_back(DirectedCorner(bb, 1, +1, +1));
            corners.push_back(DirectedCorner(bb, 2, +1, +1));
        }

        if ((!neighbor[OT::edge0PP()] && allocNeighbor[OT::edge0PP()]) && neighbor[OT::edge1PN()]
            && (!neighbor[OT::edge2NP()] && allocNeighbor[OT::edge2NP()]))
        {
            corners.push_back(DirectedCorner(bb, 1, -1, +1));
            corners.push_back(DirectedCorner(bb, 1, -1, -1));
            corners.push_back(DirectedCorner(bb, 1, +1, +1));

            corners.push_back(DirectedCorner(bb, 0, -1, +1));
            corners.push_back(DirectedCorner(bb, 2, -1, -1));
        }

        if ((!neighbor[OT::edge0PP()] && allocNeighbor[OT::edge0PP()])
            && (!neighbor[OT::edge1PN()] && allocNeighbor[OT::edge1PN()])
            && neighbor[OT::edge2NP()])
        {
            corners.push_back(DirectedCorner(bb, 2, +1, -1));
            corners.push_back(DirectedCorner(bb, 2, +1, +1));
            corners.push_back(DirectedCorner(bb, 2, -1, -1));

            corners.push_back(DirectedCorner(bb, 0, +1, -1));
            corners.push_back(DirectedCorner(bb, 1, -1, -1));
        }
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1N()] && neighbor[OT::surface2N()]
        && (!neighbor[OT::cornerPNN()] && allocNeighbor[OT::cornerPNN()]))
    {
        Box3D bb = Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0);
        if (neighbor[OT::edge0NN()] && (!neighbor[OT::edge1NP()] && allocNeighbor[OT::edge1NP()])
            && (!neighbor[OT::edge2PN()] && allocNeighbor[OT::edge2PN()]))
        {
            corners.push_back(DirectedCorner(bb, 0, +1, +1));
            corners.push_back(DirectedCorner(bb, 0, -1, +1));
            corners.push_back(DirectedCorner(bb, 0, +1, -1));

            corners.push_back(DirectedCorner(bb, 1, -1, -1));
            corners.push_back(DirectedCorner(bb, 2, -1, -1));
        }

        if ((!neighbor[OT::edge0NN()] && allocNeighbor[OT::edge0NN()]) && neighbor[OT::edge1NP()]
            && (!neighbor[OT::edge2PN()] && allocNeighbor[OT::edge2PN()]))
        {
            corners.push_back(DirectedCorner(bb, 1, +1, +1));
            corners.push_back(DirectedCorner(bb, 1, -1, -1));
            corners.push_back(DirectedCorner(bb, 1, +1, -1));

            corners.push_back(DirectedCorner(bb, 0, +1, -1));
            corners.push_back(DirectedCorner(bb, 2, +1, +1));
        }

        if ((!neighbor[OT::edge0NN()] && allocNeighbor[OT::edge0NN()])
            && (!neighbor[OT::edge1NP()] && allocNeighbor[OT::edge1NP()])
            && neighbor[OT::edge2PN()])
        {
            corners.push_back(DirectedCorner(bb, 2, -1, +1));
            corners.push_back(DirectedCorner(bb, 2, +1, +1));
            corners.push_back(DirectedCorner(bb, 2, -1, -1));

            corners.push_back(DirectedCorner(bb, 0, -1, +1));
            corners.push_back(DirectedCorner(bb, 1, +1, +1));
        }
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1N()] && neighbor[OT::surface2P()]
        && (!neighbor[OT::cornerPNP()] && allocNeighbor[OT::cornerPNP()]))
    {
        Box3D bb = Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1);
        if (neighbor[OT::edge0NP()] && (!neighbor[OT::edge1PP()] && allocNeighbor[OT::edge1PP()])
            && (!neighbor[OT::edge2PN()] && allocNeighbor[OT::edge2PN()]))
        {
            corners.push_back(DirectedCorner(bb, 0, +1, -1));
            corners.push_back(DirectedCorner(bb, 0, -1, -1));
            corners.push_back(DirectedCorner(bb, 0, +1, +1));

            corners.push_back(DirectedCorner(bb, 1, +1, -1));
            corners.push_back(DirectedCorner(bb, 2, -1, -1));
        }

        if ((!neighbor[OT::edge0NP()] && allocNeighbor[OT::edge0NP()]) && neighbor[OT::edge1PP()]
            && (!neighbor[OT::edge2PN()] && allocNeighbor[OT::edge2PN()]))
        {
            corners.push_back(DirectedCorner(bb, 1, -1, -1));
            corners.push_back(DirectedCorner(bb, 1, -1, +1));
            corners.push_back(DirectedCorner(bb, 1, +1, -1));

            corners.push_back(DirectedCorner(bb, 0, +1, +1));
            corners.push_back(DirectedCorner(bb, 2, +1, +1));
        }

        if ((!neighbor[OT::edge0NP()] && allocNeighbor[OT::edge0NP()])
            && (!neighbor[OT::edge1PP()] && allocNeighbor[OT::edge1PP()])
            && neighbor[OT::edge2PN()])
        {
            corners.push_back(DirectedCorner(bb, 2, -1, +1));
            corners.push_back(DirectedCorner(bb, 2, +1, +1));
            corners.push_back(DirectedCorner(bb, 2, -1, -1));

            corners.push_back(DirectedCorner(bb, 0, -1, -1));
            corners.push_back(DirectedCorner(bb, 1, -1, +1));
        }
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1P()] && neighbor[OT::surface2N()]
        && (!neighbor[OT::cornerPPN()] && allocNeighbor[OT::cornerPPN()]))
    {
        Box3D bb = Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0);
        if (neighbor[OT::edge0PN()] && (!neighbor[OT::edge1NP()] && allocNeighbor[OT::edge1NP()])
            && (!neighbor[OT::edge2PP()] && allocNeighbor[OT::edge2PP()]))
        {
            corners.push_back(DirectedCorner(bb, 0, -1, +1));
            corners.push_back(DirectedCorner(bb, 0, +1, +1));
            corners.push_back(DirectedCorner(bb, 0, -1, -1));

            corners.push_back(DirectedCorner(bb, 1, -1, -1));
            corners.push_back(DirectedCorner(bb, 2, -1, +1));
        }

        if ((!neighbor[OT::edge0PN()] && allocNeighbor[OT::edge0PN()]) && neighbor[OT::edge1NP()]
            && (!neighbor[OT::edge2PP()] && allocNeighbor[OT::edge2PP()]))
        {
            corners.push_back(DirectedCorner(bb, 1, +1, -1));
            corners.push_back(DirectedCorner(bb, 1, -1, -1));
            corners.push_back(DirectedCorner(bb, 1, +1, +1));

            corners.push_back(DirectedCorner(bb, 0, -1, -1));
            corners.push_back(DirectedCorner(bb, 2, +1, -1));
        }

        if ((!neighbor[OT::edge0PN()] && allocNeighbor[OT::edge0PN()])
            && (!neighbor[OT::edge1NP()] && allocNeighbor[OT::edge1NP()])
            && neighbor[OT::edge2PP()])
        {
            corners.push_back(DirectedCorner(bb, 2, -1, -1));
            corners.push_back(DirectedCorner(bb, 2, +1, -1));
            corners.push_back(DirectedCorner(bb, 2, -1, +1));

            corners.push_back(DirectedCorner(bb, 0, +1, +1));
            corners.push_back(DirectedCorner(bb, 1, +1, +1));
        }
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1P()] && neighbor[OT::surface2P()]
        && (!neighbor[OT::cornerPPP()] && allocNeighbor[OT::cornerPPP()]))
    {
        Box3D bb = Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1);
        if (neighbor[OT::edge0PP()] && (!neighbor[OT::edge1PP()] && allocNeighbor[OT::edge1PP()])
            && (!neighbor[OT::edge2PP()] && allocNeighbor[OT::edge2PP()]))
        {
            corners.push_back(DirectedCorner(bb, 0, -1, -1));
            corners.push_back(DirectedCorner(bb, 0, +1, -1));
            corners.push_back(DirectedCorner(bb, 0, -1, +1));

            corners.push_back(DirectedCorner(bb, 1, +1, -1));
            corners.push_back(DirectedCorner(bb, 2, -1, +1));
        }

        if ((!neighbor[OT::edge0PP()] && allocNeighbor[OT::edge0PP()]) && neighbor[OT::edge1PP()]
            && (!neighbor[OT::edge2PP()] && allocNeighbor[OT::edge2PP()]))
        {
            corners.push_back(DirectedCorner(bb, 1, -1, -1));
            corners.push_back(DirectedCorner(bb, 1, +1, -1));
            corners.push_back(DirectedCorner(bb, 1, -1, +1));

            corners.push_back(DirectedCorner(bb, 0, -1, +1));
            corners.push_back(DirectedCorner(bb, 2, +1, -1));
        }

        if ((!neighbor[OT::edge0PP()] && allocNeighbor[OT::edge0PP()])
            && (!neighbor[OT::edge1PP()] && allocNeighbor[OT::edge1PP()])
            && neighbor[OT::edge2PP()])
        {
            corners.push_back(DirectedCorner(bb, 2, -1, -1));
            corners.push_back(DirectedCorner(bb, 2, +1, -1));
            corners.push_back(DirectedCorner(bb, 2, -1, +1));

            corners.push_back(DirectedCorner(bb, 0, +1, -1));
            corners.push_back(DirectedCorner(bb, 1, -1, +1));
        }
    }

    // =================================================================== //
    if (neighbor[OT::surface0P()] && neighbor[OT::surface1P()] && neighbor[OT::surface2P()]
        && neighbor[OT::edge0PP()] && neighbor[OT::edge1PP()] && neighbor[OT::edge2PP()]
        && (!neighbor[OT::cornerPPP()] && allocNeighbor[OT::cornerPPP()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 0, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 1, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 2, -1, -1));
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1P()] && neighbor[OT::surface2N()]
        && neighbor[OT::edge0PN()] && neighbor[OT::edge1NP()] && neighbor[OT::edge2PP()]
        && (!neighbor[OT::cornerPPN()] && allocNeighbor[OT::cornerPPN()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 0, -1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 1, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 2, -1, -1));
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1N()] && neighbor[OT::surface2P()]
        && neighbor[OT::edge0NP()] && neighbor[OT::edge1PP()] && neighbor[OT::edge2PN()]
        && (!neighbor[OT::cornerPNP()] && allocNeighbor[OT::cornerPNP()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 0, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 1, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 2, -1, +1));
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1N()] && neighbor[OT::surface2N()]
        && neighbor[OT::edge0NN()] && neighbor[OT::edge1NP()] && neighbor[OT::edge2PN()]
        && (!neighbor[OT::cornerPNN()] && allocNeighbor[OT::cornerPNN()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 0, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 1, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 2, -1, +1));
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1P()] && neighbor[OT::surface2P()]
        && neighbor[OT::edge0PP()] && neighbor[OT::edge1PN()] && neighbor[OT::edge2NP()]
        && (!neighbor[OT::cornerNPP()] && allocNeighbor[OT::cornerNPP()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 0, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 1, -1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 2, +1, -1));
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1P()] && neighbor[OT::surface2N()]
        && neighbor[OT::edge0PN()] && neighbor[OT::edge1NN()] && neighbor[OT::edge2NP()]
        && (!neighbor[OT::cornerNPN()] && allocNeighbor[OT::cornerNPN()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 0, -1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 1, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 2, +1, -1));
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1N()] && neighbor[OT::surface2P()]
        && neighbor[OT::edge0NP()] && neighbor[OT::edge1PN()] && neighbor[OT::edge2NN()]
        && (!neighbor[OT::cornerNNP()] && allocNeighbor[OT::cornerNNP()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 0, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 1, -1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 2, +1, +1));
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1N()] && neighbor[OT::surface2N()]
        && neighbor[OT::edge0NN()] && neighbor[OT::edge1NN()] && neighbor[OT::edge2NN()]
        && (!neighbor[OT::cornerNNN()] && allocNeighbor[OT::cornerNNN()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 0, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 1, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 2, +1, +1));
    }

    // ==================================================================== //
    if (neighbor[OT::surface0P()] && neighbor[OT::surface1P()] && neighbor[OT::surface2P()]
        && neighbor[OT::cornerPPP()] && (neighbor[OT::edge0PP()] && allocNeighbor[OT::edge0PP()])
        && (!neighbor[OT::edge1PP()] && allocNeighbor[OT::edge1PP()])
        && (!neighbor[OT::edge2PP()] && allocNeighbor[OT::edge2PP()]))
    {
        // pcout << "entering new loop 3" << std::endl;
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 0, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 0, -1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 0, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 1, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 2, -1, +1));
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1P()] && neighbor[OT::surface2N()]
        && neighbor[OT::cornerPPN()] && (neighbor[OT::edge0PN()] && allocNeighbor[OT::edge0PN()])
        && (!neighbor[OT::edge1NP()] && allocNeighbor[OT::edge1NP()])
        && (!neighbor[OT::edge2PP()] && allocNeighbor[OT::edge2PP()]))
    {
        // pcout << "entering new loop 1" << std::endl;
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 0, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 0, -1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 0, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 1, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 2, -1, +1));
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1N()] && neighbor[OT::surface2P()]
        && neighbor[OT::cornerPNP()] && (neighbor[OT::edge0NP()] && allocNeighbor[OT::edge0NP()])
        && (!neighbor[OT::edge1PP()] && allocNeighbor[OT::edge1PP()])
        && (!neighbor[OT::edge2PN()] && allocNeighbor[OT::edge2PN()]))
    {
        // pcout << "entering new loop 2" << std::endl;
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 0, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 0, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 0, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 1, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 2, -1, -1));
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1N()] && neighbor[OT::surface2N()]
        && neighbor[OT::cornerPNN()] && (neighbor[OT::edge0NN()] && allocNeighbor[OT::edge0NN()])
        && (!neighbor[OT::edge1NP()] && allocNeighbor[OT::edge1NP()])
        && (!neighbor[OT::edge2PN()] && allocNeighbor[OT::edge2PN()]))
    {
        // pcout << "entering new loop 4" << std::endl;
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 0, -1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 0, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 0, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 1, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 2, -1, -1));
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1P()] && neighbor[OT::surface2P()]
        && neighbor[OT::cornerNPP()] && (neighbor[OT::edge0PP()] && allocNeighbor[OT::edge0PP()])
        && (!neighbor[OT::edge1PN()] && allocNeighbor[OT::edge1PN()])
        && (!neighbor[OT::edge2NP()] && allocNeighbor[OT::edge2NP()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 0, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 0, -1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 0, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 1, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 2, +1, +1));
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1P()] && neighbor[OT::surface2N()]
        && neighbor[OT::cornerNPN()] && (neighbor[OT::edge0PN()] && allocNeighbor[OT::edge0PN()])
        && (!neighbor[OT::edge1NN()] && allocNeighbor[OT::edge1NN()])
        && (!neighbor[OT::edge2NP()] && allocNeighbor[OT::edge2NP()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 0, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 0, -1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 0, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 1, -1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 2, +1, +1));
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1N()] && neighbor[OT::surface2P()]
        && neighbor[OT::cornerNNP()] && (neighbor[OT::edge0NP()] && allocNeighbor[OT::edge0NP()])
        && (!neighbor[OT::edge1PN()] && allocNeighbor[OT::edge1PN()])
        && (!neighbor[OT::edge2NN()] && allocNeighbor[OT::edge2NN()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 0, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 0, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 0, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 1, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 2, +1, -1));
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1N()] && neighbor[OT::surface2N()]
        && neighbor[OT::cornerNNN()] && (neighbor[OT::edge0NN()] && allocNeighbor[OT::edge0NN()])
        && (!neighbor[OT::edge1NN()] && allocNeighbor[OT::edge1NN()])
        && (!neighbor[OT::edge2NN()] && allocNeighbor[OT::edge2NN()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 0, -1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 0, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 0, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 1, -1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 2, +1, -1));
    }

    // New things added.

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1P()] && neighbor[OT::surface2N()]
        && neighbor[OT::cornerPPN()] && (!neighbor[OT::edge0PN()] && allocNeighbor[OT::edge0PN()])
        && (!neighbor[OT::edge1NP()] && allocNeighbor[OT::edge1NP()])
        && (neighbor[OT::edge2PP()] && allocNeighbor[OT::edge2PP()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 0, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 1, +1, +1));

        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 2, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 2, -1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 2, -1, -1));
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1P()] && neighbor[OT::surface2P()]
        && neighbor[OT::cornerPPP()] && (!neighbor[OT::edge0PP()] && allocNeighbor[OT::edge0PP()])
        && (!neighbor[OT::edge1PP()] && allocNeighbor[OT::edge1PP()])
        && (neighbor[OT::edge2PP()] && allocNeighbor[OT::edge2PP()]))
    {
        // Validated (because 9 is)
        // pcout << "entering new loop 10" << std::endl;
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 0, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 1, -1, +1));

        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 2, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 2, -1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 2, -1, -1));
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1P()] && neighbor[OT::surface2P()]
        && neighbor[OT::cornerNPP()] && (!neighbor[OT::edge0PP()] && allocNeighbor[OT::edge0PP()])
        && (!neighbor[OT::edge1PN()] && allocNeighbor[OT::edge1PN()])
        && (neighbor[OT::edge2NP()] && allocNeighbor[OT::edge2NP()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 0, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 1, -1, -1));

        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 2, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 2, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 2, +1, -1));
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1P()] && neighbor[OT::surface2N()]
        && neighbor[OT::cornerNPN()] && (!neighbor[OT::edge0PN()] && allocNeighbor[OT::edge0PN()])
        && (!neighbor[OT::edge1NN()] && allocNeighbor[OT::edge1NN()])
        && (neighbor[OT::edge2NP()] && allocNeighbor[OT::edge2NP()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 0, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 1, +1, -1));

        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 2, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 2, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 2, +1, -1));
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1N()] && neighbor[OT::surface2P()]
        && neighbor[OT::cornerPNP()] && (!neighbor[OT::edge0NP()] && allocNeighbor[OT::edge0NP()])
        && (!neighbor[OT::edge1PP()] && allocNeighbor[OT::edge1PP()])
        && (neighbor[OT::edge2PN()] && allocNeighbor[OT::edge2PN()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 0, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 1, -1, +1));

        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 2, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 2, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 2, -1, +1));
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1N()] && neighbor[OT::surface2N()]
        && neighbor[OT::cornerPNN()] && (!neighbor[OT::edge0NN()] && allocNeighbor[OT::edge0NN()])
        && (!neighbor[OT::edge1NP()] && allocNeighbor[OT::edge1NP()])
        && (neighbor[OT::edge2PN()] && allocNeighbor[OT::edge2PN()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 0, -1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 1, +1, +1));

        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 2, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 2, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 2, -1, +1));
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1N()] && neighbor[OT::surface2P()]
        && neighbor[OT::cornerNNP()] && (!neighbor[OT::edge0NP()] && allocNeighbor[OT::edge0NP()])
        && (!neighbor[OT::edge1PN()] && allocNeighbor[OT::edge1PN()])
        && (neighbor[OT::edge2NN()] && allocNeighbor[OT::edge2NN()]))
    {
        // Validated (because 16 is)
        // pcout << "entering new loop 15" << std::endl;
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 0, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 1, -1, -1));

        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 2, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 2, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 2, -1, +1));
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1N()] && neighbor[OT::surface2N()]
        && neighbor[OT::cornerNNN()] && (!neighbor[OT::edge0NN()] && allocNeighbor[OT::edge0NN()])
        && (!neighbor[OT::edge1NN()] && allocNeighbor[OT::edge1NN()])
        && (neighbor[OT::edge2NN()] && allocNeighbor[OT::edge2NN()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 0, -1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 1, +1, -1));

        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 2, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 2, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 2, -1, +1));
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1P()] && neighbor[OT::surface2P()]
        && neighbor[OT::cornerPPP()] && (!neighbor[OT::edge0PP()] && allocNeighbor[OT::edge0PP()])
        && (neighbor[OT::edge1PP()] && allocNeighbor[OT::edge1PP()])
        && (!neighbor[OT::edge2PP()] && allocNeighbor[OT::edge2PP()]))
    {
        // Validated (because 18 is)
        // pcout << "entering new loop 17" << std::endl;
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 0, -1, +1));

        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 1, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 1, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 1, -1, +1));

        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 2, +1, -1));
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1N()] && neighbor[OT::surface2P()]
        && neighbor[OT::cornerPNP()] && (!neighbor[OT::edge0NP()] && allocNeighbor[OT::edge0NP()])
        && (neighbor[OT::edge1PP()] && allocNeighbor[OT::edge1PP()])
        && (!neighbor[OT::edge2PN()] && allocNeighbor[OT::edge2PN()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 0, +1, +1));

        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 1, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 1, +1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 1, -1, +1));

        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 2, +1, +1));
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1P()] && neighbor[OT::surface2N()]
        && neighbor[OT::cornerNPN()] && (!neighbor[OT::edge0PN()] && allocNeighbor[OT::edge0PN()])
        && (neighbor[OT::edge1NN()] && allocNeighbor[OT::edge1NN()])
        && (!neighbor[OT::edge2NP()] && allocNeighbor[OT::edge2NP()]))
    {
        // Validated (because 20 is)
        // pcout << "entering new loop 19" << std::endl;
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 0, -1, +1));

        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 1, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 1, -1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 1, +1, -1));

        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 2, +1, -1));
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1N()] && neighbor[OT::surface2N()]
        && neighbor[OT::cornerNNN()] && (!neighbor[OT::edge0NN()] && allocNeighbor[OT::edge0NN()])
        && (neighbor[OT::edge1NN()] && allocNeighbor[OT::edge1NN()])
        && (!neighbor[OT::edge2NN()] && allocNeighbor[OT::edge2NN()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 0, +1, +1));

        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 1, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 1, -1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 1, +1, -1));

        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 2, +1, +1));
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1P()] && neighbor[OT::surface2P()]
        && neighbor[OT::cornerNPP()] && (!neighbor[OT::edge0PP()] && allocNeighbor[OT::edge0PP()])
        && (neighbor[OT::edge1PN()] && allocNeighbor[OT::edge1PN()])
        && (!neighbor[OT::edge2NP()] && allocNeighbor[OT::edge2NP()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 0, -1, +1));

        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 1, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 1, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 1, -1, +1));

        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 2, -1, -1));
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1N()] && neighbor[OT::surface2P()]
        && neighbor[OT::cornerNNP()] && (!neighbor[OT::edge0NP()] && allocNeighbor[OT::edge0NP()])
        && (neighbor[OT::edge1PN()] && allocNeighbor[OT::edge1PN()])
        && (!neighbor[OT::edge2NN()] && allocNeighbor[OT::edge2NN()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 0, +1, +1));

        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 1, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 1, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 1, -1, +1));

        corners.push_back(
            DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 2, -1, +1));
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1P()] && neighbor[OT::surface2N()]
        && neighbor[OT::cornerPPN()] && (!neighbor[OT::edge0PN()] && allocNeighbor[OT::edge0PN()])
        && (neighbor[OT::edge1NP()] && allocNeighbor[OT::edge1NP()])
        && (!neighbor[OT::edge2PP()] && allocNeighbor[OT::edge2PP()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 0, -1, -1));

        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 1, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 1, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 1, +1, -1));

        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 2, +1, -1));
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1N()] && neighbor[OT::surface2N()]
        && neighbor[OT::cornerPNN()] && (!neighbor[OT::edge0NN()] && allocNeighbor[OT::edge0NN()])
        && (neighbor[OT::edge1NP()] && allocNeighbor[OT::edge1NP()])
        && (!neighbor[OT::edge2PN()] && allocNeighbor[OT::edge2PN()]))
    {
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 0, +1, -1));

        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 1, +1, +1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 1, -1, -1));
        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 1, +1, -1));

        corners.push_back(
            DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 2, +1, +1));
    }

    // UNTIL HERE IT IS NOT VALIDATED

    // Corners must be added on the end of internal corners
    if (neighbor[OT::surface0P()] && neighbor[OT::surface1P()] && neighbor[OT::surface2P()]
        && ((!neighbor[OT::edge0PP()] && allocNeighbor[OT::edge0PP()])
            && (!neighbor[OT::edge1PP()] && allocNeighbor[OT::edge1PP()])
            && (!neighbor[OT::edge2PP()] && allocNeighbor[OT::edge2PP()])
            && (!neighbor[OT::cornerPPP()] && allocNeighbor[OT::cornerPPP()])))
    {
        Box3D bb = Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1);
        corners.push_back(DirectedCorner(bb, 0, +1, +1));
        corners.push_back(DirectedCorner(bb, 0, -1, +1));
        corners.push_back(DirectedCorner(bb, 0, +1, -1));

        corners.push_back(DirectedCorner(bb, 1, +1, +1));
        corners.push_back(DirectedCorner(bb, 1, -1, +1));
        corners.push_back(DirectedCorner(bb, 1, +1, -1));

        corners.push_back(DirectedCorner(bb, 2, +1, +1));
        corners.push_back(DirectedCorner(bb, 2, -1, +1));
        corners.push_back(DirectedCorner(bb, 2, +1, -1));
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1P()] && neighbor[OT::surface2N()]
        && ((!neighbor[OT::edge0PN()] && allocNeighbor[OT::edge0PN()])
            && (!neighbor[OT::edge1NP()] && allocNeighbor[OT::edge1NP()])
            && (!neighbor[OT::edge2PP()] && allocNeighbor[OT::edge2PP()])
            && (!neighbor[OT::cornerPPN()] && allocNeighbor[OT::cornerPPN()])))
    {
        Box3D bb = Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0);
        corners.push_back(DirectedCorner(bb, 0, +1, -1));
        corners.push_back(DirectedCorner(bb, 0, +1, +1));
        corners.push_back(DirectedCorner(bb, 0, -1, -1));

        corners.push_back(DirectedCorner(bb, 1, -1, +1));
        corners.push_back(DirectedCorner(bb, 1, +1, +1));
        corners.push_back(DirectedCorner(bb, 1, -1, -1));

        corners.push_back(DirectedCorner(bb, 2, +1, +1));
        corners.push_back(DirectedCorner(bb, 2, +1, -1));
        corners.push_back(DirectedCorner(bb, 2, -1, +1));
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1N()] && neighbor[OT::surface2P()]
        && ((!neighbor[OT::edge0NP()] && allocNeighbor[OT::edge0NP()])
            && (!neighbor[OT::edge1PP()] && allocNeighbor[OT::edge1PP()])
            && (!neighbor[OT::edge2PN()] && allocNeighbor[OT::edge2PN()])
            && (!neighbor[OT::cornerPNP()] && allocNeighbor[OT::cornerPNP()])))
    {
        Box3D bb = Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1);
        corners.push_back(DirectedCorner(bb, 0, -1, +1));
        corners.push_back(DirectedCorner(bb, 0, -1, -1));
        corners.push_back(DirectedCorner(bb, 0, +1, +1));

        corners.push_back(DirectedCorner(bb, 1, +1, +1));
        corners.push_back(DirectedCorner(bb, 1, -1, +1));
        corners.push_back(DirectedCorner(bb, 1, +1, -1));

        corners.push_back(DirectedCorner(bb, 2, +1, -1));
        corners.push_back(DirectedCorner(bb, 2, -1, -1));
        corners.push_back(DirectedCorner(bb, 2, +1, +1));
    }

    if (neighbor[OT::surface0P()] && neighbor[OT::surface1N()] && neighbor[OT::surface2N()]
        && ((!neighbor[OT::edge0NN()] && allocNeighbor[OT::edge0NN()])
            && (!neighbor[OT::edge1NP()] && allocNeighbor[OT::edge1NP()])
            && (!neighbor[OT::edge2PN()] && allocNeighbor[OT::edge2PN()])
            && (!neighbor[OT::cornerPNN()] && allocNeighbor[OT::cornerPNN()])))
    {
        Box3D bb = Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0);
        corners.push_back(DirectedCorner(bb, 0, -1, -1));
        corners.push_back(DirectedCorner(bb, 0, -1, +1));
        corners.push_back(DirectedCorner(bb, 0, +1, -1));

        corners.push_back(DirectedCorner(bb, 1, -1, +1));
        corners.push_back(DirectedCorner(bb, 1, -1, -1));
        corners.push_back(DirectedCorner(bb, 1, +1, +1));

        corners.push_back(DirectedCorner(bb, 2, +1, -1));
        corners.push_back(DirectedCorner(bb, 2, -1, -1));
        corners.push_back(DirectedCorner(bb, 2, +1, +1));
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1P()] && neighbor[OT::surface2P()]
        && ((!neighbor[OT::edge0PP()] && allocNeighbor[OT::edge0PP()])
            && (!neighbor[OT::edge1PN()] && allocNeighbor[OT::edge1PN()])
            && (!neighbor[OT::edge2NP()] && allocNeighbor[OT::edge2NP()])
            && (!neighbor[OT::cornerNPP()] && allocNeighbor[OT::cornerNPP()])))
    {
        Box3D bb = Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1);
        corners.push_back(DirectedCorner(bb, 0, +1, +1));
        corners.push_back(DirectedCorner(bb, 0, -1, +1));
        corners.push_back(DirectedCorner(bb, 0, +1, -1));

        corners.push_back(DirectedCorner(bb, 1, +1, -1));
        corners.push_back(DirectedCorner(bb, 1, -1, -1));
        corners.push_back(DirectedCorner(bb, 1, +1, +1));

        corners.push_back(DirectedCorner(bb, 2, -1, +1));
        corners.push_back(DirectedCorner(bb, 2, -1, -1));
        corners.push_back(DirectedCorner(bb, 2, +1, +1));
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1P()] && neighbor[OT::surface2N()]
        && ((!neighbor[OT::edge0PN()] && allocNeighbor[OT::edge0PN()])
            && (!neighbor[OT::edge1NN()] && allocNeighbor[OT::edge1NN()])
            && (!neighbor[OT::edge2NP()] && allocNeighbor[OT::edge2NP()])
            && (!neighbor[OT::cornerNPN()] && allocNeighbor[OT::cornerNPN()])))
    {
        Box3D bb = Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0);
        corners.push_back(DirectedCorner(bb, 0, +1, -1));
        corners.push_back(DirectedCorner(bb, 0, +1, +1));
        corners.push_back(DirectedCorner(bb, 0, -1, -1));

        corners.push_back(DirectedCorner(bb, 1, -1, -1));
        corners.push_back(DirectedCorner(bb, 1, -1, +1));
        corners.push_back(DirectedCorner(bb, 1, +1, -1));

        corners.push_back(DirectedCorner(bb, 2, -1, +1));
        corners.push_back(DirectedCorner(bb, 2, -1, -1));
        corners.push_back(DirectedCorner(bb, 2, +1, +1));
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1N()] && neighbor[OT::surface2P()]
        && ((!neighbor[OT::edge0NP()] && allocNeighbor[OT::edge0NP()])
            && (!neighbor[OT::edge1PN()] && allocNeighbor[OT::edge1PN()])
            && (!neighbor[OT::edge2NN()] && allocNeighbor[OT::edge2NN()])
            && (!neighbor[OT::cornerNNP()] && allocNeighbor[OT::cornerNNP()])))
    {
        Box3D bb = Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1);
        corners.push_back(DirectedCorner(bb, 0, -1, +1));
        corners.push_back(DirectedCorner(bb, 0, -1, -1));
        corners.push_back(DirectedCorner(bb, 0, +1, +1));

        corners.push_back(DirectedCorner(bb, 1, +1, -1));
        corners.push_back(DirectedCorner(bb, 1, +1, +1));
        corners.push_back(DirectedCorner(bb, 1, -1, -1));

        corners.push_back(DirectedCorner(bb, 2, -1, -1));
        corners.push_back(DirectedCorner(bb, 2, -1, +1));
        corners.push_back(DirectedCorner(bb, 2, +1, -1));
    }

    if (neighbor[OT::surface0N()] && neighbor[OT::surface1N()] && neighbor[OT::surface2N()]
        && ((!neighbor[OT::edge0NN()] && allocNeighbor[OT::edge0NN()])
            && (!neighbor[OT::edge1NN()] && allocNeighbor[OT::edge1NN()])
            && (!neighbor[OT::edge2NN()] && allocNeighbor[OT::edge2NN()])
            && (!neighbor[OT::cornerNNN()] && allocNeighbor[OT::cornerNNN()])))
    {
        Box3D bb = Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0);
        corners.push_back(DirectedCorner(bb, 0, -1, -1));
        corners.push_back(DirectedCorner(bb, 0, -1, +1));
        corners.push_back(DirectedCorner(bb, 0, +1, -1));

        corners.push_back(DirectedCorner(bb, 1, -1, -1));
        corners.push_back(DirectedCorner(bb, 1, -1, +1));
        corners.push_back(DirectedCorner(bb, 1, +1, -1));

        corners.push_back(DirectedCorner(bb, 2, -1, -1));
        corners.push_back(DirectedCorner(bb, 2, -1, +1));
        corners.push_back(DirectedCorner(bb, 2, +1, -1));
    }

    PLB_ASSERT(containsNoDuplicates(corners, false) && "Corners must be unique.");
    // PLB_ASSERT(corners.size() == 24 && "There are 24 directed corners in a box.");

    return corners;
}

std::vector<DirectedCorner> getBoundaryDirectedCorners(
    const Box3D &box, const Array<bool, 26> &neighbor, const Array<bool, 26> &bcNeighbor,
    const Array<bool, 26> &allocNeighbor)
{
    std::vector<DirectedCorner> corners;

    if (bcNeighbor[OT::surface0N()]) {
        // exterior corners
        if ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
            && !bcNeighbor[OT::surface1N()]
            && (!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
            && !bcNeighbor[OT::surface2N()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 1, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 2, -1, -1));
        }

        if ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
            && !bcNeighbor[OT::surface1N()]
            && (!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
            && !bcNeighbor[OT::surface2P()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 1, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 2, -1, -1));
        }

        if ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
            && !bcNeighbor[OT::surface1P()]
            && (!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
            && !bcNeighbor[OT::surface2N()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 1, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 2, -1, +1));
        }

        if ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
            && !bcNeighbor[OT::surface1P()]
            && (!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
            && !bcNeighbor[OT::surface2P()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 1, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 2, -1, +1));
        }
        // interior corners
        if (neighbor[OT::surface1N()] && !bcNeighbor[OT::surface1N()] && neighbor[OT::surface2N()]
            && !bcNeighbor[OT::surface2N()]
            && (!neighbor[OT::edge0NN()] && allocNeighbor[OT::edge0NN()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 1, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 2, -1, +1));
        }

        if (neighbor[OT::surface1N()] && !bcNeighbor[OT::surface1N()] && neighbor[OT::surface2P()]
            && !bcNeighbor[OT::surface2P()]
            && (!neighbor[OT::edge0NP()] && allocNeighbor[OT::edge0NP()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 1, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 2, -1, +1));
        }

        if (neighbor[OT::surface1P()] && !bcNeighbor[OT::surface1P()] && neighbor[OT::surface2N()]
            && !bcNeighbor[OT::surface2N()]
            && (!neighbor[OT::edge0PN()] && allocNeighbor[OT::edge0PN()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 1, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 2, -1, -1));
        }

        if (neighbor[OT::surface1P()] && !bcNeighbor[OT::surface1P()] && neighbor[OT::surface2P()]
            && !bcNeighbor[OT::surface2P()]
            && (!neighbor[OT::edge0PP()] && allocNeighbor[OT::edge0PP()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 1, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 2, -1, -1));
        }
    }

    if (bcNeighbor[OT::surface0P()]) {
        // exterior corners
        if ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
            && !bcNeighbor[OT::surface1N()]
            && (!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
            && !bcNeighbor[OT::surface2N()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 1, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 2, +1, -1));
        }

        if ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
            && !bcNeighbor[OT::surface1N()]
            && (!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
            && !bcNeighbor[OT::surface2P()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 1, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 2, +1, -1));
        }

        if ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
            && !bcNeighbor[OT::surface1P()]
            && (!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
            && !bcNeighbor[OT::surface2N()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 1, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 2, +1, +1));
        }

        if ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
            && !bcNeighbor[OT::surface1P()]
            && (!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
            && !bcNeighbor[OT::surface2P()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 1, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 2, +1, -1));
        }

        // interior corners
        if (neighbor[OT::surface1N()] && !bcNeighbor[OT::surface1N()] && neighbor[OT::surface2N()]
            && !bcNeighbor[OT::surface2N()]
            && (!neighbor[OT::edge0NN()] && allocNeighbor[OT::edge0NN()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 1, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 2, +1, +1));
        }
        if (neighbor[OT::surface1N()] && !bcNeighbor[OT::surface1N()] && neighbor[OT::surface2P()]
            && !bcNeighbor[OT::surface2P()]
            && (!neighbor[OT::edge0NP()] && allocNeighbor[OT::edge0NP()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 1, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 2, +1, +1));
        }
        if (neighbor[OT::surface1P()] && !bcNeighbor[OT::surface1P()] && neighbor[OT::surface2N()]
            && !bcNeighbor[OT::surface2N()]
            && (!neighbor[OT::edge0PN()] && allocNeighbor[OT::edge0PN()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 1, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 2, +1, -1));
        }

        if (neighbor[OT::surface1P()] && !bcNeighbor[OT::surface1P()] && neighbor[OT::surface2P()]
            && !bcNeighbor[OT::surface2P()]
            && (!neighbor[OT::edge0PP()] && allocNeighbor[OT::edge0PP()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 1, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 2, +1, -1));
        }
    }

    if (bcNeighbor[OT::surface1N()]) {
        // exterior corners

        if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
            && !bcNeighbor[OT::surface0N()]
            && (!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
            && !bcNeighbor[OT::surface2N()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 0, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 2, -1, -1));
        }

        if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
            && !bcNeighbor[OT::surface0N()]
            && (!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
            && !bcNeighbor[OT::surface2P()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 0, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 2, -1, -1));
        }

        if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
            && !bcNeighbor[OT::surface0P()]
            && (!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
            && !bcNeighbor[OT::surface2N()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 0, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 2, +1, -1));
        }

        if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
            && !bcNeighbor[OT::surface0P()]
            && (!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
            && !bcNeighbor[OT::surface2P()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 0, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 2, +1, -1));
        }

        // interior corners
        if (neighbor[OT::surface0N()] && !bcNeighbor[OT::surface0N()] && neighbor[OT::surface2N()]
            && !bcNeighbor[OT::surface2N()]
            && (!neighbor[OT::edge1NN()] && allocNeighbor[OT::edge1NN()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 0, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 2, +1, -1));
        }

        if (neighbor[OT::surface0N()] && !bcNeighbor[OT::surface0N()] && neighbor[OT::surface2P()]
            && !bcNeighbor[OT::surface2P()]
            && (!neighbor[OT::edge1PN()] && allocNeighbor[OT::edge1PN()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 0, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 2, +1, -1));
        }

        if (neighbor[OT::surface0P()] && !bcNeighbor[OT::surface0P()] && neighbor[OT::surface2N()]
            && !bcNeighbor[OT::surface2N()]
            && (!neighbor[OT::edge1NP()] && allocNeighbor[OT::edge1NP()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 0, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 2, -1, -1));
        }

        if (neighbor[OT::surface0P()] && !bcNeighbor[OT::surface0P()] && neighbor[OT::surface2P()]
            && !bcNeighbor[OT::surface2P()]
            && (!neighbor[OT::edge1PP()] && allocNeighbor[OT::edge1PP()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 0, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 2, -1, -1));
        }
    }

    if (bcNeighbor[OT::surface1P()]) {
        // exterior corners
        if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
            && !bcNeighbor[OT::surface0N()]
            && (!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
            && !bcNeighbor[OT::surface2N()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 0, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 2, -1, +1));
        }

        if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
            && !bcNeighbor[OT::surface0N()]
            && (!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
            && !bcNeighbor[OT::surface2P()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 0, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 2, -1, +1));
        }

        if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
            && !bcNeighbor[OT::surface0P()]
            && (!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])
            && !bcNeighbor[OT::surface2N()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 0, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 2, +1, +1));
        }

        if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
            && !bcNeighbor[OT::surface0P()]
            && (!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])
            && !bcNeighbor[OT::surface2P()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 0, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 2, +1, +1));
        }

        // interior corners
        if (neighbor[OT::surface0N()] && !bcNeighbor[OT::surface0N()] && neighbor[OT::surface2N()]
            && !bcNeighbor[OT::surface2N()]
            && (!neighbor[OT::edge1NN()] && allocNeighbor[OT::edge1NN()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 0, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 2, +1, +1));
        }

        if (neighbor[OT::surface0N()] && !bcNeighbor[OT::surface0N()] && neighbor[OT::surface2P()]
            && !bcNeighbor[OT::surface2P()]
            && (!neighbor[OT::edge1PN()] && allocNeighbor[OT::edge1PN()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 0, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 2, +1, +1));
        }

        if (neighbor[OT::surface0P()] && !bcNeighbor[OT::surface0P()] && neighbor[OT::surface2N()]
            && !bcNeighbor[OT::surface2N()]
            && (!neighbor[OT::edge1NP()] && allocNeighbor[OT::edge1NP()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 0, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 2, -1, +1));
        }

        if (neighbor[OT::surface0P()] && !bcNeighbor[OT::surface0P()] && neighbor[OT::surface2P()]
            && !bcNeighbor[OT::surface2P()]
            && (!neighbor[OT::edge1PP()] && allocNeighbor[OT::edge1PP()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 0, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 2, -1, +1));
        }
    }

    if (bcNeighbor[OT::surface2N()]) {
        // exterior corners

        if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
            && !bcNeighbor[OT::surface0N()]
            && (!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
            && !bcNeighbor[OT::surface1N()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 0, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 1, -1, -1));
        }

        if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
            && !bcNeighbor[OT::surface0N()]
            && (!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
            && !bcNeighbor[OT::surface1P()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 0, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 1, -1, -1));
        }

        if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
            && !bcNeighbor[OT::surface0P()]
            && (!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
            && !bcNeighbor[OT::surface1N()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 0, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 1, -1, +1));
        }

        if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
            && !bcNeighbor[OT::surface0P()]
            && (!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
            && !bcNeighbor[OT::surface1P()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 0, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 1, -1, +1));
        }
        // interior corners
        if (neighbor[OT::surface0N()] && !bcNeighbor[OT::surface0N()] && neighbor[OT::surface1N()]
            && !bcNeighbor[OT::surface1N()]
            && (!neighbor[OT::edge2NN()] && allocNeighbor[OT::edge2NN()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 0, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 1, -1, +1));
        }

        if (neighbor[OT::surface0N()] && !bcNeighbor[OT::surface0N()] && neighbor[OT::surface1P()]
            && !bcNeighbor[OT::surface1P()]
            && (!neighbor[OT::edge2NP()] && allocNeighbor[OT::edge2NP()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 0, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 1, -1, +1));
        }

        if (neighbor[OT::surface0P()] && !bcNeighbor[OT::surface0P()] && neighbor[OT::surface1N()]
            && !bcNeighbor[OT::surface1N()]
            && (!neighbor[OT::edge2PN()] && allocNeighbor[OT::edge2PN()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 0, +1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 1, -1, -1));
        }

        if (neighbor[OT::surface0P()] && !bcNeighbor[OT::surface0P()] && neighbor[OT::surface1P()]
            && !bcNeighbor[OT::surface1P()]
            && (!neighbor[OT::edge2PP()] && allocNeighbor[OT::edge2PP()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 0, -1, -1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 1, -1, -1));
        }
    }

    if (bcNeighbor[OT::surface2P()]) {
        // exterior corners
        if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
            && !bcNeighbor[OT::surface0N()]
            && (!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
            && !bcNeighbor[OT::surface1N()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 0, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 1, +1, -1));
        }

        if ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])
            && !bcNeighbor[OT::surface0N()]
            && (!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
            && !bcNeighbor[OT::surface1P()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 0, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 1, +1, -1));
        }

        if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
            && !bcNeighbor[OT::surface0P()]
            && (!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])
            && !bcNeighbor[OT::surface1N()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 0, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 1, +1, +1));
        }

        if ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])
            && !bcNeighbor[OT::surface0P()]
            && (!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])
            && !bcNeighbor[OT::surface1P()])
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 0, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 1, +1, +1));
        }
        // interior corners
        if (neighbor[OT::surface0N()] && !bcNeighbor[OT::surface0N()] && neighbor[OT::surface1N()]
            && !bcNeighbor[OT::surface1N()]
            && (!neighbor[OT::edge2NN()] && allocNeighbor[OT::edge2NN()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 0, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 1, +1, +1));
        }

        if (neighbor[OT::surface0N()] && !bcNeighbor[OT::surface0N()] && neighbor[OT::surface1P()]
            && !bcNeighbor[OT::surface1P()]
            && (!neighbor[OT::edge2NP()] && allocNeighbor[OT::edge2NP()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 0, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 1, +1, +1));
        }

        if (neighbor[OT::surface0P()] && !bcNeighbor[OT::surface0P()] && neighbor[OT::surface1N()]
            && !bcNeighbor[OT::surface1N()]
            && (!neighbor[OT::edge2PN()] && allocNeighbor[OT::edge2PN()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 0, +1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 1, +1, -1));
        }

        if (neighbor[OT::surface0P()] && !bcNeighbor[OT::surface0P()] && neighbor[OT::surface1P()]
            && !bcNeighbor[OT::surface1P()]
            && (!neighbor[OT::edge2PP()] && allocNeighbor[OT::edge2PP()]))
        {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 0, -1, +1));
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 1, +1, -1));
        }
    }

    // Boundary corners on boundary edges
    // ============ Z-EDGE ============ //
    if (bcNeighbor[OT::surface0P()] && bcNeighbor[OT::surface1P()]) {
        if (!bcNeighbor[OT::surface2P()]
            && (!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 2, +1, +1));
        }

        if (!bcNeighbor[OT::surface2N()]
            && (!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 2, +1, +1));
        }
    }
    if (bcNeighbor[OT::surface0P()] && bcNeighbor[OT::surface1N()]) {
        if (!bcNeighbor[OT::surface2P()]
            && (!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 2, +1, -1));
        }

        if (!bcNeighbor[OT::surface2N()]
            && (!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 2, +1, -1));
        }
    }

    if (bcNeighbor[OT::surface0N()] && bcNeighbor[OT::surface1P()]) {
        if (!bcNeighbor[OT::surface2P()]
            && (!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 2, -1, +1));
        }

        if (!bcNeighbor[OT::surface2N()]
            && (!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 2, -1, +1));
        }
    }

    if (bcNeighbor[OT::surface0N()] && bcNeighbor[OT::surface1N()]) {
        if (!bcNeighbor[OT::surface2P()]
            && (!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 2, -1, -1));
        }

        if (!bcNeighbor[OT::surface2N()]
            && (!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 2, -1, -1));
        }
    }

    // ============ Y-EDGE ============ //
    if (bcNeighbor[OT::surface0P()] && bcNeighbor[OT::surface2P()]) {
        if (!bcNeighbor[OT::surface1P()]
            && (!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 1, +1, +1));
        }

        if (!bcNeighbor[OT::surface1N()]
            && (!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 1, +1, +1));
        }
    }
    if (bcNeighbor[OT::surface0P()] && bcNeighbor[OT::surface2N()]) {
        if (!bcNeighbor[OT::surface1P()]
            && (!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 1, -1, +1));
        }

        if (!bcNeighbor[OT::surface1N()]
            && (!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 1, -1, +1));
        }
    }

    if (bcNeighbor[OT::surface0N()] && bcNeighbor[OT::surface2P()]) {
        if (!bcNeighbor[OT::surface1P()]
            && (!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 1, +1, -1));
        }

        if (!bcNeighbor[OT::surface1N()]
            && (!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 1, +1, -1));
        }
    }

    if (bcNeighbor[OT::surface0N()] && bcNeighbor[OT::surface2N()]) {
        if (!bcNeighbor[OT::surface1P()]
            && (!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 1, -1, -1));
        }

        if (!bcNeighbor[OT::surface1N()]
            && (!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 1, -1, -1));
        }
    }

    // ============ X-EDGE ============ //
    if (bcNeighbor[OT::surface1P()] && bcNeighbor[OT::surface2P()]) {
        if (!bcNeighbor[OT::surface0P()]
            && (!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z1, box.z1), 0, +1, +1));
        }

        if (!bcNeighbor[OT::surface0N()]
            && (!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z1, box.z1), 0, +1, +1));
        }
    }
    if (bcNeighbor[OT::surface1P()] && bcNeighbor[OT::surface2N()]) {
        if (!bcNeighbor[OT::surface0P()]
            && (!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y1, box.y1, box.z0, box.z0), 0, +1, -1));
        }

        if (!bcNeighbor[OT::surface0N()]
            && (!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y1, box.y1, box.z0, box.z0), 0, +1, -1));
        }
    }

    if (bcNeighbor[OT::surface1N()] && bcNeighbor[OT::surface2P()]) {
        if (!bcNeighbor[OT::surface0P()]
            && (!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z1, box.z1), 0, -1, +1));
        }

        if (!bcNeighbor[OT::surface0N()]
            && (!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z1, box.z1), 0, -1, +1));
        }
    }

    if (bcNeighbor[OT::surface1N()] && bcNeighbor[OT::surface2N()]) {
        if (!bcNeighbor[OT::surface0P()]
            && (!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x1, box.x1, box.y0, box.y0, box.z0, box.z0), 0, -1, -1));
        }

        if (!bcNeighbor[OT::surface0N()]
            && (!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()])) {
            corners.push_back(
                DirectedCorner(Box3D(box.x0, box.x0, box.y0, box.y0, box.z0, box.z0), 0, -1, -1));
        }
    }

    PLB_ASSERT(containsNoDuplicates(corners, false) && "Corners must be unique.");
    // PLB_ASSERT(corners.size() == 24 && "There are 24 directed corners in a box.");

    return corners;
}

bool hasFinerNeighbor(OctreeGridStructure const &ogs, plint blockId)
{
    for (plint iA = 0; iA < 26; ++iA) {
        if (ogs.neighborIsAtFinerLevel(blockId, iA)) {
            return true;
        }
    }

    return false;
}

bool hasCoarserNeighbor(OctreeGridStructure const &ogs, plint blockId)
{
    for (plint iA = 0; iA < 26; ++iA) {
        if (ogs.neighborIsAtCoarserLevel(blockId, iA)) {
            return true;
        }
    }

    return false;
}

Array<bool, 26> getNeighbors(OctreeGridStructure const &ogs, plint blockId)
{
    Array<bool, 26> neighbors;

    for (plint iA = 0; iA < 26; ++iA) {
        neighbors[iA] =
            (ogs.neighborIsAtSameLevel(blockId, iA) || ogs.neighborIsAtFinerLevel(blockId, iA))
            && ogs.neighborIsAllocated(blockId, iA);
    }

    return neighbors;
}

Array<bool, 26> getBcNeighbors(OctreeGridStructure const &ogs, plint blockId)
{
    Array<bool, 26> neighbors;

    for (plint iA = 0; iA < 26; ++iA) {
        neighbors[iA] = ogs.neighborIsBoundary(blockId, iA);
    }

    return neighbors;
}

Array<bool, 26> getAllocatedNeighbors(OctreeGridStructure const &ogs, plint blockId)
{
    Array<bool, 26> neighbors;

    for (plint iA = 0; iA < 26; ++iA) {
        neighbors[iA] = ogs.neighborIsAllocated(blockId, iA);
    }

    return neighbors;
}

std::pair<plint, plint> getDirections(plint dir)
{
    switch (dir) {
    case 0:
        return std::make_pair<plint, plint>(OT::surface0N(), OT::surface0P());
        break;
    case 1:
        return std::make_pair<plint, plint>(OT::surface1N(), OT::surface1P());
        break;
    case 2:
        return std::make_pair<plint, plint>(OT::surface2N(), OT::surface2P());
        break;
    default:
        return std::make_pair<plint, plint>(-1, -1);
        break;
    }
    return std::make_pair<plint, plint>(-1, -1);
}

std::pair<Array<plint, 2>, Array<plint, 2> > getEdgesIndices(
    plint dir, plint ori, const Array<bool, 26> &neighbors, const Array<bool, 26> &bcNeighbors,
    plint l1, plint l2)
{
    PLB_ASSERT((dir == 0 || dir == 1 || dir == 2) && "Dir must be 0,1, or 2.");
    PLB_ASSERT((ori == -1 || ori == +1) && "Ori must be -1, or +1.");

    const std::pair<plint, plint> dirOne = getDirections((dir + 1) % 3);
    const std::pair<plint, plint> dirTwo = getDirections((dir + 2) % 3);

    Array<plint, 2> dirOneExt(0, 0);
    Array<plint, 2> dirTwoExt(0, 0);

    if (dir == 0) {
        if (ori == +1) {
            if (neighbors[dirOne.first] && !bcNeighbors[dirOne.first] && neighbors[OT::edge2PN()]) {
                dirOneExt[0] = l1;
            }
            if (neighbors[dirOne.second] && !bcNeighbors[dirOne.second] && neighbors[OT::edge2PP()])
            {
                dirOneExt[1] = l2;
            }
            if (neighbors[dirTwo.first] && !bcNeighbors[dirTwo.first] && neighbors[OT::edge1NP()]) {
                dirTwoExt[0] = l1;
            }
            if (neighbors[dirTwo.second] && !bcNeighbors[dirTwo.second] && neighbors[OT::edge1PP()])
            {
                dirTwoExt[1] = l2;
            }
        } else if (ori == -1) {
            if (neighbors[dirOne.first] && !bcNeighbors[dirOne.first] && neighbors[OT::edge2NN()]) {
                dirOneExt[0] = l1;
            }
            if (neighbors[dirOne.second] && !bcNeighbors[dirOne.second] && neighbors[OT::edge2NP()])
            {
                dirOneExt[1] = l2;
            }
            if (neighbors[dirTwo.first] && !bcNeighbors[dirTwo.first] && neighbors[OT::edge1NN()]) {
                dirTwoExt[0] = l1;
            }
            if (neighbors[dirTwo.second] && !bcNeighbors[dirTwo.second] && neighbors[OT::edge1PN()])
            {
                dirTwoExt[1] = l2;
            }
        }
    } else if (dir == 1) {
        if (ori == +1) {
            if (neighbors[dirOne.first] && !bcNeighbors[dirOne.first] && neighbors[OT::edge0PN()]) {
                dirOneExt[0] = l1;
            }
            if (neighbors[dirOne.second] && !bcNeighbors[dirOne.second] && neighbors[OT::edge0PP()])
            {
                dirOneExt[1] = l2;
            }
            if (neighbors[dirTwo.first] && !bcNeighbors[dirTwo.first] && neighbors[OT::edge2NP()]) {
                dirTwoExt[0] = l1;
            }
            if (neighbors[dirTwo.second] && !bcNeighbors[dirTwo.second] && neighbors[OT::edge2PP()])
            {
                dirTwoExt[1] = l2;
            }
        } else if (ori == -1) {
            if (neighbors[dirOne.first] && !bcNeighbors[dirOne.first] && neighbors[OT::edge0NN()]) {
                dirOneExt[0] = l1;
            }
            if (neighbors[dirOne.second] && !bcNeighbors[dirOne.second] && neighbors[OT::edge0NP()])
            {
                dirOneExt[1] = l2;
            }
            if (neighbors[dirTwo.first] && !bcNeighbors[dirTwo.first] && neighbors[OT::edge2NN()]) {
                dirTwoExt[0] = l1;
            }
            if (neighbors[dirTwo.second] && !bcNeighbors[dirTwo.second] && neighbors[OT::edge2PN()])
            {
                dirTwoExt[1] = l2;
            }
        }
    } else if (dir == 2) {
        if (ori == +1) {
            if (neighbors[dirOne.first] && !bcNeighbors[dirOne.first] && neighbors[OT::edge1PN()]) {
                dirOneExt[0] = l1;
            }
            if (neighbors[dirOne.second] && !bcNeighbors[dirOne.second] && neighbors[OT::edge1PP()])
            {
                dirOneExt[1] = l2;
            }
            if (neighbors[dirTwo.first] && !bcNeighbors[dirTwo.first] && neighbors[OT::edge0NP()]) {
                dirTwoExt[0] = l1;
            }
            if (neighbors[dirTwo.second] && !bcNeighbors[dirTwo.second] && neighbors[OT::edge0PP()])
            {
                dirTwoExt[1] = l2;
            }
        } else if (ori == -1) {
            if (neighbors[dirOne.first] && !bcNeighbors[dirOne.first] && neighbors[OT::edge1NN()]) {
                dirOneExt[0] = l1;
            }
            if (neighbors[dirOne.second] && !bcNeighbors[dirOne.second] && neighbors[OT::edge1NP()])
            {
                dirOneExt[1] = l2;
            }
            if (neighbors[dirTwo.first] && !bcNeighbors[dirTwo.first] && neighbors[OT::edge0NN()]) {
                dirTwoExt[0] = l1;
            }
            if (neighbors[dirTwo.second] && !bcNeighbors[dirTwo.second] && neighbors[OT::edge0PN()])
            {
                dirTwoExt[1] = l2;
            }
        }
    }

    return std::pair<Array<plint, 2>, Array<plint, 2> >(dirOneExt, dirTwoExt);
}

void printBox(Box3D const &b, std::string name)
{
    pcout << name << " = { " << b.x0 << ", " << b.x1 << "; " << b.y0 << ", " << b.y1 << "; " << b.z0
          << ", " << b.z1 << " }" << std::endl;
}

void printBoxes(
    std::vector<Box3D> const &b, std::string fname, plint level, double dx,
    Array<double, 3> const &pos)
{
    TriangleSet<double> *triangleSet = 0;
    Array<plint, 3> nSegments = Array<plint, 3>(1, 1, 1);

    triangleSet = new TriangleSet<double>(DBL);

    for (plint iA = 0; iA < (plint)b.size(); ++iA) {
        TriangleSet<double> cuboidSet = constructCuboid<double>(
            Cuboid<double>(b[iA].multiply(util::intTwoToThePower(level))), nSegments);
        triangleSet->append(cuboidSet);
    }
    triangleSet->scale(dx);
    triangleSet->translate(pos);
    triangleSet->writeBinarySTL(fname);
    delete triangleSet;
}

void makeIntersectionsConsistent(
    const Box3D &box, const std::vector<Box3D> &boxes, std::vector<Box3D> &intBoxes)
{
    for (plint iA = 0; iA < (plint)boxes.size(); ++iA) {
        Box3D bb;
        bool doesIntersect = intersect(box, boxes[iA], bb);
        if (doesIntersect) {
            intBoxes.push_back(bb);
        }
    }
}

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
)
{
    // No neighbors of the same level and no coarser neighbors
    // therefore one must first reduce the size of the Box3D in the
    // in the x1, y1, z1 direction.

    // =========== Coarse to Fine interfaces ================= //

    // this box represents the interface Coarse to fine, in fine units
    Box3D coarseToFineBoxFineUnits = Box3D(
        box.x0, box.x1 - !(neighbor[OT::surface0P()]), box.y0,
        box.y1 - !(neighbor[OT::surface1P()]), box.z0, box.z1 - !(neighbor[OT::surface2P()]));

    // from fine to coarse one does not need directed planes/edges/corners.
    // Nevertheless this contains only the interior of the surfaces of the
    // Box3D and therefore the edges/corners must be added
    std::vector<boxLogic::DirectedPlane> coarseToFinePlanesFineUnitsTmp =
        boxLogic::getDirectedPlanes(coarseToFineBoxFineUnits, neighbor, bcNeighbor, allocNeighbor);

    // there will the interface boxes which are extended by 1 fine lattice points in
    // each direction normal to the plane because the planes are defined without the edges.
    for (plint iA = 0; iA < (plint)coarseToFinePlanesFineUnitsTmp.size(); ++iA) {
        plint dir = coarseToFinePlanesFineUnitsTmp[iA].direction;
        plint ori = coarseToFinePlanesFineUnitsTmp[iA].orientation;
        std::pair<Array<plint, 2>, Array<plint, 2> > nEdge =
            getEdgesIndices(dir, ori, neighbor, bcNeighbor, 2, 1);

        coarseToFinePlanesFineUnits.push_back(
            coarseToFinePlanesFineUnitsTmp[iA]
                .bb.enlargeInNormalPlane(1, coarseToFinePlanesFineUnitsTmp[iA].direction)
                .enlargeInNormalPlane(nEdge.first, nEdge.second, dir));
    }

    // The box must be enlarged in the direction where there is a neighbor of the same level.
    Array<plint, 6> enlargeAryCoarse(
        -neighbor[OT::surface0N()], neighbor[OT::surface0P()], -neighbor[OT::surface1N()],
        neighbor[OT::surface1P()], -neighbor[OT::surface2N()], neighbor[OT::surface2P()]);
    // this box represents the interface Coarse to fine, in coarse units
    Box3D coarseToFineBoxCoarseUnits = coarseToFineBoxFineUnits.divide(2);
    std::vector<boxLogic::DirectedPlane> coarseToFinePlanesCoarseUnitsTmp =
        boxLogic::getDirectedPlanes(
            coarseToFineBoxCoarseUnits, neighbor, bcNeighbor, allocNeighbor);
    for (plint iA = 0; iA < (plint)coarseToFinePlanesCoarseUnitsTmp.size(); ++iA) {
        plint dir = coarseToFinePlanesCoarseUnitsTmp[iA].direction;
        plint ori = coarseToFinePlanesCoarseUnitsTmp[iA].orientation;

        std::pair<plint, plint> dirOne = getDirections((dir + 1) % 3);
        Array<plint, 2> nCellsOne(neighbor[dirOne.first], neighbor[dirOne.second]);

        std::pair<plint, plint> dirTwo = getDirections((dir + 2) % 3);
        Array<plint, 2> nCellsTwo(neighbor[dirTwo.first], neighbor[dirTwo.second]);
        std::pair<Array<plint, 2>, Array<plint, 2> > nEdge =
            getEdgesIndices(dir, ori, neighbor, bcNeighbor, 2, 2);

        boxLogic::DirectedPlane dPlane(
            coarseToFinePlanesCoarseUnitsTmp[iA].bb.enlargeInNormalPlane(nCellsOne, nCellsTwo, dir),
            dir, ori);

        // enlarge to account for the eventual extra cells needed for interpolation
        coarseToFinePlanesCoarseUnits.push_back(dPlane);

        nCellsOne = Array<plint, 2>(
            2 * (!neighbor[dirOne.first] && !bcNeighbor[dirOne.first]) + bcNeighbor[dirOne.first],
            2 * (!neighbor[dirOne.second] && !bcNeighbor[dirOne.second])
                + bcNeighbor[dirOne.second]);

        nCellsTwo = Array<plint, 2>(
            2 * (!neighbor[dirTwo.first] && !bcNeighbor[dirTwo.first]) + bcNeighbor[dirTwo.first],
            2 * (!neighbor[dirTwo.second] && !bcNeighbor[dirTwo.second])
                + bcNeighbor[dirTwo.second]);
        // enlarge to account for the eventual extra cells needed for interpolation
        coarseToFinePlanesCoarseUnitsExtended.push_back(boxLogic::DirectedPlane(
            dPlane.bb.enlargeInNormalPlane(nCellsOne + nEdge.first, nCellsTwo + nEdge.second, dir),
            dir, ori));
    }

    std::vector<boxLogic::DirectedEdge> coarseToFineEdgesCoarseUnitsTmp =
        boxLogic::getDirectedEdges(
            coarseToFineBoxCoarseUnits, neighbor, bcNeighbor,
            allocNeighbor);  // used for interpolation PF
    for (plint iA = 0; iA < (plint)coarseToFineEdgesCoarseUnitsTmp.size(); ++iA) {
        // dir is the direction along the edge
        plint dir = coarseToFineEdgesCoarseUnitsTmp[iA].getDirectionAlongEdge();

        std::pair<plint, plint> dirOne = getDirections(dir);

        coarseToFineEdgesCoarseUnits.push_back(boxLogic::DirectedEdge(
            coarseToFineEdgesCoarseUnitsTmp[iA].bb.enlarge(
                Array<plint, 2>(neighbor[dirOne.first], neighbor[dirOne.second]), dir),
            coarseToFineEdgesCoarseUnitsTmp[iA].planeDir, coarseToFineEdgesCoarseUnitsTmp[iA].dir1,
            coarseToFineEdgesCoarseUnitsTmp[iA].dir2));
    }

    std::vector<boxLogic::DirectedCorner> coarseToFineCornersCoarseUnitsTmp =
        boxLogic::getDirectedCorners(
            coarseToFineBoxCoarseUnits, neighbor, bcNeighbor,
            allocNeighbor);  // used for interpolation PF;

    coarseToFineCornersCoarseUnits.insert(
        coarseToFineCornersCoarseUnits.end(), coarseToFineCornersCoarseUnitsTmp.begin(),
        coarseToFineCornersCoarseUnitsTmp.end());

    std::vector<boxLogic::DirectedEdge> coarseToFineBoundaryEdgesCoarseUnitsTmp =
        boxLogic::getBoundaryDirectedEdges(
            coarseToFineBoxCoarseUnits, neighbor, bcNeighbor,
            allocNeighbor);  // used for special interpolation of edges.

    coarseToFineBoundaryEdgesCoarseUnits.insert(
        coarseToFineBoundaryEdgesCoarseUnits.end(), coarseToFineBoundaryEdgesCoarseUnitsTmp.begin(),
        coarseToFineBoundaryEdgesCoarseUnitsTmp.end());

    std::vector<boxLogic::DirectedCorner> coarseToFineBoundaryCornersCoarseUnitsTmp =
        boxLogic::getBoundaryDirectedCorners(
            coarseToFineBoxCoarseUnits, neighbor, bcNeighbor,
            allocNeighbor);  // used for special interpolation of corners.

    coarseToFineBoundaryCornersCoarseUnits.insert(
        coarseToFineBoundaryCornersCoarseUnits.end(),
        coarseToFineBoundaryCornersCoarseUnitsTmp.begin(),
        coarseToFineBoundaryCornersCoarseUnitsTmp.end());

    // =========== Fine to Coarse interfaces ================= //

    plint dx0 = ((!neighbor[OT::surface0N()] && allocNeighbor[OT::surface0N()]))
                - bcNeighbor[OT::surface0N()];
    plint dx1 = ((!neighbor[OT::surface0P()] && allocNeighbor[OT::surface0P()]))
                - bcNeighbor[OT::surface0P()];
    plint dy0 = ((!neighbor[OT::surface1N()] && allocNeighbor[OT::surface1N()]))
                - bcNeighbor[OT::surface1N()];
    plint dy1 = ((!neighbor[OT::surface1P()] && allocNeighbor[OT::surface1P()]))
                - bcNeighbor[OT::surface1P()];
    plint dz0 = ((!neighbor[OT::surface2N()] && allocNeighbor[OT::surface2N()]))
                - bcNeighbor[OT::surface2N()];
    plint dz1 = ((!neighbor[OT::surface2P()] && allocNeighbor[OT::surface2P()]))
                - bcNeighbor[OT::surface2P()];

    Array<plint, 6> enlargeAry(2 * dx0, -2 * dx1, 2 * dy0, -2 * dy1, 2 * dz0, -2 * dz1);
    // this box represents the interface Coarse to fine, in fine units
    Box3D fineToCoarseBoxFineUnitsTmp =
        coarseToFineBoxFineUnits.enlarge(overlapWidthCoarseUnits * enlargeAry);
    // this box represents the interface Coarse to fine, in coarse units
    Box3D fineToCoarseBoxCoarseUnits = fineToCoarseBoxFineUnitsTmp.divide(2);

    // in order to be generic the fine units box must be correct by 1 node
    // in he directions where there is a neighbor
    // Array<plint,6> correctAry(+dx0,-dx1, +dy0, -dy1, +dz0, -dz1);
    // Box3D fineToCoarseBoxFineUnits = fineToCoarseBoxFineUnitsTmp.enlarge(correctAry);
    Box3D fineToCoarseBoxFineUnits = fineToCoarseBoxFineUnitsTmp;

    // from fine to coarse one does not need directed planes/edges/corners
    // where the fine lattice is decomposed (must be enlarged because)
    // the nodes will be filtered.
    std::vector<boxLogic::DirectedPlane> fineToCoarsePlanesFineUnitsTmp =
        boxLogic::getDirectedPlanes(fineToCoarseBoxFineUnits, neighbor, bcNeighbor, allocNeighbor);
    for (plint iA = 0; iA < (plint)fineToCoarsePlanesFineUnitsTmp.size(); ++iA) {
        // the boxes must be enlarged in the direction parallel to the normal
        // of the plane by 1 and 2 in the others beucase the planes are only taking
        // into account the interior of a face of the Box3D.

        plint dir = fineToCoarsePlanesFineUnitsTmp[iA].direction;
        plint ori = fineToCoarsePlanesFineUnitsTmp[iA].orientation;
        std::pair<plint, plint> dirOne = getDirections((dir + 1) % 3);
        // Array<plint,2> nCellsOne(!neighbor[dirOne.first] && bcNeighbor[dirOne.first],
        // !neighbor[dirOne.second] && bcNeighbor[dirOne.second]);
        Array<plint, 2> nCellsOne(!neighbor[dirOne.first], !neighbor[dirOne.second]);

        std::pair<plint, plint> dirTwo = getDirections((dir + 2) % 3);
        // Array<plint,2> nCellsTwo(!neighbor[dirTwo.first] && bcNeighbor[dirTwo.first],
        // !neighbor[dirTwo.second] && bcNeighbor[dirTwo.second]);
        Array<plint, 2> nCellsTwo(!neighbor[dirTwo.first], !neighbor[dirTwo.second]);
        std::pair<Array<plint, 2>, Array<plint, 2> > nEdge =
            getEdgesIndices(dir, ori, neighbor, bcNeighbor, 5, 4);

        Array<plint, 6> bcAry(
            -!bcNeighbor[OT::surface0N()], !bcNeighbor[OT::surface0P()],
            -!bcNeighbor[OT::surface1N()], !bcNeighbor[OT::surface1P()],
            -!bcNeighbor[OT::surface2N()], !bcNeighbor[OT::surface2P()]);
        fineToCoarsePlanesFineUnits.push_back(
            fineToCoarsePlanesFineUnitsTmp[iA].bb.enlarge(bcAry).enlargeInNormalPlane(
                nCellsOne + nEdge.first, nCellsTwo + nEdge.second, dir));
    }

    std::vector<boxLogic::DirectedPlane> fineToCoarsePlanesCoarseUnitsTmp =
        boxLogic::getCompleteDirectedPlanes(
            fineToCoarseBoxCoarseUnits, neighbor, bcNeighbor, allocNeighbor);
    for (plint iA = 0; iA < (plint)fineToCoarsePlanesCoarseUnitsTmp.size(); ++iA) {
        plint dir = fineToCoarsePlanesCoarseUnitsTmp[iA].direction;
        plint ori = fineToCoarsePlanesCoarseUnitsTmp[iA].orientation;

        std::pair<Array<plint, 2>, Array<plint, 2> > nEdge =
            getEdgesIndices(dir, ori, neighbor, bcNeighbor, 2, 2);

        Array<plint, 6> bcAry(
            bcNeighbor[OT::surface0N()], -bcNeighbor[OT::surface0P()], bcNeighbor[OT::surface1N()],
            -bcNeighbor[OT::surface1P()], bcNeighbor[OT::surface2N()],
            -bcNeighbor[OT::surface2P()]);

        fineToCoarsePlanesCoarseUnits.push_back(
            fineToCoarsePlanesCoarseUnitsTmp[iA]
                .bb.enlargeInNormalPlane(nEdge.first, nEdge.second, dir)
                .enlarge(bcAry));
        // fineToCoarsePlanesCoarseUnits.push_back(
        //             fineToCoarsePlanesCoarseUnitsTmp[iA].bb );
    }

    std::vector<boxLogic::Edge> fineToCoarseBoundaryEdgesCoarseUnitsTmp =
        boxLogic::getBoundaryEdges(
            fineToCoarseBoxCoarseUnits, neighbor, bcNeighbor,
            allocNeighbor);  // used for only copy of boundary edges

    fineToCoarseBoundaryEdgesCoarseUnits.insert(
        fineToCoarseBoundaryEdgesCoarseUnits.end(), fineToCoarseBoundaryEdgesCoarseUnitsTmp.begin(),
        fineToCoarseBoundaryEdgesCoarseUnitsTmp.end());
}

}  // namespace boxLogic

}  // namespace plb
