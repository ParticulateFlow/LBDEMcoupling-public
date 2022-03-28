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

/** \file
 * Coupling between grids of different refinement level -- implementation.
 */

#include "gridRefinement/couplingInterfaceGenerator3D.h"

#include <map>
#include <vector>

#include "core/globalDefs.h"
#include "gridRefinement/boxLogic3D.hh"
#include "gridRefinement/couplingInterfaceGenerator3D.hh"
#include "gridRefinement/gridRefinementFunctional3D.h"
#include "gridRefinement/gridRefinementFunctional3D.hh"
#include "gridRefinement/octreeGridStructure.h"
#include "gridRefinement/rescaleEngine.h"
#include "io/imageWriter.h"
#include "offLattice/triangleSet.hh"
#include "offLattice/triangleSetGenerator.hh"

namespace plb {

CouplingInterfaces3D::CouplingInterfaces3D(
    OctreeGridStructure const &ogs_, int level_, int overlapWidth) :
    ogs(ogs_), level(level_)
{
    bool includeOverlaps = false;
    std::vector<plint> blockIdsPerLevel = ogs.getBlockIdsAtLevel(level, includeOverlaps);

    Box3D bulk;
    plint levelTmp = 0;
    plint processId = 0;
    // TriangleSet<double> *triangleSet = 0;
    // Array<plint,3> nSegments = Array<plint,3>(1, 1, 1);

    // triangleSet = new TriangleSet<double>(DBL);
    for (plint iA = 0; iA < (plint)blockIdsPerLevel.size(); ++iA) {
        ogs.getBlock(blockIdsPerLevel[iA], bulk, levelTmp, processId);

        const Array<bool, 26> neighbors = boxLogic::getNeighbors(ogs, blockIdsPerLevel[iA]);

        const Array<bool, 26> bcNeighbors = boxLogic::getBcNeighbors(ogs, blockIdsPerLevel[iA]);

        const Array<bool, 26> allocatedNeighbors =
            boxLogic::getAllocatedNeighbors(ogs, blockIdsPerLevel[iA]);

        boxLogic::getInterfaces(
            bulk, neighbors, bcNeighbors, allocatedNeighbors, overlapWidth,
            coarseToFinePlanesFineUnits, coarseToFinePlanesCoarseUnits,
            coarseToFineEdgesCoarseUnits, coarseToFineCornersCoarseUnits,
            coarseToFinePlanesCoarseUnitsExtended, fineToCoarsePlanesFineUnits,
            fineToCoarsePlanesCoarseUnits, coarseToFineBoundaryEdgesCoarseUnits,
            coarseToFineBoundaryCornersCoarseUnits, fineToCoarseBoundaryEdgesCoarseUnits);
    }

    coarseToFinePlanesFineUnits = boxLogic::merge(coarseToFinePlanesFineUnits);
    coarseToFinePlanesCoarseUnits = boxLogic::merge(coarseToFinePlanesCoarseUnits);
    coarseToFineEdgesCoarseUnits = boxLogic::merge(coarseToFineEdgesCoarseUnits);
    // coarseToFineCornersCoarseUnits         = boxLogic::merge(coarseToFineCornersCoarseUnits);
    coarseToFinePlanesCoarseUnitsExtended = boxLogic::merge(coarseToFinePlanesCoarseUnitsExtended);
    fineToCoarsePlanesFineUnits = boxLogic::merge(fineToCoarsePlanesFineUnits);
    fineToCoarsePlanesCoarseUnits = boxLogic::merge(fineToCoarsePlanesCoarseUnits);
    coarseToFineBoundaryEdgesCoarseUnits = boxLogic::merge(coarseToFineBoundaryEdgesCoarseUnits);
    coarseToFineBoundaryCornersCoarseUnits =
        boxLogic::merge(coarseToFineBoundaryCornersCoarseUnits);
    fineToCoarseBoundaryEdgesCoarseUnits = boxLogic::merge(fineToCoarseBoundaryEdgesCoarseUnits);
}

void CouplingInterfaces3D::writeInterfaces(double dx, const Array<double, 3> &pos) const
{
    TriangleSet<double> *triangleSet = 0;
    Array<plint, 3> nSegments = Array<plint, 3>(1, 1, 1);

    triangleSet = new TriangleSet<double>(DBL);

    Box3D bulk;
    plint levelTmp = 0;
    plint processId = 0;

    bool includeOverlaps = false;
    std::vector<plint> blockIdsPerLevel = ogs.getBlockIdsAtLevel(level, includeOverlaps);
    for (plint iA = 0; iA < (plint)blockIdsPerLevel.size(); ++iA) {
        ogs.getBlock(blockIdsPerLevel[iA], bulk, levelTmp, processId);

        TriangleSet<double> cuboidSet = constructCuboid<double>(
            Cuboid<double>(bulk.multiply(util::intTwoToThePower(ogs.getNumLevels() - level - 1))),
            nSegments);
        triangleSet->append(cuboidSet);
    }
    std::string fname =
        createFileName(global::directories().getOutputDir() + "bulk", level, 2) + ".stl";
    // pcout << fname << std::endl;
    triangleSet->scale(dx);
    triangleSet->translate(pos);
    triangleSet->writeBinarySTL(fname);
    delete triangleSet;

    // pcout << "writing int" << level << std::endl;
    triangleSet = new TriangleSet<double>(DBL);

    for (plint iA = 0; iA < (plint)coarseToFinePlanesFineUnits.size(); ++iA) {
        TriangleSet<double> cuboidSet = constructCuboid<double>(
            Cuboid<double>(coarseToFinePlanesFineUnits[iA].bb.multiply(
                util::intTwoToThePower(ogs.getNumLevels() - level - 1))),
            nSegments);
        triangleSet->append(cuboidSet);
    }
    fname = createFileName(
                global::directories().getOutputDir() + "coarseToFinePlanesFineUnits", level, 2)
            + ".stl";
    // pcout << fname << std::endl;
    triangleSet->scale(dx);
    triangleSet->translate(pos);
    triangleSet->writeBinarySTL(fname);
    delete triangleSet;

    triangleSet = new TriangleSet<double>(DBL);
    for (plint iA = 0; iA < (plint)coarseToFinePlanesCoarseUnits.size(); ++iA) {
        Array<plint, 2> nCells;
        nCells.resetToZero();
        nCells[0] = (coarseToFinePlanesCoarseUnits[iA].orientation == -1) ? 1 : 0;
        nCells[1] = (coarseToFinePlanesCoarseUnits[iA].orientation == +1) ? 1 : 0;
        TriangleSet<double> cuboidSet = constructCuboid<double>(
            Cuboid<double>(coarseToFinePlanesCoarseUnits[iA]
                               .bb.multiply(util::intTwoToThePower(ogs.getNumLevels() - level))
                               .enlarge(nCells, coarseToFinePlanesCoarseUnits[iA].direction)),
            nSegments);
        triangleSet->append(cuboidSet);
    }
    fname = createFileName(
                global::directories().getOutputDir() + "coarseToFinePlanesCoarseUnits", level, 2)
            + ".stl";
    triangleSet->scale(dx);
    triangleSet->translate(pos);
    triangleSet->writeBinarySTL(fname);
    delete triangleSet;

    triangleSet = new TriangleSet<double>(DBL);
    for (plint iA = 0; iA < (plint)coarseToFinePlanesCoarseUnitsExtended.size(); ++iA) {
        TriangleSet<double> cuboidSet = constructCuboid<double>(
            Cuboid<double>(coarseToFinePlanesCoarseUnitsExtended[iA].bb.multiply(
                util::intTwoToThePower(ogs.getNumLevels() - level))),
            nSegments);
        triangleSet->append(cuboidSet);
    }
    fname = createFileName(
                global::directories().getOutputDir() + "coarseToFinePlanesCoarseUnitsExtended",
                level, 2)
            + ".stl";
    triangleSet->scale(dx);
    triangleSet->translate(pos);
    triangleSet->writeBinarySTL(fname);
    delete triangleSet;

    triangleSet = new TriangleSet<double>(DBL);
    for (plint iA = 0; iA < (plint)fineToCoarsePlanesFineUnits.size(); ++iA) {
        TriangleSet<double> cuboidSet = constructCuboid<double>(
            Cuboid<double>(fineToCoarsePlanesFineUnits[iA].multiply(
                util::intTwoToThePower(ogs.getNumLevels() - level - 1))),
            nSegments);
        triangleSet->append(cuboidSet);
    }
    fname = createFileName(
                global::directories().getOutputDir() + "fineToCoarsePlanesFineUnits", level, 2)
            + ".stl";
    triangleSet->scale(dx);
    triangleSet->translate(pos);
    triangleSet->writeBinarySTL(fname);
    delete triangleSet;

    triangleSet = new TriangleSet<double>(DBL);
    for (plint iA = 0; iA < (plint)fineToCoarsePlanesCoarseUnits.size(); ++iA) {
        TriangleSet<double> cuboidSet = constructCuboid<double>(
            Cuboid<double>(fineToCoarsePlanesCoarseUnits[iA].bb.multiply(
                util::intTwoToThePower(ogs.getNumLevels() - level))),
            nSegments);
        triangleSet->append(cuboidSet);
    }
    fname = createFileName(
                global::directories().getOutputDir() + "fineToCoarsePlanesCoarseUnits", level, 2)
            + ".stl";
    triangleSet->scale(dx);
    triangleSet->translate(pos);
    triangleSet->writeBinarySTL(fname);
    delete triangleSet;

    // pcout << dx << " " << pos[0] << ", " << pos[1] << ", " << pos[2] << std::endl;

    triangleSet = new TriangleSet<double>(DBL);
    for (plint iA = 0; iA < (plint)coarseToFineEdgesCoarseUnits.size(); ++iA) {
        Array<plint, 6> ary;
        ary.resetToZero();
        plint planeDir = coarseToFineEdgesCoarseUnits[iA].planeDir;
        plint dir1 = coarseToFineEdgesCoarseUnits[iA].dir1;
        plint dir2 = coarseToFineEdgesCoarseUnits[iA].dir2;

        if (planeDir == 0) {
            if (dir1 == 0) {
                if (dir2 == 1) {
                    ary[5] = -1;
                } else if (dir2 == -1) {
                    ary[4] = 1;
                }
            } else if (dir1 == 1) {
                ary[3] = -1;
            } else if (dir1 == -1) {
                ary[2] = 1;
            }
        } else if (planeDir == 1) {
            if (dir1 == 0) {
                if (dir2 == 1) {
                    ary[1] = -1;
                } else if (dir2 == -1) {
                    ary[0] = 1;
                }
            } else if (dir1 == 1) {
                ary[5] = -1;
            } else if (dir1 == -1) {
                ary[4] = 1;
            }
        } else if (planeDir == 2) {
            if (dir1 == 0) {
                if (dir2 == 1) {
                    ary[3] = -1;
                } else if (dir2 == -1) {
                    ary[2] = 1;
                }
            } else if (dir1 == 1) {
                ary[1] = -1;
            } else if (dir1 == -1) {
                ary[0] = 1;
            }
        }

        TriangleSet<double> cuboidSet = constructCuboid<double>(
            Cuboid<double>(coarseToFineEdgesCoarseUnits[iA]
                               .bb.multiply(util::intTwoToThePower(ogs.getNumLevels() - level))
                               .enlarge(ary)),
            nSegments);
        triangleSet->append(cuboidSet);
    }
    fname = createFileName(
                global::directories().getOutputDir() + "coarseToFineEdgesCoarseUnits", level, 2)
            + ".stl";
    triangleSet->scale(dx);
    triangleSet->translate(pos);
    triangleSet->writeBinarySTL(fname);
    delete triangleSet;

    triangleSet = new TriangleSet<double>(DBL);
    for (plint iA = 0; iA < (plint)coarseToFineCornersCoarseUnits.size(); ++iA) {
        Array<plint, 6> ary;
        ary.resetToZero();
        plint planeDir = coarseToFineCornersCoarseUnits[iA].planeDir;
        plint dir1 = coarseToFineCornersCoarseUnits[iA].dir1;
        plint dir2 = coarseToFineCornersCoarseUnits[iA].dir2;

        if (planeDir == 0) {
            if (dir2 == 1) {
                ary[5] = -dir2;
            } else if (dir2 == -1) {
                ary[4] = -dir2;
            } else {
                PLB_ASSERT(false && "dir2 must be different from zero");
            }

            if (dir1 == 1) {
                ary[3] = -dir1;
            } else if (dir1 == -1) {
                ary[2] = -dir1;
            } else {
                PLB_ASSERT(false && "dir1 must be different from zero");
            }

        } else if (planeDir == 1) {
            if (dir2 == 1) {
                ary[1] = -dir2;
            } else if (dir2 == -1) {
                ary[0] = -dir2;
            } else {
                PLB_ASSERT(false && "dir2 must be different from zero");
            }

            if (dir1 == 1) {
                ary[5] = -dir1;
            } else if (dir1 == -1) {
                ary[4] = -dir1;
            } else {
                PLB_ASSERT(false && "dir1 must be different from zero");
            }
        } else if (planeDir == 2) {
            if (dir2 == 1) {
                ary[3] = -dir2;
            } else if (dir2 == -1) {
                ary[2] = -dir2;
            }

            if (dir1 == 1) {
                ary[1] = -dir1;
            } else if (dir1 == -1) {
                ary[0] = -dir1;
            }
        }

        TriangleSet<double> cuboidSet = constructCuboid<double>(
            Cuboid<double>(coarseToFineCornersCoarseUnits[iA]
                               .bb.multiply(util::intTwoToThePower(ogs.getNumLevels() - level))
                               .enlarge(ary)),
            nSegments);
        triangleSet->append(cuboidSet);
    }

    fname = createFileName(
                global::directories().getOutputDir() + "coarseToFineCornersCoarseUnits", level, 2)
            + ".stl";
    triangleSet->scale(dx);
    triangleSet->translate(pos);
    triangleSet->writeBinarySTL(fname);
    delete triangleSet;

    triangleSet = new TriangleSet<double>(DBL);
    for (plint iA = 0; iA < (plint)coarseToFineBoundaryEdgesCoarseUnits.size(); ++iA) {
        Array<plint, 6> ary;
        ary.resetToZero();
        plint planeDir = coarseToFineBoundaryEdgesCoarseUnits[iA].planeDir;
        plint dir1 = coarseToFineBoundaryEdgesCoarseUnits[iA].dir1;
        plint dir2 = coarseToFineBoundaryEdgesCoarseUnits[iA].dir2;

        if (planeDir == 0) {
            if (dir1 == 0) {
                if (dir2 == 1) {
                    ary[5] = -1;
                } else if (dir2 == -1) {
                    ary[4] = 1;
                }
            } else if (dir1 == 1) {
                ary[3] = -1;
            } else if (dir1 == -1) {
                ary[2] = 1;
            }
        } else if (planeDir == 1) {
            if (dir1 == 0) {
                if (dir2 == 1) {
                    ary[1] = -1;
                } else if (dir2 == -1) {
                    ary[0] = 1;
                }
            } else if (dir1 == 1) {
                ary[5] = -1;
            } else if (dir1 == -1) {
                ary[4] = 1;
            }
        } else if (planeDir == 2) {
            if (dir1 == 0) {
                if (dir2 == 1) {
                    ary[3] = -1;
                } else if (dir2 == -1) {
                    ary[2] = 1;
                }
            } else if (dir1 == 1) {
                ary[1] = -1;
            } else if (dir1 == -1) {
                ary[0] = 1;
            }
        }

        TriangleSet<double> cuboidSet = constructCuboid<double>(
            Cuboid<double>(coarseToFineBoundaryEdgesCoarseUnits[iA]
                               .bb.multiply(util::intTwoToThePower(ogs.getNumLevels() - level))
                               .enlarge(Array<plint, 6>(0, 1, 0, 1, 0, 1))),
            nSegments);
        triangleSet->append(cuboidSet);
    }

    fname =
        createFileName(
            global::directories().getOutputDir() + "coarseToFineBoundaryEdgesCoarseUnits", level, 2)
        + ".stl";
    triangleSet->scale(dx);
    triangleSet->translate(pos);
    triangleSet->writeBinarySTL(fname);
    delete triangleSet;

    triangleSet = new TriangleSet<double>(DBL);
    for (plint iA = 0; iA < (plint)coarseToFineBoundaryCornersCoarseUnits.size(); ++iA) {
        Array<plint, 6> ary;
        ary.resetToZero();
        plint planeDir = coarseToFineBoundaryCornersCoarseUnits[iA].planeDir;
        plint dir1 = coarseToFineBoundaryCornersCoarseUnits[iA].dir1;
        plint dir2 = coarseToFineBoundaryCornersCoarseUnits[iA].dir2;

        if (planeDir == 0) {
            if (dir2 == 1) {
                ary[5] = -dir2;
            } else if (dir2 == -1) {
                ary[4] = -dir2;
            }

            if (dir1 == 1) {
                ary[3] = -dir1;
            } else if (dir1 == -1) {
                ary[2] = -dir1;
            }
        } else if (planeDir == 1) {
            if (dir2 == 1) {
                ary[1] = -dir2;
            } else if (dir2 == -1) {
                ary[0] = -dir2;
            }

            if (dir1 == 1) {
                ary[5] = -dir1;
            } else if (dir1 == -1) {
                ary[4] = -dir1;
            }
        } else if (planeDir == 2) {
            if (dir2 == 1) {
                ary[3] = -dir2;
            } else if (dir2 == -1) {
                ary[2] = -dir2;
            }

            if (dir1 == 1) {
                ary[1] = -dir1;
            } else if (dir1 == -1) {
                ary[0] = -dir1;
            }
        }

        TriangleSet<double> cuboidSet = constructCuboid<double>(
            Cuboid<double>(coarseToFineBoundaryCornersCoarseUnits[iA]
                               .bb.multiply(util::intTwoToThePower(ogs.getNumLevels() - level))
                               .enlarge(ary)),
            nSegments);
        triangleSet->append(cuboidSet);
    }

    fname = createFileName(
                global::directories().getOutputDir() + "coarseToFineBoundaryCornersCoarseUnits",
                level, 2)
            + ".stl";
    triangleSet->scale(dx);
    triangleSet->translate(pos);
    triangleSet->writeBinarySTL(fname);
    delete triangleSet;

    triangleSet = new TriangleSet<double>(DBL);
    for (plint iA = 0; iA < (plint)fineToCoarseBoundaryEdgesCoarseUnits.size(); ++iA) {
        TriangleSet<double> cuboidSet = constructCuboid<double>(
            Cuboid<double>(fineToCoarseBoundaryEdgesCoarseUnits[iA]
                               .bb.multiply(util::intTwoToThePower(ogs.getNumLevels() - level))
                               .enlarge(Array<plint, 6>(0, 1, 0, 1, 0, 1))),
            nSegments);
        triangleSet->append(cuboidSet);
    }

    fname =
        createFileName(
            global::directories().getOutputDir() + "fineToCoarseBoundaryEdgesCoarseUnits", level, 2)
        + ".stl";
    triangleSet->scale(dx);
    triangleSet->translate(pos);
    triangleSet->writeBinarySTL(fname);
    delete triangleSet;
}

}  // namespace plb
