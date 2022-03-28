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
 * Management of multi grid domains -- Implementation
 */

#include "multiGrid/multiGridManagement3D.h"

#include <memory>

#include "io/parallelIO.h"
#include "multiBlock/multiBlockGenerator3D.h"
#include "multiBlock/multiBlockGenerator3D.hh"
#include "multiBlock/multiBlockManagement3D.h"
#include "multiBlock/multiContainerBlock3D.h"
#include "multiBlock/multiDataField3D.hh"
#include "multiBlock/nonLocalTransfer3D.h"
#include "multiBlock/nonLocalTransfer3D.hh"
#include "multiBlock/reductiveMultiDataProcessorWrapper3D.h"
#include "offLattice/makeSparse3D.hh"
#include "parallelism/parallelMultiDataField3D.h"
#include "parallelism/parallelMultiDataField3D.hh"

namespace plb {

enum { left = 0, right, bottom, top, front, back };

MultiGridManagement3D::MultiGridManagement3D(
    plint nx, plint ny, plint nz, plint numLevels, plint referenceLevel_) :
    referenceLevel(referenceLevel_),
    boundingBoxes(numLevels),
    coarseGridInterfaces(numLevels),
    fineGridInterfaces(numLevels),
    bulks(numLevels)
{
    PLB_ASSERT(numLevels > 0);
    PLB_ASSERT(referenceLevel < numLevels);
    // mpiProcess.reserve(numLevels);
    initialize(Box3D(0, nx - 1, 0, ny - 1, 0, nz - 1));
}

MultiGridManagement3D::MultiGridManagement3D(
    Box3D coarseBoundingBox, plint numLevels, plint referenceLevel_) :
    referenceLevel(referenceLevel_),
    boundingBoxes(numLevels),
    coarseGridInterfaces(numLevels),
    fineGridInterfaces(numLevels),
    bulks(numLevels)
{
    PLB_ASSERT(numLevels > 0);
    PLB_ASSERT(referenceLevel < numLevels);
    // mpiProcess.reserve(numLevels);
    initialize(coarseBoundingBox);
}

MultiGridManagement3D::~MultiGridManagement3D()
{
    delete scaleManager;
}

void MultiGridManagement3D::initialize(Box3D const &level0_box)
{
    scaleManager = global::getDefaultMultiScaleManager().clone();

    // The parameters nx and ny refer to the coarsest lattice (lattice 0).
    boundingBoxes[0] = level0_box;
    // Multiply the values maxX and maxY by two to get, iteratively,
    //   the boundingBoxes of the next finer lattice respectively.
    for (pluint iDim = 1; iDim < boundingBoxes.size(); ++iDim) {
        boundingBoxes[iDim] = scaleManager->scaleBox(boundingBoxes[iDim - 1], 1);
    }

    // Add a block for the full domain of the first instantiated multi-block.
    bulks[referenceLevel].push_back(boundingBoxes[referenceLevel]);
}

MultiGridManagement3D::MultiGridManagement3D(MultiGridManagement3D const &rhs) :
    referenceLevel(rhs.referenceLevel),
    boundingBoxes(rhs.boundingBoxes),
    mpiProcess(rhs.mpiProcess),
    coarseGridInterfaces(rhs.coarseGridInterfaces),
    fineGridInterfaces(rhs.fineGridInterfaces),
    bulks(rhs.bulks),
    scaleManager(rhs.scaleManager->clone())
{ }

void MultiGridManagement3D::swap(MultiGridManagement3D &rhs)
{
    std::swap(referenceLevel, rhs.referenceLevel);
    boundingBoxes.swap(rhs.boundingBoxes);
    mpiProcess.swap(rhs.mpiProcess);
    coarseGridInterfaces.swap(rhs.coarseGridInterfaces);
    fineGridInterfaces.swap(rhs.fineGridInterfaces);
    bulks.swap(rhs.bulks);
    std::swap(scaleManager, rhs.scaleManager);
}

Box3D MultiGridManagement3D::getBoundingBox(plint level) const
{
    PLB_PRECONDITION(level >= 0 && level < (plint)bulks.size());
    return boundingBoxes[level];
}

std::vector<Box3D> const &MultiGridManagement3D::getBulks(plint iLevel) const
{
    PLB_PRECONDITION(iLevel >= 0 && iLevel < (plint)bulks.size());
    return bulks[iLevel];
}

std::vector<std::vector<Box3D> > const &MultiGridManagement3D::getBulks() const
{
    return bulks;
}

void MultiGridManagement3D::refine(plint coarseLevel, Box3D coarseDomain)
{
    // The finest multi-block, at level blocks.size()-1, cannot be further refined.
    PLB_PRECONDITION(coarseLevel >= 0 && coarseLevel < (plint)bulks.size() - 1);
    plint fineLevel = coarseLevel + 1;

    // container for booleans that indicate if coarseDomain touches the six walls of the
    // domain. The convention is the following: 0=left,1=right,2=bottom,3=top,4=back,5=front
    bool touches[6] = {false, false, false, false, false, false};
    trimDomain(coarseLevel, coarseDomain, touches);

    // Same as coarseDomain, but in units of the fine lattice.
    Box3D fineDomain(coarseDomain.multiply(2));

    /************ INCLUDED BY HELEN ******************************************/
    // Test whether the reduced fine domain touches the coarse domain to avoid
    // creating wrong interfaces
    bool touchesFineDomain[6] = {false, false, false, false, false, false};

    /// LEFT
    if (!touches[left]) {
        for (pluint iBlock = 0; iBlock < bulks[fineLevel].size(); ++iBlock) {
            // find bounding box of current fineLevel blocks
            if (bulks[fineLevel][iBlock].x1 - 2 == fineDomain.x0) {
                touchesFineDomain[left] = true;
            }
        }
    }
    /// BOTTOM
    if (!touches[bottom]) {
        for (pluint iBlock = 0; iBlock < bulks[fineLevel].size(); ++iBlock) {
            // find bounding box of current fineLevel blocks
            if (bulks[fineLevel][iBlock].y1 - 2 == fineDomain.y0) {
                touchesFineDomain[bottom] = true;
            }
        }
    }
    /// FRONT
    if (!touches[front]) {
        for (pluint iBlock = 0; iBlock < bulks[fineLevel].size(); ++iBlock) {
            // find bounding box of current fineLevel blocks
            if (bulks[fineLevel][iBlock].z1 - 2 == fineDomain.z0) {
                touchesFineDomain[front] = true;
            }
        }
    }
    /// RIGHT
    if (!touches[right]) {
        for (pluint iBlock = 0; iBlock < bulks[fineLevel].size(); ++iBlock) {
            // find bounding box of current fineLevel blocks
            if (bulks[fineLevel][iBlock].x0 + 2 == fineDomain.x1) {
                touchesFineDomain[right] = true;
            }
        }
    }
    /// TOP
    if (!touches[top]) {
        for (pluint iBlock = 0; iBlock < bulks[fineLevel].size(); ++iBlock) {
            // find bounding box of current fineLevel blocks
            if (bulks[fineLevel][iBlock].y0 + 2 == fineDomain.y1) {
                touchesFineDomain[top] = true;
            }
        }
    }
    /// BACK
    if (!touches[back]) {
        for (pluint iBlock = 0; iBlock < bulks[fineLevel].size(); ++iBlock) {
            // find bounding box of current fineLevel blocks
            if (bulks[fineLevel][iBlock].z0 + 2 == fineDomain.z1) {
                touchesFineDomain[back] = true;
            }
        }
    }
    /*******************************************************************/

    // The reduced coarse domain is the one which is going to be excluded from
    //   the original coarse lattice.
    Box3D reducedCoarseDomain(coarseDomain.enlarge(-1));
    // The extended coarse domain it the one which is going to be added
    //   to the original fine lattice.
    Box3D extendedCoarseDomain(coarseDomain.enlarge(1));

    // If the domain in question touches a boundary of the multi-block,
    //   both the reduced and the extended coarse domain are
    //   identified with the boundary location.
    if (touches[left]) {  // left
        reducedCoarseDomain.x0 -= 1;
        extendedCoarseDomain.x0 += 1;
    }
    if (touches[right]) {  // right
        reducedCoarseDomain.x1 += 1;
        extendedCoarseDomain.x1 -= 1;
    }
    if (touches[bottom]) {  // bottom
        reducedCoarseDomain.y0 -= 1;
        extendedCoarseDomain.y0 += 1;
    }
    if (touches[top]) {  // top
        reducedCoarseDomain.y1 += 1;
        extendedCoarseDomain.y1 -= 1;
    }
    if (touches[front]) {  // front
        reducedCoarseDomain.z0 -= 1;
        extendedCoarseDomain.z0 += 1;
    }
    if (touches[back]) {  // back
        reducedCoarseDomain.z1 += 1;
        extendedCoarseDomain.z1 -= 1;
    }

    if (touchesFineDomain[left]) {  // left
        reducedCoarseDomain.x0 -= 2;
        extendedCoarseDomain.x0 += 2;
    }
    if (touchesFineDomain[right]) {  // right
        reducedCoarseDomain.x1 += 2;
        extendedCoarseDomain.x1 -= 2;
    }
    if (touchesFineDomain[bottom]) {  // bottom
        reducedCoarseDomain.y0 -= 2;
        extendedCoarseDomain.y0 += 2;
    }
    if (touchesFineDomain[top]) {  // top
        reducedCoarseDomain.y1 += 2;
        extendedCoarseDomain.y1 -= 2;
    }
    if (touchesFineDomain[front]) {  // front
        reducedCoarseDomain.z0 -= 2;
        extendedCoarseDomain.z0 += 2;
    }
    if (touchesFineDomain[back]) {  // back
        reducedCoarseDomain.z1 += 2;
        extendedCoarseDomain.z1 -= 2;
    }

    // Extract reduced coarse domain from the original coarse multi-block, for the bulks.
    std::vector<Box3D> exceptedBulks;
    for (pluint iBlock = 0; iBlock < bulks[coarseLevel].size(); ++iBlock) {
        except(bulks[coarseLevel][iBlock], reducedCoarseDomain, exceptedBulks);
    }
    exceptedBulks.swap(bulks[coarseLevel]);

    // Convert the extended coarse domain to fine units,
    //   and add to the original fine multi-block.
    Box3D extendedFineDomain(extendedCoarseDomain.multiply(2));
    bulks[fineLevel].push_back(extendedFineDomain);

    // Define coupling interfaces for all six sides of the coarsened domain, unless they
    //   touch a boundary of the multi-block.
    std::vector<Box3D> BoxesToExceptCoarse, BoxesToExceptFine;

    if (!touches[left]) {  /// LEFT
        if (!touchesFineDomain[left]) {
            coarseGridInterfaces[coarseLevel].push_back(Box3D(
                coarseDomain.x0, coarseDomain.x0, coarseDomain.y0, coarseDomain.y1, coarseDomain.z0,
                coarseDomain.z1));
            fineGridInterfaces[fineLevel].push_back(Box3D(
                extendedCoarseDomain.x0, extendedCoarseDomain.x0, extendedCoarseDomain.y0,
                extendedCoarseDomain.y1, extendedCoarseDomain.z0, extendedCoarseDomain.z1));
        } else {
            BoxesToExceptCoarse.push_back(Box3D(
                coarseDomain.x0, coarseDomain.x0, coarseDomain.y0 + 1, coarseDomain.y1 - 1,
                coarseDomain.z0 + 1, coarseDomain.z1 - 1));
            BoxesToExceptFine.push_back(Box3D(
                extendedCoarseDomain.x0, extendedCoarseDomain.x0, extendedCoarseDomain.y0 + 1,
                extendedCoarseDomain.y1 - 1, extendedCoarseDomain.z0 + 1,
                extendedCoarseDomain.z1 - 1));
        }
    }
    if (!touches[right]) {  /// RIGHT
        if (!touchesFineDomain[right]) {
            coarseGridInterfaces[coarseLevel].push_back(Box3D(
                coarseDomain.x1, coarseDomain.x1, coarseDomain.y0, coarseDomain.y1, coarseDomain.z0,
                coarseDomain.z1));
            fineGridInterfaces[fineLevel].push_back(Box3D(
                extendedCoarseDomain.x1, extendedCoarseDomain.x1, extendedCoarseDomain.y0,
                extendedCoarseDomain.y1, extendedCoarseDomain.z0, extendedCoarseDomain.z1));
        } else {
            BoxesToExceptCoarse.push_back(Box3D(
                coarseDomain.x1, coarseDomain.x1, coarseDomain.y0 + 1, coarseDomain.y1 - 1,
                coarseDomain.z0 + 1, coarseDomain.z1 - 1));
            BoxesToExceptFine.push_back(Box3D(
                extendedCoarseDomain.x1, extendedCoarseDomain.x1, extendedCoarseDomain.y0 + 1,
                extendedCoarseDomain.y1 - 1, extendedCoarseDomain.z0 + 1,
                extendedCoarseDomain.z1 - 1));
        }
    }
    if (!touches[bottom]) {  /// BOTTOM
        if (!touchesFineDomain[bottom]) {
            coarseGridInterfaces[coarseLevel].push_back(Box3D(
                coarseDomain.x0, coarseDomain.x1, coarseDomain.y0, coarseDomain.y0, coarseDomain.z0,
                coarseDomain.z1));
            fineGridInterfaces[fineLevel].push_back(Box3D(
                extendedCoarseDomain.x0, extendedCoarseDomain.x1, extendedCoarseDomain.y0,
                extendedCoarseDomain.y0, extendedCoarseDomain.z0, extendedCoarseDomain.z1));
        } else {
            BoxesToExceptCoarse.push_back(Box3D(
                coarseDomain.x0 + 1, coarseDomain.x1 - 1, coarseDomain.y0, coarseDomain.y0,
                coarseDomain.z0 + 1, coarseDomain.z1 - 1));
            BoxesToExceptFine.push_back(Box3D(
                extendedCoarseDomain.x0 + 1, extendedCoarseDomain.x1 - 1, extendedCoarseDomain.y0,
                extendedCoarseDomain.y0, extendedCoarseDomain.z0 + 1, extendedCoarseDomain.z1 - 1));
        }
    }
    if (!touches[top]) {  /// TOP
        if (!touchesFineDomain[top]) {
            coarseGridInterfaces[coarseLevel].push_back(Box3D(
                coarseDomain.x0, coarseDomain.x1, coarseDomain.y1, coarseDomain.y1, coarseDomain.z0,
                coarseDomain.z1));
            fineGridInterfaces[fineLevel].push_back(Box3D(
                extendedCoarseDomain.x0, extendedCoarseDomain.x1, extendedCoarseDomain.y1,
                extendedCoarseDomain.y1, extendedCoarseDomain.z0, extendedCoarseDomain.z1));
        } else {
            BoxesToExceptCoarse.push_back(Box3D(
                coarseDomain.x0 + 1, coarseDomain.x1 - 1, coarseDomain.y1, coarseDomain.y1,
                coarseDomain.z0 + 1, coarseDomain.z1 - 1));
            BoxesToExceptFine.push_back(Box3D(
                extendedCoarseDomain.x0 + 1, extendedCoarseDomain.x1 - 1, extendedCoarseDomain.y1,
                extendedCoarseDomain.y1, extendedCoarseDomain.z0 + 1, extendedCoarseDomain.z1 - 1));
        }
    }
    if (!touches[front]) {  /// FRONT
        if (!touchesFineDomain[front]) {
            coarseGridInterfaces[coarseLevel].push_back(Box3D(
                coarseDomain.x0, coarseDomain.x1, coarseDomain.y0, coarseDomain.y1, coarseDomain.z0,
                coarseDomain.z0));
            fineGridInterfaces[fineLevel].push_back(Box3D(
                extendedCoarseDomain.x0, extendedCoarseDomain.x1, extendedCoarseDomain.y0,
                extendedCoarseDomain.y1, extendedCoarseDomain.z0, extendedCoarseDomain.z0));
        } else {
            BoxesToExceptCoarse.push_back(Box3D(
                coarseDomain.x0 + 1, coarseDomain.x1 - 1, coarseDomain.y0 + 1, coarseDomain.y1 - 1,
                coarseDomain.z0, coarseDomain.z0));
            BoxesToExceptFine.push_back(Box3D(
                extendedCoarseDomain.x0 + 1, extendedCoarseDomain.x1 - 1,
                extendedCoarseDomain.y0 + 1, extendedCoarseDomain.y1 - 1, extendedCoarseDomain.z0,
                extendedCoarseDomain.z0));
        }
    }
    if (!touches[back]) {  /// BACK
        if (!touchesFineDomain[back]) {
            coarseGridInterfaces[coarseLevel].push_back(Box3D(
                coarseDomain.x0, coarseDomain.x1, coarseDomain.y0, coarseDomain.y1, coarseDomain.z1,
                coarseDomain.z1));
            fineGridInterfaces[fineLevel].push_back(Box3D(
                extendedCoarseDomain.x0, extendedCoarseDomain.x1, extendedCoarseDomain.y0,
                extendedCoarseDomain.y1, extendedCoarseDomain.z1, extendedCoarseDomain.z1));
        } else {
            BoxesToExceptCoarse.push_back(Box3D(
                coarseDomain.x0 + 1, coarseDomain.x1 - 1, coarseDomain.y0 + 1, coarseDomain.y1 - 1,
                coarseDomain.z1, coarseDomain.z1));
            BoxesToExceptFine.push_back(Box3D(
                extendedCoarseDomain.x0 + 1, extendedCoarseDomain.x1 - 1,
                extendedCoarseDomain.y0 + 1, extendedCoarseDomain.y1 - 1, extendedCoarseDomain.z1,
                extendedCoarseDomain.z1));
        }
    }

    // Extract the boxes to except from the corresponding interfaces
    std::vector<Box3D> fineInterfaces;
    Box3D box, boxForExcept;
    std::vector<Box3D> newInterfaces;
    bool boxFound = false;
    for (pluint iComp = 0; iComp < fineGridInterfaces[fineLevel].size(); ++iComp) {
        box = fineGridInterfaces[fineLevel][iComp];
        // interface is left or right
        if (box.x0 == box.x1) {
            if (!touchesFineDomain[left] && !touchesFineDomain[right]) {
                fineInterfaces.push_back(box);
            }
            if (touchesFineDomain[left]) {
                // find box to except from vector
                for (pluint iBox = 0; iBox < BoxesToExceptFine.size(); ++iBox) {
                    if ((BoxesToExceptFine[iBox].x0 == box.x0)
                        && (BoxesToExceptFine[iBox].x1 == box.x0)) {
                        boxForExcept = BoxesToExceptFine[iBox];
                        boxFound = true;
                    }
                }
                if (boxFound) {
                    except(box, boxForExcept, newInterfaces);
                    for (pluint i = 0; i < newInterfaces.size(); ++i) {
                        fineInterfaces.push_back(newInterfaces[i]);
                    }
                    newInterfaces.clear();
                    boxFound = false;
                } else {
                    fineInterfaces.push_back(box);
                }
            }
            if (touchesFineDomain[right]) {
                // find box to except from vector
                for (pluint iBox = 0; iBox < BoxesToExceptFine.size(); ++iBox) {
                    if ((BoxesToExceptFine[iBox].x0 == box.x0)
                        && (BoxesToExceptFine[iBox].x1 == box.x0)) {
                        boxForExcept = BoxesToExceptFine[iBox];
                        boxFound = true;
                    }
                }
                if (boxFound) {
                    except(box, boxForExcept, newInterfaces);
                    for (pluint i = 0; i < newInterfaces.size(); ++i) {
                        fineInterfaces.push_back(newInterfaces[i]);
                    }
                    newInterfaces.clear();
                } else {
                    fineInterfaces.push_back(box);
                }
                boxFound = false;
            }
        } else if (box.y0 == box.y1) {  // interface is top or bottom
            if (!touchesFineDomain[top] && !touchesFineDomain[bottom]) {
                fineInterfaces.push_back(box);
            }
            if (touchesFineDomain[top]) {
                // find box to except from vector
                for (pluint iBox = 0; iBox < BoxesToExceptFine.size(); ++iBox) {
                    if ((BoxesToExceptFine[iBox].y0 == box.y0)
                        && (BoxesToExceptFine[iBox].y1 == box.y0)) {
                        boxForExcept = BoxesToExceptFine[iBox];
                        boxFound = true;
                    }
                }
                if (boxFound) {
                    except(box, boxForExcept, newInterfaces);
                    for (pluint i = 0; i < newInterfaces.size(); ++i) {
                        fineInterfaces.push_back(newInterfaces[i]);
                    }
                    newInterfaces.clear();
                    boxFound = false;
                } else {
                    fineInterfaces.push_back(box);
                }
            }
            if (touchesFineDomain[bottom]) {
                // find box to except from vector
                for (pluint iBox = 0; iBox < BoxesToExceptFine.size(); ++iBox) {
                    if ((BoxesToExceptFine[iBox].y0 == box.y0)
                        && (BoxesToExceptFine[iBox].y1 == box.y0)) {
                        boxForExcept = BoxesToExceptFine[iBox];
                        boxFound = true;
                    }
                }
                if (boxFound) {
                    except(box, boxForExcept, newInterfaces);
                    for (pluint i = 0; i < newInterfaces.size(); ++i) {
                        fineInterfaces.push_back(newInterfaces[i]);
                    }
                    newInterfaces.clear();
                } else {
                    fineInterfaces.push_back(box);
                }
                boxFound = false;
            }
        } else {  // interface is front or back (box.z0 == box.z1)
            if (!touchesFineDomain[front] && !touchesFineDomain[back]) {
                fineInterfaces.push_back(box);
            }
            if (touchesFineDomain[front]) {
                // find box to except from vector
                for (pluint iBox = 0; iBox < BoxesToExceptFine.size(); ++iBox) {
                    if ((BoxesToExceptFine[iBox].z0 == box.z0)
                        && (BoxesToExceptFine[iBox].z1 == box.z0)) {
                        boxForExcept = BoxesToExceptFine[iBox];
                        boxFound = true;
                    }
                }
                if (boxFound) {
                    except(box, boxForExcept, newInterfaces);
                    for (pluint i = 0; i < newInterfaces.size(); ++i) {
                        fineInterfaces.push_back(newInterfaces[i]);
                    }
                    newInterfaces.clear();
                    boxFound = false;
                } else {
                    fineInterfaces.push_back(box);
                }
            }
            if (touchesFineDomain[back]) {
                // find box to except from vector
                for (pluint iBox = 0; iBox < BoxesToExceptFine.size(); ++iBox) {
                    if ((BoxesToExceptFine[iBox].z0 == box.z0)
                        && (BoxesToExceptFine[iBox].z1 == box.z0)) {
                        boxForExcept = BoxesToExceptFine[iBox];
                        boxFound = true;
                    }
                }
                if (boxFound) {
                    except(box, boxForExcept, newInterfaces);
                    for (pluint i = 0; i < newInterfaces.size(); ++i) {
                        fineInterfaces.push_back(newInterfaces[i]);
                    }
                    newInterfaces.clear();
                    boxFound = false;
                } else {
                    fineInterfaces.push_back(box);
                }
            }
        }
    }
    fineGridInterfaces[fineLevel] = fineInterfaces;

    std::vector<Box3D> coarseInterfaces;
    newInterfaces.clear();
    for (pluint iComp = 0; iComp < coarseGridInterfaces[coarseLevel].size(); ++iComp) {
        box = coarseGridInterfaces[coarseLevel][iComp];
        // interface is left or right
        if (box.x0 == box.x1) {
            if (!touchesFineDomain[left] && !touchesFineDomain[right]) {
                coarseInterfaces.push_back(box);
            }
            if (touchesFineDomain[left]) {
                // find box to except from vector
                for (pluint iBox = 0; iBox < BoxesToExceptCoarse.size(); ++iBox) {
                    if ((BoxesToExceptCoarse[iBox].x0 == box.x0)
                        && (BoxesToExceptCoarse[iBox].x1 == box.x0)) {
                        boxForExcept = BoxesToExceptCoarse[iBox];
                        boxFound = true;
                    }
                }
                if (boxFound) {
                    except(box, boxForExcept, newInterfaces);
                    for (pluint i = 0; i < newInterfaces.size(); ++i) {
                        coarseInterfaces.push_back(newInterfaces[i]);
                    }
                    newInterfaces.clear();
                    boxFound = false;
                } else {
                    coarseInterfaces.push_back(box);
                }
            }
            if (touchesFineDomain[right]) {
                // find box to except from vector
                for (pluint iBox = 0; iBox < BoxesToExceptCoarse.size(); ++iBox) {
                    if ((BoxesToExceptCoarse[iBox].x0 == box.x0)
                        && (BoxesToExceptCoarse[iBox].x1 == box.x0)) {
                        boxForExcept = BoxesToExceptCoarse[iBox];
                        boxFound = true;
                    }
                }
                if (boxFound) {
                    except(box, boxForExcept, newInterfaces);
                    for (pluint i = 0; i < newInterfaces.size(); ++i) {
                        coarseInterfaces.push_back(newInterfaces[i]);
                    }
                    newInterfaces.clear();
                    boxFound = false;
                } else {
                    coarseInterfaces.push_back(box);
                }
            }
        } else if (box.y0 == box.y1) {  // interface is top or bottom
            if (!touchesFineDomain[top] && !touchesFineDomain[bottom]) {
                coarseInterfaces.push_back(box);
            }
            if (touchesFineDomain[top]) {
                // find box to except from vector
                for (pluint iBox = 0; iBox < BoxesToExceptCoarse.size(); ++iBox) {
                    if ((BoxesToExceptCoarse[iBox].y0 == box.y0)
                        && (BoxesToExceptCoarse[iBox].y1 == box.y0)) {
                        boxForExcept = BoxesToExceptCoarse[iBox];
                        boxFound = true;
                    }
                }
                if (boxFound) {
                    except(box, boxForExcept, newInterfaces);
                    for (pluint i = 0; i < newInterfaces.size(); ++i) {
                        coarseInterfaces.push_back(newInterfaces[i]);
                    }
                    newInterfaces.clear();
                    boxFound = false;
                } else {
                    coarseInterfaces.push_back(box);
                }
            }
            if (touchesFineDomain[bottom]) {
                // find box to except from vector
                for (pluint iBox = 0; iBox < BoxesToExceptCoarse.size(); ++iBox) {
                    if ((BoxesToExceptCoarse[iBox].y0 == box.y0)
                        && (BoxesToExceptCoarse[iBox].y1 == box.y0)) {
                        boxForExcept = BoxesToExceptCoarse[iBox];
                        boxFound = true;
                    }
                }
                if (boxFound) {
                    except(box, boxForExcept, newInterfaces);
                    for (pluint i = 0; i < newInterfaces.size(); ++i) {
                        coarseInterfaces.push_back(newInterfaces[i]);
                    }
                    newInterfaces.clear();
                    boxFound = false;
                } else {
                    coarseInterfaces.push_back(box);
                }
            }
        } else {  // interface is front or back (box.z0 == box.z1)
            if (!touchesFineDomain[front] && !touchesFineDomain[back]) {
                coarseInterfaces.push_back(box);
            }
            if (touchesFineDomain[front]) {
                // find box to except from vector
                for (pluint iBox = 0; iBox < BoxesToExceptCoarse.size(); ++iBox) {
                    if ((BoxesToExceptCoarse[iBox].z0 == box.z0)
                        && (BoxesToExceptCoarse[iBox].z1 == box.z0)) {
                        boxForExcept = BoxesToExceptCoarse[iBox];
                        boxFound = true;
                    }
                }
                if (boxFound) {
                    except(box, boxForExcept, newInterfaces);
                    for (pluint i = 0; i < newInterfaces.size(); ++i) {
                        coarseInterfaces.push_back(newInterfaces[i]);
                    }
                    newInterfaces.clear();
                    boxFound = false;
                } else {
                    coarseInterfaces.push_back(box);
                }
            }
            if (touchesFineDomain[back]) {
                // find box to except from vector
                for (pluint iBox = 0; iBox < BoxesToExceptCoarse.size(); ++iBox) {
                    if ((BoxesToExceptCoarse[iBox].z0 == box.z0)
                        && (BoxesToExceptCoarse[iBox].z1 == box.z0)) {
                        boxForExcept = BoxesToExceptCoarse[iBox];
                        boxFound = true;
                    }
                }
                if (boxFound) {
                    except(box, boxForExcept, newInterfaces);
                    for (pluint i = 0; i < newInterfaces.size(); ++i) {
                        coarseInterfaces.push_back(newInterfaces[i]);
                    }
                    newInterfaces.clear();
                    boxFound = false;
                } else {
                    coarseInterfaces.push_back(box);
                }
            }
        }
    }
    coarseGridInterfaces[coarseLevel] = coarseInterfaces;

    /*pcout << "new fine grid : " << std::endl;
    for (pluint iComp=0; iComp<bulks[fineLevel].size(); ++iComp){
        Box3D box(bulks[fineLevel][iComp]);
        pcout << " >>>>>> " << box.x0 << " " << box.x1 << " " << box.y0 << " "
              << box.y1 << " " << box.z0 << " " << box.z1 << std::endl;
    }

    pcout << "new coarse grid : " << std::endl;
    for (pluint iComp=0; iComp<bulks[coarseLevel].size(); ++iComp){
        Box3D box(bulks[coarseLevel][iComp]);
        pcout << " >>>>>> " << box.x0 << " " << box.x1 << " " << box.y0 << " "
              << box.y1 << " " << box.z0 << " " << box.z1 << std::endl;
    }

    pcout << "Fine Interfaces : " << std::endl;
    for (pluint iComp=0; iComp < fineGridInterfaces[fineLevel].size(); ++iComp){
        Box3D box(fineGridInterfaces[fineLevel][iComp]);
        pcout << " >>> " << box.x0 << " " << box.x1 << " " << box.y0 << " "
              << box.y1 << " " << box.z0 << " " << box.z1 << std::endl;
    }
    pcout << "....\n";
    pcout << "Coarse Interfaces : " << std::endl;
    for (pluint iComp=0; iComp < coarseGridInterfaces[coarseLevel].size(); ++iComp){
        Box3D box(coarseGridInterfaces[coarseLevel][iComp]);
        pcout << " >>> " << box.x0 << " " << box.x1 << " " << box.y0 << " "
              << box.y1 << " " << box.z0 << " " << box.z1 << std::endl;
    }*/
}

void MultiGridManagement3D::coarsen(plint fineLevel, Box3D coarseDomain)
{
    PLB_PRECONDITION(fineLevel >= 1 && fineLevel < (plint)bulks.size());
    plint coarseLevel = fineLevel - 1;

    // First, trim the domain coarseDomain in case it exceeds the extent of the
    //   multi-block, and determine whether coarseDomain touches one of the boundaries
    //   of the multi-block. This information is needed, because the coarse domain
    //   fully replaces the fine domain on boundaries of the multi-block, and there
    //   is therefore no need to create a coarse-fine coupling.
    bool touches[6] = {false, false, false, false, false, false};
    trimDomain(coarseLevel, coarseDomain, touches);

    // Convert the coarse domain to fine units.
    Box3D fineDomain(coarseDomain.multiply(2));

    /************ INCLUDED BY HELEN ******************************************/
    // Test whether the reduced fine domain touches the coarse domain to avoid
    // creating wrong interfaces
    bool touchesCoarseDomain[6] = {false, false, false, false, false, false};

    /// LEFT
    if (!touches[left]) {
        for (pluint iBlock = 0; iBlock < bulks[fineLevel].size(); ++iBlock) {
            // find bounding box of current fineLevel blocks
            if (bulks[fineLevel][iBlock].x0 + 2 == fineDomain.x0) {
                touchesCoarseDomain[left] = true;
            }
        }
    }
    /// BOTTOM
    if (!touches[bottom]) {
        for (pluint iBlock = 0; iBlock < bulks[fineLevel].size(); ++iBlock) {
            // find bounding box of current fineLevel blocks
            if (bulks[fineLevel][iBlock].y0 + 2 == fineDomain.y0) {
                touchesCoarseDomain[bottom] = true;
            }
        }
    }
    /// FRONT
    if (!touches[front]) {
        for (pluint iBlock = 0; iBlock < bulks[fineLevel].size(); ++iBlock) {
            // find bounding box of current fineLevel blocks
            if (bulks[fineLevel][iBlock].z0 + 2 == fineDomain.z0) {
                touchesCoarseDomain[front] = true;
            }
        }
    }
    /// RIGHT
    if (!touches[right]) {
        for (pluint iBlock = 0; iBlock < bulks[fineLevel].size(); ++iBlock) {
            // find bounding box of current fineLevel blocks
            if (bulks[fineLevel][iBlock].x1 - 2 == fineDomain.x1) {
                touchesCoarseDomain[right] = true;
            }
        }
    }
    /// TOP
    if (!touches[top]) {
        for (pluint iBlock = 0; iBlock < bulks[fineLevel].size(); ++iBlock) {
            // find bounding box of current fineLevel blocks
            if (bulks[fineLevel][iBlock].y1 - 2 == fineDomain.y1) {
                touchesCoarseDomain[top] = true;
            }
        }
    }
    /// BACK
    if (!touches[back]) {
        for (pluint iBlock = 0; iBlock < bulks[fineLevel].size(); ++iBlock) {
            // find bounding box of current fineLevel blocks
            if (bulks[fineLevel][iBlock].z1 - 2 == fineDomain.z1) {
                touchesCoarseDomain[back] = true;
            }
        }
    }

    /**************************************************************************/

    // The reduced fine domain is the one which is going to be excluded from
    //   the original fine lattice.
    Box3D reducedFineDomain(fineDomain.enlarge(-3));
    // Box3D reducedFineDomain(fineDomain.enlarge(-1));

    // The extended coarse domain is the one which is going to be added
    //   to the original coarse lattice.
    Box3D extendedCoarseDomain(coarseDomain.enlarge(0));
    // Box3D extendedCoarseDomain(coarseDomain.enlarge(1));

    // If the domain in question touches a boundary of the multi-block,
    //   both the reduced fine domain and the extended coarse domain are
    //   identified with the boundary location.
    if (touches[left]) {  // left
        extendedCoarseDomain.x0 += 0;
        reducedFineDomain.x0 -= 3;
    }
    if (touches[right]) {  // right
        extendedCoarseDomain.x1 -= 0;
        reducedFineDomain.x1 += 3;
    }
    if (touches[bottom]) {  // bottom
        extendedCoarseDomain.y0 += 0;
        reducedFineDomain.y0 -= 3;
    }
    if (touches[top]) {  // top
        extendedCoarseDomain.y1 -= 0;
        reducedFineDomain.y1 += 3;
    }
    if (touches[front]) {  // front
        extendedCoarseDomain.z0 += 0;
        reducedFineDomain.z0 -= 3;
    }
    if (touches[back]) {  // back
        extendedCoarseDomain.z1 -= 0;
        reducedFineDomain.z1 += 3;
    }
    /*
    if (touches[left]) { // left
        extendedCoarseDomain.x0 += 1;
        reducedFineDomain.x0 -= 1;
    }
    if (touches[right]) { // right
        extendedCoarseDomain.x1 -= 1;
        reducedFineDomain.x1 += 1;
    }
    if (touches[bottom]) { // bottom
        extendedCoarseDomain.y0 += 1;
        reducedFineDomain.y0 -= 1;
    }
    if (touches[top]) { // top
        extendedCoarseDomain.y1 -= 1;
        reducedFineDomain.y1 += 1;
    }
    if (touches[front]) { // front
        extendedCoarseDomain.z0 += 1;
        reducedFineDomain.z0 -= 1;
    }
    if (touches[back]) { // back
        extendedCoarseDomain.z1 -= 1;
        reducedFineDomain.z1 += 1;
    }*/

    // If the domain in question touches the coarse domain,
    //   both the reduced fine domain and the extended coarse domain are
    //   identified with the boundary location.
    if (touchesCoarseDomain[left]) {   // left
        extendedCoarseDomain.x0 += 1;  // 1 because we need to "erase" the interface
        reducedFineDomain.x0 -= 5;     // 5 because we need to "erase" the interface
    }
    if (touchesCoarseDomain[right]) {  // right
        extendedCoarseDomain.x1 -= 1;
        reducedFineDomain.x1 += 5;
    }
    if (touchesCoarseDomain[bottom]) {  // bottom
        extendedCoarseDomain.y0 += 1;
        reducedFineDomain.y0 -= 5;
    }
    if (touchesCoarseDomain[top]) {  // top
        extendedCoarseDomain.y1 -= 1;
        reducedFineDomain.y1 += 5;
    }
    if (touchesCoarseDomain[front]) {  // front
        extendedCoarseDomain.z0 += 1;
        reducedFineDomain.z0 -= 5;
    }
    if (touchesCoarseDomain[back]) {  // back
        extendedCoarseDomain.z1 -= 1;
        reducedFineDomain.z1 += 5;
    }

    // Extract reduced fine domain from the original fine multi-block.
    std::vector<Box3D> exceptedBlocks;
    for (pluint iBlock = 0; iBlock < bulks[fineLevel].size(); ++iBlock) {
        except(bulks[fineLevel][iBlock], reducedFineDomain, exceptedBlocks);
    }
    exceptedBlocks.swap(bulks[fineLevel]);

    // Add extended coarse domain to the original coarse multi-block.
    bulks[coarseLevel].push_back(extendedCoarseDomain);

    /*pcout << "new fine grid : " << std::endl;
    for (pluint iComp=0; iComp<bulks[fineLevel].size(); ++iComp){
        Box3D box(bulks[fineLevel][iComp]);
        pcout << " >>>>>> " << box.x0 << " " << box.x1 << " " << box.y0 << " "
              << box.y1 << " " << box.z0 << " " << box.z1 << std::endl;
    }*/

    /*pcout << "Added to the coarse grid : " << std::endl;
    Box3D box(extendedCoarseDomain);
    pcout << " >>> " << box.x0 << " " << box.x1 << " " << box.y0 << " "
              << box.y1 << " " << box.z0 << " " << box.z1 << std::endl;*/
    /*pcout << "new coarse grid : " << std::endl;
    for (pluint iComp=0; iComp<bulks[coarseLevel].size(); ++iComp){
        Box3D box(bulks[coarseLevel][iComp]);
        pcout << " >>>>>> " << box.x0 << " " << box.x1 << " " << box.y0 << " "
              << box.y1 << " " << box.z0 << " " << box.z1 << std::endl;
    }*/

    /**** CHANGED BY HELEN ****/
    // needed to create interfaces properly
    Box3D interfaceDomain(coarseDomain);
    if (!touches[left] && !touchesCoarseDomain[left]) {
        interfaceDomain.x0 += 1;
    } else if (touchesCoarseDomain[left])
        interfaceDomain.x0 -= 1;
    if (!touches[bottom] && !touchesCoarseDomain[bottom]) {
        interfaceDomain.y0 += 1;
    } else if (touchesCoarseDomain[bottom])
        interfaceDomain.y0 -= 1;
    if (!touches[front] && !touchesCoarseDomain[front]) {
        interfaceDomain.z0 += 1;
    } else if (touchesCoarseDomain[front])
        interfaceDomain.z0 -= 1;
    if (!touches[right] && !touchesCoarseDomain[right]) {
        interfaceDomain.x1 -= 1;
    } else if (touchesCoarseDomain[right])
        interfaceDomain.x1 += 1;
    if (!touches[top] && !touchesCoarseDomain[top]) {
        interfaceDomain.y1 -= 1;
    } else if (touchesCoarseDomain[top])
        interfaceDomain.y1 += 1;
    if (!touches[back] && !touchesCoarseDomain[back]) {
        interfaceDomain.z1 -= 1;
    } else if (touchesCoarseDomain[back])
        interfaceDomain.z1 += 1;

    Box3D coarseInterfaceDomain(extendedCoarseDomain);
    if (touchesCoarseDomain[left])
        coarseInterfaceDomain.x0 -= 1;
    if (touchesCoarseDomain[bottom])
        coarseInterfaceDomain.y0 -= 1;
    if (touchesCoarseDomain[front])
        coarseInterfaceDomain.z0 -= 1;
    if (touchesCoarseDomain[right])
        coarseInterfaceDomain.x1 += 1;
    if (touchesCoarseDomain[top])
        coarseInterfaceDomain.y1 += 1;
    if (touchesCoarseDomain[back])
        coarseInterfaceDomain.z1 += 1;

    // Define coupling interfaces for all four sides of the refined domain, unless they
    //   touch a boundary of the multi-block or the coarse grid
    if (!touches[left] && !touchesCoarseDomain[left]) {  /// LEFT
        fineGridInterfaces[fineLevel].push_back(Box3D(
            interfaceDomain.x0, interfaceDomain.x0, interfaceDomain.y0, interfaceDomain.y1,
            interfaceDomain.z0, interfaceDomain.z1));
        coarseGridInterfaces[coarseLevel].push_back(Box3D(
            coarseInterfaceDomain.x0, coarseInterfaceDomain.x0, coarseInterfaceDomain.y0,
            coarseInterfaceDomain.y1, coarseInterfaceDomain.z0, coarseInterfaceDomain.z1));
    }
    if (!touches[right] && !touchesCoarseDomain[right]) {  /// RIGHT
        fineGridInterfaces[fineLevel].push_back(Box3D(
            interfaceDomain.x1, interfaceDomain.x1, interfaceDomain.y0, interfaceDomain.y1,
            interfaceDomain.z0, interfaceDomain.z1));
        coarseGridInterfaces[coarseLevel].push_back(Box3D(
            coarseInterfaceDomain.x1, coarseInterfaceDomain.x1, coarseInterfaceDomain.y0,
            coarseInterfaceDomain.y1, coarseInterfaceDomain.z0, coarseInterfaceDomain.z1));
    }
    if (!touches[bottom] && !touchesCoarseDomain[bottom]) {  /// BOTTOM
        fineGridInterfaces[fineLevel].push_back(Box3D(
            interfaceDomain.x0, interfaceDomain.x1, interfaceDomain.y0, interfaceDomain.y0,
            interfaceDomain.z0, interfaceDomain.z1));
        coarseGridInterfaces[coarseLevel].push_back(Box3D(
            coarseInterfaceDomain.x0, coarseInterfaceDomain.x1, coarseInterfaceDomain.y0,
            coarseInterfaceDomain.y0, coarseInterfaceDomain.z0, coarseInterfaceDomain.z1));
    }
    if (!touches[top] && !touchesCoarseDomain[top]) {  /// TOP
        fineGridInterfaces[fineLevel].push_back(Box3D(
            interfaceDomain.x0, interfaceDomain.x1, interfaceDomain.y1, interfaceDomain.y1,
            interfaceDomain.z0, interfaceDomain.z1));
        coarseGridInterfaces[coarseLevel].push_back(Box3D(
            coarseInterfaceDomain.x0, coarseInterfaceDomain.x1, coarseInterfaceDomain.y1,
            coarseInterfaceDomain.y1, coarseInterfaceDomain.z0, coarseInterfaceDomain.z1));
    }
    if (!touches[front] && !touchesCoarseDomain[front]) {  /// FRONT
        fineGridInterfaces[fineLevel].push_back(Box3D(
            interfaceDomain.x0, interfaceDomain.x1, interfaceDomain.y0, interfaceDomain.y1,
            interfaceDomain.z0, interfaceDomain.z0));
        coarseGridInterfaces[coarseLevel].push_back(Box3D(
            coarseInterfaceDomain.x0, coarseInterfaceDomain.x1, coarseInterfaceDomain.y0,
            coarseInterfaceDomain.y1, coarseInterfaceDomain.z0, coarseInterfaceDomain.z0));
    }
    if (!touches[back] && !touchesCoarseDomain[back]) {  /// BACK
        fineGridInterfaces[fineLevel].push_back(Box3D(
            interfaceDomain.x0, interfaceDomain.x1, interfaceDomain.y0, interfaceDomain.y1,
            interfaceDomain.z1, interfaceDomain.z1));
        coarseGridInterfaces[coarseLevel].push_back(Box3D(
            coarseInterfaceDomain.x0, coarseInterfaceDomain.x1, coarseInterfaceDomain.y0,
            coarseInterfaceDomain.y1, coarseInterfaceDomain.z1, coarseInterfaceDomain.z1));
    }
    /**************************************************************************/
    /*// Define coupling interfaces for all four sides of the refined domain, unless they
    //   touch a boundary of the multi-block.
    if (!touches[left]) { ///LEFT
        fineGridInterfaces[fineLevel].push_back( Box3D( coarseDomain.x0, coarseDomain.x0,
                                                        coarseDomain.y0, coarseDomain.y1,
                                                        coarseDomain.z0, coarseDomain.z1) );
        coarseGridInterfaces[coarseLevel].push_back( Box3D( extendedCoarseDomain.x0,
    extendedCoarseDomain.x0, extendedCoarseDomain.y0, extendedCoarseDomain.y1,
                                                          extendedCoarseDomain.z0,
    extendedCoarseDomain.z1) );
    }
    if (!touches[right]) { ///RIGHT
        fineGridInterfaces[fineLevel].push_back( Box3D( coarseDomain.x1, coarseDomain.x1,
                                                          coarseDomain.y0, coarseDomain.y1,
                                                          coarseDomain.z0, coarseDomain.z1) );
        coarseGridInterfaces[coarseLevel].push_back( Box3D( extendedCoarseDomain.x1,
    extendedCoarseDomain.x1, extendedCoarseDomain.y0, extendedCoarseDomain.y1,
                                                        extendedCoarseDomain.z0,
    extendedCoarseDomain.z1) );
    }
    if (!touches[bottom]) { ///BOTTOM
        fineGridInterfaces[fineLevel].push_back( Box3D( coarseDomain.x0, coarseDomain.x1,
                                                          coarseDomain.y0, coarseDomain.y0,
                                                          coarseDomain.z0, coarseDomain.z1) );
        coarseGridInterfaces[coarseLevel].push_back( Box3D( extendedCoarseDomain.x0,
    extendedCoarseDomain.x1, extendedCoarseDomain.y0, extendedCoarseDomain.y0,
                                                        extendedCoarseDomain.z0,
    extendedCoarseDomain.z1) );
    }
    if (!touches[top]) { ///TOP
        fineGridInterfaces[fineLevel].push_back( Box3D( coarseDomain.x0, coarseDomain.x1,
                                                          coarseDomain.y1, coarseDomain.y1,
                                                          coarseDomain.z0, coarseDomain.z1) );
        coarseGridInterfaces[coarseLevel].push_back( Box3D( extendedCoarseDomain.x0,
    extendedCoarseDomain.x1, extendedCoarseDomain.y1, extendedCoarseDomain.y1,
                                                          extendedCoarseDomain.z0,
    extendedCoarseDomain.z1) );
    }
    if (!touches[front]) { ///FRONT
        fineGridInterfaces[fineLevel].push_back( Box3D( coarseDomain.x0, coarseDomain.x1,
                                                          coarseDomain.y0, coarseDomain.y1,
                                                          coarseDomain.z0, coarseDomain.z0) );
        coarseGridInterfaces[coarseLevel].push_back( Box3D( extendedCoarseDomain.x0,
    extendedCoarseDomain.x1, extendedCoarseDomain.y0, extendedCoarseDomain.y1,
                                                          extendedCoarseDomain.z0,
    extendedCoarseDomain.z0) );
    }
    if (!touches[back]) { ///BACK
        fineGridInterfaces[fineLevel].push_back( Box3D( coarseDomain.x0, coarseDomain.x1,
                                                          coarseDomain.y0, coarseDomain.y1,
                                                          coarseDomain.z1, coarseDomain.z1) );
        coarseGridInterfaces[coarseLevel].push_back( Box3D( extendedCoarseDomain.x0,
    extendedCoarseDomain.x1, extendedCoarseDomain.y0, extendedCoarseDomain.y1,
                                                            extendedCoarseDomain.z1,
    extendedCoarseDomain.z1) );
    }*/
    /********** ADDED BY HELEN *************************************************/
    // go through the interfaces to make sure that there are none
    // that don't need to be there!
    std::vector<Box3D> fineInterfaces;
    Box3D box;
    std::vector<Box3D> newInterfaces;
    Box3D boxForExcept;
    plint bfex0, bfex1, bfey0, bfey1, bfez0, bfez1;  // bfe=boxForExcept
    for (pluint iComp = 0; iComp < fineGridInterfaces[fineLevel].size(); ++iComp) {
        box = fineGridInterfaces[fineLevel][iComp];
        // interface is left or right
        if (box.x0 == box.x1) {
            if (!touchesCoarseDomain[left] && !touchesCoarseDomain[right]) {
                fineInterfaces.push_back(box);
            } else if (touchesCoarseDomain[left] && (box.x0 == extendedCoarseDomain.x0 - 2)) {
                bfex0 = extendedCoarseDomain.x0 - 2;

                if (touches[right] || touchesCoarseDomain[right]) {
                    bfex1 = extendedCoarseDomain.x1 + 2;
                } else
                    bfex1 = extendedCoarseDomain.x1 - 2;
                if (touches[bottom] || touchesCoarseDomain[bottom]) {
                    bfey0 = extendedCoarseDomain.y0 - 2;
                } else
                    bfey0 = extendedCoarseDomain.y0 + 1;  // used to be +2
                if (touches[top] || touchesCoarseDomain[top]) {
                    bfey1 = extendedCoarseDomain.y1 + 2;
                } else
                    bfey1 = extendedCoarseDomain.y1 - 2;
                if (touches[front] || touchesCoarseDomain[front]) {
                    bfez0 = extendedCoarseDomain.z0 - 2;
                } else
                    bfez0 = extendedCoarseDomain.z0 + 2;
                if (touches[back] || touchesCoarseDomain[back]) {
                    bfez1 = extendedCoarseDomain.z1 + 2;
                } else
                    bfez1 = extendedCoarseDomain.z1 - 2;
                boxForExcept = Box3D(bfex0, bfex1, bfey0, bfey1, bfez0, bfez1);
                except(box, boxForExcept, newInterfaces);
                for (pluint i = 0; i < newInterfaces.size(); ++i) {
                    fineInterfaces.push_back(newInterfaces[i]);
                }
                newInterfaces.clear();
            } else if (touchesCoarseDomain[right] && (box.x1 == extendedCoarseDomain.x1 + 2)) {
                if (touches[left] || touchesCoarseDomain[left]) {
                    bfex0 = extendedCoarseDomain.x0 - 2;
                } else
                    bfex0 = extendedCoarseDomain.x0 + 2;

                bfex1 = extendedCoarseDomain.x1 + 2;

                if (touches[bottom] || touchesCoarseDomain[bottom]) {
                    bfey0 = extendedCoarseDomain.y0 - 2;
                } else
                    bfey0 = extendedCoarseDomain.y0 + 1;
                if (touches[top] || touchesCoarseDomain[top]) {
                    bfey1 = extendedCoarseDomain.y1 + 2;
                } else
                    bfey1 = extendedCoarseDomain.y1 - 2;
                if (touches[front] || touchesCoarseDomain[front]) {
                    bfez0 = extendedCoarseDomain.z0 - 2;
                } else
                    bfez0 = extendedCoarseDomain.z0 + 2;
                if (touches[back] || touchesCoarseDomain[back]) {
                    bfez1 = extendedCoarseDomain.z1 + 2;
                } else
                    bfez1 = extendedCoarseDomain.z1 - 2;
                boxForExcept = Box3D(bfex0, bfex1, bfey0, bfey1, bfez0, bfez1);
                except(box, boxForExcept, newInterfaces);
                for (pluint i = 0; i < newInterfaces.size(); ++i) {
                    fineInterfaces.push_back(newInterfaces[i]);
                }
                newInterfaces.clear();
            } else
                fineInterfaces.push_back(box);
        } else if (box.y0 == box.y1) {  // interface is top or bottom
            if (!touchesCoarseDomain[top] && !touchesCoarseDomain[bottom]) {
                fineInterfaces.push_back(box);
            } else if (touchesCoarseDomain[bottom] && (box.y0 == extendedCoarseDomain.y0 - 2)) {
                if (touches[left] || touchesCoarseDomain[left]) {
                    bfex0 = extendedCoarseDomain.x0 - 2;
                } else
                    bfex0 = extendedCoarseDomain.x0 + 2;
                if (touches[right] || touchesCoarseDomain[right]) {
                    bfex1 = extendedCoarseDomain.x1 + 2;
                } else
                    bfex1 = extendedCoarseDomain.x1 - 2;

                bfey0 = extendedCoarseDomain.y0 - 2;

                if (touches[top] || touchesCoarseDomain[top]) {
                    bfey1 = extendedCoarseDomain.y1 + 2;
                } else
                    bfey1 = extendedCoarseDomain.y1 - 2;
                if (touches[front] || touchesCoarseDomain[front]) {
                    bfez0 = extendedCoarseDomain.z0 - 2;
                } else
                    bfez0 = extendedCoarseDomain.z0 + 2;
                if (touches[back] || touchesCoarseDomain[back]) {
                    bfez1 = extendedCoarseDomain.z1 + 2;
                } else
                    bfez1 = extendedCoarseDomain.z1 - 2;
                boxForExcept = Box3D(bfex0, bfex1, bfey0, bfey1, bfez0, bfez1);
                except(box, boxForExcept, newInterfaces);
                for (pluint i = 0; i < newInterfaces.size(); ++i) {
                    fineInterfaces.push_back(newInterfaces[i]);
                }
                newInterfaces.clear();
            } else if (touchesCoarseDomain[top] && (box.y1 == extendedCoarseDomain.y1 + 2)) {
                if (touches[left] || touchesCoarseDomain[left]) {
                    bfex0 = extendedCoarseDomain.x0 - 2;
                } else
                    bfex0 = extendedCoarseDomain.x0 + 2;
                if (touches[right] || touchesCoarseDomain[right]) {
                    bfex1 = extendedCoarseDomain.x1 + 2;
                } else
                    bfex1 = extendedCoarseDomain.x1 - 2;
                if (touches[bottom] || touchesCoarseDomain[bottom]) {
                    bfey0 = extendedCoarseDomain.y0 - 2;
                } else
                    bfey0 = extendedCoarseDomain.y0 + 2;

                bfey1 = extendedCoarseDomain.y1 + 2;

                if (touches[front] || touchesCoarseDomain[front]) {
                    bfez0 = extendedCoarseDomain.z0 - 2;
                } else
                    bfez0 = extendedCoarseDomain.z0 + 2;
                if (touches[back] || touchesCoarseDomain[back]) {
                    bfez1 = extendedCoarseDomain.z1 + 2;
                } else
                    bfez1 = extendedCoarseDomain.z1 - 2;
                boxForExcept = Box3D(bfex0, bfex1, bfey0, bfey1, bfez0, bfez1);
                except(box, boxForExcept, newInterfaces);
                for (pluint i = 0; i < newInterfaces.size(); ++i) {
                    fineInterfaces.push_back(newInterfaces[i]);
                }
                newInterfaces.clear();
            } else
                fineInterfaces.push_back(box);
        } else {  // interface is front or back (box.z0 == box.z1)
            if (!touchesCoarseDomain[front] && !touchesCoarseDomain[back]) {
                fineInterfaces.push_back(box);
            } else if (touchesCoarseDomain[front] && (box.z0 == extendedCoarseDomain.z0 - 2)) {
                if (touches[left] || touchesCoarseDomain[left]) {
                    bfex0 = extendedCoarseDomain.x0 - 2;
                } else
                    bfex0 = extendedCoarseDomain.x0 + 2;
                if (touches[right] || touchesCoarseDomain[right]) {
                    bfex1 = extendedCoarseDomain.x1 + 2;
                } else
                    bfex1 = extendedCoarseDomain.x1 - 2;
                if (touches[bottom] || touchesCoarseDomain[bottom]) {
                    bfey0 = extendedCoarseDomain.y0 - 2;
                } else
                    bfey0 = extendedCoarseDomain.y0;
                if (touches[top] || touchesCoarseDomain[top]) {
                    bfey1 = extendedCoarseDomain.y1 + 2;
                } else
                    bfey1 = extendedCoarseDomain.y1 - 2;

                bfez0 = extendedCoarseDomain.z0 - 2;

                if (touches[back] || touchesCoarseDomain[back]) {
                    bfez1 = extendedCoarseDomain.z1 + 2;
                } else
                    bfez1 = extendedCoarseDomain.z1 - 2;
                boxForExcept = Box3D(bfex0, bfex1, bfey0, bfey1, bfez0, bfez1);
                except(box, boxForExcept, newInterfaces);
                for (pluint i = 0; i < newInterfaces.size(); ++i) {
                    fineInterfaces.push_back(newInterfaces[i]);
                }
                newInterfaces.clear();
            } else if (touchesCoarseDomain[back] && (box.z1 == extendedCoarseDomain.z1 + 2)) {
                if (touches[left] || touchesCoarseDomain[left]) {
                    bfex0 = extendedCoarseDomain.x0 - 2;
                } else
                    bfex0 = extendedCoarseDomain.x0 + 2;
                if (touches[right] || touchesCoarseDomain[right]) {
                    bfex1 = extendedCoarseDomain.x1 + 2;
                } else
                    bfex1 = extendedCoarseDomain.x1 - 2;
                if (touches[bottom] || touchesCoarseDomain[bottom]) {
                    bfey0 = extendedCoarseDomain.y0 - 2;
                } else
                    bfey0 = extendedCoarseDomain.y0;
                if (touches[top] || touchesCoarseDomain[top]) {
                    bfey1 = extendedCoarseDomain.y1 + 2;
                } else
                    bfey1 = extendedCoarseDomain.y1 - 2;
                if (touches[front] || touchesCoarseDomain[front]) {
                    bfez0 = extendedCoarseDomain.z0 - 2;
                } else
                    bfez0 = extendedCoarseDomain.z0 + 2;

                bfez1 = extendedCoarseDomain.z1 + 2;

                boxForExcept = Box3D(bfex0, bfex1, bfey0, bfey1, bfez0, bfez1);
                except(box, boxForExcept, newInterfaces);
                for (pluint i = 0; i < newInterfaces.size(); ++i) {
                    fineInterfaces.push_back(newInterfaces[i]);
                }
                newInterfaces.clear();
            } else
                fineInterfaces.push_back(box);
        }
    }
    fineGridInterfaces[fineLevel] = fineInterfaces;

    std::vector<Box3D> coarseInterfaces;
    newInterfaces.clear();
    for (pluint iComp = 0; iComp < coarseGridInterfaces[coarseLevel].size(); ++iComp) {
        box = coarseGridInterfaces[coarseLevel][iComp];
        // interface is left or right
        if (box.x0 == box.x1) {
            if (!touchesCoarseDomain[left] && !touchesCoarseDomain[right]) {
                coarseInterfaces.push_back(box);
            } else if (touchesCoarseDomain[left] && (box.x0 == extendedCoarseDomain.x0 - 1)) {
                bfex0 = extendedCoarseDomain.x0 - 1;

                if (touches[right] || touchesCoarseDomain[right]) {
                    bfex1 = extendedCoarseDomain.x1 + 1;
                } else
                    bfex1 = extendedCoarseDomain.x1 - 1;
                if (touches[bottom] || touchesCoarseDomain[bottom]) {
                    bfey0 = extendedCoarseDomain.y0 - 1;
                } else
                    bfey0 = extendedCoarseDomain.y0;  // used to be +1
                if (touches[top] || touchesCoarseDomain[top]) {
                    bfey1 = extendedCoarseDomain.y1 + 1;
                } else
                    bfey1 = extendedCoarseDomain.y1 - 1;
                if (touches[front] || touchesCoarseDomain[front]) {
                    bfez0 = extendedCoarseDomain.z0 - 1;
                } else
                    bfez0 = extendedCoarseDomain.z0 + 1;
                if (touches[back] || touchesCoarseDomain[back]) {
                    bfez1 = extendedCoarseDomain.z1 + 1;
                } else
                    bfez1 = extendedCoarseDomain.z1 - 1;
                boxForExcept = Box3D(bfex0, bfex1, bfey0, bfey1, bfez0, bfez1);
                except(box, boxForExcept, newInterfaces);
                for (pluint i = 0; i < newInterfaces.size(); ++i) {
                    coarseInterfaces.push_back(newInterfaces[i]);
                }
                newInterfaces.clear();
            } else if (touchesCoarseDomain[right] && (box.x1 == extendedCoarseDomain.x1 + 1)) {
                if (touches[left] || touchesCoarseDomain[left]) {
                    bfex0 = extendedCoarseDomain.x0 - 1;
                } else
                    bfex0 = extendedCoarseDomain.x0 + 1;

                bfex1 = extendedCoarseDomain.x1 + 1;

                if (touches[bottom] || touchesCoarseDomain[bottom]) {
                    bfey0 = extendedCoarseDomain.y0 - 1;
                } else
                    bfey0 = extendedCoarseDomain.y0;
                if (touches[top] || touchesCoarseDomain[top]) {
                    bfey1 = extendedCoarseDomain.y1 + 1;
                } else
                    bfey1 = extendedCoarseDomain.y1 - 1;
                if (touches[front] || touchesCoarseDomain[front]) {
                    bfez0 = extendedCoarseDomain.z0 - 1;
                } else
                    bfez0 = extendedCoarseDomain.z0 + 1;
                if (touches[back] || touchesCoarseDomain[back]) {
                    bfez1 = extendedCoarseDomain.z1 + 1;
                } else
                    bfez1 = extendedCoarseDomain.z1 - 1;
                boxForExcept = Box3D(bfex0, bfex1, bfey0, bfey1, bfez0, bfez1);
                except(box, boxForExcept, newInterfaces);
                for (pluint i = 0; i < newInterfaces.size(); ++i) {
                    coarseInterfaces.push_back(newInterfaces[i]);
                }
                newInterfaces.clear();
            } else
                coarseInterfaces.push_back(box);
        } else if (box.y0 == box.y1) {  // interface is top or bottom
            if (!touchesCoarseDomain[top] && !touchesCoarseDomain[bottom]) {
                coarseInterfaces.push_back(box);
            } else if (touchesCoarseDomain[bottom] && (box.y0 == extendedCoarseDomain.y0 - 1)) {
                if (touches[left] || touchesCoarseDomain[left]) {
                    bfex0 = extendedCoarseDomain.x0 - 1;
                } else
                    bfex0 = extendedCoarseDomain.x0 + 1;
                if (touches[right] || touchesCoarseDomain[right]) {
                    bfex1 = extendedCoarseDomain.x1 + 1;
                } else
                    bfex1 = extendedCoarseDomain.x1 - 1;

                bfey0 = extendedCoarseDomain.y0 - 1;

                if (touches[top] || touchesCoarseDomain[top]) {
                    bfey1 = extendedCoarseDomain.y1 + 1;
                } else
                    bfey1 = extendedCoarseDomain.y1 - 1;
                if (touches[front] || touchesCoarseDomain[front]) {
                    bfez0 = extendedCoarseDomain.z0 - 1;
                } else
                    bfez0 = extendedCoarseDomain.z0 + 1;
                if (touches[back] || touchesCoarseDomain[back]) {
                    bfez1 = extendedCoarseDomain.z1 + 1;
                } else
                    bfez1 = extendedCoarseDomain.z1 - 1;
                boxForExcept = Box3D(bfex0, bfex1, bfey0, bfey1, bfez0, bfez1);
                except(box, boxForExcept, newInterfaces);
                for (pluint i = 0; i < newInterfaces.size(); ++i) {
                    coarseInterfaces.push_back(newInterfaces[i]);
                }
                newInterfaces.clear();
            } else if (touchesCoarseDomain[top] && (box.y1 == extendedCoarseDomain.y1 + 1)) {
                if (touches[left] || touchesCoarseDomain[left]) {
                    bfex0 = extendedCoarseDomain.x0 - 1;
                } else
                    bfex0 = extendedCoarseDomain.x0 + 1;
                if (touches[right] || touchesCoarseDomain[right]) {
                    bfex1 = extendedCoarseDomain.x1 + 1;
                } else
                    bfex1 = extendedCoarseDomain.x1 - 1;
                if (touches[bottom] || touchesCoarseDomain[bottom]) {
                    bfey0 = extendedCoarseDomain.y0 - 1;
                } else
                    bfey0 = extendedCoarseDomain.y0 + 1;

                bfey1 = extendedCoarseDomain.y1 + 1;

                if (touches[front] || touchesCoarseDomain[front]) {
                    bfez0 = extendedCoarseDomain.z0 - 1;
                } else
                    bfez0 = extendedCoarseDomain.z0 + 1;
                if (touches[back] || touchesCoarseDomain[back]) {
                    bfez1 = extendedCoarseDomain.z1 + 1;
                } else
                    bfez1 = extendedCoarseDomain.z1 - 1;
                boxForExcept = Box3D(bfex0, bfex1, bfey0, bfey1, bfez0, bfez1);
                except(box, boxForExcept, newInterfaces);
                for (pluint i = 0; i < newInterfaces.size(); ++i) {
                    coarseInterfaces.push_back(newInterfaces[i]);
                }
                newInterfaces.clear();
            }
        } else {  // interface is front or back (box.z0 == box.z1)
            if (!touchesCoarseDomain[front] && !touchesCoarseDomain[back]) {
                coarseInterfaces.push_back(box);
            } else if (touchesCoarseDomain[front] && (box.z0 == extendedCoarseDomain.z0 - 1)) {
                if (touches[left] || touchesCoarseDomain[left]) {
                    bfex0 = extendedCoarseDomain.x0 - 1;
                } else
                    bfex0 = extendedCoarseDomain.x0 + 1;
                if (touches[right] || touchesCoarseDomain[right]) {
                    bfex1 = extendedCoarseDomain.x1 + 1;
                } else
                    bfex1 = extendedCoarseDomain.x1 - 1;
                if (touches[bottom] || touchesCoarseDomain[bottom]) {
                    bfey0 = extendedCoarseDomain.y0 - 1;
                } else
                    bfey0 = extendedCoarseDomain.y0;
                if (touches[top] || touchesCoarseDomain[top]) {
                    bfey1 = extendedCoarseDomain.y1 + 1;
                } else
                    bfey1 = extendedCoarseDomain.y1 - 1;

                bfez0 = extendedCoarseDomain.z0 - 1;

                if (touches[back] || touchesCoarseDomain[back]) {
                    bfez1 = extendedCoarseDomain.z1 + 1;
                } else
                    bfez1 = extendedCoarseDomain.z1 - 1;
                boxForExcept = Box3D(bfex0, bfex1, bfey0, bfey1, bfez0, bfez1);
                except(box, boxForExcept, newInterfaces);
                for (pluint i = 0; i < newInterfaces.size(); ++i) {
                    coarseInterfaces.push_back(newInterfaces[i]);
                }
                newInterfaces.clear();
            } else if (touchesCoarseDomain[front] && (box.z1 == extendedCoarseDomain.z1 + 1)) {
                if (touches[left] || touchesCoarseDomain[left]) {
                    bfex0 = extendedCoarseDomain.x0 - 1;
                } else
                    bfex0 = extendedCoarseDomain.x0 + 1;
                if (touches[right] || touchesCoarseDomain[right]) {
                    bfex1 = extendedCoarseDomain.x1 + 1;
                } else
                    bfex1 = extendedCoarseDomain.x1 - 1;
                if (touches[bottom] || touchesCoarseDomain[bottom]) {
                    bfey0 = extendedCoarseDomain.y0 - 1;
                } else
                    bfey0 = extendedCoarseDomain.y0;
                if (touches[top] || touchesCoarseDomain[top]) {
                    bfey1 = extendedCoarseDomain.y1 + 1;
                } else
                    bfey1 = extendedCoarseDomain.y1 - 1;
                if (touches[front] || touchesCoarseDomain[front]) {
                    bfez0 = extendedCoarseDomain.z0 - 1;
                } else
                    bfez0 = extendedCoarseDomain.z0 + 1;

                bfez1 = extendedCoarseDomain.z1 + 1;

                boxForExcept = Box3D(bfex0, bfex1, bfey0, bfey1, bfez0, bfez1);
                except(box, boxForExcept, newInterfaces);
                for (pluint i = 0; i < newInterfaces.size(); ++i) {
                    coarseInterfaces.push_back(newInterfaces[i]);
                }
                newInterfaces.clear();
            } else
                coarseInterfaces.push_back(box);
        }
    }
    coarseGridInterfaces[coarseLevel] = coarseInterfaces;

    /**************************************************************************/

    /*pcout << "Fine Interfaces : " << std::endl;
    for (pluint iComp=0; iComp < fineGridInterfaces[fineLevel].size(); ++iComp){
        Box3D box(fineGridInterfaces[fineLevel][iComp]);
        pcout << " >>> " << box.x0 << " " << box.x1 << " " << box.y0 << " "
              << box.y1 << " " << box.z0 << " " << box.z1 << std::endl;
    }
    pcout << "....\n";
    pcout << "Coarse Interfaces : " << std::endl;
    for (pluint iComp=0; iComp < coarseGridInterfaces[coarseLevel].size(); ++iComp){
        Box3D box(coarseGridInterfaces[coarseLevel][iComp]);
        pcout << " >>> " << box.x0 << " " << box.x1 << " " << box.y0 << " "
              << box.y1 << " " << box.z0 << " " << box.z1 << std::endl;
    }*/
}

std::vector<std::vector<Box3D> > const &MultiGridManagement3D::getCoarseInterface() const
{
    return coarseGridInterfaces;
}

std::vector<std::vector<Box3D> > const &MultiGridManagement3D::getFineInterface() const
{
    return fineGridInterfaces;
}

void MultiGridManagement3D::trimDomain(plint whichLevel, Box3D &domain, bool *touches) const
{
    Box3D bbox = boundingBoxes[whichLevel];
    /// LEFT
    if (domain.x0 <= bbox.x0) {
        domain.x0 = bbox.x0;  // Trim in case multi-block extent is exceeded.
        touches[left] = true;
    }
    /// BOTTOM
    if (domain.y0 <= bbox.y0) {
        domain.y0 = bbox.y0;  // Trim in case multi-block extent is exceeded.
        touches[bottom] = true;
    }
    /// FRONT
    if (domain.z0 <= bbox.z0) {
        domain.z0 = bbox.z0;  // Trim in case multi-block extent is exceeded.
        touches[front] = true;
    }
    /// RIGHT
    if (domain.x1 >= bbox.x1) {
        domain.x1 = bbox.x1;  // Trim in case multi-block extent is exceeded.
        touches[right] = true;
    }
    /// TOP
    if (domain.y1 >= bbox.y1) {
        domain.y1 = bbox.y1;  // Trim in case multi-block extent is exceeded.
        touches[top] = true;
    }
    /// BACK
    if (domain.z1 >= bbox.z1) {
        domain.z1 = bbox.z1;  // Trim in case multi-block extent is exceeded.
        touches[back] = true;
    }
}

/**
 * Use one of the Parallelizer3D objects to recompute the parallelization of the
 *   domain.
 */
void MultiGridManagement3D::parallelize(Parallelizer3D *parallelizer)
{
    parallelizer->parallelize();
    bulks = parallelizer->getRecomputedBlocks();
    mpiProcess = parallelizer->getMpiDistribution();
    delete parallelizer;
}

plint MultiGridManagement3D::getNumLevels() const
{
    return (plint)bulks.size();
}

plint MultiGridManagement3D::getReferenceLevel() const
{
    return referenceLevel;
}

std::vector<std::vector<plint> > const &MultiGridManagement3D::getMpiProcesses() const
{
    return mpiProcess;
}

/// Extract a domain (in coarse coordinates) from a MultiGridManagement3D
/** This is achieved by taking all the geometric information of the MultiGridManagement3D
 *  and intersecting it with the so called coarsestDomain. All the intersecting parts
 *  are kept on the new MultiGridManagement3D.
 *
 *  Implemented by Helen Morrison, 2014 - based on the (corrected) 2D-version
 */
MultiGridManagement3D extractManagement(
    MultiGridManagement3D management, Box3D coarsestDomain, bool crop)
{
    // TODO: implement
    // return management;

    std::unique_ptr<MultiGridManagement3D> result;

    if (crop) {
        result = std::unique_ptr<MultiGridManagement3D>(new MultiGridManagement3D(
            coarsestDomain, management.getNumLevels(), management.getReferenceLevel()));
    } else {
        result = std::unique_ptr<MultiGridManagement3D>(new MultiGridManagement3D(
            management.getBoundingBox(0), management.getNumLevels(),
            management.getReferenceLevel()));
    }

    // erase the bulk added by default by the constructor
    result->bulks[0].erase(result->bulks[0].begin());

    MultiScaleManager *scaleManager = result->scaleManager->clone();

    // every level must contain the informations intersected with coarsestDomain
    for (plint iLevel = 0; iLevel < management.getNumLevels(); ++iLevel) {
        std::vector<plint> mpiProcessLevel;
        Box3D rescaledBox = scaleManager->scaleBox(coarsestDomain, iLevel);

        // boundingBoxes
        result->boundingBoxes[iLevel] = rescaledBox;

        if (iLevel < management.getNumLevels() - 1) {
            // coarseGridInterfaces
            for (pluint iInterface = 0; iInterface < management.coarseGridInterfaces[iLevel].size();
                 ++iInterface)
            {
                Box3D intersection;
                if (intersect(
                        rescaledBox, management.coarseGridInterfaces[iLevel][iInterface],
                        intersection)) {
                    result->coarseGridInterfaces[iLevel].push_back(intersection);
                }
            }
        }

        if (iLevel > 0) {
            // fineGridInterfaces
            for (pluint iInterface = 0; iInterface < management.fineGridInterfaces[iLevel].size();
                 ++iInterface) {
                Box3D intersection;
                if (intersect(
                        rescaledBox, management.fineGridInterfaces[iLevel][iInterface],
                        intersection)) {
                    result->fineGridInterfaces[iLevel].push_back(intersection);
                }
            }
        }

        // bulk
        for (pluint iBox = 0; iBox < management.bulks[iLevel].size(); ++iBox) {
            Box3D intersection;
            // verify if the bulk is contained
            if (intersect(rescaledBox, management.bulks[iLevel][iBox], intersection)) {
                result->bulks[iLevel].push_back(intersection);
                // copy also the mpi process associated
                mpiProcessLevel.push_back(management.mpiProcess[iLevel][iBox]);
            }
        }
        result->mpiProcess.push_back(mpiProcessLevel);
    }
    delete scaleManager;
    return *result;
}

void MultiGridManagement3D::eliminateBlocksOnFinestLevel(
    MultiScalarField3D<int> &domainMatrix, plint blockSizeX, plint blockSizeY, plint blockSizeZ)
{
    plint finestLevel = this->getNumLevels() - 1;

    MultiScalarField3D<int> *newDomainMatrix = new MultiScalarField3D<int>(
        *plb::reparallelize(domainMatrix, blockSizeX, blockSizeY, blockSizeZ));

    MultiContainerBlock3D multiFlagBlock(*newDomainMatrix);
    std::vector<MultiBlock3D *> args;
    args.push_back(newDomainMatrix);
    args.push_back(&multiFlagBlock);

    ComputeSparsityFunctional3D<int> sparsityFunctional;
    applyProcessingFunctional(sparsityFunctional, newDomainMatrix->getBoundingBox(), args);

    MultiBlockManagement3D const &management = multiFlagBlock.getMultiBlockManagement();
    ThreadAttribution const &threadAttribution = management.getThreadAttribution();
    SparseBlockStructure3D const &sparseBlock = management.getSparseBlockStructure();

    std::map<plint, Box3D> const &domains = sparseBlock.getBulks();
    std::vector<plint> domainIds(domains.size());
    std::vector<int> extractThisBlock(domains.size());

    std::map<plint, Box3D>::const_iterator it = domains.begin();
    plint pos = 0;
    for (; it != domains.end(); ++it) {
        plint id = it->first;
        domainIds[pos] = id;
        if (threadAttribution.isLocal(id)) {
            AtomicContainerBlock3D const &flagBlock = multiFlagBlock.getComponent(id);
            FlagData3D const *data = dynamic_cast<FlagData3D const *>(flagBlock.getData());
            PLB_ASSERT(data);
            extractThisBlock[pos] = data->keepThisBlock ? 0 : 1;
        } else {
            extractThisBlock[pos] = 0;
        }
        ++pos;
    }

#ifdef PLB_MPI_PARALLEL
    std::vector<int> tmp(extractThisBlock.size());
    global::mpi().reduceVect(extractThisBlock, tmp, MPI_SUM);
    global::mpi().bCast(&tmp[0], tmp.size());
    tmp.swap(extractThisBlock);
#endif

    //    pcout << "old fine grid : " << std::endl;
    //    for (pluint iComp=0; iComp<bulks[finestLevel].size(); ++iComp){
    //        Box3D box(bulks[finestLevel][iComp]);
    //        pcout << " >>>>>> " << box.x0 << " " << box.x1 << " " << box.y0 << " "
    //              << box.y1 << " " << box.z0 << " " << box.z1 << std::endl;
    //    }

    for (pluint iBlock = 0; iBlock < extractThisBlock.size(); ++iBlock) {
        if (extractThisBlock[iBlock]) {
            plint id = domainIds[iBlock];
            Box3D bulk;
            sparseBlock.getBulk(id, bulk);

            // Extract bulk from fine level bulks.
            std::vector<Box3D> exceptedBulks;
            for (pluint iBlock = 0; iBlock < bulks[finestLevel].size(); ++iBlock) {
                except(bulks[finestLevel][iBlock], bulk, exceptedBulks);
            }
            exceptedBulks.swap(bulks[finestLevel]);
        }
    }

    //    pcout << "new fine grid : " << std::endl;
    //    for (pluint iComp=0; iComp<bulks[finestLevel].size(); ++iComp){
    //        Box3D box(bulks[finestLevel][iComp]);
    //        pcout << " >>>>>> " << box.x0 << " " << box.x1 << " " << box.y0 << " "
    //              << box.y1 << " " << box.z0 << " " << box.z1 << std::endl;
    //    }
}

}  // namespace plb
