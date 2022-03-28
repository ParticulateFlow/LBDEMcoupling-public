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
 * Management of multi grid domains -- Header file
 */

#ifndef MULTI_GRID_MANAGEMENT_2D_H
#define MULTI_GRID_MANAGEMENT_2D_H

#include "core/geometry2D.h"
#include "core/globalDefs.h"
#include "multiGrid/multiScale.h"
#include "multiGrid/parallelizer2D.h"

namespace plb {

class MultiGridManagement2D {
public:
    MultiGridManagement2D(
        plint coarseNx, plint coarseNy, plint numLevels, plint overlapWidth_ = 1,
        plint referenceLevel = 0);
    MultiGridManagement2D(
        Box2D coarseBoundingBox, plint numLevels, plint overlapWidth_ = 1,
        plint referenceLevel = 0);

    MultiGridManagement2D(MultiGridManagement2D const &rhs);
    ~MultiGridManagement2D();

    void swap(MultiGridManagement2D &rhs);

    /// Retrieve the blocks for each level
    std::vector<Box2D> const &getBulks(plint iLevel) const;
    std::vector<std::vector<Box2D> > const &getBulks() const;

    /// Retrieve the bounding box at any given level
    Box2D getBoundingBox(plint level) const;
    std::vector<Box2D> const &getBoundingBoxes() const;

    /// Add a refinement in level at a given domain
    void refine(plint level, Box2D domain);
    /// Coarsen a given domain in certain level
    void coarsen(plint level, Box2D coarseDomain);

    /// Add a refinement patch at a certain level. Leave level unchanged
    void refineMultiGrid(plint level, Box2D domain);

    /// Store the overlap coordinates for further usage
    /// Get the Ids of the parallelized regions
    std::vector<std::vector<plint> > const &getMpiProcesses() const;

    std::vector<std::vector<Box2D> > const &getCoarseInterface() const;
    std::vector<std::vector<Box2D> > const &getFineInterface() const;
    std::vector<std::vector<Array<plint, 2> > > const &getCoarseOrientation() const;
    std::vector<std::vector<Array<plint, 2> > > const &getFineOrientation() const;

    void trimDomain(
        plint whichLevel, Box2D &domain, bool &touchLeft, bool &touchRight, bool &touchBottom,
        bool &touchTop) const;

    plint getNumLevels() const;
    plint getReferenceLevel() const;
    plint getOverlapWidth() const;

    /// Paralllelize the domain using a certain technique
    void parallelize(Parallelizer2D *parallelizer);

    /// Eliminate the redundant fine interfaces
    void eliminateUnnecessaryFineInterfaces();

    friend MultiGridManagement2D extractManagement(
        MultiGridManagement2D management, Box2D coarsestDomain, bool crop);

private:
    void initialize(Box2D const &level0_box);
    void parallelizeLevel(
        plint whichLevel, std::vector<Box2D> const &parallelRegions,
        std::vector<plint> const &regionIDs);

private:
    plint overlapWidth;
    plint referenceLevel;

    /// All the bounding boxes
    std::vector<Box2D> boundingBoxes;

    /// Cells which form an interface with a grid of next-finer level.
    /** "Interface on which I am a coarse grid." */
    std::vector<std::vector<Box2D> > coarseGridInterfaces;
    /// Cells which form an interface with a grid of next-coarser level.
    /** "Interface on which I am a fine grid." */
    std::vector<std::vector<Box2D> > fineGridInterfaces;

    /// Geometry specification: list all atomic-blocks in each multi-block.
    std::vector<std::vector<Box2D> > bulks;

    /// directions and orientations of the regions
    std::vector<std::vector<Array<plint, 2> > > coarseInterfaceOrientations;
    std::vector<std::vector<Array<plint, 2> > > fineInterfaceOrientations;

    /// IDs of the blocks, in case the code is parallelized.
    std::vector<std::vector<plint> > mpiProcess;

    MultiScaleManager *scaleManager;
};

MultiGridManagement2D extractManagement(
    MultiGridManagement2D management, Box2D coarsestDomain, bool crop = true);

}  // namespace plb

#endif  // MULTI_GRID_MANAGEMENT_2D_H
