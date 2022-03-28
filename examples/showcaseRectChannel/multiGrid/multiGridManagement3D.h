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

#ifndef MULTI_GRID_MANAGEMENT_3D_H
#define MULTI_GRID_MANAGEMENT_3D_H

#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "multiBlock/multiDataField3D.h"
#include "multiGrid/multiScale.h"
#include "multiGrid/parallelizer3D.h"

namespace plb {

class MultiGridManagement3D {
public:
    MultiGridManagement3D(
        plint coarseNx, plint coarseNy, plint coarseNz, plint numLevels, plint referenceLevel = 0);
    MultiGridManagement3D(Box3D coarseBoundingBox, plint numLevels, plint referenceLevel = 0);

    MultiGridManagement3D(MultiGridManagement3D const &rhs, Box3D coarsestDomain);
    MultiGridManagement3D(MultiGridManagement3D const &rhs);
    ~MultiGridManagement3D();

    void swap(MultiGridManagement3D &rhs);

    /// Retrieve the blocks for each level
    std::vector<Box3D> const &getBulks(plint iLevel) const;
    std::vector<std::vector<Box3D> > const &getBulks() const;

    /// Retrieve the bounding box at any given level
    Box3D getBoundingBox(plint level) const;

    /// Add a refinement in level at a given domain
    void refine(plint level, Box3D domain);
    /// Coarsen a given domain at certain level
    void coarsen(plint level, Box3D domain);

    /// Get the Ids of the parallelized regions
    std::vector<std::vector<plint> > const &getMpiProcesses() const;

    std::vector<std::vector<Box3D> > const &getCoarseInterface() const;
    std::vector<std::vector<Box3D> > const &getFineInterface() const;

    void trimDomain(plint whichLevel, Box3D &domain, bool *touches) const;

    plint getNumLevels() const;
    plint getReferenceLevel() const;

    /// Parallelize the management before the creation of the multiGrid
    void parallelize(Parallelizer3D *parallelizer);

    friend MultiGridManagement3D extractManagement(
        MultiGridManagement3D management, Box3D coarsestDomain, bool crop);

    void eliminateBlocksOnFinestLevel(
        MultiScalarField3D<int> &domainMatrix, plint blockSizeX, plint blockSizeY,
        plint blockSizeZ);

private:
    void initialize(Box3D const &level0_box);
    void parallelizeLevel(
        plint whichLevel, std::vector<Box3D> const &parallelRegions,
        std::vector<plint> const &regionIDs);

private:
    plint referenceLevel;
    std::vector<Box3D> boundingBoxes;
    /// IDs of the blocks, in case the code is parallelized.
    std::vector<std::vector<plint> > mpiProcess;
    /// Cells which form an interface with a grid of next-finer level.
    /** "Interface on which I am a coarse grid." */
    std::vector<std::vector<Box3D> > coarseGridInterfaces;
    /// Cells which form an interface with a grid of next-coarser level.
    /** "Interface on which I am a fine grid." */
    std::vector<std::vector<Box3D> > fineGridInterfaces;

    /// Geometry specification: list all atomic-blocks in each multi-block.
    std::vector<std::vector<Box3D> > bulks;

    MultiScaleManager *scaleManager;
};

MultiGridManagement3D extractManagement(
    MultiGridManagement3D management, Box3D coarsestDomain, bool crop = true);

}  // namespace plb

#endif  // MULTI_GRID_MANAGEMENT_3D_H
