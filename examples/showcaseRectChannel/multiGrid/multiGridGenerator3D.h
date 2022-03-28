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
 * Various factories that use a multiGridManagement3D -- Header file
 */

#ifndef MULTI_GRID_GENERATOR_3D_H
#define MULTI_GRID_GENERATOR_3D_H

#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"
#include "multiGrid/multiGridManagement3D.h"

namespace plb {

/// Compute the reduced bulk of the domain for fine to coarse copy
inline Box3D computeCopyReducedBulk(Box3D domain);
inline void computeCopyEdges(Box3D domain, std::vector<Box3D> &edges);
inline void computeCopyCorners(Box3D domain, std::vector<Box3D> &corners);

template <typename T>
void computeFilteringIndicesEdges(Box3D domain, std::vector<std::vector<plint> > &indices);

template <typename T, template <typename U> class Descriptor>
std::vector<MultiBlockLattice3D<T, Descriptor> *> generateLattices(
    MultiGridManagement3D management, std::vector<Dynamics<T, Descriptor> *> backgroundDynamics,
    std::vector<BlockCommunicator3D *> communicators,
    std::vector<CombinedStatistics *> combinedStatistics);

template <typename T, template <typename U> class Descriptor>
std::vector<MultiBlockLattice3D<T, Descriptor> *> generateLattices(
    MultiGridManagement3D management, std::vector<Dynamics<T, Descriptor> *> backgroundDynamics);

template <typename T, template <typename U> class Descriptor>
void createInterfaces(
    std::vector<MultiBlockLattice3D<T, Descriptor> *> &multiBlocks,
    MultiGridManagement3D management);

template <typename T, template <typename U> class Descriptor>
void createCoarseGridInterface(
    plint coarseLevel, Box3D coarseGridInterface,
    std::vector<MultiBlockLattice3D<T, Descriptor> *> &multiBlocks);

template <typename T, template <typename U> class Descriptor>
void createFineGridInterface(
    plint coarseLevel, Box3D fineGridInterface,
    std::vector<MultiBlockLattice3D<T, Descriptor> *> &multiBlocks);

template <typename T>
std::vector<MultiScalarField3D<T> *> generateScalarFields(
    MultiGridManagement3D const &management, std::vector<BlockCommunicator3D *> communicators,
    std::vector<CombinedStatistics *> combinedStatistics);

template <typename T, int nDim>
std::vector<MultiTensorField3D<T, nDim> *> generateTensorFields(
    MultiGridManagement3D const &management, std::vector<BlockCommunicator3D *> communicators,
    std::vector<CombinedStatistics *> combinedStatistics);

}  // namespace plb

#endif  // MULTI_GRID_GENERATOR_3D_H
