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
 * Helper functions for domain initialization -- header file.
 */
#ifndef META_STUFF_WRAPPER_3D_H
#define META_STUFF_WRAPPER_3D_H

#include <memory>

#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "core/globalDefs.h"
#include "dataProcessors/metaStuffFunctional3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
void extractTopMostDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<int> &dynamicsId, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<int> > extractTopMostDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<int> > extractTopMostDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
void extractBottomMostDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<int> &dynamicsId, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<int> > extractBottomMostDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<int> > extractBottomMostDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
void uniqueDynamicsChains(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain,
    std::vector<std::vector<int> > &chains, pluint &maxChainLength);

template <typename T, template <typename U> class Descriptor>
void uniqueDynamicsIds(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, std::vector<int> &ids);

template <typename T, template <typename U> class Descriptor>
void extractDynamicsChain(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<int> &dynamicsId,
    std::map<int, std::string> &nameOfDynamics, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<int> > extractDynamicsChain(
    MultiBlockLattice3D<T, Descriptor> &lattice, std::map<int, std::string> &nameOfDynamics);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<int> > extractDynamicsChain(
    MultiBlockLattice3D<T, Descriptor> &lattice, std::map<int, std::string> &nameOfDynamics,
    Box3D domain);

template <typename T, template <typename U> class Descriptor>
void copyEntireCells(
    MultiBlockLattice3D<T, Descriptor> &sourceLattice,
    MultiBlockLattice3D<T, Descriptor> &destinationLattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > copyEntireCells(
    MultiBlockLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > copyEntireCells(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain);

bool allFlagsTrue(MultiBlock3D *multiBlock);

void getThreadNum(MultiScalarField3D<int> &threadNum);

void getBlockNum(MultiScalarField3D<plint> &blockNum);

void getRandomBlockNum(MultiScalarField3D<plint> &blockNum);

}  // namespace plb

#endif  // META_STUFF_WRAPPER_3D_H
