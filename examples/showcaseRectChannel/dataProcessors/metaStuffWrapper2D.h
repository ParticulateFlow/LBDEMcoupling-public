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
#ifndef META_STUFF_WRAPPER_2D_H
#define META_STUFF_WRAPPER_2D_H

#include <memory>

#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include "core/globalDefs.h"
#include "dataProcessors/metaStuffFunctional2D.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiDataField2D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
void extractTopMostDynamics(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<int> &dynamicsId, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<int> > extractTopMostDynamics(
    MultiBlockLattice2D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<int> > extractTopMostDynamics(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
void extractBottomMostDynamics(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<int> &dynamicsId, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<int> > extractBottomMostDynamics(
    MultiBlockLattice2D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<int> > extractBottomMostDynamics(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
void uniqueDynamicsChains(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain,
    std::vector<std::vector<int> > &chains, pluint &maxChainLength);

template <typename T, template <typename U> class Descriptor>
void uniqueDynamicsIds(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, std::vector<int> &ids);

template <typename T, template <typename U> class Descriptor>
void extractDynamicsChain(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<int> &dynamicsId,
    std::map<int, std::string> &nameOfDynamics, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<int> > extractDynamicsChain(
    MultiBlockLattice2D<T, Descriptor> &lattice, std::map<int, std::string> &nameOfDynamics);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<int> > extractDynamicsChain(
    MultiBlockLattice2D<T, Descriptor> &lattice, std::map<int, std::string> &nameOfDynamics,
    Box2D domain);

template <typename T, template <typename U> class Descriptor>
void copyEntireCells(
    MultiBlockLattice2D<T, Descriptor> &sourceLattice,
    MultiBlockLattice2D<T, Descriptor> &destinationLattice, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > copyEntireCells(
    MultiBlockLattice2D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > copyEntireCells(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain);

bool allFlagsTrue(MultiBlock2D *multiBlock);

}  // namespace plb

#endif  // META_STUFF_WRAPPER_2D_H
