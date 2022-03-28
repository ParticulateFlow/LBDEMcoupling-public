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
#ifndef META_STUFF_WRAPPER_2D_HH
#define META_STUFF_WRAPPER_2D_HH

#include "core/dynamicsIdentifiers.h"
#include "dataProcessors/dataInitializerWrapper2D.h"
#include "dataProcessors/metaStuffFunctional2D.h"
#include "dataProcessors/metaStuffWrapper2D.h"
#include "multiBlock/multiContainerBlock2D.h"
#include "multiBlock/multiDataProcessorWrapper2D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper2D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
void extractTopMostDynamics(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<int> &dynamicsId, Box2D domain)
{
    applyProcessingFunctional(
        new ExtractTopMostDynamicsFunctional2D<T, Descriptor>(), domain, lattice, dynamicsId);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<int> > extractTopMostDynamics(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    MultiScalarField2D<int> *dynamicsId = new MultiScalarField2D<int>(lattice, domain);
    extractTopMostDynamics(lattice, *dynamicsId, domain);
    return std::unique_ptr<MultiScalarField2D<int> >(dynamicsId);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<int> > extractTopMostDynamics(
    MultiBlockLattice2D<T, Descriptor> &lattice)
{
    return extractTopMostDynamics(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void extractBottomMostDynamics(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<int> &dynamicsId, Box2D domain)
{
    applyProcessingFunctional(
        new ExtractBottomMostDynamicsFunctional2D<T, Descriptor>(), domain, lattice, dynamicsId);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<int> > extractBottomMostDynamics(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    MultiScalarField2D<int> *dynamicsId = new MultiScalarField2D<int>(lattice, domain);
    extractBottomMostDynamics(lattice, *dynamicsId, domain);
    return std::unique_ptr<MultiScalarField2D<int> >(dynamicsId);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<int> > extractBottomMostDynamics(
    MultiBlockLattice2D<T, Descriptor> &lattice)
{
    return extractBottomMostDynamics(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void uniqueDynamicsChains(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain,
    std::vector<std::vector<int> > &chains, pluint &maxChainLength)
{
    MultiContainerBlock2D container(lattice);
    std::vector<MultiBlock2D *> latticeAndContainer, containerVector;
    latticeAndContainer.push_back(&lattice);
    latticeAndContainer.push_back(&container);
    containerVector.push_back(&container);

    StoreDynamicsFunctional2D<T, Descriptor> SDFunctional;
    applyProcessingFunctional(SDFunctional, domain, latticeAndContainer);
    maxChainLength = SDFunctional.getMaxChainLength();

    std::vector<int> nextMaximum(maxChainLength);
    for (pluint i = 0; i < nextMaximum.size(); ++i) {
        nextMaximum[i] = -1;
    }
    std::vector<int> emptyVector(nextMaximum);
    std::vector<std::vector<int> > maxima;
    do {
        IterateDynamicsFunctional2D functional(nextMaximum);
        applyProcessingFunctional(functional, domain, containerVector);
        nextMaximum = functional.getNextMaximum();
        if (!vectorEquals(nextMaximum, emptyVector)) {
            maxima.push_back(nextMaximum);
        }
    } while (!vectorEquals(nextMaximum, emptyVector));
    chains.swap(maxima);
}

template <typename T, template <typename U> class Descriptor>
void uniqueDynamicsIds(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, std::vector<int> &ids)
{
    std::vector<std::vector<int> > chains;
    pluint maxChainLength;
    uniqueDynamicsChains(lattice, domain, chains, maxChainLength);
    std::set<int> idSet;
    for (pluint iChain = 0; iChain < chains.size(); ++iChain) {
        idSet.insert(chains[iChain].begin(), chains[iChain].end());
    }
    idSet.erase(-1);
    std::vector<int>(idSet.begin(), idSet.end()).swap(ids);
}

template <typename T, template <typename U> class Descriptor>
void extractDynamicsChain(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<int> &dynamicsId,
    std::map<int, std::string> &nameOfDynamics, Box2D domain)
{
    std::vector<std::vector<int> > chains;
    pluint maxChainLength;
    uniqueDynamicsChains(lattice, domain, chains, maxChainLength);
    nameOfDynamics.clear();
    typename ExtractDynamicsChainFunctional2D<T, Descriptor>::DMap dynamicsMap;
    for (pluint iChain = 0; iChain < chains.size(); ++iChain) {
        dynamicsMap[chains[iChain]] = iChain;
        nameOfDynamics[iChain] = meta::constructIdNameChain<T, Descriptor>(chains[iChain], " >> ");
    }
    setToConstant(dynamicsId, dynamicsId.getBoundingBox(), -1);
    applyProcessingFunctional(
        new ExtractDynamicsChainFunctional2D<T, Descriptor>(dynamicsMap, maxChainLength), domain,
        lattice, dynamicsId);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<int> > extractDynamicsChain(
    MultiBlockLattice2D<T, Descriptor> &lattice, std::map<int, std::string> &nameOfDynamics,
    Box2D domain)
{
    MultiScalarField2D<int> *dynamicsId = new MultiScalarField2D<int>(lattice, domain);
    extractDynamicsChain(lattice, *dynamicsId, nameOfDynamics, domain);
    return std::unique_ptr<MultiScalarField2D<int> >(dynamicsId);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<int> > extractDynamicsChain(
    MultiBlockLattice2D<T, Descriptor> &lattice, std::map<int, std::string> &nameOfDynamics)
{
    return extractDynamicsChain(lattice, nameOfDynamics, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void copyEntireCells(
    MultiBlockLattice2D<T, Descriptor> &sourceLattice,
    MultiBlockLattice2D<T, Descriptor> &destinationLattice, Box2D domain)
{
    applyProcessingFunctional(
        new AssignEntireCellFunctional2D<T, Descriptor>, domain, sourceLattice, destinationLattice);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > copyEntireCells(
    MultiBlockLattice2D<T, Descriptor> &lattice)
{
    return copyEntireCells(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > copyEntireCells(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain)
{
    MultiBlockLattice2D<T, Descriptor> *newLattice =
        new MultiBlockLattice2D<T, Descriptor>(lattice, domain);
    copyEntireCells(lattice, *newLattice, domain);
    return std::unique_ptr<MultiBlockLattice2D<T, Descriptor> >(newLattice);
}

}  // namespace plb

#endif  // META_STUFF_WRAPPER_2D_HH
