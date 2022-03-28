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
#ifndef META_STUFF_WRAPPER_3D_HH
#define META_STUFF_WRAPPER_3D_HH

#include <numeric>
#include <set>

#include "core/dynamicsIdentifiers.h"
#include "dataProcessors/dataInitializerWrapper3D.h"
#include "dataProcessors/metaStuffFunctional3D.h"
#include "dataProcessors/metaStuffWrapper3D.h"
#include "multiBlock/multiContainerBlock3D.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
void extractTopMostDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<int> &dynamicsId, Box3D domain)
{
    applyProcessingFunctional(
        new ExtractTopMostDynamicsFunctional3D<T, Descriptor>(), domain, lattice, dynamicsId);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<int> > extractTopMostDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    MultiScalarField3D<int> *dynamicsId = new MultiScalarField3D<int>(lattice, domain);
    extractTopMostDynamics(lattice, *dynamicsId, domain);
    return std::unique_ptr<MultiScalarField3D<int> >(dynamicsId);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<int> > extractTopMostDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return extractTopMostDynamics(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void extractBottomMostDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<int> &dynamicsId, Box3D domain)
{
    applyProcessingFunctional(
        new ExtractBottomMostDynamicsFunctional3D<T, Descriptor>(), domain, lattice, dynamicsId);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<int> > extractBottomMostDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    MultiScalarField3D<int> *dynamicsId = new MultiScalarField3D<int>(lattice, domain);
    extractBottomMostDynamics(lattice, *dynamicsId, domain);
    return std::unique_ptr<MultiScalarField3D<int> >(dynamicsId);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<int> > extractBottomMostDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return extractBottomMostDynamics(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void uniqueDynamicsChains(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain,
    std::vector<std::vector<int> > &chains, pluint &maxChainLength)
{
    MultiContainerBlock3D container(lattice);
    std::vector<MultiBlock3D *> latticeAndContainer;
    latticeAndContainer.push_back(&lattice);
    latticeAndContainer.push_back(&container);

    StoreDynamicsFunctional3D<T, Descriptor> SDFunctional;
    applyProcessingFunctional(SDFunctional, domain, latticeAndContainer);
    maxChainLength = SDFunctional.getMaxChainLength();

    MultiBlockManagement3D const &management = container.getMultiBlockManagement();
    ThreadAttribution const &threadAttribution = management.getThreadAttribution();
    SparseBlockStructure3D const &sparseBlock = management.getSparseBlockStructure();

    std::vector<int> localSerializedChains;
    std::vector<plint> localBlocks = sparseBlock.getLocalBlocks(threadAttribution);
    for (pluint i = 0; i < localBlocks.size(); ++i) {
        plint blockId = localBlocks[i];
        AtomicContainerBlock3D &atomicContainer = container.getComponent(blockId);
        StoreDynamicsID *storeId = dynamic_cast<StoreDynamicsID *>(atomicContainer.getData());
        PLB_ASSERT(storeId);
        storeId->startIterations();
        while (!storeId->empty()) {
            localSerializedChains.insert(
                localSerializedChains.end(), storeId->getCurrent().begin(),
                storeId->getCurrent().end());
            storeId->iterate();
        }
        localSerializedChains.push_back(-1);
    }

#ifdef PLB_MPI_PARALLEL
    std::vector<plint> idsPerProc(global::mpi().getSize(), 0);
    idsPerProc[global::mpi().getRank()] = (plint)localSerializedChains.size();
    global::mpi().allReduceVect(idsPerProc, MPI_SUM);

    std::vector<plint> sumIdsPerProc(idsPerProc.size() + 1);
    sumIdsPerProc[0] = 0;
    std::partial_sum(idsPerProc.begin(), idsPerProc.end(), sumIdsPerProc.begin() + 1);
    plint totalNumIds = sumIdsPerProc.back();

    std::vector<int> allIds(totalNumIds, 0);
    plint start = sumIdsPerProc[global::mpi().getRank()];
    for (pluint i = 0; i < localSerializedChains.size(); ++i) {
        allIds[start + i] = localSerializedChains[i];
    }
    global::mpi().allReduceVect(allIds, MPI_SUM);
#else
    std::vector<int> allIds(localSerializedChains);
#endif

    std::set<std::vector<int>, VectorIsLess> allChains;
    std::vector<int> nextChain;
    for (pluint i = 0; i < allIds.size(); ++i) {
        if (allIds[i] == -1) {
            allChains.insert(nextChain);
            nextChain.clear();
        } else {
            nextChain.push_back(allIds[i]);
        }
    }

    chains.assign(allChains.begin(), allChains.end());
}

template <typename T, template <typename U> class Descriptor>
void uniqueDynamicsIds(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, std::vector<int> &ids)
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
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiScalarField3D<int> &dynamicsId,
    std::map<int, std::string> &nameOfDynamics, Box3D domain)
{
    std::vector<std::vector<int> > chains;
    pluint maxChainLength;
    uniqueDynamicsChains(lattice, domain, chains, maxChainLength);
    nameOfDynamics.clear();
    typename ExtractDynamicsChainFunctional3D<T, Descriptor>::DMap dynamicsMap;
    for (pluint iChain = 0; iChain < chains.size(); ++iChain) {
        dynamicsMap[chains[iChain]] = iChain;
        nameOfDynamics[iChain] = meta::constructIdNameChain<T, Descriptor>(chains[iChain], " >> ");
    }
    setToConstant(dynamicsId, dynamicsId.getBoundingBox(), -1);
    applyProcessingFunctional(
        new ExtractDynamicsChainFunctional3D<T, Descriptor>(dynamicsMap, maxChainLength), domain,
        lattice, dynamicsId);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<int> > extractDynamicsChain(
    MultiBlockLattice3D<T, Descriptor> &lattice, std::map<int, std::string> &nameOfDynamics,
    Box3D domain)
{
    MultiScalarField3D<int> *dynamicsId = new MultiScalarField3D<int>(lattice, domain);
    extractDynamicsChain(lattice, *dynamicsId, nameOfDynamics, domain);
    return std::unique_ptr<MultiScalarField3D<int> >(dynamicsId);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<int> > extractDynamicsChain(
    MultiBlockLattice3D<T, Descriptor> &lattice, std::map<int, std::string> &nameOfDynamics)
{
    return extractDynamicsChain(lattice, nameOfDynamics, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void copyEntireCells(
    MultiBlockLattice3D<T, Descriptor> &sourceLattice,
    MultiBlockLattice3D<T, Descriptor> &destinationLattice, Box3D domain)
{
    applyProcessingFunctional(
        new AssignEntireCellFunctional3D<T, Descriptor>, domain, sourceLattice, destinationLattice);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > copyEntireCells(
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    return copyEntireCells(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > copyEntireCells(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    MultiBlockLattice3D<T, Descriptor> *newLattice =
        new MultiBlockLattice3D<T, Descriptor>(lattice, domain);
    copyEntireCells(lattice, *newLattice, domain);
    return std::unique_ptr<MultiBlockLattice3D<T, Descriptor> >(newLattice);
}

}  // namespace plb

#endif  // META_STUFF_WRAPPER_3D_HH
