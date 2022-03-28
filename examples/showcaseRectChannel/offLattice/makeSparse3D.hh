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

#ifndef MAKE_SPARSE_3D_HH
#define MAKE_SPARSE_3D_HH

#include "atomicBlock/atomicContainerBlock3D.h"
#include "atomicBlock/reductiveDataProcessingFunctional3D.h"
#include "core/globalDefs.h"
#include "offLattice/domainClustering3D.h"
#include "offLattice/makeSparse3D.h"
#include "parallelism/mpiManager.h"

namespace plb {

struct FlagData3D : public ContainerBlockData {
    bool keepThisBlock;
    virtual FlagData3D *clone() const
    {
        return new FlagData3D(*this);
    }
};

/* ******** ComputeSparsityFunctional3D ************************************ */

template <typename T>
ComputeSparsityFunctional3D<T>::ComputeSparsityFunctional3D() :
    numBlocksId(this->getStatistics().subscribeIntSum())
{ }

template <typename T>
void ComputeSparsityFunctional3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    ScalarField3D<T> *field = dynamic_cast<ScalarField3D<T> *>(blocks[0]);
    AtomicContainerBlock3D *container = dynamic_cast<AtomicContainerBlock3D *>(blocks[1]);
    PLB_ASSERT(field);
    PLB_ASSERT(container);
    bool exclusivelyEliminateCells = true;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (field->get(iX, iY, iZ) != 0) {
                    exclusivelyEliminateCells = false;
                }
            }
        }
    }
    FlagData3D *flagData = new FlagData3D;
    if (exclusivelyEliminateCells) {
        flagData->keepThisBlock = false;
        this->getStatistics().gatherIntSum(numBlocksId, 1);
    } else {
        flagData->keepThisBlock = true;
    }
    container->setData(flagData);
}

template <typename T>
ComputeSparsityFunctional3D<T> *ComputeSparsityFunctional3D<T>::clone() const
{
    return new ComputeSparsityFunctional3D<T>(*this);
}

template <typename T>
void ComputeSparsityFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Scalar Field.
    modified[1] = modif::staticVariables;  // Container Block with flag data.
}

template <typename T>
BlockDomain::DomainT ComputeSparsityFunctional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T>
pluint ComputeSparsityFunctional3D<T>::getNumBlocks() const
{
    return this->getStatistics().getIntSum(numBlocksId);
}

/* ******** computeSparseManagement ************************************ */

template <typename T>
MultiBlockManagement3D computeSparseManagement(MultiScalarField3D<T> &field, plint newEnvelopeWidth)
{
    MultiContainerBlock3D multiFlagBlock(field);
    std::vector<MultiBlock3D *> args;
    args.push_back(&field);
    args.push_back(&multiFlagBlock);
    ComputeSparsityFunctional3D<T> sparsityFunctional;
    applyProcessingFunctional(sparsityFunctional, field.getBoundingBox(), args);

    MultiBlockManagement3D const &management = multiFlagBlock.getMultiBlockManagement();
    ThreadAttribution const &threadAttribution = management.getThreadAttribution();
    SparseBlockStructure3D const &sparseBlock = management.getSparseBlockStructure();

    std::map<plint, Box3D> const &domains = sparseBlock.getBulks();
    std::vector<plint> domainIds(domains.size());
    std::vector<int> keepThisBlock(domains.size());

    std::map<plint, Box3D>::const_iterator it = domains.begin();
    plint pos = 0;
    for (; it != domains.end(); ++it) {
        plint id = it->first;
        domainIds[pos] = id;
        if (threadAttribution.isLocal(id)) {
            AtomicContainerBlock3D const &flagBlock = multiFlagBlock.getComponent(id);
            FlagData3D const *data = dynamic_cast<FlagData3D const *>(flagBlock.getData());
            PLB_ASSERT(data);
            keepThisBlock[pos] = data->keepThisBlock ? 1 : 0;
        } else {
            keepThisBlock[pos] = 0;
        }
        ++pos;
    }

#ifdef PLB_MPI_PARALLEL
    std::vector<int> tmp(keepThisBlock.size());
    global::mpi().reduceVect(keepThisBlock, tmp, MPI_SUM);
    global::mpi().bCast(&tmp[0], tmp.size());
    tmp.swap(keepThisBlock);
#endif

    SparseBlockStructure3D newSparseBlock(field.getBoundingBox());
    plint newId = 0;
    for (pluint iBlock = 0; iBlock < keepThisBlock.size(); ++iBlock) {
        if (keepThisBlock[iBlock]) {
            plint id = domainIds[iBlock];
            Box3D bulk, uniqueBulk;
            sparseBlock.getBulk(id, bulk);
            sparseBlock.getUniqueBulk(id, uniqueBulk);
            newSparseBlock.addBlock(bulk, uniqueBulk, newId++);
        }
    }
    // If this assertion fails, that means that the domain covered
    // by the sparse block-structure is empty.
    PLB_ASSERT(newId > 0);

    ExplicitThreadAttribution *newAttribution = new ExplicitThreadAttribution;
    std::vector<std::pair<plint, plint> > ranges;
    plint numRanges = std::min(newId, (plint)global::mpi().getSize());
    util::linearRepartition(0, newId - 1, numRanges, ranges);

    for (pluint iProc = 0; iProc < ranges.size(); ++iProc) {
        for (plint blockId = ranges[iProc].first; blockId <= ranges[iProc].second; ++blockId) {
            newAttribution->addBlock(blockId, iProc);
        }
    }

    MultiBlockManagement3D newManagement(
        newSparseBlock, newAttribution, newEnvelopeWidth, management.getRefinementLevel());
    return newManagement;
}

}  // namespace plb

#endif  // MAKE_SPARSE_3D_HH
