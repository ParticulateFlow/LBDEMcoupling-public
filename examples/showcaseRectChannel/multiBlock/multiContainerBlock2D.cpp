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
 * Serial implementation of scalar, vector and tensor fields for 2D data analysis.
 * -- header file
 */

#include "multiBlock/multiContainerBlock2D.h"

#include "core/blockIdentifiers.h"
#include "core/globalDefs.h"
#include "multiBlock/defaultMultiBlockPolicy2D.h"

namespace plb {

MultiContainerBlock2D::MultiContainerBlock2D(
    MultiBlockManagement2D const &multiBlockManagement_, CombinedStatistics *combinedStatistics_)

    :
    MultiBlock2D(
        multiBlockManagement_, defaultMultiBlockPolicy2D().getBlockCommunicator(),
        combinedStatistics_)
{
    allocateBlocks();
}

MultiContainerBlock2D::~MultiContainerBlock2D()
{
    deAllocateBlocks();
}

MultiContainerBlock2D::MultiContainerBlock2D(plint nx_, plint ny_) :
    MultiBlock2D(
        // Default envelope-width to 0
        defaultMultiBlockPolicy2D().getMultiBlockManagement(nx_, ny_, 0),
        defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics())
{
    allocateBlocks();
}

MultiContainerBlock2D::MultiContainerBlock2D(MultiBlock2D const &rhs) : MultiBlock2D(rhs)
{
    allocateBlocks();
}

MultiContainerBlock2D::MultiContainerBlock2D(MultiBlock2D const &rhs, Box2D subDomain, bool crop) :
    MultiBlock2D(
        intersect(rhs.getMultiBlockManagement(), subDomain, crop),
        rhs.getBlockCommunicator().clone(), rhs.getCombinedStatistics().clone())
{
    allocateBlocks();
}

MultiContainerBlock2D::MultiContainerBlock2D(MultiContainerBlock2D const &rhs) : MultiBlock2D(rhs)
{
    allocateBlocks();
}

MultiContainerBlock2D &MultiContainerBlock2D::operator=(MultiContainerBlock2D const &rhs)
{
    MultiContainerBlock2D tmp(rhs);
    swap(tmp);
    return *this;
}

void MultiContainerBlock2D::swap(MultiContainerBlock2D &rhs)
{
    blocks.swap(rhs.blocks);
    MultiBlock2D::swap(rhs);
}

MultiContainerBlock2D *MultiContainerBlock2D::clone() const
{
    return new MultiContainerBlock2D(*this);
}

MultiContainerBlock2D *MultiContainerBlock2D::clone(
    MultiBlockManagement2D const &multiBlockManagement) const
{
    // By definition, a multi container block cannot be redistributed over
    //   a different block arrangement. Consequently, this function
    //   default to plain clone().
    return clone();
}

void MultiContainerBlock2D::allocateBlocks()
{
    for (pluint iBlock = 0; iBlock < this->getLocalInfo().getBlocks().size(); ++iBlock) {
        plint blockId = this->getLocalInfo().getBlocks()[iBlock];
        SmartBulk2D bulk(this->getMultiBlockManagement(), blockId);
        Box2D envelope = bulk.computeEnvelope();
        AtomicContainerBlock2D *newBlock =
            new AtomicContainerBlock2D(envelope.getNx(), envelope.getNy());
        newBlock->setLocation(Dot2D(envelope.x0, envelope.y0));
        blocks[blockId] = newBlock;
    }
}

void MultiContainerBlock2D::allocateBlocks(MultiContainerBlock2D const &rhs)
{
    for (pluint iBlock = 0; iBlock < this->getLocalInfo().getBlocks().size(); ++iBlock) {
        plint blockId = this->getLocalInfo().getBlocks()[iBlock];
        BlockMap::const_iterator it(rhs.blocks.find(blockId));
        PLB_ASSERT(it != rhs.blocks.end());
        AtomicContainerBlock2D *newBlock = new AtomicContainerBlock2D(*it->second);
        blocks[blockId] = newBlock;
    }
}

void MultiContainerBlock2D::deAllocateBlocks()
{
    for (BlockMap::iterator it = blocks.begin(); it != blocks.end(); ++it) {
        delete it->second;
    }
}

AtomicContainerBlock2D &MultiContainerBlock2D::getComponent(plint blockId)
{
    BlockMap::iterator it = blocks.find(blockId);
    PLB_ASSERT(it != blocks.end());
    return *it->second;
}

AtomicContainerBlock2D const &MultiContainerBlock2D::getComponent(plint blockId) const
{
    BlockMap::const_iterator it = blocks.find(blockId);
    PLB_ASSERT(it != blocks.end());
    return *it->second;
}

plint MultiContainerBlock2D::sizeOfCell() const
{
    return 0;
}

plint MultiContainerBlock2D::getCellDim() const
{
    return 0;
}

int MultiContainerBlock2D::getStaticId() const
{
    return 0;
}

void MultiContainerBlock2D::copyReceive(
    MultiBlock2D const &fromBlock, Box2D const &fromDomain, Box2D const &toDomain,
    modif::ModifT whichData)
{
    PLB_ASSERT(false);
}

std::string MultiContainerBlock2D::getBlockName() const
{
    return std::string("ContainerBlock2D");
}

std::vector<std::string> MultiContainerBlock2D::getTypeInfo() const
{
    return std::vector<std::string>();
}

MultiContainerBlock2D *createContainerBlock(MultiBlock2D &templ, ContainerBlockData *data)
{
    MultiContainerBlock2D *dataContainer = new MultiContainerBlock2D(templ);

    MultiBlockManagement2D const &management = templ.getMultiBlockManagement();
    ThreadAttribution const &threadAttribution = management.getThreadAttribution();
    SparseBlockStructure2D const &sparseBlock = management.getSparseBlockStructure();
    std::map<plint, Box2D> const &domains = sparseBlock.getBulks();

    std::map<plint, Box2D>::const_iterator it = domains.begin();
    plint pos = 0;
    for (; it != domains.end(); ++it, ++pos) {
        plint id = it->first;
        if (threadAttribution.isLocal(id)) {
            AtomicContainerBlock2D &atomicDataContainer = dataContainer->getComponent(id);
            ContainerBlockData *nextData = data->clone();
            nextData->setUniqueID(pos);
            atomicDataContainer.setData(nextData);
        }
    }

    delete data;
    return dataContainer;
}

}  // namespace plb
