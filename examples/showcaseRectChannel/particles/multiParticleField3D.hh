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

#ifndef MULTI_PARTICLE_FIELD_3D_HH
#define MULTI_PARTICLE_FIELD_3D_HH

#include <memory>

#include "core/globalDefs.h"
#include "core/multiBlockIdentifiers3D.h"
#include "multiBlock/defaultMultiBlockPolicy3D.h"
#include "particles/multiParticleField3D.h"
#include "particles/particleNonLocalTransfer3D.h"

namespace plb {

// TODO:: Why this unused unnamedDummyArg?
template <class ParticleFieldT>
std::unique_ptr<MultiParticleField3D<ParticleFieldT> > defaultGenerateParticleField3D(
    MultiBlockManagement3D const &management, plint unnamedDummyArg)
{
    return std::unique_ptr<MultiParticleField3D<ParticleFieldT> >(
        new MultiParticleField3D<ParticleFieldT>(
            management, defaultMultiBlockPolicy3D().getCombinedStatistics()));
}

template <class ParticleFieldT>
const int MultiParticleField3D<ParticleFieldT>::staticId = meta::registerMultiBlock3D(
    MultiParticleField3D<ParticleFieldT>::basicType(),
    MultiParticleField3D<ParticleFieldT>::descriptorType(),
    MultiParticleField3D<ParticleFieldT>::blockName(),
    defaultGenerateParticleField3D<ParticleFieldT>);

template <class ParticleFieldT>
MultiParticleField3D<ParticleFieldT>::MultiParticleField3D(
    MultiBlockManagement3D const &multiBlockManagement_, CombinedStatistics *combinedStatistics_)

    :
    MultiBlock3D(
        multiBlockManagement_, defaultMultiBlockPolicy3D().getBlockCommunicator(),
        combinedStatistics_)
{
    allocateBlocks();
}

template <class ParticleFieldT>
MultiParticleField3D<ParticleFieldT>::~MultiParticleField3D()
{
    deAllocateBlocks();
}

template <class ParticleFieldT>
MultiParticleField3D<ParticleFieldT>::MultiParticleField3D(plint nx_, plint ny_, plint nz_) :
    MultiBlock3D(
        // Default envelope-width to 1
        defaultMultiBlockPolicy3D().getMultiBlockManagement(nx_, ny_, nz_, 1),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics())
{
    allocateBlocks();
}

template <class ParticleFieldT>
MultiParticleField3D<ParticleFieldT>::MultiParticleField3D(MultiBlock3D const &rhs) :
    MultiBlock3D(rhs)
{
    allocateBlocks();
}

template <class ParticleFieldT>
MultiParticleField3D<ParticleFieldT>::MultiParticleField3D(
    MultiBlock3D const &rhs, Box3D subDomain, bool crop) :
    MultiBlock3D(
        intersect(rhs.getMultiBlockManagement(), subDomain, crop),
        rhs.getBlockCommunicator().clone(), rhs.getCombinedStatistics().clone())
{
    allocateBlocks();
}

template <class ParticleFieldT>
MultiParticleField3D<ParticleFieldT>::MultiParticleField3D(
    MultiParticleField3D<ParticleFieldT> const &rhs) :
    MultiBlock3D(rhs)
{
    allocateBlocks();
    typename BlockMap::iterator it = blocks.begin();
    typename BlockMap::const_iterator rhsIt = rhs.blocks.begin();

    for (; it != blocks.end(); ++it, ++rhsIt) {
        *(it->second) = *(rhsIt->second);
    }
}

template <class ParticleFieldT>
MultiParticleField3D<ParticleFieldT> &MultiParticleField3D<ParticleFieldT>::operator=(
    MultiParticleField3D<ParticleFieldT> const &rhs)
{
    MultiParticleField3D<ParticleFieldT> tmp(rhs);
    swap(tmp);
    return *this;
}

template <class ParticleFieldT>
void MultiParticleField3D<ParticleFieldT>::swap(MultiParticleField3D<ParticleFieldT> &rhs)
{
    blocks.swap(rhs.blocks);
    MultiBlock3D::swap(rhs);
}

template <class ParticleFieldT>
MultiParticleField3D<ParticleFieldT> *MultiParticleField3D<ParticleFieldT>::clone() const
{
    return new MultiParticleField3D<ParticleFieldT>(*this);
}

template <class ParticleFieldT>
MultiParticleField3D<ParticleFieldT> *MultiParticleField3D<ParticleFieldT>::clone(
    MultiBlockManagement3D const &newManagement) const
{
    MultiParticleField3D<ParticleFieldT> *newField = new MultiParticleField3D<ParticleFieldT>(
        newManagement, this->getCombinedStatistics().clone());
    copy(*this, this->getBoundingBox(), *newField, newField->getBoundingBox());
    return newField;
}

template <class ParticleFieldT>
void MultiParticleField3D<ParticleFieldT>::allocateBlocks()
{
    for (pluint iBlock = 0; iBlock < this->getLocalInfo().getBlocks().size(); ++iBlock) {
        plint blockId = this->getLocalInfo().getBlocks()[iBlock];
        SmartBulk3D bulk(this->getMultiBlockManagement(), blockId);
        Box3D envelope = bulk.computeEnvelope();
        ParticleFieldT *newBlock =
            new ParticleFieldT(envelope.getNx(), envelope.getNy(), envelope.getNz());
        newBlock->setLocation(Dot3D(envelope.x0, envelope.y0, envelope.z0));
        blocks[blockId] = newBlock;
    }
}

template <class ParticleFieldT>
void MultiParticleField3D<ParticleFieldT>::deAllocateBlocks()
{
    for (typename BlockMap::iterator it = blocks.begin(); it != blocks.end(); ++it) {
        delete it->second;
    }
}

template <class ParticleFieldT>
ParticleFieldT &MultiParticleField3D<ParticleFieldT>::getComponent(plint blockId)
{
    typename BlockMap::iterator it = blocks.find(blockId);
    PLB_ASSERT(it != blocks.end());
    return *it->second;
}

template <class ParticleFieldT>
ParticleFieldT const &MultiParticleField3D<ParticleFieldT>::getComponent(plint blockId) const
{
    typename BlockMap::const_iterator it = blocks.find(blockId);
    PLB_ASSERT(it != blocks.end());
    return *it->second;
}

template <class ParticleFieldT>
plint MultiParticleField3D<ParticleFieldT>::sizeOfCell() const
{
    /// Particles are dynamic objects and have no static cell size.
    return 0;
}

template <class ParticleFieldT>
plint MultiParticleField3D<ParticleFieldT>::getCellDim() const
{
    /// Particles are dynamic objects and have no well-defined dimensionality.
    return 0;
}

template <class ParticleFieldT>
int MultiParticleField3D<ParticleFieldT>::getStaticId() const
{
    return staticId;
}

template <class ParticleFieldT>
void MultiParticleField3D<ParticleFieldT>::copyReceive(
    MultiBlock3D const &fromBlock, Box3D const &fromDomain, Box3D const &toDomain,
    [[maybe_unused]] modif::ModifT whichData)
{
    MultiParticleField3D<ParticleFieldT> const *fromField =
        dynamic_cast<MultiParticleField3D<ParticleFieldT> const *>(&fromBlock);
    PLB_ASSERT(fromField);
    copy(*fromField, fromDomain, *this, toDomain);
}

template <class ParticleFieldT>
std::string MultiParticleField3D<ParticleFieldT>::getBlockName() const
{
    return blockName();
}

template <class ParticleFieldT>
std::vector<std::string> MultiParticleField3D<ParticleFieldT>::getTypeInfo() const
{
    std::vector<std::string> answer;
    answer.push_back(basicType());
    answer.push_back(descriptorType());
    return answer;
}

template <class ParticleFieldT>
std::string MultiParticleField3D<ParticleFieldT>::blockName()
{
    return ParticleFieldT::getBlockName();
}

template <class ParticleFieldT>
std::string MultiParticleField3D<ParticleFieldT>::basicType()
{
    return ParticleFieldT::basicType();
}

template <class ParticleFieldT>
std::string MultiParticleField3D<ParticleFieldT>::descriptorType()
{
    return ParticleFieldT::descriptorType();
}

template <class ParticleFieldT>
MultiParticleField3D<ParticleFieldT> &findMultiParticleField3D(id_t id)
{
    MultiBlock3D *multiBlock = multiBlockRegistration3D().find(id);
    if (!multiBlock || multiBlock->getStaticId() != MultiParticleField3D<ParticleFieldT>::staticId)
    {
        throw PlbLogicException("Trying to access a multi block lattice that is not registered.");
    }
    return (MultiParticleField3D<ParticleFieldT> &)(*multiBlock);
}

}  // namespace plb

#endif  // MULTI_PARTICLE_FIELD_3D_H
