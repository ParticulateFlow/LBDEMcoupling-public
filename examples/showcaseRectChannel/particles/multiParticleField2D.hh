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

#ifndef MULTI_PARTICLE_FIELD_2D_HH
#define MULTI_PARTICLE_FIELD_2D_HH

#include <memory>

#include "core/globalDefs.h"
#include "particles/multiParticleField2D.h"
#include "particles/particleNonLocalTransfer2D.h"

namespace plb {

template <class ParticleFieldT>
std::unique_ptr<MultiParticleField2D<ParticleFieldT> > defaultGenerateParticleField2D(
    MultiBlockManagement2D const &management, plint unnamedDummyArg)
{
    return std::unique_ptr<MultiParticleField2D<ParticleFieldT> >(
        new MultiParticleField2D<ParticleFieldT>(
            management, defaultMultiBlockPolicy2D().getCombinedStatistics()));
}

template <class ParticleFieldT>
const int MultiParticleField2D<ParticleFieldT>::staticId = meta::registerMultiBlock2D(
    MultiParticleField2D<ParticleFieldT>::basicType(),
    MultiParticleField2D<ParticleFieldT>::descriptorType(),
    MultiParticleField2D<ParticleFieldT>::blockName(),
    defaultGenerateParticleField2D<ParticleFieldT>);

template <class ParticleFieldT>
MultiParticleField2D<ParticleFieldT>::MultiParticleField2D(
    MultiBlockManagement2D const &multiBlockManagement_, CombinedStatistics *combinedStatistics_)

    :
    MultiBlock2D(
        multiBlockManagement_, defaultMultiBlockPolicy2D().getBlockCommunicator(),
        combinedStatistics_)
{
    allocateBlocks();
}

template <class ParticleFieldT>
MultiParticleField2D<ParticleFieldT>::~MultiParticleField2D()
{
    deAllocateBlocks();
}

template <class ParticleFieldT>
MultiParticleField2D<ParticleFieldT>::MultiParticleField2D(plint nx_, plint ny_) :
    MultiBlock2D(
        // Default envelope-width to 1
        defaultMultiBlockPolicy2D().getMultiBlockManagement(nx_, ny_, 1),
        defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics())
{
    allocateBlocks();
}

template <class ParticleFieldT>
MultiParticleField2D<ParticleFieldT>::MultiParticleField2D(MultiBlock2D const &rhs) :
    MultiBlock2D(rhs)
{
    allocateBlocks();
}

template <class ParticleFieldT>
MultiParticleField2D<ParticleFieldT>::MultiParticleField2D(
    MultiBlock2D const &rhs, Box2D subDomain, bool crop) :
    MultiBlock2D(
        intersect(rhs.getMultiBlockManagement(), subDomain, crop),
        rhs.getBlockCommunicator().clone(), rhs.getCombinedStatistics().clone())
{
    allocateBlocks();
}

template <class ParticleFieldT>
MultiParticleField2D<ParticleFieldT>::MultiParticleField2D(
    MultiParticleField2D<ParticleFieldT> const &rhs) :
    MultiBlock2D(rhs)
{
    allocateBlocks();
    typename BlockMap::iterator it = blocks.begin();
    typename BlockMap::const_iterator rhsIt = rhs.blocks.begin();

    for (; it != blocks.end(); ++it, ++rhsIt) {
        *(it->second) = *(rhsIt->second);
    }
}

template <class ParticleFieldT>
MultiParticleField2D<ParticleFieldT> &MultiParticleField2D<ParticleFieldT>::operator=(
    MultiParticleField2D<ParticleFieldT> const &rhs)
{
    MultiParticleField2D<ParticleFieldT> tmp(rhs);
    swap(tmp);
    return *this;
}

template <class ParticleFieldT>
void MultiParticleField2D<ParticleFieldT>::swap(MultiParticleField2D<ParticleFieldT> &rhs)
{
    blocks.swap(rhs.blocks);
    MultiBlock2D::swap(rhs);
}

template <class ParticleFieldT>
MultiParticleField2D<ParticleFieldT> *MultiParticleField2D<ParticleFieldT>::clone() const
{
    return new MultiParticleField2D<ParticleFieldT>(*this);
}

template <class ParticleFieldT>
MultiParticleField2D<ParticleFieldT> *MultiParticleField2D<ParticleFieldT>::clone(
    MultiBlockManagement2D const &newManagement) const
{
    MultiParticleField2D<ParticleFieldT> *newField = new MultiParticleField2D<ParticleFieldT>(
        newManagement, this->getCombinedStatistics().clone());
    copy(*this, this->getBoundingBox(), *newField, newField->getBoundingBox());
    return newField;
}

template <class ParticleFieldT>
void MultiParticleField2D<ParticleFieldT>::allocateBlocks()
{
    for (pluint iBlock = 0; iBlock < this->getLocalInfo().getBlocks().size(); ++iBlock) {
        plint blockId = this->getLocalInfo().getBlocks()[iBlock];
        SmartBulk2D bulk(this->getMultiBlockManagement(), blockId);
        Box2D envelope = bulk.computeEnvelope();
        ParticleFieldT *newBlock = new ParticleFieldT(envelope.getNx(), envelope.getNy());
        newBlock->setLocation(Dot2D(envelope.x0, envelope.y0));
        blocks[blockId] = newBlock;
    }
}

template <class ParticleFieldT>
void MultiParticleField2D<ParticleFieldT>::deAllocateBlocks()
{
    for (typename BlockMap::iterator it = blocks.begin(); it != blocks.end(); ++it) {
        delete it->second;
    }
}

template <class ParticleFieldT>
ParticleFieldT &MultiParticleField2D<ParticleFieldT>::getComponent(plint blockId)
{
    typename BlockMap::iterator it = blocks.find(blockId);
    PLB_ASSERT(it != blocks.end());
    return *it->second;
}

template <class ParticleFieldT>
ParticleFieldT const &MultiParticleField2D<ParticleFieldT>::getComponent(plint blockId) const
{
    typename BlockMap::const_iterator it = blocks.find(blockId);
    PLB_ASSERT(it != blocks.end());
    return *it->second;
}

template <class ParticleFieldT>
plint MultiParticleField2D<ParticleFieldT>::sizeOfCell() const
{
    /// Particles are dynamic objects and have no static cell size.
    return 0;
}

template <class ParticleFieldT>
plint MultiParticleField2D<ParticleFieldT>::getCellDim() const
{
    /// Particles are dynamic objects and have no well-defined dimensionality.
    return 0;
}

template <class ParticleFieldT>
int MultiParticleField2D<ParticleFieldT>::getStaticId() const
{
    return staticId;
}

// TODO: whichData is not used. Should we remove it?
template <class ParticleFieldT>
void MultiParticleField2D<ParticleFieldT>::copyReceive(
    MultiBlock2D const &fromBlock, Box2D const &fromDomain, Box2D const &toDomain,
    modif::ModifT whichData)
{
    MultiParticleField2D<ParticleFieldT> const *fromField =
        dynamic_cast<MultiParticleField2D<ParticleFieldT> const *>(&fromBlock);
    PLB_ASSERT(fromField);
    copy(*fromField, fromDomain, *this, toDomain);
}

template <class ParticleFieldT>
std::string MultiParticleField2D<ParticleFieldT>::getBlockName() const
{
    return blockName();
}

template <class ParticleFieldT>
std::vector<std::string> MultiParticleField2D<ParticleFieldT>::getTypeInfo() const
{
    std::vector<std::string> answer;
    answer.push_back(basicType());
    answer.push_back(descriptorType());
    return answer;
}

template <class ParticleFieldT>
std::string MultiParticleField2D<ParticleFieldT>::blockName()
{
    return ParticleFieldT::getBlockName();
}

template <class ParticleFieldT>
std::string MultiParticleField2D<ParticleFieldT>::basicType()
{
    return ParticleFieldT::basicType();
}

template <class ParticleFieldT>
std::string MultiParticleField2D<ParticleFieldT>::descriptorType()
{
    return ParticleFieldT::descriptorType();
}

template <class ParticleFieldT>
MultiParticleField2D<ParticleFieldT> &findMultiParticleField2D(id_t id)
{
    MultiBlock2D *multiBlock = multiBlockRegistration2D().find(id);
    if (!multiBlock || multiBlock->getStaticId() != MultiParticleField2D<ParticleFieldT>::staticId)
    {
        throw PlbLogicException("Trying to access a multi block lattice that is not registered.");
    }
    return (MultiParticleField2D<ParticleFieldT> &)(*multiBlock);
}

}  // namespace plb

#endif  // MULTI_PARTICLE_FIELD_2D_H
