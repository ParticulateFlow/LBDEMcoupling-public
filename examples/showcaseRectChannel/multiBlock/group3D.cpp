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
 * The 3D group -- implementation file.
 */
#include "multiBlock/group3D.h"

#include "multiBlock/multiBlockGenerator3D.h"

namespace plb {

Group3D::Group3D(MultiBlock3D *block, std::string name)
{
    addNoCheck(block, name);
}

Group3D::~Group3D()
{
    for (pluint i = 0; i < blocks.size(); ++i) {
        delete blocks[i];
    }
}

Group3D::Group3D(Group3D const &rhs) : ids(rhs.ids), names(rhs.names)
{
    blocks.resize(rhs.blocks.size());
    for (pluint i = 0; i < blocks.size(); ++i) {
        blocks[i] = rhs.blocks[i]->clone();
    }
}

Group3D &Group3D::operator=(Group3D const &rhs)
{
    Group3D(rhs).swap(*this);
    return *this;
}

void Group3D::swap(Group3D &rhs)
{
    blocks.swap(rhs.blocks);
    ids.swap(rhs.ids);
    names.swap(rhs.names);
}

plint Group3D::addNoCheck(MultiBlock3D *block, std::string name)
{
    PLB_ASSERT(!hasBlock(name));
    plint blockId = (plint)blocks.size();
    if (name == "") {
        name = "block_" + util::val2str(blockId);
    }
    ids[name] = blockId;
    names[blockId] = name;
    block->periodicity().toggleAll(true);
    blocks.push_back(block);
    return blockId;
}

plint Group3D::add(MultiBlock3D *block, std::string name)
{
    PLB_ASSERT(!hasBlock(name));
    block->periodicity().toggleAll(true);
    plint blockId = (plint)blocks.size();
    if (name == "") {
        name = "block_" + util::val2str(blockId);
    }
    ids[name] = blockId;
    names[blockId] = name;
    if (blocks.empty()) {
        blocks.push_back(block);
    } else {
        MultiBlockManagement3D const &management = blocks[0]->getMultiBlockManagement();
        if (block->getMultiBlockManagement().equivalentTo(management)) {
            blocks.push_back(block);
        } else {
            pcout << "Adjusting parallelization of block \"" << name
                  << "\" to be equal to the one of block " << names[0] << std::endl;
            MultiBlockManagement3D newManagement(management);
            newManagement.changeEnvelopeWidth(block->getMultiBlockManagement().getEnvelopeWidth());
            newManagement.setRefinementLevel(block->getMultiBlockManagement().getRefinementLevel());
            MultiBlock3D *newBlock = block->clone(newManagement);
            delete block;
            blocks.push_back(newBlock);
        }
    }
    return blockId;
}

plint Group3D::getNumBlocks() const
{
    return (plint)blocks.size();
}

std::string Group3D::getName(plint id) const
{
    std::map<plint, std::string>::const_iterator it = names.find(id);
    PLB_ASSERT(it != names.end());
    return it->second;
}

MultiBlock3D &Group3D::get(plint id)
{
    if (id >= (plint)blocks.size()) {
        plbLogicError("Group has no block of ID " + util::val2str(id));
    }
    return *blocks[id];
}

MultiBlock3D &Group3D::get(std::string name)
{
    std::map<std::string, plint>::const_iterator it = ids.find(name);
    if (it == ids.end()) {
        plbLogicError("Group has no block of name " + name);
    }
    plint id = it->second;
    PLB_ASSERT(id < (plint)blocks.size());
    return *blocks[id];
}

bool Group3D::hasBlock(std::string name) const
{
    std::map<std::string, plint>::const_iterator it = ids.find(name);
    return it != ids.end();
}

Box3D Group3D::getBoundingBox() const
{
    // Although groups offer a default constructor for efficiency reasons, you
    // must start by adding a multi-block to the group before anything else.
    PLB_ASSERT(!blocks.empty());
    return blocks[0]->getBoundingBox();
}

MultiBlockManagement3D const &Group3D::getMultiBlockManagement() const
{
    // Although groups offer a default constructor for efficiency reasons, you
    // must start by adding a multi-block to the group before anything else.
    PLB_ASSERT(!blocks.empty());
    return blocks[0]->getMultiBlockManagement();
}

void Group3D::replaceManagement(MultiBlockManagement3D const &management)
{
    for (plint i = 0; i < (plint)blocks.size(); ++i) {
        MultiBlock3D *oldBlock = blocks[i];
        MultiBlockManagement3D const &oldManagement = oldBlock->getMultiBlockManagement();
        MultiBlockManagement3D newManagement(
            management.getSparseBlockStructure(), management.getThreadAttribution().clone(),
            oldManagement.getEnvelopeWidth(), oldManagement.getRefinementLevel());
        blocks[i] = oldBlock->clone(newManagement);
        delete oldBlock;
    }
}

void Group3D::replace(plint id, MultiBlock3D *block)
{
    if (id >= (plint)blocks.size()) {
        plbLogicError("Group has no block of ID " + util::val2str(id));
    }
    delete blocks[id];
    blocks[id] = block;
    for (plint i = 0; i < (plint)blocks.size(); ++i) {
        if (i != id) {
            MultiBlock3D *oldBlock = blocks[i];
            MultiBlockManagement3D const &oldManagement = oldBlock->getMultiBlockManagement();
            MultiBlockManagement3D newManagement(
                block->getMultiBlockManagement().getSparseBlockStructure(),
                block->getMultiBlockManagement().getThreadAttribution().clone(),
                oldManagement.getEnvelopeWidth(), oldManagement.getRefinementLevel());
            blocks[i] = oldBlock->clone(newManagement);
            delete oldBlock;
        }
    }
}

void Group3D::replace(std::string name, MultiBlock3D *block)
{
    std::map<std::string, plint>::const_iterator it = ids.find(name);
    if (it == ids.end()) {
        plbLogicError("Group has no block of name " + name);
    }
    plint id = it->second;
    PLB_ASSERT(id < (plint)blocks.size());
    replace(id, block);
}

plint Group3D::generateContainer(std::string name, plint envelopeWidth, plint gridLevel)
{
    // Although groups offer a default constructor for efficiency reasons, you
    // must start by adding a multi-block to the group before anything else.
    PLB_ASSERT(!blocks.empty());
    MultiBlock3D *block = generateMultiContainerBlock(*blocks[0], envelopeWidth).release();
    block->setRefinementLevel(gridLevel);
    return add(block, name);
}

plint Group3D::generateContainerWithData(
    std::string name, ContainerBlockData *data, plint envelopeWidth, plint gridLevel)
{
    // Although groups offer a default constructor for efficiency reasons, you
    // must start by adding a multi-block to the group before anything else.
    PLB_ASSERT(!blocks.empty());
    MultiBlock3D *block = createContainerBlock(*blocks[0], data, envelopeWidth);
    block->setRefinementLevel(gridLevel);
    return add(block, name);
}

MultiContainerBlock3D &Group3D::getContainer(plint id)
{
    if (id >= (plint)blocks.size()) {
        plbLogicError("Group has no block of ID " + util::val2str(id));
    }
    MultiContainerBlock3D *container = dynamic_cast<MultiContainerBlock3D *>(blocks[id]);
    if (!container) {
        plbLogicError("Block with ID " + util::val2str(id) + " is not a container");
    }
    return *container;
}

MultiContainerBlock3D &Group3D::getContainer(std::string name)
{
    std::map<std::string, plint>::const_iterator it = ids.find(name);
    if (it == ids.end()) {
        plbLogicError("Group has no block of name " + name);
    }
    plint id = it->second;
    PLB_ASSERT(id < (plint)blocks.size());
    MultiContainerBlock3D *container = dynamic_cast<MultiContainerBlock3D *>(blocks[id]);
    if (!container) {
        plbLogicError("Block with name \"" + name + "\" is not a container");
    }
    return *container;
}

}  // namespace plb
