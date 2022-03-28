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
 * The 2D group -- implementation file.
 */
#include "multiBlock/group2D.h"

#include "multiBlock/multiBlockGenerator2D.h"

namespace plb {

Group2D::Group2D(MultiBlock2D *block, std::string name)
{
    add(block, name);
}

Group2D::~Group2D()
{
    for (pluint i = 0; i < blocks.size(); ++i) {
        delete blocks[i];
    }
}

Group2D::Group2D(Group2D const &rhs) : ids(rhs.ids), names(rhs.names)
{
    blocks.resize(rhs.blocks.size());
    for (pluint i = 0; i < blocks.size(); ++i) {
        blocks[i] = rhs.blocks[i]->clone();
    }
}

Group2D &Group2D::operator=(Group2D const &rhs)
{
    Group2D(rhs).swap(*this);
    return *this;
}

void Group2D::swap(Group2D &rhs)
{
    blocks.swap(rhs.blocks);
    ids.swap(rhs.ids);
    names.swap(rhs.names);
}

plint Group2D::add(MultiBlock2D *block, std::string name)
{
    PLB_ASSERT(!hasBlock(name));
    plint blockId = (plint)blocks.size();
    if (name == "") {
        name = "block_" + util::val2str(blockId);
    }
    ids[name] = blockId;
    names[blockId] = name;
    if (blocks.empty()) {
        blocks.push_back(block);
    } else {
        MultiBlockManagement2D const &management = blocks[0]->getMultiBlockManagement();
        if (block->getMultiBlockManagement().equivalentTo(management)) {
            blocks.push_back(block);
        } else {
            pcout << "Adjusting parallelization of block \"" << name
                  << "\" to be equal to the one of block " << names[0] << std::endl;
            MultiBlockManagement2D newManagement(management);
            newManagement.changeEnvelopeWidth(block->getMultiBlockManagement().getEnvelopeWidth());
            newManagement.setRefinementLevel(block->getMultiBlockManagement().getRefinementLevel());
            MultiBlock2D *newBlock = block->clone(newManagement);
            delete block;
            blocks.push_back(newBlock);
        }
    }
    return blockId;
}

plint Group2D::getNumBlocks() const
{
    return (plint)blocks.size();
}

std::string Group2D::getName(plint id) const
{
    std::map<plint, std::string>::const_iterator it = names.find(id);
    PLB_ASSERT(it != names.end());
    return it->second;
}

MultiBlock2D &Group2D::get(plint id)
{
    PLB_ASSERT(id < (plint)blocks.size());
    return *blocks[id];
}

MultiBlock2D &Group2D::get(std::string name)
{
    std::map<std::string, plint>::const_iterator it = ids.find(name);
    PLB_ASSERT(it != ids.end());
    plint id = it->second;
    PLB_ASSERT(id < (plint)blocks.size());
    return *blocks[id];
}

bool Group2D::hasBlock(std::string name) const
{
    std::map<std::string, plint>::const_iterator it = ids.find(name);
    return it != ids.end();
}

Box2D Group2D::getBoundingBox() const
{
    // Although groups offer a default constructor for efficiency reasons, you
    // must start by adding a multi-block to the group before anything else.
    PLB_ASSERT(!blocks.empty());
    return blocks[0]->getBoundingBox();
}

MultiBlockManagement2D const &Group2D::getMultiBlockManagement() const
{
    // Although groups offer a default constructor for efficiency reasons, you
    // must start by adding a multi-block to the group before anything else.
    PLB_ASSERT(!blocks.empty());
    return blocks[0]->getMultiBlockManagement();
}

void Group2D::replaceManagement(MultiBlockManagement2D const &management)
{
    for (plint i = 0; i < (plint)blocks.size(); ++i) {
        MultiBlock2D *oldBlock = blocks[i];
        MultiBlockManagement2D const &oldManagement = oldBlock->getMultiBlockManagement();
        MultiBlockManagement2D newManagement(
            management.getSparseBlockStructure(), management.getThreadAttribution().clone(),
            oldManagement.getEnvelopeWidth(), oldManagement.getRefinementLevel());
        blocks[i] = oldBlock->clone(newManagement);
        delete oldBlock;
    }
}

void Group2D::replace(plint id, MultiBlock2D *block)
{
    if (id >= (plint)blocks.size()) {
        plbLogicError("Group has no block of ID " + util::val2str(id));
    }
    delete blocks[id];
    blocks[id] = block;
    for (plint i = 0; i < (plint)blocks.size(); ++i) {
        if (i != id) {
            MultiBlock2D *oldBlock = blocks[i];
            MultiBlockManagement2D const &oldManagement = oldBlock->getMultiBlockManagement();
            MultiBlockManagement2D newManagement(
                block->getMultiBlockManagement().getSparseBlockStructure(),
                block->getMultiBlockManagement().getThreadAttribution().clone(),
                oldManagement.getEnvelopeWidth(), oldManagement.getRefinementLevel());
            blocks[i] = oldBlock->clone(newManagement);
            delete oldBlock;
        }
    }
}

void Group2D::replace(std::string name, MultiBlock2D *block)
{
    std::map<std::string, plint>::const_iterator it = ids.find(name);
    if (it == ids.end()) {
        plbLogicError("Group has no block of name " + name);
    }
    plint id = it->second;
    PLB_ASSERT(id < (plint)blocks.size());
    replace(id, block);
}

plint Group2D::generateContainer(std::string name, plint envelopeWidth, plint gridLevel)
{
    // Although groups offer a default constructor for efficiency reasons, you
    // must start by adding a multi-block to the group before anything else.
    PLB_ASSERT(!blocks.empty());
    MultiBlock2D *block = generateMultiContainerBlock(*blocks[0], envelopeWidth).release();
    block->setRefinementLevel(gridLevel);
    return add(block, name);
}

MultiContainerBlock2D &Group2D::getContainer(plint id)
{
    if (id >= (plint)blocks.size()) {
        plbLogicError("Group has no block of ID " + util::val2str(id));
    }
    MultiContainerBlock2D *container = dynamic_cast<MultiContainerBlock2D *>(blocks[id]);
    if (!container) {
        plbLogicError("Block with ID " + util::val2str(id) + " is not a container");
    }
    return *container;
}

MultiContainerBlock2D &Group2D::getContainer(std::string name)
{
    std::map<std::string, plint>::const_iterator it = ids.find(name);
    if (it == ids.end()) {
        plbLogicError("Group has no block of name " + name);
    }
    plint id = it->second;
    PLB_ASSERT(id < (plint)blocks.size());
    MultiContainerBlock2D *container = dynamic_cast<MultiContainerBlock2D *>(blocks[id]);
    if (!container) {
        plbLogicError("Block with name \"" + name + "\" is not a container");
    }
    return *container;
}

}  // namespace plb
