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
 * 2D groups -- generic implementation
 */
#ifndef GROUP_2D_HH
#define GROUP_2D_HH

#include "core/globalDefs.h"
#include "core/util.h"
#include "io/parallelIO.h"
#include "multiBlock/group2D.h"
#include "multiBlock/multiBlockGenerator2D.h"
#include "multiBlock/nonLocalTransfer2D.h"

namespace plb {

template <typename T>
plint Group2D::generateScalar(std::string name, plint envelopeWidth, plint gridLevel)
{
    // Although groups offer a default constructor for efficiency reasons, you
    // must start by adding a multi-block to the group before anything else.
    PLB_ASSERT(!blocks.empty());
    MultiBlock2D *block = generateMultiScalarField<T>(*blocks[0], envelopeWidth).release();
    block->setRefinementLevel(gridLevel);
    return add(block, name);
}

template <typename T>
plint Group2D::generateNTensor(std::string name, plint nDim, plint envelopeWidth, plint gridLevel)
{
    // Although groups offer a default constructor for efficiency reasons, you
    // must start by adding a multi-block to the group before anything else.
    PLB_ASSERT(!blocks.empty());
    MultiBlock2D *block = generateMultiNTensorField2D<T>(*blocks[0], envelopeWidth, nDim);
    block->setRefinementLevel(gridLevel);
    return add(block, name);
}

template <typename T, int nDim>
plint Group2D::generateTensor(std::string name, plint envelopeWidth, plint gridLevel)
{
    // Although groups offer a default constructor for efficiency reasons, you
    // must start by adding a multi-block to the group before anything else.
    PLB_ASSERT(!blocks.empty());
    MultiBlock2D *block = generateMultiTensorField<T, nDim>(*blocks[0], envelopeWidth).release();
    block->setRefinementLevel(gridLevel);
    return add(block, name);
}

template <typename T, template <typename U> class Descriptor>
plint Group2D::generateLattice(std::string name, plint envelopeWidth, plint gridLevel)
{
    // Although groups offer a default constructor for efficiency reasons, you
    // must start by adding a multi-block to the group before anything else.
    PLB_ASSERT(!blocks.empty());
    MultiBlock2D *block =
        generateMultiBlockLattice<T, Descriptor>(*blocks[0], envelopeWidth).release();
    block->setRefinementLevel(gridLevel);
    return add(block, name);
}

template <class ParticleFieldT>
plint Group2D::generateParticleField(std::string name, plint envelopeWidth, plint gridLevel)
{
    // Although groups offer a default constructor for efficiency reasons, you
    // must start by adding a multi-block to the group before anything else.
    PLB_ASSERT(!blocks.empty());
    MultiBlock2D *block =
        generateMultiParticleField2D<ParticleFieldT>(*blocks[0], envelopeWidth).release();
    block->setRefinementLevel(gridLevel);
    return add(block, name);
}

template <typename T>
MultiScalarField2D<T> &Group2D::getScalar(plint id)
{
    if (id >= (plint)blocks.size()) {
        plbLogicError("Group has no block of ID " + util::val2str(id));
    }
    MultiScalarField2D<T> *scalar = dynamic_cast<MultiScalarField2D<T> *>(blocks[id]);
    if (!scalar) {
        plbLogicError("Block with ID " + util::val2str(id) + " is not a scalar-field");
    }
    return *scalar;
}

template <typename T>
MultiScalarField2D<T> &Group2D::getScalar(std::string name)
{
    std::map<std::string, plint>::const_iterator it = ids.find(name);
    if (it == ids.end()) {
        plbLogicError("Group has no block of name " + name);
    }
    plint id = it->second;
    PLB_ASSERT(id < (plint)blocks.size());
    MultiScalarField2D<T> *scalar = dynamic_cast<MultiScalarField2D<T> *>(blocks[id]);
    if (!scalar) {
        plbLogicError("Block with name \"" + name + "\" is not a scalar-field");
    }
    return *scalar;
}

template <typename T>
MultiNTensorField2D<T> &Group2D::getNTensor(plint id)
{
    if (id >= (plint)blocks.size()) {
        plbLogicError("Group has no block of ID " + util::val2str(id));
    }
    MultiNTensorField2D<T> *nTensor = dynamic_cast<MultiNTensorField2D<T> *>(blocks[id]);
    if (!nTensor) {
        plbLogicError("Block with ID " + util::val2str(id) + " is not a ntensor-field");
    }
    return *nTensor;
}

template <typename T>
MultiNTensorField2D<T> &Group2D::getNTensor(std::string name)
{
    std::map<std::string, plint>::const_iterator it = ids.find(name);
    if (it == ids.end()) {
        plbLogicError("Group has no block of name " + name);
    }
    plint id = it->second;
    PLB_ASSERT(id < (plint)blocks.size());
    MultiNTensorField2D<T> *nTensor = dynamic_cast<MultiNTensorField2D<T> *>(blocks[id]);
    if (!nTensor) {
        plbLogicError("Block with name \"" + name + "\" is not a ntensor-field");
    }
    return *nTensor;
}

template <typename T, int nDim>
MultiTensorField2D<T, nDim> &Group2D::getTensor(plint id)
{
    if (id >= (plint)blocks.size()) {
        plbLogicError("Group has no block of ID " + util::val2str(id));
    }
    MultiTensorField2D<T, nDim> *tensor = dynamic_cast<MultiTensorField2D<T, nDim> *>(blocks[id]);
    if (!tensor) {
        plbLogicError("Block with ID " + util::val2str(id) + " is not a scalar-field");
    }
    return *tensor;
}

template <typename T, int nDim>
MultiTensorField2D<T, nDim> &Group2D::getTensor(std::string name)
{
    std::map<std::string, plint>::const_iterator it = ids.find(name);
    if (it == ids.end()) {
        plbLogicError("Group has no block of name " + name);
    }
    plint id = it->second;
    PLB_ASSERT(id < (plint)blocks.size());
    MultiTensorField2D<T, nDim> *tensor = dynamic_cast<MultiTensorField2D<T, nDim> *>(blocks[id]);
    if (!tensor) {
        plbLogicError("Block with name \"" + name + "\" is not a tensor-field");
    }
    return *tensor;
}

template <typename T, template <typename U> class Descriptor>
MultiBlockLattice2D<T, Descriptor> &Group2D::getLattice(plint id)
{
    if (id >= (plint)blocks.size()) {
        plbLogicError("Group has no block of ID " + util::val2str(id));
    }
    MultiBlockLattice2D<T, Descriptor> *lattice =
        dynamic_cast<MultiBlockLattice2D<T, Descriptor> *>(blocks[id]);
    if (!lattice) {
        plbLogicError("Block with ID " + util::val2str(id) + " is not a block-lattice");
    }
    return *lattice;
}

template <typename T, template <typename U> class Descriptor>
MultiBlockLattice2D<T, Descriptor> &Group2D::getLattice(std::string name)
{
    std::map<std::string, plint>::const_iterator it = ids.find(name);
    if (it == ids.end()) {
        plbLogicError("Group has no block of name " + name);
    }
    plint id = it->second;
    PLB_ASSERT(id < (plint)blocks.size());
    MultiBlockLattice2D<T, Descriptor> *lattice =
        dynamic_cast<MultiBlockLattice2D<T, Descriptor> *>(blocks[id]);
    if (!lattice) {
        plbLogicError("Block with name \"" + name + "\" is not a block-lattice");
    }
    return *lattice;
}

template <class ParticleFieldT>
MultiParticleField2D<ParticleFieldT> &Group2D::getParticleField(plint id)
{
    if (id >= (plint)blocks.size()) {
        plbLogicError("Group has no block of ID " + util::val2str(id));
    }
    MultiParticleField2D<ParticleFieldT> *field =
        dynamic_cast<MultiParticleField2D<ParticleFieldT> *>(blocks[id]);
    if (!field) {
        plbLogicError("Block with ID " + util::val2str(id) + " is not a particle-field");
    }
    return *field;
}

template <class ParticleFieldT>
MultiParticleField2D<ParticleFieldT> &Group2D::getParticleField(std::string name)
{
    std::map<std::string, plint>::const_iterator it = ids.find(name);
    if (it == ids.end()) {
        plbLogicError("Group has no block of name " + name);
    }
    plint id = it->second;
    PLB_ASSERT(id < (plint)blocks.size());
    MultiParticleField2D<ParticleFieldT> *field =
        dynamic_cast<MultiParticleField2D<ParticleFieldT> *>(blocks[id]);
    if (!field) {
        plbLogicError("Block with name \"" + name + "\" is not a particle-field");
    }
    return *field;
}

}  // namespace plb

#endif  // GROUP_2D_HH
