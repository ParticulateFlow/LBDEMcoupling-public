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
 * 3D groups -- generic implementation
 */
#ifndef GROUP_3D_HH
#define GROUP_3D_HH

#include "core/globalDefs.h"
#include "core/runTimeDiagnostics.h"
#include "core/util.h"
#include "dataProcessors/dataAnalysisWrapper3D.h"
#include "io/parallelIO.h"
#include "multiBlock/group3D.h"
#include "multiBlock/multiBlockGenerator3D.h"
#include "multiBlock/nonLocalTransfer3D.h"

namespace plb {

template <typename T>
plint Group3D::generateScalar(std::string name, plint envelopeWidth, plint gridLevel)
{
    // Although groups offer a default constructor for efficiency reasons, you
    // must start by adding a multi-block to the group before anything else.
    PLB_ASSERT(!blocks.empty());
    MultiBlock3D *block = generateMultiScalarField<T>(*blocks[0], envelopeWidth).release();
    block->setRefinementLevel(gridLevel);
    return addNoCheck(block, name);
}

template <typename T>
plint Group3D::generateNTensor(std::string name, plint nDim, plint envelopeWidth, plint gridLevel)
{
    // Although groups offer a default constructor for efficiency reasons, you
    // must start by adding a multi-block to the group before anything else.
    PLB_ASSERT(!blocks.empty());
    MultiBlock3D *block = generateMultiNTensorField3D<T>(*blocks[0], envelopeWidth, nDim);
    block->setRefinementLevel(gridLevel);
    return addNoCheck(block, name);
}

template <typename T>
plint Group3D::generateNTensor(
    std::string name, plint nDim, std::vector<Box3D> const &bulks, plint envelopeWidth,
    plint gridLevel, bool cropBoundingBox)
{
    // Although groups offer a default constructor for efficiency reasons, you
    // must start by adding a multi-block to the group before anything else.
    PLB_ASSERT(!blocks.empty());
    MultiBlockManagement3D newManagement = align(
        bulks, blocks[0]->getMultiBlockManagement(), envelopeWidth, gridLevel, cropBoundingBox);
    MultiNTensorField3D<T> *field = new MultiNTensorField3D<T>(
        nDim, newManagement, defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>());
    field->setRefinementLevel(gridLevel);
    return addNoCheck(field, name);
}

template <typename T, int nDim>
plint Group3D::generateTensor(std::string name, plint envelopeWidth, plint gridLevel)
{
    // Although groups offer a default constructor for efficiency reasons, you
    // must start by adding a multi-block to the group before anything else.
    PLB_ASSERT(!blocks.empty());
    MultiBlock3D *block = generateMultiTensorField<T, nDim>(*blocks[0], envelopeWidth).release();
    block->setRefinementLevel(gridLevel);
    return addNoCheck(block, name);
}

template <typename T, template <typename U> class Descriptor>
plint Group3D::generateLattice(std::string name, plint envelopeWidth, plint gridLevel)
{
    // Although groups offer a default constructor for efficiency reasons, you
    // must start by adding a multi-block to the group before anything else.
    PLB_ASSERT(!blocks.empty());
    MultiBlock3D *block =
        generateMultiBlockLattice<T, Descriptor>(*blocks[0], envelopeWidth).release();
    block->setRefinementLevel(gridLevel);
    return addNoCheck(block, name);
}

template <typename T, template <typename U> class Descriptor>
plint Group3D::generateDenseParticles(std::string name, plint envelopeWidth, plint gridLevel)
{
    // Although groups offer a default constructor for efficiency reasons, you
    // must start by adding a multi-block to the group before anything else.
    PLB_ASSERT(!blocks.empty());
    MultiBlock3D *block =
        generateMultiParticleField3D<T, Descriptor, DenseParticleField3D<T, Descriptor> >(
            *blocks[0], envelopeWidth)
            .release();
    block->setRefinementLevel(gridLevel);
    return addNoCheck(block, name);
}

template <typename T>
MultiScalarField3D<T> &Group3D::getScalar(plint id)
{
    if (id >= (plint)blocks.size()) {
        plbLogicError("Group has no block of ID " + util::val2str(id));
    }
    MultiScalarField3D<T> *scalar = dynamic_cast<MultiScalarField3D<T> *>(blocks[id]);
    if (!scalar) {
        plbLogicError("Block with ID " + util::val2str(id) + " is not a scalar-field");
    }
    return *scalar;
}

template <typename T>
MultiScalarField3D<T> &Group3D::getScalar(std::string name)
{
    std::map<std::string, plint>::const_iterator it = ids.find(name);
    if (it == ids.end()) {
        plbLogicError("Group has no block of name " + name);
    }
    plint id = it->second;
    PLB_ASSERT(id < (plint)blocks.size());
    MultiScalarField3D<T> *scalar = dynamic_cast<MultiScalarField3D<T> *>(blocks[id]);
    if (!scalar) {
        plbLogicError("Block with name \"" + name + "\" is not a scalar-field");
    }
    return *scalar;
}

template <typename T>
MultiNTensorField3D<T> &Group3D::getNTensor(plint id)
{
    if (id >= (plint)blocks.size()) {
        plbLogicError("Group has no block of ID " + util::val2str(id));
    }
    MultiNTensorField3D<T> *nTensor = dynamic_cast<MultiNTensorField3D<T> *>(blocks[id]);
    if (!nTensor) {
        plbLogicError("Block with ID " + util::val2str(id) + " is not a ntensor-field");
    }
    return *nTensor;
}

template <typename T>
MultiNTensorField3D<T> &Group3D::getNTensor(std::string name)
{
    std::map<std::string, plint>::const_iterator it = ids.find(name);
    if (it == ids.end()) {
        plbLogicError("Group has no block of name " + name);
    }
    plint id = it->second;
    PLB_ASSERT(id < (plint)blocks.size());
    MultiNTensorField3D<T> *nTensor = dynamic_cast<MultiNTensorField3D<T> *>(blocks[id]);
    if (!nTensor) {
        plbLogicError("Block with name \"" + name + "\" is not a ntensor-field");
    }
    return *nTensor;
}

template <typename T, int nDim>
MultiTensorField3D<T, nDim> &Group3D::getTensor(plint id)
{
    if (id >= (plint)blocks.size()) {
        plbLogicError("Group has no block of ID " + util::val2str(id));
    }
    MultiTensorField3D<T, nDim> *tensor = dynamic_cast<MultiTensorField3D<T, nDim> *>(blocks[id]);
    if (!tensor) {
        plbLogicError("Block with ID " + util::val2str(id) + " is not a scalar-field");
    }
    return *tensor;
}

template <typename T, int nDim>
MultiTensorField3D<T, nDim> &Group3D::getTensor(std::string name)
{
    std::map<std::string, plint>::const_iterator it = ids.find(name);
    if (it == ids.end()) {
        plbLogicError("Group has no block of name " + name);
    }
    plint id = it->second;
    PLB_ASSERT(id < (plint)blocks.size());
    MultiTensorField3D<T, nDim> *tensor = dynamic_cast<MultiTensorField3D<T, nDim> *>(blocks[id]);
    if (!tensor) {
        plbLogicError("Block with name \"" + name + "\" is not a tensor-field");
    }
    return *tensor;
}

template <typename T, template <typename U> class Descriptor>
MultiBlockLattice3D<T, Descriptor> &Group3D::getLattice(plint id)
{
    if (id >= (plint)blocks.size()) {
        plbLogicError("Group has no block of ID " + util::val2str(id));
    }
    MultiBlockLattice3D<T, Descriptor> *lattice =
        dynamic_cast<MultiBlockLattice3D<T, Descriptor> *>(blocks[id]);
    if (!lattice) {
        plbLogicError("Block with ID " + util::val2str(id) + " is not a block-lattice");
    }
    return *lattice;
}

template <typename T, template <typename U> class Descriptor>
MultiBlockLattice3D<T, Descriptor> &Group3D::getLattice(std::string name)
{
    std::map<std::string, plint>::const_iterator it = ids.find(name);
    if (it == ids.end()) {
        plbLogicError("Group has no block of name " + name);
    }
    plint id = it->second;
    PLB_ASSERT(id < (plint)blocks.size());
    MultiBlockLattice3D<T, Descriptor> *lattice =
        dynamic_cast<MultiBlockLattice3D<T, Descriptor> *>(blocks[id]);
    if (!lattice) {
        plbLogicError("Block with name \"" + name + "\" is not a block-lattice");
    }
    return *lattice;
}

template <typename T, template <typename U> class Descriptor>
MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &Group3D::getDenseParticles(plint id)
{
    if (id >= (plint)blocks.size()) {
        plbLogicError("Group has no block of ID " + util::val2str(id));
    }
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > *particles =
        dynamic_cast<MultiParticleField3D<DenseParticleField3D<T, Descriptor> > *>(blocks[id]);
    if (!particles) {
        plbLogicError("Block with ID " + util::val2str(id) + " is not a dense particle-field");
    }
    return *particles;
}

template <typename T, template <typename U> class Descriptor>
MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &Group3D::getDenseParticles(
    std::string name)
{
    std::map<std::string, plint>::const_iterator it = ids.find(name);
    if (it == ids.end()) {
        plbLogicError("Group has no block of name " + name);
    }
    plint id = it->second;
    PLB_ASSERT(id < (plint)blocks.size());
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > *particles =
        dynamic_cast<MultiParticleField3D<DenseParticleField3D<T, Descriptor> > *>(blocks[id]);
    if (!particles) {
        plbLogicError("Block with name \"" + name + "\" is not a dense particle-field");
    }
    return *particles;
}

template <typename T, typename TConv>
void addTransform(
    Group3D &group, std::unique_ptr<MultiScalarField3D<T> > field, std::string name,
    TConv scalingFactor, TConv additiveOffset)
{
    std::unique_ptr<MultiScalarField3D<TConv> > transformedField = copyConvert<T, TConv>(*field);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    if (!util::isZero(additiveOffset)) {
        addInPlace(*transformedField, additiveOffset);
    }
    group.add(transformedField.release(), name);
}

template <typename T>
void addTransform(
    Group3D &group, std::unique_ptr<MultiScalarField3D<T> > field, std::string name,
    T scalingFactor, T additiveOffset)
{
    group.add(add(*multiply(*field, scalingFactor), additiveOffset).release(), name);
}

template <typename T, typename TConv, int nDim>
void addTransform(
    Group3D &group, std::unique_ptr<MultiTensorField3D<T, nDim> > field, std::string name,
    TConv scalingFactor)
{
    std::unique_ptr<MultiTensorField3D<TConv, nDim> > transformedField =
        copyConvert<T, TConv>(*field);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    group.add(transformedField.release(), name);
}

template <typename T, int nDim>
void addTransform(
    Group3D &group, std::unique_ptr<MultiTensorField3D<T, nDim> > field, std::string name,
    T scalingFactor)
{
    group.add(multiply(*field, scalingFactor).release(), name);
}

template <typename T, typename TConv>
void addTransform(
    Group3D &group, MultiScalarField3D<T> &field, std::string name, TConv scalingFactor,
    TConv additiveOffset)
{
    std::unique_ptr<MultiScalarField3D<TConv> > transformedField = copyConvert<T, TConv>(field);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    if (!util::isZero(additiveOffset)) {
        addInPlace(*transformedField, additiveOffset);
    }
    group.add(transformedField.release(), name);
}

template <typename T>
void addTransform(
    Group3D &group, MultiScalarField3D<T> &field, std::string name, T scalingFactor,
    T additiveOffset)
{
    group.add(add(*multiply(field, scalingFactor), additiveOffset).release(), name);
}

template <typename T, typename TConv, int nDim>
void addTransform(
    Group3D &group, MultiTensorField3D<T, nDim> &field, std::string name, TConv scalingFactor)
{
    std::unique_ptr<MultiTensorField3D<TConv, nDim> > transformedField =
        copyConvert<T, TConv>(field);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    group.add(transformedField.release(), name);
}

template <typename T, int nDim>
void addTransform(
    Group3D &group, MultiTensorField3D<T, nDim> &field, std::string name, T scalingFactor)
{
    group.add(multiply(field, scalingFactor).release(), name);
}

}  // namespace plb

#endif  // GROUP_3D_HH
