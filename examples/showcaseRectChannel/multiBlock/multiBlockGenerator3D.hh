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
 * Copy 3D multiblocks on a new parallel distribution -- generic implementation.
 */

/*
 * All functions that take as an argument a Box3D, or compute intersections, or join,
 * or crop, do not adjust periodicity of the newly created blocks by default. The user
 * is responsible to take care of this matter explicitly.
 */

#ifndef MULTI_BLOCK_GENERATOR_3D_HH
#define MULTI_BLOCK_GENERATOR_3D_HH

#include "core/globalDefs.h"
#include "dataProcessors/dataAnalysisWrapper3D.h"
#include "dataProcessors/ntensorAnalysisWrapper3D.h"
#include "multiBlock/defaultMultiBlockPolicy3D.h"
#include "multiBlock/localMultiBlockInfo3D.h"
#include "multiBlock/multiBlockOperations3D.h"
#include "multiBlock/nonLocalTransfer3D.h"
#include "multiBlock/sparseBlockStructure3D.h"

namespace plb {

inline void transferDataProcessors(MultiBlock3D const &from, MultiBlock3D &to)
{
    std::vector<MultiBlock3D::ProcessorStorage3D> newProcessors(from.getStoredProcessors());
    // Redirect all references-to-self to the new self.
    id_t oldId = from.getId();
    id_t newId = to.getId();
    for (pluint iProcessor = 0; iProcessor < newProcessors.size(); ++iProcessor) {
        newProcessors[iProcessor].replace(oldId, newId);
        DataProcessorGenerator3D *newGenerator = newProcessors[iProcessor].getGenerator().clone();
        if (newGenerator->extract(to.getBoundingBox())) {
            addInternalProcessor(
                *newGenerator, newProcessors[iProcessor].getMultiBlocks(),
                newProcessors[iProcessor].getLevel());
        }
        delete newGenerator;
    }
}

/* *************** 1. MultiScalarField ************************************** */

// TODO: Does it make sense to have a useless domain here?
template <typename T>
void transferScalarFieldLocal(
    MultiScalarField3D<T> &from, MultiScalarField3D<T> &to, Box3D const &domain)
{
    // 1. Copy all data from the old to the new field.
    plb::copy(from, to, from.getBoundingBox());
    // 2. Reconstruct the data processors.
    transferDataProcessors(from, to);
}

template <typename T>
void transferScalarFieldNonLocal(
    MultiScalarField3D<T> const &from, MultiScalarField3D<T> &to, Box3D const &domain)
{
    // 1. Copy all data from the old to the new field. This includes dynamics
    //    objects which must be fully serialized and regenerated.
    copyNonLocal(from, to, domain);
    // 2. Reconstruct the data processors.
    transferDataProcessors(from, to);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > generateMultiScalarField(
    Box3D boundingBox, plint envelopeWidth)
{
    return std::unique_ptr<MultiScalarField3D<T> >(new MultiScalarField3D<T>(
        defaultMultiBlockPolicy3D().getMultiBlockManagement(boundingBox, envelopeWidth),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiScalarAccess<T>()));
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > generateMultiScalarField(
    Box3D boundingBox, T iniVal, plint envelopeWidth)
{
    return std::unique_ptr<MultiScalarField3D<T> >(new MultiScalarField3D<T>(
        defaultMultiBlockPolicy3D().getMultiBlockManagement(boundingBox, envelopeWidth),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiScalarAccess<T>(), iniVal));
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > generateMultiScalarField(
    MultiBlock3D &multiBlock, plint envelopeWidth)
{
    MultiBlockManagement3D sparseBlockManagement(multiBlock.getMultiBlockManagement());
    MultiScalarField3D<T> *field = new MultiScalarField3D<T>(
        MultiBlockManagement3D(
            sparseBlockManagement.getSparseBlockStructure(),
            sparseBlockManagement.getThreadAttribution().clone(), envelopeWidth,
            sparseBlockManagement.getRefinementLevel()),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiScalarAccess<T>());

    field->periodicity().toggle(0, multiBlock.periodicity().get(0));
    field->periodicity().toggle(1, multiBlock.periodicity().get(1));
    field->periodicity().toggle(2, multiBlock.periodicity().get(2));

    return std::unique_ptr<MultiScalarField3D<T> >(field);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > defaultGenerateMultiScalarField3D(
    MultiBlockManagement3D const &management, T iniVal)
{
    return std::unique_ptr<MultiScalarField3D<T> >(new MultiScalarField3D<T>(
        management, defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiScalarAccess<T>(), iniVal));
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > clone(
    MultiScalarField3D<T> &originalField, Box3D const &subDomain, bool crop)
{
    std::unique_ptr<MultiScalarField3D<T> > clonedField(
        generateMultiScalarField<T>(originalField, subDomain, crop));

    transferScalarFieldLocal(originalField, *clonedField, originalField.getBoundingBox());

    return clonedField;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > generateMultiScalarField(
    MultiBlock3D const &originalField, Box3D const &intersection, bool crop)
{
    return std::unique_ptr<MultiScalarField3D<T> >(new MultiScalarField3D<T>(
        intersect(originalField.getMultiBlockManagement(), intersection, crop),
        originalField.getBlockCommunicator().clone(), originalField.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiScalarAccess<T>()));
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > generateIntersectMultiScalarField(
    MultiBlock3D const &originalField1, MultiBlock3D const &originalField2, bool crop)
{
    return std::unique_ptr<MultiScalarField3D<T> >(new MultiScalarField3D<T>(
        intersect(
            originalField1.getMultiBlockManagement(), originalField2.getMultiBlockManagement(),
            crop),
        originalField1.getBlockCommunicator().clone(),
        originalField1.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiScalarAccess<T>()));
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > generateIntersectMultiScalarField(
    MultiBlock3D const &originalField1, MultiBlock3D const &originalField2,
    Box3D const &intersection, bool crop)
{
    MultiBlockManagement3D intersectedBlocks(intersect(
        originalField1.getMultiBlockManagement(), originalField2.getMultiBlockManagement(), crop));
    MultiBlockManagement3D intersectWithDomain(intersect(intersectedBlocks, intersection, crop));
    return std::unique_ptr<MultiScalarField3D<T> >(new MultiScalarField3D<T>(
        intersectWithDomain, originalField1.getBlockCommunicator().clone(),
        originalField1.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiScalarAccess<T>()));
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > generateJoinMultiScalarField(
    MultiBlock3D const &originalField1, MultiBlock3D const &originalField2)
{
    return std::unique_ptr<MultiScalarField3D<T> >(new MultiScalarField3D<T>(
        block_union(
            originalField1.getMultiBlockManagement(), originalField2.getMultiBlockManagement()),
        originalField1.getBlockCommunicator().clone(),
        originalField1.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiScalarAccess<T>()));
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > extendEnvelopeWidth(
    MultiScalarField3D<T> &originalBlock, plint envelopeWidth)
{
    std::unique_ptr<MultiScalarField3D<T> > extendedBlock =
        generateMultiScalarField<T>(originalBlock, envelopeWidth);
    plb::copy(originalBlock, *extendedBlock, originalBlock.getBoundingBox());
    return extendedBlock;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > extend(
    MultiScalarField3D<T> &originalBlock, Box3D const &addedBlock)
{
    std::unique_ptr<MultiScalarField3D<T> > newBlock(new MultiScalarField3D<T>(
        extend(originalBlock.getMultiBlockManagement(), addedBlock, addedBlock),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiScalarAccess<T>()));

    transferScalarFieldLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > except(
    MultiScalarField3D<T> &originalBlock, Box3D const &exceptedBlock)
{
    std::unique_ptr<MultiScalarField3D<T> > newBlock(new MultiScalarField3D<T>(
        except(originalBlock.getMultiBlockManagement(), exceptedBlock),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiScalarAccess<T>()));

    transferScalarFieldLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > redistribute(
    MultiScalarField3D<T> const &originalField, SparseBlockStructure3D const &newBlockStructure,
    bool adjustPeriodicity)
{
    std::unique_ptr<MultiScalarField3D<T> > newField(new MultiScalarField3D<T>(
        MultiBlockManagement3D(
            newBlockStructure,
            originalField.getMultiBlockManagement().getThreadAttribution().clone(),
            originalField.getMultiBlockManagement().getEnvelopeWidth(),
            originalField.getMultiBlockManagement().getRefinementLevel()),
        originalField.getBlockCommunicator().clone(), originalField.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiScalarAccess<T>()));

    if (adjustPeriodicity) {
        newField->periodicity().toggle(0, originalField.periodicity().get(0));
        newField->periodicity().toggle(1, originalField.periodicity().get(1));
        newField->periodicity().toggle(2, originalField.periodicity().get(2));
    }

    transferScalarFieldNonLocal(originalField, *newField, originalField.getBoundingBox());

    return newField;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > redistribute(
    MultiScalarField3D<T> const &originalField, SparseBlockStructure3D const &newBlockStructure,
    Box3D const &intersection, bool crop)
{
    bool adjustPeriodicity = false;
    return redistribute(
        originalField, intersect(newBlockStructure, intersection, crop), adjustPeriodicity);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > align(
    MultiScalarField3D<T> const &originalBlock, MultiBlock3D const &partnerBlock)
{
    std::unique_ptr<MultiScalarField3D<T> > newBlock(new MultiScalarField3D<T>(
        align(originalBlock.getMultiBlockManagement(), partnerBlock.getMultiBlockManagement()),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiScalarAccess<T>()));

    newBlock->periodicity().toggle(0, originalBlock.periodicity().get(0));
    newBlock->periodicity().toggle(1, originalBlock.periodicity().get(1));
    newBlock->periodicity().toggle(2, originalBlock.periodicity().get(2));

    transferScalarFieldNonLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > reparallelize(MultiScalarField3D<T> const &originalBlock)
{
    return reparallelize(originalBlock, 16, 16, 16);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > reparallelize(
    MultiScalarField3D<T> const &originalBlock, plint blockLx, plint blockLy, plint blockLz)
{
    std::unique_ptr<MultiScalarField3D<T> > newBlock(new MultiScalarField3D<T>(
        reparallelize(originalBlock.getMultiBlockManagement(), blockLx, blockLy, blockLz),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiScalarAccess<T>()));

    newBlock->periodicity().toggle(0, originalBlock.periodicity().get(0));
    newBlock->periodicity().toggle(1, originalBlock.periodicity().get(1));
    newBlock->periodicity().toggle(2, originalBlock.periodicity().get(2));

    transferScalarFieldNonLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

/* *************** 2. MultiNTensorField ************************************** */

template <typename T>
std::unique_ptr<MultiNTensorField3D<T> > defaultGenerateMultiNTensorField3D(
    MultiBlockManagement3D const &management, plint nDim)
{
    MultiNTensorField3D<T> *field = new MultiNTensorField3D<T>(
        nDim, management, defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>());
    return std::unique_ptr<MultiNTensorField3D<T> >(field);
}

template <typename T>
MultiNTensorField3D<T> *generateMultiNTensorField3D(
    MultiBlock3D &multiBlock, plint envelopeWidth, plint ndim)
{
    MultiBlockManagement3D sparseBlockManagement(multiBlock.getMultiBlockManagement());
    MultiNTensorField3D<T> *field = new MultiNTensorField3D<T>(
        ndim,
        MultiBlockManagement3D(
            sparseBlockManagement.getSparseBlockStructure(),
            sparseBlockManagement.getThreadAttribution().clone(), envelopeWidth,
            sparseBlockManagement.getRefinementLevel()),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>());

    field->periodicity().toggle(0, multiBlock.periodicity().get(0));
    field->periodicity().toggle(1, multiBlock.periodicity().get(1));
    field->periodicity().toggle(2, multiBlock.periodicity().get(2));

    return field;
}

template <typename T>
MultiNTensorField3D<T> *generateMultiNTensorField3D(Box3D const &domain, plint ndim)
{
    plint defaultEnvelopeWidth = 1;
    MultiNTensorField3D<T> *field = new MultiNTensorField3D<T>(
        ndim, defaultMultiBlockPolicy3D().getMultiBlockManagement(domain, defaultEnvelopeWidth),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>());
    return field;
}

template <typename T>
MultiNTensorField3D<T> *generateMultiNTensorField3D(
    Box3D const &domain, plint ndim, T *iniVal, plint envelopeWidth)
{
    MultiNTensorField3D<T> *field = new MultiNTensorField3D<T>(
        ndim, iniVal, defaultMultiBlockPolicy3D().getMultiBlockManagement(domain, envelopeWidth),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>());
    return field;
}

template <typename T>
void transferNTensorFieldLocal(
    MultiNTensorField3D<T> &from, MultiNTensorField3D<T> &to, [[maybe_unused]] Box3D const &domain)
// TODO: Remove the domain argument. Not used and misleading.
{
    // 1. Copy all data from the old to the new field.
    plb::copy(from, to, from.getBoundingBox());
    // 2. Reconstruct the data processors.
    transferDataProcessors(from, to);
}

template <typename T>
void transferNTensorFieldNonLocal(
    MultiNTensorField3D<T> const &from, MultiNTensorField3D<T> &to, Box3D const &domain)
{
    // 1. Copy all data from the old to the new field. This includes dynamics
    //    objects which must be fully serialized and regenerated.
    copyNonLocal(from, to, domain);
    // 2. Reconstruct the data processors.
    transferDataProcessors(from, to);
}

template <typename T>
MultiNTensorField3D<T> *clone(
    MultiNTensorField3D<T> &originalField, Box3D const &subDomain, bool crop)
{
    MultiNTensorField3D<T> *clonedField =
        generateMultiNTensorField<T>(originalField, subDomain, originalField.getNdim(), crop);

    transferNTensorFieldLocal(originalField, *clonedField, originalField.getBoundingBox());

    return clonedField;
}

template <typename T>
MultiNTensorField3D<T> *generateMultiNTensorField(
    MultiBlock3D const &originalField, Box3D const &intersection, plint nDim, bool crop)
{
    MultiNTensorField3D<T> *newField = new MultiNTensorField3D<T>(
        nDim, intersect(originalField.getMultiBlockManagement(), intersection, crop),
        originalField.getBlockCommunicator().clone(), originalField.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>());
    return newField;
}

template <typename T1, typename T2>
MultiNTensorField3D<T2> *generateNTensorFieldFromNTensor3D(
    MultiNTensorField3D<T1> const &field, Box3D const &intersection, plint nDim)
{
    MultiNTensorField3D<T2> *newField = generateMultiNTensorField<T2>(field, intersection, nDim);
    return newField;
}

template <typename T1, typename T2, template <typename U> class Descriptor>
MultiNTensorField3D<T1> *generateNTensorFieldFromBlockLattice3D(
    MultiBlockLattice3D<T2, Descriptor> const &lattice, Box3D const &intersection, plint nDim)
{
    MultiNTensorField3D<T1> *newField = generateMultiNTensorField<T1>(lattice, intersection, nDim);
    return newField;
}

template <typename T>
MultiNTensorField3D<T> *generateIntersectMultiNTensorField(
    MultiBlock3D const &originalField1, MultiBlock3D const &originalField2, plint nDim, bool crop)
{
    MultiNTensorField3D<T> *newField = new MultiNTensorField3D<T>(
        nDim,
        intersect(
            originalField1.getMultiBlockManagement(), originalField2.getMultiBlockManagement(),
            crop),
        originalField1.getBlockCommunicator().clone(),
        originalField1.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>());
    return newField;
}

template <typename T>
MultiNTensorField3D<T> *generateIntersectMultiNTensorField(
    MultiBlock3D const &originalField1, MultiBlock3D const &originalField2,
    Box3D const &intersection, plint nDim, bool crop)
{
    MultiBlockManagement3D intersectedBlocks(intersect(
        originalField1.getMultiBlockManagement(), originalField2.getMultiBlockManagement(), crop));
    MultiBlockManagement3D intersectWithDomain(intersect(intersectedBlocks, intersection, crop));
    MultiNTensorField3D<T> *newField = new MultiNTensorField3D<T>(
        nDim, intersectWithDomain, originalField1.getBlockCommunicator().clone(),
        originalField1.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>());
    return newField;
}

template <typename T>
MultiNTensorField3D<T> *generateJoinMultiNTensorField(
    MultiBlock3D const &originalField1, MultiBlock3D const &originalField2, plint nDim)
{
    MultiNTensorField3D<T> *newField = new MultiNTensorField3D<T>(
        nDim,
        block_union(
            originalField1.getMultiBlockManagement(), originalField2.getMultiBlockManagement()),
        originalField1.getBlockCommunicator().clone(),
        originalField1.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>());
    return newField;
}

template <typename T>
std::unique_ptr<MultiNTensorField3D<T> > extendEnvelopeWidth(
    MultiNTensorField3D<T> &originalBlock, plint envelopeWidth)
{
    std::unique_ptr<MultiNTensorField3D<T> > extendedBlock(
        generateMultiNTensorField3D<T>(originalBlock, envelopeWidth, originalBlock.getNdim()));
    plb::copy(originalBlock, *extendedBlock, originalBlock.getBoundingBox());
    return extendedBlock;
}

template <typename T>
MultiNTensorField3D<T> *extend(MultiNTensorField3D<T> &originalBlock, Box3D const &addedBlock)
{
    MultiNTensorField3D<T> *newBlock = new MultiNTensorField3D<T>(
        originalBlock.getNdim(),
        extend(originalBlock.getMultiBlockManagement(), addedBlock, addedBlock),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>());

    transferNTensorFieldLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T>
MultiNTensorField3D<T> *except(MultiNTensorField3D<T> &originalBlock, Box3D const &exceptedBlock)
{
    MultiNTensorField3D<T> *newBlock = new MultiNTensorField3D<T>(
        originalBlock.getNdim(), except(originalBlock.getMultiBlockManagement(), exceptedBlock),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>());

    transferNTensorFieldLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T>
MultiNTensorField3D<T> *align(
    MultiNTensorField3D<T> const &originalBlock, MultiBlock3D const &partnerBlock)
{
    MultiNTensorField3D<T> *newBlock = new MultiNTensorField3D<T>(
        originalBlock.getNdim(),
        align(originalBlock.getMultiBlockManagement(), partnerBlock.getMultiBlockManagement()),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>());

    newBlock->periodicity().toggle(0, originalBlock.periodicity().get(0));
    newBlock->periodicity().toggle(1, originalBlock.periodicity().get(1));
    newBlock->periodicity().toggle(2, originalBlock.periodicity().get(2));

    transferNTensorFieldNonLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}
template <typename T>
MultiNTensorField3D<T> *reparallelize(MultiNTensorField3D<T> const &originalBlock)
{
    return reparallelize(originalBlock, 16, 16, 16);
}

template <typename T>
MultiNTensorField3D<T> *reparallelize(
    MultiNTensorField3D<T> const &originalBlock, plint blockLx, plint blockLy, plint blockLz)
{
    MultiNTensorField3D<T> *newBlock = new MultiNTensorField3D<T>(
        originalBlock.getNdim(),
        reparallelize(originalBlock.getMultiBlockManagement(), blockLx, blockLy, blockLz),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>());

    newBlock->periodicity().toggle(0, originalBlock.periodicity().get(0));
    newBlock->periodicity().toggle(1, originalBlock.periodicity().get(1));
    newBlock->periodicity().toggle(2, originalBlock.periodicity().get(2));

    transferNTensorFieldNonLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

/* *************** 3. MultiTensorField ************************************** */

template <typename T, int nDim>
void transferTensorFieldLocal(
    MultiTensorField3D<T, nDim> &from, MultiTensorField3D<T, nDim> &to,
    [[maybe_unused]] Box3D const &domain)
// TODO: Remove thedomain argument. Not used and misleading.
{
    // 1. Copy all data from the old to the new field.
    plb::copy(from, to, from.getBoundingBox());
    // 2. Reconstruct the data processors.
    transferDataProcessors(from, to);
}

template <typename T, int nDim>
void transferTensorFieldNonLocal(
    MultiTensorField3D<T, nDim> const &from, MultiTensorField3D<T, nDim> &to, Box3D const &domain)
{
    // 1. Copy all data from the old to the new field. This includes dynamics
    //    objects which must be fully serialized and regenerated.
    copyNonLocal(from, to, domain);
    // 2. Reconstruct the data processors.
    transferDataProcessors(from, to);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > generateMultiTensorField(
    Box3D boundingBox, plint envelopeWidth)
{
    return std::unique_ptr<MultiTensorField3D<T, nDim> >(new MultiTensorField3D<T, nDim>(
        defaultMultiBlockPolicy3D().getMultiBlockManagement(boundingBox, envelopeWidth),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>()));
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > generateMultiTensorField(
    Box3D boundingBox, Array<T, nDim> const &iniVal, plint envelopeWidth)
{
    return std::unique_ptr<MultiTensorField3D<T, nDim> >(new MultiTensorField3D<T, nDim>(
        defaultMultiBlockPolicy3D().getMultiBlockManagement(boundingBox, envelopeWidth),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>(), iniVal));
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > generateMultiTensorField(
    MultiBlock3D &multiBlock, plint envelopeWidth)
{
    MultiBlockManagement3D sparseBlockManagement(multiBlock.getMultiBlockManagement());
    MultiTensorField3D<T, nDim> *field = new MultiTensorField3D<T, nDim>(
        MultiBlockManagement3D(
            sparseBlockManagement.getSparseBlockStructure(),
            sparseBlockManagement.getThreadAttribution().clone(), envelopeWidth,
            sparseBlockManagement.getRefinementLevel()),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>());

    field->periodicity().toggle(0, multiBlock.periodicity().get(0));
    field->periodicity().toggle(1, multiBlock.periodicity().get(1));
    field->periodicity().toggle(2, multiBlock.periodicity().get(2));

    return std::unique_ptr<MultiTensorField3D<T, nDim> >(field);
}

// QUESTION: Why do we have this unnamedDummyArg ?
template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > defaultGenerateMultiTensorField3D(
    MultiBlockManagement3D const &management, [[maybe_unused]] plint unnamedDummyArg)
{
    Array<T, nDim> iniVal;
    iniVal.resetToZero();
    return std::unique_ptr<MultiTensorField3D<T, nDim> >(new MultiTensorField3D<T, nDim>(
        management, defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>(), iniVal));
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > clone(
    MultiTensorField3D<T, nDim> &originalField, Box3D const &subDomain, bool crop)
{
    std::unique_ptr<MultiTensorField3D<T, nDim> > clonedField(
        generateMultiTensorField<T, nDim>(originalField, subDomain, crop));

    transferTensorFieldLocal(originalField, *clonedField, originalField.getBoundingBox());

    return clonedField;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > generateMultiTensorField(
    MultiBlock3D const &originalField, Box3D const &intersection, bool crop)
{
    return std::unique_ptr<MultiTensorField3D<T, nDim> >(new MultiTensorField3D<T, nDim>(
        intersect(originalField.getMultiBlockManagement(), intersection, crop),
        originalField.getBlockCommunicator().clone(), originalField.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>()));
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > generateIntersectMultiTensorField(
    MultiBlock3D const &originalField1, MultiBlock3D const &originalField2, bool crop)
{
    return std::unique_ptr<MultiTensorField3D<T, nDim> >(new MultiTensorField3D<T, nDim>(
        intersect(
            originalField1.getMultiBlockManagement(), originalField2.getMultiBlockManagement(),
            crop),
        originalField1.getBlockCommunicator().clone(),
        originalField1.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>()));
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > generateIntersectMultiTensorField(
    MultiBlock3D const &originalField1, MultiBlock3D const &originalField2,
    Box3D const &intersection, bool crop)
{
    MultiBlockManagement3D intersectedBlocks(intersect(
        originalField1.getMultiBlockManagement(), originalField2.getMultiBlockManagement(), crop));
    MultiBlockManagement3D intersectWithDomain(intersect(intersectedBlocks, intersection, crop));
    return std::unique_ptr<MultiTensorField3D<T, nDim> >(new MultiTensorField3D<T, nDim>(
        intersectWithDomain, originalField1.getBlockCommunicator().clone(),
        originalField1.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>()));
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > generateJoinMultiTensorField(
    MultiBlock3D const &originalField1, MultiBlock3D const &originalField2)
{
    return std::unique_ptr<MultiTensorField3D<T, nDim> >(new MultiTensorField3D<T, nDim>(
        block_union(
            originalField1.getMultiBlockManagement(), originalField2.getMultiBlockManagement()),
        originalField1.getBlockCommunicator().clone(),
        originalField1.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>()));
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > extendEnvelopeWidth(
    MultiTensorField3D<T, nDim> &originalBlock, plint envelopeWidth)
{
    std::unique_ptr<MultiTensorField3D<T, nDim> > extendedBlock =
        generateMultiTensorField<T, nDim>(originalBlock, envelopeWidth);
    plb::copy(originalBlock, *extendedBlock, originalBlock.getBoundingBox());
    return extendedBlock;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > extend(
    MultiTensorField3D<T, nDim> &originalBlock, Box3D const &addedBlock)
{
    std::unique_ptr<MultiTensorField3D<T, nDim> > newBlock(new MultiTensorField3D<T, nDim>(
        extend(originalBlock.getMultiBlockManagement(), addedBlock, addedBlock),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>()));

    transferTensorFieldLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > except(
    MultiTensorField3D<T, nDim> &originalBlock, Box3D const &exceptedBlock)
{
    std::unique_ptr<MultiTensorField3D<T, nDim> > newBlock(new MultiTensorField3D<T, nDim>(
        except(originalBlock.getMultiBlockManagement(), exceptedBlock),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>()));

    transferTensorFieldLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > redistribute(
    MultiTensorField3D<T, nDim> const &originalField,
    SparseBlockStructure3D const &newBlockStructure, bool adjustPeriodicity)
{
    std::unique_ptr<MultiTensorField3D<T, nDim> > newField(new MultiTensorField3D<T, nDim>(
        MultiBlockManagement3D(
            newBlockStructure,
            originalField.getMultiBlockManagement().getThreadAttribution().clone(),
            originalField.getMultiBlockManagement().getEnvelopeWidth(),
            originalField.getMultiBlockManagement().getRefinementLevel()),
        originalField.getBlockCommunicator().clone(), originalField.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>()));

    if (adjustPeriodicity) {
        newField->periodicity().toggle(0, originalField.periodicity().get(0));
        newField->periodicity().toggle(1, originalField.periodicity().get(1));
        newField->periodicity().toggle(2, originalField.periodicity().get(2));
    }

    transferTensorFieldNonLocal(originalField, *newField, originalField.getBoundingBox());

    return newField;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > redistribute(
    MultiTensorField3D<T, nDim> const &originalField,
    SparseBlockStructure3D const &newBlockStructure, Box3D const &intersection, bool crop)
{
    bool adjustPeriodicity = false;
    return redistribute(
        originalField, intersect(newBlockStructure, intersection, crop), adjustPeriodicity);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > align(
    MultiTensorField3D<T, nDim> const &originalBlock, MultiBlock3D const &partnerBlock)
{
    std::unique_ptr<MultiTensorField3D<T, nDim> > newBlock(new MultiTensorField3D<T, nDim>(
        align(originalBlock.getMultiBlockManagement(), partnerBlock.getMultiBlockManagement()),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>()));

    newBlock->periodicity().toggle(0, originalBlock.periodicity().get(0));
    newBlock->periodicity().toggle(1, originalBlock.periodicity().get(1));
    newBlock->periodicity().toggle(2, originalBlock.periodicity().get(2));

    transferTensorFieldNonLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > reparallelize(
    MultiTensorField3D<T, nDim> const &originalBlock)
{
    return reparallelize(originalBlock, 16, 16, 16);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > reparallelize(
    MultiTensorField3D<T, nDim> const &originalBlock, plint blockLx, plint blockLy, plint blockLz)
{
    std::unique_ptr<MultiTensorField3D<T, nDim> > newBlock(new MultiTensorField3D<T, nDim>(
        reparallelize(originalBlock.getMultiBlockManagement(), blockLx, blockLy, blockLz),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>()));

    newBlock->periodicity().toggle(0, originalBlock.periodicity().get(0));
    newBlock->periodicity().toggle(1, originalBlock.periodicity().get(1));
    newBlock->periodicity().toggle(2, originalBlock.periodicity().get(2));

    transferTensorFieldNonLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

/* *************** 4. MultiBlockLattice ************************************** */

// QUESTION: Does it make sense to have a useless domain here?
template <typename T, template <typename U> class Descriptor>
void transferBlockLatticeLocal(
    MultiBlockLattice3D<T, Descriptor> &from, MultiBlockLattice3D<T, Descriptor> &to,
    [[maybe_unused]] Box3D const &domain)
{
    // 1. Copy static and dynamic data to the new block.
    copyRegenerate(from, to, from.getBoundingBox());

    // 2. Reconstruct the data processors.
    transferDataProcessors(from, to);
}

template <typename T, template <typename U> class Descriptor>
void transferBlockLatticeNonLocal(
    MultiBlockLattice3D<T, Descriptor> const &from, MultiBlockLattice3D<T, Descriptor> &to,
    Box3D const &domain)
{
    // 1. Copy all data from the old to the new field. This includes dynamics
    //    objects which must be fully serialized and regenerated.
    copyNonLocal(from, to, domain, modif::dataStructure);
    // 2. Reconstruct the data processors.
    transferDataProcessors(from, to);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > generateMultiBlockLattice(
    Box3D boundingBox, Dynamics<T, Descriptor> *backgroundDynamics, plint envelopeWidth)
{
    return std::unique_ptr<MultiBlockLattice3D<T, Descriptor> >(
        new MultiBlockLattice3D<T, Descriptor>(
            defaultMultiBlockPolicy3D().getMultiBlockManagement(boundingBox, envelopeWidth),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(), backgroundDynamics));
}

// QUESTION: Why do we have this unnamedDummyArg?
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > defaultGenerateMultiBlockLattice3D(
    MultiBlockManagement3D const &management, [[maybe_unused]] plint unnamedDummyArg)
{
    return std::unique_ptr<MultiBlockLattice3D<T, Descriptor> >(
        new MultiBlockLattice3D<T, Descriptor>(
            management, defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(),
            new NoDynamics<T, Descriptor>));
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > generateMultiBlockLattice(
    MultiBlockManagement3D const &management, Dynamics<T, Descriptor> *backgroundDynamics)
{
    if (backgroundDynamics == 0) {
        backgroundDynamics = new NoDynamics<T, Descriptor>();
    }
    return std::unique_ptr<MultiBlockLattice3D<T, Descriptor> >(
        new MultiBlockLattice3D<T, Descriptor>(
            management, defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(), backgroundDynamics));
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > clone(
    MultiBlockLattice3D<T, Descriptor> &originalLattice, Box3D const &subDomain, bool crop)
{
    std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > clonedLattice(
        generateMultiBlockLattice<T, Descriptor>(originalLattice, subDomain, crop));

    transferBlockLatticeLocal(originalLattice, *clonedLattice, originalLattice.getBoundingBox());

    return clonedLattice;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > generateMultiBlockLattice(
    MultiBlock3D &multiBlock, plint envelopeWidth, Dynamics<T, Descriptor> *backgroundDynamics)
{
    if (backgroundDynamics == 0) {
        backgroundDynamics = new NoDynamics<T, Descriptor>();
    }
    MultiBlockManagement3D sparseBlockManagement(multiBlock.getMultiBlockManagement());
    MultiBlockLattice3D<T, Descriptor> *newBlock = new MultiBlockLattice3D<T, Descriptor>(
        MultiBlockManagement3D(
            sparseBlockManagement.getSparseBlockStructure(),
            sparseBlockManagement.getThreadAttribution().clone(), envelopeWidth,
            sparseBlockManagement.getRefinementLevel()),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(), backgroundDynamics);

    newBlock->periodicity().toggle(0, multiBlock.periodicity().get(0));
    newBlock->periodicity().toggle(1, multiBlock.periodicity().get(1));
    newBlock->periodicity().toggle(2, multiBlock.periodicity().get(2));

    return std::unique_ptr<MultiBlockLattice3D<T, Descriptor> >(newBlock);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > generateMultiBlockLattice(
    MultiBlock3D const &originalBlock, Box3D const &intersection, bool crop)
{
    return std::unique_ptr<MultiBlockLattice3D<T, Descriptor> >(
        new MultiBlockLattice3D<T, Descriptor>(
            intersect(originalBlock.getMultiBlockManagement(), intersection, crop),
            originalBlock.getBlockCommunicator().clone(),
            originalBlock.getCombinedStatistics().clone(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(),
            new NoDynamics<T, Descriptor>));
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > generateIntersectMultiBlockLattice(
    MultiBlock3D const &originalBlock1, MultiBlock3D const &originalBlock2, bool crop)
{
    return std::unique_ptr<MultiBlockLattice3D<T, Descriptor> >(
        new MultiBlockLattice3D<T, Descriptor>(
            intersect(
                originalBlock1.getMultiBlockManagement(), originalBlock2.getMultiBlockManagement(),
                crop),
            originalBlock1.getBlockCommunicator().clone(),
            originalBlock1.getCombinedStatistics().clone(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(),
            new NoDynamics<T, Descriptor>));
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > generateIntersectMultiBlockLattice(
    MultiBlock3D const &originalBlock1, MultiBlock3D const &originalBlock2,
    Box3D const &intersection, bool crop)
{
    MultiBlockManagement3D intersectedBlocks(intersect(
        originalBlock1.getMultiBlockManagement(), originalBlock2.getMultiBlockManagement(), crop));
    MultiBlockManagement3D intersectWithDomain(intersect(intersectedBlocks, intersection, crop));
    return std::unique_ptr<MultiBlockLattice3D<T, Descriptor> >(
        new MultiBlockLattice3D<T, Descriptor>(
            intersectWithDomain, originalBlock1.getBlockCommunicator().clone(),
            originalBlock1.getCombinedStatistics().clone(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(),
            new NoDynamics<T, Descriptor>));
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > generateJoinMultiBlockLattice(
    MultiBlock3D const &originalBlock1, MultiBlock3D const &originalBlock2)
{
    return std::unique_ptr<MultiBlockLattice3D<T, Descriptor> >(
        new MultiBlockLattice3D<T, Descriptor>(
            block_union(
                originalBlock1.getMultiBlockManagement(), originalBlock2.getMultiBlockManagement()),
            originalBlock1.getBlockCommunicator().clone(),
            originalBlock1.getCombinedStatistics().clone(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(),
            new NoDynamics<T, Descriptor>));
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > extend(
    MultiBlockLattice3D<T, Descriptor> &originalBlock, Box3D const &addedBlock)
{
    std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > newBlock(
        new MultiBlockLattice3D<T, Descriptor>(
            extend(originalBlock.getMultiBlockManagement(), addedBlock, addedBlock),
            originalBlock.getBlockCommunicator().clone(),
            originalBlock.getCombinedStatistics().clone(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(),
            originalBlock.getBackgroundDynamics().clone()));

    transferBlockLatticeLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > except(
    MultiBlockLattice3D<T, Descriptor> &originalBlock, Box3D const &exceptedBlock)
{
    std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > newBlock(
        new MultiBlockLattice3D<T, Descriptor>(
            except(originalBlock.getMultiBlockManagement(), exceptedBlock),
            originalBlock.getBlockCommunicator().clone(),
            originalBlock.getCombinedStatistics().clone(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(),
            originalBlock.getBackgroundDynamics().clone()));

    transferBlockLatticeLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > redistribute(
    MultiBlockLattice3D<T, Descriptor> const &originalBlock,
    SparseBlockStructure3D const &newBlockStructure, bool adjustPeriodicity)
{
    std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > newBlock(
        new MultiBlockLattice3D<T, Descriptor>(
            MultiBlockManagement3D(
                newBlockStructure,
                originalBlock.getMultiBlockManagement().getThreadAttribution().clone(),
                originalBlock.getMultiBlockManagement().getEnvelopeWidth(),
                originalBlock.getMultiBlockManagement().getRefinementLevel()),
            originalBlock.getBlockCommunicator().clone(),
            originalBlock.getCombinedStatistics().clone(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(),
            originalBlock.getBackgroundDynamics().clone()));

    if (adjustPeriodicity) {
        newBlock->periodicity().toggle(0, originalBlock.periodicity().get(0));
        newBlock->periodicity().toggle(1, originalBlock.periodicity().get(1));
        newBlock->periodicity().toggle(2, originalBlock.periodicity().get(2));
    }

    transferBlockLatticeNonLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > redistribute(
    MultiBlockLattice3D<T, Descriptor> const &originalBlock,
    SparseBlockStructure3D const &newBlockStructure, Box3D const &intersection, bool crop)
{
    bool adjustPeriodicity = false;
    return redistribute(
        originalBlock, intersect(newBlockStructure, intersection, crop), adjustPeriodicity);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > align(
    MultiBlockLattice3D<T, Descriptor> const &originalBlock, MultiBlock3D const &partnerBlock)
{
    std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > newBlock(
        new MultiBlockLattice3D<T, Descriptor>(
            align(originalBlock.getMultiBlockManagement(), partnerBlock.getMultiBlockManagement()),
            originalBlock.getBlockCommunicator().clone(),
            originalBlock.getCombinedStatistics().clone(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(),
            originalBlock.getBackgroundDynamics().clone()));

    newBlock->periodicity().toggle(0, originalBlock.periodicity().get(0));
    newBlock->periodicity().toggle(1, originalBlock.periodicity().get(1));
    newBlock->periodicity().toggle(2, originalBlock.periodicity().get(2));

    transferBlockLatticeNonLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > reparallelize(
    MultiBlockLattice3D<T, Descriptor> const &originalBlock)
{
    return reparallelize(originalBlock, 16, 16, 16);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > reparallelize(
    MultiBlockLattice3D<T, Descriptor> const &originalBlock, plint blockLx, plint blockLy,
    plint blockLz)
{
    std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > newBlock(
        new MultiBlockLattice3D<T, Descriptor>(
            reparallelize(originalBlock.getMultiBlockManagement(), blockLx, blockLy, blockLz),
            originalBlock.getBlockCommunicator().clone(),
            originalBlock.getCombinedStatistics().clone(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(),
            originalBlock.getBackgroundDynamics().clone()));

    newBlock->periodicity().toggle(0, originalBlock.periodicity().get(0));
    newBlock->periodicity().toggle(1, originalBlock.periodicity().get(1));
    newBlock->periodicity().toggle(2, originalBlock.periodicity().get(2));

    transferBlockLatticeNonLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

/* *************** 4. MultiParticleField ************************************** */

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiParticleField3D<DenseParticleField3D<T, Descriptor> > >
    generateMultiDenseParticleField(Box3D boundingBox, plint envelopeWidth)
{
    return std::unique_ptr<MultiParticleField3D<DenseParticleField3D<T, Descriptor> > >(
        new MultiParticleField3D<DenseParticleField3D<T, Descriptor> >(
            defaultMultiBlockPolicy3D().getMultiBlockManagement(boundingBox, envelopeWidth),
            defaultMultiBlockPolicy3D().getCombinedStatistics()));
}

template <typename T, template <typename U> class Descriptor, class ParticleFieldT>
std::unique_ptr<MultiParticleField3D<ParticleFieldT> > generateMultiParticleField3D(
    Box3D boundingBox, plint envelopeWidth)
{
    return std::unique_ptr<MultiParticleField3D<ParticleFieldT> >(
        new MultiParticleField3D<ParticleFieldT>(
            defaultMultiBlockPolicy3D().getMultiBlockManagement(boundingBox, envelopeWidth),
            defaultMultiBlockPolicy3D().getCombinedStatistics()));
}

template <typename T, template <typename U> class Descriptor, class ParticleFieldT>
std::unique_ptr<MultiParticleField3D<ParticleFieldT> > generateMultiParticleField3D(
    MultiBlock3D &multiBlock, plint envelopeWidth)
{
    MultiBlockManagement3D sparseBlockManagement(multiBlock.getMultiBlockManagement());
    MultiParticleField3D<ParticleFieldT> *field = new MultiParticleField3D<ParticleFieldT>(
        MultiBlockManagement3D(
            sparseBlockManagement.getSparseBlockStructure(),
            sparseBlockManagement.getThreadAttribution().clone(), envelopeWidth,
            sparseBlockManagement.getRefinementLevel()),
        defaultMultiBlockPolicy3D().getCombinedStatistics());

    field->periodicity().toggle(0, multiBlock.periodicity().get(0));
    field->periodicity().toggle(1, multiBlock.periodicity().get(1));
    field->periodicity().toggle(2, multiBlock.periodicity().get(2));

    return std::unique_ptr<MultiParticleField3D<ParticleFieldT> >(field);
}

template <typename T, template <typename U> class Descriptor, class ParticleFieldT>
std::unique_ptr<MultiParticleField3D<ParticleFieldT> > generateMultiParticleField3D(
    MultiBlockManagement3D const &management, PeriodicitySwitch3D const &periodicity,
    plint envelopeWidth)
{
    MultiBlockManagement3D sparseBlockManagement(management);
    MultiParticleField3D<ParticleFieldT> *field = new MultiParticleField3D<ParticleFieldT>(
        MultiBlockManagement3D(
            sparseBlockManagement.getSparseBlockStructure(),
            sparseBlockManagement.getThreadAttribution().clone(), envelopeWidth,
            sparseBlockManagement.getRefinementLevel()),
        defaultMultiBlockPolicy3D().getCombinedStatistics());

    field->periodicity().toggle(0, periodicity.get(0));
    field->periodicity().toggle(1, periodicity.get(1));
    field->periodicity().toggle(2, periodicity.get(2));

    return std::unique_ptr<MultiParticleField3D<ParticleFieldT> >(field);
}

template <typename T, template <typename U> class Descriptor, class ParticleFieldT>
std::unique_ptr<MultiParticleField3D<ParticleFieldT> > generateMultiParticleField3D(
    MultiBlockManagement3D const &management, plint envelopeWidth)
{
    MultiBlockManagement3D sparseBlockManagement(management);
    MultiParticleField3D<ParticleFieldT> *field = new MultiParticleField3D<ParticleFieldT>(
        MultiBlockManagement3D(
            sparseBlockManagement.getSparseBlockStructure(),
            sparseBlockManagement.getThreadAttribution().clone(), envelopeWidth,
            sparseBlockManagement.getRefinementLevel()),
        defaultMultiBlockPolicy3D().getCombinedStatistics());

    return std::unique_ptr<MultiParticleField3D<ParticleFieldT> >(field);
}

}  // namespace plb

#endif  // MULTI_BLOCK_GENERATOR_3D_HH
