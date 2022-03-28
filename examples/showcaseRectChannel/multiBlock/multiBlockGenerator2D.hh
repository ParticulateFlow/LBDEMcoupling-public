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
 * Copy 2D multiblocks on a new parallel distribution -- generic implementation.
 */

/*
 * All functions that take as an argument a Box2D, or compute intersections, or join,
 * or crop, do not adjust periodicity of the newly created blocks by default. The user
 * is responsible to take care of this matter explicitly.
 */

#ifndef MULTI_BLOCK_GENERATOR_2D_HH
#define MULTI_BLOCK_GENERATOR_2D_HH

#include "core/globalDefs.h"
#include "dataProcessors/dataAnalysisWrapper2D.h"
#include "dataProcessors/ntensorAnalysisWrapper2D.h"
#include "multiBlock/defaultMultiBlockPolicy2D.h"
#include "multiBlock/localMultiBlockInfo2D.h"
#include "multiBlock/multiBlockOperations2D.h"
#include "multiBlock/nonLocalTransfer2D.h"
#include "multiBlock/sparseBlockStructure2D.h"

namespace plb {

inline void transferDataProcessors(MultiBlock2D const &from, MultiBlock2D &to)
{
    std::vector<MultiBlock2D::ProcessorStorage2D> newProcessors(from.getStoredProcessors());
    // Redirect all references-to-self to the new self.
    id_t oldId = from.getId();
    id_t newId = to.getId();
    for (pluint iProcessor = 0; iProcessor < newProcessors.size(); ++iProcessor) {
        newProcessors[iProcessor].replace(oldId, newId);
        DataProcessorGenerator2D *newGenerator = newProcessors[iProcessor].getGenerator().clone();
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
    MultiScalarField2D<T> &from, MultiScalarField2D<T> &to, Box2D const &domain)
{
    // 1. Copy all data from the old to the new field.
    plb::copy(from, to, from.getBoundingBox());
    // 2. Reconstruct the data processors.
    transferDataProcessors(from, to);
}

template <typename T>
void transferScalarFieldNonLocal(
    MultiScalarField2D<T> const &from, MultiScalarField2D<T> &to, Box2D const &domain)
{
    // 1. Copy all data from the old to the new field. This includes dynamics
    //    objects which must be fully serialized and regenerated.
    copyNonLocal(from, to, domain);
    // 2. Reconstruct the data processors.
    transferDataProcessors(from, to);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > generateMultiScalarField(
    Box2D boundingBox, plint envelopeWidth)
{
    return std::unique_ptr<MultiScalarField2D<T> >(new MultiScalarField2D<T>(
        defaultMultiBlockPolicy2D().getMultiBlockManagement(boundingBox, envelopeWidth),
        defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiScalarAccess<T>()));
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > generateMultiScalarField(
    Box2D boundingBox, T iniVal, plint envelopeWidth)
{
    return std::unique_ptr<MultiScalarField2D<T> >(new MultiScalarField2D<T>(
        defaultMultiBlockPolicy2D().getMultiBlockManagement(boundingBox, envelopeWidth),
        defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiScalarAccess<T>(), iniVal));
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > generateMultiScalarField(
    MultiBlock2D &multiBlock, plint envelopeWidth)
{
    MultiBlockManagement2D sparseBlockManagement(multiBlock.getMultiBlockManagement());
    MultiScalarField2D<T> *field = new MultiScalarField2D<T>(
        MultiBlockManagement2D(
            sparseBlockManagement.getSparseBlockStructure(),
            sparseBlockManagement.getThreadAttribution().clone(), envelopeWidth,
            sparseBlockManagement.getRefinementLevel()),
        defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiScalarAccess<T>());

    field->periodicity().toggle(0, multiBlock.periodicity().get(0));
    field->periodicity().toggle(1, multiBlock.periodicity().get(1));

    return std::unique_ptr<MultiScalarField2D<T> >(field);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > defaultGenerateMultiScalarField2D(
    MultiBlockManagement2D const &management, T iniVal)
{
    return std::unique_ptr<MultiScalarField2D<T> >(new MultiScalarField2D<T>(
        management, defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiScalarAccess<T>(), iniVal));
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > clone(
    MultiScalarField2D<T> &originalField, Box2D const &subDomain, bool crop)
{
    std::unique_ptr<MultiScalarField2D<T> > clonedField(
        generateMultiScalarField<T>(originalField, subDomain, crop));

    transferScalarFieldLocal(originalField, *clonedField, originalField.getBoundingBox());

    return clonedField;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > generateMultiScalarField(
    MultiBlock2D const &originalField, Box2D const &intersection, bool crop)
{
    return std::unique_ptr<MultiScalarField2D<T> >(new MultiScalarField2D<T>(
        intersect(originalField.getMultiBlockManagement(), intersection, crop),
        originalField.getBlockCommunicator().clone(), originalField.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiScalarAccess<T>()));
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > generateIntersectMultiScalarField(
    MultiBlock2D const &originalField1, MultiBlock2D const &originalField2, bool crop)
{
    return std::unique_ptr<MultiScalarField2D<T> >(new MultiScalarField2D<T>(
        intersect(
            originalField1.getMultiBlockManagement(), originalField2.getMultiBlockManagement(),
            crop),
        originalField1.getBlockCommunicator().clone(),
        originalField1.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiScalarAccess<T>()));
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > generateIntersectMultiScalarField(
    MultiBlock2D const &originalField1, MultiBlock2D const &originalField2,
    Box2D const &intersection, bool crop)
{
    MultiBlockManagement2D intersectedBlocks(intersect(
        originalField1.getMultiBlockManagement(), originalField2.getMultiBlockManagement(), crop));
    MultiBlockManagement2D intersectWithDomain(intersect(intersectedBlocks, intersection, crop));
    return std::unique_ptr<MultiScalarField2D<T> >(new MultiScalarField2D<T>(
        intersectWithDomain, originalField1.getBlockCommunicator().clone(),
        originalField1.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiScalarAccess<T>()));
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > generateJoinMultiScalarField(
    MultiBlock2D const &originalField1, MultiBlock2D const &originalField2)
{
    return std::unique_ptr<MultiScalarField2D<T> >(new MultiScalarField2D<T>(
        block_union(
            originalField1.getMultiBlockManagement(), originalField2.getMultiBlockManagement()),
        originalField1.getBlockCommunicator().clone(),
        originalField1.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiScalarAccess<T>()));
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > extend(
    MultiScalarField2D<T> &originalBlock, Box2D const &addedBlock)
{
    std::unique_ptr<MultiScalarField2D<T> > newBlock(new MultiScalarField2D<T>(
        extend(originalBlock.getMultiBlockManagement(), addedBlock, addedBlock),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiScalarAccess<T>()));

    transferScalarFieldLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > except(
    MultiScalarField2D<T> &originalBlock, Box2D const &exceptedBlock)
{
    std::unique_ptr<MultiScalarField2D<T> > newBlock(new MultiScalarField2D<T>(
        except(originalBlock.getMultiBlockManagement(), exceptedBlock),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiScalarAccess<T>()));

    transferScalarFieldLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > redistribute(
    MultiScalarField2D<T> const &originalField, SparseBlockStructure2D const &newBlockStructure,
    bool adjustPeriodicity)
{
    std::unique_ptr<MultiScalarField2D<T> > newField(new MultiScalarField2D<T>(
        MultiBlockManagement2D(
            newBlockStructure,
            originalField.getMultiBlockManagement().getThreadAttribution().clone(),
            originalField.getMultiBlockManagement().getEnvelopeWidth(),
            originalField.getMultiBlockManagement().getRefinementLevel()),
        originalField.getBlockCommunicator().clone(), originalField.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiScalarAccess<T>()));

    if (adjustPeriodicity) {
        newField->periodicity().toggle(0, originalField.periodicity().get(0));
        newField->periodicity().toggle(1, originalField.periodicity().get(1));
    }

    transferScalarFieldNonLocal(originalField, *newField, originalField.getBoundingBox());

    return newField;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > redistribute(
    MultiScalarField2D<T> const &originalField, SparseBlockStructure2D const &newBlockStructure,
    Box2D const &intersection, bool crop)
{
    bool adjustPeriodicity = false;
    return redistribute(
        originalField, intersect(newBlockStructure, intersection, crop), adjustPeriodicity);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > align(
    MultiScalarField2D<T> const &originalBlock, MultiBlock2D const &partnerBlock)
{
    std::unique_ptr<MultiScalarField2D<T> > newBlock(new MultiScalarField2D<T>(
        align(originalBlock.getMultiBlockManagement(), partnerBlock.getMultiBlockManagement()),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiScalarAccess<T>()));

    newBlock->periodicity().toggle(0, originalBlock.periodicity().get(0));
    newBlock->periodicity().toggle(1, originalBlock.periodicity().get(1));

    transferScalarFieldNonLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > reparallelize(MultiScalarField2D<T> const &originalBlock)
{
    return reparallelize(originalBlock, 64, 64);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > reparallelize(
    MultiScalarField2D<T> const &originalBlock, plint blockLx, plint blockLy)
{
    std::unique_ptr<MultiScalarField2D<T> > newBlock(new MultiScalarField2D<T>(
        reparallelize(originalBlock.getMultiBlockManagement(), blockLx, blockLy),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiScalarAccess<T>()));

    newBlock->periodicity().toggle(0, originalBlock.periodicity().get(0));
    newBlock->periodicity().toggle(1, originalBlock.periodicity().get(1));

    transferScalarFieldNonLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

/* *************** 2. MultiNTensorField ************************************** */

template <typename T>
std::unique_ptr<MultiNTensorField2D<T> > defaultGenerateMultiNTensorField2D(
    MultiBlockManagement2D const &management, plint nDim)
{
    MultiNTensorField2D<T> *field = new MultiNTensorField2D<T>(
        nDim, management, defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiNTensorAccess<T>());
    return std::unique_ptr<MultiNTensorField2D<T> >(field);
}

template <typename T>
MultiNTensorField2D<T> *generateMultiNTensorField2D(
    MultiBlock2D &multiBlock, plint envelopeWidth, plint ndim)
{
    MultiBlockManagement2D sparseBlockManagement(multiBlock.getMultiBlockManagement());
    MultiNTensorField2D<T> *field = new MultiNTensorField2D<T>(
        ndim,
        MultiBlockManagement2D(
            sparseBlockManagement.getSparseBlockStructure(),
            sparseBlockManagement.getThreadAttribution().clone(), envelopeWidth,
            sparseBlockManagement.getRefinementLevel()),
        defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiNTensorAccess<T>());

    field->periodicity().toggle(0, multiBlock.periodicity().get(0));
    field->periodicity().toggle(1, multiBlock.periodicity().get(1));

    return field;
}

template <typename T>
MultiNTensorField2D<T> *generateMultiNTensorField2D(Box2D const &domain, plint ndim)
{
    plint defaultEnvelopeWidth = 1;
    MultiNTensorField2D<T> *field = new MultiNTensorField2D<T>(
        ndim, defaultMultiBlockPolicy2D().getMultiBlockManagement(domain, defaultEnvelopeWidth),
        defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiNTensorAccess<T>());
    return field;
}

template <typename T>
MultiNTensorField2D<T> *generateMultiNTensorField2D(
    Box2D const &domain, plint ndim, T *iniVal, plint envelopeWidth)
{
    MultiNTensorField2D<T> *field = new MultiNTensorField2D<T>(
        ndim, iniVal, defaultMultiBlockPolicy2D().getMultiBlockManagement(domain, envelopeWidth),
        defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiNTensorAccess<T>());
    return field;
}

template <typename T>
void transferNTensorFieldLocal(
    MultiNTensorField2D<T> &from, MultiNTensorField2D<T> &to, Box2D const &domain)
// TODO: Remove the domain argument. Not used and can be misleading.
{
    // 1. Copy all data from the old to the new field.
    plb::copy(from, to, from.getBoundingBox());
    // 2. Reconstruct the data processors.
    transferDataProcessors(from, to);
}

template <typename T>
void transferNTensorFieldNonLocal(
    MultiNTensorField2D<T> const &from, MultiNTensorField2D<T> &to, Box2D const &domain)
{
    // 1. Copy all data from the old to the new field. This includes dynamics
    //    objects which must be fully serialized and regenerated.
    copyNonLocal(from, to, domain);
    // 2. Reconstruct the data processors.
    transferDataProcessors(from, to);
}

template <typename T>
MultiNTensorField2D<T> *clone(
    MultiNTensorField2D<T> &originalField, Box2D const &subDomain, bool crop)
{
    MultiNTensorField2D<T> *clonedField =
        generateMultiNTensorField<T>(originalField, subDomain, originalField.getNdim(), crop);

    transferNTensorFieldLocal(originalField, *clonedField, originalField.getBoundingBox());

    return clonedField;
}

template <typename T>
MultiNTensorField2D<T> *generateMultiNTensorField(
    MultiBlock2D const &originalField, Box2D const &intersection, plint nDim, bool crop)
{
    MultiNTensorField2D<T> *newField = new MultiNTensorField2D<T>(
        nDim, intersect(originalField.getMultiBlockManagement(), intersection, crop),
        originalField.getBlockCommunicator().clone(), originalField.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiNTensorAccess<T>());
    return newField;
}

template <typename T1, typename T2>
MultiNTensorField2D<T2> *generateNTensorFieldFromNTensor2D(
    MultiNTensorField2D<T1> const &field, Box2D const &intersection, plint nDim)
{
    MultiNTensorField2D<T2> *newField = generateMultiNTensorField<T2>(field, intersection, nDim);
    return newField;
}

template <typename T1, typename T2, template <typename U> class Descriptor>
MultiNTensorField2D<T1> *generateNTensorFieldFromBlockLattice2D(
    MultiBlockLattice2D<T2, Descriptor> const &lattice, Box2D const &intersection, plint nDim)
{
    MultiNTensorField2D<T1> *newField = generateMultiNTensorField<T1>(lattice, intersection, nDim);
    return newField;
}

template <typename T>
MultiNTensorField2D<T> *generateIntersectMultiNTensorField(
    MultiBlock2D const &originalField1, MultiBlock2D const &originalField2, plint nDim, bool crop)
{
    MultiNTensorField2D<T> *newField = new MultiNTensorField2D<T>(
        nDim,
        intersect(
            originalField1.getMultiBlockManagement(), originalField2.getMultiBlockManagement(),
            crop),
        originalField1.getBlockCommunicator().clone(),
        originalField1.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiNTensorAccess<T>());
    return newField;
}

template <typename T>
MultiNTensorField2D<T> *generateIntersectMultiNTensorField(
    MultiBlock2D const &originalField1, MultiBlock2D const &originalField2,
    Box2D const &intersection, plint nDim, bool crop)
{
    MultiBlockManagement2D intersectedBlocks(intersect(
        originalField1.getMultiBlockManagement(), originalField2.getMultiBlockManagement(), crop));
    MultiBlockManagement2D intersectWithDomain(intersect(intersectedBlocks, intersection, crop));
    MultiNTensorField2D<T> *newField = new MultiNTensorField2D<T>(
        nDim, intersectWithDomain, originalField1.getBlockCommunicator().clone(),
        originalField1.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiNTensorAccess<T>());
    return newField;
}

template <typename T>
MultiNTensorField2D<T> *generateJoinMultiNTensorField(
    MultiBlock2D const &originalField1, MultiBlock2D const &originalField2, plint nDim)
{
    MultiNTensorField2D<T> *newField = new MultiNTensorField2D<T>(
        nDim,
        block_union(
            originalField1.getMultiBlockManagement(), originalField2.getMultiBlockManagement()),
        originalField1.getBlockCommunicator().clone(),
        originalField1.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiNTensorAccess<T>());
    return newField;
}

template <typename T>
MultiNTensorField2D<T> *extend(MultiNTensorField2D<T> &originalBlock, Box2D const &addedBlock)
{
    MultiNTensorField2D<T> *newBlock = new MultiNTensorField2D<T>(
        originalBlock.getNdim(),
        extend(originalBlock.getMultiBlockManagement(), addedBlock, addedBlock),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiNTensorAccess<T>());

    transferNTensorFieldLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T>
MultiNTensorField2D<T> *except(MultiNTensorField2D<T> &originalBlock, Box2D const &exceptedBlock)
{
    MultiNTensorField2D<T> *newBlock = new MultiNTensorField2D<T>(
        originalBlock.getNdim(), except(originalBlock.getMultiBlockManagement(), exceptedBlock),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiNTensorAccess<T>());

    transferNTensorFieldLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T>
MultiNTensorField2D<T> *align(
    MultiNTensorField2D<T> const &originalBlock, MultiBlock2D const &partnerBlock)
{
    MultiNTensorField2D<T> *newBlock = new MultiNTensorField2D<T>(
        originalBlock.getNdim(),
        align(originalBlock.getMultiBlockManagement(), partnerBlock.getMultiBlockManagement()),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiNTensorAccess<T>());

    newBlock->periodicity().toggle(0, originalBlock.periodicity().get(0));
    newBlock->periodicity().toggle(1, originalBlock.periodicity().get(1));

    transferNTensorFieldNonLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T>
MultiNTensorField2D<T> *reparallelize(MultiNTensorField2D<T> const &originalBlock)
{
    MultiNTensorField2D<T> *newBlock = new MultiNTensorField2D<T>(
        originalBlock.getNdim(), reparallelize(originalBlock.getMultiBlockManagement(), 64, 64),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiNTensorAccess<T>());

    newBlock->periodicity().toggle(0, originalBlock.periodicity().get(0));
    newBlock->periodicity().toggle(1, originalBlock.periodicity().get(1));

    transferNTensorFieldNonLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

/* *************** 3. MultiTensorField ************************************** */

template <typename T, int nDim>
void transferTensorFieldLocal(
    MultiTensorField2D<T, nDim> &from, MultiTensorField2D<T, nDim> &to, Box2D const &domain)
// TODO: remove the domain argument. It is not used and can be misleading.
{
    // 1. Copy all data from the old to the new field.
    plb::copy(from, to, from.getBoundingBox());
    // 2. Reconstruct the data processors.
    transferDataProcessors(from, to);
}

template <typename T, int nDim>
void transferTensorFieldNonLocal(
    MultiTensorField2D<T, nDim> const &from, MultiTensorField2D<T, nDim> &to, Box2D const &domain)
{
    // 1. Copy all data from the old to the new field. This includes dynamics
    //    objects which must be fully serialized and regenerated.
    copyNonLocal(from, to, domain);
    // 2. Reconstruct the data processors.
    transferDataProcessors(from, to);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > generateMultiTensorField(
    Box2D boundingBox, plint envelopeWidth)
{
    return std::unique_ptr<MultiTensorField2D<T, nDim> >(new MultiTensorField2D<T, nDim>(
        defaultMultiBlockPolicy2D().getMultiBlockManagement(boundingBox, envelopeWidth),
        defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiTensorAccess<T, nDim>()));
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > generateMultiTensorField(
    Box2D boundingBox, Array<T, nDim> const &iniVal, plint envelopeWidth)
{
    return std::unique_ptr<MultiTensorField2D<T, nDim> >(new MultiTensorField2D<T, nDim>(
        defaultMultiBlockPolicy2D().getMultiBlockManagement(boundingBox, envelopeWidth),
        defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiTensorAccess<T, nDim>(), iniVal));
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > generateMultiTensorField(
    MultiBlock2D &multiBlock, plint envelopeWidth)
{
    MultiBlockManagement2D sparseBlockManagement(multiBlock.getMultiBlockManagement());
    MultiTensorField2D<T, nDim> *field = new MultiTensorField2D<T, nDim>(
        MultiBlockManagement2D(
            sparseBlockManagement.getSparseBlockStructure(),
            sparseBlockManagement.getThreadAttribution().clone(), envelopeWidth,
            sparseBlockManagement.getRefinementLevel()),
        defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiTensorAccess<T, nDim>());

    field->periodicity().toggle(0, multiBlock.periodicity().get(0));
    field->periodicity().toggle(1, multiBlock.periodicity().get(1));

    return std::unique_ptr<MultiTensorField2D<T, nDim> >(field);
}

// TODO: why is there an unnamed arg?
template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > defaultGenerateMultiTensorField2D(
    MultiBlockManagement2D const &management, plint unnamedDummyArg)
{
    Array<T, nDim> iniVal;
    iniVal.resetToZero();
    return std::unique_ptr<MultiTensorField2D<T, nDim> >(new MultiTensorField2D<T, nDim>(
        management, defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiTensorAccess<T, nDim>(), iniVal));
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > clone(
    MultiTensorField2D<T, nDim> &originalField, Box2D const &subDomain, bool crop)
{
    std::unique_ptr<MultiTensorField2D<T, nDim> > clonedField(
        generateMultiTensorField<T, nDim>(originalField, subDomain, crop));

    transferTensorFieldLocal(originalField, *clonedField, originalField.getBoundingBox());

    return clonedField;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > generateMultiTensorField(
    MultiBlock2D const &originalField, Box2D const &intersection, bool crop)
{
    return std::unique_ptr<MultiTensorField2D<T, nDim> >(new MultiTensorField2D<T, nDim>(
        intersect(originalField.getMultiBlockManagement(), intersection, crop),
        originalField.getBlockCommunicator().clone(), originalField.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiTensorAccess<T, nDim>()));
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > generateIntersectMultiTensorField(
    MultiBlock2D const &originalField1, MultiBlock2D const &originalField2, bool crop)
{
    return std::unique_ptr<MultiTensorField2D<T, nDim> >(new MultiTensorField2D<T, nDim>(
        intersect(
            originalField1.getMultiBlockManagement(), originalField2.getMultiBlockManagement(),
            crop),
        originalField1.getBlockCommunicator().clone(),
        originalField1.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiTensorAccess<T, nDim>()));
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > generateIntersectMultiTensorField(
    MultiBlock2D const &originalField1, MultiBlock2D const &originalField2,
    Box2D const &intersection, bool crop)
{
    MultiBlockManagement2D intersectedBlocks(intersect(
        originalField1.getMultiBlockManagement(), originalField2.getMultiBlockManagement(), crop));
    MultiBlockManagement2D intersectWithDomain(intersect(intersectedBlocks, intersection, crop));
    return std::unique_ptr<MultiTensorField2D<T, nDim> >(new MultiTensorField2D<T, nDim>(
        intersectWithDomain, originalField1.getBlockCommunicator().clone(),
        originalField1.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiTensorAccess<T, nDim>()));
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > generateJoinMultiTensorField(
    MultiBlock2D const &originalField1, MultiBlock2D const &originalField2)
{
    return std::unique_ptr<MultiTensorField2D<T, nDim> >(new MultiTensorField2D<T, nDim>(
        block_union(
            originalField1.getMultiBlockManagement(), originalField2.getMultiBlockManagement()),
        originalField1.getBlockCommunicator().clone(),
        originalField1.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiTensorAccess<T, nDim>()));
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > extend(
    MultiTensorField2D<T, nDim> &originalBlock, Box2D const &addedBlock)
{
    std::unique_ptr<MultiTensorField2D<T, nDim> > newBlock(new MultiTensorField2D<T, nDim>(
        extend(originalBlock.getMultiBlockManagement(), addedBlock, addedBlock),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiTensorAccess<T, nDim>()));

    transferTensorFieldLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > except(
    MultiTensorField2D<T, nDim> &originalBlock, Box2D const &exceptedBlock)
{
    std::unique_ptr<MultiTensorField2D<T, nDim> > newBlock(new MultiTensorField2D<T, nDim>(
        except(originalBlock.getMultiBlockManagement(), exceptedBlock),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiTensorAccess<T, nDim>()));

    transferTensorFieldLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > redistribute(
    MultiTensorField2D<T, nDim> const &originalField,
    SparseBlockStructure2D const &newBlockStructure, bool adjustPeriodicity)
{
    std::unique_ptr<MultiTensorField2D<T, nDim> > newField(new MultiTensorField2D<T, nDim>(
        MultiBlockManagement2D(
            newBlockStructure,
            originalField.getMultiBlockManagement().getThreadAttribution().clone(),
            originalField.getMultiBlockManagement().getEnvelopeWidth(),
            originalField.getMultiBlockManagement().getRefinementLevel()),
        originalField.getBlockCommunicator().clone(), originalField.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiTensorAccess<T, nDim>()));

    if (adjustPeriodicity) {
        newField->periodicity().toggle(0, originalField.periodicity().get(0));
        newField->periodicity().toggle(1, originalField.periodicity().get(1));
    }

    transferTensorFieldNonLocal(originalField, *newField, originalField.getBoundingBox());

    return newField;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > redistribute(
    MultiTensorField2D<T, nDim> const &originalField,
    SparseBlockStructure2D const &newBlockStructure, Box2D const &intersection, bool crop)
{
    bool adjustPeriodicity = false;
    return redistribute(
        originalField, intersect(newBlockStructure, intersection, crop), adjustPeriodicity);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > align(
    MultiTensorField2D<T, nDim> const &originalBlock, MultiBlock2D const &partnerBlock)
{
    std::unique_ptr<MultiTensorField2D<T, nDim> > newBlock(new MultiTensorField2D<T, nDim>(
        align(originalBlock.getMultiBlockManagement(), partnerBlock.getMultiBlockManagement()),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiTensorAccess<T, nDim>()));

    newBlock->periodicity().toggle(0, originalBlock.periodicity().get(0));
    newBlock->periodicity().toggle(1, originalBlock.periodicity().get(1));

    transferTensorFieldNonLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > reparallelize(
    MultiTensorField2D<T, nDim> const &originalBlock)
{
    return reparallelize(originalBlock, 64, 64);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > reparallelize(
    MultiTensorField2D<T, nDim> const &originalBlock, plint blockLx, plint blockLy)
{
    std::unique_ptr<MultiTensorField2D<T, nDim> > newBlock(new MultiTensorField2D<T, nDim>(
        reparallelize(originalBlock.getMultiBlockManagement(), blockLx, blockLy),
        originalBlock.getBlockCommunicator().clone(), originalBlock.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy2D().getMultiTensorAccess<T, nDim>()));

    newBlock->periodicity().toggle(0, originalBlock.periodicity().get(0));
    newBlock->periodicity().toggle(1, originalBlock.periodicity().get(1));

    transferTensorFieldNonLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

/* *************** 4. MultiBlockLattice ************************************** */

// QUESTION: Does it make sense to have a useless domain here?
template <typename T, template <typename U> class Descriptor>
void transferBlockLatticeLocal(
    MultiBlockLattice2D<T, Descriptor> &from, MultiBlockLattice2D<T, Descriptor> &to,
    [[maybe_unused]] Box2D const &domain)
{
    // 1. Copy static and dynamic data to the new block.
    copyRegenerate(from, to, from.getBoundingBox());

    // 2. Reconstruct the data processors.
    transferDataProcessors(from, to);
}

template <typename T, template <typename U> class Descriptor>
void transferBlockLatticeNonLocal(
    MultiBlockLattice2D<T, Descriptor> const &from, MultiBlockLattice2D<T, Descriptor> &to,
    Box2D const &domain)
{
    // 1. Copy all data from the old to the new field. This includes dynamics
    //    objects which must be fully serialized and regenerated.
    copyNonLocal(from, to, domain, modif::dataStructure);
    // 2. Reconstruct the data processors.
    transferDataProcessors(from, to);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > generateMultiBlockLattice(
    Box2D boundingBox, Dynamics<T, Descriptor> *backgroundDynamics, plint envelopeWidth)
{
    return std::unique_ptr<MultiBlockLattice2D<T, Descriptor> >(
        new MultiBlockLattice2D<T, Descriptor>(
            defaultMultiBlockPolicy2D().getMultiBlockManagement(boundingBox, envelopeWidth),
            defaultMultiBlockPolicy2D().getBlockCommunicator(),
            defaultMultiBlockPolicy2D().getCombinedStatistics(),
            defaultMultiBlockPolicy2D().getMultiCellAccess<T, Descriptor>(), backgroundDynamics));
}

// QUESTION: Why this dummy arg
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > defaultGenerateMultiBlockLattice2D(
    MultiBlockManagement2D const &management, [[maybe_unused]] plint unnamedDummyArg)
{
    return std::unique_ptr<MultiBlockLattice2D<T, Descriptor> >(
        new MultiBlockLattice2D<T, Descriptor>(
            management, defaultMultiBlockPolicy2D().getBlockCommunicator(),
            defaultMultiBlockPolicy2D().getCombinedStatistics(),
            defaultMultiBlockPolicy2D().getMultiCellAccess<T, Descriptor>(),
            new NoDynamics<T, Descriptor>));
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > clone(
    MultiBlockLattice2D<T, Descriptor> &originalLattice, Box2D const &subDomain, bool crop)
{
    std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > clonedLattice(
        generateMultiBlockLattice<T, Descriptor>(originalLattice, subDomain, crop));

    transferBlockLatticeLocal(originalLattice, *clonedLattice, originalLattice.getBoundingBox());

    return clonedLattice;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > generateMultiBlockLattice(
    MultiBlock2D const &originalBlock, Box2D const &intersection, bool crop)
{
    return std::unique_ptr<MultiBlockLattice2D<T, Descriptor> >(
        new MultiBlockLattice2D<T, Descriptor>(
            intersect(originalBlock.getMultiBlockManagement(), intersection, crop),
            originalBlock.getBlockCommunicator().clone(),
            originalBlock.getCombinedStatistics().clone(),
            defaultMultiBlockPolicy2D().getMultiCellAccess<T, Descriptor>(),
            new NoDynamics<T, Descriptor>));
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > generateIntersectMultiBlockLattice(
    MultiBlock2D const &originalBlock1, MultiBlock2D const &originalBlock2, bool crop)
{
    return std::unique_ptr<MultiBlockLattice2D<T, Descriptor> >(
        new MultiBlockLattice2D<T, Descriptor>(
            intersect(
                originalBlock1.getMultiBlockManagement(), originalBlock2.getMultiBlockManagement(),
                crop),
            originalBlock1.getBlockCommunicator().clone(),
            originalBlock1.getCombinedStatistics().clone(),
            defaultMultiBlockPolicy2D().getMultiCellAccess<T, Descriptor>(),
            new NoDynamics<T, Descriptor>));
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > generateIntersectMultiBlockLattice(
    MultiBlock2D const &originalBlock1, MultiBlock2D const &originalBlock2,
    Box2D const &intersection, bool crop)
{
    MultiBlockManagement2D intersectedBlocks(intersect(
        originalBlock1.getMultiBlockManagement(), originalBlock2.getMultiBlockManagement(), crop));
    MultiBlockManagement2D intersectWithDomain(intersect(intersectedBlocks, intersection, crop));
    return std::unique_ptr<MultiBlockLattice2D<T, Descriptor> >(
        new MultiBlockLattice2D<T, Descriptor>(
            intersectWithDomain, originalBlock1.getBlockCommunicator().clone(),
            originalBlock1.getCombinedStatistics().clone(),
            defaultMultiBlockPolicy2D().getMultiCellAccess<T, Descriptor>(),
            new NoDynamics<T, Descriptor>));
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > generateJoinMultiBlockLattice(
    MultiBlock2D const &originalBlock1, MultiBlock2D const &originalBlock2)
{
    return std::unique_ptr<MultiBlockLattice2D<T, Descriptor> >(
        new MultiBlockLattice2D<T, Descriptor>(
            block_union(
                originalBlock1.getMultiBlockManagement(), originalBlock2.getMultiBlockManagement()),
            originalBlock1.getBlockCommunicator().clone(),
            originalBlock1.getCombinedStatistics().clone(),
            defaultMultiBlockPolicy2D().getMultiCellAccess<T, Descriptor>(),
            new NoDynamics<T, Descriptor>));
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > extend(
    MultiBlockLattice2D<T, Descriptor> &originalBlock, Box2D const &addedBlock)
{
    std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > newBlock(
        new MultiBlockLattice2D<T, Descriptor>(
            extend(originalBlock.getMultiBlockManagement(), addedBlock, addedBlock),
            originalBlock.getBlockCommunicator().clone(),
            originalBlock.getCombinedStatistics().clone(),
            defaultMultiBlockPolicy2D().getMultiCellAccess<T, Descriptor>(),
            originalBlock.getBackgroundDynamics().clone()));

    transferBlockLatticeLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > except(
    MultiBlockLattice2D<T, Descriptor> &originalBlock, Box2D const &exceptedBlock)
{
    std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > newBlock(
        new MultiBlockLattice2D<T, Descriptor>(
            except(originalBlock.getMultiBlockManagement(), exceptedBlock),
            originalBlock.getBlockCommunicator().clone(),
            originalBlock.getCombinedStatistics().clone(),
            defaultMultiBlockPolicy2D().getMultiCellAccess<T, Descriptor>(),
            originalBlock.getBackgroundDynamics().clone()));

    transferBlockLatticeLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > redistribute(
    MultiBlockLattice2D<T, Descriptor> const &originalBlock,
    SparseBlockStructure2D const &newBlockStructure, bool adjustPeriodicity)
{
    std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > newBlock(
        new MultiBlockLattice2D<T, Descriptor>(
            MultiBlockManagement2D(
                newBlockStructure,
                originalBlock.getMultiBlockManagement().getThreadAttribution().clone(),
                originalBlock.getMultiBlockManagement().getEnvelopeWidth(),
                originalBlock.getMultiBlockManagement().getRefinementLevel()),
            originalBlock.getBlockCommunicator().clone(),
            originalBlock.getCombinedStatistics().clone(),
            defaultMultiBlockPolicy2D().getMultiCellAccess<T, Descriptor>(),
            originalBlock.getBackgroundDynamics().clone()));

    if (adjustPeriodicity) {
        newBlock->periodicity().toggle(0, originalBlock.periodicity().get(0));
        newBlock->periodicity().toggle(1, originalBlock.periodicity().get(1));
    }

    transferBlockLatticeNonLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > redistribute(
    MultiBlockLattice2D<T, Descriptor> const &originalBlock,
    SparseBlockStructure2D const &newBlockStructure, Box2D const &intersection, bool crop)
{
    bool adjustPeriodicity = false;
    return redistribute(
        originalBlock, intersect(newBlockStructure, intersection, crop), adjustPeriodicity);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > align(
    MultiBlockLattice2D<T, Descriptor> const &originalBlock, MultiBlock2D const &partnerBlock)
{
    std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > newBlock(
        new MultiBlockLattice2D<T, Descriptor>(
            align(originalBlock.getMultiBlockManagement(), partnerBlock.getMultiBlockManagement()),
            originalBlock.getBlockCommunicator().clone(),
            originalBlock.getCombinedStatistics().clone(),
            defaultMultiBlockPolicy2D().getMultiCellAccess<T, Descriptor>(),
            originalBlock.getBackgroundDynamics().clone()));

    newBlock->periodicity().toggle(0, originalBlock.periodicity().get(0));
    newBlock->periodicity().toggle(1, originalBlock.periodicity().get(1));

    transferBlockLatticeNonLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > reparallelize(
    MultiBlockLattice2D<T, Descriptor> const &originalBlock)
{
    return reparallelize(originalBlock, 64, 64);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > reparallelize(
    MultiBlockLattice2D<T, Descriptor> const &originalBlock, plint blockLx, plint blockLy)
{
    std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > newBlock(
        new MultiBlockLattice2D<T, Descriptor>(
            reparallelize(originalBlock.getMultiBlockManagement(), blockLx, blockLy),
            originalBlock.getBlockCommunicator().clone(),
            originalBlock.getCombinedStatistics().clone(),
            defaultMultiBlockPolicy2D().getMultiCellAccess<T, Descriptor>(),
            originalBlock.getBackgroundDynamics().clone()));

    newBlock->periodicity().toggle(0, originalBlock.periodicity().get(0));
    newBlock->periodicity().toggle(1, originalBlock.periodicity().get(1));

    transferBlockLatticeNonLocal(originalBlock, *newBlock, originalBlock.getBoundingBox());

    return newBlock;
}

/* *************** 4. MultiParticleField ************************************** */

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiParticleField2D<DenseParticleField2D<T, Descriptor> > >
    generateMultiDenseParticleField(Box2D boundingBox, plint envelopeWidth)
{
    return std::unique_ptr<MultiParticleField2D<DenseParticleField2D<T, Descriptor> > >(
        new MultiParticleField2D<DenseParticleField2D<T, Descriptor> >(
            defaultMultiBlockPolicy2D().getMultiBlockManagement(boundingBox, envelopeWidth),
            defaultMultiBlockPolicy2D().getCombinedStatistics()));
}

template <class ParticleFieldT>
std::unique_ptr<MultiParticleField2D<ParticleFieldT> > generateMultiParticleField2D(
    Box2D boundingBox, plint envelopeWidth)
{
    return std::unique_ptr<MultiParticleField2D<ParticleFieldT> >(
        new MultiParticleField2D<ParticleFieldT>(
            defaultMultiBlockPolicy2D().getMultiBlockManagement(boundingBox, envelopeWidth),
            defaultMultiBlockPolicy2D().getCombinedStatistics()));
}

template <class ParticleFieldT>
std::unique_ptr<MultiParticleField2D<ParticleFieldT> > generateMultiParticleField2D(
    MultiBlock2D &multiBlock, plint envelopeWidth)
{
    MultiBlockManagement2D sparseBlockManagement(multiBlock.getMultiBlockManagement());
    MultiParticleField2D<ParticleFieldT> *field = new MultiParticleField2D<ParticleFieldT>(
        MultiBlockManagement2D(
            sparseBlockManagement.getSparseBlockStructure(),
            sparseBlockManagement.getThreadAttribution().clone(), envelopeWidth,
            sparseBlockManagement.getRefinementLevel()),
        defaultMultiBlockPolicy2D().getCombinedStatistics());

    field->periodicity().toggle(0, multiBlock.periodicity().get(0));
    field->periodicity().toggle(1, multiBlock.periodicity().get(1));

    return std::unique_ptr<MultiParticleField2D<ParticleFieldT> >(field);
}

}  // namespace plb

#endif  // MULTI_BLOCK_GENERATOR_2D_HH
