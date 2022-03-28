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
 * Scalar, vector and tensor fields for 3D data fields -- generic implementation.
 */

#ifndef MULTI_DATA_FIELD_3D_HH
#define MULTI_DATA_FIELD_3D_HH

#include <algorithm>
#include <limits>
#include <sstream>
#include <vector>

#include "atomicBlock/dataField3D.h"
#include "atomicBlock/dataField3D.hh"
#include "core/multiBlockIdentifiers3D.h"
#include "core/plbTypenames.h"
#include "multiBlock/defaultMultiBlockPolicy3D.h"
#include "multiBlock/multiBlockGenerator3D.h"
#include "multiBlock/multiBlockManagement3D.h"
#include "multiBlock/multiDataField3D.h"
#include "multiBlock/nonLocalTransfer3D.h"

namespace plb {

/////// Class MultiScalarField3D //////////////////////////////////

template <typename T>
const int MultiScalarField3D<T>::staticId = meta::registerMultiBlock3D(
    MultiScalarField3D<T>::basicType(), "NA", MultiScalarField3D<T>::blockName(),
    defaultGenerateMultiScalarField3D<T>);

template <typename T>
MultiScalarField3D<T>::MultiScalarField3D(
    MultiBlockManagement3D const &multiBlockManagement_, BlockCommunicator3D *blockCommunicator_,
    CombinedStatistics *combinedStatistics_, MultiScalarAccess3D<T> *multiScalarAccess_, T iniVal) :
    MultiBlock3D(multiBlockManagement_, blockCommunicator_, combinedStatistics_),
    multiScalarAccess(multiScalarAccess_),
    nTensorViewBlock(0)
{
    allocateFields(iniVal);
}

template <typename T>
MultiScalarField3D<T>::MultiScalarField3D(plint nx, plint ny, plint nz, T iniVal)
    // Envelope-width defaults to 1.
    :
    MultiBlock3D(nx, ny, nz, 1),
    multiScalarAccess(defaultMultiBlockPolicy3D().getMultiScalarAccess<T>()),
    nTensorViewBlock(0)
{
    allocateFields(iniVal);
}

template <typename T>
MultiScalarField3D<T>::~MultiScalarField3D()
{
    if (nTensorViewBlock) {
        delete nTensorViewBlock;
    }
    deAllocateFields();
    delete multiScalarAccess;
}

template <typename T>
MultiScalarField3D<T>::MultiScalarField3D(MultiScalarField3D<T> const &rhs) :
    ScalarFieldBase3D<T>(rhs),
    MultiBlock3D(rhs),
    multiScalarAccess(rhs.multiScalarAccess->clone()),
    nTensorViewBlock(0)
{
    allocateFields();
    typename BlockMap::iterator it = fields.begin();
    typename BlockMap::const_iterator rhsIt = rhs.fields.begin();

    for (; it != fields.end(); ++it, ++rhsIt) {
        PLB_ASSERT(rhsIt != rhs.fields.end());
        *(it->second) = *(rhsIt->second);
    }
}

template <typename T>
MultiScalarField3D<T>::MultiScalarField3D(MultiBlock3D const &rhs)
    // Use MultiBlock's sub-domain constructor to avoid that the data-processors are copied
    :
    MultiBlock3D(rhs, rhs.getBoundingBox(), false),
    multiScalarAccess(defaultMultiBlockPolicy3D().getMultiScalarAccess<T>()),
    nTensorViewBlock(0)
{
    allocateFields();
}

template <typename T>
MultiScalarField3D<T>::MultiScalarField3D(MultiNTensorField3D<T> &rhs, bool shareMemory)
    // Use MultiBlock's sub-domain constructor to avoid that the data-processors are copied
    :
    MultiBlock3D(rhs, rhs.getBoundingBox(), false),
    multiScalarAccess(defaultMultiBlockPolicy3D().getMultiScalarAccess<T>()),
    nTensorViewBlock(0)
{
    PLB_ASSERT(rhs.getNdim() == 1);
    if (shareMemory) {
        for (typename MultiNTensorField3D<T>::BlockMap::iterator it = rhs.fields.begin();
             it != rhs.fields.end(); ++it)
        {
            fields[it->first] = new ScalarField3D<T>(*it->second);
        }
    } else {
        allocateFields();
    }
}

template <typename T>
MultiScalarField3D<T>::MultiScalarField3D(MultiBlock3D const &rhs, Box3D subDomain, bool crop) :
    MultiBlock3D(rhs, subDomain, crop),
    multiScalarAccess(defaultMultiBlockPolicy3D().getMultiScalarAccess<T>()),
    nTensorViewBlock(0)
{
    allocateFields();
}

template <typename T>
MultiScalarField3D<T> &MultiScalarField3D<T>::operator=(MultiScalarField3D<T> const &rhs)
{
    MultiScalarField3D<T> tmp(rhs);
    swap(tmp);
    return *this;
}

template <typename T>
MultiScalarField3D<T> *MultiScalarField3D<T>::clone() const
{
    return new MultiScalarField3D<T>(*this);
}

template <typename T>
MultiScalarField3D<T> *MultiScalarField3D<T>::clone(
    MultiBlockManagement3D const &newManagement) const
{
    MultiScalarField3D<T> *newField = new MultiScalarField3D<T>(
        newManagement, this->getBlockCommunicator().clone(), this->getCombinedStatistics().clone(),
        multiScalarAccess->clone(), T());
    // Use the same domain in the "from" and "to" argument, so that the data is not shifted
    // in space during the creation of the new block.
    copy(*this, newField->getBoundingBox(), *newField, newField->getBoundingBox());
    return newField;
}

template <typename T>
void MultiScalarField3D<T>::swap(MultiScalarField3D<T> &rhs)
{
    MultiBlock3D::swap(rhs);
    fields.swap(rhs.fields);
    std::swap(multiScalarAccess, rhs.multiScalarAccess);
    std::swap(nTensorViewBlock, rhs.nTensorViewBlock);
}

template <typename T>
void MultiScalarField3D<T>::reset()
{
    for (typename BlockMap::iterator it = fields.begin(); it != fields.end(); ++it) {
        it->second->reset();
    }
}

template <typename T>
void MultiScalarField3D<T>::allocateFields(T iniVal)
{
    for (pluint iBlock = 0; iBlock < this->getLocalInfo().getBlocks().size(); ++iBlock) {
        plint blockId = this->getLocalInfo().getBlocks()[iBlock];
        SmartBulk3D bulk(this->getMultiBlockManagement(), blockId);
        Box3D envelope = bulk.computeEnvelope();
        ScalarField3D<T> *newField =
            new ScalarField3D<T>(envelope.getNx(), envelope.getNy(), envelope.getNz(), iniVal);
        newField->setLocation(Dot3D(envelope.x0, envelope.y0, envelope.z0));
        fields[blockId] = newField;
    }
}

template <typename T>
void MultiScalarField3D<T>::deAllocateFields()
{
    for (typename BlockMap::iterator it = fields.begin(); it != fields.end(); ++it) {
        delete it->second;
    }
}

template <typename T>
inline T &MultiScalarField3D<T>::get(plint iX, plint iY, plint iZ)
{
    PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
    PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
    PLB_PRECONDITION(iZ >= 0 && iZ < this->getNz());
    return multiScalarAccess->getDistributedScalar(
        iX, iY, iZ, this->getMultiBlockManagement(), fields);
}

template <typename T>
inline T const &MultiScalarField3D<T>::get(plint iX, plint iY, plint iZ) const
{
    PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
    PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
    PLB_PRECONDITION(iZ >= 0 && iZ < this->getNz());
    return multiScalarAccess->getDistributedScalar(
        iX, iY, iZ, this->getMultiBlockManagement(), fields);
}

template <typename T>
ScalarField3D<T> &MultiScalarField3D<T>::getComponent(plint blockId)
{
    typename BlockMap::iterator it = fields.find(blockId);
    PLB_ASSERT(it != fields.end());
    return *it->second;
}

template <typename T>
ScalarField3D<T> const &MultiScalarField3D<T>::getComponent(plint blockId) const
{
    typename BlockMap::const_iterator it = fields.find(blockId);
    PLB_ASSERT(it != fields.end());
    return *it->second;
}

template <typename T>
plint MultiScalarField3D<T>::sizeOfCell() const
{
    return sizeof(T);
}

template <typename T>
plint MultiScalarField3D<T>::getCellDim() const
{
    return 1;
}

template <typename T>
int MultiScalarField3D<T>::getStaticId() const
{
    return staticId;
}

// TODO: whichData is unused, why?
template <typename T>
void MultiScalarField3D<T>::copyReceive(
    MultiBlock3D const &fromBlock, Box3D const &fromDomain, Box3D const &toDomain,
    modif::ModifT whichData)
{
    MultiScalarField3D<T> const *fromField =
        dynamic_cast<MultiScalarField3D<T> const *>(&fromBlock);
    PLB_ASSERT(fromField);
    copy(*fromField, fromDomain, *this, toDomain);
}

template <typename T>
std::string MultiScalarField3D<T>::getBlockName() const
{
    return blockName();
}

template <typename T>
std::vector<std::string> MultiScalarField3D<T>::getTypeInfo() const
{
    std::vector<std::string> info;
    info.push_back(basicType());
    return info;
}

template <typename T>
std::string MultiScalarField3D<T>::blockName()
{
    return std::string("ScalarField3D");
}

template <typename T>
std::string MultiScalarField3D<T>::basicType()
{
    return NativeType<T>::getName();
}

template <typename T>
MultiNTensorField3D<T> &MultiScalarField3D<T>::nTensorView()
{
    if (!nTensorViewBlock) {
        bool shareMemory = true;
        nTensorViewBlock = new MultiNTensorField3D<T>(*this, shareMemory);
    }
    return *nTensorViewBlock;
}

template <typename T>
MultiScalarField3D<T> &findMultiScalarField3D(id_t id)
{
    MultiBlock3D *multiBlock = multiBlockRegistration3D().find(id);
    if (!multiBlock || multiBlock->getStaticId() != MultiScalarField3D<T>::staticId) {
        throw PlbLogicException("Trying to access a multi scalar field that is not registered.");
    }
    return (MultiScalarField3D<T> &)(*multiBlock);
}

//////// Class MultiTensorField3D //////////////////////////////////

template <typename T, int nDim>
const int MultiTensorField3D<T, nDim>::staticId = meta::registerMultiBlock3D(
    MultiTensorField3D<T, nDim>::basicType(), "NA", MultiTensorField3D<T, nDim>::blockName(),
    defaultGenerateMultiTensorField3D<T, nDim>);

template <typename T, int nDim>
MultiTensorField3D<T, nDim>::MultiTensorField3D(
    MultiBlockManagement3D const &multiBlockManagement_, BlockCommunicator3D *blockCommunicator_,
    CombinedStatistics *combinedStatistics_, MultiTensorAccess3D<T, nDim> *multiTensorAccess_) :
    MultiBlock3D(multiBlockManagement_, blockCommunicator_, combinedStatistics_),
    multiTensorAccess(multiTensorAccess_),
    nTensorViewBlock(0)
{
    allocateFields();
}

template <typename T, int nDim>
MultiTensorField3D<T, nDim>::MultiTensorField3D(
    MultiBlockManagement3D const &multiBlockManagement_, BlockCommunicator3D *blockCommunicator_,
    CombinedStatistics *combinedStatistics_, MultiTensorAccess3D<T, nDim> *multiTensorAccess_,
    Array<T, nDim> const &iniVal) :
    MultiBlock3D(multiBlockManagement_, blockCommunicator_, combinedStatistics_),
    multiTensorAccess(multiTensorAccess_),
    nTensorViewBlock(0)
{
    allocateFields(iniVal);
}

template <typename T, int nDim>
MultiTensorField3D<T, nDim>::MultiTensorField3D(plint nx, plint ny, plint nz)
    // Envelope-width defaults to 1.
    :
    MultiBlock3D(nx, ny, nz, 1),
    multiTensorAccess(defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>()),
    nTensorViewBlock(0)
{
    allocateFields();
}

template <typename T, int nDim>
MultiTensorField3D<T, nDim>::MultiTensorField3D(
    plint nx, plint ny, plint nz, Array<T, nDim> const &iniVal)
    // Envelope-width defaults to 1.
    :
    MultiBlock3D(nx, ny, nz, 1),
    multiTensorAccess(defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>()),
    nTensorViewBlock(0)
{
    allocateFields(iniVal);
}

template <typename T, int nDim>
MultiTensorField3D<T, nDim>::~MultiTensorField3D()
{
    if (nTensorViewBlock) {
        delete nTensorViewBlock;
    }
    deAllocateFields();
    delete multiTensorAccess;
}

template <typename T, int nDim>
MultiTensorField3D<T, nDim>::MultiTensorField3D(MultiTensorField3D<T, nDim> const &rhs) :
    TensorFieldBase3D<T, nDim>(rhs),
    MultiBlock3D(rhs),
    multiTensorAccess(rhs.multiTensorAccess->clone()),
    nTensorViewBlock(0)
{
    allocateFields();
    typename BlockMap::iterator it = fields.begin();
    typename BlockMap::const_iterator rhsIt = rhs.fields.begin();

    for (; it != fields.end(); ++it, ++rhsIt) {
        PLB_ASSERT(rhsIt != rhs.fields.end());
        *(it->second) = *(rhsIt->second);
    }
}

template <typename T, int nDim>
MultiTensorField3D<T, nDim>::MultiTensorField3D(MultiBlock3D const &rhs)
    // Use MultiBlock's sub-domain constructor to avoid that the data-processors are copied
    :
    MultiBlock3D(rhs, rhs.getBoundingBox(), false),
    multiTensorAccess(defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>()),
    nTensorViewBlock(0)
{
    allocateFields();
}

template <typename T, int nDim>
MultiTensorField3D<T, nDim>::MultiTensorField3D(MultiNTensorField3D<T> &rhs, bool shareMemory)
    // Use MultiBlock's sub-domain constructor to avoid that the data-processors are copied
    :
    MultiBlock3D(rhs, rhs.getBoundingBox(), false),
    multiTensorAccess(defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>()),
    nTensorViewBlock(0)
{
    PLB_ASSERT(rhs.getNdim() == nDim);
    if (shareMemory) {
        for (typename MultiNTensorField3D<T>::BlockMap::iterator it = rhs.fields.begin();
             it != rhs.fields.end(); ++it)
        {
            fields[it->first] = new TensorField3D<T, nDim>(*it->second);
        }
    } else {
        allocateFields();
    }
}

template <typename T, int nDim>
MultiTensorField3D<T, nDim>::MultiTensorField3D(
    MultiBlock3D const &rhs, Box3D subDomain, bool crop) :
    MultiBlock3D(rhs, subDomain, crop),
    multiTensorAccess(defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>()),
    nTensorViewBlock(0)
{
    allocateFields();
}

template <typename T, int nDim>
MultiTensorField3D<T, nDim> &MultiTensorField3D<T, nDim>::operator=(
    MultiTensorField3D<T, nDim> const &rhs)
{
    MultiTensorField3D<T, nDim> tmp(rhs);
    swap(tmp);
    return *this;
}

template <typename T, int nDim>
MultiTensorField3D<T, nDim> *MultiTensorField3D<T, nDim>::clone() const
{
    return new MultiTensorField3D<T, nDim>(*this);
}

template <typename T, int nDim>
MultiTensorField3D<T, nDim> *MultiTensorField3D<T, nDim>::clone(
    MultiBlockManagement3D const &newManagement) const
{
    Array<T, nDim> iniVal;
    iniVal.resetToZero();
    MultiTensorField3D<T, nDim> *newField = new MultiTensorField3D<T, nDim>(
        newManagement, this->getBlockCommunicator().clone(), this->getCombinedStatistics().clone(),
        multiTensorAccess->clone(), iniVal);
    // Use the same domain in the "from" and "to" argument, so that the data is not shifted
    // in space during the creation of the new block.
    copy(*this, newField->getBoundingBox(), *newField, newField->getBoundingBox());
    return newField;
}

template <typename T, int nDim>
void MultiTensorField3D<T, nDim>::swap(MultiTensorField3D<T, nDim> &rhs)
{
    MultiBlock3D::swap(rhs);
    fields.swap(rhs.fields);
    std::swap(multiTensorAccess, rhs.multiTensorAccess);
    std::swap(nTensorViewBlock, rhs.nTensorViewBlock);
}

template <typename T, int nDim>
void MultiTensorField3D<T, nDim>::reset()
{
    for (typename BlockMap::iterator it = fields.begin(); it != fields.end(); ++it) {
        it->second->reset();
    }
}

template <typename T, int nDim>
void MultiTensorField3D<T, nDim>::allocateFields()
{
    Array<T, nDim> iniVal;
    iniVal.resetToZero();
    allocateFields(iniVal);
}

template <typename T, int nDim>
void MultiTensorField3D<T, nDim>::allocateFields(Array<T, nDim> const &iniVal)
{
    for (pluint iBlock = 0; iBlock < this->getLocalInfo().getBlocks().size(); ++iBlock) {
        plint blockId = this->getLocalInfo().getBlocks()[iBlock];
        SmartBulk3D bulk(this->getMultiBlockManagement(), blockId);
        Box3D envelope = bulk.computeEnvelope();
        TensorField3D<T, nDim> *newField = new TensorField3D<T, nDim>(
            envelope.getNx(), envelope.getNy(), envelope.getNz(), iniVal);
        newField->setLocation(Dot3D(envelope.x0, envelope.y0, envelope.z0));
        fields[blockId] = newField;
    }
}

template <typename T, int nDim>
void MultiTensorField3D<T, nDim>::deAllocateFields()
{
    for (typename BlockMap::iterator it = fields.begin(); it != fields.end(); ++it) {
        delete it->second;
    }
}

template <typename T, int nDim>
inline Array<T, nDim> &MultiTensorField3D<T, nDim>::get(plint iX, plint iY, plint iZ)
{
    PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
    PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
    PLB_PRECONDITION(iZ >= 0 && iZ < this->getNz());
    return multiTensorAccess->getDistributedTensor(
        iX, iY, iZ, this->getMultiBlockManagement(), fields);
}

template <typename T, int nDim>
inline Array<T, nDim> const &MultiTensorField3D<T, nDim>::get(plint iX, plint iY, plint iZ) const
{
    PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
    PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
    PLB_PRECONDITION(iZ >= 0 && iZ < this->getNz());
    return multiTensorAccess->getDistributedTensor(
        iX, iY, iZ, this->getMultiBlockManagement(), fields);
}

template <typename T, int nDim>
TensorField3D<T, nDim> &MultiTensorField3D<T, nDim>::getComponent(plint blockId)
{
    typename BlockMap::iterator it = fields.find(blockId);
    PLB_ASSERT(it != fields.end());
    return *it->second;
}

template <typename T, int nDim>
TensorField3D<T, nDim> const &MultiTensorField3D<T, nDim>::getComponent(plint blockId) const
{
    typename BlockMap::const_iterator it = fields.find(blockId);
    PLB_ASSERT(it != fields.end());
    return *it->second;
}

template <typename T, int nDim>
plint MultiTensorField3D<T, nDim>::sizeOfCell() const
{
    return nDim * sizeof(T);
}

template <typename T, int nDim>
plint MultiTensorField3D<T, nDim>::getCellDim() const
{
    return nDim;
}

template <typename T, int nDim>
int MultiTensorField3D<T, nDim>::getStaticId() const
{
    return staticId;
}
// TODO: whichData is unused, why?
template <typename T, int nDim>
void MultiTensorField3D<T, nDim>::copyReceive(
    MultiBlock3D const &fromBlock, Box3D const &fromDomain, Box3D const &toDomain,
    modif::ModifT whichData)
{
    MultiTensorField3D<T, nDim> const *fromField =
        dynamic_cast<MultiTensorField3D<T, nDim> const *>(&fromBlock);
    PLB_ASSERT(fromField);
    copy(*fromField, fromDomain, *this, toDomain);
}

template <typename T, int nDim>
std::string MultiTensorField3D<T, nDim>::getBlockName() const
{
    return blockName();
}

template <typename T, int nDim>
std::vector<std::string> MultiTensorField3D<T, nDim>::getTypeInfo() const
{
    std::vector<std::string> info;
    info.push_back(basicType());
    return info;
}

template <typename T, int nDim>
std::string MultiTensorField3D<T, nDim>::blockName()
{
    std::stringstream ss;
    ss << nDim;
    return std::string("TensorField3D_" + ss.str());
}

template <typename T, int nDim>
std::string MultiTensorField3D<T, nDim>::basicType()
{
    return NativeType<T>::getName();
}

template <typename T, int nDim>
MultiNTensorField3D<T> &MultiTensorField3D<T, nDim>::nTensorView()
{
    if (!nTensorViewBlock) {
        bool shareMemory = true;
        nTensorViewBlock = new MultiNTensorField3D<T>(*this, shareMemory);
    }
    return *nTensorViewBlock;
}

template <typename T, int nDim>
MultiTensorField3D<T, nDim> &findMultiTensorField3D(id_t id)
{
    MultiBlock3D *multiBlock = multiBlockRegistration3D().find(id);
    if (!multiBlock || multiBlock->getStaticId() != MultiTensorField3D<T, nDim>::staticId) {
        throw PlbLogicException("Trying to access a multi scalar field that is not registered.");
    }
    return (MultiTensorField3D<T, nDim> &)(*multiBlock);
}

//////// Class MultiNTensorField3D //////////////////////////////////

template <typename T>
const int MultiNTensorField3D<T>::staticId = meta::registerMultiBlock3D(
    MultiNTensorField3D<T>::basicType(), "NA", MultiNTensorField3D<T>::blockName(),
    defaultGenerateMultiNTensorField3D<T>);

template <typename T>
MultiNTensorField3D<T>::MultiNTensorField3D(
    plint ndim, MultiBlockManagement3D const &multiBlockManagement_,
    BlockCommunicator3D *blockCommunicator_, CombinedStatistics *combinedStatistics_,
    MultiNTensorAccess3D<T> *multiNTensorAccess_) :
    NTensorFieldBase3D<T>(ndim),
    MultiBlock3D(multiBlockManagement_, blockCommunicator_, combinedStatistics_),
    multiNTensorAccess(multiNTensorAccess_),
    scalarOrTensorView(0)
{
    allocateFields();
}

template <typename T>
MultiNTensorField3D<T>::MultiNTensorField3D(
    plint ndim, T const *iniVal, MultiBlockManagement3D const &multiBlockManagement_,
    BlockCommunicator3D *blockCommunicator_, CombinedStatistics *combinedStatistics_,
    MultiNTensorAccess3D<T> *multiNTensorAccess_) :
    NTensorFieldBase3D<T>(ndim),
    MultiBlock3D(multiBlockManagement_, blockCommunicator_, combinedStatistics_),
    multiNTensorAccess(multiNTensorAccess_),
    scalarOrTensorView(0)
{
    allocateFields(iniVal);
}

template <typename T>
MultiNTensorField3D<T>::MultiNTensorField3D(plint nx, plint ny, plint nz, plint ndim) :
    NTensorFieldBase3D<T>(ndim),
    // Envelope-width defaults to 1.
    MultiBlock3D(nx, ny, nz, 1),
    multiNTensorAccess(defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>()),
    scalarOrTensorView(0)
{
    allocateFields();
}

template <typename T>
MultiNTensorField3D<T>::MultiNTensorField3D(
    plint nx, plint ny, plint nz, plint ndim, T const *iniVal) :
    NTensorFieldBase3D<T>(ndim),
    // Envelope-width defaults to 1.
    MultiBlock3D(nx, ny, nz, 1),
    multiNTensorAccess(defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>()),
    scalarOrTensorView(0)
{
    allocateFields(iniVal);
}

template <typename T>
MultiNTensorField3D<T>::~MultiNTensorField3D()
{
    if (scalarOrTensorView) {
        delete scalarOrTensorView;
    }
    deAllocateFields();
    delete multiNTensorAccess;
}

template <typename T>
MultiNTensorField3D<T>::MultiNTensorField3D(MultiNTensorField3D<T> const &rhs) :
    NTensorFieldBase3D<T>(rhs),
    MultiBlock3D(rhs),
    multiNTensorAccess(rhs.multiNTensorAccess->clone()),
    scalarOrTensorView(0)
{
    allocateFields();
    typename BlockMap::iterator it = fields.begin();
    typename BlockMap::const_iterator rhsIt = rhs.fields.begin();

    for (; it != fields.end(); ++it, ++rhsIt) {
        PLB_ASSERT(rhsIt != rhs.fields.end());
        *(it->second) = *(rhsIt->second);
    }
}

template <typename T>
MultiNTensorField3D<T>::MultiNTensorField3D(MultiScalarField3D<T> &rhs, bool shareMemory) :
    NTensorFieldBase3D<T>(1),
    // Use MultiBlock's sub-domain constructor to avoid that the data-processors are copied
    MultiBlock3D(rhs, rhs.getBoundingBox(), false),
    multiNTensorAccess(defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>()),
    scalarOrTensorView(0)
{
    if (shareMemory) {
        for (typename MultiScalarField3D<T>::BlockMap::iterator it = rhs.fields.begin();
             it != rhs.fields.end(); ++it)
        {
            fields[it->first] = new NTensorField3D<T>(*it->second);
        }
    } else {
        allocateFields();
    }
}

template <typename T>
template <int nDim>
MultiNTensorField3D<T>::MultiNTensorField3D(MultiTensorField3D<T, nDim> &rhs, bool shareMemory) :
    NTensorFieldBase3D<T>(nDim),
    // Use MultiBlock's sub-domain constructor to avoid that the data-processors are copied
    MultiBlock3D(rhs, rhs.getBoundingBox(), false),
    multiNTensorAccess(defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>()),
    scalarOrTensorView(0)
{
    if (shareMemory) {
        for (typename MultiTensorField3D<T, nDim>::BlockMap::iterator it = rhs.fields.begin();
             it != rhs.fields.end(); ++it)
        {
            fields[it->first] = new NTensorField3D<T>(*it->second);
        }
    } else {
        allocateFields();
    }
}

template <typename T>
MultiNTensorField3D<T>::MultiNTensorField3D(plint ndim, MultiBlock3D const &rhs) :
    NTensorFieldBase3D<T>(ndim),
    MultiBlock3D(rhs),
    multiNTensorAccess(defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>()),
    scalarOrTensorView(0)
{
    allocateFields();
}

template <typename T>
MultiNTensorField3D<T>::MultiNTensorField3D(
    plint ndim, MultiBlock3D const &rhs, Box3D subDomain, bool crop) :
    NTensorFieldBase3D<T>(ndim),
    MultiBlock3D(rhs, subDomain, crop),
    multiNTensorAccess(defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>()),
    scalarOrTensorView(0)
{
    allocateFields();
}

template <typename T>
MultiNTensorField3D<T> &MultiNTensorField3D<T>::operator=(MultiNTensorField3D<T> const &rhs)
{
    MultiNTensorField3D<T> tmp(rhs);
    swap(tmp);
    return *this;
}

template <typename T>
MultiNTensorField3D<T> *MultiNTensorField3D<T>::clone() const
{
    return new MultiNTensorField3D<T>(*this);
}

template <typename T>
MultiNTensorField3D<T> *MultiNTensorField3D<T>::clone(
    MultiBlockManagement3D const &newManagement) const
{
    MultiNTensorField3D<T> *newField = new MultiNTensorField3D<T>(
        this->getNdim(), newManagement, this->getBlockCommunicator().clone(),
        this->getCombinedStatistics().clone(), multiNTensorAccess->clone());
    // Use the same domain in the "from" and "to" argument, so that the data is not shifted
    // in space during the creation of the new block.
    copy(*this, newField->getBoundingBox(), *newField, newField->getBoundingBox());
    return newField;
}

template <typename T>
void MultiNTensorField3D<T>::swap(MultiNTensorField3D<T> &rhs)
{
    NTensorFieldBase3D<T>::swap(rhs);
    MultiBlock3D::swap(rhs);
    fields.swap(rhs.fields);
    std::swap(multiNTensorAccess, rhs.multiNTensorAccess);
    std::swap(scalarOrTensorView, rhs.scalarOrTensorView);
}

template <typename T>
void MultiNTensorField3D<T>::reset()
{
    for (typename BlockMap::iterator it = fields.begin(); it != fields.end(); ++it) {
        it->second->reset();
    }
}

template <typename T>
void MultiNTensorField3D<T>::allocateFields()
{
    std::vector<T> iniVal(this->getNdim());
    std::fill(iniVal.begin(), iniVal.end(), T());
    allocateFields(&iniVal[0]);
}

template <typename T>
void MultiNTensorField3D<T>::allocateFields(T const *iniVal)
{
    for (pluint iBlock = 0; iBlock < this->getLocalInfo().getBlocks().size(); ++iBlock) {
        plint blockId = this->getLocalInfo().getBlocks()[iBlock];
        SmartBulk3D bulk(this->getMultiBlockManagement(), blockId);
        Box3D envelope = bulk.computeEnvelope();
        NTensorField3D<T> *newField = new NTensorField3D<T>(
            envelope.getNx(), envelope.getNy(), envelope.getNz(), this->getNdim(), iniVal);
        newField->setLocation(Dot3D(envelope.x0, envelope.y0, envelope.z0));
        fields[blockId] = newField;
    }
}

template <typename T>
void MultiNTensorField3D<T>::deAllocateFields()
{
    for (typename BlockMap::iterator it = fields.begin(); it != fields.end(); ++it) {
        delete it->second;
    }
}

template <typename T>
inline T *MultiNTensorField3D<T>::get(plint iX, plint iY, plint iZ)
{
    PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
    PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
    PLB_PRECONDITION(iZ >= 0 && iZ < this->getNz());
    return multiNTensorAccess->getDistributedNTensor(
        iX, iY, iZ, this->getMultiBlockManagement(), fields);
}

template <typename T>
inline T const *MultiNTensorField3D<T>::get(plint iX, plint iY, plint iZ) const
{
    PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
    PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
    PLB_PRECONDITION(iZ >= 0 && iZ < this->getNz());
    return multiNTensorAccess->getDistributedNTensor(
        iX, iY, iZ, this->getMultiBlockManagement(), fields);
}

template <typename T>
NTensorField3D<T> &MultiNTensorField3D<T>::getComponent(plint blockId)
{
    typename BlockMap::iterator it = fields.find(blockId);
    PLB_ASSERT(it != fields.end());
    return *it->second;
}

template <typename T>
NTensorField3D<T> const &MultiNTensorField3D<T>::getComponent(plint blockId) const
{
    typename BlockMap::const_iterator it = fields.find(blockId);
    PLB_ASSERT(it != fields.end());
    return *it->second;
}

template <typename T>
plint MultiNTensorField3D<T>::sizeOfCell() const
{
    return this->getNdim() * sizeof(T);
}

template <typename T>
plint MultiNTensorField3D<T>::getCellDim() const
{
    return this->getNdim();
}

template <typename T>
int MultiNTensorField3D<T>::getStaticId() const
{
    return staticId;
}

// TODO: whichData and toDomain is unused, why?
template <typename T>
void MultiNTensorField3D<T>::copyReceive(
    MultiBlock3D const &fromBlock, Box3D const &fromDomain, Box3D const &toDomain,
    modif::ModifT whichData)
{
    MultiNTensorField3D<T> const *fromField =
        dynamic_cast<MultiNTensorField3D<T> const *>(&fromBlock);
    PLB_ASSERT(fromField);
    copy(*fromField, fromDomain, *this, toDomain);
}

template <typename T>
std::string MultiNTensorField3D<T>::getBlockName() const
{
    return blockName();
}

template <typename T>
std::vector<std::string> MultiNTensorField3D<T>::getTypeInfo() const
{
    std::vector<std::string> info;
    info.push_back(basicType());
    return info;
}

template <typename T>
std::string MultiNTensorField3D<T>::blockName()
{
    return std::string("NTensorField3D");
}

template <typename T>
std::string MultiNTensorField3D<T>::basicType()
{
    return NativeType<T>::getName();
}

template <typename T>
MultiScalarField3D<T> &MultiNTensorField3D<T>::scalarView()
{
    PLB_ASSERT(this->getNdim() == 1);
    if (!scalarOrTensorView) {
        bool shareMemory = true;
        scalarOrTensorView = new MultiScalarField3D<T>(*this, shareMemory);
    }
    return *dynamic_cast<MultiScalarField3D<T> *>(scalarOrTensorView);
}

template <typename T>
template <int nDim>
MultiTensorField3D<T, nDim> &MultiNTensorField3D<T>::tensorView()
{
    PLB_ASSERT(this->getNdim() > 1 && this->getNdim() == nDim);
    if (!scalarOrTensorView) {
        bool shareMemory = true;
        scalarOrTensorView = new MultiTensorField3D<T, nDim>(*this, shareMemory);
    }
    return *dynamic_cast<MultiTensorField3D<T, nDim> *>(scalarOrTensorView);
}

template <typename T>
MultiNTensorField3D<T> &findMultiNTensorField3D(id_t id)
{
    MultiBlock3D *multiBlock = multiBlockRegistration3D().find(id);
    if (!multiBlock || multiBlock->getStaticId() != MultiNTensorField3D<T>::staticId) {
        throw PlbLogicException("Trying to access a multi scalar field that is not registered.");
    }
    return (MultiNTensorField3D<T> &)(*multiBlock);
}

}  // namespace plb

#endif  // MULTI_DATA_FIELD_3D_HH
