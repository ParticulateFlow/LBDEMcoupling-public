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
 * Scalar, vector and tensor fields for 2D data fields -- generic implementation.
 */

#ifndef MULTI_DATA_FIELD_2D_HH
#define MULTI_DATA_FIELD_2D_HH

#include <algorithm>
#include <limits>
#include <sstream>
#include <vector>

#include "core/multiBlockIdentifiers2D.h"
#include "core/plbTypenames.h"
#include "multiBlock/defaultMultiBlockPolicy2D.h"
#include "multiBlock/multiBlockGenerator2D.h"
#include "multiBlock/multiBlockManagement2D.h"
#include "multiBlock/multiDataField2D.h"
#include "multiBlock/nonLocalTransfer2D.h"

namespace plb {

/////// Class MultiScalarField2D //////////////////////////////////

template <typename T>
const int MultiScalarField2D<T>::staticId = meta::registerMultiBlock2D(
    MultiScalarField2D<T>::basicType(), "NA", MultiScalarField2D<T>::blockName(),
    defaultGenerateMultiScalarField2D<T>);

template <typename T>
MultiScalarField2D<T>::MultiScalarField2D(
    MultiBlockManagement2D const &multiBlockManagement_, BlockCommunicator2D *blockCommunicator_,
    CombinedStatistics *combinedStatistics_, MultiScalarAccess2D<T> *multiScalarAccess_, T iniVal) :
    MultiBlock2D(multiBlockManagement_, blockCommunicator_, combinedStatistics_),
    multiScalarAccess(multiScalarAccess_)
{
    allocateFields(iniVal);
}

template <typename T>
MultiScalarField2D<T>::MultiScalarField2D(plint nx, plint ny, T iniVal)
    // Envelope-width defaults to 1.
    :
    MultiBlock2D(nx, ny, 1),
    multiScalarAccess(defaultMultiBlockPolicy2D().getMultiScalarAccess<T>())
{
    allocateFields(iniVal);
}

template <typename T>
MultiScalarField2D<T>::~MultiScalarField2D()
{
    deAllocateFields();
    delete multiScalarAccess;
}

template <typename T>
MultiScalarField2D<T>::MultiScalarField2D(MultiScalarField2D<T> const &rhs) :
    ScalarFieldBase2D<T>(rhs), MultiBlock2D(rhs), multiScalarAccess(rhs.multiScalarAccess->clone())
{
    allocateFields();
    typename BlockMap::iterator it = fields.begin();
    typename BlockMap::const_iterator rhsIt = rhs.fields.begin();

    for (; it != fields.end(); ++it, ++rhsIt) {
        *(it->second) = *(rhsIt->second);
    }
}

template <typename T>
MultiScalarField2D<T>::MultiScalarField2D(MultiBlock2D const &rhs)
    // Use MultiBlock's sub-domain constructor to avoid that the data-processors are copied
    :
    MultiBlock2D(rhs, rhs.getBoundingBox(), false),
    multiScalarAccess(defaultMultiBlockPolicy2D().getMultiScalarAccess<T>())
{
    allocateFields();
}

template <typename T>
MultiScalarField2D<T>::MultiScalarField2D(MultiBlock2D const &rhs, Box2D subDomain, bool crop) :
    MultiBlock2D(rhs, subDomain, crop),
    multiScalarAccess(defaultMultiBlockPolicy2D().getMultiScalarAccess<T>())
{
    allocateFields();
}

template <typename T>
MultiScalarField2D<T> &MultiScalarField2D<T>::operator=(MultiScalarField2D<T> const &rhs)
{
    MultiScalarField2D<T> tmp(rhs);
    swap(tmp);
    return *this;
}

template <typename T>
MultiScalarField2D<T> *MultiScalarField2D<T>::clone() const
{
    return new MultiScalarField2D<T>(*this);
}

template <typename T>
MultiScalarField2D<T> *MultiScalarField2D<T>::clone(
    MultiBlockManagement2D const &newManagement) const
{
    MultiScalarField2D<T> *newField = new MultiScalarField2D<T>(
        newManagement, this->getBlockCommunicator().clone(), this->getCombinedStatistics().clone(),
        multiScalarAccess->clone(), T());
    // Use the same domain in the "from" and "to" argument, so that the data is not shifted
    // in space during the creation of the new block.
    copy(*this, newField->getBoundingBox(), *newField, newField->getBoundingBox());
    return newField;
}

template <typename T>
void MultiScalarField2D<T>::swap(MultiScalarField2D<T> &rhs)
{
    MultiBlock2D::swap(rhs);
    fields.swap(rhs.fields);
    std::swap(multiScalarAccess, rhs.multiScalarAccess);
}

template <typename T>
void MultiScalarField2D<T>::reset()
{
    for (typename BlockMap::iterator it = fields.begin(); it != fields.end(); ++it) {
        it->second->reset();
    }
}

template <typename T>
void MultiScalarField2D<T>::allocateFields(T iniVal)
{
    for (pluint iBlock = 0; iBlock < this->getLocalInfo().getBlocks().size(); ++iBlock) {
        plint blockId = this->getLocalInfo().getBlocks()[iBlock];
        SmartBulk2D bulk(this->getMultiBlockManagement(), blockId);
        Box2D envelope = bulk.computeEnvelope();
        ScalarField2D<T> *newField =
            new ScalarField2D<T>(envelope.getNx(), envelope.getNy(), iniVal);
        newField->setLocation(Dot2D(envelope.x0, envelope.y0));
        fields[blockId] = newField;
    }
}

template <typename T>
void MultiScalarField2D<T>::deAllocateFields()
{
    for (typename BlockMap::iterator it = fields.begin(); it != fields.end(); ++it) {
        delete it->second;
    }
}

template <typename T>
inline T &MultiScalarField2D<T>::get(plint iX, plint iY)
{
    PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
    PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
    return multiScalarAccess->getDistributedScalar(iX, iY, this->getMultiBlockManagement(), fields);
}

template <typename T>
inline T const &MultiScalarField2D<T>::get(plint iX, plint iY) const
{
    PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
    PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
    return multiScalarAccess->getDistributedScalar(iX, iY, this->getMultiBlockManagement(), fields);
}

template <typename T>
ScalarField2D<T> &MultiScalarField2D<T>::getComponent(plint blockId)
{
    typename BlockMap::iterator it = fields.find(blockId);
    PLB_ASSERT(it != fields.end());
    return *it->second;
}

template <typename T>
ScalarField2D<T> const &MultiScalarField2D<T>::getComponent(plint blockId) const
{
    typename BlockMap::const_iterator it = fields.find(blockId);
    PLB_ASSERT(it != fields.end());
    return *it->second;
}

template <typename T>
plint MultiScalarField2D<T>::sizeOfCell() const
{
    return sizeof(T);
}

template <typename T>
plint MultiScalarField2D<T>::getCellDim() const
{
    return 1;
}

template <typename T>
int MultiScalarField2D<T>::getStaticId() const
{
    return staticId;
}

template <typename T>
void MultiScalarField2D<T>::copyReceive(
    MultiBlock2D const &fromBlock, Box2D const &fromDomain, Box2D const &toDomain,
    modif::ModifT whichData)
{
    MultiScalarField2D<T> const *fromField =
        dynamic_cast<MultiScalarField2D<T> const *>(&fromBlock);
    PLB_ASSERT(fromField);
    copy(*fromField, fromDomain, *this, toDomain);
}

template <typename T>
std::string MultiScalarField2D<T>::getBlockName() const
{
    return blockName();
}

template <typename T>
std::vector<std::string> MultiScalarField2D<T>::getTypeInfo() const
{
    std::vector<std::string> info;
    info.push_back(basicType());
    return info;
}

template <typename T>
std::string MultiScalarField2D<T>::blockName()
{
    return std::string("ScalarField2D");
}

template <typename T>
std::string MultiScalarField2D<T>::basicType()
{
    return NativeType<T>::getName();
}

template <typename T>
MultiScalarField2D<T> &findMultiScalarField2D(id_t id)
{
    MultiBlock2D *multiBlock = multiBlockRegistration2D().find(id);
    if (!multiBlock || multiBlock->getStaticId() != MultiScalarField2D<T>::staticId) {
        throw PlbLogicException("Trying to access a multi scalar field that is not registered.");
    }
    return (MultiScalarField2D<T> &)(*multiBlock);
}

//////// Class MultiTensorField2D //////////////////////////////////

template <typename T, int nDim>
const int MultiTensorField2D<T, nDim>::staticId = meta::registerMultiBlock2D(
    MultiTensorField2D<T, nDim>::basicType(), "NA", MultiTensorField2D<T, nDim>::blockName(),
    defaultGenerateMultiTensorField2D<T, nDim>);

template <typename T, int nDim>
MultiTensorField2D<T, nDim>::MultiTensorField2D(
    MultiBlockManagement2D const &multiBlockManagement_, BlockCommunicator2D *blockCommunicator_,
    CombinedStatistics *combinedStatistics_, MultiTensorAccess2D<T, nDim> *multiTensorAccess_) :
    MultiBlock2D(multiBlockManagement_, blockCommunicator_, combinedStatistics_),
    multiTensorAccess(multiTensorAccess_)
{
    allocateFields();
}

template <typename T, int nDim>
MultiTensorField2D<T, nDim>::MultiTensorField2D(
    MultiBlockManagement2D const &multiBlockManagement_, BlockCommunicator2D *blockCommunicator_,
    CombinedStatistics *combinedStatistics_, MultiTensorAccess2D<T, nDim> *multiTensorAccess_,
    Array<T, nDim> const &iniVal) :
    MultiBlock2D(multiBlockManagement_, blockCommunicator_, combinedStatistics_),
    multiTensorAccess(multiTensorAccess_)
{
    allocateFields(iniVal);
}

template <typename T, int nDim>
MultiTensorField2D<T, nDim>::MultiTensorField2D(plint nx, plint ny)
    // Envelope-width defaults to 1.
    :
    MultiBlock2D(nx, ny, 1),
    multiTensorAccess(defaultMultiBlockPolicy2D().getMultiTensorAccess<T, nDim>())
{
    allocateFields();
}

template <typename T, int nDim>
MultiTensorField2D<T, nDim>::MultiTensorField2D(plint nx, plint ny, Array<T, nDim> const &iniVal)
    // Envelope-width defaults to 1.
    :
    MultiBlock2D(nx, ny, 1),
    multiTensorAccess(defaultMultiBlockPolicy2D().getMultiTensorAccess<T, nDim>())
{
    allocateFields(iniVal);
}

template <typename T, int nDim>
MultiTensorField2D<T, nDim>::~MultiTensorField2D()
{
    deAllocateFields();
    delete multiTensorAccess;
}

template <typename T, int nDim>
MultiTensorField2D<T, nDim>::MultiTensorField2D(MultiTensorField2D<T, nDim> const &rhs) :
    TensorFieldBase2D<T, nDim>(rhs),
    MultiBlock2D(rhs),
    multiTensorAccess(rhs.multiTensorAccess->clone())
{
    allocateFields();
    typename BlockMap::iterator it = fields.begin();
    typename BlockMap::const_iterator rhsIt = rhs.fields.begin();

    for (; it != fields.end(); ++it, ++rhsIt) {
        *(it->second) = *(rhsIt->second);
    }
}

template <typename T, int nDim>
MultiTensorField2D<T, nDim>::MultiTensorField2D(MultiBlock2D const &rhs)
    // Use MultiBlock's sub-domain constructor to avoid that the data-processors are copied
    :
    MultiBlock2D(rhs, rhs.getBoundingBox(), false),
    multiTensorAccess(defaultMultiBlockPolicy2D().getMultiTensorAccess<T, nDim>())
{
    allocateFields();
}

template <typename T, int nDim>
MultiTensorField2D<T, nDim>::MultiTensorField2D(
    MultiBlock2D const &rhs, Box2D subDomain, bool crop) :
    MultiBlock2D(rhs, subDomain, crop),
    multiTensorAccess(defaultMultiBlockPolicy2D().getMultiTensorAccess<T, nDim>())
{
    allocateFields();
}

template <typename T, int nDim>
MultiTensorField2D<T, nDim> &MultiTensorField2D<T, nDim>::operator=(
    MultiTensorField2D<T, nDim> const &rhs)
{
    MultiTensorField2D<T, nDim> tmp(rhs);
    swap(tmp);
    return *this;
}

template <typename T, int nDim>
MultiTensorField2D<T, nDim> *MultiTensorField2D<T, nDim>::clone() const
{
    return new MultiTensorField2D<T, nDim>(*this);
}

template <typename T, int nDim>
MultiTensorField2D<T, nDim> *MultiTensorField2D<T, nDim>::clone(
    MultiBlockManagement2D const &newManagement) const
{
    Array<T, nDim> iniVal;
    iniVal.resetToZero();
    MultiTensorField2D<T, nDim> *newField = new MultiTensorField2D<T, nDim>(
        newManagement, this->getBlockCommunicator().clone(), this->getCombinedStatistics().clone(),
        multiTensorAccess->clone(), iniVal);
    // Use the same domain in the "from" and "to" argument, so that the data is not shifted
    // in space during the creation of the new block.
    copy(*this, newField->getBoundingBox(), *newField, newField->getBoundingBox());
    return newField;
}

template <typename T, int nDim>
void MultiTensorField2D<T, nDim>::swap(MultiTensorField2D<T, nDim> &rhs)
{
    MultiBlock2D::swap(rhs);
    fields.swap(rhs.fields);
    std::swap(multiTensorAccess, rhs.multiTensorAccess);
}

template <typename T, int nDim>
void MultiTensorField2D<T, nDim>::reset()
{
    for (typename BlockMap::iterator it = fields.begin(); it != fields.end(); ++it) {
        it->second->reset();
    }
}

template <typename T, int nDim>
void MultiTensorField2D<T, nDim>::allocateFields()
{
    Array<T, nDim> iniVal;
    iniVal.resetToZero();
    allocateFields(iniVal);
}

template <typename T, int nDim>
void MultiTensorField2D<T, nDim>::allocateFields(Array<T, nDim> const &iniVal)
{
    for (pluint iBlock = 0; iBlock < this->getLocalInfo().getBlocks().size(); ++iBlock) {
        plint blockId = this->getLocalInfo().getBlocks()[iBlock];
        SmartBulk2D bulk(this->getMultiBlockManagement(), blockId);
        Box2D envelope = bulk.computeEnvelope();
        TensorField2D<T, nDim> *newField =
            new TensorField2D<T, nDim>(envelope.getNx(), envelope.getNy(), iniVal);
        newField->setLocation(Dot2D(envelope.x0, envelope.y0));
        fields[blockId] = newField;
    }
}

template <typename T, int nDim>
void MultiTensorField2D<T, nDim>::deAllocateFields()
{
    for (typename BlockMap::iterator it = fields.begin(); it != fields.end(); ++it) {
        delete it->second;
    }
}

template <typename T, int nDim>
inline Array<T, nDim> &MultiTensorField2D<T, nDim>::get(plint iX, plint iY)
{
    PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
    PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
    return multiTensorAccess->getDistributedTensor(iX, iY, this->getMultiBlockManagement(), fields);
}

template <typename T, int nDim>
inline Array<T, nDim> const &MultiTensorField2D<T, nDim>::get(plint iX, plint iY) const
{
    PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
    PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
    return multiTensorAccess->getDistributedTensor(iX, iY, this->getMultiBlockManagement(), fields);
}

template <typename T, int nDim>
TensorField2D<T, nDim> &MultiTensorField2D<T, nDim>::getComponent(plint blockId)
{
    typename BlockMap::iterator it = fields.find(blockId);
    PLB_ASSERT(it != fields.end());
    return *it->second;
}

template <typename T, int nDim>
TensorField2D<T, nDim> const &MultiTensorField2D<T, nDim>::getComponent(plint blockId) const
{
    typename BlockMap::const_iterator it = fields.find(blockId);
    PLB_ASSERT(it != fields.end());
    return *it->second;
}

template <typename T, int nDim>
plint MultiTensorField2D<T, nDim>::sizeOfCell() const
{
    return nDim * sizeof(T);
}

template <typename T, int nDim>
plint MultiTensorField2D<T, nDim>::getCellDim() const
{
    return nDim;
}

template <typename T, int nDim>
int MultiTensorField2D<T, nDim>::getStaticId() const
{
    return staticId;
}

template <typename T, int nDim>
void MultiTensorField2D<T, nDim>::copyReceive(
    MultiBlock2D const &fromBlock, Box2D const &fromDomain, Box2D const &toDomain,
    modif::ModifT whichData)
{
    MultiTensorField2D<T, nDim> const *fromField =
        dynamic_cast<MultiTensorField2D<T, nDim> const *>(&fromBlock);
    PLB_ASSERT(fromField);
    copy(*fromField, fromDomain, *this, toDomain);
}

template <typename T, int nDim>
std::string MultiTensorField2D<T, nDim>::getBlockName() const
{
    return blockName();
}

template <typename T, int nDim>
std::vector<std::string> MultiTensorField2D<T, nDim>::getTypeInfo() const
{
    std::vector<std::string> info;
    info.push_back(basicType());
    return info;
}

template <typename T, int nDim>
std::string MultiTensorField2D<T, nDim>::blockName()
{
    std::stringstream ss;
    ss << nDim;
    return std::string("TensorField2D_" + ss.str());
}

template <typename T, int nDim>
std::string MultiTensorField2D<T, nDim>::basicType()
{
    return NativeType<T>::getName();
}

template <typename T, int nDim>
MultiTensorField2D<T, nDim> &findMultiTensorField2D(id_t id)
{
    MultiBlock2D *multiBlock = multiBlockRegistration2D().find(id);
    if (!multiBlock || multiBlock->getStaticId() != MultiTensorField2D<T, nDim>::staticId) {
        throw PlbLogicException("Trying to access a multi scalar field that is not registered.");
    }
    return (MultiTensorField2D<T, nDim> &)(*multiBlock);
}

//////// Class MultiNTensorField2D //////////////////////////////////

template <typename T>
const int MultiNTensorField2D<T>::staticId = meta::registerMultiBlock2D(
    MultiNTensorField2D<T>::basicType(), "NA", MultiNTensorField2D<T>::blockName(),
    defaultGenerateMultiNTensorField2D<T>);

template <typename T>
MultiNTensorField2D<T>::MultiNTensorField2D(
    plint ndim, MultiBlockManagement2D const &multiBlockManagement_,
    BlockCommunicator2D *blockCommunicator_, CombinedStatistics *combinedStatistics_,
    MultiNTensorAccess2D<T> *multiNTensorAccess_) :
    NTensorFieldBase2D<T>(ndim),
    MultiBlock2D(multiBlockManagement_, blockCommunicator_, combinedStatistics_),
    multiNTensorAccess(multiNTensorAccess_)
{
    allocateFields();
}

template <typename T>
MultiNTensorField2D<T>::MultiNTensorField2D(
    plint ndim, T const *iniVal, MultiBlockManagement2D const &multiBlockManagement_,
    BlockCommunicator2D *blockCommunicator_, CombinedStatistics *combinedStatistics_,
    MultiNTensorAccess2D<T> *multiNTensorAccess_) :
    NTensorFieldBase2D<T>(ndim),
    MultiBlock2D(multiBlockManagement_, blockCommunicator_, combinedStatistics_),
    multiNTensorAccess(multiNTensorAccess_)
{
    allocateFields(iniVal);
}

template <typename T>
MultiNTensorField2D<T>::MultiNTensorField2D(plint nx, plint ny, plint ndim) :
    NTensorFieldBase2D<T>(ndim),
    // Envelope-width defaults to 1.
    MultiBlock2D(nx, ny, 1),
    multiNTensorAccess(defaultMultiBlockPolicy2D().getMultiNTensorAccess<T>())
{
    allocateFields();
}

template <typename T>
MultiNTensorField2D<T>::MultiNTensorField2D(plint nx, plint ny, plint ndim, T const *iniVal) :
    NTensorFieldBase2D<T>(ndim),
    // Envelope-width defaults to 1.
    MultiBlock2D(nx, ny, 1),
    multiNTensorAccess(defaultMultiBlockPolicy2D().getMultiNTensorAccess<T>())
{
    allocateFields(iniVal);
}

template <typename T>
MultiNTensorField2D<T>::~MultiNTensorField2D()
{
    deAllocateFields();
    delete multiNTensorAccess;
}

template <typename T>
MultiNTensorField2D<T>::MultiNTensorField2D(MultiNTensorField2D<T> const &rhs) :
    NTensorFieldBase2D<T>(rhs),
    MultiBlock2D(rhs),
    multiNTensorAccess(rhs.multiNTensorAccess->clone())
{
    allocateFields();
    typename BlockMap::iterator it = fields.begin();
    typename BlockMap::const_iterator rhsIt = rhs.fields.begin();

    for (; it != fields.end(); ++it, ++rhsIt) {
        *(it->second) = *(rhsIt->second);
    }
}

template <typename T>
MultiNTensorField2D<T>::MultiNTensorField2D(plint ndim, MultiBlock2D const &rhs) :
    NTensorFieldBase2D<T>(ndim),
    MultiBlock2D(rhs),
    multiNTensorAccess(defaultMultiBlockPolicy2D().getMultiNTensorAccess<T>())
{
    allocateFields();
}

template <typename T>
MultiNTensorField2D<T>::MultiNTensorField2D(
    plint ndim, MultiBlock2D const &rhs, Box2D subDomain, bool crop) :
    NTensorFieldBase2D<T>(ndim),
    MultiBlock2D(rhs, subDomain, crop),
    multiNTensorAccess(defaultMultiBlockPolicy2D().getMultiNTensorAccess<T>())
{
    allocateFields();
}

template <typename T>
MultiNTensorField2D<T> &MultiNTensorField2D<T>::operator=(MultiNTensorField2D<T> const &rhs)
{
    MultiNTensorField2D<T> tmp(rhs);
    swap(tmp);
    return *this;
}

template <typename T>
MultiNTensorField2D<T> *MultiNTensorField2D<T>::clone() const
{
    return new MultiNTensorField2D<T>(*this);
}

template <typename T>
MultiNTensorField2D<T> *MultiNTensorField2D<T>::clone(
    MultiBlockManagement2D const &newManagement) const
{
    MultiNTensorField2D<T> *newField = new MultiNTensorField2D<T>(
        this->getNdim(), newManagement, this->getBlockCommunicator().clone(),
        this->getCombinedStatistics().clone(), multiNTensorAccess->clone());
    // Use the same domain in the "from" and "to" argument, so that the data is not shifted
    // in space during the creation of the new block.
    copy(*this, newField->getBoundingBox(), *newField, newField->getBoundingBox());
    return newField;
}

template <typename T>
void MultiNTensorField2D<T>::swap(MultiNTensorField2D<T> &rhs)
{
    NTensorFieldBase2D<T>::swap(rhs);
    MultiBlock2D::swap(rhs);
    fields.swap(rhs.fields);
    std::swap(multiNTensorAccess, rhs.multiNTensorAccess);
}

template <typename T>
void MultiNTensorField2D<T>::reset()
{
    for (typename BlockMap::iterator it = fields.begin(); it != fields.end(); ++it) {
        it->second->reset();
    }
}

template <typename T>
void MultiNTensorField2D<T>::allocateFields()
{
    std::vector<T> iniVal(this->getNdim());
    std::fill(iniVal.begin(), iniVal.end(), T());
    allocateFields(&iniVal[0]);
}

template <typename T>
void MultiNTensorField2D<T>::allocateFields(T const *iniVal)
{
    for (pluint iBlock = 0; iBlock < this->getLocalInfo().getBlocks().size(); ++iBlock) {
        plint blockId = this->getLocalInfo().getBlocks()[iBlock];
        SmartBulk2D bulk(this->getMultiBlockManagement(), blockId);
        Box2D envelope = bulk.computeEnvelope();
        NTensorField2D<T> *newField =
            new NTensorField2D<T>(envelope.getNx(), envelope.getNy(), this->getNdim(), iniVal);
        newField->setLocation(Dot2D(envelope.x0, envelope.y0));
        fields[blockId] = newField;
    }
}

template <typename T>
void MultiNTensorField2D<T>::deAllocateFields()
{
    for (typename BlockMap::iterator it = fields.begin(); it != fields.end(); ++it) {
        delete it->second;
    }
}

template <typename T>
inline T *MultiNTensorField2D<T>::get(plint iX, plint iY)
{
    PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
    PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
    return multiNTensorAccess->getDistributedNTensor(
        iX, iY, this->getMultiBlockManagement(), fields);
}

template <typename T>
inline T const *MultiNTensorField2D<T>::get(plint iX, plint iY) const
{
    PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
    PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
    return multiNTensorAccess->getDistributedNTensor(
        iX, iY, this->getMultiBlockManagement(), fields);
}

template <typename T>
NTensorField2D<T> &MultiNTensorField2D<T>::getComponent(plint blockId)
{
    typename BlockMap::iterator it = fields.find(blockId);
    PLB_ASSERT(it != fields.end());
    return *it->second;
}

template <typename T>
NTensorField2D<T> const &MultiNTensorField2D<T>::getComponent(plint blockId) const
{
    typename BlockMap::const_iterator it = fields.find(blockId);
    PLB_ASSERT(it != fields.end());
    return *it->second;
}

template <typename T>
plint MultiNTensorField2D<T>::sizeOfCell() const
{
    return this->getNdim() * sizeof(T);
}

template <typename T>
plint MultiNTensorField2D<T>::getCellDim() const
{
    return this->getNdim();
}

template <typename T>
int MultiNTensorField2D<T>::getStaticId() const
{
    return staticId;
}

template <typename T>
void MultiNTensorField2D<T>::copyReceive(
    MultiBlock2D const &fromBlock, Box2D const &fromDomain, Box2D const &toDomain,
    modif::ModifT whichData)
{
    MultiNTensorField2D<T> const *fromField =
        dynamic_cast<MultiNTensorField2D<T> const *>(&fromBlock);
    PLB_ASSERT(fromField);
    copy(*fromField, fromDomain, *this, toDomain);
}

template <typename T>
std::string MultiNTensorField2D<T>::getBlockName() const
{
    return blockName();
}

template <typename T>
std::vector<std::string> MultiNTensorField2D<T>::getTypeInfo() const
{
    std::vector<std::string> info;
    info.push_back(basicType());
    return info;
}

template <typename T>
std::string MultiNTensorField2D<T>::blockName()
{
    return std::string("NTensorField2D");
}

template <typename T>
std::string MultiNTensorField2D<T>::basicType()
{
    return NativeType<T>::getName();
}

template <typename T>
MultiNTensorField2D<T> &findMultiNTensorField2D(id_t id)
{
    MultiBlock2D *multiBlock = multiBlockRegistration2D().find(id);
    if (!multiBlock || multiBlock->getStaticId() != MultiNTensorField2D<T>::staticId) {
        throw PlbLogicException("Trying to access a multi scalar field that is not registered.");
    }
    return (MultiNTensorField2D<T> &)(*multiBlock);
}

}  // namespace plb

#endif  // MULTI_DATA_FIELD_2D_HH
