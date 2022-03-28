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
 * Serial access to elements of a scalar/tensor field -- generic implementation.
 */
#ifndef SERIAL_MULTI_DATA_FIELD_3D_HH
#define SERIAL_MULTI_DATA_FIELD_3D_HH

#include "serialMultiDataField3D.h"

/* *************** Class SerialScalarAccess3D ************************ */

namespace plb {
template <typename T>
SerialScalarAccess3D<T>::SerialScalarAccess3D() : locatedBlock(0)
{ }

template <typename T>
T &SerialScalarAccess3D<T>::getDistributedScalar(
    plint iX, plint iY, plint iZ, MultiBlockManagement3D const &multiBlockManagement,
    std::map<plint, ScalarField3D<T> *> &fields)
{
    plint localX, localY, localZ;
#ifdef PLB_DEBUG
    bool ok =
#endif
        multiBlockManagement.findInLocalBulk(iX, iY, iZ, locatedBlock, localX, localY, localZ);
    PLB_PRECONDITION(ok);
    return fields[locatedBlock]->get(localX, localY, localZ);
}

template <typename T>
T const &SerialScalarAccess3D<T>::getDistributedScalar(
    plint iX, plint iY, plint iZ, MultiBlockManagement3D const &multiBlockManagement,
    std::map<plint, ScalarField3D<T> *> const &fields) const
{
    plint localX, localY, localZ;
#ifdef PLB_DEBUG
    bool ok =
#endif
        multiBlockManagement.findInLocalBulk(iX, iY, iZ, locatedBlock, localX, localY, localZ);
    PLB_PRECONDITION(ok);
    return fields.find(locatedBlock)->second->get(localX, localY, localZ);
}

template <typename T>
SerialScalarAccess3D<T> *SerialScalarAccess3D<T>::clone() const
{
    return new SerialScalarAccess3D(*this);
}

/* *************** Class SerialTensorAccess3D ************************ */

template <typename T, int nDim>
SerialTensorAccess3D<T, nDim>::SerialTensorAccess3D() : locatedBlock(0)
{ }

template <typename T, int nDim>
Array<T, nDim> &SerialTensorAccess3D<T, nDim>::getDistributedTensor(
    plint iX, plint iY, plint iZ, MultiBlockManagement3D const &multiBlockManagement,
    std::map<plint, TensorField3D<T, nDim> *> &fields)
{
    plint localX, localY, localZ;
#ifdef PLB_DEBUG
    bool ok =
#endif
        multiBlockManagement.findInLocalBulk(iX, iY, iZ, locatedBlock, localX, localY, localZ);
    PLB_PRECONDITION(ok);
    return fields[locatedBlock]->get(localX, localY, localZ);
}

template <typename T, int nDim>
Array<T, nDim> const &SerialTensorAccess3D<T, nDim>::getDistributedTensor(
    plint iX, plint iY, plint iZ, MultiBlockManagement3D const &multiBlockManagement,
    std::map<plint, TensorField3D<T, nDim> *> const &fields) const
{
    plint localX, localY, localZ;
#ifdef PLB_DEBUG
    bool ok =
#endif
        multiBlockManagement.findInLocalBulk(iX, iY, iZ, locatedBlock, localX, localY, localZ);
    PLB_PRECONDITION(ok);
    return fields.find(locatedBlock)->second->get(localX, localY, localZ);
}

template <typename T, int nDim>
SerialTensorAccess3D<T, nDim> *SerialTensorAccess3D<T, nDim>::clone() const
{
    return new SerialTensorAccess3D(*this);
}

/* *************** Class SerialNTensorAccess3D ************************ */

template <typename T>
SerialNTensorAccess3D<T>::SerialNTensorAccess3D() : locatedBlock(0)
{ }

template <typename T>
T *SerialNTensorAccess3D<T>::getDistributedNTensor(
    plint iX, plint iY, plint iZ, MultiBlockManagement3D const &multiBlockManagement,
    std::map<plint, NTensorField3D<T> *> &fields)
{
    plint localX, localY, localZ;
#ifdef PLB_DEBUG
    bool ok =
#endif
        multiBlockManagement.findInLocalBulk(iX, iY, iZ, locatedBlock, localX, localY, localZ);
    PLB_PRECONDITION(ok);
    typename std::map<plint, NTensorField3D<T> *>::const_iterator it = fields.find(locatedBlock);
    PLB_ASSERT(it != fields.end());
    return it->second->get(localX, localY, localZ);
}

template <typename T>
T const *SerialNTensorAccess3D<T>::getDistributedNTensor(
    plint iX, plint iY, plint iZ, MultiBlockManagement3D const &multiBlockManagement,
    std::map<plint, NTensorField3D<T> *> const &fields) const
{
    plint localX, localY, localZ;
#ifdef PLB_DEBUG
    bool ok =
#endif
        multiBlockManagement.findInLocalBulk(iX, iY, iZ, locatedBlock, localX, localY, localZ);
    PLB_PRECONDITION(ok);
    typename std::map<plint, NTensorField3D<T> *>::const_iterator it = fields.find(locatedBlock);
    PLB_ASSERT(it != fields.end());
    return it->second->get(localX, localY, localZ);
}

template <typename T>
SerialNTensorAccess3D<T> *SerialNTensorAccess3D<T>::clone() const
{
    return new SerialNTensorAccess3D(*this);
}

}  // namespace plb

#endif  // SERIAL_MULTI_DATA_FIELD_3D_HH
