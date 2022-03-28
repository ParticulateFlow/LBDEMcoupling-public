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
 * Parallel access to elements of a scalar/tensor field -- generic implementation.
 */
#ifndef PARALLEL_MULTI_DATA_FIELD_3D_HH
#define PARALLEL_MULTI_DATA_FIELD_3D_HH

#include "parallelMultiDataField3D.h"

/* *************** Class ParallelScalarAccess3D ************************ */

#ifdef PLB_MPI_PARALLEL

namespace plb {

template <typename T>
ParallelScalarAccess3D<T>::ParallelScalarAccess3D() : locatedBlock(0)
{ }

template <typename T>
T &ParallelScalarAccess3D<T>::getDistributedScalar(
    plint iX, plint iY, plint iZ, MultiBlockManagement3D const &multiBlockManagement,
    std::map<plint, ScalarField3D<T> *> &fields)
{
    std::vector<plint> foundId;
    std::vector<plint> foundX, foundY, foundZ;
    bool hasBulkCell = multiBlockManagement.findAllLocalRepresentations(
        iX, iY, iZ, foundId, foundX, foundY, foundZ);
    if (hasBulkCell) {
        distributedScalar = fields[foundId[0]]->get(foundX[0], foundY[0], foundZ[0]);
    }
    global::mpi().bCastThroughMaster(&distributedScalar, 1, hasBulkCell);
    return distributedScalar;
}

template <typename T>
T const &ParallelScalarAccess3D<T>::getDistributedScalar(
    plint iX, plint iY, plint iZ, MultiBlockManagement3D const &multiBlockManagement,
    std::map<plint, ScalarField3D<T> *> const &fields) const
{
    std::vector<plint> foundId;
    std::vector<plint> foundX, foundY, foundZ;
    bool hasBulkCell = multiBlockManagement.findAllLocalRepresentations(
        iX, iY, iZ, foundId, foundX, foundY, foundZ);
    if (hasBulkCell) {
        typename std::map<plint, ScalarField3D<T> *>::const_iterator it = fields.find(foundId[0]);
        distributedScalar = it->second->get(foundX[0], foundY[0], foundZ[0]);
    }
    global::mpi().bCastThroughMaster(&distributedScalar, 1, hasBulkCell);
    return distributedScalar;
}

template <typename T>
ParallelScalarAccess3D<T> *ParallelScalarAccess3D<T>::clone() const
{
    return new ParallelScalarAccess3D(*this);
}

/* *************** Class ParallelTensorAccess3D ************************ */

template <typename T, int nDim>
ParallelTensorAccess3D<T, nDim>::ParallelTensorAccess3D() : locatedBlock(0)
{ }

template <typename T, int nDim>
Array<T, nDim> &ParallelTensorAccess3D<T, nDim>::getDistributedTensor(
    plint iX, plint iY, plint iZ, MultiBlockManagement3D const &multiBlockManagement,
    std::map<plint, TensorField3D<T, nDim> *> &fields)
{
    std::vector<plint> foundId;
    std::vector<plint> foundX, foundY, foundZ;
    bool hasBulkCell = multiBlockManagement.findAllLocalRepresentations(
        iX, iY, iZ, foundId, foundX, foundY, foundZ);
    if (hasBulkCell) {
        Array<T, nDim> const &foundTensor =
            fields[foundId[0]]->get(foundX[0], foundY[0], foundZ[0]);
        for (int iD = 0; iD < nDim; ++iD) {
            distributedTensor[iD] = foundTensor[iD];
        }
    }
    global::mpi().bCastThroughMaster(&distributedTensor[0], nDim, hasBulkCell);
    return distributedTensor;
}

template <typename T, int nDim>
Array<T, nDim> const &ParallelTensorAccess3D<T, nDim>::getDistributedTensor(
    plint iX, plint iY, plint iZ, MultiBlockManagement3D const &multiBlockManagement,
    std::map<plint, TensorField3D<T, nDim> *> const &fields) const
{
    std::vector<plint> foundId;
    std::vector<plint> foundX, foundY, foundZ;
    bool hasBulkCell = multiBlockManagement.findAllLocalRepresentations(
        iX, iY, iZ, foundId, foundX, foundY, foundZ);
    if (hasBulkCell) {
        typename std::map<plint, TensorField3D<T, nDim> *>::const_iterator it =
            fields.find(foundId[0]);
        Array<T, nDim> const &foundTensor = it->second->get(foundX[0], foundY[0], foundZ[0]);
        for (int iD = 0; iD < nDim; ++iD) {
            distributedTensor[iD] = foundTensor[iD];
        }
    }
    global::mpi().bCastThroughMaster(&distributedTensor[0], nDim, hasBulkCell);
    return distributedTensor;
}

template <typename T, int nDim>
ParallelTensorAccess3D<T, nDim> *ParallelTensorAccess3D<T, nDim>::clone() const
{
    return new ParallelTensorAccess3D(*this);
}

/* *************** Class ParallelNTensorAccess3D ************************ */

template <typename T>
ParallelNTensorAccess3D<T>::ParallelNTensorAccess3D() : distributedNTensor(0), locatedBlock(0)
{ }

template <typename T>
ParallelNTensorAccess3D<T>::~ParallelNTensorAccess3D()
{
    delete[] distributedNTensor;
}

template <typename T>
ParallelNTensorAccess3D<T>::ParallelNTensorAccess3D(ParallelNTensorAccess3D<T> const &rhs) :
    distributedNTensor(0), locatedBlock(rhs.locatedBlock)
{ }

template <typename T>
T *ParallelNTensorAccess3D<T>::getDistributedNTensor(
    plint iX, plint iY, plint iZ, MultiBlockManagement3D const &multiBlockManagement,
    std::map<plint, NTensorField3D<T> *> &fields)
{
    std::vector<plint> foundId;
    std::vector<plint> foundX, foundY, foundZ;
    bool hasBulkCell = multiBlockManagement.findAllLocalRepresentations(
        iX, iY, iZ, foundId, foundX, foundY, foundZ);
    int ndim = (int)fields[foundId[0]]->getNdim();
    delete[] distributedNTensor;
    distributedNTensor = new T[ndim];
    if (hasBulkCell) {
        T const *foundNTensor = fields[foundId[0]]->get(foundX[0], foundY[0], foundZ[0]);
        for (int iD = 0; iD < ndim; ++iD) {
            distributedNTensor[iD] = foundNTensor[iD];
        }
    }
    global::mpi().bCastThroughMaster(&distributedNTensor[0], ndim, hasBulkCell);
    return distributedNTensor;
}

template <typename T>
T const *ParallelNTensorAccess3D<T>::getDistributedNTensor(
    plint iX, plint iY, plint iZ, MultiBlockManagement3D const &multiBlockManagement,
    std::map<plint, NTensorField3D<T> *> const &fields) const
{
    std::vector<plint> foundId;
    std::vector<plint> foundX, foundY, foundZ;
    bool hasBulkCell = multiBlockManagement.findAllLocalRepresentations(
        iX, iY, iZ, foundId, foundX, foundY, foundZ);
    typename std::map<plint, NTensorField3D<T> *>::const_iterator it = fields.find(foundId[0]);
    PLB_ASSERT(it != fields.end());
    int ndim = (int)it->second->getNdim();
    delete[] distributedNTensor;
    distributedNTensor = new T[ndim];
    if (hasBulkCell) {
        T const *foundNTensor = it->second->get(foundX[0], foundY[0], foundZ[0]);
        for (int iD = 0; iD < ndim; ++iD) {
            distributedNTensor[iD] = foundNTensor[iD];
        }
    }
    global::mpi().bCastThroughMaster(&distributedNTensor[0], ndim, hasBulkCell);
    return distributedNTensor;
}

template <typename T>
ParallelNTensorAccess3D<T> *ParallelNTensorAccess3D<T>::clone() const
{
    return new ParallelNTensorAccess3D(*this);
}

}  // namespace plb

#endif  // PLB_MPI_PARALLEL

#endif  // PARALLEL_MULTI_DATA_FIELD_3D_HH
