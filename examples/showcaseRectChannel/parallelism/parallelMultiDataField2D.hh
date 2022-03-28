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
#ifndef PARALLEL_MULTI_DATA_FIELD_2D_HH
#define PARALLEL_MULTI_DATA_FIELD_2D_HH

#include "parallelMultiDataField2D.h"

#ifdef PLB_MPI_PARALLEL

namespace plb {

/* *************** Class ParallelScalarAccess2D ************************ */

template <typename T>
ParallelScalarAccess2D<T>::ParallelScalarAccess2D() : locatedBlock(0)
{ }

template <typename T>
T &ParallelScalarAccess2D<T>::getDistributedScalar(
    plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
    std::map<plint, ScalarField2D<T> *> &fields)
{
    std::vector<plint> foundId;
    std::vector<plint> foundX, foundY;
    bool hasBulkCell =
        multiBlockManagement.findAllLocalRepresentations(iX, iY, foundId, foundX, foundY);
    if (hasBulkCell) {
        distributedScalar = fields[foundId[0]]->get(foundX[0], foundY[0]);
    }
    global::mpi().bCastThroughMaster(&distributedScalar, 1, hasBulkCell);
    return distributedScalar;
}

template <typename T>
T const &ParallelScalarAccess2D<T>::getDistributedScalar(
    plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
    std::map<plint, ScalarField2D<T> *> const &fields) const
{
    std::vector<plint> foundId;
    std::vector<plint> foundX, foundY;
    bool hasBulkCell =
        multiBlockManagement.findAllLocalRepresentations(iX, iY, foundId, foundX, foundY);
    if (hasBulkCell) {
        typename std::map<plint, ScalarField2D<T> *>::const_iterator it = fields.find(foundId[0]);
        distributedScalar = it->second->get(foundX[0], foundY[0]);
    }
    global::mpi().bCastThroughMaster(&distributedScalar, 1, hasBulkCell);
    return distributedScalar;
}

template <typename T>
ParallelScalarAccess2D<T> *ParallelScalarAccess2D<T>::clone() const
{
    return new ParallelScalarAccess2D(*this);
}

/* *************** Class ParallelTensorAccess2D ************************ */

template <typename T, int nDim>
ParallelTensorAccess2D<T, nDim>::ParallelTensorAccess2D() : locatedBlock(0)
{ }

template <typename T, int nDim>
Array<T, nDim> &ParallelTensorAccess2D<T, nDim>::getDistributedTensor(
    plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
    std::map<plint, TensorField2D<T, nDim> *> &fields)
{
    std::vector<plint> foundId;
    std::vector<plint> foundX, foundY;
    bool hasBulkCell =
        multiBlockManagement.findAllLocalRepresentations(iX, iY, foundId, foundX, foundY);
    if (hasBulkCell) {
        Array<T, nDim> const &foundTensor = fields[foundId[0]]->get(foundX[0], foundY[0]);
        for (int iD = 0; iD < nDim; ++iD) {
            distributedTensor[iD] = foundTensor[iD];
        }
    }
    global::mpi().bCastThroughMaster(&distributedTensor[0], nDim, hasBulkCell);
    return distributedTensor;
}

template <typename T, int nDim>
Array<T, nDim> const &ParallelTensorAccess2D<T, nDim>::getDistributedTensor(
    plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
    std::map<plint, TensorField2D<T, nDim> *> const &fields) const
{
    std::vector<plint> foundId;
    std::vector<plint> foundX, foundY;
    bool hasBulkCell =
        multiBlockManagement.findAllLocalRepresentations(iX, iY, foundId, foundX, foundY);
    if (hasBulkCell) {
        typename std::map<plint, TensorField2D<T, nDim> *>::const_iterator it =
            fields.find(foundId[0]);
        Array<T, nDim> const &foundTensor = it->second->get(foundX[0], foundY[0]);
        for (int iD = 0; iD < nDim; ++iD) {
            distributedTensor[iD] = foundTensor[iD];
        }
    }
    global::mpi().bCastThroughMaster(&distributedTensor[0], nDim, hasBulkCell);
    return distributedTensor;
}

template <typename T, int nDim>
ParallelTensorAccess2D<T, nDim> *ParallelTensorAccess2D<T, nDim>::clone() const
{
    return new ParallelTensorAccess2D(*this);
}

/* *************** Class ParallelNTensorAccess2D ************************ */

template <typename T>
ParallelNTensorAccess2D<T>::ParallelNTensorAccess2D() : distributedNTensor(0), locatedBlock(0)
{ }

template <typename T>
ParallelNTensorAccess2D<T>::~ParallelNTensorAccess2D()
{
    delete[] distributedNTensor;
}

template <typename T>
ParallelNTensorAccess2D<T>::ParallelNTensorAccess2D(ParallelNTensorAccess2D<T> const &rhs) :
    distributedNTensor(0), locatedBlock(rhs.locatedBlock)
{ }

template <typename T>
T *ParallelNTensorAccess2D<T>::getDistributedNTensor(
    plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
    std::map<plint, NTensorField2D<T> *> &fields)
{
    std::vector<plint> foundId;
    std::vector<plint> foundX, foundY;
    bool hasBulkCell =
        multiBlockManagement.findAllLocalRepresentations(iX, iY, foundId, foundX, foundY);
    typename std::map<plint, NTensorField2D<T> *>::const_iterator it = fields.find(foundId[0]);
    PLB_ASSERT(it != fields.end());
    int ndim = (int)it->second->getNdim();
    delete[] distributedNTensor;
    distributedNTensor = new T[ndim];
    if (hasBulkCell) {
        T const *foundNTensor = fields[foundId[0]]->get(foundX[0], foundY[0]);
        for (int iD = 0; iD < ndim; ++iD) {
            distributedNTensor[iD] = foundNTensor[iD];
        }
    }
    global::mpi().bCastThroughMaster(&distributedNTensor[0], ndim, hasBulkCell);
    return distributedNTensor;
}

template <typename T>
T const *ParallelNTensorAccess2D<T>::getDistributedNTensor(
    plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
    std::map<plint, NTensorField2D<T> *> const &fields) const
{
    std::vector<plint> foundId;
    std::vector<plint> foundX, foundY;
    bool hasBulkCell =
        multiBlockManagement.findAllLocalRepresentations(iX, iY, foundId, foundX, foundY);
    typename std::map<plint, NTensorField2D<T> *>::const_iterator it = fields.find(foundId[0]);
    PLB_ASSERT(it != fields.end());
    int ndim = (int)it->second->getNdim();
    delete[] distributedNTensor;
    distributedNTensor = new T[ndim];
    if (hasBulkCell) {
        T const *foundNTensor = it->second->get(foundX[0], foundY[0]);
        for (int iD = 0; iD < ndim; ++iD) {
            distributedNTensor[iD] = foundNTensor[iD];
        }
    }
    global::mpi().bCastThroughMaster(&distributedNTensor[0], ndim, hasBulkCell);
    return distributedNTensor;
}

template <typename T>
ParallelNTensorAccess2D<T> *ParallelNTensorAccess2D<T>::clone() const
{
    return new ParallelNTensorAccess2D(*this);
}

}  // namespace plb

#endif  // PLB_MPI_PARALLEL

#endif  // PARALLEL_MULTI_DATA_FIELD_2D_HH
