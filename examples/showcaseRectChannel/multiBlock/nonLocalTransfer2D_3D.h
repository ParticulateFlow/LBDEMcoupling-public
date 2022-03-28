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

#ifdef PLB_MPI_PARALLEL

#ifndef NON_LOCAL_TRANSFER_2D_3D_H
#define NON_LOCAL_TRANSFER_2D_3D_H

#include "core/globalDefs.h"
#include "multiBlock/nonLocalTransfer2D.h"
#include "multiBlock/nonLocalTransfer3D.h"
#include "parallelism/parallelBlockCommunicator3D.h"
#include "parallelism/sendRecvPool.h"

namespace plb {
/// Determine pairs of domains associated to a data transfer between domains on two blocks.
/** It is assumed that the two domains have the same extent.
 **/
std::vector<Overlap3D> copyDomainDataTransfer(
    SparseBlockStructure2D const &block1, Box2D block1Domain, SparseBlockStructure3D const &block2,
    Box3D block2Domain);

std::vector<Overlap3D> copyDomainDataTransfer(
    SparseBlockStructure3D const &block1, Box3D block1Domain, SparseBlockStructure2D const &block2,
    Box2D block2Domain);

struct CommunicationStructure2D_3D {
    CommunicationStructure2D_3D(
        std::vector<Overlap3D> const &overlaps, MultiBlockManagement2D const &originManagement,
        MultiBlockManagement3D const &destinationManagement, plint sizeOfCell);
    CommunicationPackage3D sendPackage;
    CommunicationPackage3D recvPackage;
    CommunicationPackage3D sendRecvPackage;
    SendPoolCommunicator sendComm;
    RecvPoolCommunicator recvComm;
};

struct CommunicationStructure3D_2D {
    CommunicationStructure3D_2D(
        std::vector<Overlap3D> const &overlaps, MultiBlockManagement3D const &originManagement,
        MultiBlockManagement2D const &destinationManagement, plint sizeOfCell);
    CommunicationPackage3D sendPackage;
    CommunicationPackage3D recvPackage;
    CommunicationPackage3D sendRecvPackage;
    SendPoolCommunicator sendComm;
    RecvPoolCommunicator recvComm;
};

/// Inter-domain copy between two generic fields.
/** The two blocks are not required to have same parallelization. If the two
 *  domains have different extent, the smaller domain is taken to be the valid one
 *  (this argument is reiterated independently in each space direction).
 *  The two blocks must be of same kind. This condition is enforced by an assertion.
 **/
void copy_generic(
    MultiBlock2D const &from, Box2D const &fromDomain, MultiBlock3D &to, Box3D const &toDomain,
    modif::ModifT typeOfModif);
void copy_generic(
    MultiBlock3D const &from, Box3D const &fromDomain, MultiBlock2D &to, Box2D const &toDomain,
    modif::ModifT typeOfModif);

/// Inter-domain copy between two scalar-fields.
/** The two blocks are not required to have same parallelization. If the two
 *  domains have different extent, the smaller domain is taken to be the valid one
 *  (this argument is reiterated independently in each space direction).
 **/
template <typename T>
void copy(
    MultiScalarField2D<T> const &from, Box2D const &fromDomain, MultiScalarField3D<T> &to,
    Box3D const &toDomain)
{
    copy_generic(from, fromDomain, to, toDomain, modif::staticVariables);
}

template <typename T>
void copy(
    MultiScalarField3D<T> const &from, Box3D const &fromDomain, MultiScalarField2D<T> &to,
    Box2D const &toDomain)
{
    copy_generic(from, fromDomain, to, toDomain, modif::staticVariables);
}

/// Inter-domain copy between two ntensor-fields.
/** The two blocks are not required to have same parallelization. If the two
 *  domains have different extent, the smaller domain is taken to be the valid one
 *  (this argument is reiterated independently in each space direction).
 **/
template <typename T>
void copy(
    MultiNTensorField2D<T> const &from, Box2D const &fromDomain, MultiNTensorField3D<T> &to,
    Box3D const &toDomain)
{
    copy_generic(from, fromDomain, to, toDomain, modif::staticVariables);
}
template <typename T>
void copy(
    MultiNTensorField3D<T> const &from, Box3D const &fromDomain, MultiNTensorField2D<T> &to,
    Box2D const &toDomain)
{
    copy_generic(from, fromDomain, to, toDomain, modif::staticVariables);
}

/// Inter-domain copy between two tensor-fields.
/** The two blocks are not required to have same parallelization. If the two
 *  domains have different extent, the smaller domain is taken to be the valid one
 *  (this argument is reiterated independently in each space direction).
 **/
template <typename T, int nDim>
void copy(
    MultiTensorField2D<T, nDim> const &from, Box2D const &fromDomain,
    MultiTensorField3D<T, nDim> &to, Box3D const &toDomain)
{
    copy_generic(from, fromDomain, to, toDomain, modif::staticVariables);
}
template <typename T, int nDim>
void copy(
    MultiTensorField3D<T, nDim> const &from, Box3D const &fromDomain,
    MultiTensorField2D<T, nDim> &to, Box2D const &toDomain)
{
    copy_generic(from, fromDomain, to, toDomain, modif::staticVariables);
}

}  // namespace plb

#endif  // NON_LOCAL_TRANSFER_2D_3D_H

#endif  // PLB_MPI_PARALLEL
