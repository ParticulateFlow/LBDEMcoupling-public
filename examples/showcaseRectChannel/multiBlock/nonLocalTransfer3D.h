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

#ifndef NON_LOCAL_TRANSFER_3D_H
#define NON_LOCAL_TRANSFER_3D_H

#include <vector>

#include "core/globalDefs.h"
#include "multiBlock/localMultiBlockInfo3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"

namespace plb {

/// Determine pairs of domains associated to a data transfer between two blocks.
std::vector<Overlap3D> copyAllDataTransfer(
    SparseBlockStructure3D const &block1, SparseBlockStructure3D const &block2);

/// Determine pairs of domains associated to a data transfer between domains on two blocks.
/** It is assumed that the two domains have the same extent.
 **/
std::vector<Overlap3D> copyDomainDataTransfer(
    SparseBlockStructure3D const &block1, Box3D block1Domain, SparseBlockStructure3D const &block2,
    Box3D block2Domain);

/// Equal-domain copy between two scalar-fields with different parallelization.
/// OBSOLETE: use copy instead.
template <typename T>
void copyNonLocal(
    MultiScalarField3D<T> const &from, MultiScalarField3D<T> &to, Box3D const &domain);

/// Inter-domain copy between two scalar-fields.
/** The two blocks are not required to have same parallelization. If the two
 *  domains have different extent, the smaller domain is taken to be the valid one
 *  (this argument is reiterated independently in each space direction).
 **/
template <typename T>
void copy(
    MultiScalarField3D<T> const &from, Box3D const &fromDomain, MultiScalarField3D<T> &to,
    Box3D const &toDomain);

/// Equal-domain copy between two ntensor-fields with different parallelization.
/// OBSOLETE: use copy instead.
template <typename T>
void copyNonLocal(
    MultiNTensorField3D<T> const &from, MultiNTensorField3D<T> &to, Box3D const &domain);

/// Inter-domain copy between two ntensor-fields.
/** The two blocks are not required to have same parallelization. If the two
 *  domains have different extent, the smaller domain is taken to be the valid one
 *  (this argument is reiterated independently in each space direction).
 **/
template <typename T>
void copy(
    MultiNTensorField3D<T> const &from, Box3D const &fromDomain, MultiNTensorField3D<T> &to,
    Box3D const &toDomain);

/// Equal-domain copy between two tensor-fields with different parallelization.
/// OBSOLETE: use copy instead.
template <typename T, int nDim>
void copyNonLocal(
    MultiTensorField3D<T, nDim> const &from, MultiTensorField3D<T, nDim> &to, Box3D const &domain);

/// Inter-domain copy between two tensor-fields.
/** The two blocks are not required to have same parallelization. If the two
 *  domains have different extent, the smaller domain is taken to be the valid one
 *  (this argument is reiterated independently in each space direction).
 **/
template <typename T, int nDim>
void copy(
    MultiTensorField3D<T, nDim> const &from, Box3D const &fromDomain,
    MultiTensorField3D<T, nDim> &to, Box3D const &toDomain);

/// Equal-domain copy between two block-lattices with different parallelization.
/// OBSOLETE: use copy instead.
template <typename T, template <typename U> class Descriptor>
void copyNonLocal(
    MultiBlockLattice3D<T, Descriptor> const &from, MultiBlockLattice3D<T, Descriptor> &to,
    Box3D const &domain, modif::ModifT whichContent);

/// Inter-domain copy between two block-lattices.
/** The two blocks are not required to have same parallelization. If the two
 *  domains have different extent, the smaller domain is taken to be the valid one
 *  (this argument is reiterated independently in each space direction).
 **/
template <typename T, template <typename U> class Descriptor>
void copy(
    MultiBlockLattice3D<T, Descriptor> const &from, Box3D const &fromDomain,
    MultiBlockLattice3D<T, Descriptor> &to, Box3D const &toDomain, modif::ModifT whichContent);

/// Inter-domain copy of the populations and external scalars between two block-lattices.
/** The two blocks are not required to have same parallelization. If the two
 *  domains have different extent, the smaller domain is taken to be the valid one
 *  (this argument is reiterated independently in each space direction).
 **/
template <typename T, template <typename U> class Descriptor>
void copyPopulations(
    MultiBlockLattice3D<T, Descriptor> const &from, Box3D const &fromDomain,
    MultiBlockLattice3D<T, Descriptor> &to, Box3D const &toDomain);

/// Inter-domain copy of the populations, external scalars, and dynamic variables between two
/// block-lattices.
/** The two blocks are not required to have same parallelization. If the two
 *  domains have different extent, the smaller domain is taken to be the valid one
 *  (this argument is reiterated independently in each space direction).
 **/
template <typename T, template <typename U> class Descriptor>
void copyAll(
    MultiBlockLattice3D<T, Descriptor> const &from, Box3D const &fromDomain,
    MultiBlockLattice3D<T, Descriptor> &to, Box3D const &toDomain);

/// Inter-domain copy of all variables, and reconstitution of the dynamics objects, between two
/// block-lattices.
/** The two blocks are not required to have same parallelization. If the two
 *  domains have different extent, the smaller domain is taken to be the valid one
 *  (this argument is reiterated independently in each space direction).
 **/
template <typename T, template <typename U> class Descriptor>
void copyRegenerate(
    MultiBlockLattice3D<T, Descriptor> const &from, Box3D const &fromDomain,
    MultiBlockLattice3D<T, Descriptor> &to, Box3D const &toDomain);

/// Inter-domain copy between two generic fields.
/** The two blocks are not required to have same parallelization. If the two
 *  domains have different extent, the smaller domain is taken to be the valid one
 *  (this argument is reiterated independently in each space direction).
 *  The two blocks must be of same kind. This condition is enforced by an assertion.
 **/
void copy_generic(
    MultiBlock3D const &from, Box3D const &fromDomain, MultiBlock3D &to, Box3D const &toDomain,
    modif::ModifT whichContent);

/// Equal-domain copy between two generic fields with different parallelization.
/** The two blocks must be of same kind. This condition is enforced by an assertion.  **/
void copyNonLocal_generic(
    MultiBlock3D const &from, MultiBlock3D &to, Box3D const &domain, modif::ModifT whichContent);

}  // namespace plb

#endif  // NON_LOCAL_TRANSFER_3D_H
