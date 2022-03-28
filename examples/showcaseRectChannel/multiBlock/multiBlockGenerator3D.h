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
 * Copy 3D multiblocks on a new parallel distribution -- header file.
 */

/*
 * All functions that take as an argument a Box3D, or compute intersections, or join,
 * or crop, do not adjust periodicity of the newly created blocks by default. The user
 * is responsible to take care of this matter explicitly.
 */

#ifndef MULTI_BLOCK_GENERATOR_3D_H
#define MULTI_BLOCK_GENERATOR_3D_H

#include <memory>

#include "core/globalDefs.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiContainerBlock3D.h"
#include "multiBlock/multiDataField3D.h"
#include "multiBlock/sparseBlockStructure3D.h"
#include "particles/multiParticleField3D.h"

namespace plb {

/* *************** 1. MultiScalarField ************************************** */

/// Generate a multi-scalar-field from scratch. As opposed to the standard
///   constructor, this factory function takes a full bounding-box, as well
///   as the envelope-width, as arguments.
template <typename T>
std::unique_ptr<MultiScalarField3D<T> > generateMultiScalarField(
    Box3D boundingBox, plint envelopeWidth = 1);

/// Generate a multi-scalar-field from scratch. As opposed to the standard
///   constructor, this factory function takes a full bounding-box, as well
///   as the envelope-width, as arguments.
///   The parameter envelopeWidth has no default value to avoid confusion
///   with the two-parameter version of this overloaded function.
template <typename T>
std::unique_ptr<MultiScalarField3D<T> > generateMultiScalarField(
    Box3D boundingBox, T iniVal, plint envelopeWidth);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > generateMultiScalarField(
    MultiBlock3D &multiBlock, plint envelopeWidth);

/// Generate a multi-scalar-field from scratch. As opposed to the standard
///   constructor, this factory function takes the explicit block-management
///   object, which includes stuff like block-distribution, parallelization,
///   envelope-width, etc. An optional initialization value can be provided.
template <typename T>
std::unique_ptr<MultiScalarField3D<T> > defaultGenerateMultiScalarField3D(
    MultiBlockManagement3D const &management, T iniVal = T());

/// Create a clone of a MultiScalarField (or of a sub-domain).
/** This cannot be handled through a data processor, because the internal data
 *  processors of the scalar-field must be cloned as well.
 **/
template <typename T>
std::unique_ptr<MultiScalarField3D<T> > clone(
    MultiScalarField3D<T> &originalField, Box3D const &subDomain, bool crop = true);

/// Generate a new multi-scalar-field with same distribution as original block,
///   but intersected with a given domain.
/** No data is being copied, the new scalar-field is initialized to zero values.
 *    The new boundingBox is either equal to the original one, or equal to the
 *    intersection, depending on the parameter crop.
 **/
template <typename T>
std::unique_ptr<MultiScalarField3D<T> > generateMultiScalarField(
    MultiBlock3D const &originalField, Box3D const &intersection, bool crop = true);

/// Generate a new multi-scalar-field with a distribution which is the
///   intersection of the two original distributions.
/** No data is being copied, the new scalar-field is initialized to zero values.
 *    The block identifiers and parallel distribution are taken from
 *    originalField1. The new bounding-box is equal to the intersection
 *    of the original ones if the parameter "crop" is true, and equal to
 *    the bound around the original ones if "crop" is false.
 **/
template <typename T>
std::unique_ptr<MultiScalarField3D<T> > generateIntersectMultiScalarField(
    MultiBlock3D const &originalField1, MultiBlock3D const &originalField2, bool crop = true);

/// Generate a new multi-scalar-field with a distribution which is the
///   intersection of the two original distributions and a given domain.
/** No data is being copied, the new scalar-field is initialized to zero values.
 *    The block identifiers and parallel distribution are taken from
 *    originalField1. The new bounding-box is equal to the intersection
 *    of the original fields and the domain if the parameter "crop" is
 *    true, and equal to the bound around the original fields if "crop"
 *    is false.
 **/
template <typename T>
std::unique_ptr<MultiScalarField3D<T> > generateIntersectMultiScalarField(
    MultiBlock3D const &originalField1, MultiBlock3D const &originalField2,
    Box3D const &intersection, bool crop = true);

/// Generate a new multi-scalar-field with a distribution which is the
///   union of the two original distributions.
/** No data is being copied, the new scalar-field is initialized to zero values.
 *    The block identifiers and parallel distributions are taken from
 *    originalField1 wherever originalField1 is defined. Elsewhere,
 *    new identifiers are generated, and the parallel disribution is
 *    taken from originalField2.
 **/
template <typename T>
std::unique_ptr<MultiScalarField3D<T> > generateJoinMultiScalarField(
    MultiBlock3D const &originalField1, MultiBlock3D const &originalField2);

/// Generate a new multi-scalar-field with a distribution which is the
/// the one of the orignal block, but with a new envelopeWidth
template <typename T>
std::unique_ptr<MultiScalarField3D<T> > extendEnvelopeWidth(
    MultiScalarField3D<T> &originalBlock, plint envelopeWidth);

/// Generate a new multi-scalar-field with a distribution which is the
///   union between the original one and an additional block.
template <typename T>
std::unique_ptr<MultiScalarField3D<T> > extend(
    MultiScalarField3D<T> &originalBlock, Box3D const &addedBlock);

/// Generate a new multi-scalar-field with a distribution which is the
///   union between the original one with exception of the specified block.
template <typename T>
std::unique_ptr<MultiScalarField3D<T> > except(
    MultiScalarField3D<T> &originalBlock, Box3D const &exceptedBlock);

/// Create a clone of the original field with a different block-distribution.
template <typename T>
std::unique_ptr<MultiScalarField3D<T> > redistribute(
    MultiScalarField3D<T> const &originalField, SparseBlockStructure3D const &newBlockStructure,
    bool adjustPeriodicity = true);

/// Create a clone of the original field on the domain of intersection,
///   with a different block-distribution.
/**   The new boundingBox is either equal to the original one, or equal to the
 *    intersection, depending on the parameter crop.
 **/
template <typename T>
std::unique_ptr<MultiScalarField3D<T> > redistribute(
    MultiScalarField3D<T> const &originalField, SparseBlockStructure3D const &newBlockStructure,
    Box3D const &intersection, bool crop = true);

/// Create a clone of the original field, re-distributed in such a way
///   that it is aligned with the partner-block, i.e. it can be used
///   together with partner-block from within a data processor.
template <typename T>
std::unique_ptr<MultiScalarField3D<T> > align(
    MultiScalarField3D<T> const &originalBlock, MultiBlock3D const &partnerBlock);

/// Create a clone of the original field, on a new regular distribution.
template <typename T>
std::unique_ptr<MultiScalarField3D<T> > reparallelize(MultiScalarField3D<T> const &originalBlock);

/// Create a clone of the original field, on a new regular distribution.
/** The parameters blockLx, blockLy, and blockLz indicate the approximate size of the
 *  blocks.
 **/
template <typename T>
std::unique_ptr<MultiScalarField3D<T> > reparallelize(
    MultiScalarField3D<T> const &originalBlock, plint blockLx, plint blockLy, plint blockLz);

/* *************** 2. MultiNTensorField ************************************** */

/// Generate a multi-scalar-field from scratch. As opposed to the standard
///   constructor, this factory function takes the explicit block-management
///   object, which includes stuff like block-distribution, parallelization,
///   envelope-width, etc.
template <typename T>
std::unique_ptr<MultiNTensorField3D<T> > defaultGenerateMultiNTensorField3D(
    MultiBlockManagement3D const &management, plint nDim = 1);

template <typename T>
MultiNTensorField3D<T> *generateMultiNTensorField3D(
    MultiBlock3D &multiBlock, plint envelopeWidth, plint ndim);

template <typename T>
MultiNTensorField3D<T> *generateMultiNTensorField3D(Box3D const &domain, plint ndim);

template <typename T>
MultiNTensorField3D<T> *generateMultiNTensorField3D(
    Box3D const &domain, plint ndim, T *iniVal, plint envelopeWidth);

/// Create a clone of a MultiNTensorField (or of a sub-domain).
/** This cannot be handled through a data processor, because the internal data
 *  processors of the ntensor-field must be cloned as well.
 **/
template <typename T>
MultiNTensorField3D<T> *clone(
    MultiNTensorField3D<T> &originalField, Box3D const &subDomain, bool crop = true);

/// Generate a new multi-tensor-field with same distribution as original block,
///   but intersected with a given domain.
/** No data is being copied, the new tensor-field is initialized to zero values.
 *    The new boundingBox is either equal to the original one, or equal to the
 *    intersection, depending on the parameter crop.
 **/
template <typename T>
MultiNTensorField3D<T> *generateMultiNTensorField(
    MultiBlock3D const &originalField, Box3D const &intersection, plint nDim, bool crop = true);

template <typename T1, typename T2>
MultiNTensorField3D<T2> *generateNTensorFieldFromNTensor3D(
    MultiNTensorField3D<T1> const &field, Box3D const &intersection, plint nDim);

template <typename T1, typename T2, template <typename U> class Descriptor>
MultiNTensorField3D<T1> *generateNTensorFieldFromBlockLattice3D(
    MultiBlockLattice3D<T2, Descriptor> const &lattice, Box3D const &intersection, plint nDim);

/// Generate a new multi-tensor-field with a distribution which is the
///   intersection of the two original distributions.
/** No data is being copied, the new tensor-field is initialized to zero values.
 *    The block identifiers and parallel distribution are taken from
 *    originalField1. The new bounding-box is equal to the intersection
 *    of the original ones if the parameter "crop" is true, and equal to
 *    the bound around the original ones if "crop" is false.
 **/
template <typename T>
MultiNTensorField3D<T> *generateIntersectMultiNTensorField(
    MultiBlock3D const &originalField1, MultiBlock3D const &originalField2, plint nDim,
    bool crop = true);

/// Generate a new multi-tensor-field with a distribution which is the
///   intersection of the two original distributions and a given domain.
/** No data is being copied, the new tensor-field is initialized to zero values.
 *    The block identifiers and parallel distribution are taken from
 *    originalField1. The new bounding-box is equal to the intersection
 *    of the original fields and the domain if the parameter "crop" is
 *    true, and equal to the bound around the original fields if "crop"
 *    is false.
 **/
template <typename T>
MultiNTensorField3D<T> *generateIntersectMultiNTensorField(
    MultiBlock3D const &originalField1, MultiBlock3D const &originalField2,
    Box3D const &intersection, plint nDim, bool crop = true);

/// Generate a new multi-tensor-field with a distribution which is the
///   union of the two original distributions.
/** No data is being copied, the new tensor-field is initialized to zero values.
 *    The block identifiers and parallel distributions are taken from
 *    originalField1 wherever originalField1 is defined. Elsewhere,
 *    new identifiers are generated, and the parallel disribution is
 *    taken from originalField2.
 **/
template <typename T>
MultiNTensorField3D<T> *generateJoinMultiNTensorField(
    MultiBlock3D const &originalField1, MultiBlock3D const &originalField2, plint nDim);

/// Generate a new multi-NTensor-field with a distribution which is the
/// the one of the orignal block, but with a new envelopeWidth
template <typename T>
std::unique_ptr<MultiNTensorField3D<T> > extendEnvelopeWidth(
    MultiNTensorField3D<T> &originalBlock, plint envelopeWidth);

/// Generate a new multi-ntensor-field with a distribution which is the
///   union between the original one and an additional block.
template <typename T>
MultiNTensorField3D<T> *extend(MultiNTensorField3D<T> &originalBlock, Box3D const &addedBlock);

/// Generate a new multi-ntensor-field with a distribution which is the
///   union between the original one with exception of the specified block.
template <typename T>
MultiNTensorField3D<T> *except(MultiNTensorField3D<T> &originalBlock, Box3D const &exceptedBlock);

/// Create a clone of the original field, re-distributed in such a way
///   that it is aligned with the partner-block, i.e. it can be used
///   together with partner-block from within a data processor.
template <typename T>
MultiNTensorField3D<T> *align(
    MultiNTensorField3D<T> const &originalBlock, MultiBlock3D const &partnerBlock);

/// Create a clone of the original field, on a new regular distribution.
template <typename T>
MultiNTensorField3D<T> *reparallelize(MultiNTensorField3D<T> const &originalBlock);

/// Create a clone of the original field, on a new regular distribution.
/** The parameters blockLx, blockLy, and blockLz indicate the approximate size of the
 *  blocks.
 **/
template <typename T>
MultiNTensorField3D<T> *reparallelize(
    MultiNTensorField3D<T> const &originalBlock, plint blockLx, plint blockLy, plint blockLz);

/* *************** 3. MultiTensorField ************************************** */

/// Generate a multi-tensor-field from scratch. As opposed to the standard
///   constructor, this factory function takes a full bounding-box, as well
///   as the envelope-width, as arguments.
template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > generateMultiTensorField(
    Box3D boundingBox, plint envelopeWidth = 1);

/// Generate a multi-tensor-field from scratch. As opposed to the standard
///   constructor, this factory function takes a full bounding-box, as well
///   as the envelope-width, as arguments.
template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > generateMultiTensorField(
    Box3D boundingBox, Array<T, nDim> const &iniVal, plint envelopeWidth = 1);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > generateMultiTensorField(
    MultiBlock3D &multiBlock, plint envelopeWidth);

/// Generate a multi-tensor-field from scratch. As opposed to the standard
///   constructor, this factory function takes the explicit block-management
///   object, which includes stuff like block-distribution, parallelization,
///   envelope-width, etc. A default dummy argument is used for technical reasons.
template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > defaultGenerateMultiTensorField3D(
    MultiBlockManagement3D const &management, plint unnamedDummyArg = 1);

/// Create a clone of a MultiTensorField (or of a sub-domain).
/** This cannot be handled through a data processor, because the internal data
 *  processors of the tensor-field must be cloned as well.
 **/
template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > clone(
    MultiTensorField3D<T, nDim> &originalField, Box3D const &subDomain, bool crop = true);

/// Generate a new multi-tensor-field with same distribution as original block,
///   but intersected with a given domain.
/** No data is being copied, the new tensor-field is initialized to zero values.
 *    The new boundingBox is either equal to the original one, or equal to the
 *    intersection, depending on the parameter crop.
 **/
template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > generateMultiTensorField(
    MultiBlock3D const &originalField, Box3D const &intersection, bool crop = true);

/// Generate a new multi-tensor-field with a distribution which is the
///   intersection of the two original distributions.
/** No data is being copied, the new tensor-field is initialized to zero values.
 *    The block identifiers and parallel distribution are taken from
 *    originalField1. The new bounding-box is equal to the intersection
 *    of the original ones if the parameter "crop" is true, and equal to
 *    the bound around the original ones if "crop" is false.
 **/
template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > generateIntersectMultiTensorField(
    MultiBlock3D const &originalField1, MultiBlock3D const &originalField2, bool crop = true);

/// Generate a new multi-tensor-field with a distribution which is the
///   intersection of the two original distributions and a given domain.
/** No data is being copied, the new tensor-field is initialized to zero values.
 *    The block identifiers and parallel distribution are taken from
 *    originalField1. The new bounding-box is equal to the intersection
 *    of the original fields and the domain if the parameter "crop" is
 *    true, and equal to the bound around the original fields if "crop"
 *    is false.
 **/
template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > generateIntersectMultiTensorField(
    MultiBlock3D const &originalField1, MultiBlock3D const &originalField2,
    Box3D const &intersection, bool crop = true);

/// Generate a new multi-tensor-field with a distribution which is the
///   union of the two original distributions.
/** No data is being copied, the new tensor-field is initialized to zero values.
 *    The block identifiers and parallel distributions are taken from
 *    originalField1 wherever originalField1 is defined. Elsewhere,
 *    new identifiers are generated, and the parallel disribution is
 *    taken from originalField2.
 **/
template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > generateJoinMultiTensorField(
    MultiBlock3D const &originalField1, MultiBlock3D const &originalField2);

/// Generate a new multi-tensor-field with a distribution which is the
/// the one of the orignal block, but with a new envelopeWidth
template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > extendEnvelopeWidth(
    MultiTensorField3D<T, nDim> &originalBlock, plint envelopeWidth);

/// Generate a new multi-tensor-field with a distribution which is the
///   union between the original one and an additional block.
template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > extend(
    MultiTensorField3D<T, nDim> &originalBlock, Box3D const &addedBlock);

/// Generate a new multi-tensor-field with a distribution which is the
///   union between the original one with exception of the specified block.
template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > except(
    MultiTensorField3D<T, nDim> &originalBlock, Box3D const &exceptedBlock);

/// Create a clone of the original field with a different block-distribution.
template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > redistribute(
    MultiTensorField3D<T, nDim> const &originalField,
    SparseBlockStructure3D const &newBlockStructure, bool adjustPeriodicity = true);

/// Create a clone of the original field on the domain of intersection,
///   with a different block-distribution.
/**   The new boundingBox is either equal to the original one, or equal to the
 *    intersection, depending on the parameter crop.
 **/
template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > redistribute(
    MultiTensorField3D<T, nDim> const &originalField,
    SparseBlockStructure3D const &newBlockStructure, Box3D const &intersection, bool crop = true);

/// Create a clone of the original field, re-distributed in such a way
///   that it is aligned with the partner-block, i.e. it can be used
///   together with partner-block from within a data processor.
template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > align(
    MultiTensorField3D<T, nDim> const &originalBlock, MultiBlock3D const &partnerBlock);

/// Create a clone of the original field, on a new regular distribution.
template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > reparallelize(
    MultiTensorField3D<T, nDim> const &originalBlock);

/// Create a clone of the original field, on a new regular distribution.
/** The parameters blockLx, blockLy, and blockLz indicate the approximate size of the
 *  blocks.
 **/
template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > reparallelize(
    MultiTensorField3D<T, nDim> const &originalBlock, plint blockLx, plint blockLy, plint blockLz);

/* *************** 4. MultiBlockLattice ************************************** */

/// Generate a multi-block-lattice from scratch. As opposed to the standard
///   constructor, this factory function takes a full bounding-box, as well
///   as the envelope-width, as arguments.
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > generateMultiBlockLattice(
    Box3D boundingBox, Dynamics<T, Descriptor> *backgroundDynamics, plint envelopeWidth = 1);

/// Generate a multi-block-lattice from scratch. As opposed to the standard
///   constructor, this factory function takes the explicit block-management
///   object, which includes stuff like block-distribution, parallelization,
///   envelope-width, etc. A default dummy argument is used for technical reasons.
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > defaultGenerateMultiBlockLattice3D(
    MultiBlockManagement3D const &management, plint unnamedDummyArg = 1);

/// Generate a multi-block-lattice from scratch. As opposed to the standard
///   constructor, this factory function takes the explicit block-management
///   object, which includes stuff like block-distribution, parallelization,
///   envelope-width, etc. An optional background dynamics is also provided.
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > generateMultiBlockLattice(
    MultiBlockManagement3D const &management, Dynamics<T, Descriptor> *backgroundDynamics = 0);

/// Create a clone of a MultiBlockLattice (or of a sub-domain).
/** This cannot be handled through a data processor, because the internal data
 *  processors of the block-lattice must be cloned as well.
 **/
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > clone(
    MultiBlockLattice3D<T, Descriptor> &originalLattice, Box3D const &subDomain, bool crop = true);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > generateMultiBlockLattice(
    MultiBlock3D &multiBlock, plint envelopeWidth, Dynamics<T, Descriptor> *backgroundDynamics = 0);

/// Generate a new multi-block-lattice with same distribution as original block,
///   but intersected with a given domain.
/** No data is being copied, the new block-lattice is initialized to default values.
 *    The new boundingBox is either equal to the original one, or equal to the
 *    intersection, depending on the parameter crop.
 **/
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > generateMultiBlockLattice(
    MultiBlock3D const &originalBlock, Box3D const &intersection, bool crop = true);

/// Generate a new multi-block-lattice with a distribution which is the
///   intersection of the two original distributions.
/** No data is being copied, the new block-lattice is initialized to default values.
 *    The block identifiers and parallel distribution are taken from
 *    originalBlock1. The new bounding-box is equal to the intersection
 *    of the original ones if the parameter "crop" is true, and equal to
 *    the bound around the original ones if "crop" is false.
 **/
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > generateIntersectMultiBlockLattice(
    MultiBlock3D const &originalBlock1, MultiBlock3D const &originalBlock2, bool crop = true);

/// Generate a new multi-block-lattice with a distribution which is the
///   intersection of the two original distributions and a given domain.
/** No data is being copied, the new block-lattice is initialized to default values.
 *    The block identifiers and parallel distribution are taken from
 *    originalBlock1. The new bounding-box is equal to the intersection
 *    of the original lattices and the domain if the parameter "crop" is
 *    true, and equal to the bound around the original lattices if "crop"
 *    is false.
 **/
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > generateIntersectMultiBlockLattice(
    MultiBlock3D const &originalBlock1, MultiBlock3D const &originalBlock2,
    Box3D const &intersection, bool crop = true);

/// Generate a new multi-block-lattice with a distribution which is the
///   union of the two original distributions.
/** No data is being copied, the new block-lattice is initialized to default values.
 *    The block identifiers and parallel distributions are taken from
 *    originalBlock1 wherever originalBlock1 is defined. Elsewhere,
 *    new identifiers are generated, and the parallel disribution is
 *    taken from originalBlock2.
 **/
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > generateJoinMultiBlockLattice(
    MultiBlock3D const &originalBlock1, MultiBlock3D const &originalBlock2);

/// Generate a new multi-block-lattice with a distribution which is the
///   union between the original one and an additional block.
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > extend(
    MultiBlockLattice3D<T, Descriptor> &originalBlock, Box3D const &addedBlock);

/// Generate a new multi-block-lattice with a distribution which is the
///   union between the original one with exception of the specified block.
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > except(
    MultiBlockLattice3D<T, Descriptor> &originalBlock, Box3D const &exceptedBlock);

/// Create a clone of the original lattice with a different block-distribution.
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > redistribute(
    MultiBlockLattice3D<T, Descriptor> const &originalBlock,
    SparseBlockStructure3D const &newBlockStructure, bool adjustPeriodicity = true);

/// Create a clone of the original lattice on the domain of intersection,
///   with a different block-distribution.
/**   The new boundingBox is either equal to the original one, or equal to the
 *    intersection, depending on the parameter crop.
 **/
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > redistribute(
    MultiBlockLattice3D<T, Descriptor> const &originalBlock,
    SparseBlockStructure3D const &newBlockStructure, Box3D const &intersection, bool crop = true);

/// Create a clone of the original lattice, re-distributed in such a way
///   that it is aligned with the partner-block, i.e. it can be used
///   together with partner-block from within a data processor.
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > align(
    MultiBlockLattice3D<T, Descriptor> const &originalBlock, MultiBlock3D const &partnerBlock);

/// Create a clone of the original lattice, on a new regular distribution.
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > reparallelize(
    MultiBlockLattice3D<T, Descriptor> const &originalBlock);

/// Create a clone of the original lattice, on a new regular distribution.
/** The parameters blockLx, blockLy, and blockLz indicate the approximate size of the
 *  blocks.
 **/
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > reparallelize(
    MultiBlockLattice3D<T, Descriptor> const &originalBlock, plint blockLx, plint blockLy,
    plint blockLz);

/* *************** 4. MultiParticleField ************************************ */

/// Generate a multi-particle-field from scratch. As opposed to the standard
///   constructor, this factory function takes a full bounding-box, as well
///   as the envelope-width, as arguments.
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiParticleField3D<DenseParticleField3D<T, Descriptor> > >
    generateMultiDenseParticleField(Box3D boundingBox, plint envelopeWidth = 1);

template <typename T, template <typename U> class Descriptor, class ParticleFieldT>
std::unique_ptr<MultiParticleField3D<ParticleFieldT> > generateMultiParticleField3D(
    Box3D boundingBox, plint envelopeWidth = 1);

template <typename T, template <typename U> class Descriptor, class ParticleFieldT>
std::unique_ptr<MultiParticleField3D<ParticleFieldT> > generateMultiParticleField3D(
    MultiBlock3D &multiBlock, plint envelopeWidth);

template <typename T, template <typename U> class Descriptor, class ParticleFieldT>
std::unique_ptr<MultiParticleField3D<ParticleFieldT> > generateMultiParticleField3D(
    MultiBlockManagement3D const &management, PeriodicitySwitch3D const &periodicity,
    plint envelopeWidth);

template <typename T, template <typename U> class Descriptor, class ParticleFieldT>
std::unique_ptr<MultiParticleField3D<ParticleFieldT> > generateMultiParticleField3D(
    MultiBlockManagement3D const &management, plint envelopeWidth);

/* *************** 5. MultiContainerBlock ************************************ */

std::unique_ptr<MultiContainerBlock3D> generateMultiContainerBlock(
    MultiBlock3D &multiBlock, plint envelopeWidth);

MultiContainerBlock3D *createMultiContainerBlock3D(
    MultiBlockManagement3D const &management, PeriodicitySwitch3D const &periodicity,
    plint envelopeWidth, plint gridLevel);

MultiContainerBlock3D *createMultiContainerBlock3D(
    MultiBlockManagement3D const &management, plint envelopeWidth, plint gridLevel);

/* *************** General Functions **************************************** */

/// Copy all data processors from the old to the new block, and replace
///   the "self" block in the list of the data processor to point to the new block.
void transferDataProcessors(MultiBlock3D const &from, MultiBlock3D &to);

}  // namespace plb

#endif  // MULTI_BLOCK_GENERATOR_3D_H
