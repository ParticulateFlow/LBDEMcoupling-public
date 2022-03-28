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
 * Copy 2D multiblocks on a new parallel distribution -- header file.
 */

/*
 * All functions that take as an argument a Box2D, or compute intersections, or join,
 * or crop, do not adjust periodicity of the newly created blocks by default. The user
 * is responsible to take care of this matter explicitly.
 */

#ifndef MULTI_BLOCK_GENERATOR_2D_H
#define MULTI_BLOCK_GENERATOR_2D_H

#include <memory>

#include "core/globalDefs.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiContainerBlock2D.h"
#include "multiBlock/multiDataField2D.h"
#include "multiBlock/sparseBlockStructure2D.h"
#include "particles/multiParticleField2D.h"

namespace plb {

/* *************** 1. MultiScalarField ************************************** */

/// Generate a multi-scalar-field from scratch. As opposed to the standard
///   constructor, this factory function takes a full bounding-box, as well
///   as the envelope-width, as arguments.
template <typename T>
std::unique_ptr<MultiScalarField2D<T> > generateMultiScalarField(
    Box2D boundingBox, plint envelopeWidth = 1);

/// Generate a multi-scalar-field from scratch. As opposed to the standard
///   constructor, this factory function takes a full bounding-box, as well
///   as the envelope-width, as arguments.
template <typename T>
std::unique_ptr<MultiScalarField2D<T> > generateMultiScalarField(
    Box2D boundingBox, T iniVal, plint envelopeWidth = 1);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > generateMultiScalarField(
    MultiBlock2D &multiBlock, plint envelopeWidth);

/// Generate a multi-scalar-field from scratch. As opposed to the standard
///   constructor, this factory function takes the explicit block-management
///   object, which includes stuff like block-distribution, parallelization,
///   envelope-width, etc. An optional initialization value can be provided.
template <typename T>
std::unique_ptr<MultiScalarField2D<T> > defaultGenerateMultiScalarField2D(
    MultiBlockManagement2D const &management, T iniVal = T());

/// Create a clone of a MultiScalarField (or of a sub-domain).
/** This cannot be handled through a data processor, because the internal data
 *  processors of the scalar-field must be cloned as well.
 **/
template <typename T>
std::unique_ptr<MultiScalarField2D<T> > clone(
    MultiScalarField2D<T> &originalField, Box2D const &subDomain, bool crop = true);

/// Generate a new multi-scalar-field with same distribution as original block,
///   but intersected with a given domain.
/** No data is being copied, the new scalar-field is initialized to zero values.
 *    The new boundingBox is either equal to the original one, or equal to the
 *    intersection, depending on the parameter crop.
 **/
template <typename T>
std::unique_ptr<MultiScalarField2D<T> > generateMultiScalarField(
    MultiBlock2D const &originalField, Box2D const &intersection, bool crop = true);

/// Generate a new multi-scalar-field with a distribution which is the
///   intersection of the two original distributions.
/** No data is being copied, the new scalar-field is initialized to zero values.
 *    The block identifiers and parallel distribution are taken from
 *    originalField1. The new bounding-box is equal to the intersection
 *    of the original ones if the parameter "crop" is true, and equal to
 *    the bound around the original ones if "crop" is false.
 **/
template <typename T>
std::unique_ptr<MultiScalarField2D<T> > generateIntersectMultiScalarField(
    MultiBlock2D const &originalField1, MultiBlock2D const &originalField2, bool crop = true);

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
std::unique_ptr<MultiScalarField2D<T> > generateIntersectMultiScalarField(
    MultiBlock2D const &originalField1, MultiBlock2D const &originalField2,
    Box2D const &intersection, bool crop = true);

/// Generate a new multi-scalar-field with a distribution which is the
///   union of the two original distributions.
/** No data is being copied, the new scalar-field is initialized to zero values.
 *    The block identifiers and parallel distributions are taken from
 *    originalField1 wherever originalField1 is defined. Elsewhere,
 *    new identifiers are generated, and the parallel disribution is
 *    taken from originalField2.
 **/
template <typename T>
std::unique_ptr<MultiScalarField2D<T> > generateJoinMultiScalarField(
    MultiBlock2D const &originalField1, MultiBlock2D const &originalField2);

/// Generate a new multi-scalar-field with a distribution which is the
///   union between the original one and an additional block.
template <typename T>
std::unique_ptr<MultiScalarField2D<T> > extend(
    MultiScalarField2D<T> &originalBlock, Box2D const &addedBlock);

/// Generate a new multi-scalar-field with a distribution which is the
///   union between the original one with exception of the specified block.
template <typename T>
std::unique_ptr<MultiScalarField2D<T> > except(
    MultiScalarField2D<T> &originalBlock, Box2D const &exceptedBlock);

/// Create a clone of the original field with a different block-distribution.
template <typename T>
std::unique_ptr<MultiScalarField2D<T> > redistribute(
    MultiScalarField2D<T> const &originalField, SparseBlockStructure2D const &newBlockStructure,
    bool adjustPeriodicity = true);

/// Create a clone of the original field on the domain of intersection,
///   with a different block-distribution.
/**   The new boundingBox is either equal to the original one, or equal to the
 *    intersection, depending on the parameter crop.
 **/
template <typename T>
std::unique_ptr<MultiScalarField2D<T> > redistribute(
    MultiScalarField2D<T> const &originalField, SparseBlockStructure2D const &newBlockStructure,
    Box2D const &intersection, bool crop = true);

/// Create a clone of the original field, re-distributed in such a way
///   that it is aligned with the partner-block, i.e. it can be used
///   together with partner-block from within a data processor.
template <typename T>
std::unique_ptr<MultiScalarField2D<T> > align(
    MultiScalarField2D<T> const &originalBlock, MultiBlock2D const &partnerBlock);

/// Create a clone of the original field, on a new regular distribution.
template <typename T>
std::unique_ptr<MultiScalarField2D<T> > reparallelize(MultiScalarField2D<T> const &originalBlock);

/// Create a clone of the original field, on a new regular distribution.
/** The parameters blockLx and blockLy indicate the approximate size of the blocks.
 **/
template <typename T>
std::unique_ptr<MultiScalarField2D<T> > reparallelize(
    MultiScalarField2D<T> const &originalBlock, plint blockLx, plint blockLy);

/* *************** 2. MultiNTensorField ************************************** */

/// Generate a multi-scalar-field from scratch. As opposed to the standard
///   constructor, this factory function takes the explicit block-management
///   object, which includes stuff like block-distribution, parallelization,
///   envelope-width, etc.
template <typename T>
std::unique_ptr<MultiNTensorField2D<T> > defaultGenerateMultiNTensorField2D(
    MultiBlockManagement2D const &management, plint nDim = 1);

template <typename T>
MultiNTensorField2D<T> *generateMultiNTensorField2D(
    MultiBlock2D &multiBlock, plint envelopeWidth, plint ndim);

template <typename T>
MultiNTensorField2D<T> *generateMultiNTensorField2D(Box2D const &domain, plint ndim);

template <typename T>
MultiNTensorField2D<T> *generateMultiNTensorField2D(
    Box2D const &domain, plint ndim, T *iniVal, plint envelopeWidth);

/// Create a clone of a MultiNTensorField (or of a sub-domain).
/** This cannot be handled through a data processor, because the internal data
 *  processors of the ntensor-field must be cloned as well.
 **/
template <typename T>
MultiNTensorField2D<T> *clone(
    MultiNTensorField2D<T> &originalField, Box2D const &subDomain, bool crop = true);

/// Generate a new multi-tensor-field with same distribution as original block,
///   but intersected with a given domain.
/** No data is being copied, the new tensor-field is initialized to zero values.
 *    The new boundingBox is either equal to the original one, or equal to the
 *    intersection, depending on the parameter crop.
 **/
template <typename T>
MultiNTensorField2D<T> *generateMultiNTensorField(
    MultiBlock2D const &originalField, Box2D const &intersection, plint nDim, bool crop = true);

template <typename T1, typename T2>
MultiNTensorField2D<T2> *generateNTensorFieldFromNTensor2D(
    MultiNTensorField2D<T1> const &field, Box2D const &intersection, plint nDim);

template <typename T1, typename T2, template <typename U> class Descriptor>
MultiNTensorField2D<T1> *generateNTensorFieldFromBlockLattice2D(
    MultiBlockLattice2D<T2, Descriptor> const &lattice, Box2D const &intersection, plint nDim);

/// Generate a new multi-tensor-field with a distribution which is the
///   intersection of the two original distributions.
/** No data is being copied, the new tensor-field is initialized to zero values.
 *    The block identifiers and parallel distribution are taken from
 *    originalField1. The new bounding-box is equal to the intersection
 *    of the original ones if the parameter "crop" is true, and equal to
 *    the bound around the original ones if "crop" is false.
 **/
template <typename T>
MultiNTensorField2D<T> *generateIntersectMultiNTensorField(
    MultiBlock2D const &originalField1, MultiBlock2D const &originalField2, plint nDim,
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
MultiNTensorField2D<T> *generateIntersectMultiNTensorField(
    MultiBlock2D const &originalField1, MultiBlock2D const &originalField2,
    Box2D const &intersection, plint nDim, bool crop = true);

/// Generate a new multi-tensor-field with a distribution which is the
///   union of the two original distributions.
/** No data is being copied, the new tensor-field is initialized to zero values.
 *    The block identifiers and parallel distributions are taken from
 *    originalField1 wherever originalField1 is defined. Elsewhere,
 *    new identifiers are generated, and the parallel disribution is
 *    taken from originalField2.
 **/
template <typename T>
MultiNTensorField2D<T> *generateJoinMultiNTensorField(
    MultiBlock2D const &originalField1, MultiBlock2D const &originalField2, plint nDim);

/// Generate a new multi-ntensor-field with a distribution which is the
///   union between the original one and an additional block.
template <typename T>
MultiNTensorField2D<T> *extend(MultiNTensorField2D<T> &originalBlock, Box2D const &addedBlock);

/// Generate a new multi-ntensor-field with a distribution which is the
///   union between the original one with exception of the specified block.
template <typename T>
MultiNTensorField2D<T> *except(MultiNTensorField2D<T> &originalBlock, Box2D const &exceptedBlock);

/// Create a clone of the original field, re-distributed in such a way
///   that it is aligned with the partner-block, i.e. it can be used
///   together with partner-block from within a data processor.
template <typename T>
MultiNTensorField2D<T> *align(
    MultiNTensorField2D<T> const &originalBlock, MultiBlock2D const &partnerBlock);

/// Create a clone of the original field, on a new regular distribution.
template <typename T>
MultiNTensorField2D<T> *reparallelize(MultiNTensorField2D<T> const &originalBlock);

/* *************** 3. MultiTensorField ************************************** */

/// Generate a multi-tensor-field from scratch. As opposed to the standard
///   constructor, this factory function takes a full bounding-box, as well
///   as the envelope-width, as arguments.
template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > generateMultiTensorField(
    Box2D boundingBox, plint envelopeWidth = 1);

/// Generate a multi-tensor-field from scratch. As opposed to the standard
///   constructor, this factory function takes a full bounding-box, as well
///   as the envelope-width, as arguments.
template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > generateMultiTensorField(
    Box2D boundingBox, Array<T, nDim> const &iniVal, plint envelopeWidth = 1);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > generateMultiTensorField(
    MultiBlock2D &multiBlock, plint envelopeWidth);

/// Generate a multi-tensor-field from scratch. As opposed to the standard
///   constructor, this factory function takes the explicit block-management
///   object, which includes stuff like block-distribution, parallelization,
///   envelope-width, etc. A default dummy argument is used for technical reasons.
template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > defaultGenerateMultiTensorField2D(
    MultiBlockManagement2D const &management, plint unnamedDummyArg = 1);

/// Create a clone of a MultiTensorField (or of a sub-domain).
/** This cannot be handled through a data processor, because the internal data
 *  processors of the tensor-field must be cloned as well.
 **/
template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > clone(
    MultiTensorField2D<T, nDim> &originalField, Box2D const &subDomain, bool crop = true);

/// Generate a new multi-tensor-field with same distribution as original block,
///   but intersected with a given domain.
/** No data is being copied, the new tensor-field is initialized to zero values.
 *    The new boundingBox is either equal to the original one, or equal to the
 *    intersection, depending on the parameter crop.
 **/
template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > generateMultiTensorField(
    MultiBlock2D const &originalField, Box2D const &intersection, bool crop = true);

/// Generate a new multi-tensor-field with a distribution which is the
///   intersection of the two original distributions.
/** No data is being copied, the new tensor-field is initialized to zero values.
 *    The block identifiers and parallel distribution are taken from
 *    originalField1. The new bounding-box is equal to the intersection
 *    of the original ones if the parameter "crop" is true, and equal to
 *    the bound around the original ones if "crop" is false.
 **/
template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > generateIntersectMultiTensorField(
    MultiBlock2D const &originalField1, MultiBlock2D const &originalField2, bool crop = true);

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
std::unique_ptr<MultiTensorField2D<T, nDim> > generateIntersectMultiTensorField(
    MultiBlock2D const &originalField1, MultiBlock2D const &originalField2,
    Box2D const &intersection, bool crop = true);

/// Generate a new multi-tensor-field with a distribution which is the
///   union of the two original distributions.
/** No data is being copied, the new tensor-field is initialized to zero values.
 *    The block identifiers and parallel distributions are taken from
 *    originalField1 wherever originalField1 is defined. Elsewhere,
 *    new identifiers are generated, and the parallel disribution is
 *    taken from originalField2.
 **/
template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > generateJoinMultiTensorField(
    MultiBlock2D const &originalField1, MultiBlock2D const &originalField2);

/// Generate a new multi-tensor-field with a distribution which is the
///   union between the original one and an additional block.
template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > extend(
    MultiTensorField2D<T, nDim> &originalBlock, Box2D const &addedBlock);

/// Generate a new multi-tensor-field with a distribution which is the
///   union between the original one with exception of the specified block.
template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > except(
    MultiTensorField2D<T, nDim> &originalBlock, Box2D const &exceptedBlock);

/// Create a clone of the original field with a different block-distribution.
template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > redistribute(
    MultiTensorField2D<T, nDim> const &originalField,
    SparseBlockStructure2D const &newBlockStructure, bool adjustPeriodicity = true);

/// Create a clone of the original field on the domain of intersection,
///   with a different block-distribution.
/**   The new boundingBox is either equal to the original one, or equal to the
 *    intersection, depending on the parameter crop.
 **/
template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > redistribute(
    MultiTensorField2D<T, nDim> const &originalField,
    SparseBlockStructure2D const &newBlockStructure, Box2D const &intersection, bool crop = true);

/// Create a clone of the original field, re-distributed in such a way
///   that it is aligned with the partner-block, i.e. it can be used
///   together with partner-block from within a data processor.
template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > align(
    MultiTensorField2D<T, nDim> const &originalBlock, MultiBlock2D const &partnerBlock);

/// Create a clone of the original field, on a new regular distribution.
template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > reparallelize(
    MultiTensorField2D<T, nDim> const &originalBlock);

/// Create a clone of the original field, on a new regular distribution.
/** The parameters blockLx and blockLy indicate the approximate size of the blocks.
 **/
template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > reparallelize(
    MultiTensorField2D<T, nDim> const &originalBlock, plint blockLx, plint blockLy);

/* *************** 4. MultiBlockLattice ************************************** */

/// Generate a multi-block-lattice from scratch. As opposed to the standard
///   constructor, this factory function takes a full bounding-box, as well
///   as the envelope-width, as arguments.
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > generateMultiBlockLattice(
    Box2D boundingBox, Dynamics<T, Descriptor> *backgroundDynamics, plint envelopeWidth = 1);

/// Generate a multi-block-lattice from scratch. As opposed to the standard
///   constructor, this factory function takes the explicit block-management
///   object, which includes stuff like block-distribution, parallelization,
///   envelope-width, etc. A default dummy argument is used for technical reasons.
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > defaultGenerateMultiBlockLattice2D(
    MultiBlockManagement2D const &management, plint unnamedDummyArg = 1);

/// Create a clone of a MultiBlockLattice (or of a sub-domain).
/** This cannot be handled through a data processor, because the internal data
 *  processors of the block-lattice must be cloned as well.
 **/
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > clone(
    MultiBlockLattice2D<T, Descriptor> &originalLattice, Box2D const &subDomain, bool crop = true);

/// Generate a new multi-block-lattice with same distribution as original block,
///   but intersected with a given domain.
/** No data is being copied, the new block-lattice is initialized to default values.
 *    The new boundingBox is either equal to the original one, or equal to the
 *    intersection, depending on the parameter crop.
 **/
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > generateMultiBlockLattice(
    MultiBlock2D const &originalBlock, Box2D const &intersection, bool crop = true);

/// Generate a new multi-block-lattice with a distribution which is the
///   intersection of the two original distributions.
/** No data is being copied, the new block-lattice is initialized to default values.
 *    The block identifiers and parallel distribution are taken from
 *    originalBlock1. The new bounding-box is equal to the intersection
 *    of the original ones if the parameter "crop" is true, and equal to
 *    the bound around the original ones if "crop" is false.
 **/
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > generateIntersectMultiBlockLattice(
    MultiBlock2D const &originalBlock1, MultiBlock2D const &originalBlock2, bool crop = true);

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
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > generateIntersectMultiBlockLattice(
    MultiBlock2D const &originalBlock1, MultiBlock2D const &originalBlock2,
    Box2D const &intersection, bool crop = true);

/// Generate a new multi-block-lattice with a distribution which is the
///   union of the two original distributions.
/** No data is being copied, the new block-lattice is initialized to default values.
 *    The block identifiers and parallel distributions are taken from
 *    originalBlock1 wherever originalBlock1 is defined. Elsewhere,
 *    new identifiers are generated, and the parallel disribution is
 *    taken from originalBlock2.
 **/
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > generateJoinMultiBlockLattice(
    MultiBlock2D const &originalBlock1, MultiBlock2D const &originalBlock2);

/// Generate a new multi-block-lattice with a distribution which is the
///   union between the original one and an additional block.
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > extend(
    MultiBlockLattice2D<T, Descriptor> &originalBlock, Box2D const &addedBlock);

/// Generate a new multi-block-lattice with a distribution which is the
///   union between the original one with exception of the specified block.
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > except(
    MultiBlockLattice2D<T, Descriptor> &originalBlock, Box2D const &exceptedBlock);

/// Create a clone of the original lattice with a different block-distribution.
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > redistribute(
    MultiBlockLattice2D<T, Descriptor> const &originalBlock,
    SparseBlockStructure2D const &newBlockStructure, bool adjustPeriodicity = true);

/// Create a clone of the original lattice on the domain of intersection,
///   with a different block-distribution.
/**   The new boundingBox is either equal to the original one, or equal to the
 *    intersection, depending on the parameter crop.
 **/
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > redistribute(
    MultiBlockLattice2D<T, Descriptor> const &originalBlock,
    SparseBlockStructure2D const &newBlockStructure, Box2D const &intersection, bool crop = true);

/// Create a clone of the original lattice, re-distributed in such a way
///   that it is aligned with the partner-block, i.e. it can be used
///   together with partner-block from within a data processor.
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > align(
    MultiBlockLattice2D<T, Descriptor> const &originalBlock, MultiBlock2D const &partnerBlock);

/// Create a clone of the original lattice, on a new regular distribution.
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > reparallelize(
    MultiBlockLattice2D<T, Descriptor> const &originalBlock);

/// Create a clone of the original lattice, on a new regular distribution.
/** The parameters blockLx and blockLy indicate the approximate size of the blocks.
 **/
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > reparallelize(
    MultiBlockLattice2D<T, Descriptor> const &originalBlock, plint blockLx, plint blockLy);

/* *************** 4. MultiParticleField ************************************ */

/// Generate a multi-particle-field from scratch. As opposed to the standard
///   constructor, this factory function takes a full bounding-box, as well
///   as the envelope-width, as arguments.
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiParticleField2D<DenseParticleField2D<T, Descriptor> > >
    generateMultiDenseParticleField(Box2D boundingBox, plint envelopeWidth = 1);

template <class ParticleFieldT>
std::unique_ptr<MultiParticleField2D<ParticleFieldT> > generateMultiParticleField2D(
    Box2D boundingBox, plint envelopeWidth = 1);

template <class ParticleFieldT>
std::unique_ptr<MultiParticleField2D<ParticleFieldT> > generateMultiParticleField2D(
    MultiBlock2D &multiBlock, plint envelopeWidth);

/* *************** 5. MultiContainerBlock ************************************ */

std::unique_ptr<MultiContainerBlock2D> generateMultiContainerBlock(
    MultiBlock2D &multiBlock, plint envelopeWidth);

MultiContainerBlock2D *createMultiContainerBlock2D(
    MultiBlockManagement2D const &management, PeriodicitySwitch2D &periodicity, plint envelopeWidth,
    plint gridLevel);

/* *************** General Functions **************************************** */

/// Copy all data processors from the old to the new block, and replace
///   the "self" block in the list of the data processor to point to the new block.
void transferDataProcessors(MultiBlock2D const &from, MultiBlock2D &to);

}  // namespace plb

#endif  // MULTI_BLOCK_GENERATOR_2D_H
