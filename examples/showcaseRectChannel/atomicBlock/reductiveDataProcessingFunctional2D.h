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
 * Utilities to help users handle data processors -- header file.
 */
#ifndef REDUCTIVE_DATA_PROCESSING_FUNCTIONAL_2D_H
#define REDUCTIVE_DATA_PROCESSING_FUNCTIONAL_2D_H

#include <vector>

#include "atomicBlock/dataProcessor2D.h"
#include "core/blockSurface2D.h"
#include "core/geometry2D.h"
#include "core/globalDefs.h"

namespace plb {

// Forward declarations.

template <typename T, template <typename U> class Descriptor>
class BlockLatticeBase2D;
template <typename T, template <typename U> class Descriptor>
class BlockLattice2D;
template <typename T>
class ScalarFieldBase2D;
template <typename T>
class ScalarField2D;
template <typename T, int nDim>
class TensorFieldBase2D;
template <typename T, int nDim>
class TensorField2D;
template <typename T>
class NTensorFieldBase2D;
template <typename T>
class NTensorField2D;

/* *************** All flavors of Box processing functionals ********* */

/// Easy instantiation of reductive boxed data processor (general case)
struct ReductiveBoxProcessingFunctional2D {
    virtual ~ReductiveBoxProcessingFunctional2D() { }
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) = 0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(double dxScale, double dtScale);
    virtual void getDimensionsX(std::vector<int> &dimensions) const;
    virtual void getDimensionsT(std::vector<int> &dimensions) const;
    /// Obsolete: replaced by getTypeOfModification. Don't override.
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const = 0;
    virtual ReductiveBoxProcessingFunctional2D *clone() const = 0;
    virtual BlockStatistics const &getStatistics() const = 0;
    virtual BlockStatistics &getStatistics() = 0;
    virtual void serialize(std::string &data) const;
    virtual void unserialize(std::string &data);
    virtual int getStaticId() const;
};

/// ReductiveBoxProcessingFunctional2D which instantiates its own statistics object
class PlainReductiveBoxProcessingFunctional2D : public ReductiveBoxProcessingFunctional2D {
public:
    virtual BlockStatistics const &getStatistics() const;
    virtual BlockStatistics &getStatistics();
    void swap(PlainReductiveBoxProcessingFunctional2D &rhs)
    {
        std::swap(statistics, rhs.statistics);
    }

private:
    BlockStatistics statistics;
};

/// A reductive boxed data processor, automatically generated from a
/// ReductiveBoxProcessingFunctional2D
class ReductiveBoxProcessor2D : public DataProcessor2D {
public:
    /** \param functional_ The functional is not owned by the ReductiveBoxProcessor2D,
     *                     i.e. it is not deleted in the destructor.
     */
    ReductiveBoxProcessor2D(
        ReductiveBoxProcessingFunctional2D *functional_, Box2D domain_,
        std::vector<AtomicBlock2D *> atomicBlocks_);
    Box2D getDomain() const;
    virtual void process();
    virtual ReductiveBoxProcessor2D *clone() const;
    virtual int getStaticId() const;

private:
    ReductiveBoxProcessingFunctional2D *functional;
    Box2D domain;
    std::vector<AtomicBlock2D *> atomicBlocks;
};

/// An automatically created generator for the ReductiveBoxProcessor2D
class ReductiveBoxProcessorGenerator2D : public BoxedReductiveDataProcessorGenerator2D {
public:
    /** \param functional_ The functional is not owned by the ReductiveBoxProcessorGenerator2D,
     *                     i.e. it is not deleted in the destructor.
     */
    ReductiveBoxProcessorGenerator2D(ReductiveBoxProcessingFunctional2D *functional_, Box2D domain);
    ~ReductiveBoxProcessorGenerator2D();
    ReductiveBoxProcessorGenerator2D(ReductiveBoxProcessorGenerator2D const &rhs);
    ReductiveBoxProcessorGenerator2D &operator=(ReductiveBoxProcessorGenerator2D const &rhs);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(double dxScale, double dtScale);
    virtual void getDimensionsX(std::vector<int> &dimensions) const;
    virtual void getDimensionsT(std::vector<int> &dimensions) const;
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual DataProcessor2D *generate(std::vector<AtomicBlock2D *> atomicBlocks);
    virtual ReductiveBoxProcessorGenerator2D *clone() const;
    virtual BlockStatistics const &getStatistics() const;
    virtual BlockStatistics &getStatistics();
    ReductiveBoxProcessingFunctional2D const &getFunctional() const;
    virtual void serialize(Box2D &domain, std::string &data) const;
    virtual int getStaticId() const;

private:
    ReductiveBoxProcessingFunctional2D *functional;
};

/// Easy instantiation of boxed data processor for a single lattice
template <typename T, template <typename U> class Descriptor>
struct ReductiveBoxProcessingFunctional2D_L : public PlainReductiveBoxProcessingFunctional2D {
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single scalar field
template <typename T>
struct ReductiveBoxProcessingFunctional2D_S : public PlainReductiveBoxProcessingFunctional2D {
    virtual void process(Box2D domain, ScalarField2D<T> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single tensor field
template <typename T, int nDim>
struct ReductiveBoxProcessingFunctional2D_T : public PlainReductiveBoxProcessingFunctional2D {
    virtual void process(Box2D domain, TensorField2D<T, nDim> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single n-tensor field
template <typename T>
struct ReductiveBoxProcessingFunctional2D_N : public PlainReductiveBoxProcessingFunctional2D {
    virtual void process(Box2D domain, NTensorField2D<T> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single n-tensor field
template <typename T>
struct MaskedReductiveBoxProcessingFunctional2D_N : public PlainReductiveBoxProcessingFunctional2D {
    virtual void process(Box2D domain, NTensorField2D<T> &field, NTensorField2D<int> &mask) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for lattice-lattice coupling
template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
struct ReductiveBoxProcessingFunctional2D_LL : public PlainReductiveBoxProcessingFunctional2D {
    virtual void process(
        Box2D domain, BlockLattice2D<T1, Descriptor1> &lattice1,
        BlockLattice2D<T2, Descriptor2> &lattice2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for ScalarField-ScalarField coupling
template <typename T1, typename T2>
struct ReductiveBoxProcessingFunctional2D_SS : public PlainReductiveBoxProcessingFunctional2D {
    virtual void process(Box2D domain, ScalarField2D<T1> &field1, ScalarField2D<T2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for TensorField-TensorField coupling
template <typename T1, int nDim1, typename T2, int nDim2>
struct ReductiveBoxProcessingFunctional2D_TT : public PlainReductiveBoxProcessingFunctional2D {
    virtual void process(
        Box2D domain, TensorField2D<T1, nDim1> &field1, TensorField2D<T2, nDim2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for ScalarField-TensorField coupling
template <typename T1, typename T2, int nDim>
struct ReductiveBoxProcessingFunctional2D_ST : public PlainReductiveBoxProcessingFunctional2D {
    virtual void process(
        Box2D domain, ScalarField2D<T1> &field1, TensorField2D<T2, nDim> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for NTensorField-NTensorField coupling
template <typename T1, typename T2>
struct ReductiveBoxProcessingFunctional2D_NN : public PlainReductiveBoxProcessingFunctional2D {
    virtual void process(Box2D domain, NTensorField2D<T1> &field1, NTensorField2D<T2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for NTensorField-NTensorField coupling
template <typename T1, typename T2>
struct MaskedReductiveBoxProcessingFunctional2D_NN :
    public PlainReductiveBoxProcessingFunctional2D {
    virtual void process(
        Box2D domain, NTensorField2D<T1> &field1, NTensorField2D<T2> &field2,
        NTensorField2D<int> &mask) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for ScalarField-NTensorField coupling
template <typename T1, typename T2>
struct ReductiveBoxProcessingFunctional2D_SN : public PlainReductiveBoxProcessingFunctional2D {
    virtual void process(Box2D domain, ScalarField2D<T1> &field1, NTensorField2D<T2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-ScalarField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct ReductiveBoxProcessingFunctional2D_LS : public PlainReductiveBoxProcessingFunctional2D {
    virtual void process(
        Box2D domain, BlockLattice2D<T1, Descriptor> &lattice, ScalarField2D<T2> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-NTensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
struct ReductiveBoxProcessingFunctional2D_LT : public PlainReductiveBoxProcessingFunctional2D {
    virtual void process(
        Box2D domain, BlockLattice2D<T1, Descriptor> &lattice, TensorField2D<T2, nDim> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-NTensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct ReductiveBoxProcessingFunctional2D_LN : public PlainReductiveBoxProcessingFunctional2D {
    virtual void process(
        Box2D domain, BlockLattice2D<T1, Descriptor> &lattice, NTensorField2D<T2> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-NTensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct MaskedReductiveBoxProcessingFunctional2D_LN :
    public PlainReductiveBoxProcessingFunctional2D {
    virtual void process(
        Box2D domain, BlockLattice2D<T1, Descriptor> &lattice, NTensorField2D<T2> &field,
        NTensorField2D<int> &mask) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple lattices
template <typename T, template <typename U> class Descriptor>
struct ReductiveLatticeBoxProcessingFunctional2D : public PlainReductiveBoxProcessingFunctional2D {
    virtual void process(Box2D domain, std::vector<BlockLattice2D<T, Descriptor> *> lattices) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple scalarFields
template <typename T>
struct ReductiveScalarFieldBoxProcessingFunctional2D :
    public PlainReductiveBoxProcessingFunctional2D {
    virtual void process(Box2D domain, std::vector<ScalarField2D<T> *> scalarFields) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple tensorFields
template <typename T, int nDim>
struct ReductiveTensorFieldBoxProcessingFunctional2D :
    public PlainReductiveBoxProcessingFunctional2D {
    virtual void process(Box2D domain, std::vector<TensorField2D<T, nDim> *> tensorFields) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple n-tensorFields
template <typename T>
struct ReductiveNTensorFieldBoxProcessingFunctional2D :
    public PlainReductiveBoxProcessingFunctional2D {
    virtual void process(Box2D domain, std::vector<NTensorField2D<T> *> tensorFields) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/* *************** All flavors of Dot processing functionals ********* */

/// Easy instantiation of reductive boxed data processor (general case)
class ReductiveDotProcessingFunctional2D {
public:
    virtual ~ReductiveDotProcessingFunctional2D() { }
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks) = 0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(double dxScale, double dtScale);
    virtual void getDimensionsX(std::vector<int> &dimensions) const;
    virtual void getDimensionsT(std::vector<int> &dimensions) const;
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const = 0;
    virtual ReductiveDotProcessingFunctional2D *clone() const = 0;
    virtual BlockStatistics const &getStatistics() const = 0;
    virtual BlockStatistics &getStatistics() = 0;
};

/// ReductiveDotProcessingFunctional2D which instantiates its own statistics object
class PlainReductiveDotProcessingFunctional2D : public ReductiveDotProcessingFunctional2D {
public:
    virtual BlockStatistics const &getStatistics() const;
    virtual BlockStatistics &getStatistics();

private:
    BlockStatistics statistics;
};

/// A ReductiveDotted data processor, automatically generated from a
/// ReductiveDotProcessingFunctional2D
class ReductiveDotProcessor2D : public DataProcessor2D {
public:
    /** \param functional_ The functional is not owned by the ReductiveBoxProcessor2D,
     *                     i.e. it is not deleted in the destructor.
     */
    ReductiveDotProcessor2D(
        ReductiveDotProcessingFunctional2D *functional_, DotList2D const &dotList_,
        std::vector<AtomicBlock2D *> atomicBlocks_);
    DotList2D const &getDotList() const;
    virtual void process();
    virtual ReductiveDotProcessor2D *clone() const;

private:
    ReductiveDotProcessingFunctional2D *functional;
    DotList2D dotList;
    std::vector<AtomicBlock2D *> atomicBlocks;
};

/// An automatically created generator for the ReductiveDotProcessor2D
class ReductiveDotProcessorGenerator2D : public DottedReductiveDataProcessorGenerator2D {
public:
    ReductiveDotProcessorGenerator2D(
        ReductiveDotProcessingFunctional2D *functional_, DotList2D const &dotList);
    ~ReductiveDotProcessorGenerator2D();
    ReductiveDotProcessorGenerator2D(ReductiveDotProcessorGenerator2D const &rhs);
    ReductiveDotProcessorGenerator2D &operator=(ReductiveDotProcessorGenerator2D const &rhs);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(double dxScale, double dtScale);
    virtual void getDimensionsX(std::vector<int> &dimensions) const;
    virtual void getDimensionsT(std::vector<int> &dimensions) const;
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual DataProcessor2D *generate(std::vector<AtomicBlock2D *> atomicBlocks);
    virtual ReductiveDotProcessorGenerator2D *clone() const;
    virtual BlockStatistics const &getStatistics() const;
    virtual BlockStatistics &getStatistics();
    ReductiveDotProcessingFunctional2D const &getFunctional() const;

private:
    ReductiveDotProcessingFunctional2D *functional;
};

/// Easy instantiation of dotted data processor for a single lattice
template <typename T, template <typename U> class Descriptor>
struct ReductiveDotProcessingFunctional2D_L : public PlainReductiveDotProcessingFunctional2D {
    virtual void process(DotList2D const &dotList, BlockLattice2D<T, Descriptor> &lattice) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for a single scalar field
template <typename T>
struct ReductiveDotProcessingFunctional2D_S : public PlainReductiveDotProcessingFunctional2D {
    virtual void process(DotList2D const &dotList, ScalarField2D<T> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for a single tensor field
template <typename T, int nDim>
struct ReductiveDotProcessingFunctional2D_T : public PlainReductiveDotProcessingFunctional2D {
    virtual void process(DotList2D const &dotList, TensorField2D<T, nDim> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for a single n-tensor field
template <typename T>
struct ReductiveDotProcessingFunctional2D_N : public PlainReductiveDotProcessingFunctional2D {
    virtual void process(DotList2D const &dotList, NTensorField2D<T> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for lattice-lattice coupling
template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
struct ReductiveDotProcessingFunctional2D_LL : public PlainReductiveDotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, BlockLattice2D<T1, Descriptor1> &lattice1,
        BlockLattice2D<T2, Descriptor2> &lattice2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for ScalarField-ScalarField coupling
template <typename T1, typename T2>
struct ReductiveDotProcessingFunctional2D_SS : public PlainReductiveDotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, ScalarField2D<T1> &field1, ScalarField2D<T2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for TensorField-TensorField coupling
template <typename T1, int nDim1, typename T2, int nDim2>
struct ReductiveDotProcessingFunctional2D_TT : public PlainReductiveDotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, TensorField2D<T1, nDim1> &field1,
        TensorField2D<T2, nDim2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for ScalarField-TensorField coupling
template <typename T1, typename T2, int nDim>
struct ReductiveDotProcessingFunctional2D_ST : public PlainReductiveDotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, ScalarField2D<T1> &field1, TensorField2D<T2, nDim> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for NTensorField-NTensorField coupling
template <typename T1, typename T2>
struct ReductiveDotProcessingFunctional2D_NN : public PlainReductiveDotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, NTensorField2D<T1> &field1, NTensorField2D<T2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for ScalarField-NTensorField coupling
template <typename T1, typename T2>
struct ReductiveDotProcessingFunctional2D_SN : public PlainReductiveDotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, ScalarField2D<T1> &field1, NTensorField2D<T2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for Lattice-ScalarField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct ReductiveDotProcessingFunctional2D_LS : public PlainReductiveDotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, BlockLattice2D<T1, Descriptor> &lattice,
        ScalarField2D<T2> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for Lattice-TensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
struct ReductiveDotProcessingFunctional2D_LT : public PlainReductiveDotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, BlockLattice2D<T1, Descriptor> &lattice,
        TensorField2D<T2, nDim> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for Lattice-NTensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct ReductiveDotProcessingFunctional2D_LN : public PlainReductiveDotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, BlockLattice2D<T1, Descriptor> &lattice,
        NTensorField2D<T2> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple lattices
template <typename T, template <typename U> class Descriptor>
struct ReductiveLatticeDotProcessingFunctional2D : public PlainReductiveDotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, std::vector<BlockLattice2D<T, Descriptor> *> lattices) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple scalarFields
template <typename T>
struct ReductiveScalarFieldDotProcessingFunctional2D :
    public PlainReductiveDotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, std::vector<ScalarField2D<T> *> scalarFields) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple tensorFields
template <typename T, int nDim>
struct ReductiveTensorFieldDotProcessingFunctional2D :
    public PlainReductiveDotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, std::vector<TensorField2D<T, nDim> *> tensorFields) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple n-tensorFields
template <typename T>
struct ReductiveNTensorFieldDotProcessingFunctional2D :
    public PlainReductiveDotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, std::vector<NTensorField2D<T> *> tensorFields) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/* *************** All flavors of Bounded Box processing functionals ********* */

/// Easy instantiation of boxed reductive data processor special boundary treatment (general case)
class BoundedReductiveBoxProcessingFunctional2D {
public:
    virtual ~BoundedReductiveBoxProcessingFunctional2D() { }
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) = 0;
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain,
        std::vector<AtomicBlock2D *> atomicBlocks) = 0;
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) = 0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(double dxScale, double dtScale);
    virtual void getDimensionsX(std::vector<int> &dimensions) const;
    virtual void getDimensionsT(std::vector<int> &dimensions) const;
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const = 0;
    virtual BoundedReductiveBoxProcessingFunctional2D *clone() const = 0;
    ReductiveBoxProcessingFunctional2D *getBulkProcessor() const;
    ReductiveBoxProcessingFunctional2D *getEdgeProcessor(int direction, int orientation) const;
    ReductiveBoxProcessingFunctional2D *getCornerProcessor(int normalX, int normalY) const;
    BlockStatistics const &getStatistics() const;
    BlockStatistics &getStatistics();
    void getGenerators(
        Box2D const &fullDomain, plint boundaryWidth_,
        std::vector<ReductiveBoxProcessorGenerator2D *> &generators);

private:
    BlockStatistics statistics;

public:
    class BulkWrapperFunctional : public ReductiveBoxProcessingFunctional2D {
    public:
        BulkWrapperFunctional(BoundedReductiveBoxProcessingFunctional2D *boundedFunctional_);
        BulkWrapperFunctional(BulkWrapperFunctional const &rhs);
        ~BulkWrapperFunctional();
        BulkWrapperFunctional &operator=(BulkWrapperFunctional const &rhs);
        virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        virtual void rescale(double dxScale, double dtScale);
        virtual void getDimensionsX(std::vector<int> &dimensions) const;
        virtual void getDimensionsT(std::vector<int> &dimensions) const;
        virtual void getModificationPattern(std::vector<bool> &isWritten) const;
        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
        virtual BulkWrapperFunctional *clone() const;
        virtual BlockStatistics const &getStatistics() const;
        virtual BlockStatistics &getStatistics();

    private:
        BoundedReductiveBoxProcessingFunctional2D *boundedFunctional;
    };
    class EdgeWrapperFunctional : public ReductiveBoxProcessingFunctional2D {
    public:
        EdgeWrapperFunctional(
            BoundedReductiveBoxProcessingFunctional2D *boundedFunctional_, int direction_,
            int orientation_);
        EdgeWrapperFunctional(EdgeWrapperFunctional const &rhs);
        ~EdgeWrapperFunctional();
        EdgeWrapperFunctional &operator=(EdgeWrapperFunctional const &rhs);
        virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        virtual void rescale(double dxScale, double dtScale);
        virtual void getDimensionsX(std::vector<int> &dimensions) const;
        virtual void getDimensionsT(std::vector<int> &dimensions) const;
        virtual void getModificationPattern(std::vector<bool> &isWritten) const;
        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
        virtual EdgeWrapperFunctional *clone() const;
        virtual BlockStatistics const &getStatistics() const;
        virtual BlockStatistics &getStatistics();

    private:
        BoundedReductiveBoxProcessingFunctional2D *boundedFunctional;
        int direction, orientation;
    };
    class CornerWrapperFunctional : public ReductiveBoxProcessingFunctional2D {
    public:
        CornerWrapperFunctional(
            BoundedReductiveBoxProcessingFunctional2D *boundedFunctional_, int normalX_,
            int normalY_);
        CornerWrapperFunctional(CornerWrapperFunctional const &rhs);
        ~CornerWrapperFunctional();
        CornerWrapperFunctional &operator=(CornerWrapperFunctional const &rhs);
        virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        virtual void rescale(double dxScale, double dtScale);
        virtual void getDimensionsX(std::vector<int> &dimensions) const;
        virtual void getDimensionsT(std::vector<int> &dimensions) const;
        virtual void getModificationPattern(std::vector<bool> &isWritten) const;
        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
        virtual CornerWrapperFunctional *clone() const;
        virtual BlockStatistics const &getStatistics() const;
        virtual BlockStatistics &getStatistics();

    private:
        BoundedReductiveBoxProcessingFunctional2D *boundedFunctional;
        int normalX, normalY;
    };
};

/// Easy instantiation of bounded boxed data processor for a single lattice
template <typename T, template <typename U> class Descriptor>
struct BoundedReductiveBoxProcessingFunctional2D_L :
    public BoundedReductiveBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(Box2D domain, BlockLattice2D<T, Descriptor> &lattice) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain, BlockLattice2D<T, Descriptor> &lattice) = 0;
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, BlockLattice2D<T, Descriptor> &lattice) = 0;
};

/// Easy instantiation of bounded boxed data processor for a single scalar field
template <typename T>
struct BoundedReductiveBoxProcessingFunctional2D_S :
    public BoundedReductiveBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(Box2D domain, ScalarField2D<T> &field) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain, ScalarField2D<T> &field) = 0;
    virtual void processCorner(int normalX, int normalY, Box2D domain, ScalarField2D<T> &field) = 0;
};

/// Easy instantiation of bounded boxed data processor for a single tensor field
template <typename T, int nDim>
struct BoundedReductiveBoxProcessingFunctional2D_T :
    public BoundedReductiveBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(Box2D domain, TensorField2D<T, nDim> &field) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain, TensorField2D<T, nDim> &field) = 0;
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, TensorField2D<T, nDim> &field) = 0;
};

/// Easy instantiation of bounded boxed data processor for a single n-tensor field
template <typename T>
struct BoundedReductiveBoxProcessingFunctional2D_N :
    public BoundedReductiveBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(Box2D domain, NTensorField2D<T> &field) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain, NTensorField2D<T> &field) = 0;
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, NTensorField2D<T> &field) = 0;
};

/// Easy instantiation of bounded boxed data processor for a single n-tensor field
template <typename T>
struct BoundedMaskedReductiveBoxProcessingFunctional2D_N :
    public BoundedReductiveBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(Box2D domain, NTensorField2D<T> &field, NTensorField2D<int> &mask) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain, NTensorField2D<T> &field,
        NTensorField2D<int> &mask) = 0;
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, NTensorField2D<T> &field,
        NTensorField2D<int> &mask) = 0;
};

/// Easy instantiation of bounded boxed data processor for lattice-lattice coupling
template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
struct BoundedReductiveBoxProcessingFunctional2D_LL :
    public BoundedReductiveBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(
        Box2D domain, BlockLattice2D<T1, Descriptor1> &lattice1,
        BlockLattice2D<T2, Descriptor2> &lattice2) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain, BlockLattice2D<T1, Descriptor1> &lattice1,
        BlockLattice2D<T2, Descriptor2> &lattice2) = 0;
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, BlockLattice2D<T1, Descriptor1> &lattice1,
        BlockLattice2D<T2, Descriptor2> &lattice2) = 0;
};

/// Easy instantiation of bounded boxed data processor for
/// ScalarField-ScalarField coupling.
template <typename T1, typename T2>
struct BoundedReductiveBoxProcessingFunctional2D_SS :
    public BoundedReductiveBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(
        Box2D domain, ScalarField2D<T1> &field1, ScalarField2D<T2> &field2) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain, ScalarField2D<T1> &field1,
        ScalarField2D<T2> &field2) = 0;
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, ScalarField2D<T1> &field1,
        ScalarField2D<T2> &field2) = 0;
};

/// Easy instantiation of bounded boxed data processor for
/// TensorField-TensorField coupling.
template <typename T1, int nDim1, typename T2, int nDim2>
struct BoundedReductiveBoxProcessingFunctional2D_TT :
    public BoundedReductiveBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(
        Box2D domain, TensorField2D<T1, nDim1> &field1, TensorField2D<T2, nDim2> &field2) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain, TensorField2D<T1, nDim1> &field1,
        TensorField2D<T2, nDim2> &field2) = 0;
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, TensorField2D<T1, nDim1> &field1,
        TensorField2D<T2, nDim2> &field2) = 0;
};

/// Easy instantiation of bounded boxed data processor for
/// ScalarField-TensorField coupling.
template <typename T1, typename T2, int nDim>
struct BoundedReductiveBoxProcessingFunctional2D_ST :
    public BoundedReductiveBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(
        Box2D domain, ScalarField2D<T1> &field1, TensorField2D<T2, nDim> &field2) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain, ScalarField2D<T1> &field1,
        TensorField2D<T2, nDim> &field2) = 0;
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, ScalarField2D<T1> &field1,
        TensorField2D<T2, nDim> &field2) = 0;
};

/// Easy instantiation of bounded boxed data processor for
/// NTensorField-NTensorField coupling.
template <typename T1, typename T2>
struct BoundedReductiveBoxProcessingFunctional2D_NN :
    public BoundedReductiveBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(
        Box2D domain, NTensorField2D<T1> &field1, NTensorField2D<T2> &field2) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain, NTensorField2D<T1> &field1,
        NTensorField2D<T2> &field2) = 0;
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, NTensorField2D<T1> &field1,
        NTensorField2D<T2> &field2) = 0;
};

/// Easy instantiation of bounded boxed data processor for
/// NTensorField-NTensorField coupling with mask.
template <typename T1, typename T2>
struct BoundedMaskedReductiveBoxProcessingFunctional2D_NN :
    public BoundedReductiveBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(
        Box2D domain, NTensorField2D<T1> &field1, NTensorField2D<T2> &field2,
        NTensorField2D<int> &mask) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain, NTensorField2D<T1> &field1,
        NTensorField2D<T2> &field2, NTensorField2D<int> &mask) = 0;
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, NTensorField2D<T1> &field1,
        NTensorField2D<T2> &field2, NTensorField2D<int> &mask) = 0;
};

/// Easy instantiation of bounded boxed data processor for
/// ScalarFieldN-TensorField coupling.
template <typename T1, typename T2>
struct BoundedReductiveBoxProcessingFunctional2D_SN :
    public BoundedReductiveBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(
        Box2D domain, ScalarField2D<T1> &field1, NTensorField2D<T2> &field2) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain, ScalarField2D<T1> &field1,
        NTensorField2D<T2> &field2) = 0;
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, ScalarField2D<T1> &field1,
        NTensorField2D<T2> &field2) = 0;
};

/// Easy instantiation of bounded boxed data processor for Lattice-ScalarField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct BoundedReductiveBoxProcessingFunctional2D_LS :
    public BoundedReductiveBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(
        Box2D domain, BlockLattice2D<T1, Descriptor> &lattice, ScalarField2D<T2> &field) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain, BlockLattice2D<T1, Descriptor> &lattice,
        ScalarField2D<T2> &field) = 0;
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, BlockLattice2D<T1, Descriptor> &lattice,
        ScalarField2D<T2> &field) = 0;
};

/// Easy instantiation of bounded boxed data processor for Lattice-TensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
struct BoundedReductiveBoxProcessingFunctional2D_LT :
    public BoundedReductiveBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(
        Box2D domain, BlockLattice2D<T1, Descriptor> &lattice, TensorField2D<T2, nDim> &field) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain, BlockLattice2D<T1, Descriptor> &lattice,
        TensorField2D<T2, nDim> &field) = 0;
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, BlockLattice2D<T1, Descriptor> &lattice,
        TensorField2D<T2, nDim> &field) = 0;
};

/// Easy instantiation of bounded boxed data processor for Lattice-NTensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct BoundedReductiveBoxProcessingFunctional2D_LN :
    public BoundedReductiveBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(
        Box2D domain, BlockLattice2D<T1, Descriptor> &lattice, NTensorField2D<T2> &field) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain, BlockLattice2D<T1, Descriptor> &lattice,
        NTensorField2D<T2> &field) = 0;
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, BlockLattice2D<T1, Descriptor> &lattice,
        NTensorField2D<T2> &field) = 0;
};

/// Easy instantiation of bounded boxed data processor for Lattice-NTensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct BoundedMaskedReductiveBoxProcessingFunctional2D_LN :
    public BoundedReductiveBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(
        Box2D domain, BlockLattice2D<T1, Descriptor> &lattice, NTensorField2D<T2> &field,
        NTensorField2D<int> &mask) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain, BlockLattice2D<T1, Descriptor> &lattice,
        NTensorField2D<T2> &field, NTensorField2D<int> &mask) = 0;
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, BlockLattice2D<T1, Descriptor> &lattice,
        NTensorField2D<T2> &field, NTensorField2D<int> &mask) = 0;
};

/// Easy instantiation of bounded boxed data processor for multiple lattices.
template <typename T, template <typename U> class Descriptor>
struct BoundedReductiveLatticeBoxProcessingFunctional2D :
    public BoundedReductiveBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(
        Box2D domain, std::vector<BlockLattice2D<T, Descriptor> *> lattices) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain,
        std::vector<BlockLattice2D<T, Descriptor> *> lattices) = 0;
    virtual void processCorner(
        int normalX, int normalY, Box2D domain,
        std::vector<BlockLattice2D<T, Descriptor> *> lattices) = 0;
};

/// Easy instantiation of bounded boxed data processor for multiple scalar fields.
template <typename T>
struct BoundedReductiveScalarFieldBoxProcessingFunctional2D :
    public BoundedReductiveBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(Box2D domain, std::vector<ScalarField2D<T> *> field) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain, std::vector<ScalarField2D<T> *> field) = 0;
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, std::vector<ScalarField2D<T> *> field) = 0;
};

/// Easy instantiation of bounded boxed data processor for multiple tensor fields.
template <typename T, int nDim>
struct BoundedReductiveTensorFieldBoxProcessingFunctional2D :
    public BoundedReductiveBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(Box2D domain, std::vector<TensorField2D<T, nDim> *> field) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain,
        std::vector<TensorField2D<T, nDim> *> field) = 0;
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, std::vector<TensorField2D<T, nDim> *> field) = 0;
};

/// Easy instantiation of bounded boxed data processor for multiple n-tensor fields.
template <typename T>
struct BoundedReductiveNTensorFieldBoxProcessingFunctional2D :
    public BoundedReductiveBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(Box2D domain, std::vector<NTensorField2D<T> *> field) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain, std::vector<NTensorField2D<T> *> field) = 0;
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, std::vector<NTensorField2D<T> *> field) = 0;
};

}  // namespace plb

#endif  // REDUCTIVE_DATA_PROCESSING_FUNCTIONAL_2D_H
