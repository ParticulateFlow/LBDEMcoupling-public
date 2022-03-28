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
#ifndef REDUCTIVE_DATA_PROCESSING_FUNCTIONAL_3D_H
#define REDUCTIVE_DATA_PROCESSING_FUNCTIONAL_3D_H

#include <vector>

#include "atomicBlock/dataProcessor3D.h"
#include "core/blockSurface3D.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"

namespace plb {

// Forward declarations.

template <typename T, template <typename U> class Descriptor>
class BlockLatticeBase3D;
template <typename T, template <typename U> class Descriptor>
class BlockLattice3D;
template <typename T>
class ScalarFieldBase3D;
template <typename T>
class ScalarField3D;
template <typename T, int nDim>
class TensorFieldBase3D;
template <typename T, int nDim>
class TensorField3D;
template <typename T>
class NTensorFieldBase3D;
template <typename T>
class NTensorField3D;

/* *************** All flavors of Box processing functionals ********* */

/// Easy instantiation of reductive boxed data processor (general case)
struct ReductiveBoxProcessingFunctional3D {
    virtual ~ReductiveBoxProcessingFunctional3D() { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks) = 0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(double dxScale, double dtScale);
    virtual void getDimensionsX(std::vector<int> &dimensions) const;
    virtual void getDimensionsT(std::vector<int> &dimensions) const;
    /// Obsolete: replaced by getTypeOfModification. Don't override.
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const = 0;
    virtual ReductiveBoxProcessingFunctional3D *clone() const = 0;
    virtual BlockStatistics const &getStatistics() const = 0;
    virtual BlockStatistics &getStatistics() = 0;
    virtual void serialize(std::string &data) const;
    virtual void unserialize(std::string &data);
    virtual int getStaticId() const;
};

/// ReductiveBoxProcessingFunctional3D which instantiates its own statistics object
class PlainReductiveBoxProcessingFunctional3D : public ReductiveBoxProcessingFunctional3D {
public:
    virtual BlockStatistics const &getStatistics() const;
    virtual BlockStatistics &getStatistics();
    void swap(PlainReductiveBoxProcessingFunctional3D &rhs)
    {
        std::swap(statistics, rhs.statistics);
    }

private:
    BlockStatistics statistics;
};

/// A reductive boxed data processor, automatically generated from a
/// ReductiveBoxProcessingFunctional3D
class ReductiveBoxProcessor3D : public DataProcessor3D {
public:
    /** \param functional_ The functional is not owned by the ReductiveBoxProcessor3D,
     *                     i.e. it is not deleted in the destructor.
     */
    ReductiveBoxProcessor3D(
        ReductiveBoxProcessingFunctional3D *functional_, Box3D domain_,
        std::vector<AtomicBlock3D *> atomicBlocks_);
    Box3D getDomain() const;
    virtual void process();
    virtual ReductiveBoxProcessor3D *clone() const;
    virtual int getStaticId() const;

private:
    ReductiveBoxProcessingFunctional3D *functional;
    Box3D domain;
    std::vector<AtomicBlock3D *> atomicBlocks;
};

/// An automatically created generator for the ReductiveBoxProcessor3D
class ReductiveBoxProcessorGenerator3D : public BoxedReductiveDataProcessorGenerator3D {
public:
    /** \param functional_ The functional is not owned by the ReductiveBoxProcessorGenerator3D,
     *                     i.e. it is not deleted in the destructor.
     */
    ReductiveBoxProcessorGenerator3D(ReductiveBoxProcessingFunctional3D *functional_, Box3D domain);
    ~ReductiveBoxProcessorGenerator3D();
    ReductiveBoxProcessorGenerator3D(ReductiveBoxProcessorGenerator3D const &rhs);
    ReductiveBoxProcessorGenerator3D &operator=(ReductiveBoxProcessorGenerator3D const &rhs);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(double dxScale, double dtScale);
    virtual void getDimensionsX(std::vector<int> &dimensions) const;
    virtual void getDimensionsT(std::vector<int> &dimensions) const;
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual DataProcessor3D *generate(std::vector<AtomicBlock3D *> atomicBlocks);
    virtual ReductiveBoxProcessorGenerator3D *clone() const;
    virtual BlockStatistics const &getStatistics() const;
    virtual BlockStatistics &getStatistics();
    ReductiveBoxProcessingFunctional3D const &getFunctional() const;
    virtual void serialize(Box3D &domain, std::string &data) const;
    virtual int getStaticId() const;

private:
    ReductiveBoxProcessingFunctional3D *functional;
};

/// Easy instantiation of boxed data processor for a single lattice
template <typename T, template <typename U> class Descriptor>
struct ReductiveBoxProcessingFunctional3D_L : public PlainReductiveBoxProcessingFunctional3D {
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single scalar field
template <typename T>
struct ReductiveBoxProcessingFunctional3D_S : public PlainReductiveBoxProcessingFunctional3D {
    virtual void process(Box3D domain, ScalarField3D<T> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single tensor field
template <typename T, int nDim>
struct ReductiveBoxProcessingFunctional3D_T : public PlainReductiveBoxProcessingFunctional3D {
    virtual void process(Box3D domain, TensorField3D<T, nDim> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single n-tensor field
template <typename T>
struct ReductiveBoxProcessingFunctional3D_N : public PlainReductiveBoxProcessingFunctional3D {
    virtual void process(Box3D domain, NTensorField3D<T> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single n-tensor field
template <typename T>
struct MaskedReductiveBoxProcessingFunctional3D_N : public PlainReductiveBoxProcessingFunctional3D {
    virtual void process(Box3D domain, NTensorField3D<T> &field, NTensorField3D<int> &mask) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for lattice-lattice coupling
template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
struct ReductiveBoxProcessingFunctional3D_LL : public PlainReductiveBoxProcessingFunctional3D {
    virtual void process(
        Box3D domain, BlockLattice3D<T1, Descriptor1> &lattice1,
        BlockLattice3D<T2, Descriptor2> &lattice2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for ScalarField-ScalarField coupling
template <typename T1, typename T2>
struct ReductiveBoxProcessingFunctional3D_SS : public PlainReductiveBoxProcessingFunctional3D {
    virtual void process(Box3D domain, ScalarField3D<T1> &field1, ScalarField3D<T2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for TensorField-TensorField coupling
template <typename T1, int nDim1, typename T2, int nDim2>
struct ReductiveBoxProcessingFunctional3D_TT : public PlainReductiveBoxProcessingFunctional3D {
    virtual void process(
        Box3D domain, TensorField3D<T1, nDim1> &field1, TensorField3D<T2, nDim2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for ScalarField-TensorField coupling
template <typename T1, typename T2, int nDim>
struct ReductiveBoxProcessingFunctional3D_ST : public PlainReductiveBoxProcessingFunctional3D {
    virtual void process(
        Box3D domain, ScalarField3D<T1> &field1, TensorField3D<T2, nDim> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for NTensorField-NTensorField coupling
template <typename T1, typename T2>
struct ReductiveBoxProcessingFunctional3D_NN : public PlainReductiveBoxProcessingFunctional3D {
    virtual void process(Box3D domain, NTensorField3D<T1> &field1, NTensorField3D<T2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for NTensorField-NTensorField coupling
template <typename T1, typename T2>
struct MaskedReductiveBoxProcessingFunctional3D_NN :
    public PlainReductiveBoxProcessingFunctional3D {
    virtual void process(
        Box3D domain, NTensorField3D<T1> &field1, NTensorField3D<T2> &field2,
        NTensorField3D<int> &mask) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for ScalarField-NTensorField coupling
template <typename T1, typename T2>
struct ReductiveBoxProcessingFunctional3D_SN : public PlainReductiveBoxProcessingFunctional3D {
    virtual void process(Box3D domain, ScalarField3D<T1> &field1, NTensorField3D<T2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-ScalarField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct ReductiveBoxProcessingFunctional3D_LS : public PlainReductiveBoxProcessingFunctional3D {
    virtual void process(
        Box3D domain, BlockLattice3D<T1, Descriptor> &lattice, ScalarField3D<T2> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-TensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
struct ReductiveBoxProcessingFunctional3D_LT : public PlainReductiveBoxProcessingFunctional3D {
    virtual void process(
        Box3D domain, BlockLattice3D<T1, Descriptor> &lattice, TensorField3D<T2, nDim> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-NTensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct ReductiveBoxProcessingFunctional3D_LN : public PlainReductiveBoxProcessingFunctional3D {
    virtual void process(
        Box3D domain, BlockLattice3D<T1, Descriptor> &lattice, NTensorField3D<T2> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-NTensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct MaskedReductiveBoxProcessingFunctional3D_LN :
    public PlainReductiveBoxProcessingFunctional3D {
    virtual void process(
        Box3D domain, BlockLattice3D<T1, Descriptor> &lattice, NTensorField3D<T2> &field,
        NTensorField3D<int> &mask) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple lattices
template <typename T, template <typename U> class Descriptor>
struct ReductiveLatticeBoxProcessingFunctional3D : public PlainReductiveBoxProcessingFunctional3D {
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T, Descriptor> *> lattices) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple scalarFields
template <typename T>
struct ReductiveScalarFieldBoxProcessingFunctional3D :
    public PlainReductiveBoxProcessingFunctional3D {
    virtual void process(Box3D domain, std::vector<ScalarField3D<T> *> scalarFields) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple tensorFields
template <typename T, int nDim>
struct ReductiveTensorFieldBoxProcessingFunctional3D :
    public PlainReductiveBoxProcessingFunctional3D {
    virtual void process(Box3D domain, std::vector<TensorField3D<T, nDim> *> tensorFields) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple n-tensorFields
template <typename T>
struct ReductiveNTensorFieldBoxProcessingFunctional3D :
    public PlainReductiveBoxProcessingFunctional3D {
    virtual void process(Box3D domain, std::vector<NTensorField3D<T> *> tensorFields) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/* *************** All flavors of Dot processing functionals ********* */

/// Easy instantiation of reductive boxed data processor (general case)
class ReductiveDotProcessingFunctional3D {
public:
    virtual ~ReductiveDotProcessingFunctional3D() { }
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks) = 0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(double dxScale, double dtScale);
    virtual void getDimensionsX(std::vector<int> &dimensions) const;
    virtual void getDimensionsT(std::vector<int> &dimensions) const;
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const = 0;
    virtual ReductiveDotProcessingFunctional3D *clone() const = 0;
    virtual BlockStatistics const &getStatistics() const = 0;
    virtual BlockStatistics &getStatistics() = 0;
};

/// ReductiveDotProcessingFunctional3D which instantiates its own statistics object
class PlainReductiveDotProcessingFunctional3D : public ReductiveDotProcessingFunctional3D {
public:
    virtual BlockStatistics const &getStatistics() const;
    virtual BlockStatistics &getStatistics();

private:
    BlockStatistics statistics;
};

/// A ReductiveDotted data processor, automatically generated from a
/// ReductiveDotProcessingFunctional3D
class ReductiveDotProcessor3D : public DataProcessor3D {
public:
    /** \param functional_ The functional is not owned by the ReductiveBoxProcessor3D,
     *                     i.e. it is not deleted in the destructor.
     */
    ReductiveDotProcessor3D(
        ReductiveDotProcessingFunctional3D *functional_, DotList3D const &dotList_,
        std::vector<AtomicBlock3D *> atomicBlocks_);
    DotList3D const &getDotList() const;
    virtual void process();
    virtual ReductiveDotProcessor3D *clone() const;

private:
    ReductiveDotProcessingFunctional3D *functional;
    DotList3D dotList;
    std::vector<AtomicBlock3D *> atomicBlocks;
};

/// An automatically created generator for the ReductiveDotProcessor3D
class ReductiveDotProcessorGenerator3D : public DottedReductiveDataProcessorGenerator3D {
public:
    ReductiveDotProcessorGenerator3D(
        ReductiveDotProcessingFunctional3D *functional_, DotList3D const &dotList);
    ~ReductiveDotProcessorGenerator3D();
    ReductiveDotProcessorGenerator3D(ReductiveDotProcessorGenerator3D const &rhs);
    ReductiveDotProcessorGenerator3D &operator=(ReductiveDotProcessorGenerator3D const &rhs);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(double dxScale, double dtScale);
    virtual void getDimensionsX(std::vector<int> &dimensions) const;
    virtual void getDimensionsT(std::vector<int> &dimensions) const;
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual DataProcessor3D *generate(std::vector<AtomicBlock3D *> atomicBlocks);
    virtual ReductiveDotProcessorGenerator3D *clone() const;
    virtual BlockStatistics const &getStatistics() const;
    virtual BlockStatistics &getStatistics();
    ReductiveDotProcessingFunctional3D const &getFunctional() const;

private:
    ReductiveDotProcessingFunctional3D *functional;
};

/// Easy instantiation of dotted data processor for a single lattice
template <typename T, template <typename U> class Descriptor>
struct ReductiveDotProcessingFunctional3D_L : public PlainReductiveDotProcessingFunctional3D {
    virtual void process(DotList3D const &dotList, BlockLattice3D<T, Descriptor> &lattice) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for a single scalar field
template <typename T>
struct ReductiveDotProcessingFunctional3D_S : public PlainReductiveDotProcessingFunctional3D {
    virtual void process(DotList3D const &dotList, ScalarField3D<T> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for a single tensor field
template <typename T, int nDim>
struct ReductiveDotProcessingFunctional3D_T : public PlainReductiveDotProcessingFunctional3D {
    virtual void process(DotList3D const &dotList, TensorField3D<T, nDim> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for a single n-tensor field
template <typename T>
struct ReductiveDotProcessingFunctional3D_N : public PlainReductiveDotProcessingFunctional3D {
    virtual void process(DotList3D const &dotList, NTensorField3D<T> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for lattice-lattice coupling
template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
struct ReductiveDotProcessingFunctional3D_LL : public PlainReductiveDotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, BlockLattice3D<T1, Descriptor1> &lattice1,
        BlockLattice3D<T2, Descriptor2> &lattice2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for ScalarField-ScalarField coupling
template <typename T1, typename T2>
struct ReductiveDotProcessingFunctional3D_SS : public PlainReductiveDotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, ScalarField3D<T1> &field1, ScalarField3D<T2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for TensorField-TensorField coupling
template <typename T1, int nDim1, typename T2, int nDim2>
struct ReductiveDotProcessingFunctional3D_TT : public PlainReductiveDotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, TensorField3D<T1, nDim1> &field1,
        TensorField3D<T2, nDim2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for ScalarField-TensorField coupling
template <typename T1, typename T2, int nDim>
struct ReductiveDotProcessingFunctional3D_ST : public PlainReductiveDotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, ScalarField3D<T1> &field1, TensorField3D<T2, nDim> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for NTensorField-NTensorField coupling
template <typename T1, typename T2>
struct ReductiveDotProcessingFunctional3D_NN : public PlainReductiveDotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, NTensorField3D<T1> &field1, NTensorField3D<T2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for ScalarField-NTensorField coupling
template <typename T1, typename T2>
struct ReductiveDotProcessingFunctional3D_SN : public PlainReductiveDotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, ScalarField3D<T1> &field1, NTensorField3D<T2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for Lattice-ScalarField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct ReductiveDotProcessingFunctional3D_LS : public PlainReductiveDotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, BlockLattice3D<T1, Descriptor> &lattice,
        ScalarField3D<T2> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for Lattice-TensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
struct ReductiveDotProcessingFunctional3D_LT : public PlainReductiveDotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, BlockLattice3D<T1, Descriptor> &lattice,
        TensorField3D<T2, nDim> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for Lattice-NTensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct ReductiveDotProcessingFunctional3D_LN : public PlainReductiveDotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, BlockLattice3D<T1, Descriptor> &lattice,
        NTensorField3D<T2> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple lattices
template <typename T, template <typename U> class Descriptor>
struct ReductiveLatticeDotProcessingFunctional3D : public PlainReductiveDotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, std::vector<BlockLattice3D<T, Descriptor> *> lattices) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple scalarFields
template <typename T>
struct ReductiveScalarFieldDotProcessingFunctional3D :
    public PlainReductiveDotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, std::vector<ScalarField3D<T> *> scalarFields) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple tensorFields
template <typename T, int nDim>
struct ReductiveTensorFieldDotProcessingFunctional3D :
    public PlainReductiveDotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, std::vector<TensorField3D<T, nDim> *> tensorFields) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple n-tensorFields
template <typename T>
struct ReductiveNTensorFieldDotProcessingFunctional3D :
    public PlainReductiveDotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, std::vector<NTensorField3D<T> *> tensorFields) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/* *************** All flavors of Bounded Box processing functionals ********* */

/// Easy instantiation of boxed data processor special boundary treatment (general case)
class BoundedReductiveBoxProcessingFunctional3D {
public:
    virtual ~BoundedReductiveBoxProcessingFunctional3D() { }
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks) = 0;
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks) = 0;
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks) = 0;
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks) = 0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(double dxScale, double dtScale);
    virtual void getDimensionsX(std::vector<int> &dimensions) const;
    virtual void getDimensionsT(std::vector<int> &dimensions) const;
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const = 0;
    virtual BoundedReductiveBoxProcessingFunctional3D *clone() const = 0;
    ReductiveBoxProcessingFunctional3D *getBulkProcessor() const;
    ReductiveBoxProcessingFunctional3D *getPlaneProcessor(int direction, int orientation) const;
    ReductiveBoxProcessingFunctional3D *getEdgeProcessor(int plane, int normal1, int normal2) const;
    ReductiveBoxProcessingFunctional3D *getCornerProcessor(
        int normalX, int normalY, int normalZ) const;
    BlockStatistics const &getStatistics() const;
    BlockStatistics &getStatistics();
    void getGenerators(
        Box3D const &fullDomain, plint boundaryWidth_,
        std::vector<ReductiveBoxProcessorGenerator3D *> &generators);

private:
    BlockStatistics statistics;

public:
    class BulkWrapperFunctional : public ReductiveBoxProcessingFunctional3D {
    public:
        BulkWrapperFunctional(BoundedReductiveBoxProcessingFunctional3D *boundedFunctional_);
        BulkWrapperFunctional(BulkWrapperFunctional const &rhs);
        ~BulkWrapperFunctional();
        BulkWrapperFunctional &operator=(BulkWrapperFunctional const &rhs);
        virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
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
        BoundedReductiveBoxProcessingFunctional3D *boundedFunctional;
    };
    class PlaneWrapperFunctional : public ReductiveBoxProcessingFunctional3D {
    public:
        PlaneWrapperFunctional(
            BoundedReductiveBoxProcessingFunctional3D *boundedFunctional_, int direction_,
            int orientation_);
        PlaneWrapperFunctional(PlaneWrapperFunctional const &rhs);
        ~PlaneWrapperFunctional();
        PlaneWrapperFunctional &operator=(PlaneWrapperFunctional const &rhs);
        virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        virtual void rescale(double dxScale, double dtScale);
        virtual void getDimensionsX(std::vector<int> &dimensions) const;
        virtual void getDimensionsT(std::vector<int> &dimensions) const;
        virtual void getModificationPattern(std::vector<bool> &isWritten) const;
        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
        virtual PlaneWrapperFunctional *clone() const;
        virtual BlockStatistics const &getStatistics() const;
        virtual BlockStatistics &getStatistics();

    private:
        BoundedReductiveBoxProcessingFunctional3D *boundedFunctional;
        int direction, orientation;
    };
    class EdgeWrapperFunctional : public ReductiveBoxProcessingFunctional3D {
    public:
        EdgeWrapperFunctional(
            BoundedReductiveBoxProcessingFunctional3D *boundedFunctional_, int plane_, int normal1_,
            int normal2_);
        EdgeWrapperFunctional(EdgeWrapperFunctional const &rhs);
        ~EdgeWrapperFunctional();
        EdgeWrapperFunctional &operator=(EdgeWrapperFunctional const &rhs);
        virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
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
        BoundedReductiveBoxProcessingFunctional3D *boundedFunctional;
        int plane, normal1, normal2;
    };
    class CornerWrapperFunctional : public ReductiveBoxProcessingFunctional3D {
    public:
        CornerWrapperFunctional(
            BoundedReductiveBoxProcessingFunctional3D *boundedFunctional_, int normalX_,
            int normalY_, int normalZ_);
        CornerWrapperFunctional(CornerWrapperFunctional const &rhs);
        ~CornerWrapperFunctional();
        CornerWrapperFunctional &operator=(CornerWrapperFunctional const &rhs);
        virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
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
        BoundedReductiveBoxProcessingFunctional3D *boundedFunctional;
        int normalX, normalY, normalZ;
    };
};

/// Easy instantiation of bounded boxed data processor for a single lattice
template <typename T, template <typename U> class Descriptor>
struct BoundedReductiveBoxProcessingFunctional3D_L :
    public BoundedReductiveBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual void processBulk(Box3D domain, BlockLattice3D<T, Descriptor> &lattice) = 0;
    virtual void processPlane(
        int direction, int orientation, Box3D domain, BlockLattice3D<T, Descriptor> &lattice) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain,
        BlockLattice3D<T, Descriptor> &lattice) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain,
        BlockLattice3D<T, Descriptor> &lattice) = 0;
};

/// Easy instantiation of bounded boxed data processor for a single scalar field
template <typename T>
struct BoundedReductiveBoxProcessingFunctional3D_S :
    public BoundedReductiveBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual void processBulk(Box3D domain, ScalarField3D<T> &field) = 0;
    virtual void processPlane(
        int direction, int orientation, Box3D domain, ScalarField3D<T> &field) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain, ScalarField3D<T> &field) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain, ScalarField3D<T> &field) = 0;
};

/// Easy instantiation of bounded boxed data processor for a single tensor field
template <typename T, int nDim>
struct BoundedReductiveBoxProcessingFunctional3D_T :
    public BoundedReductiveBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual void processBulk(Box3D domain, TensorField3D<T, nDim> &lattice) = 0;
    virtual void processPlane(
        int direction, int orientation, Box3D domain, TensorField3D<T, nDim> &lattice) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain, TensorField3D<T, nDim> &lattice) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain, TensorField3D<T, nDim> &lattice) = 0;
};

/// Easy instantiation of bounded boxed data processor for a single n-tensor field
template <typename T>
struct BoundedReductiveBoxProcessingFunctional3D_N :
    public BoundedReductiveBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual void processBulk(Box3D domain, NTensorField3D<T> &lattice) = 0;
    virtual void processPlane(
        int direction, int orientation, Box3D domain, NTensorField3D<T> &lattice) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain, NTensorField3D<T> &lattice) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain, NTensorField3D<T> &lattice) = 0;
};

/// Easy instantiation of bounded boxed data processor for a single n-tensor field with bool mask
template <typename T>
struct BoundedMaskedReductiveBoxProcessingFunctional3D_N :
    public BoundedReductiveBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual void processBulk(
        Box3D domain, NTensorField3D<T> &lattice, NTensorField3D<int> &mask) = 0;
    virtual void processPlane(
        int direction, int orientation, Box3D domain, NTensorField3D<T> &lattice,
        NTensorField3D<int> &mask) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain, NTensorField3D<T> &lattice,
        NTensorField3D<int> &mask) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain, NTensorField3D<T> &lattice,
        NTensorField3D<int> &mask) = 0;
};

/// Easy instantiation of bounded boxed data processor for lattice-lattice coupling
template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
struct BoundedReductiveBoxProcessingFunctional3D_LL :
    public BoundedReductiveBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual void processBulk(
        Box3D domain, BlockLattice3D<T1, Descriptor1> &lattice1,
        BlockLattice3D<T2, Descriptor2> &lattice2) = 0;
    virtual void processPlane(
        int direction, int orientation, Box3D domain, BlockLattice3D<T1, Descriptor1> &lattice1,
        BlockLattice3D<T2, Descriptor2> &lattice2) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain,
        BlockLattice3D<T1, Descriptor1> &lattice1, BlockLattice3D<T2, Descriptor2> &lattice2) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain,
        BlockLattice3D<T1, Descriptor1> &lattice1, BlockLattice3D<T2, Descriptor2> &lattice2) = 0;
};

/// Easy instantiation of bounded boxed data processor for ScalarField-ScalarField coupling
template <typename T1, typename T2>
struct BoundedReductiveBoxProcessingFunctional3D_SS :
    public BoundedReductiveBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual void processBulk(
        Box3D domain, ScalarField3D<T1> &field1, ScalarField3D<T2> &field2) = 0;
    virtual void processPlane(
        int direction, int orientation, Box3D domain, ScalarField3D<T1> &field1,
        ScalarField3D<T2> &field2) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain, ScalarField3D<T1> &field1,
        ScalarField3D<T2> &field2) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain, ScalarField3D<T1> &field1,
        ScalarField3D<T2> &field2) = 0;
};

/// Easy instantiation of bounded boxed data processor for TensorField-TensorField coupling
template <typename T1, int nDim1, typename T2, int nDim2>
struct BoundedReductiveBoxProcessingFunctional3D_TT :
    public BoundedReductiveBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual void processBulk(
        Box3D domain, TensorField3D<T1, nDim1> &field1, TensorField3D<T2, nDim2> &field2) = 0;
    virtual void processPlane(
        int direction, int orientation, Box3D domain, TensorField3D<T1, nDim1> &field1,
        TensorField3D<T2, nDim2> &field2) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain, TensorField3D<T1, nDim1> &field1,
        TensorField3D<T2, nDim2> &field2) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain, TensorField3D<T1, nDim1> &field1,
        TensorField3D<T2, nDim2> &field2) = 0;
};

/// Easy instantiation of bounded boxed data processor for ScalarField-TensorField coupling
template <typename T1, typename T2, int nDim>
struct BoundedReductiveBoxProcessingFunctional3D_ST :
    public BoundedReductiveBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual void processBulk(
        Box3D domain, ScalarField3D<T1> &field1, TensorField3D<T2, nDim> &field2) = 0;
    virtual void processPlane(
        int direction, int orientation, Box3D domain, ScalarField3D<T1> &field1,
        TensorField3D<T2, nDim> &field2) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain, ScalarField3D<T1> &field1,
        TensorField3D<T2, nDim> &field2) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain, ScalarField3D<T1> &field1,
        TensorField3D<T2, nDim> &field2) = 0;
};

/// Easy instantiation of bounded boxed data processor for NTensorField-NTensorField coupling
template <typename T1, typename T2>
struct BoundedReductiveBoxProcessingFunctional3D_NN :
    public BoundedReductiveBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual void processBulk(
        Box3D domain, NTensorField3D<T1> &field1, NTensorField3D<T2> &field2) = 0;
    virtual void processPlane(
        int direction, int orientation, Box3D domain, NTensorField3D<T1> &field1,
        NTensorField3D<T2> &field2) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain, NTensorField3D<T1> &field1,
        NTensorField3D<T2> &field2) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain, NTensorField3D<T1> &field1,
        NTensorField3D<T2> &field2) = 0;
};

/// Easy instantiation of bounded boxed data processor for NTensorField-NTensorField coupling
template <typename T1, typename T2>
struct BoundedMaskedReductiveBoxProcessingFunctional3D_NN :
    public BoundedReductiveBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual void processBulk(
        Box3D domain, NTensorField3D<T1> &field1, NTensorField3D<T2> &field2,
        NTensorField3D<int> &mask) = 0;
    virtual void processPlane(
        int direction, int orientation, Box3D domain, NTensorField3D<T1> &field1,
        NTensorField3D<T2> &field2, NTensorField3D<int> &mask) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain, NTensorField3D<T1> &field1,
        NTensorField3D<T2> &field2, NTensorField3D<int> &mask) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain, NTensorField3D<T1> &field1,
        NTensorField3D<T2> &field2, NTensorField3D<int> &mask) = 0;
};

/// Easy instantiation of bounded boxed data processor for ScalarField-NTensorField coupling
template <typename T1, typename T2>
struct BoundedReductiveBoxProcessingFunctional3D_SN :
    public BoundedReductiveBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual void processBulk(
        Box3D domain, ScalarField3D<T1> &field1, NTensorField3D<T2> &field2) = 0;
    virtual void processPlane(
        int direction, int orientation, Box3D domain, ScalarField3D<T1> &field1,
        NTensorField3D<T2> &field2) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain, ScalarField3D<T1> &field1,
        NTensorField3D<T2> &field2) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain, ScalarField3D<T1> &field1,
        NTensorField3D<T2> &field2) = 0;
};

/// Easy instantiation of bounded boxed data processor for Lattice-ScalarField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct BoundedReductiveBoxProcessingFunctional3D_LS :
    public BoundedReductiveBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual void processBulk(
        Box3D domain, BlockLattice3D<T1, Descriptor> &lattice, ScalarField3D<T2> &field) = 0;
    virtual void processPlane(
        int direction, int orientation, Box3D domain, BlockLattice3D<T1, Descriptor> &lattice,
        ScalarField3D<T2> &field) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain, BlockLattice3D<T1, Descriptor> &lattice,
        ScalarField3D<T2> &field) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain,
        BlockLattice3D<T1, Descriptor> &lattice, ScalarField3D<T2> &field) = 0;
};

/// Easy instantiation of bounded boxed data processor for Lattice-TensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
struct BoundedReductiveBoxProcessingFunctional3D_LT :
    public BoundedReductiveBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual void processBulk(
        Box3D domain, BlockLattice3D<T1, Descriptor> &lattice, TensorField3D<T2, nDim> &field) = 0;
    virtual void processPlane(
        int direction, int orientation, Box3D domain, BlockLattice3D<T1, Descriptor> &lattice,
        TensorField3D<T2, nDim> &field) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain, BlockLattice3D<T1, Descriptor> &lattice,
        TensorField3D<T2, nDim> &field) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain,
        BlockLattice3D<T1, Descriptor> &lattice, TensorField3D<T2, nDim> &field) = 0;
};

/// Easy instantiation of bounded boxed data processor for Lattice-NTensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct BoundedReductiveBoxProcessingFunctional3D_LN :
    public BoundedReductiveBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual void processBulk(
        Box3D domain, BlockLattice3D<T1, Descriptor> &lattice, NTensorField3D<T2> &field) = 0;
    virtual void processPlane(
        int direction, int orientation, Box3D domain, BlockLattice3D<T1, Descriptor> &lattice,
        NTensorField3D<T2> &field) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain, BlockLattice3D<T1, Descriptor> &lattice,
        NTensorField3D<T2> &field) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain,
        BlockLattice3D<T1, Descriptor> &lattice, NTensorField3D<T2> &field) = 0;
};

/// Easy instantiation of bounded boxed data processor for Lattice-NTensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct BoundedMaskedReductiveBoxProcessingFunctional3D_LN :
    public BoundedReductiveBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual void processBulk(
        Box3D domain, BlockLattice3D<T1, Descriptor> &lattice, NTensorField3D<T2> &field,
        NTensorField3D<int> &mask) = 0;
    virtual void processPlane(
        int direction, int orientation, Box3D domain, BlockLattice3D<T1, Descriptor> &lattice,
        NTensorField3D<T2> &field, NTensorField3D<int> &mask) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain, BlockLattice3D<T1, Descriptor> &lattice,
        NTensorField3D<T2> &field, NTensorField3D<int> &mask) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain,
        BlockLattice3D<T1, Descriptor> &lattice, NTensorField3D<T2> &field,
        NTensorField3D<int> &mask) = 0;
};

/// Easy instantiation of bounded boxed data processor for multiple lattices.
template <typename T, template <typename U> class Descriptor>
struct BoundedReductiveLatticeBoxProcessingFunctional3D :
    public BoundedReductiveBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual void processBulk(
        Box3D domain, std::vector<BlockLattice3D<T, Descriptor> *> lattices) = 0;
    virtual void processPlane(
        int direction, int orientation, Box3D domain,
        std::vector<BlockLattice3D<T, Descriptor> *> lattices) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<BlockLattice3D<T, Descriptor> *> lattices) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<BlockLattice3D<T, Descriptor> *> lattices) = 0;
};

/// Easy instantiation of bounded boxed data processor for multiple scalar fields.
template <typename T>
struct BoundedReductiveScalarFieldBoxProcessingFunctional3D :
    public BoundedReductiveBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual void processBulk(Box3D domain, std::vector<ScalarField3D<T> *> field) = 0;
    virtual void processPlane(
        int direction, int orientation, Box3D domain, std::vector<ScalarField3D<T> *> field) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<ScalarField3D<T> *> field) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<ScalarField3D<T> *> field) = 0;
};

/// Easy instantiation of bounded boxed data processor for multiple tensor fields.
template <typename T, int nDim>
struct BoundedReductiveTensorFieldBoxProcessingFunctional3D :
    public BoundedReductiveBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual void processBulk(Box3D domain, std::vector<TensorField3D<T, nDim> *> field) = 0;
    virtual void processPlane(
        int direction, int orientation, Box3D domain,
        std::vector<TensorField3D<T, nDim> *> field) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<TensorField3D<T, nDim> *> field) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<TensorField3D<T, nDim> *> field) = 0;
};

/// Easy instantiation of bounded boxed data processor for multiple n-tensor fields.
template <typename T>
struct BoundedReductiveNTensorFieldBoxProcessingFunctional3D :
    public BoundedReductiveBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual void processBulk(Box3D domain, std::vector<NTensorField3D<T> *> field) = 0;
    virtual void processPlane(
        int direction, int orientation, Box3D domain, std::vector<NTensorField3D<T> *> field) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<NTensorField3D<T> *> field) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<NTensorField3D<T> *> field) = 0;
};

}  // namespace plb

#endif  // REDUCTIVE_DATA_PROCESSING_FUNCTIONAL_3D_H
