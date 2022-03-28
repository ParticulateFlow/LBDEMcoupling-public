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
#ifndef DATA_PROCESSING_FUNCTIONAL_2D_H
#define DATA_PROCESSING_FUNCTIONAL_2D_H

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

/// Easy instantiation of boxed data processor (general case)
class BoxProcessingFunctional2D {
public:
    BoxProcessingFunctional2D();
    virtual ~BoxProcessingFunctional2D() { }
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) = 0;
    virtual BlockDomain::DomainT appliesTo() const;
    /// Obsolete: replaced by setscale.
    virtual void rescale(double dxScale, double dtScale);
    virtual void setscale(int dxScale_, int dtScale_);
    /// Obsolete: replaced by getTypeOfModification. Don't override.
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const = 0;
    virtual BoxProcessingFunctional2D *clone() const = 0;
    int getDxScale() const;
    int getDtScale() const;
    virtual void serialize(std::string &data) const;
    virtual void unserialize(std::string &data);
    virtual int getStaticId() const;

private:
    int dxScale, dtScale;
};

/// A Boxed data processor, automatically generated from a BoxProcessingFunctional2D
class BoxProcessor2D : public DataProcessor2D {
public:
    BoxProcessor2D(
        BoxProcessingFunctional2D *functional_, Box2D domain_,
        std::vector<AtomicBlock2D *> atomicBlocks_);
    BoxProcessor2D(BoxProcessor2D const &rhs);
    BoxProcessor2D &operator=(BoxProcessor2D const &rhs);
    ~BoxProcessor2D();
    Box2D getDomain() const;
    virtual void process();
    virtual BoxProcessor2D *clone() const;
    virtual int getStaticId() const;

private:
    BoxProcessingFunctional2D *functional;
    Box2D domain;
    std::vector<AtomicBlock2D *> atomicBlocks;
};

/// An automatically created generator for the BoxProcessor2D
class BoxProcessorGenerator2D : public BoxedDataProcessorGenerator2D {
public:
    BoxProcessorGenerator2D(BoxProcessingFunctional2D *functional_, Box2D domain);
    ~BoxProcessorGenerator2D();
    BoxProcessorGenerator2D(BoxProcessorGenerator2D const &rhs);
    BoxProcessorGenerator2D &operator=(BoxProcessorGenerator2D const &rhs);
    virtual BlockDomain::DomainT appliesTo() const;
    /// Obsolete: replaced by setscale.
    virtual void rescale(double dxScale, double dtScale);
    virtual void setscale(int dxScale_, int dtScale_);
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual DataProcessor2D *generate(std::vector<AtomicBlock2D *> atomicBlocks) const;
    virtual BoxProcessorGenerator2D *clone() const;
    virtual void serialize(Box2D &domain, std::string &data) const;
    virtual int getStaticId() const;

private:
    BoxProcessingFunctional2D *functional;
};

/// Easy instantiation of boxed data processor for a single lattice
template <typename T, template <typename U> class Descriptor>
struct BoxProcessingFunctional2D_L : public BoxProcessingFunctional2D {
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single lattice with bool mask.
template <typename T, template <typename U> class Descriptor>
struct MaskedBoxProcessingFunctional2D_L : public BoxProcessingFunctional2D {
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, NTensorField2D<int> &mask) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single scalar field
template <typename T>
struct BoxProcessingFunctional2D_S : public BoxProcessingFunctional2D {
    virtual void process(Box2D domain, ScalarField2D<T> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single tensor field
template <typename T, int nDim>
struct BoxProcessingFunctional2D_T : public BoxProcessingFunctional2D {
    virtual void process(Box2D domain, TensorField2D<T, nDim> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single n-tensor field
template <typename T>
struct BoxProcessingFunctional2D_N : public BoxProcessingFunctional2D {
    virtual void process(Box2D domain, NTensorField2D<T> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single n-tensor field with bool mask
template <typename T>
struct MaskedBoxProcessingFunctional2D_N : public BoxProcessingFunctional2D {
    virtual void process(Box2D domain, NTensorField2D<T> &field, NTensorField2D<int> &mask) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for lattice-lattice coupling
template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
struct BoxProcessingFunctional2D_LL : public BoxProcessingFunctional2D {
    virtual void process(
        Box2D domain, BlockLattice2D<T1, Descriptor1> &lattice1,
        BlockLattice2D<T2, Descriptor2> &lattice2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for ScalarField-ScalarField coupling
template <typename T1, typename T2>
struct BoxProcessingFunctional2D_SS : public BoxProcessingFunctional2D {
    virtual void process(Box2D domain, ScalarField2D<T1> &field1, ScalarField2D<T2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for TensorField-TensorField coupling
template <typename T1, int nDim1, typename T2, int nDim2>
struct BoxProcessingFunctional2D_TT : public BoxProcessingFunctional2D {
    virtual void process(
        Box2D domain, TensorField2D<T1, nDim1> &field1, TensorField2D<T2, nDim2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for ScalarField-TensorField coupling
template <typename T1, typename T2, int nDim>
struct BoxProcessingFunctional2D_ST : public BoxProcessingFunctional2D {
    virtual void process(
        Box2D domain, ScalarField2D<T1> &field1, TensorField2D<T2, nDim> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for NTensorField-NTensorField coupling
template <typename T1, typename T2>
struct BoxProcessingFunctional2D_NN : public BoxProcessingFunctional2D {
    virtual void process(Box2D domain, NTensorField2D<T1> &field1, NTensorField2D<T2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for NTensorField-NTensorField coupling
template <typename T1, typename T2>
struct MaskedBoxProcessingFunctional2D_NN : public BoxProcessingFunctional2D {
    virtual void process(
        Box2D domain, NTensorField2D<T1> &field1, NTensorField2D<T2> &field2,
        NTensorField2D<int> &mask) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for ScalarField-NTensorField coupling
template <typename T1, typename T2>
struct BoxProcessingFunctional2D_SN : public BoxProcessingFunctional2D {
    virtual void process(Box2D domain, ScalarField2D<T1> &field1, NTensorField2D<T2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-ScalarField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct BoxProcessingFunctional2D_LS : public BoxProcessingFunctional2D {
    virtual void process(
        Box2D domain, BlockLattice2D<T1, Descriptor> &lattice, ScalarField2D<T2> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-TensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
struct BoxProcessingFunctional2D_LT : public BoxProcessingFunctional2D {
    virtual void process(
        Box2D domain, BlockLattice2D<T1, Descriptor> &lattice, TensorField2D<T2, nDim> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-NTensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct BoxProcessingFunctional2D_LN : public BoxProcessingFunctional2D {
    virtual void process(
        Box2D domain, BlockLattice2D<T1, Descriptor> &lattice, NTensorField2D<T2> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-NTensorField coupling with mask.
template <typename T1, template <typename U> class Descriptor, typename T2>
struct MaskedBoxProcessingFunctional2D_LN : public BoxProcessingFunctional2D {
    virtual void process(
        Box2D domain, BlockLattice2D<T1, Descriptor> &lattice, NTensorField2D<T2> &field,
        NTensorField2D<int> &mask) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple lattices
template <typename T, template <typename U> class Descriptor>
struct LatticeBoxProcessingFunctional2D : public BoxProcessingFunctional2D {
    virtual void process(Box2D domain, std::vector<BlockLattice2D<T, Descriptor> *> lattices) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple scalarFields
template <typename T>
struct ScalarFieldBoxProcessingFunctional2D : public BoxProcessingFunctional2D {
    virtual void process(Box2D domain, std::vector<ScalarField2D<T> *> scalarFields) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple tensorFields
template <typename T, int nDim>
struct TensorFieldBoxProcessingFunctional2D : public BoxProcessingFunctional2D {
    virtual void process(Box2D domain, std::vector<TensorField2D<T, nDim> *> tensorFields) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple n-tensorFields
template <typename T>
struct NTensorFieldBoxProcessingFunctional2D : public BoxProcessingFunctional2D {
    virtual void process(Box2D domain, std::vector<NTensorField2D<T> *> tensorFields) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple n-tensorFields
template <typename T>
struct MaskedNTensorFieldBoxProcessingFunctional2D : public BoxProcessingFunctional2D {
    virtual void process(
        Box2D domain, std::vector<NTensorField2D<T> *> tensorFields, NTensorField2D<int> &mask) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
};

/* *************** All flavors of Dot processing functionals ********* */

/// Easy instantiation of dotted data processor (general case)
struct DotProcessingFunctional2D {
    virtual ~DotProcessingFunctional2D() { }
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks) = 0;
    virtual BlockDomain::DomainT appliesTo() const;
    /// Obsolete: replaced by setscale.
    virtual void rescale(double dxScale, double dtScale);
    virtual void setscale(int dxScale, int dtScale);
    /// Obsolete: replaced by getTypeOfModification
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const = 0;
    virtual DotProcessingFunctional2D *clone() const = 0;
};

/// A Dotted data processor, automatically generated from a DotProcessingFunctional2D
class DotProcessor2D : public DataProcessor2D {
public:
    DotProcessor2D(
        DotProcessingFunctional2D *functional_, DotList2D const &dotList_,
        std::vector<AtomicBlock2D *> atomicBlocks_);
    DotProcessor2D(DotProcessor2D const &rhs);
    DotProcessor2D &operator=(DotProcessor2D const &rhs);
    ~DotProcessor2D();
    virtual void process();
    virtual DotProcessor2D *clone() const;
    DotList2D const &getDotList() const;

private:
    DotProcessingFunctional2D *functional;
    DotList2D dotList;
    std::vector<AtomicBlock2D *> atomicBlocks;
};

/// An automatically created generator for the DotProcessor2D
class DotProcessorGenerator2D : public DottedDataProcessorGenerator2D {
public:
    DotProcessorGenerator2D(DotProcessingFunctional2D *functional_, DotList2D const &dotList);
    ~DotProcessorGenerator2D();
    DotProcessorGenerator2D(DotProcessorGenerator2D const &rhs);
    DotProcessorGenerator2D &operator=(DotProcessorGenerator2D const &rhs);
    virtual BlockDomain::DomainT appliesTo() const;
    /// Obsolete: replaced by setscale.
    virtual void rescale(double dxScale, double dtScale);
    virtual void setscale(int dxScale, int dtScale);
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual DataProcessor2D *generate(std::vector<AtomicBlock2D *> atomicBlocks) const;
    virtual DotProcessorGenerator2D *clone() const;

private:
    DotProcessingFunctional2D *functional;
};

/// Easy instantiation of dotted data processor for a single lattice
template <typename T, template <typename U> class Descriptor>
struct DotProcessingFunctional2D_L : public DotProcessingFunctional2D {
    virtual void process(DotList2D const &dotList, BlockLattice2D<T, Descriptor> &lattice) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for a single scalar field
template <typename T>
struct DotProcessingFunctional2D_S : public DotProcessingFunctional2D {
    virtual void process(DotList2D const &dotList, ScalarField2D<T> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for a single tensor field
template <typename T, int nDim>
struct DotProcessingFunctional2D_T : public DotProcessingFunctional2D {
    virtual void process(DotList2D const &dotList, TensorField2D<T, nDim> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for a single n-tensor field
template <typename T>
struct DotProcessingFunctional2D_N : public DotProcessingFunctional2D {
    virtual void process(DotList2D const &dotList, NTensorField2D<T> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for lattice-lattice coupling
template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
struct DotProcessingFunctional2D_LL : public DotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, BlockLattice2D<T1, Descriptor1> &lattice1,
        BlockLattice2D<T2, Descriptor2> &lattice2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for ScalarField-ScalarField coupling
template <typename T1, typename T2>
struct DotProcessingFunctional2D_SS : public DotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, ScalarField2D<T1> &field1, ScalarField2D<T2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for TensorField-TensorField coupling
template <typename T1, int nDim1, typename T2, int nDim2>
struct DotProcessingFunctional2D_TT : public DotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, TensorField2D<T1, nDim1> &field1,
        TensorField2D<T2, nDim2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for ScalarField-TensorField coupling
template <typename T1, typename T2, int nDim>
struct DotProcessingFunctional2D_ST : public DotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, ScalarField2D<T1> &field1, TensorField2D<T2, nDim> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for NTensorField-NTensorField coupling
template <typename T1, typename T2>
struct DotProcessingFunctional2D_NN : public DotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, NTensorField2D<T1> &field1, NTensorField2D<T2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for ScalarField-NTensorField coupling
template <typename T1, typename T2>
struct DotProcessingFunctional2D_SN : public DotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, ScalarField2D<T1> &field1, NTensorField2D<T2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for Lattice-ScalarField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct DotProcessingFunctional2D_LS : public DotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, BlockLattice2D<T1, Descriptor> &lattice,
        ScalarField2D<T2> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for Lattice-TensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
struct DotProcessingFunctional2D_LT : public DotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, BlockLattice2D<T1, Descriptor> &lattice,
        TensorField2D<T2, nDim> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for Lattice-NTensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct DotProcessingFunctional2D_LN : public DotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, BlockLattice2D<T1, Descriptor> &lattice,
        NTensorField2D<T2> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple lattices
template <typename T, template <typename U> class Descriptor>
struct LatticeDotProcessingFunctional2D : public DotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, std::vector<BlockLattice2D<T, Descriptor> *> lattices) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple scalarFields
template <typename T>
struct ScalarFieldDotProcessingFunctional2D : public DotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, std::vector<ScalarField2D<T> *> scalarFields) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple tensorFields
template <typename T, int nDim>
struct TensorFieldDotProcessingFunctional2D : public DotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, std::vector<TensorField2D<T, nDim> *> tensorFields) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple n-tensorFields
template <typename T>
struct NTensorFieldDotProcessingFunctional2D : public DotProcessingFunctional2D {
    virtual void process(
        DotList2D const &dotList, std::vector<NTensorField2D<T> *> tensorFields) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks);
};

/* *************** All flavors of Bounded Box processing functionals ********* */

/// Easy instantiation of boxed data processor special boundary treatment (general case)
class BoundedBoxProcessingFunctional2D {
public:
    BoundedBoxProcessingFunctional2D();
    virtual ~BoundedBoxProcessingFunctional2D() { }
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) = 0;
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain,
        std::vector<AtomicBlock2D *> atomicBlocks) = 0;
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) = 0;
    virtual BlockDomain::DomainT appliesTo() const;
    /// Obsolete: replaced by setscale.
    virtual void rescale(double dxScale, double dtScale);
    virtual void setscale(int dxScale_, int dtScale_);
    /// Obsolete: replace by getTypeOfModification
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const = 0;
    virtual BoundedBoxProcessingFunctional2D *clone() const = 0;
    BoxProcessingFunctional2D *getBulkProcessor() const;
    BoxProcessingFunctional2D *getEdgeProcessor(int direction, int orientation) const;
    BoxProcessingFunctional2D *getCornerProcessor(int normalX, int normalY) const;
    void getGenerators(
        Box2D const &fullDomain, plint boundaryWidth_,
        std::vector<BoxProcessorGenerator2D *> &generators);
    int getDxScale() const;
    int getDtScale() const;

private:
    int dxScale, dtScale;

public:
    class BulkWrapperFunctional : public BoxProcessingFunctional2D {
    public:
        BulkWrapperFunctional(BoundedBoxProcessingFunctional2D *boundedFunctional_);
        BulkWrapperFunctional(BulkWrapperFunctional const &rhs);
        ~BulkWrapperFunctional();
        BulkWrapperFunctional &operator=(BulkWrapperFunctional const &rhs);
        virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        /// Obsolete: replaced by setscale.
        virtual void rescale(double dxScale, double dtScale);
        virtual void setscale(int dxScale_, int dtScale_);
        virtual void getModificationPattern(std::vector<bool> &isWritten) const;
        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
        virtual BulkWrapperFunctional *clone() const;

    private:
        BoundedBoxProcessingFunctional2D *boundedFunctional;
    };
    class EdgeWrapperFunctional : public BoxProcessingFunctional2D {
    public:
        EdgeWrapperFunctional(
            BoundedBoxProcessingFunctional2D *boundedFunctional_, int direction_, int orientation_);
        EdgeWrapperFunctional(EdgeWrapperFunctional const &rhs);
        ~EdgeWrapperFunctional();
        EdgeWrapperFunctional &operator=(EdgeWrapperFunctional const &rhs);
        virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        /// Obsolete: replaced by setscale.
        virtual void rescale(double dxScale, double dtScale);
        virtual void setscale(int dxScale_, int dtScale_);
        virtual void getModificationPattern(std::vector<bool> &isWritten) const;
        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
        virtual EdgeWrapperFunctional *clone() const;

    private:
        BoundedBoxProcessingFunctional2D *boundedFunctional;
        int direction, orientation;
    };
    class CornerWrapperFunctional : public BoxProcessingFunctional2D {
    public:
        CornerWrapperFunctional(
            BoundedBoxProcessingFunctional2D *boundedFunctional_, int normalX_, int normalY_);
        CornerWrapperFunctional(CornerWrapperFunctional const &rhs);
        ~CornerWrapperFunctional();
        CornerWrapperFunctional &operator=(CornerWrapperFunctional const &rhs);
        virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        /// Obsolete: replaced by setscale.
        virtual void rescale(double dxScale, double dtScale);
        virtual void setscale(int dxScale_, int dtScale_);
        virtual void getModificationPattern(std::vector<bool> &isWritten) const;
        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
        virtual CornerWrapperFunctional *clone() const;

    private:
        BoundedBoxProcessingFunctional2D *boundedFunctional;
        int normalX, normalY;
    };
};

/// Easy instantiation of bounded boxed data processor for a single lattice
template <typename T, template <typename U> class Descriptor>
struct BoundedBoxProcessingFunctional2D_L : public BoundedBoxProcessingFunctional2D {
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
struct BoundedBoxProcessingFunctional2D_S : public BoundedBoxProcessingFunctional2D {
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
struct BoundedBoxProcessingFunctional2D_T : public BoundedBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(Box2D domain, TensorField2D<T, nDim> &lattice) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain, TensorField2D<T, nDim> &lattice) = 0;
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, TensorField2D<T, nDim> &lattice) = 0;
};

/// Easy instantiation of bounded boxed data processor for a single n-tensor field
template <typename T>
struct BoundedBoxProcessingFunctional2D_N : public BoundedBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(Box2D domain, NTensorField2D<T> &lattice) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain, NTensorField2D<T> &lattice) = 0;
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, NTensorField2D<T> &lattice) = 0;
};

/// Easy instantiation of bounded boxed data processor for a single n-tensor field with bool mask.
template <typename T>
struct BoundedMaskedBoxProcessingFunctional2D_N : public BoundedBoxProcessingFunctional2D {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks);

    virtual void processBulk(
        Box2D domain, NTensorField2D<T> &lattice, NTensorField2D<int> &mask) = 0;
    virtual void processEdge(
        int direction, int orientation, Box2D domain, NTensorField2D<T> &lattice,
        NTensorField2D<int> &mask) = 0;
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, NTensorField2D<T> &lattice,
        NTensorField2D<int> &mask) = 0;
};

/// Easy instantiation of bounded boxed data processor for lattice-lattice coupling
template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
struct BoundedBoxProcessingFunctional2D_LL : public BoundedBoxProcessingFunctional2D {
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
struct BoundedBoxProcessingFunctional2D_SS : public BoundedBoxProcessingFunctional2D {
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
struct BoundedBoxProcessingFunctional2D_TT : public BoundedBoxProcessingFunctional2D {
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
struct BoundedBoxProcessingFunctional2D_ST : public BoundedBoxProcessingFunctional2D {
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
struct BoundedBoxProcessingFunctional2D_NN : public BoundedBoxProcessingFunctional2D {
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
/// NTensorField-NTensorField coupling with bool mask.
template <typename T1, typename T2>
struct BoundedMaskedBoxProcessingFunctional2D_NN : public BoundedBoxProcessingFunctional2D {
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
/// ScalarField-TensorField coupling.
template <typename T1, typename T2>
struct BoundedBoxProcessingFunctional2D_SN : public BoundedBoxProcessingFunctional2D {
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
struct BoundedBoxProcessingFunctional2D_LS : public BoundedBoxProcessingFunctional2D {
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
struct BoundedBoxProcessingFunctional2D_LT : public BoundedBoxProcessingFunctional2D {
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
struct BoundedBoxProcessingFunctional2D_LN : public BoundedBoxProcessingFunctional2D {
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
struct BoundedMaskedBoxProcessingFunctional2D_LN : public BoundedBoxProcessingFunctional2D {
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
struct BoundedLatticeBoxProcessingFunctional2D : public BoundedBoxProcessingFunctional2D {
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
struct BoundedScalarFieldBoxProcessingFunctional2D : public BoundedBoxProcessingFunctional2D {
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
struct BoundedTensorFieldBoxProcessingFunctional2D : public BoundedBoxProcessingFunctional2D {
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
struct BoundedNTensorFieldBoxProcessingFunctional2D : public BoundedBoxProcessingFunctional2D {
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

#endif  // DATA_PROCESSING_FUNCTIONAL_2D_H
