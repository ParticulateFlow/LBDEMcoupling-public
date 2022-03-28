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
#ifndef DATA_PROCESSING_FUNCTIONAL_3D_H
#define DATA_PROCESSING_FUNCTIONAL_3D_H

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

/// Easy instantiation of boxed data processor (general case)
class BoxProcessingFunctional3D {
public:
    BoxProcessingFunctional3D();
    virtual ~BoxProcessingFunctional3D() { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks) = 0;
    virtual BlockDomain::DomainT appliesTo() const;
    /// Obsolete: replaced by setscale.
    virtual void rescale(double dxScale, double dtScale);
    virtual void setscale(int dxScale_, int dtScale_);
    /// Obsolete: replaced by getTypeOfModification.
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const = 0;
    virtual BoxProcessingFunctional3D *clone() const = 0;
    int getDxScale() const;
    int getDtScale() const;
    virtual void serialize(std::string &data) const;
    virtual void unserialize(std::string &data);
    virtual int getStaticId() const;

private:
    int dxScale, dtScale;
};

/// A Boxed data processor, automatically generated from a BoxProcessingFunctional3D
class BoxProcessor3D : public DataProcessor3D {
public:
    BoxProcessor3D(
        BoxProcessingFunctional3D *functional_, Box3D domain_,
        std::vector<AtomicBlock3D *> atomicBlocks_);
    BoxProcessor3D(BoxProcessor3D const &rhs);
    BoxProcessor3D &operator=(BoxProcessor3D const &rhs);
    ~BoxProcessor3D();
    Box3D getDomain() const;
    virtual void process();
    virtual BoxProcessor3D *clone() const;
    virtual int getStaticId() const;

private:
    BoxProcessingFunctional3D *functional;
    Box3D domain;
    std::vector<AtomicBlock3D *> atomicBlocks;
};

/// A Boxed data processor, automatically generated from a MultiBoxProcessingFunctional3D
class MultiBoxProcessor3D : public DataProcessor3D {
public:
    MultiBoxProcessor3D(
        BoxProcessingFunctional3D *functional_, std::vector<Box3D> const &domains_,
        std::vector<AtomicBlock3D *> atomicBlocks_);
    MultiBoxProcessor3D(MultiBoxProcessor3D const &rhs);
    MultiBoxProcessor3D &operator=(MultiBoxProcessor3D const &rhs);
    ~MultiBoxProcessor3D();
    std::vector<Box3D> const &getDomains() const;
    virtual void process();
    virtual MultiBoxProcessor3D *clone() const;
    virtual int getStaticId() const;

private:
    BoxProcessingFunctional3D *functional;
    std::vector<Box3D> domains;
    std::vector<AtomicBlock3D *> atomicBlocks;
};

/// An automatically created generator for the BoxProcessor3D
class BoxProcessorGenerator3D : public BoxedDataProcessorGenerator3D {
public:
    BoxProcessorGenerator3D(BoxProcessingFunctional3D *functional_, Box3D domain);
    ~BoxProcessorGenerator3D();
    BoxProcessorGenerator3D(BoxProcessorGenerator3D const &rhs);
    BoxProcessorGenerator3D &operator=(BoxProcessorGenerator3D const &rhs);
    virtual BlockDomain::DomainT appliesTo() const;
    /// Obsolete: replaced by setscale.
    virtual void rescale(double dxScale, double dtScale);
    virtual void setscale(int dxScale_, int dtScale_);
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual DataProcessor3D *generate(std::vector<AtomicBlock3D *> atomicBlocks) const;
    virtual BoxProcessorGenerator3D *clone() const;
    virtual void serialize(Box3D &domain, std::string &data) const;
    virtual int getStaticId() const;

private:
    BoxProcessingFunctional3D *functional;
};

/// An automatically created generator for the BoxProcessor3D
class MultiBoxProcessorGenerator3D : public MultiBoxedDataProcessorGenerator3D {
public:
    MultiBoxProcessorGenerator3D(
        BoxProcessingFunctional3D *functional_, std::vector<Box3D> const &domains);
    ~MultiBoxProcessorGenerator3D();
    MultiBoxProcessorGenerator3D(MultiBoxProcessorGenerator3D const &rhs);
    MultiBoxProcessorGenerator3D &operator=(MultiBoxProcessorGenerator3D const &rhs);
    virtual BlockDomain::DomainT appliesTo() const;
    /// Obsolete: replaced by setscale.
    virtual void rescale(double dxScale, double dtScale);
    virtual void setscale(int dxScale_, int dtScale_);
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual DataProcessor3D *generate(std::vector<AtomicBlock3D *> atomicBlocks) const;
    virtual MultiBoxProcessorGenerator3D *clone() const;
    virtual void serialize(Box3D &domain, std::string &data) const;
    virtual int getStaticId() const;

private:
    BoxProcessingFunctional3D *functional;
};

/// Easy instantiation of boxed data processor for a single lattice
template <typename T, template <typename U> class Descriptor>
struct BoxProcessingFunctional3D_L : public BoxProcessingFunctional3D {
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice) = 0;
    virtual BoxProcessingFunctional3D_L<T, Descriptor> *clone() const = 0;

    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single lattice with bool mask.
template <typename T, template <typename U> class Descriptor>
struct MaskedBoxProcessingFunctional3D_L : public BoxProcessingFunctional3D {
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, NTensorField3D<int> &mask) = 0;
    virtual MaskedBoxProcessingFunctional3D_L<T, Descriptor> *clone() const = 0;

    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single scalar field
template <typename T>
struct BoxProcessingFunctional3D_S : public BoxProcessingFunctional3D {
    virtual void process(Box3D domain, ScalarField3D<T> &field) = 0;
    virtual BoxProcessingFunctional3D_S<T> *clone() const = 0;

    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single tensor field
template <typename T, int nDim>
struct BoxProcessingFunctional3D_T : public BoxProcessingFunctional3D {
    virtual void process(Box3D domain, TensorField3D<T, nDim> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single n-tensor field
template <typename T>
struct BoxProcessingFunctional3D_N : public BoxProcessingFunctional3D {
    virtual void process(Box3D domain, NTensorField3D<T> &field) = 0;
    virtual BoxProcessingFunctional3D_N<T> *clone() const = 0;

    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single n-tensor field with bool mask.
template <typename T>
struct MaskedBoxProcessingFunctional3D_N : public BoxProcessingFunctional3D {
    virtual void process(Box3D domain, NTensorField3D<T> &field, NTensorField3D<int> &mask) = 0;
    virtual MaskedBoxProcessingFunctional3D_N<T> *clone() const = 0;

    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for lattice-lattice coupling
template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
struct BoxProcessingFunctional3D_LL : public BoxProcessingFunctional3D {
    virtual void process(
        Box3D domain, BlockLattice3D<T1, Descriptor1> &lattice1,
        BlockLattice3D<T2, Descriptor2> &lattice2) = 0;
    virtual BoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *clone() const = 0;

    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for ScalarField-ScalarField coupling
template <typename T1, typename T2>
struct BoxProcessingFunctional3D_SS : public BoxProcessingFunctional3D {
    virtual void process(Box3D domain, ScalarField3D<T1> &field1, ScalarField3D<T2> &field2) = 0;
    virtual BoxProcessingFunctional3D_SS<T1, T2> *clone() const = 0;

    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for TensorField-TensorField coupling
template <typename T1, int nDim1, typename T2, int nDim2>
struct BoxProcessingFunctional3D_TT : public BoxProcessingFunctional3D {
    virtual void process(
        Box3D domain, TensorField3D<T1, nDim1> &field1, TensorField3D<T2, nDim2> &field2) = 0;
    virtual BoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *clone() const = 0;

    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for ScalarField-TensorField coupling
template <typename T1, typename T2, int nDim>
struct BoxProcessingFunctional3D_ST : public BoxProcessingFunctional3D {
    virtual void process(
        Box3D domain, ScalarField3D<T1> &field1, TensorField3D<T2, nDim> &field2) = 0;
    virtual BoxProcessingFunctional3D_ST<T1, T2, nDim> *clone() const = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for NTensorField-NTensorField coupling
template <typename T1, typename T2>
struct BoxProcessingFunctional3D_NN : public BoxProcessingFunctional3D {
    virtual void process(Box3D domain, NTensorField3D<T1> &field1, NTensorField3D<T2> &field2) = 0;
    virtual BoxProcessingFunctional3D_NN<T1, T2> *clone() const = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for NTensorField-NTensorField coupling with bool
/// mask.
template <typename T1, typename T2>
struct MaskedBoxProcessingFunctional3D_NN : public BoxProcessingFunctional3D {
    virtual void process(
        Box3D domain, NTensorField3D<T1> &field1, NTensorField3D<T2> &field2,
        NTensorField3D<int> &mask) = 0;
    virtual MaskedBoxProcessingFunctional3D_NN<T1, T2> *clone() const = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for ScalarField-NTensorField coupling
template <typename T1, typename T2>
struct BoxProcessingFunctional3D_SN : public BoxProcessingFunctional3D {
    virtual void process(Box3D domain, ScalarField3D<T1> &field1, NTensorField3D<T2> &field2) = 0;
    virtual BoxProcessingFunctional3D_SN<T1, T2> *clone() const = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-ScalarField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct BoxProcessingFunctional3D_LS : public BoxProcessingFunctional3D {
    virtual void process(
        Box3D domain, BlockLattice3D<T1, Descriptor> &lattice, ScalarField3D<T2> &field) = 0;
    virtual BoxProcessingFunctional3D_LS<T1, Descriptor, T2> *clone() const = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-TensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
struct BoxProcessingFunctional3D_LT : public BoxProcessingFunctional3D {
    virtual void process(
        Box3D domain, BlockLattice3D<T1, Descriptor> &lattice, TensorField3D<T2, nDim> &field) = 0;
    virtual BoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *clone() const = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-NTensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct BoxProcessingFunctional3D_LN : public BoxProcessingFunctional3D {
    virtual void process(
        Box3D domain, BlockLattice3D<T1, Descriptor> &lattice, NTensorField3D<T2> &field) = 0;
    virtual BoxProcessingFunctional3D_LN<T1, Descriptor, T2> *clone() const = 0;

    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-NTensorField coupling with bool mask
template <typename T1, template <typename U> class Descriptor, typename T2>
struct MaskedBoxProcessingFunctional3D_LN : public BoxProcessingFunctional3D {
    virtual void process(
        Box3D domain, BlockLattice3D<T1, Descriptor> &lattice, NTensorField3D<T2> &field,
        NTensorField3D<int> &mask) = 0;
    virtual MaskedBoxProcessingFunctional3D_LN<T1, Descriptor, T2> *clone() const = 0;

    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple lattices
template <typename T, template <typename U> class Descriptor>
struct LatticeBoxProcessingFunctional3D : public BoxProcessingFunctional3D {
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T, Descriptor> *> lattices) = 0;
    virtual LatticeBoxProcessingFunctional3D<T, Descriptor> *clone() const = 0;

    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple scalarFields
template <typename T>
struct ScalarFieldBoxProcessingFunctional3D : public BoxProcessingFunctional3D {
    virtual void process(Box3D domain, std::vector<ScalarField3D<T> *> scalarFields) = 0;
    virtual ScalarFieldBoxProcessingFunctional3D<T> *clone() const = 0;

    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple tensorFields
template <typename T, int nDim>
struct TensorFieldBoxProcessingFunctional3D : public BoxProcessingFunctional3D {
    virtual void process(Box3D domain, std::vector<TensorField3D<T, nDim> *> tensorFields) = 0;
    virtual TensorFieldBoxProcessingFunctional3D<T, nDim> *clone() const = 0;

    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple n-tensorFields
template <typename T>
struct NTensorFieldBoxProcessingFunctional3D : public BoxProcessingFunctional3D {
    virtual void process(Box3D domain, std::vector<NTensorField3D<T> *> tensorFields) = 0;
    virtual NTensorFieldBoxProcessingFunctional3D<T> *clone() const = 0;

    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple n-tensorFields
template <typename T>
struct MaskedNTensorFieldBoxProcessingFunctional3D : public BoxProcessingFunctional3D {
    virtual void process(
        Box3D domain, std::vector<NTensorField3D<T> *> tensorFields, NTensorField3D<int> &mask) = 0;
    virtual MaskedNTensorFieldBoxProcessingFunctional3D<T> *clone() const = 0;

    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
};

/* *************** All flavors of Dot processing functionals ********* */

/// Easy instantiation of dotted data processor (general case)
struct DotProcessingFunctional3D {
    virtual ~DotProcessingFunctional3D() { }
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks) = 0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(double dxScale, double dtScale);
    virtual void setscale(int dxScale, int dtScale);
    /// Obsolete: replaced by getTypeOfModification
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const = 0;
    virtual DotProcessingFunctional3D *clone() const = 0;
};

/// A Dotted data processor, automatically generated from a DotProcessingFunctional3D
class DotProcessor3D : public DataProcessor3D {
public:
    DotProcessor3D(
        DotProcessingFunctional3D *functional_, DotList3D const &dotList_,
        std::vector<AtomicBlock3D *> atomicBlocks_);
    DotProcessor3D(DotProcessor3D const &rhs);
    DotProcessor3D &operator=(DotProcessor3D const &rhs);
    ~DotProcessor3D();
    virtual void process();
    virtual DotProcessor3D *clone() const;
    DotList3D const &getDotList() const;

private:
    DotProcessingFunctional3D *functional;
    DotList3D dotList;
    std::vector<AtomicBlock3D *> atomicBlocks;
};

/// An automatically created generator for the DotProcessor3D
class DotProcessorGenerator3D : public DottedDataProcessorGenerator3D {
public:
    DotProcessorGenerator3D(DotProcessingFunctional3D *functional_, DotList3D const &dotList);
    ~DotProcessorGenerator3D();
    DotProcessorGenerator3D(DotProcessorGenerator3D const &rhs);
    DotProcessorGenerator3D &operator=(DotProcessorGenerator3D const &rhs);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(double dxScale, double dtScale);
    virtual void setscale(int dxScale, int dtScale);
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual DataProcessor3D *generate(std::vector<AtomicBlock3D *> atomicBlocks) const;
    virtual DotProcessorGenerator3D *clone() const;

private:
    DotProcessingFunctional3D *functional;
};

/// Easy instantiation of dotted data processor for a single lattice
template <typename T, template <typename U> class Descriptor>
struct DotProcessingFunctional3D_L : public DotProcessingFunctional3D {
    virtual void process(DotList3D const &dotList, BlockLattice3D<T, Descriptor> &lattice) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for a single scalar field
template <typename T>
struct DotProcessingFunctional3D_S : public DotProcessingFunctional3D {
    virtual void process(DotList3D const &dotList, ScalarField3D<T> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for a single tensor field
template <typename T, int nDim>
struct DotProcessingFunctional3D_T : public DotProcessingFunctional3D {
    virtual void process(DotList3D const &dotList, TensorField3D<T, nDim> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for a single n-tensor field
template <typename T>
struct DotProcessingFunctional3D_N : public DotProcessingFunctional3D {
    virtual void process(DotList3D const &dotList, NTensorField3D<T> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for lattice-lattice coupling
template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
struct DotProcessingFunctional3D_LL : public DotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, BlockLattice3D<T1, Descriptor1> &lattice1,
        BlockLattice3D<T2, Descriptor2> &lattice2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for ScalarField-ScalarField coupling
template <typename T1, typename T2>
struct DotProcessingFunctional3D_SS : public DotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, ScalarField3D<T1> &field1, ScalarField3D<T2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for TensorField-TensorField coupling
template <typename T1, int nDim1, typename T2, int nDim2>
struct DotProcessingFunctional3D_TT : public DotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, TensorField3D<T1, nDim1> &field1,
        TensorField3D<T2, nDim2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for ScalarField-TensorField coupling
template <typename T1, typename T2, int nDim>
struct DotProcessingFunctional3D_ST : public DotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, ScalarField3D<T1> &field1, TensorField3D<T2, nDim> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for NTensorField-NTensorField coupling
template <typename T1, typename T2>
struct DotProcessingFunctional3D_NN : public DotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, NTensorField3D<T1> &field1, NTensorField3D<T2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for ScalarField-NTensorField coupling
template <typename T1, typename T2>
struct DotProcessingFunctional3D_SN : public DotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, ScalarField3D<T1> &field1, NTensorField3D<T2> &field2) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for Lattice-ScalarField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct DotProcessingFunctional3D_LS : public DotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, BlockLattice3D<T1, Descriptor> &lattice,
        ScalarField3D<T2> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for Lattice-TensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
struct DotProcessingFunctional3D_LT : public DotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, BlockLattice3D<T1, Descriptor> &lattice,
        TensorField3D<T2, nDim> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for Lattice-NTensorField coupling
template <typename T1, template <typename U> class Descriptor, typename T2>
struct DotProcessingFunctional3D_LN : public DotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, BlockLattice3D<T1, Descriptor> &lattice,
        NTensorField3D<T2> &field) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple lattices
template <typename T, template <typename U> class Descriptor>
struct LatticeDotProcessingFunctional3D : public DotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, std::vector<BlockLattice3D<T, Descriptor> *> lattices) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple scalarFields
template <typename T>
struct ScalarFieldDotProcessingFunctional3D : public DotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, std::vector<ScalarField3D<T> *> scalarFields) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple tensorFields
template <typename T, int nDim>
struct TensorFieldDotProcessingFunctional3D : public DotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, std::vector<TensorField3D<T, nDim> *> tensorFields) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple n-tensorFields
template <typename T>
struct NTensorFieldDotProcessingFunctional3D : public DotProcessingFunctional3D {
    virtual void process(
        DotList3D const &dotList, std::vector<NTensorField3D<T> *> tensorFields) = 0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(
        DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks);
};

/* *************** All flavors of Bounded Box processing functionals ********* */

/// Easy instantiation of boxed data processor special boundary treatment (general case)
class BoundedBoxProcessingFunctional3D {
public:
    BoundedBoxProcessingFunctional3D();
    virtual ~BoundedBoxProcessingFunctional3D() { }
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
    /// Obsolete: replaced by setscale.
    virtual void rescale(double dxScale, double dtScale);
    virtual void setscale(int dxScale_, int dtScale_);
    /// Obsolete: replace by getTypeOfModification
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const = 0;
    virtual BoundedBoxProcessingFunctional3D *clone() const = 0;
    BoxProcessingFunctional3D *getBulkProcessor() const;
    BoxProcessingFunctional3D *getPlaneProcessor(int direction, int orientation) const;
    BoxProcessingFunctional3D *getEdgeProcessor(int plane, int normal1, int normal2) const;
    BoxProcessingFunctional3D *getCornerProcessor(int normalX, int normalY, int normalZ) const;
    void getGenerators(
        Box3D const &fullDomain, plint boundaryWidth_,
        std::vector<BoxProcessorGenerator3D *> &generators);
    int getDxScale() const;
    int getDtScale() const;

private:
    int dxScale, dtScale;

public:
    class BulkWrapperFunctional : public BoxProcessingFunctional3D {
    public:
        BulkWrapperFunctional(BoundedBoxProcessingFunctional3D *boundedFunctional_);
        BulkWrapperFunctional(BulkWrapperFunctional const &rhs);
        ~BulkWrapperFunctional();
        BulkWrapperFunctional &operator=(BulkWrapperFunctional const &rhs);
        virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        virtual void rescale(double dxScale, double dtScale);
        virtual void setscale(int dxScale_, int dtScale_);
        virtual void getModificationPattern(std::vector<bool> &isWritten) const;
        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
        virtual BulkWrapperFunctional *clone() const;

    private:
        BoundedBoxProcessingFunctional3D *boundedFunctional;
    };
    class PlaneWrapperFunctional : public BoxProcessingFunctional3D {
    public:
        PlaneWrapperFunctional(
            BoundedBoxProcessingFunctional3D *boundedFunctional_, int direction_, int orientation_);
        PlaneWrapperFunctional(PlaneWrapperFunctional const &rhs);
        ~PlaneWrapperFunctional();
        PlaneWrapperFunctional &operator=(PlaneWrapperFunctional const &rhs);
        virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        /// Obsolete: replaced by setscale.
        virtual void rescale(double dxScale, double dtScale);
        virtual void setscale(int dxScale_, int dtScale_);
        virtual void getModificationPattern(std::vector<bool> &isWritten) const;
        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
        virtual PlaneWrapperFunctional *clone() const;

    private:
        BoundedBoxProcessingFunctional3D *boundedFunctional;
        int direction, orientation;
    };
    class EdgeWrapperFunctional : public BoxProcessingFunctional3D {
    public:
        EdgeWrapperFunctional(
            BoundedBoxProcessingFunctional3D *boundedFunctional_, int plane_, int normal1_,
            int normal2_);
        EdgeWrapperFunctional(EdgeWrapperFunctional const &rhs);
        ~EdgeWrapperFunctional();
        EdgeWrapperFunctional &operator=(EdgeWrapperFunctional const &rhs);
        virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        /// Obsolete: replaced by setscale.
        virtual void rescale(double dxScale, double dtScale);
        virtual void setscale(int dxScale_, int dtScale_);
        virtual void getModificationPattern(std::vector<bool> &isWritten) const;
        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
        virtual EdgeWrapperFunctional *clone() const;

    private:
        BoundedBoxProcessingFunctional3D *boundedFunctional;
        int plane, normal1, normal2;
    };
    class CornerWrapperFunctional : public BoxProcessingFunctional3D {
    public:
        CornerWrapperFunctional(
            BoundedBoxProcessingFunctional3D *boundedFunctional_, int normalX_, int normalY_,
            int normalZ_);
        CornerWrapperFunctional(CornerWrapperFunctional const &rhs);
        ~CornerWrapperFunctional();
        CornerWrapperFunctional &operator=(CornerWrapperFunctional const &rhs);
        virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        /// Obsolete: replaced by setscale.
        virtual void rescale(double dxScale, double dtScale);
        virtual void setscale(int dxScale_, int dtScale_);
        virtual void getModificationPattern(std::vector<bool> &isWritten) const;
        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
        virtual CornerWrapperFunctional *clone() const;

    private:
        BoundedBoxProcessingFunctional3D *boundedFunctional;
        int normalX, normalY, normalZ;
    };
};

/// Easy instantiation of bounded boxed data processor for a single lattice
template <typename T, template <typename U> class Descriptor>
struct BoundedBoxProcessingFunctional3D_L : public BoundedBoxProcessingFunctional3D {
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
struct BoundedBoxProcessingFunctional3D_S : public BoundedBoxProcessingFunctional3D {
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
struct BoundedBoxProcessingFunctional3D_T : public BoundedBoxProcessingFunctional3D {
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
        int direction, int orientation, Box3D domain, TensorField3D<T, nDim> &field) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain, TensorField3D<T, nDim> &field) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain, TensorField3D<T, nDim> &field) = 0;
};

/// Easy instantiation of bounded boxed data processor for a single n-tensor field
template <typename T>
struct BoundedBoxProcessingFunctional3D_N : public BoundedBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual void processBulk(Box3D domain, NTensorField3D<T> &field) = 0;
    virtual void processPlane(
        int direction, int orientation, Box3D domain, NTensorField3D<T> &field) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain, NTensorField3D<T> &field) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain, NTensorField3D<T> &field) = 0;
};

/// Easy instantiation of bounded boxed data processor for a single n-tensor field with bool mask.
template <typename T>
struct BoundedMaskedBoxProcessingFunctional3D_N : public BoundedBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual void processBulk(Box3D domain, NTensorField3D<T> &field, NTensorField3D<int> &mask) = 0;
    virtual void processPlane(
        int direction, int orientation, Box3D domain, NTensorField3D<T> &field,
        NTensorField3D<int> &mask) = 0;
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain, NTensorField3D<T> &field,
        NTensorField3D<int> &mask) = 0;
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain, NTensorField3D<T> &field,
        NTensorField3D<int> &mask) = 0;
};

/// Easy instantiation of bounded boxed data processor for lattice-lattice coupling
template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
struct BoundedBoxProcessingFunctional3D_LL : public BoundedBoxProcessingFunctional3D {
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
struct BoundedBoxProcessingFunctional3D_SS : public BoundedBoxProcessingFunctional3D {
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
struct BoundedBoxProcessingFunctional3D_TT : public BoundedBoxProcessingFunctional3D {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks);

    virtual BoundedBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *clone() const = 0;

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
struct BoundedBoxProcessingFunctional3D_ST : public BoundedBoxProcessingFunctional3D {
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
struct BoundedBoxProcessingFunctional3D_NN : public BoundedBoxProcessingFunctional3D {
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

/// Easy instantiation of bounded boxed data processor for NTensorField-NTensorField coupling with
/// bool mask.
template <typename T1, typename T2>
struct BoundedMaskedBoxProcessingFunctional3D_NN : public BoundedBoxProcessingFunctional3D {
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
struct BoundedBoxProcessingFunctional3D_SN : public BoundedBoxProcessingFunctional3D {
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
struct BoundedBoxProcessingFunctional3D_LS : public BoundedBoxProcessingFunctional3D {
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
struct BoundedBoxProcessingFunctional3D_LT : public BoundedBoxProcessingFunctional3D {
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
struct BoundedBoxProcessingFunctional3D_LN : public BoundedBoxProcessingFunctional3D {
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
struct BoundedMaskedBoxProcessingFunctional3D_LN : public BoundedBoxProcessingFunctional3D {
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
struct BoundedLatticeBoxProcessingFunctional3D : public BoundedBoxProcessingFunctional3D {
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
struct BoundedScalarFieldBoxProcessingFunctional3D : public BoundedBoxProcessingFunctional3D {
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
struct BoundedTensorFieldBoxProcessingFunctional3D : public BoundedBoxProcessingFunctional3D {
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
struct BoundedNTensorFieldBoxProcessingFunctional3D : public BoundedBoxProcessingFunctional3D {
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

#endif  // DATA_PROCESSING_FUNCTIONAL_3D_H
