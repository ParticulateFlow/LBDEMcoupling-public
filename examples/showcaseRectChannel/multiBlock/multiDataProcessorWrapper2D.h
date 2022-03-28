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
 *
 * This file offers explicit wrappers to execute data-processors on atomic-
 * blocks, or to integrate the data-processors into the atomic-blocks. A
 * wrapper is offered for explicit choices of data types, for one or two
 * block arguments. These wrappers are luxury: it would have been sufficient
 * to have a generic wrapper which takes a single atomic-block as argument,
 * then a generic wrapper for two atomic-blocks. The advantage of having an
 * explicit wrapper for each choice of data types is type safety. It is
 * guaranteed that the type of atomic-blocks fits the type of the data
 * processor, which avoids awkward runtime errors.
 */

#ifndef MULTI_DATA_PROCESSOR_WRAPPER_2D_H
#define MULTI_DATA_PROCESSOR_WRAPPER_2D_H

#include "atomicBlock/dataProcessingFunctional2D.h"

namespace plb {

class MultiBlock2D;
template <typename T, template <typename U> class Descriptor>
class MultiBlockLattice2D;
template <typename T>
class MultiScalarField2D;
template <typename T, int nDim>
class MultiTensorField2D;
template <typename T>
class MultiNTensorField2D;

/* *************** Generic wrappers, boxed functionals ********************** */

/// This is the most general wrapper for applying a 2D boxed data
/// functional. Use this if none of the more specific wrappers works.
void applyProcessingFunctional(
    BoxProcessingFunctional2D *functional, Box2D domain, std::vector<MultiBlock2D *> multiBlocks);

/// This is the most general wrapper for integrating a 2D boxed data
/// functional. Use this if none of the more specific wrappers works.
void integrateProcessingFunctional(
    BoxProcessingFunctional2D *functional, Box2D domain, std::vector<MultiBlock2D *> multiBlocks,
    plint level = 0);

/// Apply a functional on a sequence of block-lattices. If the number
/// of lattices is 1 or 2, you should prefer the _L and _LL version
/// of the functional.
template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    LatticeBoxProcessingFunctional2D<T, Descriptor> *functional, Box2D domain,
    std::vector<MultiBlockLattice2D<T, Descriptor> *> lattices);

/// Integrate a functional into a sequence of block-lattices. If the
/// number of lattices is 1 or 2, you should prefer the _L and _LL version
/// of the functional.
template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    LatticeBoxProcessingFunctional2D<T, Descriptor> *functional, Box2D domain,
    std::vector<MultiBlockLattice2D<T, Descriptor> *> lattices, plint level = 0);

/// Apply a functional on a sequence of scalar-fields. If the number
/// of lattices is 1 or 2, you should prefer the _S and _SS version
/// of the functional.
template <typename T>
void applyProcessingFunctional(
    ScalarFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiScalarField2D<T> *> fields);

/// Integrate a functional into a sequence of scalar-fields. If the
/// number of lattices is 1 or 2, you should prefer the _S and _SS version
/// of the functional.
template <typename T>
void integrateProcessingFunctional(
    ScalarFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiScalarField2D<T> *> fields, plint level = 0);

/// Apply a functional on a sequence of tensor-fields. If the number
/// of lattices is 1 or 2, you should prefer the _T and _TT version
/// of the functional.
template <typename T, int nDim>
void applyProcessingFunctional(
    TensorFieldBoxProcessingFunctional2D<T, nDim> *functional, Box2D domain,
    std::vector<MultiTensorField2D<T, nDim> *> fields);

/// Integrate a functional into a sequence of tensor-fields. If the
/// number of lattices is 1 or 2, you should prefer the _T and _TT version
/// of the functional.
template <typename T, int nDim>
void integrateProcessingFunctional(
    TensorFieldBoxProcessingFunctional2D<T, nDim> *functional, Box2D domain,
    std::vector<MultiTensorField2D<T, nDim> *> fields, plint level = 0);

/// Apply a functional on a sequence of tensor-fields. If the number
/// of lattices is 1 or 2, you should prefer the _N and _NN version
/// of the functional.
template <typename T>
void applyProcessingFunctional(
    NTensorFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiNTensorField2D<T> *> fields);

/// Apply a functional on a sequence of n/tensor-fields, with bool mask.
/// If the number of lattices is 1 or 2, you should prefer the _N and _NN
/// version of the functional.
template <typename T>
void applyProcessingFunctional(
    MaskedNTensorFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiNTensorField2D<T> *> fields, MultiNTensorField2D<int> &mask);

/// Integrate a functional into a sequence of tensor-fields. If the
/// number of lattices is 1 or 2, you should prefer the _N and _NN version
/// of the functional.
template <typename T>
void integrateProcessingFunctional(
    NTensorFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiNTensorField2D<T> *> fields, plint level = 0);

/// Integrate a functional into a sequence of tensor-fields, with mask. If the
/// number of lattices is 1 or 2, you should prefer the _N and _NN version
/// of the functional.
template <typename T>
void integrateProcessingFunctional(
    MaskedNTensorFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiNTensorField2D<T> *> fields, MultiNTensorField2D<int> &mask, plint level = 0);

/* *************** Typed wrappers with a single argument, boxed functionals * */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_L<T, Descriptor> *functional, Box2D domain,
    MultiBlockLattice2D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_L<T, Descriptor> *functional, Box2D domain,
    MultiBlockLattice2D<T, Descriptor> &lattice, plint level = 0);

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    MaskedBoxProcessingFunctional2D_L<T, Descriptor> *functional, Box2D domain,
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiNTensorField2D<int> &mask);

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    MaskedBoxProcessingFunctional2D_L<T, Descriptor> *functional, Box2D domain,
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiNTensorField2D<int> &mask, plint level = 0);

template <typename T>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_S<T> *functional, Box2D domain, MultiScalarField2D<T> &field);
template <typename T>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_S<T> *functional, Box2D domain, MultiScalarField2D<T> &field,
    plint level = 0);

template <typename T, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_T<T, nDim> *functional, Box2D domain,
    MultiTensorField2D<T, nDim> &field);
template <typename T, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_T<T, nDim> *functional, Box2D domain,
    MultiTensorField2D<T, nDim> &field, plint level = 0);

template <typename T>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_N<T> *functional, Box2D domain, MultiNTensorField2D<T> &field);
template <typename T>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_N<T> *functional, Box2D domain, MultiNTensorField2D<T> &field,
    plint level = 0);

template <typename T>
void applyProcessingFunctional(
    MaskedBoxProcessingFunctional2D_N<T> *functional, Box2D domain, MultiNTensorField2D<T> &field,
    MultiNTensorField2D<int> &mask);
template <typename T>
void integrateProcessingFunctional(
    MaskedBoxProcessingFunctional2D_N<T> *functional, Box2D domain, MultiNTensorField2D<T> &field,
    MultiNTensorField2D<int> &mask, plint level = 0);

/* *************** Typed wrappers with two arguments ************************ */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor1> &lattice1, MultiBlockLattice2D<T2, Descriptor2> &lattice2);

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor1> &lattice1, MultiBlockLattice2D<T2, Descriptor2> &lattice2,
    plint level = 0);

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_SS<T1, T2> *functional, Box2D domain, MultiScalarField2D<T1> &field1,
    MultiScalarField2D<T2> &field2);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_SS<T1, T2> *functional, Box2D domain, MultiScalarField2D<T1> &field1,
    MultiScalarField2D<T2> &field2, plint level = 0);

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> *functional, Box2D domain,
    MultiTensorField2D<T1, nDim1> &field1, MultiTensorField2D<T2, nDim2> &field2);

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> *functional, Box2D domain,
    MultiTensorField2D<T1, nDim1> &field1, MultiTensorField2D<T2, nDim2> &field2, plint level = 0);

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_NN<T1, T2> *functional, Box2D domain, MultiNTensorField2D<T1> &field1,
    MultiNTensorField2D<T2> &field2);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_NN<T1, T2> *functional, Box2D domain, MultiNTensorField2D<T1> &field1,
    MultiNTensorField2D<T2> &field2, plint level = 0);

template <typename T1, typename T2>
void applyProcessingFunctional(
    MaskedBoxProcessingFunctional2D_NN<T1, T2> *functional, Box2D domain,
    MultiNTensorField2D<T1> &field1, MultiNTensorField2D<T2> &field2,
    MultiNTensorField2D<int> &mask);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    MaskedBoxProcessingFunctional2D_NN<T1, T2> *functional, Box2D domain,
    MultiNTensorField2D<T1> &field1, MultiNTensorField2D<T2> &field2,
    MultiNTensorField2D<int> &mask, plint level = 0);

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_ST<T1, T2, nDim> *functional, Box2D domain,
    MultiScalarField2D<T1> &field1, MultiTensorField2D<T2, nDim> &field2);

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_ST<T1, T2, nDim> *functional, Box2D domain,
    MultiScalarField2D<T1> &field1, MultiTensorField2D<T2, nDim> &field2, plint level = 0);

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_SN<T1, T2> *functional, Box2D domain, MultiScalarField2D<T1> &field1,
    MultiNTensorField2D<T2> &field2);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_SN<T1, T2> *functional, Box2D domain, MultiScalarField2D<T1> &field1,
    MultiNTensorField2D<T2> &field2, plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_LS<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiScalarField2D<T2> &field);

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_LS<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiScalarField2D<T2> &field, plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiTensorField2D<T2, nDim> &field);

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiTensorField2D<T2, nDim> &field,
    plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_LN<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field);

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_LN<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field, plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    MaskedBoxProcessingFunctional2D_LN<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field,
    MultiNTensorField2D<int> &mask);

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    MaskedBoxProcessingFunctional2D_LN<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field,
    MultiNTensorField2D<int> &mask, plint level = 0);

/* *************** Generic wrappers, dotted functionals ********************* */

void applyProcessingFunctional(
    DotProcessingFunctional2D *functional, DotList2D const &dotList,
    std::vector<MultiBlock2D *> multiBlocks);

void integrateProcessingFunctional(
    DotProcessingFunctional2D *functional, DotList2D const &dotList,
    std::vector<MultiBlock2D *> multiBlocks, plint level = 0);

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    LatticeDotProcessingFunctional2D<T, Descriptor> *functional, DotList2D const &dotList,
    std::vector<MultiBlockLattice2D<T, Descriptor> *> lattices);

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    LatticeDotProcessingFunctional2D<T, Descriptor> *functional, DotList2D const &dotList,
    std::vector<MultiBlockLattice2D<T, Descriptor> *> lattices, plint level = 0);

template <typename T>
void applyProcessingFunctional(
    ScalarFieldDotProcessingFunctional2D<T> *functional, DotList2D const &dotList,
    std::vector<MultiScalarField2D<T> *> fields);

template <typename T>
void integrateProcessingFunctional(
    ScalarFieldDotProcessingFunctional2D<T> *functional, DotList2D const &dotList,
    std::vector<MultiScalarField2D<T> *> fields, plint level = 0);

template <typename T, int nDim>
void applyProcessingFunctional(
    TensorFieldDotProcessingFunctional2D<T, nDim> *functional, DotList2D const &dotList,
    std::vector<MultiTensorField2D<T, nDim> *> fields);

template <typename T, int nDim>
void integrateProcessingFunctional(
    TensorFieldDotProcessingFunctional2D<T, nDim> *functional, DotList2D const &dotList,
    std::vector<MultiTensorField2D<T, nDim> *> fields, plint level = 0);

template <typename T>
void applyProcessingFunctional(
    NTensorFieldDotProcessingFunctional2D<T> *functional, DotList2D const &dotList,
    std::vector<MultiNTensorField2D<T> *> fields);

template <typename T>
void integrateProcessingFunctional(
    NTensorFieldDotProcessingFunctional2D<T> *functional, DotList2D const &dotList,
    std::vector<MultiNTensorField2D<T> *> fields, plint level = 0);

/* *************** Typed wrappers with a single argument, dotted functionals* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    DotProcessingFunctional2D_L<T, Descriptor> *functional, DotList2D const &dotList,
    MultiBlockLattice2D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_L<T, Descriptor> *functional, DotList2D const &dotList,
    MultiBlockLattice2D<T, Descriptor> &lattice, plint level = 0);

template <typename T>
void applyProcessingFunctional(
    DotProcessingFunctional2D_S<T> *functional, DotList2D const &dotList,
    MultiScalarField2D<T> &field);

template <typename T>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_S<T> *functional, DotList2D const &dotList,
    MultiScalarField2D<T> &field, plint level = 0);

template <typename T, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional2D_T<T, nDim> *functional, DotList2D const &dotList,
    MultiTensorField2D<T, nDim> &field);

template <typename T, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_T<T, nDim> *functional, DotList2D const &dotList,
    MultiTensorField2D<T, nDim> &field, plint level = 0);

template <typename T>
void applyProcessingFunctional(
    DotProcessingFunctional2D_N<T> *functional, DotList2D const &dotList,
    MultiNTensorField2D<T> &field);

template <typename T>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_N<T> *functional, DotList2D const &dotList,
    MultiNTensorField2D<T> &field, plint level = 0);

/* *************** Typed wrappers with two arguments, dotted functionals* *****/

/// Easy instantiation of dotted data processor for lattice-lattice coupling
template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    DotProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> *functional,
    DotList2D const &dotList, MultiBlockLattice2D<T1, Descriptor1> &lattice1,
    MultiBlockLattice2D<T2, Descriptor2> &lattice2);

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> *functional,
    DotList2D const &dotList, MultiBlockLattice2D<T1, Descriptor1> &lattice1,
    MultiBlockLattice2D<T2, Descriptor2> &lattice2, plint level = 0);

template <typename T1, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional2D_SS<T1, T2> *functional, DotList2D const &dotList,
    MultiScalarField2D<T1> &field1, MultiScalarField2D<T2> &field2);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_SS<T1, T2> *functional, DotList2D const &dotList,
    MultiScalarField2D<T1> &field1, MultiScalarField2D<T2> &field2, plint level = 0);

/// Easy instantiation of dotted data processor for MultiTensorField-MultiTensorField coupling
template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    DotProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> *functional, DotList2D const &dotList,
    MultiTensorField2D<T1, nDim1> &field1, MultiTensorField2D<T2, nDim2> &field2);

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> *functional, DotList2D const &dotList,
    MultiTensorField2D<T1, nDim1> &field1, MultiTensorField2D<T2, nDim2> &field2, plint level = 0);

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    DotProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> *functional, DotList2D const &dotList,
    MultiTensorField2D<T1, nDim1> &field1, MultiTensorField2D<T2, nDim2> &field2);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_NN<T1, T2> *functional, DotList2D const &dotList,
    MultiNTensorField2D<T1> &field1, MultiNTensorField2D<T2> &field2, plint level = 0);

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional2D_ST<T1, T2, nDim> *functional, DotList2D const &dotList,
    MultiScalarField2D<T1> &field1, MultiTensorField2D<T2, nDim> &field2);

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_ST<T1, T2, nDim> *functional, DotList2D const &dotList,
    MultiScalarField2D<T1> &field1, MultiTensorField2D<T2, nDim> &field2, plint level = 0);

template <typename T1, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional2D_SN<T1, T2> *functional, DotList2D const &dotList,
    MultiScalarField2D<T1> &field1, MultiNTensorField2D<T2> &field2);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_SN<T1, T2> *functional, DotList2D const &dotList,
    MultiScalarField2D<T1> &field1, MultiNTensorField2D<T2> &field2, plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional2D_LS<T1, Descriptor, T2> *functional, DotList2D const &dotList,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiScalarField2D<T2> &field);

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_LS<T1, Descriptor, T2> *functional, DotList2D const &dotList,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiScalarField2D<T2> &field, plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> *functional, DotList2D const &dotList,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiTensorField2D<T2, nDim> &field);

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> *functional, DotList2D const &dotList,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiTensorField2D<T2, nDim> &field,
    plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional2D_LN<T1, Descriptor, T2> *functional, DotList2D const &dotList,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field);

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_LN<T1, Descriptor, T2> *functional, DotList2D const &dotList,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field, plint level = 0);

/* *************** Generic wrappers, bounded and boxed functionals ********** */

void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D *functional, Box2D domain,
    std::vector<MultiBlock2D *> multiBlocks, plint boundaryWidth);

void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D *functional, Box2D domain,
    std::vector<MultiBlock2D *> multiBlocks, plint boundaryWidth, plint level = 0);

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoundedLatticeBoxProcessingFunctional2D<T, Descriptor> *functional, Box2D domain,
    std::vector<MultiBlockLattice2D<T, Descriptor> *> lattices,
    plint boundaryWidth = Descriptor<T>::vicinity);

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoundedLatticeBoxProcessingFunctional2D<T, Descriptor> *functional, Box2D domain,
    std::vector<MultiBlockLattice2D<T, Descriptor> *> lattices,
    plint boundaryWidth = Descriptor<T>::vicinity, plint level = 0);

template <typename T>
void applyProcessingFunctional(
    BoundedScalarFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiScalarField2D<T> *> fields, plint boundaryWidth);

template <typename T>
void integrateProcessingFunctional(
    BoundedScalarFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiScalarField2D<T> *> fields, plint boundaryWidth);

template <typename T, int nDim>
void applyProcessingFunctional(
    BoundedTensorFieldBoxProcessingFunctional2D<T, nDim> *functional, Box2D domain,
    std::vector<MultiTensorField2D<T, nDim> *> fields, plint boundaryWidth);

template <typename T, int nDim>
void integrateProcessingFunctional(
    BoundedTensorFieldBoxProcessingFunctional2D<T, nDim> *functional, Box2D domain,
    std::vector<MultiTensorField2D<T, nDim> *> fields, plint boundaryWidth);

template <typename T>
void applyProcessingFunctional(
    BoundedNTensorFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiNTensorField2D<T> *> fields, plint boundaryWidth);

template <typename T>
void integrateProcessingFunctional(
    BoundedNTensorFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiNTensorField2D<T> *> fields, plint boundaryWidth);

/* ** Typed wrappers with a single argument, bounded and boxed functionals * */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_L<T, Descriptor> *functional, Box2D domain,
    MultiBlockLattice2D<T, Descriptor> &lattice, plint boundaryWidth = Descriptor<T>::vicinity);

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_L<T, Descriptor> *functional, Box2D domain,
    MultiBlockLattice2D<T, Descriptor> &lattice, plint boundaryWidth = Descriptor<T>::vicinity,
    plint level = 0);

template <typename T>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_S<T> *functional, Box2D domain, MultiScalarField2D<T> &field,
    plint boundaryWidth);

template <typename T>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_S<T> *functional, Box2D domain, MultiScalarField2D<T> &field,
    plint boundaryWidth, plint level = 0);

template <typename T, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_T<T, nDim> *functional, Box2D domain,
    MultiTensorField2D<T, nDim> &field, plint boundaryWidth);

template <typename T, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_T<T, nDim> *functional, Box2D domain,
    MultiTensorField2D<T, nDim> &field, plint boundaryWidth, plint level = 0);

template <typename T>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_N<T> *functional, Box2D domain, MultiNTensorField2D<T> &field,
    plint boundaryWidth);

template <typename T>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_N<T> *functional, Box2D domain, MultiNTensorField2D<T> &field,
    plint boundaryWidth, plint level = 0);

template <typename T>
void applyProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional2D_N<T> *functional, Box2D domain,
    MultiNTensorField2D<T> &field, MultiNTensorField2D<int> &mask, plint boundaryWidth);

template <typename T>
void integrateProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional2D_N<T> *functional, Box2D domain,
    MultiNTensorField2D<T> &field, MultiNTensorField2D<int> &mask, plint boundaryWidth,
    plint level = 0);

/* ** Typed wrappers with two arguments, bounded and boxed functionals *** */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor1> &lattice1, MultiBlockLattice2D<T2, Descriptor2> &lattice2,
    plint boundaryWidth = Descriptor1<T1>::vicinity);

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor1> &lattice1, MultiBlockLattice2D<T2, Descriptor2> &lattice2,
    plint boundaryWidth = Descriptor1<T1>::vicinity, plint level = 0);

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_SS<T1, T2> *functional, Box2D domain,
    MultiScalarField2D<T1> &field1, MultiScalarField2D<T2> &field2, plint boundaryWidth);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_SS<T1, T2> *functional, Box2D domain,
    MultiScalarField2D<T1> &field1, MultiScalarField2D<T2> &field2, plint boundaryWidth,
    plint level = 0);

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> *functional, Box2D domain,
    MultiTensorField2D<T1, nDim1> &field1, MultiTensorField2D<T2, nDim2> &field2,
    plint boundaryWidth);

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> *functional, Box2D domain,
    MultiTensorField2D<T1, nDim1> &field1, MultiTensorField2D<T2, nDim2> &field2,
    plint boundaryWidth, plint level = 0);

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_NN<T1, T2> *functional, Box2D domain,
    MultiNTensorField2D<T1> &field1, MultiNTensorField2D<T2> &field2, plint boundaryWidth);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_NN<T1, T2> *functional, Box2D domain,
    MultiNTensorField2D<T1> &field1, MultiNTensorField2D<T2> &field2, plint boundaryWidth,
    plint level = 0);

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional2D_NN<T1, T2> *functional, Box2D domain,
    MultiNTensorField2D<T1> &field1, MultiNTensorField2D<T2> &field2,
    MultiNTensorField2D<int> &mask, plint boundaryWidth);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional2D_NN<T1, T2> *functional, Box2D domain,
    MultiNTensorField2D<T1> &field1, MultiNTensorField2D<T2> &field2,
    MultiNTensorField2D<int> &mask, plint boundaryWidth, plint level = 0);

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_ST<T1, T2, nDim> *functional, Box2D domain,
    MultiScalarField2D<T1> &field1, MultiTensorField2D<T2, nDim> &field2, plint boundaryWidth);

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_ST<T1, T2, nDim> *functional, Box2D domain,
    MultiScalarField2D<T1> &field1, MultiTensorField2D<T2, nDim> &field2, plint boundaryWidth,
    plint level = 0);

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_SN<T1, T2> *functional, Box2D domain,
    MultiScalarField2D<T1> &field1, MultiNTensorField2D<T2> &field2, plint boundaryWidth);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_SN<T1, T2> *functional, Box2D domain,
    MultiScalarField2D<T1> &field1, MultiNTensorField2D<T2> &field2, plint boundaryWidth,
    plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_LS<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiScalarField2D<T2> &field,
    plint boundaryWidth = Descriptor<T1>::vicinity);

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_LS<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiScalarField2D<T2> &field,
    plint boundaryWidth = Descriptor<T1>::vicinity, plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiTensorField2D<T2, nDim> &field,
    plint boundaryWidth = Descriptor<T1>::vicinity);

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiTensorField2D<T2, nDim> &field,
    plint boundaryWidth = Descriptor<T1>::vicinity, plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_LN<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field,
    plint boundaryWidth = Descriptor<T1>::vicinity);

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_LN<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field,
    plint boundaryWidth = Descriptor<T1>::vicinity, plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional2D_LN<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field,
    MultiNTensorField2D<int> &mask, plint boundaryWidth = Descriptor<T1>::vicinity);

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional2D_LN<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field,
    MultiNTensorField2D<int> &mask, plint boundaryWidth = Descriptor<T1>::vicinity,
    plint level = 0);

}  // namespace plb

#endif  // MULTI_DATA_PROCESSOR_WRAPPER_2D_H
