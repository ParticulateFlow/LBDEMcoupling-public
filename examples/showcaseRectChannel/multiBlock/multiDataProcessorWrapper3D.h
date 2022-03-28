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

#ifndef MULTI_DATA_PROCESSOR_WRAPPER_3D_H
#define MULTI_DATA_PROCESSOR_WRAPPER_3D_H

#include "atomicBlock/dataProcessingFunctional3D.h"

namespace plb {

class MultiBlock3D;
template <typename T, template <typename U> class Descriptor>
class MultiBlockLattice3D;
template <typename T>
class MultiScalarField3D;
template <typename T, int nDim>
class MultiTensorField3D;
template <typename T>
class MultiNTensorField3D;

/* *************** Generic wrappers, boxed functionals ********************** */

/// This is the most general wrapper for applying a 3D boxed data
/// functional. Use this if none of the more specific wrappers works.
void applyProcessingFunctional(
    BoxProcessingFunctional3D *functional, Box3D domain, std::vector<MultiBlock3D *> multiBlocks);

/// This is the most general wrapper for integrating a 3D boxed data
/// functional. Use this if none of the more specific wrappers works.
void integrateProcessingFunctional(
    BoxProcessingFunctional3D *functional, Box3D domain, std::vector<MultiBlock3D *> multiBlocks,
    plint level = 0);

/// This is the most general wrapper for integrating a 3D boxed data
/// functional. Use this if none of the more specific wrappers works.
void integrateProcessingFunctional(
    BoxProcessingFunctional3D *functional, Box3D domain, MultiBlock3D &actor,
    std::vector<MultiBlock3D *> multiBlockArgs, plint level = 0);

/// Apply a functional on a sequence of block-lattices. If the number
/// of lattices is 1 or 2, you should prefer the _L and _LL version
/// of the functional.
template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    LatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain,
    std::vector<MultiBlockLattice3D<T, Descriptor> *> lattices);

/// Integrate a functional into a sequence of block-lattices. If the
/// number of lattices is 1 or 2, you should prefer the _L and _LL version
/// of the functional.
template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    LatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain,
    std::vector<MultiBlockLattice3D<T, Descriptor> *> lattices, plint level = 0);

/// Apply a functional on a sequence of scalar-fields. If the number
/// of lattices is 1 or 2, you should prefer the _S and _SS version
/// of the functional.
template <typename T>
void applyProcessingFunctional(
    ScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiScalarField3D<T> *> fields);

/// Integrate a functional into a sequence of scalar-fields. If the
/// number of lattices is 1 or 2, you should prefer the _S and _SS version
/// of the functional.
template <typename T>
void integrateProcessingFunctional(
    ScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiScalarField3D<T> *> fields, plint level = 0);

/// Apply a functional on a sequence of tensor-fields. If the number
/// of lattices is 1 or 2, you should prefer the _T and _TT version
/// of the functional.
template <typename T, int nDim>
void applyProcessingFunctional(
    TensorFieldBoxProcessingFunctional3D<T, nDim> *functional, Box3D domain,
    std::vector<MultiTensorField3D<T, nDim> *> fields);

/// Integrate a functional into a sequence of tensor-fields. If the
/// number of lattices is 1 or 2, you should prefer the _T and _TT version
/// of the functional.
template <typename T, int nDim>
void integrateProcessingFunctional(
    TensorFieldBoxProcessingFunctional3D<T, nDim> *functional, Box3D domain,
    std::vector<MultiTensorField3D<T, nDim> *> fields, plint level = 0);

/// Apply a functional on a sequence of n/tensor-fields. If the number
/// of lattices is 1 or 2, you should prefer the _N and _NN version
/// of the functional.
template <typename T>
void applyProcessingFunctional(
    NTensorFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiNTensorField3D<T> *> fields);

/// Apply a functional on a sequence of n/tensor-fields, with bool mask.
/// If the number of lattices is 1 or 2, you should prefer the _N and _NN
/// version of the functional.
template <typename T>
void applyProcessingFunctional(
    MaskedNTensorFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiNTensorField3D<T> *> fields, MultiNTensorField3D<int> &mask);

/// Integrate a functional into a sequence of n/tensor-fields. If the
/// number of lattices is 1 or 2, you should prefer the _N and _NN version
/// of the functional.
template <typename T>
void integrateProcessingFunctional(
    NTensorFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiNTensorField3D<T> *> fields, plint level = 0);

/// Integrate a functional into a sequence of n/tensor-fields, with mask. If the
/// number of lattices is 1 or 2, you should prefer the _N and _NN version
/// of the functional.
template <typename T>
void integrateProcessingFunctional(
    MaskedNTensorFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiNTensorField3D<T> *> fields, MultiNTensorField3D<int> &mask, plint level = 0);

/* *************** Typed wrappers with a single argument, boxed functionals * */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    MultiBlockLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    MultiBlockLattice3D<T, Descriptor> &lattice, plint level = 0);

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    MaskedBoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiNTensorField3D<int> &mask);

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiNTensorField3D<int> &mask, plint level = 0);

template <typename T>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_S<T> *functional, Box3D domain, MultiScalarField3D<T> &field);
template <typename T>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_S<T> *functional, Box3D domain, MultiScalarField3D<T> &field,
    plint level = 0);

template <typename T, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_T<T, nDim> *functional, Box3D domain,
    MultiTensorField3D<T, nDim> &field);
template <typename T, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_T<T, nDim> *functional, Box3D domain,
    MultiTensorField3D<T, nDim> &field, plint level = 0);

template <typename T>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_N<T> *functional, Box3D domain, MultiNTensorField3D<T> &field);
template <typename T>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_N<T> *functional, Box3D domain, MultiNTensorField3D<T> &field,
    plint level = 0);

template <typename T>
void applyProcessingFunctional(
    MaskedBoxProcessingFunctional3D_N<T> *functional, Box3D domain, MultiNTensorField3D<T> &field,
    MultiNTensorField3D<int> &mask);
template <typename T>
void integrateProcessingFunctional(
    MaskedBoxProcessingFunctional3D_N<T> *functional, Box3D domain, MultiNTensorField3D<T> &field,
    MultiNTensorField3D<int> &mask, plint level = 0);

/* *************** Typed wrappers with two arguments ************************ */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor1> &lattice1, MultiBlockLattice3D<T2, Descriptor2> &lattice2);

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor1> &lattice1, MultiBlockLattice3D<T2, Descriptor2> &lattice2,
    plint level = 0);

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain, MultiScalarField3D<T1> &field1,
    MultiScalarField3D<T2> &field2);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain, MultiScalarField3D<T1> &field1,
    MultiScalarField3D<T2> &field2, plint level = 0);

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    MultiTensorField3D<T1, nDim1> &field1, MultiTensorField3D<T2, nDim2> &field2);

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    MultiTensorField3D<T1, nDim1> &field1, MultiTensorField3D<T2, nDim2> &field2, plint level = 0);

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_NN<T1, T2> *functional, Box3D domain, MultiNTensorField3D<T1> &field1,
    MultiNTensorField3D<T2> &field2);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_NN<T1, T2> *functional, Box3D domain, MultiNTensorField3D<T1> &field1,
    MultiNTensorField3D<T2> &field2, plint level = 0);

template <typename T1, typename T2>
void applyProcessingFunctional(
    MaskedBoxProcessingFunctional3D_NN<T1, T2> *functional, Box3D domain,
    MultiNTensorField3D<T1> &field1, MultiNTensorField3D<T2> &field2,
    MultiNTensorField3D<int> &mask);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    MaskedBoxProcessingFunctional3D_NN<T1, T2> *functional, Box3D domain,
    MultiNTensorField3D<T1> &field1, MultiNTensorField3D<T2> &field2,
    MultiNTensorField3D<int> &mask, plint level = 0);

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain,
    MultiScalarField3D<T1> &field1, MultiTensorField3D<T2, nDim> &field2);

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain,
    MultiScalarField3D<T1> &field1, MultiTensorField3D<T2, nDim> &field2, plint level = 0);

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_SN<T1, T2> *functional, Box3D domain, MultiScalarField3D<T1> &field1,
    MultiNTensorField3D<T2> &field2);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_SN<T1, T2> *functional, Box3D domain, MultiScalarField3D<T1> &field1,
    MultiNTensorField3D<T2> &field2, plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiScalarField3D<T2> &field);

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiScalarField3D<T2> &field, plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiTensorField3D<T2, nDim> &field);

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiTensorField3D<T2, nDim> &field,
    plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiNTensorField3D<T2> &field);

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiNTensorField3D<T2> &field, plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    MaskedBoxProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiNTensorField3D<T2> &field,
    MultiNTensorField3D<int> &mask);

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    MaskedBoxProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiNTensorField3D<T2> &field,
    MultiNTensorField3D<int> &mask, plint level = 0);

/* *************** Generic wrappers, dotted functionals ********************* */

void applyProcessingFunctional(
    DotProcessingFunctional3D *functional, DotList3D const &dotList,
    std::vector<MultiBlock3D *> multiBlocks);

void integrateProcessingFunctional(
    DotProcessingFunctional3D *functional, DotList3D const &dotList,
    std::vector<MultiBlock3D *> multiBlocks, plint level = 0);

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    LatticeDotProcessingFunctional3D<T, Descriptor> *functional, DotList3D const &dotList,
    std::vector<MultiBlockLattice3D<T, Descriptor> *> lattices);

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    LatticeDotProcessingFunctional3D<T, Descriptor> *functional, DotList3D const &dotList,
    std::vector<MultiBlockLattice3D<T, Descriptor> *> lattices, plint level = 0);

template <typename T>
void applyProcessingFunctional(
    ScalarFieldDotProcessingFunctional3D<T> *functional, DotList3D const &dotList,
    std::vector<MultiScalarField3D<T> *> fields);

template <typename T>
void integrateProcessingFunctional(
    ScalarFieldDotProcessingFunctional3D<T> *functional, DotList3D const &dotList,
    std::vector<MultiScalarField3D<T> *> fields, plint level = 0);

template <typename T, int nDim>
void applyProcessingFunctional(
    TensorFieldDotProcessingFunctional3D<T, nDim> *functional, DotList3D const &dotList,
    std::vector<MultiTensorField3D<T, nDim> *> fields);

template <typename T, int nDim>
void integrateProcessingFunctional(
    TensorFieldDotProcessingFunctional3D<T, nDim> *functional, DotList3D const &dotList,
    std::vector<MultiTensorField3D<T, nDim> *> fields, plint level = 0);

/* *************** Typed wrappers with a single argument, dotted functionals* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    DotProcessingFunctional3D_L<T, Descriptor> *functional, DotList3D const &dotList,
    MultiBlockLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_L<T, Descriptor> *functional, DotList3D const &dotList,
    MultiBlockLattice3D<T, Descriptor> &lattice, plint level = 0);

template <typename T>
void applyProcessingFunctional(
    DotProcessingFunctional3D_S<T> *functional, DotList3D const &dotList,
    MultiScalarField3D<T> &field);

template <typename T>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_S<T> *functional, DotList3D const &dotList,
    MultiScalarField3D<T> &field, plint level = 0);

template <typename T, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional3D_T<T, nDim> *functional, DotList3D const &dotList,
    MultiTensorField3D<T, nDim> &field);

template <typename T, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_T<T, nDim> *functional, DotList3D const &dotList,
    MultiTensorField3D<T, nDim> &field, plint level = 0);

/* *************** Typed wrappers with two arguments, dotted functionals* *****/

/// Easy instantiation of dotted data processor for lattice-lattice coupling
template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional,
    DotList3D const &dotList, MultiBlockLattice3D<T1, Descriptor1> &lattice1,
    MultiBlockLattice3D<T2, Descriptor2> &lattice2);

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional,
    DotList3D const &dotList, MultiBlockLattice3D<T1, Descriptor1> &lattice1,
    MultiBlockLattice3D<T2, Descriptor2> &lattice2, plint level = 0);

template <typename T1, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_SS<T1, T2> *functional, DotList3D const &dotList,
    MultiScalarField3D<T1> &field1, MultiScalarField3D<T2> &field2);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_SS<T1, T2> *functional, DotList3D const &dotList,
    MultiScalarField3D<T1> &field1, MultiScalarField3D<T2> &field2, plint level = 0);

/// Easy instantiation of dotted data processor for MultiTensorField-MultiTensorField coupling
template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, DotList3D const &dotList,
    MultiTensorField3D<T1, nDim1> &field1, MultiTensorField3D<T2, nDim2> &field2);

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, DotList3D const &dotList,
    MultiTensorField3D<T1, nDim1> &field1, MultiTensorField3D<T2, nDim2> &field2, plint level = 0);

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional3D_ST<T1, T2, nDim> *functional, DotList3D const &dotList,
    MultiScalarField3D<T1> &field1, MultiTensorField3D<T2, nDim> &field2);

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_ST<T1, T2, nDim> *functional, DotList3D const &dotList,
    MultiScalarField3D<T1> &field1, MultiTensorField3D<T2, nDim> &field2, plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, DotList3D const &dotList,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiScalarField3D<T2> &field);

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, DotList3D const &dotList,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiScalarField3D<T2> &field, plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, DotList3D const &dotList,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiTensorField3D<T2, nDim> &field);

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, DotList3D const &dotList,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiTensorField3D<T2, nDim> &field,
    plint level = 0);

/* *************** Generic wrappers, bounded and boxed functionals ********** */

void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D *functional, Box3D domain,
    std::vector<MultiBlock3D *> multiBlocks, plint boundaryWidth);

void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D *functional, Box3D domain,
    std::vector<MultiBlock3D *> multiBlocks, plint boundaryWidth, plint level = 0);

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoundedLatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain,
    std::vector<MultiBlockLattice3D<T, Descriptor> *> lattices,
    plint boundaryWidth = Descriptor<T>::vicinity);

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoundedLatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain,
    std::vector<MultiBlockLattice3D<T, Descriptor> *> lattices,
    plint boundaryWidth = Descriptor<T>::vicinity, plint level = 0);

template <typename T>
void applyProcessingFunctional(
    BoundedScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiScalarField3D<T> *> fields, plint boundaryWidth);

template <typename T>
void integrateProcessingFunctional(
    BoundedScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiScalarField3D<T> *> fields, plint boundaryWidth);

template <typename T, int nDim>
void applyProcessingFunctional(
    BoundedTensorFieldBoxProcessingFunctional3D<T, nDim> *functional, Box3D domain,
    std::vector<MultiTensorField3D<T, nDim> *> fields, plint boundaryWidth);

template <typename T, int nDim>
void integrateProcessingFunctional(
    BoundedTensorFieldBoxProcessingFunctional3D<T, nDim> *functional, Box3D domain,
    std::vector<MultiTensorField3D<T, nDim> *> fields, plint boundaryWidth);

/* ** Typed wrappers with a single argument, bounded and boxed functionals * */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    MultiBlockLattice3D<T, Descriptor> &lattice, plint boundaryWidth = Descriptor<T>::vicinity);

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    MultiBlockLattice3D<T, Descriptor> &lattice, plint boundaryWidth = Descriptor<T>::vicinity,
    plint level = 0);

template <typename T>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_S<T> *functional, Box3D domain, MultiScalarField3D<T> &field,
    plint boundaryWidth);

template <typename T>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_S<T> *functional, Box3D domain, MultiScalarField3D<T> &field,
    plint boundaryWidth, plint level = 0);

template <typename T, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_T<T, nDim> *functional, Box3D domain,
    MultiTensorField3D<T, nDim> &field, plint boundaryWidth);

template <typename T, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_T<T, nDim> *functional, Box3D domain,
    MultiTensorField3D<T, nDim> &field, plint boundaryWidth, plint level = 0);

template <typename T>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_N<T> *functional, Box3D domain, MultiNTensorField3D<T> &field,
    plint boundaryWidth);

template <typename T>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_N<T> *functional, Box3D domain, MultiNTensorField3D<T> &field,
    plint boundaryWidth, plint level = 0);

template <typename T>
void applyProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional3D_N<T> *functional, Box3D domain,
    MultiNTensorField3D<T> &field, MultiNTensorField3D<int> &mask, plint boundaryWidth);

template <typename T>
void integrateProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional3D_N<T> *functional, Box3D domain,
    MultiNTensorField3D<T> &field, MultiNTensorField3D<int> &mask, plint boundaryWidth,
    plint level = 0);

/* ** Typed wrappers with two arguments, bounded and boxed functionals *** */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor1> &lattice1, MultiBlockLattice3D<T2, Descriptor2> &lattice2,
    plint boundaryWidth = Descriptor1<T1>::vicinity);

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor1> &lattice1, MultiBlockLattice3D<T2, Descriptor2> &lattice2,
    plint boundaryWidth = Descriptor1<T1>::vicinity, plint level = 0);

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain,
    MultiScalarField3D<T1> &field1, MultiScalarField3D<T2> &field2, plint boundaryWidth);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain,
    MultiScalarField3D<T1> &field1, MultiScalarField3D<T2> &field2, plint boundaryWidth,
    plint level = 0);

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    MultiTensorField3D<T1, nDim1> &field1, MultiTensorField3D<T2, nDim2> &field2,
    plint boundaryWidth);

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    MultiTensorField3D<T1, nDim1> &field1, MultiTensorField3D<T2, nDim2> &field2,
    plint boundaryWidth, plint level = 0);

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_NN<T1, T2> *functional, Box3D domain,
    MultiNTensorField3D<T1> &field1, MultiNTensorField3D<T2> &field2, plint boundaryWidth);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_NN<T1, T2> *functional, Box3D domain,
    MultiNTensorField3D<T1> &field1, MultiNTensorField3D<T2> &field2, plint boundaryWidth,
    plint level = 0);

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional3D_NN<T1, T2> *functional, Box3D domain,
    MultiNTensorField3D<T1> &field1, MultiNTensorField3D<T2> &field2,
    MultiNTensorField3D<int> &mask, plint boundaryWidth);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional3D_NN<T1, T2> *functional, Box3D domain,
    MultiNTensorField3D<T1> &field1, MultiNTensorField3D<T2> &field2,
    MultiNTensorField3D<int> &mask, plint boundaryWidth, plint level = 0);

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain,
    MultiScalarField3D<T1> &field1, MultiTensorField3D<T2, nDim> &field2, plint boundaryWidth);

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain,
    MultiScalarField3D<T1> &field1, MultiTensorField3D<T2, nDim> &field2, plint boundaryWidth,
    plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiScalarField3D<T2> &field,
    plint boundaryWidth = Descriptor<T1>::vicinity);

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiScalarField3D<T2> &field,
    plint boundaryWidth = Descriptor<T1>::vicinity, plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiTensorField3D<T2, nDim> &field,
    plint boundaryWidth = Descriptor<T1>::vicinity);

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiTensorField3D<T2, nDim> &field,
    plint boundaryWidth = Descriptor<T1>::vicinity, plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiNTensorField3D<T2> &field,
    plint boundaryWidth = Descriptor<T1>::vicinity);

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiNTensorField3D<T2> &field,
    plint boundaryWidth = Descriptor<T1>::vicinity, plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiNTensorField3D<T2> &field,
    MultiNTensorField3D<int> &mask, plint boundaryWidth = Descriptor<T1>::vicinity);

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiNTensorField3D<T2> &field,
    MultiNTensorField3D<int> &mask, plint boundaryWidth = Descriptor<T1>::vicinity,
    plint level = 0);

}  // namespace plb

#endif  // MULTI_DATA_PROCESSOR_WRAPPER_3D_H
