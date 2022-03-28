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
 *
 *
 */

#ifndef MULTI_GRID_DATA_PROCESSOR_WRAPPER_3D_H
#define MULTI_GRID_DATA_PROCESSOR_WRAPPER_3D_H

#include "atomicBlock/dataProcessingFunctional3D.h"

namespace plb {

class MultiGrid3D;
template <typename T, template <typename U> class Descriptor>
class MultiGridLattice3D;
template <typename T>
class MultiGridScalarField3D;
template <typename T, int nDim>
class MultiGridTensorField3D;

/* *************** Generic wrappers, boxed functionals ********************** */

/// This is the most general wrapper for applying a 3D boxed data
/// functional. Use this if none of the more specific wrappers works.
void applyProcessingFunctional(
    BoxProcessingFunctional3D *functional, Box3D domain, std::vector<MultiGrid3D *> multiBlocks,
    plint referenceLevel);

/// This is the most general wrapper for integrating a 3D boxed data
/// functional. Use this if none of the more specific wrappers works.
void integrateProcessingFunctional(
    BoxProcessingFunctional3D *functional, Box3D domain, std::vector<MultiGrid3D *> multiBlocks,
    plint referenceLevel, plint level = 0);

/// Apply a functional on a sequence of block-lattices. If the number
/// of lattices is 1 or 2, you should prefer the _L and _LL version
/// of the functional.
template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    LatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain,
    std::vector<MultiGridLattice3D<T, Descriptor> *> lattices, plint referenceLevel);

/// Integrate a functional into a sequence of block-lattices. If the
/// number of lattices is 1 or 2, you should prefer the _L and _LL version
/// of the functional.
template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    LatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain,
    std::vector<MultiGridLattice3D<T, Descriptor> *> lattices, plint referenceLevel,
    plint level = 0);

/// Apply a functional on a sequence of scalar-fields. If the number
/// of lattices is 1 or 2, you should prefer the _S and _SS version
/// of the functional.
template <typename T>
void applyProcessingFunctional(
    ScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiGridScalarField3D<T> *> fields, plint referenceLevel);

/// Integrate a functional into a sequence of scalar-fields. If the
/// number of lattices is 1 or 2, you should prefer the _S and _SS version
/// of the functional.
template <typename T>
void integrateProcessingFunctional(
    ScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiGridScalarField3D<T> *> fields, plint referenceLevel, plint level = 0);

/// Apply a functional on a sequence of tensor-fields. If the number
/// of lattices is 1 or 2, you should prefer the _T and _TT version
/// of the functional.
template <typename T, int nDim>
void applyProcessingFunctional(
    TensorFieldBoxProcessingFunctional3D<T, nDim> *functional, Box3D domain,
    std::vector<MultiGridTensorField3D<T, nDim> *> fields, plint referenceLevel);

/// Integrate a functional into a sequence of tensor-fields. If the
/// number of lattices is 1 or 2, you should prefer the _T and _TT version
/// of the functional.
template <typename T, int nDim>
void integrateProcessingFunctional(
    TensorFieldBoxProcessingFunctional3D<T, nDim> *functional, Box3D domain,
    std::vector<MultiGridTensorField3D<T, nDim> *> fields, plint referenceLevel, plint level = 0);

/* *************** Typed wrappers with a single argument, boxed functionals * */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    MultiGridLattice3D<T, Descriptor> &lattice, plint referenceLevel);

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    MultiGridLattice3D<T, Descriptor> &lattice, plint referenceLevel, plint level = 0);

template <typename T>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_S<T> *functional, Box3D domain, MultiGridScalarField3D<T> &field,
    plint referenceLevel);
template <typename T>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_S<T> *functional, Box3D domain, MultiGridScalarField3D<T> &field,
    plint referenceLevel, plint level = 0);

template <typename T, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_T<T, nDim> *functional, Box3D domain,
    MultiGridTensorField3D<T, nDim> &field, plint referenceLevel);
template <typename T, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_T<T, nDim> *functional, Box3D domain,
    MultiGridTensorField3D<T, nDim> &field, plint referenceLevel, plint level = 0);

/* *************** Typed wrappers with two arguments ************************ */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor1> &lattice1, MultiGridLattice3D<T2, Descriptor2> &lattice2,
    plint referenceLevel);

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor1> &lattice1, MultiGridLattice3D<T2, Descriptor2> &lattice2,
    plint referenceLevel, plint level = 0);

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain,
    MultiGridScalarField3D<T1> &field1, MultiGridScalarField3D<T2> &field2, plint referenceLevel);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain,
    MultiGridScalarField3D<T1> &field1, MultiGridScalarField3D<T2> &field2, plint referenceLevel,
    plint level = 0);

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    MultiGridTensorField3D<T1, nDim1> &field1, MultiGridTensorField3D<T2, nDim2> &field2,
    plint referenceLevel);

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    MultiGridTensorField3D<T1, nDim1> &field1, MultiGridTensorField3D<T2, nDim2> &field2,
    plint referenceLevel, plint level = 0);

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain,
    MultiGridScalarField3D<T1> &field1, MultiGridTensorField3D<T2, nDim> &field2,
    plint referenceLevel);

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain,
    MultiGridScalarField3D<T1> &field1, MultiGridTensorField3D<T2, nDim> &field2,
    plint referenceLevel, plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridScalarField3D<T2> &field,
    plint referenceLevel);

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridScalarField3D<T2> &field,
    plint referenceLevel, plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridTensorField3D<T2, nDim> &field,
    plint referenceLevel);

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridTensorField3D<T2, nDim> &field,
    plint referenceLevel, plint level = 0);

/* *************** Generic wrappers, dotted functionals ********************* */

void applyProcessingFunctional(
    DotProcessingFunctional3D *functional, DotList3D const &dotList,
    std::vector<MultiGrid3D *> multiBlocks, plint referenceLevel);

void integrateProcessingFunctional(
    DotProcessingFunctional3D *functional, DotList3D const &dotList,
    std::vector<MultiGrid3D *> multiBlocks, plint referenceLevel, plint level = 0);

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    LatticeDotProcessingFunctional3D<T, Descriptor> *functional, DotList3D const &dotList,
    std::vector<MultiGridLattice3D<T, Descriptor> *> lattices, plint referenceLevel);

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    LatticeDotProcessingFunctional3D<T, Descriptor> *functional, DotList3D const &dotList,
    std::vector<MultiGridLattice3D<T, Descriptor> *> lattices, plint referenceLevel,
    plint level = 0);

template <typename T>
void applyProcessingFunctional(
    ScalarFieldDotProcessingFunctional3D<T> *functional, DotList3D const &dotList,
    std::vector<MultiGridScalarField3D<T> *> fields, plint referenceLevel);

template <typename T>
void integrateProcessingFunctional(
    ScalarFieldDotProcessingFunctional3D<T> *functional, DotList3D const &dotList,
    std::vector<MultiGridScalarField3D<T> *> fields, plint referenceLevel, plint level = 0);

template <typename T, int nDim>
void applyProcessingFunctional(
    TensorFieldDotProcessingFunctional3D<T, nDim> *functional, DotList3D const &dotList,
    std::vector<MultiGridTensorField3D<T, nDim> *> fields, plint referenceLevel);

template <typename T, int nDim>
void integrateProcessingFunctional(
    TensorFieldDotProcessingFunctional3D<T, nDim> *functional, DotList3D const &dotList,
    std::vector<MultiGridTensorField3D<T, nDim> *> fields, plint referenceLevel, plint level = 0);

/* *************** Typed wrappers with a single argument, dotted functionals* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    DotProcessingFunctional3D_L<T, Descriptor> *functional, DotList3D const &dotList,
    MultiGridLattice3D<T, Descriptor> &lattice, plint referenceLevel);

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_L<T, Descriptor> *functional, DotList3D const &dotList,
    MultiGridLattice3D<T, Descriptor> &lattice, plint referenceLevel, plint level = 0);

template <typename T>
void applyProcessingFunctional(
    DotProcessingFunctional3D_S<T> *functional, DotList3D const &dotList,
    MultiGridScalarField3D<T> &field, plint referenceLevel);

template <typename T>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_S<T> *functional, DotList3D const &dotList,
    MultiGridScalarField3D<T> &field, plint referenceLevel, plint level = 0);

template <typename T, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional3D_T<T, nDim> *functional, DotList3D const &dotList,
    MultiGridTensorField3D<T, nDim> &field, plint referenceLevel);

template <typename T, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_T<T, nDim> *functional, DotList3D const &dotList,
    MultiGridTensorField3D<T, nDim> &field, plint referenceLevel, plint level = 0);

/* *************** Typed wrappers with two arguments, dotted functionals* *****/

/// Easy instantiation of dotted data processor for lattice-lattice coupling
template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional,
    DotList3D const &dotList, MultiGridLattice3D<T1, Descriptor1> &lattice1,
    MultiGridLattice3D<T2, Descriptor2> &lattice2, plint referenceLevel);

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional,
    DotList3D const &dotList, MultiGridLattice3D<T1, Descriptor1> &lattice1,
    MultiGridLattice3D<T2, Descriptor2> &lattice2, plint referenceLevel, plint level = 0);

template <typename T1, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_SS<T1, T2> *functional, DotList3D const &dotList,
    MultiGridScalarField3D<T1> &field1, MultiGridScalarField3D<T2> &field2, plint referenceLevel);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_SS<T1, T2> *functional, DotList3D const &dotList,
    MultiGridScalarField3D<T1> &field1, MultiGridScalarField3D<T2> &field2, plint referenceLevel,
    plint level = 0);

/// Easy instantiation of dotted data processor for MultiTensorField-MultiTensorField coupling
template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, DotList3D const &dotList,
    MultiGridTensorField3D<T1, nDim1> &field1, MultiGridTensorField3D<T2, nDim2> &field2,
    plint referenceLevel);

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, DotList3D const &dotList,
    MultiGridTensorField3D<T1, nDim1> &field1, MultiGridTensorField3D<T2, nDim2> &field2,
    plint referenceLevel, plint level = 0);

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, DotList3D const &dotList,
    MultiGridTensorField3D<T1, nDim1> &field1, MultiGridTensorField3D<T2, nDim2> &field2,
    plint referenceLevel);

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional3D_ST<T1, T2, nDim> *functional, DotList3D const &dotList,
    MultiGridScalarField3D<T1> &field1, MultiGridTensorField3D<T2, nDim> &field2,
    plint referenceLevel);

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_ST<T1, T2, nDim> *functional, DotList3D const &dotList,
    MultiGridScalarField3D<T1> &field1, MultiGridTensorField3D<T2, nDim> &field2,
    plint referenceLevel, plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, DotList3D const &dotList,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridScalarField3D<T2> &field,
    plint referenceLevel);

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, DotList3D const &dotList,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridScalarField3D<T2> &field,
    plint referenceLevel, plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, DotList3D const &dotList,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridTensorField3D<T2, nDim> &field,
    plint referenceLevel);

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, DotList3D const &dotList,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridTensorField3D<T2, nDim> &field,
    plint referenceLevel, plint level = 0);

/* *************** Generic wrappers, bounded and boxed functionals ********** */

void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D *functional, Box3D domain,
    std::vector<MultiGrid3D *> multiBlocks, plint referenceLevel, plint boundaryWidth);

void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D *functional, Box3D domain,
    std::vector<MultiGrid3D *> multiBlocks, plint boundaryWidth, plint referenceLevel,
    plint level = 0);

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoundedLatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain,
    std::vector<MultiGridLattice3D<T, Descriptor> *> lattices, plint referenceLevel,
    plint boundaryWidth = Descriptor<T>::boundaryWidth);

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoundedLatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain,
    std::vector<MultiGridLattice3D<T, Descriptor> *> lattices, plint referenceLevel,
    plint boundaryWidth = Descriptor<T>::boundaryWidth, plint level = 0);

template <typename T>
void applyProcessingFunctional(
    BoundedScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiGridScalarField3D<T> *> fields, plint referenceLevel, plint boundaryWidth);

template <typename T>
void integrateProcessingFunctional(
    BoundedScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiGridScalarField3D<T> *> fields, plint referenceLevel, plint boundaryWidth);

template <typename T, int nDim>
void applyProcessingFunctional(
    BoundedTensorFieldBoxProcessingFunctional3D<T, nDim> *functional, Box3D domain,
    std::vector<MultiGridTensorField3D<T, nDim> *> fields, plint referenceLevel,
    plint boundaryWidth);

template <typename T, int nDim>
void integrateProcessingFunctional(
    BoundedTensorFieldBoxProcessingFunctional3D<T, nDim> *functional, Box3D domain,
    std::vector<MultiGridTensorField3D<T, nDim> *> fields, plint referenceLevel,
    plint boundaryWidth);

/* ** Typed wrappers with a single argument, bounded and boxed functionals * */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    MultiGridLattice3D<T, Descriptor> &lattice, plint referenceLevel,
    plint boundaryWidth = Descriptor<T>::boundaryWidth);

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    MultiGridLattice3D<T, Descriptor> &lattice, plint referenceLevel,
    plint boundaryWidth = Descriptor<T>::boundaryWidth, plint level = 0);

template <typename T>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_S<T> *functional, Box3D domain,
    MultiGridScalarField3D<T> &field, plint referenceLevel, plint boundaryWidth);

template <typename T>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_S<T> *functional, Box3D domain,
    MultiGridScalarField3D<T> &field, plint referenceLevel, plint boundaryWidth, plint level = 0);

template <typename T, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_T<T, nDim> *functional, Box3D domain,
    MultiGridTensorField3D<T, nDim> &field, plint referenceLevel, plint boundaryWidth);

template <typename T, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_T<T, nDim> *functional, Box3D domain,
    MultiGridTensorField3D<T, nDim> &field, plint referenceLevel, plint boundaryWidth,
    plint level = 0);

/* ** Typed wrappers with two arguments, bounded and boxed functionals *** */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor1> &lattice1, MultiGridLattice3D<T2, Descriptor2> &lattice2,
    plint referenceLevel, plint boundaryWidth = Descriptor1<T1>::boundaryWidth);

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor1> &lattice1, MultiGridLattice3D<T2, Descriptor2> &lattice2,
    plint referenceLevel, plint boundaryWidth = Descriptor1<T1>::boundaryWidth, plint level = 0);

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain,
    MultiGridScalarField3D<T1> &field1, MultiGridScalarField3D<T2> &field2, plint referenceLevel,
    plint boundaryWidth);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain,
    MultiGridScalarField3D<T1> &field1, MultiGridScalarField3D<T2> &field2, plint referenceLevel,
    plint boundaryWidth, plint level = 0);

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    MultiGridTensorField3D<T1, nDim1> &field1, MultiGridTensorField3D<T2, nDim2> &field2,
    plint referenceLevel, plint boundaryWidth);

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    MultiGridTensorField3D<T1, nDim1> &field1, MultiGridTensorField3D<T2, nDim2> &field2,
    plint referenceLevel, plint boundaryWidth, plint level = 0);

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain,
    MultiGridScalarField3D<T1> &field1, MultiGridTensorField3D<T2, nDim> &field2,
    plint referenceLevel, plint boundaryWidth);

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain,
    MultiGridScalarField3D<T1> &field1, MultiGridTensorField3D<T2, nDim> &field2,
    plint referenceLevel, plint boundaryWidth, plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridScalarField3D<T2> &field,
    plint referenceLevel, plint boundaryWidth = Descriptor<T1>::boundaryWidth);

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridScalarField3D<T2> &field,
    plint referenceLevel, plint boundaryWidth = Descriptor<T1>::boundaryWidth, plint level = 0);

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridTensorField3D<T2, nDim> &field,
    plint referenceLevel, plint boundaryWidth = Descriptor<T1>::boundaryWidth);

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridTensorField3D<T2, nDim> &field,
    plint referenceLevel, plint boundaryWidth = Descriptor<T1>::boundaryWidth, plint level = 0);

}  // namespace plb

#endif  // MULTI_GRID_DATA_PROCESSOR_WRAPPER_3D_H
