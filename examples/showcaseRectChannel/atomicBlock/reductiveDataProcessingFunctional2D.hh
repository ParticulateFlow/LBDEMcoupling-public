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

#ifndef REDUCTIVE_DATA_PROCESSING_FUNCTIONAL_2D_HH
#define REDUCTIVE_DATA_PROCESSING_FUNCTIONAL_2D_HH

#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include "atomicBlock/dataProcessor2D.h"
#include "atomicBlock/reductiveDataProcessorWrapper2D.h"
#include "core/plbDebug.h"

namespace plb {

/* *************** ReductiveBoxProcessing2D_L ******************************************* */

template <typename T, template <typename U> class Descriptor>
void ReductiveBoxProcessingFunctional2D_L<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(domain, dynamic_cast<BlockLattice2D<T, Descriptor> &>(*atomicBlocks[0]));
}

/* *************** ReductiveBoxProcessing2D_S ******************************************* */

template <typename T>
void ReductiveBoxProcessingFunctional2D_S<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(domain, dynamic_cast<ScalarField2D<T> &>(*atomicBlocks[0]));
}

/* *************** ReductiveBoxProcessing2D_T ******************************************* */

template <typename T, int nDim>
void ReductiveBoxProcessingFunctional2D_T<T, nDim>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(domain, dynamic_cast<TensorField2D<T, nDim> &>(*atomicBlocks[0]));
}

/* *************** ReductiveBoxProcessing2D_N ******************************************* */

template <typename T>
void ReductiveBoxProcessingFunctional2D_N<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(domain, dynamic_cast<NTensorField2D<T> &>(*atomicBlocks[0]));
}

/* *************** MaskedReductiveBoxProcessing2D_N ******************************************* */

template <typename T>
void MaskedReductiveBoxProcessingFunctional2D_N<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(
        domain, dynamic_cast<NTensorField2D<T> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<int> &>(*atomicBlocks[1]));
}

/* *************** ReductiveBoxProcessing2D_LL******************************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void ReductiveBoxProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    PLB_PRECONDITION(atomicBlocks.size() == 2);
    process(
        domain, dynamic_cast<BlockLattice2D<T1, Descriptor1> &>(*atomicBlocks[0]),
        dynamic_cast<BlockLattice2D<T2, Descriptor2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveBoxProcessing2D_SS ****************************************** */

template <typename T1, typename T2>
void ReductiveBoxProcessingFunctional2D_SS<T1, T2>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(
        domain, dynamic_cast<ScalarField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<ScalarField2D<T2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveBoxProcessing2D_TT ****************************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void ReductiveBoxProcessingFunctional2D_TT<T1, nDim1, T2, nDim2>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(
        domain, dynamic_cast<TensorField2D<T1, nDim1> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField2D<T2, nDim2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveBoxProcessing2D_NN ****************************************** */

template <typename T1, typename T2>
void ReductiveBoxProcessingFunctional2D_NN<T1, T2>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(
        domain, dynamic_cast<NTensorField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]));
}

/* *************** MaskedReductiveBoxProcessing2D_NN ****************************************** */

template <typename T1, typename T2>
void MaskedReductiveBoxProcessingFunctional2D_NN<T1, T2>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(
        domain, dynamic_cast<NTensorField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]),
        dynamic_cast<NTensorField2D<int> &>(*atomicBlocks[2]));
}

/* *************** ReductiveBoxProcessing2D_ST ****************************************** */

template <typename T1, typename T2, int nDim>
void ReductiveBoxProcessingFunctional2D_ST<T1, T2, nDim>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(
        domain, dynamic_cast<ScalarField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField2D<T2, nDim> &>(*atomicBlocks[1]));
}

/* *************** ReductiveBoxProcessing2D_SN ****************************************** */

template <typename T1, typename T2>
void ReductiveBoxProcessingFunctional2D_SN<T1, T2>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(
        domain, dynamic_cast<ScalarField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveBoxProcessing2D_LS ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void ReductiveBoxProcessingFunctional2D_LS<T1, Descriptor, T2>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(
        domain, dynamic_cast<BlockLattice2D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<ScalarField2D<T2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveBoxProcessing2D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void ReductiveBoxProcessingFunctional2D_LT<T1, Descriptor, T2, nDim>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(
        domain, dynamic_cast<BlockLattice2D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField2D<T2, nDim> &>(*atomicBlocks[1]));
}

/* *************** ReductiveBoxProcessing2D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void ReductiveBoxProcessingFunctional2D_LN<T1, Descriptor, T2>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(
        domain, dynamic_cast<BlockLattice2D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]));
}

/* *************** MaskedReductiveBoxProcessing2D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void MaskedReductiveBoxProcessingFunctional2D_LN<T1, Descriptor, T2>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(
        domain, dynamic_cast<BlockLattice2D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]),
        dynamic_cast<NTensorField2D<int> &>(*atomicBlocks[2]));
}

/* *************** ReductiveLatticeBoxProcessing2D ************************************** */

template <typename T, template <typename U> class Descriptor>
void ReductiveLatticeBoxProcessingFunctional2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    std::vector<BlockLattice2D<T, Descriptor> *> lattices(atomicBlocks.size());
    for (pluint iLattice = 0; iLattice < atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice2D<T, Descriptor> *>(atomicBlocks[iLattice]);
    }
    process(domain, lattices);
}

/* *************** ReductiveScalarFieldBoxProcessing2D *********************************** */

template <typename T>
void ReductiveScalarFieldBoxProcessingFunctional2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    std::vector<ScalarField2D<T> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField2D<T> *>(atomicBlocks[iField]);
    }
    process(domain, fields);
}

/* *************** ReductiveTensorFieldBoxProcessing2D *********************************** */

template <typename T, int nDim>
void ReductiveTensorFieldBoxProcessingFunctional2D<T, nDim>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    std::vector<TensorField2D<T, nDim> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField2D<T, nDim> *>(atomicBlocks[iField]);
    }
    process(domain, fields);
}

/* *************** ReductiveNTensorFieldBoxProcessing2D *********************************** */

template <typename T>
void ReductiveNTensorFieldBoxProcessingFunctional2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    std::vector<NTensorField2D<T> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<NTensorField2D<T> *>(atomicBlocks[iField]);
    }
    process(domain, fields);
}

/* *************** ReductiveDotProcessing2D_L ******************************************* */

template <typename T, template <typename U> class Descriptor>
void ReductiveDotProcessingFunctional2D_L<T, Descriptor>::processGenericBlocks(
    DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(dotList, dynamic_cast<BlockLattice2D<T, Descriptor> &>(*atomicBlocks[0]));
}

/* *************** ReductiveDotProcessing2D_S ******************************************* */

template <typename T>
void ReductiveDotProcessingFunctional2D_S<T>::processGenericBlocks(
    DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(dotList, dynamic_cast<ScalarField2D<T> &>(*atomicBlocks[0]));
}

/* *************** ReductiveDotProcessing2D_T ******************************************* */

template <typename T, int nDim>
void ReductiveDotProcessingFunctional2D_T<T, nDim>::processGenericBlocks(
    DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(dotList, dynamic_cast<TensorField2D<T, nDim> &>(*atomicBlocks[0]));
}

/* *************** ReductiveDotProcessing2D_N ******************************************* */

template <typename T>
void ReductiveDotProcessingFunctional2D_N<T>::processGenericBlocks(
    DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(dotList, dynamic_cast<NTensorField2D<T> &>(*atomicBlocks[0]));
}

/* *************** ReductiveDotProcessing2D_LL******************************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void ReductiveDotProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2>::processGenericBlocks(
    DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks)
{
    PLB_PRECONDITION(atomicBlocks.size() == 2);
    process(
        dotList, dynamic_cast<BlockLattice2D<T1, Descriptor1> &>(*atomicBlocks[0]),
        dynamic_cast<BlockLattice2D<T2, Descriptor2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveDotProcessing2D_SS ****************************************** */

template <typename T1, typename T2>
void ReductiveDotProcessingFunctional2D_SS<T1, T2>::processGenericBlocks(
    DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(
        dotList, dynamic_cast<ScalarField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<ScalarField2D<T2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveDotProcessing2D_TT ****************************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void ReductiveDotProcessingFunctional2D_TT<T1, nDim1, T2, nDim2>::processGenericBlocks(
    DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(
        dotList, dynamic_cast<TensorField2D<T1, nDim1> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField2D<T2, nDim2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveDotProcessing2D_NN ****************************************** */

template <typename T1, typename T2>
void ReductiveDotProcessingFunctional2D_NN<T1, T2>::processGenericBlocks(
    DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(
        dotList, dynamic_cast<NTensorField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveDotProcessing2D_ST ****************************************** */

template <typename T1, typename T2, int nDim>
void ReductiveDotProcessingFunctional2D_ST<T1, T2, nDim>::processGenericBlocks(
    DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(
        dotList, dynamic_cast<ScalarField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField2D<T2, nDim> &>(*atomicBlocks[1]));
}

/* *************** ReductiveDotProcessing2D_SN ****************************************** */

template <typename T1, typename T2>
void ReductiveDotProcessingFunctional2D_SN<T1, T2>::processGenericBlocks(
    DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(
        dotList, dynamic_cast<ScalarField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveDotProcessing2D_LS ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void ReductiveDotProcessingFunctional2D_LS<T1, Descriptor, T2>::processGenericBlocks(
    DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(
        dotList, dynamic_cast<BlockLattice2D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<ScalarField2D<T2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveDotProcessing2D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void ReductiveDotProcessingFunctional2D_LT<T1, Descriptor, T2, nDim>::processGenericBlocks(
    DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(
        dotList, dynamic_cast<BlockLattice2D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField2D<T2, nDim> &>(*atomicBlocks[1]));
}

/* *************** ReductiveDotProcessing2D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void ReductiveDotProcessingFunctional2D_LN<T1, Descriptor, T2>::processGenericBlocks(
    DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks)
{
    process(
        dotList, dynamic_cast<BlockLattice2D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveLatticeDotProcessing2D ******************************************* */

template <typename T, template <typename U> class Descriptor>
void ReductiveLatticeDotProcessingFunctional2D<T, Descriptor>::processGenericBlocks(
    DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks)
{
    std::vector<BlockLattice2D<T, Descriptor> *> lattices(atomicBlocks.size());
    for (pluint iLattice = 0; iLattice < atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice2D<T, Descriptor> *>(atomicBlocks[iLattice]);
    }
    process(dotList, lattices);
}

/* *************** ReductiveReductiveScalarFieldDotProcessing2D ********************************* */

template <typename T>
void ReductiveScalarFieldDotProcessingFunctional2D<T>::processGenericBlocks(
    DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks)
{
    std::vector<ScalarField2D<T> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField2D<T> *>(atomicBlocks[iField]);
    }
    process(dotList, fields);
}

/* *************** ReductiveTensorFieldDotProcessing2D ********************************* */

template <typename T, int nDim>
void ReductiveTensorFieldDotProcessingFunctional2D<T, nDim>::processGenericBlocks(
    DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks)
{
    std::vector<TensorField2D<T, nDim> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField2D<T, nDim> *>(atomicBlocks[iField]);
    }
    process(dotList, fields);
}

/* *************** ReductiveNTensorFieldDotProcessing2D ********************************* */

template <typename T>
void ReductiveNTensorFieldDotProcessingFunctional2D<T>::processGenericBlocks(
    DotList2D const &dotList, std::vector<AtomicBlock2D *> atomicBlocks)
{
    std::vector<NTensorField2D<T> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<NTensorField2D<T> *>(atomicBlocks[iField]);
    }
    process(dotList, fields);
}

/* *************** BoundedReductiveBoxProcessing2D_L ********************************* */

template <typename T, template <typename U> class Descriptor>
void BoundedReductiveBoxProcessingFunctional2D_L<T, Descriptor>::processBulkGeneric(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processBulk(domain, dynamic_cast<BlockLattice2D<T, Descriptor> &>(*atomicBlocks[0]));
}

template <typename T, template <typename U> class Descriptor>
void BoundedReductiveBoxProcessingFunctional2D_L<T, Descriptor>::processEdgeGeneric(
    int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processEdge(
        direction, orientation, domain,
        dynamic_cast<BlockLattice2D<T, Descriptor> &>(*atomicBlocks[0]));
}

template <typename T, template <typename U> class Descriptor>
void BoundedReductiveBoxProcessingFunctional2D_L<T, Descriptor>::processCornerGeneric(
    int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, domain, dynamic_cast<BlockLattice2D<T, Descriptor> &>(*atomicBlocks[0]));
}

/* *************** BoundedReductiveBoxProcessing2D_S ******************************************* */

template <typename T>
void BoundedReductiveBoxProcessingFunctional2D_S<T>::processBulkGeneric(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processBulk(domain, dynamic_cast<ScalarField2D<T> &>(*atomicBlocks[0]));
}

template <typename T>
void BoundedReductiveBoxProcessingFunctional2D_S<T>::processEdgeGeneric(
    int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processEdge(direction, orientation, domain, dynamic_cast<ScalarField2D<T> &>(*atomicBlocks[0]));
}

template <typename T>
void BoundedReductiveBoxProcessingFunctional2D_S<T>::processCornerGeneric(
    int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processCorner(normalX, normalY, domain, dynamic_cast<ScalarField2D<T> &>(*atomicBlocks[0]));
}

/* *************** BoundedReductiveBoxProcessing2D_T ********************************* */

template <typename T, int nDim>
void BoundedReductiveBoxProcessingFunctional2D_T<T, nDim>::processBulkGeneric(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processBulk(domain, dynamic_cast<TensorField2D<T, nDim> &>(*atomicBlocks[0]));
}

template <typename T, int nDim>
void BoundedReductiveBoxProcessingFunctional2D_T<T, nDim>::processEdgeGeneric(
    int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processEdge(
        direction, orientation, domain, dynamic_cast<TensorField2D<T, nDim> &>(*atomicBlocks[0]));
}

template <typename T, int nDim>
void BoundedReductiveBoxProcessingFunctional2D_T<T, nDim>::processCornerGeneric(
    int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, domain, dynamic_cast<TensorField2D<T, nDim> &>(*atomicBlocks[0]));
}

/* *************** BoundedReductiveBoxProcessing2D_N ********************************* */

template <typename T>
void BoundedReductiveBoxProcessingFunctional2D_N<T>::processBulkGeneric(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processBulk(domain, dynamic_cast<NTensorField2D<T> &>(*atomicBlocks[0]));
}

template <typename T>
void BoundedReductiveBoxProcessingFunctional2D_N<T>::processEdgeGeneric(
    int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processEdge(
        direction, orientation, domain, dynamic_cast<NTensorField2D<T> &>(*atomicBlocks[0]));
}

template <typename T>
void BoundedReductiveBoxProcessingFunctional2D_N<T>::processCornerGeneric(
    int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processCorner(normalX, normalY, domain, dynamic_cast<NTensorField2D<T> &>(*atomicBlocks[0]));
}

/* *************** BoundedMaskedReductiveBoxProcessing2D_N ********************************* */

template <typename T>
void BoundedMaskedReductiveBoxProcessingFunctional2D_N<T>::processBulkGeneric(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processBulk(
        domain, dynamic_cast<NTensorField2D<T> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<int> &>(*atomicBlocks[1]));
}

template <typename T>
void BoundedMaskedReductiveBoxProcessingFunctional2D_N<T>::processEdgeGeneric(
    int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processEdge(
        direction, orientation, domain, dynamic_cast<NTensorField2D<T> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<int> &>(*atomicBlocks[1]));
}

template <typename T>
void BoundedMaskedReductiveBoxProcessingFunctional2D_N<T>::processCornerGeneric(
    int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, domain, dynamic_cast<NTensorField2D<T> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<int> &>(*atomicBlocks[1]));
}

/* *************** BoundedReductiveBoxProcessing2D_LL ****************************************** */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void BoundedReductiveBoxProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2>::
    processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    PLB_PRECONDITION(atomicBlocks.size() == 2);
    processBulk(
        domain, dynamic_cast<BlockLattice2D<T1, Descriptor1> &>(*atomicBlocks[0]),
        dynamic_cast<BlockLattice2D<T2, Descriptor2> &>(*atomicBlocks[1]));
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void BoundedReductiveBoxProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2>::
    processEdgeGeneric(
        int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    PLB_PRECONDITION(atomicBlocks.size() == 2);
    processEdge(
        direction, orientation, domain,
        dynamic_cast<BlockLattice2D<T1, Descriptor1> &>(*atomicBlocks[0]),
        dynamic_cast<BlockLattice2D<T2, Descriptor2> &>(*atomicBlocks[1]));
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void BoundedReductiveBoxProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2>::
    processCornerGeneric(
        int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    PLB_PRECONDITION(atomicBlocks.size() == 2);
    processCorner(
        normalX, normalY, domain, dynamic_cast<BlockLattice2D<T1, Descriptor1> &>(*atomicBlocks[0]),
        dynamic_cast<BlockLattice2D<T2, Descriptor2> &>(*atomicBlocks[1]));
}

/* *************** BoundedReductiveBoxProcessing2D_SS ****************************************** */

template <typename T1, typename T2>
void BoundedReductiveBoxProcessingFunctional2D_SS<T1, T2>::processBulkGeneric(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processBulk(
        domain, dynamic_cast<ScalarField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<ScalarField2D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, typename T2>
void BoundedReductiveBoxProcessingFunctional2D_SS<T1, T2>::processEdgeGeneric(
    int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processEdge(
        direction, orientation, domain, dynamic_cast<ScalarField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<ScalarField2D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, typename T2>
void BoundedReductiveBoxProcessingFunctional2D_SS<T1, T2>::processCornerGeneric(
    int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, domain, dynamic_cast<ScalarField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<ScalarField2D<T2> &>(*atomicBlocks[1]));
}

/* *************** BoundedReductiveBoxProcessing2D_TT ****************************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void BoundedReductiveBoxProcessingFunctional2D_TT<T1, nDim1, T2, nDim2>::processBulkGeneric(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processBulk(
        domain, dynamic_cast<TensorField2D<T1, nDim1> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField2D<T2, nDim2> &>(*atomicBlocks[1]));
}

template <typename T1, int nDim1, typename T2, int nDim2>
void BoundedReductiveBoxProcessingFunctional2D_TT<T1, nDim1, T2, nDim2>::processEdgeGeneric(
    int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processEdge(
        direction, orientation, domain, dynamic_cast<TensorField2D<T1, nDim1> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField2D<T2, nDim2> &>(*atomicBlocks[1]));
}

template <typename T1, int nDim1, typename T2, int nDim2>
void BoundedReductiveBoxProcessingFunctional2D_TT<T1, nDim1, T2, nDim2>::processCornerGeneric(
    int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, domain, dynamic_cast<TensorField2D<T1, nDim1> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField2D<T2, nDim2> &>(*atomicBlocks[1]));
}

/* *************** BoundedReductiveBoxProcessing2D_NN ****************************************** */

template <typename T1, typename T2>
void BoundedReductiveBoxProcessingFunctional2D_NN<T1, T2>::processBulkGeneric(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processBulk(
        domain, dynamic_cast<NTensorField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, typename T2>
void BoundedReductiveBoxProcessingFunctional2D_NN<T1, T2>::processEdgeGeneric(
    int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processEdge(
        direction, orientation, domain, dynamic_cast<NTensorField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, typename T2>
void BoundedReductiveBoxProcessingFunctional2D_NN<T1, T2>::processCornerGeneric(
    int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, domain, dynamic_cast<NTensorField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]));
}

/* *************** BoundedMaskedReductiveBoxProcessing2D_NN
 * ****************************************** */

template <typename T1, typename T2>
void BoundedMaskedReductiveBoxProcessingFunctional2D_NN<T1, T2>::processBulkGeneric(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processBulk(
        domain, dynamic_cast<NTensorField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]),
        dynamic_cast<NTensorField2D<int> &>(*atomicBlocks[2]));
}

template <typename T1, typename T2>
void BoundedMaskedReductiveBoxProcessingFunctional2D_NN<T1, T2>::processEdgeGeneric(
    int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processEdge(
        direction, orientation, domain, dynamic_cast<NTensorField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]),
        dynamic_cast<NTensorField2D<int> &>(*atomicBlocks[2]));
}

template <typename T1, typename T2>
void BoundedMaskedReductiveBoxProcessingFunctional2D_NN<T1, T2>::processCornerGeneric(
    int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, domain, dynamic_cast<NTensorField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]),
        dynamic_cast<NTensorField2D<int> &>(*atomicBlocks[2]));
}

/* *************** BoundedReductiveBoxProcessing2D_ST ****************************************** */

template <typename T1, typename T2, int nDim>
void BoundedReductiveBoxProcessingFunctional2D_ST<T1, T2, nDim>::processBulkGeneric(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processBulk(
        domain, dynamic_cast<ScalarField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField2D<T2, nDim> &>(*atomicBlocks[1]));
}

template <typename T1, typename T2, int nDim>
void BoundedReductiveBoxProcessingFunctional2D_ST<T1, T2, nDim>::processEdgeGeneric(
    int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processEdge(
        direction, orientation, domain, dynamic_cast<ScalarField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField2D<T2, nDim> &>(*atomicBlocks[1]));
}

template <typename T1, typename T2, int nDim>
void BoundedReductiveBoxProcessingFunctional2D_ST<T1, T2, nDim>::processCornerGeneric(
    int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, domain, dynamic_cast<ScalarField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField2D<T2, nDim> &>(*atomicBlocks[1]));
}

/* *************** BoundedReductiveBoxProcessing2D_SN ****************************************** */

template <typename T1, typename T2>
void BoundedReductiveBoxProcessingFunctional2D_SN<T1, T2>::processBulkGeneric(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processBulk(
        domain, dynamic_cast<ScalarField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, typename T2>
void BoundedReductiveBoxProcessingFunctional2D_SN<T1, T2>::processEdgeGeneric(
    int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processEdge(
        direction, orientation, domain, dynamic_cast<ScalarField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, typename T2>
void BoundedReductiveBoxProcessingFunctional2D_SN<T1, T2>::processCornerGeneric(
    int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, domain, dynamic_cast<ScalarField2D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]));
}

/* *************** BoundedReductiveBoxProcessing2D_LS ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void BoundedReductiveBoxProcessingFunctional2D_LS<T1, Descriptor, T2>::processBulkGeneric(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processBulk(
        domain, dynamic_cast<BlockLattice2D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<ScalarField2D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void BoundedReductiveBoxProcessingFunctional2D_LS<T1, Descriptor, T2>::processEdgeGeneric(
    int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processEdge(
        direction, orientation, domain,
        dynamic_cast<BlockLattice2D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<ScalarField2D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void BoundedReductiveBoxProcessingFunctional2D_LS<T1, Descriptor, T2>::processCornerGeneric(
    int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, domain, dynamic_cast<BlockLattice2D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<ScalarField2D<T2> &>(*atomicBlocks[1]));
}

/* *************** BoundedReductiveBoxProcessing2D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void BoundedReductiveBoxProcessingFunctional2D_LT<T1, Descriptor, T2, nDim>::processBulkGeneric(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processBulk(
        domain, dynamic_cast<BlockLattice2D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField2D<T2, nDim> &>(*atomicBlocks[1]));
}

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void BoundedReductiveBoxProcessingFunctional2D_LT<T1, Descriptor, T2, nDim>::processEdgeGeneric(
    int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processEdge(
        direction, orientation, domain,
        dynamic_cast<BlockLattice2D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField2D<T2, nDim> &>(*atomicBlocks[1]));
}

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void BoundedReductiveBoxProcessingFunctional2D_LT<T1, Descriptor, T2, nDim>::processCornerGeneric(
    int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, domain, dynamic_cast<BlockLattice2D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField2D<T2, nDim> &>(*atomicBlocks[1]));
}

/* *************** BoundedReductiveBoxProcessing2D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void BoundedReductiveBoxProcessingFunctional2D_LN<T1, Descriptor, T2>::processBulkGeneric(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processBulk(
        domain, dynamic_cast<BlockLattice2D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void BoundedReductiveBoxProcessingFunctional2D_LN<T1, Descriptor, T2>::processEdgeGeneric(
    int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processEdge(
        direction, orientation, domain,
        dynamic_cast<BlockLattice2D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void BoundedReductiveBoxProcessingFunctional2D_LN<T1, Descriptor, T2>::processCornerGeneric(
    int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, domain, dynamic_cast<BlockLattice2D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]));
}

/* *************** BoundedMaskedReductiveBoxProcessing2D_LN
 * ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void BoundedMaskedReductiveBoxProcessingFunctional2D_LN<T1, Descriptor, T2>::processBulkGeneric(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processBulk(
        domain, dynamic_cast<BlockLattice2D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]),
        dynamic_cast<NTensorField2D<int> &>(*atomicBlocks[2]));
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void BoundedMaskedReductiveBoxProcessingFunctional2D_LN<T1, Descriptor, T2>::processEdgeGeneric(
    int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processEdge(
        direction, orientation, domain,
        dynamic_cast<BlockLattice2D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]),
        dynamic_cast<NTensorField2D<int> &>(*atomicBlocks[2]));
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void BoundedMaskedReductiveBoxProcessingFunctional2D_LN<T1, Descriptor, T2>::processCornerGeneric(
    int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, domain, dynamic_cast<BlockLattice2D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField2D<T2> &>(*atomicBlocks[1]),
        dynamic_cast<NTensorField2D<int> &>(*atomicBlocks[2]));
}

/* *************** BoundedReductiveLatticeBoxProcessing2D ********************************* */

template <typename T, template <typename U> class Descriptor>
void BoundedReductiveLatticeBoxProcessingFunctional2D<T, Descriptor>::processBulkGeneric(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    std::vector<BlockLattice2D<T, Descriptor> *> lattices(atomicBlocks.size());
    for (pluint iLattice = 0; iLattice < atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice2D<T, Descriptor> *>(atomicBlocks[iLattice]);
    }
    processBulk(domain, lattices);
}

template <typename T, template <typename U> class Descriptor>
void BoundedReductiveLatticeBoxProcessingFunctional2D<T, Descriptor>::processEdgeGeneric(
    int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    std::vector<BlockLattice2D<T, Descriptor> *> lattices(atomicBlocks.size());
    for (pluint iLattice = 0; iLattice < atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice2D<T, Descriptor> *>(atomicBlocks[iLattice]);
    }
    processEdge(direction, orientation, domain, lattices);
}

template <typename T, template <typename U> class Descriptor>
void BoundedReductiveLatticeBoxProcessingFunctional2D<T, Descriptor>::processCornerGeneric(
    int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    std::vector<BlockLattice2D<T, Descriptor> *> lattices(atomicBlocks.size());
    for (pluint iLattice = 0; iLattice < atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice2D<T, Descriptor> *>(atomicBlocks[iLattice]);
    }
    processCorner(normalX, normalY, domain, lattices);
}

/* *************** BoundedReductiveScalarFieldBoxProcessing2D ************************ */

template <typename T>
void BoundedReductiveScalarFieldBoxProcessingFunctional2D<T>::processBulkGeneric(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    std::vector<ScalarField2D<T> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField2D<T> *>(atomicBlocks[iField]);
    }
    processBulk(domain, fields);
}

template <typename T>
void BoundedReductiveScalarFieldBoxProcessingFunctional2D<T>::processEdgeGeneric(
    int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    std::vector<ScalarField2D<T> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField2D<T> *>(atomicBlocks[iField]);
    }
    processEdge(direction, orientation, domain, fields);
}

template <typename T>
void BoundedReductiveScalarFieldBoxProcessingFunctional2D<T>::processCornerGeneric(
    int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    std::vector<ScalarField2D<T> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField2D<T> *>(atomicBlocks[iField]);
    }
    processCorner(normalX, normalY, domain, fields);
}

/* *************** BoundedReductiveTensorFieldBoxProcessing2D ************************ */

template <typename T, int nDim>
void BoundedReductiveTensorFieldBoxProcessingFunctional2D<T, nDim>::processBulkGeneric(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    std::vector<TensorField2D<T, nDim> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField2D<T, nDim> *>(atomicBlocks[iField]);
    }
    processBulk(domain, fields);
}

template <typename T, int nDim>
void BoundedReductiveTensorFieldBoxProcessingFunctional2D<T, nDim>::processEdgeGeneric(
    int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    std::vector<TensorField2D<T, nDim> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField2D<T, nDim> *>(atomicBlocks[iField]);
    }
    processEdge(direction, orientation, domain, fields);
}

template <typename T, int nDim>
void BoundedReductiveTensorFieldBoxProcessingFunctional2D<T, nDim>::processCornerGeneric(
    int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    std::vector<TensorField2D<T, nDim> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField2D<T, nDim> *>(atomicBlocks[iField]);
    }
    processCorner(normalX, normalY, domain, fields);
}

/* *************** BoundedReductiveNTensorFieldBoxProcessing2D ************************ */

template <typename T>
void BoundedReductiveNTensorFieldBoxProcessingFunctional2D<T>::processBulkGeneric(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    std::vector<NTensorField2D<T> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<NTensorField2D<T> *>(atomicBlocks[iField]);
    }
    processBulk(domain, fields);
}

template <typename T>
void BoundedReductiveNTensorFieldBoxProcessingFunctional2D<T>::processEdgeGeneric(
    int direction, int orientation, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    std::vector<NTensorField2D<T> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<NTensorField2D<T> *>(atomicBlocks[iField]);
    }
    processEdge(direction, orientation, domain, fields);
}

template <typename T>
void BoundedReductiveNTensorFieldBoxProcessingFunctional2D<T>::processCornerGeneric(
    int normalX, int normalY, Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    std::vector<NTensorField2D<T> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<NTensorField2D<T> *>(atomicBlocks[iField]);
    }
    processCorner(normalX, normalY, domain, fields);
}

}  // namespace plb

#endif  // REDUCTIVE_DATA_PROCESSING_FUNCTIONAL_2D_HH
