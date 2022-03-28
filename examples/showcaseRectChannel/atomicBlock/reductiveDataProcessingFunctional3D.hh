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

#ifndef REDUCTIVE_DATA_PROCESSING_FUNCTIONAL_3D_HH
#define REDUCTIVE_DATA_PROCESSING_FUNCTIONAL_3D_HH

#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "atomicBlock/reductiveDataProcessorWrapper3D.h"
#include "core/plbDebug.h"

namespace plb {

/* *************** ReductiveBoxProcessing3D_L ******************************************* */

template <typename T, template <typename U> class Descriptor>
void ReductiveBoxProcessingFunctional3D_L<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(domain, dynamic_cast<BlockLattice3D<T, Descriptor> &>(*atomicBlocks[0]));
}

/* *************** ReductiveBoxProcessing3D_S ******************************************* */

template <typename T>
void ReductiveBoxProcessingFunctional3D_S<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(domain, dynamic_cast<ScalarField3D<T> &>(*atomicBlocks[0]));
}

/* *************** ReductiveBoxProcessing3D_T ******************************************* */

template <typename T, int nDim>
void ReductiveBoxProcessingFunctional3D_T<T, nDim>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(domain, dynamic_cast<TensorField3D<T, nDim> &>(*atomicBlocks[0]));
}

/* *************** ReductiveBoxProcessing3D_N ******************************************* */

template <typename T>
void ReductiveBoxProcessingFunctional3D_N<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(domain, dynamic_cast<NTensorField3D<T> &>(*atomicBlocks[0]));
}

/* *************** ReductiveBoxProcessing3D_N ******************************************* */

template <typename T>
void MaskedReductiveBoxProcessingFunctional3D_N<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(
        domain, dynamic_cast<NTensorField3D<T> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<int> &>(*atomicBlocks[1]));
}

/* *************** ReductiveBoxProcessing3D_LL******************************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void ReductiveBoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    PLB_PRECONDITION(atomicBlocks.size() == 2);
    process(
        domain, dynamic_cast<BlockLattice3D<T1, Descriptor1> &>(*atomicBlocks[0]),
        dynamic_cast<BlockLattice3D<T2, Descriptor2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveBoxProcessing3D_SS ****************************************** */

template <typename T1, typename T2>
void ReductiveBoxProcessingFunctional3D_SS<T1, T2>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(
        domain, dynamic_cast<ScalarField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<ScalarField3D<T2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveBoxProcessing3D_TT ****************************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void ReductiveBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(
        domain, dynamic_cast<TensorField3D<T1, nDim1> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField3D<T2, nDim2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveBoxProcessing3D_NN ****************************************** */

template <typename T1, typename T2>
void ReductiveBoxProcessingFunctional3D_NN<T1, T2>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(
        domain, dynamic_cast<NTensorField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveBoxProcessing3D_NN ****************************************** */

template <typename T1, typename T2>
void MaskedReductiveBoxProcessingFunctional3D_NN<T1, T2>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(
        domain, dynamic_cast<NTensorField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]),
        dynamic_cast<NTensorField3D<int> &>(*atomicBlocks[2]));
}

/* *************** ReductiveBoxProcessing3D_ST ****************************************** */

template <typename T1, typename T2, int nDim>
void ReductiveBoxProcessingFunctional3D_ST<T1, T2, nDim>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(
        domain, dynamic_cast<ScalarField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField3D<T2, nDim> &>(*atomicBlocks[1]));
}

/* *************** ReductiveBoxProcessing3D_SN ****************************************** */

template <typename T1, typename T2>
void ReductiveBoxProcessingFunctional3D_SN<T1, T2>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(
        domain, dynamic_cast<ScalarField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveBoxProcessing3D_LS ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void ReductiveBoxProcessingFunctional3D_LS<T1, Descriptor, T2>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(
        domain, dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<ScalarField3D<T2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveBoxProcessing3D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void ReductiveBoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(
        domain, dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField3D<T2, nDim> &>(*atomicBlocks[1]));
}

/* *************** ReductiveBoxProcessing3D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void ReductiveBoxProcessingFunctional3D_LN<T1, Descriptor, T2>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(
        domain, dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]));
}

/* *************** MaskedReductiveBoxProcessing3D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void MaskedReductiveBoxProcessingFunctional3D_LN<T1, Descriptor, T2>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(
        domain, dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]),
        dynamic_cast<NTensorField3D<int> &>(*atomicBlocks[2]));
}

/* *************** ReductiveLatticeBoxProcessing3D ************************************** */

template <typename T, template <typename U> class Descriptor>
void ReductiveLatticeBoxProcessingFunctional3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<BlockLattice3D<T, Descriptor> *> lattices(atomicBlocks.size());
    for (pluint iLattice = 0; iLattice < atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice3D<T, Descriptor> *>(atomicBlocks[iLattice]);
    }
    process(domain, lattices);
}

/* *************** ReductiveScalarFieldBoxProcessing3D *********************************** */

template <typename T>
void ReductiveScalarFieldBoxProcessingFunctional3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<ScalarField3D<T> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[iField]);
    }
    process(domain, fields);
}

/* *************** ReductiveTensorFieldBoxProcessing3D *********************************** */

template <typename T, int nDim>
void ReductiveTensorFieldBoxProcessingFunctional3D<T, nDim>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<TensorField3D<T, nDim> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField3D<T, nDim> *>(atomicBlocks[iField]);
    }
    process(domain, fields);
}

/* *************** ReductiveNTensorFieldBoxProcessing3D *********************************** */

template <typename T>
void ReductiveNTensorFieldBoxProcessingFunctional3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<NTensorField3D<T> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<NTensorField3D<T> *>(atomicBlocks[iField]);
    }
    process(domain, fields);
}

/* *************** ReductiveDotProcessing3D_L ******************************************* */

template <typename T, template <typename U> class Descriptor>
void ReductiveDotProcessingFunctional3D_L<T, Descriptor>::processGenericBlocks(
    DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(dotList, dynamic_cast<BlockLattice3D<T, Descriptor> &>(*atomicBlocks[0]));
}

/* *************** ReductiveDotProcessing3D_S ******************************************* */

template <typename T>
void ReductiveDotProcessingFunctional3D_S<T>::processGenericBlocks(
    DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(dotList, dynamic_cast<ScalarField3D<T> &>(*atomicBlocks[0]));
}

/* *************** ReductiveDotProcessing3D_T ******************************************* */

template <typename T, int nDim>
void ReductiveDotProcessingFunctional3D_T<T, nDim>::processGenericBlocks(
    DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(dotList, dynamic_cast<TensorField3D<T, nDim> &>(*atomicBlocks[0]));
}

/* *************** ReductiveDotProcessing3D_N ******************************************* */

template <typename T>
void ReductiveDotProcessingFunctional3D_N<T>::processGenericBlocks(
    DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(dotList, dynamic_cast<NTensorField3D<T> &>(*atomicBlocks[0]));
}

/* *************** ReductiveDotProcessing3D_LL******************************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void ReductiveDotProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2>::processGenericBlocks(
    DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks)
{
    PLB_PRECONDITION(atomicBlocks.size() == 2);
    process(
        dotList, dynamic_cast<BlockLattice3D<T1, Descriptor1> &>(*atomicBlocks[0]),
        dynamic_cast<BlockLattice3D<T2, Descriptor2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveDotProcessing3D_SS ****************************************** */

template <typename T1, typename T2>
void ReductiveDotProcessingFunctional3D_SS<T1, T2>::processGenericBlocks(
    DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(
        dotList, dynamic_cast<ScalarField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<ScalarField3D<T2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveDotProcessing3D_TT ****************************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void ReductiveDotProcessingFunctional3D_TT<T1, nDim1, T2, nDim2>::processGenericBlocks(
    DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(
        dotList, dynamic_cast<TensorField3D<T1, nDim1> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField3D<T2, nDim2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveDotProcessing3D_NN ****************************************** */

template <typename T1, typename T2>
void ReductiveDotProcessingFunctional3D_NN<T1, T2>::processGenericBlocks(
    DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(
        dotList, dynamic_cast<NTensorField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveDotProcessing3D_ST ****************************************** */

template <typename T1, typename T2, int nDim>
void ReductiveDotProcessingFunctional3D_ST<T1, T2, nDim>::processGenericBlocks(
    DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(
        dotList, dynamic_cast<ScalarField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField3D<T2, nDim> &>(*atomicBlocks[1]));
}

/* *************** ReductiveDotProcessing3D_SN ****************************************** */

template <typename T1, typename T2>
void ReductiveDotProcessingFunctional3D_SN<T1, T2>::processGenericBlocks(
    DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(
        dotList, dynamic_cast<ScalarField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveDotProcessing3D_LS ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void ReductiveDotProcessingFunctional3D_LS<T1, Descriptor, T2>::processGenericBlocks(
    DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(
        dotList, dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<ScalarField3D<T2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveDotProcessing3D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void ReductiveDotProcessingFunctional3D_LT<T1, Descriptor, T2, nDim>::processGenericBlocks(
    DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(
        dotList, dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField3D<T2, nDim> &>(*atomicBlocks[1]));
}

/* *************** ReductiveDotProcessing3D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void ReductiveDotProcessingFunctional3D_LN<T1, Descriptor, T2>::processGenericBlocks(
    DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks)
{
    process(
        dotList, dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]));
}

/* *************** ReductiveLatticeDotProcessing3D ******************************************* */

template <typename T, template <typename U> class Descriptor>
void ReductiveLatticeDotProcessingFunctional3D<T, Descriptor>::processGenericBlocks(
    DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<BlockLattice3D<T, Descriptor> *> lattices(atomicBlocks.size());
    for (pluint iLattice = 0; iLattice < atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice3D<T, Descriptor> *>(atomicBlocks[iLattice]);
    }
    process(dotList, lattices);
}

/* *************** ReductiveReductiveScalarFieldDotProcessing3D ********************************* */

template <typename T>
void ReductiveScalarFieldDotProcessingFunctional3D<T>::processGenericBlocks(
    DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<ScalarField3D<T> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[iField]);
    }
    process(dotList, fields);
}

/* *************** ReductiveTensorFieldDotProcessing3D ********************************* */

template <typename T, int nDim>
void ReductiveTensorFieldDotProcessingFunctional3D<T, nDim>::processGenericBlocks(
    DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<TensorField3D<T, nDim> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField3D<T, nDim> *>(atomicBlocks[iField]);
    }
    process(dotList, fields);
}

/* *************** ReductiveNTensorFieldDotProcessing3D ********************************* */

template <typename T>
void ReductiveNTensorFieldDotProcessingFunctional3D<T>::processGenericBlocks(
    DotList3D const &dotList, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<NTensorField3D<T> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<NTensorField3D<T> *>(atomicBlocks[iField]);
    }
    process(dotList, fields);
}

/* *************** BoundedReductiveBoxProcessing3D_L ******************************************* */

template <typename T, template <typename U> class Descriptor>
void BoundedReductiveBoxProcessingFunctional3D_L<T, Descriptor>::processBulkGeneric(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processBulk(domain, dynamic_cast<BlockLattice3D<T, Descriptor> &>(*atomicBlocks[0]));
}

template <typename T, template <typename U> class Descriptor>
void BoundedReductiveBoxProcessingFunctional3D_L<T, Descriptor>::processPlaneGeneric(
    int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processPlane(
        direction, orientation, domain,
        dynamic_cast<BlockLattice3D<T, Descriptor> &>(*atomicBlocks[0]));
}

template <typename T, template <typename U> class Descriptor>
void BoundedReductiveBoxProcessingFunctional3D_L<T, Descriptor>::processEdgeGeneric(
    int plane, int normal1, int normal2, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processEdge(
        plane, normal1, normal2, domain,
        dynamic_cast<BlockLattice3D<T, Descriptor> &>(*atomicBlocks[0]));
}

template <typename T, template <typename U> class Descriptor>
void BoundedReductiveBoxProcessingFunctional3D_L<T, Descriptor>::processCornerGeneric(
    int normalX, int normalY, int normalZ, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, normalZ, domain,
        dynamic_cast<BlockLattice3D<T, Descriptor> &>(*atomicBlocks[0]));
}

/* *************** BoundedReductiveBoxProcessing3D_S ******************************************* */

template <typename T>
void BoundedReductiveBoxProcessingFunctional3D_S<T>::processBulkGeneric(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processBulk(domain, dynamic_cast<ScalarField3D<T> &>(*atomicBlocks[0]));
}

template <typename T>
void BoundedReductiveBoxProcessingFunctional3D_S<T>::processPlaneGeneric(
    int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processPlane(
        direction, orientation, domain, dynamic_cast<ScalarField3D<T> &>(*atomicBlocks[0]));
}

template <typename T>
void BoundedReductiveBoxProcessingFunctional3D_S<T>::processEdgeGeneric(
    int plane, int normal1, int normal2, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processEdge(
        plane, normal1, normal2, domain, dynamic_cast<ScalarField3D<T> &>(*atomicBlocks[0]));
}

template <typename T>
void BoundedReductiveBoxProcessingFunctional3D_S<T>::processCornerGeneric(
    int normalX, int normalY, int normalZ, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, normalZ, domain, dynamic_cast<ScalarField3D<T> &>(*atomicBlocks[0]));
}

/* *************** BoundedReductiveBoxProcessing3D_T ******************************************* */

template <typename T, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_T<T, nDim>::processBulkGeneric(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processBulk(domain, dynamic_cast<TensorField3D<T, nDim> &>(*atomicBlocks[0]));
}

template <typename T, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_T<T, nDim>::processPlaneGeneric(
    int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processPlane(
        direction, orientation, domain, dynamic_cast<TensorField3D<T, nDim> &>(*atomicBlocks[0]));
}

template <typename T, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_T<T, nDim>::processEdgeGeneric(
    int plane, int normal1, int normal2, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processEdge(
        plane, normal1, normal2, domain, dynamic_cast<TensorField3D<T, nDim> &>(*atomicBlocks[0]));
}

template <typename T, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_T<T, nDim>::processCornerGeneric(
    int normalX, int normalY, int normalZ, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, normalZ, domain,
        dynamic_cast<TensorField3D<T, nDim> &>(*atomicBlocks[0]));
}

/* *************** BoundedReductiveBoxProcessing3D_N ******************************************* */

template <typename T>
void BoundedReductiveBoxProcessingFunctional3D_N<T>::processBulkGeneric(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processBulk(domain, dynamic_cast<NTensorField3D<T> &>(*atomicBlocks[0]));
}

template <typename T>
void BoundedReductiveBoxProcessingFunctional3D_N<T>::processPlaneGeneric(
    int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processPlane(
        direction, orientation, domain, dynamic_cast<NTensorField3D<T> &>(*atomicBlocks[0]));
}

template <typename T>
void BoundedReductiveBoxProcessingFunctional3D_N<T>::processEdgeGeneric(
    int plane, int normal1, int normal2, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processEdge(
        plane, normal1, normal2, domain, dynamic_cast<NTensorField3D<T> &>(*atomicBlocks[0]));
}

template <typename T>
void BoundedReductiveBoxProcessingFunctional3D_N<T>::processCornerGeneric(
    int normalX, int normalY, int normalZ, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, normalZ, domain, dynamic_cast<NTensorField3D<T> &>(*atomicBlocks[0]));
}

/* *************** BoundedMaskedReductiveBoxProcessing3D_N
 * ******************************************* */

template <typename T>
void BoundedMaskedReductiveBoxProcessingFunctional3D_N<T>::processBulkGeneric(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processBulk(
        domain, dynamic_cast<NTensorField3D<T> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<int> &>(*atomicBlocks[1]));
}

template <typename T>
void BoundedMaskedReductiveBoxProcessingFunctional3D_N<T>::processPlaneGeneric(
    int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processPlane(
        direction, orientation, domain, dynamic_cast<NTensorField3D<T> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<int> &>(*atomicBlocks[1]));
}

template <typename T>
void BoundedMaskedReductiveBoxProcessingFunctional3D_N<T>::processEdgeGeneric(
    int plane, int normal1, int normal2, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processEdge(
        plane, normal1, normal2, domain, dynamic_cast<NTensorField3D<T> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<int> &>(*atomicBlocks[1]));
}

template <typename T>
void BoundedMaskedReductiveBoxProcessingFunctional3D_N<T>::processCornerGeneric(
    int normalX, int normalY, int normalZ, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, normalZ, domain, dynamic_cast<NTensorField3D<T> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<int> &>(*atomicBlocks[1]));
}

/* *************** BoundedReductiveBoxProcessing3D_LL ****************************************** */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void BoundedReductiveBoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2>::
    processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    PLB_PRECONDITION(atomicBlocks.size() == 2);
    processBulk(
        domain, dynamic_cast<BlockLattice3D<T1, Descriptor1> &>(*atomicBlocks[0]),
        dynamic_cast<BlockLattice3D<T2, Descriptor2> &>(*atomicBlocks[1]));
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void BoundedReductiveBoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2>::
    processPlaneGeneric(
        int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    PLB_PRECONDITION(atomicBlocks.size() == 2);
    processPlane(
        direction, orientation, domain,
        dynamic_cast<BlockLattice3D<T1, Descriptor1> &>(*atomicBlocks[0]),
        dynamic_cast<BlockLattice3D<T2, Descriptor2> &>(*atomicBlocks[1]));
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void BoundedReductiveBoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2>::
    processEdgeGeneric(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks)
{
    PLB_PRECONDITION(atomicBlocks.size() == 2);
    processEdge(
        plane, normal1, normal2, domain,
        dynamic_cast<BlockLattice3D<T1, Descriptor1> &>(*atomicBlocks[0]),
        dynamic_cast<BlockLattice3D<T2, Descriptor2> &>(*atomicBlocks[1]));
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void BoundedReductiveBoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2>::
    processCornerGeneric(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<AtomicBlock3D *> atomicBlocks)
{
    PLB_PRECONDITION(atomicBlocks.size() == 2);
    processCorner(
        normalX, normalY, normalZ, domain,
        dynamic_cast<BlockLattice3D<T1, Descriptor1> &>(*atomicBlocks[0]),
        dynamic_cast<BlockLattice3D<T2, Descriptor2> &>(*atomicBlocks[1]));
}

/* *************** BoundedReductiveBoxProcessing3D_SS ****************************************** */

template <typename T1, typename T2>
void BoundedReductiveBoxProcessingFunctional3D_SS<T1, T2>::processBulkGeneric(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processBulk(
        domain, dynamic_cast<ScalarField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<ScalarField3D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, typename T2>
void BoundedReductiveBoxProcessingFunctional3D_SS<T1, T2>::processPlaneGeneric(
    int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processPlane(
        direction, orientation, domain, dynamic_cast<ScalarField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<ScalarField3D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, typename T2>
void BoundedReductiveBoxProcessingFunctional3D_SS<T1, T2>::processEdgeGeneric(
    int plane, int normal1, int normal2, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processEdge(
        plane, normal1, normal2, domain, dynamic_cast<ScalarField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<ScalarField3D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, typename T2>
void BoundedReductiveBoxProcessingFunctional3D_SS<T1, T2>::processCornerGeneric(
    int normalX, int normalY, int normalZ, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, normalZ, domain, dynamic_cast<ScalarField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<ScalarField3D<T2> &>(*atomicBlocks[1]));
}

/* *************** BoundedReductiveBoxProcessing3D_TT ****************************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void BoundedReductiveBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2>::processBulkGeneric(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processBulk(
        domain, dynamic_cast<TensorField3D<T1, nDim1> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField3D<T2, nDim2> &>(*atomicBlocks[1]));
}

template <typename T1, int nDim1, typename T2, int nDim2>
void BoundedReductiveBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2>::processPlaneGeneric(
    int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processPlane(
        direction, orientation, domain, dynamic_cast<TensorField3D<T1, nDim1> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField3D<T2, nDim2> &>(*atomicBlocks[1]));
}

template <typename T1, int nDim1, typename T2, int nDim2>
void BoundedReductiveBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2>::processEdgeGeneric(
    int plane, int normal1, int normal2, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processEdge(
        plane, normal1, normal2, domain, dynamic_cast<TensorField3D<T1, nDim1> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField3D<T2, nDim2> &>(*atomicBlocks[1]));
}

template <typename T1, int nDim1, typename T2, int nDim2>
void BoundedReductiveBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2>::processCornerGeneric(
    int normalX, int normalY, int normalZ, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, normalZ, domain,
        dynamic_cast<TensorField3D<T1, nDim1> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField3D<T2, nDim2> &>(*atomicBlocks[1]));
}

/* *************** BoundedReductiveBoxProcessing3D_NN ****************************************** */

template <typename T1, typename T2>
void BoundedReductiveBoxProcessingFunctional3D_NN<T1, T2>::processBulkGeneric(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processBulk(
        domain, dynamic_cast<NTensorField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, typename T2>
void BoundedReductiveBoxProcessingFunctional3D_NN<T1, T2>::processPlaneGeneric(
    int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processPlane(
        direction, orientation, domain, dynamic_cast<NTensorField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, typename T2>
void BoundedReductiveBoxProcessingFunctional3D_NN<T1, T2>::processEdgeGeneric(
    int plane, int normal1, int normal2, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processEdge(
        plane, normal1, normal2, domain, dynamic_cast<NTensorField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, typename T2>
void BoundedReductiveBoxProcessingFunctional3D_NN<T1, T2>::processCornerGeneric(
    int normalX, int normalY, int normalZ, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, normalZ, domain, dynamic_cast<NTensorField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]));
}

/* *************** BoundedMaskedReductiveBoxProcessing3D_NN
 * ******************************************* */

template <typename T1, typename T2>
void BoundedMaskedReductiveBoxProcessingFunctional3D_NN<T1, T2>::processBulkGeneric(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processBulk(
        domain, dynamic_cast<NTensorField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]),
        dynamic_cast<NTensorField3D<int> &>(*atomicBlocks[2]));
}

template <typename T1, typename T2>
void BoundedMaskedReductiveBoxProcessingFunctional3D_NN<T1, T2>::processPlaneGeneric(
    int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processPlane(
        direction, orientation, domain, dynamic_cast<NTensorField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]),
        dynamic_cast<NTensorField3D<int> &>(*atomicBlocks[2]));
}

template <typename T1, typename T2>
void BoundedMaskedReductiveBoxProcessingFunctional3D_NN<T1, T2>::processEdgeGeneric(
    int plane, int normal1, int normal2, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processEdge(
        plane, normal1, normal2, domain, dynamic_cast<NTensorField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]),
        dynamic_cast<NTensorField3D<int> &>(*atomicBlocks[2]));
}

template <typename T1, typename T2>
void BoundedMaskedReductiveBoxProcessingFunctional3D_NN<T1, T2>::processCornerGeneric(
    int normalX, int normalY, int normalZ, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, normalZ, domain, dynamic_cast<NTensorField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]),
        dynamic_cast<NTensorField3D<int> &>(*atomicBlocks[2]));
}

/* *************** BoundedReductiveBoxProcessing3D_ST ****************************************** */

template <typename T1, typename T2, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_ST<T1, T2, nDim>::processBulkGeneric(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processBulk(
        domain, dynamic_cast<ScalarField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField3D<T2, nDim> &>(*atomicBlocks[1]));
}

template <typename T1, typename T2, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_ST<T1, T2, nDim>::processPlaneGeneric(
    int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processPlane(
        direction, orientation, domain, dynamic_cast<ScalarField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField3D<T2, nDim> &>(*atomicBlocks[1]));
}

template <typename T1, typename T2, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_ST<T1, T2, nDim>::processEdgeGeneric(
    int plane, int normal1, int normal2, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processEdge(
        plane, normal1, normal2, domain, dynamic_cast<ScalarField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField3D<T2, nDim> &>(*atomicBlocks[1]));
}

template <typename T1, typename T2, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_ST<T1, T2, nDim>::processCornerGeneric(
    int normalX, int normalY, int normalZ, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, normalZ, domain, dynamic_cast<ScalarField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField3D<T2, nDim> &>(*atomicBlocks[1]));
}

/* *************** BoundedReductiveBoxProcessing3D_SN ****************************************** */

template <typename T1, typename T2>
void BoundedReductiveBoxProcessingFunctional3D_SN<T1, T2>::processBulkGeneric(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processBulk(
        domain, dynamic_cast<ScalarField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, typename T2>
void BoundedReductiveBoxProcessingFunctional3D_SN<T1, T2>::processPlaneGeneric(
    int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processPlane(
        direction, orientation, domain, dynamic_cast<ScalarField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, typename T2>
void BoundedReductiveBoxProcessingFunctional3D_SN<T1, T2>::processEdgeGeneric(
    int plane, int normal1, int normal2, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processEdge(
        plane, normal1, normal2, domain, dynamic_cast<ScalarField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, typename T2>
void BoundedReductiveBoxProcessingFunctional3D_SN<T1, T2>::processCornerGeneric(
    int normalX, int normalY, int normalZ, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, normalZ, domain, dynamic_cast<ScalarField3D<T1> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]));
}

/* *************** BoundedReductiveBoxProcessing3D_LS ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void BoundedReductiveBoxProcessingFunctional3D_LS<T1, Descriptor, T2>::processBulkGeneric(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processBulk(
        domain, dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<ScalarField3D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void BoundedReductiveBoxProcessingFunctional3D_LS<T1, Descriptor, T2>::processPlaneGeneric(
    int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processPlane(
        direction, orientation, domain,
        dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<ScalarField3D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void BoundedReductiveBoxProcessingFunctional3D_LS<T1, Descriptor, T2>::processEdgeGeneric(
    int plane, int normal1, int normal2, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processEdge(
        plane, normal1, normal2, domain,
        dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<ScalarField3D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void BoundedReductiveBoxProcessingFunctional3D_LS<T1, Descriptor, T2>::processCornerGeneric(
    int normalX, int normalY, int normalZ, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, normalZ, domain,
        dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<ScalarField3D<T2> &>(*atomicBlocks[1]));
}

/* *************** BoundedReductiveBoxProcessing3D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim>::processBulkGeneric(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processBulk(
        domain, dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField3D<T2, nDim> &>(*atomicBlocks[1]));
}

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim>::processPlaneGeneric(
    int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processPlane(
        direction, orientation, domain,
        dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField3D<T2, nDim> &>(*atomicBlocks[1]));
}

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim>::processEdgeGeneric(
    int plane, int normal1, int normal2, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processEdge(
        plane, normal1, normal2, domain,
        dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField3D<T2, nDim> &>(*atomicBlocks[1]));
}

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim>::processCornerGeneric(
    int normalX, int normalY, int normalZ, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, normalZ, domain,
        dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<TensorField3D<T2, nDim> &>(*atomicBlocks[1]));
}

/* *************** BoundedReductiveBoxProcessing3D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void BoundedReductiveBoxProcessingFunctional3D_LN<T1, Descriptor, T2>::processBulkGeneric(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processBulk(
        domain, dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void BoundedReductiveBoxProcessingFunctional3D_LN<T1, Descriptor, T2>::processPlaneGeneric(
    int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processPlane(
        direction, orientation, domain,
        dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void BoundedReductiveBoxProcessingFunctional3D_LN<T1, Descriptor, T2>::processEdgeGeneric(
    int plane, int normal1, int normal2, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processEdge(
        plane, normal1, normal2, domain,
        dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]));
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void BoundedReductiveBoxProcessingFunctional3D_LN<T1, Descriptor, T2>::processCornerGeneric(
    int normalX, int normalY, int normalZ, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, normalZ, domain,
        dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]));
}

/* *************** BoundedMaskedReductiveBoxProcessing3D_LN
 * ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void BoundedMaskedReductiveBoxProcessingFunctional3D_LN<T1, Descriptor, T2>::processBulkGeneric(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processBulk(
        domain, dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]),
        dynamic_cast<NTensorField3D<int> &>(*atomicBlocks[2]));
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void BoundedMaskedReductiveBoxProcessingFunctional3D_LN<T1, Descriptor, T2>::processPlaneGeneric(
    int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processPlane(
        direction, orientation, domain,
        dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]),
        dynamic_cast<NTensorField3D<int> &>(*atomicBlocks[2]));
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void BoundedMaskedReductiveBoxProcessingFunctional3D_LN<T1, Descriptor, T2>::processEdgeGeneric(
    int plane, int normal1, int normal2, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processEdge(
        plane, normal1, normal2, domain,
        dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]),
        dynamic_cast<NTensorField3D<int> &>(*atomicBlocks[2]));
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void BoundedMaskedReductiveBoxProcessingFunctional3D_LN<T1, Descriptor, T2>::processCornerGeneric(
    int normalX, int normalY, int normalZ, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    processCorner(
        normalX, normalY, normalZ, domain,
        dynamic_cast<BlockLattice3D<T1, Descriptor> &>(*atomicBlocks[0]),
        dynamic_cast<NTensorField3D<T2> &>(*atomicBlocks[1]),
        dynamic_cast<NTensorField3D<int> &>(*atomicBlocks[2]));
}

/* *************** BoundedReductiveLatticeBoxProcessing3D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void BoundedReductiveLatticeBoxProcessingFunctional3D<T, Descriptor>::processBulkGeneric(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<BlockLattice3D<T, Descriptor> *> lattices(atomicBlocks.size());
    for (pluint iLattice = 0; iLattice < atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice3D<T, Descriptor> *>(atomicBlocks[iLattice]);
    }
    processBulk(domain, lattices);
}

template <typename T, template <typename U> class Descriptor>
void BoundedReductiveLatticeBoxProcessingFunctional3D<T, Descriptor>::processPlaneGeneric(
    int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<BlockLattice3D<T, Descriptor> *> lattices(atomicBlocks.size());
    for (pluint iLattice = 0; iLattice < atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice3D<T, Descriptor> *>(atomicBlocks[iLattice]);
    }
    processPlane(direction, orientation, domain, lattices);
}

template <typename T, template <typename U> class Descriptor>
void BoundedReductiveLatticeBoxProcessingFunctional3D<T, Descriptor>::processEdgeGeneric(
    int plane, int normal1, int normal2, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<BlockLattice3D<T, Descriptor> *> lattices(atomicBlocks.size());
    for (pluint iLattice = 0; iLattice < atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice3D<T, Descriptor> *>(atomicBlocks[iLattice]);
    }
    processEdge(plane, normal1, normal2, domain, lattices);
}

template <typename T, template <typename U> class Descriptor>
void BoundedReductiveLatticeBoxProcessingFunctional3D<T, Descriptor>::processCornerGeneric(
    int normalX, int normalY, int normalZ, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<BlockLattice3D<T, Descriptor> *> lattices(atomicBlocks.size());
    for (pluint iLattice = 0; iLattice < atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice3D<T, Descriptor> *>(atomicBlocks[iLattice]);
    }
    processCorner(normalX, normalY, normalZ, domain, lattices);
}

/* *************** BoundedReductiveScalarFieldBoxProcessing3D
 * ******************************************* */

template <typename T>
void BoundedReductiveScalarFieldBoxProcessingFunctional3D<T>::processBulkGeneric(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<ScalarField3D<T> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[iField]);
    }
    processBulk(domain, fields);
}

template <typename T>
void BoundedReductiveScalarFieldBoxProcessingFunctional3D<T>::processPlaneGeneric(
    int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<ScalarField3D<T> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[iField]);
    }
    processPlane(direction, orientation, domain, fields);
}

template <typename T>
void BoundedReductiveScalarFieldBoxProcessingFunctional3D<T>::processEdgeGeneric(
    int plane, int normal1, int normal2, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<ScalarField3D<T> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[iField]);
    }
    processEdge(plane, normal1, normal2, domain, fields);
}

template <typename T>
void BoundedReductiveScalarFieldBoxProcessingFunctional3D<T>::processCornerGeneric(
    int normalX, int normalY, int normalZ, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<ScalarField3D<T> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[iField]);
    }
    processCorner(normalX, normalY, normalZ, domain, fields);
}

/* *************** BoundedReductiveTensorFieldBoxProcessing3D
 * ******************************************* */

template <typename T, int nDim>
void BoundedReductiveTensorFieldBoxProcessingFunctional3D<T, nDim>::processBulkGeneric(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<TensorField3D<T, nDim> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField3D<T, nDim> *>(atomicBlocks[iField]);
    }
    processBulk(domain, fields);
}

template <typename T, int nDim>
void BoundedReductiveTensorFieldBoxProcessingFunctional3D<T, nDim>::processPlaneGeneric(
    int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<TensorField3D<T, nDim> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField3D<T, nDim> *>(atomicBlocks[iField]);
    }
    processPlane(direction, orientation, domain, fields);
}

template <typename T, int nDim>
void BoundedReductiveTensorFieldBoxProcessingFunctional3D<T, nDim>::processEdgeGeneric(
    int plane, int normal1, int normal2, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<TensorField3D<T, nDim> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField3D<T, nDim> *>(atomicBlocks[iField]);
    }
    processEdge(plane, normal1, normal2, domain, fields);
}

template <typename T, int nDim>
void BoundedReductiveTensorFieldBoxProcessingFunctional3D<T, nDim>::processCornerGeneric(
    int normalX, int normalY, int normalZ, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<TensorField3D<T, nDim> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField3D<T, nDim> *>(atomicBlocks[iField]);
    }
    processCorner(normalX, normalY, normalZ, domain, fields);
}

/* *************** BoundedReductiveNTensorFieldBoxProcessing3D
 * ******************************************* */

template <typename T>
void BoundedReductiveNTensorFieldBoxProcessingFunctional3D<T>::processBulkGeneric(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<NTensorField3D<T> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<NTensorField3D<T> *>(atomicBlocks[iField]);
    }
    processBulk(domain, fields);
}

template <typename T>
void BoundedReductiveNTensorFieldBoxProcessingFunctional3D<T>::processPlaneGeneric(
    int direction, int orientation, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<NTensorField3D<T> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<NTensorField3D<T> *>(atomicBlocks[iField]);
    }
    processPlane(direction, orientation, domain, fields);
}

template <typename T>
void BoundedReductiveNTensorFieldBoxProcessingFunctional3D<T>::processEdgeGeneric(
    int plane, int normal1, int normal2, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<NTensorField3D<T> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<NTensorField3D<T> *>(atomicBlocks[iField]);
    }
    processEdge(plane, normal1, normal2, domain, fields);
}

template <typename T>
void BoundedReductiveNTensorFieldBoxProcessingFunctional3D<T>::processCornerGeneric(
    int normalX, int normalY, int normalZ, Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    std::vector<NTensorField3D<T> *> fields(atomicBlocks.size());
    for (pluint iField = 0; iField < atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<NTensorField3D<T> *>(atomicBlocks[iField]);
    }
    processCorner(normalX, normalY, normalZ, domain, fields);
}

}  // namespace plb

#endif  // REDUCTIVE_DATA_PROCESSING_FUNCTIONAL_3D_HH
