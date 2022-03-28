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

#ifndef MULTI_DATA_PROCESSOR_WRAPPER_2D_HH
#define MULTI_DATA_PROCESSOR_WRAPPER_2D_HH

#include "atomicBlock/dataProcessor2D.h"
#include "core/plbDebug.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiBlockOperations2D.h"
#include "multiBlock/multiDataField2D.h"
#include "multiBlock/multiDataProcessorWrapper2D.h"

namespace plb {

/* *************** LatticeBoxProcessing2D ****************************************** */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    LatticeBoxProcessingFunctional2D<T, Descriptor> *functional, Box2D domain,
    std::vector<MultiBlockLattice2D<T, Descriptor> *> lattices)
{
    std::vector<MultiBlock2D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock2D *>(lattices[iLattice]);
    }
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), multiBlocks);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    LatticeBoxProcessingFunctional2D<T, Descriptor> *functional, Box2D domain,
    std::vector<MultiBlockLattice2D<T, Descriptor> *> lattices, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock2D *>(lattices[iLattice]);
    }
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), multiBlocks, level);
}

/* *************** MultiScalarFieldBoxProcessing2D ******************************* */

template <typename T>
void applyProcessingFunctional(
    ScalarFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiScalarField2D<T> *> fields)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), multiBlocks);
}

template <typename T>
void integrateProcessingFunctional(
    ScalarFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiScalarField2D<T> *> fields, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), multiBlocks, level);
}

/* *************** TensorFieldBoxProcessing2D ******************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    TensorFieldBoxProcessingFunctional2D<T, nDim> *functional, Box2D domain,
    std::vector<MultiTensorField2D<T, nDim> *> fields)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), multiBlocks);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    TensorFieldBoxProcessingFunctional2D<T, nDim> *functional, Box2D domain,
    std::vector<MultiTensorField2D<T, nDim> *> fields, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), multiBlocks, level);
}

/* *************** NTensorFieldBoxProcessing2D ******************************* */

template <typename T>
void applyProcessingFunctional(
    NTensorFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiNTensorField2D<T> *> fields)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), multiBlocks);
}

template <typename T>
void applyProcessingFunctional(
    MaskedNTensorFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiNTensorField2D<T> *> fields, MultiNTensorField2D<int> &mask)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size() + 1);
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    multiBlocks[fields.size()] = dynamic_cast<MultiBlock2D *>(&mask);
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), multiBlocks);
}

template <typename T>
void integrateProcessingFunctional(
    NTensorFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiNTensorField2D<T> *> fields, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), multiBlocks, level);
}

template <typename T>
void integrateProcessingFunctional(
    MaskedNTensorFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiNTensorField2D<T> *> fields, MultiNTensorField2D<int> &mask, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size() + 1);
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    multiBlocks[fields.size()] = dynamic_cast<MultiBlock2D *>(&mask);
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), multiBlocks, level);
}

/* *************** BoxProcessing2D_L ******************************************* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_L<T, Descriptor> *functional, Box2D domain,
    MultiBlockLattice2D<T, Descriptor> &lattice)
{
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), lattice);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_L<T, Descriptor> *functional, Box2D domain,
    MultiBlockLattice2D<T, Descriptor> &lattice, plint level)
{
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), lattice, level);
}

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    MaskedBoxProcessingFunctional2D_L<T, Descriptor> *functional, Box2D domain,
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiNTensorField2D<int> &mask)
{
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), lattice, mask);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    MaskedBoxProcessingFunctional2D_L<T, Descriptor> *functional, Box2D domain,
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiNTensorField2D<int> &mask, plint level)
{
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), lattice, mask, level);
}

/* *************** BoxProcessing2D_S ******************************************* */

template <typename T>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_S<T> *functional, Box2D domain, MultiScalarField2D<T> &field)
{
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), field);
}

template <typename T>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_S<T> *functional, Box2D domain, MultiScalarField2D<T> &field,
    plint level)
{
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), field, level);
}

/* *************** BoxProcessing2D_T ******************************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_T<T, nDim> *functional, Box2D domain,
    MultiTensorField2D<T, nDim> &field)
{
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), field);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_T<T, nDim> *functional, Box2D domain,
    MultiTensorField2D<T, nDim> &field, plint level)
{
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), field, level);
}

/* *************** BoxProcessing2D_N ******************************************* */

template <typename T>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_N<T> *functional, Box2D domain, MultiNTensorField2D<T> &field)
{
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), field);
}

template <typename T>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_N<T> *functional, Box2D domain, MultiNTensorField2D<T> &field,
    plint level)
{
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), field, level);
}

/* *************** MaskedBoxProcessing2D_N ******************************************* */

template <typename T>
void applyProcessingFunctional(
    MaskedBoxProcessingFunctional2D_N<T> *functional, Box2D domain, MultiNTensorField2D<T> &field,
    MultiNTensorField2D<int> &mask)
{
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), field, mask);
}

template <typename T>
void integrateProcessingFunctional(
    MaskedBoxProcessingFunctional2D_N<T> *functional, Box2D domain, MultiNTensorField2D<T> &field,
    MultiNTensorField2D<int> &mask, plint level)
{
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), field, mask, level);
}

/* *************** BoxProcessing2D_LL******************************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor1> &lattice1, MultiBlockLattice2D<T2, Descriptor2> &lattice2)
{
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), lattice1, lattice2);
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor1> &lattice1, MultiBlockLattice2D<T2, Descriptor2> &lattice2,
    plint level)
{
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), lattice1, lattice2, level);
}

/* *************** BoxProcessing2D_SS ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_SS<T1, T2> *functional, Box2D domain, MultiScalarField2D<T1> &field1,
    MultiScalarField2D<T2> &field2)
{
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), field1, field2);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_SS<T1, T2> *functional, Box2D domain, MultiScalarField2D<T1> &field1,
    MultiScalarField2D<T2> &field2, plint level)
{
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), field1, field2, level);
}

/* *************** BoxProcessing2D_TT ****************************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> *functional, Box2D domain,
    MultiTensorField2D<T1, nDim1> &field1, MultiTensorField2D<T2, nDim2> &field2)
{
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), field1, field2);
}

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> *functional, Box2D domain,
    MultiTensorField2D<T1, nDim1> &field1, MultiTensorField2D<T2, nDim2> &field2, plint level)
{
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), field1, field2, level);
}

/* *************** BoxProcessing2D_NN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_NN<T1, T2> *functional, Box2D domain, MultiNTensorField2D<T1> &field1,
    MultiNTensorField2D<T2> &field2)
{
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), field1, field2);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_NN<T1, T2> *functional, Box2D domain, MultiNTensorField2D<T1> &field1,
    MultiNTensorField2D<T2> &field2, plint level)
{
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), field1, field2, level);
}

/* *************** MaskedBoxProcessing2D_NN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    MaskedBoxProcessingFunctional2D_NN<T1, T2> *functional, Box2D domain,
    MultiNTensorField2D<T1> &field1, MultiNTensorField2D<T2> &field2,
    MultiNTensorField2D<int> &mask)
{
    std::vector<MultiBlock2D *> multiBlocks(3);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field2);
    multiBlocks[2] = dynamic_cast<MultiBlock2D *>(&mask);
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), multiBlocks);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    MaskedBoxProcessingFunctional2D_NN<T1, T2> *functional, Box2D domain,
    MultiNTensorField2D<T1> &field1, MultiNTensorField2D<T2> &field2,
    MultiNTensorField2D<int> &mask, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(3);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field2);
    multiBlocks[2] = dynamic_cast<MultiBlock2D *>(&mask);
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), multiBlocks, level);
}

/* *************** BoxProcessing2D_ST ****************************************** */

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_ST<T1, T2, nDim> *functional, Box2D domain,
    MultiScalarField2D<T1> &field1, MultiTensorField2D<T2, nDim> &field2)
{
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), field1, field2);
}

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_ST<T1, T2, nDim> *functional, Box2D domain,
    MultiScalarField2D<T1> &field1, MultiTensorField2D<T2, nDim> &field2, plint level)
{
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), field1, field2, level);
}

/* *************** BoxProcessing2D_SN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_SN<T1, T2> *functional, Box2D domain, MultiScalarField2D<T1> &field1,
    MultiNTensorField2D<T2> &field2)
{
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), field1, field2);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_SN<T1, T2> *functional, Box2D domain, MultiScalarField2D<T1> &field1,
    MultiNTensorField2D<T2> &field2, plint level)
{
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), field1, field2, level);
}

/* *************** BoxProcessing2D_LS ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_LS<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiScalarField2D<T2> &field)
{
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), lattice, field);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_LS<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiScalarField2D<T2> &field, plint level)
{
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), lattice, field, level);
}

/* *************** BoxProcessing2D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiTensorField2D<T2, nDim> &field)
{
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), lattice, field);
}

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiTensorField2D<T2, nDim> &field, plint level)
{
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), lattice, field, level);
}

/* *************** BoxProcessing2D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_LN<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field)
{
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), lattice, field);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_LN<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field, plint level)
{
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), lattice, field, level);
}

/* *************** BoxProcessing2D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    MaskedBoxProcessingFunctional2D_LN<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field,
    MultiNTensorField2D<int> &mask)
{
    std::vector<MultiBlock2D *> multiBlocks(3);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field);
    multiBlocks[2] = dynamic_cast<MultiBlock2D *>(&mask);
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), multiBlocks);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    MaskedBoxProcessingFunctional2D_LN<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field,
    MultiNTensorField2D<int> &mask, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(3);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field);
    multiBlocks[2] = dynamic_cast<MultiBlock2D *>(&mask);
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), multiBlocks, level);
}

/* *************** LatticeDotProcessing2D ****************************************** */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    LatticeDotProcessingFunctional2D<T, Descriptor> *functional, DotList2D const &dotList,
    std::vector<MultiBlockLattice2D<T, Descriptor> *> lattices)
{
    std::vector<MultiBlock2D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock2D *>(lattices[iLattice]);
    }
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), multiBlocks);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    LatticeDotProcessingFunctional2D<T, Descriptor> *functional, DotList2D const &dotList,
    std::vector<MultiBlockLattice2D<T, Descriptor> *> lattices, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock2D *>(lattices[iLattice]);
    }
    addInternalProcessor(DotProcessorGenerator2D(functional, dotList), multiBlocks, level);
}

/* *************** ScalarFieldDotProcessing2D ****************************************** */

template <typename T>
void applyProcessingFunctional(
    ScalarFieldDotProcessingFunctional2D<T> *functional, DotList2D const &dotList,
    std::vector<MultiScalarField2D<T> *> fields)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), multiBlocks);
}

template <typename T>
void integrateProcessingFunctional(
    ScalarFieldDotProcessingFunctional2D<T> *functional, DotList2D const &dotList,
    std::vector<MultiScalarField2D<T> *> fields, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    addInternalProcessor(DotProcessorGenerator2D(functional, dotList), multiBlocks, level);
}

/* *************** TensorFieldDotProcessing2D ****************************************** */

template <typename T, int nDim>
void applyProcessingFunctional(
    TensorFieldDotProcessingFunctional2D<T, nDim> *functional, DotList2D const &dotList,
    std::vector<MultiTensorField2D<T, nDim> *> fields)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), multiBlocks);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    TensorFieldDotProcessingFunctional2D<T, nDim> *functional, DotList2D const &dotList,
    std::vector<MultiTensorField2D<T, nDim> *> fields, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    addInternalProcessor(DotProcessorGenerator2D(functional, dotList), multiBlocks, level);
}

/* *************** NTensorFieldDotProcessing2D ****************************************** */

template <typename T>
void applyProcessingFunctional(
    NTensorFieldDotProcessingFunctional2D<T> *functional, DotList2D const &dotList,
    std::vector<MultiNTensorField2D<T> *> fields)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), multiBlocks);
}

template <typename T>
void integrateProcessingFunctional(
    NTensorFieldDotProcessingFunctional2D<T> *functional, DotList2D const &dotList,
    std::vector<MultiNTensorField2D<T> *> fields, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    addInternalProcessor(DotProcessorGenerator2D(functional, dotList), multiBlocks, level);
}

/* *************** DotProcessing2D_L ******************************************* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    DotProcessingFunctional2D_L<T, Descriptor> *functional, DotList2D const &dotList,
    MultiBlockLattice2D<T, Descriptor> &lattice)
{
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), lattice);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_L<T, Descriptor> *functional, DotList2D const &dotList,
    MultiBlockLattice2D<T, Descriptor> &lattice, plint level)
{
    addInternalProcessor(DotProcessorGenerator2D(functional, dotList), lattice, level);
}

/* *************** DotProcessing2D_S ******************************************* */

template <typename T>
void applyProcessingFunctional(
    DotProcessingFunctional2D_S<T> *functional, DotList2D const &dotList,
    MultiScalarField2D<T> &field)
{
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), field);
}

template <typename T>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_S<T> *functional, DotList2D const &dotList,
    MultiScalarField2D<T> &field, plint level)
{
    addInternalProcessor(DotProcessorGenerator2D(functional, dotList), field, level);
}

/* *************** DotProcessing2D_T ******************************************* */
template <typename T, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional2D_T<T, nDim> *functional, DotList2D const &dotList,
    MultiTensorField2D<T, nDim> &field)
{
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), field);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_T<T, nDim> *functional, DotList2D const &dotList,
    MultiTensorField2D<T, nDim> &field, plint level)
{
    addInternalProcessor(DotProcessorGenerator2D(functional, dotList), field, level);
}

/* *************** DotProcessing2D_N ******************************************* */

template <typename T>
void applyProcessingFunctional(
    DotProcessingFunctional2D_N<T> *functional, DotList2D const &dotList,
    MultiNTensorField2D<T> &field)
{
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), field);
}

template <typename T>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_N<T> *functional, DotList2D const &dotList,
    MultiNTensorField2D<T> &field, plint level)
{
    addInternalProcessor(DotProcessorGenerator2D(functional, dotList), field, level);
}

/* *************** DotProcessing2D_LL******************************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    DotProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> *functional,
    DotList2D const &dotList, MultiBlockLattice2D<T1, Descriptor1> &lattice1,
    MultiBlockLattice2D<T2, Descriptor2> &lattice2)
{
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), lattice1, lattice2);
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> *functional,
    DotList2D const &dotList, MultiBlockLattice2D<T1, Descriptor1> &lattice1,
    MultiBlockLattice2D<T2, Descriptor2> &lattice2, plint level)
{
    addInternalProcessor(DotProcessorGenerator2D(functional, dotList), lattice1, lattice2, level);
}

/* *************** DotProcessing2D_SS ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional2D_SS<T1, T2> *functional, DotList2D const &dotList,
    MultiScalarField2D<T1> &field1, MultiScalarField2D<T2> &field2)
{
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), field1, field2);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_SS<T1, T2> *functional, DotList2D const &dotList,
    MultiScalarField2D<T1> &field1, MultiScalarField2D<T2> &field2, plint level)
{
    addInternalProcessor(DotProcessorGenerator2D(functional, dotList), field1, field2, level);
}

/* *************** DotProcessing2D_TT ****************************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    DotProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> *functional, DotList2D const &dotList,
    MultiTensorField2D<T1, nDim1> &field1, MultiTensorField2D<T2, nDim2> &field2)
{
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), field1, field2);
}

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> *functional, DotList2D const &dotList,
    MultiTensorField2D<T1, nDim1> &field1, MultiTensorField2D<T2, nDim2> &field2, plint level)
{
    addInternalProcessor(DotProcessorGenerator2D(functional, dotList), field1, field2, level);
}

/* *************** DotProcessing2D_NN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional2D_NN<T1, T2> *functional, DotList2D const &dotList,
    MultiNTensorField2D<T1> &field1, MultiNTensorField2D<T2> &field2)
{
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), field1, field2);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_NN<T1, T2> *functional, DotList2D const &dotList,
    MultiNTensorField2D<T1> &field1, MultiNTensorField2D<T2> &field2, plint level)
{
    addInternalProcessor(DotProcessorGenerator2D(functional, dotList), field1, field2, level);
}

/* *************** DotProcessing2D_ST ****************************************** */

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional2D_ST<T1, T2, nDim> *functional, DotList2D const &dotList,
    MultiScalarField2D<T1> &field1, MultiTensorField2D<T2, nDim> &field2)
{
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), field1, field2);
}

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_ST<T1, T2, nDim> *functional, DotList2D const &dotList,
    MultiScalarField2D<T1> &field1, MultiTensorField2D<T2, nDim> &field2, plint level)
{
    addInternalProcessor(DotProcessorGenerator2D(functional, dotList), field1, field2, level);
}

/* *************** DotProcessing2D_SN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional2D_SN<T1, T2> *functional, DotList2D const &dotList,
    MultiScalarField2D<T1> &field1, MultiNTensorField2D<T2> &field2)
{
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), field1, field2);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_SN<T1, T2> *functional, DotList2D const &dotList,
    MultiScalarField2D<T1> &field1, MultiNTensorField2D<T2> &field2, plint level)
{
    addInternalProcessor(DotProcessorGenerator2D(functional, dotList), field1, field2, level);
}

/* *************** DotProcessing2D_LS ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional2D_LS<T1, Descriptor, T2> *functional, DotList2D const &dotList,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiScalarField2D<T2> &field)
{
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), lattice, field);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_LS<T1, Descriptor, T2> *functional, DotList2D const &dotList,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiScalarField2D<T2> &field, plint level)
{
    addInternalProcessor(DotProcessorGenerator2D(functional, dotList), lattice, field, level);
}

/* *************** DotProcessing2D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> *functional, DotList2D const &dotList,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiTensorField2D<T2, nDim> &field)
{
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), lattice, field);
}

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> *functional, DotList2D const &dotList,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiTensorField2D<T2, nDim> &field, plint level)
{
    addInternalProcessor(DotProcessorGenerator2D(functional, dotList), lattice, field, level);
}

/* *************** DotProcessing2D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional2D_LN<T1, Descriptor, T2> *functional, DotList2D const &dotList,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field)
{
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), lattice, field);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_LN<T1, Descriptor, T2> *functional, DotList2D const &dotList,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field, plint level)
{
    addInternalProcessor(DotProcessorGenerator2D(functional, dotList), lattice, field, level);
}

/* *************** BoundedLatticeBoxProcessing2D **************************** */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoundedLatticeBoxProcessingFunctional2D<T, Descriptor> *functional, Box2D domain,
    std::vector<MultiBlockLattice2D<T, Descriptor> *> lattices, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock2D *>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoundedLatticeBoxProcessingFunctional2D<T, Descriptor> *functional, Box2D domain,
    std::vector<MultiBlockLattice2D<T, Descriptor> *> lattices, plint boundaryWidth, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock2D *>(lattices[iLattice]);
    }
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedMultiScalarFieldBoxProcessing2D ************************ */

template <typename T>
void applyProcessingFunctional(
    BoundedScalarFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiScalarField2D<T> *> fields, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T>
void integrateProcessingFunctional(
    BoundedScalarFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiScalarField2D<T> *> fields, plint boundaryWidth, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedTensorFieldBoxProcessing2D ************************ */

template <typename T, int nDim>
void applyProcessingFunctional(
    BoundedTensorFieldBoxProcessingFunctional2D<T, nDim> *functional, Box2D domain,
    std::vector<MultiTensorField2D<T, nDim> *> fields, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    BoundedTensorFieldBoxProcessingFunctional2D<T, nDim> *functional, Box2D domain,
    std::vector<MultiTensorField2D<T, nDim> *> fields, plint boundaryWidth, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedNTensorFieldBoxProcessing2D ************************ */

template <typename T>
void applyProcessingFunctional(
    BoundedNTensorFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiNTensorField2D<T> *> fields, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T>
void integrateProcessingFunctional(
    BoundedNTensorFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiNTensorField2D<T> *> fields, plint boundaryWidth, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_L ********************************* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_L<T, Descriptor> *functional, Box2D domain,
    MultiBlockLattice2D<T, Descriptor> &lattice, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&lattice);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_L<T, Descriptor> *functional, Box2D domain,
    MultiBlockLattice2D<T, Descriptor> &lattice, plint boundaryWidth, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&lattice);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_S ********************************* */

template <typename T>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_S<T> *functional, Box2D domain, MultiScalarField2D<T> &field,
    plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_S<T> *functional, Box2D domain, MultiScalarField2D<T> &field,
    plint boundaryWidth, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_T ********************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_T<T, nDim> *functional, Box2D domain,
    MultiTensorField2D<T, nDim> &field, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_T<T, nDim> *functional, Box2D domain,
    MultiTensorField2D<T, nDim> &field, plint boundaryWidth, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_N ********************************* */

template <typename T>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_N<T> *functional, Box2D domain, MultiNTensorField2D<T> &field,
    plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_N<T> *functional, Box2D domain, MultiNTensorField2D<T> &field,
    plint boundaryWidth, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedMaskedBoxProcessing2D_N ********************************* */

template <typename T>
void applyProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional2D_N<T> *functional, Box2D domain,
    MultiNTensorField2D<T> &field, MultiNTensorField2D<int> &mask, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&mask);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T>
void integrateProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional2D_N<T> *functional, Box2D domain,
    MultiNTensorField2D<T> &field, MultiNTensorField2D<int> &mask, plint boundaryWidth, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&mask);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_LL********************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor1> &lattice1, MultiBlockLattice2D<T2, Descriptor2> &lattice2,
    plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&lattice2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor1> &lattice1, MultiBlockLattice2D<T2, Descriptor2> &lattice2,
    plint boundaryWidth, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&lattice2);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_SS ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_SS<T1, T2> *functional, Box2D domain,
    MultiScalarField2D<T1> &field1, MultiScalarField2D<T2> &field2, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_SS<T1, T2> *functional, Box2D domain,
    MultiScalarField2D<T1> &field1, MultiScalarField2D<T2> &field2, plint boundaryWidth,
    plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field2);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_TT ********************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> *functional, Box2D domain,
    MultiTensorField2D<T1, nDim1> &field1, MultiTensorField2D<T2, nDim2> &field2,
    plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> *functional, Box2D domain,
    MultiTensorField2D<T1, nDim1> &field1, MultiTensorField2D<T2, nDim2> &field2,
    plint boundaryWidth, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field2);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_NN ********************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_NN<T1, T2> *functional, Box2D domain,
    MultiNTensorField2D<T1> &field1, MultiNTensorField2D<T2> &field2, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_NN<T1, T2> *functional, Box2D domain,
    MultiNTensorField2D<T1> &field1, MultiNTensorField2D<T2> &field2, plint boundaryWidth,
    plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field2);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedMaskedBoxProcessing2D_NN ********************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional2D_NN<T1, T2> *functional, Box2D domain,
    MultiNTensorField2D<T1> &field1, MultiNTensorField2D<T2> &field2,
    MultiNTensorField2D<int> &mask, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(3);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field2);
    multiBlocks[2] = dynamic_cast<MultiBlock2D *>(&mask);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional2D_NN<T1, T2> *functional, Box2D domain,
    MultiNTensorField2D<T1> &field1, MultiNTensorField2D<T2> &field2,
    MultiNTensorField2D<int> &mask, plint boundaryWidth, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(3);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field2);
    multiBlocks[2] = dynamic_cast<MultiBlock2D *>(&mask);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_ST ******************************** */

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_ST<T1, T2, nDim> *functional, Box2D domain,
    MultiScalarField2D<T1> &field1, MultiTensorField2D<T2, nDim> &field2, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_ST<T1, T2, nDim> *functional, Box2D domain,
    MultiScalarField2D<T1> &field1, MultiTensorField2D<T2, nDim> &field2, plint boundaryWidth,
    plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field2);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_SN ******************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_SN<T1, T2> *functional, Box2D domain,
    MultiScalarField2D<T1> &field1, MultiNTensorField2D<T2> &field2, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_SN<T1, T2> *functional, Box2D domain,
    MultiScalarField2D<T1> &field1, MultiNTensorField2D<T2> &field2, plint boundaryWidth,
    plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field2);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_LS ******************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_LS<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiScalarField2D<T2> &field,
    plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_LS<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiScalarField2D<T2> &field,
    plint boundaryWidth, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiTensorField2D<T2, nDim> &field,
    plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiTensorField2D<T2, nDim> &field,
    plint boundaryWidth, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_LN<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field,
    plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_LN<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field,
    plint boundaryWidth, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional2D_LN<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field,
    MultiNTensorField2D<int> &mask, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(3);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field);
    multiBlocks[2] = dynamic_cast<MultiBlock2D *>(&mask);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional2D_LN<T1, Descriptor, T2> *functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field,
    MultiNTensorField2D<int> &mask, plint boundaryWidth, plint level)
{
    std::vector<MultiBlock2D *> multiBlocks(3);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field);
    multiBlocks[2] = dynamic_cast<MultiBlock2D *>(&mask);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

}  // namespace plb

#endif  // MULTI_DATA_PROCESSOR_WRAPPER_2D_HH
