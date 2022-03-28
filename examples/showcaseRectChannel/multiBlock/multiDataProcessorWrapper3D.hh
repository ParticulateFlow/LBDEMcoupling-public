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

#ifndef MULTI_DATA_PROCESSOR_WRAPPER_3D_HH
#define MULTI_DATA_PROCESSOR_WRAPPER_3D_HH

#include "atomicBlock/dataProcessor3D.h"
#include "core/plbDebug.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiBlockOperations3D.h"
#include "multiBlock/multiDataField3D.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"

namespace plb {

/* *************** LatticeBoxProcessing3D ****************************************** */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    LatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain,
    std::vector<MultiBlockLattice3D<T, Descriptor> *> lattices)
{
    std::vector<MultiBlock3D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock3D *>(lattices[iLattice]);
    }
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), multiBlocks);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    LatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain,
    std::vector<MultiBlockLattice3D<T, Descriptor> *> lattices, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock3D *>(lattices[iLattice]);
    }
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), multiBlocks, level);
}

/* *************** ScalarFieldBoxProcessing3D ******************************* */

template <typename T>
void applyProcessingFunctional(
    ScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiScalarField3D<T> *> fields)
{
    std::vector<MultiBlock3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D *>(fields[iField]);
    }
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), multiBlocks);
}

template <typename T>
void integrateProcessingFunctional(
    ScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiScalarField3D<T> *> fields, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D *>(fields[iField]);
    }
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), multiBlocks, level);
}

/* *************** MultiTensorFieldBoxProcessing3D ******************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    TensorFieldBoxProcessingFunctional3D<T, nDim> *functional, Box3D domain,
    std::vector<MultiTensorField3D<T, nDim> *> fields)
{
    std::vector<MultiBlock3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D *>(fields[iField]);
    }
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), multiBlocks);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    TensorFieldBoxProcessingFunctional3D<T, nDim> *functional, Box3D domain,
    std::vector<MultiTensorField3D<T, nDim> *> fields, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D *>(fields[iField]);
    }
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), multiBlocks, level);
}

/* *************** MultiNTensorFieldBoxProcessing3D ******************************* */

template <typename T>
void applyProcessingFunctional(
    NTensorFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiNTensorField3D<T> *> fields)
{
    std::vector<MultiBlock3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D *>(fields[iField]);
    }
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), multiBlocks);
}

template <typename T>
void applyProcessingFunctional(
    MaskedNTensorFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiNTensorField3D<T> *> fields, MultiNTensorField3D<int> &mask)
{
    std::vector<MultiBlock3D *> multiBlocks(fields.size() + 1);
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D *>(fields[iField]);
    }
    multiBlocks[fields.size()] = dynamic_cast<MultiBlock3D *>(&mask);
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), multiBlocks);
}

template <typename T>
void integrateProcessingFunctional(
    NTensorFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiNTensorField3D<T> *> fields, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D *>(fields[iField]);
    }
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), multiBlocks, level);
}

template <typename T>
void integrateProcessingFunctional(
    MaskedNTensorFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiNTensorField3D<T> *> fields, MultiNTensorField3D<int> &mask, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(fields.size() + 1);
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D *>(fields[iField]);
    }
    multiBlocks[fields.size()] = dynamic_cast<MultiBlock3D *>(&mask);
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), multiBlocks, level);
}

/* *************** BoxProcessing3D_L ******************************************* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), lattice);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    MultiBlockLattice3D<T, Descriptor> &lattice, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), lattice, level);
}

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    MaskedBoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiNTensorField3D<int> &mask)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), lattice, mask);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    MaskedBoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiNTensorField3D<int> &mask, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), lattice, mask, level);
}

/* *************** BoxProcessing3D_S ******************************************* */

template <typename T>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_S<T> *functional, Box3D domain, MultiScalarField3D<T> &field)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), field);
}

template <typename T>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_S<T> *functional, Box3D domain, MultiScalarField3D<T> &field,
    plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), field, level);
}

/* *************** BoxProcessing3D_T ******************************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_T<T, nDim> *functional, Box3D domain,
    MultiTensorField3D<T, nDim> &field)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), field);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_T<T, nDim> *functional, Box3D domain,
    MultiTensorField3D<T, nDim> &field, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), field, level);
}

/* *************** BoxProcessing3D_N ******************************************* */

template <typename T>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_N<T> *functional, Box3D domain, MultiNTensorField3D<T> &field)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), field);
}

template <typename T>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_N<T> *functional, Box3D domain, MultiNTensorField3D<T> &field,
    plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), field, level);
}

/* *************** MaskedBoxProcessing3D_N ******************************************* */

template <typename T>
void applyProcessingFunctional(
    MaskedBoxProcessingFunctional3D_N<T> *functional, Box3D domain, MultiNTensorField3D<T> &field,
    MultiNTensorField3D<int> &mask)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), field, mask);
}

template <typename T>
void integrateProcessingFunctional(
    MaskedBoxProcessingFunctional3D_N<T> *functional, Box3D domain, MultiNTensorField3D<T> &field,
    MultiNTensorField3D<int> &mask, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), field, mask, level);
}

/* *************** BoxProcessing3D_LL******************************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor1> &lattice1, MultiBlockLattice3D<T2, Descriptor2> &lattice2)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), lattice1, lattice2);
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor1> &lattice1, MultiBlockLattice3D<T2, Descriptor2> &lattice2,
    plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), lattice1, lattice2, level);
}

/* *************** BoxProcessing3D_SS ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain, MultiScalarField3D<T1> &field1,
    MultiScalarField3D<T2> &field2)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), field1, field2);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain, MultiScalarField3D<T1> &field1,
    MultiScalarField3D<T2> &field2, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), field1, field2, level);
}

/* *************** BoxProcessing3D_TT ****************************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    MultiTensorField3D<T1, nDim1> &field1, MultiTensorField3D<T2, nDim2> &field2)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), field1, field2);
}

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    MultiTensorField3D<T1, nDim1> &field1, MultiTensorField3D<T2, nDim2> &field2, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), field1, field2, level);
}

/* *************** BoxProcessing3D_NN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_NN<T1, T2> *functional, Box3D domain, MultiNTensorField3D<T1> &field1,
    MultiNTensorField3D<T2> &field2)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), field1, field2);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_NN<T1, T2> *functional, Box3D domain, MultiNTensorField3D<T1> &field1,
    MultiNTensorField3D<T2> &field2, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), field1, field2, level);
}

/* *************** MaskedBoxProcessing3D_NN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    MaskedBoxProcessingFunctional3D_NN<T1, T2> *functional, Box3D domain,
    MultiNTensorField3D<T1> &field1, MultiNTensorField3D<T2> &field2,
    MultiNTensorField3D<int> &mask)
{
    std::vector<MultiBlock3D *> multiBlocks(3);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field2);
    multiBlocks[2] = dynamic_cast<MultiBlock3D *>(&mask);
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), multiBlocks);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    MaskedBoxProcessingFunctional3D_NN<T1, T2> *functional, Box3D domain,
    MultiNTensorField3D<T1> &field1, MultiNTensorField3D<T2> &field2,
    MultiNTensorField3D<int> &mask, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(3);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field2);
    multiBlocks[2] = dynamic_cast<MultiBlock3D *>(&mask);
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), multiBlocks, level);
}

/* *************** BoxProcessing3D_ST ****************************************** */

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain,
    MultiScalarField3D<T1> &field1, MultiTensorField3D<T2, nDim> &field2)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), field1, field2);
}

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain,
    MultiScalarField3D<T1> &field1, MultiTensorField3D<T2, nDim> &field2, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), field1, field2, level);
}

/* *************** BoxProcessing3D_SN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_SN<T1, T2> *functional, Box3D domain, MultiScalarField3D<T1> &field1,
    MultiNTensorField3D<T2> &field2)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), field1, field2);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_SN<T1, T2> *functional, Box3D domain, MultiScalarField3D<T1> &field1,
    MultiNTensorField3D<T2> &field2, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), field1, field2, level);
}

/* *************** BoxProcessing3D_LS ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiScalarField3D<T2> &field)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), lattice, field);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiScalarField3D<T2> &field, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), lattice, field, level);
}

/* *************** BoxProcessing3D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiTensorField3D<T2, nDim> &field)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), lattice, field);
}

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiTensorField3D<T2, nDim> &field, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), lattice, field, level);
}

/* *************** BoxProcessing3D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiNTensorField3D<T2> &field)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), lattice, field);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiNTensorField3D<T2> &field, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), lattice, field, level);
}

/* *************** MaskedBoxProcessing3D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    MaskedBoxProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiNTensorField3D<T2> &field,
    MultiNTensorField3D<int> &mask)
{
    std::vector<MultiBlock3D *> multiBlocks(3);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field);
    multiBlocks[2] = dynamic_cast<MultiBlock3D *>(&mask);
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), multiBlocks);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    MaskedBoxProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiNTensorField3D<T2> &field,
    MultiNTensorField3D<int> &mask, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(3);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field);
    multiBlocks[2] = dynamic_cast<MultiBlock3D *>(&mask);
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), multiBlocks, level);
}

/* *************** LatticeDotProcessing3D ****************************************** */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    LatticeDotProcessingFunctional3D<T, Descriptor> *functional, DotList3D const &dotList,
    std::vector<MultiBlockLattice3D<T, Descriptor> *> lattices)
{
    std::vector<MultiBlock3D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock3D *>(lattices[iLattice]);
    }
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), multiBlocks);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    LatticeDotProcessingFunctional3D<T, Descriptor> *functional, DotList3D const &dotList,
    std::vector<MultiBlockLattice3D<T, Descriptor> *> lattices, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock3D *>(lattices[iLattice]);
    }
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), multiBlocks, level);
}

/* *************** ScalarFieldDotProcessing3D ****************************************** */

template <typename T>
void applyProcessingFunctional(
    ScalarFieldDotProcessingFunctional3D<T> *functional, DotList3D const &dotList,
    std::vector<MultiScalarField3D<T> *> fields)
{
    std::vector<MultiBlock3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D *>(fields[iField]);
    }
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), multiBlocks);
}

template <typename T>
void integrateProcessingFunctional(
    ScalarFieldDotProcessingFunctional3D<T> *functional, DotList3D const &dotList,
    std::vector<MultiScalarField3D<T> *> fields, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D *>(fields[iField]);
    }
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), multiBlocks, level);
}

/* *************** TensorFieldDotProcessing3D ****************************************** */

template <typename T, int nDim>
void applyProcessingFunctional(
    TensorFieldDotProcessingFunctional3D<T, nDim> *functional, DotList3D const &dotList,
    std::vector<MultiTensorField3D<T, nDim> *> fields)
{
    std::vector<MultiBlock3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D *>(fields[iField]);
    }
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), multiBlocks);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    TensorFieldDotProcessingFunctional3D<T, nDim> *functional, DotList3D const &dotList,
    std::vector<MultiTensorField3D<T, nDim> *> fields, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D *>(fields[iField]);
    }
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), multiBlocks, level);
}

/* *************** NTensorFieldDotProcessing3D ****************************************** */

template <typename T>
void applyProcessingFunctional(
    NTensorFieldDotProcessingFunctional3D<T> *functional, DotList3D const &dotList,
    std::vector<MultiNTensorField3D<T> *> fields)
{
    std::vector<MultiBlock3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D *>(fields[iField]);
    }
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), multiBlocks);
}

template <typename T>
void integrateProcessingFunctional(
    NTensorFieldDotProcessingFunctional3D<T> *functional, DotList3D const &dotList,
    std::vector<MultiNTensorField3D<T> *> fields, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D *>(fields[iField]);
    }
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), multiBlocks, level);
}

/* *************** DotProcessing3D_L ******************************************* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    DotProcessingFunctional3D_L<T, Descriptor> *functional, DotList3D const &dotList,
    MultiBlockLattice3D<T, Descriptor> &lattice)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), lattice);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_L<T, Descriptor> *functional, DotList3D const &dotList,
    MultiBlockLattice3D<T, Descriptor> &lattice, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), lattice, level);
}

/* *************** DotProcessing3D_S ******************************************* */

template <typename T>
void applyProcessingFunctional(
    DotProcessingFunctional3D_S<T> *functional, DotList3D const &dotList,
    MultiScalarField3D<T> &field)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), field);
}

template <typename T>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_S<T> *functional, DotList3D const &dotList,
    MultiScalarField3D<T> &field, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), field, level);
}

/* *************** DotProcessing3D_T ******************************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional3D_T<T, nDim> *functional, DotList3D const &dotList,
    MultiTensorField3D<T, nDim> &field)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), field);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_T<T, nDim> *functional, DotList3D const &dotList,
    MultiTensorField3D<T, nDim> &field, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), field, level);
}

/* *************** DotProcessing3D_N ******************************************* */

template <typename T>
void applyProcessingFunctional(
    DotProcessingFunctional3D_N<T> *functional, DotList3D const &dotList,
    MultiNTensorField3D<T> &field)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), field);
}

template <typename T>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_N<T> *functional, DotList3D const &dotList,
    MultiNTensorField3D<T> &field, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), field, level);
}

/* *************** DotProcessing3D_LL******************************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional,
    DotList3D const &dotList, MultiBlockLattice3D<T1, Descriptor1> &lattice1,
    MultiBlockLattice3D<T2, Descriptor2> &lattice2)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), lattice1, lattice2);
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional,
    DotList3D const &dotList, MultiBlockLattice3D<T1, Descriptor1> &lattice1,
    MultiBlockLattice3D<T2, Descriptor2> &lattice2, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), lattice1, lattice2, level);
}

/* *************** DotProcessing3D_SS ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_SS<T1, T2> *functional, DotList3D const &dotList,
    MultiScalarField3D<T1> &field1, MultiScalarField3D<T2> &field2)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), field1, field2);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_SS<T1, T2> *functional, DotList3D const &dotList,
    MultiScalarField3D<T1> &field1, MultiScalarField3D<T2> &field2, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), field1, field2, level);
}

/* *************** DotProcessing3D_TT ****************************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, DotList3D const &dotList,
    MultiTensorField3D<T1, nDim1> &field1, MultiTensorField3D<T2, nDim2> &field2)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), field1, field2);
}

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, DotList3D const &dotList,
    MultiTensorField3D<T1, nDim1> &field1, MultiTensorField3D<T2, nDim2> &field2, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), field1, field2, level);
}

/* *************** DotProcessing3D_NN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_NN<T1, T2> *functional, DotList3D const &dotList,
    MultiNTensorField3D<T1> &field1, MultiNTensorField3D<T2> &field2)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), field1, field2);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_NN<T1, T2> *functional, DotList3D const &dotList,
    MultiNTensorField3D<T1> &field1, MultiNTensorField3D<T2> &field2, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), field1, field2, level);
}

/* *************** DotProcessing3D_ST ****************************************** */

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional3D_ST<T1, T2, nDim> *functional, DotList3D const &dotList,
    MultiScalarField3D<T1> &field1, MultiTensorField3D<T2, nDim> &field2)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), field1, field2);
}

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_ST<T1, T2, nDim> *functional, DotList3D const &dotList,
    MultiScalarField3D<T1> &field1, MultiTensorField3D<T2, nDim> &field2, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), field1, field2, level);
}

/* *************** DotProcessing3D_SN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_SN<T1, T2> *functional, DotList3D const &dotList,
    MultiScalarField3D<T1> &field1, MultiNTensorField3D<T2> &field2)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), field1, field2);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_SN<T1, T2> *functional, DotList3D const &dotList,
    MultiScalarField3D<T1> &field1, MultiNTensorField3D<T2> &field2, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), field1, field2, level);
}

/* *************** DotProcessing3D_LS ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, DotList3D const &dotList,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiScalarField3D<T2> &field)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), lattice, field);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, DotList3D const &dotList,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiScalarField3D<T2> &field, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), lattice, field, level);
}

/* *************** DotProcessing3D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, DotList3D const &dotList,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiTensorField3D<T2, nDim> &field)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), lattice, field);
}

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, DotList3D const &dotList,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiTensorField3D<T2, nDim> &field, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), lattice, field, level);
}

/* *************** DotProcessing3D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, DotList3D const &dotList,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiNTensorField3D<T2> &field)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), lattice, field);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, DotList3D const &dotList,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiNTensorField3D<T2> &field, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), lattice, field, level);
}

/* *************** BoundedLatticeBoxProcessing3D **************************** */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoundedLatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain,
    std::vector<MultiBlockLattice3D<T, Descriptor> *> lattices, plint boundaryWidth)
{
    std::vector<MultiBlock3D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock3D *>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoundedLatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain,
    std::vector<MultiBlockLattice3D<T, Descriptor> *> lattices, plint boundaryWidth, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock3D *>(lattices[iLattice]);
    }
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedScalarFieldBoxProcessing3D ************************ */

template <typename T>
void applyProcessingFunctional(
    BoundedScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiScalarField3D<T> *> fields, plint boundaryWidth)
{
    std::vector<MultiBlock3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D *>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T>
void integrateProcessingFunctional(
    BoundedScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiScalarField3D<T> *> fields, plint boundaryWidth, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D *>(fields[iField]);
    }
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedMultiTensorFieldBoxProcessing3D ************************ */

template <typename T, int nDim>
void applyProcessingFunctional(
    BoundedTensorFieldBoxProcessingFunctional3D<T, nDim> *functional, Box3D domain,
    std::vector<MultiTensorField3D<T, nDim> *> fields, plint boundaryWidth)
{
    std::vector<MultiBlock3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D *>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    BoundedTensorFieldBoxProcessingFunctional3D<T, nDim> *functional, Box3D domain,
    std::vector<MultiTensorField3D<T, nDim> *> fields, plint boundaryWidth, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D *>(fields[iField]);
    }
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedMultiNTensorFieldBoxProcessing3D ************************ */

template <typename T>
void applyProcessingFunctional(
    BoundedNTensorFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiNTensorField3D<T> *> fields, plint boundaryWidth)
{
    std::vector<MultiBlock3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D *>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T>
void integrateProcessingFunctional(
    BoundedNTensorFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiNTensorField3D<T> *> fields, plint boundaryWidth, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D *>(fields[iField]);
    }
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_L ********************************* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    MultiBlockLattice3D<T, Descriptor> &lattice, plint boundaryWidth)
{
    std::vector<MultiBlock3D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&lattice);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    MultiBlockLattice3D<T, Descriptor> &lattice, plint boundaryWidth, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&lattice);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_S ********************************* */

template <typename T>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_S<T> *functional, Box3D domain, MultiScalarField3D<T> &field,
    plint boundaryWidth)
{
    std::vector<MultiBlock3D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_S<T> *functional, Box3D domain, MultiScalarField3D<T> &field,
    plint boundaryWidth, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&field);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_T ********************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_T<T, nDim> *functional, Box3D domain,
    MultiTensorField3D<T, nDim> &field, plint boundaryWidth)
{
    std::vector<MultiBlock3D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_T<T, nDim> *functional, Box3D domain,
    MultiTensorField3D<T, nDim> &field, plint boundaryWidth, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&field);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_N ********************************* */

template <typename T>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_N<T> *functional, Box3D domain, MultiNTensorField3D<T> &field,
    plint boundaryWidth)
{
    std::vector<MultiBlock3D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_N<T> *functional, Box3D domain, MultiNTensorField3D<T> &field,
    plint boundaryWidth, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&field);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedMaskedBoxProcessing3D_N ********************************* */

template <typename T>
void applyProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional3D_N<T> *functional, Box3D domain,
    MultiNTensorField3D<T> &field, MultiNTensorField3D<int> &mask, plint boundaryWidth)
{
    std::vector<MultiBlock3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&field);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&mask);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T>
void integrateProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional3D_N<T> *functional, Box3D domain,
    MultiNTensorField3D<T> &field, MultiNTensorField3D<int> &mask, plint boundaryWidth, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&field);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&mask);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_LL********************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor1> &lattice1, MultiBlockLattice3D<T2, Descriptor2> &lattice2,
    plint boundaryWidth)
{
    std::vector<MultiBlock3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&lattice2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor1> &lattice1, MultiBlockLattice3D<T2, Descriptor2> &lattice2,
    plint boundaryWidth, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&lattice2);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_SS ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain,
    MultiScalarField3D<T1> &field1, MultiScalarField3D<T2> &field2, plint boundaryWidth)
{
    std::vector<MultiBlock3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain,
    MultiScalarField3D<T1> &field1, MultiScalarField3D<T2> &field2, plint boundaryWidth,
    plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field2);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_TT ********************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    MultiTensorField3D<T1, nDim1> &field1, MultiTensorField3D<T2, nDim2> &field2,
    plint boundaryWidth)
{
    std::vector<MultiBlock3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    MultiTensorField3D<T1, nDim1> &field1, MultiTensorField3D<T2, nDim2> &field2,
    plint boundaryWidth, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field2);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_NN ********************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_NN<T1, T2> *functional, Box3D domain,
    MultiNTensorField3D<T1> &field1, MultiNTensorField3D<T2> &field2, plint boundaryWidth)
{
    std::vector<MultiBlock3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_NN<T1, T2> *functional, Box3D domain,
    MultiNTensorField3D<T1> &field1, MultiNTensorField3D<T2> &field2, plint boundaryWidth,
    plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field2);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedMaskedBoxProcessing3D_NN ********************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional3D_NN<T1, T2> *functional, Box3D domain,
    MultiNTensorField3D<T1> &field1, MultiNTensorField3D<T2> &field2,
    MultiNTensorField3D<int> &mask, plint boundaryWidth)
{
    std::vector<MultiBlock3D *> multiBlocks(3);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field2);
    multiBlocks[2] = dynamic_cast<MultiBlock3D *>(&mask);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional3D_NN<T1, T2> *functional, Box3D domain,
    MultiNTensorField3D<T1> &field1, MultiNTensorField3D<T2> &field2,
    MultiNTensorField3D<int> &mask, plint boundaryWidth, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(3);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field2);
    multiBlocks[2] = dynamic_cast<MultiBlock3D *>(&mask);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_ST ******************************** */

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain,
    MultiScalarField3D<T1> &field1, MultiTensorField3D<T2, nDim> &field2, plint boundaryWidth)
{
    std::vector<MultiBlock3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain,
    MultiScalarField3D<T1> &field1, MultiTensorField3D<T2, nDim> &field2, plint boundaryWidth,
    plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field2);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_SN ******************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_SN<T1, T2> *functional, Box3D domain,
    MultiScalarField3D<T1> &field1, MultiNTensorField3D<T2> &field2, plint boundaryWidth)
{
    std::vector<MultiBlock3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_SN<T1, T2> *functional, Box3D domain,
    MultiScalarField3D<T1> &field1, MultiNTensorField3D<T2> &field2, plint boundaryWidth,
    plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field2);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_LS ******************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiScalarField3D<T2> &field,
    plint boundaryWidth)
{
    std::vector<MultiBlock3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiScalarField3D<T2> &field,
    plint boundaryWidth, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiTensorField3D<T2, nDim> &field,
    plint boundaryWidth)
{
    std::vector<MultiBlock3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiTensorField3D<T2, nDim> &field,
    plint boundaryWidth, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiNTensorField3D<T2> &field,
    plint boundaryWidth)
{
    std::vector<MultiBlock3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiNTensorField3D<T2> &field,
    plint boundaryWidth, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

/* *************** BoundedMaskedBoxProcessing3D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiNTensorField3D<T2> &field,
    MultiNTensorField3D<int> &mask, plint boundaryWidth)
{
    std::vector<MultiBlock3D *> multiBlocks(3);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field);
    multiBlocks[2] = dynamic_cast<MultiBlock3D *>(&mask);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoundedMaskedBoxProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, Box3D domain,
    MultiBlockLattice3D<T1, Descriptor> &lattice, MultiNTensorField3D<T2> &field,
    MultiNTensorField3D<int> &mask, plint boundaryWidth, plint level)
{
    std::vector<MultiBlock3D *> multiBlocks(3);
    multiBlocks[0] = dynamic_cast<MultiBlock3D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D *>(&field);
    multiBlocks[2] = dynamic_cast<MultiBlock3D *>(&mask);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

}  // namespace plb

#endif  // MULTI_DATA_PROCESSOR_WRAPPER_3D_HH
