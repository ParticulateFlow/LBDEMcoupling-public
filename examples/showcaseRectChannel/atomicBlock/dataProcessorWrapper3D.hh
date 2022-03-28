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

#ifndef DATA_PROCESSOR_WRAPPER_3D_HH
#define DATA_PROCESSOR_WRAPPER_3D_HH

#include "atomicBlock/atomicBlockOperations3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "core/plbDebug.h"

namespace plb {

/* *************** LatticeBoxProcessing3D ****************************************** */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    LatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain,
    std::vector<BlockLattice3D<T, Descriptor> *> lattices)
{
    std::vector<AtomicBlock3D *> atomicBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        atomicBlocks[iLattice] = dynamic_cast<AtomicBlock3D *>(lattices[iLattice]);
    }
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), atomicBlocks);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    LatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain,
    std::vector<BlockLattice3D<T, Descriptor> *> lattices, plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        atomicBlocks[iLattice] = dynamic_cast<AtomicBlock3D *>(lattices[iLattice]);
    }
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), atomicBlocks, level);
}

/* *************** ScalarFieldBoxProcessing3D ******************************* */

template <typename T>
void applyProcessingFunctional(
    ScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<ScalarField3D<T> *> fields)
{
    std::vector<AtomicBlock3D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D *>(fields[iField]);
    }
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), atomicBlocks);
}

template <typename T>
void integrateProcessingFunctional(
    ScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<ScalarField3D<T> *> fields, plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D *>(fields[iField]);
    }
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), atomicBlocks, level);
}

/* *************** TensorFieldBoxProcessing3D ******************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    TensorFieldBoxProcessingFunctional3D<T, nDim> *functional, Box3D domain,
    std::vector<TensorField3D<T, nDim> *> fields)
{
    std::vector<AtomicBlock3D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D *>(fields[iField]);
    }
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), atomicBlocks);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    TensorFieldBoxProcessingFunctional3D<T, nDim> *functional, Box3D domain,
    std::vector<TensorField3D<T, nDim> *> fields, plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D *>(fields[iField]);
    }
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), atomicBlocks, level);
}

/* *************** NTensorFieldBoxProcessing3D ******************************* */

template <typename T>
void applyProcessingFunctional(
    NTensorFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<NTensorField3D<T> *> fields)
{
    std::vector<AtomicBlock3D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D *>(fields[iField]);
    }
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), atomicBlocks);
}

template <typename T>
void integrateProcessingFunctional(
    NTensorFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<NTensorField3D<T> *> fields, plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D *>(fields[iField]);
    }
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), atomicBlocks, level);
}

/* *************** BoxProcessing3D_L ******************************************* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    BlockLattice3D<T, Descriptor> &lattice)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), lattice);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    BlockLattice3D<T, Descriptor> &lattice, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), lattice, level);
}

/* *************** BoxProcessing3D_S ******************************************* */

template <typename T>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_S<T> *functional, Box3D domain, ScalarField3D<T> &field)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), field);
}

template <typename T>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_S<T> *functional, Box3D domain, ScalarField3D<T> &field, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), field, level);
}

/* *************** BoxProcessing3D_T ******************************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_T<T, nDim> *functional, Box3D domain, TensorField3D<T, nDim> &field)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), field);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_T<T, nDim> *functional, Box3D domain, TensorField3D<T, nDim> &field,
    plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), field, level);
}

/* *************** BoxProcessing3D_N ******************************************* */

template <typename T>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_N<T> *functional, Box3D domain, NTensorField3D<T> &field)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), field);
}

template <typename T>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_N<T> *functional, Box3D domain, NTensorField3D<T> &field, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), field, level);
}

/* *************** BoxProcessing3D_LL******************************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box3D domain,
    BlockLattice3D<T1, Descriptor1> &lattice1, BlockLattice3D<T2, Descriptor2> &lattice2)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), lattice1, lattice2);
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box3D domain,
    BlockLattice3D<T1, Descriptor1> &lattice1, BlockLattice3D<T2, Descriptor2> &lattice2,
    plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), lattice1, lattice2, level);
}

/* *************** BoxProcessing3D_SS ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain, ScalarField3D<T1> &field1,
    ScalarField3D<T2> &field2)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), field1, field2);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain, ScalarField3D<T1> &field1,
    ScalarField3D<T2> &field2, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), field1, field2, level);
}

/* *************** BoxProcessing3D_TT ****************************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    TensorField3D<T1, nDim1> &field1, TensorField3D<T2, nDim2> &field2)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), field1, field2);
}

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    TensorField3D<T1, nDim1> &field1, TensorField3D<T2, nDim2> &field2, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), field1, field2, level);
}

/* *************** BoxProcessing3D_NN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_NN<T1, T2> *functional, Box3D domain, NTensorField3D<T1> &field1,
    NTensorField3D<T2> &field2)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), field1, field2);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_NN<T1, T2> *functional, Box3D domain, NTensorField3D<T1> &field1,
    NTensorField3D<T2> &field2, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), field1, field2, level);
}

/* *************** BoxProcessing3D_ST ****************************************** */

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain, ScalarField3D<T1> &field1,
    TensorField3D<T2, nDim> &field2)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), field1, field2);
}

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain, ScalarField3D<T1> &field1,
    TensorField3D<T2, nDim> &field2, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), field1, field2, level);
}

/* *************** BoxProcessing3D_SN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_SN<T1, T2> *functional, Box3D domain, ScalarField3D<T1> &field1,
    NTensorField3D<T2> &field2)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), field1, field2);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_SN<T1, T2> *functional, Box3D domain, ScalarField3D<T1> &field1,
    NTensorField3D<T2> &field2, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), field1, field2, level);
}

/* *************** BoxProcessing3D_LS ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, Box3D domain,
    BlockLattice3D<T1, Descriptor> &lattice, ScalarField3D<T2> &field)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), lattice, field);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, Box3D domain,
    BlockLattice3D<T1, Descriptor> &lattice, ScalarField3D<T2> &field, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), lattice, field, level);
}

/* *************** BoxProcessing3D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, Box3D domain,
    BlockLattice3D<T1, Descriptor> &lattice, TensorField3D<T2, nDim> &field)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), lattice, field);
}

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, Box3D domain,
    BlockLattice3D<T1, Descriptor> &lattice, TensorField3D<T2, nDim> &field, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), lattice, field, level);
}

/* *************** BoxProcessing3D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, Box3D domain,
    BlockLattice3D<T1, Descriptor> &lattice, NTensorField3D<T2> &field)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), lattice, field);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, Box3D domain,
    BlockLattice3D<T1, Descriptor> &lattice, NTensorField3D<T2> &field, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), lattice, field, level);
}

/* *************** LatticeDotProcessing3D ****************************************** */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    LatticeDotProcessingFunctional3D<T, Descriptor> *functional, DotList3D const &dotList,
    std::vector<BlockLattice3D<T, Descriptor> *> lattices)
{
    std::vector<AtomicBlock3D *> atomicBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        atomicBlocks[iLattice] = dynamic_cast<AtomicBlock3D *>(lattices[iLattice]);
    }
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), atomicBlocks);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    LatticeDotProcessingFunctional3D<T, Descriptor> *functional, DotList3D const &dotList,
    std::vector<BlockLattice3D<T, Descriptor> *> lattices, plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        atomicBlocks[iLattice] = dynamic_cast<AtomicBlock3D *>(lattices[iLattice]);
    }
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), atomicBlocks, level);
}

/* *************** ScalarFieldDotProcessing3D ****************************************** */

template <typename T>
void applyProcessingFunctional(
    ScalarFieldDotProcessingFunctional3D<T> *functional, DotList3D const &dotList,
    std::vector<ScalarField3D<T> *> fields)
{
    std::vector<AtomicBlock3D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D *>(fields[iField]);
    }
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), atomicBlocks);
}

template <typename T>
void integrateProcessingFunctional(
    ScalarFieldDotProcessingFunctional3D<T> *functional, DotList3D const &dotList,
    std::vector<ScalarField3D<T> *> fields, plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D *>(fields[iField]);
    }
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), atomicBlocks, level);
}

/* *************** TensorFieldDotProcessing3D ****************************************** */

template <typename T, int nDim>
void applyProcessingFunctional(
    TensorFieldDotProcessingFunctional3D<T, nDim> *functional, DotList3D const &dotList,
    std::vector<TensorField3D<T, nDim> *> fields)
{
    std::vector<AtomicBlock3D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D *>(fields[iField]);
    }
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), atomicBlocks);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    TensorFieldDotProcessingFunctional3D<T, nDim> *functional, DotList3D const &dotList,
    std::vector<TensorField3D<T, nDim> *> fields, plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D *>(fields[iField]);
    }
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), atomicBlocks, level);
}

/* *************** NTensorFieldDotProcessing3D ****************************************** */

template <typename T>
void applyProcessingFunctional(
    NTensorFieldDotProcessingFunctional3D<T> *functional, DotList3D const &dotList,
    std::vector<NTensorField3D<T> *> fields)
{
    std::vector<AtomicBlock3D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D *>(fields[iField]);
    }
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), atomicBlocks);
}

template <typename T>
void integrateProcessingFunctional(
    NTensorFieldDotProcessingFunctional3D<T> *functional, DotList3D const &dotList,
    std::vector<NTensorField3D<T> *> fields, plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D *>(fields[iField]);
    }
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), atomicBlocks, level);
}

/* *************** DotProcessing3D_L ******************************************* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    DotProcessingFunctional3D_L<T, Descriptor> *functional, DotList3D const &dotList,
    BlockLattice3D<T, Descriptor> &lattice)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), lattice);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_L<T, Descriptor> *functional, DotList3D const &dotList,
    BlockLattice3D<T, Descriptor> &lattice, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), lattice, level);
}

/* *************** DotProcessing3D_S ******************************************* */

template <typename T>
void applyProcessingFunctional(
    DotProcessingFunctional3D_S<T> *functional, DotList3D const &dotList, ScalarField3D<T> &field)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), field);
}

template <typename T>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_S<T> *functional, DotList3D const &dotList, ScalarField3D<T> &field,
    plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), field, level);
}

/* *************** DotProcessing3D_T ******************************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional3D_T<T, nDim> *functional, DotList3D const &dotList,
    TensorField3D<T, nDim> &field)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), field);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_T<T, nDim> *functional, DotList3D const &dotList,
    TensorField3D<T, nDim> &field, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), field, level);
}

/* *************** DotProcessing3D_N ******************************************* */

template <typename T>
void applyProcessingFunctional(
    DotProcessingFunctional3D_N<T> *functional, DotList3D const &dotList, NTensorField3D<T> &field)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), field);
}

template <typename T>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_N<T> *functional, DotList3D const &dotList, NTensorField3D<T> &field,
    plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), field, level);
}

/* *************** DotProcessing3D_LL******************************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional,
    DotList3D const &dotList, BlockLattice3D<T1, Descriptor1> &lattice1,
    BlockLattice3D<T2, Descriptor2> &lattice2)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), lattice1, lattice2);
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional,
    DotList3D const &dotList, BlockLattice3D<T1, Descriptor1> &lattice1,
    BlockLattice3D<T2, Descriptor2> &lattice2, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), lattice1, lattice2, level);
}

/* *************** DotProcessing3D_SS ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_SS<T1, T2> *functional, DotList3D const &dotList,
    ScalarField3D<T1> &field1, ScalarField3D<T2> &field2)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), field1, field2);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_SS<T1, T2> *functional, DotList3D const &dotList,
    ScalarField3D<T1> &field1, ScalarField3D<T2> &field2, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), field1, field2, level);
}

/* *************** DotProcessing3D_TT ****************************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, DotList3D const &dotList,
    TensorField3D<T1, nDim1> &field1, TensorField3D<T2, nDim2> &field2)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), field1, field2);
}

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, DotList3D const &dotList,
    TensorField3D<T1, nDim1> &field1, TensorField3D<T2, nDim2> &field2, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), field1, field2, level);
}

/* *************** DotProcessing3D_NN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_NN<T1, T2> *functional, DotList3D const &dotList,
    NTensorField3D<T1> &field1, NTensorField3D<T2> &field2)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), field1, field2);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_NN<T1, T2> *functional, DotList3D const &dotList,
    NTensorField3D<T1> &field1, NTensorField3D<T2> &field2, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), field1, field2, level);
}

/* *************** DotProcessing3D_ST ****************************************** */

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional3D_ST<T1, T2, nDim> *functional, DotList3D const &dotList,
    ScalarField3D<T1> &field1, TensorField3D<T2, nDim> &field2)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), field1, field2);
}

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_ST<T1, T2, nDim> *functional, DotList3D const &dotList,
    ScalarField3D<T1> &field1, TensorField3D<T2, nDim> &field2, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), field1, field2, level);
}

/* *************** DotProcessing3D_SN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_SN<T1, T2> *functional, DotList3D const &dotList,
    ScalarField3D<T1> &field1, NTensorField3D<T2> &field2)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), field1, field2);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_SN<T1, T2> *functional, DotList3D const &dotList,
    ScalarField3D<T1> &field1, NTensorField3D<T2> &field2, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), field1, field2, level);
}

/* *************** DotProcessing3D_LS ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, DotList3D const &dotList,
    BlockLattice3D<T1, Descriptor> &lattice, ScalarField3D<T2> &field)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), lattice, field);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, DotList3D const &dotList,
    BlockLattice3D<T1, Descriptor> &lattice, ScalarField3D<T2> &field, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), lattice, field, level);
}

/* *************** DotProcessing3D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, DotList3D const &dotList,
    BlockLattice3D<T1, Descriptor> &lattice, TensorField3D<T2, nDim> &field)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), lattice, field);
}

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, DotList3D const &dotList,
    BlockLattice3D<T1, Descriptor> &lattice, TensorField3D<T2, nDim> &field, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), lattice, field, level);
}

/* *************** DotProcessing3D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, DotList3D const &dotList,
    BlockLattice3D<T1, Descriptor> &lattice, NTensorField3D<T2> &field)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), lattice, field);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, DotList3D const &dotList,
    BlockLattice3D<T1, Descriptor> &lattice, NTensorField3D<T2> &field, plint level)
{
    addInternalProcessor(DotProcessorGenerator3D(functional, dotList), lattice, field, level);
}

/* *************** BoundedLatticeBoxProcessing3D **************************** */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoundedLatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain,
    std::vector<BlockLattice3D<T, Descriptor> *> lattices, plint boundaryWidth)
{
    std::vector<AtomicBlock3D *> atomicBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        atomicBlocks[iLattice] = dynamic_cast<AtomicBlock3D *>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoundedLatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain,
    std::vector<BlockLattice3D<T, Descriptor> *> lattices, plint boundaryWidth, plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        atomicBlocks[iLattice] = dynamic_cast<AtomicBlock3D *>(lattices[iLattice]);
    }
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}

/* *************** BoundedScalarFieldBoxProcessing3D ************************ */

template <typename T>
void applyProcessingFunctional(
    BoundedScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<ScalarField3D<T> *> fields, plint boundaryWidth)
{
    std::vector<AtomicBlock3D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D *>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template <typename T>
void integrateProcessingFunctional(
    BoundedScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<ScalarField3D<T> *> fields, plint boundaryWidth, plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D *>(fields[iField]);
    }
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}

/* *************** BoundedTensorFieldBoxProcessing3D ************************ */

template <typename T, int nDim>
void applyProcessingFunctional(
    BoundedTensorFieldBoxProcessingFunctional3D<T, nDim> *functional, Box3D domain,
    std::vector<TensorField3D<T, nDim> *> fields, plint boundaryWidth)
{
    std::vector<AtomicBlock3D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D *>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    BoundedTensorFieldBoxProcessingFunctional3D<T, nDim> *functional, Box3D domain,
    std::vector<TensorField3D<T, nDim> *> fields, plint boundaryWidth, plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D *>(fields[iField]);
    }
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}

/* *************** BoundedNTensorFieldBoxProcessing3D ************************ */

template <typename T>
void applyProcessingFunctional(
    BoundedNTensorFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<NTensorField3D<T> *> fields, plint boundaryWidth)
{
    std::vector<AtomicBlock3D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D *>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template <typename T>
void integrateProcessingFunctional(
    BoundedNTensorFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<NTensorField3D<T> *> fields, plint boundaryWidth, plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D *>(fields[iField]);
    }
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_L ********************************* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    BlockLattice3D<T, Descriptor> &lattice, plint boundaryWidth)
{
    std::vector<AtomicBlock3D *> atomicBlocks(1);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&lattice);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    BlockLattice3D<T, Descriptor> &lattice, plint boundaryWidth, plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(1);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&lattice);
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_S ********************************* */

template <typename T>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_S<T> *functional, Box3D domain, ScalarField3D<T> &field,
    plint boundaryWidth)
{
    std::vector<AtomicBlock3D *> atomicBlocks(1);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template <typename T>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_S<T> *functional, Box3D domain, ScalarField3D<T> &field,
    plint boundaryWidth, plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(1);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&field);
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_T ********************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_T<T, nDim> *functional, Box3D domain,
    TensorField3D<T, nDim> &field, plint boundaryWidth)
{
    std::vector<AtomicBlock3D *> atomicBlocks(1);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_T<T, nDim> *functional, Box3D domain,
    TensorField3D<T, nDim> &field, plint boundaryWidth, plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(1);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&field);
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_N ********************************* */

template <typename T>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_N<T> *functional, Box3D domain, NTensorField3D<T> &field,
    plint boundaryWidth)
{
    std::vector<AtomicBlock3D *> atomicBlocks(1);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template <typename T>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_N<T> *functional, Box3D domain, NTensorField3D<T> &field,
    plint boundaryWidth, plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(1);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&field);
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_LL********************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box3D domain,
    BlockLattice3D<T1, Descriptor1> &lattice1, BlockLattice3D<T2, Descriptor2> &lattice2,
    plint boundaryWidth)
{
    std::vector<AtomicBlock3D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&lattice1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D *>(&lattice2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box3D domain,
    BlockLattice3D<T1, Descriptor1> &lattice1, BlockLattice3D<T2, Descriptor2> &lattice2,
    plint boundaryWidth, plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&lattice1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D *>(&lattice2);
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_SS ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain,
    ScalarField3D<T1> &field1, ScalarField3D<T2> &field2, plint boundaryWidth)
{
    std::vector<AtomicBlock3D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D *>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain,
    ScalarField3D<T1> &field1, ScalarField3D<T2> &field2, plint boundaryWidth, plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D *>(&field2);
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_TT ********************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    TensorField3D<T1, nDim1> &field1, TensorField3D<T2, nDim2> &field2, plint boundaryWidth)
{
    std::vector<AtomicBlock3D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D *>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    TensorField3D<T1, nDim1> &field1, TensorField3D<T2, nDim2> &field2, plint boundaryWidth,
    plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D *>(&field2);
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_NN ********************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_NN<T1, T2> *functional, Box3D domain,
    NTensorField3D<T1> &field1, NTensorField3D<T2> &field2, plint boundaryWidth)
{
    std::vector<AtomicBlock3D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D *>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_NN<T1, T2> *functional, Box3D domain,
    NTensorField3D<T1> &field1, NTensorField3D<T2> &field2, plint boundaryWidth, plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D *>(&field2);
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_ST ******************************** */

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain,
    ScalarField3D<T1> &field1, TensorField3D<T2, nDim> &field2, plint boundaryWidth)
{
    std::vector<AtomicBlock3D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D *>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain,
    ScalarField3D<T1> &field1, TensorField3D<T2, nDim> &field2, plint boundaryWidth, plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D *>(&field2);
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_SN ******************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_SN<T1, T2> *functional, Box3D domain,
    ScalarField3D<T1> &field1, NTensorField3D<T2> &field2, plint boundaryWidth)
{
    std::vector<AtomicBlock3D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D *>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_SN<T1, T2> *functional, Box3D domain,
    ScalarField3D<T1> &field1, NTensorField3D<T2> &field2, plint boundaryWidth, plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D *>(&field2);
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_LS ******************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, Box3D domain,
    BlockLattice3D<T1, Descriptor> &lattice, ScalarField3D<T2> &field, plint boundaryWidth)
{
    std::vector<AtomicBlock3D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D *>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, Box3D domain,
    BlockLattice3D<T1, Descriptor> &lattice, ScalarField3D<T2> &field, plint boundaryWidth,
    plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D *>(&field);
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, Box3D domain,
    BlockLattice3D<T1, Descriptor> &lattice, TensorField3D<T2, nDim> &field, plint boundaryWidth)
{
    std::vector<AtomicBlock3D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D *>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, Box3D domain,
    BlockLattice3D<T1, Descriptor> &lattice, TensorField3D<T2, nDim> &field, plint boundaryWidth,
    plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D *>(&field);
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, Box3D domain,
    BlockLattice3D<T1, Descriptor> &lattice, NTensorField3D<T2> &field, plint boundaryWidth)
{
    std::vector<AtomicBlock3D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D *>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LN<T1, Descriptor, T2> *functional, Box3D domain,
    BlockLattice3D<T1, Descriptor> &lattice, NTensorField3D<T2> &field, plint boundaryWidth,
    plint level)
{
    std::vector<AtomicBlock3D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D *>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D *>(&field);
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}

}  // namespace plb

#endif  // DATA_PROCESSOR_WRAPPER_3D_HH
