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

#ifndef REDUCTIVE_DATA_PROCESSOR_WRAPPER_2D_HH
#define REDUCTIVE_DATA_PROCESSOR_WRAPPER_2D_HH

#include "atomicBlock/atomicBlockOperations2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include "atomicBlock/dataProcessor2D.h"
#include "atomicBlock/reductiveDataProcessorWrapper2D.h"
#include "core/plbDebug.h"

namespace plb {

/* *************** LatticeBoxProcessing2D ****************************************** */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    ReductiveLatticeBoxProcessingFunctional2D<T, Descriptor> &functional, Box2D domain,
    std::vector<BlockLattice2D<T, Descriptor> *> lattices)
{
    std::vector<AtomicBlock2D *> atomicBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        atomicBlocks[iLattice] = dynamic_cast<AtomicBlock2D *>(lattices[iLattice]);
    }
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, atomicBlocks);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ScalarFieldBoxProcessing2D ******************************* */

template <typename T>
void applyProcessingFunctional(
    ReductiveScalarFieldBoxProcessingFunctional2D<T> &functional, Box2D domain,
    std::vector<ScalarField2D<T> *> fields)
{
    std::vector<AtomicBlock2D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock2D *>(fields[iField]);
    }
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, atomicBlocks);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** TensorFieldBoxProcessing2D ******************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    ReductiveTensorFieldBoxProcessingFunctional2D<T, nDim> &functional, Box2D domain,
    std::vector<TensorField2D<T, nDim> *> fields)
{
    std::vector<AtomicBlock2D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock2D *>(fields[iField]);
    }
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, atomicBlocks);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** NTensorFieldBoxProcessing2D ******************************* */

template <typename T>
void applyProcessingFunctional(
    ReductiveNTensorFieldBoxProcessingFunctional2D<T> &functional, Box2D domain,
    std::vector<NTensorField2D<T> *> fields)
{
    std::vector<AtomicBlock2D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock2D *>(fields[iField]);
    }
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, atomicBlocks);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_L ******************************************* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_L<T, Descriptor> &functional, Box2D domain,
    BlockLattice2D<T, Descriptor> &lattice)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, lattice);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_S ******************************************* */

template <typename T>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_S<T> &functional, Box2D domain, ScalarField2D<T> &field)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_T ******************************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_T<T, nDim> &functional, Box2D domain,
    TensorField2D<T, nDim> &field)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_N ******************************************* */

template <typename T>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_N<T> &functional, Box2D domain, NTensorField2D<T> &field)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_LL******************************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> &functional,
    Box2D domain, BlockLattice2D<T1, Descriptor1> &lattice1,
    BlockLattice2D<T2, Descriptor2> &lattice2)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, lattice1, lattice2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_SS ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_SS<T1, T2> &functional, Box2D domain,
    ScalarField2D<T1> &field1, ScalarField2D<T2> &field2)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, field1, field2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_TT ****************************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> &functional, Box2D domain,
    TensorField2D<T1, nDim1> &field1, TensorField2D<T2, nDim2> &field2)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, field1, field2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_NN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_NN<T1, T2> &functional, Box2D domain,
    NTensorField2D<T1> &field1, NTensorField2D<T2> &field2)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, field1, field2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_ST ****************************************** */

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_ST<T1, T2, nDim> &functional, Box2D domain,
    ScalarField2D<T1> &field1, TensorField2D<T2, nDim> &field2)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, field1, field2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_SN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_SN<T1, T2> &functional, Box2D domain,
    ScalarField2D<T1> &field1, NTensorField2D<T2> &field2)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, field1, field2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_LS ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_LS<T1, Descriptor, T2> &functional, Box2D domain,
    BlockLattice2D<T1, Descriptor> &lattice, ScalarField2D<T2> &field)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, lattice, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> &functional, Box2D domain,
    BlockLattice2D<T1, Descriptor> &lattice, TensorField2D<T2, nDim> &field)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, lattice, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_LN<T1, Descriptor, T2> &functional, Box2D domain,
    BlockLattice2D<T1, Descriptor> &lattice, NTensorField2D<T2> &field)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, lattice, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** LatticeDotProcessing2D ****************************************** */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    ReductiveLatticeDotProcessingFunctional2D<T, Descriptor> &functional, DotList2D const &dotList,
    std::vector<BlockLattice2D<T, Descriptor> *> lattices)
{
    std::vector<AtomicBlock2D *> atomicBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        atomicBlocks[iLattice] = dynamic_cast<AtomicBlock2D *>(lattices[iLattice]);
    }
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, atomicBlocks);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ScalarFieldDotProcessing2D ****************************************** */

template <typename T>
void applyProcessingFunctional(
    ReductiveScalarFieldDotProcessingFunctional2D<T> &functional, DotList2D const &dotList,
    std::vector<ScalarField2D<T> *> fields)
{
    std::vector<AtomicBlock2D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock2D *>(fields[iField]);
    }
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, atomicBlocks);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** TensorFieldDotProcessing2D ****************************************** */

template <typename T, int nDim>
void applyProcessingFunctional(
    ReductiveTensorFieldDotProcessingFunctional2D<T, nDim> &functional, DotList2D const &dotList,
    std::vector<TensorField2D<T, nDim> *> fields)
{
    std::vector<AtomicBlock2D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock2D *>(fields[iField]);
    }
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, atomicBlocks);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** NTensorFieldDotProcessing2D ****************************************** */

template <typename T>
void applyProcessingFunctional(
    ReductiveNTensorFieldDotProcessingFunctional2D<T> &functional, DotList2D const &dotList,
    std::vector<NTensorField2D<T> *> fields)
{
    std::vector<AtomicBlock2D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock2D *>(fields[iField]);
    }
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, atomicBlocks);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_L ******************************************* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_L<T, Descriptor> &functional, DotList2D const &dotList,
    BlockLattice2D<T, Descriptor> &lattice)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, lattice);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_S ******************************************* */

template <typename T>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_S<T> &functional, DotList2D const &dotList,
    ScalarField2D<T> &field)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_T ******************************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_T<T, nDim> &functional, DotList2D const &dotList,
    TensorField2D<T, nDim> &field)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_N ******************************************* */

template <typename T>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_N<T> &functional, DotList2D const &dotList,
    NTensorField2D<T> &field)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_LL******************************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> &functional,
    DotList2D const &dotList, BlockLattice2D<T1, Descriptor1> &lattice1,
    BlockLattice2D<T2, Descriptor2> &lattice2)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, lattice1, lattice2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_SS ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_SS<T1, T2> &functional, DotList2D const &dotList,
    ScalarField2D<T1> &field1, ScalarField2D<T2> &field2)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, field1, field2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_TT ****************************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> &functional,
    DotList2D const &dotList, TensorField2D<T1, nDim1> &field1, TensorField2D<T2, nDim2> &field2)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, field1, field2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_NN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_NN<T1, T2> &functional, DotList2D const &dotList,
    NTensorField2D<T1> &field1, NTensorField2D<T2> &field2)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, field1, field2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_ST ****************************************** */

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_ST<T1, T2, nDim> &functional, DotList2D const &dotList,
    ScalarField2D<T1> &field1, TensorField2D<T2, nDim> &field2)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, field1, field2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_SN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_SN<T1, T2> &functional, DotList2D const &dotList,
    ScalarField2D<T1> &field1, NTensorField2D<T2> &field2)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, field1, field2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_LS ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_LS<T1, Descriptor, T2> &functional, DotList2D const &dotList,
    BlockLattice2D<T1, Descriptor> &lattice, ScalarField2D<T2> &field)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, lattice, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> &functional,
    DotList2D const &dotList, BlockLattice2D<T1, Descriptor> &lattice,
    TensorField2D<T2, nDim> &field)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, lattice, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_LN<T1, Descriptor, T2> &functional, DotList2D const &dotList,
    BlockLattice2D<T1, Descriptor> &lattice, NTensorField2D<T2> &field)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, lattice, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** BoundedReductiveLatticeBoxProcessing2D **************************** */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoundedReductiveLatticeBoxProcessingFunctional2D<T, Descriptor> &functional, Box2D domain,
    std::vector<BlockLattice2D<T, Descriptor> *> lattices, plint boundaryWidth)
{
    std::vector<AtomicBlock2D *> atomicBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        atomicBlocks[iLattice] = dynamic_cast<AtomicBlock2D *>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

/* *************** BoundedReductiveScalarFieldBoxProcessing2D ************************ */

template <typename T>
void applyProcessingFunctional(
    BoundedReductiveScalarFieldBoxProcessingFunctional2D<T> &functional, Box2D domain,
    std::vector<ScalarField2D<T> *> fields, plint boundaryWidth)
{
    std::vector<AtomicBlock2D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock2D *>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

/* *************** BoundedReductiveTensorFieldBoxProcessing2D ************************ */

template <typename T, int nDim>
void applyProcessingFunctional(
    BoundedReductiveTensorFieldBoxProcessingFunctional2D<T, nDim> &functional, Box2D domain,
    std::vector<TensorField2D<T, nDim> *> fields, plint boundaryWidth)
{
    std::vector<AtomicBlock2D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock2D *>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

/* *************** BoundedReductiveNTensorFieldBoxProcessing2D ************************ */

template <typename T>
void applyProcessingFunctional(
    BoundedReductiveNTensorFieldBoxProcessingFunctional2D<T> &functional, Box2D domain,
    std::vector<NTensorField2D<T> *> fields, plint boundaryWidth)
{
    std::vector<AtomicBlock2D *> atomicBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock2D *>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_L ********************************* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_L<T, Descriptor> &functional, Box2D domain,
    BlockLattice2D<T, Descriptor> &lattice, plint boundaryWidth)
{
    std::vector<AtomicBlock2D *> atomicBlocks(1);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D *>(&lattice);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_S ********************************* */

template <typename T>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_S<T> &functional, Box2D domain,
    ScalarField2D<T> &field, plint boundaryWidth)
{
    std::vector<AtomicBlock2D *> atomicBlocks(1);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D *>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_T ********************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_T<T, nDim> &functional, Box2D domain,
    TensorField2D<T, nDim> &field, plint boundaryWidth)
{
    std::vector<AtomicBlock2D *> atomicBlocks(1);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D *>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_N ********************************* */

template <typename T>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_N<T> &functional, Box2D domain,
    NTensorField2D<T> &field, plint boundaryWidth)
{
    std::vector<AtomicBlock2D *> atomicBlocks(1);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D *>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_LL********************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> &functional,
    Box2D domain, BlockLattice2D<T1, Descriptor1> &lattice1,
    BlockLattice2D<T2, Descriptor2> &lattice2, plint boundaryWidth)
{
    std::vector<AtomicBlock2D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D *>(&lattice1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D *>(&lattice2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_SS ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_SS<T1, T2> &functional, Box2D domain,
    ScalarField2D<T1> &field1, ScalarField2D<T2> &field2, plint boundaryWidth)
{
    std::vector<AtomicBlock2D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D *>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D *>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_TT ********************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> &functional, Box2D domain,
    TensorField2D<T1, nDim1> &field1, TensorField2D<T2, nDim2> &field2, plint boundaryWidth)
{
    std::vector<AtomicBlock2D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D *>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D *>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_NN ********************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_NN<T1, T2> &functional, Box2D domain,
    NTensorField2D<T1> &field1, NTensorField2D<T2> &field2, plint boundaryWidth)
{
    std::vector<AtomicBlock2D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D *>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D *>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_ST ******************************** */

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_ST<T1, T2, nDim> &functional, Box2D domain,
    ScalarField2D<T1> &field1, TensorField2D<T2, nDim> &field2, plint boundaryWidth)
{
    std::vector<AtomicBlock2D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D *>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D *>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_SN ******************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_SN<T1, T2> &functional, Box2D domain,
    ScalarField2D<T1> &field1, NTensorField2D<T2> &field2, plint boundaryWidth)
{
    std::vector<AtomicBlock2D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D *>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D *>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_LS ******************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_LS<T1, Descriptor, T2> &functional, Box2D domain,
    BlockLattice2D<T1, Descriptor> &lattice, ScalarField2D<T2> &field, plint boundaryWidth)
{
    std::vector<AtomicBlock2D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D *>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D *>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> &functional,
    Box2D domain, BlockLattice2D<T1, Descriptor> &lattice, TensorField2D<T2, nDim> &field,
    plint boundaryWidth)
{
    std::vector<AtomicBlock2D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D *>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D *>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_LN<T1, Descriptor, T2> &functional, Box2D domain,
    BlockLattice2D<T1, Descriptor> &lattice, NTensorField2D<T2> &field, plint boundaryWidth)
{
    std::vector<AtomicBlock2D *> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D *>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D *>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

}  // namespace plb

#endif  // REDUCTIVE_DATA_PROCESSOR_WRAPPER_2D_HH
