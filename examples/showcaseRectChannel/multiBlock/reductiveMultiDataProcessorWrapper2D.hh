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

#ifndef REDUCTIVE_MULTI_DATA_PROCESSOR_WRAPPER_2D_HH
#define REDUCTIVE_MULTI_DATA_PROCESSOR_WRAPPER_2D_HH

#include "atomicBlock/dataProcessor2D.h"
#include "core/plbDebug.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiBlockOperations2D.h"
#include "multiBlock/multiDataField2D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper2D.h"

namespace plb {

/* *************** LatticeBoxProcessing2D ****************************************** */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    ReductiveLatticeBoxProcessingFunctional2D<T, Descriptor> &functional, Box2D domain,
    std::vector<MultiBlockLattice2D<T, Descriptor> *> lattices)
{
    std::vector<MultiBlock2D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock2D *>(lattices[iLattice]);
    }
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, multiBlocks);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** MultiScalarFieldBoxProcessing2D ******************************* */

template <typename T>
void applyProcessingFunctional(
    ReductiveScalarFieldBoxProcessingFunctional2D<T> &functional, Box2D domain,
    std::vector<MultiScalarField2D<T> *> fields)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, multiBlocks);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** TensorFieldBoxProcessing2D ******************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    ReductiveTensorFieldBoxProcessingFunctional2D<T, nDim> &functional, Box2D domain,
    std::vector<MultiTensorField2D<T, nDim> *> fields)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, multiBlocks);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** NTensorFieldBoxProcessing2D ******************************* */

template <typename T>
void applyProcessingFunctional(
    ReductiveNTensorFieldBoxProcessingFunctional2D<T> &functional, Box2D domain,
    std::vector<MultiNTensorField2D<T> *> fields)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, multiBlocks);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_L ******************************************* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_L<T, Descriptor> &functional, Box2D domain,
    MultiBlockLattice2D<T, Descriptor> &lattice)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, lattice);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_S ******************************************* */

template <typename T>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_S<T> &functional, Box2D domain, MultiScalarField2D<T> &field)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_T ******************************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_T<T, nDim> &functional, Box2D domain,
    MultiTensorField2D<T, nDim> &field)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_N ******************************************* */

template <typename T>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_N<T> &functional, Box2D domain,
    MultiNTensorField2D<T> &field)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_N ******************************************* */

template <typename T>
void applyProcessingFunctional(
    MaskedReductiveBoxProcessingFunctional2D_N<T> &functional, Box2D domain,
    MultiNTensorField2D<T> &field, MultiNTensorField2D<int> &mask)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, field, mask);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_LL******************************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> &functional,
    Box2D domain, MultiBlockLattice2D<T1, Descriptor1> &lattice1,
    MultiBlockLattice2D<T2, Descriptor2> &lattice2)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, lattice1, lattice2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_SS ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_SS<T1, T2> &functional, Box2D domain,
    MultiScalarField2D<T1> &field1, MultiScalarField2D<T2> &field2)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, field1, field2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_TT ****************************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> &functional, Box2D domain,
    MultiTensorField2D<T1, nDim1> &field1, MultiTensorField2D<T2, nDim2> &field2)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, field1, field2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_NN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_NN<T1, T2> &functional, Box2D domain,
    MultiNTensorField2D<T1> &field1, MultiNTensorField2D<T2> &field2)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, field1, field2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_NN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    MaskedReductiveBoxProcessingFunctional2D_NN<T1, T2> &functional, Box2D domain,
    MultiNTensorField2D<T1> &field1, MultiNTensorField2D<T2> &field2,
    MultiNTensorField2D<int> &mask)
{
    std::vector<MultiBlock2D *> multiBlocks(3);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field2);
    multiBlocks[2] = dynamic_cast<MultiBlock2D *>(&mask);
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, multiBlocks);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_ST ****************************************** */

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_ST<T1, T2, nDim> &functional, Box2D domain,
    MultiScalarField2D<T1> &field1, MultiTensorField2D<T2, nDim> &field2)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, field1, field2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_SN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_SN<T1, T2> &functional, Box2D domain,
    MultiScalarField2D<T1> &field1, MultiNTensorField2D<T2> &field2)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, field1, field2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_LS ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_LS<T1, Descriptor, T2> &functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiScalarField2D<T2> &field)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, lattice, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> &functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiTensorField2D<T2, nDim> &field)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, lattice, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    ReductiveBoxProcessingFunctional2D_LN<T1, Descriptor, T2> &functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field)
{
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, lattice, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** MaskedReductiveBoxProcessing2D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    MaskedReductiveBoxProcessingFunctional2D_LN<T1, Descriptor, T2> &functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field,
    MultiNTensorField2D<int> &mask)
{
    std::vector<MultiBlock2D *> multiBlocks(3);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field);
    multiBlocks[2] = dynamic_cast<MultiBlock2D *>(&mask);
    ReductiveBoxProcessorGenerator2D generator(functional.clone(), domain);
    executeDataProcessor(generator, multiBlocks);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** LatticeDotProcessing2D ****************************************** */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    ReductiveLatticeDotProcessingFunctional2D<T, Descriptor> &functional, DotList2D const &dotList,
    std::vector<MultiBlockLattice2D<T, Descriptor> *> lattices)
{
    std::vector<MultiBlock2D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock2D *>(lattices[iLattice]);
    }
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, multiBlocks);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ScalarFieldDotProcessing2D ****************************************** */

template <typename T>
void applyProcessingFunctional(
    ReductiveScalarFieldDotProcessingFunctional2D<T> &functional, DotList2D const &dotList,
    std::vector<MultiScalarField2D<T> *> fields)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, multiBlocks);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** MultiTensorFieldDotProcessing2D ****************************************** */

template <typename T, int nDim>
void applyProcessingFunctional(
    ReductiveTensorFieldDotProcessingFunctional2D<T, nDim> &functional, DotList2D const &dotList,
    std::vector<MultiTensorField2D<T, nDim> *> fields)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, multiBlocks);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** MultiNTensorFieldDotProcessing2D ****************************************** */

template <typename T>
void applyProcessingFunctional(
    ReductiveNTensorFieldDotProcessingFunctional2D<T> &functional, DotList2D const &dotList,
    std::vector<MultiNTensorField2D<T> *> fields)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, multiBlocks);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_L ******************************************* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_L<T, Descriptor> &functional, DotList2D const &dotList,
    MultiBlockLattice2D<T, Descriptor> &lattice)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, lattice);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_S ******************************************* */

template <typename T>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_S<T> &functional, DotList2D const &dotList,
    MultiScalarField2D<T> &field)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_T ******************************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_T<T, nDim> &functional, DotList2D const &dotList,
    MultiTensorField2D<T, nDim> &field)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_N ******************************************* */

template <typename T>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_N<T> &functional, DotList2D const &dotList,
    MultiNTensorField2D<T> &field)
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
    DotList2D const &dotList, MultiBlockLattice2D<T1, Descriptor1> &lattice1,
    MultiBlockLattice2D<T2, Descriptor2> &lattice2)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, lattice1, lattice2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_SS ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_SS<T1, T2> &functional, DotList2D const &dotList,
    MultiScalarField2D<T1> &field1, MultiScalarField2D<T2> &field2)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, field1, field2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_TT ****************************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> &functional,
    DotList2D const &dotList, MultiTensorField2D<T1, nDim1> &field1,
    MultiTensorField2D<T2, nDim2> &field2)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, field1, field2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_NN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_NN<T1, T2> &functional, DotList2D const &dotList,
    MultiNTensorField2D<T1> &field1, MultiNTensorField2D<T2> &field2)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, field1, field2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_ST ****************************************** */

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_ST<T1, T2, nDim> &functional, DotList2D const &dotList,
    MultiScalarField2D<T1> &field1, MultiTensorField2D<T2, nDim> &field2)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, field1, field2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_SN ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_SN<T1, T2> &functional, DotList2D const &dotList,
    MultiScalarField2D<T1> &field1, MultiNTensorField2D<T2> &field2)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, field1, field2);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_LS ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_LS<T1, Descriptor, T2> &functional, DotList2D const &dotList,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiScalarField2D<T2> &field)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, lattice, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> &functional,
    DotList2D const &dotList, MultiBlockLattice2D<T1, Descriptor> &lattice,
    MultiTensorField2D<T2, nDim> &field)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, lattice, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    ReductiveDotProcessingFunctional2D_LN<T1, Descriptor, T2> &functional, DotList2D const &dotList,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field)
{
    ReductiveDotProcessorGenerator2D generator(functional.clone(), dotList);
    executeDataProcessor(generator, lattice, field);
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** BoundedReductiveLatticeBoxProcessing2D **************************** */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoundedReductiveLatticeBoxProcessingFunctional2D<T, Descriptor> &functional, Box2D domain,
    std::vector<MultiBlockLattice2D<T, Descriptor> *> lattices, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock2D *>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

/* *************** BoundedReductiveMultiScalarFieldBoxProcessing2D ************************ */

template <typename T>
void applyProcessingFunctional(
    BoundedReductiveScalarFieldBoxProcessingFunctional2D<T> &functional, Box2D domain,
    std::vector<MultiScalarField2D<T> *> fields, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

/* *************** BoundedReductiveTensorFieldBoxProcessing2D ************************ */

template <typename T, int nDim>
void applyProcessingFunctional(
    BoundedReductiveTensorFieldBoxProcessingFunctional2D<T, nDim> &functional, Box2D domain,
    std::vector<MultiTensorField2D<T, nDim> *> fields, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

/* *************** BoundedReductiveNTensorFieldBoxProcessing2D ************************ */

template <typename T>
void applyProcessingFunctional(
    BoundedReductiveNTensorFieldBoxProcessingFunctional2D<T> &functional, Box2D domain,
    std::vector<MultiNTensorField2D<T> *> fields, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D *>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_L ********************************* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_L<T, Descriptor> &functional, Box2D domain,
    MultiBlockLattice2D<T, Descriptor> &lattice, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&lattice);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_S ********************************* */

template <typename T>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_S<T> &functional, Box2D domain,
    MultiScalarField2D<T> &field, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_T ********************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_T<T, nDim> &functional, Box2D domain,
    MultiTensorField2D<T, nDim> &field, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_N ********************************* */

template <typename T>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_N<T> &functional, Box2D domain,
    MultiNTensorField2D<T> &field, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

/* *************** BoundedMaskedReductiveBoxProcessing2D_N ********************************* */

template <typename T>
void applyProcessingFunctional(
    BoundedMaskedReductiveBoxProcessingFunctional2D_N<T> &functional, Box2D domain,
    MultiNTensorField2D<T> &field, MultiNTensorField2D<int> &mask, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&mask);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_LL********************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> &functional,
    Box2D domain, MultiBlockLattice2D<T1, Descriptor1> &lattice1,
    MultiBlockLattice2D<T2, Descriptor2> &lattice2, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&lattice2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_SS ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_SS<T1, T2> &functional, Box2D domain,
    MultiScalarField2D<T1> &field1, MultiScalarField2D<T2> &field2, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_TT ********************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> &functional, Box2D domain,
    MultiTensorField2D<T1, nDim1> &field1, MultiTensorField2D<T2, nDim2> &field2,
    plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_NN ********************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_NN<T1, T2> &functional, Box2D domain,
    MultiNTensorField2D<T1> &field1, MultiNTensorField2D<T2> &field2, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

/* *************** BoundedMaskedReductiveBoxProcessing2D_NN ********************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedMaskedReductiveBoxProcessingFunctional2D_NN<T1, T2> &functional, Box2D domain,
    MultiNTensorField2D<T1> &field1, MultiNTensorField2D<T2> &field2,
    MultiNTensorField2D<int> &mask, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(3);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field2);
    multiBlocks[2] = dynamic_cast<MultiBlock2D *>(&mask);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_ST ******************************** */

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_ST<T1, T2, nDim> &functional, Box2D domain,
    MultiScalarField2D<T1> &field1, MultiTensorField2D<T2, nDim> &field2, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_SN ******************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_SN<T1, T2> &functional, Box2D domain,
    MultiScalarField2D<T1> &field1, MultiNTensorField2D<T2> &field2, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_LS ******************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_LS<T1, Descriptor, T2> &functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiScalarField2D<T2> &field,
    plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> &functional,
    Box2D domain, MultiBlockLattice2D<T1, Descriptor> &lattice, MultiTensorField2D<T2, nDim> &field,
    plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

/* *************** BoundedReductiveBoxProcessing2D_LN ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoundedReductiveBoxProcessingFunctional2D_LN<T1, Descriptor, T2> &functional, Box2D domain,
    MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field,
    plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

/* *************** BoundedMaskedReductiveBoxProcessing2D_LN
 * ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoundedMaskedReductiveBoxProcessingFunctional2D_LN<T1, Descriptor, T2> &functional,
    Box2D domain, MultiBlockLattice2D<T1, Descriptor> &lattice, MultiNTensorField2D<T2> &field,
    MultiNTensorField2D<int> &mask, plint boundaryWidth)
{
    std::vector<MultiBlock2D *> multiBlocks(3);
    multiBlocks[0] = dynamic_cast<MultiBlock2D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D *>(&field);
    multiBlocks[2] = dynamic_cast<MultiBlock2D *>(&mask);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

}  // namespace plb

#endif  // REDUCTIVE_MULTI_DATA_PROCESSOR_WRAPPER_2D_HH
