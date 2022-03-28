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

#ifndef MULTI_GRID_DATA_PROCESSOR_WRAPPER_2D_HH
#define MULTI_GRID_DATA_PROCESSOR_WRAPPER_2D_HH

#include "multiGrid/multiGridDataField2D.h"
#include "multiGrid/multiGridDataProcessorWrapper2D.h"
#include "multiGrid/multiGridLattice2D.h"
#include "multiGrid/multiGridOperations2D.h"

namespace plb {

/* *************** LatticeBoxProcessing2D ****************************************** */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    LatticeBoxProcessingFunctional2D<T, Descriptor> *functional, Box2D domain,
    std::vector<MultiGridLattice2D<T, Descriptor> *> lattices, plint referenceLevel)
{
    std::vector<MultiGrid2D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiGrid2D *>(lattices[iLattice]);
    }
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), multiBlocks, referenceLevel);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    LatticeBoxProcessingFunctional2D<T, Descriptor> *functional, Box2D domain,
    std::vector<MultiGridLattice2D<T, Descriptor> *> lattices, plint referenceLevel, plint level)
{
    std::vector<MultiGrid2D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiGrid2D *>(lattices[iLattice]);
    }
    addInternalProcessor(
        BoxProcessorGenerator2D(functional, domain), multiBlocks, referenceLevel, level);
}

/* *************** MultiScalarFieldBoxProcessing2D ******************************* */

template <typename T>
void applyProcessingFunctional(
    ScalarFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiGridScalarField2D<T> *> fields, plint referenceLevel)
{
    std::vector<MultiGrid2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid2D *>(fields[iField]);
    }
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), multiBlocks, referenceLevel);
}

template <typename T>
void integrateProcessingFunctional(
    ScalarFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiGridScalarField2D<T> *> fields, plint referenceLevel, plint level)
{
    std::vector<MultiGrid2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid2D *>(fields[iField]);
    }
    addInternalProcessor(
        BoxProcessorGenerator2D(functional, domain), multiBlocks, referenceLevel, level);
}

/* *************** TensorFieldBoxProcessing2D ******************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    TensorFieldBoxProcessingFunctional2D<T, nDim> *functional, Box2D domain,
    std::vector<MultiGridTensorField2D<T, nDim> *> fields, plint referenceLevel)
{
    std::vector<MultiGrid2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid2D *>(fields[iField]);
    }
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), multiBlocks, referenceLevel);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    TensorFieldBoxProcessingFunctional2D<T, nDim> *functional, Box2D domain,
    std::vector<MultiGridTensorField2D<T, nDim> *> fields, plint referenceLevel, plint level)
{
    std::vector<MultiGrid2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid2D *>(fields[iField]);
    }
    addInternalProcessor(
        BoxProcessorGenerator2D(functional, domain), multiBlocks, referenceLevel, level);
}

/* *************** BoxProcessing2D_L ******************************************* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_L<T, Descriptor> *functional, Box2D domain,
    MultiGridLattice2D<T, Descriptor> &lattice, plint referenceLevel)
{
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), lattice, referenceLevel);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_L<T, Descriptor> *functional, Box2D domain,
    MultiGridLattice2D<T, Descriptor> &lattice, plint referenceLevel, plint level)
{
    addInternalProcessor(
        BoxProcessorGenerator2D(functional, domain), lattice, referenceLevel, level);
}

/* *************** BoxProcessing2D_S ******************************************* */

template <typename T>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_S<T> *functional, Box2D domain, MultiGridScalarField2D<T> &field,
    plint referenceLevel)
{
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), field, referenceLevel);
}

template <typename T>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_S<T> *functional, Box2D domain, MultiGridScalarField2D<T> &field,
    plint referenceLevel, plint level)
{
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), field, referenceLevel, level);
}

/* *************** BoxProcessing2D_T ******************************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_T<T, nDim> *functional, Box2D domain,
    MultiGridTensorField2D<T, nDim> &field, plint referenceLevel)
{
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), field, referenceLevel);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_T<T, nDim> *functional, Box2D domain,
    MultiGridTensorField2D<T, nDim> &field, plint referenceLevel, plint level)
{
    addInternalProcessor(BoxProcessorGenerator2D(functional, domain), field, referenceLevel, level);
}

/* *************** BoxProcessing2D_LL******************************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box2D domain,
    MultiGridLattice2D<T1, Descriptor1> &lattice1, MultiGridLattice2D<T2, Descriptor2> &lattice2,
    plint referenceLevel)
{
    executeDataProcessor(
        BoxProcessorGenerator2D(functional, domain), lattice1, lattice2, referenceLevel);
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box2D domain,
    MultiGridLattice2D<T1, Descriptor1> &lattice1, MultiGridLattice2D<T2, Descriptor2> &lattice2,
    plint referenceLevel, plint level)
{
    addInternalProcessor(
        BoxProcessorGenerator2D(functional, domain), lattice1, lattice2, referenceLevel, level);
}

/* *************** BoxProcessing2D_SS ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_SS<T1, T2> *functional, Box2D domain,
    MultiGridScalarField2D<T1> &field1, MultiGridScalarField2D<T2> &field2, plint referenceLevel)
{
    executeDataProcessor(
        BoxProcessorGenerator2D(functional, domain), field1, field2, referenceLevel);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_SS<T1, T2> *functional, Box2D domain,
    MultiGridScalarField2D<T1> &field1, MultiGridScalarField2D<T2> &field2, plint referenceLevel,
    plint level)
{
    addInternalProcessor(
        BoxProcessorGenerator2D(functional, domain), field1, field2, referenceLevel, level);
}

/* *************** BoxProcessing2D_TT ****************************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> *functional, Box2D domain,
    MultiGridTensorField2D<T1, nDim1> &field1, MultiGridTensorField2D<T2, nDim2> &field2,
    plint referenceLevel)
{
    executeDataProcessor(
        BoxProcessorGenerator2D(functional, domain), field1, field2, referenceLevel);
}

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> *functional, Box2D domain,
    MultiGridTensorField2D<T1, nDim1> &field1, MultiGridTensorField2D<T2, nDim2> &field2,
    plint referenceLevel, plint level)
{
    addInternalProcessor(
        BoxProcessorGenerator2D(functional, domain), field1, field2, referenceLevel, level);
}

/* *************** BoxProcessing2D_ST ****************************************** */

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_ST<T1, T2, nDim> *functional, Box2D domain,
    MultiGridScalarField2D<T1> &field1, MultiGridTensorField2D<T2, nDim> &field2,
    plint referenceLevel)
{
    executeDataProcessor(
        BoxProcessorGenerator2D(functional, domain), field1, field2, referenceLevel);
}

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_ST<T1, T2, nDim> *functional, Box2D domain,
    MultiGridScalarField2D<T1> &field1, MultiGridTensorField2D<T2, nDim> &field2,
    plint referenceLevel, plint level)
{
    addInternalProcessor(
        BoxProcessorGenerator2D(functional, domain), field1, field2, referenceLevel, level);
}

/* *************** BoxProcessing2D_LS ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_LS<T1, Descriptor, T2> *functional, Box2D domain,
    MultiGridLattice2D<T1, Descriptor> &lattice, MultiGridScalarField2D<T2> &field,
    plint referenceLevel)
{
    executeDataProcessor(
        BoxProcessorGenerator2D(functional, domain), lattice, field, referenceLevel);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_LS<T1, Descriptor, T2> *functional, Box2D domain,
    MultiGridLattice2D<T1, Descriptor> &lattice, MultiGridScalarField2D<T2> &field,
    plint referenceLevel, plint level)
{
    addInternalProcessor(
        BoxProcessorGenerator2D(functional, domain), lattice, field, referenceLevel, level);
}

/* *************** BoxProcessing2D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> *functional, Box2D domain,
    MultiGridLattice2D<T1, Descriptor> &lattice, MultiGridTensorField2D<T2, nDim> &field,
    plint referenceLevel)
{
    executeDataProcessor(
        BoxProcessorGenerator2D(functional, domain), lattice, field, referenceLevel);
}

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> *functional, Box2D domain,
    MultiGridLattice2D<T1, Descriptor> &lattice, MultiGridTensorField2D<T2, nDim> &field,
    plint referenceLevel, plint level)
{
    addInternalProcessor(
        BoxProcessorGenerator2D(functional, domain), lattice, field, referenceLevel, level);
}

/* *************** LatticeDotProcessing2D ****************************************** */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    LatticeDotProcessingFunctional2D<T, Descriptor> *functional, DotList2D const &dotList,
    std::vector<MultiGridLattice2D<T, Descriptor> *> lattices, plint referenceLevel)
{
    std::vector<MultiGrid2D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiGrid2D *>(lattices[iLattice]);
    }
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), multiBlocks, referenceLevel);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    LatticeDotProcessingFunctional2D<T, Descriptor> *functional, DotList2D const &dotList,
    std::vector<MultiGridLattice2D<T, Descriptor> *> lattices, plint referenceLevel, plint level)
{
    std::vector<MultiGrid2D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiGrid2D *>(lattices[iLattice]);
    }
    addInternalProcessor(
        DotProcessorGenerator2D(functional, dotList), multiBlocks, referenceLevel, level);
}

/* *************** ScalarFieldDotProcessing2D ****************************************** */

template <typename T>
void applyProcessingFunctional(
    ScalarFieldDotProcessingFunctional2D<T> *functional, DotList2D const &dotList,
    std::vector<MultiGridScalarField2D<T> *> fields, plint referenceLevel)
{
    std::vector<MultiGrid2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid2D *>(fields[iField]);
    }
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), multiBlocks, referenceLevel);
}

template <typename T>
void integrateProcessingFunctional(
    ScalarFieldDotProcessingFunctional2D<T> *functional, DotList2D const &dotList,
    std::vector<MultiGridScalarField2D<T> *> fields, plint referenceLevel, plint level)
{
    std::vector<MultiGrid2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid2D *>(fields[iField]);
    }
    addInternalProcessor(
        DotProcessorGenerator2D(functional, dotList), multiBlocks, referenceLevel, level);
}

/* *************** TensorFieldDotProcessing2D ****************************************** */

template <typename T, int nDim>
void applyProcessingFunctional(
    TensorFieldDotProcessingFunctional2D<T, nDim> *functional, DotList2D const &dotList,
    std::vector<MultiGridTensorField2D<T, nDim> *> fields, plint referenceLevel)
{
    std::vector<MultiGrid2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid2D *>(fields[iField]);
    }
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), multiBlocks, referenceLevel);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    TensorFieldDotProcessingFunctional2D<T, nDim> *functional, DotList2D const &dotList,
    std::vector<MultiGridTensorField2D<T, nDim> *> fields, plint referenceLevel, plint level)
{
    std::vector<MultiGrid2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid2D *>(fields[iField]);
    }
    addInternalProcessor(
        DotProcessorGenerator2D(functional, dotList), multiBlocks, referenceLevel, level);
}

/* *************** DotProcessing2D_L ******************************************* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    DotProcessingFunctional2D_L<T, Descriptor> *functional, DotList2D const &dotList,
    MultiGridLattice2D<T, Descriptor> &lattice, plint referenceLevel)
{
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), lattice, referenceLevel);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_L<T, Descriptor> *functional, DotList2D const &dotList,
    MultiGridLattice2D<T, Descriptor> &lattice, plint referenceLevel, plint level)
{
    addInternalProcessor(
        DotProcessorGenerator2D(functional, dotList), lattice, referenceLevel, level);
}

/* *************** DotProcessing2D_S ******************************************* */

template <typename T>
void applyProcessingFunctional(
    DotProcessingFunctional2D_S<T> *functional, DotList2D const &dotList,
    MultiGridScalarField2D<T> &field, plint referenceLevel)
{
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), field, referenceLevel);
}

template <typename T>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_S<T> *functional, DotList2D const &dotList,
    MultiGridScalarField2D<T> &field, plint referenceLevel, plint level)
{
    addInternalProcessor(
        DotProcessorGenerator2D(functional, dotList), field, referenceLevel, level);
}

/* *************** DotProcessing2D_T ******************************************* */
template <typename T, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional2D_T<T, nDim> *functional, DotList2D const &dotList,
    MultiGridTensorField2D<T, nDim> &field, plint referenceLevel)
{
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), field, referenceLevel);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_T<T, nDim> *functional, DotList2D const &dotList,
    MultiGridTensorField2D<T, nDim> &field, plint referenceLevel, plint level)
{
    addInternalProcessor(
        DotProcessorGenerator2D(functional, dotList), field, referenceLevel, level);
}

/* *************** DotProcessing2D_LL******************************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    DotProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> *functional,
    DotList2D const &dotList, MultiGridLattice2D<T1, Descriptor1> &lattice1,
    MultiGridLattice2D<T2, Descriptor2> &lattice2, plint referenceLevel)
{
    executeDataProcessor(
        DotProcessorGenerator2D(functional, dotList), lattice1, lattice2, referenceLevel);
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> *functional,
    DotList2D const &dotList, MultiGridLattice2D<T1, Descriptor1> &lattice1,
    MultiGridLattice2D<T2, Descriptor2> &lattice2, plint referenceLevel, plint level)
{
    addInternalProcessor(
        DotProcessorGenerator2D(functional, dotList), lattice1, lattice2, referenceLevel, level);
}

/* *************** DotProcessing2D_SS ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional2D_SS<T1, T2> *functional, DotList2D const &dotList,
    MultiGridScalarField2D<T1> &field1, MultiGridScalarField2D<T2> &field2, plint referenceLevel)
{
    executeDataProcessor(
        DotProcessorGenerator2D(functional, dotList), field1, field2, referenceLevel);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_SS<T1, T2> *functional, DotList2D const &dotList,
    MultiGridScalarField2D<T1> &field1, MultiGridScalarField2D<T2> &field2, plint referenceLevel,
    plint level)
{
    addInternalProcessor(
        DotProcessorGenerator2D(functional, dotList), field1, field2, referenceLevel, level);
}

/* *************** DotProcessing2D_TT ****************************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    DotProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> *functional, DotList2D const &dotList,
    MultiGridTensorField2D<T1, nDim1> &field1, MultiGridTensorField2D<T2, nDim2> &field2,
    plint referenceLevel)
{
    executeDataProcessor(
        DotProcessorGenerator2D(functional, dotList), field1, field2, referenceLevel);
}

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> *functional, DotList2D const &dotList,
    MultiGridTensorField2D<T1, nDim1> &field1, MultiGridTensorField2D<T2, nDim2> &field2,
    plint referenceLevel, plint level)
{
    addInternalProcessor(
        DotProcessorGenerator2D(functional, dotList), field1, field2, referenceLevel, level);
}

/* *************** DotProcessing2D_ST ****************************************** */

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional2D_ST<T1, T2, nDim> *functional, DotList2D const &dotList,
    MultiGridScalarField2D<T1> &field1, MultiGridTensorField2D<T2, nDim> &field2,
    plint referenceLevel)
{
    executeDataProcessor(
        DotProcessorGenerator2D(functional, dotList), field1, field2, referenceLevel);
}

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_ST<T1, T2, nDim> *functional, DotList2D const &dotList,
    MultiGridScalarField2D<T1> &field1, MultiGridTensorField2D<T2, nDim> &field2,
    plint referenceLevel, plint level)
{
    addInternalProcessor(
        DotProcessorGenerator2D(functional, dotList), field1, field2, referenceLevel, level);
}

/* *************** DotProcessing2D_LS ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional2D_LS<T1, Descriptor, T2> *functional, DotList2D const &dotList,
    MultiGridLattice2D<T1, Descriptor> &lattice, MultiGridScalarField2D<T2> &field,
    plint referenceLevel)
{
    executeDataProcessor(
        DotProcessorGenerator2D(functional, dotList), lattice, field, referenceLevel);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_LS<T1, Descriptor, T2> *functional, DotList2D const &dotList,
    MultiGridLattice2D<T1, Descriptor> &lattice, MultiGridScalarField2D<T2> &field,
    plint referenceLevel, plint level)
{
    addInternalProcessor(
        DotProcessorGenerator2D(functional, dotList), lattice, field, referenceLevel, level);
}

/* *************** DotProcessing2D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> *functional, DotList2D const &dotList,
    MultiGridLattice2D<T1, Descriptor> &lattice, MultiGridTensorField2D<T2, nDim> &field,
    plint referenceLevel)
{
    executeDataProcessor(
        DotProcessorGenerator2D(functional, dotList), lattice, field, referenceLevel);
}

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> *functional, DotList2D const &dotList,
    MultiGridLattice2D<T1, Descriptor> &lattice, MultiGridTensorField2D<T2, nDim> &field,
    plint referenceLevel, plint level)
{
    addInternalProcessor(
        DotProcessorGenerator2D(functional, dotList), lattice, field, referenceLevel, level);
}

/* *************** BoundedLatticeBoxProcessing2D **************************** */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoundedLatticeBoxProcessingFunctional2D<T, Descriptor> *functional, Box2D domain,
    std::vector<MultiGridLattice2D<T, Descriptor> *> lattices, plint referenceLevel,
    plint boundaryWidth)
{
    std::vector<MultiGrid2D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiGrid2D *>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoundedLatticeBoxProcessingFunctional2D<T, Descriptor> *functional, Box2D domain,
    std::vector<MultiGridLattice2D<T, Descriptor> *> lattices, plint referenceLevel,
    plint boundaryWidth, plint level)
{
    std::vector<MultiGrid2D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiGrid2D *>(lattices[iLattice]);
    }
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

/* *************** BoundedMultiScalarFieldBoxProcessing2D ************************ */

template <typename T>
void applyProcessingFunctional(
    BoundedScalarFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiGridScalarField2D<T> *> fields, plint referenceLevel, plint boundaryWidth)
{
    std::vector<MultiGrid2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid2D *>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <typename T>
void integrateProcessingFunctional(
    BoundedScalarFieldBoxProcessingFunctional2D<T> *functional, Box2D domain,
    std::vector<MultiGridScalarField2D<T> *> fields, plint referenceLevel, plint boundaryWidth,
    plint level)
{
    std::vector<MultiGrid2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid2D *>(fields[iField]);
    }
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

/* *************** BoundedTensorFieldBoxProcessing2D ************************ */

template <typename T, int nDim>
void applyProcessingFunctional(
    BoundedTensorFieldBoxProcessingFunctional2D<T, nDim> *functional, Box2D domain,
    std::vector<MultiGridTensorField2D<T, nDim> *> fields, plint referenceLevel,
    plint boundaryWidth)
{
    std::vector<MultiGrid2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid2D *>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    BoundedTensorFieldBoxProcessingFunctional2D<T, nDim> *functional, Box2D domain,
    std::vector<MultiGridTensorField2D<T, nDim> *> fields, plint referenceLevel,
    plint boundaryWidth, plint level)
{
    std::vector<MultiGrid2D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid2D *>(fields[iField]);
    }
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_L ********************************* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_L<T, Descriptor> *functional, Box2D domain,
    MultiGridLattice2D<T, Descriptor> &lattice, plint referenceLevel, plint boundaryWidth)
{
    std::vector<MultiGrid2D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiGrid2D *>(&lattice);
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_L<T, Descriptor> *functional, Box2D domain,
    MultiGridLattice2D<T, Descriptor> &lattice, plint referenceLevel, plint boundaryWidth,
    plint level)
{
    std::vector<MultiGrid2D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiGrid2D *>(&lattice);
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_S ********************************* */

template <typename T>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_S<T> *functional, Box2D domain,
    MultiGridScalarField2D<T> &field, plint referenceLevel, plint boundaryWidth)
{
    std::vector<MultiGrid2D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiGrid2D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <typename T>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_S<T> *functional, Box2D domain,
    MultiGridScalarField2D<T> &field, plint referenceLevel, plint boundaryWidth, plint level)
{
    std::vector<MultiGrid2D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiGrid2D *>(&field);
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_T ********************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_T<T, nDim> *functional, Box2D domain,
    MultiGridTensorField2D<T, nDim> &field, plint referenceLevel, plint boundaryWidth)
{
    std::vector<MultiGrid2D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiGrid2D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_T<T, nDim> *functional, Box2D domain,
    MultiGridTensorField2D<T, nDim> &field, plint referenceLevel, plint boundaryWidth, plint level)
{
    std::vector<MultiGrid2D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiGrid2D *>(&field);
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_LL********************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box2D domain,
    MultiGridLattice2D<T1, Descriptor1> &lattice1, MultiGridLattice2D<T2, Descriptor2> &lattice2,
    plint referenceLevel, plint boundaryWidth)
{
    std::vector<MultiGrid2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid2D *>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiGrid2D *>(&lattice2);
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box2D domain,
    MultiGridLattice2D<T1, Descriptor1> &lattice1, MultiGridLattice2D<T2, Descriptor2> &lattice2,
    plint referenceLevel, plint boundaryWidth, plint level)
{
    std::vector<MultiGrid2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid2D *>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiGrid2D *>(&lattice2);
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_SS ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_SS<T1, T2> *functional, Box2D domain,
    MultiGridScalarField2D<T1> &field1, MultiGridScalarField2D<T2> &field2, plint referenceLevel,
    plint boundaryWidth)
{
    std::vector<MultiGrid2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiGrid2D *>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_SS<T1, T2> *functional, Box2D domain,
    MultiGridScalarField2D<T1> &field1, MultiGridScalarField2D<T2> &field2, plint referenceLevel,
    plint boundaryWidth, plint level)
{
    std::vector<MultiGrid2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiGrid2D *>(&field2);
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_TT ********************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> *functional, Box2D domain,
    MultiGridTensorField2D<T1, nDim1> &field1, MultiGridTensorField2D<T2, nDim2> &field2,
    plint referenceLevel, plint boundaryWidth)
{
    std::vector<MultiGrid2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiGrid2D *>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_TT<T1, nDim1, T2, nDim2> *functional, Box2D domain,
    MultiGridTensorField2D<T1, nDim1> &field1, MultiGridTensorField2D<T2, nDim2> &field2,
    plint referenceLevel, plint boundaryWidth, plint level)
{
    std::vector<MultiGrid2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiGrid2D *>(&field2);
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_ST ******************************** */

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_ST<T1, T2, nDim> *functional, Box2D domain,
    MultiGridScalarField2D<T1> &field1, MultiGridTensorField2D<T2, nDim> &field2,
    plint referenceLevel, plint boundaryWidth)
{
    std::vector<MultiGrid2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiGrid2D *>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_ST<T1, T2, nDim> *functional, Box2D domain,
    MultiGridScalarField2D<T1> &field1, MultiGridTensorField2D<T2, nDim> &field2,
    plint referenceLevel, plint boundaryWidth, plint level)
{
    std::vector<MultiGrid2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid2D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiGrid2D *>(&field2);
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_LS ******************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_LS<T1, Descriptor, T2> *functional, Box2D domain,
    MultiGridLattice2D<T1, Descriptor> &lattice, MultiGridScalarField2D<T2> &field,
    plint referenceLevel, plint boundaryWidth)
{
    std::vector<MultiGrid2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid2D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiGrid2D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_LS<T1, Descriptor, T2> *functional, Box2D domain,
    MultiGridLattice2D<T1, Descriptor> &lattice, MultiGridScalarField2D<T2> &field,
    plint referenceLevel, plint boundaryWidth, plint level)
{
    std::vector<MultiGrid2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid2D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiGrid2D *>(&field);
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing2D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> *functional, Box2D domain,
    MultiGridLattice2D<T1, Descriptor> &lattice, MultiGridTensorField2D<T2, nDim> &field,
    plint referenceLevel, plint boundaryWidth)
{
    std::vector<MultiGrid2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid2D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiGrid2D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D_LT<T1, Descriptor, T2, nDim> *functional, Box2D domain,
    MultiGridLattice2D<T1, Descriptor> &lattice, MultiGridTensorField2D<T2, nDim> &field,
    plint referenceLevel, plint boundaryWidth, plint level)
{
    std::vector<MultiGrid2D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid2D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiGrid2D *>(&field);
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

}  // namespace plb

#endif  // MULTI_GRID_DATA_PROCESSOR_WRAPPER_2D_HH
