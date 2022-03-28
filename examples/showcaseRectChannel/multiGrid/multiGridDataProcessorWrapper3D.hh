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

#ifndef MULTI_GRID_DATA_PROCESSOR_WRAPPER_3D_HH
#define MULTI_GRID_DATA_PROCESSOR_WRAPPER_3D_HH

#include "multiGrid/multiGridDataField3D.h"
#include "multiGrid/multiGridDataProcessorWrapper3D.h"
#include "multiGrid/multiGridLattice3D.h"
#include "multiGrid/multiGridOperations3D.h"

namespace plb {

/* *************** LatticeBoxProcessing3D ****************************************** */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    LatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain,
    std::vector<MultiGridLattice3D<T, Descriptor> *> lattices, plint referenceLevel)
{
    std::vector<MultiGrid3D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiGrid3D *>(lattices[iLattice]);
    }
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), multiBlocks, referenceLevel);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    LatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain,
    std::vector<MultiGridLattice3D<T, Descriptor> *> lattices, plint referenceLevel, plint level)
{
    std::vector<MultiGrid3D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiGrid3D *>(lattices[iLattice]);
    }
    addInternalProcessor(
        BoxProcessorGenerator3D(functional, domain), multiBlocks, referenceLevel, level);
}

/* *************** MultiScalarFieldBoxProcessing3D ******************************* */

template <typename T>
void applyProcessingFunctional(
    ScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiGridScalarField3D<T> *> fields, plint referenceLevel)
{
    std::vector<MultiGrid3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid3D *>(fields[iField]);
    }
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), multiBlocks, referenceLevel);
}

template <typename T>
void integrateProcessingFunctional(
    ScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiGridScalarField3D<T> *> fields, plint referenceLevel, plint level)
{
    std::vector<MultiGrid3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid3D *>(fields[iField]);
    }
    addInternalProcessor(
        BoxProcessorGenerator3D(functional, domain), multiBlocks, referenceLevel, level);
}

/* *************** TensorFieldBoxProcessing3D ******************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    TensorFieldBoxProcessingFunctional3D<T, nDim> *functional, Box3D domain,
    std::vector<MultiGridTensorField3D<T, nDim> *> fields, plint referenceLevel)
{
    std::vector<MultiGrid3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid3D *>(fields[iField]);
    }
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), multiBlocks, referenceLevel);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    TensorFieldBoxProcessingFunctional3D<T, nDim> *functional, Box3D domain,
    std::vector<MultiGridTensorField3D<T, nDim> *> fields, plint referenceLevel, plint level)
{
    std::vector<MultiGrid3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid3D *>(fields[iField]);
    }
    addInternalProcessor(
        BoxProcessorGenerator3D(functional, domain), multiBlocks, referenceLevel, level);
}

/* *************** BoxProcessing3D_L ******************************************* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    MultiGridLattice3D<T, Descriptor> &lattice, plint referenceLevel)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), lattice, referenceLevel);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    MultiGridLattice3D<T, Descriptor> &lattice, plint referenceLevel, plint level)
{
    addInternalProcessor(
        BoxProcessorGenerator3D(functional, domain), lattice, referenceLevel, level);
}

/* *************** BoxProcessing3D_S ******************************************* */

template <typename T>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_S<T> *functional, Box3D domain, MultiGridScalarField3D<T> &field,
    plint referenceLevel)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), field, referenceLevel);
}

template <typename T>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_S<T> *functional, Box3D domain, MultiGridScalarField3D<T> &field,
    plint referenceLevel, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), field, referenceLevel, level);
}

/* *************** BoxProcessing3D_T ******************************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_T<T, nDim> *functional, Box3D domain,
    MultiGridTensorField3D<T, nDim> &field, plint referenceLevel)
{
    executeDataProcessor(BoxProcessorGenerator3D(functional, domain), field, referenceLevel);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_T<T, nDim> *functional, Box3D domain,
    MultiGridTensorField3D<T, nDim> &field, plint referenceLevel, plint level)
{
    addInternalProcessor(BoxProcessorGenerator3D(functional, domain), field, referenceLevel, level);
}

/* *************** BoxProcessing3D_LL******************************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor1> &lattice1, MultiGridLattice3D<T2, Descriptor2> &lattice2,
    plint referenceLevel)
{
    executeDataProcessor(
        BoxProcessorGenerator3D(functional, domain), lattice1, lattice2, referenceLevel);
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor1> &lattice1, MultiGridLattice3D<T2, Descriptor2> &lattice2,
    plint referenceLevel, plint level)
{
    addInternalProcessor(
        BoxProcessorGenerator3D(functional, domain), lattice1, lattice2, referenceLevel, level);
}

/* *************** BoxProcessing3D_SS ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain,
    MultiGridScalarField3D<T1> &field1, MultiGridScalarField3D<T2> &field2, plint referenceLevel)
{
    executeDataProcessor(
        BoxProcessorGenerator3D(functional, domain), field1, field2, referenceLevel);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain,
    MultiGridScalarField3D<T1> &field1, MultiGridScalarField3D<T2> &field2, plint referenceLevel,
    plint level)
{
    addInternalProcessor(
        BoxProcessorGenerator3D(functional, domain), field1, field2, referenceLevel, level);
}

/* *************** BoxProcessing3D_TT ****************************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    MultiGridTensorField3D<T1, nDim1> &field1, MultiGridTensorField3D<T2, nDim2> &field2,
    plint referenceLevel)
{
    executeDataProcessor(
        BoxProcessorGenerator3D(functional, domain), field1, field2, referenceLevel);
}

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    MultiGridTensorField3D<T1, nDim1> &field1, MultiGridTensorField3D<T2, nDim2> &field2,
    plint referenceLevel, plint level)
{
    addInternalProcessor(
        BoxProcessorGenerator3D(functional, domain), field1, field2, referenceLevel, level);
}

/* *************** BoxProcessing3D_ST ****************************************** */

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain,
    MultiGridScalarField3D<T1> &field1, MultiGridTensorField3D<T2, nDim> &field2,
    plint referenceLevel)
{
    executeDataProcessor(
        BoxProcessorGenerator3D(functional, domain), field1, field2, referenceLevel);
}

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain,
    MultiGridScalarField3D<T1> &field1, MultiGridTensorField3D<T2, nDim> &field2,
    plint referenceLevel, plint level)
{
    addInternalProcessor(
        BoxProcessorGenerator3D(functional, domain), field1, field2, referenceLevel, level);
}

/* *************** BoxProcessing3D_LS ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridScalarField3D<T2> &field,
    plint referenceLevel)
{
    executeDataProcessor(
        BoxProcessorGenerator3D(functional, domain), lattice, field, referenceLevel);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridScalarField3D<T2> &field,
    plint referenceLevel, plint level)
{
    addInternalProcessor(
        BoxProcessorGenerator3D(functional, domain), lattice, field, referenceLevel, level);
}

/* *************** BoxProcessing3D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridTensorField3D<T2, nDim> &field,
    plint referenceLevel)
{
    executeDataProcessor(
        BoxProcessorGenerator3D(functional, domain), lattice, field, referenceLevel);
}

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridTensorField3D<T2, nDim> &field,
    plint referenceLevel, plint level)
{
    addInternalProcessor(
        BoxProcessorGenerator3D(functional, domain), lattice, field, referenceLevel, level);
}

/* *************** LatticeDotProcessing3D ****************************************** */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    LatticeDotProcessingFunctional3D<T, Descriptor> *functional, DotList3D const &dotList,
    std::vector<MultiGridLattice3D<T, Descriptor> *> lattices, plint referenceLevel)
{
    std::vector<MultiGrid3D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiGrid3D *>(lattices[iLattice]);
    }
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), multiBlocks, referenceLevel);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    LatticeDotProcessingFunctional3D<T, Descriptor> *functional, DotList3D const &dotList,
    std::vector<MultiGridLattice3D<T, Descriptor> *> lattices, plint referenceLevel, plint level)
{
    std::vector<MultiGrid3D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiGrid3D *>(lattices[iLattice]);
    }
    addInternalProcessor(
        DotProcessorGenerator3D(functional, dotList), multiBlocks, referenceLevel, level);
}

/* *************** ScalarFieldDotProcessing3D ****************************************** */

template <typename T>
void applyProcessingFunctional(
    ScalarFieldDotProcessingFunctional3D<T> *functional, DotList3D const &dotList,
    std::vector<MultiGridScalarField3D<T> *> fields, plint referenceLevel)
{
    std::vector<MultiGrid3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid3D *>(fields[iField]);
    }
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), multiBlocks, referenceLevel);
}

template <typename T>
void integrateProcessingFunctional(
    ScalarFieldDotProcessingFunctional3D<T> *functional, DotList3D const &dotList,
    std::vector<MultiGridScalarField3D<T> *> fields, plint referenceLevel, plint level)
{
    std::vector<MultiGrid3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid3D *>(fields[iField]);
    }
    addInternalProcessor(
        DotProcessorGenerator3D(functional, dotList), multiBlocks, referenceLevel, level);
}

/* *************** TensorFieldDotProcessing3D ****************************************** */

template <typename T, int nDim>
void applyProcessingFunctional(
    TensorFieldDotProcessingFunctional3D<T, nDim> *functional, DotList3D const &dotList,
    std::vector<MultiGridTensorField3D<T, nDim> *> fields, plint referenceLevel)
{
    std::vector<MultiGrid3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid3D *>(fields[iField]);
    }
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), multiBlocks, referenceLevel);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    TensorFieldDotProcessingFunctional3D<T, nDim> *functional, DotList3D const &dotList,
    std::vector<MultiGridTensorField3D<T, nDim> *> fields, plint referenceLevel, plint level)
{
    std::vector<MultiGrid3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid3D *>(fields[iField]);
    }
    addInternalProcessor(
        DotProcessorGenerator3D(functional, dotList), multiBlocks, referenceLevel, level);
}

/* *************** DotProcessing3D_L ******************************************* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    DotProcessingFunctional3D_L<T, Descriptor> *functional, DotList3D const &dotList,
    MultiGridLattice3D<T, Descriptor> &lattice, plint referenceLevel)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), lattice, referenceLevel);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_L<T, Descriptor> *functional, DotList3D const &dotList,
    MultiGridLattice3D<T, Descriptor> &lattice, plint referenceLevel, plint level)
{
    addInternalProcessor(
        DotProcessorGenerator3D(functional, dotList), lattice, referenceLevel, level);
}

/* *************** DotProcessing3D_S ******************************************* */

template <typename T>
void applyProcessingFunctional(
    DotProcessingFunctional3D_S<T> *functional, DotList3D const &dotList,
    MultiGridScalarField3D<T> &field, plint referenceLevel)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), field, referenceLevel);
}

template <typename T>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_S<T> *functional, DotList3D const &dotList,
    MultiGridScalarField3D<T> &field, plint referenceLevel, plint level)
{
    addInternalProcessor(
        DotProcessorGenerator3D(functional, dotList), field, referenceLevel, level);
}

/* *************** DotProcessing3D_T ******************************************* */
template <typename T, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional3D_T<T, nDim> *functional, DotList3D const &dotList,
    MultiGridTensorField3D<T, nDim> &field, plint referenceLevel)
{
    executeDataProcessor(DotProcessorGenerator3D(functional, dotList), field, referenceLevel);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_T<T, nDim> *functional, DotList3D const &dotList,
    MultiGridTensorField3D<T, nDim> &field, plint referenceLevel, plint level)
{
    addInternalProcessor(
        DotProcessorGenerator3D(functional, dotList), field, referenceLevel, level);
}

/* *************** DotProcessing3D_LL******************************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional,
    DotList3D const &dotList, MultiGridLattice3D<T1, Descriptor1> &lattice1,
    MultiGridLattice3D<T2, Descriptor2> &lattice2, plint referenceLevel)
{
    executeDataProcessor(
        DotProcessorGenerator3D(functional, dotList), lattice1, lattice2, referenceLevel);
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional,
    DotList3D const &dotList, MultiGridLattice3D<T1, Descriptor1> &lattice1,
    MultiGridLattice3D<T2, Descriptor2> &lattice2, plint referenceLevel, plint level)
{
    addInternalProcessor(
        DotProcessorGenerator3D(functional, dotList), lattice1, lattice2, referenceLevel, level);
}

/* *************** DotProcessing3D_SS ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_SS<T1, T2> *functional, DotList3D const &dotList,
    MultiGridScalarField3D<T1> &field1, MultiGridScalarField3D<T2> &field2, plint referenceLevel)
{
    executeDataProcessor(
        DotProcessorGenerator3D(functional, dotList), field1, field2, referenceLevel);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_SS<T1, T2> *functional, DotList3D const &dotList,
    MultiGridScalarField3D<T1> &field1, MultiGridScalarField3D<T2> &field2, plint referenceLevel,
    plint level)
{
    addInternalProcessor(
        DotProcessorGenerator3D(functional, dotList), field1, field2, referenceLevel, level);
}

/* *************** DotProcessing3D_TT ****************************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, DotList3D const &dotList,
    MultiGridTensorField3D<T1, nDim1> &field1, MultiGridTensorField3D<T2, nDim2> &field2,
    plint referenceLevel)
{
    executeDataProcessor(
        DotProcessorGenerator3D(functional, dotList), field1, field2, referenceLevel);
}

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, DotList3D const &dotList,
    MultiGridTensorField3D<T1, nDim1> &field1, MultiGridTensorField3D<T2, nDim2> &field2,
    plint referenceLevel, plint level)
{
    addInternalProcessor(
        DotProcessorGenerator3D(functional, dotList), field1, field2, referenceLevel, level);
}

/* *************** DotProcessing3D_ST ****************************************** */

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional3D_ST<T1, T2, nDim> *functional, DotList3D const &dotList,
    MultiGridScalarField3D<T1> &field1, MultiGridTensorField3D<T2, nDim> &field2,
    plint referenceLevel)
{
    executeDataProcessor(
        DotProcessorGenerator3D(functional, dotList), field1, field2, referenceLevel);
}

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_ST<T1, T2, nDim> *functional, DotList3D const &dotList,
    MultiGridScalarField3D<T1> &field1, MultiGridTensorField3D<T2, nDim> &field2,
    plint referenceLevel, plint level)
{
    addInternalProcessor(
        DotProcessorGenerator3D(functional, dotList), field1, field2, referenceLevel, level);
}

/* *************** DotProcessing3D_LS ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    DotProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, DotList3D const &dotList,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridScalarField3D<T2> &field,
    plint referenceLevel)
{
    executeDataProcessor(
        DotProcessorGenerator3D(functional, dotList), lattice, field, referenceLevel);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, DotList3D const &dotList,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridScalarField3D<T2> &field,
    plint referenceLevel, plint level)
{
    addInternalProcessor(
        DotProcessorGenerator3D(functional, dotList), lattice, field, referenceLevel, level);
}

/* *************** DotProcessing3D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    DotProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, DotList3D const &dotList,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridTensorField3D<T2, nDim> &field,
    plint referenceLevel)
{
    executeDataProcessor(
        DotProcessorGenerator3D(functional, dotList), lattice, field, referenceLevel);
}

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    DotProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, DotList3D const &dotList,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridTensorField3D<T2, nDim> &field,
    plint referenceLevel, plint level)
{
    addInternalProcessor(
        DotProcessorGenerator3D(functional, dotList), lattice, field, referenceLevel, level);
}

/* *************** BoundedLatticeBoxProcessing3D **************************** */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoundedLatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain,
    std::vector<MultiGridLattice3D<T, Descriptor> *> lattices, plint referenceLevel,
    plint boundaryWidth)
{
    std::vector<MultiGrid3D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiGrid3D *>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoundedLatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain,
    std::vector<MultiGridLattice3D<T, Descriptor> *> lattices, plint referenceLevel,
    plint boundaryWidth, plint level)
{
    std::vector<MultiGrid3D *> multiBlocks(lattices.size());
    for (pluint iLattice = 0; iLattice < lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiGrid3D *>(lattices[iLattice]);
    }
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

/* *************** BoundedMultiScalarFieldBoxProcessing3D ************************ */

template <typename T>
void applyProcessingFunctional(
    BoundedScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiGridScalarField3D<T> *> fields, plint referenceLevel, plint boundaryWidth)
{
    std::vector<MultiGrid3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid3D *>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <typename T>
void integrateProcessingFunctional(
    BoundedScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain,
    std::vector<MultiGridScalarField3D<T> *> fields, plint referenceLevel, plint boundaryWidth,
    plint level)
{
    std::vector<MultiGrid3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid3D *>(fields[iField]);
    }
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

/* *************** BoundedTensorFieldBoxProcessing3D ************************ */

template <typename T, int nDim>
void applyProcessingFunctional(
    BoundedTensorFieldBoxProcessingFunctional3D<T, nDim> *functional, Box3D domain,
    std::vector<MultiGridTensorField3D<T, nDim> *> fields, plint referenceLevel,
    plint boundaryWidth)
{
    std::vector<MultiGrid3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid3D *>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    BoundedTensorFieldBoxProcessingFunctional3D<T, nDim> *functional, Box3D domain,
    std::vector<MultiGridTensorField3D<T, nDim> *> fields, plint referenceLevel,
    plint boundaryWidth, plint level)
{
    std::vector<MultiGrid3D *> multiBlocks(fields.size());
    for (pluint iField = 0; iField < fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiGrid3D *>(fields[iField]);
    }
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_L ********************************* */

template <typename T, template <typename U> class Descriptor>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    MultiGridLattice3D<T, Descriptor> &lattice, plint referenceLevel, plint boundaryWidth)
{
    std::vector<MultiGrid3D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiGrid3D *>(&lattice);
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <typename T, template <typename U> class Descriptor>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain,
    MultiGridLattice3D<T, Descriptor> &lattice, plint referenceLevel, plint boundaryWidth,
    plint level)
{
    std::vector<MultiGrid3D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiGrid3D *>(&lattice);
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_S ********************************* */

template <typename T>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_S<T> *functional, Box3D domain,
    MultiGridScalarField3D<T> &field, plint referenceLevel, plint boundaryWidth)
{
    std::vector<MultiGrid3D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiGrid3D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <typename T>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_S<T> *functional, Box3D domain,
    MultiGridScalarField3D<T> &field, plint referenceLevel, plint boundaryWidth, plint level)
{
    std::vector<MultiGrid3D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiGrid3D *>(&field);
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_T ********************************* */

template <typename T, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_T<T, nDim> *functional, Box3D domain,
    MultiGridTensorField3D<T, nDim> &field, plint referenceLevel, plint boundaryWidth)
{
    std::vector<MultiGrid3D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiGrid3D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_T<T, nDim> *functional, Box3D domain,
    MultiGridTensorField3D<T, nDim> &field, plint referenceLevel, plint boundaryWidth, plint level)
{
    std::vector<MultiGrid3D *> multiBlocks(1);
    multiBlocks[0] = dynamic_cast<MultiGrid3D *>(&field);
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_LL********************************* */

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor1> &lattice1, MultiGridLattice3D<T2, Descriptor2> &lattice2,
    plint referenceLevel, plint boundaryWidth)
{
    std::vector<MultiGrid3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid3D *>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiGrid3D *>(&lattice2);
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor1> &lattice1, MultiGridLattice3D<T2, Descriptor2> &lattice2,
    plint referenceLevel, plint boundaryWidth, plint level)
{
    std::vector<MultiGrid3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid3D *>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiGrid3D *>(&lattice2);
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_SS ****************************************** */

template <typename T1, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain,
    MultiGridScalarField3D<T1> &field1, MultiGridScalarField3D<T2> &field2, plint referenceLevel,
    plint boundaryWidth)
{
    std::vector<MultiGrid3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid3D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiGrid3D *>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain,
    MultiGridScalarField3D<T1> &field1, MultiGridScalarField3D<T2> &field2, plint referenceLevel,
    plint boundaryWidth, plint level)
{
    std::vector<MultiGrid3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid3D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiGrid3D *>(&field2);
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_TT ********************************** */

template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    MultiGridTensorField3D<T1, nDim1> &field1, MultiGridTensorField3D<T2, nDim2> &field2,
    plint referenceLevel, plint boundaryWidth)
{
    std::vector<MultiGrid3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid3D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiGrid3D *>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    MultiGridTensorField3D<T1, nDim1> &field1, MultiGridTensorField3D<T2, nDim2> &field2,
    plint referenceLevel, plint boundaryWidth, plint level)
{
    std::vector<MultiGrid3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid3D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiGrid3D *>(&field2);
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_ST ******************************** */

template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain,
    MultiGridScalarField3D<T1> &field1, MultiGridTensorField3D<T2, nDim> &field2,
    plint referenceLevel, plint boundaryWidth)
{
    std::vector<MultiGrid3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid3D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiGrid3D *>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain,
    MultiGridScalarField3D<T1> &field1, MultiGridTensorField3D<T2, nDim> &field2,
    plint referenceLevel, plint boundaryWidth, plint level)
{
    std::vector<MultiGrid3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid3D *>(&field1);
    multiBlocks[1] = dynamic_cast<MultiGrid3D *>(&field2);
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_LS ******************************** */

template <typename T1, template <typename U> class Descriptor, typename T2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridScalarField3D<T2> &field,
    plint referenceLevel, plint boundaryWidth)
{
    std::vector<MultiGrid3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid3D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiGrid3D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <typename T1, template <typename U> class Descriptor, typename T2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LS<T1, Descriptor, T2> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridScalarField3D<T2> &field,
    plint referenceLevel, plint boundaryWidth, plint level)
{
    std::vector<MultiGrid3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid3D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiGrid3D *>(&field);
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

/* *************** BoundedBoxProcessing3D_LT ****************************************** */

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridTensorField3D<T2, nDim> &field,
    plint referenceLevel, plint boundaryWidth)
{
    std::vector<MultiGrid3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid3D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiGrid3D *>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, referenceLevel, boundaryWidth);
}

template <typename T1, template <typename U> class Descriptor, typename T2, int nDim>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_LT<T1, Descriptor, T2, nDim> *functional, Box3D domain,
    MultiGridLattice3D<T1, Descriptor> &lattice, MultiGridTensorField3D<T2, nDim> &field,
    plint referenceLevel, plint boundaryWidth, plint level)
{
    std::vector<MultiGrid3D *> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiGrid3D *>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiGrid3D *>(&field);
    integrateProcessingFunctional(
        functional, domain, multiBlocks, referenceLevel, boundaryWidth, level);
}

}  // namespace plb

#endif  // MULTI_GRID_DATA_PROCESSOR_WRAPPER_3D_HH
