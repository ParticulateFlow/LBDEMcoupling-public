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
 * Coupling between grids of different refinement level -- header file.
 */

#ifndef MULTI_LEVEL_WRAPPER_3D_H
#define MULTI_LEVEL_WRAPPER_3D_H

#include <map>
#include <vector>

#include "core/globalDefs.h"
#include "gridRefinement/couplingInterfaceGenerator3D.h"
#include "gridRefinement/multiLevelNTensorField3D.h"
#include "gridRefinement/multiLevelScalarField3D.h"
#include "gridRefinement/multiLevelTensorField3D.h"

namespace plb {

// ======================================================================= //
// ==================== MULTI LEVEL LATTICE WRAPPERS ===================== //
// ======================================================================= //

// ====================== BoxProcessingFunctional3D ====================== //
inline void applyProcessingFunctional(
    BoxProcessingFunctional3D *functional, Box3D domain, plint levelOfDomain,
    std::vector<MultiLevel3D *> &levels, plint numLevels);

inline void integrateProcessingFunctional(
    BoxProcessingFunctional3D *functional, Box3D domain, plint levelOfDomain,
    std::vector<MultiLevel3D *> &levels, plint numLevels, plint levelOfIntegration = 0);

inline void integrateProcessingFunctional(
    BoxProcessingFunctional3D *functional, Box3D domain, plint levelOfDomain, MultiLevel3D &actor,
    std::vector<MultiLevel3D *> multiBlockArgs, plint numLevels, plint levelOfIntegration = 0);

// ====================== BoxProcessingFunctional3D_L ====================== //
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_L<T, Descriptor> *, Box3D domain, plint levelOfDomain,
    MultiLevelCoupling3D<T, Descriptor, Engine> &lattices);

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelCoupling3D<T, Descriptor, Engine> &lattices, plint levelOfIntegration = 0);

// ====================== BoxProcessingFunctional3D_LL ====================== //
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor3> class Engine, typename V,
    template <typename W> class Descriptor2,
    template <typename V2, template <typename W2> class Descriptor4> class Engine2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LL<T, Descriptor, V, Descriptor2> *, Box3D domain,
    plint levelOfDomain, MultiLevelCoupling3D<T, Descriptor, Engine> &lattices1,
    MultiLevelCoupling3D<V, Descriptor2, Engine2> &lattices2);

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor3> class Engine, typename V,
    template <typename W> class Descriptor2,
    template <typename V2, template <typename W2> class Descriptor4> class Engine2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LL<T, Descriptor, V, Descriptor2> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelCoupling3D<T, Descriptor, Engine> &lattices1,
    MultiLevelCoupling3D<V, Descriptor2, Engine2> &lattices2, plint levelOfIntegration = 0);

// ====================== BoxProcessingFunctional3D_LS ====================== //
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LS<T, Descriptor, T2> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelCoupling3D<T, Descriptor, Engine> &lattices,
    MultiLevelScalarField3D<T2> &scalarFields);

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LS<T, Descriptor, T2> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelCoupling3D<T, Descriptor, Engine> &lattices,
    MultiLevelScalarField3D<T2> &scalarFields, plint levelOfIntegration = 0);

// ==================== BoxProcessingFunctional3D_LS For Output ============ //
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LS<T, Descriptor, T2> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelCoupling3D<T, Descriptor, Engine> &lattices,
    MultiLevelScalarFieldForOutput3D<T2> &scalarFields);

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LS<T, Descriptor, T2> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelCoupling3D<T, Descriptor, Engine> &lattices,
    MultiLevelScalarFieldForOutput3D<T2> &scalarFields, plint levelOfIntegration = 0);

// ====================== BoxProcessingFunctional3D_LT ====================== //
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2,
    int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LT<T, Descriptor, T2, nDim> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelCoupling3D<T, Descriptor, Engine> &lattices,
    MultiLevelTensorField3D<T2, nDim> &tensorFields);

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2,
    int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LT<T, Descriptor, T2, nDim> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelCoupling3D<T, Descriptor, Engine> &lattices,
    MultiLevelTensorField3D<T2, nDim> &tensorFields, plint levelOfIntegration = 0);

// ================= BoxProcessingFunctional3D_LT for output================== //
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2,
    int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LT<T, Descriptor, T2, nDim> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelCoupling3D<T, Descriptor, Engine> &lattices,
    MultiLevelTensorFieldForOutput3D<T2, nDim> &tensorFields);

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2,
    int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LT<T, Descriptor, T2, nDim> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelCoupling3D<T, Descriptor, Engine> &lattices,
    MultiLevelTensorFieldForOutput3D<T2, nDim> &tensorFields, plint levelOfIntegration = 0);

// ====================== ScalarFieldBoxProcessingFunctional3D ====================== //
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void applyProcessingFunctional(
    LatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain, plint levelOfDomain,
    std::vector<MultiLevelCoupling3D<T, Descriptor, Engine> *> fields);

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void integrateProcessingFunctional(
    LatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain, plint levelOfDomain,
    std::vector<MultiLevelCoupling3D<T, Descriptor, Engine> *> fields,
    plint levelOfIntegration = 0);

// ======================================================================= //
// ======================== SCALAR FIELD WRAPPERS ======================== //
// ======================================================================= //

template <typename T>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_S<T> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarField3D<T> &scalarField);

template <typename T>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_S<T> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarField3D<T> &scalarField, plint levelOfIntegration);

// ====================== ForOutput ====================== //

template <typename T>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_S<T> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarFieldForOutput3D<T> &scalarField);

template <typename T>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_S<T> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarFieldForOutput3D<T> &scalarField, plint levelOfIntegration);

// ====================== BoxProcessingFunctional3D_SS ====================== //
template <typename T1, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarField3D<T1> &scalarField1, MultiLevelScalarField3D<T2> &scalarField);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarField3D<T1> &scalarField1, MultiLevelScalarField3D<T2> &scalarField2,
    plint levelOfIntegration = 0);

// ================= BoxProcessingFunctional3D_SS ForOutput ================= //
template <typename T1, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarFieldForOutput3D<T1> &scalarField1,
    MultiLevelScalarFieldForOutput3D<T2> &scalarField);

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarFieldForOutput3D<T1> &scalarField1,
    MultiLevelScalarFieldForOutput3D<T2> &scalarField2, plint levelOfIntegration = 0);

// ====================== BoxProcessingFunctional3D_ST ====================== //
template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarField3D<T1> &scalar, MultiLevelTensorField3D<T2, nDim> &tensor);

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarField3D<T1> &scalar, MultiLevelTensorField3D<T2, nDim> &tensor,
    plint levelOfIntegration = 0);

// ====================== BoxProcessingFunctional3D_ST for Output====================== //
template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarFieldForOutput3D<T1> &scalar,
    MultiLevelTensorFieldForOutput3D<T2, nDim> &tensor);

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarFieldForOutput3D<T1> &scalar,
    MultiLevelTensorFieldForOutput3D<T2, nDim> &tensor, plint levelOfIntegration = 0);

// ====================== ScalarFieldBoxProcessingFunctional3D ====================== //
template <typename T>
void applyProcessingFunctional(
    ScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain, plint levelOfDomain,
    std::vector<MultiLevelScalarField3D<T> *> fields);

template <typename T>
void integrateProcessingFunctional(
    ScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain, plint levelOfDomain,
    std::vector<MultiLevelScalarField3D<T> *> fields, plint levelOfIntegration = 0);

// ====================== ScalarFieldBoxProcessingFunctional3D For Output ====================== //
template <typename T>
void applyProcessingFunctional(
    ScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain, plint levelOfDomain,
    std::vector<MultiLevelScalarFieldForOutput3D<T> *> fields);

template <typename T>
void integrateProcessingFunctional(
    ScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain, plint levelOfDomain,
    std::vector<MultiLevelScalarFieldForOutput3D<T> *> fields, plint levelOfIntegration = 0);

// ======================================================================= //
// ======================== TENSOR FIELD WRAPPERS ======================== //
// ======================================================================= //

// ====================== BoxProcessingFunctional3D_TT ====================== //
template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelTensorField3D<T1, nDim1> &tensorFields1,
    MultiLevelTensorField3D<T2, nDim2> &tensorFields);

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelTensorField3D<T1, nDim1> &tensorFields1,
    MultiLevelTensorField3D<T2, nDim2> &tensorFields2, plint levelOfIntegration = 0);

// ============== BoxProcessingFunctional3D_TT For Output ================ //
template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelTensorFieldForOutput3D<T1, nDim1> &tensorFields1,
    MultiLevelTensorFieldForOutput3D<T2, nDim2> &tensorFields);

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelTensorFieldForOutput3D<T1, nDim1> &tensorFields1,
    MultiLevelTensorFieldForOutput3D<T2, nDim2> &tensorFields2, plint levelOfIntegration = 0);

// ====================== BoundedBoxProcessingFunctional3D_TT ====================== //
template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelTensorField3D<T1, nDim1> &tensorFields1,
    MultiLevelTensorField3D<T2, nDim2> &tensorFields, plint boundaryWidth);

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelTensorField3D<T1, nDim1> &tensorFields1,
    MultiLevelTensorField3D<T2, nDim2> &tensorFields2, plint boundaryWidth,
    plint levelOfIntegration = 0);

// =============== BoundedBoxProcessingFunctional3D_TT For Output ================ //
template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelTensorFieldForOutput3D<T1, nDim1> &tensorFields1,
    MultiLevelTensorFieldForOutput3D<T2, nDim2> &tensorFields, plint boundaryWidth);

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelTensorFieldForOutput3D<T1, nDim1> &tensorFields1,
    MultiLevelTensorFieldForOutput3D<T2, nDim2> &tensorFields2, plint boundaryWidth,
    plint levelOfIntegration = 0);

// ====================== TensorFieldBoxProcessingFunctional3D ====================== //
template <typename T, int nDim>
void applyProcessingFunctional(
    TensorFieldBoxProcessingFunctional3D<T, nDim> *functional, const Box3D &domain,
    plint levelOfDomain, std::vector<MultiLevelTensorField3D<T, nDim> *> fields);

template <typename T, int nDim>
void integrateProcessingFunctional(
    TensorFieldBoxProcessingFunctional3D<T, nDim> *functional, const Box3D &domain,
    plint levelOfDomain, std::vector<MultiLevelTensorField3D<T, nDim> *> fields,
    plint levelOfIntegration = 0);

// =============== TensorFieldBoxProcessingFunctional3D For Output ================ //
template <typename T, int nDim>
void applyProcessingFunctional(
    TensorFieldBoxProcessingFunctional3D<T, nDim> *functional, const Box3D &domain,
    plint levelOfDomain, std::vector<MultiLevelTensorFieldForOutput3D<T, nDim> *> fields);

template <typename T, int nDim>
void integrateProcessingFunctional(
    TensorFieldBoxProcessingFunctional3D<T, nDim> *functional, const Box3D &domain,
    plint levelOfDomain, std::vector<MultiLevelTensorFieldForOutput3D<T, nDim> *> fields,
    plint levelOfIntegration = 0);

// ====================== BoxProcessingFunctional3D_L ====================== //
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_L<T, Descriptor> *, Box3D domain, plint levelOfDomain,
    MultiLevelActions3D<T, Descriptor, Engine> &lattices);

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelActions3D<T, Descriptor, Engine> &lattices, plint levelOfIntegration = 0);

// ====================== BoxProcessingFunctional3D_LL ====================== //
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor3> class Engine, typename V,
    template <typename W> class Descriptor2,
    template <typename V2, template <typename W2> class Descriptor4> class Engine2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LL<T, Descriptor, V, Descriptor2> *, Box3D domain,
    plint levelOfDomain, MultiLevelActions3D<T, Descriptor, Engine> &lattices1,
    MultiLevelActions3D<V, Descriptor2, Engine2> &lattices2);

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor3> class Engine, typename V,
    template <typename W> class Descriptor2,
    template <typename V2, template <typename W2> class Descriptor4> class Engine2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LL<T, Descriptor, V, Descriptor2> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelActions3D<T, Descriptor, Engine> &lattices1,
    MultiLevelActions3D<V, Descriptor2, Engine2> &lattices2, plint levelOfIntegration = 0);

// ====================== BoxProcessingFunctional3D_LS ====================== //
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LS<T, Descriptor, T2> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelActions3D<T, Descriptor, Engine> &lattices,
    MultiLevelScalarField3D<T2> &scalarFields);

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LS<T, Descriptor, T2> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelActions3D<T, Descriptor, Engine> &lattices, MultiLevelScalarField3D<T2> &scalarFields,
    plint levelOfIntegration = 0);

// ====================== BoxProcessingFunctional3D_LT ====================== //
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2,
    int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LT<T, Descriptor, T2, nDim> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelActions3D<T, Descriptor, Engine> &lattices,
    MultiLevelTensorField3D<T2, nDim> &tensorFields);

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2,
    int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LT<T, Descriptor, T2, nDim> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelActions3D<T, Descriptor, Engine> &lattices,
    MultiLevelTensorField3D<T2, nDim> &tensorFields, plint levelOfIntegration = 0);

// ====================== ScalarFieldBoxProcessingFunctional3D ====================== //
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void applyProcessingFunctional(
    LatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain, plint levelOfDomain,
    std::vector<MultiLevelActions3D<T, Descriptor, Engine> *> fields);

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void integrateProcessingFunctional(
    LatticeBoxProcessingFunctional3D<T, Descriptor> *functional, Box3D domain, plint levelOfDomain,
    std::vector<MultiLevelActions3D<T, Descriptor, Engine> *> fields, plint levelOfIntegration = 0);

}  // namespace plb

#endif  // MULTI_LEVEL_WRAPPER_3D_H
