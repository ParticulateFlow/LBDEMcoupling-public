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
 * Helper functions for domain initialization -- header file.
 */
#ifndef GRID_REFINEMENT_FUNCTIONAL_3D_H
#define GRID_REFINEMENT_FUNCTIONAL_3D_H

#include "atomicBlock/dataProcessingFunctional3D.h"
#include "core/globalDefs.h"

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

template <typename T>
class TemporalInterpolationFunctional3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual TemporalInterpolationFunctional3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class CopyAndSpatialInterpolationPlaneFunctional3D : public BoxProcessingFunctional3D_NN<T, T> {
public:
    CopyAndSpatialInterpolationPlaneFunctional3D(plint direction_);

    virtual void process(Box3D domain, NTensorField3D<T> &cTensor, NTensorField3D<T> &fTensor);
    virtual CopyAndSpatialInterpolationPlaneFunctional3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    plint direction;
};

template <typename T>
class CopyAndSpatialInterpolationEdgeFunctional3D : public BoxProcessingFunctional3D_NN<T, T> {
public:
    CopyAndSpatialInterpolationEdgeFunctional3D(
        plint direction_, plint orientOne_, plint orientTwo_);

    virtual void process(Box3D domain, NTensorField3D<T> &cTensor, NTensorField3D<T> &fTensor);
    virtual CopyAndSpatialInterpolationEdgeFunctional3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    plint direction, orientOne, orientTwo;
};

template <typename T>
class CopyAndSpatialInterpolationCornerFunctional3D : public BoxProcessingFunctional3D_NN<T, T> {
public:
    CopyAndSpatialInterpolationCornerFunctional3D(
        plint direction_, plint orientOne_, plint orientTwo_);

    virtual void process(Box3D domain, NTensorField3D<T> &cTensor, NTensorField3D<T> &fTensor);
    virtual CopyAndSpatialInterpolationCornerFunctional3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    /* Direction is 0 for the interface parallel to the X
     * axis and 1 for parallel in the Y axis. Orientation is -1 for negative in the axis xComp
     * +1 is for positive in the axis xComp. **/
    plint direction, orientOne, orientTwo;
};

template <typename T>
class CopyAndSpatialInterpolationBoundaryEdgeFunctional3D :
    public BoxProcessingFunctional3D_NN<T, T> {
public:
    CopyAndSpatialInterpolationBoundaryEdgeFunctional3D(
        plint direction_, plint orientOne_, plint orientTwo_);

    virtual void process(Box3D domain, NTensorField3D<T> &cTensor, NTensorField3D<T> &fTensor);
    virtual CopyAndSpatialInterpolationBoundaryEdgeFunctional3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    /* Direction is 0 for the interface parallel to the X
     * axis and 1 for parallel in the Y axis. Orientation is -1 for negative in the axis xComp
     * +1 is for positive in the axis xComp. **/
    plint direction, orientOne, orientTwo;
};

template <typename T>
class CopyAndSpatialInterpolationBoundaryCornerFunctional3D :
    public BoxProcessingFunctional3D_NN<T, T> {
public:
    CopyAndSpatialInterpolationBoundaryCornerFunctional3D(
        plint direction_, plint orientOne_, plint orientTwo_);

    virtual void process(Box3D domain, NTensorField3D<T> &cTensor, NTensorField3D<T> &fTensor);
    virtual CopyAndSpatialInterpolationBoundaryCornerFunctional3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    /* Direction is 0 for the interface parallel to the X
     * axis, 1 for parallel in the Y axis, 2 for parallel in the Y axis. Orientation is -1 for
     * negative in the axis xComp +1 is for positive in the axis xComp. **/
    plint direction, orientOne, orientTwo;
};

template <typename T, template <typename U> class Descriptor>
class CopyAndFilterNonEquilibriumFunctional3D : public BoxProcessingFunctional3D_NN<T, T> {
public:
    CopyAndFilterNonEquilibriumFunctional3D(bool filterAll = true);
    virtual void process(Box3D domain, NTensorField3D<T> &cTensor, NTensorField3D<T> &fTensor);
    virtual CopyAndFilterNonEquilibriumFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    plint minIndexFilter;
};

template <typename T, template <typename U> class Descriptor>
class CopyAndSelectiveFilterNonEquilibriumFunctional3D : public BoxProcessingFunctional3D_NN<T, T> {
public:
    CopyAndSelectiveFilterNonEquilibriumFunctional3D(T sigma_, bool filterAll = true);
    virtual void process(Box3D domain, NTensorField3D<T> &cTensor, NTensorField3D<T> &fTensor);
    virtual CopyAndSelectiveFilterNonEquilibriumFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T sigma;
    plint minIndexFilter;
    Array<int, 3> c;
    Array<T, 3> d;
};

template <typename T, template <typename U> class Descriptor>
class CopyFunctional3D : public BoxProcessingFunctional3D_NN<T, T> {
public:
    virtual void process(Box3D domain, NTensorField3D<T> &cTensor, NTensorField3D<T> &fTensor);
    virtual CopyFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
class DecomposeAndRescaleFunctional3D : public BoxProcessingFunctional3D_LN<T, Descriptor, T> {
public:
    DecomposeAndRescaleFunctional3D(T xDt_, T xDxInv_, plint order_);

    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, NTensorField3D<T> &tensor);
    virtual DecomposeAndRescaleFunctional3D<T, Descriptor, Engine> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    T xDt, xDxInv;
    plint order;
};

template <typename T, template <typename U> class Descriptor>
class RecomposeFunctional3D : public BoxProcessingFunctional3D_LN<T, Descriptor, T> {
public:
    RecomposeFunctional3D(plint order_);

    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, NTensorField3D<T> &tensor);
    virtual RecomposeFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    plint order;
};

}  // namespace plb

#endif  // DATA_ANALYSIS_FUNCTIONAL_3D_H
