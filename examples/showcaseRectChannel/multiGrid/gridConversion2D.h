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

#ifndef GRID_CONVERSION_2D_H
#define GRID_CONVERSION_2D_H

#include "atomicBlock/dataProcessingFunctional2D.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "core/dynamics.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "multiBlock/multiDataProcessorWrapper2D.h"

namespace plb {

Box2D scaleBox(Box2D box, plint nLevel);

/// Copy data from a fine to a coarse ScalarField2D.
template <typename T>
class ScalarCopyFineToCoarseZerothOrder2D : public BoxProcessingFunctional2D_SS<T, T> {
public:
    ScalarCopyFineToCoarseZerothOrder2D(plint dimDx_, plint dimDt_, plint dxScale_, plint dtScale_);
    ScalarCopyFineToCoarseZerothOrder2D(ScalarCopyFineToCoarseZerothOrder2D<T> const &rhs);
    ScalarCopyFineToCoarseZerothOrder2D<T> &operator=(
        ScalarCopyFineToCoarseZerothOrder2D<T> const &rhs);
    virtual void process(
        Box2D fineDomain, ScalarField2D<T> &fineField, ScalarField2D<T> &coarseField);
    virtual ScalarCopyFineToCoarseZerothOrder2D<T> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const
    {
        return BlockDomain::bulk;
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    }

private:
    plint dimDx, dimDt;
    plint dxScale, dtScale;
};

/// Copy data from a fine to a coarse TensorField2D.
template <typename T, int nDim>
class TensorCopyFineToCoarseZerothOrder2D : public BoxProcessingFunctional2D_TT<T, nDim, T, nDim> {
public:
    TensorCopyFineToCoarseZerothOrder2D(plint dimDx_, plint dimDt_, plint dxScale_, plint dtScale_);
    TensorCopyFineToCoarseZerothOrder2D(TensorCopyFineToCoarseZerothOrder2D<T, nDim> const &rhs);
    TensorCopyFineToCoarseZerothOrder2D<T, nDim> &operator=(
        TensorCopyFineToCoarseZerothOrder2D<T, nDim> const &rhs);
    virtual void process(
        Box2D fineDomain, TensorField2D<T, nDim> &fineField, TensorField2D<T, nDim> &coarseField);
    virtual TensorCopyFineToCoarseZerothOrder2D<T, nDim> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const
    {
        return BlockDomain::bulk;
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    }

private:
    plint dimDx, dimDt;
    plint dxScale, dtScale;
};

/// Copy data from a fine to a coarse ScalarField2D.
template <typename T, template <typename U> class Descriptor>
class LatticeCopyFineToCoarseZerothOrder2D :
    public BoxProcessingFunctional2D_LL<T, Descriptor, T, Descriptor> {
public:
    LatticeCopyFineToCoarseZerothOrder2D(plint dxScale_, plint dtScale_);
    LatticeCopyFineToCoarseZerothOrder2D(
        LatticeCopyFineToCoarseZerothOrder2D<T, Descriptor> const &rhs);
    LatticeCopyFineToCoarseZerothOrder2D<T, Descriptor> &operator=(
        LatticeCopyFineToCoarseZerothOrder2D<T, Descriptor> const &rhs);
    virtual void process(
        Box2D fineDomain, BlockLattice2D<T, Descriptor> &fineLattice,
        BlockLattice2D<T, Descriptor> &coarseLattice);
    virtual LatticeCopyFineToCoarseZerothOrder2D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const
    {
        return BlockDomain::bulk;
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    }

private:
    plint dxScale, dtScale;
};

/// Copy data from a coarse to a fine ScalarField2D, using bilinear interpolation.
template <typename T>
class ScalarInterpolateCoarseToFine2D : public BoxProcessingFunctional2D_SS<T, T> {
public:
    ScalarInterpolateCoarseToFine2D(plint dimDx_, plint dimDt_, plint dxScale_, plint dtScale_);

    ScalarInterpolateCoarseToFine2D(ScalarInterpolateCoarseToFine2D<T> const &rhs);
    ScalarInterpolateCoarseToFine2D<T> &operator=(ScalarInterpolateCoarseToFine2D<T> const &rhs);

    virtual void process(
        Box2D coarseDomain, ScalarField2D<T> &coarseField, ScalarField2D<T> &fineField);

    virtual ScalarInterpolateCoarseToFine2D<T> *clone() const;

    virtual BlockDomain::DomainT appliesTo() const
    {
        return BlockDomain::bulk;
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    }

private:
    plint dimDx, dimDt;
    plint dxScale, dtScale;
};

/// Copy data from a coarse to a fine TensorField2D, using bilinear interpolation.
template <typename T, int nDim>
class TensorInterpolateCoarseToFine2D : public BoxProcessingFunctional2D_TT<T, nDim, T, nDim> {
public:
    TensorInterpolateCoarseToFine2D(plint dimDx_, plint dimDt_, plint dxScale_, plint dtScale_);

    TensorInterpolateCoarseToFine2D(TensorInterpolateCoarseToFine2D<T, nDim> const &rhs);
    TensorInterpolateCoarseToFine2D<T, nDim> &operator=(
        TensorInterpolateCoarseToFine2D<T, nDim> const &rhs);

    virtual void process(
        Box2D coarseDomain, TensorField2D<T, nDim> &coarseField, TensorField2D<T, nDim> &fineField);

    virtual TensorInterpolateCoarseToFine2D<T, nDim> *clone() const;

    virtual BlockDomain::DomainT appliesTo() const
    {
        return BlockDomain::bulk;
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    }

private:
    plint dimDx, dimDt;
    plint dxScale, dtScale;
};

/// Copy data from a coarse to a fine BlockLattice, using bilinear interpolation.
template <typename T, template <typename U> class Descriptor>
class LatticeInterpolateCoarseToFine2D :
    public BoxProcessingFunctional2D_LL<T, Descriptor, T, Descriptor> {
public:
    LatticeInterpolateCoarseToFine2D(plint dxScale_, plint dtScale_);

    LatticeInterpolateCoarseToFine2D(LatticeInterpolateCoarseToFine2D<T, Descriptor> const &rhs);
    LatticeInterpolateCoarseToFine2D<T, Descriptor> &operator=(
        LatticeInterpolateCoarseToFine2D<T, Descriptor> const &rhs);

    virtual void process(
        Box2D coarseDomain, BlockLattice2D<T, Descriptor> &coarseLattice,
        BlockLattice2D<T, Descriptor> &fineLattice);

    virtual LatticeInterpolateCoarseToFine2D<T, Descriptor> *clone() const;

    virtual BlockDomain::DomainT appliesTo() const
    {
        return BlockDomain::bulk;
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    }

private:
    plint dxScale, dtScale;
};

/// Copy data from a coarse to a fine BlockLattice, using bilinear interpolation.
template <typename T, template <typename U> class Descriptor>
class LatticeLinearInterpolateCoarseToFine2D :
    public BoxProcessingFunctional2D_LL<T, Descriptor, T, Descriptor> {
public:
    LatticeLinearInterpolateCoarseToFine2D(plint dxScale_, plint dtScale_);

    LatticeLinearInterpolateCoarseToFine2D(
        LatticeLinearInterpolateCoarseToFine2D<T, Descriptor> const &rhs);
    LatticeLinearInterpolateCoarseToFine2D<T, Descriptor> &operator=(
        LatticeLinearInterpolateCoarseToFine2D<T, Descriptor> const &rhs);

    virtual void process(
        Box2D coarseDomain, BlockLattice2D<T, Descriptor> &coarseLattice,
        BlockLattice2D<T, Descriptor> &fineLattice);

    virtual LatticeLinearInterpolateCoarseToFine2D<T, Descriptor> *clone() const;

    virtual BlockDomain::DomainT appliesTo() const
    {
        return BlockDomain::bulk;
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    }

private:
    plint dxScale, dtScale;
};

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > coarsen(
    MultiScalarField2D<T> &fineField, plint dimDx, plint dimDt, plint dxScale, plint dtScale);

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > refine(
    MultiScalarField2D<T> &coarseField, plint dimDx, plint dimDt, plint dxScale, plint dtScale);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > coarsen(
    MultiTensorField2D<T, nDim> &fineField, plint dimDx, plint dimDt, plint dxScale, plint dtScale);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > refine(
    MultiTensorField2D<T, nDim> &coarseField, plint dimDx, plint dimDt, plint dxScale,
    plint dtScale);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > coarsen(
    MultiBlockLattice2D<T, Descriptor> &fineLattice, plint dxScale, plint dtScale,
    Dynamics<T, Descriptor> *backgroundDynamics);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > refine(
    MultiBlockLattice2D<T, Descriptor> &coarseLattice, plint dxScale, plint dtScale,
    Dynamics<T, Descriptor> *backgroundDynamics);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > linearRefine(
    MultiBlockLattice2D<T, Descriptor> &coarseLattice, plint dxScale, plint dtScale,
    Dynamics<T, Descriptor> *backgroundDynamics);

}  // namespace plb

#endif  // GRID_CONVERSION_2D_H
