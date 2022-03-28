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

#ifndef GRID_CONVERSION_3D_H
#define GRID_CONVERSION_3D_H

#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "core/dynamics.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"
#include "particles/multiParticleField3D.h"

namespace plb {

Box3D scaleBox(Box3D box, plint nLevel);

/// Copy data from a fine to a coarse ScalarField3D.
template <typename T>
class ScalarCopyFineToCoarseZerothOrder3D : public BoxProcessingFunctional3D_SS<T, T> {
public:
    ScalarCopyFineToCoarseZerothOrder3D(plint dimDx_, plint dimDt_, plint dxScale_, plint dtScale_);
    ScalarCopyFineToCoarseZerothOrder3D(ScalarCopyFineToCoarseZerothOrder3D<T> const &rhs);
    ScalarCopyFineToCoarseZerothOrder3D<T> &operator=(
        ScalarCopyFineToCoarseZerothOrder3D<T> const &rhs);
    virtual void process(
        Box3D fineDomain, ScalarField3D<T> &fineField, ScalarField3D<T> &coarseField);
    virtual ScalarCopyFineToCoarseZerothOrder3D<T> *clone() const;
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

/// Copy data from a fine to a coarse TensorField3D.
template <typename T, int nDim>
class TensorCopyFineToCoarseZerothOrder3D : public BoxProcessingFunctional3D_TT<T, nDim, T, nDim> {
public:
    TensorCopyFineToCoarseZerothOrder3D(plint dimDx_, plint dimDt_, plint dxScale_, plint dtScale_);
    TensorCopyFineToCoarseZerothOrder3D(TensorCopyFineToCoarseZerothOrder3D<T, nDim> const &rhs);
    TensorCopyFineToCoarseZerothOrder3D<T, nDim> &operator=(
        TensorCopyFineToCoarseZerothOrder3D<T, nDim> const &rhs);
    virtual void process(
        Box3D fineDomain, TensorField3D<T, nDim> &fineField, TensorField3D<T, nDim> &coarseField);
    virtual TensorCopyFineToCoarseZerothOrder3D<T, nDim> *clone() const;
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

/// Copy data from a fine to a coarse ScalarField3D.
template <typename T, template <typename U> class Descriptor>
class LatticeCopyFineToCoarseZerothOrder3D :
    public BoxProcessingFunctional3D_LL<T, Descriptor, T, Descriptor> {
public:
    LatticeCopyFineToCoarseZerothOrder3D(plint dxScale_, plint dtScale_);
    LatticeCopyFineToCoarseZerothOrder3D(
        LatticeCopyFineToCoarseZerothOrder3D<T, Descriptor> const &rhs);
    LatticeCopyFineToCoarseZerothOrder3D<T, Descriptor> &operator=(
        LatticeCopyFineToCoarseZerothOrder3D<T, Descriptor> const &rhs);
    virtual void process(
        Box3D fineDomain, BlockLattice3D<T, Descriptor> &fineLattice,
        BlockLattice3D<T, Descriptor> &coarseLattice);
    virtual LatticeCopyFineToCoarseZerothOrder3D<T, Descriptor> *clone() const;
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

/// Copy data from a fine to a coarse ScalarField3D.
template <typename T, template <typename U> class Descriptor>
class LatticeCopyFineToCoarseNoScale3D :
    public BoxProcessingFunctional3D_LL<T, Descriptor, T, Descriptor> {
public:
    LatticeCopyFineToCoarseNoScale3D(plint dxScale_);
    LatticeCopyFineToCoarseNoScale3D(LatticeCopyFineToCoarseNoScale3D<T, Descriptor> const &rhs);
    LatticeCopyFineToCoarseNoScale3D<T, Descriptor> &operator=(
        LatticeCopyFineToCoarseNoScale3D<T, Descriptor> const &rhs);
    virtual void process(
        Box3D fineDomain, BlockLattice3D<T, Descriptor> &fineLattice,
        BlockLattice3D<T, Descriptor> &coarseLattice);
    virtual LatticeCopyFineToCoarseNoScale3D<T, Descriptor> *clone() const;
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
    plint dxScale;
};

/// Copy particles from a fine to a coarse field.
template <typename T, template <typename U> class Descriptor>
class CopyRescaleParticles3D : public BoxProcessingFunctional3D {
public:
    CopyRescaleParticles3D(plint dxScale_, plint dtScale_);
    virtual void processGenericBlocks(Box3D fineDomain, std::vector<AtomicBlock3D *> blocks);
    virtual CopyRescaleParticles3D<T, Descriptor> *clone() const;
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

/// Copy data from a coarse to a fine ScalarField3D, using bilinear interpolation.
template <typename T>
class ScalarInterpolateCoarseToFine3D : public BoxProcessingFunctional3D_SS<T, T> {
public:
    ScalarInterpolateCoarseToFine3D(plint dimDx_, plint dimDt_, plint dxScale_, plint dtScale_);

    ScalarInterpolateCoarseToFine3D(ScalarInterpolateCoarseToFine3D<T> const &rhs);
    ScalarInterpolateCoarseToFine3D<T> &operator=(ScalarInterpolateCoarseToFine3D<T> const &rhs);

    virtual void process(
        Box3D coarseDomain, ScalarField3D<T> &coarseField, ScalarField3D<T> &fineField);

    virtual ScalarInterpolateCoarseToFine3D<T> *clone() const;

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

/// Copy data from a coarse to a fine TensorField3D, using bilinear interpolation.
template <typename T, int nDim>
class TensorInterpolateCoarseToFine3D : public BoxProcessingFunctional3D_TT<T, nDim, T, nDim> {
public:
    TensorInterpolateCoarseToFine3D(plint dimDx_, plint dimDt_, plint dxScale_, plint dtScale_);

    TensorInterpolateCoarseToFine3D(TensorInterpolateCoarseToFine3D<T, nDim> const &rhs);
    TensorInterpolateCoarseToFine3D<T, nDim> &operator=(
        TensorInterpolateCoarseToFine3D<T, nDim> const &rhs);

    virtual void process(
        Box3D coarseDomain, TensorField3D<T, nDim> &coarseField, TensorField3D<T, nDim> &fineField);

    virtual TensorInterpolateCoarseToFine3D<T, nDim> *clone() const;

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
class LatticeInterpolateCoarseToFine3D :
    public BoxProcessingFunctional3D_LL<T, Descriptor, T, Descriptor> {
public:
    LatticeInterpolateCoarseToFine3D(plint dxScale_, plint dtScale_);

    LatticeInterpolateCoarseToFine3D(LatticeInterpolateCoarseToFine3D<T, Descriptor> const &rhs);
    LatticeInterpolateCoarseToFine3D<T, Descriptor> &operator=(
        LatticeInterpolateCoarseToFine3D<T, Descriptor> const &rhs);

    virtual void process(
        Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
        BlockLattice3D<T, Descriptor> &fineLattice);

    virtual LatticeInterpolateCoarseToFine3D<T, Descriptor> *clone() const;

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
class LatticeLinearInterpolateCoarseToFine3D :
    public BoxProcessingFunctional3D_LL<T, Descriptor, T, Descriptor> {
public:
    LatticeLinearInterpolateCoarseToFine3D(plint dxScale_, plint dtScale_);

    LatticeLinearInterpolateCoarseToFine3D(
        LatticeLinearInterpolateCoarseToFine3D<T, Descriptor> const &rhs);
    LatticeLinearInterpolateCoarseToFine3D<T, Descriptor> &operator=(
        LatticeLinearInterpolateCoarseToFine3D<T, Descriptor> const &rhs);

    virtual void process(
        Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
        BlockLattice3D<T, Descriptor> &fineLattice);

    virtual LatticeLinearInterpolateCoarseToFine3D<T, Descriptor> *clone() const;

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
std::unique_ptr<MultiScalarField3D<T> > coarsen(
    MultiScalarField3D<T> &fineField, plint dimDx, plint dimDt, plint dxScale, plint dtScale);

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > refine(
    MultiScalarField3D<T> &coarseField, plint dimDx, plint dimDt, plint dxScale, plint dtScale);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > coarsen(
    MultiTensorField3D<T, nDim> &fineField, plint dimDx, plint dimDt, plint dxScale, plint dtScale);

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > refine(
    MultiTensorField3D<T, nDim> &coarseField, plint dimDx, plint dimDt, plint dxScale,
    plint dtScale);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > coarsen(
    MultiBlockLattice3D<T, Descriptor> &fineLattice, plint dxScale, plint dtScale,
    Dynamics<T, Descriptor> *backgroundDynamics);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > coarsenNoScale(
    MultiBlockLattice3D<T, Descriptor> &fineLattice, plint dxScale,
    Dynamics<T, Descriptor> *backgroundDynamics);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > refine(
    MultiBlockLattice3D<T, Descriptor> &coarseLattice, plint dxScale, plint dtScale,
    Dynamics<T, Descriptor> *backgroundDynamics);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiParticleField3D<DenseParticleField3D<T, Descriptor> > > refine(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &coarseField, plint dxScale,
    plint dtScale);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiParticleField3D<DenseParticleField3D<T, Descriptor> > > coarsen(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &fineField, plint dxScale,
    plint dtScale);

}  // namespace plb

#endif  // GRID_CONVERSION_3D_H
