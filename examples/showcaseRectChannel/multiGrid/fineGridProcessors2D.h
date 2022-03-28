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
 * Dynamics and data processors used to implement 2D grid refinement -- header file.
 */

#ifndef FINE_GRID_PROCESSORS_2D_H
#define FINE_GRID_PROCESSORS_2D_H

#include <vector>

#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "core/globalDefs.h"
#include "finiteDifference/fdStencils1D.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiGrid/gridRefinement.h"
#include "multiGrid/gridRefinementDynamics.h"

namespace plb {

/// Coupling to be added to coarse lattice: interpolates on coarse and copies
///   to in-between fine nodes.
/** This processor interpolates to the left and right of each coarse node contained in the
 *  indicated domain. There is an error check to avoid out-of-range errors when accessing
 *  coarse and fine cells, both of which could happen, because the processor acts on bulk
 *  and boundary nodes. In an out-of-range case, interpolation is not executed.
 */
template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
class CopyCoarseToFineLinearInterp2D :
    public BoxProcessingFunctional2D_LL<T, Descriptor1, T, Descriptor2> {
public:
    CopyCoarseToFineLinearInterp2D(
        RescaleEngine<T, Descriptor1> *rescaleEngine_, plint direction_, plint orientation_);
    virtual ~CopyCoarseToFineLinearInterp2D();
    CopyCoarseToFineLinearInterp2D(
        CopyCoarseToFineLinearInterp2D<T, Descriptor1, Descriptor2> const &rhs);
    CopyCoarseToFineLinearInterp2D<T, Descriptor1, Descriptor2> &operator=(
        CopyCoarseToFineLinearInterp2D<T, Descriptor1, Descriptor2> const &rhs);
    virtual void process(
        Box2D coarseDomain, BlockLattice2D<T, Descriptor1> &coarseLattice,
        BlockLattice2D<T, Descriptor2> &fineLattice);
    virtual CopyCoarseToFineLinearInterp2D<T, Descriptor1, Descriptor2> *clone() const;

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::allVariables;
    }

private:
    RescaleEngine<T, Descriptor1> *rescaleEngine;
    plint direction;    // with respect to BC convention
    plint orientation;  // with respect to BC convention
};

/// Coupling to be added to coarse lattice: interpolates on the boundary of a grid-refined domain.
/** This processor interpolates in the desired direction from the chosen coarse node.
 *  There is an error check to avoid out-of-range errors when accessing coarse and fine cells,
 *  both of which could happen, because the processor acts on bulk and boundary nodes. In an
 *  out-of-range case, interpolation is not executed.
 */
template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
class CopyCoarseToFineBoundaryLinearInterp2D :
    public BoxProcessingFunctional2D_LL<T, Descriptor1, T, Descriptor2> {
public:
    CopyCoarseToFineBoundaryLinearInterp2D(
        RescaleEngine<T, Descriptor1> *rescaleEngine_, plint direction_, plint orientation_,
        plint whereToInterpolate_);
    virtual ~CopyCoarseToFineBoundaryLinearInterp2D();
    CopyCoarseToFineBoundaryLinearInterp2D(
        CopyCoarseToFineBoundaryLinearInterp2D<T, Descriptor1, Descriptor2> const &rhs);
    CopyCoarseToFineBoundaryLinearInterp2D<T, Descriptor1, Descriptor2> &operator=(
        CopyCoarseToFineBoundaryLinearInterp2D<T, Descriptor1, Descriptor2> const &rhs);
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor1> &coarseLattice,
        BlockLattice2D<T, Descriptor2> &fineLattice);
    virtual CopyCoarseToFineBoundaryLinearInterp2D<T, Descriptor1, Descriptor2> *clone() const;

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::allVariables;
    }

private:
    RescaleEngine<T, Descriptor1> *rescaleEngine;
    plint direction, orientation, whereToInterpolate;
};

/// Copy from coarse to fine grid and perform a third order interpolation for unknown populations
template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
class CopyCoarseToFineCubicInterp2D :
    public BoxProcessingFunctional2D_LL<T, Descriptor1, T, Descriptor2> {
public:
    CopyCoarseToFineCubicInterp2D(
        RescaleEngine<T, Descriptor1> *rescaleEngine_, plint direction_, plint orientation_);
    virtual ~CopyCoarseToFineCubicInterp2D();
    CopyCoarseToFineCubicInterp2D(
        CopyCoarseToFineCubicInterp2D<T, Descriptor1, Descriptor2> const &rhs);
    CopyCoarseToFineCubicInterp2D<T, Descriptor1, Descriptor2> &operator=(
        CopyCoarseToFineCubicInterp2D<T, Descriptor1, Descriptor2> const &rhs);
    virtual void process(
        Box2D coarseDomain, BlockLattice2D<T, Descriptor1> &coarseLattice,
        BlockLattice2D<T, Descriptor2> &fineLattice);
    virtual CopyCoarseToFineCubicInterp2D<T, Descriptor1, Descriptor2> *clone() const;

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::allVariables;
    }

private:
    RescaleEngine<T, Descriptor1> *rescaleEngine;
    plint direction;
    plint orientation;
};

/// Copy from coarse to fine grid and perform a third order interpolation for unknown populations
template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
class CopyCoarseToFineCubicInterpWithDeconvolution2D :
    public BoxProcessingFunctional2D_LL<T, Descriptor1, T, Descriptor2> {
public:
    CopyCoarseToFineCubicInterpWithDeconvolution2D(
        RescaleEngine<T, Descriptor1> *rescaleEngine_, plint direction_, plint orientation_);
    virtual ~CopyCoarseToFineCubicInterpWithDeconvolution2D();
    CopyCoarseToFineCubicInterpWithDeconvolution2D(
        CopyCoarseToFineCubicInterpWithDeconvolution2D<T, Descriptor1, Descriptor2> const &rhs);
    CopyCoarseToFineCubicInterpWithDeconvolution2D<T, Descriptor1, Descriptor2> &operator=(
        CopyCoarseToFineCubicInterpWithDeconvolution2D<T, Descriptor1, Descriptor2> const &rhs);
    virtual void process(
        Box2D coarseDomain, BlockLattice2D<T, Descriptor1> &coarseLattice,
        BlockLattice2D<T, Descriptor2> &fineLattice);
    virtual CopyCoarseToFineCubicInterpWithDeconvolution2D<T, Descriptor1, Descriptor2> *clone()
        const;

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::allVariables;
    }

private:
    std::vector<T> computeDeconvolutionExtrapolation(
        BlockLattice2D<T, Descriptor2> &fineLattice, plint iX, plint iY, plint dx, plint dy) const;

private:
    RescaleEngine<T, Descriptor1> *rescaleEngine;
    plint direction;    // 0 for x, 1 for y
    plint orientation;  // either +1 or -1
};

/// Copy over the boundary from coarse to fine with third order interpolation for unknown
/// populations
template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
class CopyCoarseToFineBoundaryCubicInterp2D :
    public BoxProcessingFunctional2D_LL<T, Descriptor1, T, Descriptor2> {
public:
    CopyCoarseToFineBoundaryCubicInterp2D(
        RescaleEngine<T, Descriptor1> *rescaleEngine_, plint direction_, plint orientation_,
        plint whereToInterpolate_);
    virtual ~CopyCoarseToFineBoundaryCubicInterp2D();
    CopyCoarseToFineBoundaryCubicInterp2D(
        CopyCoarseToFineBoundaryCubicInterp2D<T, Descriptor1, Descriptor2> const &rhs);
    CopyCoarseToFineBoundaryCubicInterp2D<T, Descriptor1, Descriptor2> &operator=(
        CopyCoarseToFineBoundaryCubicInterp2D<T, Descriptor1, Descriptor2> const &rhs);
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor1> &coarseLattice,
        BlockLattice2D<T, Descriptor2> &fineLattice);
    virtual CopyCoarseToFineBoundaryCubicInterp2D<T, Descriptor1, Descriptor2> *clone() const;

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::allVariables;
    }

private:
    RescaleEngine<T, Descriptor1> *rescaleEngine;
    plint direction;
    plint orientation;
    plint whereToInterpolate;
};

/// When performing a cubic interpolation there is one site that might not have the correct values
/// near the border
template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
class CopyCoarseToFineBoundaryHelper2D :
    public BoxProcessingFunctional2D_LL<T, Descriptor1, T, Descriptor2> {
public:
    CopyCoarseToFineBoundaryHelper2D(
        RescaleEngine<T, Descriptor1> *rescaleEngine_, plint direction_, plint orientation_,
        plint whereToInterpolate_);
    virtual ~CopyCoarseToFineBoundaryHelper2D();
    CopyCoarseToFineBoundaryHelper2D(
        CopyCoarseToFineBoundaryHelper2D<T, Descriptor1, Descriptor2> const &rhs);
    CopyCoarseToFineBoundaryHelper2D<T, Descriptor1, Descriptor2> &operator=(
        CopyCoarseToFineBoundaryHelper2D<T, Descriptor1, Descriptor2> const &rhs);
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor1> &coarseLattice,
        BlockLattice2D<T, Descriptor2> &fineLattice);
    virtual CopyCoarseToFineBoundaryHelper2D<T, Descriptor1, Descriptor2> *clone() const;

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::allVariables;
    }

private:
    RescaleEngine<T, Descriptor1> *rescaleEngine;
    plint direction;
    plint orientation;
    plint whereToInterpolate;
};

/// Processor to be added to coarse lattice: copies t1 values to t0 after each numTimeSteps
/// iteration
template <typename T, template <typename U> class Descriptor>
class Copy_t1_to_t0_2D : public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    Copy_t1_to_t0_2D(plint numTimeSteps_, plint executionTime_);
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &fineLattice);
    virtual Copy_t1_to_t0_2D<T, Descriptor> *clone() const;

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::allVariables;
    }

private:
    plint numTimeSteps;
    plint executionTime;
};

/// A helper class for grid refinement processors, used to convert a location in
///   the coarse lattice to a location in the fine lattice.
template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
class CoarseToFineConverter2D {
public:
    CoarseToFineConverter2D(
        BlockLattice2D<T, Descriptor1> const &coarseLattice,
        BlockLattice2D<T, Descriptor2> &fineLattice_) :
        fineLattice(fineLattice_),
        coarseLocation(coarseLattice.getLocation()),
        fineLocation(fineLattice.getLocation())
    { }

    /// Convert the x-component of a location from coarse lattice to fine lattice.
    plint fineX(plint coarseX) const
    {
        return (coarseX + coarseLocation.x) * 2 - fineLocation.x;
    }

    /// Convert the y-component of a location from coarse lattice to fine lattice.
    plint fineY(plint coarseY) const
    {
        return (coarseY + coarseLocation.y) * 2 - fineLocation.y;
    }

    /// Convert location from coarse to fine, add the value delta, and return a
    ///    reference to the boundary dynamics on the fine grid.
    FineGridBoundaryDynamics<T, Descriptor2> &fineDynamics(
        plint coarseX, plint coarseY, plint deltaX, plint deltaY)
    {
        return dynamic_cast<FineGridBoundaryDynamics<T, Descriptor2> &>(
            fineLattice.get(fineX(coarseX) + deltaX, fineY(coarseY) + deltaY).getDynamics());
    }

private:
    BlockLattice2D<T, Descriptor2> &fineLattice;
    Dot2D coarseLocation, fineLocation;
};

}  // namespace plb

#endif  // FINE_GRID_PROCESSORS_2D_H
