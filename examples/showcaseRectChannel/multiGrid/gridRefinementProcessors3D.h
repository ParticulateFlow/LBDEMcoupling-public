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
 * Dynamics and data processors used to implement 3D grid refinement -- header file.
 */

#ifndef GRID_REFINEMENT_PROCESSORS_3D_H
#define GRID_REFINEMENT_PROCESSORS_3D_H

#include <vector>

#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "core/globalDefs.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiGrid/gridRefinement.h"
#include "multiGrid/gridRefinementDynamics.h"

namespace plb {

/// Processor that computes the interpolations in the fine grid from values copied from
/// the coarse grid inside the bulk
template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
class InterpolateCoarseToFineDynamics3D :
    public BoxProcessingFunctional3D_LL<T, Descriptor1, T, Descriptor2> {
public:
    InterpolateCoarseToFineDynamics3D(
        RescaleEngine<T, Descriptor1> *rescaleEngine_,
        InterpolationEngine2D<T, Descriptor1> *interpolationEngine_, plint direction_,
        const Dot3D delta_[2]);
    virtual ~InterpolateCoarseToFineDynamics3D();
    InterpolateCoarseToFineDynamics3D(
        InterpolateCoarseToFineDynamics3D<T, Descriptor1, Descriptor2> const &rhs);
    InterpolateCoarseToFineDynamics3D<T, Descriptor1, Descriptor2> &operator=(
        InterpolateCoarseToFineDynamics3D<T, Descriptor1, Descriptor2> const &rhs);
    virtual void process(
        Box3D coarseDomain, BlockLattice3D<T, Descriptor1> &coarseLattice,
        BlockLattice3D<T, Descriptor2> &fineLattice);
    virtual InterpolateCoarseToFineDynamics3D<T, Descriptor1, Descriptor2> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const
    {
        return BlockDomain::bulkAndEnvelope;
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    }

private:
    RescaleEngine<T, Descriptor1> *rescaleEngine;
    InterpolationEngine2D<T, Descriptor1> *interpolationEngine;
    plint direction;
    std::vector<Dot3D> delta;
};

///
template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
class InterpolateCoarseToFineBoundaryDynamics3D :
    public BoxProcessingFunctional3D_LL<T, Descriptor1, T, Descriptor2> {
public:
    InterpolateCoarseToFineBoundaryDynamics3D(
        RescaleEngine<T, Descriptor1> *rescaleEngine_,
        InterpolationEngine2D<T, Descriptor1> *interpolationEngine_, const Dot3D delta_[2]);
    virtual ~InterpolateCoarseToFineBoundaryDynamics3D();
    InterpolateCoarseToFineBoundaryDynamics3D(
        InterpolateCoarseToFineBoundaryDynamics3D<T, Descriptor1, Descriptor2> const &rhs);
    InterpolateCoarseToFineBoundaryDynamics3D<T, Descriptor1, Descriptor2> &operator=(
        InterpolateCoarseToFineBoundaryDynamics3D<T, Descriptor1, Descriptor2> const &rhs);
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor1> &coarseLattice,
        BlockLattice3D<T, Descriptor2> &fineLattice);
    virtual InterpolateCoarseToFineBoundaryDynamics3D<T, Descriptor1, Descriptor2> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const
    {
        return BlockDomain::bulkAndEnvelope;
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    }

private:
    RescaleEngine<T, Descriptor1> *rescaleEngine;
    InterpolationEngine2D<T, Descriptor1> *interpolationEngine;
    std::vector<Dot3D> delta;
};

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
class InterpolateCoarseToFineCornerDynamics3D :
    public BoxProcessingFunctional3D_LL<T, Descriptor1, T, Descriptor2> {
public:
    InterpolateCoarseToFineCornerDynamics3D(
        RescaleEngine<T, Descriptor1> *rescaleEngine_,
        InterpolationEngine2D<T, Descriptor1> *interpolationEngine_, const Dot3D delta_[2]);
    virtual ~InterpolateCoarseToFineCornerDynamics3D();
    InterpolateCoarseToFineCornerDynamics3D(
        InterpolateCoarseToFineCornerDynamics3D<T, Descriptor1, Descriptor2> const &rhs);
    InterpolateCoarseToFineCornerDynamics3D<T, Descriptor1, Descriptor2> &operator=(
        InterpolateCoarseToFineCornerDynamics3D<T, Descriptor1, Descriptor2> const &rhs);
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor1> &coarseLattice,
        BlockLattice3D<T, Descriptor2> &fineLattice);
    virtual InterpolateCoarseToFineCornerDynamics3D<T, Descriptor1, Descriptor2> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const
    {
        return BlockDomain::bulkAndEnvelope;
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    }

private:
    RescaleEngine<T, Descriptor1> *rescaleEngine;
    InterpolationEngine2D<T, Descriptor1> *interpolationEngine;
    std::vector<Dot3D> delta;
};

/// Processor to be added to coarse lattice: copies t1 values to t0 after each numTimeSteps
/// iteration
template <typename T, template <typename U> class Descriptor>
class Copy_t1_to_t0_3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    Copy_t1_to_t0_3D(plint numTimeSteps_, plint executionTime_);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &fineLattice);
    virtual Copy_t1_to_t0_3D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const
    {
        return BlockDomain::bulkAndEnvelope;
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
    }

private:
    plint numTimeSteps;
    plint executionTime;
};

/// Coupling to be added to fine lattice: copies data to coarse lattice after numTimeSteps
/// interations
template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
class CopyFineToCoarse3D : public BoxProcessingFunctional3D_LL<T, Descriptor1, T, Descriptor2> {
public:
    CopyFineToCoarse3D(
        RescaleEngine<T, Descriptor1> *rescaleEngine_, plint numTimeSteps_, plint executionTime_);
    virtual ~CopyFineToCoarse3D();
    CopyFineToCoarse3D(CopyFineToCoarse3D<T, Descriptor1, Descriptor2> const &rhs);
    CopyFineToCoarse3D<T, Descriptor1, Descriptor2> &operator=(
        CopyFineToCoarse3D<T, Descriptor1, Descriptor2> const &rhs);
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor1> &fineLattice,
        BlockLattice3D<T, Descriptor2> &coarseLattice);
    virtual CopyFineToCoarse3D<T, Descriptor1, Descriptor2> *clone() const;
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
    RescaleEngine<T, Descriptor1> *rescaleEngine;
    plint numTimeSteps;
    plint executionTime;
};

/// A helper class for grid refinement processors, used to convert a location in
///   the coarse lattice to a location in the fine lattice.
template <typename T, template <typename U> class Descriptor>
class CoarseToFineConverter3D {
public:
    CoarseToFineConverter3D(
        BlockLattice3D<T, Descriptor> const &coarseLattice,
        BlockLattice3D<T, Descriptor> &fineLattice_) :
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
    /// Convert the z-component of a location from coarse lattice to fine lattice.
    plint fineZ(plint coarseZ) const
    {
        return (coarseZ + coarseLocation.z) * 2 - fineLocation.z;
    }

    /// Convert location from coarse to fine, add the value delta, and return a
    ///    reference to the boundary dynamics on the fine grid.
    FineGridBoundaryDynamics<T, Descriptor> &fineDynamics(Dot3D coarsePos, Dot3D delta)
    {
        return dynamic_cast<FineGridBoundaryDynamics<T, Descriptor> &>(
            fineLattice
                .get(
                    fineX(coarsePos.x) + delta.x, fineY(coarsePos.y) + delta.y,
                    fineZ(coarsePos.z) + delta.z)
                .getDynamics());
    }

    FineGridBoundaryDynamics<T, Descriptor> &fineDynamics(
        plint x, plint y, plint z, plint deltaX, plint deltaY, plint deltaZ)
    {
        return fineDynamics(Dot3D(x, y, z), Dot3D(deltaX, deltaY, deltaZ));
    }

private:
    BlockLattice3D<T, Descriptor> &fineLattice;
    Dot3D coarseLocation, fineLocation;
};

}  // namespace plb

#endif  // GRID_REFINEMENT_PROCESSORS_3D_H
