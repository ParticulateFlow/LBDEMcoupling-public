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

#ifndef COARSE_GRID_PROCESSORS_2D_H
#define COARSE_GRID_PROCESSORS_2D_H

#include <vector>

#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "core/globalDefs.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiGrid/gridRefinement.h"
#include "multiGrid/gridRefinementDynamics.h"

namespace plb {

/// Coupling to be added to fine lattice: copies data to coarse lattice after numTimeSteps
/// interations
template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
class CopyFineToCoarse2D : public BoxProcessingFunctional2D_LL<T, Descriptor1, T, Descriptor2> {
public:
    CopyFineToCoarse2D(
        RescaleEngine<T, Descriptor1> *rescaleEngine_, plint numTimeSteps_, plint executionTime_,
        plint direction_, plint orientation_);
    virtual ~CopyFineToCoarse2D();
    CopyFineToCoarse2D(CopyFineToCoarse2D<T, Descriptor1, Descriptor2> const &rhs);
    CopyFineToCoarse2D<T, Descriptor1, Descriptor2> &operator=(
        CopyFineToCoarse2D<T, Descriptor1, Descriptor2> const &rhs);
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor1> &fineLattice,
        BlockLattice2D<T, Descriptor2> &coarseLattice);
    virtual CopyFineToCoarse2D<T, Descriptor1, Descriptor2> *clone() const;
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
    plint direction;  // with respect to the Palabos boundary convention
    plint orientation;
};

/// Coupling to be added to fine lattice: copies data to coarse lattice after numTimeSteps
/// interations and does a filtering operation
template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
class CopyFineToCoarseWithFiltering2D :
    public BoxProcessingFunctional2D_LL<T, Descriptor1, T, Descriptor2> {
public:
    CopyFineToCoarseWithFiltering2D(
        RescaleEngine<T, Descriptor1> *rescaleEngine_, plint numTimeSteps_, plint executionTime_,
        plint direction_, plint orientation_);
    virtual ~CopyFineToCoarseWithFiltering2D();
    CopyFineToCoarseWithFiltering2D(
        CopyFineToCoarseWithFiltering2D<T, Descriptor1, Descriptor2> const &rhs);
    CopyFineToCoarseWithFiltering2D<T, Descriptor1, Descriptor2> &operator=(
        CopyFineToCoarseWithFiltering2D<T, Descriptor1, Descriptor2> const &rhs);
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor1> &fineLattice,
        BlockLattice2D<T, Descriptor2> &coarseLattice);
    virtual CopyFineToCoarseWithFiltering2D<T, Descriptor1, Descriptor2> *clone() const;
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
    plint direction;  // according to palabos BC convention
    plint orientation;
};

}  // namespace plb

#endif  // COARSE_GRID_PROCESSORS_2D_H
