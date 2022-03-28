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

/** \file A helper for initialising 2D boundaries -- header file.  */
#ifndef ADVECTION_DIFFUSION_BOUNDARY_CONDITION_2D_H
#define ADVECTION_DIFFUSION_BOUNDARY_CONDITION_2D_H

#include <list>
#include <vector>

#include "atomicBlock/blockLattice2D.h"
#include "complexDynamics/advectionDiffusionDynamics.h"
#include "core/globalDefs.h"
#include "multiBlock/multiBlockLattice2D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class OnLatticeAdvectionDiffusionBoundaryCondition2D {
public:
    virtual ~OnLatticeAdvectionDiffusionBoundaryCondition2D() { }

    virtual void addTemperatureBoundary0N(Box2D domain, BlockLattice2D<T, Descriptor> &lattice) = 0;
    virtual void addTemperatureBoundary0P(Box2D domain, BlockLattice2D<T, Descriptor> &lattice) = 0;
    virtual void addTemperatureBoundary1N(Box2D domain, BlockLattice2D<T, Descriptor> &lattice) = 0;
    virtual void addTemperatureBoundary1P(Box2D domain, BlockLattice2D<T, Descriptor> &lattice) = 0;

    virtual void addTemperatureCornerNN(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice) = 0;
    virtual void addTemperatureCornerNP(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice) = 0;
    virtual void addTemperatureCornerPN(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice) = 0;
    virtual void addTemperatureCornerPP(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice) = 0;

    void setTemperatureConditionOnBlockBoundaries(BlockLattice2D<T, Descriptor> &lattice);

    virtual void addTemperatureBoundary0N(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice) = 0;
    virtual void addTemperatureBoundary0P(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice) = 0;
    virtual void addTemperatureBoundary1N(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice) = 0;
    virtual void addTemperatureBoundary1P(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice) = 0;

    virtual void addTemperatureCornerNN(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice) = 0;
    virtual void addTemperatureCornerNP(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice) = 0;
    virtual void addTemperatureCornerPN(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice) = 0;
    virtual void addTemperatureCornerPP(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice) = 0;

    void setTemperatureConditionOnBlockBoundaries(MultiBlockLattice2D<T, Descriptor> &lattice);
};

//////  Factory function for Zou-He Thermal BC

template <typename T, template <typename U> class Descriptor>
OnLatticeAdvectionDiffusionBoundaryCondition2D<T, Descriptor>
    *createLocalAdvectionDiffusionBoundaryCondition2D();

//////  Factory function for Regularized Thermal BC

template <typename T, template <typename U> class Descriptor>
OnLatticeAdvectionDiffusionBoundaryCondition2D<T, Descriptor>
    *createLocalRegularizedAdvectionDiffusionBoundaryCondition2D();

//////  Factory function for Complete Regularized Thermal BC

template <typename T, template <typename U> class Descriptor>
OnLatticeAdvectionDiffusionBoundaryCondition2D<T, Descriptor>
    *createLocalCompleteRegularizedAdvectionDiffusionBoundaryCondition2D();

}  // namespace plb

#endif
