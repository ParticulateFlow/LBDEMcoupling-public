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

/** \file A helper for initialising 3D boundaries -- header file.  */

#ifndef ADVECTION_DIFFUSION_BOUNDARY_CONDITION_3D_H
#define ADVECTION_DIFFUSION_BOUNDARY_CONDITION_3D_H

#include <list>
#include <vector>

#include "atomicBlock/blockLattice3D.h"
#include "complexDynamics/advectionDiffusionDynamics.h"
#include "core/dynamics.h"
#include "core/globalDefs.h"
#include "multiBlock/multiBlockLattice3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class OnLatticeAdvectionDiffusionBoundaryCondition3D {
public:
    virtual ~OnLatticeAdvectionDiffusionBoundaryCondition3D() { }

    // 3D boundary condition for temperature:
    virtual void addTemperatureBoundary0N(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureBoundary0P(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureBoundary1N(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureBoundary1P(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureBoundary2N(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureBoundary2P(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;

    virtual void addTemperatureEdge0NN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureEdge0NP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureEdge0PN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureEdge0PP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureEdge1NN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureEdge1NP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureEdge1PN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureEdge1PP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureEdge2NN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureEdge2NP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureEdge2PN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureEdge2PP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;

    virtual void addTemperatureCornerNNN(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureCornerNNP(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureCornerNPN(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureCornerNPP(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureCornerPNN(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureCornerPNP(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureCornerPPN(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureCornerPPP(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;

    void setTemperatureConditionOnBlockBoundaries(
        BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType = boundary::dirichlet);

    void setTemperatureConditionOnBlockBoundaries(
        BlockLattice3D<T, Descriptor> &lattice, Box3D applicationDomain,
        boundary::BcType bcType = boundary::dirichlet);

    void setTemperatureConditionOnBlockBoundaries(
        BlockLattice3D<T, Descriptor> &lattice, Box3D block, Box3D applicationDomain,
        boundary::BcType bcType = boundary::dirichlet);

    // 3D boundary condition for temperature:
    virtual void addTemperatureBoundary0N(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureBoundary0P(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureBoundary1N(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureBoundary1P(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureBoundary2N(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureBoundary2P(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;

    virtual void addTemperatureEdge0NN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureEdge0NP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureEdge0PN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureEdge0PP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureEdge1NN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureEdge1NP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureEdge1PN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureEdge1PP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureEdge2NN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureEdge2NP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureEdge2PN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureEdge2PP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;

    virtual void addTemperatureCornerNNN(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureCornerNNP(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureCornerNPN(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureCornerNPP(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureCornerPNN(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureCornerPNP(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureCornerPPN(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addTemperatureCornerPPP(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;

    // void setTemperatureConditionOnBlockBoundaries(MultiBlockLattice3D<T,Descriptor>& lattice,
    //            boundary::BcType bcType=boundary::dirichlet);
    void setTemperatureConditionOnBlockBoundaries(
        MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType = boundary::dirichlet);

    void setTemperatureConditionOnBlockBoundaries(
        MultiBlockLattice3D<T, Descriptor> &lattice, Box3D applicationDomain,
        boundary::BcType bcType = boundary::dirichlet);

    void setTemperatureConditionOnBlockBoundaries(
        MultiBlockLattice3D<T, Descriptor> &lattice, Box3D block, Box3D applicationDomain,
        boundary::BcType bcType = boundary::dirichlet);
};

//////  Factory function for Zou-He Thermal BC

template <typename T, template <typename U> class Descriptor>
OnLatticeAdvectionDiffusionBoundaryCondition3D<T, Descriptor>
    *createLocalAdvectionDiffusionBoundaryCondition3D();

template <typename T, template <typename U> class Descriptor>
OnLatticeAdvectionDiffusionBoundaryCondition3D<T, Descriptor>
    *createLocalRegularizedAdvectionDiffusionBoundaryCondition3D();

template <typename T, template <typename U> class Descriptor>
OnLatticeAdvectionDiffusionBoundaryCondition3D<T, Descriptor>
    *createLocalCompleteRegularizedAdvectionDiffusionBoundaryCondition3D();

}  // namespace plb

#endif
