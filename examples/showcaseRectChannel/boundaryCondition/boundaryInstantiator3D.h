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
#ifndef BOUNDARY_INSTANTIATOR_3D_H
#define BOUNDARY_INSTANTIATOR_3D_H

#include "atomicBlock/atomicBlockOperations3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "boundaryCondition/boundaryCondition3D.h"
#include "boundaryCondition/neumannCondition3D.h"
#include "core/cell.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "dataProcessors/dataInitializerFunctional3D.h"
#include "dataProcessors/dataInitializerWrapper3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiBlockOperations3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
class BoundaryConditionInstantiator3D : public OnLatticeBoundaryCondition3D<T, Descriptor> {
public:
    virtual BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager> *clone() const;

    // PART I: Atomic-block version.

    void addVelocityBoundary0N(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addVelocityBoundary0P(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addVelocityBoundary1N(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addVelocityBoundary1P(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addVelocityBoundary2N(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addVelocityBoundary2P(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addPressureBoundary0N(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addPressureBoundary0P(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addPressureBoundary1N(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addPressureBoundary1P(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addPressureBoundary2N(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addPressureBoundary2P(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addExternalVelocityEdge0NN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityEdge0NP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityEdge0PN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityEdge0PP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityEdge1NN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityEdge1NP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityEdge1PN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityEdge1PP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityEdge2NN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityEdge2NP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityEdge2PN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityEdge2PP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addInternalVelocityEdge0NN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityEdge0NP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityEdge0PN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityEdge0PP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityEdge1NN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityEdge1NP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityEdge1PN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityEdge1PP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityEdge2NN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityEdge2NP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityEdge2PN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityEdge2PP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addExternalVelocityCornerNNN(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityCornerNNP(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityCornerNPN(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityCornerNPP(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityCornerPNN(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityCornerPNP(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityCornerPPN(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityCornerPPP(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addInternalVelocityCornerNNN(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityCornerNNP(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityCornerNPN(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityCornerNPP(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityCornerPNN(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityCornerPNP(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityCornerPPN(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityCornerPPP(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);

    // PART II: Multi-block version.

    void addVelocityBoundary0N(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addVelocityBoundary0P(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addVelocityBoundary1N(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addVelocityBoundary1P(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addVelocityBoundary2N(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addVelocityBoundary2P(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addPressureBoundary0N(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addPressureBoundary0P(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addPressureBoundary1N(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addPressureBoundary1P(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addPressureBoundary2N(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addPressureBoundary2P(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addExternalVelocityEdge0NN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityEdge0NP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityEdge0PN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityEdge0PP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityEdge1NN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityEdge1NP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityEdge1PN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityEdge1PP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityEdge2NN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityEdge2NP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityEdge2PN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityEdge2PP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addInternalVelocityEdge0NN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityEdge0NP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityEdge0PN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityEdge0PP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityEdge1NN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityEdge1NP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityEdge1PN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityEdge1PP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityEdge2NN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityEdge2NP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityEdge2PN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityEdge2PP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addExternalVelocityCornerNNN(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
    void addExternalVelocityCornerNNP(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
    void addExternalVelocityCornerNPN(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
    void addExternalVelocityCornerNPP(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
    void addExternalVelocityCornerPNN(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
    void addExternalVelocityCornerPNP(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
    void addExternalVelocityCornerPPN(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
    void addExternalVelocityCornerPPP(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);

    void addInternalVelocityCornerNNN(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
    void addInternalVelocityCornerNNP(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
    void addInternalVelocityCornerNPN(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
    void addInternalVelocityCornerNPP(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
    void addInternalVelocityCornerPNN(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
    void addInternalVelocityCornerPNP(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
    void addInternalVelocityCornerPPN(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
    void addInternalVelocityCornerPPP(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);

private:
    // PART I: Atomic-block version.

    template <int direction, int orientation>
    void addVelocityBoundary(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    template <int direction, int orientation>
    void addPressureBoundary(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    template <int plane, int normal1, int normal2>
    void addExternalVelocityEdge(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    template <int plane, int normal1, int normal2>
    void addInternalVelocityEdge(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    template <int normalX, int normalY, int normalZ>
    void addExternalVelocityCorner(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    template <int normalX, int normalY, int normalZ>
    void addInternalVelocityCorner(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);

    // PART II: Multi-block version.

    template <int direction, int orientation>
    void addVelocityBoundary(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    template <int direction, int orientation>
    void addPressureBoundary(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    template <int plane, int normal1, int normal2>
    void addExternalVelocityEdge(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    template <int plane, int normal1, int normal2>
    void addInternalVelocityEdge(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    template <int normalX, int normalY, int normalZ>
    void addExternalVelocityCorner(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
    template <int normalX, int normalY, int normalZ>
    void addInternalVelocityCorner(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
};

///////// class BoundaryConditionInstantiator3D ////////////////////////

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>
    *BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::clone() const
{
    return new BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>(*this);
}

// PART I: Atomic-block version.

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int direction, int orientation>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addVelocityBoundary(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    PLB_PRECONDITION(domain.x0 == domain.x1 || domain.y0 == domain.y1 || domain.z0 == domain.z1);

    // Convert (direction,orientation) description of the normal vector into a
    // (normalX,normalY,normalZ)
    //   description, as it is requried by the data processor for Neumann boundaries.
    enum {
        normalX = (direction == 0) ? orientation : 0,
        normalY = (direction == 1) ? orientation : 0,
        normalZ = (direction == 2) ? orientation : 0
    };

    // Instantiate the dynamics of the boundary as a composite dynamics, based on the one currently
    //   residing on the lattice.
    setCompositeDynamics(
        lattice, domain,
        BoundaryManager::template getVelocityBoundaryDynamics<direction, orientation>(
            new NoDynamics<T, Descriptor>));

    // In case an outflow condition is used, start by instantiating a data processor which copies
    //   all velocity values from the previous lattice cell.
    if (bcType == boundary::outflow || bcType == boundary::neumann) {
        integrateProcessingFunctional(
            new CopyVelocityFunctional3D<T, Descriptor, normalX, normalY, normalZ>, domain,
            lattice);
    }

    // In case a freeslip condition is used, start by instantiating a data processor which copies
    //   the tangential velocity values from the previous lattice cell.
    else if (bcType == boundary::freeslip)
    {
        integrateProcessingFunctional(
            new CopyTangentialVelocityFunctional3D<T, Descriptor, normalX, normalY, normalZ>,
            domain, lattice);
    }

    // In case a normal outflow condition is used, start by instantiating a data processor which
    // copies
    //   the normal velocity values from the previous lattice cell, and sets the other ones to zero.
    else if (bcType == boundary::normalOutflow)
    {
        integrateProcessingFunctional(
            new CopyNormalVelocityFunctional3D<T, Descriptor, normalX, normalY, normalZ>, domain,
            lattice);
    }

    // If the boundary condition has a non-local component, instantiate a corresponding data
    // processor.
    BoxProcessingFunctional3D_L<T, Descriptor> *functional =
        BoundaryManager::template getVelocityBoundaryFunctional<direction, orientation>();
    if (functional) {
        integrateProcessingFunctional(functional, domain, lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int direction, int orientation>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addPressureBoundary(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    PLB_PRECONDITION(domain.x0 == domain.x1 || domain.y0 == domain.y1 || domain.z0 == domain.z1);

    // Convert (direction,orientation) description of the normal vector into a
    // (normalX,normalY,normalZ)
    //   description, as it is requried by the data processor for Neumann boundaries.
    enum {
        normalX = (direction == 0) ? orientation : 0,
        normalY = (direction == 1) ? orientation : 0,
        normalZ = (direction == 2) ? orientation : 0
    };

    // Instantiate the dynamics of the boundary as a composite dynamics, based on the one currently
    //   residing on the lattice.
    setCompositeDynamics(
        lattice, domain,
        BoundaryManager::template getPressureBoundaryDynamics<direction, orientation>(
            new NoDynamics<T, Descriptor>));

    // In case a Neumann condition is used, start by instantiating a data processor which copies
    //   the density value from the previous lattice cell.
    if (bcType == boundary::neumann) {
        integrateProcessingFunctional(
            new CopyDensityFunctional3D<T, Descriptor, normalX, normalY, normalZ>, domain, lattice);
    }

    // If the boundary condition has a non-local component, instantiate a corresponding data
    // processor.
    BoxProcessingFunctional3D_L<T, Descriptor> *functional =
        BoundaryManager::template getPressureBoundaryFunctional<direction, orientation>();
    if (functional) {
        integrateProcessingFunctional(functional, domain, lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int plane, int normal1, int normal2>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    PLB_PRECONDITION(
        (domain.x0 == domain.x1 && domain.y0 == domain.y1)
        || (domain.x0 == domain.x1 && domain.z0 == domain.z1)
        || (domain.y0 == domain.y1 && domain.z0 == domain.z1));

    // Convert (plane,normal1,normal2) description of the normal vector into a
    // (normalX,normalY,normalZ)
    //   description, as it is requried by the data processor for Neumann boundaries.
    enum {
        normalX = plane == 0 ? 0 : (plane == 1 ? normal2 : normal1),
        normalY = plane == 1 ? 0 : (plane == 2 ? normal2 : normal1),
        normalZ = plane == 2 ? 0 : (plane == 0 ? normal2 : normal1)
    };

    setCompositeDynamics(
        lattice, domain,
        BoundaryManager::template getExternalVelocityEdgeDynamics<plane, normal1, normal2>(
            new NoDynamics<T, Descriptor>));

    if (bcType == boundary::neumann || bcType == boundary::outflow) {
        integrateProcessingFunctional(
            new CopyVelocityFunctional3D<T, Descriptor, normalX, normalY, normalZ>, domain,
            lattice);
    } else if (bcType == boundary::freeslip) {
        integrateProcessingFunctional(
            new CopyTangentialVelocityFunctional3D<T, Descriptor, normalX, normalY, normalZ>,
            domain, lattice);
    } else if (bcType == boundary::normalOutflow) {
        integrateProcessingFunctional(
            new CopyNormalVelocityFunctional3D<T, Descriptor, normalX, normalY, normalZ>, domain,
            lattice);
    }

    BoxProcessingFunctional3D_L<T, Descriptor> *functional =
        BoundaryManager::template getExternalVelocityEdgeFunctional<plane, normal1, normal2>();
    if (functional) {
        integrateProcessingFunctional(functional, domain, lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int plane, int normal1, int normal2>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    PLB_PRECONDITION(
        (domain.x0 == domain.x1 && domain.y0 == domain.y1)
        || (domain.x0 == domain.x1 && domain.z0 == domain.z1)
        || (domain.y0 == domain.y1 && domain.z0 == domain.z1));

    // Convert (plane,normal1,normal2) description of the normal vector into a
    // (normalX,normalY,normalZ)
    //   description, as it is requried by the data processor for Neumann boundaries.
    enum {
        normalX = plane == 0 ? 0 : (plane == 1 ? normal2 : normal1),
        normalY = plane == 1 ? 0 : (plane == 2 ? normal2 : normal1),
        normalZ = plane == 2 ? 0 : (plane == 0 ? normal2 : normal1)
    };

    setCompositeDynamics(
        lattice, domain,
        BoundaryManager::template getInternalVelocityEdgeDynamics<plane, normal1, normal2>(
            new NoDynamics<T, Descriptor>));

    if (bcType == boundary::neumann || bcType == boundary::outflow) {
        integrateProcessingFunctional(
            new CopyVelocityFunctional3D<T, Descriptor, normalX, normalY, normalZ>, domain,
            lattice);
    } else if (bcType == boundary::freeslip) {
        integrateProcessingFunctional(
            new CopyTangentialVelocityFunctional3D<T, Descriptor, normalX, normalY, normalZ>,
            domain, lattice);
    } else if (bcType == boundary::normalOutflow) {
        integrateProcessingFunctional(
            new CopyNormalVelocityFunctional3D<T, Descriptor, normalX, normalY, normalZ>, domain,
            lattice);
    }

    BoxProcessingFunctional3D_L<T, Descriptor> *functional =
        BoundaryManager::template getInternalVelocityEdgeFunctional<plane, normal1, normal2>();
    if (functional) {
        integrateProcessingFunctional(functional, domain, lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int xNormal, int yNormal, int zNormal>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityCorner(
    plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    setCompositeDynamics(
        lattice, Box3D(x, x, y, y, z, z),
        BoundaryManager::template getExternalVelocityCornerDynamics<xNormal, yNormal, zNormal>(
            new NoDynamics<T, Descriptor>));

    if (bcType == boundary::neumann || bcType == boundary::outflow) {
        integrateProcessingFunctional(
            new CopyVelocityFunctional3D<T, Descriptor, xNormal, yNormal, zNormal>,
            Box3D(x, x, y, y, z, z), lattice);
    } else if (bcType == boundary::freeslip) {
        integrateProcessingFunctional(
            new CopyTangentialVelocityFunctional3D<T, Descriptor, xNormal, yNormal, zNormal>,
            Box3D(x, x, y, y, z, z), lattice);
    } else if (bcType == boundary::normalOutflow) {
        integrateProcessingFunctional(
            new CopyNormalVelocityFunctional3D<T, Descriptor, xNormal, yNormal, zNormal>,
            Box3D(x, x, y, y, z, z), lattice);
    }

    BoxProcessingFunctional3D_L<T, Descriptor> *functional =
        BoundaryManager::template getExternalVelocityCornerFunctional<xNormal, yNormal, zNormal>();
    if (functional) {
        integrateProcessingFunctional(functional, Box3D(x, x, y, y, z, z), lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int xNormal, int yNormal, int zNormal>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityCorner(
    plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    setCompositeDynamics(
        lattice, Box3D(x, x, y, y, z, z),
        BoundaryManager::template getInternalVelocityCornerDynamics<xNormal, yNormal, zNormal>(
            new NoDynamics<T, Descriptor>));

    if (bcType == boundary::neumann || bcType == boundary::outflow) {
        integrateProcessingFunctional(
            new CopyVelocityFunctional3D<T, Descriptor, xNormal, yNormal, zNormal>,
            Box3D(x, x, y, y, z, z), lattice);
    }

    if (bcType == boundary::freeslip) {
        integrateProcessingFunctional(
            new CopyTangentialVelocityFunctional3D<T, Descriptor, xNormal, yNormal, zNormal>,
            Box3D(x, x, y, y, z, z), lattice);
    }

    if (bcType == boundary::normalOutflow) {
        integrateProcessingFunctional(
            new CopyNormalVelocityFunctional3D<T, Descriptor, xNormal, yNormal, zNormal>,
            Box3D(x, x, y, y, z, z), lattice);
    }

    BoxProcessingFunctional3D_L<T, Descriptor> *functional =
        BoundaryManager::template getInternalVelocityCornerFunctional<xNormal, yNormal, zNormal>();
    if (functional) {
        integrateProcessingFunctional(functional, Box3D(x, x, y, y, z, z), lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addVelocityBoundary0N(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addVelocityBoundary<0, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addVelocityBoundary0P(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addVelocityBoundary<0, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addVelocityBoundary1N(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addVelocityBoundary<1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addVelocityBoundary1P(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addVelocityBoundary<1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addVelocityBoundary2N(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addVelocityBoundary<2, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addVelocityBoundary2P(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addVelocityBoundary<2, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addPressureBoundary0N(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addPressureBoundary<0, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addPressureBoundary0P(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addPressureBoundary<0, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addPressureBoundary1N(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addPressureBoundary<1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addPressureBoundary1P(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addPressureBoundary<1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addPressureBoundary2N(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addPressureBoundary<2, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addPressureBoundary2P(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addPressureBoundary<2, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge0NN(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<0, -1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge0NP(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<0, -1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge0PN(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<0, 1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge0PP(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<0, 1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge1NN(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<1, -1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge1NP(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<1, -1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge1PN(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<1, 1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge1PP(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<1, 1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge2NN(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<2, -1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge2NP(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<2, -1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge2PN(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<2, 1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge2PP(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<2, 1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge0NN(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<0, -1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge0NP(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<0, -1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge0PN(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<0, 1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge0PP(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<0, 1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge1NN(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<1, -1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge1NP(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<1, -1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge1PN(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<1, 1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge1PP(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<1, 1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge2NN(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<2, -1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge2NP(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<2, -1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge2PN(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<2, 1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge2PP(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<2, 1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerNNN(
    plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<-1, -1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerNNP(
    plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<-1, -1, 1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerNPN(
    plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<-1, 1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerNPP(
    plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<-1, 1, 1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerPNN(
    plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<1, -1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerPNP(
    plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<1, -1, 1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerPPN(
    plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<1, 1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerPPP(
    plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<1, 1, 1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerNNN(
    plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<-1, -1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerNNP(
    plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<-1, -1, 1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerNPN(
    plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<-1, 1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerNPP(
    plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<-1, 1, 1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerPNN(
    plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<1, -1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerPNP(
    plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<1, -1, 1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerPPN(
    plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<1, 1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerPPP(
    plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<1, 1, 1>(x, y, z, lattice, bcType);
}

// PART II: Multi-block version.

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int direction, int orientation>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addVelocityBoundary(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    PLB_PRECONDITION(domain.x0 == domain.x1 || domain.y0 == domain.y1 || domain.z0 == domain.z1);

    // Convert (direction,orientation) description of the normal vector into a
    // (normalX,normalY,normalZ)
    //   description, as it is requried by the data processor for Neumann boundaries.
    enum {
        normalX = (direction == 0) ? orientation : 0,
        normalY = (direction == 1) ? orientation : 0,
        normalZ = (direction == 2) ? orientation : 0
    };

    // Instantiate the dynamics of the boundary as a composite dynamics, based on the one currently
    //   residing on the lattice.
    // (The NoDynamics is just a filler. The VelocityBoundaryDynamics actually
    // end up with the original dynamics on the lattice after this function is called!)
    setCompositeDynamics(
        lattice, domain,
        BoundaryManager::template getVelocityBoundaryDynamics<direction, orientation>(
            new NoDynamics<T, Descriptor>));

    // In case an outflow condition is used, start by instantiating a data processor which copies
    //   all velocity values from the previous lattice cell.
    if (bcType == boundary::outflow || bcType == boundary::neumann) {
        integrateProcessingFunctional(
            new CopyVelocityFunctional3D<T, Descriptor, normalX, normalY, normalZ>, domain,
            lattice);
    }

    // In case a freeslip condition is used, start by instantiating a data processor which copies
    //   the tangential velocity values from the previous lattice cell.
    else if (bcType == boundary::freeslip)
    {
        integrateProcessingFunctional(
            new CopyTangentialVelocityFunctional3D<T, Descriptor, normalX, normalY, normalZ>,
            domain, lattice);
    }

    // In case a normal outflow condition is used, start by instantiating a data processor which
    // copies
    //   the normal velocity values from the previous lattice cell, and sets the other ones to zero.
    else if (bcType == boundary::normalOutflow)
    {
        integrateProcessingFunctional(
            new CopyNormalVelocityFunctional3D<T, Descriptor, normalX, normalY, normalZ>, domain,
            lattice);
    }

    // If the boundary condition has a non-local component, instantiate a corresponding data
    // processor.
    BoxProcessingFunctional3D_L<T, Descriptor> *functional =
        BoundaryManager::template getVelocityBoundaryFunctional<direction, orientation>();
    if (functional) {
        integrateProcessingFunctional(functional, domain, lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int direction, int orientation>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addPressureBoundary(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    PLB_PRECONDITION(domain.x0 == domain.x1 || domain.y0 == domain.y1 || domain.z0 == domain.z1);

    // Convert (direction,orientation) description of the normal vector into a
    // (normalX,normalY,normalZ)
    //   description, as it is requried by the data processor for Neumann boundaries.
    enum {
        normalX = (direction == 0) ? orientation : 0,
        normalY = (direction == 1) ? orientation : 0,
        normalZ = (direction == 2) ? orientation : 0
    };

    // Instantiate the dynamics of the boundary as a composite dynamics, based on the one currently
    //   residing on the lattice.
    setCompositeDynamics(
        lattice, domain,
        BoundaryManager::template getPressureBoundaryDynamics<direction, orientation>(
            new NoDynamics<T, Descriptor>));

    // In case a Neumann condition is used, start by instantiating a data processor which copies
    //   the density value from the previous lattice cell.
    if (bcType == boundary::neumann) {
        integrateProcessingFunctional(
            new CopyDensityFunctional3D<T, Descriptor, normalX, normalY, normalZ>, domain, lattice);
    }

    // If the boundary condition has a non-local component, instantiate a corresponding data
    // processor.
    BoxProcessingFunctional3D_L<T, Descriptor> *functional =
        BoundaryManager::template getPressureBoundaryFunctional<direction, orientation>();
    if (functional) {
        integrateProcessingFunctional(functional, domain, lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int plane, int normal1, int normal2>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    PLB_PRECONDITION(
        (domain.x0 == domain.x1 && domain.y0 == domain.y1)
        || (domain.x0 == domain.x1 && domain.z0 == domain.z1)
        || (domain.y0 == domain.y1 && domain.z0 == domain.z1));

    // Convert (plane,normal1,normal2) description of the normal vector into a
    // (normalX,normalY,normalZ)
    //   description, as it is requried by the data processor for Neumann boundaries.
    enum {
        normalX = plane == 0 ? 0 : (plane == 1 ? normal2 : normal1),
        normalY = plane == 1 ? 0 : (plane == 2 ? normal2 : normal1),
        normalZ = plane == 2 ? 0 : (plane == 0 ? normal2 : normal1)
    };

    setCompositeDynamics(
        lattice, domain,
        BoundaryManager::template getExternalVelocityEdgeDynamics<plane, normal1, normal2>(
            new NoDynamics<T, Descriptor>));

    if (bcType == boundary::neumann || bcType == boundary::outflow) {
        integrateProcessingFunctional(
            new CopyVelocityFunctional3D<T, Descriptor, normalX, normalY, normalZ>, domain,
            lattice);
    } else if (bcType == boundary::freeslip) {
        integrateProcessingFunctional(
            new CopyTangentialVelocityFunctional3D<T, Descriptor, normalX, normalY, normalZ>,
            domain, lattice);
    } else if (bcType == boundary::normalOutflow) {
        integrateProcessingFunctional(
            new CopyNormalVelocityFunctional3D<T, Descriptor, normalX, normalY, normalZ>, domain,
            lattice);
    }

    BoxProcessingFunctional3D_L<T, Descriptor> *functional =
        BoundaryManager::template getExternalVelocityEdgeFunctional<plane, normal1, normal2>();
    if (functional) {
        integrateProcessingFunctional(functional, domain, lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int plane, int normal1, int normal2>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    PLB_PRECONDITION(
        (domain.x0 == domain.x1 && domain.y0 == domain.y1)
        || (domain.x0 == domain.x1 && domain.z0 == domain.z1)
        || (domain.y0 == domain.y1 && domain.z0 == domain.z1));

    // Convert (plane,normal1,normal2) description of the normal vector into a
    // (normalX,normalY,normalZ)
    //   description, as it is requried by the data processor for Neumann boundaries.
    enum {
        normalX = plane == 0 ? 0 : (plane == 1 ? normal2 : normal1),
        normalY = plane == 1 ? 0 : (plane == 2 ? normal2 : normal1),
        normalZ = plane == 2 ? 0 : (plane == 0 ? normal2 : normal1)
    };

    setCompositeDynamics(
        lattice, domain,
        BoundaryManager::template getInternalVelocityEdgeDynamics<plane, normal1, normal2>(
            new NoDynamics<T, Descriptor>));

    if (bcType == boundary::neumann || bcType == boundary::outflow) {
        integrateProcessingFunctional(
            new CopyVelocityFunctional3D<T, Descriptor, normalX, normalY, normalZ>, domain,
            lattice);
    } else if (bcType == boundary::freeslip) {
        integrateProcessingFunctional(
            new CopyTangentialVelocityFunctional3D<T, Descriptor, normalX, normalY, normalZ>,
            domain, lattice);
    } else if (bcType == boundary::normalOutflow) {
        integrateProcessingFunctional(
            new CopyNormalVelocityFunctional3D<T, Descriptor, normalX, normalY, normalZ>, domain,
            lattice);
    }

    BoxProcessingFunctional3D_L<T, Descriptor> *functional =
        BoundaryManager::template getInternalVelocityEdgeFunctional<plane, normal1, normal2>();
    if (functional) {
        integrateProcessingFunctional(functional, domain, lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int xNormal, int yNormal, int zNormal>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityCorner(
    plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    setCompositeDynamics(
        lattice, Box3D(x, x, y, y, z, z),
        BoundaryManager::template getExternalVelocityCornerDynamics<xNormal, yNormal, zNormal>(
            new NoDynamics<T, Descriptor>));

    if (bcType == boundary::neumann || bcType == boundary::outflow) {
        integrateProcessingFunctional(
            new CopyVelocityFunctional3D<T, Descriptor, xNormal, yNormal, zNormal>,
            Box3D(x, x, y, y, z, z), lattice);
    } else if (bcType == boundary::freeslip) {
        integrateProcessingFunctional(
            new CopyTangentialVelocityFunctional3D<T, Descriptor, xNormal, yNormal, zNormal>,
            Box3D(x, x, y, y, z, z), lattice);
    } else if (bcType == boundary::normalOutflow) {
        integrateProcessingFunctional(
            new CopyNormalVelocityFunctional3D<T, Descriptor, xNormal, yNormal, zNormal>,
            Box3D(x, x, y, y, z, z), lattice);
    }

    BoxProcessingFunctional3D_L<T, Descriptor> *functional =
        BoundaryManager::template getExternalVelocityCornerFunctional<xNormal, yNormal, zNormal>();
    if (functional) {
        integrateProcessingFunctional(functional, Box3D(x, x, y, y, z, z), lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int xNormal, int yNormal, int zNormal>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityCorner(
    plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    setCompositeDynamics(
        lattice, Box3D(x, x, y, y, z, z),
        BoundaryManager::template getInternalVelocityCornerDynamics<xNormal, yNormal, zNormal>(
            new NoDynamics<T, Descriptor>));

    if (bcType == boundary::neumann || bcType == boundary::outflow) {
        integrateProcessingFunctional(
            new CopyVelocityFunctional3D<T, Descriptor, xNormal, yNormal, zNormal>,
            Box3D(x, x, y, y, z, z), lattice);
    }

    if (bcType == boundary::freeslip) {
        integrateProcessingFunctional(
            new CopyTangentialVelocityFunctional3D<T, Descriptor, xNormal, yNormal, zNormal>,
            Box3D(x, x, y, y, z, z), lattice);
    }

    if (bcType == boundary::normalOutflow) {
        integrateProcessingFunctional(
            new CopyNormalVelocityFunctional3D<T, Descriptor, xNormal, yNormal, zNormal>,
            Box3D(x, x, y, y, z, z), lattice);
    }

    BoxProcessingFunctional3D_L<T, Descriptor> *functional =
        BoundaryManager::template getInternalVelocityCornerFunctional<xNormal, yNormal, zNormal>();
    if (functional) {
        integrateProcessingFunctional(functional, Box3D(x, x, y, y, z, z), lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addVelocityBoundary0N(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addVelocityBoundary<0, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addVelocityBoundary0P(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addVelocityBoundary<0, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addVelocityBoundary1N(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addVelocityBoundary<1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addVelocityBoundary1P(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addVelocityBoundary<1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addVelocityBoundary2N(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addVelocityBoundary<2, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addVelocityBoundary2P(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addVelocityBoundary<2, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addPressureBoundary0N(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addPressureBoundary<0, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addPressureBoundary0P(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addPressureBoundary<0, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addPressureBoundary1N(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addPressureBoundary<1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addPressureBoundary1P(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addPressureBoundary<1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addPressureBoundary2N(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addPressureBoundary<2, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addPressureBoundary2P(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addPressureBoundary<2, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge0NN(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<0, -1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge0NP(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<0, -1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge0PN(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<0, 1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge0PP(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<0, 1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge1NN(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<1, -1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge1NP(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<1, -1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge1PN(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<1, 1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge1PP(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<1, 1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge2NN(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<2, -1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge2NP(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<2, -1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge2PN(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<2, 1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityEdge2PP(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityEdge<2, 1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge0NN(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<0, -1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge0NP(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<0, -1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge0PN(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<0, 1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge0PP(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<0, 1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge1NN(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<1, -1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge1NP(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<1, -1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge1PN(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<1, 1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge1PP(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<1, 1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge2NN(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<2, -1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge2NP(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<2, -1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge2PN(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<2, 1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityEdge2PP(
    Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityEdge<2, 1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerNNN(
    plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<-1, -1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerNNP(
    plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<-1, -1, 1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerNPN(
    plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<-1, 1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerNPP(
    plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<-1, 1, 1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerPNN(
    plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<1, -1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerPNP(
    plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<1, -1, 1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerPPN(
    plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<1, 1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerPPP(
    plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<1, 1, 1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerNNN(
    plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<-1, -1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerNNP(
    plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<-1, -1, 1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerNPN(
    plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<-1, 1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerNPP(
    plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<-1, 1, 1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerPNN(
    plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<1, -1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerPNP(
    plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<1, -1, 1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerPPN(
    plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<1, 1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerPPP(
    plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<1, 1, 1>(x, y, z, lattice, bcType);
}

}  // namespace plb

#endif
