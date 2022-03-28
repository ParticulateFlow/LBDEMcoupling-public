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
#ifndef BOUNDARY_INSTANTIATOR_2D_H
#define BOUNDARY_INSTANTIATOR_2D_H

#include "atomicBlock/atomicBlockOperations2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "boundaryCondition/boundaryCondition2D.h"
#include "boundaryCondition/neumannCondition2D.h"
#include "core/cell.h"
#include "core/geometry2D.h"
#include "core/globalDefs.h"
#include "dataProcessors/dataInitializerFunctional2D.h"
#include "dataProcessors/dataInitializerWrapper2D.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiBlockOperations2D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
class BoundaryConditionInstantiator2D : public OnLatticeBoundaryCondition2D<T, Descriptor> {
public:
    virtual BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager> *clone() const;

    // PART I: Atomic-block version.

    void addVelocityBoundary0N(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addVelocityBoundary0P(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addVelocityBoundary1N(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addVelocityBoundary1P(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addPressureBoundary0N(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addPressureBoundary0P(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addPressureBoundary1N(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addPressureBoundary1P(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addExternalVelocityCornerNN(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityCornerNP(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityCornerPN(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityCornerPP(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addInternalVelocityCornerNN(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityCornerNP(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityCornerPN(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityCornerPP(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    // PART II: Multi-block version.

    void addVelocityBoundary0N(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addVelocityBoundary0P(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addVelocityBoundary1N(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addVelocityBoundary1P(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addPressureBoundary0N(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addPressureBoundary0P(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addPressureBoundary1N(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addPressureBoundary1P(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addExternalVelocityCornerNN(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityCornerNP(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityCornerPN(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addExternalVelocityCornerPP(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addInternalVelocityCornerNN(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityCornerNP(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityCornerPN(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addInternalVelocityCornerPP(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);

private:
    // PART I: Atomic-block version.
    template <int direction, int orientation>
    void addVelocityBoundary(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    template <int direction, int orientation>
    void addPressureBoundary(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    template <int normalX, int normalY>
    void addExternalVelocityCorner(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    template <int normalX, int normalY>
    void addInternalVelocityCorner(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);

    // PART II: Multi-block version.

    template <int direction, int orientation>
    void addVelocityBoundary(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    template <int direction, int orientation>
    void addPressureBoundary(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    template <int normalX, int normalY>
    void addExternalVelocityCorner(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    template <int normalX, int normalY>
    void addInternalVelocityCorner(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
};

///////// class BoundaryConditionInstantiator2D ////////////////////////

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>
    *BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::clone() const
{
    return new BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>(*this);
}

// PART I: Atomic-block version.

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int direction, int orientation>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addVelocityBoundary(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    PLB_PRECONDITION(domain.x0 == domain.x1 || domain.y0 == domain.y1);

    // Convert (direction,orientation) description of the normal vector into a (normalX,normalY)
    //   description, as it is requried by the data processor for Neumann boundaries.
    enum {
        normalX = (direction == 0) ? orientation : 0,
        normalY = (direction == 1) ? orientation : 0
    };

    // Instantiate the dynamics of the boundary as a composite dynamics, based on the one currently
    //   residing on the lattice.
    setCompositeDynamics(
        lattice, domain,
        BoundaryManager::template getVelocityBoundaryDynamics<direction, orientation>(
            new NoDynamics<T, Descriptor>));

    // In case an outflow condition is used, start by instantiating a data processor which copies
    //   all velocity value from the previous lattice cell.
    if (bcType == boundary::outflow || bcType == boundary::neumann) {
        integrateProcessingFunctional(
            new CopyVelocityFunctional2D<T, Descriptor, normalX, normalY>, domain, lattice);
    }
    // In case a normal outflow condition is used, start by instantiating a data processor which
    // copies
    //   the normal velocity value from the previous lattice cell, and sets the other components to
    //   zero.
    if (bcType == boundary::normalOutflow) {
        integrateProcessingFunctional(
            new CopyNormalVelocityFunctional2D<T, Descriptor, normalX, normalY>, domain, lattice);
    } else
        // In case a freeslip condition is used, start by instantiating a data processor which
        // copies
        //   the tangential velocity values from the previous lattice cell.
        if (bcType == boundary::freeslip) {
            integrateProcessingFunctional(
                new CopyTangentialVelocityFunctional2D<T, Descriptor, normalX, normalY>, domain,
                lattice);
        }

    // If the boundary condition has a non-local component, instantiate a corresponding data
    // processor.
    BoxProcessingFunctional2D_L<T, Descriptor> *functional =
        BoundaryManager::template getVelocityBoundaryFunctional<direction, orientation>();
    if (functional) {
        integrateProcessingFunctional(functional, domain, lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int direction, int orientation>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addPressureBoundary(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    PLB_PRECONDITION(domain.x0 == domain.x1 || domain.y0 == domain.y1);

    // Convert (direction,orientation) description of the normal vector into a (normalX,normalY)
    //   description, as it is requried by the data processor for Neumann boundaries.
    enum {
        normalX = (direction == 0) ? orientation : 0,
        normalY = (direction == 1) ? orientation : 0
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
            new CopyDensityFunctional2D<T, Descriptor, normalX, normalY>, domain, lattice);
    }

    // If the boundary condition has a non-local component, instantiate a corresponding data
    // processor.
    BoxProcessingFunctional2D_L<T, Descriptor> *functional =
        BoundaryManager::template getPressureBoundaryFunctional<direction, orientation>();
    if (functional) {
        integrateProcessingFunctional(functional, domain, lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int xNormal, int yNormal>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addExternalVelocityCorner(
    plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    setCompositeDynamics(
        lattice, Box2D(x, x, y, y),
        BoundaryManager::template getExternalVelocityCornerDynamics<xNormal, yNormal>(
            new NoDynamics<T, Descriptor>));

    if (bcType == boundary::neumann || bcType == boundary::outflow) {
        integrateProcessingFunctional(
            new CopyVelocityFunctional2D<T, Descriptor, xNormal, yNormal>, Box2D(x, x, y, y),
            lattice);
    } else if (bcType == boundary::freeslip) {
        integrateProcessingFunctional(
            new CopyTangentialVelocityFunctional2D<T, Descriptor, xNormal, yNormal>,
            Box2D(x, x, y, y), lattice);
    } else if (bcType == boundary::normalOutflow) {
        integrateProcessingFunctional(
            new CopyNormalVelocityFunctional2D<T, Descriptor, xNormal, yNormal>, Box2D(x, x, y, y),
            lattice);
    }

    BoxProcessingFunctional2D_L<T, Descriptor> *functional =
        BoundaryManager::template getExternalVelocityCornerFunctional<xNormal, yNormal>();
    if (functional) {
        integrateProcessingFunctional(functional, Box2D(x, x, y, y), lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int xNormal, int yNormal>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addInternalVelocityCorner(
    plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    setCompositeDynamics(
        lattice, Box2D(x, x, y, y),
        BoundaryManager::template getInternalVelocityCornerDynamics<xNormal, yNormal>(
            new NoDynamics<T, Descriptor>));

    if (bcType == boundary::neumann || bcType == boundary::outflow) {
        integrateProcessingFunctional(
            new CopyVelocityFunctional2D<T, Descriptor, xNormal, yNormal>, Box2D(x, x, y, y),
            lattice);
    } else if (bcType == boundary::freeslip) {
        integrateProcessingFunctional(
            new CopyTangentialVelocityFunctional2D<T, Descriptor, xNormal, yNormal>,
            Box2D(x, x, y, y), lattice);
    } else if (bcType == boundary::normalOutflow) {
        integrateProcessingFunctional(
            new CopyNormalVelocityFunctional2D<T, Descriptor, xNormal, yNormal>, Box2D(x, x, y, y),
            lattice);
    }

    BoxProcessingFunctional2D_L<T, Descriptor> *functional =
        BoundaryManager::template getInternalVelocityCornerFunctional<xNormal, yNormal>();
    if (functional) {
        integrateProcessingFunctional(functional, Box2D(x, x, y, y), lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addVelocityBoundary0N(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addVelocityBoundary<0, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addVelocityBoundary0P(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addVelocityBoundary<0, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addVelocityBoundary1N(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addVelocityBoundary<1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addVelocityBoundary1P(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addVelocityBoundary<1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addPressureBoundary0N(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addPressureBoundary<0, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addPressureBoundary0P(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addPressureBoundary<0, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addPressureBoundary1N(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addPressureBoundary<1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addPressureBoundary1P(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addPressureBoundary<1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerNN(
    plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<-1, -1>(x, y, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerNP(
    plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<-1, 1>(x, y, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerPN(
    plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<1, -1>(x, y, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerPP(
    plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<1, 1>(x, y, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerNN(
    plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<-1, -1>(x, y, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerNP(
    plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<-1, 1>(x, y, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerPN(
    plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<1, -1>(x, y, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerPP(
    plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<1, 1>(x, y, lattice, bcType);
}

// PART II: Multi-block version.

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int direction, int orientation>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addVelocityBoundary(
    Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    PLB_PRECONDITION(domain.x0 == domain.x1 || domain.y0 == domain.y1);

    // Convert (direction,orientation) description of the normal vector into a (normalX,normalY)
    //   description, as it is requried by the data processor for Neumann boundaries.
    enum {
        normalX = (direction == 0) ? orientation : 0,
        normalY = (direction == 1) ? orientation : 0
    };

    // Instantiate the dynamics of the boundary as a composite dynamics, based on the one currently
    //   residing on the lattice.
    setCompositeDynamics(
        lattice, domain,
        BoundaryManager::template getVelocityBoundaryDynamics<direction, orientation>(
            new NoDynamics<T, Descriptor>));

    // In case an outflow condition is used, start by instantiating a data processor which copies
    //   all velocity value from the previous lattice cell.
    if (bcType == boundary::outflow || bcType == boundary::neumann) {
        integrateProcessingFunctional(
            new CopyVelocityFunctional2D<T, Descriptor, normalX, normalY>, domain, lattice);
    }
    // In case a normal outflow condition is used, start by instantiating a data processor which
    // copies
    //   the normal velocity value from the previous lattice cell, and sets the other components to
    //   zero.
    if (bcType == boundary::normalOutflow) {
        integrateProcessingFunctional(
            new CopyNormalVelocityFunctional2D<T, Descriptor, normalX, normalY>, domain, lattice);
    } else
        // In case a freeslip condition is used, start by instantiating a data processor which
        // copies
        //   the tangential velocity values from the previous lattice cell.
        if (bcType == boundary::freeslip) {
            integrateProcessingFunctional(
                new CopyTangentialVelocityFunctional2D<T, Descriptor, normalX, normalY>, domain,
                lattice);
        }

    // If the boundary condition has a non-local component, instantiate a corresponding data
    // processor.
    BoxProcessingFunctional2D_L<T, Descriptor> *functional =
        BoundaryManager::template getVelocityBoundaryFunctional<direction, orientation>();
    if (functional) {
        integrateProcessingFunctional(functional, domain, lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int direction, int orientation>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addPressureBoundary(
    Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    PLB_PRECONDITION(domain.x0 == domain.x1 || domain.y0 == domain.y1);

    // Convert (direction,orientation) description of the normal vector into a (normalX,normalY)
    //   description, as it is requried by the data processor for Neumann boundaries.
    enum {
        normalX = (direction == 0) ? orientation : 0,
        normalY = (direction == 1) ? orientation : 0
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
            new CopyDensityFunctional2D<T, Descriptor, normalX, normalY>, domain, lattice);
    }

    // If the boundary condition has a non-local component, instantiate a corresponding data
    // processor.
    BoxProcessingFunctional2D_L<T, Descriptor> *functional =
        BoundaryManager::template getPressureBoundaryFunctional<direction, orientation>();
    if (functional) {
        integrateProcessingFunctional(functional, domain, lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int xNormal, int yNormal>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addExternalVelocityCorner(
    plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    setCompositeDynamics(
        lattice, Box2D(x, x, y, y),
        BoundaryManager::template getExternalVelocityCornerDynamics<xNormal, yNormal>(
            new NoDynamics<T, Descriptor>));

    if (bcType == boundary::neumann || bcType == boundary::outflow) {
        integrateProcessingFunctional(
            new CopyVelocityFunctional2D<T, Descriptor, xNormal, yNormal>, Box2D(x, x, y, y),
            lattice);
    } else if (bcType == boundary::freeslip) {
        integrateProcessingFunctional(
            new CopyTangentialVelocityFunctional2D<T, Descriptor, xNormal, yNormal>,
            Box2D(x, x, y, y), lattice);
    } else if (bcType == boundary::normalOutflow) {
        integrateProcessingFunctional(
            new CopyNormalVelocityFunctional2D<T, Descriptor, xNormal, yNormal>, Box2D(x, x, y, y),
            lattice);
    }

    BoxProcessingFunctional2D_L<T, Descriptor> *functional =
        BoundaryManager::template getExternalVelocityCornerFunctional<xNormal, yNormal>();
    if (functional) {
        integrateProcessingFunctional(functional, Box2D(x, x, y, y), lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int xNormal, int yNormal>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addInternalVelocityCorner(
    plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    setCompositeDynamics(
        lattice, Box2D(x, x, y, y),
        BoundaryManager::template getInternalVelocityCornerDynamics<xNormal, yNormal>(
            new NoDynamics<T, Descriptor>));

    if (bcType == boundary::neumann || bcType == boundary::outflow) {
        integrateProcessingFunctional(
            new CopyVelocityFunctional2D<T, Descriptor, xNormal, yNormal>, Box2D(x, x, y, y),
            lattice);
    } else if (bcType == boundary::freeslip) {
        integrateProcessingFunctional(
            new CopyTangentialVelocityFunctional2D<T, Descriptor, xNormal, yNormal>,
            Box2D(x, x, y, y), lattice);
    } else if (bcType == boundary::normalOutflow) {
        integrateProcessingFunctional(
            new CopyNormalVelocityFunctional2D<T, Descriptor, xNormal, yNormal>, Box2D(x, x, y, y),
            lattice);
    }

    BoxProcessingFunctional2D_L<T, Descriptor> *functional =
        BoundaryManager::template getInternalVelocityCornerFunctional<xNormal, yNormal>();
    if (functional) {
        integrateProcessingFunctional(functional, Box2D(x, x, y, y), lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addVelocityBoundary0N(
    Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addVelocityBoundary<0, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addVelocityBoundary0P(
    Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addVelocityBoundary<0, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addVelocityBoundary1N(
    Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addVelocityBoundary<1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addVelocityBoundary1P(
    Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addVelocityBoundary<1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addPressureBoundary0N(
    Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addPressureBoundary<0, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addPressureBoundary0P(
    Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addPressureBoundary<0, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addPressureBoundary1N(
    Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addPressureBoundary<1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addPressureBoundary1P(
    Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addPressureBoundary<1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerNN(
    plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<-1, -1>(x, y, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerNP(
    plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<-1, 1>(x, y, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerPN(
    plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<1, -1>(x, y, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addExternalVelocityCornerPP(
    plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addExternalVelocityCorner<1, 1>(x, y, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerNN(
    plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<-1, -1>(x, y, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerNP(
    plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<-1, 1>(x, y, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerPN(
    plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<1, -1>(x, y, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::addInternalVelocityCornerPP(
    plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addInternalVelocityCorner<1, 1>(x, y, lattice, bcType);
}

}  // namespace plb

#endif  // BOUNDARY_INSTANTIATOR_2D_H
