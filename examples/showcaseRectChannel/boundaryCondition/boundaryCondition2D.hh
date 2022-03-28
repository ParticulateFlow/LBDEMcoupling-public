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
 * A helper for initialising 2D boundaries -- generic implementation.
 */
#ifndef BOUNDARY_CONDITION_2D_HH
#define BOUNDARY_CONDITION_2D_HH

#include "atomicBlock/blockLattice2D.h"
#include "boundaryCondition/boundaryCondition2D.h"
#include "boundaryCondition/boundaryInstantiator2D.h"
#include "boundaryCondition/equilibriumBoundaryDynamics.h"
#include "boundaryCondition/regularizedBoundaryDynamics2D.h"
#include "boundaryCondition/wrappedLocalBoundaryProcessor2D.h"
#include "core/blockSurface2D.h"
#include "core/plbDebug.h"
#include "multiBlock/multiBlockLattice2D.h"

namespace plb {

// PART I: Atomic-block version.

template <typename T, template <typename U> class Descriptor>
void OnLatticeBoundaryCondition2D<T, Descriptor>::setVelocityConditionOnBlockBoundaries(
    BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    setVelocityConditionOnBlockBoundaries(lattice, lattice.getBoundingBox(), bcType);
}

template <typename T, template <typename U> class Descriptor>
void OnLatticeBoundaryCondition2D<T, Descriptor>::setVelocityConditionOnBlockBoundaries(
    BlockLattice2D<T, Descriptor> &lattice, Box2D applicationDomain, boundary::BcType bcType)
{
    setVelocityConditionOnBlockBoundaries(
        lattice, lattice.getBoundingBox(), applicationDomain, bcType);
}

template <typename T, template <typename U> class Descriptor>
void OnLatticeBoundaryCondition2D<T, Descriptor>::setVelocityConditionOnBlockBoundaries(
    BlockLattice2D<T, Descriptor> &lattice, Box2D block, Box2D applicationDomain,
    boundary::BcType bcType)
{
    plint boundaryWidth = 1;
    BlockSurface2D surf(block, boundaryWidth);
    Box2D intersection;
    if (intersect(surf.edge0N(), applicationDomain, intersection)) {
        addVelocityBoundary0N(intersection, lattice, bcType);
    }
    if (intersect(surf.edge0P(), applicationDomain, intersection)) {
        addVelocityBoundary0P(intersection, lattice, bcType);
    }
    if (intersect(surf.edge1N(), applicationDomain, intersection)) {
        addVelocityBoundary1N(intersection, lattice, bcType);
    }
    if (intersect(surf.edge1P(), applicationDomain, intersection)) {
        addVelocityBoundary1P(intersection, lattice, bcType);
    }

    if (intersect(surf.cornerNN(), applicationDomain, intersection)) {
        addExternalVelocityCornerNN(intersection.x0, intersection.y0, lattice, bcType);
    }
    if (intersect(surf.cornerNP(), applicationDomain, intersection)) {
        addExternalVelocityCornerNP(intersection.x0, intersection.y0, lattice, bcType);
    }
    if (intersect(surf.cornerPN(), applicationDomain, intersection)) {
        addExternalVelocityCornerPN(intersection.x0, intersection.y0, lattice, bcType);
    }
    if (intersect(surf.cornerPP(), applicationDomain, intersection)) {
        addExternalVelocityCornerPP(intersection.x0, intersection.y0, lattice, bcType);
    }
}

template <typename T, template <typename U> class Descriptor>
void OnLatticeBoundaryCondition2D<T, Descriptor>::setPressureConditionOnBlockBoundaries(
    BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    setPressureConditionOnBlockBoundaries(lattice, lattice.getBoundingBox(), bcType);
}

template <typename T, template <typename U> class Descriptor>
void OnLatticeBoundaryCondition2D<T, Descriptor>::setPressureConditionOnBlockBoundaries(
    BlockLattice2D<T, Descriptor> &lattice, Box2D applicationDomain, boundary::BcType bcType)
{
    setPressureConditionOnBlockBoundaries(
        lattice, lattice.getBoundingBox(), applicationDomain, bcType);
}

template <typename T, template <typename U> class Descriptor>
void OnLatticeBoundaryCondition2D<T, Descriptor>::setPressureConditionOnBlockBoundaries(
    BlockLattice2D<T, Descriptor> &lattice, Box2D block, Box2D applicationDomain,
    boundary::BcType bcType)
{
    plint boundaryWidth = 1;
    BlockSurface2D surf(block, boundaryWidth);
    Box2D intersection;
    if (intersect(surf.edge0N(), applicationDomain, intersection)) {
        addPressureBoundary0N(intersection, lattice, bcType);
    }
    if (intersect(surf.edge0P(), applicationDomain, intersection)) {
        addPressureBoundary0P(intersection, lattice, bcType);
    }
    if (intersect(surf.edge1N(), applicationDomain, intersection)) {
        addPressureBoundary1N(intersection, lattice, bcType);
    }
    if (intersect(surf.edge1P(), applicationDomain, intersection)) {
        addPressureBoundary1P(intersection, lattice, bcType);
    }

    /*
    if (intersect(surf.cornerNN(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    if (intersect(surf.cornerNP(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    if (intersect(surf.cornerPN(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    if (intersect(surf.cornerPP(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    */
}

// PART II: Multi-block version.

template <typename T, template <typename U> class Descriptor>
void OnLatticeBoundaryCondition2D<T, Descriptor>::setVelocityConditionOnBlockBoundaries(
    MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    setVelocityConditionOnBlockBoundaries(lattice, lattice.getBoundingBox(), bcType);
}

template <typename T, template <typename U> class Descriptor>
void OnLatticeBoundaryCondition2D<T, Descriptor>::setVelocityConditionOnBlockBoundaries(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D applicationDomain, boundary::BcType bcType)
{
    setVelocityConditionOnBlockBoundaries(
        lattice, lattice.getBoundingBox(), applicationDomain, bcType);
}

template <typename T, template <typename U> class Descriptor>
void OnLatticeBoundaryCondition2D<T, Descriptor>::setVelocityConditionOnBlockBoundaries(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D block, Box2D applicationDomain,
    boundary::BcType bcType)
{
    plint boundaryWidth = 1;
    BlockSurface2D surf(block, boundaryWidth);
    Box2D intersection;
    if (intersect(surf.edge0N(), applicationDomain, intersection)) {
        addVelocityBoundary0N(intersection, lattice, bcType);
    }
    if (intersect(surf.edge0P(), applicationDomain, intersection)) {
        addVelocityBoundary0P(intersection, lattice, bcType);
    }
    if (intersect(surf.edge1N(), applicationDomain, intersection)) {
        addVelocityBoundary1N(intersection, lattice, bcType);
    }
    if (intersect(surf.edge1P(), applicationDomain, intersection)) {
        addVelocityBoundary1P(intersection, lattice, bcType);
    }

    if (intersect(surf.cornerNN(), applicationDomain, intersection)) {
        addExternalVelocityCornerNN(intersection.x0, intersection.y0, lattice, bcType);
    }
    if (intersect(surf.cornerNP(), applicationDomain, intersection)) {
        addExternalVelocityCornerNP(intersection.x0, intersection.y0, lattice, bcType);
    }
    if (intersect(surf.cornerPN(), applicationDomain, intersection)) {
        addExternalVelocityCornerPN(intersection.x0, intersection.y0, lattice, bcType);
    }
    if (intersect(surf.cornerPP(), applicationDomain, intersection)) {
        addExternalVelocityCornerPP(intersection.x0, intersection.y0, lattice, bcType);
    }
}

template <typename T, template <typename U> class Descriptor>
void OnLatticeBoundaryCondition2D<T, Descriptor>::setPressureConditionOnBlockBoundaries(
    MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    setPressureConditionOnBlockBoundaries(lattice, lattice.getBoundingBox(), bcType);
}

template <typename T, template <typename U> class Descriptor>
void OnLatticeBoundaryCondition2D<T, Descriptor>::setPressureConditionOnBlockBoundaries(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D applicationDomain, boundary::BcType bcType)
{
    setPressureConditionOnBlockBoundaries(
        lattice, lattice.getBoundingBox(), applicationDomain, bcType);
}

template <typename T, template <typename U> class Descriptor>
void OnLatticeBoundaryCondition2D<T, Descriptor>::setPressureConditionOnBlockBoundaries(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D block, Box2D applicationDomain,
    boundary::BcType bcType)
{
    plint boundaryWidth = 1;
    BlockSurface2D surf(block, boundaryWidth);
    Box2D intersection;
    if (intersect(surf.edge0N(), applicationDomain, intersection)) {
        addPressureBoundary0N(intersection, lattice, bcType);
    }
    if (intersect(surf.edge0P(), applicationDomain, intersection)) {
        addPressureBoundary0P(intersection, lattice, bcType);
    }
    if (intersect(surf.edge1N(), applicationDomain, intersection)) {
        addPressureBoundary1N(intersection, lattice, bcType);
    }
    if (intersect(surf.edge1P(), applicationDomain, intersection)) {
        addPressureBoundary1P(intersection, lattice, bcType);
    }

    if (intersect(surf.cornerNN(), applicationDomain, intersection)) {
        PLB_ASSERT(false);
    }
    if (intersect(surf.cornerNP(), applicationDomain, intersection)) {
        PLB_ASSERT(false);
    }
    if (intersect(surf.cornerPN(), applicationDomain, intersection)) {
        PLB_ASSERT(false);
    }
    if (intersect(surf.cornerPP(), applicationDomain, intersection)) {
        PLB_ASSERT(false);
    }
}

template <typename T, template <typename U> class Descriptor>
void OnLatticeBoundaryCondition2D<T, Descriptor>::setVelocityConditionOnInnerBlockBoundaries(
    MultiBlockLattice2D<T, Descriptor> &lattice, Box2D block, boundary::BcType bcType)
{
    // Create Velocity boundary conditions on the domain
    addVelocityBoundary0N(Box2D(block.x1, block.x1, block.y0 + 1, block.y1 - 1), lattice, bcType);
    addVelocityBoundary0P(Box2D(block.x0, block.x0, block.y0 + 1, block.y1 - 1), lattice, bcType);

    addVelocityBoundary1N(Box2D(block.x0 + 1, block.x1 - 1, block.y1, block.y1), lattice, bcType);
    addVelocityBoundary1P(Box2D(block.x0 + 1, block.x1 - 1, block.y0, block.y0), lattice, bcType);

    addInternalVelocityCornerPP(block.x0, block.y0, lattice, bcType);
    addInternalVelocityCornerPN(block.x0, block.y1, lattice, bcType);
    addInternalVelocityCornerNP(block.x1, block.y0, lattice, bcType);
    addInternalVelocityCornerNN(block.x1, block.y1, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor>
class RegularizedBoundaryManager2D {
public:
    template <int direction, int orientation>
    static BoundaryCompositeDynamics<T, Descriptor> *getVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int direction, int orientation>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getVelocityBoundaryFunctional();

    template <int direction, int orientation>
    static BoundaryCompositeDynamics<T, Descriptor> *getPressureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int direction, int orientation>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getPressureBoundaryFunctional();

    template <int xNormal, int yNormal>
    static BoundaryCompositeDynamics<T, Descriptor> *getExternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int xNormal, int yNormal>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getExternalVelocityCornerFunctional();

    template <int xNormal, int yNormal>
    static BoundaryCompositeDynamics<T, Descriptor> *getInternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int xNormal, int yNormal>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getInternalVelocityCornerFunctional();
};

template <typename T, template <typename U> class Descriptor>
class WrappedRegularizedBoundaryManager2D {
public:
    template <int direction, int orientation>
    static BoundaryCompositeDynamics<T, Descriptor> *getVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int direction, int orientation>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getVelocityBoundaryFunctional();

    template <int direction, int orientation>
    static BoundaryCompositeDynamics<T, Descriptor> *getPressureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int direction, int orientation>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getPressureBoundaryFunctional();

    template <int xNormal, int yNormal>
    static BoundaryCompositeDynamics<T, Descriptor> *getExternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int xNormal, int yNormal>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getExternalVelocityCornerFunctional();

    template <int xNormal, int yNormal>
    static BoundaryCompositeDynamics<T, Descriptor> *getInternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int xNormal, int yNormal>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getInternalVelocityCornerFunctional();
};

template <typename T, template <typename U> class Descriptor>
class EquilibriumBoundaryManager2D {
public:
    template <int direction, int orientation>
    static BoundaryCompositeDynamics<T, Descriptor> *getVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int direction, int orientation>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getVelocityBoundaryFunctional();

    template <int direction, int orientation>
    static BoundaryCompositeDynamics<T, Descriptor> *getPressureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int direction, int orientation>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getPressureBoundaryFunctional();

    template <int xNormal, int yNormal>
    static BoundaryCompositeDynamics<T, Descriptor> *getExternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int xNormal, int yNormal>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getExternalVelocityCornerFunctional();

    template <int xNormal, int yNormal>
    static BoundaryCompositeDynamics<T, Descriptor> *getInternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int xNormal, int yNormal>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getInternalVelocityCornerFunctional();
};

template <typename T, template <typename U> class Descriptor>
class InterpolationBoundaryManager2D {
public:
    template <int direction, int orientation>
    static BoundaryCompositeDynamics<T, Descriptor> *getVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int direction, int orientation>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getVelocityBoundaryFunctional();

    template <int direction, int orientation>
    static BoundaryCompositeDynamics<T, Descriptor> *getPressureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int direction, int orientation>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getPressureBoundaryFunctional();

    template <int xNormal, int yNormal>
    static BoundaryCompositeDynamics<T, Descriptor> *getExternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int xNormal, int yNormal>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getExternalVelocityCornerFunctional();

    template <int xNormal, int yNormal>
    static BoundaryCompositeDynamics<T, Descriptor> *getInternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int xNormal, int yNormal>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getInternalVelocityCornerFunctional();
};

////////// RegularizedBoundaryManager2D /////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *RegularizedBoundaryManager2D<T, Descriptor>::getVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new RegularizedVelocityBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional2D_L<T, Descriptor>
    *RegularizedBoundaryManager2D<T, Descriptor>::getVelocityBoundaryFunctional()
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *RegularizedBoundaryManager2D<T, Descriptor>::getPressureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new RegularizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional2D_L<T, Descriptor>
    *RegularizedBoundaryManager2D<T, Descriptor>::getPressureBoundaryFunctional()
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *RegularizedBoundaryManager2D<T, Descriptor>::getExternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new StoreVelocityDynamics<T, Descriptor>(baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor>
    *RegularizedBoundaryManager2D<T, Descriptor>::getExternalVelocityCornerFunctional()
{
    return new OuterVelocityCornerFunctional2D<T, Descriptor, xNormal, yNormal>;
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *RegularizedBoundaryManager2D<T, Descriptor>::getInternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new RegularizedVelocityInnerCornerDynamics2D<T, Descriptor, xNormal, yNormal>(
        baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor>
    *RegularizedBoundaryManager2D<T, Descriptor>::getInternalVelocityCornerFunctional()
{
    return 0;
}

////////// WrappedRegularizedBoundaryManager2D /////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedRegularizedBoundaryManager2D<T, Descriptor>::getVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    bool automaticPrepareCollision = false;
    return new RegularizedVelocityBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics, automaticPrepareCollision);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional2D_L<T, Descriptor>
    *WrappedRegularizedBoundaryManager2D<T, Descriptor>::getVelocityBoundaryFunctional()
{
    return new WrappedLocalBoundaryFunctional2D<T, Descriptor>;
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedRegularizedBoundaryManager2D<T, Descriptor>::getPressureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    bool automaticPrepareCollision = false;
    return new RegularizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics, automaticPrepareCollision);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional2D_L<T, Descriptor>
    *WrappedRegularizedBoundaryManager2D<T, Descriptor>::getPressureBoundaryFunctional()
{
    return new WrappedLocalBoundaryFunctional2D<T, Descriptor>;
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedRegularizedBoundaryManager2D<T, Descriptor>::getExternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new StoreVelocityDynamics<T, Descriptor>(baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor>
    *WrappedRegularizedBoundaryManager2D<T, Descriptor>::getExternalVelocityCornerFunctional()
{
    return new OuterVelocityCornerFunctional2D<T, Descriptor, xNormal, yNormal>;
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedRegularizedBoundaryManager2D<T, Descriptor>::getInternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    bool automaticPrepareCollision = false;
    return new RegularizedVelocityInnerCornerDynamics2D<T, Descriptor, xNormal, yNormal>(
        baseDynamics, automaticPrepareCollision);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor>
    *WrappedRegularizedBoundaryManager2D<T, Descriptor>::getInternalVelocityCornerFunctional()
{
    return new WrappedLocalBoundaryFunctional2D<T, Descriptor>;
}

////////// EquilibriumBoundaryManager2D /////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *EquilibriumBoundaryManager2D<T, Descriptor>::getVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new EquilibriumVelocityBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional2D_L<T, Descriptor>
    *EquilibriumBoundaryManager2D<T, Descriptor>::getVelocityBoundaryFunctional()
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *EquilibriumBoundaryManager2D<T, Descriptor>::getPressureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new EquilibriumDensityBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional2D_L<T, Descriptor>
    *EquilibriumBoundaryManager2D<T, Descriptor>::getPressureBoundaryFunctional()
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *EquilibriumBoundaryManager2D<T, Descriptor>::getExternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new StoreVelocityDynamics<T, Descriptor>(baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor>
    *EquilibriumBoundaryManager2D<T, Descriptor>::getExternalVelocityCornerFunctional()
{
    return new OuterVelocityCornerFunctional2D<T, Descriptor, xNormal, yNormal>;
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *EquilibriumBoundaryManager2D<T, Descriptor>::getInternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new RegularizedVelocityInnerCornerDynamics2D<T, Descriptor, xNormal, yNormal>(
        baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor>
    *EquilibriumBoundaryManager2D<T, Descriptor>::getInternalVelocityCornerFunctional()
{
    return 0;
}

////////// InterpolationBoundaryManager2D /////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *InterpolationBoundaryManager2D<T, Descriptor>::getVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new VelocityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional2D_L<T, Descriptor>
    *InterpolationBoundaryManager2D<T, Descriptor>::getVelocityBoundaryFunctional()
{
    return new StraightFdBoundaryFunctional2D<T, Descriptor, direction, orientation>;
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *InterpolationBoundaryManager2D<T, Descriptor>::getPressureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional2D_L<T, Descriptor>
    *InterpolationBoundaryManager2D<T, Descriptor>::getPressureBoundaryFunctional()
{
    return new StraightFdBoundaryFunctional2D<T, Descriptor, direction, orientation>;
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *InterpolationBoundaryManager2D<T, Descriptor>::getExternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new StoreVelocityDynamics<T, Descriptor>(baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor>
    *InterpolationBoundaryManager2D<T, Descriptor>::getExternalVelocityCornerFunctional()
{
    return new OuterVelocityCornerFunctional2D<T, Descriptor, xNormal, yNormal>;
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *InterpolationBoundaryManager2D<T, Descriptor>::getInternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new RegularizedVelocityInnerCornerDynamics2D<T, Descriptor, xNormal, yNormal>(
        baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor>
    *InterpolationBoundaryManager2D<T, Descriptor>::getInternalVelocityCornerFunctional()
{
    return 0;
}

////////// Factory functions //////////////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
OnLatticeBoundaryCondition2D<T, Descriptor> *createLocalBoundaryCondition2D()
{
    // For the default boundary condition, use the wrapped version which uses
    //   data processors. This guarantees that other data processor based
    //   algorithms work properly on boundaries; example: grid refinement.
    return new BoundaryConditionInstantiator2D<
        T, Descriptor, WrappedRegularizedBoundaryManager2D<T, Descriptor> >;
}

template <typename T, template <typename U> class Descriptor>
OnLatticeBoundaryCondition2D<T, Descriptor> *createDynamicsBasedLocalBoundaryCondition2D()
{
    // Purely local, no data proessors (except in external corners).
    return new BoundaryConditionInstantiator2D<
        T, Descriptor, RegularizedBoundaryManager2D<T, Descriptor> >;
}

template <typename T, template <typename U> class Descriptor>
OnLatticeBoundaryCondition2D<T, Descriptor> *createEquilibriumBoundaryCondition2D()
{
    return new BoundaryConditionInstantiator2D<
        T, Descriptor, EquilibriumBoundaryManager2D<T, Descriptor> >;
}

template <typename T, template <typename U> class Descriptor>
OnLatticeBoundaryCondition2D<T, Descriptor> *createInterpBoundaryCondition2D()
{
    return new BoundaryConditionInstantiator2D<
        T, Descriptor, InterpolationBoundaryManager2D<T, Descriptor> >;
}

}  // namespace plb

#endif  // BOUNDARY_CONDITION_2D_HH
