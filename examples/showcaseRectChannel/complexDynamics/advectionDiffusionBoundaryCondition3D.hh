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
 * A helper for initialising 3D boundaries -- generic implementation.
 */
#ifndef ADVECTION_DIFFUSION_BOUNDARY_CONDITION_3D_HH
#define ADVECTION_DIFFUSION_BOUNDARY_CONDITION_3D_HH

#include "complexDynamics/advectionDiffusionBoundaries.h"
#include "complexDynamics/advectionDiffusionBoundaryCondition3D.h"
#include "complexDynamics/advectionDiffusionBoundaryInstantiator3D.h"

namespace plb {

// template<typename T, template<typename U> class Descriptor>
// void
// OnLatticeAdvectionDiffusionBoundaryCondition3D<T,Descriptor>::setTemperatureConditionOnBlockBoundaries
// (
//         BlockLattice3D<T,Descriptor>& lattice, boundary::BcType bcType )
//{
//     plint nx = lattice.getNx();
//     plint ny = lattice.getNy();
//     plint nz = lattice.getNz();

//    addTemperatureBoundary0N(Box3D(   0,   0,   1,ny-2,   1,nz-2), lattice, bcType);
//    addTemperatureBoundary0P(Box3D(nx-1,nx-1,   1,ny-2,   1,nz-2), lattice, bcType);
//    addTemperatureBoundary1N(Box3D(   1,nx-2,   0,   0,   1,nz-2), lattice, bcType);
//    addTemperatureBoundary1P(Box3D(   1,nx-2,ny-1,ny-1,   1,nz-2), lattice, bcType);
//    addTemperatureBoundary2N(Box3D(   1,nx-2,   1,ny-2,   0,   0), lattice, bcType);
//    addTemperatureBoundary2P(Box3D(   1,nx-2,   1,ny-2,nz-1,nz-1), lattice, bcType);

//    addTemperatureEdge0NN(Box3D(   1,nx-2,   0,   0,   0,   0), lattice, bcType);
//    addTemperatureEdge0NP(Box3D(   1,nx-2,   0,   0,nz-1,nz-1), lattice, bcType);
//    addTemperatureEdge0PN(Box3D(   1,nx-2,ny-1,ny-1,   0,   0), lattice, bcType);
//    addTemperatureEdge0PP(Box3D(   1,nx-2,ny-1,ny-1,nz-1,nz-1), lattice, bcType);

//    addTemperatureEdge1NN(Box3D(   0,   0,   1,ny-2,   0,   0), lattice, bcType);
//    addTemperatureEdge1NP(Box3D(nx-1,nx-1,   1,ny-2,   0,   0), lattice, bcType);
//    addTemperatureEdge1PN(Box3D(   0,   0,   1,ny-2,nz-1,nz-1), lattice, bcType);
//    addTemperatureEdge1PP(Box3D(nx-1,nx-1,   1,ny-2,nz-1,nz-1), lattice, bcType);

//    addTemperatureEdge2NN(Box3D(   0,   0,   0,   0,   1,nz-2), lattice, bcType);
//    addTemperatureEdge2NP(Box3D(   0,   0,ny-1,ny-1,   1,nz-2), lattice, bcType);
//    addTemperatureEdge2PN(Box3D(nx-1,nx-1,   0,   0,   1,nz-2), lattice, bcType);
//    addTemperatureEdge2PP(Box3D(nx-1,nx-1,ny-1,ny-1,   1,nz-2), lattice, bcType);

//    addTemperatureCornerNNN(   0,   0,   0, lattice, bcType);
//    addTemperatureCornerNNP(   0,   0,nz-1, lattice, bcType);
//    addTemperatureCornerNPN(   0,ny-1,   0, lattice, bcType);
//    addTemperatureCornerNPP(   0,ny-1,nz-1, lattice, bcType);
//    addTemperatureCornerPNN(nx-1,   0,   0, lattice, bcType);
//    addTemperatureCornerPNP(nx-1,   0,nz-1, lattice, bcType);
//    addTemperatureCornerPPN(nx-1,ny-1,   0, lattice, bcType);
//    addTemperatureCornerPPP(nx-1,ny-1,nz-1, lattice, bcType);
//}

// template<typename T, template<typename U> class Descriptor>
// void
// OnLatticeAdvectionDiffusionBoundaryCondition3D<T,Descriptor>::setTemperatureConditionOnBlockBoundaries
// (
//         MultiBlockLattice3D<T,Descriptor>& lattice, boundary::BcType bcType )
//{
//     plint nx = lattice.getNx();
//     plint ny = lattice.getNy();
//     plint nz = lattice.getNz();

//    addTemperatureBoundary0N(Box3D(   0,   0,   1,ny-2,   1,nz-2), lattice, bcType);
//    addTemperatureBoundary0P(Box3D(nx-1,nx-1,   1,ny-2,   1,nz-2), lattice, bcType);
//    addTemperatureBoundary1N(Box3D(   1,nx-2,   0,   0,   1,nz-2), lattice, bcType);
//    addTemperatureBoundary1P(Box3D(   1,nx-2,ny-1,ny-1,   1,nz-2), lattice, bcType);
//    addTemperatureBoundary2N(Box3D(   1,nx-2,   1,ny-2,   0,   0), lattice, bcType);
//    addTemperatureBoundary2P(Box3D(   1,nx-2,   1,ny-2,nz-1,nz-1), lattice, bcType);

//    addTemperatureEdge0NN(Box3D(   1,nx-2,   0,   0,   0,   0), lattice, bcType);
//    addTemperatureEdge0NP(Box3D(   1,nx-2,   0,   0,nz-1,nz-1), lattice, bcType);
//    addTemperatureEdge0PN(Box3D(   1,nx-2,ny-1,ny-1,   0,   0), lattice, bcType);
//    addTemperatureEdge0PP(Box3D(   1,nx-2,ny-1,ny-1,nz-1,nz-1), lattice, bcType);

//    addTemperatureEdge1NN(Box3D(   0,   0,   1,ny-2,   0,   0), lattice, bcType);
//    addTemperatureEdge1NP(Box3D(nx-1,nx-1,   1,ny-2,   0,   0), lattice, bcType);
//    addTemperatureEdge1PN(Box3D(   0,   0,   1,ny-2,nz-1,nz-1), lattice, bcType);
//    addTemperatureEdge1PP(Box3D(nx-1,nx-1,   1,ny-2,nz-1,nz-1), lattice, bcType);

//    addTemperatureEdge2NN(Box3D(   0,   0,   0,   0,   1,nz-2), lattice, bcType);
//    addTemperatureEdge2NP(Box3D(   0,   0,ny-1,ny-1,   1,nz-2), lattice, bcType);
//    addTemperatureEdge2PN(Box3D(nx-1,nx-1,   0,   0,   1,nz-2), lattice, bcType);
//    addTemperatureEdge2PP(Box3D(nx-1,nx-1,ny-1,ny-1,   1,nz-2), lattice, bcType);

//    addTemperatureCornerNNN(   0,   0,   0, lattice, bcType);
//    addTemperatureCornerNNP(   0,   0,nz-1, lattice, bcType);
//    addTemperatureCornerNPN(   0,ny-1,   0, lattice, bcType);
//    addTemperatureCornerNPP(   0,ny-1,nz-1, lattice, bcType);
//    addTemperatureCornerPNN(nx-1,   0,   0, lattice, bcType);
//    addTemperatureCornerPNP(nx-1,   0,nz-1, lattice, bcType);
//    addTemperatureCornerPPN(nx-1,ny-1,   0, lattice, bcType);
//    addTemperatureCornerPPP(nx-1,ny-1,nz-1, lattice, bcType);
//}

template <typename T, template <typename U> class Descriptor>
void OnLatticeAdvectionDiffusionBoundaryCondition3D<T, Descriptor>::
    setTemperatureConditionOnBlockBoundaries(
        BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    setTemperatureConditionOnBlockBoundaries(lattice, lattice.getBoundingBox(), bcType);
}

template <typename T, template <typename U> class Descriptor>
void OnLatticeAdvectionDiffusionBoundaryCondition3D<T, Descriptor>::
    setTemperatureConditionOnBlockBoundaries(
        BlockLattice3D<T, Descriptor> &lattice, Box3D applicationDomain, boundary::BcType bcType)
{
    setTemperatureConditionOnBlockBoundaries(
        lattice, lattice.getBoundingBox(), applicationDomain, bcType);
}

template <typename T, template <typename U> class Descriptor>
void OnLatticeAdvectionDiffusionBoundaryCondition3D<T, Descriptor>::
    setTemperatureConditionOnBlockBoundaries(
        BlockLattice3D<T, Descriptor> &lattice, Box3D block, Box3D applicationDomain,
        boundary::BcType bcType)
{
    plint boundaryWidth = 1;
    BlockSurface3D surf(block, boundaryWidth);
    Box3D intersection;
    if (intersect(surf.surface0N(), applicationDomain, intersection)) {
        addTemperatureBoundary0N(intersection, lattice, bcType);
    }
    if (intersect(surf.surface0P(), applicationDomain, intersection)) {
        addTemperatureBoundary0P(intersection, lattice, bcType);
    }
    if (intersect(surf.surface1N(), applicationDomain, intersection)) {
        addTemperatureBoundary1N(intersection, lattice, bcType);
    }
    if (intersect(surf.surface1P(), applicationDomain, intersection)) {
        addTemperatureBoundary1P(intersection, lattice, bcType);
    }
    if (intersect(surf.surface2N(), applicationDomain, intersection)) {
        addTemperatureBoundary2N(intersection, lattice, bcType);
    }
    if (intersect(surf.surface2P(), applicationDomain, intersection)) {
        addTemperatureBoundary2P(intersection, lattice, bcType);
    }

    if (intersect(surf.edge0NN(), applicationDomain, intersection)) {
        addTemperatureEdge0NN(intersection, lattice, bcType);
    }
    if (intersect(surf.edge0NP(), applicationDomain, intersection)) {
        addTemperatureEdge0NP(intersection, lattice, bcType);
    }
    if (intersect(surf.edge0PN(), applicationDomain, intersection)) {
        addTemperatureEdge0PN(intersection, lattice, bcType);
    }
    if (intersect(surf.edge0PP(), applicationDomain, intersection)) {
        addTemperatureEdge0PP(intersection, lattice, bcType);
    }

    if (intersect(surf.edge1NN(), applicationDomain, intersection)) {
        addTemperatureEdge1NN(intersection, lattice, bcType);
    }
    if (intersect(surf.edge1NP(), applicationDomain, intersection)) {
        addTemperatureEdge1NP(intersection, lattice, bcType);
    }
    if (intersect(surf.edge1PN(), applicationDomain, intersection)) {
        addTemperatureEdge1PN(intersection, lattice, bcType);
    }
    if (intersect(surf.edge1PP(), applicationDomain, intersection)) {
        addTemperatureEdge1PP(intersection, lattice, bcType);
    }

    if (intersect(surf.edge2NN(), applicationDomain, intersection)) {
        addTemperatureEdge2NN(intersection, lattice, bcType);
    }
    if (intersect(surf.edge2NP(), applicationDomain, intersection)) {
        addTemperatureEdge2NP(intersection, lattice, bcType);
    }
    if (intersect(surf.edge2PN(), applicationDomain, intersection)) {
        addTemperatureEdge2PN(intersection, lattice, bcType);
    }
    if (intersect(surf.edge2PP(), applicationDomain, intersection)) {
        addTemperatureEdge2PP(intersection, lattice, bcType);
    }

    if (intersect(surf.cornerNNN(), applicationDomain, intersection)) {
        addTemperatureCornerNNN(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (intersect(surf.cornerNNP(), applicationDomain, intersection)) {
        addTemperatureCornerNNP(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (intersect(surf.cornerNPN(), applicationDomain, intersection)) {
        addTemperatureCornerNPN(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (intersect(surf.cornerNPP(), applicationDomain, intersection)) {
        addTemperatureCornerNPP(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (intersect(surf.cornerPNN(), applicationDomain, intersection)) {
        addTemperatureCornerPNN(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (intersect(surf.cornerPNP(), applicationDomain, intersection)) {
        addTemperatureCornerPNP(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (intersect(surf.cornerPPN(), applicationDomain, intersection)) {
        addTemperatureCornerPPN(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (intersect(surf.cornerPPP(), applicationDomain, intersection)) {
        addTemperatureCornerPPP(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
}

template <typename T, template <typename U> class Descriptor>
void OnLatticeAdvectionDiffusionBoundaryCondition3D<T, Descriptor>::
    setTemperatureConditionOnBlockBoundaries(
        MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    setTemperatureConditionOnBlockBoundaries(lattice, lattice.getBoundingBox(), bcType);
}

template <typename T, template <typename U> class Descriptor>
void OnLatticeAdvectionDiffusionBoundaryCondition3D<T, Descriptor>::
    setTemperatureConditionOnBlockBoundaries(
        MultiBlockLattice3D<T, Descriptor> &lattice, Box3D applicationDomain,
        boundary::BcType bcType)
{
    setTemperatureConditionOnBlockBoundaries(
        lattice, lattice.getBoundingBox(), applicationDomain, bcType);
}

template <typename T, template <typename U> class Descriptor>
void OnLatticeAdvectionDiffusionBoundaryCondition3D<T, Descriptor>::
    setTemperatureConditionOnBlockBoundaries(
        MultiBlockLattice3D<T, Descriptor> &lattice, Box3D block, Box3D applicationDomain,
        boundary::BcType bcType)
{
    plint boundaryWidth = 1;
    BlockSurface3D surf(block, boundaryWidth);
    Box3D intersection;
    if (intersect(surf.surface0N(), applicationDomain, intersection)) {
        addTemperatureBoundary0N(intersection, lattice, bcType);
    }
    if (intersect(surf.surface0P(), applicationDomain, intersection)) {
        addTemperatureBoundary0P(intersection, lattice, bcType);
    }
    if (intersect(surf.surface1N(), applicationDomain, intersection)) {
        addTemperatureBoundary1N(intersection, lattice, bcType);
    }
    if (intersect(surf.surface1P(), applicationDomain, intersection)) {
        addTemperatureBoundary1P(intersection, lattice, bcType);
    }
    if (intersect(surf.surface2N(), applicationDomain, intersection)) {
        addTemperatureBoundary2N(intersection, lattice, bcType);
    }
    if (intersect(surf.surface2P(), applicationDomain, intersection)) {
        addTemperatureBoundary2P(intersection, lattice, bcType);
    }

    if (intersect(surf.edge0NN(), applicationDomain, intersection)) {
        addTemperatureEdge0NN(intersection, lattice, bcType);
    }
    if (intersect(surf.edge0NP(), applicationDomain, intersection)) {
        addTemperatureEdge0NP(intersection, lattice, bcType);
    }
    if (intersect(surf.edge0PN(), applicationDomain, intersection)) {
        addTemperatureEdge0PN(intersection, lattice, bcType);
    }
    if (intersect(surf.edge0PP(), applicationDomain, intersection)) {
        addTemperatureEdge0PP(intersection, lattice, bcType);
    }

    if (intersect(surf.edge1NN(), applicationDomain, intersection)) {
        addTemperatureEdge1NN(intersection, lattice, bcType);
    }
    if (intersect(surf.edge1NP(), applicationDomain, intersection)) {
        addTemperatureEdge1NP(intersection, lattice, bcType);
    }
    if (intersect(surf.edge1PN(), applicationDomain, intersection)) {
        addTemperatureEdge1PN(intersection, lattice, bcType);
    }
    if (intersect(surf.edge1PP(), applicationDomain, intersection)) {
        addTemperatureEdge1PP(intersection, lattice, bcType);
    }

    if (intersect(surf.edge2NN(), applicationDomain, intersection)) {
        addTemperatureEdge2NN(intersection, lattice, bcType);
    }
    if (intersect(surf.edge2NP(), applicationDomain, intersection)) {
        addTemperatureEdge2NP(intersection, lattice, bcType);
    }
    if (intersect(surf.edge2PN(), applicationDomain, intersection)) {
        addTemperatureEdge2PN(intersection, lattice, bcType);
    }
    if (intersect(surf.edge2PP(), applicationDomain, intersection)) {
        addTemperatureEdge2PP(intersection, lattice, bcType);
    }

    if (intersect(surf.cornerNNN(), applicationDomain, intersection)) {
        addTemperatureCornerNNN(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (intersect(surf.cornerNNP(), applicationDomain, intersection)) {
        addTemperatureCornerNNP(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (intersect(surf.cornerNPN(), applicationDomain, intersection)) {
        addTemperatureCornerNPN(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (intersect(surf.cornerNPP(), applicationDomain, intersection)) {
        addTemperatureCornerNPP(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (intersect(surf.cornerPNN(), applicationDomain, intersection)) {
        addTemperatureCornerPNN(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (intersect(surf.cornerPNP(), applicationDomain, intersection)) {
        addTemperatureCornerPNP(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (intersect(surf.cornerPPN(), applicationDomain, intersection)) {
        addTemperatureCornerPPN(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (intersect(surf.cornerPPP(), applicationDomain, intersection)) {
        addTemperatureCornerPPP(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
}

template <typename T, template <typename U> class Descriptor>
class AdvectionDiffusionBoundaryManager3D {
public:
    template <int direction, int orientation>
    static BoundaryCompositeDynamics<T, Descriptor> *getTemperatureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics);

    template <int direction, int orientation>
    static BoxProcessingFunctional3D_L<T, Descriptor> *getTemperatureBoundaryProcessor(
        Box3D domain);

    template <int plane, int normal1, int normal2>
    static BoundaryCompositeDynamics<T, Descriptor> *getTemperatureEdgeDynamics(
        Dynamics<T, Descriptor> *baseDynamics);

    template <int plane, int normal1, int normal2>
    static BoxProcessingFunctional3D_L<T, Descriptor> *getTemperatureEdgeProcessor(Box3D domain);

    template <int xNormal, int yNormal, int zNormal>
    static BoundaryCompositeDynamics<T, Descriptor> *getTemperatureCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int xNormal, int yNormal, int zNormal>
    static BoxProcessingFunctional3D_L<T, Descriptor> *getTemperatureCornerProcessor(
        plint x, plint y, plint z);
};

////////// AdvectionDiffusionBoundaryManager3D /////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *AdvectionDiffusionBoundaryManager3D<T, Descriptor>::getTemperatureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new AdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional3D_L<T, Descriptor>
    *AdvectionDiffusionBoundaryManager3D<T, Descriptor>::getTemperatureBoundaryProcessor(
        Box3D domain)
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int plane, int normal1, int normal2>
BoundaryCompositeDynamics<T, Descriptor>
    *AdvectionDiffusionBoundaryManager3D<T, Descriptor>::getTemperatureEdgeDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new AdvectionDiffusionEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>(
        baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int plane, int normal1, int normal2>
BoxProcessingFunctional3D_L<T, Descriptor>
    *AdvectionDiffusionBoundaryManager3D<T, Descriptor>::getTemperatureEdgeProcessor(Box3D domain)
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal, int zNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *AdvectionDiffusionBoundaryManager3D<T, Descriptor>::getTemperatureCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new AdvectionDiffusionCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal>(
        baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal, int zNormal>
BoxProcessingFunctional3D_L<T, Descriptor>
    *AdvectionDiffusionBoundaryManager3D<T, Descriptor>::getTemperatureCornerProcessor(
        plint x, plint y, plint z)
{
    return 0;
}

////////// RegularizedAdvectionDiffusionBoundaryManager3D /////////////

template <typename T, template <typename U> class Descriptor>
class RegularizedAdvectionDiffusionBoundaryManager3D {
public:
    template <int direction, int orientation>
    static BoundaryCompositeDynamics<T, Descriptor> *getTemperatureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics);

    template <int direction, int orientation>
    static BoxProcessingFunctional3D_L<T, Descriptor> *getTemperatureBoundaryProcessor(
        Box3D domain);

    template <int plane, int normal1, int normal2>
    static BoundaryCompositeDynamics<T, Descriptor> *getTemperatureEdgeDynamics(
        Dynamics<T, Descriptor> *baseDynamics);

    template <int plane, int normal1, int normal2>
    static BoxProcessingFunctional3D_L<T, Descriptor> *getTemperatureEdgeProcessor(Box3D domain);

    template <int xNormal, int yNormal, int zNormal>
    static BoundaryCompositeDynamics<T, Descriptor> *getTemperatureCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int xNormal, int yNormal, int zNormal>
    static BoxProcessingFunctional3D_L<T, Descriptor> *getTemperatureCornerProcessor(
        plint x, plint y, plint z);
};

////////// RegularizedAdvectionDiffusionBoundaryManager3D /////////////

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *RegularizedAdvectionDiffusionBoundaryManager3D<T, Descriptor>::getTemperatureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new RegularizedAdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional3D_L<T, Descriptor>
    *RegularizedAdvectionDiffusionBoundaryManager3D<T, Descriptor>::getTemperatureBoundaryProcessor(
        Box3D domain)
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int plane, int normal1, int normal2>
BoundaryCompositeDynamics<T, Descriptor>
    *RegularizedAdvectionDiffusionBoundaryManager3D<T, Descriptor>::getTemperatureEdgeDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new AdvectionDiffusionEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>(
        baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int plane, int normal1, int normal2>
BoxProcessingFunctional3D_L<T, Descriptor>
    *RegularizedAdvectionDiffusionBoundaryManager3D<T, Descriptor>::getTemperatureEdgeProcessor(
        Box3D domain)
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal, int zNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *RegularizedAdvectionDiffusionBoundaryManager3D<T, Descriptor>::getTemperatureCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new AdvectionDiffusionCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal>(
        baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal, int zNormal>
BoxProcessingFunctional3D_L<T, Descriptor>
    *RegularizedAdvectionDiffusionBoundaryManager3D<T, Descriptor>::getTemperatureCornerProcessor(
        plint x, plint y, plint z)
{
    return 0;
}

////////// CompleteRegularizedAdvectionDiffusionBoundaryManager3D /////////////

template <typename T, template <typename U> class Descriptor>
class CompleteRegularizedAdvectionDiffusionBoundaryManager3D {
public:
    template <int direction, int orientation>
    static BoundaryCompositeDynamics<T, Descriptor> *getTemperatureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics);

    template <int direction, int orientation>
    static BoxProcessingFunctional3D_L<T, Descriptor> *getTemperatureBoundaryProcessor(
        Box3D domain);

    template <int plane, int normal1, int normal2>
    static BoundaryCompositeDynamics<T, Descriptor> *getTemperatureEdgeDynamics(
        Dynamics<T, Descriptor> *baseDynamics);

    template <int plane, int normal1, int normal2>
    static BoxProcessingFunctional3D_L<T, Descriptor> *getTemperatureEdgeProcessor(Box3D domain);

    template <int xNormal, int yNormal, int zNormal>
    static BoundaryCompositeDynamics<T, Descriptor> *getTemperatureCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int xNormal, int yNormal, int zNormal>
    static BoxProcessingFunctional3D_L<T, Descriptor> *getTemperatureCornerProcessor(
        plint x, plint y, plint z);
};

////////// CompleteRegularizedAdvectionDiffusionBoundaryManager3D /////////////

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor> *CompleteRegularizedAdvectionDiffusionBoundaryManager3D<
    T, Descriptor>::getTemperatureBoundaryDynamics(Dynamics<T, Descriptor> *baseDynamics)
{
    return new CompleteRegularizedAdvectionDiffusionBoundaryDynamics<
        T, Descriptor, direction, orientation>(baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional3D_L<T, Descriptor> *CompleteRegularizedAdvectionDiffusionBoundaryManager3D<
    T, Descriptor>::getTemperatureBoundaryProcessor(Box3D domain)
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int plane, int normal1, int normal2>
BoundaryCompositeDynamics<T, Descriptor> *CompleteRegularizedAdvectionDiffusionBoundaryManager3D<
    T, Descriptor>::getTemperatureEdgeDynamics(Dynamics<T, Descriptor> *baseDynamics)
{
    return new StoreDensityDynamics<T, Descriptor>(baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int plane, int normal1, int normal2>
BoxProcessingFunctional3D_L<T, Descriptor> *CompleteRegularizedAdvectionDiffusionBoundaryManager3D<
    T, Descriptor>::getTemperatureEdgeProcessor(Box3D domain)
{
    return new CompleteAdvectionDiffusionEdgeBoundaryFunctional3D<
        T, Descriptor, plane, normal1, normal2>();
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal, int zNormal>
BoundaryCompositeDynamics<T, Descriptor> *CompleteRegularizedAdvectionDiffusionBoundaryManager3D<
    T, Descriptor>::getTemperatureCornerDynamics(Dynamics<T, Descriptor> *baseDynamics)
{
    return new StoreDensityDynamics<T, Descriptor>(baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal, int zNormal>
BoxProcessingFunctional3D_L<T, Descriptor> *CompleteRegularizedAdvectionDiffusionBoundaryManager3D<
    T, Descriptor>::getTemperatureCornerProcessor(plint x, plint y, plint z)
{
    return new CompleteAdvectionDiffusionCornerBoundaryFunctional3D<
        T, Descriptor, xNormal, yNormal, zNormal>();
}

////////// Factory functions //////////////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
OnLatticeAdvectionDiffusionBoundaryCondition3D<T, Descriptor>
    *createLocalAdvectionDiffusionBoundaryCondition3D()
{
    return new AdvectionDiffusionBoundaryConditionInstantiator3D<
        T, Descriptor, AdvectionDiffusionBoundaryManager3D<T, Descriptor> >();
}

template <typename T, template <typename U> class Descriptor>
OnLatticeAdvectionDiffusionBoundaryCondition3D<T, Descriptor>
    *createLocalRegularizedAdvectionDiffusionBoundaryCondition3D()
{
    return new AdvectionDiffusionBoundaryConditionInstantiator3D<
        T, Descriptor, RegularizedAdvectionDiffusionBoundaryManager3D<T, Descriptor> >();
}

template <typename T, template <typename U> class Descriptor>
OnLatticeAdvectionDiffusionBoundaryCondition3D<T, Descriptor>
    *createLocalCompleteRegularizedAdvectionDiffusionBoundaryCondition3D()
{
    return new AdvectionDiffusionBoundaryConditionInstantiator3D<
        T, Descriptor, CompleteRegularizedAdvectionDiffusionBoundaryManager3D<T, Descriptor> >();
}

}  // namespace plb

#endif
