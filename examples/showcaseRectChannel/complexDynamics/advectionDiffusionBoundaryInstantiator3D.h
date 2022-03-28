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
#ifndef ADVECTION_DIFFUSION_BOUNDARY_INSTANTIATOR_3D_H
#define ADVECTION_DIFFUSION_BOUNDARY_INSTANTIATOR_3D_H

#include "complexDynamics/advectionDiffusionBoundaryCondition3D.h"
#include "core/cell.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "dataProcessors/dataInitializerFunctional3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
class AdvectionDiffusionBoundaryConditionInstantiator3D :
    public OnLatticeAdvectionDiffusionBoundaryCondition3D<T, Descriptor> {
public:
    AdvectionDiffusionBoundaryConditionInstantiator3D();

    void addTemperatureBoundary0N(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureBoundary0P(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet);
    void addTemperatureBoundary1N(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureBoundary1P(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureBoundary2N(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureBoundary2P(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addTemperatureEdge0NN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureEdge0NP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureEdge0PN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureEdge0PP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureEdge1NN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureEdge1NP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureEdge1PN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureEdge1PP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureEdge2NN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureEdge2NP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureEdge2PN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureEdge2PP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addTemperatureCornerNNN(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureCornerNNP(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureCornerNPN(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureCornerNPP(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureCornerPNN(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureCornerPNP(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureCornerPPN(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureCornerPPP(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addTemperatureBoundary0N(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureBoundary0P(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureBoundary1N(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureBoundary1P(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureBoundary2N(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureBoundary2P(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addTemperatureEdge0NN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureEdge0NP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureEdge0PN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureEdge0PP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureEdge1NN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureEdge1NP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureEdge1PN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureEdge1PP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureEdge2NN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureEdge2NP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureEdge2PN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureEdge2PP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addTemperatureCornerNNN(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
    void addTemperatureCornerNNP(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
    void addTemperatureCornerNPN(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
    void addTemperatureCornerNPP(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
    void addTemperatureCornerPNN(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
    void addTemperatureCornerPNP(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
    void addTemperatureCornerPPN(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
    void addTemperatureCornerPPP(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);

private:
    template <int direction, int orientation>
    void addTemperatureBoundary(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    template <int plane, int normal1, int normal2>
    void addTemperatureEdge(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    template <int normalX, int normalY, int normalZ>
    void addTemperatureCorner(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);

    template <int direction, int orientation>
    void addTemperatureBoundary(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    template <int plane, int normal1, int normal2>
    void addTemperatureEdge(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType);
    template <int normalX, int normalY, int normalZ>
    void addTemperatureCorner(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType);
};

///////// class AdvectionDiffusionBoundaryConditionInstantiator3D ////////////////////////

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
AdvectionDiffusionBoundaryConditionInstantiator3D<
    T, Descriptor, BoundaryManager>::AdvectionDiffusionBoundaryConditionInstantiator3D()
{ }

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int direction, int orientation>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    PLB_PRECONDITION(domain.x0 == domain.x1 || domain.y0 == domain.y1 || domain.z0 == domain.z1);

    setCompositeDynamics(
        lattice, domain,
        BoundaryManager::template getTemperatureBoundaryDynamics<direction, orientation>(
            new NoDynamics<T, Descriptor>));

    // In case an outflow condition is used, start by instantiating a data processor which copies
    //   all velocity values from the previous lattice cell.
    if (bcType == boundary::neumann) {
        integrateProcessingFunctional(
            new FlatAdiabaticBoundaryFunctional3D<T, Descriptor, direction, orientation>, domain,
            lattice);
    }

    BoxProcessingFunctional3D_L<T, Descriptor> *functional =
        BoundaryManager::template getTemperatureBoundaryProcessor<direction, orientation>(domain);
    if (functional) {
        integrateProcessingFunctional(functional, domain, lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int plane, int normal1, int normal2>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    PLB_PRECONDITION(
        (domain.x0 == domain.x1 && domain.y0 == domain.y1)
        || (domain.x0 == domain.x1 && domain.z0 == domain.z1)
        || (domain.y0 == domain.y1 && domain.z0 == domain.z1));

    setCompositeDynamics(
        lattice, domain,
        BoundaryManager::template getTemperatureEdgeDynamics<plane, normal1, normal2>(
            new NoDynamics<T, Descriptor>));

    BoxProcessingFunctional3D_L<T, Descriptor> *functional =
        BoundaryManager::template getTemperatureEdgeProcessor<plane, normal1, normal2>(domain);
    if (functional) {
        integrateProcessingFunctional(functional, domain, lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int xNormal, int yNormal, int zNormal>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureCorner(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    setCompositeDynamics(
        lattice, Box3D(x, x, y, y, z, z),
        BoundaryManager::template getTemperatureCornerDynamics<xNormal, yNormal, zNormal>(
            new NoDynamics<T, Descriptor>));

    BoxProcessingFunctional3D_L<T, Descriptor> *functional =
        BoundaryManager::template getTemperatureCornerProcessor<xNormal, yNormal, zNormal>(x, y, z);
    if (functional) {
        integrateProcessingFunctional(functional, Box3D(x, x, y, y, z, z), lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary0N(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureBoundary<0, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary0P(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureBoundary<0, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary1N(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureBoundary<1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary1P(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureBoundary<1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary2N(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureBoundary<2, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary2P(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureBoundary<2, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge0NN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<0, -1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge0NP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<0, -1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge0PN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<0, 1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge0PP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<0, 1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge1NN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<1, -1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge1NP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<1, -1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge1PN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<1, 1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge1PP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<1, 1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge2NN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<2, -1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge2NP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<2, -1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge2PN(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<2, 1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge2PP(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<2, 1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerNNN(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureCorner<-1, -1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerNNP(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureCorner<-1, -1, 1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerNPN(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureCorner<-1, 1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerNPP(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureCorner<-1, 1, 1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerPNN(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureCorner<1, -1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerPNP(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureCorner<1, -1, 1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerPPN(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureCorner<1, 1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerPPP(
        plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureCorner<1, 1, 1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int direction, int orientation>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    PLB_PRECONDITION(domain.x0 == domain.x1 || domain.y0 == domain.y1 || domain.z0 == domain.z1);

    setCompositeDynamics(
        lattice, domain,
        BoundaryManager::template getTemperatureBoundaryDynamics<direction, orientation>(
            new NoDynamics<T, Descriptor>));

    // In case an outflow condition is used, start by instantiating a data processor which copies
    //   all velocity values from the previous lattice cell.
    if (bcType == boundary::neumann) {
        integrateProcessingFunctional(
            new FlatAdiabaticBoundaryFunctional3D<T, Descriptor, direction, orientation>, domain,
            lattice);
    }

    BoxProcessingFunctional3D_L<T, Descriptor> *functional =
        BoundaryManager::template getTemperatureBoundaryProcessor<direction, orientation>(domain);

    if (functional) {
        integrateProcessingFunctional(functional, domain, lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int plane, int normal1, int normal2>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    PLB_PRECONDITION(
        (domain.x0 == domain.x1 && domain.y0 == domain.y1)
        || (domain.x0 == domain.x1 && domain.z0 == domain.z1)
        || (domain.y0 == domain.y1 && domain.z0 == domain.z1));

    setCompositeDynamics(
        lattice, domain,
        BoundaryManager::template getTemperatureEdgeDynamics<plane, normal1, normal2>(
            new NoDynamics<T, Descriptor>));

    BoxProcessingFunctional3D_L<T, Descriptor> *functional =
        BoundaryManager::template getTemperatureEdgeProcessor<plane, normal1, normal2>(domain);
    if (functional) {
        integrateProcessingFunctional(functional, domain, lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int xNormal, int yNormal, int zNormal>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureCorner(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType)
{
    setCompositeDynamics(
        lattice, Box3D(x, x, y, y, z, z),
        BoundaryManager::template getTemperatureCornerDynamics<xNormal, yNormal, zNormal>(
            new NoDynamics<T, Descriptor>));

    BoxProcessingFunctional3D_L<T, Descriptor> *functional =
        BoundaryManager::template getTemperatureCornerProcessor<xNormal, yNormal, zNormal>(x, y, z);
    if (functional) {
        integrateProcessingFunctional(functional, Box3D(x, x, y, y, z, z), lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary0N(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureBoundary<0, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary0P(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureBoundary<0, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary1N(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureBoundary<1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary1P(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureBoundary<1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary2N(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureBoundary<2, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary2P(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureBoundary<2, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge0NN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<0, -1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge0NP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<0, -1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge0PN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<0, 1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge0PP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<0, 1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge1NN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<1, -1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge1NP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<1, -1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge1PN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<1, 1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge1PP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<1, 1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge2NN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<2, -1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge2NP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<2, -1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge2PN(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<2, 1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureEdge2PP(
        Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureEdge<2, 1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerNNN(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType)
{
    addTemperatureCorner<-1, -1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerNNP(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType)
{
    addTemperatureCorner<-1, -1, 1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerNPN(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType)
{
    addTemperatureCorner<-1, 1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerNPP(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType)
{
    addTemperatureCorner<-1, 1, 1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerPNN(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType)
{
    addTemperatureCorner<1, -1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerPNP(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType)
{
    addTemperatureCorner<1, -1, 1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerPPN(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType)
{
    addTemperatureCorner<1, 1, -1>(x, y, z, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerPPP(
        plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
        boundary::BcType bcType)
{
    addTemperatureCorner<1, 1, 1>(x, y, z, lattice, bcType);
}

}  // namespace plb

#endif
