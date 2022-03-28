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

#ifndef GUO_ADV_DIFF_OFF_LATTICE_MODEL_3D_HH
#define GUO_ADV_DIFF_OFF_LATTICE_MODEL_3D_HH

#include <algorithm>
#include <cmath>

#include "core/plbProfiler.h"
#include "latticeBoltzmann/advectionDiffusionMomentTemplates.h"
#include "latticeBoltzmann/externalFieldAccess.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/indexTemplates.h"
#include "offLattice/guoAdvDiffOffLatticeModel3D.h"
#include "offLattice/nextNeighbors3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
GuoAdvDiffOffLatticeModel3D<T, Descriptor>::GuoAdvDiffOffLatticeModel3D(
    BoundaryShape3D<T, Array<T, 2> > *shape_, int flowType_) :
    OffLatticeModel3D<T, Array<T, 2> >(shape_, flowType_), secondOrderFlag(true)
{ }

template <typename T, template <typename U> class Descriptor>
GuoAdvDiffOffLatticeModel3D<T, Descriptor>::GuoAdvDiffOffLatticeModel3D(
    GuoAdvDiffOffLatticeModel3D<T, Descriptor> const &rhs) :
    OffLatticeModel3D<T, Array<T, 2> >(rhs), secondOrderFlag(rhs.secondOrderFlag)
{ }

template <typename T, template <typename U> class Descriptor>
GuoAdvDiffOffLatticeModel3D<T, Descriptor> &GuoAdvDiffOffLatticeModel3D<T, Descriptor>::operator=(
    GuoAdvDiffOffLatticeModel3D<T, Descriptor> const &rhs)
{
    OffLatticeModel3D<T, Array<T, 2> >::operator=(rhs);
    return *this;
}

template <typename T, template <typename U> class Descriptor>
GuoAdvDiffOffLatticeModel3D<T, Descriptor> *GuoAdvDiffOffLatticeModel3D<T, Descriptor>::clone()
    const
{
    return new GuoAdvDiffOffLatticeModel3D(*this);
}

template <typename T, template <typename U> class Descriptor>
plint GuoAdvDiffOffLatticeModel3D<T, Descriptor>::getNumNeighbors() const
{
    if (usesSecondOrder()) {
        return 2;
    } else {
        return 1;
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoAdvDiffOffLatticeModel3D<T, Descriptor>::prepareCell(
    Dot3D const &cellLocation, AtomicContainerBlock3D &container)
{
    Dot3D offset = container.getLocation();
    GuoAdvDiffOffLatticeInfo3D *info =
        dynamic_cast<GuoAdvDiffOffLatticeInfo3D *>(container.getData());
    PLB_ASSERT(info);
    if (this->isSolid(cellLocation + offset)) {
        std::vector<std::pair<int, int> > liquidNeighbors;
        std::vector<plint> ids;
        for (int iNeighbor = 0; iNeighbor < NextNeighbor<T>::numNeighbors; ++iNeighbor) {
            int const *c = NextNeighbor<T>::c[iNeighbor];
            Dot3D neighbor(cellLocation.x + c[0], cellLocation.y + c[1], cellLocation.z + c[2]);
            // If the non-fluid node has a fluid neighbor ...
            if (this->isFluid(neighbor + offset)) {
                // ... check how many fluid nodes it has ahead of it ...
                int depth = 1;
                for (int iDepth = 2; iDepth <= getNumNeighbors(); ++iDepth) {
                    Dot3D nextNeighbor(
                        cellLocation.x + iDepth * c[0], cellLocation.y + iDepth * c[1],
                        cellLocation.z + iDepth * c[2]);
                    if (this->isFluid(nextNeighbor + offset)) {
                        depth = iDepth;
                    } else {
                        break;
                    }
                }
                // ... then add this node to the list.
                liquidNeighbors.push_back(std::make_pair(iNeighbor, depth));
                plint iTriangle = -1;
                global::timer("intersect").start();
                Array<T, 3> locatedPoint;
                T distance;
                Array<T, 3> wallNormal;
                Array<T, 2> surfaceData;
                OffBoundary::Type bdType;
                bool ok = this->pointOnSurface(
                    cellLocation + offset, Dot3D(c[0], c[1], c[2]), locatedPoint, distance,
                    wallNormal, surfaceData, bdType, iTriangle);
                global::timer("intersect").stop();
                if (!ok) {
                    global::plbErrors().registerError(
                        "Guo AD off-lattice model could not find an intersection with a triangle.");
                }
                ids.push_back(iTriangle);
                // normals.push_back(wallNormal);
            }
        }
        if (!liquidNeighbors.empty()) {
            info->getDryNodes().push_back(cellLocation);
            info->getDryNodeFluidDirections().push_back(liquidNeighbors);
            info->getDryNodeIds().push_back(ids);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ContainerBlockData *GuoAdvDiffOffLatticeModel3D<T, Descriptor>::generateOffLatticeInfo() const
{
    return new GuoAdvDiffOffLatticeInfo3D;
}

template <typename T, template <typename U> class Descriptor>
void GuoAdvDiffOffLatticeModel3D<T, Descriptor>::boundaryCompletion(
    AtomicBlock3D &nonTypeLattice, AtomicContainerBlock3D &container,
    std::vector<AtomicBlock3D *> const &args)
{
    BlockLattice3D<T, Descriptor> &lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> &>(nonTypeLattice);
    GuoAdvDiffOffLatticeInfo3D *info =
        dynamic_cast<GuoAdvDiffOffLatticeInfo3D *>(container.getData());
    PLB_ASSERT(info);
    std::vector<Dot3D> const &dryNodes = info->getDryNodes();
    std::vector<std::vector<std::pair<int, int> > > const &dryNodeFluidDirections =
        info->getDryNodeFluidDirections();
    std::vector<std::vector<plint> > const &dryNodeIds = info->getDryNodeIds();
    if (dryNodes.size() != dryNodeFluidDirections.size()) {
        global::plbErrors().registerError(
            "Error in the Guo AD off-lattice model boundary completion.");
    }

    Dot3D absoluteOffset = container.getLocation();

    for (pluint iDry = 0; iDry < dryNodes.size(); ++iDry) {
        cellCompletion(
            lattice, dryNodes[iDry], dryNodeFluidDirections[iDry], dryNodeIds[iDry],
            absoluteOffset);
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoAdvDiffOffLatticeModel3D<T, Descriptor>::cellCompletion(
    BlockLattice3D<T, Descriptor> &lattice, Dot3D const &guoNode,
    std::vector<std::pair<int, int> > const &dryNodeFluidDirections,
    std::vector<plint> const &dryNodeIds, Dot3D const &absoluteOffset)
{
    typedef Descriptor<T> D;
    Cell<T, Descriptor> &cell = lattice.get(guoNode.x, guoNode.y, guoNode.z);
    plint numDirections = (plint)dryNodeFluidDirections.size();
    std::vector<T> weights(numDirections);
    std::vector<T> rhoBar_vect(numDirections);
    std::vector<Array<T, Descriptor<T>::d> > jNeq_vect(numDirections);
    T sumWeights = T();
    Array<T, 3> wallNormal;
    for (plint iDirection = 0; iDirection < numDirections; ++iDirection) {
        int iNeighbor = dryNodeFluidDirections[iDirection].first;
        int const *c = NextNeighbor<T>::c[iNeighbor];
        Dot3D fluidDirection(c[0], c[1], c[2]);
        plint dryNodeId = dryNodeIds[iDirection];
        int depth = dryNodeFluidDirections[iDirection].second;

        Array<T, 3> wallNode;
        Array<T, 2> wallData;
        T wallDistance;
        OffBoundary::Type bdType;
        bool ok = this->pointOnSurface(
            guoNode + absoluteOffset, fluidDirection, wallNode, wallDistance, wallNormal, wallData,
            bdType, dryNodeId);
        if (!ok) {
            global::plbErrors().registerError(
                "Guo AD off-lattice model could not find an intersection with a triangle.");
        }
        T invDistanceToNeighbor = NextNeighbor<T>::invD[iNeighbor];
        if (!(wallDistance <= NextNeighbor<T>::d[iNeighbor])) {
            global::plbErrors().registerError(
                "Error treating the geometry in the Guo AD off-lattice model.");
        }
        T delta = (T)1. - wallDistance * invDistanceToNeighbor;
        Array<T, 3> normalFluidDirection(
            (T)fluidDirection.x, (T)fluidDirection.y, (T)fluidDirection.z);
        normalFluidDirection *= invDistanceToNeighbor;
        weights[iDirection] = util::sqr(std::fabs(dot(normalFluidDirection, wallNormal)));
        sumWeights += weights[iDirection];

        computeRhoBarJNeq(
            lattice, guoNode, fluidDirection, depth, wallNode, delta, wallData, bdType, wallNormal,

            rhoBar_vect[iDirection], jNeq_vect[iDirection]);
    }

    T rhoBar = T();
    Array<T, D::d> jNeq;
    jNeq.resetToZero();
    for (plint iDirection = 0; iDirection < numDirections; ++iDirection) {
        rhoBar += rhoBar_vect[iDirection] * weights[iDirection];
        jNeq += jNeq_vect[iDirection] * weights[iDirection];
    }
    rhoBar /= sumWeights;
    jNeq /= sumWeights;

    T dummyJsqr = T();
    Array<T, SymmetricTensor<T, Descriptor>::n> dummyPiNeq;

    Array<T, Descriptor<T>::d> j;
    j.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt));
    j *= Descriptor<T>::fullRho(rhoBar);  // At this point we got jEq.
    j += jNeq;                            // And at this point we add off-equilibrium.

    Dynamics<T, Descriptor> const &dynamics = cell.getDynamics();
    if (this->getPartialReplace()) {
        Cell<T, Descriptor> saveCell(cell);
        dynamics.regularize(cell, rhoBar, j, dummyJsqr, dummyPiNeq);
        for (plint iDirection = 0; iDirection < numDirections; ++iDirection) {
            int iNeighbor = dryNodeFluidDirections[iDirection].first;
            plint iPop = nextNeighborPop<T, Descriptor>(iNeighbor);
            plint oppPop = indexTemplates::opposite<D>(iPop);
            cell[oppPop] = saveCell[oppPop];
        }
    } else {
        dynamics.regularize(cell, rhoBar, j, dummyJsqr, dummyPiNeq);
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoAdvDiffOffLatticeModel3D<T, Descriptor>::computeRhoBarJNeq(
    BlockLattice3D<T, Descriptor> const &lattice, Dot3D const &guoNode, Dot3D const &fluidDirection,
    int depth, [[maybe_unused]] Array<T, 3> const &wallNode, T delta, Array<T, 2> wallData,
    OffBoundary::Type bdType, [[maybe_unused]] Array<T, 3> const &wallNormal, T &rhoBar,
    Array<T, Descriptor<T>::d> &jNeq) const
{
    if (depth > getNumNeighbors()) {
        global::plbErrors().registerError(
            "The number of fluid neighbors is wrong in the Guo AD off-lattice model.");
    }
    T wall_scalar = wallData[0];
    T rhoBar1 = T(), rhoBar2 = T();
    Array<T, Descriptor<T>::d> jNeq1(Array<T, Descriptor<T>::d>::zero()),
        jNeq2(Array<T, Descriptor<T>::d>::zero()), jEqDummy(Array<T, Descriptor<T>::d>::zero());
    Cell<T, Descriptor> const &cell1 = lattice.get(
        guoNode.x + fluidDirection.x, guoNode.y + fluidDirection.y, guoNode.z + fluidDirection.z);

    advectionDiffusionMomentTemplates<T, Descriptor>::get_rhoBar_jEq_jNeq(
        cell1, rhoBar1, jEqDummy, jNeq1);

    if (usesSecondOrder()) {
        Cell<T, Descriptor> const &cell2 = lattice.get(
            guoNode.x + 2 * fluidDirection.x, guoNode.y + 2 * fluidDirection.y,
            guoNode.z + 2 * fluidDirection.z);
        advectionDiffusionMomentTemplates<T, Descriptor>::get_rhoBar_jEq_jNeq(
            cell2, rhoBar2, jEqDummy, jNeq2);
    }
    T wall_rhoBar = Descriptor<T>::rhoBar(wall_scalar);

    if (depth < 2) {
        jNeq = jNeq1;
        if (delta < (T)0.25) {
            rhoBar = wall_rhoBar;
        } else {
            rhoBar = (T)1. / delta * (wall_rhoBar + (delta - (T)1.) * rhoBar1);
        }
    } else {  // depth >= 2
        if (delta < (T)0.75) {
            jNeq = delta * jNeq1 + ((T)1. - delta) * jNeq2;
            rhoBar = wall_rhoBar + (delta - (T)1.) * rhoBar1
                     + ((T)1. - delta) / ((T)1. + delta)
                           * ((T)2. * wall_rhoBar + (delta - (T)1.) * rhoBar2);
        } else {
            jNeq = jNeq1;
            rhoBar = (T)1. / delta * (wall_rhoBar + (delta - (T)1.) * rhoBar1);
        }
    }
    if (bdType == OffBoundary::neumann) {
        if (depth < 2) {
            rhoBar = rhoBar1;
        } else {
            rhoBar = ((T)4.0 * delta * rhoBar1 + ((T)1.0 - (T)2.0 * delta) * rhoBar2)
                     / ((T)1.0 + (T)2.0 * delta);
        }
    } else if (bdType == OffBoundary::flux) {
        T dx = std::sqrt(
            util::sqr((T)fluidDirection.x) + util::sqr((T)fluidDirection.y)
            + util::sqr((T)fluidDirection.z));
        T gradRho = wallData[0];
        rhoBar = rhoBar1 + dx * gradRho;
    } else if (bdType == OffBoundary::isolation) {
        T dx = std::sqrt(
            util::sqr((T)fluidDirection.x) + util::sqr((T)fluidDirection.y)
            + util::sqr((T)fluidDirection.z));
        T rhoBarOut = Descriptor<T>::rhoBar(wallData[0]);
        T kappa = wallData[1];
        rhoBar = (rhoBar1 + kappa * dx * rhoBarOut) / ((T)1 + kappa * dx);
    }
}

}  // namespace plb

#endif  // GUO_ADV_DIFF_OFF_LATTICE_MODEL_3D_HH
