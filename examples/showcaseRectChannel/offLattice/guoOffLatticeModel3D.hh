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

#ifndef GUO_OFF_LATTICE_MODEL_3D_HH
#define GUO_OFF_LATTICE_MODEL_3D_HH

#include <algorithm>
#include <cmath>
#include <vector>

#include "core/plbTimer.h"
#include "core/runTimeDiagnostics.h"
#include "latticeBoltzmann/externalFieldAccess.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/indexTemplates.h"
#include "offLattice/guoOffLatticeModel3D.h"
#include "offLattice/nextNeighbors3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
GuoOffLatticeModel3D<T, Descriptor>::LiquidNeighbor::LiquidNeighbor(
    plint iNeighbor_, plint depth_, plint iTriangle_, Array<T, 3> wallNormal) :
    iNeighbor(iNeighbor_), depth(depth_), iTriangle(iTriangle_)
{
    int const *c = NextNeighbor<T>::c[iNeighbor];
    Array<T, 3> neighborVect(c[0], c[1], c[2]);
    cosAngle = std::fabs(dot(neighborVect, wallNormal)) * NextNeighbor<T>::invD[iNeighbor];
}

template <typename T, template <typename U> class Descriptor>
bool GuoOffLatticeModel3D<T, Descriptor>::LiquidNeighbor::operator<(LiquidNeighbor const &rhs) const
{
    return cosAngle < rhs.cosAngle;
}

/**
 * This class implements the Guo (GZS,2002) boundary condition on a BoundaryShape.
 * The BoundaryShape determines whether the points of the discrete lattice are "inside"
 * or "outside" some geometry.
 *
 * It can handle moving boundaries using the momentum correction of ladd (LADD, 1994).
 * The wall velocity is recovered from SurfaceData stored in BoundaryShape3D<T,SurfaceData>*
 *
 * (GZS, 2002) Z. Guo, C. Zheng, and B. Shi, “An extrapolation method for boundary conditions in
 lattice Boltzmann method,”
 * Physics of Fluids, vol. 14, no. 6, pp. 2007–2010, Jun. 2002, doi: 10.1063/1.1471914.

 * (LADD, 1994) A. J. C. Ladd, “Numerical simulations of particulate suspensions via a discretized
 Boltzmann equation. Part 1. Theoretical foundation,”
 *              Journal of Fluid Mechanics, vol. 271, pp. 285–309, Jul. 1994,
 doi: 10.1017/S0022112094001771.
 *
 * @tparam T
 * @tparam Descriptor
 */
template <typename T, template <typename U> class Descriptor>
GuoOffLatticeModel3D<T, Descriptor>::GuoOffLatticeModel3D(
    BoundaryShape3D<T, Array<T, 3> > *shape_, int flowType_, bool useAllDirections_) :
    OffLatticeModel3D<T, Array<T, 3> >(shape_, flowType_), useAllDirections(useAllDirections_)
{ }

template <typename T, template <typename U> class Descriptor>
GuoOffLatticeModel3D<T, Descriptor> *GuoOffLatticeModel3D<T, Descriptor>::clone() const
{
    return new GuoOffLatticeModel3D(*this);
}

template <typename T, template <typename U> class Descriptor>
plint GuoOffLatticeModel3D<T, Descriptor>::getNumNeighbors() const
{
    return 2;
}

template <typename T, template <typename U> class Descriptor>
bool GuoOffLatticeModel3D<T, Descriptor>::isExtrapolated() const
{
    // Guo is a completion scheme for a layer of cells on the
    // "solid" side of the boundary.
    return true;
}

template <typename T, template <typename U> class Descriptor>
void GuoOffLatticeModel3D<T, Descriptor>::prepareCell(
    Dot3D const &cellLocation, AtomicContainerBlock3D &container)
{
    Dot3D offset = container.getLocation();
    GuoOffLatticeInfo3D *info = dynamic_cast<GuoOffLatticeInfo3D *>(container.getData());
    PLB_ASSERT(info);
    std::vector<LiquidNeighbor> liquidNeighbors;
    if (this->isSolid(cellLocation + offset)) {
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
                plint iTriangle = -1;
                global::timer("intersect").start();
                Array<T, 3> locatedPoint;
                T distance;
                Array<T, 3> wallNormal;
                Array<T, 3> surfaceData;
                OffBoundary::Type bdType;
                bool ok = this->pointOnSurface(
                    cellLocation + offset, Dot3D(c[0], c[1], c[2]), locatedPoint, distance,
                    wallNormal, surfaceData, bdType, iTriangle);
                // In the following, the importance of directions is sorted wrt. how well they
                //   are aligned with the wall normal. It is better to take the continuous normal,
                //   because it is not sensitive to the choice of the triangle when we shoot at
                //   an edge.
                // wallNormal = this->computeContinuousNormal(locatedPoint, iTriangle);
                global::timer("intersect").stop();
                if (!ok) {
                    global::plbErrors().registerError(
                        "Guo off-lattice model could not find an intersection with a triangle.");
                }
                // ... then add this node to the list.
                liquidNeighbors.push_back(LiquidNeighbor(iNeighbor, depth, iTriangle, wallNormal));
            }
        }
        if (!liquidNeighbors.empty()) {
            info->getDryNodes().push_back(cellLocation);
            std::sort(liquidNeighbors.begin(), liquidNeighbors.end());
            std::vector<std::pair<int, int> > neighborDepthPairs;
            std::vector<plint> ids;
            if (useAllDirections) {
                for (pluint i = 0; i < liquidNeighbors.size(); ++i) {
                    neighborDepthPairs.push_back(
                        std::make_pair(liquidNeighbors[i].iNeighbor, liquidNeighbors[i].depth));
                    ids.push_back(liquidNeighbors[i].iTriangle);
                }
            } else {
                plint i = liquidNeighbors.size() - 1;
                neighborDepthPairs.push_back(
                    std::make_pair(liquidNeighbors[i].iNeighbor, liquidNeighbors[i].depth));
                ids.push_back(liquidNeighbors[i].iTriangle);
            }
            info->getDryNodeFluidDirections().push_back(neighborDepthPairs);
            info->getDryNodeIds().push_back(ids);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ContainerBlockData *GuoOffLatticeModel3D<T, Descriptor>::generateOffLatticeInfo() const
{
    return new GuoOffLatticeInfo3D;
}

template <typename T, template <typename U> class Descriptor>
Array<T, 3> GuoOffLatticeModel3D<T, Descriptor>::getLocalForce(
    AtomicContainerBlock3D &container) const
{
    GuoOffLatticeInfo3D *info = dynamic_cast<GuoOffLatticeInfo3D *>(container.getData());
    PLB_ASSERT(info);
    return info->getLocalForce();
}

template <typename T, template <typename U> class Descriptor>
void GuoOffLatticeModel3D<T, Descriptor>::boundaryCompletion(
    AtomicBlock3D &nonTypeLattice, AtomicContainerBlock3D &container,
    std::vector<AtomicBlock3D *> const &args)
{
    BlockLattice3D<T, Descriptor> &lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> &>(nonTypeLattice);
    GuoOffLatticeInfo3D *info = dynamic_cast<GuoOffLatticeInfo3D *>(container.getData());
    PLB_ASSERT(info);
    std::vector<Dot3D> const &dryNodes = info->getDryNodes();
    std::vector<std::vector<std::pair<int, int> > > const &dryNodeFluidDirections =
        info->getDryNodeFluidDirections();
    std::vector<std::vector<plint> > const &dryNodeIds = info->getDryNodeIds();
    if (dryNodes.size() != dryNodeFluidDirections.size()) {
        global::plbErrors().registerError(
            "Error in the Guo off-lattice model boundary completion.");
    }

    Dot3D absoluteOffset = container.getLocation();

    Array<T, 3> &localForce = info->getLocalForce();
    localForce.resetToZero();
    for (pluint iDry = 0; iDry < dryNodes.size(); ++iDry) {
        cellCompletion(
            lattice, dryNodes[iDry], dryNodeFluidDirections[iDry], dryNodeIds[iDry], absoluteOffset,
            localForce, args);
    }
}

template <typename T, template <typename U> class Descriptor>
class GuoAlgorithm3D {
public:
    typedef Descriptor<T> D;
    GuoAlgorithm3D(
        OffLatticeModel3D<T, Array<T, 3> > &model_, BlockLattice3D<T, Descriptor> &lattice_,
        Dot3D const &guoNode_, std::vector<std::pair<int, int> > const &dryNodeFluidDirections_,
        std::vector<plint> const &dryNodeIds_, Dot3D const &absoluteOffset_,
        Array<T, 3> &localForce_, std::vector<AtomicBlock3D *> const &args_, bool computeStat_,
        bool secondOrder_);
    virtual ~GuoAlgorithm3D() { }
    virtual bool computeNeighborData();
    void finalize();

    virtual void extrapolateVariables(
        Dot3D const &fluidDirection, int depth, Array<T, 3> const &wallNode, T delta,
        Array<T, 3> const &wall_vel, OffBoundary::Type bdType, Array<T, 3> const &wallNormal,
        plint triangleId, plint iDirection) = 0;
    virtual void reduceVariables(T sumWeights) = 0;
    virtual void complete() = 0;

protected:
    OffLatticeModel3D<T, Array<T, 3> > &model;
    BlockLattice3D<T, Descriptor> &lattice;
    Dot3D const &guoNode;
    Cell<T, Descriptor> &cell;
    std::vector<std::pair<int, int> > const &dryNodeFluidDirections;
    std::vector<plint> const &dryNodeIds;
    Dot3D absoluteOffset;
    Array<T, 3> &localForce;
    std::vector<AtomicBlock3D *> const &args;

    plint numDirections;
    std::vector<T> weights;
    std::vector<T> rhoBarVect;
    std::vector<Array<T, Descriptor<T>::d> > jVect;

    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    bool computeStat, secondOrder;
};

template <typename T, template <typename U> class Descriptor>
GuoAlgorithm3D<T, Descriptor>::GuoAlgorithm3D(
    OffLatticeModel3D<T, Array<T, 3> > &model_, BlockLattice3D<T, Descriptor> &lattice_,
    Dot3D const &guoNode_, std::vector<std::pair<int, int> > const &dryNodeFluidDirections_,
    std::vector<plint> const &dryNodeIds_, Dot3D const &absoluteOffset_, Array<T, 3> &localForce_,
    std::vector<AtomicBlock3D *> const &args_, bool computeStat_, bool secondOrder_) :
    model(model_),
    lattice(lattice_),
    guoNode(guoNode_),
    cell(lattice.get(guoNode.x, guoNode.y, guoNode.z)),
    dryNodeFluidDirections(dryNodeFluidDirections_),
    dryNodeIds(dryNodeIds_),
    absoluteOffset(absoluteOffset_),
    localForce(localForce_),
    args(args_),
    computeStat(computeStat_),
    secondOrder(secondOrder_)
{
    numDirections = (plint)dryNodeFluidDirections.size();
    weights.resize(numDirections);
    rhoBarVect.resize(numDirections);
    jVect.resize(numDirections);
}

template <typename T, template <typename U> class Descriptor>
bool GuoAlgorithm3D<T, Descriptor>::computeNeighborData()
{
    T sumWeights = T();
    Array<T, 3> wallNormal;
    for (plint iDirection = 0; iDirection < numDirections; ++iDirection) {
        int iNeighbor = dryNodeFluidDirections[iDirection].first;
        int const *c = NextNeighbor<T>::c[iNeighbor];
        Dot3D fluidDirection(c[0], c[1], c[2]);
        plint dryNodeId = dryNodeIds[iDirection];
        int depth = dryNodeFluidDirections[iDirection].second;

        Array<T, 3> wallNode, wall_vel;
        T wallDistance;
        OffBoundary::Type bdType;
        bool ok = this->model.pointOnSurface(
            guoNode + absoluteOffset, fluidDirection, wallNode, wallDistance, wallNormal, wall_vel,
            bdType, dryNodeId);
        if (!ok) {
            global::plbErrors().registerError(
                "Guo off-lattice model could not find an intersection with a triangle.");
        }
        if (!(bdType == OffBoundary::dirichlet || bdType == OffBoundary::neumann
              || bdType == OffBoundary::freeSlip || bdType == OffBoundary::constRhoInlet
              || bdType == OffBoundary::densityNeumann))
        {
            return false;
        }
        if (bdType == OffBoundary::dirichlet) {
            for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
                // Use the formula uLB = uP - 1/2 g. If there is no external force,
                //   the force term automatically evaluates to zero.
                wall_vel[iD] -= (T)0.5 * getExternalForceComponent(cell, iD);
            }
        }
        T invDistanceToNeighbor = NextNeighbor<T>::invD[iNeighbor];
        if (!(wallDistance <= NextNeighbor<T>::d[iNeighbor])) {
            global::plbErrors().registerError(
                "Error treating the geometry in the Guo off-lattice model.");
        }
        T delta = (T)1. - wallDistance * invDistanceToNeighbor;
        Array<T, 3> normalFluidDirection(
            (T)fluidDirection.x, (T)fluidDirection.y, (T)fluidDirection.z);
        normalFluidDirection *= invDistanceToNeighbor;
        weights[iDirection] = std::fabs(dot(normalFluidDirection, wallNormal));
        sumWeights += weights[iDirection];
        this->extrapolateVariables(
            fluidDirection, depth, wallNode, delta, wall_vel, bdType, wallNormal, dryNodeId,
            iDirection);
    }
    this->reduceVariables(sumWeights);

    return true;
}

template <typename T, template <typename U> class Descriptor>
void GuoAlgorithm3D<T, Descriptor>::finalize()
{
    Array<T, D::d> deltaJ;
    deltaJ.resetToZero();
    if (computeStat) {
        for (plint iDirection = 0; iDirection < numDirections; ++iDirection) {
            int iNeighbor = dryNodeFluidDirections[iDirection].first;
            int iPop = nextNeighborPop<T, Descriptor>(iNeighbor);
            if (iPop >= 0) {
                plint oppPop = indexTemplates::opposite<D>(iPop);
                deltaJ[0] += D::c[oppPop][0] * cell[oppPop];
                deltaJ[1] += D::c[oppPop][1] * cell[oppPop];
                deltaJ[2] += D::c[oppPop][2] * cell[oppPop];
            }
        }
    }

    this->complete();

    if (computeStat) {
        Cell<T, Descriptor> collidedCell(cell);
        BlockStatistics statsCopy(lattice.getInternalStatistics());
        collidedCell.collide(statsCopy);

        for (plint iDirection = 0; iDirection < numDirections; ++iDirection) {
            int iNeighbor = dryNodeFluidDirections[iDirection].first;
            plint iPop = nextNeighborPop<T, Descriptor>(iNeighbor);
            if (iPop >= 0) {
                deltaJ[0] -= D::c[iPop][0] * collidedCell[iPop];
                deltaJ[1] -= D::c[iPop][1] * collidedCell[iPop];
                deltaJ[2] -= D::c[iPop][2] * collidedCell[iPop];
            }
        }
        // Don't divide by rho. Here we just divide by rho0=1. Remember that,
        //   while j is conserved along a channel, rho is not due to the
        //   pressure drop. Dividing by rho instead of rho0 would yield a
        //   larger result upstream than downstream, independent of whether
        //   the velocity is u or j.
        localForce += deltaJ;
    }
}

template <typename T, template <typename U> class Descriptor>
class GuoPiNeqAlgorithm3D : public GuoAlgorithm3D<T, Descriptor> {
public:
    typedef Descriptor<T> D;
    GuoPiNeqAlgorithm3D(
        OffLatticeModel3D<T, Array<T, 3> > &model_, BlockLattice3D<T, Descriptor> &lattice_,
        Dot3D const &guoNode_, std::vector<std::pair<int, int> > const &dryNodeFluidDirections_,
        std::vector<plint> const &dryNodeIds_, Dot3D const &absoluteOffset_,
        Array<T, 3> &localForce_, std::vector<AtomicBlock3D *> const &args_, bool computeStat_,
        bool secondOrder_);
    virtual void extrapolateVariables(
        Dot3D const &fluidDirection, int depth, Array<T, 3> const &wallNode, T delta,
        Array<T, 3> const &wall_vel, OffBoundary::Type bdType, Array<T, 3> const &wallNormal,
        plint triangleId, plint iDirection);
    virtual void reduceVariables(T sumWeights);
    virtual void complete();

private:
    std::vector<Array<T, SymmetricTensor<T, Descriptor>::n> > PiNeqVect;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
};

template <typename T, template <typename U> class Descriptor>
GuoPiNeqAlgorithm3D<T, Descriptor>::GuoPiNeqAlgorithm3D(
    OffLatticeModel3D<T, Array<T, 3> > &model_, BlockLattice3D<T, Descriptor> &lattice_,
    Dot3D const &guoNode_, std::vector<std::pair<int, int> > const &dryNodeFluidDirections_,
    std::vector<plint> const &dryNodeIds_, Dot3D const &absoluteOffset_, Array<T, 3> &localForce_,
    std::vector<AtomicBlock3D *> const &args_, bool computeStat_, bool secondOrder_) :
    GuoAlgorithm3D<T, Descriptor>(
        model_, lattice_, guoNode_, dryNodeFluidDirections_, dryNodeIds_, absoluteOffset_,
        localForce_, args_, computeStat_, secondOrder_)
{
    PiNeqVect.resize(this->numDirections);
    PiNeq.resetToZero();
}

template <typename T, template <typename U> class Descriptor>
void GuoPiNeqAlgorithm3D<T, Descriptor>::extrapolateVariables(
    Dot3D const &fluidDirection, int depth, Array<T, 3> const &wallNode, T delta,
    Array<T, 3> const &wall_vel, OffBoundary::Type bdType, Array<T, 3> const &wallNormal,
    plint triangleId, plint iDirection)
{
    if (!this->secondOrder) {
        depth = 1;
    }
    T rhoBar1;
    Array<T, Descriptor<T>::d> j1, j2;
    Cell<T, Descriptor> const &cell1 = this->lattice.get(
        this->guoNode.x + fluidDirection.x, this->guoNode.y + fluidDirection.y,
        this->guoNode.z + fluidDirection.z);
    cell1.getDynamics().computeRhoBarJPiNeq(cell1, rhoBar1, j1, this->PiNeqVect[iDirection]);
    if (this->args.empty()) {
        Cell<T, Descriptor> const &cell2 = this->lattice.get(
            this->guoNode.x + 2 * fluidDirection.x, this->guoNode.y + 2 * fluidDirection.y,
            this->guoNode.z + 2 * fluidDirection.z);

        T tmpRhoBar;
        cell2.getDynamics().computeRhoBarJ(cell2, tmpRhoBar, j2);
    } else {
        if ((plint)this->args.size() == 1) {
            NTensorField3D<T> const *macroField =
                dynamic_cast<NTensorField3D<T> const *>(this->args[0]);
            PLB_ASSERT(macroField);
            // 1 Variable for rhoBar, 3 variables for j.
            PLB_ASSERT(macroField->getNdim() == 4);
            Dot3D offset = computeRelativeDisplacement(this->lattice, *macroField);
            T const *macroscopic = macroField->get(
                this->guoNode.x + 2 * fluidDirection.x + offset.x,
                this->guoNode.y + 2 * fluidDirection.y + offset.y,
                this->guoNode.z + 2 * fluidDirection.z + offset.z);
            j2.from_cArray(macroscopic + 1);
        } else if ((plint)this->args.size() == 2) {
            // 1 field for rhoBar, 1 field for j.
#ifdef PLB_DEBUG
            ScalarField3D<T> const *rhoBarField =
                dynamic_cast<ScalarField3D<T> const *>(this->args[0]);
#endif
            TensorField3D<T, 3> const *jField =
                dynamic_cast<TensorField3D<T, 3> const *>(this->args[1]);
            PLB_ASSERT(rhoBarField);
            PLB_ASSERT(jField);
            Dot3D offset = computeRelativeDisplacement(this->lattice, *jField);
            j2 = jField->get(
                this->guoNode.x + 2 * fluidDirection.x + offset.x,
                this->guoNode.y + 2 * fluidDirection.y + offset.y,
                this->guoNode.z + 2 * fluidDirection.z + offset.z);
        } else {
            PLB_ASSERT(false);  // Not implemented for 3 args.
        }
    }

    if (bdType == OffBoundary::constRhoInlet || bdType == OffBoundary::densityNeumann) {
        this->rhoBarVect[iDirection] = Descriptor<T>::rhoBar(wall_vel[0]);
    } else {
        this->rhoBarVect[iDirection] = rhoBar1;
    }
    Array<T, 3> wall_j(
        this->model.velIsJ() ? wall_vel
                             : (Descriptor<T>::fullRho(this->rhoBarVect[iDirection]) * wall_vel));
    if (depth < 2) {
        if (delta < (T)0.25) {
            this->jVect[iDirection] = wall_j;
        } else {
            // j = (T)1./delta * (wall_j+(delta-(T)1.)*j1);
            this->jVect[iDirection] = wall_j;  // Temporary fix for complex geometries.
        }
    } else {                    // depth >= 2           d=1   d=0.5   d=0
        if (delta < (T)0.75) {  // x---|=========o-------------o
            this->jVect[iDirection] =
                wall_j + (delta - (T)1.) * j1
                + ((T)1. - delta) / ((T)1. + delta) * ((T)2. * wall_j + (delta - (T)1.) * j2);
        }  //       d=1   d=0.5   d=0
        else
        {  // x===|---------o-------------o
            this->jVect[iDirection] = (T)1. / delta * (wall_j + (delta - (T)1.) * j1);
        }
    }
    if (bdType == OffBoundary::neumann) {
        this->jVect[iDirection] = j1;
    } else if (bdType == OffBoundary::densityNeumann) {
        this->jVect[iDirection] = dot(j1, wallNormal) * wallNormal;
    } else if (bdType == OffBoundary::freeSlip) {
        Array<T, 3> continuousNormal = this->model.computeContinuousNormal(wallNode, triangleId);
        this->jVect[iDirection] = j1 - dot(j1, continuousNormal) * continuousNormal;
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoPiNeqAlgorithm3D<T, Descriptor>::reduceVariables(T sumWeights)
{
    this->rhoBar = T();
    this->j.resetToZero();
    PiNeq.resetToZero();
    for (plint iDirection = 0; iDirection < this->numDirections; ++iDirection) {
        this->rhoBar += this->rhoBarVect[iDirection] * this->weights[iDirection];
        this->j += this->jVect[iDirection] * this->weights[iDirection];
        PiNeq += PiNeqVect[iDirection] * this->weights[iDirection];
    }
    this->rhoBar /= sumWeights;
    this->j /= sumWeights;
    PiNeq /= sumWeights;
}

template <typename T, template <typename U> class Descriptor>
void GuoPiNeqAlgorithm3D<T, Descriptor>::complete()
{
    Dynamics<T, Descriptor> const &dynamics = this->cell.getDynamics();
    T jSqr = normSqr(this->j);
    if (this->model.getPartialReplace()) {
        Cell<T, Descriptor> saveCell(this->cell);
        dynamics.regularize(this->cell, this->rhoBar, this->j, jSqr, PiNeq);
        for (plint iDirection = 0; iDirection < this->numDirections; ++iDirection) {
            int iNeighbor = this->dryNodeFluidDirections[iDirection].first;
            plint iPop = nextNeighborPop<T, Descriptor>(iNeighbor);
            if (iPop >= 0) {
                plint oppPop = indexTemplates::opposite<D>(iPop);
                this->cell[oppPop] = saveCell[oppPop];
            }
        }
    } else {
        dynamics.regularize(this->cell, this->rhoBar, this->j, jSqr, PiNeq);
    }
}

template <typename T, template <typename U> class Descriptor>
class GuoOffPopAlgorithm3D : public GuoAlgorithm3D<T, Descriptor> {
public:
    typedef Descriptor<T> D;
    GuoOffPopAlgorithm3D(
        OffLatticeModel3D<T, Array<T, 3> > &model_, BlockLattice3D<T, Descriptor> &lattice_,
        Dot3D const &guoNode_, std::vector<std::pair<int, int> > const &dryNodeFluidDirections_,
        std::vector<plint> const &dryNodeIds_, Dot3D const &absoluteOffset_,
        Array<T, 3> &localForce_, std::vector<AtomicBlock3D *> const &args_, bool computeStat_,
        bool secondOrder_);
    virtual void extrapolateVariables(
        Dot3D const &fluidDirection, int depth, Array<T, 3> const &wallNode, T delta,
        Array<T, 3> const &wall_vel, OffBoundary::Type bdType, Array<T, 3> const &wallNormal,
        plint triangleId, plint iDirection);
    virtual void reduceVariables(T sumWeights);
    virtual void complete();

private:
    std::vector<Array<T, Descriptor<T>::q> > fNeqVect;
    Array<T, Descriptor<T>::q> fNeq;
};

template <typename T, template <typename U> class Descriptor>
GuoOffPopAlgorithm3D<T, Descriptor>::GuoOffPopAlgorithm3D(
    OffLatticeModel3D<T, Array<T, 3> > &model_, BlockLattice3D<T, Descriptor> &lattice_,
    Dot3D const &guoNode_, std::vector<std::pair<int, int> > const &dryNodeFluidDirections_,
    std::vector<plint> const &dryNodeIds_, Dot3D const &absoluteOffset_, Array<T, 3> &localForce_,
    std::vector<AtomicBlock3D *> const &args_, bool computeStat_, bool secondOrder_) :
    GuoAlgorithm3D<T, Descriptor>(
        model_, lattice_, guoNode_, dryNodeFluidDirections_, dryNodeIds_, absoluteOffset_,
        localForce_, args_, computeStat_, secondOrder_)
{
    fNeqVect.resize(this->numDirections);
    fNeq.resetToZero();
}

template <typename T, template <typename U> class Descriptor>
void GuoOffPopAlgorithm3D<T, Descriptor>::extrapolateVariables(
    Dot3D const &fluidDirection, int depth, Array<T, 3> const &wallNode, T delta,
    Array<T, 3> const &wall_vel, OffBoundary::Type bdType, Array<T, 3> const &wallNormal,
    plint triangleId, plint iDirection)
{
    T rhoBar;
    Array<T, 3> j;
    Cell<T, Descriptor> const &cell =
        this->lattice.get(this->guoNode.x, this->guoNode.y, this->guoNode.z);
    cell.getDynamics().computeRhoBarJ(cell, rhoBar, j);

    if (!this->secondOrder) {
        depth = 1;
    }
    T rhoBar1;
    Array<T, Descriptor<T>::d> j1, j2;
    Array<T, Descriptor<T>::q> fNeq1;
    Cell<T, Descriptor> const &cell1 = this->lattice.get(
        this->guoNode.x + fluidDirection.x, this->guoNode.y + fluidDirection.y,
        this->guoNode.z + fluidDirection.z);
    cell1.getDynamics().computeRhoBarJ(cell1, rhoBar1, j1);

    if (this->args.empty()) {
        Cell<T, Descriptor> const &cell2 = this->lattice.get(
            this->guoNode.x + 2 * fluidDirection.x, this->guoNode.y + 2 * fluidDirection.y,
            this->guoNode.z + 2 * fluidDirection.z);

        T tmpRhoBar;
        cell2.getDynamics().computeRhoBarJ(cell2, tmpRhoBar, j2);
    } else {
        if ((plint)this->args.size() == 1) {
            NTensorField3D<T> const *macroField =
                dynamic_cast<NTensorField3D<T> const *>(this->args[0]);
            PLB_ASSERT(macroField);
            // 1 Variable for rhoBar, 3 variables for j.
            PLB_ASSERT(macroField->getNdim() == 4);
            Dot3D offset = computeRelativeDisplacement(this->lattice, *macroField);
            T const *macroscopic = macroField->get(
                this->guoNode.x + 2 * fluidDirection.x + offset.x,
                this->guoNode.y + 2 * fluidDirection.y + offset.y,
                this->guoNode.z + 2 * fluidDirection.z + offset.z);
            j2.from_cArray(macroscopic + 1);
        } else if ((plint)this->args.size() == 2) {
            // 1 field for rhoBar, 1 field for j.
#ifdef PLB_DEBUG
            ScalarField3D<T> const *rhoBarField =
                dynamic_cast<ScalarField3D<T> const *>(this->args[0]);
#endif
            TensorField3D<T, 3> const *jField =
                dynamic_cast<TensorField3D<T, 3> const *>(this->args[1]);
            PLB_ASSERT(rhoBarField);
            PLB_ASSERT(jField);
            Dot3D offset = computeRelativeDisplacement(this->lattice, *jField);
            j2 = jField->get(
                this->guoNode.x + 2 * fluidDirection.x + offset.x,
                this->guoNode.y + 2 * fluidDirection.y + offset.y,
                this->guoNode.z + 2 * fluidDirection.z + offset.z);
        } else {
            PLB_ASSERT(false);  // Not implemented for 3 args.
        }
    }

    T j1sqr = normSqr(j1);
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        fNeq1[iPop] =
            cell1[iPop] - cell1.getDynamics().computeEquilibrium(iPop, rhoBar1, j1, j1sqr);
    }

    if (bdType == OffBoundary::constRhoInlet || bdType == OffBoundary::densityNeumann) {
        this->rhoBarVect[iDirection] = Descriptor<T>::rhoBar(wall_vel[0]);
    } else {
        this->rhoBarVect[iDirection] = rhoBar1;
    }
    Array<T, 3> wall_j(
        this->model.velIsJ() ? wall_vel
                             : (Descriptor<T>::fullRho(this->rhoBarVect[iDirection]) * wall_vel));
    Array<T, 3> jNeumann;
    this->fNeqVect[iDirection] = fNeq1;
    jNeumann = j1;
    if (depth < 2) {
        if (delta < (T)0.25) {
            this->jVect[iDirection] = wall_j;
        } else {
            // this->jVect[iDirection] = (T)1./delta * (wall_j+(delta-(T)1.)*j1);
            this->jVect[iDirection] = wall_j;  // Temporary fix for complex geometries.
        }
    } else {  // depth >= 2
        if (delta < (T)0.75) {
            this->jVect[iDirection] =
                wall_j + (delta - (T)1.) * j1
                + ((T)1. - delta) / ((T)1. + delta) * ((T)2. * wall_j + (delta - (T)1.) * j2);
        } else {
            this->jVect[iDirection] = (T)1. / delta * (wall_j + (delta - (T)1.) * j1);
        }
    }
    if (bdType == OffBoundary::neumann) {
        this->jVect[iDirection] = jNeumann;
    } else if (bdType == OffBoundary::densityNeumann) {
        this->jVect[iDirection] = dot(jNeumann, wallNormal) * wallNormal;
    } else if (bdType == OffBoundary::freeSlip) {
        Array<T, 3> continuousNormal = this->model.computeContinuousNormal(wallNode, triangleId);
        this->jVect[iDirection] = jNeumann - dot(jNeumann, continuousNormal) * continuousNormal;
    }
    rhoBar1 = rhoBar;
}

template <typename T, template <typename U> class Descriptor>
void GuoOffPopAlgorithm3D<T, Descriptor>::reduceVariables(T sumWeights)
{
    this->rhoBar = T();
    this->j.resetToZero();
    fNeq.resetToZero();
    for (plint iDirection = 0; iDirection < this->numDirections; ++iDirection) {
        this->rhoBar += this->rhoBarVect[iDirection] * this->weights[iDirection];
        this->j += this->jVect[iDirection] * this->weights[iDirection];
        fNeq += fNeqVect[iDirection] * this->weights[iDirection];
    }
    this->rhoBar /= sumWeights;
    this->j /= sumWeights;
    fNeq /= sumWeights;
}

template <typename T, template <typename U> class Descriptor>
void GuoOffPopAlgorithm3D<T, Descriptor>::complete()
{
    T jSqr = normSqr(this->j);
    if (this->model.getPartialReplace()) {
        for (plint iDirection = 0; iDirection < this->numDirections; ++iDirection) {
            int iNeighbor = this->dryNodeFluidDirections[iDirection].first;
            plint iPop = nextNeighborPop<T, Descriptor>(iNeighbor);
            if (iPop >= 0) {
                this->cell[iPop] =
                    this->cell.computeEquilibrium(iPop, this->rhoBar, this->j, jSqr) + fNeq[iPop];
            }
        }
    } else {
        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            this->cell[iPop] =
                this->cell.computeEquilibrium(iPop, this->rhoBar, this->j, jSqr) + fNeq[iPop];
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoOffLatticeModel3D<T, Descriptor>::cellCompletion(
    BlockLattice3D<T, Descriptor> &lattice, Dot3D const &guoNode,
    std::vector<std::pair<int, int> > const &dryNodeFluidDirections,
    std::vector<plint> const &dryNodeIds, Dot3D const &absoluteOffset, Array<T, 3> &localForce,
    std::vector<AtomicBlock3D *> const &args)
{
    GuoAlgorithm3D<T, Descriptor> *algorithm = 0;
    if (this->usesRegularizedModel()) {
        algorithm = new GuoPiNeqAlgorithm3D<T, Descriptor>(
            *this, lattice, guoNode, dryNodeFluidDirections, dryNodeIds, absoluteOffset, localForce,
            args, this->computesStat(), this->usesSecondOrder());
    } else {
        algorithm = new GuoOffPopAlgorithm3D<T, Descriptor>(
            *this, lattice, guoNode, dryNodeFluidDirections, dryNodeIds, absoluteOffset, localForce,
            args, this->computesStat(), this->usesSecondOrder());
    }
    bool ok = algorithm->computeNeighborData();
    if (!ok) {
        global::plbErrors().registerError(
            "Error treating the geometry in the Guo off-lattice model.");
    }
    algorithm->finalize();
    delete algorithm;
}

template <typename T, template <typename U> class Descriptor>
class GuoDefineVelocityAlgorithm3D : public GuoAlgorithm3D<T, Descriptor> {
public:
    typedef Descriptor<T> D;
    GuoDefineVelocityAlgorithm3D(
        OffLatticeModel3D<T, Array<T, 3> > &model_, BlockLattice3D<T, Descriptor> &lattice_,
        Dot3D const &guoNode_, std::vector<std::pair<int, int> > const &dryNodeFluidDirections_,
        std::vector<plint> const &dryNodeIds_, Dot3D const &absoluteOffset_,
        Array<T, 3> &localForce_, std::vector<AtomicBlock3D *> const &args_, bool computeStat_,
        bool secondOrder_);
    virtual void extrapolateVariables(
        Dot3D const &fluidDirection, int depth, Array<T, 3> const &wallNode, T delta,
        Array<T, 3> const &wall_vel, OffBoundary::Type bdType, Array<T, 3> const &wallNormal,
        plint triangleId, plint iDirection);
    virtual void reduceVariables(T sumWeights);
    virtual void complete();
};

template <typename T, template <typename U> class Descriptor>
GuoDefineVelocityAlgorithm3D<T, Descriptor>::GuoDefineVelocityAlgorithm3D(
    OffLatticeModel3D<T, Array<T, 3> > &model_, BlockLattice3D<T, Descriptor> &lattice_,
    Dot3D const &guoNode_, std::vector<std::pair<int, int> > const &dryNodeFluidDirections_,
    std::vector<plint> const &dryNodeIds_, Dot3D const &absoluteOffset_, Array<T, 3> &localForce_,
    std::vector<AtomicBlock3D *> const &args_, bool computeStat_, bool secondOrder_) :
    GuoAlgorithm3D<T, Descriptor>(
        model_, lattice_, guoNode_, dryNodeFluidDirections_, dryNodeIds_, absoluteOffset_,
        localForce_, args_, computeStat_, secondOrder_)
{ }

template <typename T, template <typename U> class Descriptor>
void GuoDefineVelocityAlgorithm3D<T, Descriptor>::extrapolateVariables(
    Dot3D const &fluidDirection, int depth, Array<T, 3> const &wallNode, T delta,
    Array<T, 3> const &wall_vel, OffBoundary::Type bdType, Array<T, 3> const &wallNormal,
    plint triangleId, plint iDirection)
{
    if (!this->secondOrder) {
        depth = 1;
    }
    T rhoBar1;
    Array<T, Descriptor<T>::d> j1, j2;
    Cell<T, Descriptor> const &cell1 = this->lattice.get(
        this->guoNode.x + fluidDirection.x, this->guoNode.y + fluidDirection.y,
        this->guoNode.z + fluidDirection.z);
    cell1.getDynamics().computeRhoBarJ(cell1, rhoBar1, j1);
    if (this->args.empty()) {
        Cell<T, Descriptor> const &cell2 = this->lattice.get(
            this->guoNode.x + 2 * fluidDirection.x, this->guoNode.y + 2 * fluidDirection.y,
            this->guoNode.z + 2 * fluidDirection.z);

        T tmpRhoBar;
        cell2.getDynamics().computeRhoBarJ(cell2, tmpRhoBar, j2);
    } else {
        if ((plint)this->args.size() == 1) {
            NTensorField3D<T> const *macroField =
                dynamic_cast<NTensorField3D<T> const *>(this->args[0]);
            PLB_ASSERT(macroField);
            // 1 Variable for rhoBar, 3 variables for j.
            PLB_ASSERT(macroField->getNdim() == 4);
            Dot3D offset = computeRelativeDisplacement(this->lattice, *macroField);
            T const *macroscopic = macroField->get(
                this->guoNode.x + 2 * fluidDirection.x + offset.x,
                this->guoNode.y + 2 * fluidDirection.y + offset.y,
                this->guoNode.z + 2 * fluidDirection.z + offset.z);
            j2.from_cArray(macroscopic + 1);
        } else if ((plint)this->args.size() == 2) {
            // 1 field for rhoBar, 1 field for j.
            ScalarField3D<T> const *rhoBarField =
                dynamic_cast<ScalarField3D<T> const *>(this->args[0]);
            TensorField3D<T, 3> const *jField =
                dynamic_cast<TensorField3D<T, 3> const *>(this->args[1]);
            PLB_ASSERT(rhoBarField);
            PLB_ASSERT(jField);
            Dot3D offset = computeRelativeDisplacement(this->lattice, *jField);
            j2 = jField->get(
                this->guoNode.x + 2 * fluidDirection.x + offset.x,
                this->guoNode.y + 2 * fluidDirection.y + offset.y,
                this->guoNode.z + 2 * fluidDirection.z + offset.z);
        } else {
            PLB_ASSERT(false);  // Not implemented for 3 args.
        }
    }

    if (bdType == OffBoundary::constRhoInlet || bdType == OffBoundary::densityNeumann) {
        this->rhoBarVect[iDirection] = Descriptor<T>::rhoBar(wall_vel[0]);
    } else {
        this->rhoBarVect[iDirection] = rhoBar1;
    }
    Array<T, 3> wall_j(
        this->model.velIsJ() ? wall_vel
                             : (Descriptor<T>::fullRho(this->rhoBarVect[iDirection]) * wall_vel));
    if (depth < 2) {
        if (delta < (T)0.25) {
            this->jVect[iDirection] = wall_j;
        } else {
            // j = (T)1./delta * (wall_j+(delta-(T)1.)*j1);
            this->jVect[iDirection] = wall_j;  // Temporary fix for complex geometries.
        }
    } else {                    // depth >= 2           d=1   d=0.5   d=0
        if (delta < (T)0.75) {  // x---|=========o-------------o
            this->jVect[iDirection] =
                wall_j + (delta - (T)1.) * j1
                + ((T)1. - delta) / ((T)1. + delta) * ((T)2. * wall_j + (delta - (T)1.) * j2);
        }  //       d=1   d=0.5   d=0
        else
        {  // x===|---------o-------------o
            this->jVect[iDirection] = (T)1. / delta * (wall_j + (delta - (T)1.) * j1);
        }
    }
    if (bdType == OffBoundary::neumann) {
        this->jVect[iDirection] = j1;
    } else if (bdType == OffBoundary::densityNeumann) {
        this->jVect[iDirection] = dot(j1, wallNormal) * wallNormal;
    } else if (bdType == OffBoundary::freeSlip) {
        Array<T, 3> continuousNormal = this->model.computeContinuousNormal(wallNode, triangleId);
        this->jVect[iDirection] = j1 - dot(j1, continuousNormal) * continuousNormal;
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoDefineVelocityAlgorithm3D<T, Descriptor>::reduceVariables(T sumWeights)
{
    this->rhoBar = T();
    this->j.resetToZero();
    for (plint iDirection = 0; iDirection < this->numDirections; ++iDirection) {
        this->rhoBar += this->rhoBarVect[iDirection] * this->weights[iDirection];
        this->j += this->jVect[iDirection] * this->weights[iDirection];
    }
    this->rhoBar /= sumWeights;
    this->j /= sumWeights;
}

template <typename T, template <typename U> class Descriptor>
void GuoDefineVelocityAlgorithm3D<T, Descriptor>::complete()
{
    if (this->args.empty()) {
        Cell<T, Descriptor> &cell =
            this->lattice.get(this->guoNode.x, this->guoNode.y, this->guoNode.z);

        cell.getDynamics().computeEquilibria(
            cell.getRawPopulations(), this->rhoBar, this->j, normSqr(this->j));
    } else {
        if ((plint)this->args.size() == 1) {
            NTensorField3D<T> *macroField = dynamic_cast<NTensorField3D<T> *>(this->args[0]);
            PLB_ASSERT(macroField);
            // 1 Variable for rhoBar, 3 variables for j.
            PLB_ASSERT(macroField->getNdim() == 4);
            Dot3D offset = computeRelativeDisplacement(this->lattice, *macroField);
            T *macroscopic = macroField->get(
                this->guoNode.x + offset.x, this->guoNode.y + offset.y, this->guoNode.z + offset.z);
            // pcout << "rhoBar = " << this->rhoBar << ", ux = " << this->j[0] << ", uy = " <<
            // this->j[1] << ", uz = " << this->j[2] << std::endl;
            macroscopic[0] = this->rhoBar;
            this->j.to_cArray(macroscopic + 1);
        } else if ((plint)this->args.size() == 2) {
            // 1 field for rhoBar, 1 field for j.
            ScalarField3D<T> *rhoBarField = dynamic_cast<ScalarField3D<T> *>(this->args[0]);

            TensorField3D<T, 3> *jField = dynamic_cast<TensorField3D<T, 3> *>(this->args[1]);
            PLB_ASSERT(rhoBarField);
            PLB_ASSERT(jField);

            Dot3D offsetScalar = computeRelativeDisplacement(this->lattice, *rhoBarField);
            rhoBarField->get(
                this->guoNode.x + offsetScalar.x, this->guoNode.y + offsetScalar.y,
                this->guoNode.z + offsetScalar.z) = this->rhoBar;
            Dot3D offset = computeRelativeDisplacement(this->lattice, *jField);
            jField->get(
                this->guoNode.x + offset.x, this->guoNode.y + offset.y,
                this->guoNode.z + offset.z) = this->j;
        } else {
            PLB_ASSERT(false);  // Not implemented for 3 args.
        }
    }

    // Dynamics<T,Descriptor> const& dynamics = this->cell.getDynamics();
    // T jSqr = normSqr(this->j);
    // if (this->model.getPartialReplace()) {
    //     Cell<T,Descriptor> saveCell(this->cell);
    //     dynamics.regularize(this->cell, this->rhoBar, this->j, jSqr, PiNeq);
    //     for (plint iDirection=0; iDirection<this->numDirections; ++iDirection) {
    //         int iNeighbor = this->dryNodeFluidDirections[iDirection].first;
    //         plint iPop = nextNeighborPop<T,Descriptor>(iNeighbor);
    //         if (iPop>=0) {
    //             plint oppPop = indexTemplates::opposite<D>(iPop);
    //             this->cell[oppPop] = saveCell[oppPop];
    //         }
    //     }
    // }
    // else {
    //     dynamics.regularize(this->cell, this->rhoBar, this->j, jSqr, PiNeq);
    // }
}

template <typename T, template <typename U> class Descriptor>
class GuoFdCompletionAlgorithm3D : public GuoAlgorithm3D<T, Descriptor> {
public:
    typedef Descriptor<T> D;
    GuoFdCompletionAlgorithm3D(
        OffLatticeModel3D<T, Array<T, 3> > &model_, BlockLattice3D<T, Descriptor> &lattice_,
        Dot3D const &guoNode_, std::vector<std::pair<int, int> > const &dryNodeFluidDirections_,
        std::vector<plint> const &dryNodeIds_, Dot3D const &absoluteOffset_,
        Array<T, 3> &localForce_, std::vector<AtomicBlock3D *> const &args_, bool computeStat_,
        bool secondOrder_, std::pair<int, int> const &xDerivDirAndOrder,
        std::pair<int, int> const &yDerivDirAndOrder, std::pair<int, int> const &zDerivDirAndOrder);
    virtual void extrapolateVariables(
        Dot3D const &fluidDirection, int depth, Array<T, 3> const &wallNode, T delta,
        Array<T, 3> const &wall_vel, OffBoundary::Type bdType, Array<T, 3> const &wallNormal,
        plint triangleId, plint iDirection);
    virtual void reduceVariables(T sumWeights);
    virtual void complete();
    virtual bool computeNeighborData();

private:
    void computeCentralGradient(
        const BlockLattice3D<T, Descriptor> &lattice, const Dot3D &pos, const Dot3D &dx,
        Array<T, Descriptor<T>::d> &d_u) const;
    void computeCentralGradient(
        NTensorField3D<T> const *macroField, const Dot3D &pos, const Dot3D &dx,
        Array<T, Descriptor<T>::d> &d_u);
    void computeCentralGradient(
        ScalarField3D<T> const *rhoBarField, TensorField3D<T, 3> const *jField, const Dot3D &pos,
        const Dot3D &pos2, const Dot3D &dx, Array<T, Descriptor<T>::d> &d_u);
    void computeFwdGradient(
        const BlockLattice3D<T, Descriptor> &lattice, const Array<T, Descriptor<T>::d> &u0,
        const Dot3D &pos, const Dot3D &dx, int orient, Array<T, Descriptor<T>::d> &d_u) const;
    void computeFwdGradient(
        NTensorField3D<T> const *macroField, const Array<T, Descriptor<T>::d> &u0, const Dot3D &pos,
        const Dot3D &dx, int orient, Array<T, Descriptor<T>::d> &d_u);
    void computeFwdGradient(
        ScalarField3D<T> const *rhoBarField, TensorField3D<T, 3> const *jField,
        const Array<T, Descriptor<T>::d> &u0, const Dot3D &pos, const Dot3D &pos2, const Dot3D &dx,
        int orient, Array<T, Descriptor<T>::d> &d_u);

private:
    std::pair<int, int> xDerivDirAndOrder, yDerivDirAndOrder, zDerivDirAndOrder;
};

template <typename T, template <typename U> class Descriptor>
GuoFdCompletionAlgorithm3D<T, Descriptor>::GuoFdCompletionAlgorithm3D(
    OffLatticeModel3D<T, Array<T, 3> > &model_, BlockLattice3D<T, Descriptor> &lattice_,
    Dot3D const &guoNode_, std::vector<std::pair<int, int> > const &dryNodeFluidDirections_,
    std::vector<plint> const &dryNodeIds_, Dot3D const &absoluteOffset_, Array<T, 3> &localForce_,
    std::vector<AtomicBlock3D *> const &args_, bool computeStat_, bool secondOrder_,
    std::pair<int, int> const &xDerivDirAndOrder_, std::pair<int, int> const &yDerivDirAndOrder_,
    std::pair<int, int> const &zDerivDirAndOrder_) :
    GuoAlgorithm3D<T, Descriptor>(
        model_, lattice_, guoNode_, dryNodeFluidDirections_, dryNodeIds_, absoluteOffset_,
        localForce_, args_, computeStat_, secondOrder_),
    xDerivDirAndOrder(xDerivDirAndOrder_),
    yDerivDirAndOrder(yDerivDirAndOrder_),
    zDerivDirAndOrder(zDerivDirAndOrder_)
{ }

template <typename T, template <typename U> class Descriptor>
void GuoFdCompletionAlgorithm3D<T, Descriptor>::extrapolateVariables(
    Dot3D const &fluidDirection, int depth, Array<T, 3> const &wallNode, T delta,
    Array<T, 3> const &wall_vel, OffBoundary::Type bdType, Array<T, 3> const &wallNormal,
    plint triangleId, plint iDirection)
{
    PLB_ASSERT(false && "Nothing to extrapolate in fd guo completion.");
}

template <typename T, template <typename U> class Descriptor>
bool GuoFdCompletionAlgorithm3D<T, Descriptor>::computeNeighborData()
{
    reduceVariables((T)0);
    return true;
}

template <typename T, template <typename U> class Descriptor>
void GuoFdCompletionAlgorithm3D<T, Descriptor>::reduceVariables([[maybe_unused]] T sumWeights)
{
    if (this->args.empty()) {
        this->cell.getDynamics().computeRhoBarJ(this->cell, this->rhoBar, this->j);
    } else {
        if ((plint)this->args.size() == 1) {
            NTensorField3D<T> const *macroField =
                dynamic_cast<NTensorField3D<T> const *>(this->args[0]);
            PLB_ASSERT(macroField);
            // 1 Variable for rhoBar, 3 variables for j.
            PLB_ASSERT(macroField->getNdim() == 4);
            Dot3D offset = computeRelativeDisplacement(this->lattice, *macroField);
            T const *macro_0 = macroField->get(
                this->guoNode.x + offset.x, this->guoNode.y + offset.y, this->guoNode.z + offset.z);
            this->rhoBar = macro_0[0];
            this->j.from_cArray(macro_0 + 1);
        } else if ((plint)this->args.size() == 2) {
            // 1 field for rhoBar, 1 field for j.
            ScalarField3D<T> const *rhoBarField =
                dynamic_cast<ScalarField3D<T> const *>(this->args[0]);
            TensorField3D<T, 3> const *jField =
                dynamic_cast<TensorField3D<T, 3> const *>(this->args[1]);
            PLB_ASSERT(rhoBarField);
            PLB_ASSERT(jField);

            Dot3D offsetScalar = computeRelativeDisplacement(this->lattice, *rhoBarField);
            Dot3D offset = computeRelativeDisplacement(this->lattice, *jField);

            this->j = jField->get(
                this->guoNode.x + offset.x, this->guoNode.y + offset.y, this->guoNode.z + offset.z);
            this->rhoBar = rhoBarField->get(
                this->guoNode.x + offsetScalar.x, this->guoNode.y + offsetScalar.y,
                this->guoNode.z + offsetScalar.z);
        } else {
            PLB_ASSERT(false);  // Not implemented for 3 args.
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoFdCompletionAlgorithm3D<T, Descriptor>::computeCentralGradient(
    const BlockLattice3D<T, Descriptor> &lattice, const Dot3D &pos, const Dot3D &dx,
    Array<T, Descriptor<T>::d> &d_u) const
{
    Dot3D p1 = pos + dx;
    Dot3D m1 = pos - dx;

    Array<T, Descriptor<T>::d> u_p1;
    lattice.get(p1.x, p1.y, p1.z).computeVelocity(u_p1);
    Array<T, Descriptor<T>::d> u_m1;
    lattice.get(m1.x, m1.y, m1.z).computeVelocity(u_m1);

    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        d_u[iD] = fd::ctl_diff(u_p1[iD], u_m1[iD]);
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoFdCompletionAlgorithm3D<T, Descriptor>::computeCentralGradient(
    NTensorField3D<T> const *macroField, const Dot3D &pos, const Dot3D &dx,
    Array<T, Descriptor<T>::d> &d_u)
{
    Dot3D p1 = pos + dx;
    Dot3D m1 = pos - dx;

    T const *macro_p1 = macroField->get(p1.x, p1.y, p1.z);
    T const *macro_m1 = macroField->get(m1.x, m1.y, m1.z);

    Array<T, Descriptor<T>::d> u_p1;
    u_p1.from_cArray(macro_p1 + 1);
    Array<T, Descriptor<T>::d> u_m1;
    u_m1.from_cArray(macro_m1 + 1);

    if (!this->model.velIsJ()) {
        u_p1 *= Descriptor<T>::invRho(macro_p1[0]);
        u_m1 *= Descriptor<T>::invRho(macro_m1[0]);
    }

    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        d_u[iD] = fd::ctl_diff(u_p1[iD], u_m1[iD]);
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoFdCompletionAlgorithm3D<T, Descriptor>::computeCentralGradient(
    ScalarField3D<T> const *rhoBarField, TensorField3D<T, 3> const *jField, const Dot3D &pos,
    const Dot3D &pos2, const Dot3D &dx, Array<T, Descriptor<T>::d> &d_u)
{
    Dot3D p1 = pos + dx;  // position for the tensor field
    Dot3D m1 = pos - dx;

    Array<T, 3> u_p1 = jField->get(p1.x, p1.y, p1.z);
    Array<T, 3> u_m1 = jField->get(m1.x, m1.y, m1.z);

    p1 = pos2 + dx;  // position for the scalar field
    m1 = pos2 - dx;
    if (!this->model.velIsJ()) {
        u_p1 *= Descriptor<T>::invRho(rhoBarField->get(p1.x, p1.y, p1.z));
        u_m1 *= Descriptor<T>::invRho(rhoBarField->get(m1.x, m1.y, m1.z));
    }

    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        d_u[iD] = fd::ctl_diff(u_p1[iD], u_m1[iD]);
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoFdCompletionAlgorithm3D<T, Descriptor>::computeFwdGradient(
    const BlockLattice3D<T, Descriptor> &lattice, const Array<T, Descriptor<T>::d> &u0,
    const Dot3D &pos, const Dot3D &dx, int orient, Array<T, Descriptor<T>::d> &d_u) const
{
    Dot3D p1 = pos + dx;
    Array<T, Descriptor<T>::d> u_p1;
    lattice.get(p1.x, p1.y, p1.z).computeVelocity(u_p1);

    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        d_u[iD] = (T)orient * fd::o1_fwd_diff(u0[iD], u_p1[iD]);
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoFdCompletionAlgorithm3D<T, Descriptor>::computeFwdGradient(
    NTensorField3D<T> const *macroField, const Array<T, Descriptor<T>::d> &u0, const Dot3D &pos,
    const Dot3D &dx, int orient, Array<T, Descriptor<T>::d> &d_u)
{
    Dot3D p1 = pos + dx;
    T const *macro_p1 = macroField->get(p1.x, p1.y, p1.z);
    Array<T, Descriptor<T>::d> u_p1;
    u_p1.from_cArray(macro_p1 + 1);

    if (!this->model.velIsJ()) {
        u_p1 *= Descriptor<T>::invRho(macro_p1[0]);
    }

    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        d_u[iD] = (T)orient * fd::o1_fwd_diff(u0[iD], u_p1[iD]);
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoFdCompletionAlgorithm3D<T, Descriptor>::computeFwdGradient(
    ScalarField3D<T> const *rhoBarField, TensorField3D<T, 3> const *jField,
    const Array<T, Descriptor<T>::d> &u0, const Dot3D &pos, const Dot3D &pos2, const Dot3D &dx,
    int orient, Array<T, Descriptor<T>::d> &d_u)
{
    Dot3D p1 = pos + dx;  // position for the tensor field
    Array<T, 3> u_p1 = jField->get(p1.x, p1.y, p1.z);

    p1 = pos2 + dx;  // position for the scalar field
    if (!this->model.velIsJ()) {
        u_p1 *= Descriptor<T>::invRho(rhoBarField->get(p1.x, p1.y, p1.z));
    }

    for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
        d_u[iD] = (T)orient * fd::o1_fwd_diff(u0[iD], u_p1[iD]);
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoFdCompletionAlgorithm3D<T, Descriptor>::complete()
{
    Array<T, 3> u0 = this->j;
    if (!this->model.velIsJ()) {
        u0 *= Descriptor<T>::invRho(this->rhoBar);
    }

    Array<T, Descriptor<T>::d> dx_u, dy_u, dz_u;
    if (this->args.empty()) {
        if (xDerivDirAndOrder.second == 2) {
            computeCentralGradient(this->lattice, this->guoNode, Dot3D(1, 0, 0), dx_u);
        } else if (xDerivDirAndOrder.second == 1) {
            PLB_ASSERT(
                (xDerivDirAndOrder.first == 1 || xDerivDirAndOrder.first == -1)
                && "Invalid xDerivDirAndOrder.first value (must be = 1 or = -1");
            computeFwdGradient(
                this->lattice, u0, this->guoNode, Dot3D(xDerivDirAndOrder.first, 0, 0),
                xDerivDirAndOrder.first, dx_u);
        } else if (xDerivDirAndOrder.second == 0) {
            dx_u.resetToZero();
        } else {
            PLB_ASSERT(false && "xDerivDirAndOrder has invalid order");
        }
        if (yDerivDirAndOrder.second == 2) {
            computeCentralGradient(this->lattice, this->guoNode, Dot3D(0, 1, 0), dy_u);
        } else if (yDerivDirAndOrder.second == 1) {
            PLB_ASSERT(
                (yDerivDirAndOrder.first == 1 || yDerivDirAndOrder.first == -1)
                && "Invalid yDerivDirAndOrder.first value (must be = 1 or = -1");
            computeFwdGradient(
                this->lattice, u0, this->guoNode, Dot3D(0, yDerivDirAndOrder.first, 0),
                yDerivDirAndOrder.first, dy_u);
        } else if (yDerivDirAndOrder.second == 0) {
            dy_u.resetToZero();
        } else {
            PLB_ASSERT(false && "yDerivDirAndOrder has invalid order");
        }
        if (zDerivDirAndOrder.second == 2) {
            computeCentralGradient(this->lattice, this->guoNode, Dot3D(0, 0, 1), dz_u);
        } else if (zDerivDirAndOrder.second == 1) {
            PLB_ASSERT(
                (zDerivDirAndOrder.first == 1 || zDerivDirAndOrder.first == -1)
                && "Invalid zDerivDirAndOrder.first value (must be = 1 or = -1");
            computeFwdGradient(
                this->lattice, u0, this->guoNode, Dot3D(0, 0, zDerivDirAndOrder.first),
                zDerivDirAndOrder.first, dz_u);
        } else if (zDerivDirAndOrder.second == 0) {
            dz_u.resetToZero();
        } else {
            PLB_ASSERT(false && "zDerivDirAndOrder has invalid order");
        }

    } else {
        if ((plint)this->args.size() == 1) {
            NTensorField3D<T> const *macroField =
                dynamic_cast<NTensorField3D<T> const *>(this->args[0]);
            PLB_ASSERT(macroField);
            // 1 Variable for rhoBar, 3 variables for j.
            PLB_ASSERT(macroField->getNdim() == 4);
            Dot3D offset = computeRelativeDisplacement(this->lattice, *macroField);
            Array<T, 3> u0 = this->j;
            if (!this->model.velIsJ()) {
                u0 *= Descriptor<T>::invRho(this->rhoBar);
            }
            if (xDerivDirAndOrder.second == 2) {
                computeCentralGradient(macroField, this->guoNode + offset, Dot3D(1, 0, 0), dx_u);
            } else if (xDerivDirAndOrder.second == 1) {
                PLB_ASSERT(
                    (xDerivDirAndOrder.first == 1 || xDerivDirAndOrder.first == -1)
                    && "Invalid xDerivDirAndOrder.first value (must be = 1 or = -1");
                computeFwdGradient(
                    macroField, u0, this->guoNode + offset, Dot3D(xDerivDirAndOrder.first, 0, 0),
                    xDerivDirAndOrder.first, dx_u);
            } else if (xDerivDirAndOrder.second == 0) {
                dx_u.resetToZero();
            } else {
                PLB_ASSERT(false && "xDerivDirAndOrder has invalid order");
            }
            if (yDerivDirAndOrder.second == 2) {
                computeCentralGradient(macroField, this->guoNode + offset, Dot3D(0, 1, 0), dy_u);
            } else if (yDerivDirAndOrder.second == 1) {
                PLB_ASSERT(
                    (yDerivDirAndOrder.first == 1 || yDerivDirAndOrder.first == -1)
                    && "Invalid yDerivDirAndOrder.first value (must be = 1 or = -1");
                computeFwdGradient(
                    macroField, u0, this->guoNode + offset, Dot3D(0, yDerivDirAndOrder.first, 0),
                    yDerivDirAndOrder.first, dy_u);
            } else if (yDerivDirAndOrder.second == 0) {
                dy_u.resetToZero();
            } else {
                PLB_ASSERT(false && "yDerivDirAndOrder has invalid order");
            }
            if (zDerivDirAndOrder.second == 2) {
                computeCentralGradient(macroField, this->guoNode + offset, Dot3D(0, 0, 1), dz_u);
            } else if (zDerivDirAndOrder.second == 1) {
                PLB_ASSERT(
                    (zDerivDirAndOrder.first == 1 || zDerivDirAndOrder.first == -1)
                    && "Invalid zDerivDirAndOrder.first value (must be = 1 or = -1");
                computeFwdGradient(
                    macroField, u0, this->guoNode + offset, Dot3D(0, 0, zDerivDirAndOrder.first),
                    zDerivDirAndOrder.first, dz_u);
            } else if (zDerivDirAndOrder.second == 0) {
                dz_u.resetToZero();
            } else {
                PLB_ASSERT(false && "zDerivDirAndOrder has invalid order");
            }

        } else if ((plint)this->args.size() == 2) {
            // 1 field for rhoBar, 1 field for j.
            ScalarField3D<T> const *rhoBarField =
                dynamic_cast<ScalarField3D<T> const *>(this->args[0]);

            TensorField3D<T, 3> const *jField =
                dynamic_cast<TensorField3D<T, 3> const *>(this->args[1]);
            PLB_ASSERT(rhoBarField);
            PLB_ASSERT(jField);

            Dot3D offsetScalar = computeRelativeDisplacement(this->lattice, *rhoBarField);
            Dot3D offset = computeRelativeDisplacement(this->lattice, *jField);

            if (xDerivDirAndOrder.second == 2) {
                computeCentralGradient(
                    rhoBarField, jField, this->guoNode + offset, this->guoNode + offsetScalar,
                    Dot3D(1, 0, 0), dx_u);
            } else if (xDerivDirAndOrder.second == 1) {
                PLB_ASSERT(
                    (xDerivDirAndOrder.first == 1 || xDerivDirAndOrder.first == -1)
                    && "Invalid xDerivDirAndOrder.first value (must be = 1 or = -1");
                computeFwdGradient(
                    rhoBarField, jField, u0, this->guoNode + offset, this->guoNode + offsetScalar,
                    Dot3D(xDerivDirAndOrder.first, 0, 0), xDerivDirAndOrder.first, dx_u);
            } else if (xDerivDirAndOrder.second == 0) {
                dx_u.resetToZero();
            } else {
                PLB_ASSERT(false && "xDerivDirAndOrder has invalid order");
            }
            if (yDerivDirAndOrder.second == 2) {
                computeCentralGradient(
                    rhoBarField, jField, this->guoNode + offset, this->guoNode + offsetScalar,
                    Dot3D(0, 1, 0), dy_u);
            } else if (yDerivDirAndOrder.second == 1) {
                PLB_ASSERT(
                    (yDerivDirAndOrder.first == 1 || yDerivDirAndOrder.first == -1)
                    && "Invalid yDerivDirAndOrder.first value (must be = 1 or = -1");
                computeFwdGradient(
                    rhoBarField, jField, u0, this->guoNode + offset, this->guoNode + offsetScalar,
                    Dot3D(0, yDerivDirAndOrder.first, 0), yDerivDirAndOrder.first, dy_u);
            } else if (yDerivDirAndOrder.second == 0) {
                dy_u.resetToZero();
            } else {
                PLB_ASSERT(false && "yDerivDirAndOrder has invalid order");
            }
            if (zDerivDirAndOrder.second == 2) {
                computeCentralGradient(
                    rhoBarField, jField, this->guoNode + offset, this->guoNode + offsetScalar,
                    Dot3D(0, 0, 1), dz_u);
            } else if (zDerivDirAndOrder.second == 1) {
                PLB_ASSERT(
                    (zDerivDirAndOrder.first == 1 || zDerivDirAndOrder.first == -1)
                    && "Invalid zDerivDirAndOrder.first value (must be = 1 or = -1");
                computeFwdGradient(
                    rhoBarField, jField, u0, this->guoNode + offset, this->guoNode + offsetScalar,
                    Dot3D(0, 0, zDerivDirAndOrder.first), zDerivDirAndOrder.first, dz_u);
            } else if (zDerivDirAndOrder.second == 0) {
                dz_u.resetToZero();
            } else {
                PLB_ASSERT(false && "zDerivDirAndOrder has invalid order");
            }

        } else {
            PLB_ASSERT(false);  // Not implemented for 3 args.
        }
    }

    T dx_ux = dx_u[0];
    T dy_ux = dy_u[0];
    T dz_ux = dz_u[0];
    T dx_uy = dx_u[1];
    T dy_uy = dy_u[1];
    T dz_uy = dz_u[1];
    T dx_uz = dx_u[2];
    T dy_uz = dy_u[2];
    T dz_uz = dz_u[2];
    T omega = this->cell.getDynamics().getOmega();
    T rho = Descriptor<T>::fullRho(this->rhoBar);
    T sToPi = -rho * Descriptor<T>::cs2 / omega;

    typedef SymmetricTensorImpl<T, 3> S;
    Array<T, SymmetricTensor<T, Descriptor>::n> pi;
    pi[S::xx] = (T)2 * dx_ux * sToPi;
    pi[S::yy] = (T)2 * dy_uy * sToPi;
    pi[S::zz] = (T)2 * dz_uz * sToPi;
    pi[S::xy] = (dx_uy + dy_ux) * sToPi;
    pi[S::xz] = (dx_uz + dz_ux) * sToPi;
    pi[S::yz] = (dy_uz + dz_uy) * sToPi;

    T jSqr = VectorTemplate<T, Descriptor>::normSqr(this->j);

    this->cell.getDynamics().regularize(this->cell, this->rhoBar, this->j, jSqr, pi);
}

// ============================================================================ //
// ==================== GuoOffLatticeFdModel3D ==================== //
// ============================================================================ //

template <typename T, template <typename U> class Descriptor>
GuoOffLatticeFdModel3D<T, Descriptor>::LiquidNeighbor::LiquidNeighbor(
    plint iNeighbor_, plint depth_, plint iTriangle_, Array<T, 3> wallNormal) :
    iNeighbor(iNeighbor_), depth(depth_), iTriangle(iTriangle_)
{
    int const *c = NextNeighbor<T>::c[iNeighbor];
    Array<T, 3> neighborVect(c[0], c[1], c[2]);
    cosAngle = std::fabs(dot(neighborVect, wallNormal)) * NextNeighbor<T>::invD[iNeighbor];
}

template <typename T, template <typename U> class Descriptor>
bool GuoOffLatticeFdModel3D<T, Descriptor>::LiquidNeighbor::operator<(
    LiquidNeighbor const &rhs) const
{
    return cosAngle < rhs.cosAngle;
}

template <typename T, template <typename U> class Descriptor>
GuoOffLatticeFdModel3D<T, Descriptor>::GuoOffLatticeFdModel3D(
    BoundaryShape3D<T, Array<T, 3> > *shape_, int flowType_, bool useAllDirections_) :
    OffLatticeModel3D<T, Array<T, 3> >(shape_, flowType_), useAllDirections(useAllDirections_)
{ }

template <typename T, template <typename U> class Descriptor>
GuoOffLatticeFdModel3D<T, Descriptor> *GuoOffLatticeFdModel3D<T, Descriptor>::clone() const
{
    return new GuoOffLatticeFdModel3D(*this);
}

template <typename T, template <typename U> class Descriptor>
plint GuoOffLatticeFdModel3D<T, Descriptor>::getNumNeighbors() const
{
    return 2;
}

template <typename T, template <typename U> class Descriptor>
bool GuoOffLatticeFdModel3D<T, Descriptor>::isExtrapolated() const
{
    return true;
}

template <typename T, template <typename U> class Descriptor>
bool GuoOffLatticeFdModel3D<T, Descriptor>::isUsable(const Dot3D &pos) const
{
    if (this->isFluid(pos)) {
        return true;
    } else if (this->isSolid(pos)) {
        for (int iNeighbor = 0; iNeighbor < NextNeighbor<T>::numNeighbors; ++iNeighbor) {
            int const *c = NextNeighbor<T>::c[iNeighbor];
            Dot3D neighborPos(pos.x + c[0], pos.y + c[1], pos.z + c[2]);
            if (this->isFluid(neighborPos)) {
                return true;
            }
        }
    }

    return false;
}

template <typename T, template <typename U> class Descriptor>
std::pair<int, int> GuoOffLatticeFdModel3D<T, Descriptor>::computeOrderAndDirection(
    const Dot3D &pos, const Dot3D &dx) const
{
    if (isUsable(pos + dx) && isUsable(pos - dx)) {
        return std::make_pair(0, 2);
    } else if (isUsable(pos + dx)) {
        return std::make_pair(1, 1);
    } else if (isUsable(pos - dx)) {
        return std::make_pair(-1, 1);
    }
    return std::make_pair(0, 0);
}

template <typename T, template <typename U> class Descriptor>
void GuoOffLatticeFdModel3D<T, Descriptor>::prepareCell(
    Dot3D const &cellLocation, AtomicContainerBlock3D &container)
{
    Dot3D offset = container.getLocation();
    GuoOffLatticeInfo3D *info = dynamic_cast<GuoOffLatticeInfo3D *>(container.getData());
    PLB_ASSERT(info);
    std::vector<LiquidNeighbor> liquidNeighbors;
    if (this->isSolid(cellLocation + offset)) {
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
                plint iTriangle = -1;
                global::timer("intersect").start();
                Array<T, 3> locatedPoint;
                T distance;
                Array<T, 3> wallNormal;
                Array<T, 3> surfaceData;
                OffBoundary::Type bdType;
                bool ok = this->pointOnSurface(
                    cellLocation + offset, Dot3D(c[0], c[1], c[2]), locatedPoint, distance,
                    wallNormal, surfaceData, bdType, iTriangle);
                // In the following, the importance of directions is sorted wrt. how well they
                //   are aligned with the wall normal. It is better to take the continuous normal,
                //   because it is not sensitive to the choice of the triangle when we shoot at
                //   an edge.
                // wallNormal = this->computeContinuousNormal(locatedPoint, iTriangle);
                global::timer("intersect").stop();
                if (!ok) {
                    global::plbErrors().registerError(
                        "Guo off-lattice model could not find an intersection with a triangle.");
                }
                // ... then add this node to the list.
                liquidNeighbors.push_back(LiquidNeighbor(iNeighbor, depth, iTriangle, wallNormal));
            }
        }

        if (!liquidNeighbors.empty()) {
            info->getDryNodes().push_back(cellLocation);
            std::sort(liquidNeighbors.begin(), liquidNeighbors.end());
            std::vector<std::pair<int, int> > neighborDepthPairs;
            std::vector<plint> ids;

            Dot3D pos = cellLocation + offset;
            Dot3D dx(1, 0, 0), dy(0, 1, 0), dz(0, 0, 1);
            info->getXderivDirAndOrder().push_back(computeOrderAndDirection(pos, dx));
            info->getYderivDirAndOrder().push_back(computeOrderAndDirection(pos, dy));
            info->getZderivDirAndOrder().push_back(computeOrderAndDirection(pos, dz));

            if (useAllDirections) {
                for (pluint i = 0; i < liquidNeighbors.size(); ++i) {
                    neighborDepthPairs.push_back(
                        std::make_pair(liquidNeighbors[i].iNeighbor, liquidNeighbors[i].depth));
                    ids.push_back(liquidNeighbors[i].iTriangle);
                }
            } else {
                plint i = liquidNeighbors.size() - 1;
                neighborDepthPairs.push_back(
                    std::make_pair(liquidNeighbors[i].iNeighbor, liquidNeighbors[i].depth));
                ids.push_back(liquidNeighbors[i].iTriangle);
            }
            info->getDryNodeFluidDirections().push_back(neighborDepthPairs);
            info->getDryNodeIds().push_back(ids);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ContainerBlockData *GuoOffLatticeFdModel3D<T, Descriptor>::generateOffLatticeInfo() const
{
    return new GuoOffLatticeInfo3D;
}

template <typename T, template <typename U> class Descriptor>
Array<T, 3> GuoOffLatticeFdModel3D<T, Descriptor>::getLocalForce(
    AtomicContainerBlock3D &container) const
{
    GuoOffLatticeInfo3D *info = dynamic_cast<GuoOffLatticeInfo3D *>(container.getData());
    PLB_ASSERT(info);
    return info->getLocalForce();
}

template <typename T, template <typename U> class Descriptor>
void GuoOffLatticeFdModel3D<T, Descriptor>::boundaryCompletion(
    AtomicBlock3D &nonTypeLattice, AtomicContainerBlock3D &container,
    std::vector<AtomicBlock3D *> const &args)
{
    BlockLattice3D<T, Descriptor> &lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> &>(nonTypeLattice);
    GuoOffLatticeInfo3D *info = dynamic_cast<GuoOffLatticeInfo3D *>(container.getData());
    PLB_ASSERT(info);
    std::vector<Dot3D> const &dryNodes = info->getDryNodes();
    std::vector<std::vector<std::pair<int, int> > > const &dryNodeFluidDirections =
        info->getDryNodeFluidDirections();
    std::vector<std::vector<plint> > const &dryNodeIds = info->getDryNodeIds();
    std::vector<std::pair<int, int> > const &xDerivDirAndOrder = info->getXderivDirAndOrder();
    std::vector<std::pair<int, int> > const &yDerivDirAndOrder = info->getYderivDirAndOrder();
    std::vector<std::pair<int, int> > const &zDerivDirAndOrder = info->getZderivDirAndOrder();
    if (dryNodes.size() != dryNodeFluidDirections.size()) {
        global::plbErrors().registerError(
            "Error in the Guo off-lattice Fd model boundary completion.");
    }

    Dot3D absoluteOffset = container.getLocation();

    Array<T, 3> &localForce = info->getLocalForce();
    localForce.resetToZero();
    for (pluint iDry = 0; iDry < dryNodes.size(); ++iDry) {
        cellCompletion(
            lattice, dryNodes[iDry], dryNodeFluidDirections[iDry], dryNodeIds[iDry], absoluteOffset,
            xDerivDirAndOrder[iDry], yDerivDirAndOrder[iDry], zDerivDirAndOrder[iDry], localForce,
            args);
    }
}

template <typename T, template <typename U> class Descriptor>
void GuoOffLatticeFdModel3D<T, Descriptor>::cellCompletion(
    BlockLattice3D<T, Descriptor> &lattice, Dot3D const &guoNode,
    std::vector<std::pair<int, int> > const &dryNodeFluidDirections,
    std::vector<plint> const &dryNodeIds, Dot3D const &absoluteOffset,
    const std::pair<int, int> &xDerivDirAndOrder, const std::pair<int, int> &yDerivDirAndOrder,
    const std::pair<int, int> &zDerivDirAndOrder, Array<T, 3> &localForce,
    std::vector<AtomicBlock3D *> const &args)
{
    GuoAlgorithm3D<T, Descriptor> *algorithm = 0;
    if (this->getDefineVelocity()) {
        algorithm = new GuoDefineVelocityAlgorithm3D<T, Descriptor>(
            *this, lattice, guoNode, dryNodeFluidDirections, dryNodeIds, absoluteOffset, localForce,
            args, this->computesStat(), this->usesSecondOrder());
    } else {
        algorithm = new GuoFdCompletionAlgorithm3D<T, Descriptor>(
            *this, lattice, guoNode, dryNodeFluidDirections, dryNodeIds, absoluteOffset, localForce,
            args, this->computesStat(), this->usesSecondOrder(), xDerivDirAndOrder,
            yDerivDirAndOrder, zDerivDirAndOrder);
    }
    bool ok = algorithm->computeNeighborData();
    if (!ok) {
        global::plbErrors().registerError(
            "Error treating the geometry in the Guo off-lattice model.");
    }
    algorithm->finalize();
    delete algorithm;
}

}  // namespace plb

#endif  // GUO_OFF_LATTICE_MODEL_3D_HH
