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

#ifndef GENERALIZED_OFF_LATTICE_MODEL_3D_HH
#define GENERALIZED_OFF_LATTICE_MODEL_3D_HH

#include <algorithm>
#include <cmath>

#include "core/runTimeDiagnostics.h"
#include "guoOffLatticeModel3D.hh"
#include "latticeBoltzmann/externalFieldAccess.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "offLattice/generalizedOffLatticeModel3D.h"
#include "offLattice/nextNeighbors3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class ExtrapolationGeneralizedAlgorithm3D : public GuoAlgorithm3D<T, Descriptor> {
public:
    typedef Descriptor<T> D;
    ExtrapolationGeneralizedAlgorithm3D(
        OffLatticeModel3D<T, Array<T, 3> > &model_, BlockLattice3D<T, Descriptor> &lattice_,
        Dot3D const &guoNode_, std::vector<std::pair<int, int> > const &dryNodeFluidDirections_,
        std::vector<int> dryNodeFluidNoSolidDirections_, std::vector<plint> const &dryNodeIds_,
        Dot3D const &absoluteOffset_, Array<T, 3> &localForce_,
        std::vector<AtomicBlock3D *> const &args_, bool computeStat_, bool secondOrder_);
    virtual void extrapolateVariables(
        Dot3D const &fluidDirection, int depth, Array<T, 3> const &wallNode, T delta,
        Array<T, 3> const &wall_vel, OffBoundary::Type bdType, Array<T, 3> const &wallNormal,
        plint triangleId, plint iDirection);
    virtual void reduceVariables(T sumWeights);
    virtual void complete();

private:
    std::vector<int> dryNodeFluidNoSolidDirections;
};

template <typename T, template <typename U> class Descriptor>
ExtrapolationGeneralizedAlgorithm3D<T, Descriptor>::ExtrapolationGeneralizedAlgorithm3D(
    OffLatticeModel3D<T, Array<T, 3> > &model_, BlockLattice3D<T, Descriptor> &lattice_,
    Dot3D const &guoNode_, std::vector<std::pair<int, int> > const &dryNodeFluidDirections_,
    std::vector<int> dryNodeFluidNoSolidDirections_, std::vector<plint> const &dryNodeIds_,
    Dot3D const &absoluteOffset_, Array<T, 3> &localForce_,
    std::vector<AtomicBlock3D *> const &args_, bool computeStat_, bool secondOrder_) :
    GuoAlgorithm3D<T, Descriptor>(
        model_, lattice_, guoNode_, dryNodeFluidDirections_, dryNodeIds_, absoluteOffset_,
        localForce_, args_, computeStat_, secondOrder_),
    dryNodeFluidNoSolidDirections(dryNodeFluidNoSolidDirections_)
{ }

template <typename T, template <typename U> class Descriptor>
void ExtrapolationGeneralizedAlgorithm3D<T, Descriptor>::extrapolateVariables(
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
void ExtrapolationGeneralizedAlgorithm3D<T, Descriptor>::reduceVariables(T sumWeights)
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
void ExtrapolationGeneralizedAlgorithm3D<T, Descriptor>::complete()
{
    std::vector<plint> knownIndices;
    knownIndices.push_back(0);
    for (pluint iDirection = 0; iDirection < dryNodeFluidNoSolidDirections.size(); ++iDirection) {
        int iNeighbor = dryNodeFluidNoSolidDirections[iDirection];
        plint iPop = nextNeighborPop<T, Descriptor>(iNeighbor);
        if (iPop >= 0) {
            plint index = indexTemplates::opposite<Descriptor<T> >(iPop);
            knownIndices.push_back(index);
        }
    }
    PLB_ASSERT(knownIndices.size() >= 6);

    std::vector<plint> missingIndices =
        indexTemplates::remainingIndexes<Descriptor<T> >(knownIndices);
    Array<T, Descriptor<T>::d> u = this->j;
    if (!this->model.velIsJ())
        u *= Descriptor<T>::invRho(this->rhoBar);
    Dynamics<T, Descriptor> const &dynamics = this->cell.getDynamics();
    DirichletVelocityBoundarySolver<T, Descriptor> bc(missingIndices, knownIndices, u);
    bc.apply(this->cell, dynamics, !this->model.getPartialReplace());
}

template <typename T, template <typename U> class Descriptor>
ExtrapolatedGeneralizedOffLatticeModel3D<T, Descriptor>::ExtrapolatedGeneralizedOffLatticeModel3D(
    BoundaryShape3D<T, Array<T, 3> > *shape_, int flowType_) :
    OffLatticeModel3D<T, Array<T, 3> >(shape_, flowType_)
{ }

template <typename T, template <typename U> class Descriptor>
ExtrapolatedGeneralizedOffLatticeModel3D<T, Descriptor>
    *ExtrapolatedGeneralizedOffLatticeModel3D<T, Descriptor>::clone() const
{
    return new ExtrapolatedGeneralizedOffLatticeModel3D(*this);
}

template <typename T, template <typename U> class Descriptor>
plint ExtrapolatedGeneralizedOffLatticeModel3D<T, Descriptor>::getNumNeighbors() const
{
    return 2;
}

template <typename T, template <typename U> class Descriptor>
bool ExtrapolatedGeneralizedOffLatticeModel3D<T, Descriptor>::isExtrapolated() const
{
    return true;
}

template <typename T, template <typename U> class Descriptor>
void ExtrapolatedGeneralizedOffLatticeModel3D<T, Descriptor>::prepareCell(
    Dot3D const &cellLocation, AtomicContainerBlock3D &container)
{
    Dot3D offset = container.getLocation();
    ExtrapolatedGeneralizedOffLatticeInfo3D *info =
        dynamic_cast<ExtrapolatedGeneralizedOffLatticeInfo3D *>(container.getData());
    PLB_ASSERT(info);
    if (this->isSolid(cellLocation + offset)) {
        std::vector<std::pair<int, int> > liquidNeighbors;
        std::vector<int> liquidNeighborsNoSolid;
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
                Array<T, 3> locatedPoint;
                T distance;
                Array<T, 3> wallNormal;
                Array<T, 3> surfaceData;
                plint iTriangle = -1;
                OffBoundary::Type bdType;
#ifdef PLB_DEBUG
                bool ok =
#endif
                    this->pointOnSurface(
                        cellLocation + offset, Dot3D(c[0], c[1], c[2]), locatedPoint, distance,
                        wallNormal, surfaceData, bdType, iTriangle);
                PLB_ASSERT(ok);
                ids.push_back(iTriangle);
                liquidNeighborsNoSolid.push_back(iNeighbor);
            } else {
                bool fluidNeighbor = false;
                for (int jNeighbor = 1; jNeighbor < Descriptor<T>::q; ++jNeighbor) {
                    Dot3D next_neighbor(
                        neighbor.x + Descriptor<T>::c[jNeighbor][0],
                        neighbor.y + Descriptor<T>::c[jNeighbor][1],
                        neighbor.z + Descriptor<T>::c[jNeighbor][2]);
                    if ((this->isFluid(next_neighbor + offset))) {
                        fluidNeighbor = true;
                        break;
                    }
                }
                if (fluidNeighbor) {
                    liquidNeighborsNoSolid.push_back(iNeighbor);
                }
            }
        }
        if (!liquidNeighbors.empty()) {
            // selecting only deepest depth available (for a higher order interpolation order)
            //             std::vector<std::vector<std::pair<int,int> > >
            //             selectionLiquid(getNumNeighbors()+1); std::vector<std::vector<plint > >
            //             selectionIds(getNumNeighbors()+1); for (pluint iA = 0; iA <
            //             liquidNeighbors.size(); ++iA) {
            //                 selectionLiquid[liquidNeighbors[iA].second].push_back(liquidNeighbors[iA]);
            //                 selectionIds[liquidNeighbors[iA].second].push_back(ids[iA]);
            //             }

            info->getDryNodes().push_back(cellLocation);
            info->getDryNodeFluidWithFluidDirections().push_back(liquidNeighborsNoSolid);
            info->getDryNodeFluidDirections().push_back(liquidNeighbors);
            info->getDryNodeIds().push_back(ids);

            // add only biggest depth to the list of directions and triangles
            //             for (plint iA = getNumNeighbors(); iA >= 0; --iA) {
            //                 if (!selectionLiquid[iA].empty()) {
            //                     info->getDryNodeFluidDirections().push_back(selectionLiquid[iA]);
            //                     info->getDryNodeIds().push_back(selectionIds[iA]);
            //                     break;
            //                 }
            //             }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ContainerBlockData *
    ExtrapolatedGeneralizedOffLatticeModel3D<T, Descriptor>::generateOffLatticeInfo() const
{
    return new ExtrapolatedGeneralizedOffLatticeInfo3D;
}

template <typename T, template <typename U> class Descriptor>
Array<T, 3> ExtrapolatedGeneralizedOffLatticeModel3D<T, Descriptor>::getLocalForce(
    AtomicContainerBlock3D &container) const
{
    ExtrapolatedGeneralizedOffLatticeInfo3D *info =
        dynamic_cast<ExtrapolatedGeneralizedOffLatticeInfo3D *>(container.getData());
    PLB_ASSERT(info);
    return info->getLocalForce();
}

template <typename T, template <typename U> class Descriptor>
void ExtrapolatedGeneralizedOffLatticeModel3D<T, Descriptor>::boundaryCompletion(
    AtomicBlock3D &nonTypeLattice, AtomicContainerBlock3D &container,
    std::vector<AtomicBlock3D *> const &args)
{
    BlockLattice3D<T, Descriptor> &lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> &>(nonTypeLattice);
    ExtrapolatedGeneralizedOffLatticeInfo3D *info =
        dynamic_cast<ExtrapolatedGeneralizedOffLatticeInfo3D *>(container.getData());
    PLB_ASSERT(info);
    std::vector<Dot3D> const &dryNodes = info->getDryNodes();
    std::vector<std::vector<std::pair<int, int> > > const &dryNodeFluidDirections =
        info->getDryNodeFluidDirections();
    std::vector<std::vector<int> > const &dryNodeFluidNoSolidDirections =
        info->getDryNodeFluidWithFluidDirections();
    std::vector<std::vector<plint> > const &dryNodeIds = info->getDryNodeIds();
    PLB_ASSERT(dryNodes.size() == dryNodeFluidDirections.size());

    Dot3D absoluteOffset = container.getLocation();

    Array<T, 3> &localForce = info->getLocalForce();
    localForce.resetToZero();
    for (pluint iDry = 0; iDry < dryNodes.size(); ++iDry) {
        cellCompletion(
            lattice, dryNodes[iDry], dryNodeFluidDirections[iDry],
            dryNodeFluidNoSolidDirections[iDry], dryNodeIds[iDry], absoluteOffset, localForce,
            args);
    }
}

template <typename T, template <typename U> class Descriptor>
void ExtrapolatedGeneralizedOffLatticeModel3D<T, Descriptor>::cellCompletion(
    BlockLattice3D<T, Descriptor> &lattice, Dot3D const &guoNode,
    std::vector<std::pair<int, int> > const &dryNodeFluidDirections,
    std::vector<int> const &dryNodeFluidNoSolidDirections, std::vector<plint> const &dryNodeIds,
    Dot3D const &absoluteOffset, Array<T, 3> &localForce, std::vector<AtomicBlock3D *> const &args)
{
    ExtrapolationGeneralizedAlgorithm3D<T, Descriptor> *algorithm =
        new ExtrapolationGeneralizedAlgorithm3D<T, Descriptor>(
            *this, lattice, guoNode, dryNodeFluidDirections, dryNodeFluidNoSolidDirections,
            dryNodeIds, absoluteOffset, localForce, args, this->computesStat(),
            this->usesSecondOrder());

    bool ok = algorithm->computeNeighborData();
    if (!ok) {
        global::plbErrors().registerError(
            "Error treating the geometry in the extrapolated generalized off-lattice model.");
    }
    algorithm->finalize();
    delete algorithm;
}

// ========================================================================= //
// ====Generalized BC with the velocity ocmputed with an interpolation====== //
// ========================================================================= //

template <typename T, template <typename U> class Descriptor>
class InterpolationAlgorithm3D {
public:
    typedef Descriptor<T> D;
    InterpolationAlgorithm3D(
        OffLatticeModel3D<T, Array<T, 3> > &model_, BlockLattice3D<T, Descriptor> &lattice_,
        Dot3D const &intNode_,
        std::vector<std::pair<int, int> > const &wetNodeSolidUsableDirections_,
        std::vector<plint> const &wetNodeIds_, std::vector<int> const &wetNodeFluidDirections_,
        std::vector<int> const &solidDirections_, Dot3D const &absoluteOffset_,
        Array<T, 3> &localForce_, std::vector<AtomicBlock3D *> const &args_, bool computeStat_,
        bool secondOrder_);
    virtual ~InterpolationAlgorithm3D() { }
    virtual bool computeNeighborData();
    void finalize();

    virtual void interpolateVariables(
        Dot3D const &solidDirection, int depth, Array<T, 3> const &wallNode, T wallDistance,
        Array<T, 3> const &wall_vel, OffBoundary::Type bdType, Array<T, 3> const &wallNormal,
        plint triangleId, plint iDirection) = 0;
    virtual void reduceVariables(T sumWeights) = 0;
    virtual void complete() = 0;

protected:
    OffLatticeModel3D<T, Array<T, 3> > &model;
    BlockLattice3D<T, Descriptor> &lattice;
    Dot3D const &intNode;
    Cell<T, Descriptor> &cell;
    std::vector<std::pair<int, int> > const &wetNodeSolidUsableDirections;
    std::vector<plint> const &wetNodeIds;
    std::vector<int> const &wetNodeFluidDirections;
    std::vector<int> const &solidDirections;
    Dot3D absoluteOffset;
    Array<T, 3> &localForce;
    std::vector<AtomicBlock3D *> const &args;

    plint numDirections;
    std::vector<T> weights;
    std::vector<Array<T, Descriptor<T>::d> > uVect;
    std::vector<Array<T, Descriptor<T>::d> > normalVect;

    Array<T, Descriptor<T>::d> u;
    bool computeStat, secondOrder;
    Array<T, Descriptor<T>::d> normal;
};

template <typename T, template <typename U> class Descriptor>
InterpolationAlgorithm3D<T, Descriptor>::InterpolationAlgorithm3D(
    OffLatticeModel3D<T, Array<T, 3> > &model_, BlockLattice3D<T, Descriptor> &lattice_,
    Dot3D const &intNode_, std::vector<std::pair<int, int> > const &wetNodeSolidUsableDirections_,
    std::vector<plint> const &wetNodeIds_, std::vector<int> const &wetNodeFluidDirections_,
    std::vector<int> const &solidDirections_, Dot3D const &absoluteOffset_,
    Array<T, 3> &localForce_, std::vector<AtomicBlock3D *> const &args_, bool computeStat_,
    bool secondOrder_) :
    model(model_),
    lattice(lattice_),
    intNode(intNode_),
    cell(lattice.get(intNode.x, intNode.y, intNode.z)),
    wetNodeSolidUsableDirections(wetNodeSolidUsableDirections_),
    wetNodeIds(wetNodeIds_),
    wetNodeFluidDirections(wetNodeFluidDirections_),
    solidDirections(solidDirections_),
    absoluteOffset(absoluteOffset_),
    localForce(localForce_),
    args(args_),
    computeStat(computeStat_),
    secondOrder(secondOrder_)
{
    numDirections = (plint)wetNodeSolidUsableDirections.size();
    weights.resize(numDirections);
    uVect.resize(numDirections);
    normalVect.resize(numDirections);
}

template <typename T, template <typename U> class Descriptor>
bool InterpolationAlgorithm3D<T, Descriptor>::computeNeighborData()
{
    using namespace indexTemplates;

    T sumWeights = T();
    this->u.resetToZero();
    Array<T, 3> wallNormal;
    for (plint iDirection = 0; iDirection < numDirections; ++iDirection) {
        int iNeighbor = wetNodeSolidUsableDirections[iDirection].first;
        int depth = wetNodeSolidUsableDirections[iDirection].second;
        int const *c = NextNeighbor<T>::c[iNeighbor];

        Dot3D solidDirection(c[0], c[1], c[2]);
        plint wetNodeId = wetNodeIds[iDirection];

        Array<T, 3> wallNode, wall_vel;
        T wallDistance;
        OffBoundary::Type bdType;
        bool ok = this->model.pointOnSurface(
            intNode + absoluteOffset, solidDirection, wallNode, wallDistance, wallNormal, wall_vel,
            bdType, wetNodeId);

        if (!ok) {
            global::plbErrors().registerError(
                "Interpolated off-lattice model could not find an intersection with a triangle.");
        }
        if (!(bdType == OffBoundary::dirichlet || bdType == OffBoundary::neumann
              || bdType == OffBoundary::freeSlip || bdType == OffBoundary::constRhoInlet
              || bdType == OffBoundary::densityNeumann))
        {
            return false;
        }

        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            // Use the formula uLB = uP - 1/2 g. If there is no external force,
            //   the force term automatically evaluates to zero.
            wall_vel[iD] -= (T)0.5 * getExternalForceComponent(cell, iD);
        }
        if (!(wallDistance <= NextNeighbor<T>::d[iNeighbor])) {
            global::plbErrors().registerError(
                "Error treating the geometry in the interpolated off-lattice model.");
        }

        this->interpolateVariables(
            solidDirection, depth, wallNode, wallDistance, wall_vel, bdType, wallNormal, wetNodeId,
            iDirection);
    }

    this->reduceVariables(sumWeights);

    return true;
}

template <typename T, template <typename U> class Descriptor>
void InterpolationAlgorithm3D<T, Descriptor>::finalize()
{
    Array<T, D::d> deltaJ;
    deltaJ.resetToZero();
    if (computeStat) {
        for (plint iDirection = 0; iDirection < (plint)solidDirections.size(); ++iDirection) {
            int iNeighbor = solidDirections[iDirection];
            int iPop = nextNeighborPop<T, Descriptor>(iNeighbor);
            if (iPop >= 0) {
                deltaJ[0] += D::c[iPop][0] * cell[iPop];
                deltaJ[1] += D::c[iPop][1] * cell[iPop];
                deltaJ[2] += D::c[iPop][2] * cell[iPop];
            }
        }
    }

    this->complete();

    if (computeStat) {
        Cell<T, Descriptor> collidedCell(cell);
        BlockStatistics statsCopy(lattice.getInternalStatistics());
        collidedCell.collide(statsCopy);

        for (plint iDirection = 0; iDirection < (plint)solidDirections.size(); ++iDirection) {
            int iNeighbor = solidDirections[iDirection];
            plint iPop = nextNeighborPop<T, Descriptor>(iNeighbor);
            if (iPop >= 0) {
                plint oppPop = indexTemplates::opposite<D>(iPop);

                deltaJ[0] -= D::c[oppPop][0] * collidedCell[oppPop];
                deltaJ[1] -= D::c[oppPop][1] * collidedCell[oppPop];
                deltaJ[2] -= D::c[oppPop][2] * collidedCell[oppPop];
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

// ================ InterpolationGeneralizedAlgorithm3D ============================= //

template <typename T, template <typename U> class Descriptor>
class InterpolationGeneralizedAlgorithm3D : public InterpolationAlgorithm3D<T, Descriptor> {
public:
    typedef Descriptor<T> D;
    InterpolationGeneralizedAlgorithm3D(
        OffLatticeModel3D<T, Array<T, 3> > &model_, BlockLattice3D<T, Descriptor> &lattice_,
        Dot3D const &intNode_,
        std::vector<std::pair<int, int> > const &wetNodeSolidUsableDirections_,
        std::vector<int> const &wetNodeFluidDirections_, std::vector<plint> const &wetNodeIds_,
        std::vector<plint> const &allWetNodeIds_, std::vector<int> const &solidDirections_,
        Dot3D const &absoluteOffset_, Array<T, 3> &localForce_,
        std::vector<AtomicBlock3D *> const &args_, bool computeStat_, bool secondOrder_);
    virtual void interpolateVariables(
        Dot3D const &fluidDirection, int depth, Array<T, 3> const &wallNode, T wallDistance,
        Array<T, 3> const &wall_vel, OffBoundary::Type bdType, Array<T, 3> const &wallNormal,
        plint triangleId, plint iDirection);
    virtual void reduceVariables(T sumWeights);
    virtual void complete();

private:
    std::vector<int> wetNodeFluidDirections;
    std::vector<plint> allWetNodeIds;
    T totWeights;
};

template <typename T, template <typename U> class Descriptor>
InterpolationGeneralizedAlgorithm3D<T, Descriptor>::InterpolationGeneralizedAlgorithm3D(
    OffLatticeModel3D<T, Array<T, 3> > &model_, BlockLattice3D<T, Descriptor> &lattice_,
    Dot3D const &intNode_, std::vector<std::pair<int, int> > const &wetNodeSolidUsableDirections_,
    std::vector<int> const &wetNodeFluidDirections_, std::vector<plint> const &wetNodeIds_,
    std::vector<plint> const &allWetNodeIds_, std::vector<int> const &solidDirections_,
    Dot3D const &absoluteOffset_, Array<T, 3> &localForce_,
    std::vector<AtomicBlock3D *> const &args_, bool computeStat_, bool secondOrder_) :
    InterpolationAlgorithm3D<T, Descriptor>(
        model_, lattice_, intNode_, wetNodeSolidUsableDirections_, wetNodeIds_,
        wetNodeFluidDirections_, solidDirections_, absoluteOffset_, localForce_, args_,
        computeStat_, secondOrder_),
    wetNodeFluidDirections(wetNodeFluidDirections_),
    allWetNodeIds(allWetNodeIds_),
    totWeights(T())
{
    // pcout << this->solidDirections.size() << ", " << this->wetNodeSolidUsableDirections.size() <<
    // std::endl; pcout << this->wetNodeIds.size() << std::endl;
}

template <typename T, template <typename U> class Descriptor>
void InterpolationGeneralizedAlgorithm3D<T, Descriptor>::interpolateVariables(
    Dot3D const &solidDirection, int depth, Array<T, 3> const &wallNode, T wallDistance,
    Array<T, 3> const &wall_vel, OffBoundary::Type bdType, Array<T, 3> const &wallNormal,
    plint triangleId, plint iDirection)
{
    Array<T, Descriptor<T>::d> u1, u2;
    if (depth == 0) {
        plbWarning(true, "Warning depth 0 detected in one direction.");
        return;
    } else if (depth == 1 || !this->secondOrder) {
        Dot3D pos = this->intNode - solidDirection;
        if (this->args.empty()) {
            Cell<T, Descriptor> const &cell1 = this->lattice.get(pos.x, pos.y, pos.z);

            if (!this->model.velIsJ()) {
                cell1.getDynamics().computeVelocity(cell1, u1);
            } else {
                T rhoBar;
                cell1.getDynamics().computeRhoBarJ(cell1, rhoBar, u1);
            }
        } else {
            if ((plint)this->args.size() == 1) {
                NTensorField3D<T> const *macroField =
                    dynamic_cast<NTensorField3D<T> const *>(this->args[0]);
                PLB_ASSERT(macroField);
                // 1 Variable for rhoBar, 3 variables for j.
                PLB_ASSERT(macroField->getNdim() == 4);
                Dot3D macroPos = pos + computeRelativeDisplacement(this->lattice, *macroField);
                T const *macroscopic = macroField->get(macroPos.x, macroPos.y, macroPos.z);
                u1.from_cArray(macroscopic + 1);
                if (!this->model.velIsJ()) {
                    u1 *= Descriptor<T>::invRho(macroscopic[0]);
                }
            } else if ((plint)this->args.size() == 2) {
                // 1 field for rhoBar, 1 field for j.
                // #ifdef PLB_DEBUG
                ScalarField3D<T> const *rhoBarField =
                    dynamic_cast<ScalarField3D<T> const *>(this->args[0]);
                // #endif
                TensorField3D<T, 3> const *jField =
                    dynamic_cast<TensorField3D<T, 3> const *>(this->args[1]);
                PLB_ASSERT(rhoBarField);
                PLB_ASSERT(jField);
                Dot3D posJ = pos + computeRelativeDisplacement(this->lattice, *jField);
                Dot3D posRho = pos + computeRelativeDisplacement(this->lattice, *rhoBarField);
                u1 = jField->get(posJ.x, posJ.y, posJ.z);
                if (!this->model.velIsJ()) {
                    u1 *= Descriptor<T>::invRho(rhoBarField->get(posRho.x, posRho.y, posRho.z));
                }
            } else {
                PLB_ASSERT(false);  // Not implemented for 3 args.
            }
        }
    } else if (depth >= 2) {  // depth >= 2
        Dot3D pos = this->intNode - solidDirection;
        Dot3D pos2 = pos - solidDirection;

        if (this->args.empty()) {
            Cell<T, Descriptor> const &cell1 = this->lattice.get(pos.x, pos.y, pos.z);
            Cell<T, Descriptor> const &cell2 = this->lattice.get(pos2.x, pos2.y, pos2.z);
            if (!this->model.velIsJ()) {
                cell1.getDynamics().computeVelocity(cell1, u1);
                cell2.getDynamics().computeVelocity(cell2, u2);
            } else {
                T rhoBar;
                cell1.getDynamics().computeRhoBarJ(cell1, rhoBar, u1);
                cell2.getDynamics().computeRhoBarJ(cell2, rhoBar, u2);
            }
            Array<T, 3> r(
                -2 * (T)solidDirection.x, -2 * (T)solidDirection.y, -2 * (T)solidDirection.z);
            T W = inamuroDeltaFunction<T>().W(r);
            this->u += W * u2;
            totWeights += W;

        } else {
            if ((plint)this->args.size() == 1) {
                NTensorField3D<T> const *macroField =
                    dynamic_cast<NTensorField3D<T> const *>(this->args[0]);
                PLB_ASSERT(macroField);
                // 1 Variable for rhoBar, 3 variables for j.
                PLB_ASSERT(macroField->getNdim() == 4);
                Dot3D off = computeRelativeDisplacement(this->lattice, *macroField);
                Dot3D macroPos = pos + off;
                Dot3D macroPos2 = pos2 + off;
                T const *macroscopic = macroField->get(macroPos.x, macroPos.y, macroPos.z);
                u1.from_cArray(macroscopic + 1);
                T const *macroscopic2 = macroField->get(macroPos2.x, macroPos2.y, macroPos2.z);
                u2.from_cArray(macroscopic2 + 1);
                if (!this->model.velIsJ()) {
                    u1 *= Descriptor<T>::invRho(macroscopic[0]);
                    u2 *= Descriptor<T>::invRho(macroscopic2[0]);
                }
                Array<T, 3> r(
                    -2 * (T)solidDirection.x, -2 * (T)solidDirection.y, -2 * (T)solidDirection.z);
                T W = inamuroDeltaFunction<T>().W(r);
                this->u += W * u2;
                totWeights += W;
            } else if ((plint)this->args.size() == 2) {
                // 1 field for rhoBar, 1 field for j.
                // #ifdef PLB_DEBUG
                ScalarField3D<T> const *rhoBarField =
                    dynamic_cast<ScalarField3D<T> const *>(this->args[0]);
                // #endif
                TensorField3D<T, 3> const *jField =
                    dynamic_cast<TensorField3D<T, 3> const *>(this->args[1]);
                PLB_ASSERT(rhoBarField);
                PLB_ASSERT(jField);
                Dot3D posJ = pos + computeRelativeDisplacement(this->lattice, *jField);
                Dot3D posRho = pos + computeRelativeDisplacement(this->lattice, *rhoBarField);
                u1 = jField->get(posJ.x, posJ.y, posJ.z);
                Dot3D posJ2 = pos2 + computeRelativeDisplacement(this->lattice, *jField);
                Dot3D posRho2 = pos2 + computeRelativeDisplacement(this->lattice, *rhoBarField);
                u2 = jField->get(posJ2.x, posJ2.y, posJ2.z);
                if (!this->model.velIsJ()) {
                    u1 *= Descriptor<T>::invRho(rhoBarField->get(posRho.x, posRho.y, posRho.z));
                    u2 *= Descriptor<T>::invRho(rhoBarField->get(posRho2.x, posRho2.y, posRho2.z));
                }
                Array<T, 3> r(
                    -2 * (T)solidDirection.x, -2 * (T)solidDirection.y, -2 * (T)solidDirection.z);
                T W = inamuroDeltaFunction<T>().W(r);
                this->u += W * u2;
                totWeights += W;

            } else {
                PLB_ASSERT(false);  // Not implemented for 3 args.
            }
        }
    } else {
        PLB_ASSERT(false && "In Interpolated Generalized BC ");
    }

    Array<T, 3> r(-(T)solidDirection.x, -(T)solidDirection.y, -(T)solidDirection.z);
    T W = inamuroDeltaFunction<T>().W(r);
    this->u += W * u1;
    totWeights += W;
}

template <typename T, template <typename U> class Descriptor>
void InterpolationGeneralizedAlgorithm3D<T, Descriptor>::reduceVariables(T sumWeights)
{
    Array<T, 3> wallNormal;
    for (plint iDirection = 0; iDirection < (plint)this->solidDirections.size(); ++iDirection) {
        int const *c = NextNeighbor<T>::c[this->solidDirections[iDirection]];

        Dot3D solidDirection(c[0], c[1], c[2]);
        plint wetNodeId = allWetNodeIds[iDirection];

        Array<T, 3> wallNode, wall_vel;
        T wallDistance;
        OffBoundary::Type bdType;
        bool ok = this->model.pointOnSurface(
            this->intNode + this->absoluteOffset, solidDirection, wallNode, wallDistance,
            wallNormal, wall_vel, bdType, wetNodeId);
        if (!ok) {
            global::plbErrors().registerError(
                "Interpolated generalized off-lattice model could not find an intersection with a "
                "triangle.");
        }

        Array<T, 3> r((T)solidDirection.x, (T)solidDirection.y, (T)solidDirection.z);
        r /= norm(r);
        r *= wallDistance;
        T W = inamuroDeltaFunction<T>().W(r);
        this->u += W * wall_vel;
        totWeights += W;
    }
    this->u /= totWeights;
}

template <typename T, template <typename U> class Descriptor>
void InterpolationGeneralizedAlgorithm3D<T, Descriptor>::complete()
{
    std::vector<plint> knownIndices;
    knownIndices.push_back(0);
    for (pluint iDirection = 0; iDirection < wetNodeFluidDirections.size(); ++iDirection) {
        int iPop = wetNodeFluidDirections[iDirection];
        plint index = indexTemplates::opposite<Descriptor<T> >(iPop);
        knownIndices.push_back(index);
    }
    PLB_ASSERT(knownIndices.size() >= 6);

    std::vector<plint> missingIndices =
        indexTemplates::remainingIndexes<Descriptor<T> >(knownIndices);

    Dynamics<T, Descriptor> const &dynamics = this->cell.getDynamics();
    DirichletVelocityBoundarySolver<T, Descriptor> bc(missingIndices, knownIndices, this->u);
    bc.apply(this->cell, dynamics, !this->model.getPartialReplace());
}

// ================ InterpolationDefineVelocityOffLatticeModel3D ============================= //

template <typename T, template <typename U> class Descriptor>
class InterpolationDefineVelocityOffLatticeModel3D :
    public InterpolationAlgorithm3D<T, Descriptor> {
public:
    typedef Descriptor<T> D;
    InterpolationDefineVelocityOffLatticeModel3D(
        OffLatticeModel3D<T, Array<T, 3> > &model_, BlockLattice3D<T, Descriptor> &lattice_,
        Dot3D const &intNode_,
        std::vector<std::pair<int, int> > const &wetNodeSolidUsableDirections_,
        std::vector<int> const &wetNodeFluidDirections_, std::vector<plint> const &wetNodeIds_,
        std::vector<int> const &solidDirections_, Dot3D const &absoluteOffset_,
        Array<T, 3> &localForce_, std::vector<AtomicBlock3D *> const &args_, bool computeStat_,
        bool secondOrder_);
    virtual void interpolateVariables(
        Dot3D const &fluidDirection, int depth, Array<T, 3> const &wallNode, T wallDistance,
        Array<T, 3> const &wall_vel, OffBoundary::Type bdType, Array<T, 3> const &wallNormal,
        plint triangleId, plint iDirection);
    virtual void reduceVariables(T sumWeights);
    virtual void complete();

private:
    std::vector<int> wetNodeFluidDirections;

    T rhoBar;
    T weightsU;
    T weightsRho;
};

template <typename T, template <typename U> class Descriptor>
InterpolationDefineVelocityOffLatticeModel3D<T, Descriptor>::
    InterpolationDefineVelocityOffLatticeModel3D(
        OffLatticeModel3D<T, Array<T, 3> > &model_, BlockLattice3D<T, Descriptor> &lattice_,
        Dot3D const &intNode_,
        std::vector<std::pair<int, int> > const &wetNodeSolidUsableDirections_,
        std::vector<int> const &wetNodeFluidDirections_, std::vector<plint> const &wetNodeIds_,
        std::vector<int> const &solidDirections_, Dot3D const &absoluteOffset_,
        Array<T, 3> &localForce_, std::vector<AtomicBlock3D *> const &args_, bool computeStat_,
        bool secondOrder_) :
    InterpolationAlgorithm3D<T, Descriptor>(
        model_, lattice_, intNode_, wetNodeSolidUsableDirections_, wetNodeIds_,
        wetNodeFluidDirections_, solidDirections_, absoluteOffset_, localForce_, args_,
        computeStat_, secondOrder_),
    wetNodeFluidDirections(wetNodeFluidDirections_)
{
    rhoBar = T();
    weightsU = T();
    weightsRho = T();
}

template <typename T, template <typename U> class Descriptor>
void InterpolationDefineVelocityOffLatticeModel3D<T, Descriptor>::interpolateVariables(
    Dot3D const &solidDirection, int depth, Array<T, 3> const &wallNode, T wallDistance,
    Array<T, 3> const &wall_vel, OffBoundary::Type bdType, Array<T, 3> const &wallNormal,
    plint triangleId, plint iDirection)
{
    Array<T, Descriptor<T>::d> u1, u2;
    T rhoBar1, rhoBar2;
    if (depth == 0) {
        plbWarning(true, "Warning depth 0 detected in one direction.");
        return;
    } else if (depth == 1 || !this->secondOrder) {
        Dot3D pos = this->intNode - solidDirection;
        if (this->args.empty()) {
            Cell<T, Descriptor> const &cell1 = this->lattice.get(pos.x, pos.y, pos.z);

            if (!this->model.velIsJ()) {
                cell1.getDynamics().computeRhoBarJ(cell1, rhoBar1, u1);
                u1 *= Descriptor<T>::invRho(rhoBar1);
            } else {
                cell1.getDynamics().computeRhoBarJ(cell1, rhoBar1, u1);
            }
        } else {
            if ((plint)this->args.size() == 1) {
                NTensorField3D<T> const *macroField =
                    dynamic_cast<NTensorField3D<T> const *>(this->args[0]);
                PLB_ASSERT(macroField);
                // 1 Variable for rhoBar, 3 variables for j.
                PLB_ASSERT(macroField->getNdim() == 4);
                Dot3D macroPos = pos + computeRelativeDisplacement(this->lattice, *macroField);
                T const *macroscopic = macroField->get(macroPos.x, macroPos.y, macroPos.z);
                rhoBar1 = macroscopic[0];
                u1.from_cArray(macroscopic + 1);
                if (!this->model.velIsJ()) {
                    u1 *= Descriptor<T>::invRho(rhoBar1);
                }
            } else if ((plint)this->args.size() == 2) {
                // 1 field for rhoBar, 1 field for j.
                // #ifdef PLB_DEBUG
                ScalarField3D<T> const *rhoBarField =
                    dynamic_cast<ScalarField3D<T> const *>(this->args[0]);
                // #endif
                TensorField3D<T, 3> const *jField =
                    dynamic_cast<TensorField3D<T, 3> const *>(this->args[1]);
                PLB_ASSERT(rhoBarField);
                PLB_ASSERT(jField);
                Dot3D posJ = pos + computeRelativeDisplacement(this->lattice, *jField);
                Dot3D posRho = pos + computeRelativeDisplacement(this->lattice, *rhoBarField);
                rhoBar1 = rhoBarField->get(posRho.x, posRho.y, posRho.z);
                u1 = jField->get(posJ.x, posJ.y, posJ.z);
                if (!this->model.velIsJ()) {
                    u1 *= Descriptor<T>::invRho(rhoBar1);
                }
            } else {
                PLB_ASSERT(false);  // Not implemented for 3 args.
            }
        }
    } else if (depth >= 2) {  // depth >= 2
        Dot3D pos = this->intNode - solidDirection;
        Dot3D pos2 = pos - solidDirection;

        if (this->args.empty()) {
            Cell<T, Descriptor> const &cell1 = this->lattice.get(pos.x, pos.y, pos.z);
            Cell<T, Descriptor> const &cell2 = this->lattice.get(pos2.x, pos2.y, pos2.z);
            if (!this->model.velIsJ()) {
                cell1.getDynamics().computeRhoBarJ(cell1, rhoBar1, u1);
                cell2.getDynamics().computeRhoBarJ(cell2, rhoBar2, u2);
                u1 *= Descriptor<T>::invRho(rhoBar1);
                u2 *= Descriptor<T>::invRho(rhoBar2);
            } else {
                cell1.getDynamics().computeRhoBarJ(cell1, rhoBar1, u1);
                cell2.getDynamics().computeRhoBarJ(cell2, rhoBar2, u2);
            }
            Array<T, 3> r(
                -2 * (T)solidDirection.x, -2 * (T)solidDirection.y, -2 * (T)solidDirection.z);
            T W = inamuroDeltaFunction<T>().W(r);
            this->u += W * u2;
            rhoBar += W * rhoBar2;
            weightsU += W;
            weightsRho += W;

        } else {
            if ((plint)this->args.size() == 1) {
                NTensorField3D<T> const *macroField =
                    dynamic_cast<NTensorField3D<T> const *>(this->args[0]);
                PLB_ASSERT(macroField);
                // 1 Variable for rhoBar, 3 variables for j.
                PLB_ASSERT(macroField->getNdim() == 4);
                Dot3D off = computeRelativeDisplacement(this->lattice, *macroField);
                Dot3D macroPos = pos + off;
                Dot3D macroPos2 = pos2 + off;
                T const *macroscopic = macroField->get(macroPos.x, macroPos.y, macroPos.z);
                rhoBar1 = macroscopic[0];
                u1.from_cArray(macroscopic + 1);
                T const *macroscopic2 = macroField->get(macroPos2.x, macroPos2.y, macroPos2.z);
                rhoBar2 = macroscopic2[0];
                u2.from_cArray(macroscopic2 + 1);
                if (!this->model.velIsJ()) {
                    u1 *= Descriptor<T>::invRho(rhoBar1);
                    u2 *= Descriptor<T>::invRho(rhoBar2);
                }
                Array<T, 3> r(
                    -2 * (T)solidDirection.x, -2 * (T)solidDirection.y, -2 * (T)solidDirection.z);
                T W = inamuroDeltaFunction<T>().W(r);
                this->u += W * u2;
                rhoBar += W * rhoBar2;
                weightsU += W;
                weightsRho += W;
            } else if ((plint)this->args.size() == 2) {
                // 1 field for rhoBar, 1 field for j.
                // #ifdef PLB_DEBUG
                ScalarField3D<T> const *rhoBarField =
                    dynamic_cast<ScalarField3D<T> const *>(this->args[0]);
                // #endif
                TensorField3D<T, 3> const *jField =
                    dynamic_cast<TensorField3D<T, 3> const *>(this->args[1]);
                PLB_ASSERT(rhoBarField);
                PLB_ASSERT(jField);
                Dot3D posJ = pos + computeRelativeDisplacement(this->lattice, *jField);
                Dot3D posRho = pos + computeRelativeDisplacement(this->lattice, *rhoBarField);
                rhoBar1 = rhoBarField->get(posRho.x, posRho.y, posRho.z);
                u1 = jField->get(posJ.x, posJ.y, posJ.z);
                Dot3D posJ2 = pos2 + computeRelativeDisplacement(this->lattice, *jField);
                Dot3D posRho2 = pos2 + computeRelativeDisplacement(this->lattice, *rhoBarField);
                rhoBar2 = rhoBarField->get(posRho2.x, posRho2.y, posRho2.z);
                u2 = jField->get(posJ2.x, posJ2.y, posJ2.z);
                if (!this->model.velIsJ()) {
                    u1 *= Descriptor<T>::invRho(rhoBar1);
                    u2 *= Descriptor<T>::invRho(rhoBar2);
                }
                Array<T, 3> r(
                    -2 * (T)solidDirection.x, -2 * (T)solidDirection.y, -2 * (T)solidDirection.z);
                T W = inamuroDeltaFunction<T>().W(r);
                rhoBar += W * rhoBar2;
                this->u += W * u2;
                weightsU += W;
                weightsRho += W;
            } else {
                PLB_ASSERT(false);  // Not implemented for 3 args.
            }
        }
    } else {
        PLB_ASSERT(false && "In Interpolated Generalized BC ");
    }

    Array<T, 3> r(-(T)solidDirection.x, -(T)solidDirection.y, -(T)solidDirection.z);
    T W = inamuroDeltaFunction<T>().W(r);
    this->u += W * u1;
    rhoBar += W * rhoBar1;
    weightsU += W;
    weightsRho += W;

    r = Array<T, 3>((T)solidDirection.x, (T)solidDirection.y, (T)solidDirection.z);
    r /= norm(r);
    r *= wallDistance;

    W = inamuroDeltaFunction<T>().W(r);
    this->u += W * wall_vel;
    weightsU += W;
}

template <typename T, template <typename U> class Descriptor>
void InterpolationDefineVelocityOffLatticeModel3D<T, Descriptor>::reduceVariables(T sumWeights)
{
    this->u /= weightsRho;
    rhoBar /= weightsRho;
}

template <typename T, template <typename U> class Descriptor>
void InterpolationDefineVelocityOffLatticeModel3D<T, Descriptor>::complete()
{
    Array<T, Descriptor<T>::d> j = Descriptor<T>::fullRho(rhoBar) * this->u;
    if (this->args.empty()) {
        Cell<T, Descriptor> &cell =
            this->lattice.get(this->intNode.x, this->intNode.y, this->intNode.z);

        cell.getDynamics().computeEquilibria(cell.getRawPopulations(), rhoBar, j, normSqr(j));
    } else {
        if ((plint)this->args.size() == 1) {
            NTensorField3D<T> *macroField = dynamic_cast<NTensorField3D<T> *>(this->args[0]);
            PLB_ASSERT(macroField);
            // 1 Variable for rhoBar, 3 variables for j.
            PLB_ASSERT(macroField->getNdim() == 4);
            Dot3D offset = computeRelativeDisplacement(this->lattice, *macroField);
            T *macroscopic = macroField->get(
                this->intNode.x + offset.x, this->intNode.y + offset.y, this->intNode.z + offset.z);
            // pcout << "rhoBar = " << this->rhoBar << ", ux = " << this->j[0] << ", uy = " <<
            // this->j[1] << ", uz = " << this->j[2] << std::endl;
            macroscopic[0] = rhoBar;
            j.to_cArray(macroscopic + 1);
        } else if ((plint)this->args.size() == 2) {
            // 1 field for rhoBar, 1 field for j.
            ScalarField3D<T> *rhoBarField = dynamic_cast<ScalarField3D<T> *>(this->args[0]);

            TensorField3D<T, 3> *jField = dynamic_cast<TensorField3D<T, 3> *>(this->args[1]);
            PLB_ASSERT(rhoBarField);
            PLB_ASSERT(jField);

            Dot3D offsetScalar = computeRelativeDisplacement(this->lattice, *rhoBarField);
            rhoBarField->get(
                this->intNode.x + offsetScalar.x, this->intNode.y + offsetScalar.y,
                this->intNode.z + offsetScalar.z) = rhoBar;
            Dot3D offset = computeRelativeDisplacement(this->lattice, *jField);
            jField->get(
                this->intNode.x + offset.x, this->intNode.y + offset.y,
                this->intNode.z + offset.z) = j;
        } else {
            PLB_ASSERT(false);  // Not implemented for 3 args.
        }
    }
}

// ================ InterpolationFdCompletionAlgorithm3D ============================= //

template <typename T, template <typename U> class Descriptor>
class InterpolationFdCompletionAlgorithm3D : public InterpolationAlgorithm3D<T, Descriptor> {
public:
    typedef Descriptor<T> D;
    InterpolationFdCompletionAlgorithm3D(
        OffLatticeModel3D<T, Array<T, 3> > &model_, BlockLattice3D<T, Descriptor> &lattice_,
        Dot3D const &intNode_,
        std::vector<std::pair<int, int> > const &wetNodeSolidUsableDirections_,
        std::vector<int> const &wetNodeFluidDirections_, std::vector<plint> const &wetNodeIds_,
        std::vector<int> const &solidDirections_, Dot3D const &absoluteOffset_,
        Array<T, 3> &localForce_, std::vector<AtomicBlock3D *> const &args_, bool computeStat_,
        bool secondOrder_, std::pair<int, int> xDerivDirAndOrder_,
        std::pair<int, int> yDerivDirAndOrder_, std::pair<int, int> zDerivDirAndOrder_);
    virtual void interpolateVariables(
        Dot3D const &fluidDirection, int depth, Array<T, 3> const &wallNode, T wallDistance,
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
        Array<T, Descriptor<T>::d> &d_u) const;
    void computeCentralGradient(
        ScalarField3D<T> const *rhoBarField, TensorField3D<T, 3> const *jField, const Dot3D &pos,
        const Dot3D &pos2, const Dot3D &dx, Array<T, Descriptor<T>::d> &d_u) const;
    void computeFwdGradient(
        const BlockLattice3D<T, Descriptor> &lattice, const Array<T, Descriptor<T>::d> &u0,
        const Dot3D &pos, const Dot3D &dx, int orient, Array<T, Descriptor<T>::d> &d_u) const;
    void computeFwdGradient(
        NTensorField3D<T> const *macroField, const Array<T, Descriptor<T>::d> &u0, const Dot3D &pos,
        const Dot3D &dx, int orient, Array<T, Descriptor<T>::d> &d_u) const;
    void computeFwdGradient(
        ScalarField3D<T> const *rhoBarField, TensorField3D<T, 3> const *jField,
        const Array<T, Descriptor<T>::d> &u0, const Dot3D &pos, const Dot3D &pos2, const Dot3D &dx,
        int orient, Array<T, Descriptor<T>::d> &d_u) const;

private:
    std::vector<int> wetNodeFluidDirections;

    std::pair<int, int> xDerivDirAndOrder;
    std::pair<int, int> yDerivDirAndOrder;
    std::pair<int, int> zDerivDirAndOrder;

    T rhoBar;
};

template <typename T, template <typename U> class Descriptor>
InterpolationFdCompletionAlgorithm3D<T, Descriptor>::InterpolationFdCompletionAlgorithm3D(
    OffLatticeModel3D<T, Array<T, 3> > &model_, BlockLattice3D<T, Descriptor> &lattice_,
    Dot3D const &intNode_, std::vector<std::pair<int, int> > const &wetNodeSolidUsableDirections_,
    std::vector<int> const &wetNodeFluidDirections_, std::vector<plint> const &wetNodeIds_,
    std::vector<int> const &solidDirections_, Dot3D const &absoluteOffset_,
    Array<T, 3> &localForce_, std::vector<AtomicBlock3D *> const &args_, bool computeStat_,
    bool secondOrder_, std::pair<int, int> xDerivDirAndOrder_,
    std::pair<int, int> yDerivDirAndOrder_, std::pair<int, int> zDerivDirAndOrder_) :
    InterpolationAlgorithm3D<T, Descriptor>(
        model_, lattice_, intNode_, wetNodeSolidUsableDirections_, wetNodeIds_,
        wetNodeFluidDirections_, solidDirections_, absoluteOffset_, localForce_, args_,
        computeStat_, secondOrder_),
    wetNodeFluidDirections(wetNodeFluidDirections_),
    xDerivDirAndOrder(xDerivDirAndOrder_),
    yDerivDirAndOrder(yDerivDirAndOrder_),
    zDerivDirAndOrder(zDerivDirAndOrder_)

{ }

template <typename T, template <typename U> class Descriptor>
void InterpolationFdCompletionAlgorithm3D<T, Descriptor>::computeCentralGradient(
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
void InterpolationFdCompletionAlgorithm3D<T, Descriptor>::computeCentralGradient(
    NTensorField3D<T> const *macroField, const Dot3D &pos, const Dot3D &dx,
    Array<T, Descriptor<T>::d> &d_u) const
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
void InterpolationFdCompletionAlgorithm3D<T, Descriptor>::computeCentralGradient(
    ScalarField3D<T> const *rhoBarField, TensorField3D<T, 3> const *jField, const Dot3D &pos,
    const Dot3D &pos2, const Dot3D &dx, Array<T, Descriptor<T>::d> &d_u) const
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
void InterpolationFdCompletionAlgorithm3D<T, Descriptor>::computeFwdGradient(
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
void InterpolationFdCompletionAlgorithm3D<T, Descriptor>::computeFwdGradient(
    NTensorField3D<T> const *macroField, const Array<T, Descriptor<T>::d> &u0, const Dot3D &pos,
    const Dot3D &dx, int orient, Array<T, Descriptor<T>::d> &d_u) const
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
void InterpolationFdCompletionAlgorithm3D<T, Descriptor>::computeFwdGradient(
    ScalarField3D<T> const *rhoBarField, TensorField3D<T, 3> const *jField,
    const Array<T, Descriptor<T>::d> &u0, const Dot3D &pos, const Dot3D &pos2, const Dot3D &dx,
    int orient, Array<T, Descriptor<T>::d> &d_u) const
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
bool InterpolationFdCompletionAlgorithm3D<T, Descriptor>::computeNeighborData()
{
    reduceVariables((T)0);
    return true;
}

template <typename T, template <typename U> class Descriptor>
void InterpolationFdCompletionAlgorithm3D<T, Descriptor>::interpolateVariables(
    Dot3D const &solidDirection, int depth, Array<T, 3> const &wallNode, T wallDistance,
    Array<T, 3> const &wall_vel, OffBoundary::Type bdType, Array<T, 3> const &wallNormal,
    plint triangleId, plint iDirection)
{ }

template <typename T, template <typename U> class Descriptor>
void InterpolationFdCompletionAlgorithm3D<T, Descriptor>::reduceVariables(T sumWeights)
{
    if (this->args.empty()) {
        this->cell.getDynamics().computeRhoBarJ(this->cell, rhoBar, this->u);
    } else {
        if ((plint)this->args.size() == 1) {
            NTensorField3D<T> const *macroField =
                dynamic_cast<NTensorField3D<T> const *>(this->args[0]);
            PLB_ASSERT(macroField);
            // 1 Variable for rhoBar, 3 variables for j.
            PLB_ASSERT(macroField->getNdim() == 4);
            Dot3D offset = computeRelativeDisplacement(this->lattice, *macroField);
            T const *macro_0 = macroField->get(
                this->intNode.x + offset.x, this->intNode.y + offset.y, this->intNode.z + offset.z);
            rhoBar = macro_0[0];
            this->u.from_cArray(macro_0 + 1);
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
            rhoBar = rhoBarField->get(
                this->intNode.x + offsetScalar.x, this->intNode.y + offsetScalar.y,
                this->intNode.z + offsetScalar.z);

            this->u = jField->get(
                this->intNode.x + offset.x, this->intNode.y + offset.y, this->intNode.z + offset.z);
        } else {
            PLB_ASSERT(false);  // Not implemented for 3 args.
        }
    }
    if (!this->model.velIsJ()) {
        this->u *= Descriptor<T>::invRho(rhoBar);
    }
}

template <typename T, template <typename U> class Descriptor>
void InterpolationFdCompletionAlgorithm3D<T, Descriptor>::complete()
{
    Array<T, Descriptor<T>::d> dx_u, dy_u, dz_u;
    if (this->args.empty()) {
        if (xDerivDirAndOrder.second == 2) {
            computeCentralGradient(this->lattice, this->intNode, Dot3D(1, 0, 0), dx_u);
        } else if (xDerivDirAndOrder.second == 1) {
            PLB_ASSERT(
                (xDerivDirAndOrder.first == 1 || xDerivDirAndOrder.first == -1)
                && "Invalid xDerivDirAndOrder.first value (must be = 1 or = -1");
            computeFwdGradient(
                this->lattice, this->u, this->intNode, Dot3D(xDerivDirAndOrder.first, 0, 0),
                xDerivDirAndOrder.first, dx_u);
        } else if (xDerivDirAndOrder.second == 0) {
            dx_u.resetToZero();
        } else {
            PLB_ASSERT(false && "xDerivDirAndOrder has invalid order");
        }
        if (yDerivDirAndOrder.second == 2) {
            computeCentralGradient(this->lattice, this->intNode, Dot3D(0, 1, 0), dy_u);
        } else if (yDerivDirAndOrder.second == 1) {
            PLB_ASSERT(
                (yDerivDirAndOrder.first == 1 || yDerivDirAndOrder.first == -1)
                && "Invalid yDerivDirAndOrder.first value (must be = 1 or = -1");
            computeFwdGradient(
                this->lattice, this->u, this->intNode, Dot3D(0, yDerivDirAndOrder.first, 0),
                yDerivDirAndOrder.first, dy_u);
        } else if (yDerivDirAndOrder.second == 0) {
            dy_u.resetToZero();
        } else {
            PLB_ASSERT(false && "yDerivDirAndOrder has invalid order");
        }
        if (zDerivDirAndOrder.second == 2) {
            computeCentralGradient(this->lattice, this->intNode, Dot3D(0, 0, 1), dz_u);
        } else if (zDerivDirAndOrder.second == 1) {
            PLB_ASSERT(
                (zDerivDirAndOrder.first == 1 || zDerivDirAndOrder.first == -1)
                && "Invalid zDerivDirAndOrder.first value (must be = 1 or = -1");
            computeFwdGradient(
                this->lattice, this->u, this->intNode, Dot3D(0, 0, zDerivDirAndOrder.first),
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
            if (xDerivDirAndOrder.second == 2) {
                computeCentralGradient(macroField, this->intNode + offset, Dot3D(1, 0, 0), dx_u);
            } else if (xDerivDirAndOrder.second == 1) {
                PLB_ASSERT(
                    (xDerivDirAndOrder.first == 1 || xDerivDirAndOrder.first == -1)
                    && "Invalid xDerivDirAndOrder.first value (must be = 1 or = -1");
                computeFwdGradient(
                    macroField, this->u, this->intNode + offset,
                    Dot3D(xDerivDirAndOrder.first, 0, 0), xDerivDirAndOrder.first, dx_u);
            } else if (xDerivDirAndOrder.second == 0) {
                dx_u.resetToZero();
            } else {
                PLB_ASSERT(false && "xDerivDirAndOrder has invalid order");
            }
            if (yDerivDirAndOrder.second == 2) {
                computeCentralGradient(macroField, this->intNode + offset, Dot3D(0, 1, 0), dy_u);
            } else if (yDerivDirAndOrder.second == 1) {
                PLB_ASSERT(
                    (yDerivDirAndOrder.first == 1 || yDerivDirAndOrder.first == -1)
                    && "Invalid yDerivDirAndOrder.first value (must be = 1 or = -1");
                computeFwdGradient(
                    macroField, this->u, this->intNode + offset,
                    Dot3D(0, yDerivDirAndOrder.first, 0), yDerivDirAndOrder.first, dy_u);
            } else if (yDerivDirAndOrder.second == 0) {
                dy_u.resetToZero();
            } else {
                PLB_ASSERT(false && "yDerivDirAndOrder has invalid order");
            }
            if (zDerivDirAndOrder.second == 2) {
                computeCentralGradient(macroField, this->intNode + offset, Dot3D(0, 0, 1), dz_u);
            } else if (zDerivDirAndOrder.second == 1) {
                PLB_ASSERT(
                    (zDerivDirAndOrder.first == 1 || zDerivDirAndOrder.first == -1)
                    && "Invalid zDerivDirAndOrder.first value (must be = 1 or = -1");
                computeFwdGradient(
                    macroField, this->u, this->intNode + offset,
                    Dot3D(0, 0, zDerivDirAndOrder.first), zDerivDirAndOrder.first, dz_u);
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
                    rhoBarField, jField, this->intNode + offset, this->intNode + offsetScalar,
                    Dot3D(1, 0, 0), dx_u);
            } else if (xDerivDirAndOrder.second == 1) {
                PLB_ASSERT(
                    (xDerivDirAndOrder.first == 1 || xDerivDirAndOrder.first == -1)
                    && "Invalid xDerivDirAndOrder.first value (must be = 1 or = -1");
                computeFwdGradient(
                    rhoBarField, jField, this->u, this->intNode + offset,
                    this->intNode + offsetScalar, Dot3D(xDerivDirAndOrder.first, 0, 0),
                    xDerivDirAndOrder.first, dx_u);
            } else if (xDerivDirAndOrder.second == 0) {
                dx_u.resetToZero();
            } else {
                PLB_ASSERT(false && "xDerivDirAndOrder has invalid order");
            }
            if (yDerivDirAndOrder.second == 2) {
                computeCentralGradient(
                    rhoBarField, jField, this->intNode + offset, this->intNode + offsetScalar,
                    Dot3D(0, 1, 0), dy_u);
            } else if (yDerivDirAndOrder.second == 1) {
                PLB_ASSERT(
                    (yDerivDirAndOrder.first == 1 || yDerivDirAndOrder.first == -1)
                    && "Invalid yDerivDirAndOrder.first value (must be = 1 or = -1");
                computeFwdGradient(
                    rhoBarField, jField, this->u, this->intNode + offset,
                    this->intNode + offsetScalar, Dot3D(0, yDerivDirAndOrder.first, 0),
                    yDerivDirAndOrder.first, dy_u);
            } else if (yDerivDirAndOrder.second == 0) {
                dy_u.resetToZero();
            } else {
                PLB_ASSERT(false && "yDerivDirAndOrder has invalid order");
            }
            if (zDerivDirAndOrder.second == 2) {
                computeCentralGradient(
                    rhoBarField, jField, this->intNode + offset, this->intNode + offsetScalar,
                    Dot3D(0, 0, 1), dz_u);
            } else if (zDerivDirAndOrder.second == 1) {
                PLB_ASSERT(
                    (zDerivDirAndOrder.first == 1 || zDerivDirAndOrder.first == -1)
                    && "Invalid zDerivDirAndOrder.first value (must be = 1 or = -1");
                computeFwdGradient(
                    rhoBarField, jField, this->u, this->intNode + offset,
                    this->intNode + offsetScalar, Dot3D(0, 0, zDerivDirAndOrder.first),
                    zDerivDirAndOrder.first, dz_u);
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
    T rho = Descriptor<T>::fullRho(rhoBar);
    T sToPi = -rho * Descriptor<T>::cs2 / omega;
    Array<T, Descriptor<T>::d> j = rho * this->u;

    typedef SymmetricTensorImpl<T, 3> S;
    Array<T, SymmetricTensor<T, Descriptor>::n> pi;
    pi[S::xx] = (T)2 * dx_ux * sToPi;
    pi[S::yy] = (T)2 * dy_uy * sToPi;
    pi[S::zz] = (T)2 * dz_uz * sToPi;
    pi[S::xy] = (dx_uy + dy_ux) * sToPi;
    pi[S::xz] = (dx_uz + dz_ux) * sToPi;
    pi[S::yz] = (dy_uz + dz_uy) * sToPi;

    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    this->cell.getDynamics().regularize(this->cell, rhoBar, j, jSqr, pi);
}

// ====================================================================== //
// ============== InterpolatedGeneralizedOffLatticeModel3D ============== //
// ====================================================================== //

template <typename T, template <typename U> class Descriptor>
InterpolatedGeneralizedOffLatticeModel3D<T, Descriptor>::InterpolatedGeneralizedOffLatticeModel3D(
    BoundaryShape3D<T, Array<T, 3> > *shape_, int flowType_) :
    OffLatticeModel3D<T, Array<T, 3> >(shape_, flowType_)
{ }

template <typename T, template <typename U> class Descriptor>
InterpolatedGeneralizedOffLatticeModel3D<T, Descriptor>
    *InterpolatedGeneralizedOffLatticeModel3D<T, Descriptor>::clone() const
{
    return new InterpolatedGeneralizedOffLatticeModel3D(*this);
}

template <typename T, template <typename U> class Descriptor>
plint InterpolatedGeneralizedOffLatticeModel3D<T, Descriptor>::getNumNeighbors() const
{
    return 2;
}

template <typename T, template <typename U> class Descriptor>
bool InterpolatedGeneralizedOffLatticeModel3D<T, Descriptor>::isExtrapolated() const
{
    return false;
}

template <typename T, template <typename U> class Descriptor>
void InterpolatedGeneralizedOffLatticeModel3D<T, Descriptor>::prepareCell(
    Dot3D const &cellLocation, AtomicContainerBlock3D &container)
{
    Dot3D offset = container.getLocation();
    InterpolatedGeneralizedOffLatticeInfo3D *info =
        dynamic_cast<InterpolatedGeneralizedOffLatticeInfo3D *>(container.getData());
    PLB_ASSERT(info);
    if (this->isFluid(cellLocation + offset)) {
        std::vector<std::pair<int, int> > solidNeighborsNoSolidNeighbors;
        std::vector<int> solidNeighbors;
        std::vector<int> wetNodeFluidDirections;
        std::vector<plint> ids;
        for (int iNeighbor = 0; iNeighbor < NextNeighbor<T>::numNeighbors; ++iNeighbor) {
            int const *c = NextNeighbor<T>::c[iNeighbor];
            Dot3D neighbor(cellLocation.x + c[0], cellLocation.y + c[1], cellLocation.z + c[2]);
            // TODO use only depth two if possible (add a check for only depth
            // two selection or depth one if depth two not possible.

            // If the fluid node has a non-fluid neighbor ...
            if (this->isSolid(neighbor + offset)) {
                // used for the drag/ift computation
                solidNeighbors.push_back(iNeighbor);
                // ... check how many fluid nodes without solid neighbors
                // in the direction opposite to the solid.
                int depth = 0;
                for (int iDepth = 1; iDepth <= getNumNeighbors(); ++iDepth) {
                    Dot3D nextNeighbor(
                        cellLocation.x - iDepth * c[0], cellLocation.y - iDepth * c[1],
                        cellLocation.z - iDepth * c[2]);

                    bool solidNeighbor = false;
                    for (int iPop = 1; iPop < Descriptor<T>::q; ++iPop) {
                        Dot3D nextNeighbor_neighbor(
                            nextNeighbor.x + Descriptor<T>::c[iPop][0],
                            nextNeighbor.y + Descriptor<T>::c[iPop][1],
                            nextNeighbor.z + Descriptor<T>::c[iPop][2]);

                        if (this->isSolid(nextNeighbor_neighbor + offset)) {
                            solidNeighbor = true;
                            break;
                        }
                    }

                    if (this->isFluid(nextNeighbor + offset) && !solidNeighbor) {
                        depth = iDepth;
                    } else {
                        break;
                    }
                }
                // ... then add this node to the list.
                solidNeighborsNoSolidNeighbors.push_back(std::make_pair(iNeighbor, depth));
                Array<T, 3> locatedPoint;
                T distance;
                Array<T, 3> wallNormal;
                Array<T, 3> surfaceData;
                plint iTriangle = -1;
                OffBoundary::Type bdType;
#ifdef PLB_DEBUG
                bool ok =
#endif
                    this->pointOnSurface(
                        cellLocation + offset, Dot3D(c[0], c[1], c[2]), locatedPoint, distance,
                        wallNormal, surfaceData, bdType, iTriangle);
                PLB_ASSERT(ok);
                ids.push_back(iTriangle);
            }
        }

        if (!solidNeighborsNoSolidNeighbors.empty()) {
            // selecting only deepest depth available (for a higher order interpolation order)
            std::vector<std::vector<std::pair<int, int> > > selectionSolid(getNumNeighbors() + 1);
            std::vector<std::vector<plint> > selectionIds(getNumNeighbors() + 1);
            for (pluint iA = 0; iA < solidNeighborsNoSolidNeighbors.size(); ++iA) {
                selectionSolid[solidNeighborsNoSolidNeighbors[iA].second].push_back(
                    solidNeighborsNoSolidNeighbors[iA]);
                selectionIds[solidNeighborsNoSolidNeighbors[iA].second].push_back(ids[iA]);
            }
            // adding fluid boundary wet node's fluid directions (know populations directions).
            for (int iPop = 1; iPop < Descriptor<T>::q; ++iPop) {
                int const *c = Descriptor<T>::c[iPop];
                Dot3D potFluidNeighbor(
                    cellLocation.x + c[0], cellLocation.y + c[1], cellLocation.z + c[2]);

                if (this->isFluid(potFluidNeighbor + offset)) {
                    wetNodeFluidDirections.push_back(iPop);
                }
            }
            // pcout << "solidNeighbors.size() = " << solidNeighbors.size() << std::endl;
            // pcout << "ids.size() = " << ids.size() << std::endl;
            info->getSolidNeighbors().push_back(solidNeighbors);
            info->getWetNodes().push_back(cellLocation);
            info->getAllWetNodeIds().push_back(ids);
            // add only biggest depth to the list of directions and triangles
            for (plint iA = getNumNeighbors(); iA >= 0; --iA) {
                if (!selectionSolid[iA].empty()) {
                    info->getWetNodeSolidUsableDirections().push_back(selectionSolid[iA]);
                    info->getWetNodeIds().push_back(selectionIds[iA]);
                    break;
                }
            }
            info->getWetNodeFluidDirections().push_back(wetNodeFluidDirections);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ContainerBlockData *
    InterpolatedGeneralizedOffLatticeModel3D<T, Descriptor>::generateOffLatticeInfo() const
{
    return new InterpolatedGeneralizedOffLatticeInfo3D;
}

template <typename T, template <typename U> class Descriptor>
Array<T, 3> InterpolatedGeneralizedOffLatticeModel3D<T, Descriptor>::getLocalForce(
    AtomicContainerBlock3D &container) const
{
    InterpolatedGeneralizedOffLatticeInfo3D *info =
        dynamic_cast<InterpolatedGeneralizedOffLatticeInfo3D *>(container.getData());
    PLB_ASSERT(info);
    return info->getLocalForce();
}

template <typename T, template <typename U> class Descriptor>
void InterpolatedGeneralizedOffLatticeModel3D<T, Descriptor>::boundaryCompletion(
    AtomicBlock3D &nonTypeLattice, AtomicContainerBlock3D &container,
    std::vector<AtomicBlock3D *> const &args)
{
    BlockLattice3D<T, Descriptor> &lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> &>(nonTypeLattice);
    InterpolatedGeneralizedOffLatticeInfo3D *info =
        dynamic_cast<InterpolatedGeneralizedOffLatticeInfo3D *>(container.getData());
    PLB_ASSERT(info);
    std::vector<Dot3D> const &wetNodes = info->getWetNodes();
    std::vector<std::vector<std::pair<int, int> > > const &wetNodeSolidUsableDirections =
        info->getWetNodeSolidUsableDirections();
    std::vector<std::vector<int> > const &wetNodeFluidDirections =
        info->getWetNodeFluidDirections();
    std::vector<std::vector<plint> > const &wetNodeIds = info->getWetNodeIds();
    std::vector<std::vector<plint> > const &allWetNodeIds = info->getAllWetNodeIds();
    std::vector<std::vector<int> > const &solidDirections = info->getSolidNeighbors();
    PLB_ASSERT(wetNodes.size() == wetNodeSolidUsableDirections.size());

    Dot3D absoluteOffset = container.getLocation();

    Array<T, 3> &localForce = info->getLocalForce();
    localForce.resetToZero();
    for (pluint iWet = 0; iWet < wetNodes.size(); ++iWet) {
        cellCompletion(
            lattice, wetNodes[iWet], wetNodeSolidUsableDirections[iWet],
            wetNodeFluidDirections[iWet], wetNodeIds[iWet], allWetNodeIds[iWet],
            solidDirections[iWet], absoluteOffset, localForce, args);
    }
}

template <typename T, template <typename U> class Descriptor>
void InterpolatedGeneralizedOffLatticeModel3D<T, Descriptor>::cellCompletion(
    BlockLattice3D<T, Descriptor> &lattice, Dot3D const &intNode,
    std::vector<std::pair<int, int> > const &wetNodeSolidUsableDirections,
    std::vector<int> const &wetNodeFluidDirections, std::vector<plint> const &wetNodeIds,
    std::vector<plint> const &allWetNodeIds, std::vector<int> const &solidDirections,
    Dot3D const &absoluteOffset, Array<T, 3> &localForce, std::vector<AtomicBlock3D *> const &args)
{
    InterpolationAlgorithm3D<T, Descriptor> *algorithm =
        new InterpolationGeneralizedAlgorithm3D<T, Descriptor>(
            *this, lattice, intNode, wetNodeSolidUsableDirections, wetNodeFluidDirections,
            wetNodeIds, allWetNodeIds, solidDirections, absoluteOffset, localForce, args,
            this->computesStat(), this->usesSecondOrder());

    bool ok = algorithm->computeNeighborData();
    if (!ok) {
        global::plbErrors().registerError(
            "Error treating the geometry in the interpolated generalized off-lattice model.");
    }
    algorithm->finalize();
    delete algorithm;
}

// ====================================================================== //
// ================== InterpolatedFdOffLatticeModel3D =================== //
// ====================================================================== //

template <typename T, template <typename U> class Descriptor>
InterpolatedFdOffLatticeModel3D<T, Descriptor>::InterpolatedFdOffLatticeModel3D(
    BoundaryShape3D<T, Array<T, 3> > *shape_, int flowType_) :
    OffLatticeModel3D<T, Array<T, 3> >(shape_, flowType_)
{ }

template <typename T, template <typename U> class Descriptor>
InterpolatedFdOffLatticeModel3D<T, Descriptor>
    *InterpolatedFdOffLatticeModel3D<T, Descriptor>::clone() const
{
    return new InterpolatedFdOffLatticeModel3D(*this);
}

template <typename T, template <typename U> class Descriptor>
plint InterpolatedFdOffLatticeModel3D<T, Descriptor>::getNumNeighbors() const
{
    return 2;
}

template <typename T, template <typename U> class Descriptor>
bool InterpolatedFdOffLatticeModel3D<T, Descriptor>::isExtrapolated() const
{
    return false;
}

template <typename T, template <typename U> class Descriptor>
bool InterpolatedFdOffLatticeModel3D<T, Descriptor>::isUsable(const Dot3D &pos) const
{
    if (this->isFluid(pos)) {
        return true;
    }

    return false;
}

template <typename T, template <typename U> class Descriptor>
std::pair<int, int> InterpolatedFdOffLatticeModel3D<T, Descriptor>::computeOrderAndDirection(
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
void InterpolatedFdOffLatticeModel3D<T, Descriptor>::prepareCell(
    Dot3D const &cellLocation, AtomicContainerBlock3D &container)
{
    Dot3D offset = container.getLocation();
    InterpolatedFdOffLatticeInfo3D *info =
        dynamic_cast<InterpolatedFdOffLatticeInfo3D *>(container.getData());
    PLB_ASSERT(info);
    if (this->isFluid(cellLocation + offset)) {
        std::vector<std::pair<int, int> > solidNeighborsNoSolidNeighbors;
        std::vector<int> solidNeighbors;
        std::vector<int> wetNodeFluidDirections;
        std::vector<plint> ids;
        for (int iNeighbor = 0; iNeighbor < NextNeighbor<T>::numNeighbors; ++iNeighbor) {
            int const *c = NextNeighbor<T>::c[iNeighbor];
            Dot3D neighbor(cellLocation.x + c[0], cellLocation.y + c[1], cellLocation.z + c[2]);
            // TODO use only depth two if possible (add a check for only depth
            // two selection or depth one if depth two not possible.

            // If the fluid node has a non-fluid neighbor ...
            if (this->isSolid(neighbor + offset)) {
                solidNeighbors.push_back(iNeighbor);
                // ... check how many fluid nodes without solid neighbors
                // in the direction opposite to the solid.
                int depth = 0;
                for (int iDepth = 1; iDepth <= getNumNeighbors(); ++iDepth) {
                    Dot3D nextNeighbor(
                        cellLocation.x - iDepth * c[0], cellLocation.y - iDepth * c[1],
                        cellLocation.z - iDepth * c[2]);

                    bool solidNeighbor = false;
                    for (int iPop = 1; iPop < Descriptor<T>::q; ++iPop) {
                        Dot3D nextNeighbor_neighbor(
                            nextNeighbor.x + Descriptor<T>::c[iPop][0],
                            nextNeighbor.y + Descriptor<T>::c[iPop][1],
                            nextNeighbor.z + Descriptor<T>::c[iPop][2]);

                        if ((this->isSolid(nextNeighbor_neighbor + offset))) {
                            solidNeighbor = true;
                            break;
                        }
                    }

                    if (this->isFluid(nextNeighbor + offset) && !solidNeighbor) {
                        depth = iDepth;
                    } else {
                        break;
                    }
                }
                // ... then add this node to the list.
                solidNeighborsNoSolidNeighbors.push_back(std::make_pair(iNeighbor, depth));
                Array<T, 3> locatedPoint;
                T distance;
                Array<T, 3> wallNormal;
                Array<T, 3> surfaceData;
                plint iTriangle = -1;
                OffBoundary::Type bdType;
#ifdef PLB_DEBUG
                bool ok =
#endif
                    this->pointOnSurface(
                        cellLocation + offset, Dot3D(c[0], c[1], c[2]), locatedPoint, distance,
                        wallNormal, surfaceData, bdType, iTriangle);
                PLB_ASSERT(ok);
                ids.push_back(iTriangle);
            }
        }

        if (!solidNeighborsNoSolidNeighbors.empty()) {
            // selecting only deepest depth available (for a higher order interpolation order)
            std::vector<std::vector<std::pair<int, int> > > selectionSolid(getNumNeighbors() + 1);
            std::vector<std::vector<plint> > selectionIds(getNumNeighbors() + 1);
            for (pluint iA = 0; iA < solidNeighborsNoSolidNeighbors.size(); ++iA) {
                selectionSolid[solidNeighborsNoSolidNeighbors[iA].second].push_back(
                    solidNeighborsNoSolidNeighbors[iA]);
                selectionIds[solidNeighborsNoSolidNeighbors[iA].second].push_back(ids[iA]);
            }
            // adding fluid boundary wet node's fluid directions (know populations directions).
            for (int iPop = 1; iPop < Descriptor<T>::q; ++iPop) {
                int const *c = Descriptor<T>::c[iPop];
                Dot3D potFluidNeighbor(
                    cellLocation.x + c[0], cellLocation.y + c[1], cellLocation.z + c[2]);

                if (this->isFluid(potFluidNeighbor + offset)) {
                    wetNodeFluidDirections.push_back(iPop);
                }
            }

            Dot3D pos = cellLocation + offset;
            Dot3D dx(1, 0, 0), dy(0, 1, 0), dz(0, 0, 1);
            info->getXderivDirAndOrder().push_back(computeOrderAndDirection(pos, dx));
            info->getYderivDirAndOrder().push_back(computeOrderAndDirection(pos, dy));
            info->getZderivDirAndOrder().push_back(computeOrderAndDirection(pos, dz));

            info->getSolidNeighbors().push_back(solidNeighbors);
            info->getWetNodes().push_back(cellLocation);
            // add only biggest depth to the list of directions and triangles
            for (plint iA = getNumNeighbors(); iA >= 0; --iA) {
                if (!selectionSolid[iA].empty()) {
                    info->getWetNodeSolidUsableDirections().push_back(selectionSolid[iA]);
                    info->getWetNodeIds().push_back(selectionIds[iA]);
                    break;
                }
            }
            info->getWetNodeFluidDirections().push_back(wetNodeFluidDirections);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ContainerBlockData *InterpolatedFdOffLatticeModel3D<T, Descriptor>::generateOffLatticeInfo() const
{
    return new InterpolatedFdOffLatticeInfo3D;
}

template <typename T, template <typename U> class Descriptor>
Array<T, 3> InterpolatedFdOffLatticeModel3D<T, Descriptor>::getLocalForce(
    AtomicContainerBlock3D &container) const
{
    InterpolatedFdOffLatticeInfo3D *info =
        dynamic_cast<InterpolatedFdOffLatticeInfo3D *>(container.getData());
    PLB_ASSERT(info);
    return info->getLocalForce();
}

template <typename T, template <typename U> class Descriptor>
void InterpolatedFdOffLatticeModel3D<T, Descriptor>::boundaryCompletion(
    AtomicBlock3D &nonTypeLattice, AtomicContainerBlock3D &container,
    std::vector<AtomicBlock3D *> const &args)
{
    BlockLattice3D<T, Descriptor> &lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> &>(nonTypeLattice);
    InterpolatedFdOffLatticeInfo3D *info =
        dynamic_cast<InterpolatedFdOffLatticeInfo3D *>(container.getData());
    PLB_ASSERT(info);
    std::vector<Dot3D> const &wetNodes = info->getWetNodes();
    std::vector<std::vector<std::pair<int, int> > > const &wetNodeSolidUsableDirections =
        info->getWetNodeSolidUsableDirections();
    std::vector<std::vector<int> > const &wetNodeFluidDirections =
        info->getWetNodeFluidDirections();
    std::vector<std::vector<int> > const &solidDirections = info->getSolidNeighbors();
    std::vector<std::vector<plint> > const &wetNodeIds = info->getWetNodeIds();
    std::vector<std::pair<int, int> > const &xDerivDirAndOrder = info->getXderivDirAndOrder();
    std::vector<std::pair<int, int> > const &yDerivDirAndOrder = info->getYderivDirAndOrder();
    std::vector<std::pair<int, int> > const &zDerivDirAndOrder = info->getZderivDirAndOrder();
    PLB_ASSERT(wetNodes.size() == wetNodeSolidUsableDirections.size());

    Dot3D absoluteOffset = container.getLocation();
    Array<T, 3> &localForce = info->getLocalForce();
    localForce.resetToZero();
    for (pluint iWet = 0; iWet < wetNodes.size(); ++iWet) {
        cellCompletion(
            lattice, wetNodes[iWet], wetNodeSolidUsableDirections[iWet],
            wetNodeFluidDirections[iWet], wetNodeIds[iWet], solidDirections[iWet], absoluteOffset,
            xDerivDirAndOrder[iWet], yDerivDirAndOrder[iWet], zDerivDirAndOrder[iWet], localForce,
            args);
    }
}

template <typename T, template <typename U> class Descriptor>
void InterpolatedFdOffLatticeModel3D<T, Descriptor>::cellCompletion(
    BlockLattice3D<T, Descriptor> &lattice, Dot3D const &intNode,
    std::vector<std::pair<int, int> > const &wetNodeSolidUsableDirections,
    std::vector<int> const &wetNodeFluidDirections, std::vector<plint> const &wetNodeIds,
    std::vector<int> const &solidDirections, Dot3D const &absoluteOffset,
    const std::pair<int, int> &xDerivDirAndOrder, const std::pair<int, int> &yDerivDirAndOrder,
    const std::pair<int, int> &zDerivDirAndOrder, Array<T, 3> &localForce,
    std::vector<AtomicBlock3D *> const &args)
{
    InterpolationAlgorithm3D<T, Descriptor> *algorithm = 0;
    if (this->getDefineVelocity()) {
        algorithm = new InterpolationDefineVelocityOffLatticeModel3D<T, Descriptor>(
            *this, lattice, intNode, wetNodeSolidUsableDirections, wetNodeFluidDirections,
            wetNodeIds, solidDirections, absoluteOffset, localForce, args, this->computesStat(),
            this->usesSecondOrder());
    } else {
        algorithm = new InterpolationFdCompletionAlgorithm3D<T, Descriptor>(
            *this, lattice, intNode, wetNodeSolidUsableDirections, wetNodeFluidDirections,
            wetNodeIds, solidDirections, absoluteOffset, localForce, args, this->computesStat(),
            this->usesSecondOrder(), xDerivDirAndOrder, yDerivDirAndOrder, zDerivDirAndOrder);
    }
    bool ok = algorithm->computeNeighborData();
    if (!ok) {
        global::plbErrors().registerError(
            "Error treating the geometry in the interpolated Fd off-lattice model.");
    }
    algorithm->finalize();
    delete algorithm;
}

}  // namespace plb

#endif  // GUO_OFF_LATTICE_MODEL_3D_HH
