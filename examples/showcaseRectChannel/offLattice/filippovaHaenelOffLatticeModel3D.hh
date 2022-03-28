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

#ifndef FILIPPOVA_HAENEL_OFF_LATTICE_MODEL_3D_HH
#define FILIPPOVA_HAENEL_OFF_LATTICE_MODEL_3D_HH

#include <algorithm>
#include <cmath>
#include <vector>

#include "core/dynamics.h"
#include "latticeBoltzmann/externalFieldAccess.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "offLattice/filippovaHaenelOffLatticeModel3D.h"
#include "offLattice/nextNeighbors3D.h"

namespace plb {

/**
 * This class implements the Filippova-Haenel (FH,1998) boundary condition on a BoundaryShape.
 * The BoundaryShape determines whether the points of the discrete lattice are "inside"
 * or "outside" some geometry.
 *
 * It can handle moving boundaries using the momentum correction of ladd (LADD, 1994).
 * The wall velocity is recovered from SurfaceData stored in BoundaryShape3D<T,SurfaceData>*
 *
 * IMPORTANT NOTE: in Palabos versions before June 2020 the FilippovaHaenelLocalModel3D boundary
 * condition refers to the MeiLuoShyy (MLS,1999) variant than now has a independent implementation
 * in offLattice/meiLuoShyyOffLatticeModel3D.h. The name of the class has been changed to have a
 * BREAKING CHANGE in order to encourage the user to choose and check the implementation they need
 * in their application.
 *
 * (FH,1998) O. Filippova and D. Hänel, “Grid Refinement for Lattice-BGK Models,”
 *     Journal of Computational Physics, vol. 147, no. 1, pp. 219–228, Nov. 1998,
 * doi: 10.1006/jcph.1998.6089.
 *
 * (MLS,1999) R. Mei, L.-S. Luo, and W. Shyy, “An Accurate Curved Boundary Treatment in the Lattice
 * Boltzmann Method,” Journal of Computational Physics, vol. 155, no. 2, pp. 307–330, Nov. 1999,
 * doi: 10.1006/jcph.1999.6334.
 *
 * (LADD, 1994) A. J. C. Ladd, “Numerical simulations of particulate suspensions via a discretized
 * Boltzmann equation. Part 1. Theoretical foundation,” Journal of Fluid Mechanics, vol. 271, pp.
 * 285–309, Jul. 1994, doi: 10.1017/S0022112094001771.
 *
 * @tparam T
 * @tparam Descriptor
 */
template <typename T, template <typename U> class Descriptor>
FilippovaHaenelLocalModel3D<T, Descriptor>::FilippovaHaenelLocalModel3D(
    BoundaryShape3D<T, Array<T, 3> > *shape_, int flowType_) :
    OffLatticeModel3D<T, Array<T, 3> >(shape_, flowType_)
{ }

template <typename T, template <typename U> class Descriptor>
FilippovaHaenelLocalModel3D<T, Descriptor> *FilippovaHaenelLocalModel3D<T, Descriptor>::clone()
    const
{
    return new FilippovaHaenelLocalModel3D(*this);
}

template <typename T, template <typename U> class Descriptor>
plint FilippovaHaenelLocalModel3D<T, Descriptor>::getNumNeighbors() const
{
    return 1;
}

/**
 *  Filippova-Haenel is a completion scheme for a layer of cells on the
 *  "solid" side of the boundary, like Guo.
 * @tparam T
 * @tparam Descriptor
 * @return true
 */
template <typename T, template <typename U> class Descriptor>
bool FilippovaHaenelLocalModel3D<T, Descriptor>::isExtrapolated() const
{
    return true;
}

template <typename T, template <typename U> class Descriptor>
void FilippovaHaenelLocalModel3D<T, Descriptor>::prepareCell(
    Dot3D const &cellLocation, AtomicContainerBlock3D &container)
{
    typedef Descriptor<T> D;
    Dot3D offset = container.getLocation();
    OffLatticeInfo3D *info = dynamic_cast<OffLatticeInfo3D *>(container.getData());
    PLB_ASSERT(info);
    std::vector<int> liquidNeighbors;
    std::vector<plint> ids;
    Dot3D absLoc = cellLocation + offset;
    if (this->isSolid(absLoc)) {
        for (int iPop = 0; iPop < D::q; ++iPop) {
            Dot3D neighbor(
                cellLocation.x + D::c[iPop][0], cellLocation.y + D::c[iPop][1],
                cellLocation.z + D::c[iPop][2]);
            Dot3D neighborLoc = neighbor + offset;
            // If the non-fluid node has a fluid neighbor ...
            if (this->isFluid(neighborLoc)) {
                // ... check how many fluid nodes it has ahead of it ...
                plint iTriangle = -1;
                global::timer("intersect").start();
                Array<T, 3> locatedPoint;
                T distance;
                Array<T, 3> wallNormal;
                Array<T, 3> surfaceData;
                OffBoundary::Type bdType;
#ifdef PLB_DEBUG
                bool ok =
#endif
                    this->pointOnSurface(
                        cellLocation + offset, Dot3D(D::c[iPop][0], D::c[iPop][1], D::c[iPop][2]),
                        locatedPoint, distance, wallNormal, surfaceData, bdType, iTriangle);

                // In the following, the importance of directions is sorted wrt. how well they
                //   are aligned with the wall normal. It is better to take the continuous normal,
                //   because it is not sensitive to the choice of the triangle when we shoot at
                //   an edge.
                // wallNormal = this->computeContinuousNormal(locatedPoint, iTriangle);
                global::timer("intersect").stop();
                PLB_ASSERT(ok);
                // ... then add this node to the list.
                liquidNeighbors.push_back(iPop);
                ids.push_back(iTriangle);
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
ContainerBlockData *FilippovaHaenelLocalModel3D<T, Descriptor>::generateOffLatticeInfo() const
{
    return new OffLatticeInfo3D;
}

template <typename T, template <typename U> class Descriptor>
Array<T, 3> FilippovaHaenelLocalModel3D<T, Descriptor>::getLocalForce(
    AtomicContainerBlock3D &container) const
{
    OffLatticeInfo3D *info = dynamic_cast<OffLatticeInfo3D *>(container.getData());
    PLB_ASSERT(info);
    return info->getLocalForce();
}

template <typename T, template <typename U> class Descriptor>
void FilippovaHaenelLocalModel3D<T, Descriptor>::boundaryCompletion(
    AtomicBlock3D &nonTypeLattice, AtomicContainerBlock3D &container,
    std::vector<AtomicBlock3D *> const &args)
{
    BlockLattice3D<T, Descriptor> &lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> &>(nonTypeLattice);
    OffLatticeInfo3D *info = dynamic_cast<OffLatticeInfo3D *>(container.getData());
    PLB_ASSERT(info);
    std::vector<Dot3D> const &dryNodes = info->getDryNodes();
    std::vector<std::vector<int> > const &dryNodeFluidDirections =
        info->getDryNodeFluidDirections();
    std::vector<std::vector<plint> > const &dryNodeIds = info->getDryNodeIds();
    PLB_ASSERT(dryNodes.size() == dryNodeFluidDirections.size());

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
void FilippovaHaenelLocalModel3D<T, Descriptor>::cellCompletion(
    BlockLattice3D<T, Descriptor> &lattice, Dot3D const &guoNode,
    std::vector<int> const &dryNodeFluidDirections, std::vector<plint> const &dryNodeIds,
    Dot3D const &absoluteOffset, Array<T, 3> &localForce, std::vector<AtomicBlock3D *> const &args)
{
    typedef Descriptor<T> D;
    Cell<T, Descriptor> &s_cell = lattice.get(guoNode.x, guoNode.y, guoNode.z);
#ifdef PLB_DEBUG
    int noDynId = NoDynamics<T, Descriptor>().getId();
#endif
    PLB_ASSERT(
        s_cell.getDynamics().getId() == noDynId
        && "Filippova-Haenel BC needs the dynamics to be set to NoDynamics.");
    for (plint iDirection = 0; iDirection < (plint)dryNodeFluidDirections.size(); ++iDirection) {
        int iOpp = dryNodeFluidDirections[iDirection];
        int iPop = indexTemplates::opposite<Descriptor<T> >(iOpp);
        Dot3D fluidDirection(D::c[iOpp][0], D::c[iOpp][1], D::c[iOpp][2]);
        plint dryNodeId = dryNodeIds[iDirection];

        Array<T, 3> wallNode, wall_vel;
        T wallDistance;
        OffBoundary::Type bdType;
        Cell<T, Descriptor> const &f_cell = lattice.get(
            guoNode.x + fluidDirection.x, guoNode.y + fluidDirection.y,
            guoNode.z + fluidDirection.z);

        Cell<T, Descriptor> collidedCell(f_cell);
        BlockStatistics statsCopy(lattice.getInternalStatistics());
        collidedCell.collide(statsCopy);

        T f_rhoBar, ff_rhoBar;
        Array<T, 3> f_j, ff_j;
        Array<T, 3> wallNormal;

        if (args.empty()) {
            Cell<T, Descriptor> const &ff_cell = lattice.get(
                guoNode.x + 2 * fluidDirection.x, guoNode.y + 2 * fluidDirection.y,
                guoNode.z + 2 * fluidDirection.z);

            f_cell.getDynamics().computeRhoBarJ(f_cell, f_rhoBar, f_j);
            ff_cell.getDynamics().computeRhoBarJ(ff_cell, ff_rhoBar, ff_j);
        } else {
            if ((plint)args.size() == 1) {
                NTensorField3D<T> const *macroField =
                    dynamic_cast<NTensorField3D<T> const *>(args[0]);
                PLB_ASSERT(macroField);
                // 1 Variable for rhoBar, 3 variables for j.
                PLB_ASSERT(macroField->getNdim() == 4);
                Dot3D off = computeRelativeDisplacement(lattice, *macroField);

                Dot3D macroPos = guoNode + fluidDirection + off;
                Dot3D macroPos2 = macroPos + fluidDirection;
                T const *macroscopic = macroField->get(macroPos.x, macroPos.y, macroPos.z);
                f_rhoBar = macroscopic[0];
                f_j.from_cArray(macroscopic + 1);
                T const *macroscopic2 = macroField->get(macroPos2.x, macroPos2.y, macroPos2.z);
                ff_j.from_cArray(macroscopic2 + 1);
            } else if ((plint)args.size() == 2) {
                // 1 field for rhoBar, 1 field for j.
                // #ifdef PLB_DEBUG
                ScalarField3D<T> const *rhoBarField =
                    dynamic_cast<ScalarField3D<T> const *>(args[0]);
                // #endif
                TensorField3D<T, 3> const *jField =
                    dynamic_cast<TensorField3D<T, 3> const *>(args[1]);
                PLB_ASSERT(rhoBarField);
                PLB_ASSERT(jField);

                Dot3D posJ =
                    guoNode + fluidDirection + computeRelativeDisplacement(lattice, *jField);
                Dot3D posRho =
                    guoNode + fluidDirection + computeRelativeDisplacement(lattice, *rhoBarField);

                f_rhoBar = rhoBarField->get(posRho.x, posRho.y, posRho.z);
                f_j = jField->get(posJ.x, posJ.y, posJ.z);
                Dot3D posJ2 = posJ + fluidDirection;

                ff_j = jField->get(posJ2.x, posJ2.y, posJ2.z);
            } else {
                PLB_ASSERT(false);  // Not implemented for 3 args.
            }
        }

        T f_rho = D::fullRho(f_rhoBar);
        T f_jSqr = normSqr(f_j);

#ifdef PLB_DEBUG
        bool ok =
#endif
            this->pointOnSurface(
                guoNode + absoluteOffset, fluidDirection, wallNode, wallDistance, wallNormal,
                wall_vel, bdType, dryNodeId);
        PLB_ASSERT(ok);

        Array<T, 3> w_j = wall_vel * f_rho;
        T d = std::sqrt(D::cNormSqr[iOpp]);
        PLB_ASSERT(wallDistance <= d);
        T delta = 1.0 - wallDistance / d;

        T kappa = 0.;
        Array<T, 3> wf_j;
        wf_j.resetToZero();
        T omega = f_cell.getDynamics().getOmega();

        if (delta < 0.5) {
            wf_j = f_j;
            kappa = (omega * (2.0 * delta - 1.0)) / (1.0 - omega);

        } else {
            wf_j = f_j * ((delta - 1.0) / delta) + w_j / delta;
            kappa = omega * (2.0 * delta - 1.0);
        }

        T c_i_wf_j_f_j = D::c[iPop][0] * (wf_j[0] - f_j[0]) + D::c[iPop][1] * (wf_j[1] - f_j[1])
                         + D::c[iPop][2] * (wf_j[2] - f_j[2]);

        T c_i_w_j = D::c[iPop][0] * w_j[0] + D::c[iPop][1] * w_j[1] + D::c[iPop][2] * w_j[2];

        T f_ieq = f_cell.getDynamics().computeEquilibrium(iPop, f_rhoBar, f_j, f_jSqr)
                  + D::t[iPop] * D::invCs2 * c_i_wf_j_f_j;
        s_cell[iOpp] = (1.0 - kappa) * collidedCell[iPop] + kappa * f_ieq
                       + 2.0 * D::t[iPop] * D::invCs2 * c_i_w_j;

        localForce[0] += D::c[iPop][0] * (s_cell[iPop] + s_cell[iOpp]);
        localForce[1] += D::c[iPop][1] * (s_cell[iPop] + s_cell[iOpp]);
        localForce[2] += D::c[iPop][2] * (s_cell[iPop] + s_cell[iOpp]);
    }
}

}  // namespace plb

#endif  // FILIPPOVA_HAENEL_OFF_LATTICE_MODEL_3D_HH
