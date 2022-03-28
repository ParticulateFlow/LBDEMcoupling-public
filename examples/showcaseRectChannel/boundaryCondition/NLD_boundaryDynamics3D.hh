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
 * A collection of dynamics classes that implement the algorithms for
 * non-local dynamics (NLD) boundaries -- generic implementation.
 */
#ifndef NLD_BOUNDARY_DYNAMICS_3D_HH
#define NLD_BOUNDARY_DYNAMICS_3D_HH

#include "boundaryCondition/NLD_boundaryDynamics3D.h"
#include "boundaryCondition/boundaryTemplates.h"
#include "core/cell.h"
#include "core/dynamicsIdentifiers.h"
#include "core/plbProfiler.h"
#include "finiteDifference/fdStencils1D.h"
#include "latticeBoltzmann/advectionDiffusionDynamicsTemplates.h"
#include "latticeBoltzmann/advectionDiffusionMomentTemplates.h"
#include "latticeBoltzmann/externalFieldAccess.h"

namespace plb {

/* *************** Class NonLocalBoundaryDynamics3D ****** */

template <typename T, template <typename U> class Descriptor>
NonLocalBoundaryDynamics3D<T, Descriptor>::NonLocalBoundaryDynamics3D(
    Dynamics<T, Descriptor> *baseDynamics_) :
    CompositeDynamics<T, Descriptor>(baseDynamics_, false)
{ }

template <typename T, template <typename U> class Descriptor>
bool NonLocalBoundaryDynamics3D<T, Descriptor>::isBoundary() const
{
    return true;
}

template <typename T, template <typename U> class Descriptor>
bool NonLocalBoundaryDynamics3D<T, Descriptor>::isNonLocal() const
{
    return true;
}

/* *************** Class NLD_VelocityPlaneUtil3D ****** */

template <typename T, template <typename U> class Descriptor>
class NLD_VelocityPlaneUtil3D {
public:
    NLD_VelocityPlaneUtil3D(
        int direction_, int orientation_, plint iX_, plint iY_, plint iZ_,
        BlockLattice3D<T, Descriptor> &lattice_,
        NonLocalBoundaryDynamics3D<T, Descriptor> &dynamics_) :
        direction(direction_),
        orientation(orientation_),
        iX(iX_),
        iY(iY_),
        iZ(iZ_),
        lattice(lattice_),
        dynamics(dynamics_)
    {
        direction1 = (direction + 1) % 3;
        direction2 = (direction + 2) % 3;
    }
    // Given a velocity, compute rhoBar and j locally, using the Zou/he trick.
    void computeMacroscopic(Array<T, 3> const &u, T &rhoBar, Array<T, 3> &j) const
    {
        Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
        T uNormal = u[direction];
        T fNormal = getExternalForceComponent(cell, direction);
        rhoBar =
            flatWall::computeRhoBar<T, Descriptor>(cell, direction, orientation, uNormal, fNormal);
        j = u;
        if (!dynamics.velIsJ()) {
            j *= Descriptor<T>::fullRho(rhoBar);
        }
    }
    // Extrapolate the velocity from nearest neighbors with first-order accuracy.
    // Use all nine neighbors and weight them, to increase the numerical stability
    // and avoid checkerboard patterns.
    Array<T, 3> extrapolateVelocity() const
    {
        Array<T, 3> uTot;
        uTot.resetToZero();
        static T weight0 = 1. / 9.;
        static T weight1 = 1. / 36.;
        static T weight2 = 4. / 9.;
        static T weights[] = {weight1, weight0, weight1, weight0, weight2,
                              weight0, weight1, weight0, weight1};
        plint i = 0;
        plint numCells = 0;
        T sumWeights = T();
        for (plint i1 = -1; i1 <= +1; ++i1) {
            for (plint i2 = -1; i2 <= +1; ++i2) {
                Array<plint, 3> neighb(iX, iY, iZ);
                neighb[direction] -= orientation;
                neighb[direction1] += i1;
                neighb[direction2] += i2;
                Cell<T, Descriptor> &neighborCell = lattice.get(neighb[0], neighb[1], neighb[2]);
                if (!neighborCell.getDynamics().isBoundary()) {
                    Array<T, 3> u;
                    neighborCell.computeVelocity(u);
                    uTot += weights[i] * u;
                    ++numCells;
                    sumWeights += weights[i];
                }
                ++i;
            }
        }
        // Compute appropriate weight if not all cells were usable.
        // Avoid division by zero.
        if (numCells < 9 && numCells > 0) {
            uTot /= sumWeights;
        }
        return uTot;
    }
    // Extrapolate the density from nearest neighbors with first-order accuracy.
    // Use all nine neighbors and weight them, to increase the numerical stability
    // and avoid checkerboard patterns.
    T extrapolateDensity() const
    {
        T rhoTot = T();
        static T weight0 = 1. / 9.;
        static T weight1 = 1. / 36.;
        static T weight2 = 4. / 9.;
        static T weights[] = {weight1, weight0, weight1, weight0, weight2,
                              weight0, weight1, weight0, weight1};
        plint i = 0;
        plint numCells = 0;
        T sumWeights = T();
        for (plint i1 = -1; i1 <= +1; ++i1) {
            for (plint i2 = -1; i2 <= +1; ++i2) {
                Array<plint, 3> neighb(iX, iY, iZ);
                neighb[direction] -= orientation;
                neighb[direction1] += i1;
                neighb[direction2] += i2;
                Cell<T, Descriptor> &neighborCell = lattice.get(neighb[0], neighb[1], neighb[2]);
                if (!neighborCell.getDynamics().isBoundary()) {
                    T rho = neighborCell.computeDensity();
                    rhoTot += weights[i] * rho;
                    ++numCells;
                    sumWeights += weights[i];
                }
                ++i;
            }
        }
        // Compute appropriate weight if not all cells were usable.
        // Avoid division by zero.
        if (numCells < 9 && numCells > 0) {
            rhoTot /= sumWeights;
        }
        return rhoTot;
    }

private:
    int direction, orientation;
    int direction1, direction2;
    int iX, iY, iZ;
    BlockLattice3D<T, Descriptor> &lattice;
    NonLocalBoundaryDynamics3D<T, Descriptor> &dynamics;
};

/* *************** Class NLD_VelocityEdgeUtil3D ****** */

template <typename T, template <typename U> class Descriptor>
class NLD_VelocityEdgeUtil3D {
public:
    NLD_VelocityEdgeUtil3D(
        int plane_, int normal1_, int normal2_, plint iX_, plint iY_, plint iZ_,
        BlockLattice3D<T, Descriptor> &lattice_,
        NonLocalBoundaryDynamics3D<T, Descriptor> &dynamics_) :
        plane(plane_),
        normal1(normal1_),
        normal2(normal2_),
        iX(iX_),
        iY(iY_),
        iZ(iZ_),
        lattice(lattice_),
        dynamics(dynamics_)
    {
        direction1 = (plane + 1) % 3;
        direction2 = (plane + 2) % 3;
    }
    // Computes rhoBar on a neighbor of an edge cell, by calling the appropriate
    // NLD method of the neighbor cell.
    // step1 indicates how many cells to advance along direction1, and step2
    // how many cells to advance along direction2.
    // It is assumed that either step1 or step2 is zero: diagonal neighbors
    // are not computed here.
    T getNeighborRhoBar(plint step1, plint step2)
    {
        PLB_ASSERT(step1 == 0 || step2 == 0);
        Array<T, 3> dummyJ;
        T dummyTheta;
        T rhoBar;
        Array<int, 3> coords(iX, iY, iZ);
        coords[direction1] += -normal1 * step1;
        coords[direction2] += -normal2 * step2;
        // The neighbor of an edge along direction1 or direction2 is necessarily a plane.
        // Therefore, the computation of rhoBar is forwarded to the plane NLD dynamics of this cell.
        if (step1 == 0) {
            NonLocalBoundaryDynamics3D<T, Descriptor>::staticPlaneComputeMacroscopic(
                direction1, normal1, coords[0], coords[1], coords[2], lattice, rhoBar, dummyJ,
                dummyTheta);
            PLB_ASSERT(step2 != 0);
        } else {  // step2==0.
            NonLocalBoundaryDynamics3D<T, Descriptor>::staticPlaneComputeMacroscopic(
                direction2, normal2, coords[0], coords[1], coords[2], lattice, rhoBar, dummyJ,
                dummyTheta);
        }
        return rhoBar;
    }
    // Given u, compute rhoBar and j on the edge. For this, rhoBar is extrapolated
    // from the plane NLD neighbors with second-order accuracy.
    void computeMacroscopic(Array<T, 3> const &u, T &rhoBar, Array<T, 3> &j)
    {
        T rho10 = getNeighborRhoBar(1, 0);
        T rho01 = getNeighborRhoBar(0, 1);
        T rho20 = getNeighborRhoBar(2, 0);
        T rho02 = getNeighborRhoBar(0, 2);
        rhoBar = (T)2 / (T)3 * (rho01 + rho10) - (T)1 / (T)6 * (rho02 + rho20);
        j = u;
        if (dynamics.velIsJ()) {
            j *= Descriptor<T>::fullRho(rhoBar);
        }
    }
    // Computes the velocity-derivative along the edge direction. This works only if
    // the two neighbors on the edge are velocity NLD nodes. Otherwise, zero-values
    // are returned.
    void computeEdgeGradient(Array<T, 3> &gradient)
    {
        Array<T, 3> u_p1, u_m1;
        Array<plint, 3> posPlus(iX, iY, iZ), posMinus(iX, iY, iZ);
        ++posPlus[plane];
        --posMinus[plane];

        NLD_VelocityBoundaryDynamics3D<T, Descriptor> *dynPlus =
            dynamic_cast<NLD_VelocityBoundaryDynamics3D<T, Descriptor> *>(
                &lattice.get(posPlus[0], posPlus[1], posPlus[2]).getDynamics());
        NLD_VelocityBoundaryDynamics3D<T, Descriptor> *dynMinus =
            dynamic_cast<NLD_VelocityBoundaryDynamics3D<T, Descriptor> *>(
                &lattice.get(posMinus[0], posMinus[1], posMinus[2]).getDynamics());
        if (dynPlus && dynMinus) {
            u_p1 = dynPlus->get_u();
            u_m1 = dynMinus->get_u();
            for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
                gradient[iD] = fd::ctl_diff(u_p1[iD], u_m1[iD]);
            }
        } else {
            gradient.resetToZero();
        }
    }
    // Given the velocity on an edge node, compute its derivative into one of the two plane
    // directions (indicated by deriveDirection/deriveOrientation). You also need to indicate
    // the orientation of the plane in the given direction through the values normalDirection/
    // normalOrientation (remember that a NLD dynamics doesn't know its direction/orientation:
    // this information is provided from outside).
    // The derivative is computed with first-order accuracy.
    void computeGradientIntoPlane(
        Array<T, 3> const &u0, int deriveDirection, int deriveOrientation, int normalDirection,
        int normalOrientation, Array<T, 3> &gradient)
    {
        Array<T, 3> u1;
        Array<plint, 3> pos1(iX, iY, iZ);
        pos1[deriveDirection] -= deriveOrientation;

        T rhoBar, dummyTheta;
        NonLocalBoundaryDynamics3D<T, Descriptor>::staticPlaneComputeMacroscopic(
            normalDirection, normalOrientation, pos1[0], pos1[1], pos1[2], lattice, rhoBar, u1,
            dummyTheta);
        if (!dynamics.velIsJ()) {
            u1 /= Descriptor<T>::fullRho(rhoBar);
        }
        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            gradient[iD] = -deriveOrientation * fd::o1_fwd_diff(u0[iD], u1[iD]);
        }
    }
    // Given the velocity on an edge node, execute a completion scheme for the unknown populations.
    // The density is being extrapolated with second-order accuracy from the two plane directions,
    // and the strain-rate computed by finite differences, with second-order accuracy on the edge
    // direction and first-order accuracy on the two plane directions. Finally, all populations are
    // replaced by their regularized value.
    void boundaryCompletion(Array<T, 3> const &u)
    {
        typedef SymmetricTensorImpl<T, Descriptor<T>::d> S;
        T rhoBar;
        Array<T, 3> j;
        computeMacroscopic(u, rhoBar, j);
        T rho = Descriptor<T>::fullRho(rhoBar);
        std::vector<Array<T, 3> > dA_uB_(3);
        computeEdgeGradient(dA_uB_[0]);
        computeGradientIntoPlane(u, direction1, normal1, direction2, normal2, dA_uB_[1]);
        computeGradientIntoPlane(u, direction2, normal2, direction1, normal1, dA_uB_[2]);

        std::vector<Array<T, 3> > dA_uB(3);
        for (int iBeta = 0; iBeta < 3; ++iBeta) {
            dA_uB[plane][iBeta] = dA_uB_[0][iBeta];
            dA_uB[direction1][iBeta] = dA_uB_[1][iBeta];
            dA_uB[direction2][iBeta] = dA_uB_[2][iBeta];
        }

        T sToPi = -rho / Descriptor<T>::invCs2 / dynamics.getOmega();
        Array<T, SymmetricTensor<T, Descriptor>::n> pi;
        pi[S::xx] = (T)2 * dA_uB[0][0] * sToPi;
        pi[S::yy] = (T)2 * dA_uB[1][1] * sToPi;
        pi[S::zz] = (T)2 * dA_uB[2][2] * sToPi;
        pi[S::xy] = (dA_uB[0][1] + dA_uB[1][0]) * sToPi;
        pi[S::xz] = (dA_uB[0][2] + dA_uB[2][0]) * sToPi;
        pi[S::yz] = (dA_uB[1][2] + dA_uB[2][1]) * sToPi;

        T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

        Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            cell[iPop] = dynamics.computeEquilibrium(iPop, rhoBar, j, jSqr)
                         + offEquilibriumTemplates<T, Descriptor>::fromPiToFneq(iPop, pi);
        }
    }

private:
    int plane, normal1, normal2, direction1, direction2;
    plint iX, iY, iZ;
    BlockLattice3D<T, Descriptor> &lattice;
    NonLocalBoundaryDynamics3D<T, Descriptor> &dynamics;
};

/* *************** Class NLD_VelocityCornerUtil3D ****** */

template <typename T, template <typename U> class Descriptor>
class NLD_VelocityCornerUtil3D {
public:
    NLD_VelocityCornerUtil3D(
        int xNormal_, int yNormal_, int zNormal_, plint iX_, plint iY_, plint iZ_,
        BlockLattice3D<T, Descriptor> &lattice_,
        NonLocalBoundaryDynamics3D<T, Descriptor> &dynamics_) :
        xNormal(xNormal_),
        yNormal(yNormal_),
        zNormal(zNormal_),
        iX(iX_),
        iY(iY_),
        iZ(iZ_),
        lattice(lattice_),
        dynamics(dynamics_)
    { }
    // Given the velocity on the current node, compute velocity-derivative
    // in X-direction with first-order accuracy. It is assumed that the neighbor
    // in X-direction is a NLD edge.
    void computeVelocityGradientX(Array<T, 3> const &u0, Array<T, 3> &gradient)
    {
        Array<T, 3> u1;
        T rhoBar, dummyTheta;
        NonLocalBoundaryDynamics3D<T, Descriptor>::staticEdgeComputeMacroscopic(
            0, yNormal, zNormal, iX - 1 * xNormal, iY - 0 * yNormal, iZ - 0 * zNormal, lattice,
            rhoBar, u1, dummyTheta);
        if (!dynamics.velIsJ()) {
            u1 /= Descriptor<T>::fullRho(rhoBar);
        }
        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            gradient[iD] = -xNormal * fd::o1_fwd_diff(u0[iD], u1[iD]);
        }
    }
    // Given the velocity on the current node, compute velocity-derivative
    // in Y-direction with first-order accuracy. It is assumed that the neighbor
    // in Y-direction is a NLD edge.
    void computeVelocityGradientY(Array<T, 3> const &u0, Array<T, 3> &gradient)
    {
        Array<T, 3> u1;
        T rhoBar, dummyTheta;
        NonLocalBoundaryDynamics3D<T, Descriptor>::staticEdgeComputeMacroscopic(
            1, zNormal, xNormal, iX - 0 * xNormal, iY - 1 * yNormal, iZ - 0 * zNormal, lattice,
            rhoBar, u1, dummyTheta);
        if (!dynamics.velIsJ()) {
            u1 /= Descriptor<T>::fullRho(rhoBar);
        }
        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            gradient[iD] = -yNormal * fd::o1_fwd_diff(u0[iD], u1[iD]);
        }
    }
    // Given the velocity on the current node, compute velocity-derivative
    // in Z-direction with first-order accuracy. It is assumed that the neighbor
    // in Z-direction is a NLD edge.
    void computeVelocityGradientZ(Array<T, 3> const &u0, Array<T, 3> &gradient)
    {
        Array<T, 3> u1;
        T rhoBar, dummyTheta;
        NonLocalBoundaryDynamics3D<T, Descriptor>::staticEdgeComputeMacroscopic(
            2, xNormal, yNormal, iX - 0 * xNormal, iY - 0 * yNormal, iZ - 1 * zNormal, lattice,
            rhoBar, u1, dummyTheta);
        if (!dynamics.velIsJ()) {
            u1 /= Descriptor<T>::fullRho(rhoBar);
        }
        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            gradient[iD] = -zNormal * fd::o1_fwd_diff(u0[iD], u1[iD]);
        }
    }
    // Given the velocity on the current node, compute rhoBar and j. The value of
    // rhoBar is extrapolated with first-order accuracy from the three edge directions.
    void computeMacroscopic(Array<T, 3> const &u, T &rhoBar, Array<T, 3> &j)
    {
        T rho100, rho010, rho001;
        T dummyTheta;
        Array<T, 3> dummyJ;
        NonLocalBoundaryDynamics3D<T, Descriptor>::staticEdgeComputeMacroscopic(
            0, yNormal, zNormal, iX - 1 * xNormal, iY - 0 * yNormal, iZ - 0 * zNormal, lattice,
            rho100, dummyJ, dummyTheta);
        NonLocalBoundaryDynamics3D<T, Descriptor>::staticEdgeComputeMacroscopic(
            1, zNormal, xNormal, iX - 0 * xNormal, iY - 1 * yNormal, iZ - 0 * zNormal, lattice,
            rho010, dummyJ, dummyTheta);
        NonLocalBoundaryDynamics3D<T, Descriptor>::staticEdgeComputeMacroscopic(
            2, xNormal, yNormal, iX - 0 * xNormal, iY - 0 * yNormal, iZ - 1 * zNormal, lattice,
            rho001, dummyJ, dummyTheta);
        rhoBar = (T)1 / (T)3 * (rho001 + rho010 + rho100);
        j = u;
        if (dynamics.velIsJ()) {
            j *= Descriptor<T>::fullRho(rhoBar);
        }
    }
    // Given the velocity on a corner node, execute a completion scheme for the unknown populations.
    // The density is being extrapolated with first-order accuracy from the three edge directions,
    // and the strain-rate computed by finite differences, with first-order accuracy along the three
    // edge directions. Finally, all populations are replaced by their regularized value.
    void boundaryCompletion(Array<T, 3> const &u)
    {
        typedef SymmetricTensorImpl<T, Descriptor<T>::d> S;
        T rhoBar;
        Array<T, 3> j;
        computeMacroscopic(u, rhoBar, j);
        T rho = Descriptor<T>::fullRho(rhoBar);
        Array<T, 3> dx_u, dy_u, dz_u;
        computeVelocityGradientX(u, dx_u);
        computeVelocityGradientY(u, dy_u);
        computeVelocityGradientZ(u, dz_u);

        T dx_ux = dx_u[0];
        T dy_ux = dy_u[0];
        T dz_ux = dz_u[0];
        T dx_uy = dx_u[1];
        T dy_uy = dy_u[1];
        T dz_uy = dz_u[1];
        T dx_uz = dx_u[2];
        T dy_uz = dy_u[2];
        T dz_uz = dz_u[2];
        T omega = dynamics.getOmega();
        T sToPi = -rho / Descriptor<T>::invCs2 / omega;
        Array<T, SymmetricTensor<T, Descriptor>::n> pi;
        pi[S::xx] = (T)2 * dx_ux * sToPi;
        pi[S::yy] = (T)2 * dy_uy * sToPi;
        pi[S::zz] = (T)2 * dz_uz * sToPi;
        pi[S::xy] = (dx_uy + dy_ux) * sToPi;
        pi[S::xz] = (dx_uz + dz_ux) * sToPi;
        pi[S::yz] = (dy_uz + dz_uy) * sToPi;

        T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

        Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            cell[iPop] = dynamics.computeEquilibrium(iPop, rhoBar, j, jSqr)
                         + offEquilibriumTemplates<T, Descriptor>::fromPiToFneq(iPop, pi);
        }
    }

private:
    int xNormal, yNormal, zNormal;
    plint iX, iY, iZ;
    BlockLattice3D<T, Descriptor> &lattice;
    NonLocalBoundaryDynamics3D<T, Descriptor> &dynamics;
};

/* *************** Class NLD_VelocityBoundaryDynamics3D ****** */

template <typename T, template <typename U> class Descriptor>
int NLD_VelocityBoundaryDynamics3D<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, NLD_VelocityBoundaryDynamics3D<T, Descriptor> >(
        std::string("NLD_Velocity"));

template <typename T, template <typename U> class Descriptor>
NLD_VelocityBoundaryDynamics3D<T, Descriptor>::NLD_VelocityBoundaryDynamics3D(
    Dynamics<T, Descriptor> *baseDynamics_) :
    NonLocalBoundaryDynamics3D<T, Descriptor>(baseDynamics_)
{ }

template <typename T, template <typename U> class Descriptor>
NLD_VelocityBoundaryDynamics3D<T, Descriptor>::NLD_VelocityBoundaryDynamics3D(
    HierarchicUnserializer &unserializer) :
    NonLocalBoundaryDynamics3D<T, Descriptor>(0)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
void NLD_VelocityBoundaryDynamics3D<T, Descriptor>::defineVelocity(
    Cell<T, Descriptor> &cell, Array<T, Descriptor<T>::d> const &velocity_)
{
    u = velocity_;
}

template <typename T, template <typename U> class Descriptor>
void NLD_VelocityBoundaryDynamics3D<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &velocity_) const
{
    velocity_ = u;
}

template <typename T, template <typename U> class Descriptor>
void NLD_VelocityBoundaryDynamics3D<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    NonLocalBoundaryDynamics3D<T, Descriptor>::serialize(serializer);
    serializer.addValues<T, 3>(u);
}

template <typename T, template <typename U> class Descriptor>
void NLD_VelocityBoundaryDynamics3D<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    NonLocalBoundaryDynamics3D<T, Descriptor>::unserialize(unserializer);
    unserializer.readValues<T, 3>(u);
}

template <typename T, template <typename U> class Descriptor>
int NLD_VelocityBoundaryDynamics3D<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void NLD_VelocityBoundaryDynamics3D<T, Descriptor>::prepareCollision(Cell<T, Descriptor> &cell)
{ }

template <typename T, template <typename U> class Descriptor>
NLD_VelocityBoundaryDynamics3D<T, Descriptor>
    *NLD_VelocityBoundaryDynamics3D<T, Descriptor>::clone() const
{
    return new NLD_VelocityBoundaryDynamics3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void NLD_VelocityBoundaryDynamics3D<T, Descriptor>::planeComputeMacroscopic(
    int direction, int orientation, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar)
{
    NLD_VelocityPlaneUtil3D<T, Descriptor>(direction, orientation, iX, iY, iZ, lattice, *this)
        .computeMacroscopic(this->u, rhoBar, j);
    thetaBar = T();
}

template <typename T, template <typename U> class Descriptor>
void NLD_VelocityBoundaryDynamics3D<T, Descriptor>::edgeComputeMacroscopic(
    int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar)
{
    NLD_VelocityEdgeUtil3D<T, Descriptor> edgeUtil(
        plane, normal1, normal2, iX, iY, iZ, lattice, *this);
    edgeUtil.computeMacroscopic(this->u, rhoBar, j);
    thetaBar = T();
}

template <typename T, template <typename U> class Descriptor>
void NLD_VelocityBoundaryDynamics3D<T, Descriptor>::cornerComputeMacroscopic(
    int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar)
{
    NLD_VelocityCornerUtil3D<T, Descriptor> cornerUtil(
        xNormal, yNormal, zNormal, iX, iY, iZ, lattice, *this);
    cornerUtil.computeMacroscopic(this->u, rhoBar, j);
    thetaBar = T();
}

template <typename T, template <typename U> class Descriptor>
void NLD_VelocityBoundaryDynamics3D<T, Descriptor>::planeBoundaryCompletion(
    int direction, int orientation, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice)
{
    Array<T, 3> j;
    T rhoBar;
    rhoBar = T();
    j.resetToZero();
    NLD_VelocityPlaneUtil3D<T, Descriptor>(direction, orientation, iX, iY, iZ, lattice, *this)
        .computeMacroscopic(this->u, rhoBar, j);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    flatWall::compute_PiNeq<T, Descriptor>(
        this->getBaseDynamics(), cell, direction, orientation, rhoBar, j, jSqr, PiNeq);

    this->getBaseDynamics().regularize(cell, rhoBar, j, jSqr, PiNeq);
}

template <typename T, template <typename U> class Descriptor>
void NLD_VelocityBoundaryDynamics3D<T, Descriptor>::edgeBoundaryCompletion(
    int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice)
{
    //  x
    //  o - -
    //  |
    //  |   A few comments on periodicity. On the scheme on the left, assume that the x-axis is
    //  |   direction1, and the y-axis is direction2. Further, assume that direction2 is periodic
    //      and direction1 is not periodic. Then, the cell "x" must be in the envelope, and it
    //      necessarily has NLD dynamics (because envelopes are automatically updated with dynamic
    //      content when the domain is made periodic). Also, remember that all NLD dynamics objects
    //      are the same, no matter if they are plane, corner or edge, so there is no possible
    //      confusion here. To obtain periodicity, the user will set the "-" cells to background-
    //      dynamics, and not to NLD dynamics. That's how, in the code below, periodicity is
    //      detected.
    int direction1 = (plane + 1) % 3;
    int direction2 = (plane + 2) % 3;
    Array<plint, 3> nbCell1(iX, iY, iZ), nbCell2(iX, iY, iZ);
    nbCell1[direction1] -= normal1;
    nbCell2[direction2] -= normal2;

    // 1. Detect case where a direction is periodic.
    bool isWall1 = lattice.get(nbCell1[0], nbCell1[1], nbCell1[2]).getDynamics().isBoundary();
    bool isWall2 = lattice.get(nbCell2[0], nbCell2[1], nbCell2[2]).getDynamics().isBoundary();

    // If both directions are periodic, there is no need for a boundary condition.
    if (!isWall1 && !isWall2)
        return;
    // If exactly one direction is periodic, then the edge behaves like a plane.
    if (!isWall1) {
        planeBoundaryCompletion(direction1, normal1, iX, iY, iZ, lattice);
        return;
    } else if (!isWall2) {
        planeBoundaryCompletion(direction2, normal2, iX, iY, iZ, lattice);
        return;
    }

    // 2. Implement general edge algorithm.
    NLD_VelocityEdgeUtil3D<T, Descriptor> edgeUtil(
        plane, normal1, normal2, iX, iY, iZ, lattice, *this);
    edgeUtil.boundaryCompletion(this->u);
}

template <typename T, template <typename U> class Descriptor>
void NLD_VelocityBoundaryDynamics3D<T, Descriptor>::cornerBoundaryCompletion(
    int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice)
{
    // Important: to understand periodicity, look at the comment in the edgeBoundaryCompletion
    // function right above.

    Array<plint, 3> middleNbCell(iX, iY, iZ), nbCellX(iX, iY, iZ), nbCellY(iX, iY, iZ),
        nbCellZ(iX, iY, iZ);
    middleNbCell -= Array<plint, 3>(xNormal, yNormal, zNormal);
    nbCellX[0] -= xNormal;
    nbCellY[1] -= yNormal;
    nbCellZ[2] -= zNormal;

    // 1. To start width, detect periodic directions where the corner actually
    //    becomes edge, plane, or non-boundary.

    bool xPeriodic = !lattice.get(iX, iY - yNormal, iZ - zNormal).getDynamics().isBoundary();
    bool yPeriodic = !lattice.get(iX - xNormal, iY, iZ - zNormal).getDynamics().isBoundary();
    bool zPeriodic = !lattice.get(iX - xNormal, iY - yNormal, iZ).getDynamics().isBoundary();

    // Note that we only need to consider cases where corners become edges,
    //   because further degenerate cases (two or three periodic directions)
    //   are treated in the edge algorithm.
    if (xPeriodic) {
        edgeBoundaryCompletion(0, yNormal, zNormal, iX, iY, iZ, lattice);
        return;
    }
    if (yPeriodic) {
        edgeBoundaryCompletion(1, zNormal, xNormal, iX, iY, iZ, lattice);
        return;
    }
    if (zPeriodic) {
        edgeBoundaryCompletion(2, xNormal, yNormal, iX, iY, iZ, lattice);
        return;
    }

    // 2. Implement generic case.
    NLD_VelocityCornerUtil3D<T, Descriptor> cornerUtil(
        xNormal, yNormal, zNormal, iX, iY, iZ, lattice, *this);
    cornerUtil.boundaryCompletion(this->u);
}

/* *************** Class NLD_VelocityNeumannBoundaryDynamics3D ****** */

template <typename T, template <typename U> class Descriptor>
int NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>::id = meta::registerGeneralDynamics<
    T, Descriptor, NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor> >(
    std::string("NLD_VelocityNeumann"));

template <typename T, template <typename U> class Descriptor>
NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>::NLD_VelocityNeumannBoundaryDynamics3D(
    Dynamics<T, Descriptor> *baseDynamics_, bool noPenetration_) :
    NonLocalBoundaryDynamics3D<T, Descriptor>(baseDynamics_), noPenetration(noPenetration_)
{ }

template <typename T, template <typename U> class Descriptor>
NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>::NLD_VelocityNeumannBoundaryDynamics3D(
    HierarchicUnserializer &unserializer) :
    NonLocalBoundaryDynamics3D<T, Descriptor>(0)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
int NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>::prepareCollision(
    Cell<T, Descriptor> &cell)
{ }

template <typename T, template <typename U> class Descriptor>
void NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    NonLocalBoundaryDynamics3D<T, Descriptor>::serialize(serializer);
    serializer.addValue<bool>(noPenetration);
}

template <typename T, template <typename U> class Descriptor>
void NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    NonLocalBoundaryDynamics3D<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue<bool>(noPenetration);
}

template <typename T, template <typename U> class Descriptor>
NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>
    *NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>::clone() const
{
    return new NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>::planeComputeMacroscopic(
    int direction, int orientation, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar)
{
    NLD_VelocityPlaneUtil3D<T, Descriptor> planeUtil(
        direction, orientation, iX, iY, iZ, lattice, *this);
    Array<T, 3> u(planeUtil.extrapolateVelocity());
    planeUtil.computeMacroscopic(u, rhoBar, j);
    if (noPenetration) {
        u[direction] = T();
    }
    thetaBar = T();
}

template <typename T, template <typename U> class Descriptor>
void NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>::edgeComputeMacroscopic(
    int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar)
{
    // General philosophy: The cell being treated here implements a Neumann condition on an
    // edge. Now, if one of the two neighbors in the plane directions is Dirichlet, it is
    // better to treat the current cell as Dirichlet as well, because this is numerically
    // more stable, and this is probably what the user wanted to happen anyway.
    //
    // If none of the neighbors is Dirichlet, then it could be that this Neumann edge node
    // is part of an outlet, and that it is in contact with, say, a free-slip wall. In this
    // case, it is important that the extrapolation is done along the normal of the outlet
    // plane, i.e. in direction of this free-slip wall node.
    //
    // Finally, if all neighbors are Neumann, there exists no preferred direction. In this
    // case, extrapolation is performed from the bulk, along a diagonal direction.

    // This code is non-thermal. Default-assign the temperature.
    thetaBar = T();

    int direction1 = (plane + 1) % 3;
    int direction2 = (plane + 2) % 3;

    Array<plint, 3> nbCell1(iX, iY, iZ), nbCell2(iX, iY, iZ);
    nbCell1[direction1] -= normal1;
    nbCell2[direction2] -= normal2;
    Cell<T, Descriptor> &neighbor1 = lattice.get(nbCell1[0], nbCell1[1], nbCell1[2]);
    Cell<T, Descriptor> &neighbor2 = lattice.get(nbCell2[0], nbCell2[1], nbCell2[2]);

    // Check if the plane neighbors are Dirichlet.
    NLD_VelocityBoundaryDynamics3D<T, Descriptor> *dyn1 =
        dynamic_cast<NLD_VelocityBoundaryDynamics3D<T, Descriptor> *>(&neighbor1.getDynamics());
    NLD_VelocityBoundaryDynamics3D<T, Descriptor> *dyn2 =
        dynamic_cast<NLD_VelocityBoundaryDynamics3D<T, Descriptor> *>(&neighbor2.getDynamics());
    // General case: none of them is Dirichlet. Implement a general Neumann condition
    // by extrapolating along the normal of the outlet plane, or, if no preferred
    // direction exists, from the bulk node in diagonal direction.
    if (!(dyn1 || dyn2)) {
        NLD_VelocityEdgeUtil3D<T, Descriptor> edgeUtil(
            plane, normal1, normal2, iX, iY, iZ, lattice, *this);

        // Now, let's check if one of the plane neighbors is Neumann, while the other
        // is not. In this case it will be concluded that we're in presence of an
        // outlet plane, and that all extrapolations are normal to this plane.
        bool same1 = isSameTypeAsI(nbCell1[0], nbCell1[1], nbCell1[2], lattice);
        bool same2 = isSameTypeAsI(nbCell2[0], nbCell2[1], nbCell2[2], lattice);
        // Neighbor1 is Neumann, while neighbor2 is something else --> the outlet normal goes like
        // normal2.
        if (same1 && !same2) {
            T rhoBar_, thetaBar_;
            Array<T, 3> j_;
            // Extrapolate the velocity from neighbor2, which is on a plane with normal normal1.
            this->staticPlaneComputeMacroscopic(
                direction1, normal1, nbCell2[0], nbCell2[1], nbCell2[2], lattice, rhoBar_, j_,
                thetaBar_);
            if (!lattice.get(iX, iY, iZ).getDynamics().velIsJ()) {
                j_ *= Descriptor<T>::invRho(rhoBar_);
            }
            edgeUtil.computeMacroscopic(j_, rhoBar, j);
            if (noPenetration) {
                j[direction2] = T();
            }
            return;
        }
        // Neighbor2 is Neumann, while neighbor1 is something else --> the outlet normal goes like
        // normal1.
        if (same2 && !same1) {
            T rhoBar_, thetaBar_;
            Array<T, 3> j_;
            // Extrapolate the velocity from neighbor1, which is on a plane with normal normal2.
            this->staticPlaneComputeMacroscopic(
                direction2, normal2, nbCell1[0], nbCell1[1], nbCell1[2], lattice, rhoBar_, j_,
                thetaBar_);
            if (!lattice.get(iX, iY, iZ).getDynamics().velIsJ()) {
                j_ *= Descriptor<T>::invRho(rhoBar_);
            }
            edgeUtil.computeMacroscopic(j_, rhoBar, j);
            if (noPenetration) {
                j[direction1] = T();
            }
            return;
        }
        // Either both or none of them is Neumann. In this case, there is no preferred direction,
        // and extrapolation is taken from the diagonal bulk cell.
        Array<plint, 3> middleNbCell(iX, iY, iZ);
        middleNbCell[direction1] -= normal1;
        middleNbCell[direction2] -= normal2;
        Cell<T, Descriptor> &middleNeighbor =
            lattice.get(middleNbCell[0], middleNbCell[1], middleNbCell[2]);
        Array<T, 3> u;
        middleNeighbor.computeVelocity(u);
        edgeUtil.computeMacroscopic(u, rhoBar, j);
        if (noPenetration) {
            j[direction1] = T();
            j[direction2] = T();
        }
    }
    // At least one of the two plane neighbors is Dirichlet. Treat the current cell
    // as Dirichlet as well, inheriting the velocity value from the neighbor(s).
    else
    {
        T rhoBar1, rhoBar2;
        Array<T, 3> j1, j2;
        // Inherit Dirichlet condition from neighbor1, if neighbor1 is Dirichlet.
        if (dyn1) {
            dyn1->edgeComputeMacroscopic(
                plane, normal1, normal2, iX, iY, iZ, lattice, rhoBar1, j1, thetaBar);
            rhoBar = rhoBar1;
            j = j1;
            // Note: no-penetration condition is not applied when Dirichlet is inherited.
        }
        // Inherit Dirichlet condition from neighbor2, if neighbor2 is Dirichlet.
        if (dyn2) {
            dyn2->edgeComputeMacroscopic(
                plane, normal1, normal2, iX, iY, iZ, lattice, rhoBar2, j2, thetaBar);
            rhoBar = rhoBar2;
            j = j2;
            // Note: no-penetration condition is not applied when Dirichlet is inherited.
        }
        // If both neighbors are Dirichlet, take an average.
        if (dyn1 && dyn2) {
            j = (j1 + j2) / (T)2;
            // Note: no-penetration condition is not applied when Dirichlet is inherited.
            rhoBar = (rhoBar1 + rhoBar2) / (T)2;
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>::cornerComputeMacroscopic(
    int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar)
{
    // General philosophy: The cell being treated here implements a Neumann condition on an
    // corner. Now, if one of the three neighbors in the edge directions is Dirichlet, it is
    // better to treat the current cell as Dirichlet as well, because this is numerically
    // more stable, and this is probably what the user wanted to happen anyway.
    //
    // If none of the neighbors is Dirichlet, then it could be that this Neumann corner node
    // is part of an outlet, and that it is in contact with, say, a free-slip wall. In this
    // case, it is important that the extrapolation is done along the normal of the outlet
    // plane, i.e. in direction of this free-slip wall node.
    //
    // Finally, if either all neighbors or only one neighbor are Neumann, there exists no
    // preferred direction. In this case, extrapolation is performed from the bulk, along
    // a diagonal direction.

    // This code is non-thermal. Default-assign the temperature.
    thetaBar = T();

    Array<plint, 3> nbCell1(iX - xNormal, iY, iZ), nbCell2(iX, iY - yNormal, iZ),
        nbCell3(iX, iY, iZ - zNormal);
    Cell<T, Descriptor> &neighbor1 = lattice.get(nbCell1[0], nbCell1[1], nbCell1[2]);
    Cell<T, Descriptor> &neighbor2 = lattice.get(nbCell2[0], nbCell2[1], nbCell2[2]);
    Cell<T, Descriptor> &neighbor3 = lattice.get(nbCell3[0], nbCell3[1], nbCell3[2]);

    NLD_VelocityBoundaryDynamics3D<T, Descriptor> *dyn1 =
        dynamic_cast<NLD_VelocityBoundaryDynamics3D<T, Descriptor> *>(&neighbor1.getDynamics());
    NLD_VelocityBoundaryDynamics3D<T, Descriptor> *dyn2 =
        dynamic_cast<NLD_VelocityBoundaryDynamics3D<T, Descriptor> *>(&neighbor2.getDynamics());
    NLD_VelocityBoundaryDynamics3D<T, Descriptor> *dyn3 =
        dynamic_cast<NLD_VelocityBoundaryDynamics3D<T, Descriptor> *>(&neighbor3.getDynamics());

    // Worst-case: none of them is Dirichlet. Implement a general Neumann condition
    // by extrapolating from the bulk node in diagonal direction.
    if (!(dyn1 || dyn2 || dyn3)) {
        NLD_VelocityCornerUtil3D<T, Descriptor> cornerUtil(
            xNormal, yNormal, zNormal, iX, iY, iZ, lattice, *this);

        // Now, let's check if one of the edge neighbors is Neumann, while the other
        // two are not. In this case it will be concluded that we're in presence of an
        // outlet plane, and that all extrapolations are normal to this plane.
        bool same1 = isSameTypeAsI(nbCell1[0], nbCell1[1], nbCell1[2], lattice);
        bool same2 = isSameTypeAsI(nbCell2[0], nbCell2[1], nbCell2[2], lattice);
        bool same3 = isSameTypeAsI(nbCell3[0], nbCell3[1], nbCell3[2], lattice);
        // Neighbor 1 and 2 are Neumann, while neighbor 3 is something else --> the outlet normal
        // goes like normal3.
        if (same1 && same2 && !same3) {
            T rhoBar_, thetaBar_;
            Array<T, 3> j_;
            // Extrapolate the velocity from neighbor3, which is on an edge in z-direction.
            this->staticEdgeComputeMacroscopic(
                2, xNormal, yNormal, nbCell3[0], nbCell3[1], nbCell3[2], lattice, rhoBar_, j_,
                thetaBar_);
            if (!lattice.get(iX, iY, iZ).getDynamics().velIsJ()) {
                j_ *= Descriptor<T>::invRho(rhoBar_);
            }
            cornerUtil.computeMacroscopic(j_, rhoBar, j);
            if (noPenetration) {
                j[2] = T();
            }
            return;
        }
        // Neighbor 1 and 3 are Neumann, while neighbor 2 is something else --> the outlet normal
        // goes like normal2.
        if (same1 && same3 && !same2) {
            T rhoBar_, thetaBar_;
            Array<T, 3> j_;
            // Extrapolate the velocity from neighbor2, which is on an edge in y-direction.
            this->staticEdgeComputeMacroscopic(
                1, zNormal, xNormal, nbCell2[0], nbCell2[1], nbCell2[2], lattice, rhoBar_, j_,
                thetaBar_);
            if (!lattice.get(iX, iY, iZ).getDynamics().velIsJ()) {
                j_ *= Descriptor<T>::invRho(rhoBar_);
            }
            cornerUtil.computeMacroscopic(j_, rhoBar, j);
            if (noPenetration) {
                j[1] = T();
            }
            return;
        }
        // Neighbor 2 and 3 are Neumann, while neighbor 1 is something else --> the outlet normal
        // goes like normal1.
        if (same2 && same3 && !same1) {
            T rhoBar_, thetaBar_;
            Array<T, 3> j_;
            // Extrapolate the velocity from neighbor1, which is on an edge in x-direction.
            this->staticEdgeComputeMacroscopic(
                0, yNormal, zNormal, nbCell1[0], nbCell1[1], nbCell1[2], lattice, rhoBar_, j_,
                thetaBar_);
            if (!lattice.get(iX, iY, iZ).getDynamics().velIsJ()) {
                j_ *= Descriptor<T>::invRho(rhoBar_);
            }
            cornerUtil.computeMacroscopic(j_, rhoBar, j);
            if (noPenetration) {
                j[0] = T();
            }
            return;
        }

        // The number of Neumann neighbors is 0 or 3. In this case, there is no preferred direction,
        // and extrapolation is taken from the diagonal bulk cell.
        Array<plint, 3> nbCell(iX, iY, iZ);
        nbCell -= Array<plint, 3>(xNormal, yNormal, zNormal);
        Cell<T, Descriptor> &middleNeighbor = lattice.get(nbCell[0], nbCell[1], nbCell[2]);
        Array<T, 3> u;
        middleNeighbor.computeVelocity(u);

        cornerUtil.computeMacroscopic(u, rhoBar, j);
        if (noPenetration) {
            j.resetToZero();
        }
    }
    // At least one of the three edge neighbors is Dirichlet. Treat the current cell
    // as Dirichlet as well, inheriting the velocity value from the neighbor(s).
    else
    {
        rhoBar = T();
        j.resetToZero();
        plint numNeighbors = 0;
        // Inherit Dirichlet condition from neighbor1, if neighbor1 is Dirichlet.
        if (dyn1) {
            T rhoBar_;
            Array<T, 3> j_;
            dyn1->cornerComputeMacroscopic(
                xNormal, yNormal, zNormal, iX, iY, iZ, lattice, rhoBar_, j_, thetaBar);
            rhoBar += rhoBar_;
            // Note: no-penetration condition is not applied when Dirichlet is inherited.
            j += j_;
            ++numNeighbors;
        }
        // Inherit Dirichlet condition from neighbor2, if neighbor2 is Dirichlet.
        if (dyn2) {
            T rhoBar_;
            Array<T, 3> j_;
            dyn2->cornerComputeMacroscopic(
                xNormal, yNormal, zNormal, iX, iY, iZ, lattice, rhoBar_, j_, thetaBar);
            rhoBar += rhoBar_;
            // Note: no-penetration condition is not applied when Dirichlet is inherited.
            j += j_;
            ++numNeighbors;
        }
        // Inherit Dirichlet condition from neighbor3, if neighbor3 is Dirichlet.
        if (dyn3) {
            T rhoBar_;
            Array<T, 3> j_;
            dyn3->cornerComputeMacroscopic(
                xNormal, yNormal, zNormal, iX, iY, iZ, lattice, rhoBar_, j_, thetaBar);
            rhoBar += rhoBar_;
            // Note: no-penetration condition is not applied when Dirichlet is inherited.
            j += j_;
            ++numNeighbors;
        }
        // If more than one neighbor is Dirichlet, take an average.
        if (numNeighbors > 1) {
            rhoBar /= (T)numNeighbors;
            j /= (T)numNeighbors;
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>::planeBoundaryCompletion(
    int direction, int orientation, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice)
{
    NLD_VelocityPlaneUtil3D<T, Descriptor> planeUtil(
        direction, orientation, iX, iY, iZ, lattice, *this);
    Array<T, 3> u(planeUtil.extrapolateVelocity());
    if (noPenetration) {
        u[direction] = T();
    }
    T rhoBar;
    Array<T, 3> j;
    planeUtil.computeMacroscopic(u, rhoBar, j);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    flatWall::compute_PiNeq<T, Descriptor>(
        this->getBaseDynamics(), cell, direction, orientation, rhoBar, j, jSqr, PiNeq);

    this->getBaseDynamics().regularize(cell, rhoBar, j, jSqr, PiNeq);
}

template <typename T, template <typename U> class Descriptor>
void NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>::edgeBoundaryCompletion(
    int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice)
{
    // Important: To understand periodicity, have a look at the comments in the
    // corresponding method in NLD_VelocityBoundaryDynamics3D.

    int direction1 = (plane + 1) % 3;
    int direction2 = (plane + 2) % 3;

    Array<plint, 3> nbCell1(iX, iY, iZ), nbCell2(iX, iY, iZ), middleNbCell(iX, iY, iZ);
    nbCell1[direction1] -= normal1;
    nbCell2[direction2] -= normal2;
    middleNbCell[direction1] -= normal1;
    middleNbCell[direction2] -= normal2;

    // 1. Detect case where a direction is periodic.
    bool isWall1 = lattice.get(nbCell1[0], nbCell1[1], nbCell1[2]).getDynamics().isBoundary();
    bool isWall2 = lattice.get(nbCell2[0], nbCell2[1], nbCell2[2]).getDynamics().isBoundary();

    // If both directions are periodic, there is no need for a boundary condition.
    if (!isWall1 && !isWall2)
        return;
    // If exactly one direction is periodic, then the edge behaves like a plane.
    if (!isWall1) {
        planeBoundaryCompletion(direction1, normal1, iX, iY, iZ, lattice);
        return;
    } else if (!isWall2) {
        planeBoundaryCompletion(direction2, normal2, iX, iY, iZ, lattice);
        return;
    }

    // General case: no direction is periodic.

    Cell<T, Descriptor> &neighbor1 = lattice.get(nbCell1[0], nbCell1[1], nbCell1[2]);
    Cell<T, Descriptor> &neighbor2 = lattice.get(nbCell2[0], nbCell2[1], nbCell2[2]);
    // Check is one of the two neighbors in plane directions is Dirichlet.
    NLD_VelocityBoundaryDynamics3D<T, Descriptor> *dyn1 =
        dynamic_cast<NLD_VelocityBoundaryDynamics3D<T, Descriptor> *>(&neighbor1.getDynamics());
    NLD_VelocityBoundaryDynamics3D<T, Descriptor> *dyn2 =
        dynamic_cast<NLD_VelocityBoundaryDynamics3D<T, Descriptor> *>(&neighbor2.getDynamics());

    // General case: none of them is Dirichlet. Implement a general Neumann condition
    // by extrapolating along the normal of the outlet plane, or, if no preferred
    // direction exists, from the bulk node in diagonal direction.
    if (!(dyn1 || dyn2)) {
        NLD_VelocityEdgeUtil3D<T, Descriptor> edgeUtil(
            plane, normal1, normal2, iX, iY, iZ, lattice, *this);

        // Now, let's check if one of the plane neighbors is Neumann, while the other
        // is not. In this case it will be concluded that we're in presence of an
        // outlet plane, and that all extrapolations are normal to this plane.
        bool same1 = isSameTypeAsI(nbCell1[0], nbCell1[1], nbCell1[2], lattice);
        bool same2 = isSameTypeAsI(nbCell2[0], nbCell2[1], nbCell2[2], lattice);
        // Neighbor1 is Neumann, while neighbor2 is something else --> the outlet normal goes like
        // normal2.
        if (same1 && !same2) {
            T rhoBar_, thetaBar_;
            Array<T, 3> j_;
            // Extrapolate the velocity from neighbor2, which is on a plane with normal normal1.
            NonLocalBoundaryDynamics3D<T, Descriptor>::staticPlaneComputeMacroscopic(
                direction1, normal1, nbCell2[0], nbCell2[1], nbCell2[2], lattice, rhoBar_, j_,
                thetaBar_);
            if (!lattice.get(iX, iY, iZ).getDynamics().velIsJ()) {
                j_ *= Descriptor<T>::invRho(rhoBar_);
            }
            if (noPenetration) {
                j_[direction2] = T();
            }
            edgeUtil.boundaryCompletion(j_);
            return;
        }
        // Neighbor2 is Neumann, while neighbor1 is something else --> the outlet normal goes like
        // normal1.
        if (same2 && !same1) {
            T rhoBar_, thetaBar_;
            Array<T, 3> j_;
            // Extrapolate the velocity from neighbor1, which is on a plane with normal normal2.
            NonLocalBoundaryDynamics3D<T, Descriptor>::staticPlaneComputeMacroscopic(
                direction2, normal2, nbCell1[0], nbCell1[1], nbCell1[2], lattice, rhoBar_, j_,
                thetaBar_);
            if (!lattice.get(iX, iY, iZ).getDynamics().velIsJ()) {
                j_ *= Descriptor<T>::invRho(rhoBar_);
            }
            if (noPenetration) {
                j_[direction1] = T();
            }
            edgeUtil.boundaryCompletion(j_);
            return;
        }

        // None of the neighbors is Neumann: there is no preferred direction, and extrapolation goes
        // from the diagonal bulk cell.
        Cell<T, Descriptor> &middleNeighbor =
            lattice.get(middleNbCell[0], middleNbCell[1], middleNbCell[2]);
        Array<T, 3> u;
        middleNeighbor.computeVelocity(u);
        if (noPenetration) {
            u[direction1] = T();
            u[direction2] = T();
        }
        edgeUtil.boundaryCompletion(u);
    }
    // At least one of them is Dirichlet: implement Dirichlet on the current cell, and inherit
    // the velocity from the Dirichlet neighbor.
    else
    {
        if (dyn1) {
            dyn1->edgeBoundaryCompletion(plane, normal1, normal2, iX, iY, iZ, lattice);
        }
        if (dyn2) {
            Cell<T, Descriptor> saveCell(lattice.get(iX, iY, iZ));
            dyn2->edgeBoundaryCompletion(plane, normal1, normal2, iX, iY, iZ, lattice);
            // If both neighbors are Dirichlet, take an average.
            if (dyn1) {
                Cell<T, Descriptor> &cell(lattice.get(iX, iY, iZ));
                for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                    cell[iPop] = (cell[iPop] + saveCell[iPop]) / 2.;
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>::cornerBoundaryCompletion(
    int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice)
{
    // Important: To understand periodicity, have a look at the comments in the
    // corresponding method in NLD_VelocityBoundaryDynamics3D.

    Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
    Array<plint, 3> nbCell(iX, iY, iZ);
    nbCell -= Array<plint, 3>(xNormal, yNormal, zNormal);

    // 1. To start with, detect periodic directions where the corner actually
    //    becomes edge, plane, or non-boundary.

    bool xPeriodic = !lattice.get(iX, iY - yNormal, iZ - zNormal).getDynamics().isBoundary();
    bool yPeriodic = !lattice.get(iX - xNormal, iY, iZ - zNormal).getDynamics().isBoundary();
    bool zPeriodic = !lattice.get(iX - xNormal, iY - yNormal, iZ).getDynamics().isBoundary();

    // Note that we only need to consider cases where corners become edges,
    //   because further degenerate cases (two or three periodic directions)
    //   are treated in the edge algorithm.
    if (xPeriodic) {
        edgeBoundaryCompletion(0, yNormal, zNormal, iX, iY, iZ, lattice);
        return;
    }
    if (yPeriodic) {
        edgeBoundaryCompletion(1, zNormal, xNormal, iX, iY, iZ, lattice);
        return;
    }
    if (zPeriodic) {
        edgeBoundaryCompletion(2, xNormal, yNormal, iX, iY, iZ, lattice);
        return;
    }

    // General case: none of the directions is periodic.

    Array<plint, 3> nbCell1(iX - xNormal, iY, iZ), nbCell2(iX, iY - yNormal, iZ),
        nbCell3(iX, iY, iZ - zNormal);
    Cell<T, Descriptor> &neighbor1 = lattice.get(nbCell1[0], nbCell1[1], nbCell1[2]);
    Cell<T, Descriptor> &neighbor2 = lattice.get(nbCell2[0], nbCell2[1], nbCell2[2]);
    Cell<T, Descriptor> &neighbor3 = lattice.get(nbCell3[0], nbCell3[1], nbCell3[2]);

    // Check if one of the three neighbors in plane directions is Dirichlet.
    NLD_VelocityBoundaryDynamics3D<T, Descriptor> *dyn1 =
        dynamic_cast<NLD_VelocityBoundaryDynamics3D<T, Descriptor> *>(&neighbor1.getDynamics());
    NLD_VelocityBoundaryDynamics3D<T, Descriptor> *dyn2 =
        dynamic_cast<NLD_VelocityBoundaryDynamics3D<T, Descriptor> *>(&neighbor2.getDynamics());
    NLD_VelocityBoundaryDynamics3D<T, Descriptor> *dyn3 =
        dynamic_cast<NLD_VelocityBoundaryDynamics3D<T, Descriptor> *>(&neighbor3.getDynamics());

    // General case: none of them is Dirichlet.
    if (!(dyn1 || dyn2 || dyn3)) {
        NLD_VelocityCornerUtil3D<T, Descriptor> cornerUtil(
            xNormal, yNormal, zNormal, iX, iY, iZ, lattice, *this);

        // Now, let's check if one of the edge neighbors is Neumann, while the other
        // two are not. In this case it will be concluded that we're in presence of an
        // outlet plane, and that all extrapolations are normal to this plane.
        bool same1 = isSameTypeAsI(nbCell1[0], nbCell1[1], nbCell1[2], lattice);
        bool same2 = isSameTypeAsI(nbCell2[0], nbCell2[1], nbCell2[2], lattice);
        bool same3 = isSameTypeAsI(nbCell3[0], nbCell3[1], nbCell3[2], lattice);
        // Neighbor 1 and 2 are Neumann, while neighbor 3 is something else --> the outlet normal
        // goes like normal3.
        if (same1 && same2 && !same3) {
            T rhoBar_, thetaBar_;
            Array<T, 3> j_;
            // Extrapolate the velocity from neighbor3, which is on an edge in z-direction.
            this->staticEdgeComputeMacroscopic(
                2, xNormal, yNormal, nbCell3[0], nbCell3[1], nbCell3[2], lattice, rhoBar_, j_,
                thetaBar_);
            if (!cell.getDynamics().velIsJ()) {
                j_ *= Descriptor<T>::invRho(rhoBar_);
            }
            if (noPenetration) {
                j_[2] = T();
            }
            cornerUtil.boundaryCompletion(j_);
            return;
        }
        // Neighbor 1 and 3 are Neumann, while neighbor 2 is something else --> the outlet normal
        // goes like normal2.
        if (same1 && same3 && !same2) {
            T rhoBar_, thetaBar_;
            Array<T, 3> j_;
            // Extrapolate the velocity from neighbor2, which is on an edge in y-direction.
            this->staticEdgeComputeMacroscopic(
                1, zNormal, xNormal, nbCell2[0], nbCell2[1], nbCell2[2], lattice, rhoBar_, j_,
                thetaBar_);
            if (!cell.getDynamics().velIsJ()) {
                j_ *= Descriptor<T>::invRho(rhoBar_);
            }
            if (noPenetration) {
                j_[1] = T();
            }
            cornerUtil.boundaryCompletion(j_);
            return;
        }
        // Neighbor 2 and 3 are Neumann, while neighbor 1 is something else --> the outlet normal
        // goes like normal1.
        if (same2 && same3 && !same1) {
            T rhoBar_, thetaBar_;
            Array<T, 3> j_;
            // Extrapolate the velocity from neighbor1, which is on an edge in x-direction.
            this->staticEdgeComputeMacroscopic(
                0, yNormal, zNormal, nbCell1[0], nbCell1[1], nbCell1[2], lattice, rhoBar_, j_,
                thetaBar_);
            if (!cell.getDynamics().velIsJ()) {
                j_ *= Descriptor<T>::invRho(rhoBar_);
            }
            if (noPenetration) {
                j_[0] = T();
            }
            cornerUtil.boundaryCompletion(j_);
            return;
        }

        // The number of Neumann neighbors is 0 or 3. In this case, there is no preferred direction,
        // and extrapolation is taken from the diagonal bulk cell.
        Cell<T, Descriptor> &middleNeighbor = lattice.get(nbCell[0], nbCell[1], nbCell[2]);
        Array<T, 3> u;
        middleNeighbor.computeVelocity(u);
        if (noPenetration) {
            u.resetToZero();
        }
        cornerUtil.boundaryCompletion(u);
    }
    // At least one of them is Dirichlet: implement Dirichlet on the current cell, and inherit
    // the velocity from the Dirichlet neighbor.
    else
    {
        Cell<T, Descriptor> tmpCell(lattice.get(iX, iY, iZ));
        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            tmpCell[iPop] = T();
        }
        plint numNeighbors = 0;
        if (dyn1) {
            dyn1->cornerBoundaryCompletion(xNormal, yNormal, zNormal, iX, iY, iZ, lattice);
            for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                tmpCell[iPop] += cell[iPop];
            }
            ++numNeighbors;
        }
        if (dyn2) {
            dyn2->cornerBoundaryCompletion(xNormal, yNormal, zNormal, iX, iY, iZ, lattice);
            for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                tmpCell[iPop] += cell[iPop];
            }
            ++numNeighbors;
        }
        if (dyn3) {
            dyn3->cornerBoundaryCompletion(xNormal, yNormal, zNormal, iX, iY, iZ, lattice);
            for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                tmpCell[iPop] += cell[iPop];
            }
            ++numNeighbors;
        }
        // If more than one neighbor is Dirichlet, take an average.
        if (numNeighbors > 1) {
            for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                cell[iPop] = tmpCell[iPop] / (T)numNeighbors;
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
bool NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>::isSameTypeAsI(
    plint iX, plint iY, plint iZ, BlockLattice3D<T, Descriptor> &lattice)
{
    NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor> *neumDyn =
        dynamic_cast<NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor> *>(
            &lattice.get(iX, iY, iZ).getDynamics());
    return neumDyn;
}

/* *************** Class NLD_OutflowBoundaryDynamics3D ****** */

template <typename T, template <typename U> class Descriptor>
int NLD_OutflowBoundaryDynamics3D<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, NLD_OutflowBoundaryDynamics3D<T, Descriptor> >(
        std::string("NLD_Outflow"));

template <typename T, template <typename U> class Descriptor>
NLD_OutflowBoundaryDynamics3D<T, Descriptor>::NLD_OutflowBoundaryDynamics3D(
    Dynamics<T, Descriptor> *baseDynamics_) :
    NonLocalBoundaryDynamics3D<T, Descriptor>(baseDynamics_)
{ }

template <typename T, template <typename U> class Descriptor>
NLD_OutflowBoundaryDynamics3D<T, Descriptor>::NLD_OutflowBoundaryDynamics3D(
    HierarchicUnserializer &unserializer) :
    NonLocalBoundaryDynamics3D<T, Descriptor>(0)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
int NLD_OutflowBoundaryDynamics3D<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void NLD_OutflowBoundaryDynamics3D<T, Descriptor>::prepareCollision(Cell<T, Descriptor> &cell)
{ }

template <typename T, template <typename U> class Descriptor>
NLD_OutflowBoundaryDynamics3D<T, Descriptor> *NLD_OutflowBoundaryDynamics3D<T, Descriptor>::clone()
    const
{
    return new NLD_OutflowBoundaryDynamics3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void NLD_OutflowBoundaryDynamics3D<T, Descriptor>::planeComputeMacroscopic(
    int direction, int orientation, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar)
{
    Array<plint, 3> neighb(iX, iY, iZ);
    neighb[direction] -= orientation;
    Cell<T, Descriptor> &cellFrom = lattice.get(neighb[0], neighb[1], neighb[2]);
    Cell<T, Descriptor> cellTo(lattice.get(iX, iY, iZ));
    flatWall::extrapolatePopulationsFixPressure(cellFrom, cellTo, direction, orientation);
    cellTo.getDynamics().computeRhoBarJ(cellTo, rhoBar, j);
    thetaBar = T();
}

template <typename T, template <typename U> class Descriptor>
void NLD_OutflowBoundaryDynamics3D<T, Descriptor>::edgeComputeMacroscopic(
    int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar)
{
    NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>(this->getBaseDynamics().clone())
        .edgeComputeMacroscopic(plane, normal1, normal2, iX, iY, iZ, lattice, rhoBar, j, thetaBar);
}

template <typename T, template <typename U> class Descriptor>
void NLD_OutflowBoundaryDynamics3D<T, Descriptor>::cornerComputeMacroscopic(
    int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar)
{
    NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>(this->getBaseDynamics().clone())
        .cornerComputeMacroscopic(
            xNormal, yNormal, zNormal, iX, iY, iZ, lattice, rhoBar, j, thetaBar);
}

template <typename T, template <typename U> class Descriptor>
void NLD_OutflowBoundaryDynamics3D<T, Descriptor>::planeBoundaryCompletion(
    int direction, int orientation, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice)
{
    flatWall::extrapolatePopulationsFixPressure(lattice, iX, iY, iZ, direction, orientation);
}

template <typename T, template <typename U> class Descriptor>
void NLD_OutflowBoundaryDynamics3D<T, Descriptor>::edgeBoundaryCompletion(
    int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice)
{
    NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>(this->getBaseDynamics().clone())
        .edgeBoundaryCompletion(plane, normal1, normal2, iX, iY, iZ, lattice);
}

template <typename T, template <typename U> class Descriptor>
void NLD_OutflowBoundaryDynamics3D<T, Descriptor>::cornerBoundaryCompletion(
    int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice)
{
    NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>(this->getBaseDynamics().clone())
        .cornerBoundaryCompletion(xNormal, yNormal, zNormal, iX, iY, iZ, lattice);
}

template <typename T, template <typename U> class Descriptor>
Dynamics<T, Descriptor> *cloneAndInsertNLDatTop(
    Dynamics<T, Descriptor> const &dynamics, NonLocalBoundaryDynamics3D<T, Descriptor> *newTop)
{
    Dynamics<T, Descriptor> *clonedDynamics = dynamics.clone();
    Dynamics<T, Descriptor> *current = clonedDynamics;
    Dynamics<T, Descriptor> *topMost = 0;
    bool isComposite = false;
    do {
        if (current->isNonLocal()) {
            topMost = 0;
        } else if (!topMost) {
            topMost = current;
        }
        isComposite = current->isComposite();
        if (isComposite) {
            current =
                &(dynamic_cast<CompositeDynamics<T, Descriptor> *>(current)->getBaseDynamics());
        }
    } while (isComposite);
    PLB_ASSERT(topMost);
    Dynamics<T, Descriptor> *result = cloneAndInsertAtTopDynamics(*topMost, newTop);
    delete clonedDynamics;
    return result;
}

/* *************** Class NLD_AD_DirichletDynamics3D ****** */

template <typename T, template <typename U> class Descriptor>
int NLD_AD_DirichletDynamics3D<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, NLD_AD_DirichletDynamics3D<T, Descriptor> >(
        std::string("NLD_AD_Dirichlet3D"));

template <typename T, template <typename U> class Descriptor>
NLD_AD_DirichletDynamics3D<T, Descriptor>::NLD_AD_DirichletDynamics3D(
    Dynamics<T, Descriptor> *baseDynamics_) :
    NonLocalBoundaryDynamics3D<T, Descriptor>(baseDynamics_)
{ }

template <typename T, template <typename U> class Descriptor>
NLD_AD_DirichletDynamics3D<T, Descriptor>::NLD_AD_DirichletDynamics3D(
    HierarchicUnserializer &unserializer) :
    NonLocalBoundaryDynamics3D<T, Descriptor>(0)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
void NLD_AD_DirichletDynamics3D<T, Descriptor>::defineDensity(Cell<T, Descriptor> &cell, T density)
{
    rho = density;
}

template <typename T, template <typename U> class Descriptor>
T NLD_AD_DirichletDynamics3D<T, Descriptor>::computeDensity(Cell<T, Descriptor> const &cell) const
{
    return rho;
}

template <typename T, template <typename U> class Descriptor>
void NLD_AD_DirichletDynamics3D<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    NonLocalBoundaryDynamics3D<T, Descriptor>::serialize(serializer);
    serializer.addValue(rho);
}

template <typename T, template <typename U> class Descriptor>
void NLD_AD_DirichletDynamics3D<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    NonLocalBoundaryDynamics3D<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(rho);
}

template <typename T, template <typename U> class Descriptor>
int NLD_AD_DirichletDynamics3D<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void NLD_AD_DirichletDynamics3D<T, Descriptor>::prepareCollision(Cell<T, Descriptor> &cell)
{ }

template <typename T, template <typename U> class Descriptor>
NLD_AD_DirichletDynamics3D<T, Descriptor> *NLD_AD_DirichletDynamics3D<T, Descriptor>::clone() const
{
    return new NLD_AD_DirichletDynamics3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void NLD_AD_DirichletDynamics3D<T, Descriptor>::planeComputeMacroscopic(
    int direction, int orientation, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar)
{
    rhoBar = Descriptor<T>::rhoBar(rho);
    // Here we decide to return j^eq instead of j.
    advectionDiffusionMomentTemplates<T, Descriptor>::get_jEq(lattice.get(iX, iY, iZ), rhoBar, j);
    thetaBar = T();
}

template <typename T, template <typename U> class Descriptor>
void NLD_AD_DirichletDynamics3D<T, Descriptor>::edgeComputeMacroscopic(
    int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar)
{
    rhoBar = Descriptor<T>::rhoBar(rho);
    // Here we decide to return j^eq instead of j.
    advectionDiffusionMomentTemplates<T, Descriptor>::get_jEq(lattice.get(iX, iY, iZ), rhoBar, j);
    thetaBar = T();
}

template <typename T, template <typename U> class Descriptor>
void NLD_AD_DirichletDynamics3D<T, Descriptor>::cornerComputeMacroscopic(
    int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar)
{
    rhoBar = Descriptor<T>::rhoBar(rho);
    // Here we decide to return j^eq instead of j.
    advectionDiffusionMomentTemplates<T, Descriptor>::get_jEq(lattice.get(iX, iY, iZ), rhoBar, j);
    thetaBar = T();
}

template <typename T, template <typename U> class Descriptor>
void NLD_AD_DirichletDynamics3D<T, Descriptor>::planeBoundaryCompletion(
    int direction, int orientation, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice)
{
    typedef advectionDiffusionDynamicsTemplates<T, Descriptor> dynTempl;
    Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
    T rhoBar;
    Array<T, 3> jEq;
    T dummyThetaBar;
    planeComputeMacroscopic(
        direction, orientation, iX, iY, iZ, lattice, rhoBar, jEq, dummyThetaBar);

    for (int iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        if (Descriptor<T>::c[iPop][direction] == -orientation) {
            plint opp = indexTemplates::opposite<Descriptor<T> >(iPop);
            cell[iPop] = dynTempl::bgk_ma1_equilibrium(iPop, rhoBar, jEq)
                         + dynTempl::bgk_ma1_equilibrium(opp, rhoBar, jEq) - cell[opp];
        }
    }

    Array<T, 3> j;
    this->getBaseDynamics().computeRhoBarJ(cell, rhoBar, j);
    dynTempl::regularize(cell, rhoBar, j, jEq);
}

template <typename T, template <typename U> class Descriptor>
void NLD_AD_DirichletDynamics3D<T, Descriptor>::edgeBoundaryCompletion(
    int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice)
{
    //  x
    //  o - -
    //  |
    //  |   A few comments on periodicity. On the scheme on the left, assume that the x-axis is
    //  |   direction1, and the y-axis is direction2. Further, assume that direction2 is periodic
    //      and direction1 is not periodic. Then, the cell "x" must be in the envelope, and it
    //      necessarily has NLD dynamics (because envelopes are automatically updated with dynamic
    //      content when the domain is made periodic). Also, remember that all NLD dynamics objects
    //      are the same, no matter if they are plane, corner or edge, so there is no possible
    //      confusion here. To obtain periodicity, the user will set the "-" cells to background-
    //      dynamics, and not to NLD dynamics. That's how, in the code below, periodicity is
    //      detected.
    int direction1 = (plane + 1) % 3;
    int direction2 = (plane + 2) % 3;
    Array<plint, 3> nbCell1(iX, iY, iZ), nbCell2(iX, iY, iZ);
    nbCell1[direction1] -= normal1;
    nbCell2[direction2] -= normal2;

    // 1. Detect case where a direction is periodic.
    bool isWall1 = lattice.get(nbCell1[0], nbCell1[1], nbCell1[2]).getDynamics().isBoundary();
    bool isWall2 = lattice.get(nbCell2[0], nbCell2[1], nbCell2[2]).getDynamics().isBoundary();

    // If both directions are periodic, there is no need for a boundary condition.
    if (!isWall1 && !isWall2)
        return;
    // If exactly one direction is periodic, then the edge behaves like a plane.
    if (!isWall1) {
        planeBoundaryCompletion(direction1, normal1, iX, iY, iZ, lattice);
        return;
    } else if (!isWall2) {
        planeBoundaryCompletion(direction2, normal2, iX, iY, iZ, lattice);
        return;
    }

    // 2. Implement general edge algorithm.
    typedef advectionDiffusionDynamicsTemplates<T, Descriptor> dynTempl;
    typedef Descriptor<T> D;
    Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
    T rhoBar;
    Array<T, 3> jEq;
    T dummyThetaBar;
    edgeComputeMacroscopic(
        plane, normal1, normal2, iX, iY, iZ, lattice, rhoBar, jEq, dummyThetaBar);

    for (int iPop = 0; iPop < D::q; ++iPop) {
        if ((D::c[iPop][direction1] == -normal1 || D::c[iPop][direction1] == 0)
            && (D::c[iPop][direction2] == -normal2 || D::c[iPop][direction2] == 0)
            && (!(D::c[iPop][direction1] == 0 && D::c[iPop][direction2] == 0)))
        {
            plint opp = indexTemplates::opposite<Descriptor<T> >(iPop);
            cell[iPop] = dynTempl::bgk_ma1_equilibrium(iPop, rhoBar, jEq)
                         + dynTempl::bgk_ma1_equilibrium(opp, rhoBar, jEq) - cell[opp];
        }
    }

    Array<T, 3> j;
    this->getBaseDynamics().computeRhoBarJ(cell, rhoBar, j);
    dynTempl::regularize(cell, rhoBar, j, jEq);
}

template <typename T, template <typename U> class Descriptor>
void NLD_AD_DirichletDynamics3D<T, Descriptor>::cornerBoundaryCompletion(
    int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice)
{
    // Important: to understand periodicity, look at the comment in the edgeBoundaryCompletion
    // function right above.

    Array<plint, 3> middleNbCell(iX, iY, iZ), nbCellX(iX, iY, iZ), nbCellY(iX, iY, iZ),
        nbCellZ(iX, iY, iZ);
    middleNbCell -= Array<plint, 3>(xNormal, yNormal, zNormal);
    nbCellX[0] -= xNormal;
    nbCellY[1] -= yNormal;
    nbCellZ[2] -= zNormal;

    // 1. To start width, detect periodic directions where the corner actually
    //    becomes edge, plane, or non-boundary.

    bool xPeriodic = !lattice.get(iX, iY - yNormal, iZ - zNormal).getDynamics().isBoundary();
    bool yPeriodic = !lattice.get(iX - xNormal, iY, iZ - zNormal).getDynamics().isBoundary();
    bool zPeriodic = !lattice.get(iX - xNormal, iY - yNormal, iZ).getDynamics().isBoundary();

    // Note that we only need to consider cases where corners become edges,
    //   because further degenerate cases (two or three periodic directions)
    //   are treated in the edge algorithm.
    if (xPeriodic) {
        edgeBoundaryCompletion(0, yNormal, zNormal, iX, iY, iZ, lattice);
        return;
    }
    if (yPeriodic) {
        edgeBoundaryCompletion(1, zNormal, xNormal, iX, iY, iZ, lattice);
        return;
    }
    if (zPeriodic) {
        edgeBoundaryCompletion(2, xNormal, yNormal, iX, iY, iZ, lattice);
        return;
    }

    // 2. Implement generic case.
    typedef advectionDiffusionDynamicsTemplates<T, Descriptor> dynTempl;
    typedef Descriptor<T> D;
    Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
    T rhoBar;
    Array<T, 3> jEq;
    T dummyThetaBar;
    cornerComputeMacroscopic(
        xNormal, yNormal, zNormal, iX, iY, iZ, lattice, rhoBar, jEq, dummyThetaBar);

    for (int iPop = 0; iPop < D::q; ++iPop) {
        if ((D::c[iPop][0] == -xNormal || D::c[iPop][0] == 0)
            && (D::c[iPop][1] == -yNormal || D::c[iPop][1] == 0)
            && (D::c[iPop][2] == -zNormal || D::c[iPop][2] == 0)
            && (!(D::c[iPop][0] == 0 && D::c[iPop][1] == 0) && D::c[iPop][2] == 0))
        {
            plint opp = indexTemplates::opposite<Descriptor<T> >(iPop);
            cell[iPop] = dynTempl::bgk_ma1_equilibrium(iPop, rhoBar, jEq)
                         + dynTempl::bgk_ma1_equilibrium(opp, rhoBar, jEq) - cell[opp];
        }
    }

    Array<T, 3> j;
    this->getBaseDynamics().computeRhoBarJ(cell, rhoBar, j);
    dynTempl::regularize(cell, rhoBar, j, jEq);
}

/* *************** Class NLD_AD_NeumannBoundaryDynamics3D ****** */

template <typename T, template <typename U> class Descriptor>
int NLD_AD_NeumannBoundaryDynamics3D<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, NLD_AD_NeumannBoundaryDynamics3D<T, Descriptor> >(
        std::string("NLD_AD_Neumann"));

template <typename T, template <typename U> class Descriptor>
NLD_AD_NeumannBoundaryDynamics3D<T, Descriptor>::NLD_AD_NeumannBoundaryDynamics3D(
    Dynamics<T, Descriptor> *baseDynamics_) :
    NonLocalBoundaryDynamics3D<T, Descriptor>(baseDynamics_)
{ }

template <typename T, template <typename U> class Descriptor>
NLD_AD_NeumannBoundaryDynamics3D<T, Descriptor>::NLD_AD_NeumannBoundaryDynamics3D(
    HierarchicUnserializer &unserializer) :
    NonLocalBoundaryDynamics3D<T, Descriptor>(0)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
int NLD_AD_NeumannBoundaryDynamics3D<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void NLD_AD_NeumannBoundaryDynamics3D<T, Descriptor>::prepareCollision(Cell<T, Descriptor> &cell)
{ }

template <typename T, template <typename U> class Descriptor>
NLD_AD_NeumannBoundaryDynamics3D<T, Descriptor>
    *NLD_AD_NeumannBoundaryDynamics3D<T, Descriptor>::clone() const
{
    return new NLD_AD_NeumannBoundaryDynamics3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void NLD_AD_NeumannBoundaryDynamics3D<T, Descriptor>::planeComputeMacroscopic(
    int direction, int orientation, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar)
{
    NLD_VelocityPlaneUtil3D<T, Descriptor> planeUtil(
        direction, orientation, iX, iY, iZ, lattice, *this);
    rhoBar = Descriptor<T>::rhoBar(planeUtil.extrapolateDensity());
    // Here we decide to return j^eq instead of j.
    advectionDiffusionMomentTemplates<T, Descriptor>::get_jEq(lattice.get(iX, iY, iZ), rhoBar, j);
    thetaBar = T();
}

template <typename T, template <typename U> class Descriptor>
void NLD_AD_NeumannBoundaryDynamics3D<T, Descriptor>::edgeComputeMacroscopic(
    int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar)
{
    // General philosophy: The cell being treated here implements a Neumann condition on an
    // edge. Now, if one of the two neighbors in the plane directions is Dirichlet, it is
    // better to treat the current cell as Dirichlet as well, because this is numerically
    // more stable, and this is probably what the user wanted to happen anyway.
    //
    // If none of the neighbors is Dirichlet, then it could be that this Neumann edge node
    // is part of an outlet, and that it is in contact with, say, a free-slip wall. In this
    // case, it is important that the extrapolation is done along the normal of the outlet
    // plane, i.e. in direction of this free-slip wall node.
    //
    // Finally, if all neighbors are Neumann, there exists no preferred direction. In this
    // case, extrapolation is performed from the bulk, along a diagonal direction.

    // This code is non-thermal. Default-assign the temperature.
    thetaBar = T();

    int direction1 = (plane + 1) % 3;
    int direction2 = (plane + 2) % 3;

    Array<plint, 3> nbCell1(iX, iY, iZ), nbCell2(iX, iY, iZ);
    nbCell1[direction1] -= normal1;
    nbCell2[direction2] -= normal2;
    Cell<T, Descriptor> &neighbor1 = lattice.get(nbCell1[0], nbCell1[1], nbCell1[2]);
    Cell<T, Descriptor> &neighbor2 = lattice.get(nbCell2[0], nbCell2[1], nbCell2[2]);

    // Check if the plane neighbors are Dirichlet.
    NLD_VelocityBoundaryDynamics3D<T, Descriptor> *dyn1 =
        dynamic_cast<NLD_VelocityBoundaryDynamics3D<T, Descriptor> *>(&neighbor1.getDynamics());
    NLD_VelocityBoundaryDynamics3D<T, Descriptor> *dyn2 =
        dynamic_cast<NLD_VelocityBoundaryDynamics3D<T, Descriptor> *>(&neighbor2.getDynamics());
    // General case: none of them is Dirichlet. Implement a general Neumann condition
    // by extrapolating along the normal of the outlet plane, or, if no preferred
    // direction exists, from the bulk node in diagonal direction.
    if (!(dyn1 || dyn2)) {
        // Now, let's check if one of the plane neighbors is Neumann, while the other
        // is not. In this case it will be concluded that we're in presence of an
        // outlet plane, and that all extrapolations are normal to this plane.
        bool same1 = isSameTypeAsI(nbCell1[0], nbCell1[1], nbCell1[2], lattice);
        bool same2 = isSameTypeAsI(nbCell2[0], nbCell2[1], nbCell2[2], lattice);
        // Neighbor1 is Neumann, while neighbor2 is something else --> the outlet normal goes like
        // normal2.
        if (same1 && !same2) {
            // Extrapolate the density from neighbor2, which is on a plane with normal normal1.
            this->staticPlaneComputeMacroscopic(
                direction1, normal1, nbCell2[0], nbCell2[1], nbCell2[2], lattice, rhoBar, j,
                thetaBar);
            return;
        }
        // Neighbor2 is Neumann, while neighbor1 is something else --> the outlet normal goes like
        // normal1.
        if (same2 && !same1) {
            // Extrapolate the density from neighbor1, which is on a plane with normal normal2.
            this->staticPlaneComputeMacroscopic(
                direction2, normal2, nbCell1[0], nbCell1[1], nbCell1[2], lattice, rhoBar, j,
                thetaBar);
            return;
        }
        // Either both or none of them is Neumann. In this case, there is no preferred direction,
        // and extrapolation is taken from the diagonal bulk cell.
        Array<plint, 3> middleNbCell(iX, iY, iZ);
        middleNbCell[direction1] -= normal1;
        middleNbCell[direction2] -= normal2;
        Cell<T, Descriptor> &middleNeighbor =
            lattice.get(middleNbCell[0], middleNbCell[1], middleNbCell[2]);
        rhoBar = Descriptor<T>::rhoBar(middleNeighbor.computeDensity());
    }
    // At least one of the two plane neighbors is Dirichlet. Treat the current cell
    // as Dirichlet as well, inheriting the density value from the neighbor(s).
    else
    {
        T rhoBar1, rhoBar2;
        Array<T, 3> j1, j2;
        // Inherit Dirichlet condition from neighbor1, if neighbor1 is Dirichlet.
        if (dyn1) {
            dyn1->edgeComputeMacroscopic(
                plane, normal1, normal2, iX, iY, iZ, lattice, rhoBar1, j1, thetaBar);
            rhoBar = rhoBar1;
        }
        // Inherit Dirichlet condition from neighbor2, if neighbor2 is Dirichlet.
        if (dyn2) {
            dyn2->edgeComputeMacroscopic(
                plane, normal1, normal2, iX, iY, iZ, lattice, rhoBar2, j2, thetaBar);
            rhoBar = rhoBar2;
        }
        // If both neighbors are Dirichlet, take an average.
        if (dyn1 && dyn2) {
            rhoBar = (rhoBar1 + rhoBar2) / (T)2;
        }
    }
    advectionDiffusionMomentTemplates<T, Descriptor>::get_jEq(lattice.get(iX, iY, iZ), rhoBar, j);
}

template <typename T, template <typename U> class Descriptor>
void NLD_AD_NeumannBoundaryDynamics3D<T, Descriptor>::cornerComputeMacroscopic(
    int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar)
{
    // General philosophy: The cell being treated here implements a Neumann condition on an
    // corner. Now, if one of the three neighbors in the edge directions is Dirichlet, it is
    // better to treat the current cell as Dirichlet as well, because this is numerically
    // more stable, and this is probably what the user wanted to happen anyway.
    //
    // If none of the neighbors is Dirichlet, then it could be that this Neumann corner node
    // is part of an outlet, and that it is in contact with, say, a free-slip wall. In this
    // case, it is important that the extrapolation is done along the normal of the outlet
    // plane, i.e. in direction of this free-slip wall node.
    //
    // Finally, if either all neighbors or only one neighbor are Neumann, there exists no
    // preferred direction. In this case, extrapolation is performed from the bulk, along
    // a diagonal direction.

    // This code is non-thermal. Default-assign the temperature.
    thetaBar = T();

    Array<plint, 3> nbCell1(iX - xNormal, iY, iZ), nbCell2(iX, iY - yNormal, iZ),
        nbCell3(iX, iY, iZ - zNormal);
    Cell<T, Descriptor> &neighbor1 = lattice.get(nbCell1[0], nbCell1[1], nbCell1[2]);
    Cell<T, Descriptor> &neighbor2 = lattice.get(nbCell2[0], nbCell2[1], nbCell2[2]);
    Cell<T, Descriptor> &neighbor3 = lattice.get(nbCell3[0], nbCell3[1], nbCell3[2]);

    NLD_VelocityBoundaryDynamics3D<T, Descriptor> *dyn1 =
        dynamic_cast<NLD_VelocityBoundaryDynamics3D<T, Descriptor> *>(&neighbor1.getDynamics());
    NLD_VelocityBoundaryDynamics3D<T, Descriptor> *dyn2 =
        dynamic_cast<NLD_VelocityBoundaryDynamics3D<T, Descriptor> *>(&neighbor2.getDynamics());
    NLD_VelocityBoundaryDynamics3D<T, Descriptor> *dyn3 =
        dynamic_cast<NLD_VelocityBoundaryDynamics3D<T, Descriptor> *>(&neighbor3.getDynamics());

    // Worst-case: none of them is Dirichlet. Implement a general Neumann condition
    // by extrapolating from the bulk node in diagonal direction.
    if (!(dyn1 || dyn2 || dyn3)) {
        // Now, let's check if one of the edge neighbors is Neumann, while the other
        // two are not. In this case it will be concluded that we're in presence of an
        // outlet plane, and that all extrapolations are normal to this plane.
        bool same1 = isSameTypeAsI(nbCell1[0], nbCell1[1], nbCell1[2], lattice);
        bool same2 = isSameTypeAsI(nbCell2[0], nbCell2[1], nbCell2[2], lattice);
        bool same3 = isSameTypeAsI(nbCell3[0], nbCell3[1], nbCell3[2], lattice);
        // Neighbor 1 and 2 are Neumann, while neighbor 3 is something else --> the outlet normal
        // goes like normal3.
        if (same1 && same2 && !same3) {
            // Extrapolate the density from neighbor3, which is on an edge in z-direction.
            this->staticEdgeComputeMacroscopic(
                2, xNormal, yNormal, nbCell3[0], nbCell3[1], nbCell3[2], lattice, rhoBar, j,
                thetaBar);
            return;
        }
        // Neighbor 1 and 3 are Neumann, while neighbor 2 is something else --> the outlet normal
        // goes like normal2.
        if (same1 && same3 && !same2) {
            // Extrapolate the density from neighbor2, which is on an edge in y-direction.
            this->staticEdgeComputeMacroscopic(
                1, zNormal, xNormal, nbCell2[0], nbCell2[1], nbCell2[2], lattice, rhoBar, j,
                thetaBar);
            return;
        }
        // Neighbor 2 and 3 are Neumann, while neighbor 1 is something else --> the outlet normal
        // goes like normal1.
        if (same2 && same3 && !same1) {
            // Extrapolate the density from neighbor1, which is on an edge in x-direction.
            this->staticEdgeComputeMacroscopic(
                0, yNormal, zNormal, nbCell1[0], nbCell1[1], nbCell1[2], lattice, rhoBar, j,
                thetaBar);
            return;
        }

        // The number of Neumann neighbors is 0 or 3. In this case, there is no preferred direction,
        // and extrapolation is taken from the diagonal bulk cell.
        Array<plint, 3> nbCell(iX, iY, iZ);
        nbCell -= Array<plint, 3>(xNormal, yNormal, zNormal);
        Cell<T, Descriptor> &middleNeighbor = lattice.get(nbCell[0], nbCell[1], nbCell[2]);
        rhoBar = Descriptor<T>::rhoBar(middleNeighbor.computeDensity());
    }
    // At least one of the three edge neighbors is Dirichlet. Treat the current cell
    // as Dirichlet as well, inheriting the velocity value from the neighbor(s).
    else
    {
        rhoBar = T();
        plint numNeighbors = 0;
        // Inherit Dirichlet condition from neighbor1, if neighbor1 is Dirichlet.
        if (dyn1) {
            T rhoBar_;
            Array<T, 3> j_;
            dyn1->cornerComputeMacroscopic(
                xNormal, yNormal, zNormal, iX, iY, iZ, lattice, rhoBar_, j_, thetaBar);
            rhoBar += rhoBar_;
            ++numNeighbors;
        }
        // Inherit Dirichlet condition from neighbor2, if neighbor2 is Dirichlet.
        if (dyn2) {
            T rhoBar_;
            Array<T, 3> j_;
            dyn2->cornerComputeMacroscopic(
                xNormal, yNormal, zNormal, iX, iY, iZ, lattice, rhoBar_, j_, thetaBar);
            rhoBar += rhoBar_;
            ++numNeighbors;
        }
        // Inherit Dirichlet condition from neighbor3, if neighbor3 is Dirichlet.
        if (dyn3) {
            T rhoBar_;
            Array<T, 3> j_;
            dyn3->cornerComputeMacroscopic(
                xNormal, yNormal, zNormal, iX, iY, iZ, lattice, rhoBar_, j_, thetaBar);
            rhoBar += rhoBar_;
            ++numNeighbors;
        }
        // If more than one neighbor is Dirichlet, take an average.
        if (numNeighbors > 1) {
            rhoBar /= (T)numNeighbors;
        }
    }
    advectionDiffusionMomentTemplates<T, Descriptor>::get_jEq(lattice.get(iX, iY, iZ), rhoBar, j);
}

template <typename T, template <typename U> class Descriptor>
void NLD_AD_NeumannBoundaryDynamics3D<T, Descriptor>::planeBoundaryCompletion(
    int direction, int orientation, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice)
{
    typedef advectionDiffusionDynamicsTemplates<T, Descriptor> dynTempl;
    Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
    T rhoBar;
    Array<T, 3> jEq;
    T dummyThetaBar;
    planeComputeMacroscopic(
        direction, orientation, iX, iY, iZ, lattice, rhoBar, jEq, dummyThetaBar);

    for (int iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        if (Descriptor<T>::c[iPop][direction] == -orientation) {
            plint opp = indexTemplates::opposite<Descriptor<T> >(iPop);
            cell[iPop] = dynTempl::bgk_ma1_equilibrium(iPop, rhoBar, jEq)
                         + dynTempl::bgk_ma1_equilibrium(opp, rhoBar, jEq) - cell[opp];
        }
    }

    Array<T, 3> j;
    this->getBaseDynamics().computeRhoBarJ(cell, rhoBar, j);
    dynTempl::regularize(cell, rhoBar, j, jEq);
}

template <typename T, template <typename U> class Descriptor>
void NLD_AD_NeumannBoundaryDynamics3D<T, Descriptor>::edgeBoundaryCompletion(
    int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice)
{
    //  x
    //  o - -
    //  |
    //  |   A few comments on periodicity. On the scheme on the left, assume that the x-axis is
    //  |   direction1, and the y-axis is direction2. Further, assume that direction2 is periodic
    //      and direction1 is not periodic. Then, the cell "x" must be in the envelope, and it
    //      necessarily has NLD dynamics (because envelopes are automatically updated with dynamic
    //      content when the domain is made periodic). Also, remember that all NLD dynamics objects
    //      are the same, no matter if they are plane, corner or edge, so there is no possible
    //      confusion here. To obtain periodicity, the user will set the "-" cells to background-
    //      dynamics, and not to NLD dynamics. That's how, in the code below, periodicity is
    //      detected.
    int direction1 = (plane + 1) % 3;
    int direction2 = (plane + 2) % 3;
    Array<plint, 3> nbCell1(iX, iY, iZ), nbCell2(iX, iY, iZ);
    nbCell1[direction1] -= normal1;
    nbCell2[direction2] -= normal2;

    // 1. Detect case where a direction is periodic.
    bool isWall1 = lattice.get(nbCell1[0], nbCell1[1], nbCell1[2]).getDynamics().isBoundary();
    bool isWall2 = lattice.get(nbCell2[0], nbCell2[1], nbCell2[2]).getDynamics().isBoundary();

    // If both directions are periodic, there is no need for a boundary condition.
    if (!isWall1 && !isWall2)
        return;
    // If exactly one direction is periodic, then the edge behaves like a plane.
    if (!isWall1) {
        planeBoundaryCompletion(direction1, normal1, iX, iY, iZ, lattice);
        return;
    } else if (!isWall2) {
        planeBoundaryCompletion(direction2, normal2, iX, iY, iZ, lattice);
        return;
    }

    // 2. Implement general edge algorithm.
    typedef advectionDiffusionDynamicsTemplates<T, Descriptor> dynTempl;
    typedef Descriptor<T> D;
    Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
    T rhoBar;
    Array<T, 3> jEq;
    T dummyThetaBar;
    edgeComputeMacroscopic(
        plane, normal1, normal2, iX, iY, iZ, lattice, rhoBar, jEq, dummyThetaBar);

    for (int iPop = 0; iPop < D::q; ++iPop) {
        if ((D::c[iPop][direction1] == -normal1 || D::c[iPop][direction1] == 0)
            && (D::c[iPop][direction2] == -normal2 || D::c[iPop][direction2] == 0)
            && (!(D::c[iPop][direction1] == 0 && D::c[iPop][direction2] == 0)))
        {
            plint opp = indexTemplates::opposite<Descriptor<T> >(iPop);
            cell[iPop] = dynTempl::bgk_ma1_equilibrium(iPop, rhoBar, jEq)
                         + dynTempl::bgk_ma1_equilibrium(opp, rhoBar, jEq) - cell[opp];
        }
    }

    Array<T, 3> j;
    this->getBaseDynamics().computeRhoBarJ(cell, rhoBar, j);
    dynTempl::regularize(cell, rhoBar, j, jEq);
}

template <typename T, template <typename U> class Descriptor>
void NLD_AD_NeumannBoundaryDynamics3D<T, Descriptor>::cornerBoundaryCompletion(
    int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
    BlockLattice3D<T, Descriptor> &lattice)
{
    // Important: to understand periodicity, look at the comment in the edgeBoundaryCompletion
    // function right above.

    Array<plint, 3> middleNbCell(iX, iY, iZ), nbCellX(iX, iY, iZ), nbCellY(iX, iY, iZ),
        nbCellZ(iX, iY, iZ);
    middleNbCell -= Array<plint, 3>(xNormal, yNormal, zNormal);
    nbCellX[0] -= xNormal;
    nbCellY[1] -= yNormal;
    nbCellZ[2] -= zNormal;

    // 1. To start width, detect periodic directions where the corner actually
    //    becomes edge, plane, or non-boundary.

    bool xPeriodic = !lattice.get(iX, iY - yNormal, iZ - zNormal).getDynamics().isBoundary();
    bool yPeriodic = !lattice.get(iX - xNormal, iY, iZ - zNormal).getDynamics().isBoundary();
    bool zPeriodic = !lattice.get(iX - xNormal, iY - yNormal, iZ).getDynamics().isBoundary();

    // Note that we only need to consider cases where corners become edges,
    //   because further degenerate cases (two or three periodic directions)
    //   are treated in the edge algorithm.
    if (xPeriodic) {
        edgeBoundaryCompletion(0, yNormal, zNormal, iX, iY, iZ, lattice);
        return;
    }
    if (yPeriodic) {
        edgeBoundaryCompletion(1, zNormal, xNormal, iX, iY, iZ, lattice);
        return;
    }
    if (zPeriodic) {
        edgeBoundaryCompletion(2, xNormal, yNormal, iX, iY, iZ, lattice);
        return;
    }

    // 2. Implement generic case.
    typedef advectionDiffusionDynamicsTemplates<T, Descriptor> dynTempl;
    typedef Descriptor<T> D;
    Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
    T rhoBar;
    Array<T, 3> jEq;
    T dummyThetaBar;
    cornerComputeMacroscopic(
        xNormal, yNormal, zNormal, iX, iY, iZ, lattice, rhoBar, jEq, dummyThetaBar);

    for (int iPop = 0; iPop < D::q; ++iPop) {
        if ((D::c[iPop][0] == -xNormal || D::c[iPop][0] == 0)
            && (D::c[iPop][1] == -yNormal || D::c[iPop][1] == 0)
            && (D::c[iPop][2] == -zNormal || D::c[iPop][2] == 0)
            && (!(D::c[iPop][0] == 0 && D::c[iPop][1] == 0) && D::c[iPop][2] == 0))
        {
            plint opp = indexTemplates::opposite<Descriptor<T> >(iPop);
            cell[iPop] = dynTempl::bgk_ma1_equilibrium(iPop, rhoBar, jEq)
                         + dynTempl::bgk_ma1_equilibrium(opp, rhoBar, jEq) - cell[opp];
        }
    }

    Array<T, 3> j;
    this->getBaseDynamics().computeRhoBarJ(cell, rhoBar, j);
    dynTempl::regularize(cell, rhoBar, j, jEq);
}

template <typename T, template <typename U> class Descriptor>
bool NLD_AD_NeumannBoundaryDynamics3D<T, Descriptor>::isSameTypeAsI(
    plint iX, plint iY, plint iZ, BlockLattice3D<T, Descriptor> &lattice)
{
    NLD_AD_NeumannBoundaryDynamics3D<T, Descriptor> *neumDyn =
        dynamic_cast<NLD_AD_NeumannBoundaryDynamics3D<T, Descriptor> *>(
            &lattice.get(iX, iY, iZ).getDynamics());
    return neumDyn;
}

}  // namespace plb

// namespace plb

#endif  // NLD_BOUNDARY_DYNAMICS_3D_HH
