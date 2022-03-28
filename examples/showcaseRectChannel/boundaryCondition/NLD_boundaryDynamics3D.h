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
 * non-local dynamics (NLD) boundaries -- header file.
 */
#ifndef NLD_BOUNDARY_DYNAMICS_3D_H
#define NLD_BOUNDARY_DYNAMICS_3D_H

#include "boundaryCondition/regularizedBoundaryDynamics.h"
#include "core/array.h"
#include "core/globalDefs.h"

namespace plb {

/// Generic interface to any kind of non-local dynamics (NLD) boundary algorithms
/// (may be dirichlet, neumann, etc.).
template <typename T, template <typename U> class Descriptor>
class NonLocalBoundaryDynamics3D : public CompositeDynamics<T, Descriptor> {
public:
    NonLocalBoundaryDynamics3D(Dynamics<T, Descriptor> *baseDynamics_);
    virtual bool isBoundary() const;
    virtual bool isNonLocal() const;
    virtual NonLocalBoundaryDynamics3D<T, Descriptor> *clone() const = 0;

    virtual void planeBoundaryCompletion(
        int direction, int orientation, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice) = 0;
    virtual void edgeBoundaryCompletion(
        int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice) = 0;
    virtual void cornerBoundaryCompletion(
        int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice) = 0;

    virtual void planeComputeMacroscopic(
        int direction, int orientation, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar) = 0;
    virtual void edgeComputeMacroscopic(
        int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar) = 0;
    virtual void cornerComputeMacroscopic(
        int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar) = 0;

public:
    /// Convenience function to compute macroscopic variables on a plane NLD dynamics
    /// object. Returns default value if the dynamics is not plane NLD.
    static void staticPlaneComputeMacroscopic(
        int direction, int orientation, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar)
    {
        NonLocalBoundaryDynamics3D<T, Descriptor> *dynamics =
            dynamic_cast<NonLocalBoundaryDynamics3D<T, Descriptor> *>(
                &lattice.get(iX, iY, iZ).getDynamics());
        if (dynamics) {
            dynamics->planeComputeMacroscopic(
                direction, orientation, iX, iY, iZ, lattice, rhoBar, j, thetaBar);
        } else {
            rhoBar = T();
            j.resetToZero();
            thetaBar = T();
        }
    }
    /// Convenience function to compute macroscopic variables on a edge NLD dynamics
    /// object. Returns default value if the dynamics is not edge NLD.
    static void staticEdgeComputeMacroscopic(
        int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar)
    {
        NonLocalBoundaryDynamics3D<T, Descriptor> *dynamics =
            dynamic_cast<NonLocalBoundaryDynamics3D<T, Descriptor> *>(
                &lattice.get(iX, iY, iZ).getDynamics());
        if (dynamics) {
            dynamics->edgeComputeMacroscopic(
                plane, normal1, normal2, iX, iY, iZ, lattice, rhoBar, j, thetaBar);
        } else {
            rhoBar = T();
            j.resetToZero();
            thetaBar = T();
        }
    }
    /// Convenience function to compute macroscopic variables on a corner NLD dynamics
    /// object. Returns default value if the dynamics is not corner NLD.
    static void staticCornerComputeMacroscopic(
        int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar)
    {
        NonLocalBoundaryDynamics3D<T, Descriptor> *dynamics =
            dynamic_cast<NonLocalBoundaryDynamics3D<T, Descriptor> *>(
                &lattice.get(iX, iY, iZ).getDynamics());
        if (dynamics) {
            dynamics->cornerComputeMacroscopic(
                xNormal, yNormal, zNormal, iX, iY, iZ, lattice, rhoBar, j, thetaBar);
        } else {
            rhoBar = T();
            j.resetToZero();
            thetaBar = T();
        }
    }
};

/// Dirichlet version of the generic NLD boundary class.
template <typename T, template <typename U> class Descriptor>
class NLD_VelocityBoundaryDynamics3D : public NonLocalBoundaryDynamics3D<T, Descriptor> {
public:
    NLD_VelocityBoundaryDynamics3D(Dynamics<T, Descriptor> *baseDynamics_);
    NLD_VelocityBoundaryDynamics3D(HierarchicUnserializer &unserializer);
    virtual void defineVelocity(
        Cell<T, Descriptor> &cell, Array<T, Descriptor<T>::d> const &velocity_);
    virtual void computeVelocity(
        Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &velocity_) const;
    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);
    virtual NLD_VelocityBoundaryDynamics3D<T, Descriptor> *clone() const;
    virtual int getId() const;
    virtual void prepareCollision(Cell<T, Descriptor> &cell);
    Array<T, 3> const &get_u() const
    {
        return u;
    }

public:
    virtual void planeBoundaryCompletion(
        int direction, int orientation, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice);
    virtual void edgeBoundaryCompletion(
        int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice);
    virtual void cornerBoundaryCompletion(
        int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice);

    virtual void planeComputeMacroscopic(
        int direction, int orientation, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar);
    virtual void edgeComputeMacroscopic(
        int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar);
    virtual void cornerComputeMacroscopic(
        int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar);

private:
    Array<T, 3> u;
    static int id;
};

/// Neumann version of the generic NLD boundary class.
template <typename T, template <typename U> class Descriptor>
class NLD_VelocityNeumannBoundaryDynamics3D : public NonLocalBoundaryDynamics3D<T, Descriptor> {
public:
    NLD_VelocityNeumannBoundaryDynamics3D(
        Dynamics<T, Descriptor> *baseDynamics_, bool noPenetration_ = false);
    NLD_VelocityNeumannBoundaryDynamics3D(HierarchicUnserializer &unserializer);
    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);
    virtual NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor> *clone() const;
    virtual int getId() const;
    virtual void prepareCollision(Cell<T, Descriptor> &cell);

public:
    virtual void planeBoundaryCompletion(
        int direction, int orientation, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice);
    virtual void edgeBoundaryCompletion(
        int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice);
    virtual void cornerBoundaryCompletion(
        int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice);

    virtual void planeComputeMacroscopic(
        int direction, int orientation, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar);
    virtual void edgeComputeMacroscopic(
        int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar);
    virtual void cornerComputeMacroscopic(
        int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar);

private:
    bool isSameTypeAsI(plint iX, plint iY, plint iZ, BlockLattice3D<T, Descriptor> &lattice);

private:
    bool noPenetration;
    static int id;
};

/// NLD boundary class for outflows: fixed density, extrapolate unknown populations.
template <typename T, template <typename U> class Descriptor>
class NLD_OutflowBoundaryDynamics3D : public NonLocalBoundaryDynamics3D<T, Descriptor> {
public:
    NLD_OutflowBoundaryDynamics3D(Dynamics<T, Descriptor> *baseDynamics_);
    NLD_OutflowBoundaryDynamics3D(HierarchicUnserializer &unserializer);
    virtual NLD_OutflowBoundaryDynamics3D<T, Descriptor> *clone() const;
    virtual int getId() const;
    virtual void prepareCollision(Cell<T, Descriptor> &cell);

public:
    virtual void planeBoundaryCompletion(
        int direction, int orientation, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice);
    virtual void edgeBoundaryCompletion(
        int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice);
    virtual void cornerBoundaryCompletion(
        int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice);

    virtual void planeComputeMacroscopic(
        int direction, int orientation, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar);
    virtual void edgeComputeMacroscopic(
        int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar);
    virtual void cornerComputeMacroscopic(
        int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar);

private:
    static int id;
};

/// Dirichlet version of the generic NLD boundary class for advection-diffusion.
template <typename T, template <typename U> class Descriptor>
class NLD_AD_DirichletDynamics3D : public NonLocalBoundaryDynamics3D<T, Descriptor> {
public:
    NLD_AD_DirichletDynamics3D(Dynamics<T, Descriptor> *baseDynamics_);
    NLD_AD_DirichletDynamics3D(HierarchicUnserializer &unserializer);
    void defineDensity(Cell<T, Descriptor> &cell, T density);
    T computeDensity(Cell<T, Descriptor> const &cell) const;
    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);
    virtual NLD_AD_DirichletDynamics3D<T, Descriptor> *clone() const;
    virtual int getId() const;
    virtual void prepareCollision(Cell<T, Descriptor> &cell);

public:
    virtual void planeBoundaryCompletion(
        int direction, int orientation, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice);
    virtual void edgeBoundaryCompletion(
        int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice);
    virtual void cornerBoundaryCompletion(
        int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice);

    virtual void planeComputeMacroscopic(
        int direction, int orientation, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar);
    virtual void edgeComputeMacroscopic(
        int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar);
    virtual void cornerComputeMacroscopic(
        int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar);

private:
    T rho;
    static int id;
};

/// Neumann version of the generic NLD boundary class.
template <typename T, template <typename U> class Descriptor>
class NLD_AD_NeumannBoundaryDynamics3D : public NonLocalBoundaryDynamics3D<T, Descriptor> {
public:
    NLD_AD_NeumannBoundaryDynamics3D(Dynamics<T, Descriptor> *baseDynamics_);
    NLD_AD_NeumannBoundaryDynamics3D(HierarchicUnserializer &unserializer);
    virtual NLD_AD_NeumannBoundaryDynamics3D<T, Descriptor> *clone() const;
    virtual int getId() const;
    virtual void prepareCollision(Cell<T, Descriptor> &cell);

public:
    virtual void planeBoundaryCompletion(
        int direction, int orientation, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice);
    virtual void edgeBoundaryCompletion(
        int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice);
    virtual void cornerBoundaryCompletion(
        int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice);

    virtual void planeComputeMacroscopic(
        int direction, int orientation, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar);
    virtual void edgeComputeMacroscopic(
        int plane, int normal1, int normal2, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar);
    virtual void cornerComputeMacroscopic(
        int xNormal, int yNormal, int zNormal, plint iX, plint iY, plint iZ,
        BlockLattice3D<T, Descriptor> &lattice, T &rhoBar, Array<T, 3> &j, T &thetaBar);

private:
    bool isSameTypeAsI(plint iX, plint iY, plint iZ, BlockLattice3D<T, Descriptor> &lattice);

private:
    static int id;
};

/// Make a clone of the full composite dynamics, remove all NLD dynamics and whatever stands
///   above it. Then, insert the new NLS dynamics at the top, under the same conditions as
///   for the function cloneAndInsertAtTopDynamics.
template <typename T, template <typename U> class Descriptor>
Dynamics<T, Descriptor> *cloneAndInsertNLDatTop(
    Dynamics<T, Descriptor> const &dynamics, NonLocalBoundaryDynamics3D<T, Descriptor> *newTop);

}  // namespace plb

#endif  // NLD_BOUNDARY_DYNAMICS_3D_H
