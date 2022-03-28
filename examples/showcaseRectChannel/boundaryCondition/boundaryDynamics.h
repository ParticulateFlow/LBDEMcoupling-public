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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- header file.
 */
#ifndef BOUNDARY_DYNAMICS_H
#define BOUNDARY_DYNAMICS_H

#include "core/dynamics.h"
#include "core/globalDefs.h"
#include "core/hierarchicSerializer.h"

namespace plb {

/// Computation of the macroscopic variables is obtained after invoking completion scheme
/** You can instantiate this non-abstract class. This is not interesting, though,
 *  as it is semantically identical with the base dynamics. Classes which inherit
 *  from this one are more interesting, as they customize some aspects. Examples:
 *  computation of the velocity on a velocity boundary, or, completion scheme for
 *  a boundary node with missing particle populations (and no data processor to
 *  complete them by non-local means).
 */
template <typename T, template <typename U> class Descriptor>
class BoundaryCompositeDynamics : public PreparePopulationsDynamics<T, Descriptor> {
public:
    /* *************** Construction and Destruction ********************** */

    BoundaryCompositeDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_ = true);

    /// Clone the object, based on its dynamic type
    virtual BoundaryCompositeDynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    virtual bool isBoundary() const;

    /* *************** Computation of macroscopic variables ************** */

    /// Compute the local particle density in lattice units
    virtual T computeDensity(Cell<T, Descriptor> const &cell) const;

    /// Compute the local pressure in lattice units
    virtual T computePressure(Cell<T, Descriptor> const &cell) const;

    /// Compute the local fluid velocity in lattice units
    virtual void computeVelocity(
        Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &velocity) const;

    /// Compute the temperature in lattice units
    virtual T computeTemperature(Cell<T, Descriptor> const &cell) const;

    /// Compute the "off-equilibrium part of Pi"
    virtual void computePiNeq(
        Cell<T, Descriptor> const &cell, Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const;

    /// Compute the deviatoric stress tensor
    virtual void computeShearStress(
        Cell<T, Descriptor> const &cell, Array<T, SymmetricTensor<T, Descriptor>::n> &stress) const;

    /// Compute the heat flux in lattice units
    virtual void computeHeatFlux(
        Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &q) const;

    /// Compute additional user-defined moments
    virtual void computeMoment(Cell<T, Descriptor> const &cell, plint momentId, T *moment) const;

    /* *************** Additional moments, intended for internal use ***** */

    /// Compute order-0 moment rho-bar
    virtual T computeRhoBar(Cell<T, Descriptor> const &cell) const;

    /// Compute order-0 moment rho-bar and order-1 moment j
    virtual void computeRhoBarJ(
        Cell<T, Descriptor> const &cell, T &rhoBar_, Array<T, Descriptor<T>::d> &j) const;

    /// Compute order-0 moment rho-bar, order-1 moment j, and order-2
    ///   off-equilibrium moment PiNeq.
    virtual void computeRhoBarJPiNeq(
        Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const;

    /// Compute e-bar, which is related to the internal energy
    virtual T computeEbar(Cell<T, Descriptor> const &cell) const;

    /* *************** Default completion scheme ************************* */

    /// Default completion scheme, does nothing
    virtual void completePopulations(Cell<T, Descriptor> &cell) const;

    /// Decompose from population representation into moment representation.
    virtual void decompose(
        Cell<T, Descriptor> const &cell, std::vector<T> &rawData, plint order) const;

private:
    static int id;
};

/// Value of the density is stored inside Dynamics
template <typename T, template <typename U> class Descriptor>
class StoreDensityDynamics : public BoundaryCompositeDynamics<T, Descriptor> {
public:
    /* *************** Construction and Destruction ********************** */

    StoreDensityDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_ = true);
    StoreDensityDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object, based on its dynamic type
    virtual StoreDensityDynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);

    /* *************** Computation of macroscopic variables ************** */

    /// Compute the local particle density in lattice units
    virtual T computeDensity(Cell<T, Descriptor> const &cell) const;

    /// Define density. Stores the value inside the Dynamics object.
    virtual void defineDensity(Cell<T, Descriptor> &cell, T rho_);

    /* *************** Additional moments, intended for internal use ***** */

    /// Compute order-0 moment rho-bar
    virtual T computeRhoBar(Cell<T, Descriptor> const &cell) const;

    /// Compute order-0 moment rho-bar and order-1 moment j
    virtual void computeRhoBarJ(
        Cell<T, Descriptor> const &cell, T &rhoBar_, Array<T, Descriptor<T>::d> &j) const;

protected:
    T rhoBar;

private:
    static int id;
};

/// Value of the velocity is stored inside dynamics
template <typename T, template <typename U> class Descriptor>
class StoreVelocityDynamics : public BoundaryCompositeDynamics<T, Descriptor> {
public:
    /* *************** Construction and Destruction ********************** */

    StoreVelocityDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_ = true);
    StoreVelocityDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object, based on its dynamic type
    virtual StoreVelocityDynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);

    /* *************** Computation of macroscopic variables ************** */

    /// Compute the local fluid velocity in lattice units
    virtual void computeVelocity(
        Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &velocity_) const;

    /// Define velocity. Stores value inside Dynamics object.
    virtual void defineVelocity(
        Cell<T, Descriptor> &cell, Array<T, Descriptor<T>::d> const &velocity_);

    /* *************** Additional moments, intended for internal use ***** */

    /// Compute order-0 moment rho-bar
    virtual T computeRhoBar(Cell<T, Descriptor> const &cell) const;

    /// Compute order-0 moment rho-bar and order-1 moment j
    virtual void computeRhoBarJ(
        Cell<T, Descriptor> const &cell, T &rhoBar_, Array<T, Descriptor<T>::d> &j) const;

protected:
    Array<T, Descriptor<T>::d> velocity;

private:
    static int id;
};

/// Density and Velocity are stored inside dynamics
template <typename T, template <typename U> class Descriptor>
class StoreDensityAndVelocityDynamics : public BoundaryCompositeDynamics<T, Descriptor> {
public:
    /* *************** Construction and Destruction ********************** */

    StoreDensityAndVelocityDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_ = true);
    StoreDensityAndVelocityDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object, based on its dynamic type
    virtual StoreDensityAndVelocityDynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);

    /* *************** Computation of macroscopic variables ************** */

    /// Compute the local particle density in lattice units
    virtual T computeDensity(Cell<T, Descriptor> const &cell) const;

    /// Compute the local fluid velocity in lattice units
    virtual void computeVelocity(
        Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &velocity_) const;

    /// Define density. Stores the value inside the Dynamics object.
    virtual void defineDensity(Cell<T, Descriptor> &cell, T rho_);

    /// Define velocity. Stores value inside Dynamics object.
    virtual void defineVelocity(
        Cell<T, Descriptor> &cell, Array<T, Descriptor<T>::d> const &velocity_);

    /* *************** Additional moments, intended for internal use ***** */

    /// Compute order-0 moment rho-bar
    virtual T computeRhoBar(Cell<T, Descriptor> const &cell) const;

    /// Compute order-0 moment rho-bar and order-1 moment j
    virtual void computeRhoBarJ(
        Cell<T, Descriptor> const &cell, T &rhoBar_, Array<T, Descriptor<T>::d> &j) const;

protected:
    T rhoBar;
    Array<T, Descriptor<T>::d> velocity;

private:
    static int id;
};

/// Temperature and Velocity are stored inside dynamics
template <typename T, template <typename U> class Descriptor>
class StoreTemperatureAndVelocityDynamics : public BoundaryCompositeDynamics<T, Descriptor> {
public:
    /* *************** Construction and Destruction ********************** */

    StoreTemperatureAndVelocityDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_ = true);
    StoreTemperatureAndVelocityDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object, based on its dynamic type
    virtual StoreTemperatureAndVelocityDynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);

    /* *************** Computation of macroscopic variables ************** */

    /// Compute the local fluid velocity in lattice units
    virtual void computeVelocity(
        Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &velocity_) const;

    /// Compute the temperature in lattice units
    virtual T computeTemperature(Cell<T, Descriptor> const &cell) const;

    /// Define velocity. Stores value inside Dynamics object.
    virtual void defineVelocity(
        Cell<T, Descriptor> &cell, Array<T, Descriptor<T>::d> const &velocity_);

    /// Define density. Stores the value inside the Dynamics object.
    virtual void defineTemperature(Cell<T, Descriptor> &cell, T theta_);

    /* *************** Additional moments, intended for internal use ***** */

    /// Compute order-0 moment rho-bar
    virtual T computeRhoBar(Cell<T, Descriptor> const &cell) const;

    /// Compute order-0 moment rho-bar and order-1 moment j
    virtual void computeRhoBarJ(
        Cell<T, Descriptor> const &cell, T &rhoBar_, Array<T, Descriptor<T>::d> &j) const;

protected:
    T thetaBar;
    Array<T, Descriptor<T>::d> velocity;

private:
    static int id;
};

/// Velocity Dirichlet boundary dynamics for a straight wall
template <typename T, template <typename U> class Descriptor, int direction, int orientation>
class VelocityDirichletBoundaryDynamics : public StoreVelocityDynamics<T, Descriptor> {
public:
    VelocityDirichletBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_ = true);
    VelocityDirichletBoundaryDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object, based on its dynamic type
    virtual VelocityDirichletBoundaryDynamics<T, Descriptor, direction, orientation> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Compute density from incoming particle populations
    virtual T computeDensity(Cell<T, Descriptor> const &cell) const;
    /// Compute order-0 moment rho-bar
    virtual T computeRhoBar(Cell<T, Descriptor> const &cell) const;
    /// Compute order-0 moment rho-bar and order-1 moment j
    virtual void computeRhoBarJ(
        Cell<T, Descriptor> const &cell, T &rhoBar_, Array<T, Descriptor<T>::d> &j) const;

private:
    static int id;
};

/// Velocity Dirichlet boundary dynamics for a straight wall
template <typename T, template <typename U> class Descriptor, int direction, int orientation>
class VelocityDirichletConstRhoBoundaryDynamics : public StoreVelocityDynamics<T, Descriptor> {
public:
    VelocityDirichletConstRhoBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_ = true);
    VelocityDirichletConstRhoBoundaryDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object, based on its dynamic type
    virtual VelocityDirichletConstRhoBoundaryDynamics<T, Descriptor, direction, orientation>
        *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Compute density from incoming particle populations
    virtual T computeDensity(Cell<T, Descriptor> const &cell) const;
    /// Compute order-0 moment rho-bar
    virtual T computeRhoBar(Cell<T, Descriptor> const &cell) const;
    /// Compute order-0 moment rho-bar and order-1 moment j
    virtual void computeRhoBarJ(
        Cell<T, Descriptor> const &cell, T &rhoBar_, Array<T, Descriptor<T>::d> &j) const;

private:
    static int id;
};

/// Density Dirichlet boundary dynamics for a straight wall
template <typename T, template <typename U> class Descriptor, int direction, int orientation>
class DensityDirichletBoundaryDynamics : public StoreDensityDynamics<T, Descriptor> {
public:
    DensityDirichletBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_ = true);
    DensityDirichletBoundaryDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object, based on its dynamic type
    virtual DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Compute the local fluid velocity in lattice units
    virtual void computeVelocity(
        Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &velocity) const;
    /// Compute order-0 moment rho-bar and order-1 moment j
    virtual void computeRhoBarJ(
        Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j) const;

public:
    void computeJ(Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &j_) const;

private:
    static int id;
};

}  // namespace plb

#endif  // BOUNDARY_DYNAMICS_H
