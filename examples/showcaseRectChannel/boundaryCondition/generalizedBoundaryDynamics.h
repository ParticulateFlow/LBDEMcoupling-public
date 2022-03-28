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
#ifndef GENERALIZED_BOUNDARY_DYNAMICS_H
#define GENERALIZED_BOUNDARY_DYNAMICS_H

#include "boundaryCondition/boundaryDynamics.h"
#include "core/globalDefs.h"

namespace plb {

/// Generic velocity boundary dynamics for a straight wall
template <typename T, template <typename U> class Descriptor>
class GeneralizedVelocityBoundaryDynamics : public StoreVelocityDynamics<T, Descriptor> {
public:
    GeneralizedVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, std::vector<plint> missingIndices_,
        bool automaticPrepareCollision = true);
    GeneralizedVelocityBoundaryDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object, based on its dynamic type
    virtual GeneralizedVelocityBoundaryDynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);

    void computeUlb(Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &uLb) const;

    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T, Descriptor> &cell) const;

private:
    static int id;
    std::vector<plint> missingIndices, knownIndices;
};

/// Mass Conserving Generic velocity boundary dynamics for a straight wall
template <typename T, template <typename U> class Descriptor>
class GeneralizedMassConservingVelocityBoundaryDynamics :
    public StoreVelocityDynamics<T, Descriptor> {
public:
    GeneralizedMassConservingVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, std::vector<plint> missingIndices_,
        std::vector<plint> knownIndices_, std::vector<plint> inGoingIndices_,
        bool automaticPrepareCollision = true);
    GeneralizedMassConservingVelocityBoundaryDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object, based on its dynamic type
    virtual GeneralizedMassConservingVelocityBoundaryDynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);

    void computeUlb(Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &uLb) const;

    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T, Descriptor> &cell) const;

private:
    static int id;
    std::vector<plint> missingIndices, knownIndices, inGoingIndices;
};

/// Generic density Dirichlet boundary dynamics for a straight wall
template <typename T, template <typename U> class Descriptor, int direction, int orientation>
class GeneralizedDensityBoundaryDynamics :
    public DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation> {
public:
    GeneralizedDensityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, std::vector<plint> missingIndices_,
        bool automaticPrepareCollision = true);
    GeneralizedDensityBoundaryDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object, based on its dynamic type
    virtual GeneralizedDensityBoundaryDynamics<T, Descriptor, direction, orientation> *clone()
        const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);

    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T, Descriptor> &cell) const;

private:
    static int id;
    std::vector<plint> missingIndices, knownIndices;
};

// ========================================================================== //
// ============= Generalized Temperature and Velocity imposed BC ============ //
// ========== On wall node (in opposition to bulk nodes with unknowns======== //
// ========================================================================== //

/// Generic velocity and temperature boundary dynamics for a straight wall
template <typename T, template <typename U> class Descriptor>
class GeneralizedVelocityTemperatureBoundaryDynamics :
    public StoreTemperatureAndVelocityDynamics<T, Descriptor> {
public:
    GeneralizedVelocityTemperatureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, std::vector<plint> missingIndices_,
        bool massConserving, bool automaticPrepareCollision = true);
    GeneralizedVelocityTemperatureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, std::vector<plint> missingIndices_,
        std::vector<plint> knownIndices_, bool massConserving,
        bool automaticPrepareCollision = true);

    GeneralizedVelocityTemperatureBoundaryDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object, based on its dynamic type
    virtual GeneralizedVelocityTemperatureBoundaryDynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);

    void computeUlb(Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &uLb) const;

    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T, Descriptor> &cell) const;

private:
    static int id;
    std::vector<plint> missingIndices, knownIndices;
    bool massConserving;
};

// ========================================================================== //
//  Generalized bulk node (nodes where there are unknown f_is close to walls  //
// ========================================================================== //

/// Generic velocity and temperature boundary dynamics for a straight wall
template <typename T, template <typename U> class Descriptor>
class GeneralizedNextToBoundaryDynamics : public BoundaryCompositeDynamics<T, Descriptor> {
public:
    GeneralizedNextToBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, std::vector<plint> missingIndices_,
        bool automaticPrepareCollision = true);

    GeneralizedNextToBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, std::vector<plint> missingIndices_,
        std::vector<plint> knownIndices_, bool automaticPrepareCollision = true);

    GeneralizedNextToBoundaryDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object, based on its dynamic type
    virtual GeneralizedNextToBoundaryDynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);

    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T, Descriptor> &cell);

private:
    static int id;
    std::vector<plint> missingIndices, knownIndices;
    T rho, thetaBar;
    Array<T, Descriptor<T>::d> u;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> Qneq;
};

}  // namespace plb

#endif  // GENERALIZED_BOUNDARY_DYNAMICS_H
