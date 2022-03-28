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
#ifndef WAVE_PROPAGATION_H
#define WAVE_PROPAGATION_H

#include "core/dynamics.h"
#include "core/globalDefs.h"

namespace plb {

/// Implementation of O(Ma^2) BGK dynamics with adjustable speed of sound
template <typename T, template <typename U> class Descriptor>
class WaveDynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    WaveDynamics(T vs2_);
    WaveDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual WaveDynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;

    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);

    /* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);

    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;

    /* *************** Configurable parameters *************************** */

    /// Set local value of any generic parameter
    virtual void setParameter(plint whichParameter, T value);
    /// Get local value of any generic parameter
    virtual T getParameter(plint whichParameter) const;
    /// Set local speed of sound
    void setVs2(T vs2_);
    /// Get local speed of sound
    T getVs2() const;

private:
    /* *************** Static implementation methods********************** */

    /// Implementation of collision operator
    static T waveCollision(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T vs2);
    /// Implementation of equilibrium
    static T waveEquilibrium(
        plint iPop, T rhoBar, T invRho, Array<T, Descriptor<T>::d> const &j, T jSqr, T vs2);

private:
    T vs2;  ///< speed of sound
private:
    static int id;
};

/// This class implements the absorbing condition of H. Xu
template <typename T, template <typename U> class Descriptor>
class WaveAbsorptionDynamics : public CompositeDynamics<T, Descriptor> {
public:
    WaveAbsorptionDynamics(Dynamics<T, Descriptor> *baseDynamics_);
    WaveAbsorptionDynamics(HierarchicUnserializer &unserializer);
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);
    virtual WaveAbsorptionDynamics<T, Descriptor> *clone() const
    {
        return new WaveAbsorptionDynamics<T, Descriptor>(*this);
    }
    virtual void prepareCollision(Cell<T, Descriptor> &cell);
    /// Return a unique ID for this class.

    /// Recompose from moment representation to population representation.
    virtual void recompose(
        Cell<T, Descriptor> &cell, std::vector<T> const &rawData, plint order) const;
    virtual void recomposeOrder0(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;
    virtual void recomposeOrder1(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;

    virtual int getId() const;

private:
    static int id;
};

// Declaration of a specific "sigma" function for WaveAbsorptionDynamics.

template <typename T>
class WaveAbsorptionSigmaFunction3D {
public:
    WaveAbsorptionSigmaFunction3D(Box3D domain_, Array<plint, 6> const &numCells_, T omega_);
    T operator()(plint iX, plint iY, plint iZ) const;

private:
    void addDistance(plint from, plint pos, std::vector<plint> &distances, plint i) const;
    T sigma(T x0, T x1, T x) const;

private:
    Box3D domain;
    Array<plint, 6> numCells;
    T xi;
};

// Declaration of a specific 2D "sigma" function for WaveAbsorptionDynamics

template <typename T>
class WaveAbsorptionSigmaFunction2D {
public:
    WaveAbsorptionSigmaFunction2D(Box2D domain_, Array<plint, 4> const &numCells_, T omega_);
    T operator()(plint iX, plint iY) const;

private:
    void addDistance(plint from, plint pos, std::vector<plint> &distances, plint i) const;
    T sigma(T x0, T x1, T x) const;

private:
    Box2D domain;
    Array<plint, 4> numCells;
    T xi;
};

}  // namespace plb

#endif  // WAVE_PROPAGATION_H
