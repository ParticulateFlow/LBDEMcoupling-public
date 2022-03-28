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
 * Collision terms with external force -- header file.
 */
#ifndef EXTERNAL_FORCE_DYNAMICS_H
#define EXTERNAL_FORCE_DYNAMICS_H

#include "basicDynamics/isoThermalDynamics.h"
#include "core/dynamics.h"
#include "core/globalDefs.h"

namespace plb {

/** Implementation of the computation of the velocity in the
 * presence of an external force for compressible models
 */
template <typename T, template <typename U> class Descriptor>
class ExternalForceDynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    ExternalForceDynamics(T omega_);

    /* *************** Velocity computation ************************* */

    /// Implementation of velocity computation
    virtual void computeVelocity(
        Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const;

    virtual void computeVelocityExternal(
        Cell<T, Descriptor> const &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j,
        Array<T, Descriptor<T>::d> &u) const;
};

/** Implementation of the computation of the velocity in the
 * presence of an external force for incompressible models
 */
template <typename T, template <typename U> class Descriptor>
class IncExternalForceDynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    IncExternalForceDynamics(T omega_);

    /* *************** Velocity computation ************************* */

    /// Implementation of velocity computation
    virtual void computeVelocity(
        Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const;

    virtual void computeVelocityExternal(
        Cell<T, Descriptor> const &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j,
        Array<T, Descriptor<T>::d> &u) const;

    /// For PiNeq, subtract equilibrium term jj instead of invRho*jj.
    virtual void computeRhoBarJPiNeq(
        Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const;
};

/// Implementation BGK dynamics with a simple, linear external force
template <typename T, template <typename U> class Descriptor>
class NaiveExternalForceBGKdynamics : public ExternalForceDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    NaiveExternalForceBGKdynamics(T omega_);
    NaiveExternalForceBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual NaiveExternalForceBGKdynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

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

private:
    static int id;
};

/// Implementation BGK dynamics with a simple, linear external force
template <typename T, template <typename U> class Descriptor>
class NaiveExternalForcePrecondBGKdynamics : public ExternalForceDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    NaiveExternalForcePrecondBGKdynamics(T omega_, T invGamma_);
    NaiveExternalForcePrecondBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual NaiveExternalForcePrecondBGKdynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

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
    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);

private:
    static int id;
    T invGamma;
};

/// Implementation of O(Ma^2) BGK dynamics with an external force (Guo approach)
template <typename T, template <typename U> class Descriptor>
class GuoExternalForceBGKdynamics : public ExternalForceDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    GuoExternalForceBGKdynamics(T omega_);
    GuoExternalForceBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual GuoExternalForceBGKdynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

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

private:
    static int id;
};

/// Implementation of O(Ma^2) BGK dynamics with an external force (Guo approach)
template <typename T, template <typename U> class Descriptor>
class GuoExternalForceCompleteRegularizedBGKdynamics : public ExternalForceDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    GuoExternalForceCompleteRegularizedBGKdynamics(T omega_);
    GuoExternalForceCompleteRegularizedBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual GuoExternalForceCompleteRegularizedBGKdynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

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

    virtual void computeEquilibria(
        Array<T, Descriptor<T>::q> &fEq, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
        T thetaBar = T()) const;

    virtual void regularize(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
        Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar = T()) const;

    virtual void decomposeOrder0(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;

    virtual void recomposeOrder0(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;

private:
    static int id;
};

/// Implementation of O(Ma^2) BGK dynamics with an external force (Guo approach)
template <typename T, template <typename U> class Descriptor>
class GuoExternalForceConsistentSmagorinskyCompleteRegularizedBGKdynamics :
    public GuoExternalForceCompleteRegularizedBGKdynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    GuoExternalForceConsistentSmagorinskyCompleteRegularizedBGKdynamics(T omega_, T cSmago_);
    GuoExternalForceConsistentSmagorinskyCompleteRegularizedBGKdynamics(
        HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual GuoExternalForceConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor>
        *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    virtual void serialize(HierarchicSerializer &serializer) const;

    virtual void unserialize(HierarchicUnserializer &unserializer);

    virtual T getDynamicParameter(plint whichParameter, Cell<T, Descriptor> const &cell) const;

    /* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);

    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);

private:
    static int id;
    T cSmago;
};

/// Implementation of O(Ma^2) BGK dynamics with an external force (Shan/Chen approach)
template <typename T, template <typename U> class Descriptor>
class ShanChenExternalForceBGKdynamics : public ExternalForceDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    ShanChenExternalForceBGKdynamics(T omega_);
    ShanChenExternalForceBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual ShanChenExternalForceBGKdynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

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

private:
    static int id;
};

/// Implementation of O(Ma^2) BGK dynamics with an external force (He et al. approach)
template <typename T, template <typename U> class Descriptor>
class HeExternalForceBGKdynamics : public ExternalForceDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    HeExternalForceBGKdynamics(T omega_);
    HeExternalForceBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual HeExternalForceBGKdynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

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

private:
    static int id;
};

/// Incompressible version of the Guo external forcing
template <typename T, template <typename U> class Descriptor>
class IncGuoExternalForceBGKdynamics : public IncExternalForceDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    IncGuoExternalForceBGKdynamics(T omega_);
    IncGuoExternalForceBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual IncGuoExternalForceBGKdynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

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

    ///  Is it in the incompressible BGK model
    virtual bool velIsJ() const;

private:
    static int id;
};

/// Implementation of O(Ma^2) BGK dynamics with an external force (Guo approach)
template <typename T, template <typename U> class Descriptor>
class ShanChenExternalForceRegularizedBGKdynamics : public ExternalForceDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    ShanChenExternalForceRegularizedBGKdynamics(T omega_);
    ShanChenExternalForceRegularizedBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual ShanChenExternalForceRegularizedBGKdynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

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

private:
    static int id;
};

}  // namespace plb

#endif  // EXTERNAL_FORCE_DYNAMICS_H
