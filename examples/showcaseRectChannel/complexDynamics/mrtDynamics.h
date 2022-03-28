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

/* Orestis Malaspinas contributed this code.
 */

/** \file
 * This object is a MRT LB dynamics as described in D.Yu et al. in
 * Progress in Aerospace Sciences 39 (2003) 329-367
 *
 * The full list of implemented MRTdynamics is:
 *
 * In mrtDynamics.h:
 *
 * - MRTdynamics: "The one with rho in front of equilibrium".
 * - IncMRTdynamics: "The one with rho0 in front of equilibrium".
 *
 * In smagorinskyDynamics.h:
 *
 * - SmagorinskyMrtDynamics: "The classical, naive implementation".
 * - ConsistentSmagorinskyMRTdynamics: Consistent according to
 *       Orestis' model: doesn't mix the stress modes with the
 *       other ones when you change the relaxation time.
 * - ConsistentSmagorinskyIncMRTdynamics
 *
 * In externalForceMRTdynamics.h:
 *
 * - GuoExternalForceMRTdynamics
 * - GuoExternalForceSmagorinskyMRTdynamics
 * - GuoExternalForceSmagorinskyIncMRTdynamics
 * - GuoExternalForceConsistentSmagorinskyMRTdynamics
 * - GuoExternalForceConsistentSmagorinskyIncMRTdynamics
 *
 */
#ifndef MRT_DYNAMICS_H
#define MRT_DYNAMICS_H

#include "basicDynamics/isoThermalDynamics.h"
#include "core/globalDefs.h"

namespace plb {

/// Implementation of the MRT collision step
template <typename T, template <typename U> class Descriptor>
class MRTdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    MRTdynamics(T omega_);
    MRTdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual MRTdynamics<T, Descriptor> *clone() const;

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

/// Implementation of incompressible MRT dynamics.
/** This is the MRT equivalent of IncBGKdynamics: the "rho" moment of the
 *  populations appears only as a pressure term in the equilibrium, while
 *  the other terms are multiplied by the constant rho0.
 **/
template <typename T, template <typename U> class Descriptor>
class IncMRTdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    IncMRTdynamics(T omega_);
    IncMRTdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual IncMRTdynamics<T, Descriptor> *clone() const;

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

    virtual bool velIsJ() const;

    /* *************** Macroscopic variables ***************************** */

    /// Velocity is equal to j, not u.
    virtual void computeVelocity(
        Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const;
    /// For PiNeq, subtract equilibrium term jj instead of invRho*jj.
    virtual void computeRhoBarJPiNeq(
        Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const;

private:
    static int id;
};

template <typename T, template <typename U> class Descriptor>
class D3Q13Dynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    D3Q13Dynamics(T omega_);
    D3Q13Dynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual D3Q13Dynamics<T, Descriptor> *clone() const;

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

    virtual bool velIsJ() const;

    /* *************** Macroscopic variables ***************************** */
private:
    static int id;
};

}  // namespace plb

#endif  // MRT_DYNAMICS_H
