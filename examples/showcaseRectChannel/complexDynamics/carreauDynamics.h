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

#ifndef CARREAU_DYNAMICS_H
#define CARREAU_DYNAMICS_H

#include "complexDynamics/variableOmegaDynamics.h"
#include "core/dynamics.h"
#include "core/globalDefs.h"

namespace plb {

/// This class recomputes omega for a generalied newtonian fluid with a carreau constitutive
/// equation.
/** The constitutive equation is nu=nu0*(1+(lambda*|gamma|)^2)^((n-1)/2).
 * Note that in order to be numerically efficient there is no resolution of the implicit equation
 * rather we try to get to the fixed point by auto-replacement in the solution.(??? need to be
 * rewritten)
 */
template <typename T, template <typename U> class Descriptor, int N>
class CarreauDynamics : public OmegaFromPiDynamics<T, Descriptor> {
public:
    CarreauDynamics(Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision = true);
    virtual T getOmegaFromPiAndRhoBar(
        Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T rhoBar) const;
    CarreauDynamics<T, Descriptor, N> *clone() const;
    virtual int getId() const;

private:
    static int id;
};

/// Implementation of O(Ma^2) BGK dynamics with constant average density
/** Semantically, this class is equivalent to RLBdynamics< . , . , BGKdynamics<.,.> >,
 *  but the implementation is more efficient.
 */
template <typename T, template <typename U> class Descriptor, int N>
class BGKCarreauDynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    BGKCarreauDynamics(T omega = (T)1);
    BGKCarreauDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual BGKCarreauDynamics<T, Descriptor, N> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;

private:
    static int id;
};

/// Implementation of O(Ma^2) BGK dynamics with constant average density
/** Semantically, this class is equivalent to RLBdynamics< . , . , BGKdynamics<.,.> >,
 *  but the implementation is more efficient.
 */
template <typename T, template <typename U> class Descriptor, int N>
class RegularizedBGKCarreauDynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    RegularizedBGKCarreauDynamics(T omega = (T)1);
    RegularizedBGKCarreauDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual RegularizedBGKCarreauDynamics<T, Descriptor, N> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;

private:
    static int id;
};

}  // namespace plb

#endif  // VARIABLE_OMEGA_DYNAMICS_H
