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
#ifndef THERMAL_DYNAMICS_H
#define THERMAL_DYNAMICS_H

#include "core/dynamics.h"
#include "core/globalDefs.h"

namespace plb {

/// Common base iso-thermal (or athermal) bulk dynamics
template <typename T, template <typename U> class Descriptor>
class ThermalBulkDynamics : public BasicBulkDynamics<T, Descriptor> {
public:
    ThermalBulkDynamics(T omega_);

    /* *************** Collision, Equilibrium, and Non-equilibrium ******* */

    /// Re-compute particle populations from the leading moments
    virtual void regularize(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
        Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar = T()) const;

    /* *************** Computation of macroscopic variables ************** */

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

    /* *************** Switch between population and moment representation ****** */

    /// Number of variables required to decompose a population representation into moments.
    /** In the present implementation, the decomposition is carried out up to order-1 in the
     *    Chapman-Enskog expansion. Example: Take the D2Q9 lattice. A decomposition means:
     *    - At order 0: Decompose into rho, u, and fNeq (1+2+9=12 variables)
     *    - At order 1: Decompose into rho, u, and PiNeq (1+2+3=6 variables)
     *    - At higher order: Decompose according to order 1.
     */
    virtual plint numDecomposedVariables(plint order) const;

    /// Decompose from population representation into moment representation.
    /**   \sa numDecomposedVariables()
     */
    virtual void decompose(
        Cell<T, Descriptor> const &cell, std::vector<T> &rawData, plint order) const;

    /// Recompose from moment representation to population representation.
    /**   \sa numDecomposedVariables()
     *    This process is also known as "regularization step", and this function is therefore
     *    equivalent to regularize(), although one or the other function may be more useful
     *    in a specific context, due to the form of the parameters.
     */
    virtual void recompose(
        Cell<T, Descriptor> &cell, std::vector<T> const &rawData, plint order) const;

    /// Change the space and time scales of the variables in moment representation.
    /**   \sa numDecomposedVariables()
     *    \param xDxInv Inverse of the factor by which space scale is multiplied.
     *    \param xDt Factor by which time scale is multiplied.
     */
    virtual void rescale(std::vector<T> &rawData, T xDxInv, T xDt, plint order) const;

    /* *************** Additional moments, intended for internal use ***** */

    /// Returns 0, as a default value for isothermal flow.
    virtual T computeEbar(Cell<T, Descriptor> const &cell) const;

private:
    virtual void decomposeOrder0(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void decomposeOrder1(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void recomposeOrder0(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;
    virtual void recomposeOrder1(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;
    virtual void rescaleOrder0(std::vector<T> &rawData, T xDxInv, T xDt) const;
    virtual void rescaleOrder1(std::vector<T> &rawData, T xDxInv, T xDt) const;
};

/// Implementation of O(Ma^2) BGK dynamics
template <typename T, template <typename U> class Descriptor>
class IsoThermalBGKdynamics : public ThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    IsoThermalBGKdynamics(T omega_);
    IsoThermalBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual IsoThermalBGKdynamics<T, Descriptor> *clone() const;

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

/// Implementation of O(Ma^4) BGK dynamics
template <typename T, template <typename U> class Descriptor>
class ThermalBGKdynamics : public ThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    ThermalBGKdynamics(T omega_);
    ThermalBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual ThermalBGKdynamics<T, Descriptor> *clone() const;

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

/// Implementation of O(Ma^2) BGK dynamics
template <typename T, template <typename U> class Descriptor>
class ThermalRLBdynamics : public ThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    ThermalRLBdynamics(T omega_);
    ThermalRLBdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual ThermalRLBdynamics<T, Descriptor> *clone() const;

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

#endif
