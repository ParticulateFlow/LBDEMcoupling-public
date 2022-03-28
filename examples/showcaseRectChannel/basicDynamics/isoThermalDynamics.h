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
#ifndef ISO_THERMAL_DYNAMICS_H
#define ISO_THERMAL_DYNAMICS_H

#include "core/dynamics.h"
#include "core/globalDefs.h"

namespace plb {

/// Common base iso-thermal (or athermal) bulk dynamics
template <typename T, template <typename U> class Descriptor>
class IsoThermalBulkDynamics : public BasicBulkDynamics<T, Descriptor> {
public:
    IsoThermalBulkDynamics(T omega_);

    /* *************** Collision, Equilibrium, and Non-equilibrium ******* */

    /// Re-compute particle populations from the leading moments
    virtual void regularize(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
        Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar = T()) const;

    /* *************** Computation of macroscopic variables ************** */

    /// Returns 1, as a default value for isothermal flow.
    virtual T computeTemperature(Cell<T, Descriptor> const &cell) const;

    /// Compute the "off-equilibrium part of Pi"
    virtual void computePiNeq(
        Cell<T, Descriptor> const &cell, Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const;

    /// Compute the shear stress tensor
    virtual void computeShearStress(
        Cell<T, Descriptor> const &cell, Array<T, SymmetricTensor<T, Descriptor>::n> &stress) const;

    /// Returns 0, as a default value for isothermal flow.
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
    virtual void rescale(int dxScale, int dtScale)
    {
        Dynamics<T, Descriptor>::rescale(dxScale, dtScale);
    }

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
class BGKdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    BGKdynamics(T omega_);
    BGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual BGKdynamics<T, Descriptor> *clone() const;

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
    virtual void decomposeOrder0(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void recomposeOrder0(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;

private:
    static int id;
};

/// Implementation of O(Ma^2) BGK dynamics
template <typename T, template <typename U> class Descriptor>
class StochasticDynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    StochasticDynamics(T omega_);
    StochasticDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual StochasticDynamics<T, Descriptor> *clone() const;

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
    virtual void decomposeOrder0(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void recomposeOrder0(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;

private:
    static int id;
};

/// Implementation of O(Ma^2) BGK dynamics
template <typename T, template <typename U> class Descriptor>
class CompleteBGKdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    CompleteBGKdynamics(T omega_);
    CompleteBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual CompleteBGKdynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /* *************** Collision and Equilibrium and Regularize************************* */
    /// Re-compute particle populations from the leading moments
    virtual void regularize(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
        Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar = T()) const;

    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);

    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);

    virtual void computeEquilibria(
        Array<T, Descriptor<T>::q> &fEq, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
        T thetaBar = T()) const;

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;

private:
    virtual void decomposeOrder0(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void recomposeOrder0(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;

private:
    static int id;
};

/// Implementation of O(Ma^2) TRT (Orestis Malaspinas' style and not Irina Ginzburg's) dynamics
template <typename T, template <typename U> class Descriptor>
class CompleteRegularizedBGKdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    CompleteRegularizedBGKdynamics(T omega_);
    CompleteRegularizedBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual CompleteRegularizedBGKdynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;
    /* *************** Collision and Equilibrium and Regularize************************* */
    /// Re-compute particle populations from the leading moments
    virtual void regularize(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
        Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar = T()) const;

    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);

    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);

    virtual void computeEquilibria(
        Array<T, Descriptor<T>::q> &fEq, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
        T thetaBar = T()) const;

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;

private:
    virtual void decomposeOrder0(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void recomposeOrder0(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;

private:
    static int id;
};

/// Implementation of O(Ma^2) TRT (Orestis Malaspinas' style and not Irina Ginzburg's) dynamics
template <typename T, template <typename U> class Descriptor>
class CompleteTRTdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    CompleteTRTdynamics(T omega_, T psi_, int order_ = 6);
    CompleteTRTdynamics(T omega_, int order_ = 6);
    CompleteTRTdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual CompleteTRTdynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;

    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);

    /* *************** Collision and Equilibrium and Regularize************************* */
    /// Re-compute particle populations from the leading moments
    virtual void regularize(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
        Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar = T()) const;

    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);

    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);

    virtual void computeEquilibria(
        Array<T, Descriptor<T>::q> &fEq, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
        T thetaBar = T()) const;

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;

    /* *************** Configurable parameters *************************** */

    /// Get all relaxation frequencies (in case of an SRT model they are all equal to omega).
    virtual Array<T, Descriptor<T>::q> getRelaxationFrequencies() const;
    /// Set all relaxation frequencies (in case of an SRT model they should all equal omega).
    virtual void setRelaxationFrequencies(Array<T, Descriptor<T>::q> const &frequencies);
    /// Set local value of any generic parameter
    virtual void setParameter(plint whichParameter, T value);
    /// Get local value of any generic parameter
    virtual T getParameter(plint whichParameter) const;
    /// Set second relaxation time
    void setPsi(T psi_);
    /// Get second relaxation time
    T getPsi() const;

private:
    template <typename U>
    CompleteTRTdynamics(T omega_, T psi_, U order_);

    template <typename U>
    CompleteTRTdynamics(T omega_, U order_);

private:
    virtual void decomposeOrder0(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void recomposeOrder0(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;

private:
    static int id;
    T psi;
    int order;
};

/// Implementation of O(Ma^2) TRT (Orestis Malaspinas' style and not Irina Ginzburg's) dynamics
template <typename T, template <typename U> class Descriptor>
class CompleteRegularizedTRTdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    CompleteRegularizedTRTdynamics(T omega_, T psi_, int order_ = 6);
    CompleteRegularizedTRTdynamics(T omega_, int order_ = 6);
    CompleteRegularizedTRTdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual CompleteRegularizedTRTdynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;

    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);

    /* *************** Collision and Equilibrium and Regularize************************* */
    /// Re-compute particle populations from the leading moments
    virtual void regularize(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
        Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar = T()) const;

    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);

    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);

    virtual void computeEquilibria(
        Array<T, Descriptor<T>::q> &fEq, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
        T thetaBar = T()) const;

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;

    /* *************** Configurable parameters *************************** */

    /// Get all relaxation frequencies (in case of an SRT model they are all equal to omega).
    virtual Array<T, Descriptor<T>::q> getRelaxationFrequencies() const;
    /// Set all relaxation frequencies (in case of an SRT model they should all equal omega).
    virtual void setRelaxationFrequencies(Array<T, Descriptor<T>::q> const &frequencies);
    /// Set local value of any generic parameter
    virtual void setParameter(plint whichParameter, T value);
    /// Get local value of any generic parameter
    virtual T getParameter(plint whichParameter) const;
    /// Set second relaxation time
    void setPsi(T psi_);
    /// Get second relaxation time
    T getPsi() const;

private:
    template <typename U>
    CompleteRegularizedTRTdynamics(T omega_, T psi_, U order_);

    template <typename U>
    CompleteRegularizedTRTdynamics(T omega_, U order_);

private:
    virtual void decomposeOrder0(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void recomposeOrder0(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;

private:
    static int id;
    T psi;
    int order;
};

/// Implementation of O(Ma^2) TRT (Orestis Malaspinas' style and not Irina Ginzburg's) dynamics
/// with only second order equilibrium terms.
template <typename T, template <typename U> class Descriptor>
class TruncatedTRTdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    TruncatedTRTdynamics(T omega_, T psi_);

    TruncatedTRTdynamics(T omega_);
    TruncatedTRTdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual TruncatedTRTdynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;

    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);

    /* *************** Collision and Equilibrium and Regularize************************* */
    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);

    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);

    virtual void computeEquilibria(
        Array<T, Descriptor<T>::q> &fEq, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
        T thetaBar = T()) const;

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;

    /* *************** Configurable parameters *************************** */

    /// Get all relaxation frequencies (in case of an SRT model they are all equal to omega).
    virtual Array<T, Descriptor<T>::q> getRelaxationFrequencies() const;
    /// Set all relaxation frequencies (in case of an SRT model they should all equal omega).
    virtual void setRelaxationFrequencies(Array<T, Descriptor<T>::q> const &frequencies);
    /// Set local value of any generic parameter
    virtual void setParameter(plint whichParameter, T value);
    /// Get local value of any generic parameter
    virtual T getParameter(plint whichParameter) const;
    /// Set second relaxation time
    void setPsi(T psi_);
    /// Get second relaxation time
    T getPsi() const;

private:
    virtual void decomposeOrder0(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void recomposeOrder0(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;

private:
    static int id;
    T psi;
};

template <typename T, template <typename U> class Descriptor>
class StoreRhoBarJBGKdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    StoreRhoBarJBGKdynamics(T omega_);
    StoreRhoBarJBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual StoreRhoBarJBGKdynamics<T, Descriptor> *clone() const;

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

/// Implementation of O(Ma^2) BGK dynamics, density and momentum taken from external scalars
template <typename T, template <typename U> class Descriptor>
class ExternalMomentBGKdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    ExternalMomentBGKdynamics(T omega_);
    ExternalMomentBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual ExternalMomentBGKdynamics<T, Descriptor> *clone() const;

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

    /* *************** Moments ******************************************* */

    // The function computeDensity() is not overridden, and the default
    //   implementation is kept, for two reasons. First, it is equivalent to access
    //   rho from the external scalar or to recompute it from the bulk (this
    //   is not the case for the velocity, in a Shan/Chen multicomponent model).
    //   Second, the Shan/Chen data-processor needs computeDensity() to be
    //   default implemented, because it uses this function to treat walls with
    //   a virtual-density mechanism.

    /** Accesses velocity from external scalar. **/
    virtual void computeVelocity(
        Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const;
    /// Compute order-0 moment rho-bar
    /** Accesses rhoBar from external scalar. **/
    virtual T computeRhoBar(Cell<T, Descriptor> const &cell) const;

    /// Compute order-0 moment rho-bar and order-1 moment j
    /** Accesses rhoBar and j from external scalar. **/
    virtual void computeRhoBarJ(
        Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j) const;

private:
    static int id;
};

/// Implementation of O(Ma^2) BGK dynamics, velocity taken from external scalar
template <typename T, template <typename U> class Descriptor>
class ExternalVelocityBGKdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    ExternalVelocityBGKdynamics(T omega_);
    ExternalVelocityBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual ExternalVelocityBGKdynamics<T, Descriptor> *clone() const;

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
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const;

    /* *************** Moments ******************************************* */

    // The function computeDensity() is not overridden, and the default
    //   implementation is kept, for two reasons. First, it is equivalent to access
    //   rho from the external scalar or to recompute it from the bulk (this
    //   is not the case for the velocity, in a Shan/Chen multicomponent model).
    //   Second, the Shan/Chen data-processor needs computeDensity() to be
    //   default implemented, because it uses this function to treat walls with
    //   a virtual-density mechanism.

    /** Accesses velocity from external scalar. **/
    virtual void computeVelocity(
        Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const;
    /// Compute order-0 moment rho-bar
    /** Accesses rhoBar from external scalar. **/
    virtual T computeRhoBar(Cell<T, Descriptor> const &cell) const;

    /// Compute order-0 moment rho-bar and order-1 moment j
    /** Accesses rhoBar and j from external scalar. **/
    virtual void computeRhoBarJ(
        Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j) const;

private:
    static int id;
};

/// Implementation of quasi-incompressible BGK dynamics.
/** This is exactly the same as BGK dynamics, but computeVelocity() returns
 *  the variable j instead of the variable u.
 **/
template <typename T, template <typename U> class Descriptor>
class QuasiIncBGKdynamics : public BGKdynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    QuasiIncBGKdynamics(T omega_);
    QuasiIncBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual QuasiIncBGKdynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Say if velocity in this dynamics is computed as "j" (the order-1 moment
    ///   of the populations) or as "j/rho".
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

/// Implementation of incompressible BGK dynamics
template <typename T, template <typename U> class Descriptor>
class IncBGKdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    IncBGKdynamics(T omega_);
    IncBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual IncBGKdynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Say if velocity in this dynamics is computed as "j" (the order-1 moment
    ///   of the populations) or as "j/rho".
    virtual bool velIsJ() const;

    /* *************** Collision and Equilibrium ************************* */

    /// Velocity is equal to j, not u.
    virtual void computeVelocity(
        Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const;

    /// For PiNeq, subtract equilibrium term jj instead of invRho*jj.
    virtual void computeRhoBarJPiNeq(
        Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const;

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

/// Implementation of O(Ma^2) BGK dynamics with constant average density
template <typename T, template <typename U> class Descriptor>
class ConstRhoBGKdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    ConstRhoBGKdynamics(T omega_);
    ConstRhoBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual ConstRhoBGKdynamics<T, Descriptor> *clone() const;

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

/// Generic implementation of the Regularized BGK dynamics
/** This implementation is valid for isothermal models only.
 *  This model is substantially more stable than plain BGK, and has roughly
 *  the same efficiency. However, it cuts out the modes at higher Knudsen
 *  numbers and therefore cannot be used in the regime of rarefied gases.
 */
template <typename T, template <typename U> class Descriptor>
class RLBdynamics : public BulkCompositeDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    RLBdynamics(Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision = true);

    /// Clone the object on its dynamic type.
    virtual RLBdynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /* *************** Completion algorithm and base dynamics ************ */

    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T, Descriptor> &cell) const;

private:
    static int id;
};

/// Implementation of O(Ma^2) regularized BGK dynamics.
/** Semantically, this class is equivalent to RLBdynamics< . , . , BGKdynamics<.,.> >,
 *  but the implementation is more efficient.
 */
template <typename T, template <typename U> class Descriptor>
class RegularizedBGKdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    RegularizedBGKdynamics(T omega_);
    RegularizedBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual RegularizedBGKdynamics<T, Descriptor> *clone() const;

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

/// Implementation of O(Ma^2) stabilized, regularized BGK dynamics.
template <typename T, template <typename U> class Descriptor>
class SecuredRegularizedBGKdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    SecuredRegularizedBGKdynamics(T omega_);
    SecuredRegularizedBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual SecuredRegularizedBGKdynamics<T, Descriptor> *clone() const;

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
    inline static void constrainValue(T &value, T softLimit, T hardLimit);

private:
    static int id;
};

/// Implementation of O(Ma^2) Regularized BGK dynamics without the 1/rho in front of the non-linear
/// term.
template <typename T, template <typename U> class Descriptor>
class IncRegularizedBGKdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    IncRegularizedBGKdynamics(T omega_);
    IncRegularizedBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual IncRegularizedBGKdynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /* *************** Collision and Equilibrium ************************* */

    /// Velocity is equal to j, not u.
    virtual void computeVelocity(
        Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const;

    /// For PiNeq, subtract equilibrium term jj instead of invRho*jj.
    virtual void computeRhoBarJPiNeq(
        Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const;

    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);

    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;

    /// Say if velocity in this dynamics is computed as "j" (the order-1 moment
    ///   of the populations) or as "j/rho".
    bool velIsJ() const;

private:
    static int id;
};

/// Implementation of O(Ma^2) BGK dynamics, density and momentum taken from external scalars
template <typename T, template <typename U> class Descriptor>
class ExternalMomentRegularizedBGKdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    ExternalMomentRegularizedBGKdynamics(T omega_);
    ExternalMomentRegularizedBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual ExternalMomentRegularizedBGKdynamics<T, Descriptor> *clone() const;

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

    /* *************** Moments ******************************************* */

    // The function computeDensity() is not overridden, and the default
    //   implementation is kept, for two reasons. First, it is equivalent to access
    //   rho from the external scalar or to recompute it from the bulk (this
    //   is not the case for the velocity, in a Shan/Chen multicomponent model).
    //   Second, the Shan/Chen data-processor needs computeDensity() to be
    //   default implemented, because it uses this function to treat walls with
    //   a virtual-density mechanism.

    /** Accesses velocity from external scalar. **/
    virtual void computeVelocity(
        Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const;
    /// Compute order-0 moment rho-bar
    /** Accesses rhoBar from external scalar. **/
    virtual T computeRhoBar(Cell<T, Descriptor> const &cell) const;

    /// Compute order-0 moment rho-bar and order-1 moment j
    /** Accesses rhoBar and j from external scalar. **/
    virtual void computeRhoBarJ(
        Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j) const;

private:
    static int id;
};

/// Implementation of O(Ma^2) BGK dynamics with adjustable speed of sound
template <typename T, template <typename U> class Descriptor>
class ChopardDynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    ChopardDynamics(T vs2_, T omega_);
    ChopardDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual ChopardDynamics<T, Descriptor> *clone() const;

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
    static T chopardBgkCollision(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T vs2, T omega);
    /// Implementation of equilibrium
    static T chopardEquilibrium(
        plint iPop, T rhoBar, T invRho, Array<T, Descriptor<T>::d> const &j, T jSqr, T vs2);

private:
    T vs2;  ///< speed of sound
private:
    static int id;
};

/// Implementation of O(Ma^2) BGK dynamics with preconditioning
template <typename T, template <typename U> class Descriptor>
class PrecondBGKdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    PrecondBGKdynamics(T omega_, T invGamma_);
    PrecondBGKdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual PrecondBGKdynamics<T, Descriptor> *clone() const;

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

private:
    static int id;
    T invGamma;
};

}  // namespace plb

#endif  // ISO_THERMAL_DYNAMICS_H
