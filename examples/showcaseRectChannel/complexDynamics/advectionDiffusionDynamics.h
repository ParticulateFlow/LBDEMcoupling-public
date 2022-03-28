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
#ifndef ADVECTION_DIFFUSION_DYNAMICS_H
#define ADVECTION_DIFFUSION_DYNAMICS_H

#include "core/dynamics.h"
#include "core/globalDefs.h"

namespace plb {

/// Common base iso-thermal (or athermal) bulk dynamics
template <typename T, template <typename U> class Descriptor>
class AdvectionDiffusionDynamics : public BasicBulkDynamics<T, Descriptor> {
public:
    AdvectionDiffusionDynamics(T omega_);

    /* *************** Collision, Equilibrium, and Non-equilibrium ******* */

    /// Re-compute particle populations from the leading moments
    virtual void regularize(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
        Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar = T()) const;

    /* *************** Additional moments, intended for internal use ***** */

    /// Returns 0, as a default value for isothermal flow.
    virtual T computeEbar(Cell<T, Descriptor> const &cell) const;

    /* *************** Switch between population and moment representation ****** */

    /// Number of variables required to decompose a population representation into moments.
    virtual plint numDecomposedVariables(plint order) const;

    /// Decompose from population representation into moment representation.
    virtual void decompose(
        Cell<T, Descriptor> const &cell, std::vector<T> &rawData, plint order) const;

    /// Recompose from moment representation to population representation.
    virtual void recompose(
        Cell<T, Descriptor> &cell, std::vector<T> const &rawData, plint order) const;

    /// Change the space and time scales of the variables in moment representation.
    virtual void rescale(std::vector<T> &rawData, T xDxInv, T xDt, plint order) const;

    virtual bool isAdvectionDiffusion() const
    {
        return true;
    }

    virtual void computeVelocity(
        Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const;
};

/// Regularized Advection-Diffusion dynamics
/** It uses the regularized approximation that can be found in
 *   the thesis of J. Latt (2007).
 */
template <typename T, template <typename U> class Descriptor>
class AdvectionDiffusionRLBdynamics : public AdvectionDiffusionDynamics<T, Descriptor> {
public:
    /// Constructor
    AdvectionDiffusionRLBdynamics(T omega_);
    AdvectionDiffusionRLBdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual AdvectionDiffusionRLBdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics);
    /// Implementation of the collision step, with imposed macroscopic variables
    /// The arguments:
    /// - rhoBar: the "rhoBar" version of the scalar rho.
    /// - jEq: the equilibrium part of the second-order moment. jEq = u*rho, where u is the external
    /// convective term.
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T thetaBar,
        BlockStatistics &stat);
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T jSqr,
        T thetaBar = T()) const;

private:
    static int id;
};

/// SRT Advection-Diffusion dynamics by Perko & Patel PRE 89, 053309 (2014).
template <typename T, template <typename U> class Descriptor>
class AdvectionDiffusionPerkoDynamics : public AdvectionDiffusionDynamics<T, Descriptor> {
public:
    /// Constructor
    AdvectionDiffusionPerkoDynamics(T omega_, T omegaRef_);
    AdvectionDiffusionPerkoDynamics(HierarchicUnserializer &unserializer);
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual AdvectionDiffusionPerkoDynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics);
    /// Implementation of the collision step, with imposed macroscopic variables
    /// The arguments:
    /// - rhoBar: the "rhoBar" version of the scalar rho.
    /// - jEq: the equilibrium part of the second-order moment. jEq = u*rho, where u is the external
    /// convective term.
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T thetaBar,
        BlockStatistics &stat);
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T jSqr,
        T thetaBar = T()) const;

private:
    static int id;
    T omegaRef;
};

template <typename T, template <typename U> class Descriptor>
class AdvectionDiffusionWithSourceLinearRLBdynamics :
    public AdvectionDiffusionDynamics<T, Descriptor> {
public:
    /// Constructor
    AdvectionDiffusionWithSourceLinearRLBdynamics(T omega_);
    AdvectionDiffusionWithSourceLinearRLBdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual AdvectionDiffusionWithSourceLinearRLBdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics);
    /// Implementation of the collision step, with imposed macroscopic variables
    /// The arguments:
    /// - rhoBar: the "rhoBar" version of the scalar rho.
    /// - jEq: the equilibrium part of the second-order moment. jEq = u*rho, where u is the external
    /// convective term.
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T thetaBar,
        BlockStatistics &stat);
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T jSqr,
        T thetaBar = T()) const;

private:
    static int id;
};

template <typename T, template <typename U> class Descriptor>
class AdvectionDiffusionWithSourceRLBdynamics : public AdvectionDiffusionDynamics<T, Descriptor> {
public:
    /// Constructor
    AdvectionDiffusionWithSourceRLBdynamics(T omega_);
    AdvectionDiffusionWithSourceRLBdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual AdvectionDiffusionWithSourceRLBdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics);
    /// Implementation of the collision step, with imposed macroscopic variables
    /// The arguments:
    /// - rhoBar: the "rhoBar" version of the scalar rho.
    /// - jEq: the equilibrium part of the second-order moment. jEq = u*rho, where u is the external
    /// convective term.
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T thetaBar,
        BlockStatistics &stat);
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T jSqr,
        T thetaBar = T()) const;

private:
    static int id;
};

/// SRT Advection-Diffusion dynamics by Perko & Patel PRE 89, 053309 (2014),
/// modified to account for a source term.
template <typename T, template <typename U> class Descriptor>
class AdvectionDiffusionWithSourcePerkoDynamics : public AdvectionDiffusionDynamics<T, Descriptor> {
public:
    /// Constructor
    AdvectionDiffusionWithSourcePerkoDynamics(T omega_, T omegaRef_);
    AdvectionDiffusionWithSourcePerkoDynamics(HierarchicUnserializer &unserializer);
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual AdvectionDiffusionWithSourcePerkoDynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics);
    /// Implementation of the collision step, with imposed macroscopic variables
    /// The arguments:
    /// - rhoBar: the "rhoBar" version of the scalar rho.
    /// - jEq: the equilibrium part of the second-order moment. jEq = u*rho, where u is the external
    /// convective term.
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T thetaBar,
        BlockStatistics &stat);
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T jSqr,
        T thetaBar = T()) const;

private:
    static int id;
    T omegaRef;
};

/// Regularized Advection-Diffusion dynamics with artificial diffusivity as in the Smagorinsky
/// model.
template <typename T, template <typename U> class Descriptor>
class SmagorinskyAdvectionDiffusionRLBdynamics : public AdvectionDiffusionDynamics<T, Descriptor> {
public:
    /// Constructor
    SmagorinskyAdvectionDiffusionRLBdynamics(T omega_, T T0_, T cSmago_);
    /// Constructor from a serialized object.
    SmagorinskyAdvectionDiffusionRLBdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual SmagorinskyAdvectionDiffusionRLBdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);
    /// Collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics);
    /// Implementation of the collision step, with imposed macroscopic variables
    /// The arguments:
    /// - rhoBar: the "rhoBar" version of the scalar rho.
    /// - jEq: the equilibrium part of the second-order moment. jEq = u*rho, where u is the external
    /// convective term.
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T thetaBar,
        BlockStatistics &stat);
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T jSqr,
        T thetaBar = T()) const;

private:
    T invT0;
    T cSmago;
    static int id;
};

/// BGK Advection-Diffusion dynamics
/** This approach contains a slight error in the diffusion
 *  term.
 */
template <typename T, template <typename U> class Descriptor>
class AdvectionDiffusionBGKdynamics : public AdvectionDiffusionDynamics<T, Descriptor> {
public:
    /// Constructor
    AdvectionDiffusionBGKdynamics(T omega_);
    AdvectionDiffusionBGKdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual AdvectionDiffusionBGKdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics);
    /// Implementation of the collision step, with imposed macroscopic variables
    /// The arguments:
    /// - rhoBar: the "rhoBar" version of the scalar rho.
    /// - jEq: the equilibrium part of the second-order moment. jEq = u*rho, where u is the external
    /// convective term.
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T thetaBar,
        BlockStatistics &stat);
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T jSqr,
        T thetaBar = T()) const;

private:
    static int id;
};

/// Complete BGK Advection-Diffusion dynamics
/** This approach contains a slight error in the diffusion
 *  term. We tried to reduce it with the extended exquilibrium distribution
 */
template <typename T, template <typename U> class Descriptor>
class CompleteAdvectionDiffusionBGKdynamics : public AdvectionDiffusionDynamics<T, Descriptor> {
public:
    /// Constructor
    CompleteAdvectionDiffusionBGKdynamics(T omega_);
    CompleteAdvectionDiffusionBGKdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual CompleteAdvectionDiffusionBGKdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Computation of the density field (sum_i f_i = rho*phi), phi is the advected diffused field
    /// rho the density of the fluid
    virtual T computeDensity(Cell<T, Descriptor> const &cell) const;

    /// Collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics);
    /// Implementation of the collision step, with imposed macroscopic variables
    /// The arguments:
    /// - rhoBar: the "rhoBar" version of the scalar rho.
    /// - j: the equilibrium part of the second-order moment. j = u*rho, where u is the external
    /// convective term.
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;

    /* *************** Switch between population and moment representation ****** */
    /*
        /// Number of variables required to decompose a population representation into moments.
        virtual plint numDecomposedVariables(plint order) const;

        /// Decompose from population representation into moment representation.
        virtual void decompose(Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order)
       const;

        /// Recompose from moment representation to population representation.
        virtual void recompose(Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order)
       const;*/

    /// Re-compute particle populations from the leading moments
    virtual void regularize(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
        Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar = T()) const;

private:
    static int id;
};

/// Complete BGK Advection-Diffusion dynamics
/** This approach contains a slight error in the diffusion
 *  term. We tried to reduce it with the extended exquilibrium distribution
 */
template <typename T, template <typename U> class Descriptor>
class CompleteRegularizedAdvectionDiffusionBGKdynamics :
    public CompleteAdvectionDiffusionBGKdynamics<T, Descriptor> {
public:
    /// Constructor
    CompleteRegularizedAdvectionDiffusionBGKdynamics(T omega_);
    CompleteRegularizedAdvectionDiffusionBGKdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual CompleteRegularizedAdvectionDiffusionBGKdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Computation of the density field (sum_i f_i = rho*phi), phi is the advected diffused field
    /// rho the density of the fluid
    virtual T computeDensity(Cell<T, Descriptor> const &cell) const;

    /// Collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics);

    /* *************** Switch between population and moment representation ****** */
    /*
        /// Number of variables required to decompose a population representation into moments.
        virtual plint numDecomposedVariables(plint order) const;

        /// Decompose from population representation into moment representation.
        virtual void decompose(Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order)
       const;

        /// Recompose from moment representation to population representation.
        virtual void recompose(Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order)
       const;*/
private:
    static int id;
};

/// Complete BGK Advection-Diffusion dynamics with source
/** This approach contains a slight error in the diffusion
 *  term. We tried to reduce it with the extended exquilibrium distribution
 */
template <typename T, template <typename U> class Descriptor>
class CompleteRegularizedAdvectionDiffusionWithSourceBGKdynamics :
    public CompleteAdvectionDiffusionBGKdynamics<T, Descriptor> {
public:
    /// Constructor
    CompleteRegularizedAdvectionDiffusionWithSourceBGKdynamics(T omega_);
    CompleteRegularizedAdvectionDiffusionWithSourceBGKdynamics(
        HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual CompleteRegularizedAdvectionDiffusionWithSourceBGKdynamics<T, Descriptor> *clone()
        const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Computation of the density field (sum_i f_i = rho*phi), phi is the advected diffused field
    /// rho the density of the fluid
    virtual T computeDensity(Cell<T, Descriptor> const &cell) const;

    /// Collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics);

    /* *************** Switch between population and moment representation ****** */
    /*
        /// Number of variables required to decompose a population representation into moments.
        virtual plint numDecomposedVariables(plint order) const;

        /// Decompose from population representation into moment representation.
        virtual void decompose(Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order)
       const;

        /// Recompose from moment representation to population representation.
        virtual void recompose(Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order)
       const;*/
private:
    static int id;
};

/// Complete TRT Advection-Diffusion dynamics
/** This approach contains a slight error in the diffusion
 *  term. We tried to reduce it with the extended exquilibrium distribution
 */
template <typename T, template <typename U> class Descriptor>
class CompleteAdvectionDiffusionTRTdynamics : public AdvectionDiffusionDynamics<T, Descriptor> {
public:
    /// Constructor
    CompleteAdvectionDiffusionTRTdynamics(T omega_, T psi_);
    CompleteAdvectionDiffusionTRTdynamics(T omega_);
    CompleteAdvectionDiffusionTRTdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual CompleteAdvectionDiffusionTRTdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;

    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);

    /// Computation of the density field (sum_i f_i = rho*phi), phi is the advected diffused field
    /// rho the density of the fluid
    virtual T computeDensity(Cell<T, Descriptor> const &cell) const;

    /// Collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics);
    /// Implementation of the collision step, with imposed macroscopic variables
    /// The arguments:
    /// - rhoBar: the "rhoBar" version of the scalar rho.
    /// - j: the equilibrium part of the second-order moment. j = u*rho, where u is the external
    /// convective term.
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
    void setPsi(T psi_);
    /// Get local speed of sound
    T getPsi() const;

    /* *************** Switch between population and moment representation ****** */
    /*
        /// Number of variables required to decompose a population representation into moments.
        virtual plint numDecomposedVariables(plint order) const;

        /// Decompose from population representation into moment representation.
        virtual void decompose(Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order)
       const;

        /// Recompose from moment representation to population representation.
        virtual void recompose(Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order)
       const;*/

    /// Re-compute particle populations from the leading moments
    virtual void regularize(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
        Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar = T()) const;

private:
    static int id;
    T psi;
};

/// BGK Advection-Diffusion dynamics
/** This approach contains a slight error in the diffusion
 *  term.
 */
template <typename T, template <typename U> class Descriptor>
class AdvectionDiffusionWithSourceBGKdynamics : public AdvectionDiffusionDynamics<T, Descriptor> {
public:
    /// Constructor
    AdvectionDiffusionWithSourceBGKdynamics(T omega_);
    AdvectionDiffusionWithSourceBGKdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual AdvectionDiffusionWithSourceBGKdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics);
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &jEq, T jSqr,
        T thetaBar = T()) const;

private:
    static int id;
};

}  // namespace plb

#endif  // ADVECTION_DIFFUSION_DYNAMICS_H
