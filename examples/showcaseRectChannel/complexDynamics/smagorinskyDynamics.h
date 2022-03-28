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

/* Orestis Malaspinas designed some of the classes and concepts contained
 * in this file. */

#ifndef SMAGORINSKY_DYNAMICS_H
#define SMAGORINSKY_DYNAMICS_H

#include "basicDynamics/externalForceDynamics.h"
#include "basicDynamics/isoThermalDynamics.h"
#include "complexDynamics/mrtDynamics.h"
#include "complexDynamics/variableOmegaDynamics.h"
#include "core/globalDefs.h"
#include "core/hierarchicSerializer.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class SmagorinskyDynamics : public OmegaFromPiDynamics<T, Descriptor> {
public:
    SmagorinskyDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, T omega0_, T cSmago_,
        bool automaticPrepareCollision = true);
    SmagorinskyDynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    SmagorinskyDynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);
    /// Modify the value of omega, using the Smagorinsky algorithm based on omega0.
    virtual T getOmegaFromPiAndRhoBar(
        Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T rhoBar) const;
    /// With this method, you can modify the constant value omega0 (not the actual value of omega,
    ///  which is computed during run-time from omega0 and the local strain-rate).
    virtual void setOmega(T omega_);
    /// Returns omega0.
    virtual T getOmega() const;
    /// Return dynamic value of omega for whichParameter=dynamicOmega or smagorinskyConstant.
    virtual T getDynamicParameter(plint whichParameter, Cell<T, Descriptor> const &cell) const;

private:
    T omega0;     //< "Laminar" relaxation parameter, used when the strain-rate is zero.
    T cSmago;     //< Smagorinsky constant.
    T preFactor;  //< A factor depending on the Smagorinky constant, used to compute the effective
                  // viscosity.
    static int id;
};

template <typename T, template <typename U> class Descriptor>
class SmagorinskyBGKdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    SmagorinskyBGKdynamics(T omega0_, T cSmago_);
    SmagorinskyBGKdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual SmagorinskyBGKdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);
    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);
    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;
    /// With this method, you can modify the constant value omega0 (not the actual value of omega,
    ///  which is computed during run-time from omega0 and the local strain-rate).
    virtual void setOmega(T omega_);
    /// Returns omega0.
    virtual T getOmega() const;
    /// Return dynamic value of omega for whichParameter=dynamicOmega.
    virtual T getDynamicParameter(plint whichParameter, Cell<T, Descriptor> const &cell) const;

private:
    T omega0;     //< "Laminar" relaxation parameter, used when the strain-rate is zero.
    T cSmago;     //< Smagorinsky constant.
    T preFactor;  //< A factor depending on the Smagorinky constant, used to compute the effective
                  // viscosity.
    static int id;
};

/// Consistent according to Orestis' model: doesn't mix the stress
/// modes with the other ones when you change the relaxation time.
template <typename T, template <typename U> class Descriptor>
class ConsistentSmagorinskyBGKdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    ConsistentSmagorinskyBGKdynamics(T omega0_, T cSmago_);
    ConsistentSmagorinskyBGKdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual ConsistentSmagorinskyBGKdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);
    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);
    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;
    /// With this method, you can modify the constant value omega0 (not the actual value of omega,
    ///  which is computed during run-time from omega0 and the local strain-rate).
    virtual void setOmega(T omega_);
    /// Returns omega0.
    virtual T getOmega() const;
    /// Return dynamic value of omega for whichParameter=dynamicOmega or smagorinskyConstant.
    virtual T getDynamicParameter(plint whichParameter, Cell<T, Descriptor> const &cell) const;
    /// Get local value of any generic parameter.
    /// For the Smagorinsky constant cSmago, use dynamicParams::smagorinskyConstant.
    virtual T getParameter(plint whichParameter) const;

private:
    T omega0;  //< "Laminar" relaxation parameter, used when the strain-rate is zero.
    T cSmago;  //< Smagorinsky constant.
    static int id;
};

/// Consistent according to Orestis' model: doesn't mix the stress
/// modes with the other ones when you change the relaxation time.
template <typename T, template <typename U> class Descriptor>
class ConsistentSmagorinskyCompleteTRTdynamics : public CompleteTRTdynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    ConsistentSmagorinskyCompleteTRTdynamics(T omega_, T psi_, T cSmago_);
    ConsistentSmagorinskyCompleteTRTdynamics(T omega_, T cSmago_);
    ConsistentSmagorinskyCompleteTRTdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual ConsistentSmagorinskyCompleteTRTdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);
    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);
    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);
    /// Return dynamic value of omega for whichParameter=dynamicOmega or smagorinskyConstant.
    virtual T getDynamicParameter(plint whichParameter, Cell<T, Descriptor> const &cell) const;

private:
    T cSmago;  //< Smagorinsky constant.
    static int id;
};

/// Consistent according to Orestis' model: doesn't mix the stress
/// modes with the other ones when you change the relaxation time.
template <typename T, template <typename U> class Descriptor>
class ConsistentSmagorinskyCompleteRegularizedTRTdynamics :
    public CompleteRegularizedTRTdynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    ConsistentSmagorinskyCompleteRegularizedTRTdynamics(T omega_, T psi_, int order_, T cSmago_);
    ConsistentSmagorinskyCompleteRegularizedTRTdynamics(T omega_, int order_, T cSmago_);
    ConsistentSmagorinskyCompleteRegularizedTRTdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual ConsistentSmagorinskyCompleteRegularizedTRTdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);
    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);
    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);
    /// Return dynamic value of omega for whichParameter=dynamicOmega or smagorinskyConstant.
    virtual T getDynamicParameter(plint whichParameter, Cell<T, Descriptor> const &cell) const;

private:
    template <typename U>
    ConsistentSmagorinskyCompleteRegularizedTRTdynamics(T omega_, T psi_, U order_, T cSmago_);

    template <typename U>
    ConsistentSmagorinskyCompleteRegularizedTRTdynamics(T omega_, U order_, T cSmago_);

private:
    int order;
    T cSmago;  //< Smagorinsky constant.
    static int id;
};

/// Consistent according to Orestis' model: doesn't mix the stress
/// modes with the other ones when you change the relaxation time.
template <typename T, template <typename U> class Descriptor>
class ConsistentSmagorinskyCompleteRegularizedBGKdynamics :
    public CompleteRegularizedBGKdynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    ConsistentSmagorinskyCompleteRegularizedBGKdynamics(T omega_, T cSmago_);
    ConsistentSmagorinskyCompleteRegularizedBGKdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual ConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);
    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);
    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);
    /// Return dynamic value of omega for whichParameter=dynamicOmega or smagorinskyConstant.
    virtual T getDynamicParameter(plint whichParameter, Cell<T, Descriptor> const &cell) const;

private:
    T cSmago;  //< Smagorinsky constant.
    static int id;
};

/// Consistent according to Orestis' model: doesn't mix the stress
/// modes with the other ones when you change the relaxation time.
template <typename T, template <typename U> class Descriptor>
class ConsistentSmagorinskyTruncatedTRTdynamics : public TruncatedTRTdynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    ConsistentSmagorinskyTruncatedTRTdynamics(T omega_, T psi_, T cSmago_);
    ConsistentSmagorinskyTruncatedTRTdynamics(T omega_, T cSmago_);
    ConsistentSmagorinskyTruncatedTRTdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual ConsistentSmagorinskyTruncatedTRTdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);
    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);
    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);
    /// Return dynamic value of omega for whichParameter=dynamicOmega or smagorinskyConstant.
    virtual T getDynamicParameter(plint whichParameter, Cell<T, Descriptor> const &cell) const;

private:
    T cSmago;  //< Smagorinsky constant.
    static int id;
};

/// Consistent according to Orestis' model: doesn't mix the stress
/// modes with the other ones when you change the relaxation time.
/// Allows to put any subgridstress tensor from an external field.
template <typename T, template <typename U> class Descriptor>
class ConsistentSgsBGKdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    ConsistentSgsBGKdynamics(T omega0_);
    ConsistentSgsBGKdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual ConsistentSgsBGKdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);
    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);
    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;
    /// With this method, you can modify the constant value omega0 (not the actual value of omega,
    ///  which is computed during run-time from omega0 and the local strain-rate).
    virtual void setOmega(T omega_);
    /// Returns omega0.
    virtual T getOmega() const;
    /// Return dynamic value of omega for whichParameter=dynamicOmega or smagorinskyConstant.
    virtual T getDynamicParameter(plint whichParameter, Cell<T, Descriptor> const &cell) const;
    /// Get local value of any generic parameter.
    /// For the Smagorinsky constant cSmago, use dynamicParams::smagorinskyConstant.
    virtual T getParameter(plint whichParameter) const;

private:
    T omega0;  //< "Laminar" relaxation parameter, used when the strain-rate is zero.
    static int id;
};

template <typename T, template <typename U> class Descriptor>
class GuoExternalForceSmagorinskyBGKdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    GuoExternalForceSmagorinskyBGKdynamics(T omega0_, T cSmago_);
    GuoExternalForceSmagorinskyBGKdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual GuoExternalForceSmagorinskyBGKdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);
    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);
    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);
    /// Compute the physical velocity
    virtual void computeVelocity(
        Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const;
    /// Compute the physical velocity given rhoBar and j
    virtual void computeVelocityExternal(
        Cell<T, Descriptor> const &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j,
        Array<T, Descriptor<T>::d> &u) const;
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;
    /// With this method, you can modify the constant value omega0 (not the actual value of omega,
    ///  which is computed during run-time from omega0 and the local strain-rate).
    virtual void setOmega(T omega_);
    /// Returns omega0.
    virtual T getOmega() const;
    /// Set local value of any generic parameter.
    /// For the Smagorinsky constant cSmago, use dynamicParams::smagorinskyConstant.
    virtual void setParameter(plint whichParameter, T value);

private:
    T omega0;     //< "Laminar" relaxation parameter, used when the strain-rate is zero.
    T cSmago;     //< Smagorinsky constant.
    T preFactor;  //< A factor depending on the Smagorinky constant, used to compute the effective
                  // viscosity.
    static int id;
};

template <typename T, template <typename U> class Descriptor>
class SmagorinskyIncBGKdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    SmagorinskyIncBGKdynamics(T omega0_, T cSmago_);
    SmagorinskyIncBGKdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual SmagorinskyIncBGKdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);
    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);
    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);
    /// Compute the physical velocity
    virtual void computeVelocity(
        Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const;
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;
    /// For PiNeq, subtract equilibrium term jj instead of invRho*jj.
    virtual void computeRhoBarJPiNeq(
        Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const;
    /// With this method, you can modify the constant value omega0 (not the actual value of omega,
    ///  which is computed during run-time from omega0 and the local strain-rate).
    virtual void setOmega(T omega_);
    /// Returns omega0.
    virtual T getOmega() const;

    /// Get local value of any generic parameter.
    /// For the Smagorinsky constant cSmago, use dynamicParams::smagorinskyConstant.
    virtual T getParameter(plint whichParameter) const;

    /// Set local value of any generic parameter.
    /// For the Smagorinsky constant cSmago, use dynamicParams::smagorinskyConstant.
    virtual void setParameter(plint whichParameter, T value);
    ///  Is it in the incompressible BGK model
    virtual bool velIsJ() const;

private:
    T omega0;     //< "Laminar" relaxation parameter, used when the strain-rate is zero.
    T cSmago;     //< Smagorinsky constant.
    T preFactor;  //< A factor depending on the Smagorinky constant, used to compute the effective
                  // viscosity.
    static int id;
};

template <typename T, template <typename U> class Descriptor>
class SmagorinskyRegularizedDynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    SmagorinskyRegularizedDynamics(T omega0_, T cSmago_);
    SmagorinskyRegularizedDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual SmagorinskyRegularizedDynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);
    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);
    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;
    /// With this method, you can modify the constant value omega0 (not the actual value of omega,
    ///  which is computed during run-time from omega0 and the local strain-rate).
    virtual void setOmega(T omega_);
    /// Returns omega0.
    virtual T getOmega() const;

private:
    T omega0;     //< "Laminar" relaxation parameter, used when the strain-rate is zero.
    T cSmago;     //< Smagorinsky constant.
    T preFactor;  //< A factor depending on the Smagorinky constant, used to compute the effective
                  // viscosity.
    static int id;
};

template <typename T, template <typename U> class Descriptor>
class SecuredSmagorinskyRegularizedDynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    SecuredSmagorinskyRegularizedDynamics(T omega0_, T cSmago_);
    SecuredSmagorinskyRegularizedDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual SecuredSmagorinskyRegularizedDynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);
    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);
    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;
    /// With this method, you can modify the constant value omega0 (not the actual value of omega,
    ///  which is computed during run-time from omega0 and the local strain-rate).
    virtual void setOmega(T omega_);
    /// Returns omega0.
    virtual T getOmega() const;

private:
    inline static void constrainValue(T &value, T softLimit, T hardLimit);

private:
    T omega0;     //< "Laminar" relaxation parameter, used when the strain-rate is zero.
    T cSmago;     //< Smagorinsky constant.
    T preFactor;  //< A factor depending on the Smagorinky constant, used to compute the effective
                  // viscosity.
    static int id;
};

/// Smagorinky model where only the shear viscosity omega is modified
/// in order to obtain a Smagorinsky model
template <typename T, template <typename U> class Descriptor>
class SmagorinskyMRTdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    SmagorinskyMRTdynamics(T omega0_, T cSmago_);
    SmagorinskyMRTdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual SmagorinskyMRTdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);
    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);
    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;
    /// With this method, you can modify the constant value omega0 (not the actual value of omega,
    ///  which is computed during run-time from omega0 and the local strain-rate).
    virtual void setOmega(T omega_);
    /// Returns omega0.
    virtual T getOmega() const;
    /// Return dynamic value of omega for whichParameter=dynamicOmega.
    virtual T getDynamicParameter(plint whichParameter, Cell<T, Descriptor> const &cell) const;

private:
    T omega0;     //< "Laminar" relaxation parameter, used when the strain-rate is zero.
    T cSmago;     //< Smagorinsky constant.
    T preFactor;  //< A factor depending on the Smagorinky constant, used to compute the effective
                  // viscosity.
    static int id;
};

template <typename T, template <typename U> class Descriptor>
class ConsistentSmagorinskyMRTdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    ConsistentSmagorinskyMRTdynamics(T omega0_, T cSmago_);
    ConsistentSmagorinskyMRTdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual ConsistentSmagorinskyMRTdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);
    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);
    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;
    /// Return dynamic value of omega for whichParameter=dynamicOmega or smagorinskyConstant.
    virtual T getDynamicParameter(plint whichParameter, Cell<T, Descriptor> const &cell) const;
    /// Get local value of any generic parameter.
    /// For the Smagorinsky constant cSmago, use dynamicParams::smagorinskyConstant.
    virtual T getParameter(plint whichParameter) const;

private:
    T cSmago;  //< Smagorinsky constant.
    static int id;
};

template <typename T, template <typename U> class Descriptor>
class ConsistentSgsMRTdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    ConsistentSgsMRTdynamics(T omega0_);
    ConsistentSgsMRTdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual ConsistentSgsMRTdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);
    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;
    /// Return dynamic value of omega for whichParameter=dynamicOmega or smagorinskyConstant.
    virtual T getDynamicParameter(plint whichParameter, Cell<T, Descriptor> const &cell) const;
    /// Get local value of any generic parameter.
    /// For the Sgs constant cSmago, use dynamicParams::smagorinskyConstant.
    virtual T getParameter(plint whichParameter) const;

private:
    static int id;
};

template <typename T, template <typename U> class Descriptor>
class GuoExternalForceConsistentSgsMRTdynamics : public ExternalForceDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    GuoExternalForceConsistentSgsMRTdynamics(T omega0_);
    GuoExternalForceConsistentSgsMRTdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual GuoExternalForceConsistentSgsMRTdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;

private:
    static int id;
};

template <typename T, template <typename U> class Descriptor>
class ConsistentSmagorinskyIncMRTdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    ConsistentSmagorinskyIncMRTdynamics(T omega0_, T cSmago_);
    ConsistentSmagorinskyIncMRTdynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual ConsistentSmagorinskyIncMRTdynamics<T, Descriptor> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);

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
    /// Return dynamic value of omega for whichParameter=dynamicOmega or smagorinskyConstant.
    virtual T getDynamicParameter(plint whichParameter, Cell<T, Descriptor> const &cell) const;
    /// Get local value of any generic parameter.
    /// For the Smagorinsky constant cSmago, use dynamicParams::smagorinskyConstant.
    virtual T getParameter(plint whichParameter) const;

private:
    T cSmago;  //< Smagorinsky constant.
    static int id;
};

}  // namespace plb

#endif  // SMAGORINSKY_DYNAMICS_H
