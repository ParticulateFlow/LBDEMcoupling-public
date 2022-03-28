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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- header file.
 */
#ifndef ENTROPIC_LB_DYNAMICS_H
#define ENTROPIC_LB_DYNAMICS_H

#include "basicDynamics/isoThermalDynamics.h"
#include "complexDynamics/variableOmegaDynamics.h"
#include "core/globalDefs.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class Cell;

/// Implementation of the entropic collision step
template <typename T, template <typename U> class Descriptor>
class EntropicDynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    EntropicDynamics(T omega_);

    EntropicDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual EntropicDynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    virtual bool isEntropic() const
    {
        return true;
    }

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

/// Implementation of the forced entropic collision step
template <typename T, template <typename U> class Descriptor>
class ForcedEntropicDynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    ForcedEntropicDynamics(T omega_);

    ForcedEntropicDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual ForcedEntropicDynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    virtual bool isEntropic() const
    {
        return true;
    }

    /* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;

private:
    /// computes the entropy function H(f)=sum_i f_i*ln(f_i/t_i)
    T computeEntropy(Array<T, Descriptor<T>::q> const &f);
    /// computes the entropy growth H(f)-H(f-alpha*fNeq)
    T computeEntropyGrowth(
        Array<T, Descriptor<T>::q> const &f, Array<T, Descriptor<T>::q> const &fNeq, T alpha);
    /// computes the entropy growth derivative
    /// dH/dalpha=-sum_i fNeq_i*ln((f_i-alpha*fNeq_i)/t_i)
    T computeEntropyGrowthDerivative(
        Array<T, Descriptor<T>::q> const &f, Array<T, Descriptor<T>::q> const &fNeq, T alpha);
    /// Get the alpha parameter
    bool getAlpha(
        T &alpha, Array<T, Descriptor<T>::q> const &f, Array<T, Descriptor<T>::q> const &fNeq);

    static const int forceBeginsAt = Descriptor<T>::ExternalField::forceBeginsAt;
    static const int sizeOfForce = Descriptor<T>::ExternalField::sizeOfForce;
    static int id;
};

/*******************************************************************************
 *
 * A version of ELBM, which simply modifies the local omega by alpha
 * (may be needed for grid refinement?)
 *
 * ****************************************************************************/
/// Implementation of the entropic collision step
template <typename T, template <typename U> class Descriptor>
class VariableOmegaELBMDynamics : public VariableOmegaDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    VariableOmegaELBMDynamics(T omega0_, bool automaticPrepareCollision = true);

    VariableOmegaELBMDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual VariableOmegaELBMDynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    virtual bool isEntropic() const
    {
        return true;
    }

    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);

    /// With this method, you can modify the constant value omega0 (not the actual value of omega,
    ///  which is computed during run-time from omega0 and the local strain-rate).
    virtual void setOmega(T omega_);
    /// Returns omega0.
    virtual T getOmega() const;

    virtual T getOmegaFromCell(Cell<T, Descriptor> const &cell) const;

private:
    T omega0;
    static int id;
};

}  // namespace plb

#endif  // ENTROPIC_LB_DYNAMICS_H
