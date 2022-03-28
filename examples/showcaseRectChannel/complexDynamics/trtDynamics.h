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

#ifndef TRT_DYNAMICS_H
#define TRT_DYNAMICS_H

#include "basicDynamics/isoThermalDynamics.h"
#include "core/globalDefs.h"

namespace plb {

/// Implementation of the TRT collision step
template <typename T, template <typename U> class Descriptor>
class BaseTRTdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /// If a value of 0 is provided for omegaMinus, this parameter will default
    /// to a value depending on omega, which guarantees that the wall is
    /// exactly half-way between nodes for parallel walls, according to Eq. (49) in Ginzburg and
    /// Adler 1994 [1] [https://doi.org/10.1051/jp2:1994123], (where lambda_{2c}->-omegaMinus and
    /// lambda_{psi}->-omegaPlus). It corresponds to a "magic parameter" of Lambda = 3/16 (see Eq.
    /// (10.43) pag 426 of [2]). See also Eq. 11 in Pan et al 2006 [3] for application to porous
    /// media \param omegaPlus_       viscosity relaxation rate \param omegaMinus_      relaxation
    /// time associated to the anti-symmetric part. \param constant_magic   if set to true (default)
    /// it recomputes omegaMinus every time setOmega() is called in order
    ///                         to keep the original magic parameter
    /// \examples 1. set omegaPlus and omegaMinus and disable auto-recalculation of omegaMinus:
    /// BaseTRTdynamics(omegaPlus_,omegaMinus_, false) \examples 2. set omegaPlus and omegaMinus
    /// with auto-recalculation of omegaMinus: BaseTRTdynamics(omegaPlus_,omegaMinus_) \examples 3.
    /// set omegaPlus and magicParameter: auto trt = new BaseTRTdynamics<T,Descriptor>(omegaPlus_);
    /// trt->setMagicParam(magic); \references [1] I. Ginzbourg, P. M. Adler, “Boundary flow
    /// condition analysis for the three-dimensional lattice Boltzmann model,” J. Phys. II France,
    /// vol. 4, no. 2, pp. 191–214, Feb. 1994, doi: 10.1051/jp2:1994123.    /// \references
    /// [2]Kruger et al., The lattice Boltzmann method: principles and practice. New York, NY:
    /// Springer Berlin Heidelberg, 2016.    /// \references [3] C. Pan, L.S. Luo, C. T. Miller,“An
    /// evaluation of lattice Boltzmann schemes for porous medium flow simulation,”
    ///      Computers & Fluids, vol. 35, no. 8, pp. 898–909, Sep. 2006,
    ///      doi: 10.1016/j.compfluid.2005.03.008. https://doi.org/10.1016/j.compfluid.2005.03.008
    explicit BaseTRTdynamics(T omegaPlus_, T omegaMinus_ = T(), bool constant_magic = true);
    explicit BaseTRTdynamics(HierarchicUnserializer &unserializer);

    /// Serialize the dynamics object.
    void serialize(HierarchicSerializer &serializer) const override;
    /// Un-Serialize the dynamics object.
    void unserialize(HierarchicUnserializer &unserializer) override;
    /// Clone the object on its dynamic type.
    virtual BaseTRTdynamics<T, Descriptor> *clone() const = 0;

    /* *************** Access to Dynamics variables, e.g. omega ***************** */
    /// Get local relaxation parameter of the dynamics
    virtual T getOmegaMinus() const;

    /// Set local relaxation parameter of the dynamics
    void setOmega(T omega_) override;

    T getParameter(plint whichParameter) const override;
    /// Set local value of any generic parameter.
    void setParameter(plint whichParameter, T value) override;

    /// Set local relaxation parameter of the dynamics
    virtual void setOmegaMinus(T omegaMinus_);

    /// Get the value of the magic parameter from the value of omega and omegaMinus
    virtual T getMagicParam() const;

    /// Set omegaMinus according to the value of the provided magic parameter (see Eq. (10.43) pag
    /// 426 of [2])
    virtual void setMagicParam(T magic_);

    /* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_) = 0;

    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat) = 0;

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
        T thetaBar = T()) const = 0;

private:
    T omegaMinus;
    bool keep_magic_constant_when_setting_omega;
};

template <typename T, template <typename U> class Descriptor>
class TRTdynamics : public BaseTRTdynamics<T, Descriptor> {
public:
    // inherit constructors
    using BaseTRTdynamics<T, Descriptor>::BaseTRTdynamics;
    /// Clone the object on its dynamic type.
    TRTdynamics<T, Descriptor> *clone() const override;

    /// Return a unique ID for this class.
    int getId() const override;
    /* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_) override;

    /// Implementation of the collision step, with imposed macroscopic variables
    void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat) override;

    /// Compute equilibrium distribution function
    T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
        T thetaBar = T()) const override;

private:
    static int id;
};

template <typename T, template <typename U> class Descriptor>
class Ma1TRTdynamics : public BaseTRTdynamics<T, Descriptor> {
public:
    // inherit constructors
    using BaseTRTdynamics<T, Descriptor>::BaseTRTdynamics;
    /// Clone the object on its dynamic type.
    Ma1TRTdynamics<T, Descriptor> *clone() const override;
    ;

    /// Return a unique ID for this class.
    int getId() const override;
    /* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_) override;

    /// Implementation of the collision step, with imposed macroscopic variables
    void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat) override;

    /// Compute equilibrium distribution function
    T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
        T thetaBar = T()) const override;

private:
    static int id;
};

/// Implementation of incompressible TRT dynamics.
/** This is the TRT equivalent of IncBGKdynamics: the "rho" moment of the
 *  populations appears only as a pressure term in the equilibrium, while
 *  the other terms are multiplied by the constant rho0.
 *  breaking change omegaMinus != 1.1, now follows the baseClass default
 **/
template <typename T, template <typename U> class Descriptor>
class IncTRTdynamics : public BaseTRTdynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    // inherit constructors
    using BaseTRTdynamics<T, Descriptor>::BaseTRTdynamics;

    /// Clone the object on its dynamic type.
    IncTRTdynamics<T, Descriptor> *clone() const override;

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

private:
    static int id;
};

}  // namespace plb

#endif  // TRT_DYNAMICS_H
