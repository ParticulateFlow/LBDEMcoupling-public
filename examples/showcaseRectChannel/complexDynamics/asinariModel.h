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
 * Implementation of the LW-ACM model by Pietro Asinari and others -- header file.
 * More details in http://arxiv.org/abs/1111.2142
 */
#ifndef ASINARI_MODEL_H
#define ASINARI_MODEL_H

#include "atomicBlock/dataProcessingFunctional2D.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "basicDynamics/isoThermalDynamics.h"

namespace plb {

/// First part of collision in Asinari's LW-ACM model.
/** Here, density and momentum are computed from the populations
 *  and stored in external scalars. Then, the first part of the
 *  collision is executed: equation 3a in http://arxiv.org/abs/1111.2142 .
 *  You must combine this with an execution of the data processor
 *  AsinariPostCollide3D, right after the streaming step.
 */
template <typename T, template <typename U> class Descriptor>
class AsinariDynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    AsinariDynamics(T omega_);
    AsinariDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual AsinariDynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Set local relaxation parameter of the dynamics.
    virtual void setOmega(T omega_);

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
    /// Get local value of any generic parameter.
    /** For Asinari-Dynamics parameter with ID 1000 is the prefactor
     *  2*(omega-1)/omega.
     **/
    virtual T getParameter(plint whichParameter) const;

private:
    void computePrefactor();

private:
    static int id;
    T prefactor;
};

/// First part of collision in Asinari's LW-ACM model, with incompressible equilibrium.
template <typename T, template <typename U> class Descriptor>
class IncAsinariDynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    IncAsinariDynamics(T omega_);
    IncAsinariDynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual IncAsinariDynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Set local relaxation parameter of the dynamics.
    virtual void setOmega(T omega_);

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
    /// Get local value of any generic parameter.
    /** For Asinari-Dynamics parameter with ID 1000 is the prefactor
     *  2*(omega-1)/omega.
     **/
    virtual T getParameter(plint whichParameter) const;

private:
    void computePrefactor();

private:
    static int id;
    T prefactor;
};

/* ************* Class AsinariPostCollide3D ******************* */

/// Second part of collision in Asinari's LW-ACM model, for 3D.
/** This implements equation 3c in http://arxiv.org/abs/1111.2142 .
 */
template <typename T, template <typename U> class Descriptor>
class AsinariPostCollide3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual AsinariPostCollide3D<T, Descriptor> *clone() const
    {
        return new AsinariPostCollide3D(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
    }
};

/// Second part of collision in Asinari's LW-ACM model, for 2D.
/** This implements equation 3c in http://arxiv.org/abs/1111.2142 .
 */
template <typename T, template <typename U> class Descriptor>
class AsinariPostCollide2D : public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual AsinariPostCollide2D<T, Descriptor> *clone() const
    {
        return new AsinariPostCollide2D(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
    }
};

}  // namespace plb

#endif  // ASINARI_MODEL_H
