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
 * A collection of dynamics classes (collision models) with which a Cell object
 * can be instantiated -- header file.
 * Theoretical background about these collision models can be found in
 * Coreixas et al. 'Comprehensive comparison of collision models in the
 * lattice Boltzmann framework: Theoretical investigations', PRE, 2019.
 */
#ifndef COMPREHENSIVE_ISO_THERMAL_DYNAMICS_H
#define COMPREHENSIVE_ISO_THERMAL_DYNAMICS_H

#include "basicDynamics/isoThermalDynamics.h"
#include "core/dynamics.h"
#include "core/globalDefs.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class RMdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    RMdynamics(T omegaVisc);
    RMdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual RMdynamics<T, Descriptor> *clone() const;

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
    static Array<T, Descriptor<T>::numRelaxationTimes> allOmega;

private:
    virtual void decomposeOrder0(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void decomposeOrder1(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void recomposeOrder0(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;
    virtual void recomposeOrder1(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;

private:
    static int id;
};

template <typename T, template <typename U> class Descriptor>
class HMdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    HMdynamics(T omegaVisc);
    HMdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual HMdynamics<T, Descriptor> *clone() const;

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

    static Array<T, Descriptor<T>::numRelaxationTimes> allOmega;

private:
    virtual void decomposeOrder0(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void decomposeOrder1(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void recomposeOrder0(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;
    virtual void recomposeOrder1(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;

private:
    static int id;
};

template <typename T, template <typename U> class Descriptor>
class CMdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    CMdynamics(T omegaVisc);
    CMdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual CMdynamics<T, Descriptor> *clone() const;

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

    static Array<T, Descriptor<T>::numRelaxationTimes> allOmega;

private:
    virtual void decomposeOrder0(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void decomposeOrder1(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void recomposeOrder0(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;
    virtual void recomposeOrder1(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;

private:
    static int id;
};

template <typename T, template <typename U> class Descriptor>
class CHMdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    CHMdynamics(T omegaVisc);
    CHMdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual CHMdynamics<T, Descriptor> *clone() const;

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

    static Array<T, Descriptor<T>::numRelaxationTimes> allOmega;

private:
    virtual void decomposeOrder0(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void decomposeOrder1(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void recomposeOrder0(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;
    virtual void recomposeOrder1(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;

private:
    static int id;
};

template <typename T, template <typename U> class Descriptor>
class Kdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    Kdynamics(T omegaVisc);
    Kdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual Kdynamics<T, Descriptor> *clone() const;

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

    static Array<T, Descriptor<T>::numRelaxationTimes> allOmega;

private:
    virtual void decomposeOrder0(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void decomposeOrder1(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void recomposeOrder0(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;
    virtual void recomposeOrder1(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;

private:
    static int id;
};

template <typename T, template <typename U> class Descriptor>
class GHdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    GHdynamics(T omegaVisc);
    GHdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual GHdynamics<T, Descriptor> *clone() const;

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

    static Array<T, Descriptor<T>::numRelaxationTimes> allOmega;

private:
    virtual void decomposeOrder0(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void decomposeOrder1(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void recomposeOrder0(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;
    virtual void recomposeOrder1(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;

private:
    static int id;
};

template <typename T, template <typename U> class Descriptor>
class RRdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    RRdynamics(T omegaVisc);
    RRdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual RRdynamics<T, Descriptor> *clone() const;

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

    static Array<T, Descriptor<T>::numRelaxationTimes> allOmega;

private:
    virtual void decomposeOrder0(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void decomposeOrder1(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void recomposeOrder0(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;
    virtual void recomposeOrder1(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;

private:
    static int id;
};

}  // namespace plb

#endif  // COMPREHENSIVE_ISO_THERMAL_DYNAMICS_H
