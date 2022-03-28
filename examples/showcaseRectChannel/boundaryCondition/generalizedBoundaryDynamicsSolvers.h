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
#ifndef GENERALIZED_BOUNDARY_DYNAMICS_SOLVERS_H
#define GENERALIZED_BOUNDARY_DYNAMICS_SOLVERS_H

#include <Eigen3/Core>

#include "core/cell.h"
#include "core/globalDefs.h"

namespace plb {

// ================ GeneralizedBoundarySolver base class =================== //
template <typename T, template <typename U> class Descriptor>
class GeneralizedBoundarySolver {
public:
    GeneralizedBoundarySolver(const std::vector<plint> &mInd_, const std::vector<plint> &kInd_);
    virtual ~GeneralizedBoundarySolver() {};

    virtual void apply(
        Cell<T, Descriptor> &cell, const Dynamics<T, Descriptor> &dyn, bool replaceAll = true) = 0;

protected:
    // missing idexes (mInd) and known indexes (kInnd)
    std::vector<plint> mInd, kInd;
};

// ============== GeneralizedLinearBoundarySolver base class =============== //
template <typename T, template <typename U> class Descriptor>
class GeneralizedLinearBoundarySolver : public GeneralizedBoundarySolver<T, Descriptor> {
public:
    GeneralizedLinearBoundarySolver(
        const std::vector<plint> &mInd_, const std::vector<plint> &kInd_);

    virtual void apply(
        Cell<T, Descriptor> &cell, const Dynamics<T, Descriptor> &dyn, bool replaceAll = true);

    virtual void createLinearSystem(
        Cell<T, Descriptor> &cell, Eigen::MatrixXd &A, Eigen::VectorXd &b) = 0;
    virtual void regularizePopulations(
        Cell<T, Descriptor> &cell, const Eigen::VectorXd &x, const Dynamics<T, Descriptor> &dyn,
        bool replaceAll) = 0;

    void solveLinearSytem(
        Cell<T, Descriptor> &cell, Eigen::MatrixXd &A, Eigen::VectorXd &b, Eigen::VectorXd &x);
};

// ============ GeneralizedNonLinearBoundarySolver base class ============== //
template <typename T, template <typename U> class Descriptor>
class GeneralizedNonLinearBoundarySolver : public GeneralizedBoundarySolver<T, Descriptor> {
public:
    GeneralizedNonLinearBoundarySolver(
        const std::vector<plint> &mInd_, const std::vector<plint> &kInd_, T epsilon_);

    virtual void apply(
        Cell<T, Descriptor> &cell, const Dynamics<T, Descriptor> &dyn, bool replaceAll = true);

    bool converge(const Eigen::VectorXd &x, const Eigen::VectorXd &dx);

    virtual void fromXtoMacro(const Eigen::VectorXd &x) = 0;
    virtual void fromMacroToX(Eigen::VectorXd &x) = 0;

    virtual void iniSystem(
        Eigen::MatrixXd &Jac, Eigen::VectorXd &f, Eigen::VectorXd &x, Eigen::VectorXd &dx) = 0;
    virtual void createNonLinearSystem(
        const Cell<T, Descriptor> &cell, Eigen::MatrixXd &Jac, Eigen::VectorXd &f) = 0;
    virtual void iterateNonLinearSystem(
        const Eigen::MatrixXd &Jac, const Eigen::VectorXd &f, Eigen::VectorXd &x,
        Eigen::VectorXd &dx);
    virtual void regularizePopulations(
        Cell<T, Descriptor> &cell, const Eigen::VectorXd &x, const Dynamics<T, Descriptor> &dyn,
        bool replaceAll) = 0;

private:
    T epsilon;

protected:
    std::vector<T> macro;
};

// ========================================================================= //
// ============ Different kind of Linear BCs for incompressible fluids ===== //
// ========================================================================= //

// ========================= Dirichlet Velocity BC ========================= //
// ======================= without mass conservation ======================= //
template <typename T, template <typename U> class Descriptor>
class DirichletVelocityBoundarySolver : public GeneralizedLinearBoundarySolver<T, Descriptor> {
public:
    DirichletVelocityBoundarySolver(
        const std::vector<plint> &mInd_, const std::vector<plint> &kInd_,
        const Array<T, Descriptor<T>::d> &u_);

    virtual void createLinearSystem(
        Cell<T, Descriptor> &cell, Eigen::MatrixXd &A, Eigen::VectorXd &b);
    virtual void regularizePopulations(
        Cell<T, Descriptor> &cell, const Eigen::VectorXd &x, const Dynamics<T, Descriptor> &dyn,
        bool replaceAll);

private:
    Array<T, Descriptor<T>::d> u;
    T uSqr;
    plint sysX, sysY;
};

// ========================= Dirichlet Velocity BC ========================= //
// ======================= with mass conservation ======================= //
template <typename T, template <typename U> class Descriptor>
class DirichletMassConservingVelocityBoundarySolver :
    public GeneralizedLinearBoundarySolver<T, Descriptor> {
public:
    DirichletMassConservingVelocityBoundarySolver(
        const std::vector<plint> &mInd_, const std::vector<plint> &kInd_,
        const std::vector<plint> &inGoingIndices_, const Array<T, Descriptor<T>::d> &u_);

    virtual void createLinearSystem(
        Cell<T, Descriptor> &cell, Eigen::MatrixXd &A, Eigen::VectorXd &b);
    virtual void regularizePopulations(
        Cell<T, Descriptor> &cell, const Eigen::VectorXd &x, const Dynamics<T, Descriptor> &dyn,
        bool replaceAll);

private:
    std::vector<plint> inGoingInd;
    Array<T, Descriptor<T>::d> u;
    T uSqr;
    plint sysX, sysY;
};

// ========================================================================= //
// ======== Different kind of Non-Linear BCs for incompressible fluids ===== //
// ========================================================================= //

// ========================= Dirichlet Density BC ========================== //
// ======================= without mass conservation ======================= //
template <typename T, template <typename U> class Descriptor, int dir>
class DirichletDensityBoundarySolver : public GeneralizedNonLinearBoundarySolver<T, Descriptor> {
public:
    DirichletDensityBoundarySolver(
        const std::vector<plint> &mInd_, const std::vector<plint> &kInd_, T &rho_,
        Array<T, Descriptor<T>::d> &u_, Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq_,
        T epsilon_);

    virtual void fromXtoMacro(const Eigen::VectorXd &x);
    virtual void fromMacroToX(Eigen::VectorXd &x);
    virtual void iniSystem(
        Eigen::MatrixXd &Jac, Eigen::VectorXd &f, Eigen::VectorXd &x, Eigen::VectorXd &dx);
    virtual void createNonLinearSystem(
        const Cell<T, Descriptor> &cell, Eigen::MatrixXd &Jac, Eigen::VectorXd &f);
    virtual void regularizePopulations(
        Cell<T, Descriptor> &cell, const Eigen::VectorXd &x, const Dynamics<T, Descriptor> &dyn,
        bool replaceAll);

private:
    plint sysX, sysY;
};

}  // namespace plb

#endif  // GENERALIZED_BOUNDARY_DYNAMICS_H
