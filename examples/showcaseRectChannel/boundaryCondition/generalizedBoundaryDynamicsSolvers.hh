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
 * can be instantiated -- generic implementation.
 */
#ifndef GENERALIZED_BOUNDARY_DYNAMICS_SOLVER_HH
#define GENERALIZED_BOUNDARY_DYNAMICS_SOLVER_HH

#include <Eigen3/Cholesky>
#include <Eigen3/LU>

#include "core/cell.h"
#include "generalizedBoundaryDynamicsSolvers.h"
#include "generalizedIncompressibleBoundaryTemplates.h"
#include "latticeBoltzmann/indexTemplates.h"

namespace plb {

// ================ GeneralizedBoundarySolver base class =================== //
template <typename T, template <typename U> class Descriptor>
GeneralizedBoundarySolver<T, Descriptor>::GeneralizedBoundarySolver(
    const std::vector<plint> &mInd_, const std::vector<plint> &kInd_) :
    mInd(mInd_), kInd(kInd_)
{ }

// ============ GeneralizedLinearBoundarySolver base class ================= //
template <typename T, template <typename U> class Descriptor>
GeneralizedLinearBoundarySolver<T, Descriptor>::GeneralizedLinearBoundarySolver(
    const std::vector<plint> &mInd_, const std::vector<plint> &kInd_) :
    GeneralizedBoundarySolver<T, Descriptor>(mInd_, kInd_)
{ }

template <typename T, template <typename U> class Descriptor>
void GeneralizedLinearBoundarySolver<T, Descriptor>::apply(
    Cell<T, Descriptor> &cell, const Dynamics<T, Descriptor> &dyn, bool replaceAll)
{
    Eigen::MatrixXd A;     // lhs matrix
    Eigen::VectorXd b, x;  // rhs of system of equations, unkown of the system

    createLinearSystem(cell, A, b);
    solveLinearSytem(cell, A, b, x);
    regularizePopulations(cell, x, dyn, replaceAll);
}

template <typename T, template <typename U> class Descriptor>
void GeneralizedLinearBoundarySolver<T, Descriptor>::solveLinearSytem(
    Cell<T, Descriptor> &cell, Eigen::MatrixXd &A, Eigen::VectorXd &b, Eigen::VectorXd &x)
{
    Eigen::MatrixXd AT = A.transpose();
    A = AT * A;
    b = AT * b;

#ifdef PLB_DEBUG
    //     bool solutionExists = A.lu().solve(b,&x);   // using a LU factorization
    //     pcout << "Error in the solution = " << (A*x-b).norm()/x.norm() << std::endl;
    //     PLB_ASSERT(solutionExists);

    x = A.ldlt().solve(b);
    T relError = (A * x - b).norm() / b.norm();
    PLB_ASSERT(relError < 1.0e-12);

#else
    x = A.ldlt().solve(b);
    //     A.lu().solve(b,&x);
#endif
}

// =========== GeneralizedNonLinearBoundarySolver base class =============== //

template <typename T, template <typename U> class Descriptor>
GeneralizedNonLinearBoundarySolver<T, Descriptor>::GeneralizedNonLinearBoundarySolver(
    const std::vector<plint> &mInd_, const std::vector<plint> &kInd_, T epsilon_) :
    GeneralizedBoundarySolver<T, Descriptor>(mInd_, kInd_), epsilon(epsilon_)
{ }

template <typename T, template <typename U> class Descriptor>
bool GeneralizedNonLinearBoundarySolver<T, Descriptor>::converge(
    const Eigen::VectorXd &x, const Eigen::VectorXd &dx)
{
    for (plint iPi = 0; iPi < x.rows(); ++iPi) {
        T res = (std::fabs(x[iPi]) > 1.0e-14 ? std::fabs(dx(iPi) / x(iPi)) : std::fabs(x(iPi)));

        if (res > epsilon)
            return false;
    }
    return true;
}

template <typename T, template <typename U> class Descriptor>
void GeneralizedNonLinearBoundarySolver<T, Descriptor>::iterateNonLinearSystem(
    const Eigen::MatrixXd &Jac, const Eigen::VectorXd &f, Eigen::VectorXd &x, Eigen::VectorXd &dx)
{
    Eigen::MatrixXd JacT = Jac.transpose();
    Eigen::MatrixXd JacSqr = JacT * Jac;
    Eigen::VectorXd JacTf = JacT * f;

#ifdef PLB_DEBUG
    //     bool solutionExists = JacSqr.lu().solve(JacTf,&dx);   // using a LU factorization
    //     PLB_ASSERT(solutionExists);
    x = JacSqr.ldlt().solve(JacT);
    T relError = (JacSqr * x - JacT).norm() / JacT.norm();
    PLB_ASSERT(relError < 1.0e-12);
#else
    x = JacSqr.ldlt().solve(JacT);
    //     JacSqr.lu().solve(JacTf,&dx);
#endif

    T stepMult = (T)1;   // step size (step mult can only be <= 1 (usually = 1).
    x += stepMult * dx;  // increment solution
}

template <typename T, template <typename U> class Descriptor>
void GeneralizedNonLinearBoundarySolver<T, Descriptor>::apply(
    Cell<T, Descriptor> &cell, const Dynamics<T, Descriptor> &dyn, bool replaceAll)
{
    Eigen::MatrixXd Jacobian;  // lhs matrix
    Eigen::VectorXd f, x, dx;  // rhs of system of equations, unkown of the system

    plint maxT = 10000;
    iniSystem(Jacobian, f, x, dx);
    fromMacroToX(x);
    for (plint iT = 0; iT < maxT; ++iT) {
        createNonLinearSystem(cell, Jacobian, f);
        iterateNonLinearSystem(Jacobian, f, x, dx);
        fromXtoMacro(x);
        if (converge(x, dx)) {
            break;
        }
    }

    regularizePopulations(cell, x, dyn, replaceAll);
}

// ========================================================================= //
// ============ Different kind of Linear BCs for incompressible fluids ===== //
// ========================================================================= //

// ========================= Dirichlet Velocity BC ========================= //
// ======================= without mass conservation ======================= //

template <typename T, template <typename U> class Descriptor>
DirichletVelocityBoundarySolver<T, Descriptor>::DirichletVelocityBoundarySolver(
    const std::vector<plint> &mInd_, const std::vector<plint> &kInd_,
    const Array<T, Descriptor<T>::d> &u_) :
    GeneralizedLinearBoundarySolver<T, Descriptor>(mInd_, kInd_), u(u_)
{
    sysX = SymmetricTensor<T, Descriptor>::n + 1;
    sysY = this->kInd.size() + 1;
}

template <typename T, template <typename U> class Descriptor>
void DirichletVelocityBoundarySolver<T, Descriptor>::createLinearSystem(
    Cell<T, Descriptor> &cell, Eigen::MatrixXd &A, Eigen::VectorXd &b)
{
    // matrix of the system Ax=b
    A = Eigen::MatrixXd::Zero(sysY, sysX);
    // rhs of the equation Ax=b
    b = Eigen::VectorXd::Zero(sysY);

    T uSqr = VectorTemplate<T, Descriptor>::normSqr(u);
    // f^k = A * x
    // A = g, 1/(2c_s^4) H^2
    // with g being feq/rho and H^2 the second order Hermite polynomial
    generalizedIncomprBoundaryTemplates<T, Descriptor>::f_to_A_ma2_contrib(this->kInd, u, uSqr, A);
    generalizedIncomprBoundaryTemplates<T, Descriptor>::f_to_b_contrib(cell, this->kInd, b);

    T rhoTmp = T();
    for (pluint fInd = 0; fInd < this->kInd.size(); ++fInd) {
        plint iPop = this->kInd[fInd];
        rhoTmp += fullF<T, Descriptor>(cell[iPop], iPop);
    }
    // rhoTtmp = sum_i->known f_i.
    b[sysY - 1] = rhoTmp;

    // first row of the A matrix. imposing sum_i f_i = rho.
    Eigen::RowVectorXd e0 = Eigen::RowVectorXd::Zero(sysX);
    e0[0] = 1.0;

    Eigen::RowVectorXd sumA = Eigen::RowVectorXd::Zero(sysX);
    for (pluint fInd = 0; fInd < this->mInd.size(); ++fInd) {
        plint iPop = this->mInd[fInd];
        Eigen::RowVectorXd lineA = Eigen::RowVectorXd::Zero(sysX);
        generalizedIncomprBoundaryTemplates<T, Descriptor>::f_ma2_linear(iPop, u, uSqr, lineA);
        for (plint iVec = 0; iVec < sysX; ++iVec)
            sumA(iVec) += lineA(iVec);
    }

    A.row(sysY - 1) = e0 - sumA;
}

template <typename T, template <typename U> class Descriptor>
void DirichletVelocityBoundarySolver<T, Descriptor>::regularizePopulations(
    Cell<T, Descriptor> &cell, const Eigen::VectorXd &x, const Dynamics<T, Descriptor> &dyn,
    bool replaceAll)
{
    T rho;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    generalizedIncomprBoundaryTemplates<T, Descriptor>::fromXtoRhoAndPiNeq(x, rho, PiNeq);
    T rhoBar = Descriptor<T>::rhoBar(rho);
    Array<T, Descriptor<T>::d> j = rho * u;
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    if (replaceAll) {
        dyn.regularize(cell, rhoBar, j, jSqr, PiNeq);
    } else {
        Cell<T, Descriptor> tmpCell;
        dyn.regularize(tmpCell, rhoBar, j, jSqr, PiNeq);
        for (pluint fInd = 0; fInd < this->mInd.size(); ++fInd) {
            plint iPop = this->mInd[fInd];
            cell[iPop] = tmpCell[iPop];
        }
    }
}

// ========================= Dirichlet Velocity BC ========================= //
// ======================= with  mass conservation ======================= //

template <typename T, template <typename U> class Descriptor>
DirichletMassConservingVelocityBoundarySolver<T, Descriptor>::
    DirichletMassConservingVelocityBoundarySolver(
        const std::vector<plint> &mInd_, const std::vector<plint> &kInd_,
        const std::vector<plint> &inGoingInd_, const Array<T, Descriptor<T>::d> &u_) :
    GeneralizedLinearBoundarySolver<T, Descriptor>(mInd_, kInd_), inGoingInd(inGoingInd_), u(u_)
{
    sysX = SymmetricTensor<T, Descriptor>::n + 1;
    sysY = this->kInd.size() + 1;
}

template <typename T, template <typename U> class Descriptor>
void DirichletMassConservingVelocityBoundarySolver<T, Descriptor>::createLinearSystem(
    Cell<T, Descriptor> &cell, Eigen::MatrixXd &A, Eigen::VectorXd &b)
{
    // matrix of the system Ax=b
    A = Eigen::MatrixXd::Zero(sysY, sysX);
    // rhs of the equation Ax=b
    b = Eigen::VectorXd::Zero(sysY);

    T uSqr = VectorTemplate<T, Descriptor>::normSqr(u);
    // f^k = A * x
    // A = g, 1/(2c_s^4) H^2
    // with g being feq/rho and H^2 the second order Hermite polynomial

    generalizedIncomprBoundaryTemplates<T, Descriptor>::f_to_A_ma2_contrib(this->kInd, u, uSqr, A);
    generalizedIncomprBoundaryTemplates<T, Descriptor>::f_to_b_contrib(cell, this->kInd, b);

    // computing mass incoming in the wall
    T rhoTmp = T();
    for (pluint fInd = 0; fInd < this->inGoingInd.size(); ++fInd) {
        plint iPop = indexTemplates::opposite<Descriptor<T> >(inGoingInd[fInd]);
        rhoTmp += fullF<T, Descriptor>(cell[iPop], iPop);
    }
    // rhoTtmp = sum_i->in_wall f_i.
    b[sysY - 1] = rhoTmp;

    const T omega = cell.getDynamics().getOmega();
    Eigen::RowVectorXd sumA = Eigen::RowVectorXd::Zero(sysX);
    for (pluint fInd = 0; fInd < inGoingInd.size(); ++fInd) {
        plint iPop = inGoingInd[fInd];
        Eigen::RowVectorXd lineA = Eigen::RowVectorXd::Zero(sysX);
        generalizedIncomprBoundaryTemplates<T, Descriptor>::f_ma2_linear(
            iPop, u, uSqr, lineA, omega);
        for (plint iVec = 0; iVec < sysX; ++iVec)
            sumA(iVec) += lineA(iVec);
    }

    A.row(sysY - 1) = sumA;
}

template <typename T, template <typename U> class Descriptor>
void DirichletMassConservingVelocityBoundarySolver<T, Descriptor>::regularizePopulations(
    Cell<T, Descriptor> &cell, const Eigen::VectorXd &x, const Dynamics<T, Descriptor> &dyn,
    bool replaceAll)
{
    T rho;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    generalizedIncomprBoundaryTemplates<T, Descriptor>::fromXtoRhoAndPiNeq(x, rho, PiNeq);
    T rhoBar = Descriptor<T>::rhoBar(rho);
    Array<T, Descriptor<T>::d> j = rho * u;
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    dyn.regularize(cell, rhoBar, j, jSqr, PiNeq);
}

// ========================================================================= //
// ======== Different kind of Non Linear BCs for incompressible fluids ===== //
// ========================================================================= //

// ========================= Dirichlet Density BC ========================= //
template <typename T, template <typename U> class Descriptor, int dir>
DirichletDensityBoundarySolver<T, Descriptor, dir>::DirichletDensityBoundarySolver(
    const std::vector<plint> &mInd_, const std::vector<plint> &kInd_, T &rho_,
    Array<T, Descriptor<T>::d> &u_, Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq_,
    T epsilon_) :
    GeneralizedNonLinearBoundarySolver<T, Descriptor>(mInd_, kInd_, epsilon_)
{
    sysX = SymmetricTensor<T, Descriptor>::n + 1;
    sysY = this->kInd.size();

    this->macro.push_back(rho_);
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
        this->macro.push_back(u_[iD]);
    for (plint iPi = 0; iPi < SymmetricTensor<T, Descriptor>::n; ++iPi)
        this->macro.push_back(PiNeq_[iPi]);
}

template <typename T, template <typename U> class Descriptor, int dir>
void DirichletDensityBoundarySolver<T, Descriptor, dir>::fromMacroToX(Eigen::VectorXd &x)
{
    x[0] = this->macro[1 + dir];
    for (plint iPi = 0; iPi < SymmetricTensor<T, Descriptor>::n; ++iPi) {
        x[iPi + 1] = this->macro[1 + Descriptor<T>::d + iPi];
    }
}

template <typename T, template <typename U> class Descriptor, int dir>
void DirichletDensityBoundarySolver<T, Descriptor, dir>::fromXtoMacro(const Eigen::VectorXd &x)
{
    this->macro[1 + dir] = x[0];
    for (plint iPi = 0; iPi < SymmetricTensor<T, Descriptor>::n; ++iPi) {
        this->macro[1 + Descriptor<T>::d + iPi] = x[iPi + 1];
    }
}

template <typename T, template <typename U> class Descriptor, int dir>
void DirichletDensityBoundarySolver<T, Descriptor, dir>::iniSystem(
    Eigen::MatrixXd &Jac, Eigen::VectorXd &f, Eigen::VectorXd &x, Eigen::VectorXd &dx)
{
    Jac = Eigen::MatrixXd::Zero(sysY, sysX);
    f = Eigen::VectorXd::Zero(sysY);  // stores the non-linear function
    x = Eigen::VectorXd::Zero(sysX);  // stores the variables (u[dir] and PiNeq)
    dx = Eigen::VectorXd::Zero(
        sysX);  // contains delta_u[dir], delta_PiNeq (the increments towards the solution)
}

template <typename T, template <typename U> class Descriptor, int dir>
void DirichletDensityBoundarySolver<T, Descriptor, dir>::createNonLinearSystem(
    const Cell<T, Descriptor> &cell, Eigen::MatrixXd &Jac, Eigen::VectorXd &f)
{
    // computes the non linear function F(macro) = f_i-(f^eq+f^1)
    T rho;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    generalizedIncomprBoundaryTemplates<T, Descriptor>::from_macro_to_rho_j_pineq(
        this->macro, rho, j, PiNeq);

    T rhoBar = Descriptor<T>::rhoBar(rho);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    for (pluint iPop = 0; iPop < this->kInd.size(); ++iPop) {
        f[iPop] =
            cell[this->kInd[iPop]]
            - (dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(
                   this->kInd[iPop], rhoBar, invRho, j, jSqr)
               + offEquilibriumTemplates<T, Descriptor>::fromPiToFneq(this->kInd[iPop], PiNeq));
    }

    // computes the Jacobian of F(macro)
    Eigen::RowVectorXd df = Eigen::RowVectorXd::Zero(sysX);
    for (pluint iPop = 0; iPop < this->kInd.size(); ++iPop) {
        generalizedIncomprBoundaryTemplates<T, Descriptor>::compute_f_diff_u_dir_and_PiNeq(
            this->kInd[iPop], rho, j * invRho, df, dir);
        Jac.row(iPop) = df;
    }
}

template <typename T, template <typename U> class Descriptor, int dir>
void DirichletDensityBoundarySolver<T, Descriptor, dir>::regularizePopulations(
    Cell<T, Descriptor> &cell, const Eigen::VectorXd &x, const Dynamics<T, Descriptor> &dyn,
    bool replaceAll)
{
    T rho;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    generalizedIncomprBoundaryTemplates<T, Descriptor>::from_macro_to_rho_j_pineq(
        this->macro, rho, j, PiNeq);

    T rhoBar = Descriptor<T>::rhoBar(rho);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    if (replaceAll) {
        dyn.regularize(cell, rhoBar, j, jSqr, PiNeq);
    } else {
        Cell<T, Descriptor> tmpCell;
        dyn.regularize(tmpCell, rhoBar, j, jSqr, PiNeq);
        for (pluint fInd = 0; fInd < this->mInd.size(); ++fInd) {
            plint iPop = this->mInd[fInd];
            cell[iPop] = tmpCell[iPop];
        }
    }
}

}  // namespace plb

#endif  // GENERALIZED_BOUNDARY_DYNAMICS_HH
