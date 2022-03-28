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
#ifndef GENERALIZED_INCOMPRESSIBLE_BOUNDARY_TEMPLATES_H
#define GENERALIZED_INCOMPRESSIBLE_BOUNDARY_TEMPLATES_H

#include <Eigen3/Cholesky>
#include <Eigen3/Core>
#include <Eigen3/Dense>
#include <Eigen3/LU>
#include <Eigen3/QR>
#include <Eigen3/SVD>

#include "core/cell.h"
#include "core/dynamicsIdentifiers.h"
#include "generalizedBoundaryDynamics.h"
#include "latticeBoltzmann/hermitePolynomialsTemplates.h"
#include "latticeBoltzmann/indexTemplates.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
struct generalizedIncomprBoundaryTemplates {
    static T equilibrium_ma2_over_rho(plint iPop, const Array<T, Descriptor<T>::d> &u, T uSqr)
    {
        T c_u = Descriptor<T>::c[iPop][0] * u[0];
        for (int iD = 1; iD < Descriptor<T>::d; ++iD) {
            c_u += Descriptor<T>::c[iPop][iD] * u[iD];
        }
        return Descriptor<T>::t[iPop]
               * ((T)1 + Descriptor<T>::invCs2 * c_u
                  + Descriptor<T>::invCs2 / (T)2 * (Descriptor<T>::invCs2 * c_u * c_u - uSqr));
    }

    // f = w_i*rho*g_i+H2/(2*cs^4):PiNeq (rho and PiNeq unknowns)
    static void f_ma2_linear(
        plint iPop, const Array<T, Descriptor<T>::d> &u, T uSqr, Eigen::RowVectorXd &a)
    {
        T eqOverRho = equilibrium_ma2_over_rho(iPop, u, uSqr);
        a[0] = eqOverRho;

        T factor = 0.5 * Descriptor<T>::t[iPop] * Descriptor<T>::invCs2 * Descriptor<T>::invCs2;
        Array<T, SymmetricTensor<T, Descriptor>::n> H2 =
            HermiteTemplate<T, Descriptor>::contractedOrder2(iPop);

        for (plint iPi = 1; iPi <= SymmetricTensor<T, Descriptor>::n; ++iPi)
            a[iPi] = H2[iPi - 1] * factor;
    }

    static void f_ma2_linear(
        plint iPop, const Array<T, Descriptor<T>::d> &u, T uSqr, Eigen::RowVectorXd &a, T omega)
    {
        T eqOverRho = equilibrium_ma2_over_rho(iPop, u, uSqr);
        a(0) = eqOverRho;

        T factor = 0.5 * Descriptor<T>::t[iPop] * Descriptor<T>::invCs2 * Descriptor<T>::invCs2
                   * ((T)1 - omega);
        Array<T, SymmetricTensor<T, Descriptor>::n> H2 =
            HermiteTemplate<T, Descriptor>::contractedOrder2(iPop);

        for (plint iPi = 1; iPi <= SymmetricTensor<T, Descriptor>::n; ++iPi)
            a(iPi) = H2[iPi - 1] * factor;
    }

    static void f_to_A_ma2_contrib(
        const std::vector<plint> kInd, const Array<T, Descriptor<T>::d> &u, T uSqr,
        Eigen::MatrixXd &A)
    {
        for (pluint fInd = 0; fInd < kInd.size(); ++fInd) {
            plint iPop = kInd[fInd];
            Eigen::RowVectorXd lineA = Eigen::RowVectorXd::Zero(A.cols());
            generalizedIncomprBoundaryTemplates<T, Descriptor>::f_ma2_linear(iPop, u, uSqr, lineA);
            A.row(fInd) = lineA;
        }
    }

    static void f_to_b_contrib(
        Cell<T, Descriptor> &cell, const std::vector<plint> kInd, Eigen::VectorXd &b)
    {
        for (pluint fInd = 0; fInd < kInd.size(); ++fInd) {
            plint iPop = kInd[fInd];
            b[fInd] = fullF<T, Descriptor>(cell[iPop], iPop);
        }
    }

    static void from_macro_to_rho_j_pineq(
        const std::vector<T> &macro, T &rho, Array<T, Descriptor<T>::d> &j,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq)
    {
        rho = macro[0];
        for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
            j[iD] = rho * (macro[1 + iD]);
        for (plint iPi = 0; iPi < SymmetricTensor<T, Descriptor>::n; ++iPi)
            PiNeq[iPi] = macro[1 + Descriptor<T>::d + iPi];
    }

    // generic tranformation methods
    static void fromRhoAndPiNeqToX(
        T rho, const Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq, Eigen::VectorXd &x)
    {
        x(0) = rho;
        for (plint iPi = 0; iPi < SymmetricTensor<T, Descriptor>::n; ++iPi)
            x(iPi + 1) = PiNeq[iPi];
    }

    static void fromXtoRhoAndPiNeq(
        const Eigen::VectorXd &x, T &rho, Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq)
    {
        rho = x(0);
        for (plint iPi = 0; iPi < SymmetricTensor<T, Descriptor>::n; ++iPi)
            PiNeq[iPi] = x(iPi + 1);
    }

    static void fromUandPiNeqToX(
        const Array<T, Descriptor<T>::d> &u,
        const Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq, Eigen::VectorXd &x, plint dir)
    {
        x(0) = u[dir];
        for (plint iPi = 0; iPi < SymmetricTensor<T, Descriptor>::n; ++iPi)
            x(iPi + 1) = PiNeq[iPi];
    }

    static void fromXtoUandPiNeq(
        const Eigen::VectorXd &x, Array<T, Descriptor<T>::d> &u,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq, plint dir)
    {
        u[dir] = x(0);
        for (plint iPi = 0; iPi < SymmetricTensor<T, Descriptor>::n; ++iPi)
            PiNeq[iPi] = x(iPi + 1);
    }

    // creation of the over-determined (usually) linear system
    //     static void createLinearSystem(const Cell<T,Descriptor>& cell, const
    //     Array<T,Descriptor<T>::d> &u,
    //                                    const std::vector<plint> &missingIndices,
    //                                    const std::vector<plint> &knownIndices,
    //                                    Eigen::MatrixXd &A,Eigen::VectorXd &b) {
    //
    //         T uSqr = VectorTemplate<T,Descriptor>::normSqr(u);
    //
    //         plint systSizeX = SymmetricTensor<T,Descriptor>::n+1;
    //         plint systSizeY = knownIndices.size()+1;
    //
    //         // matrix of the system Ax=b
    //         A = Eigen::MatrixXd::Zero(systSizeY,systSizeX);
    //         // rhs of the equation Ax=b
    //         b = Eigen::VectorXd::Zero(systSizeY);
    //
    //         // f^k = A * x
    //         // A = g, 1/(2c_s^4) H^2
    //         // with g being feq/rho and H^2 the second order Hermite polynomial
    //         for (pluint fInd = 0; fInd < knownIndices.size(); ++fInd) {
    //             plint iPop = knownIndices[fInd];
    //             Eigen::RowVectorXd lineA = Eigen::RowVectorXd::Zero(systSizeX);
    //             computeMatrixRow(iPop, u, uSqr, lineA);
    //             A.row(fInd) = lineA;
    //         }
    //
    //         T rhoTmp = T();
    //         for (pluint kInd = 0; kInd < knownIndices.size(); ++kInd) {
    //             plint iPop = knownIndices[kInd];
    //             rhoTmp += fullF<T,Descriptor>(cell[iPop], iPop);
    //             b(kInd) = fullF<T,Descriptor>(cell[iPop], iPop);
    //         }
    //         // rhoTtmp = sum_i->known f_i.
    //         b(knownIndices.size()) = rhoTmp;
    //
    //         // first row of the A matrix. imposing sum_i f_i = rho.
    //         Eigen::RowVectorXd e0 = Eigen::RowVectorXd::Zero(systSizeX);
    //         e0(0) = 1.0;
    //
    //         Eigen::RowVectorXd sumA = Eigen::RowVectorXd::Zero(systSizeX);
    //         for (pluint fInd = 0; fInd < missingIndices.size(); ++fInd) {
    //             plint iPop = missingIndices[fInd];
    //             Eigen::RowVectorXd lineA = Eigen::RowVectorXd::Zero(systSizeX);
    //             computeMatrixRow(iPop, u, uSqr, lineA);
    //             for (plint iVec = 0; iVec < systSizeX; ++iVec) sumA(iVec) += lineA(iVec);
    //         }
    //
    //         A.row(knownIndices.size()) = e0-sumA;
    //
    //     }
    //
    //     static void solveLinearSystemEigen(const Cell<T,Descriptor>& cell,
    //                                       const Array<T,Descriptor<T>::d> &u,
    //                                       const std::vector<plint> &missingIndices,
    //                                       const std::vector<plint> &knownIndices,
    //                                       T &rho,
    //                                       Array<T,SymmetricTensor<T,Descriptor>::n> &PiNeq) {
    //         Eigen::MatrixXd A;
    //         Eigen::VectorXd b;
    //
    //         createLinearSystem(cell, u, missingIndices, knownIndices, A, b);
    //
    //         Eigen::VectorXd x;
    //
    //         Eigen::MatrixXd AT = A.transpose();
    //         A = AT * A;
    //         b = AT * b;
    //
    //         #ifdef PLB_DEBUG
    //         bool solutionExists = A.lu().solve(b,&x);   // using a LU factorization
    //         PLB_ASSERT(solutionExists);
    //         #else
    //         A.lu().solve(b,&x);
    //         #endif
    //
    //         fromXtoRhoAndPiNeq(x,rho,PiNeq);
    //     }

    // ========= Methods used for the density BCs ============== //
    static void compute_f_diff_u_dir_and_PiNeq(
        plint iPop, T rho, const Array<T, Descriptor<T>::d> &u, Eigen::RowVectorXd &df, plint dir)
    {
        T tcs2 = Descriptor<T>::invCs2 * Descriptor<T>::t[iPop];
        T factor = 0.5 * tcs2 * Descriptor<T>::invCs2;
        Array<T, SymmetricTensor<T, Descriptor>::n> H2 =
            HermiteTemplate<T, Descriptor>::contractedOrder2(iPop);

        T diffFeqUmissing = T();
        plint iPi = 0;
        for (plint iA = 0; iA < Descriptor<T>::d; ++iA) {
            for (plint iB = iA; iB < Descriptor<T>::d; ++iB) {
                if (iA == dir || iB == dir) {
                    if (iA == iB) {
                        diffFeqUmissing += (T)2 * H2[iPi] * u[dir];
                    } else if (iA != dir) {
                        diffFeqUmissing += H2[iPi] * u[iA];
                    } else if (iB != dir) {
                        diffFeqUmissing += H2[iPi] * u[iB];
                    }
                }
                df(iPi + 1) = factor * H2[iPi];
                ++iPi;
            }
        }
        diffFeqUmissing *= factor;
        diffFeqUmissing += Descriptor<T>::c[iPop][dir] * tcs2;
        diffFeqUmissing *= rho;

        df(0) = diffFeqUmissing;
    }

    static void computeDiffF(
        plint iPop, T rho, const Array<T, Descriptor<T>::d> &u, plint dir, Eigen::RowVectorXd &df)
    {
        T tcs2 = Descriptor<T>::invCs2 * Descriptor<T>::t[iPop];
        T factor = 0.5 * tcs2 * Descriptor<T>::invCs2;
        Array<T, SymmetricTensor<T, Descriptor>::n> H2 =
            HermiteTemplate<T, Descriptor>::contractedOrder2(iPop);

        T diffFeqUmissing = T();
        plint iPi = 0;
        for (plint iA = 0; iA < Descriptor<T>::d; ++iA) {
            for (plint iB = iA; iB < Descriptor<T>::d; ++iB) {
                if (iA == dir || iB == dir) {
                    if (iA == iB) {
                        diffFeqUmissing += (T)2 * H2[iPi] * u[dir];
                    } else if (iA != dir) {
                        diffFeqUmissing += H2[iPi] * u[iA];
                    } else if (iB != dir) {
                        diffFeqUmissing += H2[iPi] * u[iB];
                    }
                }
                df(iPi + 1) = factor * H2[iPi];
                ++iPi;
            }
        }
        diffFeqUmissing *= factor;
        diffFeqUmissing += Descriptor<T>::c[iPop][dir] * tcs2;
        diffFeqUmissing *= rho;

        df(0) = diffFeqUmissing;
    }

    static void computeJacobian(
        T rho, const Array<T, Descriptor<T>::d> &u,
        const Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq,
        const plint dir,  // direction of the unknown velocity
        const std::vector<plint> &knownIndices, Eigen::MatrixXd &Jac)
    {
        plint systSizeX = SymmetricTensor<T, Descriptor>::n + 1;

        Eigen::RowVectorXd df = Eigen::RowVectorXd::Zero(systSizeX);
        for (pluint iPop = 0; iPop < knownIndices.size(); ++iPop) {
            computeDiffF(knownIndices[iPop], rho, u, dir, df);
            Jac.row(iPop) = df;
        }
    }

    static void computeNonLinearFunction(
        const Cell<T, Descriptor> &cell, T rho, const Array<T, Descriptor<T>::d> &u, const T uSqr,
        const Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq,
        const plint dir,  // direction of the unknown velocity
        const std::vector<plint> &knownIndices, Eigen::VectorXd &f)
    {
        T rhoBar = Descriptor<T>::rhoBar(rho);
        T invRho = Descriptor<T>::invRho(rhoBar);
        Array<T, Descriptor<T>::d> j = rho * u;
        T jSqr = rho * rho * uSqr;
        for (pluint iPop = 0; iPop < knownIndices.size(); ++iPop) {
            f(iPop) = cell[knownIndices[iPop]]
                      - (dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(
                             knownIndices[iPop], rhoBar, invRho, j, jSqr)
                         + offEquilibriumTemplates<T, Descriptor>::fromPiToFneq(
                             knownIndices[iPop], PiNeq));
        }
    }

    static bool converge(Eigen::VectorXd &x, Eigen::VectorXd &dx, T epsilon)
    {
        for (plint iPi = 0; iPi < x.rows(); ++iPi) {
            T res = (std::fabs(x[iPi]) > 1.0e-14 ? std::fabs(dx(iPi) / x(iPi)) : std::fabs(x(iPi)));

            if (res > epsilon)
                return false;
        }
        return true;
    }

    static void iterativelySolveSystem(
        const Cell<T, Descriptor> &cell, T rho, Array<T, Descriptor<T>::d> &u,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq,
        const int dir,  // direction of the unknown velocity
        const std::vector<plint> &knownIndices, T epsilon)
    {
        // u and PiNeq contain the initial guess for the solution of the system
        plint maxT = 10000;

        plint systSizeX = SymmetricTensor<T, Descriptor>::n + 1;
        plint systSizeY = knownIndices.size();

        Eigen::VectorXd f = Eigen::VectorXd::Zero(systSizeY);  // stores the non-linear function
        Eigen::VectorXd x =
            Eigen::VectorXd::Zero(systSizeX);  // stores the variables (u[dir] and PiNeq)
        Eigen::VectorXd dx = Eigen::VectorXd::Zero(
            systSizeX);  // contains delta_u[dir], delta_PiNeq (the increments towards the solution)

        fromUandPiNeqToX(u, PiNeq, x, dir);
        Eigen::MatrixXd Jac = Eigen::MatrixXd::Zero(systSizeY, systSizeX);
        Eigen::MatrixXd JacT = Jac.transpose();
        for (plint iT = 0; iT < maxT; ++iT) {
            T uSqr = VectorTemplate<T, Descriptor>::normSqr(u);
            computeNonLinearFunction(cell, rho, u, uSqr, PiNeq, dir, knownIndices, f);
            //             std::cout << iT << " " << f << std::endl << std::endl;
            computeJacobian(rho, u, PiNeq, dir, knownIndices, Jac);
            JacT = Jac.transpose();

            Eigen::MatrixXd JacSqr = JacT * Jac;
            Eigen::VectorXd JacTf = JacT * f;

#ifdef PLB_DEBUG
            //             bool solutionExists = JacSqr.lu().solve(JacTf,&dx);   // using a LU
            //             factorization PLB_ASSERT(solutionExists);
            dx = JacSqr.fullPivLu().solve(JacTf);
            T relError = (JacSqr * dx - JacTf).norm() / JacTf.norm();
            PLB_ASSERT(relError < 1.0e-12);
#else
            dx = JacSqr.fullPivLu().solve(JacTf);
            //             JacSqr.lu().solve(JacTf,&dx);
            //             dx = JacSqr.fullPivLu().solve(JacTf);
#endif

            T stepMult = (T)1;   // step size (step mult can only be <= 1 (usually = 1).
            x += stepMult * dx;  // increment solution
            fromXtoUandPiNeq(x, u, PiNeq, dir);
            if (converge(x, dx, epsilon)) {
                //                 pcout << "Converged after " << iT << " iterations." << std::endl;
                //                 std::cout << x << std::endl << std::endl;
                break;
            }
        }
        //         pcout << "NEVER CONVERGED!!!." << std::endl;
    }

};  // struct generalizedIncomprBoundaryTemplates

}  // namespace plb

#endif  // GENERALIZED_BOUNDARY_DYNAMICS_HH
