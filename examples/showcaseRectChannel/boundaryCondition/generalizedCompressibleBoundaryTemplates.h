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
#ifndef GENERALIZED_COMPRESSIBLE_BOUNDARY_TEMPLATES_H
#define GENERALIZED_COMPRESSIBLE_BOUNDARY_TEMPLATES_H

#include <Eigen3/Cholesky>
#include <Eigen3/Core>
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
struct generalizedComprTempBoundaryTemplates {
    static T equilibriumOverRho(plint iPop, const Array<T, Descriptor<T>::d> &u, T uSqr, T thetaBar)
    {
        typedef Descriptor<T> L;
        T c_u = Descriptor<T>::c[iPop][0] * u[0];
        for (int iD = 1; iD < Descriptor<T>::d; ++iD) {
            c_u += Descriptor<T>::c[iPop][iD] * u[iD];
        }

        return L::t[iPop]
               * ((T)1 + L::invCs2 * c_u
                  + L::invCs2 / (T)2
                        * ((L::invCs2 * c_u * c_u - uSqr)
                           + thetaBar * (L::cNormSqr[iPop] - L::cs2 * L::d))
                  + L::invCs2 * L::invCs2 * L::invCs2 / (T)6 * c_u
                        * (c_u * c_u - (T)3 * L::cs2 * uSqr
                           + (T)3 * L::cs2 * thetaBar * (L::cNormSqr[iPop] - L::cs2 * (L::d + 2)))
                  + L::invCs2 * L::invCs2 * L::invCs2 * L::invCs2 / (T)24
                        * (c_u * c_u * c_u * c_u - (T)6 * L::cs2 * uSqr * c_u * c_u
                           + (T)3 * L::cs2 * L::cs2 * uSqr * uSqr
                           + (T)6 * L::cs2 * thetaBar
                                 * (c_u * c_u * (L::cNormSqr[iPop] - L::cs2 * (L::d + 4))
                                    + L::cs2 * uSqr * (L::cs2 * (L::d + 2) - L::cNormSqr[iPop]))
                           + (T)3 * L::cs2 * L::cs2 * thetaBar * thetaBar
                                 * (L::cNormSqr[iPop] * L::cNormSqr[iPop]
                                    - (T)2 * L::cs2 * (L::d + 2) * L::cNormSqr[iPop]
                                    + L::cs2 * L::cs2 * L::d * (L::d + 2))));
    }

    // ================================================================================== //
    // ================================================================================== //
    // = Methods used for boundary celles located "far" from the wall (not on the wall) = //
    // ================================= Non-Linear case ================================ //
    // ================================================================================== //
    // ================================================================================== //
    static void computeDiffF(
        plint iPop, T rho, const Array<T, Descriptor<T>::d> &u, T uSqr, T thetaBar,
        Eigen::RowVectorXd &df)
    {
        df(0) = equilibriumOverRho(iPop, u, uSqr, thetaBar);  // df/drho

        Array<T, SymmetricTensor<T, Descriptor>::n> H2 =
            HermiteTemplate<T, Descriptor>::order2(iPop);
        Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> H3 =
            HermiteTemplate<T, Descriptor>::order3(iPop);
        Array<T, SymmetricRankFourTensor<T, Descriptor>::n> H4 =
            HermiteTemplate<T, Descriptor>::order4(iPop);

        Array<T, Descriptor<T>::d> H3contracted, H2_u, H3_u2, H4_u, H4_u3;
        SymmetricRankThreeTensor<T, Descriptor>::contractLastTwoIndexes(H3, H3contracted);

        Array<T, SymmetricTensor<T, Descriptor>::n> H4contracted;
        SymmetricRankFourTensor<T, Descriptor>::contractLastTwoIndexes(H4, H4contracted);

        Array<T, SymmetricTensor<T, Descriptor>::n> u2;
        VectorTemplate<T, Descriptor>::symTensorProduct(u, u, u2);
        SymmetricTensor<T, Descriptor>::matVectMult(H2, u, H2_u);
        SymmetricRankThreeTensor<T, Descriptor>::multWithRankTwoSymTensor(H3, u2, H3_u2);

        Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> u3;
        VectorTemplate<T, Descriptor>::symTensorProduct(u, u, u, u3);
        SymmetricRankFourTensor<T, Descriptor>::multWithRankThreeSymTensor(H4, u3, H4_u3);
        SymmetricTensor<T, Descriptor>::matVectMult(H4contracted, u, H4_u);

        T tcs2 = Descriptor<T>::invCs2 * Descriptor<T>::t[iPop];
        T tcs4 = tcs2 * Descriptor<T>::invCs2;
        T tcs6_2 = 0.5 * tcs4 * Descriptor<T>::invCs2;
        T cs2_thetaBar = Descriptor<T>::cs2 * thetaBar;
        T tcs8_6 = tcs6_2 * Descriptor<T>::invCs2 / (T)3;
        T cs2_thetaBar_3 = cs2_thetaBar * (T)3;
        for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
            df(1 + iD) = tcs2 * Descriptor<T>::c[iPop][iD];
            df(1 + iD) += tcs4 * H2_u[iD];
            df(1 + iD) += tcs6_2 * (H3_u2[iD] + cs2_thetaBar * H3contracted[iD]);
            df(1 + iD) += tcs8_6 * (H4_u3[iD] + cs2_thetaBar_3 * H4_u[iD]);

            df(1 + iD) *= rho;
        }
        df(1 + Descriptor<T>::d) = tcs2 * SymmetricTensor<T, Descriptor>::trace(H2);
        df(1 + Descriptor<T>::d) +=
            tcs4 / (T)6 * VectorTemplate<T, Descriptor>::scalarProduct(H3contracted, u);
        df(1 + Descriptor<T>::d) +=
            tcs6_2 / (T)2
            * (SymmetricTensor<T, Descriptor>::contractIndexes(H4contracted, u2)
               + cs2_thetaBar * SymmetricTensor<T, Descriptor>::trace(H4contracted));

        H2 = HermiteTemplate<T, Descriptor>::contractedOrder2(iPop);
        for (plint iPi = 0; iPi < SymmetricTensor<T, Descriptor>::n; ++iPi) {
            df(1 + Descriptor<T>::d + 1 + iPi) = H2[iPi] * tcs4 / (T)2;
        }
        H3 = HermiteTemplate<T, Descriptor>::contractedOrder3(iPop);
        for (plint iPi = 0; iPi < SymmetricRankThreeTensor<T, Descriptor>::n; ++iPi) {
            df(1 + Descriptor<T>::d + 1 + SymmetricTensor<T, Descriptor>::n + iPi) =
                H3[iPi] * tcs6_2 / (T)3;
        }
    }

    static void computeDiffFtrLessPiNeq(
        plint iPop, T rho, const Array<T, Descriptor<T>::d> &u, T uSqr, T thetaBar,
        Eigen::RowVectorXd &df)
    {
        df(0) = equilibriumOverRho(iPop, u, uSqr, thetaBar);  // df/drho

        Array<T, SymmetricTensor<T, Descriptor>::n> H2 =
            HermiteTemplate<T, Descriptor>::order2(iPop);
        Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> H3 =
            HermiteTemplate<T, Descriptor>::order3(iPop);
        Array<T, SymmetricRankFourTensor<T, Descriptor>::n> H4 =
            HermiteTemplate<T, Descriptor>::order4(iPop);

        Array<T, Descriptor<T>::d> H3contracted, H2_u, H3_u2, H4_u, H4_u3;
        SymmetricRankThreeTensor<T, Descriptor>::contractLastTwoIndexes(H3, H3contracted);

        Array<T, SymmetricTensor<T, Descriptor>::n> H4contracted;
        SymmetricRankFourTensor<T, Descriptor>::contractLastTwoIndexes(H4, H4contracted);

        Array<T, SymmetricTensor<T, Descriptor>::n> u2;
        VectorTemplate<T, Descriptor>::symTensorProduct(u, u, u2);
        SymmetricTensor<T, Descriptor>::matVectMult(H2, u, H2_u);
        SymmetricRankThreeTensor<T, Descriptor>::multWithRankTwoSymTensor(H3, u2, H3_u2);

        Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> u3;
        VectorTemplate<T, Descriptor>::symTensorProduct(u, u, u, u3);
        SymmetricRankFourTensor<T, Descriptor>::multWithRankThreeSymTensor(H4, u3, H4_u3);
        SymmetricTensor<T, Descriptor>::matVectMult(H4contracted, u, H4_u);

        T tcs2 = Descriptor<T>::invCs2 * Descriptor<T>::t[iPop];
        T tcs4 = tcs2 * Descriptor<T>::invCs2;
        T tcs6_2 = 0.5 * tcs4 * Descriptor<T>::invCs2;
        T cs2_thetaBar = Descriptor<T>::cs2 * thetaBar;
        T tcs8_6 = tcs6_2 * Descriptor<T>::invCs2 / (T)3;
        T cs2_thetaBar_3 = cs2_thetaBar * (T)3;
        for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
            df(1 + iD) = tcs2 * Descriptor<T>::c[iPop][iD];
            df(1 + iD) += tcs4 * H2_u[iD];
            df(1 + iD) += tcs6_2 * (H3_u2[iD] + cs2_thetaBar * H3contracted[iD]);
            df(1 + iD) += tcs8_6 * (H4_u3[iD] + cs2_thetaBar_3 * H4_u[iD]);

            df(1 + iD) *= rho;
        }
        df(1 + Descriptor<T>::d) = tcs2 * SymmetricTensor<T, Descriptor>::trace(H2);
        df(1 + Descriptor<T>::d) +=
            tcs4 / (T)6 * VectorTemplate<T, Descriptor>::scalarProduct(H3contracted, u);
        df(1 + Descriptor<T>::d) +=
            tcs6_2 / (T)2
            * (SymmetricTensor<T, Descriptor>::contractIndexes(H4contracted, u2)
               + cs2_thetaBar * SymmetricTensor<T, Descriptor>::trace(H4contracted));

        H2 = HermiteTemplate<T, Descriptor>::contractedOrder2(iPop);
        df(2 + Descriptor<T>::d) =
            (H2[SymmetricTensorImpl<T, 2>::xx] - H2[SymmetricTensorImpl<T, 2>::yy]) * tcs4 / (T)2;
        df(2 + Descriptor<T>::d + 1) = H2[SymmetricTensorImpl<T, 2>::xy] * tcs4 / (T)2;
        H3 = HermiteTemplate<T, Descriptor>::contractedOrder3(iPop);
        for (plint iPi = 0; iPi < SymmetricRankThreeTensor<T, Descriptor>::n; ++iPi) {
            df(1 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n + iPi) =
                H3[iPi] * tcs6_2 / (T)3;
        }
    }

    static void computeJacobian(
        T rho, const Array<T, Descriptor<T>::d> &u, T thetaBar,
        const std::vector<plint> &knownIndices, Eigen::MatrixXd &Jac)
    {
        plint systSizeX = SymmetricRankThreeTensor<T, Descriptor>::n
                          + SymmetricTensor<T, Descriptor>::n + 1 + Descriptor<T>::d + 1;

        plint systSizeY = knownIndices.size();

        Jac = Eigen::MatrixXd::Zero(systSizeY, systSizeX);
        T uSqr = VectorTemplate<T, Descriptor>::normSqr(u);

        Eigen::RowVectorXd df = Eigen::RowVectorXd::Zero(systSizeX);
        for (pluint iPop = 0; iPop < knownIndices.size(); ++iPop) {
            computeDiffF(knownIndices[iPop], rho, u, uSqr, thetaBar, df);
            //             Jac.row(iPop+1) = df;
            Jac.row(iPop) = df;
        }
    }

    static void computeJacobianTrLessPiNeq(
        T rho, const Array<T, Descriptor<T>::d> &u, T thetaBar,
        const std::vector<plint> &knownIndices, Eigen::MatrixXd &Jac)
    {
        plint systSizeX = SymmetricRankThreeTensor<T, Descriptor>::n
                          + SymmetricTensor<T, Descriptor>::n + 1 + Descriptor<T>::d;

        plint systSizeY = knownIndices.size();

        Jac = Eigen::MatrixXd::Zero(systSizeY, systSizeX);
        T uSqr = VectorTemplate<T, Descriptor>::normSqr(u);

        Eigen::RowVectorXd df = Eigen::RowVectorXd::Zero(systSizeX);
        for (pluint iPop = 0; iPop < knownIndices.size(); ++iPop) {
            computeDiffFtrLessPiNeq(knownIndices[iPop], rho, u, uSqr, thetaBar, df);
            //             Jac.row(iPop+1) = df;
            Jac.row(iPop) = df;
        }
    }

    static void computeNonLinearFunction(
        const Cell<T, Descriptor> &cell, T rho, const Array<T, Descriptor<T>::d> &u, const T uSqr,
        T thetaBar, const Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq,
        const Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> &Qneq,
        const std::vector<plint> &knownIndices, const std::vector<plint> &missingIndices,
        Eigen::VectorXd &f)
    {
        plint systSizeY = knownIndices.size();

        f = Eigen::VectorXd::Zero(systSizeY);

        T rhoBar = Descriptor<T>::rhoBar(rho);
        T invRho = Descriptor<T>::invRho(rhoBar);
        Array<T, Descriptor<T>::d> j = rho * u;
        T jSqr = rho * rho * uSqr;
        for (pluint iPop = 0; iPop < knownIndices.size(); ++iPop) {
            f(iPop) = cell[knownIndices[iPop]]
                      - (dynamicsTemplates<T, Descriptor>::bgk_ma4_equilibrium(
                             knownIndices[iPop], rhoBar, invRho, j, jSqr, thetaBar)
                         + offEquilibriumTemplates<T, Descriptor>::fromPiAndQtoFneq(
                             knownIndices[iPop], PiNeq, Qneq));
            //             pcout << ", f[" << iPop << "] = " << f(iPop+1);
        }
        //         pcout << std::endl;
    }

    static bool converge(
        Eigen::VectorXd &x, const T rho, const Array<T, Descriptor<T>::d> &u, const T thetaBar,
        const Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq,
        const Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> &Qneq, T epsilon, T &sumRes)
    {
        //         pcout << "x = ";
        T threshold = 1.0e-13;
        Array<
            T, SymmetricTensor<T, Descriptor>::n + 2 + Descriptor<T>::d
                   + SymmetricRankThreeTensor<T, Descriptor>::n>
            res;
        if (std::fabs(rho) > threshold) {
            res[0] = std::fabs(x(0) / rho);
        } else
            res[0] = std::fabs(x(0));
        //         pcout << x(0) << ", ";

        for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
            //             pcout << x(iD+1) << ", ";
            if (std::fabs(u[iD]) > threshold) {
                res[iD + 1] = std::fabs(x(iD + 1) / u[iD]);
            } else
                res[iD + 1] = std::fabs(x(iD + 1));
        }
        //         pcout << std::endl;
        if (std::fabs(thetaBar) > threshold) {
            res[1 + Descriptor<T>::d] = std::fabs(x(1 + Descriptor<T>::d) / thetaBar);
        } else
            res[1 + Descriptor<T>::d] = std::fabs(x(1 + Descriptor<T>::d));

        for (plint iPi = 0; iPi < SymmetricTensor<T, Descriptor>::n; ++iPi) {
            //             pcout << x(iPi+1) << ", ";
            if (std::fabs(PiNeq[iPi]) > threshold) {
                res[iPi + 2 + Descriptor<T>::d] =
                    std::fabs(x(iPi + 2 + Descriptor<T>::d) / PiNeq[iPi]);
            } else
                res[iPi + 2 + Descriptor<T>::d] = std::fabs(x(iPi + 2 + Descriptor<T>::d));
        }
        //         pcout << std::endl;

        //         pcout << u[dir] << ", " << PiNeq[0] << ", " << PiNeq[1] << ", " << PiNeq[2] <<
        //         std::endl;
        for (plint iPi = 0; iPi < SymmetricRankThreeTensor<T, Descriptor>::n; ++iPi) {
            //             pcout << x(iPi+1) << ", ";
            if (std::fabs(Qneq[iPi]) > threshold) {
                res[iPi + 2 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n] = std::fabs(
                    x(iPi + 2 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n) / Qneq[iPi]);
            } else
                res[iPi + 2 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n] =
                    std::fabs(x(iPi + 2 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n));
        }

        //         pcout << "res = ";
        sumRes = T();
        for (plint iPi = 0; iPi < SymmetricRankThreeTensor<T, Descriptor>::n
                                      + SymmetricTensor<T, Descriptor>::n + 2 + Descriptor<T>::d;
             ++iPi)
        {
            sumRes += std::fabs(res[iPi]);
            //             pcout << res[iPi] << ", ";
            if (res[iPi] > epsilon) {
                return false;
            }
        }
        //         pcout << std::endl;
        //         pcout << "num = " << u[dir] << ", " << PiNeq[0] << ", " << PiNeq[1] << ", " <<
        //         PiNeq[2] << std::endl;
        return true;
    }

    static bool convergeTrLessPiNeq(
        Eigen::VectorXd &x, const T rho, const Array<T, Descriptor<T>::d> &u, const T thetaBar,
        const Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq,
        const Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> &Qneq, T epsilon, T &sumRes)
    {
        //         pcout << "x = ";
        T threshold = 1.0e-13;
        Array<
            T, SymmetricTensor<T, Descriptor>::n + 2 + Descriptor<T>::d
                   + SymmetricRankThreeTensor<T, Descriptor>::n>
            res;
        if (std::fabs(rho) > threshold) {
            res[0] = std::fabs(x(0) / rho);
        } else
            res[0] = std::fabs(x(0));
        //         pcout << x(0) << ", ";

        for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
            //             pcout << x(iD+1) << ", ";
            if (std::fabs(u[iD]) > threshold) {
                res[iD + 1] = std::fabs(x(iD + 1) / u[iD]);
            } else
                res[iD + 1] = std::fabs(x(iD + 1));
        }
        //         pcout << std::endl;
        if (std::fabs(thetaBar) > threshold) {
            res[1 + Descriptor<T>::d] = std::fabs(x(1 + Descriptor<T>::d) / thetaBar);
        } else
            res[1 + Descriptor<T>::d] = std::fabs(x(1 + Descriptor<T>::d));

        for (plint iPi = 0; iPi < SymmetricTensor<T, Descriptor>::n - 1; ++iPi) {
            //             pcout << x(iPi+1) << ", ";
            if (std::fabs(PiNeq[iPi]) > threshold) {
                res[iPi + 2 + Descriptor<T>::d] =
                    std::fabs(x(iPi + 2 + Descriptor<T>::d) / PiNeq[iPi]);
            } else
                res[iPi + 2 + Descriptor<T>::d] = std::fabs(x(iPi + 2 + Descriptor<T>::d));
        }
        //         pcout << std::endl;

        //         pcout << u[dir] << ", " << PiNeq[0] << ", " << PiNeq[1] << ", " << PiNeq[2] <<
        //         std::endl;
        for (plint iPi = 0; iPi < SymmetricRankThreeTensor<T, Descriptor>::n; ++iPi) {
            //             pcout << x(iPi+1) << ", ";
            if (std::fabs(Qneq[iPi]) > threshold) {
                res[iPi + 1 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n] = std::fabs(
                    x(iPi + 1 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n) / Qneq[iPi]);
            } else
                res[iPi + 1 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n] =
                    std::fabs(x(iPi + 1 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n));
        }

        //         pcout << "res = ";
        sumRes = T();
        for (plint iPi = 0; iPi < SymmetricRankThreeTensor<T, Descriptor>::n
                                      + SymmetricTensor<T, Descriptor>::n + 1 + Descriptor<T>::d;
             ++iPi)
        {
            sumRes += std::fabs(res[iPi]);
            //             pcout << res[iPi] << ", ";
            if (res[iPi] > epsilon) {
                return false;
            }
        }
        //         pcout << std::endl;
        //         pcout << "num = " << u[dir] << ", " << PiNeq[0] << ", " << PiNeq[1] << ", " <<
        //         PiNeq[2] << std::endl;
        return true;
    }

    static bool iterativelySolveSystem(
        const Cell<T, Descriptor> &cell, T &rho, Array<T, Descriptor<T>::d> &u, T &thetaBar,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq,
        Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> &Qneq,
        const std::vector<plint> &knownIndices, const std::vector<plint> &missingIndices, T epsilon,
        T &resSum)
    {
        // rho, u, PiNeq contain the initial guess for the solution of the system
        plint maxT = 10000;
        Eigen::VectorXd f = Eigen::VectorXd::Zero(knownIndices.size());
        Eigen::VectorXd x;

        Eigen::MatrixXd Jac;
        for (plint iT = 0; iT < maxT; ++iT) {
            T uSqr = VectorTemplate<T, Descriptor>::normSqr(u);
            computeNonLinearFunction(
                cell, rho, u, uSqr, thetaBar, PiNeq, Qneq, knownIndices, missingIndices, f);
            computeJacobian(rho, u, thetaBar, knownIndices, Jac);
            Eigen::MatrixXd JacT = Jac.transpose();

            Jac = JacT * Jac;
            f = JacT * f;

#ifdef PLB_DEBUG
            //             bool solutionExists = Jac.lu().solve(f,&x);   // using a LU factorization
            //             PLB_ASSERT(solutionExists);
            x = Jac.fullPivLu().solve(f);
            T relError = (Jac * x - f).norm() / f.norm();
            PLB_ASSERT(relError < 1.0e-12);
#else
            //             Jac.lu().solve(f,&x);
            x = Jac.fullPivLu().solve(f);
#endif

            T stepMult = (T)1;
            rho += stepMult * x(0);
            for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
                u[iD] += stepMult * x(iD + 1);
            }
            thetaBar += stepMult * x(Descriptor<T>::d + 1);
            for (plint iPi = 0; iPi < SymmetricTensor<T, Descriptor>::n; ++iPi) {
                PiNeq[iPi] += stepMult * x(iPi + 2 + Descriptor<T>::d);
            }
            for (plint iPi = 0; iPi < SymmetricRankThreeTensor<T, Descriptor>::n; ++iPi) {
                Qneq[iPi] +=
                    stepMult * x(iPi + 2 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n);
            }
            if (converge(x, rho, u, thetaBar, PiNeq, Qneq, epsilon, resSum)) {
                //                 pcout << "Res Sum = " << resSum << std::endl;
                //                 pcout << "Converged after " << iT << " iterations." << std::endl;
                return true;
            }
            //             pcout << "Res Sum = " << resSum << std::endl;
        }
        //         pcout << "NEVER CONVERGED!!!." << std::endl;
        return false;
    }

    static bool iterativelySolveSystemTrLessPiNeq(
        const Cell<T, Descriptor> &cell, T &rho, Array<T, Descriptor<T>::d> &u, T &thetaBar,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq,
        Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> &Qneq,
        const std::vector<plint> &knownIndices, const std::vector<plint> &missingIndices, T epsilon,
        T &resSum)
    {
        // rho, u, PiNeq contain the initial guess for the solution of the system
        plint maxT = 10000;
        Eigen::VectorXd f = Eigen::VectorXd::Zero(knownIndices.size());
        Eigen::VectorXd x;

        Eigen::MatrixXd Jac;
        for (plint iT = 0; iT < maxT; ++iT) {
            T uSqr = VectorTemplate<T, Descriptor>::normSqr(u);
            computeNonLinearFunction(
                cell, rho, u, uSqr, thetaBar, PiNeq, Qneq, knownIndices, missingIndices, f);
            computeJacobianTrLessPiNeq(rho, u, thetaBar, knownIndices, Jac);
            Eigen::MatrixXd JacT = Jac.transpose();

            Jac = JacT * Jac;
            f = JacT * f;

#ifdef PLB_DEBUG
            //             bool solutionExists = Jac.lu().solve(f,&x);   // using a LU factorization
            //             PLB_ASSERT(solutionExists);
            x = Jac.fullPivLu().solve(f);
            T relError = (Jac * x - f).norm() / f.norm();
            PLB_ASSERT(relError < 1.0e-12);
#else
            //             Jac.lu().solve(f,&x);
            x = Jac.fullPivLu().solve(f);
#endif

            T stepMult = (T)1;
            rho += stepMult * x(0);
            for (plint iD = 0; iD < Descriptor<T>::d; ++iD) {
                u[iD] += stepMult * x(iD + 1);
            }
            thetaBar += stepMult * x(Descriptor<T>::d + 1);
            for (plint iPi = 0; iPi < SymmetricTensor<T, Descriptor>::n - 1; ++iPi) {
                PiNeq[iPi] += stepMult * x(iPi + 2 + Descriptor<T>::d);
            }
            PiNeq[2] += -stepMult * x(2 + Descriptor<T>::d);
            for (plint iPi = 0; iPi < SymmetricRankThreeTensor<T, Descriptor>::n; ++iPi) {
                Qneq[iPi] +=
                    stepMult * x(iPi + 1 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n);
            }
            if (convergeTrLessPiNeq(x, rho, u, thetaBar, PiNeq, Qneq, epsilon, resSum)) {
                //                 pcout << "Res Sum = " << resSum << std::endl;
                //                 pcout << "Converged after " << iT << " iterations." << std::endl;
                return true;
            }
            //             pcout << "Res Sum = " << resSum << std::endl;
        }
        //         pcout << "NEVER CONVERGED!!!." << std::endl;
        return false;
    }

    // ==================================================================== //
    // ==================================================================== //
    // ============= Methods used for ON-WALL boundary nodes ============================ //
    // ==================================================================== //
    // ==================================================================== //

    // = matrix row computation for the decomposition of f_i=rho*t_i*g_i+H2:piNeq+H3:qNeq= //
    // = it also imposes tha mass conservation============================================ //
    static void computeMatrixRow(
        plint iPop, const Array<T, Descriptor<T>::d> &u, T uSqr, T thetaBar, T omega,
        Eigen::RowVectorXd &a)
    {
        T eqOverRho = equilibriumOverRho(iPop, u, uSqr, thetaBar);
        a(0) = eqOverRho;

        T oneMinusOmega = (T)1 - omega;
        T factor = oneMinusOmega * Descriptor<T>::t[iPop] * Descriptor<T>::invCs2
                   * Descriptor<T>::invCs2 / (T)2;
        Array<T, SymmetricTensor<T, Descriptor>::n> H2 =
            HermiteTemplate<T, Descriptor>::contractedOrder2(iPop);
        for (plint iPi = 0; iPi < SymmetricTensor<T, Descriptor>::n; ++iPi) {
            a(iPi + 1) = H2[iPi] * factor;
        }

        factor *= Descriptor<T>::invCs2 / (T)3;
        Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> H3 =
            HermiteTemplate<T, Descriptor>::contractedOrder3(iPop);
        for (plint iPi = 0; iPi < SymmetricRankThreeTensor<T, Descriptor>::n; ++iPi) {
            a(iPi + 1 + SymmetricTensor<T, Descriptor>::n) = H3[iPi] * factor;
        }
    }

    // = matrix row computation for the decomposition of f_i=rho*t_i*g_i+H2:piNeq+H3:qNeq =//
    // = it also imposes tha mass conservation============================================ //
    // In this case the assumption that PiNeq is traceless is done.========================//
    static void computeMatrixRowTrLessPiNeq(
        plint iPop, const Array<T, Descriptor<T>::d> &u, T uSqr, T thetaBar, T omega,
        Eigen::RowVectorXd &a)
    {
        typedef SymmetricTensorImpl<T, Descriptor<T>::d> S;
        T eqOverRho = equilibriumOverRho(iPop, u, uSqr, thetaBar);
        a(0) = eqOverRho;

        T oneMinusOmega = (T)1 - omega;
        T factor = oneMinusOmega * Descriptor<T>::t[iPop] * Descriptor<T>::invCs2
                   * Descriptor<T>::invCs2 / (T)2;
        Array<T, SymmetricTensor<T, Descriptor>::n> H2 =
            HermiteTemplate<T, Descriptor>::contractedOrder2(iPop);

        PLB_ASSERT(Descriptor<T>::d == 2);
        // TODO generalize for 2 and 3d cases....
        if (Descriptor<T>::d == 2) {
            a(1) = (H2[S::xx] - H2[S::yy]) * factor;
            a(2) = H2[S::xy] * factor;
        }

        factor *= Descriptor<T>::invCs2 / (T)3;
        Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> H3 =
            HermiteTemplate<T, Descriptor>::contractedOrder3(iPop);
        for (plint iPi = 0; iPi < SymmetricRankThreeTensor<T, Descriptor>::n; ++iPi) {
            a(iPi + SymmetricTensor<T, Descriptor>::n) = H3[iPi] * factor;
        }
    }

    // = matrix row computation for the decomposition of f_i=rho*t_i*g_i+H2:piNeq+H3:qNeq= //
    static void computeMatrixRow(
        plint iPop, const Array<T, Descriptor<T>::d> &u, T uSqr, T thetaBar, Eigen::RowVectorXd &a)
    {
        T eqOverRho = equilibriumOverRho(iPop, u, uSqr, thetaBar);
        a(0) = eqOverRho;
        T factor = Descriptor<T>::t[iPop] * Descriptor<T>::invCs2 * Descriptor<T>::invCs2 / (T)2;
        Array<T, SymmetricTensor<T, Descriptor>::n> H2 =
            HermiteTemplate<T, Descriptor>::contractedOrder2(iPop);

        for (plint iPi = 0; iPi < SymmetricTensor<T, Descriptor>::n; ++iPi) {
            a(iPi + 1) = H2[iPi] * factor;
        }

        factor *= Descriptor<T>::invCs2 / (T)3;
        Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> H3 =
            HermiteTemplate<T, Descriptor>::contractedOrder3(iPop);

        for (plint iPi = 0; iPi < SymmetricRankThreeTensor<T, Descriptor>::n; ++iPi) {
            a(iPi + 1 + SymmetricTensor<T, Descriptor>::n) = H3[iPi] * factor;
        }
    }

    // = matrix row computation for the decomposition of f_i=rho*t_i*g_i+H2:piNeq+H3:qNeq =//
    // In this case the assumption that PiNeq is traceless is done.========================//
    static void computeMatrixRowTrLessPiNeq(
        plint iPop, const Array<T, Descriptor<T>::d> &u, T uSqr, T thetaBar, Eigen::RowVectorXd &a)
    {
        typedef SymmetricTensorImpl<T, Descriptor<T>::d> S;

        T eqOverRho = equilibriumOverRho(iPop, u, uSqr, thetaBar);
        a(0) = eqOverRho;
        T factor = Descriptor<T>::t[iPop] * Descriptor<T>::invCs2 * Descriptor<T>::invCs2 / (T)2;
        Array<T, SymmetricTensor<T, Descriptor>::n> H2 =
            HermiteTemplate<T, Descriptor>::contractedOrder2(iPop);

        PLB_ASSERT(Descriptor<T>::d == 2);
        // TODO generalize for 2 and 3d cases....
        if (Descriptor<T>::d == 2) {
            a(1) = (H2[S::xx] - H2[S::yy]) * factor;
            a(2) = H2[S::xy] * factor;
        }

        factor *= Descriptor<T>::invCs2 / (T)3;
        Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> H3 =
            HermiteTemplate<T, Descriptor>::contractedOrder3(iPop);

        for (plint iPi = 0; iPi < SymmetricRankThreeTensor<T, Descriptor>::n; ++iPi) {
            a(iPi + SymmetricTensor<T, Descriptor>::n) = H3[iPi] * factor;
        }
    }

    // creation of the over-determined (usually) linear system
    // for an on-wall boundary node
    static void createLinearSystem(
        const Cell<T, Descriptor> &cell, const Array<T, Descriptor<T>::d> &u, T thetaBar,
        const std::vector<plint> &missingIndices, const std::vector<plint> &knownIndices,
        Eigen::MatrixXd &A, Eigen::VectorXd &b, bool massConservation)
    {
        T uSqr = VectorTemplate<T, Descriptor>::normSqr(u);
        T omega = cell.getDynamics().getOmega();

        plint systSizeX =
            SymmetricRankThreeTensor<T, Descriptor>::n + SymmetricTensor<T, Descriptor>::n + 1;
        //         plint systSizeX =
        //         SymmetricRankThreeTensor<T,Descriptor>::n+SymmetricTensor<T,Descriptor>::n; //
        //         PiNeq is traceless...
        plint systSizeY = knownIndices.size() + 1;

        // matrix of the system Ax=b
        A = Eigen::MatrixXd::Zero(systSizeY, systSizeX);
        // rhs of the equation Ax=b
        b = Eigen::VectorXd::Zero(systSizeY);

        if (massConservation) {
            T rhoIn = T();
            //             T rhoTest = T();
            for (pluint mInd = 0; mInd < missingIndices.size(); ++mInd) {
                plint iPop = missingIndices[mInd];
                plint opp = indexTemplates::opposite<Descriptor<T> >(iPop);
                rhoIn += fullF<T, Descriptor>(cell[opp], opp);

                //                 rhoTest += Descriptor<T>::t[iPop];
            }
            //             pcout << ", rhoTest = " << rhoTest << ", ";
            b(knownIndices.size()) = rhoIn;
            Eigen::RowVectorXd sumA = Eigen::RowVectorXd::Zero(systSizeX);
            //             pcout << ", uSqr = " << uSqr << ", thetaBar = " << thetaBar << ", omega =
            //             " << omega;
            for (pluint fInd = 0; fInd < missingIndices.size(); ++fInd) {
                plint iPop = missingIndices[fInd];
                Eigen::RowVectorXd lineA = Eigen::RowVectorXd::Zero(systSizeX);
                computeMatrixRow(iPop, u, uSqr, thetaBar, omega, lineA);
                for (plint iVec = 0; iVec < systSizeX; ++iVec)
                    sumA(iVec) += lineA(iVec);
            }
            A.row(knownIndices.size()) = sumA;
            //             pcout << ", rhoTest2 = " << sumA(0) << ", ";

        } else {
            T rhoTmp = T();
            for (pluint kInd = 0; kInd < knownIndices.size(); ++kInd) {
                plint iPop = knownIndices[kInd];
                rhoTmp += fullF<T, Descriptor>(cell[iPop], iPop);
            }
            // rhoTtmp = sum_i->known f_i.
            b(knownIndices.size()) = rhoTmp;

            // first row of the A matrix. imposing sum_i f_i = rho.
            Eigen::RowVectorXd e0 = Eigen::RowVectorXd::Zero(systSizeX);
            e0(0) = 1.0;

            Eigen::RowVectorXd sumA = Eigen::RowVectorXd::Zero(systSizeX);
            for (pluint fInd = 0; fInd < missingIndices.size(); ++fInd) {
                plint iPop = missingIndices[fInd];
                Eigen::RowVectorXd lineA = Eigen::RowVectorXd::Zero(systSizeX);
                computeMatrixRow(iPop, u, uSqr, thetaBar, lineA);
                for (plint iVec = 0; iVec < systSizeX; ++iVec)
                    sumA(iVec) += lineA(iVec);
            }

            A.row(knownIndices.size()) = e0 - sumA;
        }
        //         pcout << "rhoIn = " << b(0) << ", ";

        for (pluint fInd = 0; fInd < knownIndices.size(); ++fInd) {
            plint iPop = knownIndices[fInd];
            b(fInd) = fullF<T, Descriptor>(cell[iPop], iPop);
            Eigen::RowVectorXd lineA = Eigen::RowVectorXd::Zero(systSizeX);
            computeMatrixRow(iPop, u, uSqr, thetaBar, lineA);
            A.row(fInd) = lineA;
        }
    }

    // We impose tr(PiNeq) = 0 in this case
    // creation of the over-determined (usually) linear system
    static void createLinearSystemTrLessPiNeq(
        const Cell<T, Descriptor> &cell, const Array<T, Descriptor<T>::d> &u, T thetaBar,
        const std::vector<plint> &missingIndices, const std::vector<plint> &knownIndices,
        Eigen::MatrixXd &A, Eigen::VectorXd &b, bool massConservation)
    {
        T uSqr = VectorTemplate<T, Descriptor>::normSqr(u);
        T omega = cell.getDynamics().getOmega();

        plint systSizeX = SymmetricRankThreeTensor<T, Descriptor>::n
                          + SymmetricTensor<T, Descriptor>::n;  // PiNeq is traceless...
        plint systSizeY = knownIndices.size() + 1;

        // matrix of the system Ax=b
        A = Eigen::MatrixXd::Zero(systSizeY, systSizeX);
        // rhs of the equation Ax=b
        b = Eigen::VectorXd::Zero(systSizeY);

        if (massConservation) {
            T rhoIn = T();
            //             T rhoTest = T();
            for (pluint mInd = 0; mInd < missingIndices.size(); ++mInd) {
                plint iPop = missingIndices[mInd];
                plint opp = indexTemplates::opposite<Descriptor<T> >(iPop);
                rhoIn += fullF<T, Descriptor>(cell[opp], opp);

                //                 rhoTest += Descriptor<T>::t[iPop];
            }
            //             pcout << ", rhoTest = " << rhoTest << ", ";
            b(knownIndices.size()) = rhoIn;
            Eigen::RowVectorXd sumA = Eigen::RowVectorXd::Zero(systSizeX);
            //             pcout << ", uSqr = " << uSqr << ", thetaBar = " << thetaBar << ", omega =
            //             " << omega;
            for (pluint fInd = 0; fInd < missingIndices.size(); ++fInd) {
                plint iPop = missingIndices[fInd];
                Eigen::RowVectorXd lineA = Eigen::RowVectorXd::Zero(systSizeX);
                computeMatrixRowTrLessPiNeq(iPop, u, uSqr, thetaBar, omega, lineA);
                for (plint iVec = 0; iVec < systSizeX; ++iVec)
                    sumA(iVec) += lineA(iVec);
            }
            A.row(knownIndices.size()) = sumA;
            //             pcout << ", rhoTest2 = " << sumA(0) << ", ";

        } else {
            T rhoTmp = T();
            for (pluint kInd = 0; kInd < knownIndices.size(); ++kInd) {
                plint iPop = knownIndices[kInd];
                rhoTmp += fullF<T, Descriptor>(cell[iPop], iPop);
            }
            // rhoTtmp = sum_i->known f_i.
            b(knownIndices.size()) = rhoTmp;

            // first row of the A matrix. imposing sum_i f_i = rho.
            Eigen::RowVectorXd e0 = Eigen::RowVectorXd::Zero(systSizeX);
            e0(0) = 1.0;

            Eigen::RowVectorXd sumA = Eigen::RowVectorXd::Zero(systSizeX);
            for (pluint fInd = 0; fInd < missingIndices.size(); ++fInd) {
                plint iPop = missingIndices[fInd];
                Eigen::RowVectorXd lineA = Eigen::RowVectorXd::Zero(systSizeX);
                computeMatrixRowTrLessPiNeq(iPop, u, uSqr, thetaBar, lineA);
                for (plint iVec = 0; iVec < systSizeX; ++iVec)
                    sumA(iVec) += lineA(iVec);
            }

            A.row(knownIndices.size()) = e0 - sumA;
        }
        //         pcout << "rhoIn = " << b(0) << ", ";

        //         pcout << ", rhoIn = " << rhoIn << std::endl;
        for (pluint fInd = 0; fInd < knownIndices.size(); ++fInd) {
            plint iPop = knownIndices[fInd];
            b(fInd) = fullF<T, Descriptor>(cell[iPop], iPop);
            Eigen::RowVectorXd lineA = Eigen::RowVectorXd::Zero(systSizeX);
            computeMatrixRowTrLessPiNeq(iPop, u, uSqr, thetaBar, lineA);
            A.row(fInd) = lineA;
        }
    }

    static void solveLinearSystem(
        const Cell<T, Descriptor> &cell, const Array<T, Descriptor<T>::d> &u, T thetaBar,
        const std::vector<plint> &missingIndices, const std::vector<plint> &knownIndices, T &rho,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq,
        Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> &qNeq, bool massConservation)
    {
        Eigen::MatrixXd A;
        Eigen::VectorXd b;

        createLinearSystem(cell, u, thetaBar, missingIndices, knownIndices, A, b, massConservation);

        Eigen::VectorXd x;

        Eigen::MatrixXd AT = A.transpose();
        A = AT * A;
        b = AT * b;

#ifdef PLB_DEBUG
        //         bool solutionExists = A.lu().solve(b,&x);   // using a LU factorization
        //         PLB_ASSERT(solutionExists);
        x = A.fullPivLu().solve(b);
        T relError = (A * x - b).norm() / b.norm();
        PLB_ASSERT(relError < 1.0e-12);
#else
        //         A.lu().solve(b,&x);
        x = A.fullPivLu().solve(b);
#endif

        rho = x(0);
        //         pcout << "rho = " << rho << ", ";
        for (plint iPi = 0; iPi < SymmetricTensor<T, Descriptor>::n; ++iPi) {
            PiNeq[iPi] = x(iPi + 1);
            //             pcout << "piNeq = " << x(iPi+1) << ", ";
        }
        //         PiNeq[2] = -PiNeq[0];
        for (plint iPi = 0; iPi < SymmetricRankThreeTensor<T, Descriptor>::n; ++iPi) {
            qNeq[iPi] = x(iPi + 1 + SymmetricTensor<T, Descriptor>::n);
            //             pcout << "qNeq = " << x(iPi+1+SymmetricTensor<T,Descriptor>::n) << ", ";
        }

        //         pcout << "TrPiNeq = " << PiNeq[0] - PiNeq[2] << ", PiNeq[xx] = " << PiNeq[0] <<
        //         ", PiNeq[xy] = " << PiNeq[1] << ", PiNeq[yy] = " << PiNeq[2] << std::endl;
    }

    // ======================================================================= //
    // We impose tr(PiNeq) = 0 in this case
    // creation of the over-determined (usually) linear system
    // ======================================================================= //
    static void solveLinearSystemTrLessPiNeq(
        const Cell<T, Descriptor> &cell, const Array<T, Descriptor<T>::d> &u, T thetaBar,
        const std::vector<plint> &missingIndices, const std::vector<plint> &knownIndices, T &rho,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq,
        Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> &qNeq, bool massConservation)
    {
        Eigen::MatrixXd A;
        Eigen::VectorXd b;

        createLinearSystemTrLessPiNeq(
            cell, u, thetaBar, missingIndices, knownIndices, A, b, massConservation);

        Eigen::VectorXd x;

        Eigen::MatrixXd AT = A.transpose();
        A = AT * A;
        b = AT * b;

#ifdef PLB_DEBUG
        //         bool solutionExists = A.lu().solve(b,&x);   // using a LU factorization
        //         PLB_ASSERT(solutionExists);
        x = A.fullPivLu().solve(b);
        T relError = (A * x - b).norm() / b.norm();
        PLB_ASSERT(relError < 1.0e-12);
#else
        x = A.fullPivLu().solve(b);
        //         A.lu().solve(b,&x);
#endif

        rho = x(0);
        //         pcout << "rho = " << rho << ", ";
        for (plint iPi = 0; iPi < SymmetricTensor<T, Descriptor>::n - 1; ++iPi) {
            PiNeq[iPi] = x(iPi + 1);
            //             pcout << "piNeq = " << x(iPi+1) << ", ";
        }
        PiNeq[2] = -PiNeq[0];
        for (plint iPi = 0; iPi < SymmetricRankThreeTensor<T, Descriptor>::n; ++iPi) {
            qNeq[iPi] = x(iPi + SymmetricTensor<T, Descriptor>::n);
            //             pcout << "qNeq = " << x(iPi+1+SymmetricTensor<T,Descriptor>::n) << ", ";
        }

        //         pcout << "TrPiNeq = " << PiNeq[0] - PiNeq[2] << ", PiNeq[xx] = " << PiNeq[0] <<
        //         ", PiNeq[xy] = " << PiNeq[1] << ", PiNeq[yy] = " << PiNeq[2] << std::endl;
    }

    // ================================================================================== //
    // ================================================================================== //
    // = Methods used for boundary celles located "far" from the wall (not on the wall) = //
    // ================================= Linear case ==================================== //
    // ================================================================================== //
    // ================================================================================== //

    // Matrix row corespond to the factors in front of each Hermite coefficient (the Hermite
    // polynomials times the weights). a(0) = w_i H0, a(1,..,d) = w_i/c_s^ 2*H1_ii,
    // a(d+1,..,SymRank2::n)=w_i/(2c_s^4)Q_i a(SymRank2::n+1,..,SymRank3::n) = w_i, a(1,..,d) =
    // w_i/(6c_s^6)*H3_i a(SmyRank3::n+1,...,SmyRank4::n)=w_i/(24c_s^8)H4_i
    static void computeMatrixRow(plint iPop, Eigen::RowVectorXd &a)
    {
        T factor = Descriptor<T>::t[iPop];
        a(0) = (T)1 * factor;

        factor *= Descriptor<T>::invCs2;
        for (plint iPi = 0; iPi < Descriptor<T>::d; ++iPi)
            a(iPi + 1) = Descriptor<T>::c[iPop][iPi] * factor;

        factor *= Descriptor<T>::invCs2 / (T)2;
        Array<T, SymmetricTensor<T, Descriptor>::n> H2 =
            HermiteTemplate<T, Descriptor>::contractedOrder2(iPop);
        for (plint iPi = 0; iPi < SymmetricTensor<T, Descriptor>::n; ++iPi)
            a(iPi + 1 + Descriptor<T>::d) = H2[iPi] * factor;

        factor *= Descriptor<T>::invCs2 / (T)3;
        Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> H3 =
            HermiteTemplate<T, Descriptor>::contractedOrder3(iPop);
        for (plint iPi = 0; iPi < SymmetricRankThreeTensor<T, Descriptor>::n; ++iPi) {
            a(iPi + 1 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n) = H3[iPi] * factor;
        }

        factor *= Descriptor<T>::invCs2 / (T)4;
        Array<T, SymmetricRankFourTensor<T, Descriptor>::n> H4 =
            HermiteTemplate<T, Descriptor>::contractedOrder4(iPop);

        for (plint iPi = 0; iPi < SymmetricRankFourTensor<T, Descriptor>::n; ++iPi) {
            a(iPi + 1 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n
              + SymmetricRankThreeTensor<T, Descriptor>::n) = H4[iPi] * factor;
        }
    }

    // creation of the over-determined (usually) linear system
    // In this case the matrix is composed by the Hermite polynomials and the corresponding weights
    // while the unknowns are the Hermite coefficients.
    static void createLinearSystem(
        const Cell<T, Descriptor> &cell, const std::vector<plint> &missingIndices,
        const std::vector<plint> &knownIndices, Eigen::MatrixXd &A, Eigen::VectorXd &b)
    {
        plint systSizeX = SymmetricRankFourTensor<T, Descriptor>::n
                          + SymmetricRankThreeTensor<T, Descriptor>::n
                          + SymmetricTensor<T, Descriptor>::n + Descriptor<T>::d + 1;

        plint systSizeY = knownIndices.size();
        // matrix of the system Ax=b
        A = Eigen::MatrixXd::Zero(systSizeY, systSizeX);
        // rhs of the equation Ax=b
        b = Eigen::VectorXd::Zero(systSizeY);

        for (pluint fInd = 0; fInd < knownIndices.size(); ++fInd) {
            plint iPop = knownIndices[fInd];
            Eigen::RowVectorXd lineA = Eigen::RowVectorXd::Zero(systSizeX);
            computeMatrixRow(iPop, lineA);
            b(fInd) = fullF<T, Descriptor>(cell[iPop], iPop);
            //             std::cout << lineA << std::endl << std::endl;
            A.row(fInd) = lineA;
        }
    }

    // The linear system composed used to compute explicitely the Hermite coefficients (not directly
    // the macroscopic moments : density, velocity , temperature, ...)
    static void solveLinearSystem(
        const Cell<T, Descriptor> &cell, const std::vector<plint> &missingIndices,
        const std::vector<plint> &knownIndices, std::vector<T> &coeffs)
    {
        Eigen::MatrixXd A;
        Eigen::VectorXd b;

        createLinearSystem(cell, missingIndices, knownIndices, A, b);

        Eigen::VectorXd x;

        Eigen::MatrixXd AT = A.transpose();

        A = AT * A;
        b = AT * b;

#ifdef PLB_DEBUG
        //         bool solutionExists = A.lu().solve(b,&x);   // using a LU factorization
        //         PLB_ASSERT(solutionExists);
        x = A.fullPivLu().solve(b);
        T relError = (A * x - b).norm() / b.norm();
        PLB_ASSERT(relError < 1.0e-12);
#else
        x = A.fullPivLu().solve(b);
        //         A.lu().solve(b,&x);
#endif

        coeffs.push_back(x(0));
        //         pcout << "a0 = " << coeffs[0] << ", ";
        for (plint iPi = 0; iPi < Descriptor<T>::d; ++iPi) {
            coeffs.push_back(x(iPi + 1));
            //             pcout << "a1 = " << coeffs[iPi+1] << ", ";
        }
        for (plint iPi = 0; iPi < SymmetricTensor<T, Descriptor>::n; ++iPi) {
            coeffs.push_back(x(iPi + 1 + Descriptor<T>::d));
            //             pcout << "a2 = " << coeffs[iPi+1+Descriptor<T>::d] << ", ";
        }
        for (plint iPi = 0; iPi < SymmetricRankThreeTensor<T, Descriptor>::n; ++iPi) {
            coeffs.push_back(x(iPi + 1 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n));
            //             pcout << "a3 = " <<
            //             coeffs[iPi+1+Descriptor<T>::d+SymmetricTensor<T,Descriptor>::n] << ", ";
        }
        for (plint iPi = 0; iPi < SymmetricRankFourTensor<T, Descriptor>::n; ++iPi) {
            coeffs.push_back(
                x(iPi + 1 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n
                  + SymmetricRankThreeTensor<T, Descriptor>::n));
            //             pcout << "a4 = " <<
            //             coeffs[iPi+1+Descriptor<T>::d+SymmetricTensor<T,Descriptor>::n+SymmetricRankThreeTensor<T,Descriptor>::n]
            //             << ", ";
        }
        //         pcout << std::endl;
    }

};  // struct generalizedIncomprBoundaryTemplates

}  // namespace plb

#endif  // GENERALIZED_BOUNDARY_DYNAMICS_HH
