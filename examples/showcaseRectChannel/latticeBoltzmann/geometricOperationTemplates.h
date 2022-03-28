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
 * Templates for common geometric operations (scalar product, vector-
 * matric operations etc.).
 *  -- header file
 */
#ifndef GEOMETRIC_OPERATION_TEMPLATES_H
#define GEOMETRIC_OPERATION_TEMPLATES_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "core/array.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "io/parallelIO.h"

namespace plb {

template <typename T, int d>
struct VectorTemplateImpl;
template <typename T, template <typename U> class Descriptor>
struct SymmetricTensor;
template <typename T, template <typename U> class Descriptor>
struct SymmetricRankThreeTensor;
template <typename T, int d>
struct SymmetricTensorImpl;
template <typename T, int d>
struct SymmetricRankThreeTensorImpl;

template <typename T, template <typename U> class Descriptor>
struct VectorTemplate {
    /// Number of dimensions for current lattice
    static const int d = Descriptor<T>::d;
    /// Compute scalar product between two vectors
    static T scalarProduct(Array<T, d> const &u1, Array<T, d> const &u2)
    {
        return VectorTemplateImpl<T, d>::scalarProduct(u1, u2);
    }
    /// Compute scalar product between two a c-array and a plb-array
    static T scalarProduct(const T u1[], Array<T, d> const &u2)
    {
        return VectorTemplateImpl<T, d>::scalarProduct(u1, u2);
    }
    /// Compute scalar product between two a c-array and a plb-array
    static T intScalarProduct(const int u1[], Array<T, d> const &u2)
    {
        return VectorTemplateImpl<T, d>::intScalarProduct(u1, u2);
    }
    /// Compute norm-square of a vector
    static T normSqr(Array<T, d> const &u)
    {
        return VectorTemplateImpl<T, d>::normSqr(u);
    }
    /// Multiply vector elements component-wise by a scalar
    static void multiplyByScalar(Array<T, d> &u, T scalar)
    {
        VectorTemplateImpl<T, d>::multiplyByScalar(u, scalar);
    }
    /// Multiply vector elements component-wise by a scalar and store in second vector
    static void multiplyByScalar(Array<T, d> const &u, T scalar, Array<T, d> &result)
    {
        VectorTemplateImpl<T, d>::multiplyByScalar(u, scalar, result);
    }
    /// Symmetric tensor product of two vectors
    static void symTensorProduct(
        Array<T, d> const &u, Array<T, d> const &v,
        Array<T, SymmetricTensor<T, Descriptor>::n> &result)
    {
        VectorTemplateImpl<T, d>::symTensorProduct(u, v, result);
    }
    /// Symmetric tensor product of three vectors
    static void symTensorProduct(
        Array<T, d> const &u, Array<T, d> const &v, Array<T, d> const &w,
        Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> &result)
    {
        VectorTemplateImpl<T, d>::symTensorProduct(u, v, w, result);
    }
};

template <typename T, int d>
struct VectorTemplateImpl {
    static T scalarProduct(Array<T, d> const &u1, Array<T, d> const &u2)
    {
        T result = T();
        for (int iD = 0; iD < d; ++iD) {
            result += u1[iD] * u2[iD];
        }
        return result;
    }
    static T scalarProduct(const T u1[], Array<T, d> const &u2)
    {
        T result = T();
        for (int iD = 0; iD < d; ++iD) {
            result += u1[iD] * u2[iD];
        }
        return result;
    }
    static T intScalarProduct(const int u1[], Array<T, d> const &u2)
    {
        T result = T();
        for (int iD = 0; iD < d; ++iD) {
            result += u1[iD] * u2[iD];
        }
        return result;
    }
    static T normSqr(Array<T, d> const &u)
    {
        return scalarProduct(u, u);
    }
    static void multiplyByScalar(Array<T, d> &u, T scalar)
    {
        for (int iD = 0; iD < d; ++iD) {
            u[iD] *= scalar;
        }
    }

    static void multiplyByScalar(Array<T, d> const &u, T scalar, Array<T, d> &result)
    {
        for (int iD = 0; iD < d; ++iD) {
            result[iD] = u[iD] * scalar;
        }
    }

    static void symTensorProduct(
        Array<T, d> const &u, Array<T, d> const &v, Array<T, SymmetricTensorImpl<T, d>::n> &result)
    {
        int iPi = 0;
        for (int iA = 0; iA < d; ++iA) {
            for (int iB = iA; iB < d; ++iB) {
                result[iPi] = u[iA] * v[iB];
                ++iPi;
            }
        }
    }

    static void symTensorProduct(
        Array<T, d> const &u, Array<T, d> const &v, Array<T, d> const &w,
        Array<T, SymmetricRankThreeTensorImpl<T, d>::n> &result)
    {
        int iPi = 0;
        for (int iA = 0; iA < d; ++iA) {
            for (int iB = iA; iB < d; ++iB) {
                for (int iC = iB; iC < d; ++iC) {
                    result[iPi] = u[iA] * v[iB] * w[iC];
                    ++iPi;
                }
            }
        }
    }
};

template <typename T>
struct VectorTemplateImpl<T, 2> {
    static T scalarProduct(Array<T, 2> const &u1, Array<T, 2> const &u2)
    {
        return u1[0] * u2[0] + u1[1] * u2[1];
    }
    static T scalarProduct(const T u1[], Array<T, 2> const &u2)
    {
        return u1[0] * u2[0] + u1[1] * u2[1];
    }
    static T intScalarProduct(const int u1[], Array<T, 2> const &u2)
    {
        return u1[0] * u2[0] + u1[1] * u2[1];
    }
    static T normSqr(Array<T, 2> const &u)
    {
        return u[0] * u[0] + u[1] * u[1];
    }
    static void multiplyByScalar(Array<T, 2> &u, T scalar)
    {
        u[0] *= scalar;
        u[1] *= scalar;
    }
    static void multiplyByScalar(Array<T, 2> const &u, T scalar, Array<T, 2> &result)
    {
        result[0] = u[0] * scalar;
        result[1] = u[1] * scalar;
    }
    static void symTensorProduct(Array<T, 2> const &u, Array<T, 2> const &v, Array<T, 3> &result)
    {
        result[0] = u[0] * v[0];
        result[1] = u[0] * v[1];
        result[2] = u[1] * v[1];
    }
    static void symTensorProduct(
        Array<T, 2> const &u, Array<T, 2> const &v, Array<T, 2> const &w, Array<T, 4> &result)
    {
        result[0] = u[0] * v[0] * w[0];
        result[1] = u[0] * v[0] * w[1];
        result[2] = u[0] * v[1] * w[1];
        result[3] = u[1] * v[1] * w[1];
    }
};

template <typename T>
struct VectorTemplateImpl<T, 3> {
    static T scalarProduct(Array<T, 3> const &u1, Array<T, 3> const &u2)
    {
        return u1[0] * u2[0] + u1[1] * u2[1] + u1[2] * u2[2];
    }
    static T scalarProduct(const T u1[], Array<T, 3> const &u2)
    {
        return u1[0] * u2[0] + u1[1] * u2[1] + u1[2] * u2[2];
    }
    static T intScalarProduct(const int u1[], Array<T, 3> const &u2)
    {
        return u1[0] * u2[0] + u1[1] * u2[1] + u1[2] * u2[2];
    }
    static T normSqr(Array<T, 3> const &u)
    {
        return u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    }
    static void multiplyByScalar(Array<T, 3> &u, T scalar)
    {
        u[0] *= scalar;
        u[1] *= scalar;
        u[2] *= scalar;
    }
    static void multiplyByScalar(Array<T, 3> const &u, T scalar, Array<T, 3> &result)
    {
        result[0] = u[0] * scalar;
        result[1] = u[1] * scalar;
        result[2] = u[2] * scalar;
    }
    static void symTensorProduct(Array<T, 3> const &u, Array<T, 3> const &v, Array<T, 6> &result)
    {
        result[0] = u[0] * v[0];
        result[1] = u[0] * v[1];
        result[2] = u[0] * v[2];
        result[3] = u[1] * v[1];
        result[4] = u[1] * v[2];
        result[5] = u[2] * v[2];
    }
    static void symTensorProduct(
        Array<T, 3> const &u, Array<T, 3> const &v, Array<T, 3> const &w, Array<T, 10> &result)
    {
        result[0] = u[0] * v[0] * w[0];
        result[1] = u[0] * v[0] * w[1];
        result[2] = u[0] * v[0] * w[2];
        result[3] = u[0] * v[1] * w[1];
        result[4] = u[0] * v[1] * w[2];
        result[5] = u[0] * v[2] * w[2];
        result[6] = u[1] * v[1] * w[1];
        result[7] = u[1] * v[1] * w[2];
        result[8] = u[1] * v[2] * w[2];
        result[9] = u[2] * v[2] * w[2];
    }
};

template <typename T, int d>
struct SymmetricTensorImpl {
    static const int n = d * (d + 1) / 2;
};

template <typename T>
struct SymmetricTensorImpl<T, 0> {
    static const int n = 0;
    static void matVectMult(Array<T, n> const &mat, Array<T, 0> const &vect, Array<T, 0> &result)
    {
        result.resetToZero();
    }
    static T tensorNormSqr(Array<T, n> const &mat)
    {
        return T();
    }
    static T contractIndexes(Array<T, n> const &A, Array<T, n> const &B)
    {
        return T();
    }
    static T trace(Array<T, n> const &A)
    {
        return T();
    }
    static Array<T, n> id()
    {
        return Array<T, n>::zero();
    }
    static T det(Array<T, n> const &A)
    {
        return T();
    }
    static void inv(Array<T, n> const &A, Array<T, n> &result)
    {
        result.resetToZero();
    }
};

template <typename T>
struct SymmetricTensorImpl<T, 2> {
    static const int n = 3;
    enum Indices { xx = 0, xy = 1, yy = 2 };
    static void matVectMult(Array<T, n> const &mat, Array<T, 2> const &vect, Array<T, 2> &result)
    {
        result[0] = mat[xx] * vect[0] + mat[xy] * vect[1];
        result[1] = mat[xy] * vect[0] + mat[yy] * vect[1];
    }
    static T tensorNormSqr(Array<T, n> const &mat)
    {
        return mat[xx] * mat[xx] + mat[yy] * mat[yy] + (T)2 * mat[xy] * mat[xy];
    }
    static T contractIndexes(Array<T, n> const &A, Array<T, n> const &B)
    {
        return A[xx] * B[xx] + A[yy] * B[yy] + (T)2 * A[xy] * B[xy];
    }
    static T trace(Array<T, n> const &A)
    {
        return A[xx] + A[yy];
    }
    static Array<T, n> id()
    {
        return Array<T, n>((T)1, T(), (T)1);
    }
    static T det(Array<T, n> const &A)
    {
        return A[xx] * A[yy] - (T)2 * A[xy];
    }
    static void inv(Array<T, n> const &A, Array<T, n> &result)
    {
        T iDet = (T)1 / det(A);
        result[xx] = iDet * A[yy];
        result[yy] = iDet * A[xx];
        result[xy] = -iDet * A[xy];
    }
};

template <typename T>
struct SymmetricTensorImpl<T, 3> {
    static const int n = 6;
    enum Indices { xx = 0, xy = 1, xz = 2, yy = 3, yz = 4, zz = 5 };
    static void matVectMult(Array<T, n> const &mat, Array<T, 3> const &vect, Array<T, 3> &result)
    {
        result[0] = mat[xx] * vect[0] + mat[xy] * vect[1] + mat[xz] * vect[2];
        result[1] = mat[xy] * vect[0] + mat[yy] * vect[1] + mat[yz] * vect[2];
        result[2] = mat[xz] * vect[0] + mat[yz] * vect[1] + mat[zz] * vect[2];
    }
    static T tensorNormSqr(Array<T, n> const &mat)
    {
        return mat[xx] * mat[xx] + mat[yy] * mat[yy] + mat[zz] * mat[zz] + (T)2 * mat[xy] * mat[xy]
               + (T)2 * mat[xz] * mat[xz] + (T)2 * mat[yz] * mat[yz];
    }
    static T contractIndexes(Array<T, n> const &A, Array<T, n> const &B)
    {
        return A[xx] * B[xx] + A[yy] * B[yy] + A[zz] * B[zz]
               + (T)2 * (A[xy] * B[xy] + A[xz] * B[xz] + A[yz] * B[yz]);
    }
    static T trace(Array<T, n> const &A)
    {
        return A[xx] + A[yy] + A[zz];
    }
    static Array<T, n> id()
    {
        Array<T, n> I;

        I[0] = (T)1;
        I[1] = T();
        I[2] = T();
        I[3] = (T)1;
        I[4] = T();
        I[5] = (T)1;

        return I;
    }
    static T det(Array<T, n> const &A)
    {
        return A[xx] * (A[yy] * A[zz] - A[yz] * A[yz]) - A[xy] * (A[xy] * A[zz] - A[yz] * A[xz])
               + A[xz] * (A[xy] * A[yz] - A[yy] * A[xz]);
    }
    static void inv(Array<T, n> const &A, Array<T, n> &result)
    {
        T iDet = (T)1 / det(A);
        result[xx] = iDet * (A[yy] * A[zz] - A[yz] * A[yz]);
        result[xy] = iDet * (A[xz] * A[yz] - A[xx] * A[xy]);
        result[xz] = iDet * (A[xy] * A[yz] - A[yy] * A[xz]);
        result[yy] = iDet * (A[xx] * A[zz] - A[xz] * A[xz]);
        result[yz] = iDet * (A[xz] * A[xy] - A[yz] * A[xx]);
        result[zz] = iDet * (A[xx] * A[yy] - A[xy] * A[xy]);
    }
};

/// Operations on a symmetric tensor which stores only above-or-on-diagonal values
template <typename T, template <typename U> class Descriptor>
struct SymmetricTensor {
    /// Number of dimensions for current lattice
    static const int d = Descriptor<T>::d;
    /// Number of elements (reduced by symmetry)
    static const int n = SymmetricTensorImpl<T, d>::n;
    static void matVectMult(Array<T, n> const &mat, Array<T, d> const &vect, Array<T, d> &result)
    {
        SymmetricTensorImpl<T, d>::matVectMult(mat, vect, result);
    }
    static T tensorNormSqr(Array<T, n> const &mat)
    {
        return SymmetricTensorImpl<T, d>::tensorNormSqr(mat);
    }
    static T contractIndexes(Array<T, n> const &A, Array<T, n> const &B)
    {
        return SymmetricTensorImpl<T, d>::contractIndexes(A, B);
    }
    static T trace(Array<T, n> const &A)
    {
        return SymmetricTensorImpl<T, d>::trace(A);
    }
    static Array<T, n> id()
    {
        return SymmetricTensorImpl<T, d>::id();
    }
    static T det(Array<T, n> const &A)
    {
        return SymmetricTensorImpl<T, d>::det(A);
    }
    static void inv(Array<T, n> const &A, Array<T, n> &result)
    {
        SymmetricTensorImpl<T, d>::inv(A, result);
    }
};

template <typename T, int d>
struct SymmetricRankThreeTensorImpl {
    static const int n = d * (d + 1) * (d + 2) / 6;
};

template <typename T>
struct SymmetricRankThreeTensorImpl<T, 0> {
    static const int n = 0;
    static void contractLastTwoIndexes(Array<T, n> const &tens, Array<T, 0> &res)
    {
        res.resetToZero();
    }

    static T contractIndexes(Array<T, n> const &A, Array<T, n> const &B)
    {
        return T();
    }

    static void multWithRankTwoSymTensor(
        Array<T, n> const &A, Array<T, SymmetricTensorImpl<T, 0>::n> const &B, Array<T, 0> &x)
    {
        x.resetToZero();
    }
};

template <typename T>
struct SymmetricRankThreeTensorImpl<T, 2> {
    static const int n = 4;
    enum Indices { xxx = 0, xxy = 1, xyy = 2, yyy = 3 };
    static void contractLastTwoIndexes(Array<T, n> const &tens, Array<T, 2> &res)
    {
        res[0] = tens[xxx] + tens[xyy];
        res[1] = tens[xxy] + tens[yyy];
    }

    static T contractIndexes(Array<T, n> const &A, Array<T, n> const &B)
    {
        return A[xxx] * B[xxx] + A[yyy] * B[yyy] + (T)3 * (A[xxy] * B[xxy] + A[xyy] * B[xyy]);
    }

    static void multWithRankTwoSymTensor(
        Array<T, n> const &A, Array<T, SymmetricTensorImpl<T, 2>::n> const &B, Array<T, 2> &x)
    {
        x[0] = A[xxx] * B[SymmetricTensorImpl<T, 2>::xx]
               + (T)2 * A[xxy] * B[SymmetricTensorImpl<T, 2>::xy]
               + A[xyy] * B[SymmetricTensorImpl<T, 2>::yy];
        x[1] = A[xxy] * B[SymmetricTensorImpl<T, 2>::xx]
               + (T)2 * A[xyy] * B[SymmetricTensorImpl<T, 2>::xy]
               + A[yyy] * B[SymmetricTensorImpl<T, 2>::yy];
    }
};

template <typename T>
struct SymmetricRankThreeTensorImpl<T, 3> {
    static const int n = 10;
    enum Indices {
        xxx = 0,
        xxy = 1,
        xxz = 2,
        xyy = 3,
        xyz = 4,
        xzz = 5,
        yyy = 6,
        yyz = 7,
        yzz = 8,
        zzz = 9
    };
    static void contractLastTwoIndexes(Array<T, n> const &tens, Array<T, 3> &res)
    {
        res[0] = tens[xxx] + tens[xyy] + tens[xzz];
        res[1] = tens[xxy] + tens[yyy] + tens[yzz];
        res[2] = tens[xxz] + tens[yyz] + tens[zzz];
    }

    static T contractIndexes(Array<T, n> const &A, Array<T, n> const &B)
    {
        return A[xxx] * B[xxx] + A[yyy] * B[yyy] + A[zzz] * B[zzz]
               + (T)3
                     * (A[xxy] * B[xxy] + A[xxz] * B[xxz] + A[xyy] * B[xyy] + A[xzz] * B[xzz]
                        + A[yzz] * B[yzz] + A[yyz] * B[yyz])
               + (T)6 * A[xyz] * B[xyz];
    }

    static void multWithRankTwoSymTensor(
        Array<T, n> const &A, Array<T, SymmetricTensorImpl<T, 3>::n> const &B, Array<T, 3> &x)
    {
        x[0] = A[xxx] * B[SymmetricTensorImpl<T, 3>::xx] + A[xyy] * B[SymmetricTensorImpl<T, 3>::yy]
               + A[xzz] * B[SymmetricTensorImpl<T, 3>::zz]
               + (T)2
                     * (A[xxy] * B[SymmetricTensorImpl<T, 3>::xy]
                        + A[xxz] * B[SymmetricTensorImpl<T, 3>::xz]
                        + A[xyz] * B[SymmetricTensorImpl<T, 3>::yz]);
        x[1] = A[yyy] * B[SymmetricTensorImpl<T, 3>::yy] + A[xxy] * B[SymmetricTensorImpl<T, 3>::xx]
               + A[yzz] * B[SymmetricTensorImpl<T, 3>::zz]
               + (T)2
                     * (A[xyy] * B[SymmetricTensorImpl<T, 3>::xy]
                        + A[yyz] * B[SymmetricTensorImpl<T, 3>::yz]
                        + A[xyz] * B[SymmetricTensorImpl<T, 3>::xz]);
        x[2] = A[zzz] * B[SymmetricTensorImpl<T, 3>::zz] + A[xxz] * B[SymmetricTensorImpl<T, 3>::xx]
               + A[yyz] * B[SymmetricTensorImpl<T, 3>::yy]
               + (T)2
                     * (A[xzz] * B[SymmetricTensorImpl<T, 3>::xz]
                        + A[yzz] * B[SymmetricTensorImpl<T, 3>::yz]
                        + A[xyz] * B[SymmetricTensorImpl<T, 3>::xy]);
    }
};

/// Operations on a symmetric tensor which stores only above-or-on-diagonal values
template <typename T, template <typename U> class Descriptor>
struct SymmetricRankThreeTensor {
    /// Number of dimensions for current lattice
    static const int d = Descriptor<T>::d;
    /// Number of elements (reduced by symmetry)
    static const int n = SymmetricRankThreeTensorImpl<T, d>::n;
    /// computes the contraction of the last two indexes of the tensor
    static void contractLastTwoIndexes(Array<T, n> const &tens, Array<T, d> &res)
    {
        SymmetricRankThreeTensorImpl<T, d>::contractLastTwoIndexes(tens, res);
    }
    /// computes the contraction of the last two indexes of the tensor
    static T contractIndexes(Array<T, n> const &A, Array<T, n> const &B)
    {
        return SymmetricRankThreeTensorImpl<T, d>::contractIndexes(A, B);
    }
    static void multWithRankTwoSymTensor(
        Array<T, n> const &A, Array<T, SymmetricTensor<T, Descriptor>::n> const &B, Array<T, d> &x)
    {
        SymmetricRankThreeTensorImpl<T, d>::multWithRankTwoSymTensor(A, B, x);
    }
};

template <typename T, int d>
struct SymmetricRankFourTensorImpl {
};

template <typename T>
struct SymmetricRankFourTensorImpl<T, 2> {
    static const int n = 5;
    enum Indices { xxxx = 0, xxxy = 1, xxyy = 2, xyyy = 3, yyyy = 4 };

    static T contractIndexes(Array<T, n> const &A, Array<T, n> const &B)
    {
        return A[xxxx] * B[xxxx] + A[yyyy] * B[yyyy]
               + (T)4 * (A[xxxy] * B[xxxy] + A[xyyy] * B[xyyy]) + (T)6 * A[xxyy] * B[xxyy];
    }
    /// computes the contraction of the last two indexes of the tensor
    static void contractLastTwoIndexes(Array<T, n> const &tens, Array<T, 3> &res)
    {
        res[SymmetricTensorImpl<T, 2>::xx] = tens[xxxx] + tens[xxyy];
        res[SymmetricTensorImpl<T, 2>::xy] = tens[xxxy] + tens[xyyy];
        res[SymmetricTensorImpl<T, 2>::yy] = tens[xxyy] + tens[yyyy];
    }
    /// computes the interior product of symmetric rank 4 with rank three tensor
    static void multWithRankThreeSymTensor(
        Array<T, n> const &A, Array<T, SymmetricRankThreeTensorImpl<T, 2>::n> const &B,
        Array<T, 2> &x)
    {
        typedef SymmetricRankThreeTensorImpl<T, 2> SRT;

        x[0] = A[xxxx] * B[SRT::xxx] + (T)3 * (A[xxxy] * B[SRT::xxy] + A[xxyy] * B[SRT::xyy])
               + A[xyyy] * B[SRT::yyy];
        x[1] = A[xxxy] * B[SRT::xxx] + (T)3 * (A[xxyy] * B[SRT::xxy] + A[xyyy] * B[SRT::xyy])
               + A[yyyy] * B[SRT::yyy];
    }
};

template <typename T>
struct SymmetricRankFourTensorImpl<T, 3> {
    static const int n = 15;
    enum Indices {
        xxxx = 0,
        xxxy = 1,
        xxxz = 2,
        xxyy = 3,
        xxyz = 4,
        xxzz = 5,
        xyyy = 6,
        xyyz = 7,
        xyzz = 8,
        xzzz = 9,
        yyyy = 10,
        yyyz = 11,
        yyzz = 12,
        yzzz = 13,
        zzzz = 14
    };

    static T contractIndexes(Array<T, n> const &A, Array<T, n> const &B)
    {
        return A[xxxx] * B[xxxx] + A[yyyy] * B[yyyy] + A[zzzz] * B[zzzz]
               + (T)4
                     * (A[xxxy] * B[xxxy] + A[xxxz] * B[xxxz] + A[xyyy] * B[xyyy]
                        + A[xzzz] * B[xzzz] + A[yyyz] * B[yyyz] + A[yzzz] * B[yzzz])
               + (T)6 * (A[xxyy] * B[xxyy] + A[xxzz] * B[xxzz] + A[yyzz] * B[yyzz])
               + (T)12 * (A[xxyz] * B[xxyz] + A[xyyz] * B[xyyz] + A[xyzz] * B[xyzz]);
    }
    /// computes the contraction of the last two indexes of the tensor
    static void contractLastTwoIndexes(Array<T, n> const &tens, Array<T, 6> &res)
    {
        res[SymmetricTensorImpl<T, 3>::xx] = tens[xxxx] + tens[xxyy] + tens[xxzz];
        res[SymmetricTensorImpl<T, 3>::xy] = tens[xxxy] + tens[xyyy] + tens[xyzz];
        res[SymmetricTensorImpl<T, 3>::xz] = tens[xxxz] + tens[xyyz] + tens[xzzz];
        res[SymmetricTensorImpl<T, 3>::yy] = tens[xxyy] + tens[yyyy] + tens[yyzz];
        res[SymmetricTensorImpl<T, 3>::yz] = tens[xxyz] + tens[yyyz] + tens[yzzz];
        res[SymmetricTensorImpl<T, 3>::zz] = tens[xxzz] + tens[yyzz] + tens[zzzz];
    }
    static void multWithRankThreeSymTensor(
        Array<T, n> const &A, Array<T, SymmetricRankThreeTensorImpl<T, 3>::n> const &B,
        Array<T, 3> &x)
    {
        typedef SymmetricRankThreeTensorImpl<T, 3> SRT;

        x[0] = A[xxxx] * B[SRT::xxx] + A[xyyy] * B[SRT::yyy] + A[xzzz] * B[SRT::zzz]
               + (T)3
                     * (A[xxxy] * B[SRT::xxy] + A[xxyy] * B[SRT::xyy] + A[xxxz] * B[SRT::xxz]
                        + A[xyyz] * B[SRT::yyz] + A[xxzz] * B[SRT::xzz] + A[xyzz] * B[SRT::yzz])
               + (T)6 * A[xxyz] * B[SRT::xyz];

        x[1] = A[yyyy] * B[SRT::yyy] + A[xxxy] * B[SRT::xxx] + A[yzzz] * B[SRT::zzz]
               + (T)3
                     * (A[xxyy] * B[SRT::xxy] + A[xyyy] * B[SRT::xyy] + A[yyyz] * B[SRT::yyz]
                        + A[yyzz] * B[SRT::yzz] + A[xyzz] * B[SRT::xzz] + A[xxyz] * B[SRT::xxz])
               + (T)6 * A[xyyz] * B[SRT::xyz];

        x[2] = A[zzzz] * B[SRT::zzz] + A[xxxz] * B[SRT::xxx] + A[yyyz] * B[SRT::yyy]
               + (T)3
                     * (A[xxyz] * B[SRT::xxy] + A[xyyz] * B[SRT::xyy] + A[xxzz] * B[SRT::xxz]
                        + A[yyzz] * B[SRT::yyz] + A[xzzz] * B[SRT::xzz] + A[yzzz] * B[SRT::yzz])
               + (T)6 * A[xyzz] * B[SRT::xyz];
    }
};

/// Operations on a symmetric tensor which stores only above-or-on-diagonal values
template <typename T, template <typename U> class Descriptor>
struct SymmetricRankFourTensor {
    /// Number of dimensions for current lattice
    static const int d = Descriptor<T>::d;
    /// Number of elements (reduced by symmetry)
    static const int n = SymmetricRankFourTensorImpl<T, d>::n;

    /// computes the contraction of the last two indexes of the tensor
    static T contractIndexes(Array<T, n> const &A, Array<T, n> const &B)
    {
        return SymmetricRankFourTensorImpl<T, d>::contractIndexes(A, B);
    }
    /// computes the contraction of the last two indexes of the tensor
    static void contractLastTwoIndexes(
        Array<T, n> const &tens, Array<T, SymmetricTensor<T, Descriptor>::n> &res)
    {
        SymmetricRankFourTensorImpl<T, d>::contractLastTwoIndexes(tens, res);
    }
    static void multWithRankThreeSymTensor(
        Array<T, n> const &A, Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> const &B,
        Array<T, d> &x)
    {
        SymmetricRankFourTensorImpl<T, d>::multWithRankThreeSymTensor(A, B, x);
    }
};

template <typename T, int d>
struct SymmetricRankFiveTensorImpl {
};

template <typename T>
struct SymmetricRankFiveTensorImpl<T, 2> {
    static const int n = 6;
    enum Indices { xxxxx = 0, xxxxy = 1, xxxyy = 2, xxyyy = 3, xyyyy = 4, yyyyy = 5 };

    // static T contractIndexes(Array<T,n> const& A, Array<T,n> const& B) {
    //     return A[xxxxx]*B[xxxxx] + A[yyyy]*B[yyyy] + (T)4*(A[xxxy]*B[xxxy] + A[xyyy]*B[xyyy]) +
    //     (T)6*A[xxyy]*B[xxyy];
    // }
};

template <typename T>
struct SymmetricRankFiveTensorImpl<T, 3> {
    static const int n = 21;
    enum Indices {
        xxxxx = 0,
        xxxxy = 1,
        xxxxz = 2,
        xxxyy = 3,
        xxxyz = 4,
        xxxzz = 5,
        xxyyy = 6,
        xxyyz = 7,
        xxyzz = 8,
        xxzzz = 9,
        xyyyy = 10,
        xyyyz = 11,
        xyyzz = 12,
        xyzzz = 13,
        xzzzz = 14,
        yyyyy = 15,
        yyyyz = 16,
        yyyzz = 17,
        yyzzz = 18,
        yzzzz = 19,
        zzzzz = 20
    };

    // static T contractIndexes(Array<T,n> const& A, Array<T,n> const& B) {
    //     return A[xxxx]*B[xxxx] + A[yyyy]*B[yyyy] + A[zzzz]*B[zzzz]
    //         + (T)4*(A[xxxy]*B[xxxy] + A[xxxz]*B[xxxz] + A[xyyy]*B[xyyy] + A[xzzz]*B[xzzz] +
    //         A[yyyz]*B[yyyz] + A[yzzz]*B[yzzz])
    //         + (T)6*(A[xxyy]*B[xxyy] + A[xxzz]*B[xxzz] + A[yyzz]*B[yyzz])
    //         + (T)12*(A[xxyz]*B[xxyz] + A[xyyz]*B[xyyz] + A[xyzz]*B[xyzz]);
    // }
};

/// Operations on a symmetric tensor which stores only above-or-on-diagonal values
template <typename T, template <typename U> class Descriptor>
struct SymmetricRankFiveTensor {
    /// Number of dimensions for current lattice
    static const int d = Descriptor<T>::d;
    /// Number of elements (reduced by symmetry)
    static const int n = SymmetricRankFiveTensorImpl<T, d>::n;

    // /// computes the contraction of the last two indexes of the tensor
    // static T contractIndexes(Array<T,n> const& A, Array<T,n> const& B) {
    //     return SymmetricRankFiveTensorImpl<T,d>::contractIndexes(A, B);
    // }
};

template <typename T, int d>
struct SymmetricRankSixTensorImpl {
};

template <typename T>
struct SymmetricRankSixTensorImpl<T, 2> {
    static const int n = 7;
    enum Indices {
        xxxxxx = 0,
        xxxxxy = 1,
        xxxxyy = 2,
        xxxyyy = 3,
        xxyyyy = 4,
        xyyyyy = 5,
        yyyyyy = 6
    };

    // static T contractIndexes(Array<T,n> const& A, Array<T,n> const& B) {
    //     return A[xxxxx]*B[xxxxx] + A[yyyy]*B[yyyy] + (T)4*(A[xxxy]*B[xxxy] + A[xyyy]*B[xyyy]) +
    //     (T)6*A[xxyy]*B[xxyy];
    // }
};

template <typename T>
struct SymmetricRankSixTensorImpl<T, 3> {
    static const int n = 28;
    enum Indices {
        xxxxxx = 0,
        xxxxxy = 1,
        xxxxxz = 2,
        xxxxyy = 3,
        xxxxyz = 4,
        xxxxzz = 5,
        xxxyyy = 6,
        xxxyyz = 7,
        xxxyzz = 8,
        xxxzzz = 9,
        xxyyyy = 10,
        xxyyyz = 11,
        xxyyzz = 12,
        xxyzzz = 13,
        xxzzzz = 14,
        xyyyyy = 15,
        xyyyyz = 16,
        xyyyzz = 17,
        xyyzzz = 18,
        xyzzzz = 19,
        xzzzzz = 20,
        yyyyyy = 21,
        yyyyyz = 22,
        yyyyzz = 23,
        yyyzzz = 24,
        yyzzzz = 25,
        yzzzzz = 26,
        zzzzzz = 27
    };

    // static T contractIndexes(Array<T,n> const& A, Array<T,n> const& B) {
    //     return A[xxxx]*B[xxxx] + A[yyyy]*B[yyyy] + A[zzzz]*B[zzzz]
    //         + (T)4*(A[xxxy]*B[xxxy] + A[xxxz]*B[xxxz] + A[xyyy]*B[xyyy] + A[xzzz]*B[xzzz] +
    //         A[yyyz]*B[yyyz] + A[yzzz]*B[yzzz])
    //         + (T)6*(A[xxyy]*B[xxyy] + A[xxzz]*B[xxzz] + A[yyzz]*B[yyzz])
    //         + (T)12*(A[xxyz]*B[xxyz] + A[xyyz]*B[xyyz] + A[xyzz]*B[xyzz]);
    // }
};

/// Operations on a symmetric tensor which stores only above-or-on-diagonal values
template <typename T, template <typename U> class Descriptor>
struct SymmetricRankSixTensor {
    /// Number of dimensions for current lattice
    static const int d = Descriptor<T>::d;
    /// Number of elements (reduced by symmetry)
    static const int n = SymmetricRankSixTensorImpl<T, d>::n;

    // /// computes the contraction of the last two indexes of the tensor
    // static T contractIndexes(Array<T,n> const& A, Array<T,n> const& B) {
    //     return SymmetricRankSixTensorImpl<T,d>::contractIndexes(A, B);
    // }
};

template <typename T>
void crossProduct(Array<T, 3> const &u1, Array<T, 3> const &u2, Array<T, 3> &result)
{
    result[0] = u1[1] * u2[2] - u1[2] * u2[1];
    result[1] = u1[2] * u2[0] - u1[0] * u2[2];
    result[2] = u1[0] * u2[1] - u1[1] * u2[0];
}

template <typename T>
Array<T, 3> crossProduct(Array<T, 3> const &u1, Array<T, 3> const &u2)
{
    Array<T, 3> result;
    crossProduct(u1, u2, result);

    return result;
}

/// Scalar product between two vectors.
template <typename T, pluint n>
T dot(Array<T, n> const &v1, Array<T, n> const &v2)
{
    return VectorTemplateImpl<T, n>::scalarProduct(v1, v2);
}

/// Scalar product between a c-array and a vector.
template <typename T, pluint n>
T dot(const T v1[], Array<T, n> const &v2)
{
    return VectorTemplateImpl<T, n>::scalarProduct(v1, v2);
}
/// Scalar product between an explicitly int-based c-array and a vector.
template <typename T, pluint n>
T intDot(const int v1[], Array<T, n> const &v2)
{
    return VectorTemplateImpl<T, n>::intScalarProduct(v1, v2);
}

template <typename T, pluint n>
T normSqr(Array<T, n> const &v)
{
    return VectorTemplateImpl<T, n>::normSqr(v);
}

template <typename T, pluint n>
T norm(Array<T, n> const &v)
{
    return std::sqrt(normSqr<T, n>(v));
}

/// Compute the normal vector for a given triangle. If "isAreaWeighted" is false,
///   then the normal has length equal to one, or equal to zero if the triangle
///   has a zero area. If "isAreaWeighted" is true, then the normal has length
///   equal to twice the area of the triangle.
template <typename T>
Array<T, 3> computeTriangleNormal(
    Array<T, 3> const &v0, Array<T, 3> const &v1, Array<T, 3> const &v2, bool isAreaWeighted)
{
    Array<T, 3> e01 = v1 - v0;
    Array<T, 3> e02 = v2 - v0;

    Array<T, 3> n;
    crossProduct(e01, e02, n);
    if (!isAreaWeighted) {
        T normN = norm(n);
        if (!util::isZero(normN)) {
            n /= normN;
        } else {
            n.resetToZero();
        }
    }

    return n;
}

template <typename T>
T computeTriangleArea(Array<T, 3> const &v0, Array<T, 3> const &v1, Array<T, 3> const &v2)
{
    Array<T, 3> e01 = v1 - v0;
    Array<T, 3> e02 = v2 - v0;
    Array<T, 3> cross;
    crossProduct(e01, e02, cross);

    return (T)0.5 * norm(cross);
}

template <typename T>
T computeTriangleArea(T *v0, T *v1, T *v2)
{
    Array<T, 3> a0;
    a0.from_cArray(v0);
    Array<T, 3> a1;
    a1.from_cArray(v1);
    Array<T, 3> a2;
    a2.from_cArray(v2);
    return computeTriangleArea(a0, a1, a2);
}

template <typename T>
void computeTriangleAreaAndUnitNormal(
    Array<T, 3> const &v0, Array<T, 3> const &v1, Array<T, 3> const &v2, T &area,
    Array<T, 3> &unitNormal)
{
    Array<T, 3> e01 = v1 - v0;
    Array<T, 3> e02 = v2 - v0;

    crossProduct(e01, e02, unitNormal);
    T normN = norm(unitNormal);
    if (!util::isZero(normN)) {
        area = (T)0.5 * normN;
        unitNormal /= normN;
    } else {
        area = (T)0;
        unitNormal.resetToZero();
    }
}

template <typename T>
void computeTriangleNormalTimesArea(
    Array<T, 3> const &v0, Array<T, 3> const &v1, Array<T, 3> const &v2,
    Array<T, 3> &normalTimesArea)
{
    Array<T, 3> e01 = v1 - v0;
    Array<T, 3> e02 = v2 - v0;
    crossProduct(e01, e02, normalTimesArea);
    normalTimesArea *= (T)0.5;
}

template <typename T>
T computeTetrahedronSignedVolume(
    Array<T, 3> const &v0, Array<T, 3> const &v1, Array<T, 3> const &v2, Array<T, 3> const &v3)
{
    // The position arrays specify two oriented triangles, namely
    // (0, 1, 2) and (3, 2, 1). These triangles share
    // the common edge 1-2.

    Array<T, 3> e10 = v0 - v1;
    Array<T, 3> e12 = v2 - v1;
    Array<T, 3> e13 = v3 - v1;
    Array<T, 3> cross;
    crossProduct(e12, e13, cross);

    return dot(e10, cross) / 6.0;
}

/// Given two vectors this function computes the angle between them
///   assuming that they originate at the same point in space
template <typename T>
T angleBetweenVectors(Array<T, 3> const &v1, Array<T, 3> const &v2)
{
    Array<T, 3> cross;
    crossProduct(v1, v2, cross);
    return std::atan2(norm(cross), dot(v1, v2));
}

template <typename T>
Array<T, 3> rotateAtOrigin(Array<T, 3> const &p, Array<T, 3> const &normedAxis, T theta)
{
    // Standard (initial) implementation.
    // (This version is the slowest).

    /*
    Array<T,3> const& u = normedAxis;
    T d = std::sqrt(u[1]*u[1] + u[2]*u[2]);

    Array<T,3> q1 = p;
    Array<T,3> q2 = q1;

    T eps = 1.e-2;
    //T eps = (T)100.*std::numeric_limits<T>::epsilon();
    // Rotate about the x-axis to be in xz-plane.
    if (!util::isZero(d, eps)) {
        q2[0] = q1[0];
        q2[1] = q1[1] * u[2] / d - q1[2] * u[1] / d;
        q2[2] = q1[1] * u[1] / d + q1[2] * u[2] / d;
    }

    // Rotate about the y-axis to fall on z-axis.
    q1[0] = q2[0] * d - q2[2] * u[0];
    q1[1] = q2[1];
    q1[2] = q2[0] * u[0] + q2[2] * d;

    // Perform desired rotation.
    T ct = std::cos(theta);
    T st = std::sin(theta);
    q2[0] = q1[0] * ct - q1[1] * st;
    q2[1] = q1[0] * st + q1[1] * ct;
    q2[2] = q1[2];

    // Rotate backward around y-axis.
    q1[0] =   q2[0] * d + q2[2] * u[0];
    q1[1] =   q2[1];
    q1[2] = - q2[0] * u[0] + q2[2] * d;

    q2 = q1;
    // Rotate backward around x-axis.
    if (!util::isZero(d, eps)) {
        q2[0] =   q1[0];
        q2[1] =   q1[1] * u[2] / d + q1[2] * u[1] / d;
        q2[2] = - q1[1] * u[1] / d + q1[2] * u[2] / d;
    }

    return q2;
    */

    // Optimized implementation.
    // (At the time of writing, this implementation is faster than the
    //  next one when maximum compiler optimizations are used).

    Array<T, 3> const &e = normedAxis;
    T c = std::cos(theta);
    T s = std::sin(theta);
    T omc = (T)1 - c;

    T R[3][3];
    R[0][0] = omc * e[0] * e[0] + c;
    R[0][1] = omc * e[0] * e[1] - s * e[2];
    R[0][2] = omc * e[0] * e[2] + s * e[1];
    R[1][0] = omc * e[1] * e[0] + s * e[2];
    R[1][1] = omc * e[1] * e[1] + c;
    R[1][2] = omc * e[1] * e[2] - s * e[0];
    R[2][0] = omc * e[2] * e[0] - s * e[1];
    R[2][1] = omc * e[2] * e[1] + s * e[0];
    R[2][2] = omc * e[2] * e[2] + c;

    Array<T, 3> q;
    q[0] = R[0][0] * p[0] + R[0][1] * p[1] + R[0][2] * p[2];
    q[1] = R[1][0] * p[0] + R[1][1] * p[1] + R[1][2] * p[2];
    q[2] = R[2][0] * p[0] + R[2][1] * p[1] + R[2][2] * p[2];

    return q;

    // Another optimized implementation.
    // (This version is the fastest of the three when no compiler optimizations are used).

    /*
    Array<T,3> const& e = normedAxis;
    T c = std::cos(theta);
    T s = std::sin(theta);
    T f = ((T) 1 - c) * (e[0] * p[0] + e[1] * p[1] + e[2] * p[2]);

    Array<T,3> q;
    q[0] = c * p[0] + s * (e[1] * p[2] - e[2] * p[1]) + f * e[0];
    q[1] = c * p[1] + s * (e[2] * p[0] - e[0] * p[2]) + f * e[1];
    q[2] = c * p[2] + s * (e[0] * p[1] - e[1] * p[0]) + f * e[2];

    return q;
    */
}

// This function rotates p with "proper Euler angles" and "intrinsic" rotation zx'z''
// * phi: represents a rotation around the z axis,
// * theta: represents a rotation around the x′ axis,
// * psi: represents a rotation around the z″ axis.
template <typename T>
Array<T, 3> rotateWithEulerAngles(Array<T, 3> const &p, T phi, T theta, T psi)
{
#ifdef PLB_DEBUG
    T pi = std::acos((T)-1.0);
#endif
    PLB_ASSERT(util::greaterEqual(theta, (T)0) && util::lessEqual(theta, pi));

    T cosPhi = std::cos(phi);
    T sinPhi = std::sin(phi);
    T cosTheta = std::cos(theta);
    T sinTheta = std::sin(theta);
    T cosPsi = std::cos(psi);
    T sinPsi = std::sin(psi);

    T a[3][3];
    a[0][0] = (T)1.0;
    a[0][1] = (T)0.0;
    a[0][2] = (T)0.0;
    a[1][0] = (T)0.0;
    a[1][1] = cosTheta;
    a[1][2] = -sinTheta;
    a[2][0] = (T)0.0;
    a[2][1] = sinTheta;
    a[2][2] = cosTheta;

    T b[3][3];
    b[0][0] = cosPhi;
    b[0][1] = -sinPhi;
    b[0][2] = (T)0.0;
    b[1][0] = sinPhi;
    b[1][1] = cosPhi;
    b[1][2] = (T)0.0;
    b[2][0] = (T)0.0;
    b[2][1] = (T)0.0;
    b[2][2] = (T)1.0;

    T c[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            c[i][j] = (T)0.0;
            for (int k = 0; k < 3; k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }

    b[0][0] = cosPsi;
    b[0][1] = -sinPsi;
    b[0][2] = (T)0.0;
    b[1][0] = sinPsi;
    b[1][1] = cosPsi;
    b[1][2] = (T)0.0;
    b[2][0] = (T)0.0;
    b[2][1] = (T)0.0;
    b[2][2] = (T)1.0;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            a[i][j] = (T)0.0;
            for (int k = 0; k < 3; k++) {
                a[i][j] += b[i][k] * c[k][j];
            }
        }
    }

    Array<T, 3> x;
    for (int i = 0; i < 3; i++) {
        x[i] = (T)0.0;
        for (int j = 0; j < 3; j++) {
            x[i] += a[i][j] * p[j];
        }
    }

    return x;
}

// This function returns the new position of the point "oldPosition", after a rotation with with an
// angle "theta". It needs the angle of the rotation, a normalized vector parallel to the axis of
// rotation, and a point on the axis of rotation.
template <typename T>
Array<T, 3> getRotatedPosition(
    Array<T, 3> const &oldPosition, T theta, Array<T, 3> const &rotationAxisUnitVector,
    Array<T, 3> const &pointOnRotationAxis)
{
    // Standard optimized implementation.

    return (
        pointOnRotationAxis
        + rotateAtOrigin(oldPosition - pointOnRotationAxis, rotationAxisUnitVector, theta));
}

// This function returns the new position of the point "oldPosition", after a discrete rotation
// with time step equal to 1.0. It needs the angular velocity of the rotation, a normalized vector
// parallel to the axis of rotation, and a point on the axis of rotation.
// (The "rotationAxisUnitVector" is provided for optimization purposes).
template <typename T>
Array<T, 3> getRotatedPosition(
    Array<T, 3> const &oldPosition, Array<T, 3> const &angularVelocity,
    Array<T, 3> const &rotationAxisUnitVector, Array<T, 3> const &pointOnRotationAxis)
{
    return (getRotatedPosition(
        oldPosition, dot(angularVelocity, rotationAxisUnitVector), rotationAxisUnitVector,
        pointOnRotationAxis));
}

// This function returns the "discrete" rotational velocity of the point "oldPosition", considering
// a rotation with an angle "theta". It needs the angle of the rotation, a normalized vector
// parallel to the axis of rotation, and a point on the axis of rotation.
template <typename T>
Array<T, 3> getDiscreteRotationalVelocity(
    Array<T, 3> const &oldPosition, T theta, Array<T, 3> const &rotationAxisUnitVector,
    Array<T, 3> const &pointOnRotationAxis)
{
    Array<T, 3> newPosition =
        getRotatedPosition(oldPosition, theta, rotationAxisUnitVector, pointOnRotationAxis);
    Array<T, 3> velocity = newPosition - oldPosition;
    return (velocity);
}

// This function returns the "discrete" rotational velocity of the point "oldPosition", considering
// a time step equal to 1.0. It needs the angular velocity of the rotation, a normalized vector
// parallel to the axis of rotation, and a point on the axis of rotation.
// (The rotationAxisUnitVector is provided for optimization purposes).
template <typename T>
Array<T, 3> getDiscreteRotationalVelocity(
    Array<T, 3> const &oldPosition, Array<T, 3> const &angularVelocity,
    Array<T, 3> const &rotationAxisUnitVector, Array<T, 3> const &pointOnRotationAxis)
{
    return (getDiscreteRotationalVelocity(
        oldPosition, dot(angularVelocity, rotationAxisUnitVector), rotationAxisUnitVector,
        pointOnRotationAxis));
}

// This function returns the "exact" rotational velocity of the point "position". It needs the
// angular velocity of the rotation, and a point on the axis of rotation.
template <typename T>
Array<T, 3> getExactRotationalVelocity(
    Array<T, 3> const &position, Array<T, 3> const &angularVelocity,
    Array<T, 3> const &pointOnRotationAxis)
{
    return (crossProduct(angularVelocity, position - pointOnRotationAxis));
}

}  // namespace plb

#endif
