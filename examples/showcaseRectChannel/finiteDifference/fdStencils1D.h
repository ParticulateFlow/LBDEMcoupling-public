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

#ifndef FD_STENCILS_1D_H
#define FD_STENCILS_1D_H

#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField2D.h"
#include "atomicBlock/dataField3D.h"
#include "core/globalDefs.h"

namespace plb {

namespace fd {

/// Second-order central gradient (u_p1 = u(x+1))
template <typename T>
T ctl_diff(T u_p1, T u_m1)
{
    return (u_p1 - u_m1) / (T)2;
}

/// Fourth-order central gradient (u_p1 = u(x+1))
template <typename T>
T ctl_diff(T u_p2, T u_p1, T u_m1, T u_m2)
{
    return ((T)1 / (T)12) * (u_m2 - u_p2) + ((T)2 / (T)3) * (u_p1 - u_m1);
}

/// Sixth-order central gradient (u_p1 = u(x+1))
template <typename T>
T ctl_diff(T u_p3, T u_p2, T u_p1, T u_m1, T u_m2, T u_m3)
{
    return -(-u_p3 + (T)9 * u_p2 - (T)45 * u_p1 + (T)45 * u_m1 - (T)9 * u_m2 + u_m3) / (T)60;
}

/// Sixth-order central gradient (u_p1 = u(x+1))
template <typename T>
T ctl_diff(T u_p4, T u_p3, T u_p2, T u_p1, T u_m1, T u_m2, T u_m3, T u_m4)
{
    return -((T)3 * u_p4 - (T)32 * u_p3 + (T)168 * u_p2 - (T)672 * u_p1 + (T)672 * u_m1
             - (T)168 * u_m2 + (T)32 * u_m3 - (T)3 * u_m4)
           / (T)840;
}

/// Second-order forward gradient (u_1 = u(x+1))
template <typename T>
T fwd_diff(T u_0, T u_1, T u_2)
{
    return (-(T)3 * u_0 + (T)4 * u_1 - (T)1 * u_2) / (T)2;
}

/// Second-order backward gradient (u_1 = u(x-1))
template <typename T>
T bwd_diff(T u_0, T u_1, T u_2)
{
    return -fwd_diff(u_0, u_1, u_2);
}

/// First-order forward gradient (u_1 = u(x+1))
template <typename T>
T o1_fwd_diff(T u_0, T u_1)
{
    return (-u_0 + u_1);
}

/// Value at u_0 for which asymmetric gradient is zero (u_1 = u(x+1))
template <typename T>
T boundaryZeroGradient(T u_1, T u_2)
{
    return (T)4 / (T)3 * u_1 - (T)1 / (T)3 * u_2;
}

/// Linear interpolation (yields u0 at pos=0 and u1 at pos=1)
template <typename T>
T linearInterpolate(T u_0, T u_1, T pos)
{
    return ((T)1 - pos) * u_0 + pos * u_1;
}

/// Third order corner gradient (a_p1=a(x+1/2), a_p2=a(x+3/2)), returns a(x)
template <typename T>
inline T cornerThreePointsInterp(T a_m1, T a_p1, T a_p2)
{
    //   o -x- o --- o //
    // a_m1    a_p1  a_p2, with x the point that is interpolated
    return 0.375 * a_m1 + (T)0.75 * a_p1 - (T)0.125 * a_p2;
}

/// Fourth order central gradient (a_p1=a(x+1/2), a_p2=a(x+3/2)), returns a(x)
template <typename T>
inline T centralFourPointsInterp(T a_m2, T a_m1, T a_p1, T a_p2)
{
    //   o --- o -x- o --- o //
    // a_m1    a    a_p1  a_p2, with x the point that is interpolated
    return -(T)0.0625 * (a_m2 + a_p2) + (T)0.5625 * (a_m1 + a_p1);
}

}  // namespace fd

/// Finite Difference operations on data-fields
namespace fdDataField {

template <typename T>
inline T bulkXderiv(ScalarField3D<T> const &field, plint iX, plint iY, plint iZ)
{
    T dxu = fd::ctl_diff(field.get(iX + 1, iY, iZ), field.get(iX - 1, iY, iZ));
    return dxu;
}

template <typename T>
inline T bulkYderiv(ScalarField3D<T> const &field, plint iX, plint iY, plint iZ)
{
    T dyu = fd::ctl_diff(field.get(iX, iY + 1, iZ), field.get(iX, iY - 1, iZ));
    return dyu;
}

template <typename T>
inline T bulkZderiv(ScalarField3D<T> const &field, plint iX, plint iY, plint iZ)
{
    T dzu = fd::ctl_diff(field.get(iX, iY, iZ + 1), field.get(iX, iY, iZ - 1));
    return dzu;
}

template <typename T>
inline T planeXderiv(
    ScalarField3D<T> const &field, int direction, int orientation, plint iX, plint iY, plint iZ)
{
    if (direction == 0) {
        return -orientation
               * fd::o1_fwd_diff(field.get(iX, iY, iZ), field.get(iX - 1 * orientation, iY, iZ));
    } else {
        return bulkXderiv(field, iX, iY, iZ);
    }
}

template <typename T>
inline T planeYderiv(
    ScalarField3D<T> const &field, int direction, int orientation, plint iX, plint iY, plint iZ)
{
    if (direction == 1) {
        return -orientation
               * fd::o1_fwd_diff(field.get(iX, iY, iZ), field.get(iX, iY - 1 * orientation, iZ));
    } else {
        return bulkYderiv(field, iX, iY, iZ);
    }
}

template <typename T>
inline T planeZderiv(
    ScalarField3D<T> const &field, int direction, int orientation, plint iX, plint iY, plint iZ)
{
    if (direction == 2) {
        return -orientation
               * fd::o1_fwd_diff(field.get(iX, iY, iZ), field.get(iX, iY, iZ - 1 * orientation));
    } else {
        return bulkZderiv(field, iX, iY, iZ);
    }
}

template <typename T>
inline T edgeXderiv(
    ScalarField3D<T> const &field, int plane, int direction1, int direction2, plint iX, plint iY,
    plint iZ)
{
    if (plane == 0) {
        return bulkXderiv(field, iX, iY, iZ);
    } else {
        int orientation = plane == 1 ? direction2 : direction1;
        return -orientation
               * fd::o1_fwd_diff(field.get(iX, iY, iZ), field.get(iX - 1 * orientation, iY, iZ));
    }
}

template <typename T>
inline T edgeYderiv(
    ScalarField3D<T> const &field, int plane, int direction1, int direction2, plint iX, plint iY,
    plint iZ)
{
    if (plane == 1) {
        return bulkYderiv(field, iX, iY, iZ);
    } else {
        int orientation = plane == 0 ? direction1 : direction2;
        return -orientation
               * fd::o1_fwd_diff(field.get(iX, iY, iZ), field.get(iX, iY - 1 * orientation, iZ));
    }
}

template <typename T>
inline T edgeZderiv(
    ScalarField3D<T> const &field, int plane, int direction1, int direction2, plint iX, plint iY,
    plint iZ)
{
    if (plane == 2) {
        return bulkZderiv(field, iX, iY, iZ);
    } else {
        int orientation = plane == 0 ? direction2 : direction1;
        return -orientation
               * fd::o1_fwd_diff(field.get(iX, iY, iZ), field.get(iX, iY, iZ - 1 * orientation));
    }
}

template <typename T>
inline T cornerXderiv(
    ScalarField3D<T> const &field, int normalX, int normalY, int normalZ, plint iX, plint iY,
    plint iZ)
{
    int orientation = normalX;
    return -orientation
           * fd::o1_fwd_diff(field.get(iX, iY, iZ), field.get(iX - 1 * orientation, iY, iZ));
}

template <typename T>
inline T cornerYderiv(
    ScalarField3D<T> const &field, int normalX, int normalY, int normalZ, plint iX, plint iY,
    plint iZ)
{
    int orientation = normalY;
    return -orientation
           * fd::o1_fwd_diff(field.get(iX, iY, iZ), field.get(iX, iY - 1 * orientation, iZ));
}

template <typename T>
inline T cornerZderiv(
    ScalarField3D<T> const &field, int normalX, int normalY, int normalZ, plint iX, plint iY,
    plint iZ)
{
    int orientation = normalZ;
    return -orientation
           * fd::o1_fwd_diff(field.get(iX, iY, iZ), field.get(iX, iY, iZ - 1 * orientation));
}

template <typename T, int nDim>
inline T bulkXderiv(TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ, int iD)
{
    T dxu = fd::ctl_diff(velocity.get(iX + 1, iY, iZ)[iD], velocity.get(iX - 1, iY, iZ)[iD]);
    return dxu;
}

template <typename T, int nDim>
inline T bulkYderiv(TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ, int iD)
{
    T dyu = fd::ctl_diff(velocity.get(iX, iY + 1, iZ)[iD], velocity.get(iX, iY - 1, iZ)[iD]);
    return dyu;
}

template <typename T, int nDim>
inline T bulkZderiv(TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ, int iD)
{
    T dzu = fd::ctl_diff(velocity.get(iX, iY, iZ + 1)[iD], velocity.get(iX, iY, iZ - 1)[iD]);
    return dzu;
}

template <typename T, int nDim>
inline T bulkXderivOrderFour(
    TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ, int iD)
{
    T dxu = fd::ctl_diff(
        velocity.get(iX + 2, iY, iZ)[iD], velocity.get(iX + 1, iY, iZ)[iD],
        velocity.get(iX - 1, iY, iZ)[iD], velocity.get(iX - 2, iY, iZ)[iD]);
    return dxu;
}

template <typename T, int nDim>
inline T bulkYderivOrderFour(
    TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ, int iD)
{
    T dyu = fd::ctl_diff(
        velocity.get(iX, iY + 2, iZ)[iD], velocity.get(iX, iY + 1, iZ)[iD],
        velocity.get(iX, iY - 1, iZ)[iD], velocity.get(iX, iY - 2, iZ)[iD]);
    return dyu;
}

template <typename T, int nDim>
inline T bulkZderivOrderFour(
    TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ, int iD)
{
    T dzu = fd::ctl_diff(
        velocity.get(iX, iY, iZ + 2)[iD], velocity.get(iX, iY, iZ + 1)[iD],
        velocity.get(iX, iY, iZ - 1)[iD], velocity.get(iX, iY, iZ - 2)[iD]);
    return dzu;
}

template <typename T, int nDim>
inline T bulkXderivOrderSix(
    TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ, int iD)
{
    T dxu = fd::ctl_diff(
        velocity.get(iX + 3, iY, iZ)[iD], velocity.get(iX + 2, iY, iZ)[iD],
        velocity.get(iX + 1, iY, iZ)[iD], velocity.get(iX - 1, iY, iZ)[iD],
        velocity.get(iX - 2, iY, iZ)[iD], velocity.get(iX - 3, iY, iZ)[iD]);
    return dxu;
}

template <typename T, int nDim>
inline T bulkYderivOrderSix(
    TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ, int iD)
{
    T dyu = fd::ctl_diff(
        velocity.get(iX, iY + 3, iZ)[iD], velocity.get(iX, iY + 2, iZ)[iD],
        velocity.get(iX, iY + 1, iZ)[iD], velocity.get(iX, iY - 1, iZ)[iD],
        velocity.get(iX, iY - 2, iZ)[iD], velocity.get(iX, iY - 3, iZ)[iD]);
    return dyu;
}

template <typename T, int nDim>
inline T bulkZderivOrderSix(
    TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ, int iD)
{
    T dzu = fd::ctl_diff(
        velocity.get(iX, iY, iZ + 3)[iD], velocity.get(iX, iY, iZ + 2)[iD],
        velocity.get(iX, iY, iZ + 1)[iD], velocity.get(iX, iY, iZ - 1)[iD],
        velocity.get(iX, iY, iZ - 2)[iD], velocity.get(iX, iY, iZ - 3)[iD]);
    return dzu;
}

template <typename T, int nDim>
inline T bulkXderivOrderEight(
    TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ, int iD)
{
    T dxu = fd::ctl_diff(
        velocity.get(iX + 4, iY, iZ)[iD], velocity.get(iX + 3, iY, iZ)[iD],
        velocity.get(iX + 2, iY, iZ)[iD], velocity.get(iX + 1, iY, iZ)[iD],
        velocity.get(iX - 1, iY, iZ)[iD], velocity.get(iX - 2, iY, iZ)[iD],
        velocity.get(iX - 3, iY, iZ)[iD], velocity.get(iX - 4, iY, iZ)[iD]);
    return dxu;
}

template <typename T, int nDim>
inline T bulkYderivOrderEight(
    TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ, int iD)
{
    T dyu = fd::ctl_diff(
        velocity.get(iX, iY + 4, iZ)[iD], velocity.get(iX, iY + 3, iZ)[iD],
        velocity.get(iX, iY + 2, iZ)[iD], velocity.get(iX, iY + 1, iZ)[iD],
        velocity.get(iX, iY - 1, iZ)[iD], velocity.get(iX, iY - 2, iZ)[iD],
        velocity.get(iX, iY - 3, iZ)[iD], velocity.get(iX, iY - 4, iZ)[iD]);
    return dyu;
}

template <typename T, int nDim>
inline T bulkZderivOrderEight(
    TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ, int iD)
{
    T dzu = fd::ctl_diff(
        velocity.get(iX, iY, iZ + 4)[iD], velocity.get(iX, iY, iZ + 3)[iD],
        velocity.get(iX, iY, iZ + 2)[iD], velocity.get(iX, iY, iZ + 1)[iD],
        velocity.get(iX, iY, iZ - 1)[iD], velocity.get(iX, iY, iZ - 2)[iD],
        velocity.get(iX, iY, iZ - 3)[iD], velocity.get(iX, iY, iZ - 4)[iD]);
    return dzu;
}

template <typename T, int nDim>
inline T planeXderiv(
    TensorField3D<T, nDim> const &velocity, int direction, int orientation, plint iX, plint iY,
    plint iZ, int iD)
{
    if (direction == 0) {
        return -orientation
               * fd::o1_fwd_diff(
                   velocity.get(iX, iY, iZ)[iD], velocity.get(iX - 1 * orientation, iY, iZ)[iD]);
    } else {
        return bulkXderiv(velocity, iX, iY, iZ, iD);
    }
}

template <typename T, int nDim>
inline T planeYderiv(
    TensorField3D<T, nDim> const &velocity, int direction, int orientation, plint iX, plint iY,
    plint iZ, int iD)
{
    if (direction == 1) {
        return -orientation
               * fd::o1_fwd_diff(
                   velocity.get(iX, iY, iZ)[iD], velocity.get(iX, iY - 1 * orientation, iZ)[iD]);
    } else {
        return bulkYderiv(velocity, iX, iY, iZ, iD);
    }
}

template <typename T, int nDim>
inline T planeZderiv(
    TensorField3D<T, nDim> const &velocity, int direction, int orientation, plint iX, plint iY,
    plint iZ, int iD)
{
    if (direction == 2) {
        return -orientation
               * fd::o1_fwd_diff(
                   velocity.get(iX, iY, iZ)[iD], velocity.get(iX, iY, iZ - 1 * orientation)[iD]);
    } else {
        return bulkZderiv(velocity, iX, iY, iZ, iD);
    }
}

template <typename T, int nDim>
inline T edgeXderiv(
    TensorField3D<T, nDim> const &velocity, int plane, int direction1, int direction2, plint iX,
    plint iY, plint iZ, int iD)
{
    if (plane == 0) {
        return bulkXderiv(velocity, iX, iY, iZ, iD);
    } else {
        int orientation = plane == 1 ? direction2 : direction1;
        return -orientation
               * fd::o1_fwd_diff(
                   velocity.get(iX, iY, iZ)[iD], velocity.get(iX - 1 * orientation, iY, iZ)[iD]);
    }
}

template <typename T, int nDim>
inline T edgeYderiv(
    TensorField3D<T, nDim> const &velocity, int plane, int direction1, int direction2, plint iX,
    plint iY, plint iZ, int iD)
{
    if (plane == 1) {
        return bulkYderiv(velocity, iX, iY, iZ, iD);
    } else {
        int orientation = plane == 0 ? direction1 : direction2;
        return -orientation
               * fd::o1_fwd_diff(
                   velocity.get(iX, iY, iZ)[iD], velocity.get(iX, iY - 1 * orientation, iZ)[iD]);
    }
}

template <typename T, int nDim>
inline T edgeZderiv(
    TensorField3D<T, nDim> const &velocity, int plane, int direction1, int direction2, plint iX,
    plint iY, plint iZ, int iD)
{
    if (plane == 2) {
        return bulkZderiv(velocity, iX, iY, iZ, iD);
    } else {
        int orientation = plane == 0 ? direction2 : direction1;
        return -orientation
               * fd::o1_fwd_diff(
                   velocity.get(iX, iY, iZ)[iD], velocity.get(iX, iY, iZ - 1 * orientation)[iD]);
    }
}

template <typename T, int nDim>
inline T cornerXderiv(
    TensorField3D<T, nDim> const &velocity, int normalX, int normalY, int normalZ, plint iX,
    plint iY, plint iZ, int iD)
{
    int orientation = normalX;
    return -orientation
           * fd::o1_fwd_diff(
               velocity.get(iX, iY, iZ)[iD], velocity.get(iX - 1 * orientation, iY, iZ)[iD]);
}

template <typename T, int nDim>
inline T cornerYderiv(
    TensorField3D<T, nDim> const &velocity, int normalX, int normalY, int normalZ, plint iX,
    plint iY, plint iZ, int iD)
{
    int orientation = normalY;
    return -orientation
           * fd::o1_fwd_diff(
               velocity.get(iX, iY, iZ)[iD], velocity.get(iX, iY - 1 * orientation, iZ)[iD]);
}

template <typename T, int nDim>
inline T cornerZderiv(
    TensorField3D<T, nDim> const &velocity, int normalX, int normalY, int normalZ, plint iX,
    plint iY, plint iZ, int iD)
{
    int orientation = normalZ;
    return -orientation
           * fd::o1_fwd_diff(
               velocity.get(iX, iY, iZ)[iD], velocity.get(iX, iY, iZ - 1 * orientation)[iD]);
}

template <typename T, int nDim>
inline T bulkVorticityX(TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ)
{
    T dyuz = fdDataField::bulkYderiv(velocity, iX, iY, iZ, 2);
    T dzuy = fdDataField::bulkZderiv(velocity, iX, iY, iZ, 1);

    return dyuz - dzuy;
}

template <typename T, int nDim>
inline T bulkVorticityY(TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ)
{
    T dzux = fdDataField::bulkZderiv(velocity, iX, iY, iZ, 0);
    T dxuz = fdDataField::bulkXderiv(velocity, iX, iY, iZ, 2);

    return dzux - dxuz;
}

template <typename T, int nDim>
inline T bulkVorticityZ(TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ)
{
    T dxuy = fdDataField::bulkXderiv(velocity, iX, iY, iZ, 1);
    T dyux = fdDataField::bulkYderiv(velocity, iX, iY, iZ, 0);

    return dxuy - dyux;
}

template <typename T, int nDim>
inline T bulkVorticityXOrderFour(
    TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ)
{
    T dyuz = fdDataField::bulkYderivOrderFour(velocity, iX, iY, iZ, 2);
    T dzuy = fdDataField::bulkZderivOrderFour(velocity, iX, iY, iZ, 1);

    return dyuz - dzuy;
}

template <typename T, int nDim>
inline T bulkVorticityYOrderFour(
    TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ)
{
    T dzux = fdDataField::bulkZderivOrderFour(velocity, iX, iY, iZ, 0);
    T dxuz = fdDataField::bulkXderivOrderFour(velocity, iX, iY, iZ, 2);

    return dzux - dxuz;
}

template <typename T, int nDim>
inline T bulkVorticityZOrderFour(
    TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ)
{
    T dxuy = fdDataField::bulkXderivOrderFour(velocity, iX, iY, iZ, 1);
    T dyux = fdDataField::bulkYderivOrderFour(velocity, iX, iY, iZ, 0);

    return dxuy - dyux;
}

template <typename T, int nDim>
inline T bulkVorticityXOrderSix(
    TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ)
{
    T dyuz = fdDataField::bulkYderivOrderSix(velocity, iX, iY, iZ, 2);
    T dzuy = fdDataField::bulkZderivOrderSix(velocity, iX, iY, iZ, 1);

    return dyuz - dzuy;
}

template <typename T, int nDim>
inline T bulkVorticityYOrderSix(
    TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ)
{
    T dzux = fdDataField::bulkZderivOrderSix(velocity, iX, iY, iZ, 0);
    T dxuz = fdDataField::bulkXderivOrderSix(velocity, iX, iY, iZ, 2);

    return dzux - dxuz;
}

template <typename T, int nDim>
inline T bulkVorticityZOrderSix(
    TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ)
{
    T dxuy = fdDataField::bulkXderivOrderSix(velocity, iX, iY, iZ, 1);
    T dyux = fdDataField::bulkYderivOrderSix(velocity, iX, iY, iZ, 0);

    return dxuy - dyux;
}

template <typename T, int nDim>
inline T bulkVorticityXOrderEight(
    TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ)
{
    T dyuz = fdDataField::bulkYderivOrderEight(velocity, iX, iY, iZ, 2);
    T dzuy = fdDataField::bulkZderivOrderEight(velocity, iX, iY, iZ, 1);

    return dyuz - dzuy;
}

template <typename T, int nDim>
inline T bulkVorticityYOrderEight(
    TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ)
{
    T dzux = fdDataField::bulkZderivOrderEight(velocity, iX, iY, iZ, 0);
    T dxuz = fdDataField::bulkXderivOrderEight(velocity, iX, iY, iZ, 2);

    return dzux - dxuz;
}

template <typename T, int nDim>
inline T bulkVorticityZOrderEight(
    TensorField3D<T, nDim> const &velocity, plint iX, plint iY, plint iZ)
{
    T dxuy = fdDataField::bulkXderivOrderEight(velocity, iX, iY, iZ, 1);
    T dyux = fdDataField::bulkYderivOrderEight(velocity, iX, iY, iZ, 0);

    return dxuy - dyux;
}

template <typename T, int nDim>
inline T planeVorticityX(
    TensorField3D<T, nDim> const &velocity, int direction, int orientation, plint iX, plint iY,
    plint iZ)
{
    T dyuz = fdDataField::planeYderiv(velocity, direction, orientation, iX, iY, iZ, 2);
    T dzuy = fdDataField::planeZderiv(velocity, direction, orientation, iX, iY, iZ, 1);

    return dyuz - dzuy;
}

template <typename T, int nDim>
inline T planeVorticityY(
    TensorField3D<T, nDim> const &velocity, int direction, int orientation, plint iX, plint iY,
    plint iZ)
{
    T dzux = fdDataField::planeZderiv(velocity, direction, orientation, iX, iY, iZ, 0);
    T dxuz = fdDataField::planeXderiv(velocity, direction, orientation, iX, iY, iZ, 2);

    return dzux - dxuz;
}

template <typename T, int nDim>
inline T planeVorticityZ(
    TensorField3D<T, nDim> const &velocity, int direction, int orientation, plint iX, plint iY,
    plint iZ)
{
    T dxuy = fdDataField::planeXderiv(velocity, direction, orientation, iX, iY, iZ, 1);
    T dyux = fdDataField::planeYderiv(velocity, direction, orientation, iX, iY, iZ, 0);

    return dxuy - dyux;
}

template <typename T, int nDim>
inline T edgeVorticityX(
    TensorField3D<T, nDim> const &velocity, int plane, int normal1, int normal2, plint iX, plint iY,
    plint iZ)
{
    T dyuz = fdDataField::edgeYderiv(velocity, plane, normal1, normal2, iX, iY, iZ, 2);
    T dzuy = fdDataField::edgeZderiv(velocity, plane, normal1, normal2, iX, iY, iZ, 1);

    return dyuz - dzuy;
}

template <typename T, int nDim>
inline T edgeVorticityY(
    TensorField3D<T, nDim> const &velocity, int plane, int normal1, int normal2, plint iX, plint iY,
    plint iZ)
{
    T dzux = fdDataField::edgeZderiv(velocity, plane, normal1, normal2, iX, iY, iZ, 0);
    T dxuz = fdDataField::edgeXderiv(velocity, plane, normal1, normal2, iX, iY, iZ, 2);

    return dzux - dxuz;
}

template <typename T, int nDim>
inline T edgeVorticityZ(
    TensorField3D<T, nDim> const &velocity, int plane, int normal1, int normal2, plint iX, plint iY,
    plint iZ)
{
    T dxuy = fdDataField::edgeXderiv(velocity, plane, normal1, normal2, iX, iY, iZ, 1);
    T dyux = fdDataField::edgeYderiv(velocity, plane, normal1, normal2, iX, iY, iZ, 0);

    return dxuy - dyux;
}

template <typename T, int nDim>
inline T cornerVorticityX(
    TensorField3D<T, nDim> const &velocity, int normalX, int normalY, int normalZ, plint iX,
    plint iY, plint iZ)
{
    T dyuz = fdDataField::cornerYderiv(velocity, normalX, normalY, normalZ, iX, iY, iZ, 2);
    T dzuy = fdDataField::cornerZderiv(velocity, normalX, normalY, normalZ, iX, iY, iZ, 1);

    return dyuz - dzuy;
}

template <typename T, int nDim>
inline T cornerVorticityY(
    TensorField3D<T, nDim> const &velocity, int normalX, int normalY, int normalZ, plint iX,
    plint iY, plint iZ)
{
    T dzux = fdDataField::cornerZderiv(velocity, normalX, normalY, normalZ, iX, iY, iZ, 0);
    T dxuz = fdDataField::cornerXderiv(velocity, normalX, normalY, normalZ, iX, iY, iZ, 2);

    return dzux - dxuz;
}

template <typename T, int nDim>
inline T cornerVorticityZ(
    TensorField3D<T, nDim> const &velocity, int normalX, int normalY, int normalZ, plint iX,
    plint iY, plint iZ)
{
    T dxuy = fdDataField::cornerXderiv(velocity, normalX, normalY, normalZ, iX, iY, iZ, 1);
    T dyux = fdDataField::cornerYderiv(velocity, normalX, normalY, normalZ, iX, iY, iZ, 0);

    return dxuy - dyux;
}
template <typename T>
inline T bulkXderiv(ScalarField2D<T> const &field, plint iX, plint iY)
{
    T dxu = fd::ctl_diff(field.get(iX + 1, iY), field.get(iX - 1, iY));
    return dxu;
}

template <typename T>
inline T bulkYderiv(ScalarField2D<T> const &field, plint iX, plint iY)
{
    T dyu = fd::ctl_diff(field.get(iX, iY + 1), field.get(iX, iY - 1));
    return dyu;
}

template <typename T>
inline T edgeXderiv(
    ScalarField2D<T> const &field, int direction, int orientation, plint iX, plint iY)
{
    if (direction == 0) {
        return -orientation
               * fd::o1_fwd_diff(field.get(iX, iY), field.get(iX - 1 * orientation, iY));
    } else {
        return bulkXderiv(field, iX, iY);
    }
}

template <typename T>
inline T edgeYderiv(
    ScalarField2D<T> const &field, int direction, int orientation, plint iX, plint iY)
{
    if (direction == 1) {
        return -orientation
               * fd::o1_fwd_diff(field.get(iX, iY), field.get(iX, iY - 1 * orientation));
    } else {
        return bulkYderiv(field, iX, iY);
    }
}

template <typename T>
inline T cornerXderiv(ScalarField2D<T> const &field, int normalX, int normalY, plint iX, plint iY)
{
    int orientation = normalX;
    return -orientation * fd::o1_fwd_diff(field.get(iX, iY), field.get(iX - 1 * orientation, iY));
}

template <typename T>
inline T cornerYderiv(ScalarField2D<T> const &field, int normalX, int normalY, plint iX, plint iY)
{
    int orientation = normalY;
    return -orientation * fd::o1_fwd_diff(field.get(iX, iY), field.get(iX, iY - 1 * orientation));
}

template <typename T, int nDim>
inline T bulkXderiv(TensorField2D<T, nDim> const &velocity, plint iX, plint iY, int iD)
{
    T dxu = fd::ctl_diff(velocity.get(iX + 1, iY)[iD], velocity.get(iX - 1, iY)[iD]);
    return dxu;
}

template <typename T, int nDim>
inline T bulkYderiv(TensorField2D<T, nDim> const &velocity, plint iX, plint iY, int iD)
{
    T dyu = fd::ctl_diff(velocity.get(iX, iY + 1)[iD], velocity.get(iX, iY - 1)[iD]);
    return dyu;
}

template <typename T, int nDim>
inline T edgeXderiv(
    TensorField2D<T, nDim> const &velocity, int direction, int orientation, plint iX, plint iY,
    int iD)
{
    if (direction == 0) {
        return -orientation
               * fd::o1_fwd_diff(
                   velocity.get(iX, iY)[iD], velocity.get(iX - 1 * orientation, iY)[iD]);
    } else {
        return bulkXderiv(velocity, iX, iY, iD);
    }
}

template <typename T, int nDim>
inline T edgeYderiv(
    TensorField2D<T, nDim> const &velocity, int direction, int orientation, plint iX, plint iY,
    int iD)
{
    if (direction == 1) {
        return -orientation
               * fd::o1_fwd_diff(
                   velocity.get(iX, iY)[iD], velocity.get(iX, iY - 1 * orientation)[iD]);
    } else {
        return bulkYderiv(velocity, iX, iY, iD);
    }
}

template <typename T, int nDim>
inline T cornerXderiv(
    TensorField2D<T, nDim> const &velocity, int normalX, int normalY, plint iX, plint iY, int iD)
{
    int orientation = normalX;
    return -orientation
           * fd::o1_fwd_diff(velocity.get(iX, iY)[iD], velocity.get(iX - 1 * orientation, iY)[iD]);
}

template <typename T, int nDim>
inline T cornerYderiv(
    TensorField2D<T, nDim> const &velocity, int normalX, int normalY, plint iX, plint iY, int iD)
{
    int orientation = normalY;
    return -orientation
           * fd::o1_fwd_diff(velocity.get(iX, iY)[iD], velocity.get(iX, iY - 1 * orientation)[iD]);
}

template <typename T, int nDim>
inline T bulkVorticity(TensorField2D<T, nDim> const &velocity, plint iX, plint iY)
{
    T dxuy = bulkXderiv(velocity, iX, iY, 1);
    T dyux = bulkYderiv(velocity, iX, iY, 0);
    return dxuy - dyux;
}

template <typename T, int nDim>
inline T edgeVorticity(
    TensorField2D<T, nDim> const &velocity, int direction, int orientation, plint iX, plint iY)
{
    T dxuy = edgeXderiv(velocity, direction, orientation, iX, iY, 1);
    T dyux = edgeYderiv(velocity, direction, orientation, iX, iY, 0);
    return dxuy - dyux;
}

template <typename T, int nDim>
inline T cornerVorticity(
    TensorField2D<T, nDim> const &velocity, int normalX, int normalY, plint iX, plint iY)
{
    T dxuy = cornerXderiv(velocity, normalX, normalY, iX, iY, 1);
    T dyux = cornerYderiv(velocity, normalX, normalY, iX, iY, 0);
    return dxuy - dyux;
}

}  // namespace fdDataField

/// Finite Difference operations on data-fields
namespace fdLattice  // fdLattice
{

template <typename T, template <typename U> class Descriptor>
inline Array<T, 3> firstOrderBulkVorticity(
    BlockLattice3D<T, Descriptor> const &lattice, plint iX, plint iY, plint iZ)
{
    Array<T, 3> u000;
    lattice.get(iX, iY, iZ).computeVelocity(u000);
    Array<T, 3> u001;
    lattice.get(iX, iY, iZ + 1).computeVelocity(u001);
    Array<T, 3> u010;
    lattice.get(iX, iY + 1, iZ).computeVelocity(u010);
    Array<T, 3> u100;
    lattice.get(iX + 1, iY, iZ).computeVelocity(u100);

    Array<T, 3> vorticity;
    T dyuz = u010[2] - u000[2];
    T dzuy = u001[1] - u000[1];
    vorticity[0] = dyuz - dzuy;

    T dzux = u001[0] - u000[0];
    T dxuz = u100[2] - u000[2];
    vorticity[1] = dzux - dxuz;

    T dxuy = u100[1] - u000[1];
    T dyux = u010[0] - u000[0];
    vorticity[2] = dxuy - dyux;

    return vorticity;
}

}  // namespace fdLattice
}  // namespace plb

#endif  // FD_STENCILS_1D_H
