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
 * Helper functions for domain initialization -- header file.
 */
#ifndef GRID_REFINEMENT_FUNCTIONAL_3D_HH
#define GRID_REFINEMENT_FUNCTIONAL_3D_HH

#include <cmath>
#include <limits>

#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "core/plbDebug.h"
#include "core/util.h"
#include "finiteDifference/fdStencils1D.h"
#include "finiteDifference/fdStencils2D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

template <typename T>
class DirectedInterpolation2D {
public:
    virtual ~DirectedInterpolation2D() {};

    virtual void interp(
        const NTensorField3D<T> &cTensor, plint iX, plint iY, plint iZ, NTensorField3D<T> &fTensor,
        plint fX, plint fY, plint fZ) = 0;
};

// Interpolation is made on a X-directed plane
template <typename T>
class Xinterpolation2D : public DirectedInterpolation2D<T> {
public:
    Xinterpolation2D<T>() : DirectedInterpolation2D<T>()
    {
        dOne = 1;
        dTwo = 1;
    }

    Xinterpolation2D(plint orientationOne_, plint orientationTwo_) :
        DirectedInterpolation2D<T>(), dOne(-orientationOne_), dTwo(-orientationTwo_)
    {
        dOne = (dOne == 0) ? +1 : dOne;
        dTwo = (dTwo == 0) ? +1 : dTwo;
        PLB_PRECONDITION(dOne == -1 || dOne == +1);
        PLB_PRECONDITION(dTwo == -1 || dTwo == +1);
    }

    virtual void interp(
        const NTensorField3D<T> &cTensor, plint iX, plint iY, plint iZ, NTensorField3D<T> &fTensor,
        plint fX, plint fY, plint fZ)
    {
        //         pcout << "IN INTERPOLATOR = " << iX << " " << iY << " " << iZ << " " << fX << " "
        //         << fY << " " << fZ << std::endl; pcout << "POSITIONS = (" << fY+dOne << " " <<
        //         fZ+dTwo << ")." << std::endl;

        plint nDim = cTensor.getNdim();
        for (plint iA = 0; iA < nDim; ++iA) {
            T a_m2m2 = cTensor.get(iX, iY - dOne, iZ - dTwo)[iA];
            T a_m2m1 = cTensor.get(iX, iY - dOne, iZ)[iA];
            T a_m2p1 = cTensor.get(iX, iY - dOne, iZ + dTwo)[iA];
            T a_m2p2 = cTensor.get(iX, iY - dOne, iZ + 2 * dTwo)[iA];

            T a_m1m2 = cTensor.get(iX, iY, iZ - dTwo)[iA];
            T a_m1m1 = cTensor.get(iX, iY, iZ)[iA];
            T a_m1p1 = cTensor.get(iX, iY, iZ + dTwo)[iA];
            T a_m1p2 = cTensor.get(iX, iY, iZ + 2 * dTwo)[iA];

            T a_p1m2 = cTensor.get(iX, iY + dOne, iZ - dTwo)[iA];
            T a_p1m1 = cTensor.get(iX, iY + dOne, iZ)[iA];
            T a_p1p1 = cTensor.get(iX, iY + dOne, iZ + dTwo)[iA];
            T a_p1p2 = cTensor.get(iX, iY + dOne, iZ + 2 * dTwo)[iA];

            T a_p2m2 = cTensor.get(iX, iY + 2 * dOne, iZ - dTwo)[iA];
            T a_p2m1 = cTensor.get(iX, iY + 2 * dOne, iZ)[iA];
            T a_p2p1 = cTensor.get(iX, iY + 2 * dOne, iZ + dTwo)[iA];
            T a_p2p2 = cTensor.get(iX, iY + 2 * dOne, iZ + 2 * dTwo)[iA];

            fTensor.get(fX, fY, fZ)[iA] = cTensor.get(iX, iY, iZ)[iA];

            fTensor.get(fX, fY, fZ + dTwo)[iA] =
                fd::centralFourPointsInterp(a_m1m2, a_m1m1, a_m1p1, a_m1p2);
            fTensor.get(fX, fY + dOne, fZ)[iA] =
                fd::centralFourPointsInterp(a_m2m1, a_m1m1, a_p1m1, a_p2m1);

            fTensor.get(fX, fY + dOne, fZ + dTwo)[iA] = fd::centralSixteenPointsInterp(
                a_m2m2, a_m2m1, a_m2p1, a_m2p2, a_m1m2, a_m1m1, a_m1p1, a_m1p2, a_p1m2, a_p1m1,
                a_p1p1, a_p1p2, a_p2m2, a_p2m1, a_p2p1, a_p2p2);
        }
    }

private:
    plint dOne, dTwo;
};

// Interpolation is made on a X-directed plane, Y direction is interpolated with 3 points
template <typename T>
class XinterpolationEdgeY2D : public DirectedInterpolation2D<T> {
public:
    XinterpolationEdgeY2D() : DirectedInterpolation2D<T>()
    {
        dOne = 1;
        dTwo = 1;
    }

    XinterpolationEdgeY2D(plint orientationOne_, plint orientationTwo_) :
        DirectedInterpolation2D<T>(), dOne(-orientationOne_), dTwo(-orientationTwo_)
    {
        if (dOne == 0)
            dOne = 1;
        if (dTwo == 0)
            dTwo = 1;
    }

    virtual void interp(
        const NTensorField3D<T> &cTensor, plint iX, plint iY, plint iZ, NTensorField3D<T> &fTensor,
        plint fX, plint fY, plint fZ)
    {
        plint nDim = cTensor.getNdim();
        for (plint iA = 0; iA < nDim; ++iA) {
            T a_m1m2 = cTensor.get(iX, iY, iZ - dTwo)[iA];
            T a_m1m1 = cTensor.get(iX, iY, iZ)[iA];
            T a_m1p1 = cTensor.get(iX, iY, iZ + dTwo)[iA];
            T a_m1p2 = cTensor.get(iX, iY, iZ + 2 * dTwo)[iA];

            T a_p1m2 = cTensor.get(iX, iY + dOne, iZ - dTwo)[iA];
            T a_p1m1 = cTensor.get(iX, iY + dOne, iZ)[iA];
            T a_p1p1 = cTensor.get(iX, iY + dOne, iZ + dTwo)[iA];
            T a_p1p2 = cTensor.get(iX, iY + dOne, iZ + 2 * dTwo)[iA];

            T a_p2m2 = cTensor.get(iX, iY + 2 * dOne, iZ - dTwo)[iA];
            T a_p2m1 = cTensor.get(iX, iY + 2 * dOne, iZ)[iA];
            T a_p2p1 = cTensor.get(iX, iY + 2 * dOne, iZ + dTwo)[iA];
            T a_p2p2 = cTensor.get(iX, iY + 2 * dOne, iZ + 2 * dTwo)[iA];

            fTensor.get(fX, fY, fZ)[iA] = cTensor.get(iX, iY, iZ)[iA];

            fTensor.get(fX, fY, fZ + dTwo)[iA] =
                fd::centralFourPointsInterp(a_m1m2, a_m1m1, a_m1p1, a_m1p2);
            fTensor.get(fX, fY + dOne, fZ)[iA] =
                fd::cornerThreePointsInterp(a_m1m1, a_p1m1, a_p2m1);

            fTensor.get(fX, fY + dOne, fZ + dTwo)[iA] = fd::asymTwelvePointsInterp(
                a_m1m2, a_p1m2, a_p2m2, a_m1m1, a_p1m1, a_p2m1, a_m1p1, a_p1p1, a_p2p1, a_m1p2,
                a_p1p2, a_p2p2);
        }
    }

private:
    plint dOne, dTwo;
};

// Interpolation is made on a X-directed plane, Z direction is interpolated with 3 points
template <typename T>
class XinterpolationEdgeZ2D : public DirectedInterpolation2D<T> {
public:
    XinterpolationEdgeZ2D() : DirectedInterpolation2D<T>()
    {
        dOne = 1;
        dTwo = 1;
    }

    XinterpolationEdgeZ2D(plint orientationOne_, plint orientationTwo_) :
        DirectedInterpolation2D<T>(), dOne(-orientationOne_), dTwo(-orientationTwo_)
    {
        if (dOne == 0)
            dOne = 1;
        if (dTwo == 0)
            dTwo = 1;
    }

    virtual void interp(
        const NTensorField3D<T> &cTensor, plint iX, plint iY, plint iZ, NTensorField3D<T> &fTensor,
        plint fX, plint fY, plint fZ)
    {
        plint nDim = cTensor.getNdim();
        for (plint iA = 0; iA < nDim; ++iA) {
            T a_m2m1 = cTensor.get(iX, iY - dOne, iZ)[iA];
            T a_m2p1 = cTensor.get(iX, iY - dOne, iZ + dTwo)[iA];
            T a_m2p2 = cTensor.get(iX, iY - dOne, iZ + 2 * dTwo)[iA];

            T a_m1m1 = cTensor.get(iX, iY, iZ)[iA];
            T a_m1p1 = cTensor.get(iX, iY, iZ + dTwo)[iA];
            T a_m1p2 = cTensor.get(iX, iY, iZ + 2 * dTwo)[iA];

            T a_p1m1 = cTensor.get(iX, iY + dOne, iZ)[iA];
            T a_p1p1 = cTensor.get(iX, iY + dOne, iZ + dTwo)[iA];
            T a_p1p2 = cTensor.get(iX, iY + dOne, iZ + 2 * dTwo)[iA];

            T a_p2m1 = cTensor.get(iX, iY + 2 * dOne, iZ)[iA];
            T a_p2p1 = cTensor.get(iX, iY + 2 * dOne, iZ + dTwo)[iA];
            T a_p2p2 = cTensor.get(iX, iY + 2 * dOne, iZ + 2 * dTwo)[iA];

            fTensor.get(fX, fY, fZ)[iA] = cTensor.get(iX, iY, iZ)[iA];

            fTensor.get(fX, fY, fZ + dTwo)[iA] =
                fd::cornerThreePointsInterp(a_m1m1, a_m1p1, a_m1p2);
            fTensor.get(fX, fY + dOne, fZ)[iA] =
                fd::centralFourPointsInterp(a_m2m1, a_m1m1, a_p1m1, a_p2m1);

            fTensor.get(fX, fY + dOne, fZ + dTwo)[iA] = fd::asymTwelvePointsInterp(
                a_m2m1, a_m2p1, a_m2p2, a_m1m1, a_m1p1, a_m1p2, a_p1m1, a_p1p1, a_p1p2, a_p2m1,
                a_p2p1, a_p2p2);
        }
    }

private:
    plint dOne, dTwo;
};

// Interpolation is made on a X-directed plane, Y and Z direction is interpolated with 3 points
template <typename T>
class XinterpolationCorner2D : public DirectedInterpolation2D<T> {
public:
    XinterpolationCorner2D() : DirectedInterpolation2D<T>()
    {
        dOne = 1;
        dTwo = 1;
    }

    XinterpolationCorner2D(plint orientationOne_, plint orientationTwo_) :
        DirectedInterpolation2D<T>(), dOne(-orientationOne_), dTwo(-orientationTwo_)
    {
        if (dOne == 0)
            dOne = 1;
        if (dTwo == 0)
            dTwo = 1;
    }

    virtual void interp(
        const NTensorField3D<T> &cTensor, plint iX, plint iY, plint iZ, NTensorField3D<T> &fTensor,
        plint fX, plint fY, plint fZ)
    {
        plint nDim = cTensor.getNdim();
        for (plint iA = 0; iA < nDim; ++iA) {
            T a_m1m1 = cTensor.get(iX, iY, iZ)[iA];
            T a_m1p1 = cTensor.get(iX, iY, iZ + dTwo)[iA];
            T a_m1p2 = cTensor.get(iX, iY, iZ + 2 * dTwo)[iA];

            T a_p1m1 = cTensor.get(iX, iY + dOne, iZ)[iA];
            T a_p1p1 = cTensor.get(iX, iY + dOne, iZ + dTwo)[iA];
            T a_p1p2 = cTensor.get(iX, iY + dOne, iZ + 2 * dTwo)[iA];

            T a_p2m1 = cTensor.get(iX, iY + 2 * dOne, iZ)[iA];
            T a_p2p1 = cTensor.get(iX, iY + 2 * dOne, iZ + dTwo)[iA];
            T a_p2p2 = cTensor.get(iX, iY + 2 * dOne, iZ + 2 * dTwo)[iA];

            fTensor.get(fX, fY, fZ)[iA] = cTensor.get(iX, iY, iZ)[iA];

            fTensor.get(fX, fY, fZ + dTwo)[iA] =
                fd::cornerThreePointsInterp(a_m1m1, a_m1p1, a_m1p2);
            fTensor.get(fX, fY + dOne, fZ)[iA] =
                fd::cornerThreePointsInterp(a_m1m1, a_p1m1, a_p2m1);

            fTensor.get(fX, fY + dOne, fZ + dTwo)[iA] = fd::asymNinePointsInterp(
                a_m1m1, a_m1p1, a_m1p2, a_p1m1, a_p1p1, a_p1p2, a_p2m1, a_p2p1, a_p2p2);
        }
    }

private:
    plint dOne, dTwo;
};

// Interpolation is made on a Y-directed plane
template <typename T>
class Yinterpolation2D : public DirectedInterpolation2D<T> {
public:
    Yinterpolation2D<T>() : DirectedInterpolation2D<T>()
    {
        dOne = 1;
        dTwo = 1;
    }

    Yinterpolation2D(plint orientationOne_, plint orientationTwo_) :
        DirectedInterpolation2D<T>(), dOne(-orientationOne_), dTwo(-orientationTwo_)
    {
        dOne = (dOne == 0) ? +1 : dOne;
        dTwo = (dTwo == 0) ? +1 : dTwo;

        PLB_PRECONDITION(dOne == -1 || dOne == +1);
        PLB_PRECONDITION(dTwo == -1 || dTwo == +1);
    }

    virtual void interp(
        const NTensorField3D<T> &cTensor, plint iX, plint iY, plint iZ, NTensorField3D<T> &fTensor,
        plint fX, plint fY, plint fZ)
    {
        plint nDim = cTensor.getNdim();
        for (plint iA = 0; iA < nDim; ++iA) {
            T a_m2m2 = cTensor.get(iX - dTwo, iY, iZ - dOne)[iA];
            T a_m2m1 = cTensor.get(iX - dTwo, iY, iZ)[iA];
            T a_m2p1 = cTensor.get(iX - dTwo, iY, iZ + dOne)[iA];
            T a_m2p2 = cTensor.get(iX - dTwo, iY, iZ + 2 * dOne)[iA];

            T a_m1m2 = cTensor.get(iX, iY, iZ - dOne)[iA];
            T a_m1m1 = cTensor.get(iX, iY, iZ)[iA];
            T a_m1p1 = cTensor.get(iX, iY, iZ + dOne)[iA];
            T a_m1p2 = cTensor.get(iX, iY, iZ + 2 * dOne)[iA];

            T a_p1m2 = cTensor.get(iX + dTwo, iY, iZ - dOne)[iA];
            T a_p1m1 = cTensor.get(iX + dTwo, iY, iZ)[iA];
            T a_p1p1 = cTensor.get(iX + dTwo, iY, iZ + dOne)[iA];
            T a_p1p2 = cTensor.get(iX + dTwo, iY, iZ + 2 * dOne)[iA];

            T a_p2m2 = cTensor.get(iX + 2 * dTwo, iY, iZ - dOne)[iA];
            T a_p2m1 = cTensor.get(iX + 2 * dTwo, iY, iZ)[iA];
            T a_p2p1 = cTensor.get(iX + 2 * dTwo, iY, iZ + dOne)[iA];
            T a_p2p2 = cTensor.get(iX + 2 * dTwo, iY, iZ + 2 * dOne)[iA];

            fTensor.get(fX, fY, fZ)[iA] = cTensor.get(iX, iY, iZ)[iA];

            fTensor.get(fX, fY, fZ + dOne)[iA] =
                fd::centralFourPointsInterp(a_m1m2, a_m1m1, a_m1p1, a_m1p2);
            fTensor.get(fX + dTwo, fY, fZ)[iA] =
                fd::centralFourPointsInterp(a_m2m1, a_m1m1, a_p1m1, a_p2m1);

            fTensor.get(fX + dTwo, fY, fZ + dOne)[iA] = fd::centralSixteenPointsInterp(
                a_m2m2, a_m2m1, a_m2p1, a_m2p2, a_m1m2, a_m1m1, a_m1p1, a_m1p2, a_p1m2, a_p1m1,
                a_p1p1, a_p1p2, a_p2m2, a_p2m1, a_p2p1, a_p2p2);
        }
    }

private:
    plint dOne, dTwo;
};

// Interpolation is made on a X-directed plane, Y direction is interpolated with 3 points
template <typename T>
class YinterpolationEdgeX2D : public DirectedInterpolation2D<T> {
public:
    YinterpolationEdgeX2D() : DirectedInterpolation2D<T>()
    {
        dOne = 1;
        dTwo = 1;
    }

    YinterpolationEdgeX2D(plint orientationOne_, plint orientationTwo_) :
        DirectedInterpolation2D<T>(), dOne(-orientationOne_), dTwo(-orientationTwo_)
    {
        if (dOne == 0)
            dOne = 1;
        if (dTwo == 0)
            dTwo = 1;
    }

    virtual void interp(
        const NTensorField3D<T> &cTensor, plint iX, plint iY, plint iZ, NTensorField3D<T> &fTensor,
        plint fX, plint fY, plint fZ)
    {
        plint nDim = cTensor.getNdim();
        for (plint iA = 0; iA < nDim; ++iA) {
            T a_m1m2 = cTensor.get(iX, iY, iZ - dOne)[iA];
            T a_m1m1 = cTensor.get(iX, iY, iZ)[iA];
            T a_m1p1 = cTensor.get(iX, iY, iZ + dOne)[iA];
            T a_m1p2 = cTensor.get(iX, iY, iZ + 2 * dOne)[iA];

            T a_p1m2 = cTensor.get(iX + dTwo, iY, iZ - dOne)[iA];
            T a_p1m1 = cTensor.get(iX + dTwo, iY, iZ)[iA];
            T a_p1p1 = cTensor.get(iX + dTwo, iY, iZ + dOne)[iA];
            T a_p1p2 = cTensor.get(iX + dTwo, iY, iZ + 2 * dOne)[iA];

            T a_p2m2 = cTensor.get(iX + 2 * dTwo, iY, iZ - dOne)[iA];
            T a_p2m1 = cTensor.get(iX + 2 * dTwo, iY, iZ)[iA];
            T a_p2p1 = cTensor.get(iX + 2 * dTwo, iY, iZ + dOne)[iA];
            T a_p2p2 = cTensor.get(iX + 2 * dTwo, iY, iZ + 2 * dOne)[iA];

            fTensor.get(fX, fY, fZ)[iA] = cTensor.get(iX, iY, iZ)[iA];

            fTensor.get(fX, fY, fZ + dOne)[iA] =
                fd::centralFourPointsInterp(a_m1m2, a_m1m1, a_m1p1, a_m1p2);
            fTensor.get(fX + dTwo, fY, fZ)[iA] =
                fd::cornerThreePointsInterp(a_m1m1, a_p1m1, a_p2m1);

            fTensor.get(fX + dTwo, fY, fZ + dOne)[iA] = fd::asymTwelvePointsInterp(
                a_m1m2, a_p1m2, a_p2m2, a_m1m1, a_p1m1, a_p2m1, a_m1p1, a_p1p1, a_p2p1, a_m1p2,
                a_p1p2, a_p2p2);
        }
    }

private:
    plint dOne, dTwo;
};

// Interpolation is made on a X-directed plane, Z direction is interpolated with 3 points
template <typename T>
class YinterpolationEdgeZ2D : public DirectedInterpolation2D<T> {
public:
    YinterpolationEdgeZ2D() : DirectedInterpolation2D<T>()
    {
        dOne = 1;
        dTwo = 1;
    }

    YinterpolationEdgeZ2D(plint orientationOne_, plint orientationTwo_) :
        DirectedInterpolation2D<T>(), dOne(-orientationOne_), dTwo(-orientationTwo_)
    {
        if (dOne == 0)
            dOne = 1;
        if (dTwo == 0)
            dTwo = 1;
    }

    virtual void interp(
        const NTensorField3D<T> &cTensor, plint iX, plint iY, plint iZ, NTensorField3D<T> &fTensor,
        plint fX, plint fY, plint fZ)
    {
        plint nDim = cTensor.getNdim();
        for (plint iA = 0; iA < nDim; ++iA) {
            T a_m2m1 = cTensor.get(iX - dTwo, iY, iZ)[iA];
            T a_m2p1 = cTensor.get(iX - dTwo, iY, iZ + dOne)[iA];
            T a_m2p2 = cTensor.get(iX - dTwo, iY, iZ + 2 * dOne)[iA];

            T a_m1m1 = cTensor.get(iX, iY, iZ)[iA];
            T a_m1p1 = cTensor.get(iX, iY, iZ + dOne)[iA];
            T a_m1p2 = cTensor.get(iX, iY, iZ + 2 * dOne)[iA];

            T a_p1m1 = cTensor.get(iX + dTwo, iY, iZ)[iA];
            T a_p1p1 = cTensor.get(iX + dTwo, iY, iZ + dOne)[iA];
            T a_p1p2 = cTensor.get(iX + dTwo, iY, iZ + 2 * dOne)[iA];

            T a_p2m1 = cTensor.get(iX + 2 * dTwo, iY, iZ)[iA];
            T a_p2p1 = cTensor.get(iX + 2 * dTwo, iY, iZ + dOne)[iA];
            T a_p2p2 = cTensor.get(iX + 2 * dTwo, iY, iZ + 2 * dOne)[iA];

            fTensor.get(fX, fY, fZ)[iA] = cTensor.get(iX, iY, iZ)[iA];

            fTensor.get(fX, fY, fZ + dOne)[iA] =
                fd::cornerThreePointsInterp(a_m1m1, a_m1p1, a_m1p2);
            fTensor.get(fX + dTwo, fY, fZ)[iA] =
                fd::centralFourPointsInterp(a_m2m1, a_m1m1, a_p1m1, a_p2m1);

            fTensor.get(fX + dTwo, fY, fZ + dOne)[iA] = fd::asymTwelvePointsInterp(
                a_m2m1, a_m2p1, a_m2p2, a_m1m1, a_m1p1, a_m1p2, a_p1m1, a_p1p1, a_p1p2, a_p2m1,
                a_p2p1, a_p2p2);
        }
    }

private:
    plint dOne, dTwo;
};

// Interpolation is made on a X-directed plane, Y and Z direction is interpolated with 3 points
template <typename T>
class YinterpolationCorner2D : public DirectedInterpolation2D<T> {
public:
    YinterpolationCorner2D() : DirectedInterpolation2D<T>()
    {
        dOne = 1;
        dTwo = 1;
    }

    YinterpolationCorner2D(plint orientationOne_, plint orientationTwo_) :
        DirectedInterpolation2D<T>(), dOne(-orientationOne_), dTwo(-orientationTwo_)
    {
        if (dOne == 0)
            dOne = 1;
        if (dTwo == 0)
            dTwo = 1;
    }

    virtual void interp(
        const NTensorField3D<T> &cTensor, plint iX, plint iY, plint iZ, NTensorField3D<T> &fTensor,
        plint fX, plint fY, plint fZ)
    {
        plint nDim = cTensor.getNdim();
        for (plint iA = 0; iA < nDim; ++iA) {
            T a_m1m1 = cTensor.get(iX, iY, iZ)[iA];
            T a_m1p1 = cTensor.get(iX, iY, iZ + dOne)[iA];
            T a_m1p2 = cTensor.get(iX, iY, iZ + 2 * dOne)[iA];

            T a_p1m1 = cTensor.get(iX + dTwo, iY, iZ)[iA];
            T a_p1p1 = cTensor.get(iX + dTwo, iY, iZ + dOne)[iA];
            T a_p1p2 = cTensor.get(iX + dTwo, iY, iZ + 2 * dOne)[iA];

            T a_p2m1 = cTensor.get(iX + 2 * dTwo, iY, iZ)[iA];
            T a_p2p1 = cTensor.get(iX + 2 * dTwo, iY, iZ + dOne)[iA];
            T a_p2p2 = cTensor.get(iX + 2 * dTwo, iY, iZ + 2 * dOne)[iA];

            fTensor.get(fX, fY, fZ)[iA] = cTensor.get(iX, iY, iZ)[iA];

            fTensor.get(fX, fY, fZ + dOne)[iA] =
                fd::cornerThreePointsInterp(a_m1m1, a_m1p1, a_m1p2);
            fTensor.get(fX + dTwo, fY, fZ)[iA] =
                fd::cornerThreePointsInterp(a_m1m1, a_p1m1, a_p2m1);

            fTensor.get(fX + dTwo, fY, fZ + dOne)[iA] = fd::asymNinePointsInterp(
                a_m1m1, a_m1p1, a_m1p2, a_p1m1, a_p1p1, a_p1p2, a_p2m1, a_p2p1, a_p2p2);
        }
    }

private:
    plint dOne, dTwo;
};

// Interpolation is made on a Z-directed plane
template <typename T>
class Zinterpolation2D : public DirectedInterpolation2D<T> {
public:
    Zinterpolation2D<T>() : DirectedInterpolation2D<T>()
    {
        dOne = 1;
        dTwo = 1;
    }

    Zinterpolation2D(plint orientationOne_, plint orientationTwo_) :
        DirectedInterpolation2D<T>(), dOne(-orientationOne_), dTwo(-orientationTwo_)
    {
        dOne = (dOne == 0) ? +1 : dOne;
        dTwo = (dTwo == 0) ? +1 : dTwo;
        PLB_PRECONDITION(dOne == -1 || dOne == +1);
        PLB_PRECONDITION(dTwo == -1 || dTwo == +1);
    }

    virtual void interp(
        const NTensorField3D<T> &cTensor, plint iX, plint iY, plint iZ, NTensorField3D<T> &fTensor,
        plint fX, plint fY, plint fZ)
    {
        plint nDim = cTensor.getNdim();
        for (plint iA = 0; iA < nDim; ++iA) {
            T a_m2m2 = cTensor.get(iX - dOne, iY - dTwo, iZ)[iA];
            T a_m2m1 = cTensor.get(iX - dOne, iY, iZ)[iA];
            T a_m2p1 = cTensor.get(iX - dOne, iY + dTwo, iZ)[iA];
            T a_m2p2 = cTensor.get(iX - dOne, iY + 2 * dTwo, iZ)[iA];

            T a_m1m2 = cTensor.get(iX, iY - dTwo, iZ)[iA];
            T a_m1m1 = cTensor.get(iX, iY, iZ)[iA];
            T a_m1p1 = cTensor.get(iX, iY + dTwo, iZ)[iA];
            T a_m1p2 = cTensor.get(iX, iY + 2 * dTwo, iZ)[iA];

            T a_p1m2 = cTensor.get(iX + dOne, iY - dTwo, iZ)[iA];
            T a_p1m1 = cTensor.get(iX + dOne, iY, iZ)[iA];
            T a_p1p1 = cTensor.get(iX + dOne, iY + dTwo, iZ)[iA];
            T a_p1p2 = cTensor.get(iX + dOne, iY + 2 * dTwo, iZ)[iA];

            T a_p2m2 = cTensor.get(iX + 2 * dOne, iY - dTwo, iZ)[iA];
            T a_p2m1 = cTensor.get(iX + 2 * dOne, iY, iZ)[iA];
            T a_p2p1 = cTensor.get(iX + 2 * dOne, iY + dTwo, iZ)[iA];
            T a_p2p2 = cTensor.get(iX + 2 * dOne, iY + 2 * dTwo, iZ)[iA];

            fTensor.get(fX, fY, fZ)[iA] = cTensor.get(iX, iY, iZ)[iA];

            fTensor.get(fX, fY + dTwo, fZ)[iA] =
                fd::centralFourPointsInterp(a_m1m2, a_m1m1, a_m1p1, a_m1p2);
            fTensor.get(fX + dOne, fY, fZ)[iA] =
                fd::centralFourPointsInterp(a_m2m1, a_m1m1, a_p1m1, a_p2m1);

            fTensor.get(fX + dOne, fY + dTwo, fZ)[iA] = fd::centralSixteenPointsInterp(
                a_m2m2, a_m2m1, a_m2p1, a_m2p2, a_m1m2, a_m1m1, a_m1p1, a_m1p2, a_p1m2, a_p1m1,
                a_p1p1, a_p1p2, a_p2m2, a_p2m1, a_p2p1, a_p2p2);
        }
    }

private:
    plint dOne, dTwo;
};

// Interpolation is made on a X-directed plane, Y direction is interpolated with 3 points
template <typename T>
class ZinterpolationEdgeY2D : public DirectedInterpolation2D<T> {
public:
    ZinterpolationEdgeY2D() : DirectedInterpolation2D<T>()
    {
        dOne = 1;
        dTwo = 1;
    }

    ZinterpolationEdgeY2D(plint orientationOne_, plint orientationTwo_) :
        DirectedInterpolation2D<T>(), dOne(-orientationOne_), dTwo(-orientationTwo_)
    {
        if (dOne == 0)
            dOne = 1;
        if (dTwo == 0)
            dTwo = 1;
    }

    virtual void interp(
        const NTensorField3D<T> &cTensor, plint iX, plint iY, plint iZ, NTensorField3D<T> &fTensor,
        plint fX, plint fY, plint fZ)
    {
        plint nDim = cTensor.getNdim();
        for (plint iA = 0; iA < nDim; ++iA) {
            T a_m1m2 = cTensor.get(iX - dOne, iY, iZ)[iA];
            T a_m1m1 = cTensor.get(iX, iY, iZ)[iA];
            T a_m1p1 = cTensor.get(iX + dOne, iY, iZ)[iA];
            T a_m1p2 = cTensor.get(iX + 2 * dOne, iY, iZ)[iA];

            T a_p1m2 = cTensor.get(iX - dOne, iY + dTwo, iZ)[iA];
            T a_p1m1 = cTensor.get(iX, iY + dTwo, iZ)[iA];
            T a_p1p1 = cTensor.get(iX + dOne, iY + dTwo, iZ)[iA];
            T a_p1p2 = cTensor.get(iX + 2 * dOne, iY + dTwo, iZ)[iA];

            T a_p2m2 = cTensor.get(iX - dOne, iY + 2 * dTwo, iZ)[iA];
            T a_p2m1 = cTensor.get(iX, iY + 2 * dTwo, iZ)[iA];
            T a_p2p1 = cTensor.get(iX + dOne, iY + 2 * dTwo, iZ)[iA];
            T a_p2p2 = cTensor.get(iX + 2 * dOne, iY + 2 * dTwo, iZ)[iA];

            fTensor.get(fX, fY, fZ)[iA] = cTensor.get(iX, iY, iZ)[iA];

            fTensor.get(fX + dOne, fY, fZ)[iA] =
                fd::centralFourPointsInterp(a_m1m2, a_m1m1, a_m1p1, a_m1p2);
            fTensor.get(fX, fY + dTwo, fZ)[iA] =
                fd::cornerThreePointsInterp(a_m1m1, a_p1m1, a_p2m1);

            fTensor.get(fX + dOne, fY + dTwo, fZ)[iA] = fd::asymTwelvePointsInterp(
                a_m1m2, a_p1m2, a_p2m2, a_m1m1, a_p1m1, a_p2m1, a_m1p1, a_p1p1, a_p2p1, a_m1p2,
                a_p1p2, a_p2p2);
        }
    }

private:
    plint dOne, dTwo;
};

// Interpolation is made on a X-directed plane, Z direction is interpolated with 3 points
template <typename T>
class ZinterpolationEdgeX2D : public DirectedInterpolation2D<T> {
public:
    ZinterpolationEdgeX2D() : DirectedInterpolation2D<T>()
    {
        dOne = 1;
        dTwo = 1;
    }

    ZinterpolationEdgeX2D(plint orientationOne_, plint orientationTwo_) :
        DirectedInterpolation2D<T>(), dOne(-orientationOne_), dTwo(-orientationTwo_)
    {
        if (dOne == 0)
            dOne = 1;
        if (dTwo == 0)
            dTwo = 1;
    }

    virtual void interp(
        const NTensorField3D<T> &cTensor, plint iX, plint iY, plint iZ, NTensorField3D<T> &fTensor,
        plint fX, plint fY, plint fZ)
    {
        plint nDim = cTensor.getNdim();
        for (plint iA = 0; iA < nDim; ++iA) {
            T a_m2m1 = cTensor.get(iX, iY - dTwo, iZ)[iA];
            T a_m2p1 = cTensor.get(iX + dOne, iY - dTwo, iZ)[iA];
            T a_m2p2 = cTensor.get(iX + 2 * dOne, iY - dTwo, iZ)[iA];

            T a_m1m1 = cTensor.get(iX, iY, iZ)[iA];
            T a_m1p1 = cTensor.get(iX + dOne, iY, iZ)[iA];
            T a_m1p2 = cTensor.get(iX + 2 * dOne, iY, iZ)[iA];

            T a_p1m1 = cTensor.get(iX, iY + dTwo, iZ)[iA];
            T a_p1p1 = cTensor.get(iX + dOne, iY + dTwo, iZ)[iA];
            T a_p1p2 = cTensor.get(iX + 2 * dOne, iY + dTwo, iZ)[iA];

            T a_p2m1 = cTensor.get(iX, iY + 2 * dTwo, iZ)[iA];
            T a_p2p1 = cTensor.get(iX + dOne, iY + 2 * dTwo, iZ)[iA];
            T a_p2p2 = cTensor.get(iX + 2 * dOne, iY + 2 * dTwo, iZ)[iA];

            fTensor.get(fX, fY, fZ)[iA] = cTensor.get(iX, iY, iZ)[iA];

            fTensor.get(fX + dOne, fY, fZ)[iA] =
                fd::cornerThreePointsInterp(a_m1m1, a_m1p1, a_m1p2);
            fTensor.get(fX, fY + dTwo, fZ)[iA] =
                fd::centralFourPointsInterp(a_m2m1, a_m1m1, a_p1m1, a_p2m1);

            fTensor.get(fX + dOne, fY + dTwo, fZ)[iA] = fd::asymTwelvePointsInterp(
                a_m2m1, a_m2p1, a_m2p2, a_m1m1, a_m1p1, a_m1p2, a_p1m1, a_p1p1, a_p1p2, a_p2m1,
                a_p2p1, a_p2p2);
        }
    }

private:
    plint dOne, dTwo;
};

// Interpolation is made on a X-directed plane, Y and Z direction is interpolated with 3 points
template <typename T>
class ZinterpolationCorner2D : public DirectedInterpolation2D<T> {
public:
    ZinterpolationCorner2D() : DirectedInterpolation2D<T>()
    {
        dOne = 1;
        dTwo = 1;
    }

    ZinterpolationCorner2D(plint orientationOne_, plint orientationTwo_) :
        DirectedInterpolation2D<T>(), dOne(-orientationOne_), dTwo(-orientationTwo_)
    {
        if (dOne == 0)
            dOne = 1;
        if (dTwo == 0)
            dTwo = 1;
    }

    virtual void interp(
        const NTensorField3D<T> &cTensor, plint iX, plint iY, plint iZ, NTensorField3D<T> &fTensor,
        plint fX, plint fY, plint fZ)
    {
        plint nDim = cTensor.getNdim();
        for (plint iA = 0; iA < nDim; ++iA) {
            T a_m1m1 = cTensor.get(iX, iY, iZ)[iA];
            T a_m1p1 = cTensor.get(iX + dOne, iY, iZ)[iA];
            T a_m1p2 = cTensor.get(iX + 2 * dOne, iY, iZ)[iA];

            T a_p1m1 = cTensor.get(iX, iY + dTwo, iZ)[iA];
            T a_p1p1 = cTensor.get(iX + dOne, iY + dTwo, iZ)[iA];
            T a_p1p2 = cTensor.get(iX + 2 * dOne, iY + dTwo, iZ)[iA];

            T a_p2m1 = cTensor.get(iX, iY + 2 * dTwo, iZ)[iA];
            T a_p2p1 = cTensor.get(iX + dOne, iY + 2 * dTwo, iZ)[iA];
            T a_p2p2 = cTensor.get(iX + 2 * dOne, iY + 2 * dTwo, iZ)[iA];

            fTensor.get(fX, fY, fZ)[iA] = cTensor.get(iX, iY, iZ)[iA];

            fTensor.get(fX + dOne, fY, fZ)[iA] =
                fd::cornerThreePointsInterp(a_m1m1, a_m1p1, a_m1p2);
            fTensor.get(fX, fY + dTwo, fZ)[iA] =
                fd::cornerThreePointsInterp(a_m1m1, a_p1m1, a_p2m1);

            fTensor.get(fX + dOne, fY + dTwo, fZ)[iA] = fd::asymNinePointsInterp(
                a_m1m1, a_m1p1, a_m1p2, a_p1m1, a_p1p1, a_p1p2, a_p2m1, a_p2p1, a_p2p2);
        }
    }

private:
    plint dOne, dTwo;
};

// =================== TemporalInterpolationFunctional3D =================== //
template <typename T>
void TemporalInterpolationFunctional3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    global::timer("gr_only_time").start();

    NTensorField3D<T> &tensor_t0 = *dynamic_cast<NTensorField3D<T> *>(fields[0]);
    NTensorField3D<T> &tensor_t1 = *dynamic_cast<NTensorField3D<T> *>(fields[1]);
    NTensorField3D<T> &tensor = *dynamic_cast<NTensorField3D<T> *>(fields[2]);

    PLB_PRECONDITION(tensor_t0.getNdim() == tensor_t1.getNdim());
    PLB_PRECONDITION(tensor_t0.getNdim() == tensor.getNdim());
    PLB_PRECONDITION(tensor_t1.getNdim() == tensor.getNdim());

    plint nDim = tensor_t0.getNdim();
    Dot3D offset_t1 = computeRelativeDisplacement(tensor_t0, tensor_t1);
    Dot3D offset = computeRelativeDisplacement(tensor_t0, tensor);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = iX + offset.x;
        plint tX = iX + offset_t1.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY = iY + offset.y;
            plint tY = iY + offset_t1.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint oZ = iZ + offset.z;
                plint tZ = iZ + offset_t1.z;

                // Linear interpolation in time ((t0+t1)/2)
                for (plint iA = 0; iA < nDim; ++iA) {
                    tensor.get(oX, oY, oZ)[iA] =
                        (tensor_t0.get(iX, iY, iZ)[iA] + tensor_t1.get(tX, tY, tZ)[iA]) * (T)0.5;
                }
            }
        }
    }
    global::timer("gr_only_time").stop();
}

template <typename T>
TemporalInterpolationFunctional3D<T> *TemporalInterpolationFunctional3D<T>::clone() const
{
    return new TemporalInterpolationFunctional3D<T>(*this);
}

template <typename T>
void TemporalInterpolationFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT TemporalInterpolationFunctional3D<T>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

// ============== CopyAndSpatialInterpolationPlaneFunctional3D ============== //

template <typename T>
CopyAndSpatialInterpolationPlaneFunctional3D<T>::CopyAndSpatialInterpolationPlaneFunctional3D(
    plint direction_) :
    direction(direction_)
{
    PLB_PRECONDITION(direction == 0 || direction == 1 || direction == 2);
}

template <typename T>
void CopyAndSpatialInterpolationPlaneFunctional3D<T>::process(
    Box3D domain, NTensorField3D<T> &cTensor, NTensorField3D<T> &fTensor)
{
    global::timer("gr_only_interp").start();

    PLB_PRECONDITION(domain.x0 == domain.x1 || domain.y0 == domain.y1 || domain.z0 == domain.z1);
    PLB_PRECONDITION(cTensor.getNdim() == fTensor.getNdim());

    DirectedInterpolation2D<T> *interpolator = 0;
    if (direction == 0) {
        interpolator = new Xinterpolation2D<T>();
    } else if (direction == 1) {
        interpolator = new Yinterpolation2D<T>();
    } else if (direction == 2) {
        interpolator = new Zinterpolation2D<T>();
    } else {
        PLB_ASSERT(false);
    }

    Dot3D coarsePos = cTensor.getLocation();
    Dot3D finePos = fTensor.getLocation();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint fX = 2 * (iX + coarsePos.x) - finePos.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint fY = 2 * (iY + coarsePos.y) - finePos.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint fZ = 2 * (iZ + coarsePos.z) - finePos.z;

                interpolator->interp(cTensor, iX, iY, iZ, fTensor, fX, fY, fZ);
            }
        }
    }

    delete interpolator;
    global::timer("gr_only_interp").start();
}

template <typename T>
CopyAndSpatialInterpolationPlaneFunctional3D<T>
    *CopyAndSpatialInterpolationPlaneFunctional3D<T>::clone() const
{
    return new CopyAndSpatialInterpolationPlaneFunctional3D<T>(*this);
}

template <typename T>
void CopyAndSpatialInterpolationPlaneFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

// ============== CopyAndSpatialInterpolationEdgeFunctional3D ============== //

template <typename T>
CopyAndSpatialInterpolationEdgeFunctional3D<T>::CopyAndSpatialInterpolationEdgeFunctional3D(
    plint direction_, plint orientOne_, plint orientTwo_) :
    direction(direction_), orientOne(orientOne_), orientTwo(orientTwo_)
{
    PLB_PRECONDITION(direction == 0 || direction == 1 || direction == 2);
}

template <typename T>
void CopyAndSpatialInterpolationEdgeFunctional3D<T>::process(
    Box3D domain, NTensorField3D<T> &cTensor, NTensorField3D<T> &fTensor)
{
    global::timer("gr_only_interp").start();

    PLB_PRECONDITION(
        (domain.x0 == domain.x1 && domain.y0 == domain.y1)
        || (domain.x0 == domain.x1 && domain.z0 == domain.z1)
        || (domain.y0 == domain.y1 && domain.z0 == domain.z1));

    PLB_PRECONDITION(cTensor.getNdim() == fTensor.getNdim());

    Dot3D coarsePos = cTensor.getLocation();
    Dot3D finePos = fTensor.getLocation();

    DirectedInterpolation2D<T> *interpolator = 0;
    if (direction == 0) {
        interpolator = new Xinterpolation2D<T>(orientOne, orientTwo);
    } else if (direction == 1) {
        interpolator = new Yinterpolation2D<T>(orientOne, orientTwo);
    } else if (direction == 2) {
        interpolator = new Zinterpolation2D<T>(orientOne, orientTwo);
    } else {
        PLB_ASSERT(false);
    }

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint fX = 2 * (iX + coarsePos.x) - finePos.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint fY = 2 * (iY + coarsePos.y) - finePos.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint fZ = 2 * (iZ + coarsePos.z) - finePos.z;

                interpolator->interp(cTensor, iX, iY, iZ, fTensor, fX, fY, fZ);
            }
        }
    }

    delete interpolator;
    global::timer("gr_only_interp").stop();
}

template <typename T>
CopyAndSpatialInterpolationEdgeFunctional3D<T>
    *CopyAndSpatialInterpolationEdgeFunctional3D<T>::clone() const
{
    return new CopyAndSpatialInterpolationEdgeFunctional3D<T>(*this);
}

template <typename T>
void CopyAndSpatialInterpolationEdgeFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

// ============== CopyAndSpatialInterpolationCornerFunctional3D ============== //

template <typename T>
CopyAndSpatialInterpolationCornerFunctional3D<T>::CopyAndSpatialInterpolationCornerFunctional3D(
    plint direction_, plint orientOne_, plint orientTwo_) :
    direction(direction_), orientOne(orientOne_), orientTwo(orientTwo_)
{
    PLB_PRECONDITION(direction == 0 || direction == 1 || direction == 2);
    PLB_PRECONDITION(
        (orientOne == +1 || (orientOne == -1)) && ((orientTwo == 1 || orientTwo == -1)));
}

template <typename T>
void CopyAndSpatialInterpolationCornerFunctional3D<T>::process(
    Box3D domain, NTensorField3D<T> &cTensor, NTensorField3D<T> &fTensor)
{
    PLB_PRECONDITION(domain.x0 == domain.x1 && domain.y0 == domain.y1 && domain.z0 == domain.z1);
    PLB_PRECONDITION(cTensor.getNdim() == fTensor.getNdim());

    global::timer("gr_only_interp").start();

    Dot3D coarsePos = cTensor.getLocation();
    Dot3D finePos = fTensor.getLocation();

    DirectedInterpolation2D<T> *interpolator = 0;
    if (direction == 0) {
        interpolator = new Xinterpolation2D<T>(orientOne, orientTwo);
    } else if (direction == 1) {
        interpolator = new Yinterpolation2D<T>(orientOne, orientTwo);
    } else if (direction == 2) {
        interpolator = new Zinterpolation2D<T>(orientOne, orientTwo);
    } else {
        PLB_ASSERT(false);
    }

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint fX = 2 * (iX + coarsePos.x) - finePos.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint fY = 2 * (iY + coarsePos.y) - finePos.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint fZ = 2 * (iZ + coarsePos.z) - finePos.z;

                interpolator->interp(cTensor, iX, iY, iZ, fTensor, fX, fY, fZ);
            }
        }
    }

    delete interpolator;

    global::timer("gr_only_interp").stop();
}

template <typename T>
CopyAndSpatialInterpolationCornerFunctional3D<T>
    *CopyAndSpatialInterpolationCornerFunctional3D<T>::clone() const
{
    return new CopyAndSpatialInterpolationCornerFunctional3D<T>(*this);
}

template <typename T>
void CopyAndSpatialInterpolationCornerFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

// ============== CopyAndSpatialInterpolationBoundaryEdgeFunctional3D ============== //

template <typename T>
CopyAndSpatialInterpolationBoundaryEdgeFunctional3D<T>::
    CopyAndSpatialInterpolationBoundaryEdgeFunctional3D(
        plint direction_, plint orientOne_, plint orientTwo_) :
    direction(direction_), orientOne(orientOne_), orientTwo(orientTwo_)
{
    PLB_ASSERT(
        (direction == 0 || direction == 1 || direction == 2) && "Direction must be 1,2, or 3.");
    PLB_ASSERT(
        ((orientOne == 0 && orientTwo != 0) || (orientOne != 0 && orientTwo == 0))
        && "orientOne == 0 and orientTwo != 0 || orientOne != 0 and orientTwo == 0");
}

template <typename T>
void CopyAndSpatialInterpolationBoundaryEdgeFunctional3D<T>::process(
    Box3D domain, NTensorField3D<T> &cTensor, NTensorField3D<T> &fTensor)
{
    global::timer("gr_only_interp").start();

    PLB_PRECONDITION(domain.x0 == domain.x1 || domain.y0 == domain.y1 || domain.z0 == domain.z1);
    PLB_PRECONDITION(cTensor.getNdim() == fTensor.getNdim());

    Dot3D coarsePos = cTensor.getLocation();
    Dot3D finePos = fTensor.getLocation();

    DirectedInterpolation2D<T> *interpolator = 0;
    if (direction == 0) {
        if (orientOne == 0) {
            interpolator = new XinterpolationEdgeZ2D<T>(orientOne, orientTwo);
        } else {
            interpolator = new XinterpolationEdgeY2D<T>(orientOne, orientTwo);
        }
    } else if (direction == 1) {
        if (orientOne == 0) {
            interpolator = new YinterpolationEdgeX2D<T>(orientOne, orientTwo);
        } else {
            interpolator = new YinterpolationEdgeZ2D<T>(orientOne, orientTwo);
        }
    } else if (direction == 2) {
        if (orientOne == 0) {
            interpolator = new ZinterpolationEdgeY2D<T>(orientOne, orientTwo);
        } else {
            interpolator = new ZinterpolationEdgeX2D<T>(orientOne, orientTwo);
        }
    } else {
        PLB_ASSERT(false);
    }

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint fX = 2 * (iX + coarsePos.x) - finePos.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint fY = 2 * (iY + coarsePos.y) - finePos.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint fZ = 2 * (iZ + coarsePos.z) - finePos.z;

                interpolator->interp(cTensor, iX, iY, iZ, fTensor, fX, fY, fZ);
            }
        }
    }
    delete interpolator;

    global::timer("gr_only_interp").stop();
}

template <typename T>
CopyAndSpatialInterpolationBoundaryEdgeFunctional3D<T>
    *CopyAndSpatialInterpolationBoundaryEdgeFunctional3D<T>::clone() const
{
    return new CopyAndSpatialInterpolationBoundaryEdgeFunctional3D<T>(*this);
}

template <typename T>
void CopyAndSpatialInterpolationBoundaryEdgeFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

// ============== CopyAndSpatialInterpolationBoundaryCornerFunctional3D ============== //

template <typename T>
CopyAndSpatialInterpolationBoundaryCornerFunctional3D<T>::
    CopyAndSpatialInterpolationBoundaryCornerFunctional3D(
        plint direction_, plint orientOne_, plint orientTwo_) :
    direction(direction_), orientOne(orientOne_), orientTwo(orientTwo_)
{ }

template <typename T>
void CopyAndSpatialInterpolationBoundaryCornerFunctional3D<T>::process(
    Box3D domain, NTensorField3D<T> &cTensor, NTensorField3D<T> &fTensor)
{
    global::timer("gr_only_interp").start();
    PLB_PRECONDITION(domain.x0 == domain.x1 && domain.y0 == domain.y1 && domain.z0 == domain.z1);
    PLB_PRECONDITION(cTensor.getNdim() == fTensor.getNdim());

    Dot3D coarsePos = cTensor.getLocation();
    Dot3D finePos = fTensor.getLocation();

    DirectedInterpolation2D<T> *interpolator = 0;
    if (direction == 0) {
        interpolator = new XinterpolationCorner2D<T>(orientOne, orientTwo);
    } else if (direction == 1) {
        interpolator = new YinterpolationCorner2D<T>(orientOne, orientTwo);
    } else if (direction == 2) {
        interpolator = new ZinterpolationCorner2D<T>(orientOne, orientTwo);
    } else {
        PLB_ASSERT(false);
    }

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint fX = 2 * (iX + coarsePos.x) - finePos.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint fY = 2 * (iY + coarsePos.y) - finePos.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint fZ = 2 * (iZ + coarsePos.z) - finePos.z;

                interpolator->interp(cTensor, iX, iY, iZ, fTensor, fX, fY, fZ);
            }
        }
    }
    delete interpolator;
    global::timer("gr_only_interp").stop();
}

template <typename T>
CopyAndSpatialInterpolationBoundaryCornerFunctional3D<T>
    *CopyAndSpatialInterpolationBoundaryCornerFunctional3D<T>::clone() const
{
    return new CopyAndSpatialInterpolationBoundaryCornerFunctional3D<T>(*this);
}

template <typename T>
void CopyAndSpatialInterpolationBoundaryCornerFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

// ============== CopyAndFilterNonEquilibriumFunctional3D ============== //

template <typename T, template <typename U> class Descriptor>
CopyAndFilterNonEquilibriumFunctional3D<T, Descriptor>::CopyAndFilterNonEquilibriumFunctional3D(
    bool filterAll)
{
    minIndexFilter = filterAll ? 0 : 1 + Descriptor<T>::d;
}

template <typename T, template <typename U> class Descriptor>
void CopyAndFilterNonEquilibriumFunctional3D<T, Descriptor>::process(
    Box3D domain, NTensorField3D<T> &cTensor, NTensorField3D<T> &fTensor)
{
    global::timer("gr_only_filter").start();
    using namespace descriptors;
    PLB_PRECONDITION(cTensor.getNdim() == fTensor.getNdim());

    plint nDim = cTensor.getNdim();

    Dot3D coarsePos = cTensor.getLocation();
    Dot3D finePos = fTensor.getLocation();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint fX = 2 * (iX + coarsePos.x) - finePos.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint fY = 2 * (iY + coarsePos.y) - finePos.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint fZ = 2 * (iZ + coarsePos.z) - finePos.z;

                for (plint iA = 0; iA < minIndexFilter; ++iA) {
                    cTensor.get(iX, iY, iZ)[iA] =
                        fTensor.get(fX, fY, fZ)[iA];  // rho and u may not be filtered (only copied)
                }

                for (plint iA = minIndexFilter;
                     iA < nDim - Descriptor<T>::ExternalField::numScalars; ++iA) {
                    cTensor.get(iX, iY, iZ)[iA] = fTensor.get(fX, fY, fZ)[iA];
                    for (plint iPop = 1; iPop < D3Q27Descriptor<T>::q; ++iPop) {
                        plint nextX = fX + D3Q27Descriptor<T>::c[iPop][0];
                        plint nextY = fY + D3Q27Descriptor<T>::c[iPop][1];
                        plint nextZ = fZ + D3Q27Descriptor<T>::c[iPop][2];

                        cTensor.get(iX, iY, iZ)[iA] += fTensor.get(nextX, nextY, nextZ)[iA];
                    }
                    cTensor.get(iX, iY, iZ)[iA] /= (T)D3Q27Descriptor<T>::q;
                }

                for (plint iA = nDim - Descriptor<T>::ExternalField::numScalars; iA < nDim; ++iA) {
                    cTensor.get(iX, iY, iZ)[iA] =
                        fTensor.get(fX, fY, fZ)[iA];  // rho and u may not be filtered (only copied)
                }
            }
        }
    }
    global::timer("gr_only_filter").stop();
}

template <typename T, template <typename U> class Descriptor>
CopyAndFilterNonEquilibriumFunctional3D<T, Descriptor>
    *CopyAndFilterNonEquilibriumFunctional3D<T, Descriptor>::clone() const
{
    return new CopyAndFilterNonEquilibriumFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void CopyAndFilterNonEquilibriumFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

// ============== CopyAndSelectiveFilterNonEquilibriumFunctional3D ============== //

template <typename T, template <typename U> class Descriptor>
CopyAndSelectiveFilterNonEquilibriumFunctional3D<
    T, Descriptor>::CopyAndSelectiveFilterNonEquilibriumFunctional3D(T sigma_, bool filterAll) :
    sigma(sigma_)
{
    minIndexFilter = filterAll ? 0 : 1 + Descriptor<T>::d;
    c = Array<int, 3>(-1, 0, 1);
    d = Array<T, 3>(-0.25, 0.5, 0.25);
}

template <typename T, template <typename U> class Descriptor>
void CopyAndSelectiveFilterNonEquilibriumFunctional3D<T, Descriptor>::process(
    Box3D domain, NTensorField3D<T> &cTensor, NTensorField3D<T> &fTensor)
{
    using namespace descriptors;
    PLB_PRECONDITION(cTensor.getNdim() == fTensor.getNdim());

    plint nDim = cTensor.getNdim();

    Dot3D coarsePos = cTensor.getLocation();
    Dot3D finePos = fTensor.getLocation();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint fX = 2 * (iX + coarsePos.x) - finePos.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint fY = 2 * (iY + coarsePos.y) - finePos.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint fZ = 2 * (iZ + coarsePos.z) - finePos.z;

                for (plint iA = 0; iA < minIndexFilter; ++iA) {
                    cTensor.get(iX, iY, iZ)[iA] =
                        fTensor.get(fX, fY, fZ)[iA];  // rho and u may not be filtered (only copied)
                }

                for (plint iA = minIndexFilter;
                     iA < nDim - Descriptor<T>::ExternalField::numScalars; ++iA) {
                    T fl = T();
                    for (plint iPop = 0; iPop < 3; ++iPop) {
                        fl += fTensor.get(fX + c[iPop], fY, fZ)[iA] * d[iPop];
                        fl += fTensor.get(fX, fY + c[iPop], fZ)[iA] * d[iPop];
                        fl += fTensor.get(fX, fY, fZ + c[iPop])[iA] * d[iPop];
                    }
                    cTensor.get(iX, iY, iZ)[iA] = fTensor.get(fX, fY, fZ)[iA] - sigma * fl;
                }

                for (plint iA = nDim - Descriptor<T>::ExternalField::numScalars; iA < nDim; ++iA) {
                    cTensor.get(iX, iY, iZ)[iA] =
                        fTensor.get(fX, fY, fZ)[iA];  // rho and u may not be filtered (only copied)
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
CopyAndSelectiveFilterNonEquilibriumFunctional3D<T, Descriptor>
    *CopyAndSelectiveFilterNonEquilibriumFunctional3D<T, Descriptor>::clone() const
{
    return new CopyAndSelectiveFilterNonEquilibriumFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void CopyAndSelectiveFilterNonEquilibriumFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

// ============== CopyFunctional3D ============== //

template <typename T, template <typename U> class Descriptor>
void CopyFunctional3D<T, Descriptor>::process(
    Box3D domain, NTensorField3D<T> &cTensor, NTensorField3D<T> &fTensor)
{
    global::timer("gr_only_copy").start();

    using namespace descriptors;
    PLB_PRECONDITION(cTensor.getNdim() == fTensor.getNdim());

    plint nDim = cTensor.getNdim();

    Dot3D coarsePos = cTensor.getLocation();
    Dot3D finePos = fTensor.getLocation();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint fX = 2 * (iX + coarsePos.x) - finePos.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint fY = 2 * (iY + coarsePos.y) - finePos.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint fZ = 2 * (iZ + coarsePos.z) - finePos.z;

                for (plint iA = 0; iA < nDim; ++iA) {
                    cTensor.get(iX, iY, iZ)[iA] = fTensor.get(fX, fY, fZ)[iA];
                }
            }
        }
    }
    global::timer("gr_only_copy").stop();
}

template <typename T, template <typename U> class Descriptor>
CopyFunctional3D<T, Descriptor> *CopyFunctional3D<T, Descriptor>::clone() const
{
    return new CopyFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void CopyFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

// ============== DecomposeAndRescaleFunctional3D ============== //

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
DecomposeAndRescaleFunctional3D<T, Descriptor, Engine>::DecomposeAndRescaleFunctional3D(
    T xDt_, T xDxInv_, plint order_) :
    xDt(xDt_), xDxInv(xDxInv_), order(order_)
{ }

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void DecomposeAndRescaleFunctional3D<T, Descriptor, Engine>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, NTensorField3D<T> &tensor)
{
    global::timer("gr_only_decomp").start();
    Engine<T, Descriptor> engine;

    std::vector<T> decompAndRescaled;
    Dot3D offset = computeRelativeDisplacement(lattice, tensor);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint oZ = iZ + offset.z;

                Cell<T, Descriptor> const &cell = lattice.get(iX, iY, iZ);

                // these two last line may need to be changed. Mayb by a global change
                // in the dynamics class in order to return all the relevant relaxation times
                // etc...

                engine.decomposeAndRescale(cell, xDt, order, decompAndRescaled);

                // if ((plint)decompAndRescaled.size() != (plint)tensor.getNdim()) {
                //     pcout << "decomp = " << decompAndRescaled.size() << ", tensir = " <<
                //     tensor.getNdim() << std::endl;
                // }
                PLB_ASSERT(decompAndRescaled.size() == (pluint)tensor.getNdim());

                // std::copy(decompAndRescaled.begin(), decompAndRescaled.end(),
                // tensor.get(oX,oY,oZ));

                for (pluint iA = 0; iA < decompAndRescaled.size(); ++iA) {
                    tensor.get(oX, oY, oZ)[iA] = decompAndRescaled[iA];
                }
            }
        }
    }
    global::timer("gr_only_decomp").stop();
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
DecomposeAndRescaleFunctional3D<T, Descriptor, Engine>
    *DecomposeAndRescaleFunctional3D<T, Descriptor, Engine>::clone() const
{
    return new DecomposeAndRescaleFunctional3D<T, Descriptor, Engine>(*this);
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void DecomposeAndRescaleFunctional3D<T, Descriptor, Engine>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
BlockDomain::DomainT DecomposeAndRescaleFunctional3D<T, Descriptor, Engine>::appliesTo() const
{
    return BlockDomain::bulk;
}

// =================== RecomposeFunctional3D ==================== //

template <typename T, template <typename U> class Descriptor>
RecomposeFunctional3D<T, Descriptor>::RecomposeFunctional3D(plint order_) : order(order_)
{ }

template <typename T, template <typename U> class Descriptor>
void RecomposeFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, NTensorField3D<T> &tensor)
{
    global::timer("gr_only_recomp").start();
    PLB_ASSERT(domain.x0 == domain.x1 || domain.y0 == domain.y1 || domain.z0 == domain.z1)
    plint nDim = tensor.getNdim();

    std::vector<T> decomposed(nDim);
    Dot3D offset = computeRelativeDisplacement(lattice, tensor);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint oZ = iZ + offset.z;

                Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);

                for (plint iA = 0; iA < nDim; ++iA) {
                    decomposed[iA] = tensor.get(oX, oY, oZ)[iA];
                }

                cell.getDynamics().recompose(cell, decomposed, order);
            }
        }
    }
    global::timer("gr_only_recomp").stop();
}

template <typename T, template <typename U> class Descriptor>
RecomposeFunctional3D<T, Descriptor> *RecomposeFunctional3D<T, Descriptor>::clone() const
{
    return new RecomposeFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void RecomposeFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT RecomposeFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

}  // namespace plb

#endif  // GRID_REFINEMENT_FUNCTIONAL_3D_HH
