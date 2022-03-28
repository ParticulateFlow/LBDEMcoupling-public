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
 * Dynamics and data processors used to implement 2D grid refinement -- generic implementation.
 */
#ifndef FINE_GRID_PROCESSORS_2D_HH
#define FINE_GRID_PROCESSORS_2D_HH

#include <iomanip>

#include "core/geometry2D.h"
#include "finiteDifference/fdStencils1D.h"
#include "multiBlock/multiDataProcessorWrapper2D.h"
#include "multiGrid/fineGridProcessors2D.h"

namespace plb {

/* *************** Class CopyCoarseToFineLinearInterp2D ****************** */

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineLinearInterp2D<T, Descriptor1, Descriptor2>::CopyCoarseToFineLinearInterp2D(
    RescaleEngine<T, Descriptor1> *rescaleEngine_, plint direction_, plint orientation_) :
    rescaleEngine(rescaleEngine_), direction(direction_), orientation(orientation_)
{ }

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineLinearInterp2D<T, Descriptor1, Descriptor2>::~CopyCoarseToFineLinearInterp2D()
{
    delete rescaleEngine;
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineLinearInterp2D<T, Descriptor1, Descriptor2>::CopyCoarseToFineLinearInterp2D(
    CopyCoarseToFineLinearInterp2D<T, Descriptor1, Descriptor2> const &rhs) :
    rescaleEngine(rhs.rescaleEngine->clone()),
    direction(rhs.direction),
    orientation(rhs.orientation)
{ }

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineLinearInterp2D<T, Descriptor1, Descriptor2>
    &CopyCoarseToFineLinearInterp2D<T, Descriptor1, Descriptor2>::operator=(
        CopyCoarseToFineLinearInterp2D<T, Descriptor1, Descriptor2> const &rhs)
{
    delete rescaleEngine;
    rescaleEngine = rhs.rescaleEngine->clone();
    direction = rhs.direction;
    orientation = rhs.orientation;
    return *this;
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
void CopyCoarseToFineLinearInterp2D<T, Descriptor1, Descriptor2>::process(
    Box2D coarseDomain, BlockLattice2D<T, Descriptor1> &coarseLattice,
    BlockLattice2D<T, Descriptor2> &fineLattice)
{
    const plint x0 = coarseDomain.x0;
    const plint y0 = coarseDomain.y0;
    const plint x1 = coarseDomain.x1;
    const plint y1 = coarseDomain.y1;

    PLB_PRECONDITION(x0 == x1 || y0 == y1);

    Array<plint, Descriptor1<T>::d> delta;
    delta.resetToZero();
    if (direction == 0)
        delta[1] = 1;
    else
        delta[0] = 1;
    plint whichTime = 1;

    std::vector<T> pop1;
    std::vector<T> pop2;

    CoarseToFineConverter2D<T, Descriptor1, Descriptor2> converter(coarseLattice, fineLattice);

    // Check if the extended domain (over which the interpolation is executed) exceeds either
    //   the coarse or the fine lattice, in which case interpolation is not executed in these
    //   areas.
    Dot2D coarseLeftExcess(x0 - delta[0], y0 - delta[1]);
    Dot2D fineLeftExcess(converter.fineX(x0) - delta[0], converter.fineY(y0) - delta[1]);
    Dot2D coarseRightExcess(x1 + delta[0], y1 + delta[1]);
    Dot2D fineRightExcess(converter.fineX(x1) + delta[0], converter.fineY(y1) + delta[1]);
    // leftOutOfRange is true if interpolations are out-of-range on the left side of the domain.
    bool leftOutOfRange =
        (!contained(coarseLeftExcess, coarseLattice.getBoundingBox())
         || !contained(fineLeftExcess, fineLattice.getBoundingBox()));
    // rightOutOfRange is true if interpolations are out-of-range on the right side of the domain.
    bool rightOutOfRange =
        (!contained(coarseRightExcess, coarseLattice.getBoundingBox())
         || !contained(fineRightExcess, fineLattice.getBoundingBox()));

    // STEP 1: Add left-most cell to the interpolation engine, and perform left-most interpolation
    //   unless interpolations are out-of-range on the left side of the domain.

    if (!leftOutOfRange) {
        rescaleEngine->scaleCoarseFine(coarseLattice.get(x0 - delta[0], y0 - delta[1]), pop1);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(x0, y0), pop2);
        linearInterpolation(
            pop1, pop2,
            converter.fineDynamics(x0, y0, -delta[0], -delta[1]).getDecomposedValues(whichTime));
    }

    // STEP 2: Perform copies and right-directed interpolations on each cell, except the last one
    //   (which must be treated individually to avoid out-of-range problems).

    for (plint iX = x0; iX <= x1 - delta[0]; ++iX) {
        for (plint iY = y0; iY <= y1 - delta[1]; ++iY) {
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY), pop1);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + delta[0], iY + delta[1]), pop2);
            converter.fineDynamics(iX, iY, 0, 0)
                .getDecomposedValues(whichTime)
                .assign(pop1.begin(), pop1.end());

            // linearInterpolation the unknown population
            linearInterpolation(
                pop1, pop2,
                converter.fineDynamics(iX, iY, delta[0], delta[1]).getDecomposedValues(whichTime));
        }
    }

    // STEP 3: Copy from the last cell, and perform right-most interpolation unless interpolations
    //   are out-of-range on the right side of the domain.
    rescaleEngine->scaleCoarseFine(coarseLattice.get(x1, y1), pop1);
    converter.fineDynamics(x1, y1, 0, 0)
        .getDecomposedValues(whichTime)
        .assign(pop1.begin(), pop1.end());

    if (!rightOutOfRange) {
        rescaleEngine->scaleCoarseFine(coarseLattice.get(x1, y1), pop1);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(x1 + delta[0], y1 + delta[1]), pop2);
        linearInterpolation(
            pop1, pop2,
            converter.fineDynamics(x1, y1, delta[0], delta[1]).getDecomposedValues(whichTime));
    }
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineLinearInterp2D<T, Descriptor1, Descriptor2>
    *CopyCoarseToFineLinearInterp2D<T, Descriptor1, Descriptor2>::clone() const
{
    return new CopyCoarseToFineLinearInterp2D(*this);
}

/* *************** Class CopyCoarseToFineBoundaryLinearInterp2D ********** */

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineBoundaryLinearInterp2D<T, Descriptor1, Descriptor2>::
    CopyCoarseToFineBoundaryLinearInterp2D(
        RescaleEngine<T, Descriptor1> *rescaleEngine_, plint direction_, plint orientation_,
        plint whereToInterpolate_) :
    rescaleEngine(rescaleEngine_),
    direction(direction_),
    orientation(orientation_),
    whereToInterpolate(whereToInterpolate_)
{ }

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineBoundaryLinearInterp2D<T, Descriptor1, Descriptor2>::
    CopyCoarseToFineBoundaryLinearInterp2D(
        CopyCoarseToFineBoundaryLinearInterp2D<T, Descriptor1, Descriptor2> const &rhs) :
    rescaleEngine(rhs.rescaleEngine->clone()),
    direction(rhs.direction),
    orientation(rhs.orientation),
    whereToInterpolate(rhs.whereToInterpolate)
{ }

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineBoundaryLinearInterp2D<T, Descriptor1, Descriptor2>
    &CopyCoarseToFineBoundaryLinearInterp2D<T, Descriptor1, Descriptor2>::operator=(
        CopyCoarseToFineBoundaryLinearInterp2D<T, Descriptor1, Descriptor2> const &rhs)
{
    delete rescaleEngine;
    rescaleEngine = rhs.rescaleEngine->clone();
    direction = rhs.direction;
    orientation = rhs.orientation;
    whereToInterpolate = rhs.whereToInterpolate;
    return *this;
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineBoundaryLinearInterp2D<
    T, Descriptor1, Descriptor2>::~CopyCoarseToFineBoundaryLinearInterp2D()
{
    delete rescaleEngine;
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
void CopyCoarseToFineBoundaryLinearInterp2D<T, Descriptor1, Descriptor2>::process(
    Box2D domain, BlockLattice2D<T, Descriptor1> &coarseLattice,
    BlockLattice2D<T, Descriptor2> &fineLattice)
{
    PLB_PRECONDITION(domain.x0 == domain.x1);
    PLB_PRECONDITION(domain.y0 == domain.y1);
    PLB_PRECONDITION(whereToInterpolate == 1 || whereToInterpolate == -1);

    plint iX = domain.x0;
    plint iY = domain.y0;

    Array<plint, Descriptor1<T>::d> delta;
    delta.resetToZero();
    if (direction == 0)
        delta[1] = whereToInterpolate;
    else
        delta[0] = whereToInterpolate;
    plint whichTime = 1;

    std::vector<T> pop1;
    std::vector<T> pop2;
    CoarseToFineConverter2D<T, Descriptor1, Descriptor2> converter(coarseLattice, fineLattice);

    // Check if the extended domain (over which the interpolation is executed) exceeds either
    //   the coarse or the fine lattice, in which case interpolation is not executed.
    Dot2D coarseExcess(iX + delta[0], iY + delta[1]);
    Dot2D fineExcess(converter.fineX(iX) + delta[0], converter.fineY(iY) + delta[1]);
    bool outOfRange =
        (!contained(coarseExcess, coarseLattice.getBoundingBox())
         || !contained(fineExcess, fineLattice.getBoundingBox()));

    // retrieve the values at x0,y0
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY), pop1);
    // copy them to the corresponding fine site
    converter.fineDynamics(iX, iY, 0, 0)
        .getDecomposedValues(whichTime)
        .assign(pop1.begin(), pop1.end());

    // check if we are not out of range
    if (!outOfRange) {
        // if it is not the case, linearInterpolation
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + delta[0], iY + delta[1]), pop2);
        linearInterpolation(
            pop1, pop2,
            converter.fineDynamics(iX, iY, delta[0], delta[1]).getDecomposedValues(whichTime));
    }
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineBoundaryLinearInterp2D<T, Descriptor1, Descriptor2>
    *CopyCoarseToFineBoundaryLinearInterp2D<T, Descriptor1, Descriptor2>::clone() const
{
    return new CopyCoarseToFineBoundaryLinearInterp2D(*this);
}

/* *************** Class CopyCoarseToFineCubicInterp2D ****************** */

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineCubicInterp2D<T, Descriptor1, Descriptor2>::CopyCoarseToFineCubicInterp2D(
    RescaleEngine<T, Descriptor1> *rescaleEngine_, plint direction_, plint orientation_) :
    rescaleEngine(rescaleEngine_), direction(direction_), orientation(orientation_)
{ }

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineCubicInterp2D<T, Descriptor1, Descriptor2>::~CopyCoarseToFineCubicInterp2D()
{
    delete rescaleEngine;
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineCubicInterp2D<T, Descriptor1, Descriptor2>::CopyCoarseToFineCubicInterp2D(
    CopyCoarseToFineCubicInterp2D<T, Descriptor1, Descriptor2> const &rhs) :
    rescaleEngine(rhs.rescaleEngine->clone()),
    direction(rhs.direction),
    orientation(rhs.orientation)
{ }

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineCubicInterp2D<T, Descriptor1, Descriptor2>
    &CopyCoarseToFineCubicInterp2D<T, Descriptor1, Descriptor2>::operator=(
        CopyCoarseToFineCubicInterp2D<T, Descriptor1, Descriptor2> const &rhs)
{
    delete rescaleEngine;
    rescaleEngine = rhs.rescaleEngine->clone();
    direction = rhs.direction;
    orientation = rhs.orientation;
    return *this;
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
void CopyCoarseToFineCubicInterp2D<T, Descriptor1, Descriptor2>::process(
    Box2D coarseDomain, BlockLattice2D<T, Descriptor1> &coarseLattice,
    BlockLattice2D<T, Descriptor2> &fineLattice)
{
    const plint x0 = coarseDomain.x0;
    const plint y0 = coarseDomain.y0;
    const plint x1 = coarseDomain.x1;
    const plint y1 = coarseDomain.y1;

    CoarseToFineConverter2D<T, Descriptor1, Descriptor2> converter(coarseLattice, fineLattice);

    // prepare the delta, keep in mind that it has Palabos BC conventiion
    PLB_PRECONDITION(x0 == x1 || y0 == y1);
    Array<plint, Descriptor1<T>::d> delta;
    delta.resetToZero();
    // recomputing the correct delta in the direction of refinement
    if (direction == 0)
        delta[1] = 1;
    else
        delta[0] = 1;
    plint whichTime = 1;

    std::vector<T> pop1;  // the leftmost
    std::vector<T> pop2;  // the current
    std::vector<T> pop3;  // the next
    std::vector<T> pop4;  // the rightmost

    // Check if the extended domain (over which the interpolation is executed) exceeds either
    //   the coarse or the fine lattice, in which case interpolation is not executed in these
    //   areas.
    Dot2D coarseLeftExcess(x0 - 2 * delta[0], y0 - 2 * delta[1]);  // 2 neighbors for the coarse
    Dot2D fineLeftExcess(
        converter.fineX(x0) - 2 * delta[0],
        converter.fineY(y0) - 2 * delta[1]);  // only one for the fine

    Dot2D coarseRightExcess(x1 + 2 * delta[0], y1 + 2 * delta[1]);
    Dot2D fineRightExcess(converter.fineX(x1) + 2 * delta[0], converter.fineY(y1) + 2 * delta[1]);

    // leftOutOfRange is true if interpolations are out-of-range on the left side of the domain.
    bool leftOutOfRange =
        (!contained(coarseLeftExcess, coarseLattice.getBoundingBox())
         || !contained(fineLeftExcess, fineLattice.getBoundingBox()));
    // rightOutOfRange is true if interpolations are out-of-range on the right side of the domain.
    bool rightOutOfRange =
        (!contained(coarseRightExcess, coarseLattice.getBoundingBox())
         || !contained(fineRightExcess, fineLattice.getBoundingBox()));

    // STEP 1: Check if we need to interpolate left of x0
    if (!leftOutOfRange) {
        rescaleEngine->scaleCoarseFine(
            coarseLattice.get(x0 - 2 * delta[0], y0 - 2 * delta[1]), pop1);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(x0 - delta[0], y0 - delta[1]), pop2);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(x0, y0), pop3);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(x0 + delta[0], y0 + delta[1]), pop4);

        cubicCenteredInterpolation(
            pop1, pop2, pop3, pop4,
            converter.fineDynamics(x0, y0, -delta[0], -delta[1]).getDecomposedValues(whichTime));
    }

    // STEP 2: Perform copies and right-directed interpolations on each cell, except the last one
    //   (which must be treated individually to avoid out-of-range problems).
    for (plint iX = x0; iX <= x1 - delta[0]; ++iX) {
        for (plint iY = y0; iY <= y1 - delta[1]; ++iY) {
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - delta[0], iY - delta[1]), pop1);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY), pop2);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + delta[0], iY + delta[1]), pop3);
            rescaleEngine->scaleCoarseFine(
                coarseLattice.get(iX + 2 * delta[0], iY + 2 * delta[1]), pop4);
            // copy the known value over iX,iY in fine coordinates
            converter.fineDynamics(iX, iY, 0, 0)
                .getDecomposedValues(whichTime)
                .assign(pop2.begin(), pop2.end());
            // linearInterpolation the unknown population
            cubicCenteredInterpolation(
                pop1, pop2, pop3, pop4,
                converter.fineDynamics(iX, iY, delta[0], delta[1]).getDecomposedValues(whichTime));
        }
    }
    // copy over the last site
    rescaleEngine->scaleCoarseFine(coarseLattice.get(x1, y1), pop1);
    converter.fineDynamics(x1, y1, 0, 0)
        .getDecomposedValues(whichTime)
        .assign(pop1.begin(), pop1.end());

    // STEP 3: Check if we need to interpolate at the rightmost cells
    if (!rightOutOfRange) {
        rescaleEngine->scaleCoarseFine(coarseLattice.get(x1 - delta[0], y1 - delta[1]), pop1);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(x1, y1), pop2);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(x1 + delta[0], y1 + delta[1]), pop3);
        rescaleEngine->scaleCoarseFine(
            coarseLattice.get(x1 + 2 * delta[0], y1 + 2 * delta[1]), pop4);

        cubicCenteredInterpolation(
            pop1, pop2, pop3, pop4,
            converter.fineDynamics(x1, y1, delta[0], delta[1]).getDecomposedValues(whichTime));
    }
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineCubicInterp2D<T, Descriptor1, Descriptor2>
    *CopyCoarseToFineCubicInterp2D<T, Descriptor1, Descriptor2>::clone() const
{
    return new CopyCoarseToFineCubicInterp2D(*this);
}

/* *************** Class CopyCoarseToFineCubicInterpWithDeconvolution2D ****************** */

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineCubicInterpWithDeconvolution2D<T, Descriptor1, Descriptor2>::
    CopyCoarseToFineCubicInterpWithDeconvolution2D(
        RescaleEngine<T, Descriptor1> *rescaleEngine_, plint direction_, plint orientation_) :
    rescaleEngine(rescaleEngine_), direction(direction_), orientation(orientation_)
{ }

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineCubicInterpWithDeconvolution2D<
    T, Descriptor1, Descriptor2>::~CopyCoarseToFineCubicInterpWithDeconvolution2D()
{
    delete rescaleEngine;
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineCubicInterpWithDeconvolution2D<T, Descriptor1, Descriptor2>::
    CopyCoarseToFineCubicInterpWithDeconvolution2D(
        CopyCoarseToFineCubicInterpWithDeconvolution2D<T, Descriptor1, Descriptor2> const &rhs) :
    rescaleEngine(rhs.rescaleEngine->clone()),
    direction(rhs.direction),
    orientation(rhs.orientation)
{ }

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineCubicInterpWithDeconvolution2D<T, Descriptor1, Descriptor2>
    &CopyCoarseToFineCubicInterpWithDeconvolution2D<T, Descriptor1, Descriptor2>::operator=(
        CopyCoarseToFineCubicInterpWithDeconvolution2D<T, Descriptor1, Descriptor2> const &rhs)
{
    delete rescaleEngine;
    rescaleEngine = rhs.rescaleEngine->clone();
    direction = rhs.direction;
    orientation = rhs.orientation;
    return *this;
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
void CopyCoarseToFineCubicInterpWithDeconvolution2D<T, Descriptor1, Descriptor2>::process(
    Box2D coarseDomain, BlockLattice2D<T, Descriptor1> &coarseLattice,
    BlockLattice2D<T, Descriptor2> &fineLattice)
{
    const plint x0 = coarseDomain.x0;
    const plint y0 = coarseDomain.y0;
    const plint x1 = coarseDomain.x1;
    const plint y1 = coarseDomain.y1;

    CoarseToFineConverter2D<T, Descriptor1, Descriptor2> converter(coarseLattice, fineLattice);

    // prepare the delta, keep in mind that it has Palabos BC conventiion
    PLB_PRECONDITION(x0 == x1 || y0 == y1);
    Array<plint, Descriptor1<T>::d> delta;
    delta.resetToZero();
    // recomputing the correct delta in the direction of refinement
    if (direction == 0)
        delta[1] = 1;
    else
        delta[0] = 1;
    plint whichTime = 1;

    std::vector<T> pop1;  // the leftmost
    std::vector<T> pop2;  // the current
    std::vector<T> pop3;  // the next
    std::vector<T> pop4;  // the rightmost

    // Check if the extended domain (over which the interpolation is executed) exceeds either
    //   the coarse or the fine lattice, in which case interpolation is not executed in these
    //   areas.
    Dot2D coarseLeftExcess(x0 - 2 * delta[0], y0 - 2 * delta[1]);  // 2 neighbors for the coarse
    Dot2D fineLeftExcess(
        converter.fineX(x0) - 2 * delta[0],
        converter.fineY(y0) - 2 * delta[1]);  // only one for the fine

    Dot2D coarseRightExcess(x1 + 2 * delta[0], y1 + 2 * delta[1]);
    Dot2D fineRightExcess(converter.fineX(x1) + 2 * delta[0], converter.fineY(y1) + 2 * delta[1]);

    // leftOutOfRange is true if interpolations are out-of-range on the left side of the domain.
    bool leftOutOfRange =
        (!contained(coarseLeftExcess, coarseLattice.getBoundingBox())
         || !contained(fineLeftExcess, fineLattice.getBoundingBox()));
    // rightOutOfRange is true if interpolations are out-of-range on the right side of the domain.
    bool rightOutOfRange =
        (!contained(coarseRightExcess, coarseLattice.getBoundingBox())
         || !contained(fineRightExcess, fineLattice.getBoundingBox()));

    // STEP 1: Check if we need to interpolate left of x0
    if (!leftOutOfRange) {
        rescaleEngine->scaleCoarseFine(
            coarseLattice.get(x0 - 2 * delta[0], y0 - 2 * delta[1]), pop1);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(x0 - delta[0], y0 - delta[1]), pop2);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(x0, y0), pop3);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(x0 + delta[0], y0 + delta[1]), pop4);

        cubicCenteredInterpolation(
            pop1, pop2, pop3, pop4,
            converter.fineDynamics(x0, y0, -delta[0], -delta[1]).getDecomposedValues(whichTime));

        //         plint fX = converter.fineX(x0) - delta[0] - orientation*(1-delta[0]);
        //         plint fY = converter.fineY(y0) - delta[1] - orientation*(1-delta[1]);
        //         std::vector<T> correction = computeDeconvolutionExtrapolation(fineLattice, fX,
        //         fY, delta[0], delta[1]); std::vector<T> original   =
        //         converter.fineDynamics(x0,y0,-delta[0],-delta[1]).getDecomposedValues(whichTime);
        //
        //         for (pluint iA = 0; iA < original.size(); ++iA) {
        //             original[iA] -= correction[iA];
        //         }
        //
        //         converter.fineDynamics(x0,y0,-delta[0],-delta[1]).getDecomposedValues(whichTime).assign(original.begin(),original.end());
    }

    // STEP 2: Perform copies and right-directed interpolations on each cell, except the last one
    //   (which must be treated individually to avoid out-of-range problems).
    for (plint iX = x0; iX <= x1 - delta[0]; ++iX) {
        for (plint iY = y0; iY <= y1 - delta[1]; ++iY) {
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - delta[0], iY - delta[1]), pop1);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY), pop2);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + delta[0], iY + delta[1]), pop3);
            rescaleEngine->scaleCoarseFine(
                coarseLattice.get(iX + 2 * delta[0], iY + 2 * delta[1]), pop4);
            // copy the known value over iX,iY in fine coordinates
            converter.fineDynamics(iX, iY, 0, 0)
                .getDecomposedValues(whichTime)
                .assign(pop2.begin(), pop2.end());

            //             plint fX = converter.fineX(iX) - orientation*(1-delta[0]);
            //             plint fY = converter.fineY(iY) - orientation*(1-delta[1]);
            //
            //             std::vector<T> correction =
            //             computeDeconvolutionExtrapolation(fineLattice, fX, fY, delta[0],
            //             delta[1]); std::vector<T> original   =
            //             converter.fineDynamics(iX,iY,0,0).getDecomposedValues(whichTime);
            //
            //             for (pluint iA = 0; iA < original.size(); ++iA) {
            //                 original[iA] -= correction[iA];
            //             }
            //             converter.fineDynamics(iX,iY,0,0).getDecomposedValues(whichTime).assign(original.begin(),original.end());

            // linearInterpolation the unknown population
            cubicCenteredInterpolation(
                pop1, pop2, pop3, pop4,
                converter.fineDynamics(iX, iY, delta[0], delta[1]).getDecomposedValues(whichTime));

            //             fX +=  delta[0]; fY += delta[1];
            //             correction = computeDeconvolutionExtrapolation(fineLattice, fX, fY,
            //             delta[0], delta[1]); original   =
            //             converter.fineDynamics(iX,iY,delta[0],delta[1]).getDecomposedValues(whichTime);
            //
            //             for (pluint iA = 0; iA < original.size(); ++iA) {
            //                 original[iA] -= correction[iA];
            //             }
            //             converter.fineDynamics(iX,iY,delta[0],delta[1]).getDecomposedValues(whichTime).assign(original.begin(),original.end());
        }
    }
    // copy over the last site
    rescaleEngine->scaleCoarseFine(coarseLattice.get(x1, y1), pop1);
    converter.fineDynamics(x1, y1, 0, 0)
        .getDecomposedValues(whichTime)
        .assign(pop1.begin(), pop1.end());

    //     plint fX = converter.fineX(x1) - orientation*(1-delta[0]);
    //     plint fY = converter.fineY(y1) - orientation*(1-delta[1]);
    //
    //     std::vector<T> correction = computeDeconvolutionExtrapolation(fineLattice, fX, fY,
    //     delta[0], delta[1]); std::vector<T> original   =
    //     converter.fineDynamics(x1,y1,0,0).getDecomposedValues(whichTime);
    //
    //     for (pluint iA = 0; iA < original.size(); ++iA) {
    //         original[iA] -= correction[iA];
    //     }
    //
    //     converter.fineDynamics(x1,y1,0,0).getDecomposedValues(whichTime).assign(original.begin(),original.end());

    // STEP 3: Check if we need to interpolate at the rightmost cells
    if (!rightOutOfRange) {
        rescaleEngine->scaleCoarseFine(coarseLattice.get(x1 - delta[0], y1 - delta[1]), pop1);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(x1, y1), pop2);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(x1 + delta[0], y1 + delta[1]), pop3);
        rescaleEngine->scaleCoarseFine(
            coarseLattice.get(x1 + 2 * delta[0], y1 + 2 * delta[1]), pop4);

        cubicCenteredInterpolation(
            pop1, pop2, pop3, pop4,
            converter.fineDynamics(x1, y1, delta[0], delta[1]).getDecomposedValues(whichTime));

        //         plint fX = converter.fineX(x1) - orientation*(1-delta[0]);
        //         plint fY = converter.fineY(y1) - orientation*(1-delta[1]);
        //         std::vector<T> correction = computeDeconvolutionExtrapolation(fineLattice, fX,
        //         fY, delta[0], delta[1]); std::vector<T> original   =
        //         converter.fineDynamics(x1,y1,delta[0],delta[1]).getDecomposedValues(whichTime);
        //
        //         for (pluint iA = 0; iA < original.size(); ++iA) {
        //             original[iA] -= correction[iA];
        //         }
        //         converter.fineDynamics(x1,y1,delta[0],delta[1]).getDecomposedValues(whichTime).assign(original.begin(),original.end());
    }
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineCubicInterpWithDeconvolution2D<T, Descriptor1, Descriptor2>
    *CopyCoarseToFineCubicInterpWithDeconvolution2D<T, Descriptor1, Descriptor2>::clone() const
{
    return new CopyCoarseToFineCubicInterpWithDeconvolution2D(*this);
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
std::vector<T> CopyCoarseToFineCubicInterpWithDeconvolution2D<T, Descriptor1, Descriptor2>::
    computeDeconvolutionExtrapolation(
        BlockLattice2D<T, Descriptor2> &fineLattice, plint iX, plint iY, plint dx, plint dy) const
{
    std::vector<T> currentDecomp;
    std::vector<std::vector<T> > filteredDecomp;

    fineLattice.get(iX + dx, iY + dy)
        .getDynamics()
        .decompose(
            fineLattice.get(iX + dx, iY + dy), currentDecomp,
            rescaleEngine->getDecompositionOrder());
    filteredDecomp.push_back(currentDecomp);
    fineLattice.get(iX - dx, iY - dy)
        .getDynamics()
        .decompose(
            fineLattice.get(iX - dx, iY - dy), currentDecomp,
            rescaleEngine->getDecompositionOrder());
    filteredDecomp.push_back(currentDecomp);

    fineLattice.get(iX, iY).getDynamics().decompose(
        fineLattice.get(iX, iY), currentDecomp, rescaleEngine->getDecompositionOrder());

    std::vector<T> refDecomp = currentDecomp;

    for (pluint iA = 0; iA < filteredDecomp.size(); ++iA) {
        for (pluint iB = 0; iB < filteredDecomp[iA].size(); ++iB) {
            currentDecomp[iB] += filteredDecomp[iA][iB];
        }
    }

    for (pluint iA = 0; iA < refDecomp.size(); ++iA) {
        currentDecomp[iA] /= ((T)filteredDecomp.size() + (T)1);

        refDecomp[iA] -= currentDecomp[iA];
        refDecomp[iA] *= (T)0.8;
    }

    return refDecomp;
}

/* *************** Class CopyCoarseToFineBoundaryCubicInterp2D ************* */

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineBoundaryCubicInterp2D<T, Descriptor1, Descriptor2>::
    CopyCoarseToFineBoundaryCubicInterp2D(
        RescaleEngine<T, Descriptor1> *rescaleEngine_, plint direction_, plint orientation_,
        plint whereToInterpolate_) :
    rescaleEngine(rescaleEngine_),
    direction(direction_),
    orientation(orientation_),
    whereToInterpolate(whereToInterpolate_)
{ }

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineBoundaryCubicInterp2D<T, Descriptor1, Descriptor2>::
    CopyCoarseToFineBoundaryCubicInterp2D(
        CopyCoarseToFineBoundaryCubicInterp2D<T, Descriptor1, Descriptor2> const &rhs) :
    rescaleEngine(rhs.rescaleEngine->clone()),
    direction(rhs.direction),
    orientation(rhs.orientation),
    whereToInterpolate(rhs.whereToInterpolate)
{ }

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineBoundaryCubicInterp2D<T, Descriptor1, Descriptor2>
    &CopyCoarseToFineBoundaryCubicInterp2D<T, Descriptor1, Descriptor2>::operator=(
        CopyCoarseToFineBoundaryCubicInterp2D<T, Descriptor1, Descriptor2> const &rhs)
{
    delete rescaleEngine;
    rescaleEngine = rhs.rescaleEngine->clone();
    direction = rhs.direction;
    orientation = rhs.orientation;
    whereToInterpolate = rhs.whereToInterpolate;
    return *this;
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineBoundaryCubicInterp2D<
    T, Descriptor1, Descriptor2>::~CopyCoarseToFineBoundaryCubicInterp2D()
{
    delete rescaleEngine;
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
void CopyCoarseToFineBoundaryCubicInterp2D<T, Descriptor1, Descriptor2>::process(
    Box2D domain, BlockLattice2D<T, Descriptor1> &coarseLattice,
    BlockLattice2D<T, Descriptor2> &fineLattice)
{
    PLB_PRECONDITION(domain.x0 == domain.x1 || domain.y0 == domain.y1);
    PLB_PRECONDITION(whereToInterpolate == 1 || whereToInterpolate == -1);

    plint iX = domain.x0;
    plint iY = domain.y0;

    Array<plint, Descriptor1<T>::d> delta;
    delta.resetToZero();
    if (direction == 0)
        delta[1] = whereToInterpolate;
    else
        delta[0] = whereToInterpolate;
    plint whichTime = 1;

    std::vector<T> pop1;
    std::vector<T> pop2;
    std::vector<T> pop3;

    CoarseToFineConverter2D<T, Descriptor1, Descriptor2> converter(coarseLattice, fineLattice);

    // retrieve the values at x0,y0 and neighbors
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY), pop1);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + delta[0], iY + delta[1]), pop2);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2 * delta[0], iY + 2 * delta[1]), pop3);

    // copy (x0,y0) values from coarse to the corresponding fine site
    converter.fineDynamics(iX, iY, 0, 0)
        .getDecomposedValues(whichTime)
        .assign(pop1.begin(), pop1.end());
    // cubic decentered interpolation for the fine site next to (x0,y0)
    quadraticNonCenteredInterpolation(
        pop1, pop2, pop3,
        converter.fineDynamics(iX, iY, delta[0], delta[1]).getDecomposedValues(whichTime));
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineBoundaryCubicInterp2D<T, Descriptor1, Descriptor2>
    *CopyCoarseToFineBoundaryCubicInterp2D<T, Descriptor1, Descriptor2>::clone() const
{
    return new CopyCoarseToFineBoundaryCubicInterp2D(*this);
}

/* *************** Class CopyCoarseToFineBoundaryHelper2D ************* */

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineBoundaryHelper2D<T, Descriptor1, Descriptor2>::CopyCoarseToFineBoundaryHelper2D(
    RescaleEngine<T, Descriptor1> *rescaleEngine_, plint direction_, plint orientation_,
    plint whereToInterpolate_) :
    rescaleEngine(rescaleEngine_),
    direction(direction_),
    orientation(orientation_),
    whereToInterpolate(whereToInterpolate_)
{ }

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineBoundaryHelper2D<T, Descriptor1, Descriptor2>::CopyCoarseToFineBoundaryHelper2D(
    CopyCoarseToFineBoundaryHelper2D<T, Descriptor1, Descriptor2> const &rhs) :
    rescaleEngine(rhs.rescaleEngine->clone()),
    direction(rhs.direction),
    orientation(rhs.orientation),
    whereToInterpolate(rhs.whereToInterpolate)
{ }

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineBoundaryHelper2D<T, Descriptor1, Descriptor2>
    &CopyCoarseToFineBoundaryHelper2D<T, Descriptor1, Descriptor2>::operator=(
        CopyCoarseToFineBoundaryHelper2D<T, Descriptor1, Descriptor2> const &rhs)
{
    delete rescaleEngine;
    rescaleEngine = rhs.rescaleEngine->clone();
    direction = rhs.direction;
    orientation = rhs.orientation;
    whereToInterpolate = rhs.whereToInterpolate;
    return *this;
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineBoundaryHelper2D<T, Descriptor1, Descriptor2>::~CopyCoarseToFineBoundaryHelper2D()
{
    delete rescaleEngine;
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
void CopyCoarseToFineBoundaryHelper2D<T, Descriptor1, Descriptor2>::process(
    Box2D domain, BlockLattice2D<T, Descriptor1> &coarseLattice,
    BlockLattice2D<T, Descriptor2> &fineLattice)
{
    PLB_PRECONDITION(domain.x0 == domain.x1 || domain.y0 == domain.y1);
    PLB_PRECONDITION(whereToInterpolate == 1 || whereToInterpolate == -1);

    plint iX = domain.x0;
    plint iY = domain.y0;

    Array<plint, Descriptor1<T>::d> delta;
    delta.resetToZero();
    if (direction == 0)
        delta[1] = whereToInterpolate;
    else
        delta[0] = whereToInterpolate;
    plint whichTime = 1;

    std::vector<T> pop1;
    std::vector<T> pop2;
    std::vector<T> pop3;

    CoarseToFineConverter2D<T, Descriptor1, Descriptor2> converter(coarseLattice, fineLattice);

    // retrieve the values at x0,y0 and neighbors
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - delta[0], iY - delta[1]), pop1);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY), pop2);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + delta[0], iY + delta[1]), pop3);

    // cubic decentered interpolation for the fine site next to (x0,y0)
    quadraticNonCenteredInterpolation(
        pop1, pop2, pop3,
        converter.fineDynamics(iX, iY, -delta[0], -delta[1]).getDecomposedValues(whichTime));
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyCoarseToFineBoundaryHelper2D<T, Descriptor1, Descriptor2>
    *CopyCoarseToFineBoundaryHelper2D<T, Descriptor1, Descriptor2>::clone() const
{
    return new CopyCoarseToFineBoundaryHelper2D(*this);
}

/* *************** Class Copy_t1_to_t0_2D *********************************** */
template <typename T, template <typename U> class Descriptor>
Copy_t1_to_t0_2D<T, Descriptor>::Copy_t1_to_t0_2D(plint numTimeSteps_, plint executionTime_) :
    numTimeSteps(numTimeSteps_), executionTime(executionTime_)
{ }

template <typename T, template <typename U> class Descriptor>
void Copy_t1_to_t0_2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &fineLattice)
{
    // Determine current relative value of time steps.
    size_t relativeTime = fineLattice.getTimeCounter().getTime() % numTimeSteps;
    // Execute data processor only if one is at the end of a cycle.
    if ((plint)relativeTime == executionTime) {
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                FineGridBoundaryDynamics<T, Descriptor> &fineDynamics =
                    dynamic_cast<FineGridBoundaryDynamics<T, Descriptor> &>(
                        fineLattice.get(iX, iY).getDynamics());
                std::vector<T> &t0 = fineDynamics.getDecomposedValues(0);
                std::vector<T> &t1 = fineDynamics.getDecomposedValues(1);
                t0.assign(t1.begin(), t1.end());
                // complete populations
                // TODO: erase this and create a data processor that executes this only when needed
                fineDynamics.completePopulations(fineLattice.get(iX, iY));
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
Copy_t1_to_t0_2D<T, Descriptor> *Copy_t1_to_t0_2D<T, Descriptor>::clone() const
{
    return new Copy_t1_to_t0_2D(*this);
}

}  // namespace plb

#endif  // FINE_GRID_PROCESSORS_2D_HH
