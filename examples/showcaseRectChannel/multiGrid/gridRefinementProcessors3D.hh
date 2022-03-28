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
 * Dynamics and data processors used to implement 3D grid refinement -- generic implementation.
 */
#ifndef GRID_REFINEMENT_PROCESSORS_3D_HH
#define GRID_REFINEMENT_PROCESSORS_3D_HH

#include "core/geometry3D.h"
#include "finiteDifference/fdStencils1D.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"
#include "multiGrid/gridRefinementProcessors3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
void InterpolateOverLineAndExcess(
    Box3D domain, BlockLattice3D<T, Descriptor> &coarseLattice,
    BlockLattice3D<T, Descriptor> &fineLattice,
    InterpolationEngine2D<T, Descriptor> *interpolationEngine,
    RescaleEngine<T, Descriptor> *rescaleEngine, std::vector<Dot3D> &delta)
{
    const plint x0 = domain.x0;
    const plint y0 = domain.y0;
    const plint z0 = domain.z0;
    const plint x1 = domain.x1;
    const plint y1 = domain.y1;
    const plint z1 = domain.z1;

    plint whichTime = 1;
    CoarseToFineConverter3D<T, Descriptor> converter(coarseLattice, fineLattice);

    // delta[0] especifies the direction of the interpolation
    // delta[1] especifies the direction perpendicular to the interpolation
    // delta[2] is the sum of the other two deltas
    delta.push_back(
        Dot3D(delta[0].x + delta[1].x, delta[0].y + delta[1].y, delta[0].z + delta[1].z));

    // Check if the extended domain (over which the interpolation is executed) exceeds either
    //   the coarse or the fine lattice, in which case interpolation is not executed.
    Dot3D coarseExcess(x0 + delta[1].x, y0 + delta[1].y, z0 + delta[1].z);
    Dot3D fineExcess(
        converter.fineX(x0) + delta[1].x, converter.fineY(y0) + delta[1].y,
        converter.fineZ(z0) + delta[1].z);
    bool outOfRange =
        (!contained(coarseExcess, coarseLattice.getBoundingBox())
         || !contained(fineExcess, fineLattice.getBoundingBox()));

    Dot3D startExcess(x0 - delta[0].x, y0 - delta[0].y, z0 - delta[0].z);
    Dot3D startFineExcess(
        converter.fineX(x0) - delta[0].x, converter.fineY(y0) - delta[0].y,
        converter.fineZ(z0) - delta[0].z);
    bool startOutOfRange =
        (!contained(startExcess, coarseLattice.getBoundingBox())
         || !contained(startFineExcess, fineLattice.getBoundingBox()));
    Dot3D endExcess(x1 + delta[0].x, y1 + delta[0].y, z1 + delta[0].z);
    Dot3D endFineExcess(
        converter.fineX(x1) + delta[0].x, converter.fineY(y1) + delta[0].y,
        converter.fineZ(z1) + delta[0].z);
    bool endOutOfRange =
        (!contained(endExcess, coarseLattice.getBoundingBox())
         || !contained(endFineExcess, fineLattice.getBoundingBox()));

    if (!startOutOfRange) {
        interpolationEngine->pushToCouple(
            Array<Cell<T, Descriptor>, 2>(
                coarseLattice.get(x0 - delta[0].x, y0 - delta[0].y, z0 - delta[0].z),
                coarseLattice.get(x0, y0, z0)),
            *rescaleEngine);
    }

    if (outOfRange && !startOutOfRange) {
        interpolationEngine->pushToCouple();
    } else {
        if (!startOutOfRange && !outOfRange) {
            interpolationEngine->pushToCouple(
                Array<Cell<T, Descriptor>, 2>(
                    coarseLattice.get(
                        x0 - delta[0].x + delta[1].x, y0 - delta[0].y + delta[1].y,
                        z0 - delta[0].z + delta[1].z),
                    coarseLattice.get(x0 + delta[1].x, y0 + delta[1].y, z0 + delta[1].z)),
                *rescaleEngine);

            if (startOutOfRange) {
                interpolationEngine->pushToCouple();
            }
        }
    }

    if (!startOutOfRange) {
        interpolationEngine->copyFrom(
            converter.fineDynamics(x0, y0, z0, 0, 0, 0).getDecomposedValues(whichTime), 0);
        interpolationEngine->interpolateFromCouple(
            converter.fineDynamics(x0, y0, z0, -delta[0].x, -delta[0].y, -delta[0].z)
                .getDecomposedValues(whichTime),
            0, 1);
    }

    if (!outOfRange) {
        if (!startOutOfRange) {
            interpolationEngine->interpolateFromCouple(
                converter
                    .fineDynamics(
                        x0, y0, z0, -delta[0].x + delta[1].x, -delta[0].y + delta[1].y,
                        -delta[0].z + delta[1].z)
                    .getDecomposedValues(whichTime),
                0, 2);
            interpolationEngine->interpolateFromAll(
                converter.fineDynamics(x0, y0, z0, delta[1].x, delta[1].y, delta[1].z)
                    .getDecomposedValues(whichTime));
        } else {
            interpolationEngine->interpolateFromCouple(
                converter.fineDynamics(x0, y0, z0, delta[1].x, delta[1].y, delta[1].z)
                    .getDecomposedValues(whichTime),
                0, 1);
        }
    }

    for (plint iX = x0; iX <= x1 - delta[0].x; ++iX) {
        for (plint iY = y0; iY <= y1 - delta[0].y; ++iY) {
            for (plint iZ = z0; iZ <= z1 - delta[0].z; ++iZ) {
                interpolationEngine->pushToCouple(
                    Array<Cell<T, Descriptor>, 2>(
                        coarseLattice.get(iX, iY, iZ),
                        coarseLattice.get(iX + delta[0].x, iY + delta[0].y, iZ + delta[0].z)),
                    *rescaleEngine);

                if (outOfRange) {
                    interpolationEngine->pushToCouple();
                } else {
                    interpolationEngine->pushToCouple(
                        Array<Cell<T, Descriptor>, 2>(
                            coarseLattice.get(iX + delta[1].x, iY + delta[1].y, iZ + delta[1].z),
                            coarseLattice.get(iX + delta[2].x, iY + delta[2].y, iZ + delta[2].z)),
                        *rescaleEngine);
                }

                interpolationEngine->copyFrom(
                    converter.fineDynamics(iX, iY, iZ, 0, 0, 0).getDecomposedValues(whichTime), 0);
                interpolationEngine->copyFrom(
                    converter
                        .fineDynamics(iX, iY, iZ, 2 * delta[0].x, 2 * delta[0].y, 2 * delta[0].z)
                        .getDecomposedValues(whichTime),
                    1);
                interpolationEngine->interpolateFromCouple(
                    converter.fineDynamics(iX, iY, iZ, delta[0].x, delta[0].y, delta[0].z)
                        .getDecomposedValues(whichTime),
                    0, 1);
                if (!outOfRange) {
                    interpolationEngine->interpolateFromCouple(
                        converter.fineDynamics(iX, iY, iZ, delta[1].x, delta[1].y, delta[1].z)
                            .getDecomposedValues(whichTime),
                        0, 2);
                    interpolationEngine->interpolateFromAll(
                        converter.fineDynamics(iX, iY, iZ, delta[2].x, delta[2].y, delta[2].z)
                            .getDecomposedValues(whichTime));
                }
            }
        }
    }

    if (!endOutOfRange) {
        interpolationEngine->pushToCouple(
            Array<Cell<T, Descriptor>, 2>(
                coarseLattice.get(x1, y1, z1),
                coarseLattice.get(x1 + delta[0].x, y1 + delta[0].y, z1 + delta[0].z)),
            *rescaleEngine);
    }

    if (outOfRange && !endOutOfRange) {
        interpolationEngine->pushToCouple();
    } else {
        if (!endOutOfRange && !outOfRange) {
            interpolationEngine->pushToCouple(
                Array<Cell<T, Descriptor>, 2>(
                    coarseLattice.get(x1 + delta[1].x, y1 + delta[1].y, z1 + delta[1].z),
                    coarseLattice.get(
                        x1 + delta[0].x + delta[1].x, y1 + delta[0].y + delta[1].y,
                        z1 + delta[0].z + delta[1].z)),
                *rescaleEngine);

            if (endOutOfRange) {
                interpolationEngine->pushToCouple();
            }
        }
    }

    if (!endOutOfRange) {
        interpolationEngine->copyFrom(
            converter.fineDynamics(x1, y1, z1, 0, 0, 0).getDecomposedValues(whichTime), 0);
        interpolationEngine->interpolateFromCouple(
            converter.fineDynamics(x1, y1, z1, delta[0].x, delta[0].y, delta[0].z)
                .getDecomposedValues(whichTime),
            0, 1);
    }

    if (!outOfRange) {
        if (!endOutOfRange) {
            interpolationEngine->interpolateFromCouple(
                converter.fineDynamics(x1, y1, z1, delta[1].x, delta[1].y, delta[1].z)
                    .getDecomposedValues(whichTime),
                0, 2);
            interpolationEngine->interpolateFromAll(
                converter.fineDynamics(x1, y1, z1, delta[2].x, delta[2].y, delta[2].z)
                    .getDecomposedValues(whichTime));
        } else {
            interpolationEngine->interpolateFromCouple(
                converter.fineDynamics(x1, y1, z1, delta[1].x, delta[1].y, delta[1].z)
                    .getDecomposedValues(whichTime),
                0, 1);
        }
    }

    delete rescaleEngine;
    delete interpolationEngine;
}

std::vector<Box3D> createCorners(Box3D domain, plint direction)
{
    std::vector<Box3D> result;
    if (direction == 0) {  // the YZ plane
        result.push_back(
            Box3D(domain.x0, domain.x1, domain.y0, domain.y0, domain.z1, domain.z1));  // CUL
        result.push_back(
            Box3D(domain.x0, domain.x1, domain.y0, domain.y0, domain.z0, domain.z0));  // CLL
        result.push_back(
            Box3D(domain.x0, domain.x1, domain.y1, domain.y1, domain.z0, domain.z0));  // CLR
        result.push_back(
            Box3D(domain.x0, domain.x1, domain.y1, domain.y1, domain.z1, domain.z1));  // CUR
    }
    if (direction == 1) {  // the XZ plane
        result.push_back(Box3D(domain.x0, domain.x0, domain.y0, domain.y1, domain.z1, domain.z1));
        result.push_back(Box3D(domain.x0, domain.x0, domain.y0, domain.y1, domain.z0, domain.z0));
        result.push_back(Box3D(domain.x1, domain.x1, domain.y0, domain.y1, domain.z0, domain.z0));
        result.push_back(Box3D(domain.x1, domain.x1, domain.y0, domain.y1, domain.z1, domain.z1));
    }
    if (direction == 2) {  // the XY plane
        result.push_back(Box3D(domain.x0, domain.x0, domain.y1, domain.y1, domain.z0, domain.z1));
        result.push_back(Box3D(domain.x0, domain.x0, domain.y0, domain.y0, domain.z0, domain.z1));
        result.push_back(Box3D(domain.x1, domain.x1, domain.y0, domain.y0, domain.z0, domain.z1));
        result.push_back(Box3D(domain.x1, domain.x1, domain.y1, domain.y1, domain.z0, domain.z1));
    }

    return result;
}

Box3D createLeftDomain(Box3D originalDomain, plint direction)
{
    Box3D result;
    if (direction == 0) {  // the YZ plane
        result = Box3D(
            originalDomain.x0, originalDomain.x1, originalDomain.y0, originalDomain.y0,
            originalDomain.z0, originalDomain.z1);
    }
    if (direction == 1) {  // the XZ plane
        result = Box3D(
            originalDomain.x0, originalDomain.x0, originalDomain.y0, originalDomain.y1,
            originalDomain.z0, originalDomain.z1);
    }
    if (direction == 2) {  // the XY plane
        result = Box3D(
            originalDomain.x0, originalDomain.x0, originalDomain.y0, originalDomain.y1,
            originalDomain.z0, originalDomain.z1);
    }
    return result;
}

Box3D createRightDomain(Box3D originalDomain, plint direction)
{
    Box3D result;
    if (direction == 0) {  // the YZ plane
        result = Box3D(
            originalDomain.x0, originalDomain.x1, originalDomain.y1, originalDomain.y1,
            originalDomain.z0, originalDomain.z1);
    }
    if (direction == 1) {  // the XZ plane
        result = Box3D(
            originalDomain.x1, originalDomain.x1, originalDomain.y0, originalDomain.y1,
            originalDomain.z0, originalDomain.z1);
    }
    if (direction == 2) {  // the XY plane
        result = Box3D(
            originalDomain.x1, originalDomain.x1, originalDomain.y0, originalDomain.y1,
            originalDomain.z0, originalDomain.z1);
    }
    return result;
}

Box3D createBottomDomain(Box3D originalDomain, plint direction)
{
    Box3D result;
    if (direction == 0) {  // the YZ plane
        result = Box3D(
            originalDomain.x0, originalDomain.x1, originalDomain.y0, originalDomain.y1,
            originalDomain.z0, originalDomain.z0);
    }
    if (direction == 1) {  // the XZ plane
        result = Box3D(
            originalDomain.x0, originalDomain.x1, originalDomain.y0, originalDomain.y1,
            originalDomain.z0, originalDomain.z0);
    }
    if (direction == 2) {  // the XY plane
        result = Box3D(
            originalDomain.x0, originalDomain.x1, originalDomain.y0, originalDomain.y0,
            originalDomain.z0, originalDomain.z1);
    }
    return result;
}

Box3D createTopDomain(Box3D originalDomain, plint direction)
{
    Box3D result;
    if (direction == 0) {  // the YZ plane
        result = Box3D(
            originalDomain.x0, originalDomain.x1, originalDomain.y0, originalDomain.y1,
            originalDomain.z1, originalDomain.z1);
    }
    if (direction == 1) {  // the XZ plane
        result = Box3D(
            originalDomain.x0, originalDomain.x1, originalDomain.y0, originalDomain.y1,
            originalDomain.z1, originalDomain.z1);
    }
    if (direction == 2) {  // the XY plane
        result = Box3D(
            originalDomain.x0, originalDomain.x1, originalDomain.y1, originalDomain.y1,
            originalDomain.z0, originalDomain.z1);
    }
    return result;
}

/* *************** Class InterpolateCoarseToFineDynamics3D ****************** */

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
InterpolateCoarseToFineDynamics3D<T, Descriptor1, Descriptor2>::InterpolateCoarseToFineDynamics3D(
    RescaleEngine<T, Descriptor1> *rescaleEngine_,
    InterpolationEngine2D<T, Descriptor1> *interpolationEngine_, plint direction_,
    const Dot3D delta_[2]) :
    rescaleEngine(rescaleEngine_), interpolationEngine(interpolationEngine_), direction(direction_)
{
    delta.push_back(delta_[0]);
    delta.push_back(delta_[1]);
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
InterpolateCoarseToFineDynamics3D<T, Descriptor1, Descriptor2>::~InterpolateCoarseToFineDynamics3D()
{
    delete rescaleEngine;
    delete interpolationEngine;
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
InterpolateCoarseToFineDynamics3D<T, Descriptor1, Descriptor2>::InterpolateCoarseToFineDynamics3D(
    InterpolateCoarseToFineDynamics3D<T, Descriptor1, Descriptor2> const &rhs) :
    rescaleEngine(rhs.rescaleEngine->clone()),
    interpolationEngine(rhs.interpolationEngine->clone()),
    direction(rhs.direction),
    delta(rhs.delta)
{ }

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
InterpolateCoarseToFineDynamics3D<T, Descriptor1, Descriptor2>
    &InterpolateCoarseToFineDynamics3D<T, Descriptor1, Descriptor2>::operator=(
        InterpolateCoarseToFineDynamics3D<T, Descriptor1, Descriptor2> const &rhs)
{
    delete rescaleEngine;
    rescaleEngine = rhs.rescaleEngine->clone();
    delete interpolationEngine;
    interpolationEngine = rhs.interpolationEngine->clone();
    direction = rhs.direction;
    delta.assign(rhs.delta.begin(), rhs.delta.end());
    return *this;
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
void InterpolateCoarseToFineDynamics3D<T, Descriptor1, Descriptor2>::process(
    Box3D coarseDomain, BlockLattice3D<T, Descriptor1> &coarseLattice,
    BlockLattice3D<T, Descriptor2> &fineLattice)
{
    const plint x0 = coarseDomain.x0;
    const plint y0 = coarseDomain.y0;
    const plint z0 = coarseDomain.z0;
    const plint x1 = coarseDomain.x1;
    const plint y1 = coarseDomain.y1;
    const plint z1 = coarseDomain.z1;

    PLB_PRECONDITION(x0 == x1 || y0 == y1 || z0 == z1);
    PLB_PRECONDITION(
        delta[0].x + delta[0].y + delta[0].z == 1 && delta[1].x + delta[1].y + delta[1].z == 1);
    // creation of delta[2]
    delta.push_back(
        Dot3D(delta[0].x + delta[1].x, delta[0].y + delta[1].y, delta[0].z + delta[1].z));

    CoarseToFineConverter3D<T, Descriptor1> converter(coarseLattice, fineLattice);
    plint whichTime = 1;

    // Check if the extended domain (over which the interpolation is executed) exceeds either
    //   the coarse or the fine lattice, in which case interpolation is not executed in these
    //   areas.
    Dot3D coarseLeftExcess(x0 - delta[0].x, y0 - delta[0].y, z0 - delta[0].z);
    Dot3D fineLeftExcess(
        converter.fineX(x0) - delta[0].x, converter.fineY(y0) - delta[0].y,
        converter.fineZ(z0) - delta[0].z);
    Dot3D coarseRightExcess(x1 + delta[0].x, y1 + delta[0].y, z1 + delta[0].z);
    Dot3D fineRightExcess(
        converter.fineX(x1) + delta[0].x, converter.fineY(y1) + delta[0].y,
        converter.fineZ(z1) + delta[0].z);
    Dot3D coarseBottomExcess(x0 - delta[1].x, y0 - delta[1].y, z0 - delta[1].z);
    Dot3D fineBottomExcess(
        converter.fineX(x0) - delta[1].x, converter.fineY(y0) - delta[1].y,
        converter.fineZ(z0) - delta[1].z);
    Dot3D coarseTopExcess(x1 + delta[1].x, y1 + delta[1].y, z1 + delta[1].z);
    Dot3D fineTopExcess(
        converter.fineX(x1) + delta[1].x, converter.fineY(y1) + delta[1].y,
        converter.fineZ(z1) + delta[1].z);

    // leftOutOfRange is true if interpolations are out-of-range on the left side of the domain.
    bool leftOutOfRange =
        (!contained(coarseLeftExcess, coarseLattice.getBoundingBox())
         || !contained(fineLeftExcess, fineLattice.getBoundingBox()));
    // rightOutOfRange is true if interpolations are out-of-range on the right side of the domain.
    bool rightOutOfRange =
        (!contained(coarseRightExcess, coarseLattice.getBoundingBox())
         || !contained(fineRightExcess, fineLattice.getBoundingBox()));
    // bottomOutOfRange is true if interpolations are out-of-range on the bottom of the domain.
    bool bottomOutOfRange =
        (!contained(coarseBottomExcess, coarseLattice.getBoundingBox())
         || !contained(fineBottomExcess, fineLattice.getBoundingBox()));
    // topOutOfRange is true if interpolations are out-of-range on the top of the domain.
    bool topOutOfRange =
        (!contained(coarseTopExcess, coarseLattice.getBoundingBox())
         || !contained(fineTopExcess, fineLattice.getBoundingBox()));

    // STEP 1: We take 4 sites of the coarse grid and we interpolate all the fine values inside
    for (plint iX = x0; iX <= x1 - delta[2].x; ++iX) {
        for (plint iY = y0; iY <= y1 - delta[2].y; ++iY) {
            for (plint iZ = z0; iZ <= z1 - delta[2].z; ++iZ) {
                // inserting the corresponding coarse values inside the 2D
                // interpolation engine
                interpolationEngine->pushToCouple(
                    Array<Cell<T, Descriptor1>, 2>(
                        coarseLattice.get(iX, iY, iZ),
                        coarseLattice.get(iX + delta[0].x, iY + delta[0].y, iZ + delta[0].z)),
                    *rescaleEngine);
                interpolationEngine->pushToCouple(
                    Array<Cell<T, Descriptor1>, 2>(
                        coarseLattice.get(iX + delta[1].x, iY + delta[1].y, iZ + delta[1].z),
                        coarseLattice.get(iX + delta[2].x, iY + delta[2].y, iZ + delta[2].z)),
                    *rescaleEngine);

                // copy of the values from the coarse to the fine
                interpolationEngine->copyFrom(
                    converter.fineDynamics(iX, iY, iZ, 0, 0, 0).getDecomposedValues(whichTime), 0);
                interpolationEngine->copyFrom(
                    converter
                        .fineDynamics(iX + delta[0].x, iY + delta[0].y, iZ + delta[0].z, 0, 0, 0)
                        .getDecomposedValues(whichTime),
                    1);
                interpolationEngine->copyFrom(
                    converter
                        .fineDynamics(iX + delta[1].x, iY + delta[1].y, iZ + delta[1].z, 0, 0, 0)
                        .getDecomposedValues(whichTime),
                    2);
                interpolationEngine->copyFrom(
                    converter
                        .fineDynamics(iX + delta[2].x, iY + delta[2].y, iZ + delta[2].z, 0, 0, 0)
                        .getDecomposedValues(whichTime),
                    3);

                // interpolate the 5 missing values
                interpolationEngine->interpolateFromCouple(
                    converter.fineDynamics(iX, iY, iZ, delta[0].x, delta[0].y, delta[0].z)
                        .getDecomposedValues(whichTime),
                    0, 1);
                interpolationEngine->interpolateFromCouple(
                    converter.fineDynamics(iX, iY, iZ, delta[1].x, delta[1].y, delta[1].z)
                        .getDecomposedValues(whichTime),
                    0, 2);
                interpolationEngine->interpolateFromCouple(
                    converter
                        .fineDynamics(
                            iX + delta[0].x, iY + delta[0].y, iZ + delta[0].z, delta[1].x,
                            delta[1].y, delta[1].z)
                        .getDecomposedValues(whichTime),
                    1, 3);
                interpolationEngine->interpolateFromCouple(
                    converter
                        .fineDynamics(
                            iX + delta[1].x, iY + delta[1].y, iZ + delta[1].z, delta[0].x,
                            delta[0].y, delta[0].z)
                        .getDecomposedValues(whichTime),
                    2, 3);
                // this interpolation uses the four coarse values stored
                interpolationEngine->interpolateFromAll(
                    converter.fineDynamics(iX, iY, iZ, delta[2].x, delta[2].y, delta[2].z)
                        .getDecomposedValues(whichTime));
            }
        }
    }

    // STEP 2: We perform interpolations over the excedent regions of the grid if needed
    if (!leftOutOfRange) {
        std::vector<Dot3D> deltas;
        deltas.push_back(Dot3D(delta[1].x, delta[1].y, delta[1].z));
        deltas.push_back(Dot3D(-delta[0].x, -delta[0].y, -delta[0].z));
        InterpolateOverLineAndExcess<T, Descriptor1>(
            createLeftDomain(coarseDomain, direction), coarseLattice, fineLattice,
            interpolationEngine->clone(), rescaleEngine->clone(), deltas);
    }

    if (!bottomOutOfRange) {
        std::vector<Dot3D> deltas;
        deltas.push_back(Dot3D(delta[0].x, delta[0].y, delta[0].z));
        deltas.push_back(Dot3D(-delta[1].x, -delta[1].y, -delta[1].z));
        InterpolateOverLineAndExcess<T, Descriptor1>(
            createBottomDomain(coarseDomain, direction), coarseLattice, fineLattice,
            interpolationEngine->clone(), rescaleEngine->clone(), deltas);
    }

    if (!topOutOfRange) {
        std::vector<Dot3D> deltas;
        deltas.push_back(Dot3D(delta[0].x, delta[0].y, delta[0].z));
        deltas.push_back(Dot3D(delta[1].x, delta[1].y, delta[1].z));
        InterpolateOverLineAndExcess<T, Descriptor1>(
            createTopDomain(coarseDomain, direction), coarseLattice, fineLattice,
            interpolationEngine->clone(), rescaleEngine->clone(), deltas);
    }

    if (!rightOutOfRange) {
        std::vector<Dot3D> deltas;
        deltas.push_back(Dot3D(delta[1].x, delta[1].y, delta[1].z));
        deltas.push_back(Dot3D(delta[0].x, delta[0].y, delta[0].z));
        InterpolateOverLineAndExcess<T, Descriptor1>(
            createRightDomain(coarseDomain, direction), coarseLattice, fineLattice,
            interpolationEngine->clone(), rescaleEngine->clone(), deltas);
    }

    // STEP 3: We perform interpolations over the excedent corners if needed
    std::vector<Box3D> corners = createCorners(coarseDomain, direction);

    if (!leftOutOfRange && !topOutOfRange) {
        const Dot3D deltas[2] = {
            Dot3D(delta[1].x, delta[1].y, delta[1].z),
            Dot3D(-delta[0].x, -delta[0].y, -delta[0].z)};
        applyProcessingFunctional(
            new InterpolateCoarseToFineCornerDynamics3D<T, Descriptor1, Descriptor2>(
                rescaleEngine->clone(), interpolationEngine->clone(), deltas),
            corners[0], coarseLattice, fineLattice);
    }

    if (!leftOutOfRange && !bottomOutOfRange) {
        const Dot3D deltas[2] = {
            Dot3D(-delta[1].x, -delta[1].y, -delta[1].z),
            Dot3D(-delta[0].x, -delta[0].y, -delta[0].z)};
        applyProcessingFunctional(
            new InterpolateCoarseToFineCornerDynamics3D<T, Descriptor1, Descriptor2>(
                rescaleEngine->clone(), interpolationEngine->clone(), deltas),
            corners[1], coarseLattice, fineLattice);
    }
    if (!rightOutOfRange && !bottomOutOfRange) {
        const Dot3D deltas[2] = {
            Dot3D(-delta[1].x, -delta[1].y, -delta[1].z),
            Dot3D(delta[0].x, delta[0].y, delta[0].z)};
        applyProcessingFunctional(
            new InterpolateCoarseToFineCornerDynamics3D<T, Descriptor1, Descriptor2>(
                rescaleEngine->clone(), interpolationEngine->clone(), deltas),
            corners[2], coarseLattice, fineLattice);
    }
    if (!rightOutOfRange && !topOutOfRange) {
        const Dot3D deltas[2] = {
            Dot3D(delta[1].x, delta[1].y, delta[1].z), Dot3D(delta[0].x, delta[0].y, delta[0].z)};
        applyProcessingFunctional(
            new InterpolateCoarseToFineCornerDynamics3D<T, Descriptor1, Descriptor2>(
                rescaleEngine->clone(), interpolationEngine->clone(), deltas),
            corners[3], coarseLattice, fineLattice);
    }
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
InterpolateCoarseToFineDynamics3D<T, Descriptor1, Descriptor2>
    *InterpolateCoarseToFineDynamics3D<T, Descriptor1, Descriptor2>::clone() const
{
    return new InterpolateCoarseToFineDynamics3D(*this);
}

/* *************** Class InterpolateCoarseToFineBoundaryDynamics3D ********** */

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
InterpolateCoarseToFineBoundaryDynamics3D<T, Descriptor1, Descriptor2>::
    InterpolateCoarseToFineBoundaryDynamics3D(
        RescaleEngine<T, Descriptor1> *rescaleEngine_,
        InterpolationEngine2D<T, Descriptor1> *interpolationEngine_, const Dot3D delta_[2]) :
    rescaleEngine(rescaleEngine_), interpolationEngine(interpolationEngine_)

{
    delta.push_back(delta_[0]);
    delta.push_back(delta_[1]);
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
InterpolateCoarseToFineBoundaryDynamics3D<T, Descriptor1, Descriptor2>::
    InterpolateCoarseToFineBoundaryDynamics3D(
        InterpolateCoarseToFineBoundaryDynamics3D<T, Descriptor1, Descriptor2> const &rhs) :
    rescaleEngine(rhs.rescaleEngine->clone()),
    interpolationEngine(rhs.interpolationEngine->clone()),
    delta(rhs.delta)

{ }

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
InterpolateCoarseToFineBoundaryDynamics3D<T, Descriptor1, Descriptor2>
    &InterpolateCoarseToFineBoundaryDynamics3D<T, Descriptor1, Descriptor2>::operator=(
        InterpolateCoarseToFineBoundaryDynamics3D<T, Descriptor1, Descriptor2> const &rhs)
{
    delete rescaleEngine;
    rescaleEngine = rhs.rescaleEngine->clone();
    delete interpolationEngine;
    interpolationEngine = rhs.interpolationEngine->clone();
    delta.resize(0);
    delta.push_back(rhs.delta[0]);
    delta.push_back(rhs.delta[1]);
    return *this;
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
InterpolateCoarseToFineBoundaryDynamics3D<
    T, Descriptor1, Descriptor2>::~InterpolateCoarseToFineBoundaryDynamics3D()
{
    delete rescaleEngine;
    delete interpolationEngine;
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
void InterpolateCoarseToFineBoundaryDynamics3D<T, Descriptor1, Descriptor2>::process(
    Box3D domain, BlockLattice3D<T, Descriptor1> &coarseLattice,
    BlockLattice3D<T, Descriptor2> &fineLattice)
{
    InterpolateOverLineAndExcess<T, Descriptor1>(
        domain, coarseLattice, fineLattice, interpolationEngine->clone(), rescaleEngine->clone(),
        delta);
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
InterpolateCoarseToFineBoundaryDynamics3D<T, Descriptor1, Descriptor2>
    *InterpolateCoarseToFineBoundaryDynamics3D<T, Descriptor1, Descriptor2>::clone() const
{
    return new InterpolateCoarseToFineBoundaryDynamics3D(*this);
}

/* *************** Class InterpolateCoarseToFineCornerDynamics3D ********** */

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
InterpolateCoarseToFineCornerDynamics3D<T, Descriptor1, Descriptor2>::
    InterpolateCoarseToFineCornerDynamics3D(
        RescaleEngine<T, Descriptor1> *rescaleEngine_,
        InterpolationEngine2D<T, Descriptor1> *interpolationEngine_, const Dot3D delta_[2]) :
    rescaleEngine(rescaleEngine_), interpolationEngine(interpolationEngine_)
{
    delta.push_back(delta_[0]);
    delta.push_back(delta_[1]);
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
InterpolateCoarseToFineCornerDynamics3D<T, Descriptor1, Descriptor2>::
    InterpolateCoarseToFineCornerDynamics3D(
        InterpolateCoarseToFineCornerDynamics3D<T, Descriptor1, Descriptor2> const &rhs) :
    rescaleEngine(rhs.rescaleEngine->clone()),
    interpolationEngine(rhs.interpolationEngine->clone()),
    delta(rhs.delta)
{ }

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
InterpolateCoarseToFineCornerDynamics3D<T, Descriptor1, Descriptor2>
    &InterpolateCoarseToFineCornerDynamics3D<T, Descriptor1, Descriptor2>::operator=(
        InterpolateCoarseToFineCornerDynamics3D<T, Descriptor1, Descriptor2> const &rhs)
{
    delete rescaleEngine;
    rescaleEngine = rhs.rescaleEngine->clone();
    delete interpolationEngine;
    interpolationEngine = rhs.interpolationEngine->clone();
    delta.resize(0);
    delta.push_back(rhs.delta[0]);
    delta.push_back(rhs.delta[1]);
    return *this;
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
InterpolateCoarseToFineCornerDynamics3D<
    T, Descriptor1, Descriptor2>::~InterpolateCoarseToFineCornerDynamics3D()
{
    delete rescaleEngine;
    delete interpolationEngine;
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
void InterpolateCoarseToFineCornerDynamics3D<T, Descriptor1, Descriptor2>::process(
    Box3D domain, BlockLattice3D<T, Descriptor1> &coarseLattice,
    BlockLattice3D<T, Descriptor2> &fineLattice)
{
    PLB_PRECONDITION(domain.x0 == domain.x1);
    PLB_PRECONDITION(domain.y0 == domain.y1);
    PLB_PRECONDITION(domain.z0 == domain.z1);

    plint x0 = domain.x0;
    plint y0 = domain.y0;
    plint z0 = domain.z0;

    plint whichTime = 1;
    CoarseToFineConverter3D<T, Descriptor1> converter(coarseLattice, fineLattice);

    // Check if the extended domain (over which the interpolation is executed) exceeds either
    //   the coarse or the fine lattice, in which case interpolation is not executed.
    Dot3D coarseExcess1(x0 + delta[0].x, y0 + delta[0].y, z0 + delta[0].z);
    Dot3D fineExcess1(
        converter.fineX(x0) + delta[0].x, converter.fineY(y0) + delta[0].y,
        converter.fineZ(z0) + delta[0].z);
    Dot3D coarseExcess2(x0 + delta[1].x, y0 + delta[1].y, z0 + delta[1].z);
    Dot3D fineExcess2(
        converter.fineX(x0) + delta[1].x, converter.fineY(y0) + delta[1].y,
        converter.fineZ(z0) + delta[1].z);

    bool outOfRange1 =
        (!contained(coarseExcess1, coarseLattice.getBoundingBox())
         || !contained(fineExcess1, fineLattice.getBoundingBox()));
    bool outOfRange2 =
        (!contained(coarseExcess2, coarseLattice.getBoundingBox())
         || !contained(fineExcess2, fineLattice.getBoundingBox()));

    // Copy coarse->left on the chosen cell, and execute interpolation, unless the
    //   interpolation is out-of-range.
    if (!outOfRange1) {
        interpolationEngine->pushToCouple(
            Array<Cell<T, Descriptor1>, 2>(
                coarseLattice.get(x0, y0, z0),
                coarseLattice.get(x0 + delta[0].x, y0 + delta[0].y, z0 + delta[0].z)),
            *rescaleEngine);
    }

    if (outOfRange2) {
        interpolationEngine->pushToCouple();
    } else {
        interpolationEngine->pushToCouple(
            Array<Cell<T, Descriptor1>, 2>(
                coarseLattice.get(x0 + delta[1].x, y0 + delta[1].y, z0 + delta[1].z),
                coarseLattice.get(
                    x0 + delta[0].x + delta[1].x, y0 + delta[0].y + delta[1].y,
                    z0 + delta[0].z + delta[1].z)),
            *rescaleEngine);
    }

    interpolationEngine->copyFrom(
        converter.fineDynamics(x0, y0, z0, 0, 0, 0).getDecomposedValues(whichTime), 0);

    if (!outOfRange1) {
        interpolationEngine->interpolateFromCouple(
            converter.fineDynamics(x0, y0, z0, delta[0].x, delta[0].y, delta[0].z)
                .getDecomposedValues(whichTime),
            0, 1);
    }

    if (!outOfRange2) {
        interpolationEngine->interpolateFromCouple(
            converter.fineDynamics(x0, y0, z0, delta[1].x, delta[1].y, delta[1].z)
                .getDecomposedValues(whichTime),
            0, 2);
    }

    if (!outOfRange1 && !outOfRange2) {
        interpolationEngine->interpolateFromAll(converter
                                                    .fineDynamics(
                                                        x0, y0, z0, delta[1].x + delta[0].x,
                                                        delta[1].y + delta[0].y,
                                                        delta[1].z + delta[0].z)
                                                    .getDecomposedValues(whichTime));
    }
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
InterpolateCoarseToFineCornerDynamics3D<T, Descriptor1, Descriptor2>
    *InterpolateCoarseToFineCornerDynamics3D<T, Descriptor1, Descriptor2>::clone() const
{
    return new InterpolateCoarseToFineCornerDynamics3D(*this);
}

/* *************** Class Copy_t1_to_t0_3D *********************************** */

template <typename T, template <typename U> class Descriptor>
Copy_t1_to_t0_3D<T, Descriptor>::Copy_t1_to_t0_3D(plint numTimeSteps_, plint executionTime_) :
    numTimeSteps(numTimeSteps_), executionTime(executionTime_)
{ }

template <typename T, template <typename U> class Descriptor>
void Copy_t1_to_t0_3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &fineLattice)
{
    // Determine current relative value of time steps.
    size_t relativeTime = fineLattice.getTimeCounter().getTime() % numTimeSteps;
    // Execute data processor only if one is at the end of a cycle.
    if ((plint)relativeTime == executionTime) {
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    FineGridBoundaryDynamics<T, Descriptor> &fineDynamics =
                        dynamic_cast<FineGridBoundaryDynamics<T, Descriptor> &>(
                            fineLattice.get(iX, iY, iZ).getDynamics());
                    std::vector<T> &t0 = fineDynamics.getDecomposedValues(0);
                    std::vector<T> &t1 = fineDynamics.getDecomposedValues(1);
                    t0.assign(t1.begin(), t1.end());
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
Copy_t1_to_t0_3D<T, Descriptor> *Copy_t1_to_t0_3D<T, Descriptor>::clone() const
{
    return new Copy_t1_to_t0_3D(*this);
}

/* *************** Class CopyFineToCoarse3D ********************************* */

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyFineToCoarse3D<T, Descriptor1, Descriptor2>::CopyFineToCoarse3D(
    RescaleEngine<T, Descriptor1> *rescaleEngine_, plint numTimeSteps_, plint executionTime_) :
    rescaleEngine(rescaleEngine_), numTimeSteps(numTimeSteps_), executionTime(executionTime_)
{ }

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyFineToCoarse3D<T, Descriptor1, Descriptor2>::~CopyFineToCoarse3D()
{
    delete rescaleEngine;
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyFineToCoarse3D<T, Descriptor1, Descriptor2>::CopyFineToCoarse3D(
    CopyFineToCoarse3D<T, Descriptor1, Descriptor2> const &rhs) :
    rescaleEngine(rhs.rescaleEngine->clone()),
    numTimeSteps(rhs.numTimeSteps),
    executionTime(rhs.executionTime)
{ }

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyFineToCoarse3D<T, Descriptor1, Descriptor2>
    &CopyFineToCoarse3D<T, Descriptor1, Descriptor2>::operator=(
        CopyFineToCoarse3D<T, Descriptor1, Descriptor2> const &rhs)
{
    delete rescaleEngine;
    rescaleEngine = rhs.rescaleEngine->clone();
    numTimeSteps = rhs.numTimeSteps;
    executionTime = rhs.executionTime;
    return *this;
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
void CopyFineToCoarse3D<T, Descriptor1, Descriptor2>::process(
    Box3D fineDomain, BlockLattice3D<T, Descriptor1> &fineLattice,
    BlockLattice3D<T, Descriptor2> &coarseLattice)
{
    PLB_PRECONDITION(
        fineDomain.x0 == fineDomain.x1 || fineDomain.y0 == fineDomain.y1
        || fineDomain.z0 == fineDomain.z1);

    // Determine current relative value of time steps
    size_t relativeTime = fineLattice.getTimeCounter().getTime() % numTimeSteps;
    // Execute data processor only if one is at the end of a cycle (iT=1)
    if ((plint)relativeTime == executionTime) {
        Dot3D posFine = fineLattice.getLocation();      // Position of fine grid in multi block.
        Dot3D posCoarse = coarseLattice.getLocation();  // Position of coarse grid in multi block.
        Box3D coarseDomain(fineDomain.shift(posFine.x, posFine.y, posFine.z)
                               .  // Convert to absolute fine coordinates.
                           divideAndFitSmaller(2)
                               .  // Rescale, but don't exceed original domain.
                           shift(
                               -posCoarse.x, -posCoarse.y,
                               -posCoarse.z));  // Convert to relative coarse coordinates.
        std::vector<T> decomposedCoarseValues;
        // Loop over coarse lattice
        for (plint iX = coarseDomain.x0; iX <= coarseDomain.x1; ++iX) {
            for (plint iY = coarseDomain.y0; iY <= coarseDomain.y1; ++iY) {
                for (plint iZ = coarseDomain.z0; iZ <= coarseDomain.z1; ++iZ) {
                    // Determine corresonding coordinates on fine lattice
                    plint fineX = (iX + posCoarse.x) * 2 - posFine.x;
                    plint fineY = (iY + posCoarse.y) * 2 - posFine.y;
                    plint fineZ = (iZ + posCoarse.z) * 2 - posFine.z;

                    Cell<T, Descriptor1> &coarseCell = coarseLattice.get(iX, iY, iZ);
                    Cell<T, Descriptor2> const &fineCell = fineLattice.get(fineX, fineY, fineZ);

                    rescaleEngine->scaleFineCoarse(fineCell, decomposedCoarseValues);
                    rescaleEngine->recompose(coarseCell, decomposedCoarseValues);
                }
            }
        }
    }
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyFineToCoarse3D<T, Descriptor1, Descriptor2>
    *CopyFineToCoarse3D<T, Descriptor1, Descriptor2>::clone() const
{
    return new CopyFineToCoarse3D(*this);
}

}  // namespace plb

#endif  // GRID_REFINEMENT_PROCESSORS_3D_HH
