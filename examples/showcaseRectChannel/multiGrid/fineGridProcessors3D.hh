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
 * Dynamics and data processors used to implement 3D grid refinement -- implementation.
 */

#ifndef FINE_GRID_PROCESSORS_3D_HH
#define FINE_GRID_PROCESSORS_3D_HH

#include "multiGrid/fineGridProcessors3D.h"

namespace plb {

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

                    // TODO: erase this and create a data processor that executes this only when
                    // needed
                    fineDynamics.completePopulations(fineLattice.get(iX, iY, iZ));
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

/* ***************  ScalarCubicInterpolationYZ  *************** */

template <typename T, template <typename U> class Descriptor>
void ScalarCubicInterpolationYZ<T, Descriptor>::process(
    Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
    BlockLattice3D<T, Descriptor> &fineLattice)
{
    Dot3D posCoarse = coarseLattice.getLocation();
    Dot3D posFine = fineLattice.getLocation();

    const plint x0 = coarseDomain.x0;
    const plint y0 = coarseDomain.y0;
    const plint z0 = coarseDomain.z0;

    const plint y1 = coarseDomain.y1;
    const plint z1 = coarseDomain.z1;
    PLB_PRECONDITION(x0 == coarseDomain.x1);

    plint iX = x0;
    Array<Array<T, 2>, 5> centeredPoints;
    centeredPoints[0] = Array<T, 2>(1.0, 1.5);
    centeredPoints[1] = Array<T, 2>(1.5, 1.0);
    centeredPoints[2] = Array<T, 2>(1.5, 1.5);
    centeredPoints[3] = Array<T, 2>(1.5, 2.0);
    centeredPoints[4] = Array<T, 2>(2.0, 1.5);

    T neighbors[4][4];
    std::vector<std::vector<T> > pop(16);

    for (plint iY = y0 - 1; iY <= y1; ++iY) {
        for (plint iZ = z0 - 1; iZ <= z1; ++iZ) {
            // extracting and rescaling the known 16 coarse values
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY - 1, iZ - 1), pop[0]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ - 1), pop[1]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 1, iZ - 1), pop[2]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 2, iZ - 1), pop[3]);

            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY - 1, iZ), pop[4]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop[5]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 1, iZ), pop[6]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 2, iZ), pop[7]);

            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY - 1, iZ + 1), pop[8]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + 1), pop[9]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 1, iZ + 1), pop[10]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 2, iZ + 1), pop[11]);

            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY - 1, iZ + 2), pop[12]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + 2), pop[13]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 1, iZ + 2), pop[14]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 2, iZ + 2), pop[15]);

            plint cellDim = pop[0].size();  // TODO get the size of the cells more generally
            // the containers for the interpolated values
            std::vector<T> interpValues1(cellDim);
            std::vector<T> interpValues2(cellDim);
            std::vector<T> interpValues3(cellDim);
            std::vector<T> interpValues4(cellDim);
            std::vector<T> interpValues5(cellDim);

            for (plint iComp = 0; iComp < cellDim; ++iComp) {
                // get the 16 neighbors and interpolate where it is needed
                neighbors[0][0] = pop[0][iComp];
                neighbors[1][0] = pop[1][iComp];
                neighbors[2][0] = pop[2][iComp];
                neighbors[3][0] = pop[3][iComp];

                neighbors[0][1] = pop[4][iComp];
                neighbors[1][1] = pop[5][iComp];
                neighbors[2][1] = pop[6][iComp];
                neighbors[3][1] = pop[7][iComp];

                neighbors[0][2] = pop[8][iComp];
                neighbors[1][2] = pop[9][iComp];
                neighbors[2][2] = pop[10][iComp];
                neighbors[3][2] = pop[11][iComp];

                neighbors[0][3] = pop[12][iComp];
                neighbors[1][3] = pop[13][iComp];
                neighbors[2][3] = pop[14][iComp];
                neighbors[3][3] = pop[15][iComp];

                // interpolate the values according to the neighboring values, use centered schema
                std::vector<T> interpolatedValues = symetricCubicInterpolation<T>(neighbors);

                interpValues1[iComp] = interpolatedValues[0];
                interpValues2[iComp] = interpolatedValues[1];
                interpValues3[iComp] = interpolatedValues[2];
                interpValues4[iComp] = interpolatedValues[3];
                interpValues5[iComp] = interpolatedValues[4];
            }

            // convert the coarse coordinates to fine coordinates
            plint fineX = (iX + posCoarse.x) * 2 - posFine.x;
            plint fineY = (iY + posCoarse.y) * 2 - posFine.y;
            plint fineZ = (iZ + posCoarse.z) * 2 - posFine.z;

            // assigning the known values
            copyPopulations(pop[5], fineLattice.get(fineX, fineY, fineZ));
            copyPopulations(pop[6], fineLattice.get(fineX, fineY + 2, fineZ));
            copyPopulations(pop[9], fineLattice.get(fineX, fineY, fineZ + 2));
            copyPopulations(pop[10], fineLattice.get(fineX, fineY + 2, fineZ + 2));

            // assigning the interpolated values
            copyPopulations(interpValues1, fineLattice.get(fineX, fineY, fineZ + 1));
            copyPopulations(interpValues2, fineLattice.get(fineX, fineY + 1, fineZ));
            copyPopulations(interpValues3, fineLattice.get(fineX, fineY + 1, fineZ + 1));
            copyPopulations(interpValues4, fineLattice.get(fineX, fineY + 1, fineZ + 2));
            copyPopulations(interpValues5, fineLattice.get(fineX, fineY + 2, fineZ + 1));
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ScalarCubicInterpolationYZ<T, Descriptor> *ScalarCubicInterpolationYZ<T, Descriptor>::clone() const
{
    return new ScalarCubicInterpolationYZ<T, Descriptor>(rescaleEngine->clone());
}

/* ******************* ScalarCubicInterpolationYZLineY3D ******************* */

template <typename T, template <typename U> class Descriptor>
ScalarCubicInterpolationYZLineY3D<T, Descriptor>::ScalarCubicInterpolationYZLineY3D(
    plint delta_, RescaleEngine<T, Descriptor> *rescaleEngine_) :
    delta(delta_), rescaleEngine(rescaleEngine_)
{
    PLB_PRECONDITION(delta == 1 || delta == -1);
}

template <typename T, template <typename U> class Descriptor>
void ScalarCubicInterpolationYZLineY3D<T, Descriptor>::process(
    Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
    BlockLattice3D<T, Descriptor> &fineLattice)
{
    Dot3D posCoarse = coarseLattice.getLocation();
    Dot3D posFine = fineLattice.getLocation();

    const plint x0 = coarseDomain.x0;
    const plint y0 = coarseDomain.y0;
    const plint z0 = coarseDomain.z0;

    const plint y1 = coarseDomain.y1;

    PLB_PRECONDITION(x0 == coarseDomain.x1);
    PLB_PRECONDITION(z0 == coarseDomain.z1);

    plint iX = x0;
    plint iZ = z0;

    T neighbors[4][3];
    std::vector<std::vector<T> > pop(12);

    for (plint iY = y0 - 1; iY <= y1; ++iY) {
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY - 1, iZ), pop[0]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop[1]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 1, iZ), pop[2]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 2, iZ), pop[3]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY - 1, iZ + delta), pop[4]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + delta), pop[5]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 1, iZ + delta), pop[6]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 2, iZ + delta), pop[7]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY - 1, iZ + 2 * delta), pop[8]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + 2 * delta), pop[9]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 1, iZ + 2 * delta), pop[10]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 2, iZ + 2 * delta), pop[11]);

        plint cellDim = pop[0].size();  // TODO get the size of the cells more generally

        std::vector<T> interpValues1(cellDim);
        std::vector<T> interpValues2(cellDim);
        std::vector<T> interpValues3(cellDim);
        std::vector<T> interpValues4(cellDim);

        for (plint iComp = 0; iComp < cellDim; ++iComp) {
            // get the 12 neighbors and interpolate where it is needed
            neighbors[0][0] = pop[0][iComp];
            neighbors[1][0] = pop[1][iComp];
            neighbors[2][0] = pop[2][iComp];
            neighbors[3][0] = pop[3][iComp];

            neighbors[0][1] = pop[4][iComp];
            neighbors[1][1] = pop[5][iComp];
            neighbors[2][1] = pop[6][iComp];
            neighbors[3][1] = pop[7][iComp];

            neighbors[0][2] = pop[8][iComp];
            neighbors[1][2] = pop[9][iComp];
            neighbors[2][2] = pop[10][iComp];
            neighbors[3][2] = pop[11][iComp];

            // interpolate the values according to the neighboring values, use centered schema
            std::vector<T> interpolatedValues = asymetricCubicInterpolation<T>(neighbors);

            interpValues1[iComp] = interpolatedValues[0];
            interpValues2[iComp] = interpolatedValues[1];
            interpValues3[iComp] = interpolatedValues[2];
            interpValues4[iComp] = interpolatedValues[3];
        }

        // convert the coarse coordinates to fine coordinates
        plint fineX = (iX + posCoarse.x) * 2 - posFine.x;
        plint fineY = (iY + posCoarse.y) * 2 - posFine.y;
        plint fineZ = (iZ + posCoarse.z) * 2 - posFine.z;

        // copy the 4 known values
        copyPopulations(pop[1], fineLattice.get(fineX, fineY, fineZ));
        copyPopulations(pop[2], fineLattice.get(fineX, fineY + 2, fineZ));
        copyPopulations(pop[5], fineLattice.get(fineX, fineY, fineZ + 2 * delta));
        copyPopulations(pop[6], fineLattice.get(fineX, fineY + 2, fineZ + 2 * delta));

        // assigning the computed interpolated values
        copyPopulations(interpValues1, fineLattice.get(fineX, fineY, fineZ + delta));
        copyPopulations(interpValues2, fineLattice.get(fineX, fineY + 1, fineZ));
        copyPopulations(interpValues3, fineLattice.get(fineX, fineY + 1, fineZ + delta));
        copyPopulations(interpValues4, fineLattice.get(fineX, fineY + 2, fineZ + delta));
    }
}

template <typename T, template <typename U> class Descriptor>
ScalarCubicInterpolationYZLineY3D<T, Descriptor>
    *ScalarCubicInterpolationYZLineY3D<T, Descriptor>::clone() const
{
    return new ScalarCubicInterpolationYZLineY3D<T, Descriptor>(delta, rescaleEngine->clone());
}

/* ******************* ScalarCubicInterpolationYZLineZ3D ******************* */

template <typename T, template <typename U> class Descriptor>
ScalarCubicInterpolationYZLineZ3D<T, Descriptor>::ScalarCubicInterpolationYZLineZ3D(
    plint delta_, RescaleEngine<T, Descriptor> *rescaleEngine_) :
    delta(delta_), rescaleEngine(rescaleEngine_)
{ }

template <typename T, template <typename U> class Descriptor>
void ScalarCubicInterpolationYZLineZ3D<T, Descriptor>::process(
    Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
    BlockLattice3D<T, Descriptor> &fineLattice)
{
    Dot3D posCoarse = coarseLattice.getLocation();
    Dot3D posFine = fineLattice.getLocation();

    const plint x0 = coarseDomain.x0;
    const plint y0 = coarseDomain.y0;
    const plint z0 = coarseDomain.z0;

    const plint z1 = coarseDomain.z1;

    PLB_PRECONDITION(x0 == coarseDomain.x1);
    PLB_PRECONDITION(y0 == coarseDomain.y1);

    plint iX = x0;
    plint iY = y0;

    T neighbors[4][3];
    std::vector<std::vector<T> > pop(12);

    for (plint iZ = z0 - 1; iZ <= z1; ++iZ) {
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ - 1), pop[0]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop[1]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + 1), pop[2]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + 2), pop[3]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + delta, iZ - 1), pop[4]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + delta, iZ), pop[5]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + delta, iZ + 1), pop[6]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + delta, iZ + 2), pop[7]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 2 * delta, iZ - 1), pop[8]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 2 * delta, iZ), pop[9]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 2 * delta, iZ + 1), pop[10]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 2 * delta, iZ + 2), pop[11]);

        plint cellDim = pop[0].size();  // TODO get the size of the cells more generally

        std::vector<T> interpValues1(cellDim);
        std::vector<T> interpValues2(cellDim);
        std::vector<T> interpValues3(cellDim);
        std::vector<T> interpValues4(cellDim);

        for (plint iComp = 0; iComp < cellDim; ++iComp) {
            // get the 12 neighbors and interpolate where it is needed
            neighbors[0][0] = pop[0][iComp];
            neighbors[1][0] = pop[1][iComp];
            neighbors[2][0] = pop[2][iComp];
            neighbors[3][0] = pop[3][iComp];

            neighbors[0][1] = pop[4][iComp];
            neighbors[1][1] = pop[5][iComp];
            neighbors[2][1] = pop[6][iComp];
            neighbors[3][1] = pop[7][iComp];

            neighbors[0][2] = pop[8][iComp];
            neighbors[1][2] = pop[9][iComp];
            neighbors[2][2] = pop[10][iComp];
            neighbors[3][2] = pop[11][iComp];

            // interpolate the values according to the neighboring values, use centered schema
            std::vector<T> interpolatedValues = asymetricCubicInterpolation<T>(neighbors);

            interpValues1[iComp] = interpolatedValues[0];
            interpValues2[iComp] = interpolatedValues[1];
            interpValues3[iComp] = interpolatedValues[2];
            interpValues4[iComp] = interpolatedValues[3];
        }

        // convert the coarse coordinates to fine coordinates
        plint fineX = (iX + posCoarse.x) * 2 - posFine.x;
        plint fineY = (iY + posCoarse.y) * 2 - posFine.y;
        plint fineZ = (iZ + posCoarse.z) * 2 - posFine.z;

        // copy the 4 known values
        copyPopulations(pop[1], fineLattice.get(fineX, fineY, fineZ));
        copyPopulations(pop[5], fineLattice.get(fineX, fineY + 2 * delta, fineZ));
        copyPopulations(pop[2], fineLattice.get(fineX, fineY, fineZ + 2));
        copyPopulations(pop[6], fineLattice.get(fineX, fineY + 2 * delta, fineZ + 2));

        // assigning the computed interpolated values
        copyPopulations(interpValues1, fineLattice.get(fineX, fineY + delta, fineZ));
        copyPopulations(interpValues2, fineLattice.get(fineX, fineY, fineZ + 1));
        copyPopulations(interpValues3, fineLattice.get(fineX, fineY + delta, fineZ + 1));
        copyPopulations(interpValues4, fineLattice.get(fineX, fineY + delta, fineZ + 2));
    }
}

template <typename T, template <typename U> class Descriptor>
ScalarCubicInterpolationYZLineZ3D<T, Descriptor>
    *ScalarCubicInterpolationYZLineZ3D<T, Descriptor>::clone() const
{
    return new ScalarCubicInterpolationYZLineZ3D<T, Descriptor>(delta, rescaleEngine->clone());
}

/* ******************* ScalarCubicInterpolationYZCorner3D ******************* */
template <typename T, template <typename U> class Descriptor>
ScalarCubicInterpolationYZCorner3D<T, Descriptor>::ScalarCubicInterpolationYZCorner3D(
    plint deltaY_, plint deltaZ_, RescaleEngine<T, Descriptor> *rescaleEngine_) :
    deltaY(deltaY_), deltaZ(deltaZ_), rescaleEngine(rescaleEngine_)
{ }

template <typename T, template <typename U> class Descriptor>
void ScalarCubicInterpolationYZCorner3D<T, Descriptor>::process(
    Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
    BlockLattice3D<T, Descriptor> &fineLattice)
{
    Dot3D posCoarse = coarseLattice.getLocation();
    Dot3D posFine = fineLattice.getLocation();

    PLB_PRECONDITION(coarseDomain.x0 == coarseDomain.x1);
    PLB_PRECONDITION(coarseDomain.y0 == coarseDomain.y1);
    PLB_PRECONDITION(coarseDomain.z0 == coarseDomain.z1);

    const plint iX = coarseDomain.x0;
    const plint iY = coarseDomain.y0;
    const plint iZ = coarseDomain.z0;

    plint fineX = (iX + posCoarse.x) * 2 - posFine.x;
    plint fineY = (iY + posCoarse.y) * 2 - posFine.y;
    plint fineZ = (iZ + posCoarse.z) * 2 - posFine.z;

    T neighbors[3][3];
    std::vector<std::vector<T> > pop(9);

    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop[0]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + deltaY, iZ), pop[1]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 2 * deltaY, iZ), pop[2]);

    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + deltaZ), pop[3]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + deltaY, iZ + deltaZ), pop[4]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 2 * deltaY, iZ + deltaZ), pop[5]);

    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + 2 * deltaZ), pop[6]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + deltaY, iZ + 2 * deltaZ), pop[7]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 2 * deltaY, iZ + 2 * deltaZ), pop[8]);

    plint cellDim = pop[0].size();  // TODO get the size of the cells more generally

    std::vector<T> interpValues1(cellDim);
    std::vector<T> interpValues2(cellDim);
    std::vector<T> interpValues3(cellDim);

    // for every component in the decomposed values
    for (plint iComp = 0; iComp < cellDim; ++iComp) {
        // get the 9 neighbors and interpolate where it is needed
        neighbors[0][0] = pop[0][iComp];
        neighbors[1][0] = pop[1][iComp];
        neighbors[2][0] = pop[2][iComp];

        neighbors[0][1] = pop[3][iComp];
        neighbors[1][1] = pop[4][iComp];
        neighbors[2][1] = pop[5][iComp];

        neighbors[0][2] = pop[6][iComp];
        neighbors[1][2] = pop[7][iComp];
        neighbors[2][2] = pop[8][iComp];

        // interpolate the values according to the neighboring values, use centered schema
        std::vector<T> interpolatedValues = cornerInterpolation<T>(neighbors);

        // pack the values in vectors to use latter
        interpValues1[iComp] = interpolatedValues[0];
        interpValues2[iComp] = interpolatedValues[1];
        interpValues3[iComp] = interpolatedValues[2];
    }

    // copy the 4 known values
    copyPopulations(pop[0], fineLattice.get(fineX, fineY, fineZ));
    copyPopulations(pop[1], fineLattice.get(fineX, fineY + 2 * deltaY, fineZ));
    copyPopulations(pop[3], fineLattice.get(fineX, fineY, fineZ + 2 * deltaZ));
    copyPopulations(pop[4], fineLattice.get(fineX, fineY + 2 * deltaY, fineZ + 2 * deltaZ));

    // assigning the interpolated values
    copyPopulations(interpValues1, fineLattice.get(fineX, fineY, fineZ + deltaZ));
    copyPopulations(interpValues2, fineLattice.get(fineX, fineY + deltaY, fineZ));
    copyPopulations(interpValues3, fineLattice.get(fineX, fineY + deltaY, fineZ + deltaZ));
}

template <typename T, template <typename U> class Descriptor>
ScalarCubicInterpolationYZCorner3D<T, Descriptor>
    *ScalarCubicInterpolationYZCorner3D<T, Descriptor>::clone() const
{
    return new ScalarCubicInterpolationYZCorner3D<T, Descriptor>(
        deltaY, deltaZ, rescaleEngine->clone());
}

/* ************* ScalarCubicInterpolationXY ******************* */
template <typename T, template <typename U> class Descriptor>
void ScalarCubicInterpolationXY<T, Descriptor>::process(
    Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
    BlockLattice3D<T, Descriptor> &fineLattice)
{
    Dot3D posCoarse = coarseLattice.getLocation();
    Dot3D posFine = fineLattice.getLocation();

    const plint x0 = coarseDomain.x0;
    const plint y0 = coarseDomain.y0;
    const plint z0 = coarseDomain.z0;

    const plint x1 = coarseDomain.x1;
    const plint y1 = coarseDomain.y1;

    PLB_PRECONDITION(z0 == coarseDomain.z1);

    plint iZ = z0;
    Array<Array<T, 2>, 5> centeredPoints;
    centeredPoints[0] = Array<T, 2>(1.0, 1.5);
    centeredPoints[1] = Array<T, 2>(1.5, 1.0);
    centeredPoints[2] = Array<T, 2>(1.5, 1.5);
    centeredPoints[3] = Array<T, 2>(1.5, 2.0);
    centeredPoints[4] = Array<T, 2>(2.0, 1.5);

    T neighbors[4][4];
    std::vector<std::vector<T> > pop(16);

    for (plint iX = x0 - 1; iX <= x1; ++iX) {
        for (plint iY = y0 - 1; iY <= y1; ++iY) {
            // extracting and rescaling the known 16 coarse values
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - 1, iY - 1, iZ), pop[0]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY - 1, iZ), pop[1]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 1, iY - 1, iZ), pop[2]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2, iY - 1, iZ), pop[3]);

            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - 1, iY, iZ), pop[4]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop[5]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 1, iY, iZ), pop[6]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2, iY, iZ), pop[7]);

            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - 1, iY + 1, iZ), pop[8]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 1, iZ), pop[9]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 1, iY + 1, iZ), pop[10]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2, iY + 1, iZ), pop[11]);

            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - 1, iY + 2, iZ), pop[12]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 2, iZ), pop[13]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 1, iY + 2, iZ), pop[14]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2, iY + 2, iZ), pop[15]);

            plint cellDim = pop[0].size();  // TODO get the size of the cells more generally
            // the containers for the interpolated values
            std::vector<T> interpValues1(cellDim);
            std::vector<T> interpValues2(cellDim);
            std::vector<T> interpValues3(cellDim);
            std::vector<T> interpValues4(cellDim);
            std::vector<T> interpValues5(cellDim);

            for (plint iComp = 0; iComp < cellDim; ++iComp) {
                // get the 16 neighbors and interpolate where it is needed
                neighbors[0][0] = pop[0][iComp];
                neighbors[1][0] = pop[1][iComp];
                neighbors[2][0] = pop[2][iComp];
                neighbors[3][0] = pop[3][iComp];

                neighbors[0][1] = pop[4][iComp];
                neighbors[1][1] = pop[5][iComp];
                neighbors[2][1] = pop[6][iComp];
                neighbors[3][1] = pop[7][iComp];

                neighbors[0][2] = pop[8][iComp];
                neighbors[1][2] = pop[9][iComp];
                neighbors[2][2] = pop[10][iComp];
                neighbors[3][2] = pop[11][iComp];

                neighbors[0][3] = pop[12][iComp];
                neighbors[1][3] = pop[13][iComp];
                neighbors[2][3] = pop[14][iComp];
                neighbors[3][3] = pop[15][iComp];

                // interpolate the values according to the neighboring values, use centered schema
                std::vector<T> interpolatedValues = symetricCubicInterpolation<T>(neighbors);

                interpValues1[iComp] = interpolatedValues[0];
                interpValues2[iComp] = interpolatedValues[1];
                interpValues3[iComp] = interpolatedValues[2];
                interpValues4[iComp] = interpolatedValues[3];
                interpValues5[iComp] = interpolatedValues[4];
            }

            // convert the coarse coordinates to fine coordinates
            plint fineX = (iX + posCoarse.x) * 2 - posFine.x;
            plint fineY = (iY + posCoarse.y) * 2 - posFine.y;
            plint fineZ = (iZ + posCoarse.z) * 2 - posFine.z;

            // assigning the known values
            copyPopulations(pop[5], fineLattice.get(fineX, fineY, fineZ));
            copyPopulations(pop[6], fineLattice.get(fineX + 2, fineY, fineZ));
            copyPopulations(pop[9], fineLattice.get(fineX, fineY + 2, fineZ));
            copyPopulations(pop[10], fineLattice.get(fineX + 2, fineY + 2, fineZ));

            // assigning the interpolated values
            copyPopulations(interpValues1, fineLattice.get(fineX, fineY + 1, fineZ));
            copyPopulations(interpValues2, fineLattice.get(fineX + 1, fineY, fineZ));
            copyPopulations(interpValues3, fineLattice.get(fineX + 1, fineY + 1, fineZ));
            copyPopulations(interpValues4, fineLattice.get(fineX + 1, fineY + 2, fineZ));
            copyPopulations(interpValues5, fineLattice.get(fineX + 2, fineY + 1, fineZ));
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ScalarCubicInterpolationXY<T, Descriptor> *ScalarCubicInterpolationXY<T, Descriptor>::clone() const
{
    return new ScalarCubicInterpolationXY<T, Descriptor>(rescaleEngine->clone());
}

/* **************** ScalarCubicInterpolationXYLineX3D ********************** */

template <typename T, template <typename U> class Descriptor>
ScalarCubicInterpolationXYLineX3D<T, Descriptor>::ScalarCubicInterpolationXYLineX3D(
    plint delta_, RescaleEngine<T, Descriptor> *rescaleEngine_) :
    delta(delta_), rescaleEngine(rescaleEngine_)
{
    PLB_PRECONDITION(delta == 1 || delta == -1);
}

template <typename T, template <typename U> class Descriptor>
void ScalarCubicInterpolationXYLineX3D<T, Descriptor>::process(
    Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
    BlockLattice3D<T, Descriptor> &fineLattice)
{
    Dot3D posCoarse = coarseLattice.getLocation();
    Dot3D posFine = fineLattice.getLocation();

    const plint x0 = coarseDomain.x0;
    const plint y0 = coarseDomain.y0;
    const plint z0 = coarseDomain.z0;

    const plint x1 = coarseDomain.x1;

    PLB_PRECONDITION(y0 == coarseDomain.y1);
    PLB_PRECONDITION(z0 == coarseDomain.z1);

    plint iY = y0;
    plint iZ = z0;

    T neighbors[4][3];
    std::vector<std::vector<T> > pop(12);

    for (plint iX = x0 - 1; iX <= x1; ++iX) {
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - 1, iY, iZ), pop[0]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop[1]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 1, iY, iZ), pop[2]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2, iY, iZ), pop[3]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - 1, iY + delta, iZ), pop[4]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + delta, iZ), pop[5]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 1, iY + delta, iZ), pop[6]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2, iY + delta, iZ), pop[7]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - 1, iY + 2 * delta, iZ), pop[8]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 2 * delta, iZ), pop[9]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 1, iY + 2 * delta, iZ), pop[10]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2, iY + 2 * delta, iZ), pop[11]);

        plint cellDim = pop[0].size();  // TODO get the size of the cells more generally

        std::vector<T> interpValues1(cellDim);
        std::vector<T> interpValues2(cellDim);
        std::vector<T> interpValues3(cellDim);
        std::vector<T> interpValues4(cellDim);

        for (plint iComp = 0; iComp < cellDim; ++iComp) {
            // get the 12 neighbors and interpolate where it is needed
            neighbors[0][0] = pop[0][iComp];
            neighbors[1][0] = pop[1][iComp];
            neighbors[2][0] = pop[2][iComp];
            neighbors[3][0] = pop[3][iComp];

            neighbors[0][1] = pop[4][iComp];
            neighbors[1][1] = pop[5][iComp];
            neighbors[2][1] = pop[6][iComp];
            neighbors[3][1] = pop[7][iComp];

            neighbors[0][2] = pop[8][iComp];
            neighbors[1][2] = pop[9][iComp];
            neighbors[2][2] = pop[10][iComp];
            neighbors[3][2] = pop[11][iComp];

            // interpolate the values according to the neighboring values, use centered schema
            std::vector<T> interpolatedValues = asymetricCubicInterpolation<T>(neighbors);

            interpValues1[iComp] = interpolatedValues[0];
            interpValues2[iComp] = interpolatedValues[1];
            interpValues3[iComp] = interpolatedValues[2];
            interpValues4[iComp] = interpolatedValues[3];
        }

        // convert the coarse coordinates to fine coordinates
        plint fineX = (iX + posCoarse.x) * 2 - posFine.x;
        plint fineY = (iY + posCoarse.y) * 2 - posFine.y;
        plint fineZ = (iZ + posCoarse.z) * 2 - posFine.z;

        // copy the 4 known values
        copyPopulations(pop[1], fineLattice.get(fineX, fineY, fineZ));
        copyPopulations(pop[2], fineLattice.get(fineX + 2, fineY, fineZ));
        copyPopulations(pop[5], fineLattice.get(fineX, fineY + 2 * delta, fineZ));
        copyPopulations(pop[6], fineLattice.get(fineX + 2, fineY + 2 * delta, fineZ));

        // assigning the computed interpolated values
        copyPopulations(interpValues1, fineLattice.get(fineX, fineY + delta, fineZ));
        copyPopulations(interpValues2, fineLattice.get(fineX + 1, fineY, fineZ));
        copyPopulations(interpValues3, fineLattice.get(fineX + 1, fineY + delta, fineZ));
        copyPopulations(interpValues4, fineLattice.get(fineX + 2, fineY + delta, fineZ));
    }
}

template <typename T, template <typename U> class Descriptor>
ScalarCubicInterpolationXYLineX3D<T, Descriptor>
    *ScalarCubicInterpolationXYLineX3D<T, Descriptor>::clone() const
{
    return new ScalarCubicInterpolationXYLineX3D<T, Descriptor>(delta, rescaleEngine->clone());
}

/* ************ ScalarCubicInterpolationXYLineY3D *************** */

template <typename T, template <typename U> class Descriptor>
ScalarCubicInterpolationXYLineY3D<T, Descriptor>::ScalarCubicInterpolationXYLineY3D(
    plint delta_, RescaleEngine<T, Descriptor> *rescaleEngine_) :
    delta(delta_), rescaleEngine(rescaleEngine_)
{
    PLB_PRECONDITION(delta == 1 || delta == -1);
}

template <typename T, template <typename U> class Descriptor>
void ScalarCubicInterpolationXYLineY3D<T, Descriptor>::process(
    Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
    BlockLattice3D<T, Descriptor> &fineLattice)
{
    Dot3D posCoarse = coarseLattice.getLocation();
    Dot3D posFine = fineLattice.getLocation();

    const plint x0 = coarseDomain.x0;
    const plint y0 = coarseDomain.y0;
    const plint z0 = coarseDomain.z0;

    const plint y1 = coarseDomain.y1;

    PLB_PRECONDITION(x0 == coarseDomain.x1);
    PLB_PRECONDITION(z0 == coarseDomain.z1);

    plint iX = x0;
    plint iZ = z0;

    T neighbors[4][3];
    std::vector<std::vector<T> > pop(12);

    for (plint iY = y0 - 1; iY <= y1; ++iY) {
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY - 1, iZ), pop[0]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop[1]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 1, iZ), pop[2]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 2, iZ), pop[3]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + delta, iY - 1, iZ), pop[4]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + delta, iY, iZ), pop[5]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + delta, iY + 1, iZ), pop[6]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + delta, iY + 2, iZ), pop[7]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2 * delta, iY - 1, iZ), pop[8]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2 * delta, iY, iZ), pop[9]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2 * delta, iY + 1, iZ), pop[10]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2 * delta, iY + 2, iZ), pop[11]);

        plint cellDim = pop[0].size();  // TODO get the size of the cells more generally

        std::vector<T> interpValues1(cellDim);
        std::vector<T> interpValues2(cellDim);
        std::vector<T> interpValues3(cellDim);
        std::vector<T> interpValues4(cellDim);

        for (plint iComp = 0; iComp < cellDim; ++iComp) {
            // get the 12 neighbors and interpolate where it is needed
            neighbors[0][0] = pop[0][iComp];
            neighbors[1][0] = pop[1][iComp];
            neighbors[2][0] = pop[2][iComp];
            neighbors[3][0] = pop[3][iComp];

            neighbors[0][1] = pop[4][iComp];
            neighbors[1][1] = pop[5][iComp];
            neighbors[2][1] = pop[6][iComp];
            neighbors[3][1] = pop[7][iComp];

            neighbors[0][2] = pop[8][iComp];
            neighbors[1][2] = pop[9][iComp];
            neighbors[2][2] = pop[10][iComp];
            neighbors[3][2] = pop[11][iComp];

            // interpolate the values according to the neighboring values, use centered schema
            std::vector<T> interpolatedValues = asymetricCubicInterpolation<T>(neighbors);

            interpValues1[iComp] = interpolatedValues[0];
            interpValues2[iComp] = interpolatedValues[1];
            interpValues3[iComp] = interpolatedValues[2];
            interpValues4[iComp] = interpolatedValues[3];
        }

        // convert the coarse coordinates to fine coordinates
        plint fineX = (iX + posCoarse.x) * 2 - posFine.x;
        plint fineY = (iY + posCoarse.y) * 2 - posFine.y;
        plint fineZ = (iZ + posCoarse.z) * 2 - posFine.z;

        // copy the 4 known values
        copyPopulations(pop[1], fineLattice.get(fineX, fineY, fineZ));
        copyPopulations(pop[5], fineLattice.get(fineX + 2 * delta, fineY, fineZ));
        copyPopulations(pop[2], fineLattice.get(fineX, fineY + 2, fineZ));
        copyPopulations(pop[6], fineLattice.get(fineX + 2 * delta, fineY + 2, fineZ));

        // assigning the computed interpolated values
        copyPopulations(interpValues1, fineLattice.get(fineX + delta, fineY, fineZ));
        copyPopulations(interpValues2, fineLattice.get(fineX, fineY + 1, fineZ));
        copyPopulations(interpValues3, fineLattice.get(fineX + delta, fineY + 1, fineZ));
        copyPopulations(interpValues4, fineLattice.get(fineX + delta, fineY + 2, fineZ));
    }
}

template <typename T, template <typename U> class Descriptor>
ScalarCubicInterpolationXYLineY3D<T, Descriptor>
    *ScalarCubicInterpolationXYLineY3D<T, Descriptor>::clone() const
{
    return new ScalarCubicInterpolationXYLineY3D<T, Descriptor>(delta, rescaleEngine->clone());
}

/* ********** ScalarCubicInterpolationXYCorner3D ************* */

template <typename T, template <typename U> class Descriptor>
ScalarCubicInterpolationXYCorner3D<T, Descriptor>::ScalarCubicInterpolationXYCorner3D(
    plint deltaX_, plint deltaY_, RescaleEngine<T, Descriptor> *rescaleEngine_) :
    deltaX(deltaX_), deltaY(deltaY_), rescaleEngine(rescaleEngine_)

{ }

template <typename T, template <typename U> class Descriptor>
void ScalarCubicInterpolationXYCorner3D<T, Descriptor>::process(
    Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
    BlockLattice3D<T, Descriptor> &fineLattice)
{
    Dot3D posCoarse = coarseLattice.getLocation();
    Dot3D posFine = fineLattice.getLocation();

    PLB_PRECONDITION(coarseDomain.x0 == coarseDomain.x1);
    PLB_PRECONDITION(coarseDomain.y0 == coarseDomain.y1);
    PLB_PRECONDITION(coarseDomain.z0 == coarseDomain.z1);

    const plint iX = coarseDomain.x0;
    const plint iY = coarseDomain.y0;
    const plint iZ = coarseDomain.z0;

    plint fineX = (iX + posCoarse.x) * 2 - posFine.x;
    plint fineY = (iY + posCoarse.y) * 2 - posFine.y;
    plint fineZ = (iZ + posCoarse.z) * 2 - posFine.z;

    T neighbors[3][3];
    std::vector<std::vector<T> > pop(9);

    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop[0]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + deltaX, iY, iZ), pop[1]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2 * deltaX, iY, iZ), pop[2]);

    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + deltaY, iZ), pop[3]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + deltaX, iY + deltaY, iZ), pop[4]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2 * deltaX, iY + deltaY, iZ), pop[5]);

    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 2 * deltaY, iZ), pop[6]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + deltaX, iY + 2 * deltaY, iZ), pop[7]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2 * deltaX, iY + 2 * deltaY, iZ), pop[8]);

    plint cellDim = pop[0].size();  // TODO get the size of the cells more generally

    std::vector<T> interpValues1(cellDim);
    std::vector<T> interpValues2(cellDim);
    std::vector<T> interpValues3(cellDim);

    // for every component in the decomposed values
    for (plint iComp = 0; iComp < cellDim; ++iComp) {
        // get the 9 neighbors and interpolate where it is needed
        neighbors[0][0] = pop[0][iComp];
        neighbors[1][0] = pop[1][iComp];
        neighbors[2][0] = pop[2][iComp];

        neighbors[0][1] = pop[3][iComp];
        neighbors[1][1] = pop[4][iComp];
        neighbors[2][1] = pop[5][iComp];

        neighbors[0][2] = pop[6][iComp];
        neighbors[1][2] = pop[7][iComp];
        neighbors[2][2] = pop[8][iComp];

        // interpolate the values according to the neighboring values, use centered schema
        std::vector<T> interpolatedValues = cornerInterpolation<T>(neighbors);

        // pack the values in vectors to use latter
        interpValues1[iComp] = interpolatedValues[0];
        interpValues2[iComp] = interpolatedValues[1];
        interpValues3[iComp] = interpolatedValues[2];
    }

    // copy the 4 known values
    copyPopulations(pop[0], fineLattice.get(fineX, fineY, fineZ));
    copyPopulations(pop[1], fineLattice.get(fineX + 2 * deltaX, fineY, fineZ));
    copyPopulations(pop[3], fineLattice.get(fineX, fineY + 2 * deltaY, fineZ));
    copyPopulations(pop[4], fineLattice.get(fineX + 2 * deltaX, fineY + 2 * deltaY, fineZ));

    // assigning the interpolated values
    copyPopulations(interpValues1, fineLattice.get(fineX, fineY + deltaY, fineZ));
    copyPopulations(interpValues2, fineLattice.get(fineX + deltaX, fineY, fineZ));
    copyPopulations(interpValues3, fineLattice.get(fineX + deltaX, fineY + deltaY, fineZ));
}

template <typename T, template <typename U> class Descriptor>
ScalarCubicInterpolationXYCorner3D<T, Descriptor>
    *ScalarCubicInterpolationXYCorner3D<T, Descriptor>::clone() const
{
    return new ScalarCubicInterpolationXYCorner3D<T, Descriptor>(
        deltaX, deltaY, rescaleEngine->clone());
}

/* ************* ScalarCubicInterpolationXZ ******************* */
template <typename T, template <typename U> class Descriptor>
void ScalarCubicInterpolationXZ<T, Descriptor>::process(
    Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
    BlockLattice3D<T, Descriptor> &fineLattice)
{
    Dot3D posCoarse = coarseLattice.getLocation();
    Dot3D posFine = fineLattice.getLocation();

    const plint x0 = coarseDomain.x0;
    const plint y0 = coarseDomain.y0;
    const plint z0 = coarseDomain.z0;

    const plint x1 = coarseDomain.x1;
    const plint z1 = coarseDomain.z1;

    PLB_PRECONDITION(y0 == coarseDomain.y1);

    plint iY = y0;
    Array<Array<T, 2>, 5> centeredPoints;
    centeredPoints[0] = Array<T, 2>(1.0, 1.5);
    centeredPoints[1] = Array<T, 2>(1.5, 1.0);
    centeredPoints[2] = Array<T, 2>(1.5, 1.5);
    centeredPoints[3] = Array<T, 2>(1.5, 2.0);
    centeredPoints[4] = Array<T, 2>(2.0, 1.5);

    T neighbors[4][4];
    std::vector<std::vector<T> > pop(16);

    for (plint iX = x0 - 1; iX <= x1; ++iX) {
        for (plint iZ = z0 - 1; iZ <= z1; ++iZ) {
            // extracting and rescaling the known 16 coarse values
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - 1, iY, iZ - 1), pop[0]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ - 1), pop[1]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 1, iY, iZ - 1), pop[2]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2, iY, iZ - 1), pop[3]);

            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - 1, iY, iZ), pop[4]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop[5]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 1, iY, iZ), pop[6]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2, iY, iZ), pop[7]);

            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - 1, iY, iZ + 1), pop[8]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + 1), pop[9]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 1, iY, iZ + 1), pop[10]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2, iY, iZ + 1), pop[11]);

            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - 1, iY, iZ + 2), pop[12]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + 2), pop[13]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 1, iY, iZ + 2), pop[14]);
            rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2, iY, iZ + 2), pop[15]);

            plint cellDim = pop[0].size();  // TODO get the size of the cells more generally
            // the containers for the interpolated values
            std::vector<T> interpValues1(cellDim);
            std::vector<T> interpValues2(cellDim);
            std::vector<T> interpValues3(cellDim);
            std::vector<T> interpValues4(cellDim);
            std::vector<T> interpValues5(cellDim);

            for (plint iComp = 0; iComp < cellDim; ++iComp) {
                // get the 16 neighbors and interpolate where it is needed
                neighbors[0][0] = pop[0][iComp];
                neighbors[1][0] = pop[1][iComp];
                neighbors[2][0] = pop[2][iComp];
                neighbors[3][0] = pop[3][iComp];

                neighbors[0][1] = pop[4][iComp];
                neighbors[1][1] = pop[5][iComp];
                neighbors[2][1] = pop[6][iComp];
                neighbors[3][1] = pop[7][iComp];

                neighbors[0][2] = pop[8][iComp];
                neighbors[1][2] = pop[9][iComp];
                neighbors[2][2] = pop[10][iComp];
                neighbors[3][2] = pop[11][iComp];

                neighbors[0][3] = pop[12][iComp];
                neighbors[1][3] = pop[13][iComp];
                neighbors[2][3] = pop[14][iComp];
                neighbors[3][3] = pop[15][iComp];

                // interpolate the values according to the neighboring values, use centered schema
                std::vector<T> interpolatedValues = symetricCubicInterpolation<T>(neighbors);

                interpValues1[iComp] = interpolatedValues[0];
                interpValues2[iComp] = interpolatedValues[1];
                interpValues3[iComp] = interpolatedValues[2];
                interpValues4[iComp] = interpolatedValues[3];
                interpValues5[iComp] = interpolatedValues[4];
            }

            // convert the coarse coordinates to fine coordinates
            plint fineX = (iX + posCoarse.x) * 2 - posFine.x;
            plint fineY = (iY + posCoarse.y) * 2 - posFine.y;
            plint fineZ = (iZ + posCoarse.z) * 2 - posFine.z;

            // assigning the known values
            copyPopulations(pop[5], fineLattice.get(fineX, fineY, fineZ));
            copyPopulations(pop[6], fineLattice.get(fineX + 2, fineY, fineZ));
            copyPopulations(pop[9], fineLattice.get(fineX, fineY, fineZ + 2));
            copyPopulations(pop[10], fineLattice.get(fineX + 2, fineY, fineZ + 2));

            // assigning the interpolated values
            copyPopulations(interpValues1, fineLattice.get(fineX, fineY, fineZ + 1));
            copyPopulations(interpValues2, fineLattice.get(fineX + 1, fineY, fineZ));
            copyPopulations(interpValues3, fineLattice.get(fineX + 1, fineY, fineZ + 1));
            copyPopulations(interpValues4, fineLattice.get(fineX + 1, fineY, fineZ + 2));
            copyPopulations(interpValues5, fineLattice.get(fineX + 2, fineY, fineZ + 1));
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ScalarCubicInterpolationXZ<T, Descriptor> *ScalarCubicInterpolationXZ<T, Descriptor>::clone() const
{
    return new ScalarCubicInterpolationXZ<T, Descriptor>(rescaleEngine->clone());
}

/* **************** ScalarCubicInterpolationXYLineX3D ********************** */

template <typename T, template <typename U> class Descriptor>
ScalarCubicInterpolationXZLineX3D<T, Descriptor>::ScalarCubicInterpolationXZLineX3D(
    plint delta_, RescaleEngine<T, Descriptor> *rescaleEngine_) :
    delta(delta_), rescaleEngine(rescaleEngine_)
{
    PLB_PRECONDITION(delta == 1 || delta == -1);
    points[0] = Array<T, 2>(1.0, 0.5);
    points[1] = Array<T, 2>(1.5, 0.0);
    points[2] = Array<T, 2>(1.5, 0.5);
    points[3] = Array<T, 2>(2.0, 0.5);
}

template <typename T, template <typename U> class Descriptor>
void ScalarCubicInterpolationXZLineX3D<T, Descriptor>::process(
    Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
    BlockLattice3D<T, Descriptor> &fineLattice)
{
    Dot3D posCoarse = coarseLattice.getLocation();
    Dot3D posFine = fineLattice.getLocation();

    const plint x0 = coarseDomain.x0;
    const plint y0 = coarseDomain.y0;
    const plint z0 = coarseDomain.z0;

    const plint x1 = coarseDomain.x1;

    PLB_PRECONDITION(y0 == coarseDomain.y1);
    PLB_PRECONDITION(z0 == coarseDomain.z1);

    plint iY = y0;
    plint iZ = z0;

    T neighbors[4][3];
    std::vector<std::vector<T> > pop(12);

    for (plint iX = x0 - 1; iX <= x1; ++iX) {
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - 1, iY, iZ), pop[0]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop[1]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 1, iY, iZ), pop[2]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2, iY, iZ), pop[3]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - 1, iY, iZ + delta), pop[4]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + delta), pop[5]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 1, iY, iZ + delta), pop[6]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2, iY, iZ + delta), pop[7]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - 1, iY, iZ + 2 * delta), pop[8]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + 2 * delta), pop[9]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 1, iY, iZ + 2 * delta), pop[10]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2, iY, iZ + 2 * delta), pop[11]);

        plint cellDim = pop[0].size();  // TODO get the size of the cells more generally

        std::vector<T> interpValues1(cellDim);
        std::vector<T> interpValues2(cellDim);
        std::vector<T> interpValues3(cellDim);
        std::vector<T> interpValues4(cellDim);

        for (plint iComp = 0; iComp < cellDim; ++iComp) {
            // get the 12 neighbors and interpolate where it is needed
            neighbors[0][0] = pop[0][iComp];
            neighbors[1][0] = pop[1][iComp];
            neighbors[2][0] = pop[2][iComp];
            neighbors[3][0] = pop[3][iComp];

            neighbors[0][1] = pop[4][iComp];
            neighbors[1][1] = pop[5][iComp];
            neighbors[2][1] = pop[6][iComp];
            neighbors[3][1] = pop[7][iComp];

            neighbors[0][2] = pop[8][iComp];
            neighbors[1][2] = pop[9][iComp];
            neighbors[2][2] = pop[10][iComp];
            neighbors[3][2] = pop[11][iComp];

            // interpolate the values according to the neighboring values, use centered schema
            std::vector<T> interpolatedValues = asymetricCubicInterpolation<T>(neighbors);

            interpValues1[iComp] = interpolatedValues[0];
            interpValues2[iComp] = interpolatedValues[1];
            interpValues3[iComp] = interpolatedValues[2];
            interpValues4[iComp] = interpolatedValues[3];
        }

        // convert the coarse coordinates to fine coordinates
        plint fineX = (iX + posCoarse.x) * 2 - posFine.x;
        plint fineY = (iY + posCoarse.y) * 2 - posFine.y;
        plint fineZ = (iZ + posCoarse.z) * 2 - posFine.z;

        // copy the 4 known values
        copyPopulations(pop[1], fineLattice.get(fineX, fineY, fineZ));
        copyPopulations(pop[2], fineLattice.get(fineX + 2, fineY, fineZ));
        copyPopulations(pop[5], fineLattice.get(fineX, fineY, fineZ + 2 * delta));
        copyPopulations(pop[6], fineLattice.get(fineX + 2, fineY, fineZ + 2 * delta));

        // assigning the computed interpolated values
        copyPopulations(interpValues1, fineLattice.get(fineX, fineY, fineZ + delta));
        copyPopulations(interpValues2, fineLattice.get(fineX + 1, fineY, fineZ));
        copyPopulations(interpValues3, fineLattice.get(fineX + 1, fineY, fineZ + delta));
        copyPopulations(interpValues4, fineLattice.get(fineX + 2, fineY, fineZ + delta));
    }
}

template <typename T, template <typename U> class Descriptor>
ScalarCubicInterpolationXZLineX3D<T, Descriptor>
    *ScalarCubicInterpolationXZLineX3D<T, Descriptor>::clone() const
{
    return new ScalarCubicInterpolationXZLineX3D<T, Descriptor>(delta, rescaleEngine->clone());
}

/* ************ ScalarCubicInterpolationXZLineZ3D *************** */

template <typename T, template <typename U> class Descriptor>
ScalarCubicInterpolationXZLineZ3D<T, Descriptor>::ScalarCubicInterpolationXZLineZ3D(
    plint delta_, RescaleEngine<T, Descriptor> *rescaleEngine_) :
    delta(delta_), rescaleEngine(rescaleEngine_)
{
    PLB_PRECONDITION(delta == 1 || delta == -1);
}

template <typename T, template <typename U> class Descriptor>
void ScalarCubicInterpolationXZLineZ3D<T, Descriptor>::process(
    Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
    BlockLattice3D<T, Descriptor> &fineLattice)
{
    Dot3D posCoarse = coarseLattice.getLocation();
    Dot3D posFine = fineLattice.getLocation();

    const plint x0 = coarseDomain.x0;
    const plint y0 = coarseDomain.y0;
    const plint z0 = coarseDomain.z0;

    const plint z1 = coarseDomain.z1;

    PLB_PRECONDITION(x0 == coarseDomain.x1);
    PLB_PRECONDITION(y0 == coarseDomain.y1);

    plint iX = x0;
    plint iY = y0;

    T neighbors[4][3];
    std::vector<std::vector<T> > pop(12);

    for (plint iZ = z0 - 1; iZ <= z1; ++iZ) {
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ - 1), pop[0]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop[1]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + 1), pop[2]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + 2), pop[3]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + delta, iY, iZ - 1), pop[4]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + delta, iY, iZ), pop[5]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + delta, iY, iZ + 1), pop[6]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + delta, iY, iZ + 2), pop[7]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2 * delta, iY, iZ - 1), pop[8]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2 * delta, iY, iZ), pop[9]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2 * delta, iY, iZ + 1), pop[10]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2 * delta, iY, iZ + 2), pop[11]);

        plint cellDim = pop[0].size();  // TODO get the size of the cells more generally

        std::vector<T> interpValues1(cellDim);
        std::vector<T> interpValues2(cellDim);
        std::vector<T> interpValues3(cellDim);
        std::vector<T> interpValues4(cellDim);

        for (plint iComp = 0; iComp < cellDim; ++iComp) {
            // get the 12 neighbors and interpolate where it is needed
            neighbors[0][0] = pop[0][iComp];
            neighbors[1][0] = pop[1][iComp];
            neighbors[2][0] = pop[2][iComp];
            neighbors[3][0] = pop[3][iComp];

            neighbors[0][1] = pop[4][iComp];
            neighbors[1][1] = pop[5][iComp];
            neighbors[2][1] = pop[6][iComp];
            neighbors[3][1] = pop[7][iComp];

            neighbors[0][2] = pop[8][iComp];
            neighbors[1][2] = pop[9][iComp];
            neighbors[2][2] = pop[10][iComp];
            neighbors[3][2] = pop[11][iComp];

            // interpolate the values according to the neighboring values, use centered schema
            std::vector<T> interpolatedValues = asymetricCubicInterpolation<T>(neighbors);

            interpValues1[iComp] = interpolatedValues[0];
            interpValues2[iComp] = interpolatedValues[1];
            interpValues3[iComp] = interpolatedValues[2];
            interpValues4[iComp] = interpolatedValues[3];
        }

        // convert the coarse coordinates to fine coordinates
        plint fineX = (iX + posCoarse.x) * 2 - posFine.x;
        plint fineY = (iY + posCoarse.y) * 2 - posFine.y;
        plint fineZ = (iZ + posCoarse.z) * 2 - posFine.z;

        // copy the 4 known values
        copyPopulations(pop[1], fineLattice.get(fineX, fineY, fineZ));
        copyPopulations(pop[5], fineLattice.get(fineX + 2 * delta, fineY, fineZ));
        copyPopulations(pop[2], fineLattice.get(fineX, fineY, fineZ + 2));
        copyPopulations(pop[6], fineLattice.get(fineX + 2 * delta, fineY, fineZ + 2));

        // assigning the computed interpolated values
        copyPopulations(interpValues1, fineLattice.get(fineX + delta, fineY, fineZ));
        copyPopulations(interpValues2, fineLattice.get(fineX, fineY, fineZ + 1));
        copyPopulations(interpValues3, fineLattice.get(fineX + delta, fineY, fineZ + 1));
        copyPopulations(interpValues4, fineLattice.get(fineX + delta, fineY, fineZ + 2));
    }
}

template <typename T, template <typename U> class Descriptor>
ScalarCubicInterpolationXZLineZ3D<T, Descriptor>
    *ScalarCubicInterpolationXZLineZ3D<T, Descriptor>::clone() const
{
    return new ScalarCubicInterpolationXZLineZ3D<T, Descriptor>(delta, rescaleEngine->clone());
}

/* ********** ScalarCubicInterpolationXZCorner3D ************* */

template <typename T, template <typename U> class Descriptor>
ScalarCubicInterpolationXZCorner3D<T, Descriptor>::ScalarCubicInterpolationXZCorner3D(
    plint deltaX_, plint deltaZ_, RescaleEngine<T, Descriptor> *rescaleEngine_) :
    deltaX(deltaX_), deltaZ(deltaZ_), rescaleEngine(rescaleEngine_)
{ }

template <typename T, template <typename U> class Descriptor>
void ScalarCubicInterpolationXZCorner3D<T, Descriptor>::process(
    Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
    BlockLattice3D<T, Descriptor> &fineLattice)
{
    Dot3D posCoarse = coarseLattice.getLocation();
    Dot3D posFine = fineLattice.getLocation();

    PLB_PRECONDITION(coarseDomain.x0 == coarseDomain.x1);
    PLB_PRECONDITION(coarseDomain.y0 == coarseDomain.y1);
    PLB_PRECONDITION(coarseDomain.z0 == coarseDomain.z1);

    const plint iX = coarseDomain.x0;
    const plint iY = coarseDomain.y0;
    const plint iZ = coarseDomain.z0;

    plint fineX = (iX + posCoarse.x) * 2 - posFine.x;
    plint fineY = (iY + posCoarse.y) * 2 - posFine.y;
    plint fineZ = (iZ + posCoarse.z) * 2 - posFine.z;

    T neighbors[3][3];
    std::vector<std::vector<T> > pop(9);

    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop[0]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + deltaX, iY, iZ), pop[1]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2 * deltaX, iY, iZ), pop[2]);

    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + deltaZ), pop[3]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + deltaX, iY, iZ + deltaZ), pop[4]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2 * deltaX, iY, iZ + deltaZ), pop[5]);

    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + 2 * deltaZ), pop[6]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + deltaX, iY, iZ + 2 * deltaZ), pop[7]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2 * deltaX, iY, iZ + 2 * deltaZ), pop[8]);

    plint cellDim = pop[0].size();  // TODO get the size of the cells more generally

    std::vector<T> interpValues1(cellDim);
    std::vector<T> interpValues2(cellDim);
    std::vector<T> interpValues3(cellDim);

    // for every component in the decomposed values
    for (plint iComp = 0; iComp < cellDim; ++iComp) {
        // get the 9 neighbors and interpolate where it is needed
        neighbors[0][0] = pop[0][iComp];
        neighbors[1][0] = pop[1][iComp];
        neighbors[2][0] = pop[2][iComp];

        neighbors[0][1] = pop[3][iComp];
        neighbors[1][1] = pop[4][iComp];
        neighbors[2][1] = pop[5][iComp];

        neighbors[0][2] = pop[6][iComp];
        neighbors[1][2] = pop[7][iComp];
        neighbors[2][2] = pop[8][iComp];

        // interpolate the values according to the neighboring values, use centered schema
        std::vector<T> interpolatedValues = cornerInterpolation<T>(neighbors);

        // pack the values in vectors to use latter
        interpValues1[iComp] = interpolatedValues[0];
        interpValues2[iComp] = interpolatedValues[1];
        interpValues3[iComp] = interpolatedValues[2];
    }

    // copy the 4 known values
    copyPopulations(pop[0], fineLattice.get(fineX, fineY, fineZ));
    copyPopulations(pop[1], fineLattice.get(fineX + 2 * deltaX, fineY, fineZ));
    copyPopulations(pop[3], fineLattice.get(fineX, fineY, fineZ + 2 * deltaZ));
    copyPopulations(pop[4], fineLattice.get(fineX + 2 * deltaX, fineY, fineZ + 2 * deltaZ));

    // assigning the interpolated values
    copyPopulations(interpValues1, fineLattice.get(fineX, fineY, fineZ + deltaZ));
    copyPopulations(interpValues2, fineLattice.get(fineX + deltaX, fineY, fineZ));
    copyPopulations(interpValues3, fineLattice.get(fineX + deltaX, fineY, fineZ + deltaZ));
}

template <typename T, template <typename U> class Descriptor>
ScalarCubicInterpolationXZCorner3D<T, Descriptor>
    *ScalarCubicInterpolationXZCorner3D<T, Descriptor>::clone() const
{
    return new ScalarCubicInterpolationXZCorner3D<T, Descriptor>(
        deltaX, deltaZ, rescaleEngine->clone());
}

}  // namespace plb

#endif  // FINE_GRID_PROCESSORS_3D_HH
