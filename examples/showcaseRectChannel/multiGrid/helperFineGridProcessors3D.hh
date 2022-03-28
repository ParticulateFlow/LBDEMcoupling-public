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
 * Dynamics and data processors used to implement 3D grid refinement -- header file.
 */

#ifndef HELPER_FINE_GRID_PROCESSORS_3D_HH
#define HELPER_FINE_GRID_PROCESSORS_3D_HH

#include "multiGrid/helperFineGridProcessors3D.h"

namespace plb {

//////////////////////////////////////////////////////////////////////////////
/////////////////////          YZ plane           ////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
CubicInterpolationYZLineYHelper3D<T, Descriptor>::CubicInterpolationYZLineYHelper3D(
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
void CubicInterpolationYZLineYHelper3D<T, Descriptor>::process(
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
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY - 1, iZ - delta), pop[0]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ - delta), pop[1]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 1, iZ - delta), pop[2]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 2, iZ - delta), pop[3]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY - 1, iZ), pop[4]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop[5]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 1, iZ), pop[6]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 2, iZ), pop[7]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY - 1, iZ + delta), pop[8]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + delta), pop[9]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 1, iZ + delta), pop[10]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 2, iZ + delta), pop[11]);

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

        // assigning the computed interpolated values
        copyPopulations(interpValues1, fineLattice.get(fineX, fineY, fineZ - delta));
        copyPopulations(interpValues2, fineLattice.get(fineX, fineY + 1, fineZ - 2 * delta));
        copyPopulations(interpValues3, fineLattice.get(fineX, fineY + 1, fineZ - delta));
        copyPopulations(interpValues4, fineLattice.get(fineX, fineY + 2, fineZ - delta));
    }
}

template <typename T, template <typename U> class Descriptor>
CubicInterpolationYZLineYHelper3D<T, Descriptor>
    *CubicInterpolationYZLineYHelper3D<T, Descriptor>::clone() const
{
    return new CubicInterpolationYZLineYHelper3D<T, Descriptor>(delta, rescaleEngine->clone());
}

template <typename T, template <typename U> class Descriptor>
CubicInterpolationYZLineZHelper3D<T, Descriptor>::CubicInterpolationYZLineZHelper3D(
    plint delta_, RescaleEngine<T, Descriptor> *rescaleEngine_) :
    delta(delta_), rescaleEngine(rescaleEngine_)
{
    points[0] = Array<T, 2>(0.0, 1.5);
    points[1] = Array<T, 2>(0.5, 1.0);
    points[2] = Array<T, 2>(0.5, 1.5);
    points[3] = Array<T, 2>(0.5, 2.0);
}

template <typename T, template <typename U> class Descriptor>
void CubicInterpolationYZLineZHelper3D<T, Descriptor>::process(
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
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY - delta, iZ - 1), pop[0]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY - delta, iZ), pop[1]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY - delta, iZ + 1), pop[2]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY - delta, iZ + 2), pop[3]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ - 1), pop[4]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop[5]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + 1), pop[6]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + 2), pop[7]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + delta, iZ - 1), pop[8]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + delta, iZ), pop[9]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + delta, iZ + 1), pop[10]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + delta, iZ + 2), pop[11]);

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

        // assigning the computed interpolated values
        copyPopulations(interpValues1, fineLattice.get(fineX, fineY - delta, fineZ));
        copyPopulations(interpValues2, fineLattice.get(fineX, fineY - 2 * delta, fineZ + 1));
        copyPopulations(interpValues3, fineLattice.get(fineX, fineY - delta, fineZ + 1));
        copyPopulations(interpValues4, fineLattice.get(fineX, fineY - delta, fineZ + 2));
    }
}

template <typename T, template <typename U> class Descriptor>
CubicInterpolationYZLineZHelper3D<T, Descriptor>
    *CubicInterpolationYZLineZHelper3D<T, Descriptor>::clone() const
{
    return new CubicInterpolationYZLineZHelper3D<T, Descriptor>(delta, rescaleEngine->clone());
}

/* ************** CubicInterpolationYZCornerHelper3D **************** */

template <typename T, template <typename U> class Descriptor>
CubicInterpolationYZCornerHelper3D<T, Descriptor>::CubicInterpolationYZCornerHelper3D(
    plint deltaY_, plint deltaZ_, RescaleEngine<T, Descriptor> *rescaleEngine_) :
    deltaY(deltaY_), deltaZ(deltaZ_), rescaleEngine(rescaleEngine_)
{ }

template <typename T, template <typename U> class Descriptor>
void CubicInterpolationYZCornerHelper3D<T, Descriptor>::process(
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

    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY - deltaY, iZ - deltaZ), pop[0]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ - deltaZ), pop[1]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + deltaY, iZ - deltaZ), pop[2]);

    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY - deltaY, iZ), pop[3]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop[4]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + deltaY, iZ), pop[5]);

    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY - deltaY, iZ + deltaZ), pop[6]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + deltaZ), pop[7]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + deltaY, iZ + deltaZ), pop[8]);

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

    // assigning the interpolated values
    copyPopulations(interpValues1, fineLattice.get(fineX, fineY - 2 * deltaY, fineZ - deltaZ));
    copyPopulations(interpValues2, fineLattice.get(fineX, fineY - deltaY, fineZ - 2 * deltaZ));
    copyPopulations(interpValues3, fineLattice.get(fineX, fineY - deltaY, fineZ - deltaZ));
}

template <typename T, template <typename U> class Descriptor>
CubicInterpolationYZCornerHelper3D<T, Descriptor>
    *CubicInterpolationYZCornerHelper3D<T, Descriptor>::clone() const
{
    return new CubicInterpolationYZCornerHelper3D<T, Descriptor>(
        deltaY, deltaZ, rescaleEngine->clone());
}

//////////////////////////////////////////////////////////////////////////////
/////////////////////          XY plane           ////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
CubicInterpolationXYLineXHelper3D<T, Descriptor>::CubicInterpolationXYLineXHelper3D(
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
void CubicInterpolationXYLineXHelper3D<T, Descriptor>::process(
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
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - 1, iY - delta, iZ), pop[0]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY - delta, iZ), pop[1]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 1, iY - delta, iZ), pop[2]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2, iY - delta, iZ), pop[3]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - 1, iY, iZ), pop[4]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop[5]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 1, iY, iZ), pop[6]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2, iY, iZ), pop[7]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - 1, iY + delta, iZ), pop[8]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + delta, iZ), pop[9]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 1, iY + delta, iZ), pop[10]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2, iY + delta, iZ), pop[11]);

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

        // assigning the computed interpolated values
        copyPopulations(interpValues1, fineLattice.get(fineX, fineY - delta, fineZ));
        copyPopulations(interpValues2, fineLattice.get(fineX + 1, fineY - 2 * delta, fineZ));
        copyPopulations(interpValues3, fineLattice.get(fineX + 1, fineY - delta, fineZ));
        copyPopulations(interpValues4, fineLattice.get(fineX + 2, fineY - delta, fineZ));
    }
}

template <typename T, template <typename U> class Descriptor>
CubicInterpolationXYLineXHelper3D<T, Descriptor>
    *CubicInterpolationXYLineXHelper3D<T, Descriptor>::clone() const
{
    return new CubicInterpolationXYLineXHelper3D<T, Descriptor>(delta, rescaleEngine->clone());
}

template <typename T, template <typename U> class Descriptor>
CubicInterpolationXYLineYHelper3D<T, Descriptor>::CubicInterpolationXYLineYHelper3D(
    plint delta_, RescaleEngine<T, Descriptor> *rescaleEngine_) :
    delta(delta_), rescaleEngine(rescaleEngine_)
{
    PLB_PRECONDITION(delta == 1 || delta == -1);
    points[0] = Array<T, 2>(0.0, 1.5);
    points[1] = Array<T, 2>(0.5, 1.0);
    points[2] = Array<T, 2>(0.5, 1.5);
    points[3] = Array<T, 2>(0.5, 2.0);
}

template <typename T, template <typename U> class Descriptor>
void CubicInterpolationXYLineYHelper3D<T, Descriptor>::process(
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
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - delta, iY - 1, iZ), pop[0]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - delta, iY, iZ), pop[1]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - delta, iY + 1, iZ), pop[2]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - delta, iY + 2, iZ), pop[3]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY - 1, iZ), pop[4]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop[5]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 1, iZ), pop[6]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + 2, iZ), pop[7]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + delta, iY - 1, iZ), pop[8]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + delta, iY, iZ), pop[9]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + delta, iY + 1, iZ), pop[10]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + delta, iY + 2, iZ), pop[11]);

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

        // assigning the computed interpolated values
        copyPopulations(interpValues1, fineLattice.get(fineX - delta, fineY, fineZ));
        copyPopulations(interpValues2, fineLattice.get(fineX - 2 * delta, fineY + 1, fineZ));
        copyPopulations(interpValues3, fineLattice.get(fineX - delta, fineY + 1, fineZ));
        copyPopulations(interpValues4, fineLattice.get(fineX - delta, fineY + 2, fineZ));
    }
}

template <typename T, template <typename U> class Descriptor>
CubicInterpolationXYLineYHelper3D<T, Descriptor>
    *CubicInterpolationXYLineYHelper3D<T, Descriptor>::clone() const
{
    return new CubicInterpolationXYLineYHelper3D<T, Descriptor>(delta, rescaleEngine->clone());
}

/* ************* CubicInterpolationXYCornerHelper3D ************* */

template <typename T, template <typename U> class Descriptor>
CubicInterpolationXYCornerHelper3D<T, Descriptor>::CubicInterpolationXYCornerHelper3D(
    plint deltaX_, plint deltaY_, RescaleEngine<T, Descriptor> *rescaleEngine_) :
    deltaX(deltaX_), deltaY(deltaY_), rescaleEngine(rescaleEngine_)
{ }

template <typename T, template <typename U> class Descriptor>
void CubicInterpolationXYCornerHelper3D<T, Descriptor>::process(
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

    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - deltaX, iY - deltaY, iZ), pop[0]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY - deltaY, iZ), pop[1]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + deltaX, iY - deltaY, iZ), pop[2]);

    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - deltaX, iY, iZ), pop[3]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop[4]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + deltaX, iY, iZ), pop[5]);

    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - deltaX, iY + deltaY, iZ), pop[6]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY + deltaY, iZ), pop[7]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + deltaX, iY + deltaY, iZ), pop[8]);

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

    // assigning the interpolated values
    copyPopulations(interpValues1, fineLattice.get(fineX - 2 * deltaX, fineY - deltaY, fineZ));
    copyPopulations(interpValues2, fineLattice.get(fineX - deltaX, fineY - 2 * deltaY, fineZ));
    copyPopulations(interpValues3, fineLattice.get(fineX - deltaX, fineY - deltaY, fineZ));
}

template <typename T, template <typename U> class Descriptor>
CubicInterpolationXYCornerHelper3D<T, Descriptor>
    *CubicInterpolationXYCornerHelper3D<T, Descriptor>::clone() const
{
    return new CubicInterpolationXYCornerHelper3D<T, Descriptor>(
        deltaX, deltaY, rescaleEngine->clone());
}

//////////////////////////////////////////////////////////////////////////////
/////////////////////          XZ plane           ////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
CubicInterpolationXZLineXHelper3D<T, Descriptor>::CubicInterpolationXZLineXHelper3D(
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
void CubicInterpolationXZLineXHelper3D<T, Descriptor>::process(
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
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - 1, iY, iZ - delta), pop[0]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ - delta), pop[1]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 1, iY, iZ - delta), pop[2]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2, iY, iZ - delta), pop[3]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - 1, iY, iZ), pop[4]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop[5]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 1, iY, iZ), pop[6]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2, iY, iZ), pop[7]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - 1, iY, iZ + delta), pop[8]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + delta), pop[9]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 1, iY, iZ + delta), pop[10]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + 2, iY, iZ + delta), pop[11]);

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

        // assigning the computed interpolated values
        copyPopulations(interpValues1, fineLattice.get(fineX, fineY, fineZ - delta));
        copyPopulations(interpValues2, fineLattice.get(fineX + 1, fineY, fineZ - 2 * delta));
        copyPopulations(interpValues3, fineLattice.get(fineX + 1, fineY, fineZ - delta));
        copyPopulations(interpValues4, fineLattice.get(fineX + 2, fineY, fineZ - delta));
    }
}

template <typename T, template <typename U> class Descriptor>
CubicInterpolationXZLineXHelper3D<T, Descriptor>
    *CubicInterpolationXZLineXHelper3D<T, Descriptor>::clone() const
{
    return new CubicInterpolationXZLineXHelper3D<T, Descriptor>(delta, rescaleEngine->clone());
}

/* ************* CubicInterpolationXZLineZHelper3D ************* */
template <typename T, template <typename U> class Descriptor>
CubicInterpolationXZLineZHelper3D<T, Descriptor>::CubicInterpolationXZLineZHelper3D(
    plint delta_, RescaleEngine<T, Descriptor> *rescaleEngine_) :
    delta(delta_), rescaleEngine(rescaleEngine_)
{
    PLB_PRECONDITION(delta == 1 || delta == -1);
    points[0] = Array<T, 2>(0.0, 1.5);
    points[1] = Array<T, 2>(0.5, 1.0);
    points[2] = Array<T, 2>(0.5, 1.5);
    points[3] = Array<T, 2>(0.5, 2.0);
}

template <typename T, template <typename U> class Descriptor>
void CubicInterpolationXZLineZHelper3D<T, Descriptor>::process(
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
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - delta, iY, iZ - 1), pop[0]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - delta, iY, iZ), pop[1]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - delta, iY, iZ + 1), pop[2]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - delta, iY, iZ + 2), pop[3]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ - 1), pop[4]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop[5]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + 1), pop[6]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + 2), pop[7]);

        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + delta, iY, iZ - 1), pop[8]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + delta, iY, iZ), pop[9]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + delta, iY, iZ + 1), pop[10]);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + delta, iY, iZ + 2), pop[11]);

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

        // assigning the computed interpolated values
        copyPopulations(interpValues1, fineLattice.get(fineX - delta, fineY, fineZ));
        copyPopulations(interpValues2, fineLattice.get(fineX - 2 * delta, fineY, fineZ + 1));
        copyPopulations(interpValues3, fineLattice.get(fineX - delta, fineY, fineZ + 1));
        copyPopulations(interpValues4, fineLattice.get(fineX - delta, fineY, fineZ + 2));
    }
}

template <typename T, template <typename U> class Descriptor>
CubicInterpolationXZLineZHelper3D<T, Descriptor>
    *CubicInterpolationXZLineZHelper3D<T, Descriptor>::clone() const
{
    return new CubicInterpolationXZLineZHelper3D<T, Descriptor>(delta, rescaleEngine->clone());
}

/* ************* CubicInterpolationXZCornerHelper3D *************** */

template <typename T, template <typename U> class Descriptor>
CubicInterpolationXZCornerHelper3D<T, Descriptor>::CubicInterpolationXZCornerHelper3D(
    plint deltaX_, plint deltaZ_, RescaleEngine<T, Descriptor> *rescaleEngine_) :
    deltaX(deltaX_), deltaZ(deltaZ_), rescaleEngine(rescaleEngine_)
{ }

template <typename T, template <typename U> class Descriptor>
void CubicInterpolationXZCornerHelper3D<T, Descriptor>::process(
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

    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - deltaX, iY, iZ - deltaZ), pop[0]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ - deltaZ), pop[1]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + deltaX, iY, iZ - deltaZ), pop[2]);

    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - deltaX, iY, iZ), pop[3]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop[4]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + deltaX, iY, iZ), pop[5]);

    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX - deltaX, iY, iZ + deltaZ), pop[6]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ + deltaZ), pop[7]);
    rescaleEngine->scaleCoarseFine(coarseLattice.get(iX + deltaX, iY, iZ + deltaZ), pop[8]);

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

    // assigning the interpolated values
    copyPopulations(interpValues1, fineLattice.get(fineX - 2 * deltaX, fineY, fineZ - deltaZ));
    copyPopulations(interpValues2, fineLattice.get(fineX - deltaX, fineY, fineZ - 2 * deltaZ));
    copyPopulations(interpValues3, fineLattice.get(fineX - deltaX, fineY, fineZ - deltaZ));
}

template <typename T, template <typename U> class Descriptor>
CubicInterpolationXZCornerHelper3D<T, Descriptor>
    *CubicInterpolationXZCornerHelper3D<T, Descriptor>::clone() const
{
    return new CubicInterpolationXZCornerHelper3D<T, Descriptor>(
        deltaX, deltaZ, rescaleEngine->clone());
}

}  // namespace plb

#endif  // HELPER_FINE_GRID_PROCESSORS_3D_HH
