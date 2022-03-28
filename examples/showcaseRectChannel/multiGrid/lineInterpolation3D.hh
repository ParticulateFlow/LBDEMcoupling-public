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
 * Interpolations over a line in coarse and fine coordinates in 3D -- implementation file.
 */
#ifndef LINE_INTERPOLATION_3D_HH
#define LINE_INTERPOLATION_3D_HH

#include "multiGrid/lineInterpolation3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
void cubicInterpolationOverCoarseLineX(
    Box3D domain, BlockLattice3D<T, Descriptor> &coarseLattice,
    BlockLattice3D<T, Descriptor> &fineLattice, RescaleEngine<T, Descriptor> *rescaleEngine)
{
    // domain is given in coarse coordinates
    PLB_ASSERT(domain.y0 == domain.y1 && domain.z0 == domain.z1);
    plint iY = domain.y0;
    plint iZ = domain.z0;
    plint x0 = domain.x0;
    plint x1 = domain.x1;

    CoarseToFineConverter3D<T, Descriptor> converter(coarseLattice, fineLattice);
    Dot3D delta(1, 0, 0);
    plint whichTime = 1;
    std::vector<T> pop1;
    std::vector<T> pop2;
    std::vector<T> pop3;
    std::vector<T> pop4;

    Dot3D startExcess(x0 - 2 * delta.x, iY - 2 * delta.y, iZ - delta.z);
    Dot3D startFineExcess(
        converter.fineX(x0) - 2 * delta.x, converter.fineY(iY) - 2 * delta.y,
        converter.fineZ(iZ) - 2 * delta.z);
    bool startOutOfRange =
        (!contained(startExcess, coarseLattice.getBoundingBox())
         || !contained(startFineExcess, fineLattice.getBoundingBox()));

    Dot3D endExcess(x1 + 2 * delta.x, iY + 2 * delta.y, iZ + delta.z);
    Dot3D endFineExcess(
        converter.fineX(x1) + 2 * delta.x, converter.fineY(iY) + 2 * delta.y,
        converter.fineZ(iZ) + 2 * delta.z);
    bool endOutOfRange =
        (!contained(endExcess, coarseLattice.getBoundingBox())
         || !contained(endFineExcess, fineLattice.getBoundingBox()));

    // check if interpolation is needed in the leftmost part
    if (!startOutOfRange) {
        rescaleEngine->scaleCoarseFine(
            coarseLattice.get(x0 - 2 * delta.x, iY - 2 * delta.y, iZ - 2 * delta.z), pop1);
        rescaleEngine->scaleCoarseFine(
            coarseLattice.get(x0 - delta.x, iY - delta.y, iZ - delta.z), pop2);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(x0, iY, iZ), pop3);
        rescaleEngine->scaleCoarseFine(
            coarseLattice.get(x0 + delta.x, iY + delta.y, iZ + delta.z), pop4);

        cubicCenteredInterpolation<T>(
            pop1, pop2, pop3, pop4,
            converter.fineDynamics(Dot3D(x0, iY, iZ), Dot3D(-1, 0, 0))
                .getDecomposedValues(whichTime));
    }

    // only variations in the x-coordinate
    for (plint iX = x0; iX <= x1 - delta.x; ++iX) {
        rescaleEngine->scaleCoarseFine(
            coarseLattice.get(iX - delta.x, iY - delta.y, iZ - delta.z), pop1);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop2);
        rescaleEngine->scaleCoarseFine(
            coarseLattice.get(iX + delta.x, iY + delta.y, iZ + delta.z), pop3);
        rescaleEngine->scaleCoarseFine(
            coarseLattice.get(iX + 2 * delta.x, iY + 2 * delta.y, iZ + 2 * delta.z), pop4);

        // cubicInterpolation the unknown population
        cubicCenteredInterpolation<T>(
            pop1, pop2, pop3, pop4,
            converter.fineDynamics(Dot3D(iX, iY, iZ), delta).getDecomposedValues(whichTime));
    }

    // check if interpolation is needed in the rightmost part
    if (!endOutOfRange) {
        rescaleEngine->scaleCoarseFine(
            coarseLattice.get(x1 - delta.x, iY - delta.y, iZ - delta.z), pop1);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(x1, iY, iZ), pop2);
        rescaleEngine->scaleCoarseFine(
            coarseLattice.get(x1 + delta.x, iY + delta.y, iZ + delta.z), pop3);
        rescaleEngine->scaleCoarseFine(
            coarseLattice.get(x1 + 2 * delta.x, iY + 2 * delta.y, iZ + 2 * delta.z), pop4);

        cubicCenteredInterpolation<T>(
            pop1, pop2, pop3, pop4,
            converter.fineDynamics(Dot3D(x1, iY, iZ), delta).getDecomposedValues(whichTime));
    }
}

template <typename T, template <typename U> class Descriptor>
void cubicInterpolationOverFineLineX(Box3D domain, BlockLattice3D<T, Descriptor> &fineLattice)
{
    // attention, everything is in fine coordinates
    // however, the first site always corresponds to a coarse site!
    PLB_ASSERT(domain.y0 == domain.y1 && domain.z0 == domain.z1);
    plint iY = domain.y0;
    plint iZ = domain.z0;
    plint x0 = domain.x0;
    plint x1 = domain.x1;

    Dot3D delta(1, 0, 0);
    plint whichTime = 1;
    std::vector<T> pop1;
    std::vector<T> pop2;
    std::vector<T> pop3;
    std::vector<T> pop4;

    Dot3D startFineExcess(x0 - 4 * delta.x, iY - 4 * delta.y, iZ - 4 * delta.z);
    bool startOutOfRange = !contained(startFineExcess, fineLattice.getBoundingBox());

    Dot3D endFineExcess(x1 + 4 * delta.x, iY + 4 * delta.y, iZ + 4 * delta.z);
    bool endOutOfRange = !contained(endFineExcess, fineLattice.getBoundingBox());

    // check if interpolation is needed in the leftmost part
    if (!startOutOfRange) {
        //         pop1 = dynamic_cast<FineGridBoundaryDynamics<T,Descriptor>&> (
        //                    fineLattice.get(x0-4*delta.x,iY-4*delta.y,iZ-4*delta.z).getDynamics()
        //                    ).getDecomposedValues(whichTime);
        //         pop2 = dynamic_cast<FineGridBoundaryDynamics<T,Descriptor>&> (
        //                    fineLattice.get(x0-2*delta.x,iY-2*delta.y,iZ-2*delta.z).getDynamics()
        //                    ).getDecomposedValues(whichTime);
        //         pop3 = dynamic_cast<FineGridBoundaryDynamics<T,Descriptor>&> (
        //                    fineLattice.get(x0,iY,iZ).getDynamics()
        //                    ).getDecomposedValues(whichTime);
        //         pop4 = dynamic_cast<FineGridBoundaryDynamics<T,Descriptor>&> (
        //                    fineLattice.get(x0+2*delta.x,iY+2*delta.y,iZ+2*delta.z).getDynamics()
        //                    ).getDecomposedValues(whichTime);
        //
        //         cubicCenteredInterpolation<T>( pop1,pop2,pop3,pop4,
        //                         dynamic_cast<FineGridBoundaryDynamics<T,Descriptor>&> (
        //                             fineLattice.get(x0-delta.x,iY-delta.y,iZ-delta.z).getDynamics()
        //                                                     ).getDecomposedValues(whichTime) );
    }

    // only variations in the x-coordinate
    for (plint iX = x0; iX <= x1 - 2 * delta.x; iX += 2) {
        pop1 =
            dynamic_cast<FineGridBoundaryDynamics<T, Descriptor> &>(
                fineLattice.get(iX - 2 * delta.x, iY - 2 * delta.y, iZ - 2 * delta.z).getDynamics())
                .getDecomposedValues(whichTime);
        pop2 = dynamic_cast<FineGridBoundaryDynamics<T, Descriptor> &>(
                   fineLattice.get(iX, iY, iZ).getDynamics())
                   .getDecomposedValues(whichTime);
        pop3 =
            dynamic_cast<FineGridBoundaryDynamics<T, Descriptor> &>(
                fineLattice.get(iX + 2 * delta.x, iY + 2 * delta.y, iZ + 2 * delta.z).getDynamics())
                .getDecomposedValues(whichTime);
        pop4 =
            dynamic_cast<FineGridBoundaryDynamics<T, Descriptor> &>(
                fineLattice.get(iX + 4 * delta.x, iY + 4 * delta.y, iZ + 4 * delta.z).getDynamics())
                .getDecomposedValues(whichTime);

        cubicCenteredInterpolation<T>(
            pop1, pop2, pop3, pop4,
            dynamic_cast<FineGridBoundaryDynamics<T, Descriptor> &>(
                fineLattice.get(iX + delta.x, iY + delta.y, iZ + delta.z).getDynamics())
                .getDecomposedValues(whichTime));
    }

    // check if interpolation is needed in the rightmost part
    if (!endOutOfRange) {
        //         pop1 = dynamic_cast<FineGridBoundaryDynamics<T,Descriptor>&> (
        //                    fineLattice.get(x1-2*delta.x,iY-2*delta.y,iZ-2*delta.z).getDynamics()
        //                    ).getDecomposedValues(whichTime);
        //         pop2 = dynamic_cast<FineGridBoundaryDynamics<T,Descriptor>&> (
        //                    fineLattice.get(x1,iY,iZ).getDynamics()
        //                    ).getDecomposedValues(whichTime);
        //         pop3 = dynamic_cast<FineGridBoundaryDynamics<T,Descriptor>&> (
        //                    fineLattice.get(x1+2*delta.x,iY+2*delta.y,iZ+2*delta.z).getDynamics()
        //                    ).getDecomposedValues(whichTime);
        //         pop4 = dynamic_cast<FineGridBoundaryDynamics<T,Descriptor>&> (
        //                    fineLattice.get(x1+4*delta.x,iY+4*delta.y,iZ+4*delta.z).getDynamics()
        //                    ).getDecomposedValues(whichTime);
        //
        //         cubicCenteredInterpolation<T>( pop1,pop2,pop3,pop4,
        //                         dynamic_cast<FineGridBoundaryDynamics<T,Descriptor>&> (
        //                             fineLattice.get(x1+delta.x,iY+delta.y,iZ+delta.z).getDynamics()
        //                                                     ).getDecomposedValues(whichTime) );
    }
}

template <typename T, template <typename U> class Descriptor>
void cubicInterpolationForEdgeX(
    Box3D domain, BlockLattice3D<T, Descriptor> &coarseLattice,
    BlockLattice3D<T, Descriptor> &fineLattice, RescaleEngine<T, Descriptor> *rescaleEngine,
    plint orientation)
{
    // domain is in coarse coordinates
    PLB_ASSERT(domain.y0 == domain.y1 && domain.z0 == domain.z1);
    plint iY = domain.y0;
    plint iZ = domain.z0;
    plint x0 = domain.x0;
    plint x1 = domain.x1;

    pcout << "Perpendicular X with orientation " << orientation << std::endl;
    CoarseToFineConverter3D<T, Descriptor> converter(coarseLattice, fineLattice);

    Box3D line(
        converter.fineX(x0), converter.fineX(x1), converter.fineY(iY) + orientation,
        converter.fineY(iY) + orientation, converter.fineZ(iZ), converter.fineZ(iZ));

    plint whichTime = 1;
    for (plint iX = line.x0; iX <= line.x1; ++iX) {
        // non centered interpolation over the edges
    }
}

template <typename T, template <typename U> class Descriptor>
void cubicInterpolationOverCoarseLineY(
    Box3D domain, BlockLattice3D<T, Descriptor> &coarseLattice,
    BlockLattice3D<T, Descriptor> &fineLattice, RescaleEngine<T, Descriptor> *rescaleEngine)
{
    // domain is given in coarse coordinates
    PLB_ASSERT(domain.x0 == domain.x1 && domain.z0 == domain.z1);
    plint iX = domain.x0;
    plint iZ = domain.z0;
    plint y0 = domain.y0;
    plint y1 = domain.y1;

    CoarseToFineConverter3D<T, Descriptor> converter(coarseLattice, fineLattice);
    Dot3D delta(0, 1, 0);
    plint whichTime = 1;
    std::vector<T> pop1;
    std::vector<T> pop2;
    std::vector<T> pop3;
    std::vector<T> pop4;

    Dot3D startExcess(iX - 2 * delta.x, y0 - 2 * delta.y, iZ - delta.z);
    Dot3D startFineExcess(
        converter.fineX(iX) - 2 * delta.x, converter.fineY(y0) - 2 * delta.y,
        converter.fineZ(iZ) - 2 * delta.z);
    bool startOutOfRange =
        (!contained(startExcess, coarseLattice.getBoundingBox())
         || !contained(startFineExcess, fineLattice.getBoundingBox()));

    Dot3D endExcess(iX + 2 * delta.x, y1 + 2 * delta.y, iZ + delta.z);
    Dot3D endFineExcess(
        converter.fineX(iX) + 2 * delta.x, converter.fineY(y1) + 2 * delta.y,
        converter.fineZ(iZ) + 2 * delta.z);
    bool endOutOfRange =
        (!contained(endExcess, coarseLattice.getBoundingBox())
         || !contained(endFineExcess, fineLattice.getBoundingBox()));

    // check if interpolation is needed in the leftmost part
    if (!startOutOfRange) {
        rescaleEngine->scaleCoarseFine(
            coarseLattice.get(iX - 2 * delta.x, y0 - 2 * delta.y, iZ - 2 * delta.z), pop1);
        rescaleEngine->scaleCoarseFine(
            coarseLattice.get(iX - delta.x, y0 - delta.y, iZ - delta.z), pop2);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, y0, iZ), pop3);
        rescaleEngine->scaleCoarseFine(
            coarseLattice.get(iX + delta.x, y0 + delta.y, iZ + delta.z), pop4);

        cubicCenteredInterpolation<T>(
            pop1, pop2, pop3, pop4,
            converter.fineDynamics(Dot3D(iX, y0, iZ), Dot3D(0, -1, 0))
                .getDecomposedValues(whichTime));
    }

    // only variations in the x-coordinate
    for (plint iY = y0; iY <= y1 - delta.y; ++iY) {
        rescaleEngine->scaleCoarseFine(
            coarseLattice.get(iX - delta.x, iY - delta.y, iZ - delta.z), pop1);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, iY, iZ), pop2);
        rescaleEngine->scaleCoarseFine(
            coarseLattice.get(iX + delta.x, iY + delta.y, iZ + delta.z), pop3);
        rescaleEngine->scaleCoarseFine(
            coarseLattice.get(iX + 2 * delta.x, iY + 2 * delta.y, iZ + 2 * delta.z), pop4);

        // cubicInterpolation the unknown population
        cubicCenteredInterpolation<T>(
            pop1, pop2, pop3, pop4,
            converter.fineDynamics(Dot3D(iX, iY, iZ), delta).getDecomposedValues(whichTime));
    }

    // check if interpolation is needed in the rightmost part
    if (!endOutOfRange) {
        rescaleEngine->scaleCoarseFine(
            coarseLattice.get(iX - delta.x, y1 - delta.y, iZ - delta.z), pop1);
        rescaleEngine->scaleCoarseFine(coarseLattice.get(iX, y1, iZ), pop2);
        rescaleEngine->scaleCoarseFine(
            coarseLattice.get(iX + delta.x, y1 + delta.y, iZ + delta.z), pop3);
        rescaleEngine->scaleCoarseFine(
            coarseLattice.get(iX + 2 * delta.x, y1 + 2 * delta.y, iZ + 2 * delta.z), pop4);

        cubicCenteredInterpolation<T>(
            pop1, pop2, pop3, pop4,
            converter.fineDynamics(Dot3D(iX, y1, iZ), delta).getDecomposedValues(whichTime));
    }
}

template <typename T, template <typename U> class Descriptor>
void cubicInterpolationOverFineLineY(
    Box3D domain, BlockLattice3D<T, Descriptor> &coarseLattice,
    BlockLattice3D<T, Descriptor> &fineLattice, RescaleEngine<T, Descriptor> *rescaleEngine)
{
    // domain is given in fine coordinates
    PLB_ASSERT(domain.x0 == domain.x1 && domain.z0 == domain.z1);
    plint iX = domain.x0;
    plint iZ = domain.z0;
    plint y0 = domain.y0;
    plint y1 = domain.y1;

    CoarseToFineConverter3D<T, Descriptor> converter(coarseLattice, fineLattice);

    // compute the coarse coordinates
    Dot3D posFine = fineLattice.getLocation();      // Position of fine grid in multi block.
    Dot3D posCoarse = coarseLattice.getLocation();  // Position of coarse grid in multi block.

    Box3D coarseDomain(
        domain.shift(posFine.x, posFine.y, posFine.z)
            .  // Convert to absolute fine coordinates.
        divideAndFitSmaller(2)
            .  // Rescale, but don't exceed original domain.
        shift(
            -posCoarse.x, -posCoarse.y, -posCoarse.z));  // Convert to relative coarse coordinates.

    plint coarseY0 = coarseDomain.y0;
    plint coarseY1 = coarseDomain.y1;
}

template <typename T, template <typename U> class Descriptor>
void cubicInterpolationForEdgeY(
    Box3D domain, BlockLattice3D<T, Descriptor> &coarseLattice,
    BlockLattice3D<T, Descriptor> &fineLattice, RescaleEngine<T, Descriptor> *rescaleEngine,
    plint orientation)
{ }

template <typename T, template <typename U> class Descriptor>
void cubicInterpolationOverCoarseLineZ(
    Box3D domain, BlockLattice3D<T, Descriptor> &coarseLattice,
    BlockLattice3D<T, Descriptor> &fineLattice, RescaleEngine<T, Descriptor> *rescaleEngine)
{ }

template <typename T, template <typename U> class Descriptor>
void cubicInterpolationOverFineLineZ(
    Box3D domain, BlockLattice3D<T, Descriptor> &coarseLattice,
    BlockLattice3D<T, Descriptor> &fineLattice, RescaleEngine<T, Descriptor> *rescaleEngine)
{ }

}  // namespace plb

#endif  // LINE_INTERPOLATION_3D_HH
