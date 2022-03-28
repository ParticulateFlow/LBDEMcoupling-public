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
 * Dynamics and data processors added to the coarse grid to implement 3D grid refinement --
 * implementation file.
 */

#ifndef COARSE_GRID_PROCESSORS_3D_HH
#define COARSE_GRID_PROCESSORS_3D_HH

#include "multiGrid/coarseGridProcessors3D.h"

namespace plb {

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

                    // if dynamics = advection diffusion dynamics, simply copy
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

/* ************* CopyFineToCoarseWithFiltering3D *************** */

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyFineToCoarseWithFiltering3D<T, Descriptor1, Descriptor2>::CopyFineToCoarseWithFiltering3D(
    RescaleEngine<T, Descriptor1> *rescaleEngine_, plint numTimeSteps_, plint executionTime_,
    std::vector<plint> indices_) :
    rescaleEngine(rescaleEngine_),
    numTimeSteps(numTimeSteps_),
    executionTime(executionTime_),
    indices(indices_)
{ }

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyFineToCoarseWithFiltering3D<T, Descriptor1, Descriptor2>::~CopyFineToCoarseWithFiltering3D()
{
    delete rescaleEngine;
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyFineToCoarseWithFiltering3D<T, Descriptor1, Descriptor2>::CopyFineToCoarseWithFiltering3D(
    CopyFineToCoarseWithFiltering3D<T, Descriptor1, Descriptor2> const &rhs) :
    rescaleEngine(rhs.rescaleEngine->clone()),
    numTimeSteps(rhs.numTimeSteps),
    executionTime(rhs.executionTime),
    indices(rhs.indices)
{ }

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyFineToCoarseWithFiltering3D<T, Descriptor1, Descriptor2>
    &CopyFineToCoarseWithFiltering3D<T, Descriptor1, Descriptor2>::operator=(
        CopyFineToCoarseWithFiltering3D<T, Descriptor1, Descriptor2> const &rhs)
{
    delete rescaleEngine;
    rescaleEngine = rhs.rescaleEngine->clone();
    numTimeSteps = rhs.numTimeSteps;
    executionTime = rhs.executionTime;
    indices.resize(0);
    indices.assign(rhs.indices.begin(), rhs.indices.end());
    return *this;
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
void CopyFineToCoarseWithFiltering3D<T, Descriptor1, Descriptor2>::process(
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

        plint start = (plint)1 + Descriptor1<T>::d;
        // Loop over coarse lattice
        for (plint iX = coarseDomain.x0; iX <= coarseDomain.x1; ++iX) {
            for (plint iY = coarseDomain.y0; iY <= coarseDomain.y1; ++iY) {
                for (plint iZ = coarseDomain.z0; iZ <= coarseDomain.z1; ++iZ) {
                    // Determine corresponding coordinates on fine lattice
                    plint fineX = (iX + posCoarse.x) * 2 - posFine.x;
                    plint fineY = (iY + posCoarse.y) * 2 - posFine.y;
                    plint fineZ = (iZ + posCoarse.z) * 2 - posFine.z;

                    Cell<T, Descriptor1> &coarseCell = coarseLattice.get(iX, iY, iZ);
                    Cell<T, Descriptor2> const &fineCell = fineLattice.get(fineX, fineY, fineZ);

                    rescaleEngine->scaleFineCoarse(fineCell, decomposedCoarseValues);

                    for (plint iPop = 0; iPop < (plint)indices.size(); ++iPop) {
                        PLB_ASSERT(indices[iPop] != 0);
                        std::vector<T> tmpDec;
                        Cell<T, Descriptor2> const &nextCell = fineLattice.get(
                            fineX + Descriptor2<T>::c[indices[iPop]][0],
                            fineY + Descriptor2<T>::c[indices[iPop]][1],
                            fineZ + Descriptor2<T>::c[indices[iPop]][2]);
                        rescaleEngine->scaleFineCoarse(nextCell, tmpDec);
                        for (pluint iA = start; iA < tmpDec.size(); ++iA) {
                            decomposedCoarseValues[iA] += tmpDec[iA];
                        }
                    }
                    for (pluint iA = start; iA < decomposedCoarseValues.size(); ++iA) {
                        decomposedCoarseValues[iA] /= (T)(indices.size() + 1);
                    }

                    rescaleEngine->recompose(coarseCell, decomposedCoarseValues);
                }
            }
        }
    }
}

template <
    typename T, template <typename U> class Descriptor1, template <typename U> class Descriptor2>
CopyFineToCoarseWithFiltering3D<T, Descriptor1, Descriptor2>
    *CopyFineToCoarseWithFiltering3D<T, Descriptor1, Descriptor2>::clone() const
{
    return new CopyFineToCoarseWithFiltering3D(*this);
}

}  // namespace plb

#endif  // COARSE_GRID_PROCESSORS_3D_HH
