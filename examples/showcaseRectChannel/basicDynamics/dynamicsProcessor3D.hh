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

#ifndef DYNAMICS_PROCESSOR_3D_HH
#define DYNAMICS_PROCESSOR_3D_HH

#include "atomicBlock/blockLattice3D.h"
#include "basicDynamics/dynamicsProcessor3D.h"
#include "core/cell.h"
#include "core/plbProfiler.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/latticeTemplates.h"
#include "multiGrid/multiGridUtil.h"

namespace plb {

/* ************* Class ExternalRhoJcollideAndStream3D ******************* */

template <typename T, template <typename U> class Descriptor>
void ExternalRhoJcollideAndStream3D<T, Descriptor>::collide(
    BlockLattice3D<T, Descriptor> &lattice, Box3D const &domain,
    ScalarField3D<T> const &rhoBarField, Dot3D const &offset1, TensorField3D<T, 3> const &jField,
    Dot3D const &offset2, BlockStatistics &stat)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
                T rhoBar = rhoBarField.get(iX + offset1.x, iY + offset1.y, iZ + offset1.z);
                Array<T, 3> const &j = jField.get(iX + offset2.x, iY + offset2.y, iZ + offset2.z);
                cell.getDynamics().collideExternal(cell, rhoBar, j, T(), stat);
                cell.revert();
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void ExternalRhoJcollideAndStream3D<T, Descriptor>::bulkCollideAndStream(
    BlockLattice3D<T, Descriptor> &lattice, Box3D const &domain,
    ScalarField3D<T> const &rhoBarField, Dot3D const &offset1, TensorField3D<T, 3> const &jField,
    Dot3D const &offset2, BlockStatistics &stat)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
                T rhoBar = rhoBarField.get(iX + offset1.x, iY + offset1.y, iZ + offset1.z);
                Array<T, 3> const &j = jField.get(iX + offset2.x, iY + offset2.y, iZ + offset2.z);
                cell.getDynamics().collideExternal(cell, rhoBar, j, T(), stat);
                latticeTemplates<T, Descriptor>::swapAndStream3D(lattice.grid, iX, iY, iZ);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void ExternalRhoJcollideAndStream3D<T, Descriptor>::boundaryStream(
    BlockLattice3D<T, Descriptor> &lattice, Box3D const &bound, Box3D const &domain)
{
    // Make sure domain is contained within bound
    PLB_PRECONDITION(contained(domain, bound));

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                for (plint iPop = 1; iPop <= Descriptor<T>::q / 2; ++iPop) {
                    plint nextX = iX + Descriptor<T>::c[iPop][0];
                    plint nextY = iY + Descriptor<T>::c[iPop][1];
                    plint nextZ = iZ + Descriptor<T>::c[iPop][2];
                    if (nextX >= bound.x0 && nextX <= bound.x1 && nextY >= bound.y0
                        && nextY <= bound.y1 && nextZ >= bound.z0 && nextZ <= bound.z1)
                    {
                        std::swap(
                            lattice.grid[iX][iY][iZ][iPop + Descriptor<T>::q / 2],
                            lattice.grid[nextX][nextY][nextZ][iPop]);
                    }
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void ExternalRhoJcollideAndStream3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    global::timer("collideAndStream").start();
    BlockLattice3D<T, Descriptor> &lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> &>(*atomicBlocks[0]);
    ScalarField3D<T> const &rhoBarField = dynamic_cast<ScalarField3D<T> const &>(*atomicBlocks[1]);
    TensorField3D<T, 3> const &jField = dynamic_cast<TensorField3D<T, 3> const &>(*atomicBlocks[2]);

    BlockStatistics &stat = lattice.getInternalStatistics();

    static const plint vicinity = Descriptor<T>::vicinity;
    Box3D extDomain(domain.enlarge(vicinity));

    Dot3D offset1 = computeRelativeDisplacement(lattice, rhoBarField);
    Dot3D offset2 = computeRelativeDisplacement(lattice, jField);

    global::profiler().start("collStream");
    global::profiler().increment("collStreamCells", extDomain.nCells());

    // First, do the collision on cells within a boundary envelope of width
    // equal to the range of the lattice vectors (e.g. 1 for D2Q9)
    collide(
        lattice,
        Box3D(
            extDomain.x0, extDomain.x0 + vicinity - 1, extDomain.y0, extDomain.y1, extDomain.z0,
            extDomain.z1),
        rhoBarField, offset1, jField, offset2, stat);
    collide(
        lattice,
        Box3D(
            extDomain.x1 - vicinity + 1, extDomain.x1, extDomain.y0, extDomain.y1, extDomain.z0,
            extDomain.z1),
        rhoBarField, offset1, jField, offset2, stat);
    collide(
        lattice,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0,
            extDomain.y0 + vicinity - 1, extDomain.z0, extDomain.z1),
        rhoBarField, offset1, jField, offset2, stat);
    collide(
        lattice,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y1 - vicinity + 1,
            extDomain.y1, extDomain.z0, extDomain.z1),
        rhoBarField, offset1, jField, offset2, stat);
    collide(
        lattice,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0 + vicinity,
            extDomain.y1 - vicinity, extDomain.z0, extDomain.z0 + vicinity - 1),
        rhoBarField, offset1, jField, offset2, stat);
    collide(
        lattice,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0 + vicinity,
            extDomain.y1 - vicinity, extDomain.z1 - vicinity + 1, extDomain.z1),
        rhoBarField, offset1, jField, offset2, stat);

    // Then, do the efficient collideAndStream algorithm in the bulk,
    // excluding the envelope (this is efficient because there is no
    // if-then-else statement within the loop, given that the boundary
    // region is excluded)
    bulkCollideAndStream(
        lattice,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0 + vicinity,
            extDomain.y1 - vicinity, extDomain.z0 + vicinity, extDomain.z1 - vicinity),
        rhoBarField, offset1, jField, offset2, stat);

    // Finally, do streaming in the boundary envelope to conclude the
    // collision-stream cycle
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x0, extDomain.x0 + vicinity - 1, extDomain.y0, extDomain.y1, extDomain.z0,
            extDomain.z1));
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x1 - vicinity + 1, extDomain.x1, extDomain.y0, extDomain.y1, extDomain.z0,
            extDomain.z1));
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0,
            extDomain.y0 + vicinity - 1, extDomain.z0, extDomain.z1));
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y1 - vicinity + 1,
            extDomain.y1, extDomain.z0, extDomain.z1));
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0 + vicinity,
            extDomain.y1 - vicinity, extDomain.z0, extDomain.z0 + vicinity - 1));
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0 + vicinity,
            extDomain.y1 - vicinity, extDomain.z1 - vicinity + 1, extDomain.z1));
    global::profiler().stop("collStream");
    global::timer("collideAndStream").stop();
}

template <typename T, template <typename U> class Descriptor>
ExternalRhoJcollideAndStream3D<T, Descriptor>
    *ExternalRhoJcollideAndStream3D<T, Descriptor>::clone() const
{
    return new ExternalRhoJcollideAndStream3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void ExternalRhoJcollideAndStream3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

/* ************* Class PackedExternalRhoJcollideAndStream3D ******************* */

template <typename T, template <typename U> class Descriptor>
void PackedExternalRhoJcollideAndStream3D<T, Descriptor>::collide(
    BlockLattice3D<T, Descriptor> &lattice, Box3D const &domain,
    NTensorField3D<T> const &rhoBarJfield, Dot3D const &offset, BlockStatistics &stat)
{
    Array<T, 3> j;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
                T const *macroscopic =
                    rhoBarJfield.get(iX + offset.x, iY + offset.y, iZ + offset.z);
                T rhoBar = *macroscopic;
                j.from_cArray(macroscopic + 1);

                cell.getDynamics().collideExternal(cell, rhoBar, j, T(), stat);
                cell.revert();
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void PackedExternalRhoJcollideAndStream3D<T, Descriptor>::bulkCollideAndStream(
    BlockLattice3D<T, Descriptor> &lattice, Box3D const &domain,
    NTensorField3D<T> const &rhoBarJfield, Dot3D const &offset, BlockStatistics &stat)
{
    Array<T, 3> j;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
                T const *macroscopic =
                    rhoBarJfield.get(iX + offset.x, iY + offset.y, iZ + offset.z);
                T rhoBar = *macroscopic;
                j.from_cArray(macroscopic + 1);

                cell.getDynamics().collideExternal(cell, rhoBar, j, T(), stat);
                latticeTemplates<T, Descriptor>::swapAndStream3D(lattice.grid, iX, iY, iZ);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void PackedExternalRhoJcollideAndStream3D<T, Descriptor>::boundaryStream(
    BlockLattice3D<T, Descriptor> &lattice, Box3D const &bound, Box3D const &domain)
{
    // Make sure domain is contained within bound
    PLB_PRECONDITION(contained(domain, bound));

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                for (plint iPop = 1; iPop <= Descriptor<T>::q / 2; ++iPop) {
                    plint nextX = iX + Descriptor<T>::c[iPop][0];
                    plint nextY = iY + Descriptor<T>::c[iPop][1];
                    plint nextZ = iZ + Descriptor<T>::c[iPop][2];
                    if (nextX >= bound.x0 && nextX <= bound.x1 && nextY >= bound.y0
                        && nextY <= bound.y1 && nextZ >= bound.z0 && nextZ <= bound.z1)
                    {
                        std::swap(
                            lattice.grid[iX][iY][iZ][iPop + Descriptor<T>::q / 2],
                            lattice.grid[nextX][nextY][nextZ][iPop]);
                    }
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void PackedExternalRhoJcollideAndStream3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, NTensorField3D<T> &rhoBarJfield)
{
    global::timer("collideAndStream").start();

    PLB_ASSERT(rhoBarJfield.getNdim() == 4);

    BlockStatistics &stat = lattice.getInternalStatistics();

    static const plint vicinity = Descriptor<T>::vicinity;
    Box3D extDomain(domain.enlarge(vicinity));

    Dot3D offset = computeRelativeDisplacement(lattice, rhoBarJfield);

    global::profiler().start("collStream");
    global::profiler().increment("collStreamCells", extDomain.nCells());

    // First, do the collision on cells within a boundary envelope of width
    // equal to the range of the lattice vectors (e.g. 1 for D2Q9)
    collide(
        lattice,
        Box3D(
            extDomain.x0, extDomain.x0 + vicinity - 1, extDomain.y0, extDomain.y1, extDomain.z0,
            extDomain.z1),
        rhoBarJfield, offset, stat);
    collide(
        lattice,
        Box3D(
            extDomain.x1 - vicinity + 1, extDomain.x1, extDomain.y0, extDomain.y1, extDomain.z0,
            extDomain.z1),
        rhoBarJfield, offset, stat);
    collide(
        lattice,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0,
            extDomain.y0 + vicinity - 1, extDomain.z0, extDomain.z1),
        rhoBarJfield, offset, stat);
    collide(
        lattice,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y1 - vicinity + 1,
            extDomain.y1, extDomain.z0, extDomain.z1),
        rhoBarJfield, offset, stat);
    collide(
        lattice,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0 + vicinity,
            extDomain.y1 - vicinity, extDomain.z0, extDomain.z0 + vicinity - 1),
        rhoBarJfield, offset, stat);
    collide(
        lattice,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0 + vicinity,
            extDomain.y1 - vicinity, extDomain.z1 - vicinity + 1, extDomain.z1),
        rhoBarJfield, offset, stat);

    // Then, do the efficient collideAndStream algorithm in the bulk,
    // excluding the envelope (this is efficient because there is no
    // if-then-else statement within the loop, given that the boundary
    // region is excluded)
    bulkCollideAndStream(
        lattice,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0 + vicinity,
            extDomain.y1 - vicinity, extDomain.z0 + vicinity, extDomain.z1 - vicinity),
        rhoBarJfield, offset, stat);

    // Finally, do streaming in the boundary envelope to conclude the
    // collision-stream cycle
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x0, extDomain.x0 + vicinity - 1, extDomain.y0, extDomain.y1, extDomain.z0,
            extDomain.z1));
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x1 - vicinity + 1, extDomain.x1, extDomain.y0, extDomain.y1, extDomain.z0,
            extDomain.z1));
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0,
            extDomain.y0 + vicinity - 1, extDomain.z0, extDomain.z1));
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y1 - vicinity + 1,
            extDomain.y1, extDomain.z0, extDomain.z1));
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0 + vicinity,
            extDomain.y1 - vicinity, extDomain.z0, extDomain.z0 + vicinity - 1));
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0 + vicinity,
            extDomain.y1 - vicinity, extDomain.z1 - vicinity + 1, extDomain.z1));
    global::profiler().stop("collStream");
    global::timer("collideAndStream").stop();
}

template <typename T, template <typename U> class Descriptor>
PackedExternalRhoJcollideAndStream3D<T, Descriptor>
    *PackedExternalRhoJcollideAndStream3D<T, Descriptor>::clone() const
{
    return new PackedExternalRhoJcollideAndStream3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void PackedExternalRhoJcollideAndStream3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

/* ************* Class WaveAbsorptionExternalRhoJcollideAndStream3D ******************* */

template <typename T, template <typename U> class Descriptor>
WaveAbsorptionExternalRhoJcollideAndStream3D<T, Descriptor>::
    WaveAbsorptionExternalRhoJcollideAndStream3D(
        T rhoBarF_, Array<T, Descriptor<T>::d> const &uF_) :
    rhoBarF(rhoBarF_)
{
    jF = Descriptor<T>::fullRho(rhoBarF) * uF_;
    jFsqr = normSqr(jF);
}

template <typename T, template <typename U> class Descriptor>
void WaveAbsorptionExternalRhoJcollideAndStream3D<T, Descriptor>::collide(
    BlockLattice3D<T, Descriptor> &lattice, Box3D const &domain,
    ScalarField3D<T> const &rhoBarField, Dot3D const &offset1,
    TensorField3D<T, Descriptor<T>::d> const &jField, Dot3D const &offset2,
    ScalarField3D<T> const &sigmaField, Dot3D const &offset3, BlockStatistics &stat)
{
    static const T epsilon = 1.e4 * std::numeric_limits<T>::epsilon();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
                T rhoBar = rhoBarField.get(iX + offset1.x, iY + offset1.y, iZ + offset1.z);
                Array<T, Descriptor<T>::d> const &j =
                    jField.get(iX + offset2.x, iY + offset2.y, iZ + offset2.z);
                T sigma = sigmaField.get(iX + offset3.x, iY + offset3.y, iZ + offset3.z);

                if (sigma < epsilon) {
                    cell.getDynamics().collideExternal(cell, rhoBar, j, T(), stat);
                } else {
                    T jSqr = normSqr(j);
                    Array<T, Descriptor<T>::q> fEq;
                    cell.getDynamics().computeEquilibria(fEq, rhoBar, j, jSqr);

                    Array<T, Descriptor<T>::q> fEqF;
                    cell.getDynamics().computeEquilibria(fEqF, rhoBarF, jF, jFsqr);

                    cell.getDynamics().collideExternal(cell, rhoBar, j, T(), stat);

                    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                        cell[iPop] -= sigma * (fEq[iPop] - fEqF[iPop]);
                    }
                }

                cell.revert();
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void WaveAbsorptionExternalRhoJcollideAndStream3D<T, Descriptor>::bulkCollideAndStream(
    BlockLattice3D<T, Descriptor> &lattice, Box3D const &domain,
    ScalarField3D<T> const &rhoBarField, Dot3D const &offset1,
    TensorField3D<T, Descriptor<T>::d> const &jField, Dot3D const &offset2,
    ScalarField3D<T> const &sigmaField, Dot3D const &offset3, BlockStatistics &stat)
{
    static const T epsilon = 1.e4 * std::numeric_limits<T>::epsilon();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
                T rhoBar = rhoBarField.get(iX + offset1.x, iY + offset1.y, iZ + offset1.z);
                Array<T, Descriptor<T>::d> const &j =
                    jField.get(iX + offset2.x, iY + offset2.y, iZ + offset2.z);
                T sigma = sigmaField.get(iX + offset3.x, iY + offset3.y, iZ + offset3.z);

                if (sigma < epsilon) {
                    cell.getDynamics().collideExternal(cell, rhoBar, j, T(), stat);
                } else {
                    T jSqr = normSqr(j);
                    Array<T, Descriptor<T>::q> fEq;
                    cell.getDynamics().computeEquilibria(fEq, rhoBar, j, jSqr);

                    Array<T, Descriptor<T>::q> fEqF;
                    cell.getDynamics().computeEquilibria(fEqF, rhoBarF, jF, jFsqr);

                    cell.getDynamics().collideExternal(cell, rhoBar, j, T(), stat);

                    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                        cell[iPop] -= sigma * (fEq[iPop] - fEqF[iPop]);
                    }
                }

                latticeTemplates<T, Descriptor>::swapAndStream3D(lattice.grid, iX, iY, iZ);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void WaveAbsorptionExternalRhoJcollideAndStream3D<T, Descriptor>::boundaryStream(
    BlockLattice3D<T, Descriptor> &lattice, Box3D const &bound, Box3D const &domain)
{
    // Make sure domain is contained within bound
    PLB_PRECONDITION(contained(domain, bound));

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                for (plint iPop = 1; iPop <= Descriptor<T>::q / 2; ++iPop) {
                    plint nextX = iX + Descriptor<T>::c[iPop][0];
                    plint nextY = iY + Descriptor<T>::c[iPop][1];
                    plint nextZ = iZ + Descriptor<T>::c[iPop][2];
                    if (nextX >= bound.x0 && nextX <= bound.x1 && nextY >= bound.y0
                        && nextY <= bound.y1 && nextZ >= bound.z0 && nextZ <= bound.z1)
                    {
                        std::swap(
                            lattice.grid[iX][iY][iZ][iPop + Descriptor<T>::q / 2],
                            lattice.grid[nextX][nextY][nextZ][iPop]);
                    }
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void WaveAbsorptionExternalRhoJcollideAndStream3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    global::timer("collideAndStream").start();
    BlockLattice3D<T, Descriptor> &lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> &>(*atomicBlocks[0]);
    ScalarField3D<T> const &rhoBarField = dynamic_cast<ScalarField3D<T> const &>(*atomicBlocks[1]);
    TensorField3D<T, Descriptor<T>::d> const &jField =
        dynamic_cast<TensorField3D<T, Descriptor<T>::d> const &>(*atomicBlocks[2]);
    ScalarField3D<T> const &sigmaField = dynamic_cast<ScalarField3D<T> const &>(*atomicBlocks[3]);

    BlockStatistics &stat = lattice.getInternalStatistics();

    static const plint vicinity = Descriptor<T>::vicinity;
    Box3D extDomain(domain.enlarge(vicinity));

    Dot3D offset1 = computeRelativeDisplacement(lattice, rhoBarField);
    Dot3D offset2 = computeRelativeDisplacement(lattice, jField);
    Dot3D offset3 = computeRelativeDisplacement(lattice, sigmaField);

    global::profiler().start("collStream");
    global::profiler().increment("collStreamCells", extDomain.nCells());

    // First, do the collision on cells within a boundary envelope of width
    // equal to the range of the lattice vectors (e.g. 1 for D2Q9)
    collide(
        lattice,
        Box3D(
            extDomain.x0, extDomain.x0 + vicinity - 1, extDomain.y0, extDomain.y1, extDomain.z0,
            extDomain.z1),
        rhoBarField, offset1, jField, offset2, sigmaField, offset3, stat);
    collide(
        lattice,
        Box3D(
            extDomain.x1 - vicinity + 1, extDomain.x1, extDomain.y0, extDomain.y1, extDomain.z0,
            extDomain.z1),
        rhoBarField, offset1, jField, offset2, sigmaField, offset3, stat);
    collide(
        lattice,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0,
            extDomain.y0 + vicinity - 1, extDomain.z0, extDomain.z1),
        rhoBarField, offset1, jField, offset2, sigmaField, offset3, stat);
    collide(
        lattice,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y1 - vicinity + 1,
            extDomain.y1, extDomain.z0, extDomain.z1),
        rhoBarField, offset1, jField, offset2, sigmaField, offset3, stat);
    collide(
        lattice,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0 + vicinity,
            extDomain.y1 - vicinity, extDomain.z0, extDomain.z0 + vicinity - 1),
        rhoBarField, offset1, jField, offset2, sigmaField, offset3, stat);
    collide(
        lattice,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0 + vicinity,
            extDomain.y1 - vicinity, extDomain.z1 - vicinity + 1, extDomain.z1),
        rhoBarField, offset1, jField, offset2, sigmaField, offset3, stat);

    // Then, do the efficient collideAndStream algorithm in the bulk,
    // excluding the envelope (this is efficient because there is no
    // if-then-else statement within the loop, given that the boundary
    // region is excluded)
    bulkCollideAndStream(
        lattice,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0 + vicinity,
            extDomain.y1 - vicinity, extDomain.z0 + vicinity, extDomain.z1 - vicinity),
        rhoBarField, offset1, jField, offset2, sigmaField, offset3, stat);

    // Finally, do streaming in the boundary envelope to conclude the
    // collision-stream cycle
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x0, extDomain.x0 + vicinity - 1, extDomain.y0, extDomain.y1, extDomain.z0,
            extDomain.z1));
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x1 - vicinity + 1, extDomain.x1, extDomain.y0, extDomain.y1, extDomain.z0,
            extDomain.z1));
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0,
            extDomain.y0 + vicinity - 1, extDomain.z0, extDomain.z1));
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y1 - vicinity + 1,
            extDomain.y1, extDomain.z0, extDomain.z1));
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0 + vicinity,
            extDomain.y1 - vicinity, extDomain.z0, extDomain.z0 + vicinity - 1));
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0 + vicinity,
            extDomain.y1 - vicinity, extDomain.z1 - vicinity + 1, extDomain.z1));
    global::profiler().stop("collStream");
    global::timer("collideAndStream").stop();
}

template <typename T, template <typename U> class Descriptor>
WaveAbsorptionExternalRhoJcollideAndStream3D<T, Descriptor>
    *WaveAbsorptionExternalRhoJcollideAndStream3D<T, Descriptor>::clone() const
{
    return new WaveAbsorptionExternalRhoJcollideAndStream3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void WaveAbsorptionExternalRhoJcollideAndStream3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
    modified[3] = modif::nothing;
}

/* ************* Class OnLinkExternalRhoJcollideAndStream3D ******************* */

template <typename T, template <typename U> class Descriptor>
void OnLinkExternalRhoJcollideAndStream3D<T, Descriptor>::collide(
    BlockLattice3D<T, Descriptor> &lattice, Box3D const &domain,
    ScalarField3D<T> const &rhoBarField, Dot3D const &offset1, TensorField3D<T, 3> const &jField,
    Dot3D const &offset2, BlockStatistics &stat)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
                T rhoBar = rhoBarField.get(iX + offset1.x, iY + offset1.y, iZ + offset1.z);
                Array<T, 3> const &j = jField.get(iX + offset2.x, iY + offset2.y, iZ + offset2.z);
                cell.getDynamics().collideExternal(cell, rhoBar, j, T(), stat);
                cell.revert();
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void onLinkSwapAndStream3D(Cell<T, Descriptor> ***grid, plint iX, plint iY, plint iZ)
{
    static int bbId = BounceBack<T, Descriptor>().getId();
    const plint half = Descriptor<T>::q / 2;
    for (plint iPop = 1; iPop <= half; ++iPop) {
        plint nextX = iX + Descriptor<T>::c[iPop][0];
        plint nextY = iY + Descriptor<T>::c[iPop][1];
        plint nextZ = iZ + Descriptor<T>::c[iPop][2];
        if (grid[iX][iY][iZ].getDynamics().getId() == bbId
            || grid[nextX][nextY][nextZ].getDynamics().getId() == bbId)
        {
            std::swap(grid[iX][iY][iZ][iPop], grid[iX][iY][iZ][iPop + half]);
        } else {
            T fTmp = grid[iX][iY][iZ][iPop];
            grid[iX][iY][iZ][iPop] = grid[iX][iY][iZ][iPop + half];
            grid[iX][iY][iZ][iPop + half] = grid[nextX][nextY][nextZ][iPop];
            grid[nextX][nextY][nextZ][iPop] = fTmp;
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void OnLinkExternalRhoJcollideAndStream3D<T, Descriptor>::bulkCollideAndStream(
    BlockLattice3D<T, Descriptor> &lattice, Box3D const &domain,
    ScalarField3D<T> const &rhoBarField, Dot3D const &offset1, TensorField3D<T, 3> const &jField,
    Dot3D const &offset2, BlockStatistics &stat)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
                T rhoBar = rhoBarField.get(iX + offset1.x, iY + offset1.y, iZ + offset1.z);
                Array<T, 3> const &j = jField.get(iX + offset2.x, iY + offset2.y, iZ + offset2.z);
                cell.getDynamics().collideExternal(cell, rhoBar, j, T(), stat);
                onLinkSwapAndStream3D<T, Descriptor>(lattice.grid, iX, iY, iZ);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void OnLinkExternalRhoJcollideAndStream3D<T, Descriptor>::boundaryStream(
    BlockLattice3D<T, Descriptor> &lattice, Box3D const &bound, Box3D const &domain)
{
    // Make sure domain is contained within bound
    PLB_PRECONDITION(contained(domain, bound));

    int bbId = BounceBack<T, Descriptor>().getId();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                for (plint iPop = 1; iPop <= Descriptor<T>::q / 2; ++iPop) {
                    plint nextX = iX + Descriptor<T>::c[iPop][0];
                    plint nextY = iY + Descriptor<T>::c[iPop][1];
                    plint nextZ = iZ + Descriptor<T>::c[iPop][2];
                    if (nextX >= bound.x0 && nextX <= bound.x1 && nextY >= bound.y0
                        && nextY <= bound.y1 && nextZ >= bound.z0 && nextZ <= bound.z1)
                    {
                        if (lattice.grid[iX][iY][iZ].getDynamics().getId() != bbId
                            && lattice.grid[nextX][nextY][nextZ].getDynamics().getId() != bbId)
                        {
                            std::swap(
                                lattice.grid[iX][iY][iZ][iPop + Descriptor<T>::q / 2],
                                lattice.grid[nextX][nextY][nextZ][iPop]);
                        }
                    }
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void OnLinkExternalRhoJcollideAndStream3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    global::timer("collideAndStream").start();
    BlockLattice3D<T, Descriptor> &lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> &>(*atomicBlocks[0]);
    ScalarField3D<T> const &rhoBarField = dynamic_cast<ScalarField3D<T> const &>(*atomicBlocks[1]);
    TensorField3D<T, 3> const &jField = dynamic_cast<TensorField3D<T, 3> const &>(*atomicBlocks[2]);

    BlockStatistics &stat = lattice.getInternalStatistics();

    static const plint vicinity = Descriptor<T>::vicinity;
    Box3D extDomain(domain.enlarge(vicinity));

    Dot3D offset1 = computeRelativeDisplacement(lattice, rhoBarField);
    Dot3D offset2 = computeRelativeDisplacement(lattice, jField);

    global::profiler().start("collStream");
    global::profiler().increment("collStreamCells", extDomain.nCells());

    // First, do the collision on cells within a boundary envelope of width
    // equal to the range of the lattice vectors (e.g. 1 for D2Q9)
    collide(
        lattice,
        Box3D(
            extDomain.x0, extDomain.x0 + vicinity - 1, extDomain.y0, extDomain.y1, extDomain.z0,
            extDomain.z1),
        rhoBarField, offset1, jField, offset2, stat);
    collide(
        lattice,
        Box3D(
            extDomain.x1 - vicinity + 1, extDomain.x1, extDomain.y0, extDomain.y1, extDomain.z0,
            extDomain.z1),
        rhoBarField, offset1, jField, offset2, stat);
    collide(
        lattice,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0,
            extDomain.y0 + vicinity - 1, extDomain.z0, extDomain.z1),
        rhoBarField, offset1, jField, offset2, stat);
    collide(
        lattice,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y1 - vicinity + 1,
            extDomain.y1, extDomain.z0, extDomain.z1),
        rhoBarField, offset1, jField, offset2, stat);
    collide(
        lattice,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0 + vicinity,
            extDomain.y1 - vicinity, extDomain.z0, extDomain.z0 + vicinity - 1),
        rhoBarField, offset1, jField, offset2, stat);
    collide(
        lattice,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0 + vicinity,
            extDomain.y1 - vicinity, extDomain.z1 - vicinity + 1, extDomain.z1),
        rhoBarField, offset1, jField, offset2, stat);

    // Then, do the efficient collideAndStream algorithm in the bulk,
    // excluding the envelope (this is efficient because there is no
    // if-then-else statement within the loop, given that the boundary
    // region is excluded)
    bulkCollideAndStream(
        lattice,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0 + vicinity,
            extDomain.y1 - vicinity, extDomain.z0 + vicinity, extDomain.z1 - vicinity),
        rhoBarField, offset1, jField, offset2, stat);

    // Finally, do streaming in the boundary envelope to conclude the
    // collision-stream cycle
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x0, extDomain.x0 + vicinity - 1, extDomain.y0, extDomain.y1, extDomain.z0,
            extDomain.z1));
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x1 - vicinity + 1, extDomain.x1, extDomain.y0, extDomain.y1, extDomain.z0,
            extDomain.z1));
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0,
            extDomain.y0 + vicinity - 1, extDomain.z0, extDomain.z1));
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y1 - vicinity + 1,
            extDomain.y1, extDomain.z0, extDomain.z1));
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0 + vicinity,
            extDomain.y1 - vicinity, extDomain.z0, extDomain.z0 + vicinity - 1));
    boundaryStream(
        lattice, extDomain,
        Box3D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0 + vicinity,
            extDomain.y1 - vicinity, extDomain.z1 - vicinity + 1, extDomain.z1));
    global::profiler().stop("collStream");
    global::timer("collideAndStream").stop();
}

template <typename T, template <typename U> class Descriptor>
OnLinkExternalRhoJcollideAndStream3D<T, Descriptor>
    *OnLinkExternalRhoJcollideAndStream3D<T, Descriptor>::clone() const
{
    return new OnLinkExternalRhoJcollideAndStream3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void OnLinkExternalRhoJcollideAndStream3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

/* ************* Class MaskedCollide3D ******************* */

template <typename T, template <typename U> class Descriptor>
MaskedCollide3D<T, Descriptor>::MaskedCollide3D(
    Dynamics<T, Descriptor> const &dynamics_, int flag_) :
    dynamics(dynamics_.clone()), flag(flag_)
{ }

template <typename T, template <typename U> class Descriptor>
MaskedCollide3D<T, Descriptor>::MaskedCollide3D(MaskedCollide3D<T, Descriptor> const &rhs) :
    dynamics(rhs.dynamics->clone()), flag(rhs.flag)
{ }

template <typename T, template <typename U> class Descriptor>
MaskedCollide3D<T, Descriptor> &MaskedCollide3D<T, Descriptor>::operator=(
    MaskedCollide3D<T, Descriptor> const &rhs)
{
    MaskedCollide3D<T, Descriptor>(rhs).swap(*this);
}

template <typename T, template <typename U> class Descriptor>
void MaskedCollide3D<T, Descriptor>::swap(MaskedCollide3D<T, Descriptor> &rhs)
{
    std::swap(dynamics, rhs.dynamics);
    std::swap(flag, rhs.flag);
}

template <typename T, template <typename U> class Descriptor>
void MaskedCollide3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &mask)
{
    BlockStatistics dummyStatistics;
    dummyStatistics.subscribeAverage();  // Subscribe average rho-bar
    dummyStatistics.subscribeAverage();  // Subscribe average uSqr
    dummyStatistics.subscribeMax();      // Subscribe max uSqr

    Dot3D offset = computeRelativeDisplacement(lattice, mask);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (mask.get(iX + offset.x, iY + offset.y, iZ + offset.z) == flag) {
                    dynamics->collide(lattice.get(iX, iY, iZ), dummyStatistics);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
MaskedCollide3D<T, Descriptor> *MaskedCollide3D<T, Descriptor>::clone() const
{
    return new MaskedCollide3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void MaskedCollide3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, template <typename U> class Descriptor>
void maskedCollide(
    MultiBlockLattice3D<T, Descriptor> &lattice, Dynamics<T, Descriptor> const &dynamics,
    MultiScalarField3D<int> &mask, int flag, Box3D domain)
{
    applyProcessingFunctional(
        new MaskedCollide3D<T, Descriptor>(dynamics, flag), domain, lattice, mask);
}

template <typename T, template <typename U> class Descriptor>
void maskedCollide(
    MultiBlockLattice3D<T, Descriptor> &lattice, Dynamics<T, Descriptor> const &dynamics,
    MultiScalarField3D<int> &mask, int flag)
{
    applyProcessingFunctional(
        new MaskedCollide3D<T, Descriptor>(dynamics, flag), lattice.getBoundingBox(), lattice,
        mask);
}

}  // namespace plb

#endif  // DYNAMICS_PROCESSOR_3D_HH
