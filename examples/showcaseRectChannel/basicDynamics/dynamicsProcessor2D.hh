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

#ifndef DYNAMICS_PROCESSOR_2D_HH
#define DYNAMICS_PROCESSOR_2D_HH

#include "atomicBlock/blockLattice2D.h"
#include "basicDynamics/dynamicsProcessor2D.h"
#include "core/cell.h"
#include "core/plbProfiler.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/latticeTemplates.h"
#include "multiGrid/multiGridUtil.h"

namespace plb {

/* ************* Class ExternalRhoJcollideAndStream2D ******************* */

template <typename T, template <typename U> class Descriptor>
void ExternalRhoJcollideAndStream2D<T, Descriptor>::collide(
    BlockLattice2D<T, Descriptor> &lattice, Box2D const &domain,
    ScalarField2D<T> const &rhoBarField, Dot2D const &offset1, TensorField2D<T, 2> const &jField,
    Dot2D const &offset2, BlockStatistics &stat)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Cell<T, Descriptor> &cell = lattice.get(iX, iY);
            T rhoBar = rhoBarField.get(iX + offset1.x, iY + offset1.y);
            Array<T, 2> const &j = jField.get(iX + offset2.x, iY + offset2.y);
            cell.getDynamics().collideExternal(cell, rhoBar, j, T(), stat);
            cell.revert();
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void ExternalRhoJcollideAndStream2D<T, Descriptor>::bulkCollideAndStream(
    BlockLattice2D<T, Descriptor> &lattice, Box2D const &domain,
    ScalarField2D<T> const &rhoBarField, Dot2D const &offset1, TensorField2D<T, 2> const &jField,
    Dot2D const &offset2, BlockStatistics &stat)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Cell<T, Descriptor> &cell = lattice.get(iX, iY);
            T rhoBar = rhoBarField.get(iX + offset1.x, iY + offset1.y);
            Array<T, 2> const &j = jField.get(iX + offset2.x, iY + offset2.y);
            cell.getDynamics().collideExternal(cell, rhoBar, j, T(), stat);
            latticeTemplates<T, Descriptor>::swapAndStream2D(lattice.grid, iX, iY);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void ExternalRhoJcollideAndStream2D<T, Descriptor>::boundaryStream(
    BlockLattice2D<T, Descriptor> &lattice, Box2D const &bound, Box2D const &domain)
{
    // Make sure domain is contained within bound
    PLB_PRECONDITION(contained(domain, bound));

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iPop = 1; iPop <= Descriptor<T>::q / 2; ++iPop) {
                plint nextX = iX + Descriptor<T>::c[iPop][0];
                plint nextY = iY + Descriptor<T>::c[iPop][1];
                if (nextX >= bound.x0 && nextX <= bound.x1 && nextY >= bound.y0
                    && nextY <= bound.y1) {
                    std::swap(
                        lattice.grid[iX][iY][iPop + Descriptor<T>::q / 2],
                        lattice.grid[nextX][nextY][iPop]);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void ExternalRhoJcollideAndStream2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    BlockLattice2D<T, Descriptor> &lattice =
        dynamic_cast<BlockLattice2D<T, Descriptor> &>(*atomicBlocks[0]);
    ScalarField2D<T> const &rhoBarField = dynamic_cast<ScalarField2D<T> const &>(*atomicBlocks[1]);
    TensorField2D<T, 2> const &jField = dynamic_cast<TensorField2D<T, 2> const &>(*atomicBlocks[2]);

    BlockStatistics &stat = lattice.getInternalStatistics();

    static const plint vicinity = Descriptor<T>::vicinity;
    Box2D extDomain(domain.enlarge(vicinity));

    Dot2D offset1 = computeRelativeDisplacement(lattice, rhoBarField);
    Dot2D offset2 = computeRelativeDisplacement(lattice, jField);

    global::profiler().start("collStream");
    global::profiler().increment("collStreamCells", extDomain.nCells());

    // First, do the collision on cells within a boundary envelope of width
    // equal to the range of the lattice vectors (e.g. 1 for D2Q9)
    collide(
        lattice, Box2D(extDomain.x0, extDomain.x0 + vicinity - 1, extDomain.y0, extDomain.y1),
        rhoBarField, offset1, jField, offset2, stat);
    collide(
        lattice, Box2D(extDomain.x1 - vicinity + 1, extDomain.x1, extDomain.y0, extDomain.y1),
        rhoBarField, offset1, jField, offset2, stat);
    collide(
        lattice,
        Box2D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0,
            extDomain.y0 + vicinity - 1),
        rhoBarField, offset1, jField, offset2, stat);
    collide(
        lattice,
        Box2D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y1 - vicinity + 1,
            extDomain.y1),
        rhoBarField, offset1, jField, offset2, stat);

    // Then, do the efficient collideAndStream algorithm in the bulk,
    // excluding the envelope (this is efficient because there is no
    // if-then-else statement within the loop, given that the boundary
    // region is excluded)
    bulkCollideAndStream(
        lattice,
        Box2D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0 + vicinity,
            extDomain.y1 - vicinity),
        rhoBarField, offset1, jField, offset2, stat);

    // Finally, do streaming in the boundary envelope to conclude the
    // collision-stream cycle
    boundaryStream(
        lattice, extDomain,
        Box2D(extDomain.x0, extDomain.x0 + vicinity - 1, extDomain.y0, extDomain.y1));
    boundaryStream(
        lattice, extDomain,
        Box2D(extDomain.x1 - vicinity + 1, extDomain.x1, extDomain.y0, extDomain.y1));
    boundaryStream(
        lattice, extDomain,
        Box2D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y0,
            extDomain.y0 + vicinity - 1));
    boundaryStream(
        lattice, extDomain,
        Box2D(
            extDomain.x0 + vicinity, extDomain.x1 - vicinity, extDomain.y1 - vicinity + 1,
            extDomain.y1));

    global::profiler().stop("collStream");
}

template <typename T, template <typename U> class Descriptor>
ExternalRhoJcollideAndStream2D<T, Descriptor>
    *ExternalRhoJcollideAndStream2D<T, Descriptor>::clone() const
{
    return new ExternalRhoJcollideAndStream2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void ExternalRhoJcollideAndStream2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

/* ************* Class Tau1CollideAndStream2D ******************* */

template <typename T, template <typename U> class Descriptor>
void Tau1CollideAndStream2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    std::vector<T> rhoBarCol(domain.getNy());
    std::vector<Array<T, 2> > jCol1(domain.getNy());
    Array<T, Descriptor<T>::q> fEq1(domain.getNy());
    std::vector<Array<T, 2> > jCol2(domain.getNy());
    Array<T, Descriptor<T>::q> fEq2(domain.getNy());

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Cell<T, Descriptor> &cell = lattice.get(iX, iY);
        }
    }
}

}  // namespace plb

#endif  // DYNAMICS_PROCESSOR_2D_HH
