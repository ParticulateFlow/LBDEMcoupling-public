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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- generic implementation.
 */
#ifndef BOUNDARY_TEMPLATES_H
#define BOUNDARY_TEMPLATES_H

#include "atomicBlock/blockLattice3D.h"
#include "boundaryCondition/regularizedBoundaryDynamics.h"
#include "core/cell.h"
#include "latticeBoltzmann/indexTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"

namespace plb {

/// Computation of flat-wall quantities with static genericity: direction and orientation
/// are template parameters.
template <typename T, template <typename U> class Descriptor, int direction, int orientation>
struct boundaryTemplates {
    static void compute_PiNeq(
        Dynamics<T, Descriptor> const &dynamics, Cell<T, Descriptor> const &cell, T rhoBar,
        Array<T, Descriptor<T>::d> const &j, T jSqr,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq)
    {
        typedef Descriptor<T> L;

        std::vector<plint> const &onWallIndices = indexTemplates::subIndex<L, direction, 0>();
        std::vector<plint> const &normalIndices =
            indexTemplates::subIndex<L, direction, orientation>();

        // Compute off-equilibrium for known particle populations.
        Array<T, Descriptor<T>::q> fNeq;
        for (pluint fIndex = 0; fIndex < onWallIndices.size(); ++fIndex) {
            plint iPop = onWallIndices[fIndex];
            fNeq[iPop] = cell[iPop] - dynamics.computeEquilibrium(iPop, rhoBar, j, jSqr);
        }
        for (pluint fIndex = 0; fIndex < normalIndices.size(); ++fIndex) {
            plint iPop = normalIndices[fIndex];
            if (iPop == 0) {
                fNeq[iPop] = T();  // fNeq[0] will not be used anyway
            } else {
                fNeq[iPop] = cell[iPop] - dynamics.computeEquilibrium(iPop, rhoBar, j, jSqr);
            }
        }

        // Compute PiNeq from fNeq, by using "bounce-back of off-equilibrium part" rule.
        int iPi = 0;
        for (int iAlpha = 0; iAlpha < L::d; ++iAlpha) {
            for (int iBeta = iAlpha; iBeta < L::d; ++iBeta) {
                PiNeq[iPi] = T();
                for (pluint fIndex = 0; fIndex < onWallIndices.size(); ++fIndex) {
                    const plint iPop = onWallIndices[fIndex];
                    PiNeq[iPi] += L::c[iPop][iAlpha] * L::c[iPop][iBeta] * fNeq[iPop];
                }
                for (pluint fIndex = 0; fIndex < normalIndices.size(); ++fIndex) {
                    const plint iPop = normalIndices[fIndex];
                    PiNeq[iPi] += (T)2 * L::c[iPop][iAlpha] * L::c[iPop][iBeta] * fNeq[iPop];
                }
                ++iPi;
            }
        }
    }

    static void compute_jNeq(
        Dynamics<T, Descriptor> const &dynamics, Cell<T, Descriptor> const &cell, T rhoBar,
        Array<T, Descriptor<T>::d> const &j, T jSqr, Array<T, Descriptor<T>::d> &jNeq)
    {
        typedef Descriptor<T> D;

        std::vector<plint> const &onWallIndices = indexTemplates::subIndex<D, direction, 0>();
        std::vector<plint> const &normalIndices =
            indexTemplates::subIndex<D, direction, orientation>();

        // Compute off-equilibrium for known particle populations.
        Array<T, Descriptor<T>::q> fNeq;
        for (pluint fIndex = 0; fIndex < onWallIndices.size(); ++fIndex) {
            plint iPop = onWallIndices[fIndex];
            fNeq[iPop] = cell[iPop] - dynamics.computeEquilibrium(iPop, rhoBar, j, jSqr);
        }
        for (pluint fIndex = 0; fIndex < normalIndices.size(); ++fIndex) {
            plint iPop = normalIndices[fIndex];
            plint iOpp = indexTemplates::opposite<D>(iPop);
            if (iPop == 0) {
                fNeq[iPop] = T();  // fNeq[0] will not be used anyway
            } else {
                fNeq[iPop] = cell[iPop] - dynamics.computeEquilibrium(iPop, rhoBar, j, jSqr);
                fNeq[iOpp] = -fNeq[iPop];
            }
        }

        momentTemplatesImpl<T, D>::get_j(fNeq, jNeq);
    }

};  // struct boundaryTemplates

/// Computation of flat-wall quantities with dynamic genericity: direction and
/// orientation are dynamic arguments.
namespace flatWall {

template <typename T, template <typename U> class Descriptor>
T computeRhoBar(Cell<T, Descriptor> &cell, int direction, int orientation, T u, T f)
{
    T rhoOnWall = T();
    T rhoNormal = T();
    for (int iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        if (Descriptor<T>::c[iPop][direction] == 0) {
            rhoOnWall += cell[iPop];
        } else if (Descriptor<T>::c[iPop][direction] == orientation) {
            rhoNormal += cell[iPop];
        }
    }

    T velNormal = (T)orientation * (u - (T)0.5 * f);
    T rhoBar = T();
    if (cell.getDynamics().velIsJ()) {
        rhoBar = (T)2 * rhoNormal + rhoOnWall - velNormal;
    } else {
        rhoBar = ((T)2 * rhoNormal + rhoOnWall - Descriptor<T>::SkordosFactor() * velNormal)
                 / ((T)1 + velNormal);
    }
    return rhoBar;
}

template <typename T, template <typename U> class Descriptor>
Array<T, 3> computeJ(
    Cell<T, Descriptor> const &cell, int direction, int orientation, Array<T, 3> j, T rhoBar)
{
    T rhoOnWall = T(), rhoNormal = T();
    for (int iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        if (Descriptor<T>::c[iPop][direction] == 0) {
            rhoOnWall += cell[iPop];
        } else if (Descriptor<T>::c[iPop][direction] == orientation) {
            rhoNormal += cell[iPop];
        }
    }

    j[direction] = (T)orientation * ((T)2 * rhoNormal + rhoOnWall - rhoBar);
    return j;
}

template <typename T, template <typename U> class Descriptor>
void compute_PiNeq(
    Dynamics<T, Descriptor> const &dynamics, Cell<T, Descriptor> const &cell, int direction,
    int orientation, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq)
{
    typedef Descriptor<T> L;

    // Compute off-equilibrium for known particle populations.
    Array<T, Descriptor<T>::q> fNeq;

    for (int iPop = 0; iPop < L::q; ++iPop) {
        if (iPop == 0) {
            fNeq[iPop] = T();  // fNeq[0] will not be used anyway
        } else {
            int cNormal = L::c[iPop][direction];
            if (cNormal == 0 || cNormal == orientation) {
                fNeq[iPop] = cell[iPop] - dynamics.computeEquilibrium(iPop, rhoBar, j, jSqr);
            }
        }
    }

    // Compute PiNeq from fNeq, by using "bounce-back of off-equilibrium part" rule.
    int iPi = 0;
    for (int iAlpha = 0; iAlpha < L::d; ++iAlpha) {
        for (int iBeta = iAlpha; iBeta < L::d; ++iBeta) {
            PiNeq[iPi] = T();
            for (int iPop = 0; iPop < L::q; ++iPop) {
                int cNormal = L::c[iPop][direction];
                if (cNormal == 0) {
                    PiNeq[iPi] += L::c[iPop][iAlpha] * L::c[iPop][iBeta] * fNeq[iPop];
                } else if (cNormal == orientation) {
                    PiNeq[iPi] += (T)2 * L::c[iPop][iAlpha] * L::c[iPop][iBeta] * fNeq[iPop];
                }
            }
            ++iPi;
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void extrapolatePopulationsFixPressure(
    Cell<T, Descriptor> const &cellFrom, Cell<T, Descriptor> &cellTo, int direction,
    int orientation)
{
    for (pluint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        if (Descriptor<T>::c[iPop][direction] == -orientation) {
            cellTo[iPop] = cellFrom[iPop];
        }
    }
    T rhoBar;
    Array<T, 3> j;
    cellTo.getDynamics().computeRhoBarJ(cellTo, rhoBar, j);
    Array<T, Descriptor<T>::q> oldFeq, newFeq;
    T jSqr = normSqr(j);
    cellTo.getDynamics().computeEquilibria(oldFeq, rhoBar, j, jSqr);
    cellTo.getDynamics().computeEquilibria(newFeq, T(), j, jSqr);
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cellTo[iPop] += newFeq[iPop] - oldFeq[iPop];
    }
}

template <typename T, template <typename U> class Descriptor>
void extrapolatePopulationsFixPressure(
    BlockLattice3D<T, Descriptor> &lattice, plint iX, plint iY, plint iZ, int direction,
    int orientation)
{
    Array<plint, 3> neighb(iX, iY, iZ);
    neighb[direction] -= orientation;
    Cell<T, Descriptor> &cellTo = lattice.get(iX, iY, iZ);
    Cell<T, Descriptor> &cellFrom = lattice.get(neighb[0], neighb[1], neighb[2]);
    extrapolatePopulationsFixPressure(cellFrom, cellTo, direction, orientation);
}

}  // namespace flatWall

}  // namespace plb

#endif  // BOUNDARY_TEMPLATES_H
