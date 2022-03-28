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
 * Helper functions for the implementation of lattice operations. This file is all
 * about efficiency. The generic template code is specialized for commonly
 * used Lattices, so that a maximum performance can be taken out of each
 * case.
 */
#ifndef LATTICE_TEMPLATES_H
#define LATTICE_TEMPLATES_H

#include "core/cell.h"
#include "core/globalDefs.h"
#include "core/util.h"

namespace plb {

/// Helper functions with full-lattice access
template <typename T, template <typename U> class Descriptor>
struct latticeTemplates {
    /// Swap ("bounce-back") values of a cell (2D), and apply streaming step
    static void swapAndStream2D(Cell<T, Descriptor> **grid, plint iX, plint iY)
    {
        const plint half = Descriptor<T>::q / 2;
        for (plint iPop = 1; iPop <= half; ++iPop) {
            plint nextX = iX + Descriptor<T>::c[iPop][0];
            plint nextY = iY + Descriptor<T>::c[iPop][1];
            T fTmp = grid[iX][iY][iPop];
            grid[iX][iY][iPop] = grid[iX][iY][iPop + half];
            grid[iX][iY][iPop + half] = grid[nextX][nextY][iPop];
            grid[nextX][nextY][iPop] = fTmp;
        }
    }

    /// Swap ("bounce-back") values of a cell (3D), and apply streaming step
    static void swapAndStream3D(Cell<T, Descriptor> ***grid, plint iX, plint iY, plint iZ)
    {
        const plint half = Descriptor<T>::q / 2;
        for (plint iPop = 1; iPop <= half; ++iPop) {
            plint nextX = iX + Descriptor<T>::c[iPop][0];
            plint nextY = iY + Descriptor<T>::c[iPop][1];
            plint nextZ = iZ + Descriptor<T>::c[iPop][2];
            T fTmp = grid[iX][iY][iZ][iPop];
            grid[iX][iY][iZ][iPop] = grid[iX][iY][iZ][iPop + half];
            grid[iX][iY][iZ][iPop + half] = grid[nextX][nextY][nextZ][iPop];
            grid[nextX][nextY][nextZ][iPop] = fTmp;
        }
    }
};

}  // namespace plb

#include "latticeBoltzmann/latticeTemplates2D.h"
#include "latticeBoltzmann/latticeTemplates3D.h"

#endif
