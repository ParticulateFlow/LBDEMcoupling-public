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

#ifndef BOX_LOGIC_3D_HH
#define BOX_LOGIC_3D_HH

#include <algorithm>
#include <string>
#include <vector>

#include "core/array.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "gridRefinement/boxLogic3D.h"
#include "gridRefinement/octree.h"
#include "offLattice/triangleSetGenerator.h"

namespace plb {

namespace boxLogic {

typedef OctreeTables OT;

// does a vector of boxes have two identical boxes?
// if the size of the box is of 1, then an exception can be made
template <class Boundary>
bool containsDuplicates(const std::vector<Boundary> &boxes, bool exceptSizeOne)
{
    for (plint iA = 0; iA < (plint)boxes.size(); ++iA) {
        Boundary localBox = boxes[iA];
        for (plint iB = iA + 1; iB < (plint)boxes.size(); ++iB) {
            if (localBox == boxes[iB]) {
                if (!(exceptSizeOne && (localBox.nCells() == 1))) {
                    std::cerr << "localBox = " << localBox.toStr() << std::endl;
                    std::cerr << "boxes[" << iB << "] = " << boxes[iB].toStr() << std::endl;
                    return true;
                }
            }
        }
    }
    return false;
}

// does a vector of boxes have no two identical boxes?
template <class Box>
bool containsNoDuplicates(const std::vector<Box> &boxes, bool exceptSizeOne)
{
    return !containsDuplicates<Box>(boxes, exceptSizeOne);
}

template <class Surf>
std::vector<Surf> merge(const std::vector<Surf> &origBoxes)
{
    std::vector<Surf> result = origBoxes;

    plint size = 0;
    while (size != (plint)result.size()) {
        size = (plint)result.size();
        for (plint iA = 0; iA < (plint)result.size(); ++iA) {
            for (plint iB = iA + 1; iB < (plint)result.size(); ++iB) {
                bool merged = merge(result[iA], result[iB]);
                if (merged) {
                    result.erase(result.begin() + iB);
                    iB -= 1;
                }
            }
        }
    }

    return result;
}

}  // namespace boxLogic

}  // namespace plb

#endif  // BOX_LOGIC_3D_HH
