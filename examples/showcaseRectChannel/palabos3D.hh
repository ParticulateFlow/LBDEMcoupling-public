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

#include <algorithm/headers3D.hh>
#include <atomicBlock/headers3D.hh>
#include <basicDynamics/headers3D.hh>
#include <boundaryCondition/headers3D.hh>
#include <coProcessors/headers3D.hh>
#include <complexDynamics/headers3D.hh>
#include <core/headers3D.hh>
#include <dataProcessors/headers3D.hh>
#include <finiteDifference/headers3D.hh>
#include <gridRefinement/headers3D.hh>
#include <io/headers3D.hh>
#include <latticeBoltzmann/headers3D.hh>
#include <libraryInterfaces/headers3D.hh>
#include <multiBlock/headers3D.hh>
#include <multiGrid/headers3D.hh>
#include <multiPhysics/headers3D.hh>
#include <offLattice/headers3D.hh>
#include <parallelism/headers3D.hh>
#include <particles/headers3D.hh>

// Include also 2D data structures. They are for
// example required to write 2D images.
#include <atomicBlock/headers2D.hh>
#include <core/headers2D.hh>
#include <multiBlock/headers2D.hh>
#include <multiGrid/headers2D.hh>
