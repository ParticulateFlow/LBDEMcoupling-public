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
 * Groups all the template files for the 3D off-lattice directory.
 */

#include "offLattice/asciiDataFileIO.hh"
#include "offLattice/boundaryShapes3D.hh"
#include "offLattice/bouzidiOffLatticeModel3D.hh"
#include "offLattice/connectedTriangleGenerator.hh"
#include "offLattice/connectedTriangleMesh.hh"
#include "offLattice/connectedTriangleSet.hh"
#include "offLattice/connectedTriangleUtil.hh"
#include "offLattice/filippovaHaenelOffLatticeModel3D.hh"
#include "offLattice/guoAdvDiffOffLatticeModel3D.hh"
#include "offLattice/guoOffLatticeModel3D.hh"
#include "offLattice/immersedAdvectionDiffusionWalls3D.hh"
#include "offLattice/immersedWalls3D.hh"
#include "offLattice/makeSparse3D.hh"
#include "offLattice/marchingCube.hh"
#include "offLattice/meiLuoShyyOffLatticeModel3D.hh"
#include "offLattice/meshAlgorithm.hh"
#include "offLattice/nextNeighbors3D.hh"
#include "offLattice/obsoleteFormatWrapper.hh"
#include "offLattice/offFileIO.hh"
#include "offLattice/offLatticeBoundaryCondition3D.hh"
#include "offLattice/offLatticeBoundaryProcessor3D.hh"
#include "offLattice/offLatticeBoundaryProfiles3D.hh"
#include "offLattice/offLatticeModel3D.hh"
#include "offLattice/rawTriangleMesh.hh"
#include "offLattice/stlFileIO.hh"
#include "offLattice/triangleBoundary3D.hh"
#include "offLattice/triangleHash.hh"
#include "offLattice/triangleSet.hh"
#include "offLattice/triangleSetGenerator.hh"
#include "offLattice/triangleToDef.hh"
#include "offLattice/triangleUtil.hh"
#include "offLattice/triangularSurfaceMesh.hh"
#include "offLattice/voxelizer.hh"
#include "offLattice/vtkFileIO.hh"

#ifndef PLB_BGP
#ifdef PLB_USE_EIGEN
#include "offLattice/generalizedOffLatticeModel3D.hh"
#endif
#endif
