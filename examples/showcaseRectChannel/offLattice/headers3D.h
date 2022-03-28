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
 * Groups all the include files for the 3D off-lattice directory.
 */

#include "offLattice/asciiDataFileIO.h"
#include "offLattice/boundaryShapes3D.h"
#include "offLattice/bouzidiOffLatticeModel3D.h"
#include "offLattice/connectedTriangleGenerator.h"
#include "offLattice/connectedTriangleMesh.h"
#include "offLattice/connectedTriangleSet.h"
#include "offLattice/connectedTriangleUtil.h"
#include "offLattice/filippovaHaenelOffLatticeModel3D.h"
#include "offLattice/guoAdvDiffOffLatticeModel3D.h"
#include "offLattice/guoOffLatticeModel3D.h"
#include "offLattice/immersedAdvectionDiffusionWalls3D.h"
#include "offLattice/immersedWalls3D.h"
#include "offLattice/makeSparse3D.h"
#include "offLattice/marchingCube.h"
#include "offLattice/meiLuoShyyOffLatticeModel3D.h"
#include "offLattice/meshAlgorithm.h"
#include "offLattice/nextNeighbors3D.h"
#include "offLattice/obsoleteFormatWrapper.h"
#include "offLattice/offFileIO.h"
#include "offLattice/offLatticeBoundaryCondition3D.h"
#include "offLattice/offLatticeBoundaryProcessor3D.h"
#include "offLattice/offLatticeBoundaryProfiles3D.h"
#include "offLattice/offLatticeModel3D.h"
#include "offLattice/rawTriangleMesh.h"
#include "offLattice/stlFileIO.h"
#include "offLattice/triangleBoundary3D.h"
#include "offLattice/triangleHash.h"
#include "offLattice/triangleMesh.h"
#include "offLattice/triangleSelector.h"
#include "offLattice/triangleSet.h"
#include "offLattice/triangleSetGenerator.h"
#include "offLattice/triangleToDef.h"
#include "offLattice/triangleUtil.h"
#include "offLattice/triangularSurfaceMesh.h"
#include "offLattice/voxelizer.h"
#include "offLattice/vtkFileIO.h"

#ifndef PLB_BGP
#ifdef PLB_USE_EIGEN
#include "offLattice/generalizedOffLatticeModel3D.h"
#endif
#endif
