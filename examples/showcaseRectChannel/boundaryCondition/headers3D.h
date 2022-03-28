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
 * #include "core/globalDefs.h"
 * Groups all the 3D include files in the boundaryConditions directory.
 */

#include "boundaryCondition/NLD_boundaries3D.h"
#include "boundaryCondition/NLD_boundaryDynamics3D.h"
#include "boundaryCondition/bounceBackModels.h"
#include "boundaryCondition/bounceBackModels3D.h"
#include "boundaryCondition/boundaryCondition3D.h"
#include "boundaryCondition/boundaryDynamics.h"
#include "boundaryCondition/equilibriumBoundaryDynamics.h"
#include "boundaryCondition/finiteDifferenceBoundaryProcessor3D.h"
#include "boundaryCondition/inamuroAnalyticalDynamics.h"
#include "boundaryCondition/neumannCondition3D.h"
#include "boundaryCondition/partialBBdynamics.h"
#include "boundaryCondition/regularizedBoundaryDynamics.h"
#include "boundaryCondition/regularizedBoundaryDynamics3D.h"
#include "boundaryCondition/spongeZones3D.h"
#include "boundaryCondition/wrappedLocalBoundaryProcessor3D.h"
#include "boundaryCondition/zouHeBoundary3D.h"
#include "boundaryCondition/zouHeDynamics.h"

#ifndef PLB_BGP
#ifdef PLB_USE_EIGEN
#include "boundaryCondition/generalizedBoundaryCondition3D.h"
#include "boundaryCondition/generalizedBoundaryDynamics.h"
#endif
#endif
