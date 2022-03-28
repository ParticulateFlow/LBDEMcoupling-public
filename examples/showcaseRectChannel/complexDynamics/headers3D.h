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
 * Groups all the 3D include files in the complexDynamics directory.
 */

#include "complexDynamics/adiabaticBoundaryProcessor3D.h"
#include "complexDynamics/advectionDiffusionBoundaries.h"
#include "complexDynamics/advectionDiffusionBoundaryCondition3D.h"
#include "complexDynamics/advectionDiffusionDynamics.h"
#include "complexDynamics/advectionDiffusionProcessor3D.h"
#include "complexDynamics/advectionDiffusionUnits.h"
#include "complexDynamics/asinariModel.h"
#include "complexDynamics/carreauDynamics.h"
#include "complexDynamics/carreauDynamicsTemplates.h"
#include "complexDynamics/carreauGlobalDefs.h"
#include "complexDynamics/carreauUnits.h"
#include "complexDynamics/dynamicSmagorinskyDynamics.h"
#include "complexDynamics/dynamicSmagorinskyLattices3D.h"
#include "complexDynamics/dynamicSmagorinskyProcessor3D.h"
#include "complexDynamics/entropicDynamics.h"
#include "complexDynamics/externalForceMrtDynamics.h"
#include "complexDynamics/mrtDynamics.h"
#include "complexDynamics/smagorinskyDynamics.h"
#include "complexDynamics/smagorinskyDynamics3D.h"
#include "complexDynamics/trtDynamics.h"
#include "complexDynamics/variableOmegaDynamics.h"
#include "complexDynamics/wavePropagation.h"
