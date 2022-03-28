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
 * Groups all the include files for the 2D multiBlock.
 */

#include "multiBlock/combinedStatistics.h"
#include "multiBlock/coupling2D.h"
#include "multiBlock/defaultMultiBlockPolicy2D.h"
#include "multiBlock/group2D.h"
#include "multiBlock/localMultiBlockInfo2D.h"
#include "multiBlock/multiBlock2D.h"
#include "multiBlock/multiBlockGenerator2D.h"
#include "multiBlock/multiBlockInfo2D.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiBlockManagement2D.h"
#include "multiBlock/multiBlockOperations2D.h"
#include "multiBlock/multiBlockSerializer2D.h"
#include "multiBlock/multiContainerBlock2D.h"
#include "multiBlock/multiDataField2D.h"
#include "multiBlock/multiDataProcessorWrapper2D.h"
#include "multiBlock/nonLocalTransfer2D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper2D.h"
#include "multiBlock/serialBlockCommunicator2D.h"
#include "multiBlock/serialMultiBlockLattice2D.h"
#include "multiBlock/serialMultiDataField2D.h"
#include "multiBlock/staticRepartitions2D.h"
