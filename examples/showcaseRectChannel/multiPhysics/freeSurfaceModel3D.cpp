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

#include "multiPhysics/freeSurfaceModel3D.h"

#include "atomicBlock/atomicContainerBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "core/block3D.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "multiPhysics/freeSurfaceUtil3D.h"

#ifdef PLB_MPI_PARALLEL
// DISABLE_WARNING_PUSH
// DISABLE_WARNING_CAST_FUNCTION_TYPE
#include <mpi.h>
// DISABLE_WARNING_POP
#endif

namespace plb {

#ifdef PLB_MPI_PARALLEL

/* *************** Class FreeSurfaceComputeReductions3D *******************************************
 */

// CAUTION: This data processor will not work (will hang) if the multi-block passed does
//          not have at least one atomic block on every MPI process.
void FreeSurfaceComputeReductions3D::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 1);

    if (!reductionData.used) {
        double sendbuf[3] = {
            reductionData.localLostMass, reductionData.localTotalMass,
            (double)reductionData.localNumInterfaceCells};
        double recvbuf[3] = {0.0, 0.0, 0.0};

        (void)MPI_Allreduce(
            static_cast<void *>(sendbuf), static_cast<void *>(recvbuf), 3, MPI_DOUBLE, MPI_SUM,
            reductionCommunicator);

        reductionData.lostMass = recvbuf[0];
        reductionData.totalMass = recvbuf[1];
        reductionData.numInterfaceCells = util::roundToInt(recvbuf[2]);

        reductionData.used = true;
    }
}

#else

/* *************** Class FreeSurfaceComputeSerialReductions3D
 * ******************************************* */

void FreeSurfaceComputeSerialReductions3D::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 1);

    if (!reductionData.used) {
        reductionData.lostMass = reductionData.localLostMass;
        reductionData.totalMass = reductionData.localTotalMass;
        reductionData.numInterfaceCells = reductionData.localNumInterfaceCells;

        reductionData.used = true;
    }
}

#endif

/* *************** Class FreeSurfaceResetReductionData3D *******************************************
 */

void FreeSurfaceResetReductionData3D::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 1);
    reductionData.reset();
}

}  // namespace plb
