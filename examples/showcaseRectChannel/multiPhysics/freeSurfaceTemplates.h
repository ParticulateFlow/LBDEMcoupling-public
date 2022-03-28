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

#ifndef FREE_SURFACE_TEMPLATES_3D_HH
#define FREE_SURFACE_TEMPLATES_3D_HH

#include "atomicBlock/atomicContainerBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "core/block3D.h"
#include "core/globalDefs.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "multiPhysics/freeSurfaceModel3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
struct freeSurfaceTemplates {
    static void massExchangeFluidCell(
        FreeSurfaceProcessorParam3D<T, Descriptor> &param, plint iX, plint iY, plint iZ)
    {
        typedef Descriptor<T> D;
        using namespace freeSurfaceFlag;
        // Calculate mass at time t+1 --> eq 6 Thurey's paper.
        for (plint iPop = 1; iPop < D::q; ++iPop) {
            plint nextX = iX + D::c[iPop][0];
            plint nextY = iY + D::c[iPop][1];
            plint nextZ = iZ + D::c[iPop][2];
            int nextFlag = param.flag(nextX, nextY, nextZ);
            if (isWet(nextFlag)) {
                // In Thuerey's paper, the mass balance is computed locally on one cell, but
                // N. Thuerey uses outgoing populations. Palabos works with incoming populations
                // and uses the relation f_out_i(x,t) = f_in_opp(i)(x+c_i,t+1).
                plint opp = indexTemplates::opposite<D>(iPop);
                param.mass(iX, iY, iZ) +=
                    param.cell(iX, iY, iZ)[opp] - param.cell(nextX, nextY, nextZ)[iPop];
            }
        }
    }
};

template <typename T>
struct freeSurfaceTemplates<T, descriptors::ForcedD3Q19Descriptor> {
    static void massExchangeFluidCell(
        FreeSurfaceProcessorParam3D<T, descriptors::ForcedD3Q19Descriptor> &param, plint iX,
        plint iY, plint iZ)
    {
        using namespace freeSurfaceFlag;
        // Calculate mass at time t+1 --> eq 6 Thurey's paper.
        plint nextX, nextY, nextZ, nextFlag;

        nextX = iX + -1;
        nextY = iY + 0;
        nextZ = iZ + 0;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[10] - param.cell(nextX, nextY, nextZ)[1];
        }

        nextX = iX + 0;
        nextY = iY + -1;
        nextZ = iZ + 0;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[11] - param.cell(nextX, nextY, nextZ)[2];
        }

        nextX = iX + 0;
        nextY = iY + 0;
        nextZ = iZ + -1;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[12] - param.cell(nextX, nextY, nextZ)[3];
        }

        nextX = iX + -1;
        nextY = iY + -1;
        nextZ = iZ + 0;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[13] - param.cell(nextX, nextY, nextZ)[4];
        }

        nextX = iX + -1;
        nextY = iY + 1;
        nextZ = iZ + 0;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[14] - param.cell(nextX, nextY, nextZ)[5];
        }

        nextX = iX + -1;
        nextY = iY + 0;
        nextZ = iZ + -1;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[15] - param.cell(nextX, nextY, nextZ)[6];
        }

        nextX = iX + -1;
        nextY = iY + 0;
        nextZ = iZ + 1;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[16] - param.cell(nextX, nextY, nextZ)[7];
        }

        nextX = iX + 0;
        nextY = iY + -1;
        nextZ = iZ + -1;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[17] - param.cell(nextX, nextY, nextZ)[8];
        }

        nextX = iX + 0;
        nextY = iY + -1;
        nextZ = iZ + 1;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[18] - param.cell(nextX, nextY, nextZ)[9];
        }

        nextX = iX + 1;
        nextY = iY + 0;
        nextZ = iZ + 0;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[1] - param.cell(nextX, nextY, nextZ)[10];
        }

        nextX = iX + 0;
        nextY = iY + 1;
        nextZ = iZ + 0;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[2] - param.cell(nextX, nextY, nextZ)[11];
        }

        nextX = iX + 0;
        nextY = iY + 0;
        nextZ = iZ + 1;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[3] - param.cell(nextX, nextY, nextZ)[12];
        }

        nextX = iX + 1;
        nextY = iY + 1;
        nextZ = iZ + 0;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[4] - param.cell(nextX, nextY, nextZ)[13];
        }

        nextX = iX + 1;
        nextY = iY + -1;
        nextZ = iZ + 0;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[5] - param.cell(nextX, nextY, nextZ)[14];
        }

        nextX = iX + 1;
        nextY = iY + 0;
        nextZ = iZ + 1;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[6] - param.cell(nextX, nextY, nextZ)[15];
        }

        nextX = iX + 1;
        nextY = iY + 0;
        nextZ = iZ + -1;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[7] - param.cell(nextX, nextY, nextZ)[16];
        }

        nextX = iX + 0;
        nextY = iY + 1;
        nextZ = iZ + 1;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[8] - param.cell(nextX, nextY, nextZ)[17];
        }

        nextX = iX + 0;
        nextY = iY + 1;
        nextZ = iZ + -1;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[9] - param.cell(nextX, nextY, nextZ)[18];
        }
    }
};

template <typename T>
struct freeSurfaceTemplates<T, descriptors::D3Q19Descriptor> {
    static void massExchangeFluidCell(
        FreeSurfaceProcessorParam3D<T, descriptors::D3Q19Descriptor> &param, plint iX, plint iY,
        plint iZ)
    {
        using namespace freeSurfaceFlag;
        // Calculate mass at time t+1 --> eq 6 Thurey's paper.
        plint nextX, nextY, nextZ, nextFlag;

        nextX = iX + -1;
        nextY = iY + 0;
        nextZ = iZ + 0;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[10] - param.cell(nextX, nextY, nextZ)[1];
        }

        nextX = iX + 0;
        nextY = iY + -1;
        nextZ = iZ + 0;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[11] - param.cell(nextX, nextY, nextZ)[2];
        }

        nextX = iX + 0;
        nextY = iY + 0;
        nextZ = iZ + -1;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[12] - param.cell(nextX, nextY, nextZ)[3];
        }

        nextX = iX + -1;
        nextY = iY + -1;
        nextZ = iZ + 0;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[13] - param.cell(nextX, nextY, nextZ)[4];
        }

        nextX = iX + -1;
        nextY = iY + 1;
        nextZ = iZ + 0;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[14] - param.cell(nextX, nextY, nextZ)[5];
        }

        nextX = iX + -1;
        nextY = iY + 0;
        nextZ = iZ + -1;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[15] - param.cell(nextX, nextY, nextZ)[6];
        }

        nextX = iX + -1;
        nextY = iY + 0;
        nextZ = iZ + 1;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[16] - param.cell(nextX, nextY, nextZ)[7];
        }

        nextX = iX + 0;
        nextY = iY + -1;
        nextZ = iZ + -1;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[17] - param.cell(nextX, nextY, nextZ)[8];
        }

        nextX = iX + 0;
        nextY = iY + -1;
        nextZ = iZ + 1;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[18] - param.cell(nextX, nextY, nextZ)[9];
        }

        nextX = iX + 1;
        nextY = iY + 0;
        nextZ = iZ + 0;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[1] - param.cell(nextX, nextY, nextZ)[10];
        }

        nextX = iX + 0;
        nextY = iY + 1;
        nextZ = iZ + 0;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[2] - param.cell(nextX, nextY, nextZ)[11];
        }

        nextX = iX + 0;
        nextY = iY + 0;
        nextZ = iZ + 1;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[3] - param.cell(nextX, nextY, nextZ)[12];
        }

        nextX = iX + 1;
        nextY = iY + 1;
        nextZ = iZ + 0;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[4] - param.cell(nextX, nextY, nextZ)[13];
        }

        nextX = iX + 1;
        nextY = iY + -1;
        nextZ = iZ + 0;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[5] - param.cell(nextX, nextY, nextZ)[14];
        }

        nextX = iX + 1;
        nextY = iY + 0;
        nextZ = iZ + 1;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[6] - param.cell(nextX, nextY, nextZ)[15];
        }

        nextX = iX + 1;
        nextY = iY + 0;
        nextZ = iZ + -1;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[7] - param.cell(nextX, nextY, nextZ)[16];
        }

        nextX = iX + 0;
        nextY = iY + 1;
        nextZ = iZ + 1;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[8] - param.cell(nextX, nextY, nextZ)[17];
        }

        nextX = iX + 0;
        nextY = iY + 1;
        nextZ = iZ + -1;
        nextFlag = param.flag(nextX, nextY, nextZ);
        if (isWet(nextFlag)) {
            param.mass(iX, iY, iZ) +=
                param.cell(iX, iY, iZ)[9] - param.cell(nextX, nextY, nextZ)[18];
        }
    }
};

}  // namespace plb

#endif  // FREE_SURFACE_TEMPLATES_3D_HH
