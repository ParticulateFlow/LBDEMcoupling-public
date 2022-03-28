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

#include "multiPhysics/freeSurfaceAnalysis3D.h"

#include "atomicBlock/reductiveDataProcessingFunctional3D.h"
#include "atomicBlock/reductiveDataProcessingFunctional3D.hh"
#include "multiBlock/reductiveMultiDataProcessorWrapper3D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper3D.hh"
#include "multiPhysics/freeSurfaceUtil3D.h"

namespace plb {

/* ************** class FreeSurfaceCountWallCells3D ********************************** */

FreeSurfaceCountWallCells3D::FreeSurfaceCountWallCells3D() :
    sumScalarId(this->getStatistics().subscribeIntSum())
{ }

void FreeSurfaceCountWallCells3D::process(Box3D domain, ScalarField3D<int> &flag)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (freeSurfaceFlag::isAnyWall(flag.get(iX, iY, iZ))) {
                    statistics.gatherIntSum(sumScalarId, (plint)1);
                }
            }
        }
    }
}

FreeSurfaceCountWallCells3D *FreeSurfaceCountWallCells3D::clone() const
{
    return new FreeSurfaceCountWallCells3D(*this);
}

plint FreeSurfaceCountWallCells3D::getNumWallCells() const
{
    return this->getStatistics().getIntSum(sumScalarId);
}

plint freeSurfaceCountWallCells(MultiScalarField3D<int> &flag, Box3D domain)
{
    FreeSurfaceCountWallCells3D functional;
    applyProcessingFunctional(functional, domain, flag);
    return functional.getNumWallCells();
}

/* ************** class FreeSurfaceComputeFluidBoundingBox3D ********************************** */

FreeSurfaceComputeFluidBoundingBox3D::FreeSurfaceComputeFluidBoundingBox3D() :
    minIdX(this->getStatistics().subscribeMax()),
    minIdY(this->getStatistics().subscribeMax()),
    minIdZ(this->getStatistics().subscribeMax()),
    maxIdX(this->getStatistics().subscribeMax()),
    maxIdY(this->getStatistics().subscribeMax()),
    maxIdZ(this->getStatistics().subscribeMax())
{ }

void FreeSurfaceComputeFluidBoundingBox3D::process(Box3D domain, ScalarField3D<int> &flag)
{
    Dot3D absOfs = flag.getLocation();
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint x = iX + absOfs.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint y = iY + absOfs.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint z = iZ + absOfs.z;
                if (freeSurfaceFlag::isWet(flag.get(iX, iY, iZ))) {
                    // BlockStatistics computes only maximum, no minimum. Therefore,
                    //   the relation min(x) = -max(-x) is used.
                    statistics.gatherMax(minIdX, (double)-x);
                    statistics.gatherMax(minIdY, (double)-y);
                    statistics.gatherMax(minIdZ, (double)-z);
                    statistics.gatherMax(maxIdX, (double)x);
                    statistics.gatherMax(maxIdY, (double)y);
                    statistics.gatherMax(maxIdZ, (double)z);
                }
            }
        }
    }
}

FreeSurfaceComputeFluidBoundingBox3D *FreeSurfaceComputeFluidBoundingBox3D::clone() const
{
    return new FreeSurfaceComputeFluidBoundingBox3D(*this);
}

Box3D FreeSurfaceComputeFluidBoundingBox3D::getFluidBoundingBox() const
{
    // The minus sign accounts for the relation min(x) = -max(-x).
    double doubleMinX = -this->getStatistics().getMax(minIdX);
    double doubleMinY = -this->getStatistics().getMax(minIdY);
    double doubleMinZ = -this->getStatistics().getMax(minIdZ);
    double doubleMaxX = this->getStatistics().getMax(maxIdX);
    double doubleMaxY = this->getStatistics().getMax(maxIdY);
    double doubleMaxZ = this->getStatistics().getMax(maxIdZ);

    plint minX = util::roundToInt(doubleMinX);
    plint minY = util::roundToInt(doubleMinY);
    plint minZ = util::roundToInt(doubleMinZ);
    plint maxX = util::roundToInt(doubleMaxX);
    plint maxY = util::roundToInt(doubleMaxY);
    plint maxZ = util::roundToInt(doubleMaxZ);

    return Box3D(minX, maxX, minY, maxY, minZ, maxZ);
}

Box3D freeSurfaceComputeFluidBoundingBox(MultiScalarField3D<int> &flag, Box3D domain)
{
    FreeSurfaceComputeFluidBoundingBox3D functional;
    applyProcessingFunctional(functional, domain, flag);
    return functional.getFluidBoundingBox();
}

Box3D freeSurfaceComputeFluidBoundingBox(MultiScalarField3D<int> &flag)
{
    return freeSurfaceComputeFluidBoundingBox(flag, flag.getBoundingBox());
}

}  // namespace plb
