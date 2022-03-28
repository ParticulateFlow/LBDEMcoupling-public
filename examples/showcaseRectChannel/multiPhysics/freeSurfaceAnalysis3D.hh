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

#ifndef FREE_SURFACE_ANALYSIS_3D_HH
#define FREE_SURFACE_ANALYSIS_3D_HH

#include "multiPhysics/freeSurfaceAnalysis3D.h"
#include "multiPhysics/freeSurfaceUtil3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
FS_AverageMassFunctional3D<T, Descriptor>::FS_AverageMassFunctional3D() :
    averageMassId(this->getStatistics().subscribeAverage())
{ }

template <typename T, template <typename U> class Descriptor>
void FS_AverageMassFunctional3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);
    BlockStatistics &statistics = this->getStatistics();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                statistics.gatherAverage(averageMassId, param.mass(iX, iY, iZ));
                statistics.incrementStats();
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
FS_AverageMassFunctional3D<T, Descriptor> *FS_AverageMassFunctional3D<T, Descriptor>::clone() const
{
    return new FS_AverageMassFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
T FS_AverageMassFunctional3D<T, Descriptor>::getAverageMass() const
{
    return this->getStatistics().getAverage(averageMassId);
}

template <typename T, template <typename U> class Descriptor>
T freeSurfaceAverageMass(std::vector<MultiBlock3D *> freeSurfaceArgs, Box3D domain)
{
    FS_AverageMassFunctional3D<T, Descriptor> functional;
    applyProcessingFunctional(functional, domain, freeSurfaceArgs);
    return functional.getAverageMass();
}

template <typename T, template <typename U> class Descriptor>
FS_TotalMassFunctional3D<T, Descriptor>::FS_TotalMassFunctional3D() :
    totalMassId(this->getStatistics().subscribeSum())
{ }

template <typename T, template <typename U> class Descriptor>
void FS_TotalMassFunctional3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);
    BlockStatistics &statistics = this->getStatistics();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                statistics.gatherSum(totalMassId, param.mass(iX, iY, iZ));
                statistics.incrementStats();
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
FS_TotalMassFunctional3D<T, Descriptor> *FS_TotalMassFunctional3D<T, Descriptor>::clone() const
{
    return new FS_TotalMassFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
T FS_TotalMassFunctional3D<T, Descriptor>::getTotalMass() const
{
    return this->getStatistics().getSum(totalMassId);
}

template <typename T, template <typename U> class Descriptor>
T freeSurfaceTotalMass(std::vector<MultiBlock3D *> freeSurfaceArgs, Box3D domain)
{
    FS_TotalMassFunctional3D<T, Descriptor> functional;
    applyProcessingFunctional(functional, domain, freeSurfaceArgs);
    return functional.getTotalMass();
}

template <typename T, template <typename U> class Descriptor>
FS_AverageDensityFunctional3D<T, Descriptor>::FS_AverageDensityFunctional3D() :
    averageDensityId(this->getStatistics().subscribeAverage())
{ }

template <typename T, template <typename U> class Descriptor>
void FS_AverageDensityFunctional3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);
    BlockStatistics &statistics = this->getStatistics();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (param.cell(iX, iY, iZ).getDynamics().hasMoments()) {
                    statistics.gatherAverage(averageDensityId, param.getDensity(iX, iY, iZ));
                    statistics.incrementStats();
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
FS_AverageDensityFunctional3D<T, Descriptor> *FS_AverageDensityFunctional3D<T, Descriptor>::clone()
    const
{
    return new FS_AverageDensityFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
T FS_AverageDensityFunctional3D<T, Descriptor>::getAverageDensity() const
{
    return this->getStatistics().getAverage(averageDensityId);
}

template <typename T, template <typename U> class Descriptor>
T freeSurfaceAverageDensity(std::vector<MultiBlock3D *> freeSurfaceArgs, Box3D domain)
{
    FS_AverageDensityFunctional3D<T, Descriptor> functional;
    applyProcessingFunctional(functional, domain, freeSurfaceArgs);
    return functional.getAverageDensity();
}

template <typename T, template <typename U> class Descriptor>
FS_AverageVolumeFractionFunctional3D<T, Descriptor>::FS_AverageVolumeFractionFunctional3D() :
    averageVfId(this->getStatistics().subscribeAverage())
{ }

template <typename T, template <typename U> class Descriptor>
void FS_AverageVolumeFractionFunctional3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);
    BlockStatistics &statistics = this->getStatistics();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (param.cell(iX, iY, iZ).getDynamics().hasMoments()) {
                    statistics.gatherAverage(averageVfId, param.volumeFraction(iX, iY, iZ));
                    statistics.incrementStats();
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
FS_AverageVolumeFractionFunctional3D<T, Descriptor>
    *FS_AverageVolumeFractionFunctional3D<T, Descriptor>::clone() const
{
    return new FS_AverageVolumeFractionFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
T FS_AverageVolumeFractionFunctional3D<T, Descriptor>::getAverageVolumeFraction() const
{
    return this->getStatistics().getAverage(averageVfId);
}

template <typename T, template <typename U> class Descriptor>
T freeSurfaceAverageVolumeFraction(std::vector<MultiBlock3D *> freeSurfaceArgs, Box3D domain)
{
    FS_AverageVolumeFractionFunctional3D<T, Descriptor> functional;
    applyProcessingFunctional(functional, domain, freeSurfaceArgs);
    return functional.getAverageVolumeFraction();
}

template <typename T, template <typename U> class Descriptor>
CountFreeSurfaceElementsFunctional3D<T, Descriptor>::CountFreeSurfaceElementsFunctional3D(
    plint flagToLookFor_) :
    numCellsId(this->getStatistics().subscribeIntSum()), flagToLookFor(flagToLookFor_)
{ }

template <typename T, template <typename U> class Descriptor>
void CountFreeSurfaceElementsFunctional3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);
    BlockStatistics &statistics = this->getStatistics();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                int materialIndex = param.flag(iX, iY, iZ);
                if (materialIndex == flagToLookFor) {  // Fluid Cell
                    statistics.gatherIntSum(numCellsId, 1);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
CountFreeSurfaceElementsFunctional3D<T, Descriptor>
    *CountFreeSurfaceElementsFunctional3D<T, Descriptor>::clone() const
{
    return new CountFreeSurfaceElementsFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
plint CountFreeSurfaceElementsFunctional3D<T, Descriptor>::getNumElements() const
{
    return this->getStatistics().getIntSum(numCellsId);
}

template <typename T, template <typename U> class Descriptor>
plint countFreeSurfaceElements(
    std::vector<MultiBlock3D *> freeSurfaceArgs, plint flagToLookFor, Box3D domain)
{
    CountFreeSurfaceElementsFunctional3D<T, Descriptor> functional(flagToLookFor);
    applyProcessingFunctional(functional, domain, freeSurfaceArgs);
    return functional.getNumElements();
}

template <typename T, template <typename U> class Descriptor>
FS_AverageMomentumFunctional3D<T, Descriptor>::FS_AverageMomentumFunctional3D() :
    averageMomentumId(
        this->getStatistics().subscribeAverage(), this->getStatistics().subscribeAverage(),
        this->getStatistics().subscribeAverage())
{ }

template <typename T, template <typename U> class Descriptor>
void FS_AverageMomentumFunctional3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);
    BlockStatistics &statistics = this->getStatistics();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (param.cell(iX, iY, iZ).getDynamics().hasMoments()) {
                    Array<T, Descriptor<T>::d> j = param.getMomentum(iX, iY, iZ);
                    statistics.gatherAverage(averageMomentumId[0], j[0]);
                    statistics.gatherAverage(averageMomentumId[1], j[1]);
                    statistics.gatherAverage(averageMomentumId[2], j[2]);
                    statistics.incrementStats();
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
FS_AverageMomentumFunctional3D<T, Descriptor>
    *FS_AverageMomentumFunctional3D<T, Descriptor>::clone() const
{
    return new FS_AverageMomentumFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
Array<T, 3> FS_AverageMomentumFunctional3D<T, Descriptor>::getAverageMomentum() const
{
    return Array<T, 3>(
        this->getStatistics().getAverage(averageMomentumId[0]),
        this->getStatistics().getAverage(averageMomentumId[1]),
        this->getStatistics().getAverage(averageMomentumId[2]));
}

template <typename T, template <typename U> class Descriptor>
Array<T, 3> freeSurfaceAverageMomentum(std::vector<MultiBlock3D *> freeSurfaceArgs, Box3D domain)
{
    FS_AverageMomentumFunctional3D<T, Descriptor> functional;
    applyProcessingFunctional(functional, domain, freeSurfaceArgs);
    return functional.getAverageMomentum();
}

template <typename T, template <typename U> class Descriptor>
FS_AverageHeightFunctional3D<T, Descriptor>::FS_AverageHeightFunctional3D() :
    averageHeightId(this->getStatistics().subscribeAverage())
{ }

template <typename T, template <typename U> class Descriptor>
void FS_AverageHeightFunctional3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);
    BlockStatistics &statistics = this->getStatistics();
    Dot3D absOffset = param.absOffset();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (param.cell(iX, iY, iZ).getDynamics().hasMoments()) {
                    T localHeight = T(0);
                    if (util::greaterEqual(param.volumeFraction(iX, iY, iZ), (T)1)) {
                        localHeight = T(absOffset.z);
                    }
                    statistics.gatherAverage(averageHeightId, localHeight);
                    statistics.incrementStats();
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
FS_AverageHeightFunctional3D<T, Descriptor> *FS_AverageHeightFunctional3D<T, Descriptor>::clone()
    const
{
    return new FS_AverageHeightFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
T FS_AverageHeightFunctional3D<T, Descriptor>::getAverageHeight() const
{
    return this->getStatistics().getAverage(averageHeightId);
}

template <typename T, template <typename U> class Descriptor>
T freeSurfaceAverageHeight(std::vector<MultiBlock3D *> freeSurfaceArgs, Box3D domain)
{
    FS_AverageHeightFunctional3D<T, Descriptor> functional;
    applyProcessingFunctional(functional, domain, freeSurfaceArgs);
    return functional.getAverageHeight();
}

template <typename T, template <typename U> class Descriptor>
GetWaterLevelAtxyFunctional3D<T, Descriptor>::GetWaterLevelAtxyFunctional3D() :
    numFluidOccupiedCellId(this->getStatistics().subscribeIntSum())
{ }

template <typename T, template <typename U> class Descriptor>
void GetWaterLevelAtxyFunctional3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);
    BlockStatistics &statistics = this->getStatistics();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (param.cell(iX, iY, iZ).getDynamics().hasMoments()) {
                    if (param.volumeFraction(iX, iY, iZ) >= 0.5) {
                        statistics.gatherIntSum(numFluidOccupiedCellId, 1);
                    }
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
GetWaterLevelAtxyFunctional3D<T, Descriptor> *GetWaterLevelAtxyFunctional3D<T, Descriptor>::clone()
    const
{
    return new GetWaterLevelAtxyFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
plint GetWaterLevelAtxyFunctional3D<T, Descriptor>::getNumFluidCellsAtXY() const
{
    return this->getStatistics().getIntSum(numFluidOccupiedCellId);
}

template <typename T, template <typename U> class Descriptor>
T getAverageHeightAtXY(std::vector<MultiBlock3D *> freeSurfaceArgs, plint N, Box3D domain)
{
    GetWaterLevelAtxyFunctional3D<T, Descriptor> functional;
    applyProcessingFunctional(functional, domain, freeSurfaceArgs);
    plint length_domain = domain.x1 - domain.x0;  // number of cell along y direction
    if (length_domain == 0)
        length_domain = 1;
    plint width_domain = domain.y1 - domain.y0;  // number of cell along y direction
    if (width_domain == 0)
        width_domain = 1;
    T heightAtXY = functional.getNumFluidCellsAtXY() / (T(N) * length_domain * width_domain);
    return heightAtXY;
}

template <typename T>
T freeSurfaceComputePorosity(MultiScalarField3D<int> &flag, Box3D domain)
{
    plint numWallCells = freeSurfaceCountWallCells(flag, domain);
    T porosity = (T)1 - (T)numWallCells / (T)domain.nCells();
    return porosity;
}

template <typename T>
FreeSurfaceComputeFluidVolume3D<T>::FreeSurfaceComputeFluidVolume3D() :
    sumScalarId(this->getStatistics().subscribeSum())
{ }

template <typename T>
void FreeSurfaceComputeFluidVolume3D<T>::process(
    Box3D domain, ScalarField3D<T> &volumeFraction, ScalarField3D<int> &flag)
{
    Dot3D offset = computeRelativeDisplacement(volumeFraction, flag);
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                int localFlag = flag.get(iX + offset.x, iY + offset.y, iZ + offset.z);
                if (freeSurfaceFlag::isFullWet(localFlag)) {
                    statistics.gatherSum(sumScalarId, (double)1.0);
                } else if (localFlag == freeSurfaceFlag::interface) {
                    statistics.gatherSum(sumScalarId, (double)volumeFraction.get(iX, iY, iZ));
                }
            }
        }
    }
}

template <typename T>
FreeSurfaceComputeFluidVolume3D<T> *FreeSurfaceComputeFluidVolume3D<T>::clone() const
{
    return new FreeSurfaceComputeFluidVolume3D<T>(*this);
}

template <typename T>
T FreeSurfaceComputeFluidVolume3D<T>::getFluidVolume() const
{
    double fluidVolume = this->getStatistics().getSum(sumScalarId);
    // The sum is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    if (std::numeric_limits<T>::is_integer) {
        return (T)util::roundToInt(fluidVolume);
    }
    return (T)fluidVolume;
}

template <typename T>
T freeSurfaceComputeFluidVolume(
    MultiScalarField3D<T> &volumeFraction, MultiScalarField3D<int> &flag, Box3D domain)
{
    FreeSurfaceComputeFluidVolume3D<T> functional;
    applyProcessingFunctional(functional, domain, volumeFraction, flag);
    return functional.getFluidVolume();
}

template <typename T>
MaskedFreeSurfaceComputeFluidVolume3D<T>::MaskedFreeSurfaceComputeFluidVolume3D() :
    sumScalarId(this->getStatistics().subscribeSum())
{ }

template <typename T>
void MaskedFreeSurfaceComputeFluidVolume3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    ScalarField3D<T> &volumeFraction = *dynamic_cast<ScalarField3D<T> *>(fields[0]);
    ScalarField3D<int> &flag = *dynamic_cast<ScalarField3D<int> *>(fields[1]);
    ScalarField3D<int> &mask = *dynamic_cast<ScalarField3D<int> *>(fields[2]);
    Dot3D offsetFlag = computeRelativeDisplacement(volumeFraction, flag);
    Dot3D offsetMask = computeRelativeDisplacement(volumeFraction, mask);
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (mask.get(iX + offsetMask.x, iY + offsetMask.y, iZ + offsetMask.z)) {
                    int localFlag =
                        flag.get(iX + offsetFlag.x, iY + offsetFlag.y, iZ + offsetFlag.z);
                    if (freeSurfaceFlag::isFullWet(localFlag)) {
                        statistics.gatherSum(sumScalarId, (double)1.0);
                    } else if (localFlag == freeSurfaceFlag::interface) {
                        statistics.gatherSum(sumScalarId, (double)volumeFraction.get(iX, iY, iZ));
                    }
                }
            }
        }
    }
}

template <typename T>
MaskedFreeSurfaceComputeFluidVolume3D<T> *MaskedFreeSurfaceComputeFluidVolume3D<T>::clone() const
{
    return new MaskedFreeSurfaceComputeFluidVolume3D<T>(*this);
}

template <typename T>
T MaskedFreeSurfaceComputeFluidVolume3D<T>::getFluidVolume() const
{
    double fluidVolume = this->getStatistics().getSum(sumScalarId);
    // The sum is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    if (std::numeric_limits<T>::is_integer) {
        return (T)util::roundToInt(fluidVolume);
    }
    return (T)fluidVolume;
}

template <typename T>
T freeSurfaceComputeFluidVolume(
    MultiScalarField3D<T> &volumeFraction, MultiScalarField3D<int> &flag,
    MultiScalarField3D<int> &mask)
{
    MaskedFreeSurfaceComputeFluidVolume3D<T> functional;
    std::vector<MultiBlock3D *> args;
    args.push_back(&volumeFraction);
    args.push_back(&flag);
    args.push_back(&mask);
    applyProcessingFunctional(functional, volumeFraction.getBoundingBox(), args);
    return functional.getFluidVolume();
}

template <typename T>
T freeSurfaceComputeSaturation(
    T porosity, MultiScalarField3D<T> &volumeFraction, MultiScalarField3D<int> &flag, Box3D domain)
{
    T fluidVolume = freeSurfaceComputeFluidVolume(volumeFraction, flag, domain);
    T totalVolume = domain.nCells();
    T saturation = fluidVolume / (porosity * totalVolume);
    return saturation;
}

template <typename T, template <typename U> class Descriptor>
FreeSurfaceBoxSumForcedEnergyFunctional3D<
    T, Descriptor>::FreeSurfaceBoxSumForcedEnergyFunctional3D() :
    sumEnergyId(this->getStatistics().subscribeSum()),
    sumCellsId(this->getStatistics().subscribeSum())
{ }

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceBoxSumForcedEnergyFunctional3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 3);
    BlockLattice3D<T, Descriptor> *lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> *>(blocks[0]);
    PLB_ASSERT(lattice);
    TensorField3D<T, Descriptor<T>::d> *force =
        dynamic_cast<TensorField3D<T, Descriptor<T>::d> *>(blocks[1]);
    PLB_ASSERT(force);
    ScalarField3D<int> *flag = dynamic_cast<ScalarField3D<int> *>(blocks[2]);
    PLB_ASSERT(flag);

    Dot3D ofsForce = computeRelativeDisplacement(*lattice, *force);
    Dot3D ofsFlag = computeRelativeDisplacement(*lattice, *flag);

    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (freeSurfaceFlag::isWet(
                        flag->get(iX + ofsFlag.x, iY + ofsFlag.y, iZ + ofsFlag.z))) {
                    Array<T, Descriptor<T>::d> velocity;
                    lattice->get(iX, iY, iZ).computeVelocity(velocity);
                    Array<T, Descriptor<T>::d> const &f =
                        force->get(iX + ofsForce.x, iY + ofsForce.y, iZ + ofsForce.z);
                    velocity[0] += (T)0.5 * f[0];
                    velocity[1] += (T)0.5 * f[1];
                    velocity[2] += (T)0.5 * f[2];
                    T uNormSqr = VectorTemplate<T, Descriptor>::normSqr(velocity);
                    statistics.gatherSum(sumEnergyId, uNormSqr);
                    statistics.gatherSum(sumCellsId, (T)1);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
FreeSurfaceBoxSumForcedEnergyFunctional3D<T, Descriptor>
    *FreeSurfaceBoxSumForcedEnergyFunctional3D<T, Descriptor>::clone() const
{
    return new FreeSurfaceBoxSumForcedEnergyFunctional3D(*this);
}

template <typename T, template <typename U> class Descriptor>
T FreeSurfaceBoxSumForcedEnergyFunctional3D<T, Descriptor>::getAverageEnergy() const
{
    T sumEnergy = this->getStatistics().getSum(sumEnergyId) / (T)2;
    T sumCells = this->getStatistics().getSum(sumCellsId);
    T averageEnergy = (T)0;
    if (!util::isZero(sumCells)) {
        averageEnergy = sumEnergy / sumCells;
    }
    return averageEnergy;
}

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceBoxSumForcedEnergyFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;  // Lattice.
    modified[1] = modif::nothing;  // Force.
    modified[2] = modif::nothing;  // Flag.
}

template <typename T, template <typename U> class Descriptor>
T freeSurfaceComputeAverageForcedEnergy(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<int> &flag, Box3D domain)
{
    FreeSurfaceBoxSumForcedEnergyFunctional3D<T, Descriptor> functional;
    std::vector<MultiBlock3D *> args;
    args.push_back(&lattice);
    args.push_back(&force);
    args.push_back(&flag);
    applyProcessingFunctional(functional, domain, args);
    return functional.getAverageEnergy();
}

template <typename T, template <typename U> class Descriptor>
T freeSurfaceComputeAverageForcedEnergy(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<int> &flag)
{
    return freeSurfaceComputeAverageForcedEnergy(lattice, force, flag, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void freeSurfaceComputeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<int> &flag, MultiScalarField3D<T> &velocityNorm, Box3D domain)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&lattice);
    args.push_back(&force);
    args.push_back(&velocityNorm);
    applyProcessingFunctional(new BoxForcedVelocityNormFunctional3D<T, Descriptor>, domain, args);

    setToConstant<T>(velocityNorm, flag, (int)freeSurfaceFlag::empty, domain, (T)0);
    setToConstant<T>(velocityNorm, flag, (int)freeSurfaceFlag::protectEmpty, domain, (T)0);
    setToConstant<T>(velocityNorm, flag, (int)freeSurfaceFlag::wall, domain, (T)0);
    setToConstant<T>(velocityNorm, flag, (int)freeSurfaceFlag::slipWall, domain, (T)0);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > freeSurfaceComputeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<int> &flag, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > velocityNorm =
        generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxForcedVelocityNormFunctional3D() acts on both bulk and envelope, you would expect the
    // envelope layer around the domain, on the velocityNorm multi-block, to be assigned some proper
    // values too. By default, this is however not what happens, because the physical space occupied
    // by these envelopes does not intersect with the domain "domain". We work around this issue by
    // extending the domain. There's no problem if the enlarged domain gets beyond the actual extent
    // of the lattice, because Palabos handles these situations properly.

    freeSurfaceComputeForcedVelocityNorm(
        lattice, force, flag, *velocityNorm,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return velocityNorm;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > freeSurfaceComputeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<int> &flag)
{
    return freeSurfaceComputeForcedVelocityNorm(lattice, force, flag, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void freeSurfaceComputeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<int> &flag, MultiScalarField3D<T> &velocityComponent, Box3D domain,
    plint iComponent)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&lattice);
    args.push_back(&force);
    args.push_back(&velocityComponent);
    applyProcessingFunctional(
        new BoxForcedVelocityComponentFunctional3D<T, Descriptor>(iComponent), domain, args);

    setToConstant<T>(velocityComponent, flag, (int)freeSurfaceFlag::empty, domain, (T)0);
    setToConstant<T>(velocityComponent, flag, (int)freeSurfaceFlag::protectEmpty, domain, (T)0);
    setToConstant<T>(velocityComponent, flag, (int)freeSurfaceFlag::wall, domain, (T)0);
    setToConstant<T>(velocityComponent, flag, (int)freeSurfaceFlag::slipWall, domain, (T)0);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > freeSurfaceComputeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<int> &flag, Box3D domain, plint iComponent)
{
    std::unique_ptr<MultiScalarField3D<T> > velocityComponent =
        generateMultiScalarField<T>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxForcedVelocityComponentFunctional3D() acts on both bulk and envelope, you would expect the
    // envelope layer around the domain, on the velocityComponent multi-block, to be assigned some
    // proper values too. By default, this is however not what happens, because the physical space
    // occupied by these envelopes does not intersect with the domain "domain". We work around this
    // issue by extending the domain. There's no problem if the enlarged domain gets beyond the
    // actual extent of the lattice, because Palabos handles these situations properly.

    freeSurfaceComputeForcedVelocityComponent(
        lattice, force, flag, *velocityComponent,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()), iComponent);
    return velocityComponent;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > freeSurfaceComputeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<int> &flag, plint iComponent)
{
    return freeSurfaceComputeForcedVelocityComponent(
        lattice, force, flag, lattice.getBoundingBox(), iComponent);
}

template <typename T, template <typename U> class Descriptor>
void freeSurfaceComputeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<int> &flag, MultiTensorField3D<T, Descriptor<T>::d> &velocity, Box3D domain)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&lattice);
    args.push_back(&force);
    args.push_back(&velocity);
    applyProcessingFunctional(new BoxForcedVelocityFunctional3D<T, Descriptor>, domain, args);

    setToConstant<T, Descriptor<T>::d>(
        velocity, flag, (int)freeSurfaceFlag::empty, domain, Array<T, Descriptor<T>::d>::zero());
    setToConstant<T, Descriptor<T>::d>(
        velocity, flag, (int)freeSurfaceFlag::protectEmpty, domain,
        Array<T, Descriptor<T>::d>::zero());
    setToConstant<T, Descriptor<T>::d>(
        velocity, flag, (int)freeSurfaceFlag::wall, domain, Array<T, Descriptor<T>::d>::zero());
    setToConstant<T, Descriptor<T>::d>(
        velocity, flag, (int)freeSurfaceFlag::slipWall, domain, Array<T, Descriptor<T>::d>::zero());
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > freeSurfaceComputeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<int> &flag, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > velocity =
        generateMultiTensorField<T, Descriptor<T>::d>(lattice, domain);

    // The domain needs to be extended to the outer envelopes, for the following reason. Imagine
    // that the domain is smaller than the bounding-box. Given that the
    // BoxForcedVelocityFunctional3D() acts on both bulk and envelope, you would expect the envelope
    // layer around the domain, on the velocity multi-block, to be assigned some proper values too.
    // By default, this is however not what happens, because the physical space occupied by these
    // envelopes does not intersect with the domain "domain". We work around this issue by extending
    // the domain. There's no problem if the enlarged domain gets beyond the actual extent of the
    // lattice, because Palabos handles these situations properly.

    freeSurfaceComputeForcedVelocity(
        lattice, force, flag, *velocity,
        domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
    return velocity;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > freeSurfaceComputeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<int> &flag)
{
    return freeSurfaceComputeForcedVelocity(lattice, force, flag, lattice.getBoundingBox());
}

}  // namespace plb

#endif  // FREE_SURFACE_ANALYSIS_3D_HH
