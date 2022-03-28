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

#ifndef FREE_SURFACE_ANALYSIS_3D_H
#define FREE_SURFACE_ANALYSIS_3D_H

#include <memory>

#include "core/globalDefs.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
T freeSurfaceAverageMass(std::vector<MultiBlock3D *> freeSurfaceArgs, Box3D domain);

template <typename T, template <typename U> class Descriptor>
class FS_AverageMassFunctional3D : public PlainReductiveBoxProcessingFunctional3D {
public:
    FS_AverageMassFunctional3D();
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual FS_AverageMassFunctional3D<T, Descriptor> *clone() const;
    T getAverageMass() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        for (pluint i = 0; i < modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }

private:
    plint averageMassId;
};

template <typename T, template <typename U> class Descriptor>
T freeSurfaceTotalMass(std::vector<MultiBlock3D *> freeSurfaceArgs, Box3D domain);

template <typename T, template <typename U> class Descriptor>
class FS_TotalMassFunctional3D : public PlainReductiveBoxProcessingFunctional3D {
public:
    FS_TotalMassFunctional3D();
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual FS_TotalMassFunctional3D<T, Descriptor> *clone() const;
    T getTotalMass() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        for (pluint i = 0; i < modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }

private:
    plint totalMassId;
};

template <typename T, template <typename U> class Descriptor>
T freeSurfaceAverageDensity(std::vector<MultiBlock3D *> freeSurfaceArgs, Box3D domain);

template <typename T, template <typename U> class Descriptor>
class FS_AverageDensityFunctional3D : public PlainReductiveBoxProcessingFunctional3D {
public:
    FS_AverageDensityFunctional3D();
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual FS_AverageDensityFunctional3D<T, Descriptor> *clone() const;
    T getAverageDensity() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        for (pluint i = 0; i < modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }

private:
    plint averageDensityId;
};

template <typename T, template <typename U> class Descriptor>
T freeSurfaceAverageVolumeFraction(std::vector<MultiBlock3D *> freeSurfaceArgs, Box3D domain);

template <typename T, template <typename U> class Descriptor>
class FS_AverageVolumeFractionFunctional3D : public PlainReductiveBoxProcessingFunctional3D {
public:
    FS_AverageVolumeFractionFunctional3D();
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual FS_AverageVolumeFractionFunctional3D<T, Descriptor> *clone() const;
    T getAverageVolumeFraction() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        for (pluint i = 0; i < modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }

private:
    plint averageVfId;
};

template <typename T, template <typename U> class Descriptor>
plint countFreeSurfaceElements(
    std::vector<MultiBlock3D *> freeSurfaceArgs, plint flagToLookFor, Box3D domain);

template <typename T, template <typename U> class Descriptor>
class CountFreeSurfaceElementsFunctional3D : public PlainReductiveBoxProcessingFunctional3D {
public:
    CountFreeSurfaceElementsFunctional3D(plint flagToLookFor_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual CountFreeSurfaceElementsFunctional3D<T, Descriptor> *clone() const;
    plint getNumElements() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        for (pluint i = 0; i < modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }

private:
    plint numCellsId;
    plint flagToLookFor;
};

template <typename T, template <typename U> class Descriptor>
Array<T, 3> freeSurfaceAverageMomentum(std::vector<MultiBlock3D *> freeSurfaceArgs, Box3D domain);

template <typename T, template <typename U> class Descriptor>
class FS_AverageMomentumFunctional3D : public PlainReductiveBoxProcessingFunctional3D {
public:
    FS_AverageMomentumFunctional3D();
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual FS_AverageMomentumFunctional3D<T, Descriptor> *clone() const;
    Array<T, 3> getAverageMomentum() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        for (pluint i = 0; i < modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }

private:
    Array<plint, 3> averageMomentumId;
};

template <typename T, template <typename U> class Descriptor>
T freeSurfaceAverageHeight(std::vector<MultiBlock3D *> freeSurfaceArgs, Box3D domain);

template <typename T, template <typename U> class Descriptor>
class FS_AverageHeightFunctional3D : public PlainReductiveBoxProcessingFunctional3D {
public:
    FS_AverageHeightFunctional3D();
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual FS_AverageHeightFunctional3D<T, Descriptor> *clone() const;
    T getAverageHeight() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        for (pluint i = 0; i < modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }

private:
    plint averageHeightId;
};

template <typename T, template <typename U> class Descriptor>
T getAverageHeightAtXY(std::vector<MultiBlock3D *> freeSurfaceArgs, plint N, Box3D domain);

template <typename T, template <typename U> class Descriptor>
class GetWaterLevelAtxyFunctional3D : public PlainReductiveBoxProcessingFunctional3D {
public:
    GetWaterLevelAtxyFunctional3D();
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual GetWaterLevelAtxyFunctional3D<T, Descriptor> *clone() const;
    plint getNumFluidCellsAtXY() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        for (pluint i = 0; i < modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }

private:
    plint numFluidOccupiedCellId;
};

// Here we count both the no-slip and the free-slip wall cells.
class FreeSurfaceCountWallCells3D : public ReductiveBoxProcessingFunctional3D_S<int> {
public:
    FreeSurfaceCountWallCells3D();
    virtual void process(Box3D domain, ScalarField3D<int> &flag);
    virtual FreeSurfaceCountWallCells3D *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;  // Flags.
    }
    plint getNumWallCells() const;

private:
    plint sumScalarId;
};

plint freeSurfaceCountWallCells(MultiScalarField3D<int> &flag, Box3D domain);

template <typename T>
T freeSurfaceComputePorosity(MultiScalarField3D<int> &flag, Box3D domain);

template <typename T>
class FreeSurfaceComputeFluidVolume3D : public ReductiveBoxProcessingFunctional3D_SS<T, int> {
public:
    FreeSurfaceComputeFluidVolume3D();
    virtual void process(Box3D domain, ScalarField3D<T> &volumeFraction, ScalarField3D<int> &flag);
    virtual FreeSurfaceComputeFluidVolume3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;  // Volume fraction.
        modified[1] = modif::nothing;  // Flags.
    }
    T getFluidVolume() const;

private:
    plint sumScalarId;
};

template <typename T>
class MaskedFreeSurfaceComputeFluidVolume3D : public PlainReductiveBoxProcessingFunctional3D {
public:
    MaskedFreeSurfaceComputeFluidVolume3D();
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual MaskedFreeSurfaceComputeFluidVolume3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;  // Volume fraction.
        modified[1] = modif::nothing;  // Flags.
        modified[2] = modif::nothing;  // Mask.
    }
    T getFluidVolume() const;

private:
    plint sumScalarId;
};

template <typename T>
T freeSurfaceComputeFluidVolume(
    MultiScalarField3D<T> &volumeFraction, MultiScalarField3D<int> &flag, Box3D domain);

template <typename T>
T freeSurfaceComputeFluidVolume(
    MultiScalarField3D<T> &volumeFraction, MultiScalarField3D<int> &flag,
    MultiScalarField3D<int> &mask);

// The porosity must have been computed inside the exact same domain as the argument of this
// function. Usually the porosity is constant in time, this is why we pass it as an argument, and do
// not recompute it inside this function.
template <typename T>
T freeSurfaceComputeSaturation(
    T porosity, MultiScalarField3D<T> &volumeFraction, MultiScalarField3D<int> &flag, Box3D domain);

// Here we find the smallest Box3D that includes all fluid and interface cells
class FreeSurfaceComputeFluidBoundingBox3D : public ReductiveBoxProcessingFunctional3D_S<int> {
public:
    FreeSurfaceComputeFluidBoundingBox3D();
    virtual void process(Box3D domain, ScalarField3D<int> &flag);
    virtual FreeSurfaceComputeFluidBoundingBox3D *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;  // Flags.
    }
    Box3D getFluidBoundingBox() const;

private:
    plint minIdX, minIdY, minIdZ;
    plint maxIdX, maxIdY, maxIdZ;
};

Box3D freeSurfaceComputeFluidBoundingBox(MultiScalarField3D<int> &flag, Box3D domain);
Box3D freeSurfaceComputeFluidBoundingBox(MultiScalarField3D<int> &flag);

template <typename T, template <typename U> class Descriptor>
class FreeSurfaceBoxSumForcedEnergyFunctional3D : public PlainReductiveBoxProcessingFunctional3D {
public:
    FreeSurfaceBoxSumForcedEnergyFunctional3D();
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual FreeSurfaceBoxSumForcedEnergyFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    T getAverageEnergy() const;

private:
    plint sumEnergyId;
    plint sumCellsId;
};

template <typename T, template <typename U> class Descriptor>
T freeSurfaceComputeAverageForcedEnergy(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<int> &flag, Box3D domain);

template <typename T, template <typename U> class Descriptor>
T freeSurfaceComputeAverageForcedEnergy(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<int> &flag);

template <typename T, template <typename U> class Descriptor>
void freeSurfaceComputeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<int> &flag, MultiScalarField3D<T> &velocityNorm, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > freeSurfaceComputeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<int> &flag, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > freeSurfaceComputeForcedVelocityNorm(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<int> &flag);

template <typename T, template <typename U> class Descriptor>
void freeSurfaceComputeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<int> &flag, MultiScalarField3D<T> &velocityComponent, Box3D domain,
    plint iComponent);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > freeSurfaceComputeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<int> &flag, Box3D domain, plint iComponent);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > freeSurfaceComputeForcedVelocityComponent(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<int> &flag, plint iComponent);

template <typename T, template <typename U> class Descriptor>
void freeSurfaceComputeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<int> &flag, MultiTensorField3D<T, Descriptor<T>::d> &velocity, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > freeSurfaceComputeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<int> &flag, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > freeSurfaceComputeForcedVelocity(
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &force,
    MultiScalarField3D<int> &flag);

}  // namespace plb

#endif  // FREE_SURFACE_ANALYSIS_3D_H
