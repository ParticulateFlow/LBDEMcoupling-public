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

#ifndef OFF_LATTICE_BOUNDARY_CONDITION_3D_H
#define OFF_LATTICE_BOUNDARY_CONDITION_3D_H

#include "core/globalDefs.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "offLattice/guoOffLatticeModel3D.h"
#include "offLattice/offLatticeBoundaryProfiles3D.h"
#include "offLattice/triangleBoundary3D.h"
#include "offLattice/triangleToDef.h"
#include "offLattice/triangularSurfaceMesh.h"

namespace plb {

template <typename T>
struct OffLatticeBoundaryConditionBase3D {
    virtual ~OffLatticeBoundaryConditionBase3D() { }
    virtual void apply() = 0;
    virtual void insert(plint processorLevel) = 0;
    virtual void apply(std::vector<MultiBlock3D *> const &completionArg) = 0;
    virtual void insert(std::vector<MultiBlock3D *> const &completionArg, plint processorLevel) = 0;
    virtual Array<T, 3> getForceOnObject() = 0;
};

template <typename T, template <typename U> class Descriptor, class BoundaryType>
class OffLatticeBoundaryCondition3D : public OffLatticeBoundaryConditionBase3D<T> {
public:
    OffLatticeBoundaryCondition3D(
        OffLatticeModel3D<T, BoundaryType> *offLatticeModel_,
        VoxelizedDomain3D<T> &voxelizedDomain_, MultiBlockLattice3D<T, Descriptor> &lattice_);
    OffLatticeBoundaryCondition3D(
        OffLatticeModel3D<T, BoundaryType> *offLatticeModel_,
        VoxelizedDomain3D<T> &voxelizedDomain_, MultiBlockLattice3D<T, Descriptor> &lattice_,
        MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particleField_);
    OffLatticeBoundaryCondition3D(
        OffLatticeBoundaryCondition3D<T, Descriptor, BoundaryType> const &rhs);
    ~OffLatticeBoundaryCondition3D();
    MultiBlockLattice3D<T, Descriptor> const &getLattice() const
    {
        return lattice;
    }
    VoxelizedDomain3D<T> const &getVoxelizedDomain() const
    {
        return voxelizedDomain;
    }
    VoxelizedDomain3D<T> &getVoxelizedDomain()
    {
        return voxelizedDomain;
    }
    virtual void apply();
    virtual void insert(plint processorLevel = 1);
    virtual void apply(std::vector<MultiBlock3D *> const &completionArg);
    virtual void insert(std::vector<MultiBlock3D *> const &completionArg, plint processorLevel = 1);
    virtual Array<T, 3> getForceOnObject();
    std::unique_ptr<MultiTensorField3D<T, 3> > computeVelocity(Box3D domain);
    std::unique_ptr<MultiTensorField3D<T, 3> > computeVelocity();
    std::unique_ptr<MultiTensorField3D<T, 3> > computeVorticity(Box3D domain);
    std::unique_ptr<MultiTensorField3D<T, 3> > computeVorticity();
    std::unique_ptr<MultiScalarField3D<T> > computeVelocityNorm(Box3D domain);
    std::unique_ptr<MultiScalarField3D<T> > computeVelocityNorm();
    std::unique_ptr<MultiScalarField3D<T> > computeVorticityNorm(Box3D domain);
    std::unique_ptr<MultiScalarField3D<T> > computeVorticityNorm();
    std::unique_ptr<MultiScalarField3D<T> > computeVelocityComponent(Box3D domain, plint iComp);
    std::unique_ptr<MultiScalarField3D<T> > computeVelocityComponent(plint iComp);
    std::unique_ptr<MultiScalarField3D<T> > computePressure(Box3D domain);
    std::unique_ptr<MultiScalarField3D<T> > computePressure();
    std::unique_ptr<MultiScalarField3D<T> > computeDensity(Box3D domain, T solidDensity = T());
    std::unique_ptr<MultiScalarField3D<T> > computeDensity(T solidDensity = T());
    std::unique_ptr<MultiScalarField3D<T> > computeStrainRateNorm();
    std::unique_ptr<MultiScalarField3D<T> > computeStrainRateNorm(Box3D domain);
    std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computeStrainRate();
    std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computeStrainRate(
        Box3D domain);
    std::unique_ptr<MultiScalarField3D<T> > computeShearStressNorm();
    std::unique_ptr<MultiScalarField3D<T> > computeShearStressNorm(Box3D domain);
    std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computeShearStress();
    std::unique_ptr<MultiTensorField3D<T, SymmetricTensor<T, Descriptor>::n> > computeShearStress(
        Box3D domain);
    T computeAverageVelocityComponent(Box3D domain, plint iComponent);
    Array<T, 3> computeAverageVelocity(Box3D domain);
    T computeAverageDensity(Box3D domain);
    T computeAverageDensity();
    T computeAverageEnergy(Box3D domain);
    T computeAverageEnergy();
    T computeRMSvorticity(Box3D domain);
    T computeRMSvorticity();
    T computeAverageShearStressNorm(Box3D domain);
    T computeAverageShearStressNorm();
    T computeRMSshearStressNorm(Box3D domain);
    T computeRMSshearStressNorm();

private:
    VoxelizedDomain3D<T> &voxelizedDomain;
    MultiBlockLattice3D<T, Descriptor> &lattice;
    MultiBlock3D &boundaryShapeArg;
    OffLatticeModel3D<T, BoundaryType> *offLatticeModel;
    MultiContainerBlock3D offLatticePattern;
};

}  // namespace plb

#endif  // OFF_LATTICE_BOUNDARY_CONDITION_3D_H
