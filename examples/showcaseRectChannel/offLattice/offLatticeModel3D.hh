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

#ifndef OFF_LATTICE_MODEL_3D_HH
#define OFF_LATTICE_MODEL_3D_HH

#include <algorithm>
#include <cmath>

#include "latticeBoltzmann/externalFieldAccess.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "offLattice/offLatticeModel3D.h"
#include "offLattice/voxelizer.h"

namespace plb {

template <typename T, class SurfaceData>
OffLatticeModel3D<T, SurfaceData>::OffLatticeModel3D(
    BoundaryShape3D<T, SurfaceData> *shape_, int flowType_) :
    shape(shape_),
    flowType(flowType_),
    velIsJflag(false),
    partialReplaceFlag(false),
    secondOrderFlag(true),
    regularizedModel(true),
    computeStat(true),
    defineVelocity(true)
{
    PLB_ASSERT(flowType == voxelFlag::inside || flowType == voxelFlag::outside);
}

template <typename T, class SurfaceData>
OffLatticeModel3D<T, SurfaceData>::OffLatticeModel3D(OffLatticeModel3D<T, SurfaceData> const &rhs) :
    shape(rhs.shape->clone()),
    flowType(rhs.flowType),
    velIsJflag(rhs.velIsJflag),
    partialReplaceFlag(rhs.partialReplaceFlag),
    secondOrderFlag(rhs.secondOrderFlag),
    regularizedModel(rhs.regularizedModel),
    computeStat(rhs.computeStat),
    defineVelocity(rhs.defineVelocity)
{ }

template <typename T, class SurfaceData>
OffLatticeModel3D<T, SurfaceData> &OffLatticeModel3D<T, SurfaceData>::operator=(
    OffLatticeModel3D<T, SurfaceData> const &rhs)
{
    delete shape;
    shape = rhs.shape->clone();
    flowType = rhs.flowType;
    velIsJflag = rhs.velIsJflag;
    partialReplaceFlag = rhs.partialReplaceFlag;
    secondOrderFlag = rhs.secondOrderFlag;
    regularizedModel = rhs.regularizedModel;
    computeStat = rhs.computeStat;
    defineVelocity = rhs.defineVelocity;
    return *this;
}

template <typename T, class SurfaceData>
OffLatticeModel3D<T, SurfaceData>::~OffLatticeModel3D()
{
    delete shape;
}

template <typename T, class SurfaceData>
void OffLatticeModel3D<T, SurfaceData>::provideShapeArguments(std::vector<AtomicBlock3D *> args)
{
    BoundaryShape3D<T, SurfaceData> *newShape = shape->clone(args);
    std::swap(shape, newShape);
    delete newShape;
}

template <typename T, class SurfaceData>
plint OffLatticeModel3D<T, SurfaceData>::getTag(plint id) const
{
    return shape->getTag(id);
}

template <typename T, class SurfaceData>
bool OffLatticeModel3D<T, SurfaceData>::pointOnSurface(
    Dot3D const &fromPoint, Dot3D const &direction, Array<T, 3> &locatedPoint, T &distance,
    Array<T, 3> &wallNormal, SurfaceData &surfaceData, OffBoundary::Type &bdType, plint &id) const
{
    return shape->gridPointOnSurface(
        fromPoint, direction, locatedPoint, distance, wallNormal, surfaceData, bdType, id);
}

template <typename T, class SurfaceData>
Array<T, 3> OffLatticeModel3D<T, SurfaceData>::computeContinuousNormal(
    Array<T, 3> const &p, plint id, bool isAreaWeighted) const
{
    return shape->computeContinuousNormal(p, id, isAreaWeighted);
}

template <typename T, class SurfaceData>
bool OffLatticeModel3D<T, SurfaceData>::intersectsSurface(
    Dot3D const &p1, Dot3D const &p2, plint &id) const
{
    return shape->intersectsSurface(p1, p2, id);
}

template <typename T, class SurfaceData>
bool OffLatticeModel3D<T, SurfaceData>::isFluid(Dot3D const &location) const
{
    if (flowType == voxelFlag::inside) {
        return shape->isInside(location);
    } else {
        return shape->isOutside(location);
    }
}

template <typename T, class SurfaceData>
bool OffLatticeModel3D<T, SurfaceData>::isSolid(Dot3D const &location) const
{
    if (flowType == voxelFlag::inside) {
        return shape->isOutside(location);
    } else {
        return shape->isInside(location);
    }
}

template <typename T, template <typename U> class Descriptor, class SurfaceData>
OffLatticeCompletionFunctional3D<T, Descriptor, SurfaceData>::OffLatticeCompletionFunctional3D(
    OffLatticeModel3D<T, SurfaceData> *offLatticeModel_, plint numShapeArgs_,
    plint numCompletionArgs_) :
    offLatticeModel(offLatticeModel_),
    numShapeArgs(numShapeArgs_),
    numCompletionArgs(numCompletionArgs_)
{ }

template <typename T, template <typename U> class Descriptor, class SurfaceData>
OffLatticeCompletionFunctional3D<T, Descriptor, SurfaceData>::~OffLatticeCompletionFunctional3D()
{
    delete offLatticeModel;
}

template <typename T, template <typename U> class Descriptor, class SurfaceData>
OffLatticeCompletionFunctional3D<T, Descriptor, SurfaceData>::OffLatticeCompletionFunctional3D(
    OffLatticeCompletionFunctional3D<T, Descriptor, SurfaceData> const &rhs) :
    offLatticeModel(rhs.offLatticeModel->clone()),
    numShapeArgs(rhs.numShapeArgs),
    numCompletionArgs(rhs.numCompletionArgs)
{ }

template <typename T, template <typename U> class Descriptor, class SurfaceData>
OffLatticeCompletionFunctional3D<T, Descriptor, SurfaceData>
    &OffLatticeCompletionFunctional3D<T, Descriptor, SurfaceData>::operator=(
        OffLatticeCompletionFunctional3D<T, Descriptor, SurfaceData> const &rhs)
{
    OffLatticeCompletionFunctional3D<T, Descriptor, SurfaceData>(rhs).swap(*this);
    return *this;
}

template <typename T, template <typename U> class Descriptor, class SurfaceData>
void OffLatticeCompletionFunctional3D<T, Descriptor, SurfaceData>::swap(
    OffLatticeCompletionFunctional3D<T, Descriptor, SurfaceData> &rhs)
{
    std::swap(offLatticeModel, rhs.offLatticeModel);
    std::swap(numShapeArgs, rhs.numShapeArgs);
    std::swap(numCompletionArgs, rhs.numCompletionArgs);
}

template <typename T, template <typename U> class Descriptor, class SurfaceData>
OffLatticeCompletionFunctional3D<T, Descriptor, SurfaceData>
    *OffLatticeCompletionFunctional3D<T, Descriptor, SurfaceData>::clone() const
{
    return new OffLatticeCompletionFunctional3D<T, Descriptor, SurfaceData>(*this);
}

template <typename T, template <typename U> class Descriptor, class SurfaceData>
void OffLatticeCompletionFunctional3D<T, Descriptor, SurfaceData>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    PLB_ASSERT((plint)modified.size() == 2 + numShapeArgs + numCompletionArgs);
    modified[0] = modif::staticVariables;  // Lattice.
    // It is very important that the "offLatticePattern" container block
    // which is passed as the second atomic-block with the off-lattice info
    // has the same multi-block management as the lattice used in the simulation.
    modified[1] = modif::nothing;  // Container for wet/dry nodes.
    // Possible additional parameters for the shape function and
    //   for the completion algorithm are read-only.
    for (pluint i = 2; i < modified.size(); ++i) {
        modified[i] = modif::nothing;
    }
}

template <typename T, template <typename U> class Descriptor, class SurfaceData>
BlockDomain::DomainT OffLatticeCompletionFunctional3D<T, Descriptor, SurfaceData>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor, class SurfaceData>
void OffLatticeCompletionFunctional3D<T, Descriptor, SurfaceData>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_PRECONDITION((plint)fields.size() == 2 + numShapeArgs + numCompletionArgs);
    AtomicBlock3D *lattice = fields[0];
    PLB_ASSERT(lattice);

    // It is very important that the "offLatticePattern" container block
    // which is passed as the second atomic-block with the off-lattice info
    // has the same multi-block management as the lattice used in the simulation.
    AtomicContainerBlock3D *container =  // Container for wet/dry nodes.
        dynamic_cast<AtomicContainerBlock3D *>(fields[1]);
    PLB_ASSERT(container);

    if (numShapeArgs > 0) {
        std::vector<AtomicBlock3D *> shapeParameters(numShapeArgs);
        for (plint i = 0; i < numShapeArgs; ++i) {
            shapeParameters[i] = fields[i + 2];
        }
        offLatticeModel->provideShapeArguments(shapeParameters);
    }
    std::vector<AtomicBlock3D *> completionParameters(numCompletionArgs);
    for (plint i = 0; i < numCompletionArgs; ++i) {
        completionParameters[i] = fields[i + 2 + numShapeArgs];
    }

    offLatticeModel->boundaryCompletion(*lattice, *container, completionParameters);
}

template <typename T, class SurfaceData>
OffLatticePatternFunctional3D<T, SurfaceData>::OffLatticePatternFunctional3D(
    OffLatticeModel3D<T, SurfaceData> *offLatticeModel_) :
    offLatticeModel(offLatticeModel_)
{ }

template <typename T, class SurfaceData>
OffLatticePatternFunctional3D<T, SurfaceData>::~OffLatticePatternFunctional3D()
{
    delete offLatticeModel;
}

template <typename T, class SurfaceData>
OffLatticePatternFunctional3D<T, SurfaceData>::OffLatticePatternFunctional3D(
    OffLatticePatternFunctional3D<T, SurfaceData> const &rhs) :
    offLatticeModel(rhs.offLatticeModel->clone())
{ }

template <typename T, class SurfaceData>
OffLatticePatternFunctional3D<T, SurfaceData>
    &OffLatticePatternFunctional3D<T, SurfaceData>::operator=(
        OffLatticePatternFunctional3D<T, SurfaceData> const &rhs)
{
    OffLatticePatternFunctional3D<T, SurfaceData>(rhs).swap(*this);
    return *this;
}

template <typename T, class SurfaceData>
void OffLatticePatternFunctional3D<T, SurfaceData>::swap(
    OffLatticePatternFunctional3D<T, SurfaceData> &rhs)
{
    std::swap(offLatticeModel, rhs.offLatticeModel);
}

template <typename T, class SurfaceData>
OffLatticePatternFunctional3D<T, SurfaceData>
    *OffLatticePatternFunctional3D<T, SurfaceData>::clone() const
{
    return new OffLatticePatternFunctional3D<T, SurfaceData>(*this);
}

template <typename T, class SurfaceData>
void OffLatticePatternFunctional3D<T, SurfaceData>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    // It is very important that the "offLatticePattern" container block
    // (the first atomic-block passed) has the same multi-block management
    // as the lattice used in the simulation.
    modified[0] = modif::staticVariables;  // Container.
    // Possible additional parameters for the shape function are read-only.
    for (pluint i = 1; i < modified.size(); ++i) {
        modified[i] = modif::nothing;
    }
}

template <typename T, class SurfaceData>
BlockDomain::DomainT OffLatticePatternFunctional3D<T, SurfaceData>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, class SurfaceData>
void OffLatticePatternFunctional3D<T, SurfaceData>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_PRECONDITION(fields.size() >= 1);
    AtomicContainerBlock3D *container = dynamic_cast<AtomicContainerBlock3D *>(fields[0]);
    PLB_ASSERT(container);
    ContainerBlockData *storeInfo = offLatticeModel->generateOffLatticeInfo();
    container->setData(storeInfo);

    if (fields.size() > 1) {
        std::vector<AtomicBlock3D *> shapeParameters(fields.size() - 1);
        for (pluint i = 0; i < shapeParameters.size(); ++i) {
            shapeParameters[i] = fields[i + 1];
        }
        offLatticeModel->provideShapeArguments(shapeParameters);
    }

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                offLatticeModel->prepareCell(Dot3D(iX, iY, iZ), *container);
            }
        }
    }
}

template <typename T, class SurfaceData>
GetForceOnObjectFunctional3D<T, SurfaceData>::GetForceOnObjectFunctional3D(
    OffLatticeModel3D<T, SurfaceData> *offLatticeModel_) :
    offLatticeModel(offLatticeModel_),
    forceId(
        this->getStatistics().subscribeSum(), this->getStatistics().subscribeSum(),
        this->getStatistics().subscribeSum())
{ }

template <typename T, class SurfaceData>
GetForceOnObjectFunctional3D<T, SurfaceData>::~GetForceOnObjectFunctional3D()
{
    delete offLatticeModel;
}

template <typename T, class SurfaceData>
GetForceOnObjectFunctional3D<T, SurfaceData>::GetForceOnObjectFunctional3D(
    GetForceOnObjectFunctional3D<T, SurfaceData> const &rhs) :
    PlainReductiveBoxProcessingFunctional3D(rhs),
    offLatticeModel(rhs.offLatticeModel->clone()),
    forceId(rhs.forceId)
{ }

template <typename T, class SurfaceData>
GetForceOnObjectFunctional3D<T, SurfaceData>
    &GetForceOnObjectFunctional3D<T, SurfaceData>::operator=(
        GetForceOnObjectFunctional3D<T, SurfaceData> const &rhs)
{
    delete offLatticeModel;
    offLatticeModel = rhs.offLatticeModel->clone();
    forceId = rhs.forceId;
    PlainReductiveBoxProcessingFunctional3D::operator=(rhs);
    return *this;
}

template <typename T, class SurfaceData>
void GetForceOnObjectFunctional3D<T, SurfaceData>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_PRECONDITION(fields.size() == 1);
    AtomicContainerBlock3D *offLatticeInfo = dynamic_cast<AtomicContainerBlock3D *>(fields[0]);
    PLB_ASSERT(offLatticeInfo);

    Array<T, 3> force = offLatticeModel->getLocalForce(*offLatticeInfo);
    this->getStatistics().gatherSum(forceId[0], force[0]);
    this->getStatistics().gatherSum(forceId[1], force[1]);
    this->getStatistics().gatherSum(forceId[2], force[2]);
}

template <typename T, class SurfaceData>
GetForceOnObjectFunctional3D<T, SurfaceData> *GetForceOnObjectFunctional3D<T, SurfaceData>::clone()
    const
{
    return new GetForceOnObjectFunctional3D<T, SurfaceData>(*this);
}

template <typename T, class SurfaceData>
void GetForceOnObjectFunctional3D<T, SurfaceData>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;  // Off-lattice info.
}

template <typename T, class SurfaceData>
BlockDomain::DomainT GetForceOnObjectFunctional3D<T, SurfaceData>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, class SurfaceData>
Array<T, 3> GetForceOnObjectFunctional3D<T, SurfaceData>::getForce() const
{
    return Array<T, 3>(
        this->getStatistics().getSum(forceId[0]), this->getStatistics().getSum(forceId[1]),
        this->getStatistics().getSum(forceId[2]));
}

template <typename T, class BoundaryType>
Array<T, 3> getForceOnObject(
    MultiBlock3D &offLatticePattern, OffLatticeModel3D<T, BoundaryType> const &offLatticeModel)
{
    std::vector<MultiBlock3D *> arg;
    arg.push_back(&offLatticePattern);
    GetForceOnObjectFunctional3D<T, BoundaryType> functional(offLatticeModel.clone());
    applyProcessingFunctional(functional, offLatticePattern.getBoundingBox(), arg);
    return functional.getForce();
}

}  // namespace plb

#endif  // OFF_LATTICE_MODEL_3D_HH
