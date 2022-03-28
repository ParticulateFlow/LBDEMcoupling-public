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

#ifndef OFF_LATTICE_BOUNDARY_PROCESSOR_3D_HH
#define OFF_LATTICE_BOUNDARY_PROCESSOR_3D_HH

#include <algorithm>
#include <cmath>

#include "offLattice/nextNeighbors3D.h"
#include "offLattice/offLatticeBoundaryProcessor3D.h"

namespace plb {

template <typename T>
CheckVoxelizationFunctional3D<T>::CheckVoxelizationFunctional3D(
    BoundaryShape3D<T, Array<T, 3> > *shape_, int flowType_) :
    numErrorsId(this->getStatistics().subscribeIntSum()), shape(shape_), flowType(flowType_)
{ }

template <typename T>
CheckVoxelizationFunctional3D<T>::~CheckVoxelizationFunctional3D()
{
    delete shape;
}

template <typename T>
CheckVoxelizationFunctional3D<T>::CheckVoxelizationFunctional3D(
    CheckVoxelizationFunctional3D<T> const &rhs) :
    PlainReductiveBoxProcessingFunctional3D(rhs),
    numErrorsId(rhs.numErrorsId),
    shape(rhs.shape->clone()),
    flowType(rhs.flowType)
{ }

template <typename T>
CheckVoxelizationFunctional3D<T> &CheckVoxelizationFunctional3D<T>::operator=(
    CheckVoxelizationFunctional3D<T> const &rhs)
{
    PlainReductiveBoxProcessingFunctional3D::operator=(rhs);
    numErrorsId = rhs.numErrorsId;
    shape = rhs.shape->clone();
    flowType = rhs.flowType;
    return *this;
}

template <typename T>
CheckVoxelizationFunctional3D<T> *CheckVoxelizationFunctional3D<T>::clone() const
{
    return new CheckVoxelizationFunctional3D<T>(*this);
}

template <typename T>
void CheckVoxelizationFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Flag matrix.
    // Possible additional parameters for the shape function are read-only.
    for (pluint i = 1; i < modified.size(); ++i) {
        modified[i] = modif::nothing;
    }
}

template <typename T>
BlockDomain::DomainT CheckVoxelizationFunctional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T>
bool CheckVoxelizationFunctional3D<T>::isFluid(Dot3D const &location) const
{
    if (flowType == voxelFlag::inside) {
        return shape->isInside(location);
    } else {
        return !shape->isInside(location);
    }
}

template <typename T>
void CheckVoxelizationFunctional3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_PRECONDITION(fields.size() >= 1);
    ScalarField3D<int> *flagMatrix = dynamic_cast<ScalarField3D<int> *>(fields[0]);
    PLB_ASSERT(flagMatrix);
    Dot3D absoluteOffset = flagMatrix->getLocation();

    BoundaryShape3D<T, Array<T, 3> > *newShape = 0;
    if (fields.size() > 1) {
        std::vector<AtomicBlock3D *> shapeParameters(fields.size() - 1);
        for (pluint i = 0; i < shapeParameters.size(); ++i) {
            shapeParameters[i] = fields[i + 1];
        }
        newShape = shape->clone(shapeParameters);
        std::swap(shape, newShape);
    }

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                computeCell(Dot3D(iX, iY, iZ), *flagMatrix, absoluteOffset);
            }
        }
    }

    if (newShape) {
        std::swap(shape, newShape);
        delete newShape;
    }
}

template <typename T>
void CheckVoxelizationFunctional3D<T>::computeCell(
    Dot3D const &cellLocation, ScalarField3D<int> &flagMatrix, Dot3D const &offset)
{
    //  Non-Fluid nodes.
    if (!isFluid(cellLocation + offset)) {
        plint numNeighbors = 0;
        plint numShallow = 0;
        plint numFailures = 0;
        for (int iNeighbor = 0; iNeighbor < NextNeighbor<T>::numNeighbors; ++iNeighbor) {
            int const *c = NextNeighbor<T>::c[iNeighbor];
            Dot3D neighbor(cellLocation.x + c[0], cellLocation.y + c[1], cellLocation.z + c[2]);
            Dot3D nextNeighbor(
                cellLocation.x + 2 * c[0], cellLocation.y + 2 * c[1], cellLocation.z + 2 * c[2]);
            // If the non-fluid node has a fluid neighbor ...
            if (isFluid(neighbor + offset)) {
                ++numNeighbors;
                Array<T, 3> wallNode;
                T wallDistance;
                Array<T, 3> wall_u, wallNormal;
                OffBoundary::Type bdType;
                plint id = -1;  // No optimization.
                bool ok = shape->pointOnSurface(
                    cellLocation + offset, Dot3D(c[0], c[1], c[2]), wallNode, wallDistance,
                    wallNormal, wall_u, bdType, id);
                if (!ok) {
                    ++numFailures;
                }
                if (!isFluid(nextNeighbor + offset)) {
                    ++numShallow;
                }
            }
        }
        if (numFailures > 0) {
            this->getStatistics().gatherIntSum(numErrorsId, 1);
            flagMatrix.get(cellLocation.x, cellLocation.y, cellLocation.z) = 1;
        }
        if (numNeighbors > 0 && numNeighbors == numShallow) {
            flagMatrix.get(cellLocation.x, cellLocation.y, cellLocation.z) = 2;
        }
    }
}

template <typename T>
plint CheckVoxelizationFunctional3D<T>::getNumErrors() const
{
    return this->getStatistics().getIntSum(numErrorsId);
}

}  // namespace plb

#endif  // OFF_LATTICE_BOUNDARY_PROCESSOR_3D_H
