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

#ifndef VOXELIZER_HH
#define VOXELIZER_HH

#include "atomicBlock/dataField3D.h"
#include "core/globalDefs.h"
#include "core/plbTimer.h"
#include "dataProcessors/dataAnalysisWrapper3D.h"
#include "dataProcessors/dataInitializerWrapper3D.h"
#include "dataProcessors/metaStuffWrapper3D.h"
#include "multiBlock/defaultMultiBlockPolicy3D.h"
#include "multiBlock/multiBlockGenerator3D.h"
#include "offLattice/voxelizer.h"

namespace plb {

namespace voxelFlag {
inline int invert(int arg)
{
    switch (arg) {
    case inside:
        return outside;
    case outside:
        return inside;
    case innerBorder:
        return outerBorder;
    case outerBorder:
        return innerBorder;
    case undetermined:
        return undetermined;
    default:
        PLB_ASSERT(false);
    }
    return undetermined;
}
inline int bulkFlag(int arg)
{
    if (arg == innerBorder || arg == inside) {
        return inside;
    } else if (arg == outerBorder || arg == outside) {
        return outside;
    } else {
        return undetermined;
    }
}
inline int borderFlag(int arg)
{
    if (arg == inside || arg == innerBorder) {
        return innerBorder;
    } else if (arg == outside || arg == outerBorder) {
        return outerBorder;
    } else {
        return undetermined;
    }
}
inline bool insideFlag(int arg)
{
    return arg == inside || arg == innerBorder;
}
inline bool outsideFlag(int arg)
{
    return arg == outside || arg == outerBorder;
}

}  // namespace voxelFlag

template <typename T>
std::unique_ptr<MultiScalarField3D<int> > voxelize(
    TriangularSurfaceMesh<T> const &mesh, plint symmetricLayer, plint borderWidth)
{
    Array<T, 2> xRange, yRange, zRange;
    mesh.computeBoundingBox(xRange, yRange, zRange);
    // Creation of the multi-scalar field. The +1 is because if the resolution is N,
    //   the number of nodes is N+1.
    plint nx = (plint)(xRange[1] - xRange[0]) + 1 + 2 * symmetricLayer;
    plint ny = (plint)(yRange[1] - yRange[0]) + 1 + 2 * symmetricLayer;
    plint nz = (plint)(zRange[1] - zRange[0]) + 1 + 2 * symmetricLayer;

    return voxelize(mesh, Box3D(0, nx - 1, 0, ny - 1, 0, nz - 1), borderWidth);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<int> > voxelize(
    TriangularSurfaceMesh<T> const &mesh, Box3D const &domain, plint borderWidth)
{
    // As initial seed, a one-cell layer around the outer boundary is tagged
    //   as ouside cells.
    plint envelopeWidth = 1;
    std::unique_ptr<MultiScalarField3D<int> > voxelMatrix =
        generateMultiScalarField<int>(domain, voxelFlag::outside, envelopeWidth);
    setToConstant(*voxelMatrix, voxelMatrix->getBoundingBox().enlarge(-1), voxelFlag::undetermined);

    MultiContainerBlock3D hashContainer(*voxelMatrix);
    std::vector<MultiBlock3D *> container_arg;
    container_arg.push_back(&hashContainer);
    applyProcessingFunctional(
        new CreateTriangleHash<T>(mesh), hashContainer.getBoundingBox(), container_arg);

    std::vector<MultiBlock3D *> flag_hash_arg;
    flag_hash_arg.push_back(voxelMatrix.get());
    flag_hash_arg.push_back(&hashContainer);

    voxelMatrix->resetFlags();  // Flags are used internally by VoxelizeMeshFunctional3D.
    plint maxIteration = 5000;
    plint i = 0;
    while (!allFlagsTrue(voxelMatrix.get()) && i < maxIteration) {
        applyProcessingFunctional(
            new VoxelizeMeshFunctional3D<T>(mesh), voxelMatrix->getBoundingBox(), flag_hash_arg);
        ++i;
    }
    if (i == maxIteration) {
        pcout << "Warning: Voxelization failed." << std::endl;
    }

    detectBorderLine(*voxelMatrix, voxelMatrix->getBoundingBox(), borderWidth);

    return voxelMatrix;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<int> > voxelize(
    TriangularSurfaceMesh<T> const &mesh, Box3D const &domain, plint borderWidth, Box3D seed)
{
    // As initial seed, a one-cell layer around the outer boundary is tagged
    //   as ouside cells.
    plint envelopeWidth = 1;

    std::unique_ptr<MultiScalarField3D<int> > voxelMatrix =
        generateMultiScalarField<int>(domain, voxelFlag::undetermined, envelopeWidth);
    setToConstant(*voxelMatrix, seed, voxelFlag::outside);

    MultiContainerBlock3D hashContainer(*voxelMatrix);
    std::vector<MultiBlock3D *> container_arg;
    container_arg.push_back(&hashContainer);
    applyProcessingFunctional(
        new CreateTriangleHash<T>(mesh), hashContainer.getBoundingBox(), container_arg);

    std::vector<MultiBlock3D *> flag_hash_arg;
    flag_hash_arg.push_back(voxelMatrix.get());
    flag_hash_arg.push_back(&hashContainer);

    voxelMatrix->resetFlags();  // Flags are used internally by VoxelizeMeshFunctional3D.
    plint maxIteration = 5000;
    plint i = 0;
    while (!allFlagsTrue(voxelMatrix.get()) && i < maxIteration) {
        bool useFullVoxelizationRange = (i == 0 ? true : false);
        applyProcessingFunctional(
            new VoxelizeMeshFunctional3D<T>(mesh, useFullVoxelizationRange),
            voxelMatrix->getBoundingBox(), flag_hash_arg);
        ++i;
    }
    if (i == maxIteration) {
        pcout << "Warning: Voxelization failed." << std::endl;
    }

    detectBorderLine(*voxelMatrix, voxelMatrix->getBoundingBox(), borderWidth);

    return voxelMatrix;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<int> > voxelize(
    TriangularSurfaceMesh<T> const &mesh, MultiBlockManagement3D const &management,
    plint borderWidth, Box3D seed)
{
    // As initial seed, a one-cell layer around the outer boundary is tagged
    //   as ouside cells.
    std::unique_ptr<MultiScalarField3D<int> > voxelMatrix =
        defaultGenerateMultiScalarField3D<int>(management, voxelFlag::undetermined);
    setToConstant(*voxelMatrix, seed, voxelFlag::outside);

    MultiContainerBlock3D hashContainer(*voxelMatrix);
    std::vector<MultiBlock3D *> container_arg;
    container_arg.push_back(&hashContainer);
    applyProcessingFunctional(
        new CreateTriangleHash<T>(mesh), hashContainer.getBoundingBox(), container_arg);

    std::vector<MultiBlock3D *> flag_hash_arg;
    flag_hash_arg.push_back(voxelMatrix.get());
    flag_hash_arg.push_back(&hashContainer);

    voxelMatrix->resetFlags();  // Flags are used internally by VoxelizeMeshFunctional3D.
    plint maxIteration = 5000;
    plint i = 0;
    while (!allFlagsTrue(voxelMatrix.get()) && i < maxIteration) {
        bool useFullVoxelizationRange = (i == 0 ? true : false);
        applyProcessingFunctional(
            new VoxelizeMeshFunctional3D<T>(mesh, useFullVoxelizationRange),
            voxelMatrix->getBoundingBox(), flag_hash_arg);
        ++i;
    }
    if (i == maxIteration) {
        pcout << "Warning: Voxelization failed." << std::endl;
    }

    detectBorderLine(*voxelMatrix, voxelMatrix->getBoundingBox(), borderWidth);

    return voxelMatrix;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<int> > revoxelize(
    TriangularSurfaceMesh<T> const &mesh, MultiScalarField3D<int> &oldVoxelMatrix,
    MultiContainerBlock3D &hashContainer, plint borderWidth)
{
    // As initial seed, a one-cell layer around the outer boundary is tagged
    //   as ouside cells.
    Box3D domain(oldVoxelMatrix.getBoundingBox());
    std::unique_ptr<MultiScalarField3D<int> > voxelMatrix(
        new MultiScalarField3D<int>((MultiBlock3D &)oldVoxelMatrix));
    setToConstant(*voxelMatrix, domain, voxelFlag::outside);
    setToConstant(*voxelMatrix, voxelMatrix->getBoundingBox().enlarge(-1), voxelFlag::undetermined);

    std::vector<MultiBlock3D *> flag_hash_arg;
    flag_hash_arg.push_back(voxelMatrix.get());
    flag_hash_arg.push_back(&hashContainer);

    voxelMatrix->resetFlags();  // Flags are used internally by VoxelizeMeshFunctional3D.
    plint maxIteration = 5000;
    plint i = 0;
    while (!allFlagsTrue(voxelMatrix.get()) && i < maxIteration) {
        applyProcessingFunctional(
            new VoxelizeMeshFunctional3D<T>(mesh), voxelMatrix->getBoundingBox(), flag_hash_arg);
        ++i;
    }
    if (i == maxIteration) {
        pcout << "Warning: Voxelization failed." << std::endl;
    }

    detectBorderLine(*voxelMatrix, voxelMatrix->getBoundingBox(), borderWidth);

    return voxelMatrix;
}

/* ******** VoxelizeMeshFunctional3D ************************************* */

template <typename T>
VoxelizeMeshFunctional3D<T>::VoxelizeMeshFunctional3D(
    TriangularSurfaceMesh<T> const &mesh_, bool useFullVoxelizationRange_) :
    mesh(mesh_), useFullVoxelizationRange(useFullVoxelizationRange_)
{ }

template <typename T>
bool VoxelizeMeshFunctional3D<T>::distanceToSurface(
    AtomicContainerBlock3D &hashContainer, Array<T, 3> const &point, T &distance,
    bool &isBehind) const
{
    T maxDistance = std::sqrt((T)3);
    Array<T, 2> xRange(point[0] - maxDistance, point[0] + maxDistance);
    Array<T, 2> yRange(point[1] - maxDistance, point[1] + maxDistance);
    Array<T, 2> zRange(point[2] - maxDistance, point[2] + maxDistance);
    TriangleHash<T> triangleHash(hashContainer);
    std::vector<plint> possibleTriangles;
    triangleHash.getTriangles(xRange, yRange, zRange, possibleTriangles);

    T tmpDistance;
    bool tmpIsBehind;
    bool triangleFound = false;

    for (pluint iPossible = 0; iPossible < possibleTriangles.size(); ++iPossible) {
        plint iTriangle = possibleTriangles[iPossible];
        mesh.distanceToTriangle(point, iTriangle, tmpDistance, tmpIsBehind);
        if (!triangleFound || tmpDistance < distance) {
            distance = tmpDistance;
            isBehind = tmpIsBehind;
            triangleFound = true;
        }
    }
    return triangleFound;
}

template <typename T>
bool VoxelizeMeshFunctional3D<T>::checkIfFacetsCrossed(
    AtomicContainerBlock3D &hashContainer, Array<T, 3> const &point1, Array<T, 3> const &point2,
    T &distance, plint &whichTriangle)
{
    Array<T, 2> xRange(std::min(point1[0], point2[0]), std::max(point1[0], point2[0]));
    Array<T, 2> yRange(std::min(point1[1], point2[1]), std::max(point1[1], point2[1]));
    Array<T, 2> zRange(std::min(point1[2], point2[2]), std::max(point1[2], point2[2]));
    TriangleHash<T> triangleHash(hashContainer);
    std::vector<plint> possibleTriangles;
    triangleHash.getTriangles(xRange, yRange, zRange, possibleTriangles);

    int flag = 0;              // Check for crossings inside the point1-point2 segment.
    Array<T, 3> intersection;  // Dummy variable.
    Array<T, 3> normal;        // Dummy variable.
    T tmpDistance;             // Dummy variable.

    if (global::counter("voxelizer-debug").getCount() == 1) {
        std::cout << "{";
    }
    std::vector<T> crossings;
    for (pluint iPossible = 0; iPossible < possibleTriangles.size(); ++iPossible) {
        plint iTriangle = possibleTriangles[iPossible];
        if (mesh.pointOnTriangle(point1, point2, flag, iTriangle, intersection, normal, tmpDistance)
            == 1) {
            if (global::counter("voxelizer-debug").getCount() == 1) {
                std::cout << "(" << iTriangle << ";" << tmpDistance << ")";
            }
            if (!util::fpequal_abs(tmpDistance, T(), mesh.eps1)) {
                crossings.push_back(tmpDistance);
            }
            if (crossings.size() == 1 || tmpDistance < distance) {
                distance = tmpDistance;
                whichTriangle = iTriangle;
            }
        }
    }
    if (global::counter("voxelizer-debug").getCount() == 1) {
        std::cout << "}";
    }

    if (crossings.size() == 0) {
        return false;
    } else {
        bool hasCrossed = true;
        for (pluint iCrossing = 1; iCrossing < crossings.size(); ++iCrossing) {
            // const T eps1 = std::numeric_limits<double>::epsilon()*1.e2;
            // if ( !util::fpequal(crossings[iCrossing], crossings[iCrossing-1], eps1) )

            // const T eps1 = std::numeric_limits<double>::epsilon()*1.e4;

            const T eps1 = std::numeric_limits<double>::epsilon() * 1.e4;
            if (std::fabs(crossings[iCrossing] - crossings[iCrossing - 1]) > eps1) {
                hasCrossed = !hasCrossed;
            }
        }
        return hasCrossed;
    }
}

template <typename T>
bool VoxelizeMeshFunctional3D<T>::createVoxelizationRange(
    Box3D const &domain, ScalarField3D<int> &voxels, Array<plint, 2> &xRange,
    Array<plint, 2> &yRange, Array<plint, 2> &zRange)
{
    // The purpose of the first three loops is to locate the eight
    //   corners of the cube. One voxel per corner would be insufficient
    //   because a potential seed is situated differently, depending on
    //   whether it is on the boundary of the multi-block or somewhere inside.
    for (plint dx = 0; dx <= +1; ++dx) {
        plint xMin = domain.x0 + dx * domain.getNx() - 1;
        plint xMax = domain.x0 + dx * domain.getNx();
        for (plint dy = 0; dy <= +1; ++dy) {
            plint yMin = domain.y0 + dy * domain.getNy() - 1;
            plint yMax = domain.y0 + dy * domain.getNy();
            for (plint dz = 0; dz <= +1; ++dz) {
                plint zMin = domain.z0 + dz * domain.getNz() - 1;
                plint zMax = domain.z0 + dz * domain.getNz();

                // Locate a potential seed in one of the corners.
                for (plint iX = xMin; iX <= xMax; ++iX) {
                    for (plint iY = yMin; iY <= yMax; ++iY) {
                        for (plint iZ = zMin; iZ <= zMax; ++iZ) {
                            if (voxels.get(iX, iY, iZ) != voxelFlag::undetermined) {
                                xRange[0] = domain.x0 + dx * (domain.getNx() - 1);
                                xRange[1] = domain.x0 + (1 - dx) * (domain.getNx() - 1);
                                yRange[0] = domain.y0 + dy * (domain.getNy() - 1);
                                yRange[1] = domain.y0 + (1 - dy) * (domain.getNy() - 1);
                                zRange[0] = domain.z0 + dz * (domain.getNz() - 1);
                                zRange[1] = domain.z0 + (1 - dz) * (domain.getNz() - 1);
                                return true;
                            }
                        }
                    }
                }
            }
        }
    }
    return false;
}

template <typename T>
void VoxelizeMeshFunctional3D<T>::printOffender(
    ScalarField3D<int> const &voxels, AtomicContainerBlock3D &hashContainer, Dot3D pos)
{
    std::set<plint> triangles;
    Dot3D offset = voxels.getLocation();
    Dot3D pos_ = pos + offset;
    std::cout << "Position (" << pos_.x << "," << pos_.y << "," << pos_.z << ")" << std::endl;
    for (plint dx = -1; dx <= +1; ++dx) {
        for (plint dy = -1; dy <= +1; ++dy) {
            for (plint dz = -1; dz <= +1; ++dz) {
                if (!(dx == 0 && dy == 0 && dz == 0)) {
                    Dot3D neigh = pos + offset + Dot3D(dx, dy, dz);
                    int typeOfNeighbor = voxels.get(pos.x + dx, pos.y + dy, pos.z + dz);
                    if (typeOfNeighbor != voxelFlag::undetermined) {
                        T distance;
                        plint whichTriangle;
                        Array<T, 3> p1(pos_.x, pos_.y, pos_.z);
                        Array<T, 3> p2(neigh.x, neigh.y, neigh.z);
                        global::counter("voxelizer-debug").increment(1);
                        bool crossed =
                            checkIfFacetsCrossed(hashContainer, p1, p2, distance, whichTriangle);
                        global::counter("voxelizer-debug").reset();
                        std::cout << "Neighbor (" << dx << "," << dy << "," << dz << "); is "
                                  << (voxelFlag::insideFlag(typeOfNeighbor) ? "inside" : "outside");
                        if (crossed) {
                            triangles.insert(whichTriangle);
                            std::cout << " inters. at distance " << distance << " with triangle "
                                      << whichTriangle << std::endl;
                        } else {
                            std::cout << " no inters." << std::endl;
                        }
                    }
                }
            }
        }
    }
    std::set<plint>::iterator it = triangles.begin();
    for (; it != triangles.end(); ++it) {
        std::cout << "Triangle " << *it << " [" << std::flush;
        Array<T, 3> p0 = mesh.getVertex(*it, 0);
        Array<T, 3> p1 = mesh.getVertex(*it, 1);
        Array<T, 3> p2 = mesh.getVertex(*it, 2);
        std::cout << p0[0] << " " << p1[0] << " " << p2[0] << " " << p0[0] << "], [" << p0[1] << " "
                  << p1[1] << " " << p2[1] << " " << p0[1] << "], [" << p0[2] << " " << p1[2] << " "
                  << p2[2] << " " << p0[2] << "]" << std::endl;
    }
}

template <typename T>
bool VoxelizeMeshFunctional3D<T>::voxelizeFromNeighbor(
    ScalarField3D<int> const &voxels, AtomicContainerBlock3D &hashContainer, Dot3D pos,
    Dot3D neighbor, int &voxelType)
{
    int verificationLevel = 0;
    Dot3D offset = voxels.getLocation();
    int typeOfNeighbor = voxels.get(neighbor.x, neighbor.y, neighbor.z);
    if (typeOfNeighbor == voxelFlag::undetermined) {
        return true;
    }
    // If there is no verification and the voxel has already been voxelized,
    //   it is not being re-voxelized here.
    if (verificationLevel == 0) {
        if (voxelType != voxelFlag::undetermined) {
            return true;
        }
    }
    Dot3D pos_ = pos + offset;
    Dot3D neighbor_ = neighbor + offset;
    Array<T, 3> point1((T)pos_.x, (T)pos_.y, (T)pos_.z);
    Array<T, 3> point2((T)neighbor_.x, (T)neighbor_.y, (T)neighbor_.z);
    int newVoxelType = voxelFlag::undetermined;
    T distance1, distance2, distance3, distance4;
    bool isBehind1, isBehind2;
    plint whichTriangle1, whichTriangle2;
    if (checkIfFacetsCrossed(hashContainer, point1, point2, distance1, whichTriangle1)) {
        // make sure that there aren't two triangles between the considered points!
        checkIfFacetsCrossed(hashContainer, point2, point1, distance2, whichTriangle2);

        if (whichTriangle1 == whichTriangle2)
            newVoxelType = voxelFlag::invert(typeOfNeighbor);

        // Additional consistency checks only at the ultimate level of verification.
        if (verificationLevel == 2) {
            PLB_ASSERT(distance1 < std::sqrt((T)3) + (T)0.0001);
#ifdef PLB_DEBUG
            bool ok =
                checkIfFacetsCrossed(hashContainer, point2, point1, distance2, whichTriangle2);
#else
            (void)checkIfFacetsCrossed(hashContainer, point2, point1, distance2, whichTriangle2);
#endif
            PLB_ASSERT(ok);
            PLB_ASSERT(distance2 < std::sqrt((T)3) + (T)0.0001);

#ifdef PLB_DEBUG
            bool ok1 = distanceToSurface(hashContainer, point1, distance3, isBehind1);
#else
            (void)distanceToSurface(hashContainer, point1, distance3, isBehind1);
#endif

            PLB_ASSERT(ok1);
            PLB_ASSERT(distance1 < std::sqrt((T)3) + (T)0.0001);
            // Attention: At this moment, the following consistency check fails sometimes,
            //   god knows why. It might be that there is a bug in the method
            //   mesh.distanceToSurface.
            PLB_ASSERT(
                (voxelFlag::insideFlag(newVoxelType) && isBehind1)
                || (voxelFlag::outsideFlag(newVoxelType) && !isBehind1));

#ifdef PLB_DEBUG
            bool ok2 = distanceToSurface(hashContainer, point2, distance4, isBehind2);
#else
            (void)distanceToSurface(hashContainer, point2, distance4, isBehind2);
#endif
            PLB_ASSERT(ok2);
            PLB_ASSERT(distance2 < std::sqrt((T)3) + (T)0.0001);
            PLB_ASSERT(
                (voxelFlag::insideFlag(typeOfNeighbor) && isBehind2)
                || (voxelFlag::outsideFlag(typeOfNeighbor) && !isBehind2));
        }
    } else {
        newVoxelType = typeOfNeighbor;
    }
    int oldVoxelType = voxelType;
    voxelType = newVoxelType;
    if (oldVoxelType == voxelFlag::undetermined) {
        return true;
    } else {
        return oldVoxelType == newVoxelType;
    }
}

template <typename T>
void VoxelizeMeshFunctional3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    ScalarField3D<int> *voxels = dynamic_cast<ScalarField3D<int> *>(blocks[0]);
    PLB_ASSERT(voxels);
    AtomicContainerBlock3D *container = dynamic_cast<AtomicContainerBlock3D *>(blocks[1]);
    PLB_ASSERT(container);

    // Return if this block is already voxelized.
    if (voxels->getFlag()) {
        return;
    }

    Array<plint, 2> xRange, yRange, zRange;
    if (useFullVoxelizationRange) {
        xRange[0] = domain.x0;
        xRange[1] = domain.x1;
        yRange[0] = domain.y0;
        yRange[1] = domain.y1;
        zRange[0] = domain.z0;
        zRange[1] = domain.z1;
    } else {
        if (!createVoxelizationRange(domain, *voxels, xRange, yRange, zRange)) {
            // If no seed has been found in the envelope, just return and wait
            //   for the next round.
            return;
        }
    }

    bool voxelizedEverything = true;

    // Specify if the loops go in positive or negative direction.
    plint xIncr = xRange[1] > xRange[0] ? 1 : -1;
    plint yIncr = yRange[1] > yRange[0] ? 1 : -1;
    plint zIncr = zRange[1] > zRange[0] ? 1 : -1;
    // The ranges are closed on both ends. Here, the range[1] end
    //   is converted to an open one so we can use != checks in the loops.
    xRange[1] += xIncr;
    yRange[1] += yIncr;
    zRange[1] += zIncr;
    for (plint iX = xRange[0]; iX != xRange[1]; iX += xIncr) {
        for (plint iY = yRange[0]; iY != yRange[1]; iY += yIncr) {
            for (plint iZ = zRange[0]; iZ != zRange[1]; iZ += zIncr) {
                Dot3D pos(iX, iY, iZ);
                int voxelType = voxels->get(iX, iY, iZ);
                if (voxelType == voxelFlag::undetermined) {
                    for (plint dx = -1; dx <= +1; ++dx) {
                        for (plint dy = -1; dy <= +1; ++dy) {
                            for (plint dz = -1; dz <= +1; ++dz) {
                                if (!(dx == 0 && dy == 0 && dz == 0)) {
                                    Dot3D neighbor(iX + dx, iY + dy, iZ + dz);
                                    bool ok = voxelizeFromNeighbor(
                                        *voxels, *container, pos, neighbor, voxelType);
                                    if (!ok) {
                                        printOffender(*voxels, *container, pos);
                                    }
                                    PLB_ASSERT(ok);
                                }
                            }
                        }
                    }
                    voxels->get(iX, iY, iZ) = voxelType;
                    if (voxelType == voxelFlag::undetermined) {
                        voxelizedEverything = false;
                    }
                }
            }
        }
    }
    // Indicate that this atomic-block has been voxelized.
    if (voxelizedEverything) {
        voxels->setFlag(true);
    }
}

template <typename T>
VoxelizeMeshFunctional3D<T> *VoxelizeMeshFunctional3D<T>::clone() const
{
    return new VoxelizeMeshFunctional3D<T>(*this);
}

template <typename T>
void VoxelizeMeshFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Voxels
    modified[1] = modif::nothing;          // Hash Container
}

template <typename T>
BlockDomain::DomainT VoxelizeMeshFunctional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** DetectBorderLineFunctional3D ************************************* */

template <typename T>
void detectBorderLine(MultiScalarField3D<T> &voxelMatrix, Box3D const &domain, plint borderWidth)
{
    applyProcessingFunctional(
        new DetectBorderLineFunctional3D<T>(borderWidth), domain, voxelMatrix);
}

template <typename T>
DetectBorderLineFunctional3D<T>::DetectBorderLineFunctional3D(plint borderWidth_) :
    borderWidth(borderWidth_)
{ }

template <typename T>
void DetectBorderLineFunctional3D<T>::process(Box3D domain, ScalarField3D<T> &voxels)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                for (plint dx = -borderWidth; dx <= borderWidth; ++dx)
                    for (plint dy = -borderWidth; dy <= borderWidth; ++dy)
                        for (plint dz = -borderWidth; dz <= borderWidth; ++dz)
                            if (!(dx == 0 && dy == 0 && dz == 0)) {
                                plint nextX = iX + dx;
                                plint nextY = iY + dy;
                                plint nextZ = iZ + dz;
                                if (contained(Dot3D(nextX, nextY, nextZ), voxels.getBoundingBox()))
                                {
                                    if (voxelFlag::outsideFlag(voxels.get(iX, iY, iZ))
                                        && voxelFlag::insideFlag(voxels.get(nextX, nextY, nextZ)))
                                    {
                                        voxels.get(iX, iY, iZ) = voxelFlag::outerBorder;
                                    }
                                    if (voxelFlag::insideFlag(voxels.get(iX, iY, iZ))
                                        && voxelFlag::outsideFlag(voxels.get(nextX, nextY, nextZ)))
                                    {
                                        voxels.get(iX, iY, iZ) = voxelFlag::innerBorder;
                                    }
                                }
                            }
            }
        }
    }
}

template <typename T>
DetectBorderLineFunctional3D<T> *DetectBorderLineFunctional3D<T>::clone() const
{
    return new DetectBorderLineFunctional3D<T>(*this);
}

template <typename T>
void DetectBorderLineFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT DetectBorderLineFunctional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

///* ******** ResetBorderLineFunctional3D ************************************* */

template <typename T>
void resetBorderFlags(MultiScalarField3D<T> &voxelMatrix, Box3D const &domain)
{
    applyProcessingFunctional(new ResetBorderLineFunctional3D<T>(), domain, voxelMatrix);
}

template <typename T>
ResetBorderLineFunctional3D<T>::ResetBorderLineFunctional3D()
{ }

template <typename T>
void ResetBorderLineFunctional3D<T>::process(Box3D domain, ScalarField3D<T> &voxels)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (voxels.get(iX, iY, iZ) == voxelFlag::innerBorder) {
                    voxels.get(iX, iY, iZ) = voxelFlag::inside;
                }
                if (voxels.get(iX, iY, iZ) == voxelFlag::outerBorder) {
                    voxels.get(iX, iY, iZ) = voxelFlag::outside;
                }
            }
        }
    }
}

template <typename T>
ResetBorderLineFunctional3D<T> *ResetBorderLineFunctional3D<T>::clone() const
{
    return new ResetBorderLineFunctional3D<T>(*this);
}

template <typename T>
void ResetBorderLineFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT ResetBorderLineFunctional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

}  // namespace plb

#endif  // VOXELIZER_HH
