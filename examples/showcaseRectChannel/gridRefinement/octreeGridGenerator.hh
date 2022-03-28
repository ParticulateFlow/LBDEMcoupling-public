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

#ifndef OCTREE_GRID_GENERATOR_HH
#define OCTREE_GRID_GENERATOR_HH

#include <cmath>
#include <cstdio>
#include <limits>

#include "gridRefinement/octree.h"
#include "gridRefinement/octree.hh"
#include "gridRefinement/octreeGridGenerator.h"
#include "io/imageWriter.h"
#include "io/parallelIO.h"
#include "offLattice/triangleSet.h"
#include "offLattice/triangleSet.hh"
#include "offLattice/triangleSetGenerator.h"
#include "offLattice/triangleSetGenerator.hh"

namespace plb {

template <typename T>
OctreeData<T>::OctreeData(Cuboid<T> const &cuboid_) :
    cuboid(cuboid_),
    touch(true),
    bulk(-1, -1, -1, -1, -1, -1),
    blockId(-10),  // Value used for debugging.
    isOverlap(false),
    allocate(true),
    processId(-1)
{ }

namespace octreeDataUtil {

template <typename T>
bool isAllocatableLeaf(OctreeNode<OctreeData<T> > *node)
{
    return (node->isLeaf && node->data->allocate);
}

template <typename T>
bool isOverlap(OctreeNode<OctreeData<T> > *node)
{
    return (node->data->isOverlap);
}

template <typename T>
bool isAllocatable(OctreeNode<OctreeData<T> > *node)
{
    return (isAllocatableLeaf(node) || isOverlap(node));
}

template <typename T>
Cuboid<T> getOctant(Cuboid<T> const &cuboid, int i)
{
    PLB_ASSERT(i >= 0 && i < 8);

    int z = i % 2;
    int y = (i / 2) % 2;
    int x = ((i / 2) / 2) % 2;

    T dx = (T)0.5 * (cuboid.x1() - cuboid.x0());
    T dy = (T)0.5 * (cuboid.y1() - cuboid.y0());
    T dz = (T)0.5 * (cuboid.z1() - cuboid.z0());

    Cuboid<T> octant;
    octant.x0() = cuboid.x0() + (T)x * dx;
    octant.y0() = cuboid.y0() + (T)y * dy;
    octant.z0() = cuboid.z0() + (T)z * dz;
    octant.x1() = octant.x0() + dx;
    octant.y1() = octant.y0() + dy;
    octant.z1() = octant.z0() + dz;

    return (octant);
}

template <typename T>
void getAllocatedMinMaxOctreeLeafNodeLevels(
    OctreeNode<OctreeData<T> > *root, int &minLevel, int &maxLevel)
{
    // minLevel and maxLevel must be properly initialized before calling this function.
    if (root == 0) {
        return;
    }

    if (isAllocatableLeaf(root)) {
        minLevel = std::min(minLevel, root->level);
        maxLevel = std::max(maxLevel, root->level);
    } else {
        for (int i = 0; i < 8; i++) {
            OctreeNode<OctreeData<T> > *child = root->child[i];
            if (child != 0) {
                getAllocatedMinMaxOctreeLeafNodeLevels(child, minLevel, maxLevel);
            }
        }
    }
}

}  // namespace octreeDataUtil

template <typename T>
void AssignOctreeData<T>::operator()(OctreeNode<OctreeData<T> > *node)
{
    using namespace octreeDataUtil;

    if (node->data == 0 && node->parent != 0) {
        Cuboid<T> const &parentCuboid = node->parent->data->cuboid;
        Cuboid<T> childCuboid(getOctant(parentCuboid, octreeChildType(node)));
        OctreeData<T> *childData = new OctreeData<T>(childCuboid);
        node->data = childData;
    }
}

template <typename T>
GridDensityFunction<T>::GridDensityFunction(std::string fileName, bool boundFromBelow)
{
    plbIOErrorIfCannotOpenFileForReading(fileName);
    std::ifstream file(fileName.c_str());

    file >> cuboid.x0();
    file >> cuboid.x1();
    file >> cuboid.y0();
    file >> cuboid.y1();
    file >> cuboid.z0();
    file >> cuboid.z1();

    file >> dx;

    file >> nx;
    file >> ny;
    file >> nz;

    gridDensity = new ScalarField3D<T>(nx, ny, nz);
    for (plint iX = 0; iX < nx; iX++) {
        for (plint iY = 0; iY < ny; iY++) {
            for (plint iZ = 0; iZ < nz; iZ++) {
                file >> gridDensity->get(iX, iY, iZ);
            }
        }
    }

    file.close();

    // Bound the grid density function.
    existNegativeValues = false;
    for (plint iX = 0; iX < nx; iX++) {
        for (plint iY = 0; iY < ny; iY++) {
            for (plint iZ = 0; iZ < nz; iZ++) {
                T &value = gridDensity->get(iX, iY, iZ);
                if (value > (T)1) {
                    value = (T)1;
                } else if (value < (T)0) {
                    existNegativeValues = true;
                    if (boundFromBelow) {
                        value = (T)0;
                    }
                }
            }
        }
    }
}

template <typename T>
GridDensityFunction<T>::GridDensityFunction(
    Cuboid<T> const &cuboid_, T dx_, ScalarField3D<T> const &gridDensity_, bool boundFromBelow) :
    cuboid(cuboid_), dx(dx_)
{
    gridDensity = new ScalarField3D<T>(gridDensity_);
    nx = gridDensity->getNx();
    ny = gridDensity->getNy();
    nz = gridDensity->getNz();

    // Bound the grid density function.
    existNegativeValues = false;
    for (plint iX = 0; iX < nx; iX++) {
        for (plint iY = 0; iY < ny; iY++) {
            for (plint iZ = 0; iZ < nz; iZ++) {
                T &value = gridDensity->get(iX, iY, iZ);
                if (value > (T)1) {
                    value = (T)1;
                } else if (value < (T)0) {
                    existNegativeValues = true;
                    if (boundFromBelow) {
                        value = (T)0;
                    }
                }
            }
        }
    }
}

template <typename T>
GridDensityFunction<T>::GridDensityFunction(GridDensityFunction<T> const &rhs) :
    cuboid(rhs.cuboid),
    dx(rhs.dx),
    nx(rhs.nx),
    ny(rhs.ny),
    nz(rhs.nz),
    existNegativeValues(rhs.existNegativeValues)
{
    gridDensity = new ScalarField3D<T>(*rhs.gridDensity);
}

template <typename T>
void GridDensityFunction<T>::swap(GridDensityFunction<T> &rhs)
{
    std::swap(cuboid, rhs.cuboid);
    std::swap(dx, rhs.dx);
    std::swap(nx, rhs.nx);
    std::swap(ny, rhs.ny);
    std::swap(nz, rhs.nz);
    std::swap(gridDensity, rhs.gridDensity);
    std::swap(existNegativeValues, rhs.existNegativeValues);
}

template <typename T>
GridDensityFunction<T> &GridDensityFunction<T>::operator=(GridDensityFunction<T> const &rhs)
{
    GridDensityFunction<T>(rhs).swap(*this);
    return (*this);
}

template <typename T>
GridDensityFunction<T> *GridDensityFunction<T>::clone() const
{
    return (new GridDensityFunction<T>(*this));
}

template <typename T>
GridDensityFunction<T>::~GridDensityFunction()
{
    delete gridDensity;
}

template <typename T>
bool GridDensityFunction<T>::hasNegativeValues() const
{
    return (existNegativeValues);
}

// The grid density function takes values in [0, 1], or negative values in regions where the user
// explicitly does not want blocks to be allocated (like inside a solid object for example).
template <typename T>
T GridDensityFunction<T>::operator()(Array<T, 3> const &pos) const
{
    T x = (pos[0] - cuboid.x0()) / dx;
    T y = (pos[1] - cuboid.y0()) / dx;
    T z = (pos[2] - cuboid.z0()) / dx;

    if (x <= (T)0 || x >= (T)(nx - 1) || y <= (T)0 || y >= (T)(ny - 1) || z <= (T)0
        || z >= (T)(nz - 1)) {
        return (0);
    }

    Array<T, 3> position(x, y, z);

    return (linearInterpolateScalarField(*gridDensity, position));
}

// The grid density function takes values in [0, 1], or negative values in regions where the user
// explicitly does not want blocks to be allocated (like inside a solid object for example).
template <typename T>
T GridDensityFunction<T>::operator()(Cuboid<T> const &c) const
{
    plint x0 = (plint)((c.x0() - cuboid.x0()) / dx);
    plint y0 = (plint)((c.y0() - cuboid.y0()) / dx);
    plint z0 = (plint)((c.z0() - cuboid.z0()) / dx);

    plint x1 = (plint)((c.x1() - cuboid.x0()) / dx) + 1;
    plint y1 = (plint)((c.y1() - cuboid.y0()) / dx) + 1;
    plint z1 = (plint)((c.z1() - cuboid.z0()) / dx) + 1;

    /*
    plint x0 = util::roundToInt((c.x0() - cuboid.x0()) / dx);
    plint y0 = util::roundToInt((c.y0() - cuboid.y0()) / dx);
    plint z0 = util::roundToInt((c.z0() - cuboid.z0()) / dx);

    plint x1 = util::roundToInt((c.x1() - cuboid.x0()) / dx);
    plint y1 = util::roundToInt((c.y1() - cuboid.y0()) / dx);
    plint z1 = util::roundToInt((c.z1() - cuboid.z0()) / dx);
    */

    Box3D domain(x0, x1, y0, y1, z0, z1);
    Box3D boundingBox(gridDensity->getBoundingBox());
    Box3D intersection;
    if (intersect(domain, boundingBox, intersection)) {
        T maxGridDensity = -std::numeric_limits<T>::max();
        for (plint iX = intersection.x0; iX <= intersection.x1; iX++) {
            for (plint iY = intersection.y0; iY <= intersection.y1; iY++) {
                for (plint iZ = intersection.z0; iZ <= intersection.z1; iZ++) {
                    maxGridDensity = std::max(maxGridDensity, gridDensity->get(iX, iY, iZ));
                }
            }
        }
        return (maxGridDensity);
    }

    return (0);
}

template <typename T>
RefineOctree<T>::RefineOctree(
    GridDensityFunction<T> const &gridDensityFunction_, T gridDensityScaleFactor_,
    int minLeafLevel_, int maxLeafLevel_, bool useSamples_, plint numSamples_) :
    gridDensityFunction(gridDensityFunction_),
    gridDensityScaleFactor(gridDensityScaleFactor_),
    minLeafLevel(minLeafLevel_),
    maxLeafLevel(maxLeafLevel_),
    useSamples(useSamples_),
    numSamples(numSamples_),
    refined(false)
{ }

template <typename T>
void RefineOctree<T>::reset()
{
    refined = false;
}

template <typename T>
bool RefineOctree<T>::refinedOctree() const
{
    return (refined);
}

template <typename T>
void RefineOctree<T>::operator()(OctreeNode<OctreeData<T> > *node)
{
    using namespace octreeDataUtil;

    if (isAllocatableLeaf(node) && node->data->touch) {
        if (node->level < minLeafLevel) {
            splitOctreeNode(node, assignOctreeData);
            refined = true;
        } else {
            Cuboid<T> const &cuboid = node->data->cuboid;
            T continuousLevel = -std::numeric_limits<T>::max();
            if (useSamples) {
                plint numSamplesAtLevel = numSamples;
                if (numSamplesAtLevel < 2) {
                    numSamplesAtLevel = util::intTwoToThePower(maxLeafLevel - node->level) + 1;
                    // numSamplesAtLevel = util::intTwoToThePower(maxLeafLevel + 1 - node->level) +
                    // 1;
                }
                continuousLevel = (maxLeafLevel - minLeafLevel)
                                      * getMaxScaledGridDensityFromPoints(cuboid, numSamplesAtLevel)
                                  + minLeafLevel;
            } else {
                continuousLevel =
                    (maxLeafLevel - minLeafLevel) * getMaxScaledGridDensityFromVolumes(cuboid)
                    + minLeafLevel;
            }

            // plint discreteLevel = (int) continuousLevel;
            // plint discreteLevel = (int) continuousLevel + 1;
            plint discreteLevel = util::roundToInt(continuousLevel);
            if (node->level < discreteLevel) {
                splitOctreeNode(node, assignOctreeData);
                refined = true;
            } else {
                node->data->touch = false;
            }
        }
    }
}

template <typename T>
T RefineOctree<T>::getMaxScaledGridDensityFromPoints(Cuboid<T> const &cuboid, plint n) const
{
    T a = (T)1 / (T)(n - 1);
    T dx = a * (cuboid.x1() - cuboid.x0());
    T dy = a * (cuboid.y1() - cuboid.y0());
    T dz = a * (cuboid.z1() - cuboid.z0());
    T maxVal = -std::numeric_limits<T>::max();
    for (plint i = 0; i < n; i++) {
        T x = cuboid.x0() + (T)i * dx;
        for (plint j = 0; j < n; j++) {
            T y = cuboid.y0() + (T)j * dy;
            for (plint k = 0; k < n; k++) {
                T z = cuboid.z0() + (T)k * dz;
                Array<T, 3> pos(x, y, z);
                maxVal = std::max(maxVal, gridDensityFunction(pos));
            }
        }
    }
    maxVal *= gridDensityScaleFactor;
    maxVal = maxVal > (T)1 ? (T)1 : maxVal;
    return (maxVal);
}

template <typename T>
T RefineOctree<T>::getMaxScaledGridDensityFromVolumes(Cuboid<T> const &cuboid) const
{
    T maxVal = gridDensityFunction(cuboid);
    maxVal *= gridDensityScaleFactor;
    maxVal = maxVal > (T)1 ? (T)1 : maxVal;
    return (maxVal);
}

template <typename T>
CalibrateOctree<T>::CalibrateOctree() : calibrated(false)
{ }

template <typename T>
void CalibrateOctree<T>::reset()
{
    calibrated = false;
}

template <typename T>
bool CalibrateOctree<T>::calibratedOctree() const
{
    return (calibrated);
}

template <typename T>
void CalibrateOctree<T>::operator()(OctreeNode<OctreeData<T> > *node)
{
    using namespace octreeDataUtil;

    if (isAllocatableLeaf(node)) {
        std::vector<OctreeNode<OctreeData<T> > *> neighbors =
            gatherOctreeGreaterEqualNeighbors(node);
        for (int i = 0; i < (int)neighbors.size(); i++) {
            if (neighbors[i] != 0 && isAllocatableLeaf(neighbors[i])
                && node->level > neighbors[i]->level + 1) {
                splitOctreeNode(neighbors[i], assignOctreeData);
                calibrated = true;
            }
        }
    }
}

template <typename T>
DeallocateOctreeLeafNodes<T>::DeallocateOctreeLeafNodes(
    GridDensityFunction<T> const &gridDensityFunction_, int maxLeafLevel_, bool useSamples_,
    plint numSamples_) :
    gridDensityFunction(gridDensityFunction_),
    maxLeafLevel(maxLeafLevel_),
    useSamples(useSamples_),
    numSamples(numSamples_)
{ }

template <typename T>
void DeallocateOctreeLeafNodes<T>::operator()(OctreeNode<OctreeData<T> > *node)
{
    using namespace octreeDataUtil;

    if (isAllocatableLeaf(node)) {
        Cuboid<T> const &cuboid = node->data->cuboid;
        bool isMaxGridDensityNegative = false;
        if (useSamples) {
            plint numSamplesAtLevel = numSamples;
            if (numSamplesAtLevel < 2) {
                numSamplesAtLevel = util::intTwoToThePower(maxLeafLevel - node->level) + 1;
                // numSamplesAtLevel = util::intTwoToThePower(maxLeafLevel + 1 - node->level) + 1;
            }
            isMaxGridDensityNegative = isMaxGridDensityNegativeOnPoints(cuboid, numSamplesAtLevel);
        } else {
            isMaxGridDensityNegative = isMaxGridDensityNegativeInVolumes(cuboid);
        }
        if (isMaxGridDensityNegative) {
            node->data->allocate = false;
        }
    }
}

template <typename T>
bool DeallocateOctreeLeafNodes<T>::isMaxGridDensityNegativeOnPoints(
    Cuboid<T> const &cuboid, plint n) const
{
    T a = (T)1 / (T)(n - 1);
    T dx = a * (cuboid.x1() - cuboid.x0());
    T dy = a * (cuboid.y1() - cuboid.y0());
    T dz = a * (cuboid.z1() - cuboid.z0());
    T maxVal = -std::numeric_limits<T>::max();
    for (plint i = 0; i < n; i++) {
        T x = cuboid.x0() + (T)i * dx;
        for (plint j = 0; j < n; j++) {
            T y = cuboid.y0() + (T)j * dy;
            for (plint k = 0; k < n; k++) {
                T z = cuboid.z0() + (T)k * dz;
                Array<T, 3> pos(x, y, z);
                maxVal = std::max(maxVal, gridDensityFunction(pos));
            }
        }
    }
    return (maxVal < (T)0);
}

template <typename T>
bool DeallocateOctreeLeafNodes<T>::isMaxGridDensityNegativeInVolumes(Cuboid<T> const &cuboid) const
{
    T maxVal = gridDensityFunction(cuboid);
    return (maxVal < (T)0);
}

template <typename T>
CountOctreeLeafNodes<T>::CountOctreeLeafNodes() : numLeaves(0)
{ }

template <typename T>
void CountOctreeLeafNodes<T>::reset()
{
    numLeaves = 0;
}

template <typename T>
plint CountOctreeLeafNodes<T>::getNumLeafNodes() const
{
    return (numLeaves);
}

template <typename T>
void CountOctreeLeafNodes<T>::operator()(OctreeNode<OctreeData<T> > *node)
{
    using namespace octreeDataUtil;

    if (isAllocatableLeaf(node)) {
        numLeaves++;
    }
}

template <typename T>
void GatherOctreeLeafCuboids<T>::reset()
{
    cuboids.clear();
}

template <typename T>
std::vector<Cuboid<T> > const &GatherOctreeLeafCuboids<T>::getCuboids() const
{
    return (cuboids);
}

template <typename T>
void GatherOctreeLeafCuboids<T>::operator()(OctreeNode<OctreeData<T> > *node)
{
    using namespace octreeDataUtil;

    if (isAllocatableLeaf(node)) {
        cuboids.push_back(node->data->cuboid);
    }
}

template <typename T>
void TouchAllOctreeNodes<T>::operator()(OctreeNode<OctreeData<T> > *node)
{
    node->data->touch = true;
}

template <typename T>
void AssignOctreeOverlapNodes<T>::operator()(OctreeNode<OctreeData<T> > *node)
{
    using namespace octreeDataUtil;

    if (isAllocatableLeaf(node)) {
        OctreeNode<OctreeData<T> > *parent = node->parent;
        if (parent != 0 && parent->data->touch) {
            parent->data->isOverlap = false;
            std::vector<OctreeNode<OctreeData<T> > *> neighbors =
                gatherOctreeEqualNeighbors(parent);
            for (int i = 0; i < (int)neighbors.size(); i++) {
                if (neighbors[i] != 0 && isAllocatableLeaf(neighbors[i])) {
                    parent->data->isOverlap = true;
                    break;
                }
            }
            parent->data->touch = false;
        }
    }
}

template <typename T>
void CheckOctreeOverlapNodes<T>::operator()(OctreeNode<OctreeData<T> > *node)
{
    typedef OctreeTables OT;
    using namespace octreeDataUtil;

    if (isOverlap(node)) {
        /*
        PLB_ASSERT(!node->isLeaf);
        int numCoupledChildren = 0;
        std::vector<OctreeNode<OctreeData<T> >*> neighbors = gatherOctreeEqualNeighbors(node);
        for (int i = 0; i < (int) neighbors.size(); i++) {
            if (neighbors[i] != 0 && isAllocatableLeaf(neighbors[i])) {
                for (int j = 0; j < 8; j++) {
                    if (OT::adj(i, j) && isAllocatable(node->child[j])) {
                        numCoupledChildren++;
                        PLB_ASSERT(!isOverlap(node->child[j]));
                    }
                }
            }
        }
        PLB_ASSERT(numCoupledChildren > 0);
        */

        PLB_ASSERT(!node->isLeaf);
        int numCoupledChildren = 0;
        std::vector<OctreeNode<OctreeData<T> > *> neighbors = gatherOctreeEqualNeighbors(node);
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < (int)neighbors.size(); j++) {
                if (OT::adj(j, i) && neighbors[j] != 0 && isAllocatableLeaf(neighbors[j])
                    && isAllocatable(node->child[i]))
                {
                    numCoupledChildren++;
                    PLB_ASSERT(!isOverlap(node->child[i]));
                }
            }
        }
        PLB_ASSERT(numCoupledChildren > 0);
    }
}

template <typename T>
ComputeOctreeBulks<T>::ComputeOctreeBulks(
    T dxFinest_, int maxLeafLevel_, Cuboid<T> const &fullDomain_, plint nx_, plint ny_, plint nz_) :
    dxFinest(dxFinest_),
    maxLeafLevel(maxLeafLevel_),
    fullDomain(fullDomain_),
    nx(nx_),
    ny(ny_),
    nz(nz_)
{ }

template <typename T>
void ComputeOctreeBulks<T>::operator()(OctreeNode<OctreeData<T> > *node)
{
    using namespace octreeDataUtil;

    // Bulks are computed only for the nodes to be allocated (leaves and overlaps).
    if (isAllocatable(node)) {
        Cuboid<T> const &cuboid = node->data->cuboid;
        Box3D &bulk = node->data->bulk;
        int level = node->level;

        T dx = dxFinest * (T)util::intTwoToThePower(maxLeafLevel - level);
        bulk.x0 = util::roundToInt((cuboid.x0() - fullDomain.x0()) / dx);
        bulk.x1 = bulk.x0 + nx - 1;
        bulk.y0 = util::roundToInt((cuboid.y0() - fullDomain.y0()) / dx);
        bulk.y1 = bulk.y0 + ny - 1;
        bulk.z0 = util::roundToInt((cuboid.z0() - fullDomain.z0()) / dx);
        bulk.z1 = bulk.z0 + nz - 1;
    }
}

template <typename T>
AssignOctreeLeafAndOverlapIds<T>::AssignOctreeLeafAndOverlapIds() : blockId(0)
{ }

template <typename T>
void AssignOctreeLeafAndOverlapIds<T>::operator()(OctreeNode<OctreeData<T> > *node)
{
    using namespace octreeDataUtil;

    if (isAllocatable(node)) {
        node->data->blockId = blockId;
        blockId++;
    }
}

template <typename T>
AssignOctreeOverlapProcessIds<T>::AssignOctreeOverlapProcessIds(
    int level_, int numLevelsToGroupOverlaps_, OctreeProcessLoads &octreeProcessLoads_) :
    level(level_),
    numLevelsToGroupOverlaps(numLevelsToGroupOverlaps_),
    octreeProcessLoads(octreeProcessLoads_),
    singleProcessId(-1)
{
    PLB_ASSERT(level >= 0);
    PLB_ASSERT(numLevelsToGroupOverlaps >= 0);

    groupingLevel = level - numLevelsToGroupOverlaps;
    if (groupingLevel < 0) {
        singleProcessId = octreeProcessLoads.getProcessWithMinLoad();
    }
}

template <typename T>
void AssignOctreeOverlapProcessIds<T>::operator()(OctreeNode<OctreeData<T> > *node)
{
    if (groupingLevel < 0) {
        if (needsProcessId(node)) {
            node->data->processId = singleProcessId;
            octreeProcessLoads.addLoad(singleProcessId, level);
        }
    } else {
        if (node->level == groupingLevel) {
            singleProcessId = -1;
            assignProcessIdsPreOrder(node);
        }
    }
}

template <typename T>
bool AssignOctreeOverlapProcessIds<T>::needsProcessId(OctreeNode<OctreeData<T> > *node)
{
    using namespace octreeDataUtil;
    return (node->level == level && isOverlap(node) && node->data->processId < 0);
}

template <typename T>
void AssignOctreeOverlapProcessIds<T>::assignProcessIdsPreOrder(OctreeNode<OctreeData<T> > *root)
{
    // Top-Down processing: Start processing at the root.
    if (needsProcessId(root)) {
        if (singleProcessId < 0) {
            singleProcessId = octreeProcessLoads.getProcessWithMinLoad();
        }
        root->data->processId = singleProcessId;
        octreeProcessLoads.addLoad(singleProcessId, level);
    }
    for (int i = 0; i < 8; i++) {
        OctreeNode<OctreeData<T> > *child = root->child[i];
        if (child != 0 && child->level <= level) {
            assignProcessIdsPreOrder(child);
        }
    }
}

template <typename T>
EnforceOctreeProcessAssignmentConstraints<T>::EnforceOctreeProcessAssignmentConstraints(
    OctreeProcessLoads &octreeProcessLoads_) :
    octreeProcessLoads(octreeProcessLoads_)
{ }

template <typename T>
void EnforceOctreeProcessAssignmentConstraints<T>::operator()(OctreeNode<OctreeData<T> > *node)
{
    typedef OctreeTables OT;
    using namespace octreeDataUtil;

    if (isOverlap(node)) {
        PLB_ASSERT(node->data->processId >= 0);
        plint overlapProcessId = node->data->processId;

        /*
        std::vector<OctreeNode<OctreeData<T> >*> neighbors = gatherOctreeEqualNeighbors(node);
        for (int i = 0; i < (int) neighbors.size(); i++) {
            if (neighbors[i] != 0 && isAllocatableLeaf(neighbors[i])) {
                for (int j = 0; j < 8; j++) {
                    if (OT::adj(i, j) && isAllocatable(node->child[j])) {
                        plint& childProcessId = node->child[j]->data->processId;
                        PLB_ASSERT(childProcessId < 0 || childProcessId == overlapProcessId);
                        if (childProcessId < 0) {
                            childProcessId = overlapProcessId;
                            octreeProcessLoads.addLoad(childProcessId, node->child[j]->level);
                        }
                    }
                }
            }
        }
        */

        std::vector<OctreeNode<OctreeData<T> > *> neighbors = gatherOctreeEqualNeighbors(node);
        for (int i = 0; i < 8; i++) {
            plint &childProcessId = node->child[i]->data->processId;
            if (childProcessId < 0) {
                for (int j = 0; j < (int)neighbors.size(); j++) {
                    if (OT::adj(j, i) && neighbors[j] != 0 && isAllocatableLeaf(neighbors[j])
                        && isAllocatable(node->child[i]))
                    {
                        childProcessId = overlapProcessId;
                        octreeProcessLoads.addLoad(childProcessId, node->child[i]->level);
                        break;
                    }
                }
            } else {
                PLB_ASSERT(isOverlap(node->child[i]));
            }
        }
    }
}

template <typename T>
AssignOctreeLeafProcessIds<T>::AssignOctreeLeafProcessIds(
    int level_, int numLevelsToGroupBlocks_, bool strongGrouping_,
    OctreeProcessLoads &octreeProcessLoads_) :
    level(level_),
    numLevelsToGroupBlocks(numLevelsToGroupBlocks_),
    strongGrouping(strongGrouping_),
    octreeProcessLoads(octreeProcessLoads_),
    singleProcessId(-1)
{
    PLB_ASSERT(level >= 0);
    PLB_ASSERT(numLevelsToGroupBlocks >= 0);

    overlapLevel = (level > 0 ? level - 1 : 0);

    groupingLevel = level - numLevelsToGroupBlocks;
    if (groupingLevel < 0) {
        singleProcessId = octreeProcessLoads.getProcessWithMinLoad();
    }
}

template <typename T>
void AssignOctreeLeafProcessIds<T>::operator()(OctreeNode<OctreeData<T> > *node)
{
    if (groupingLevel < 0) {
        if (needsProcessId(node)) {
            node->data->processId = singleProcessId;
            octreeProcessLoads.addLoad(singleProcessId, level);
        }
    } else {
        if (node->level == groupingLevel) {
            singleProcessId = -1;
            if (strongGrouping) {
                // The single overlap processId determination is only possible
                // when the overlaps are grouped at a level greater than or
                // equal to numLevelsToGroupBlocks-1. If they are grouped at a
                // level less than numLevelsToGroupBlocks-1 then we are not
                // be able to find a single overlap processId here.
                determineSingleOverlapProcessIdPreOrder(node);
            }
            assignProcessIdsPreOrder(node);
        }
    }
}

template <typename T>
bool AssignOctreeLeafProcessIds<T>::needsProcessId(OctreeNode<OctreeData<T> > *node)
{
    using namespace octreeDataUtil;
    return (node->level == level && isAllocatableLeaf(node) && node->data->processId < 0);
}

template <typename T>
void AssignOctreeLeafProcessIds<T>::determineSingleOverlapProcessIdPreOrder(
    OctreeNode<OctreeData<T> > *root)
{
    using namespace octreeDataUtil;

    // Top-Down processing: Start processing at the root.
    if (isOverlap(root) && root->level == overlapLevel) {
        if (singleProcessId < 0) {
            singleProcessId = root->data->processId;
        } else {
            PLB_ASSERT(singleProcessId == root->data->processId);
        }
    }
    for (int i = 0; i < 8; i++) {
        OctreeNode<OctreeData<T> > *child = root->child[i];
        if (child != 0 && child->level <= overlapLevel) {
            determineSingleOverlapProcessIdPreOrder(child);
        }
    }
}

template <typename T>
void AssignOctreeLeafProcessIds<T>::assignProcessIdsPreOrder(OctreeNode<OctreeData<T> > *root)
{
    // Top-Down processing: Start processing at the root.
    if (needsProcessId(root)) {
        if (singleProcessId < 0) {
            singleProcessId = octreeProcessLoads.getProcessWithMinLoad();
        }
        root->data->processId = singleProcessId;
        octreeProcessLoads.addLoad(singleProcessId, level);
    }
    for (int i = 0; i < 8; i++) {
        OctreeNode<OctreeData<T> > *child = root->child[i];
        if (child != 0 && child->level <= level) {
            assignProcessIdsPreOrder(child);
        }
    }
}

template <typename T>
GatherOctreeCuboidsInLatticeUnits<T>::GatherOctreeCuboidsInLatticeUnits(
    bool gatherOverlaps_, int maxLeafLevel_, int whichLevel_, plint whichProcess_) :
    gatherOverlaps(gatherOverlaps_),
    maxLeafLevel(maxLeafLevel_),
    whichLevel(whichLevel_),
    whichProcess(whichProcess_)
{
    PLB_ASSERT(whichLevel <= maxLeafLevel);
}

template <typename T>
std::vector<Cuboid<T> > const &GatherOctreeCuboidsInLatticeUnits<T>::getCuboids() const
{
    return (cuboids);
}

template <typename T>
void GatherOctreeCuboidsInLatticeUnits<T>::operator()(OctreeNode<OctreeData<T> > *node)
{
    using namespace octreeDataUtil;

    if ((!gatherOverlaps && isAllocatableLeaf(node)) || (gatherOverlaps && isOverlap(node))) {
        if ((whichLevel < 0 || node->level == whichLevel)
            && (whichProcess < 0 || node->data->processId == whichProcess))
        {
            plint stretch = util::intTwoToThePower(maxLeafLevel - node->level);
            Box3D bulk = node->data->bulk.multiply(stretch);

            Array<T, 3> llc((T)bulk.x0, (T)bulk.y0, (T)bulk.z0);
            Array<T, 3> urc((T)bulk.x1, (T)bulk.y1, (T)bulk.z1);
            cuboids.push_back(Cuboid<T>(llc, urc));
        }
    }
}

template <typename T>
ComputeOctreeGridStatistics<T>::ComputeOctreeGridStatistics(
    int minLeafLevel_, int maxLeafLevel_, plint numProcesses_, plint numCellsPerBlock_) :
    minLeafLevel(minLeafLevel_),
    maxLeafLevel(maxLeafLevel_),
    numProcesses(numProcesses_),
    numCellsPerBlock(numCellsPerBlock_),
    totalLoad(0)
{
    numLeafBlocksOnLevel.resize(maxLeafLevel + 1, 0);
    numOverlapBlocksOnLevel.resize(maxLeafLevel + 1, 0);

    numBlocksOnLevelOnProcess.resize(maxLeafLevel + 1);
    for (int level = 0; level <= maxLeafLevel; level++) {
        numBlocksOnLevelOnProcess[level].resize(numProcesses, 0);
    }

    loadOnLevelOnProcess.resize(maxLeafLevel + 1);
    for (int level = 0; level <= maxLeafLevel; level++) {
        loadOnLevelOnProcess[level].resize(numProcesses, 0);
    }

    totalLoadOnLevel.resize(maxLeafLevel + 1, 0);
}

template <typename T>
plint ComputeOctreeGridStatistics<T>::getNumLeafBlocksOnLevel(int level) const
{
    PLB_ASSERT(level >= 0 && level <= maxLeafLevel);
    return (numLeafBlocksOnLevel[level]);
}

template <typename T>
plint ComputeOctreeGridStatistics<T>::getNumOverlapBlocksOnLevel(int level) const
{
    PLB_ASSERT(level >= 0 && level <= maxLeafLevel);
    return (numOverlapBlocksOnLevel[level]);
}

template <typename T>
plint ComputeOctreeGridStatistics<T>::getNumBlocksOnLevelOnProcess(int level, plint processId) const
{
    PLB_ASSERT(level >= 0 && level <= maxLeafLevel);
    PLB_ASSERT(processId >= 0 && processId < numProcesses);
    return (numBlocksOnLevelOnProcess[level][processId]);
}

template <typename T>
plint ComputeOctreeGridStatistics<T>::getNumBlocksOnProcess(plint processId) const
{
    PLB_ASSERT(processId >= 0 && processId < numProcesses);
    plint numBlocksOnProcess = 0;
    for (int level = 0; level <= maxLeafLevel; level++) {
        numBlocksOnProcess += numBlocksOnLevelOnProcess[level][processId];
    }
    return (numBlocksOnProcess);
}

template <typename T>
T ComputeOctreeGridStatistics<T>::getRelativeLoadOnLevelOnProcess(int level, plint processId) const
{
    PLB_ASSERT(level >= 0 && level <= maxLeafLevel);
    PLB_ASSERT(processId >= 0 && processId < numProcesses);
    if (totalLoadOnLevel[level] != 0) {
        return (
            (T)numProcesses * (T)loadOnLevelOnProcess[level][processId]
            / (T)totalLoadOnLevel[level]);
    } else {
        return (1);
    }
}

template <typename T>
T ComputeOctreeGridStatistics<T>::getRelativeLoadOnProcess(plint processId) const
{
    PLB_ASSERT(processId >= 0 && processId < numProcesses);
    plint totalLoadOnProcess = 0;
    for (int level = 0; level <= maxLeafLevel; level++) {
        totalLoadOnProcess += loadOnLevelOnProcess[level][processId];
    }
    return ((T)numProcesses * (T)totalLoadOnProcess / (T)totalLoad);
}

template <typename T>
T ComputeOctreeGridStatistics<T>::getStandardDeviationOfRelativeLoadOnLevel(int level) const
{
    T sigma = 0.0;
    for (plint iProcess = 0; iProcess < numProcesses; iProcess++) {
        sigma += util::sqr(getRelativeLoadOnLevelOnProcess(level, iProcess) - (T)1);
    }
    sigma = std::sqrt(sigma / (T)numProcesses);
    return (sigma);
}

template <typename T>
T ComputeOctreeGridStatistics<T>::getStandardDeviationOfRelativeLoad() const
{
    T sigma = 0.0;
    for (plint iProcess = 0; iProcess < numProcesses; iProcess++) {
        sigma += util::sqr(getRelativeLoadOnProcess(iProcess) - (T)1);
    }
    sigma = std::sqrt(sigma / (T)numProcesses);
    return (sigma);
}

template <typename T>
T ComputeOctreeGridStatistics<T>::getRangeOfRelativeLoadOnLevel(int level) const
{
    PLB_ASSERT(level >= 0 && level <= maxLeafLevel);
    if (totalLoadOnLevel[level] == 0) {
        return (0);
    }

    plint minLoad = std::numeric_limits<plint>::max();
    plint maxLoad = 0;
    for (plint iProcess = 0; iProcess < numProcesses; iProcess++) {
        plint load = loadOnLevelOnProcess[level][iProcess];
        if (load < minLoad) {
            minLoad = load;
        }
        if (load > maxLoad) {
            maxLoad = load;
        }
    }
    T range = (T)numProcesses * (T)(maxLoad - minLoad) / (T)totalLoadOnLevel[level];
    return (range);
}

template <typename T>
T ComputeOctreeGridStatistics<T>::getRangeOfRelativeLoad() const
{
    plint minLoad = std::numeric_limits<plint>::max();
    plint maxLoad = 0;
    for (plint iProcess = 0; iProcess < numProcesses; iProcess++) {
        plint load = 0;
        for (int level = 0; level <= maxLeafLevel; level++) {
            load += loadOnLevelOnProcess[level][iProcess];
        }
        if (load < minLoad) {
            minLoad = load;
        }
        if (load > maxLoad) {
            maxLoad = load;
        }
    }
    T range = (T)numProcesses * (T)(maxLoad - minLoad) / (T)totalLoad;
    return (range);
}

template <typename T>
void ComputeOctreeGridStatistics<T>::operator()(OctreeNode<OctreeData<T> > *node)
{
    using namespace octreeDataUtil;

    if (isAllocatableLeaf(node)) {
        numLeafBlocksOnLevel[node->level]++;
        numBlocksOnLevelOnProcess[node->level][node->data->processId]++;

        plint blockLoad = numCellsPerBlock * util::intTwoToThePower(node->level - minLeafLevel);
        loadOnLevelOnProcess[node->level][node->data->processId] += blockLoad;

        totalLoadOnLevel[node->level] += blockLoad;
        totalLoad += blockLoad;
    } else if (isOverlap(node)) {
        numOverlapBlocksOnLevel[node->level]++;
        numBlocksOnLevelOnProcess[node->level][node->data->processId]++;

        plint blockLoad = numCellsPerBlock * util::intTwoToThePower(node->level - minLeafLevel);
        loadOnLevelOnProcess[node->level][node->data->processId] += blockLoad;

        totalLoadOnLevel[node->level] += blockLoad;
        totalLoad += blockLoad;
    }
}

template <typename T>
ConstructOctreeGridStructure<T>::ConstructOctreeGridStructure(
    OctreeGridStructure &octreeGridStructure_, int minLeafLevel_) :
    octreeGridStructure(octreeGridStructure_), minLeafLevel(minLeafLevel_)
{ }

template <typename T>
void ConstructOctreeGridStructure<T>::operator()(OctreeNode<OctreeData<T> > *node)
{
    typedef OctreeTables OT;
    using namespace octreeDataUtil;

    if (isAllocatable(node)) {
        bool isOverlap = node->data->isOverlap;

        plint blockId = node->data->blockId;
        Box3D const &bulk = node->data->bulk;
        plint level = node->level - minLeafLevel;
        plint processId = node->data->processId;
        octreeGridStructure.addBlock(blockId, bulk, level, processId, isOverlap);

        // Set the connectivity only for the non-overlapping nodes.
        if (!isOverlap) {
            std::vector<OctreeNode<OctreeData<T> > *> neighbors =
                gatherOctreeGreaterEqualNeighbors(node);
            PLB_ASSERT(neighbors.size() == 26);
            Array<plint, 26> neighborLeafIds;
            for (int i = 0; i < 26; i++) {
                if (neighbors[i] == 0) {
                    neighborLeafIds[i] = OT::border();  // The node is touching the boundary.
                } else {
                    if (neighbors[i]->isLeaf) {
                        if (neighbors[i]->data->allocate) {
                            neighborLeafIds[i] = neighbors[i]->data->blockId;
                        } else {
                            neighborLeafIds[i] = OT::undef();
                        }
                    } else {
                        neighborLeafIds[i] =
                            OT::smaller();  // The node is touching a finer level of refinement.
                    }
                }
            }
            octreeGridStructure.addNeighborBlockIds(blockId, neighborLeafIds);
        }
    }
}

template <typename T>
OctreeGridGenerator<T>::OctreeGridGenerator(std::string xmlFileName) : outS(&pcout), errS(&pcerr)
{
    abortIfCannotOpenFileForReading(xmlFileName);
    XMLreader document(xmlFileName);

    Array<T, 2> x, y, z;
    document["simulationDomain"]["x"].read<T, 2>(x);
    plbIOError(util::lessEqual(x[1], x[0]), "The x-extent of the simulation domain is wrong.");
    document["simulationDomain"]["y"].read<T, 2>(y);
    plbIOError(util::lessEqual(y[1], y[0]), "The y-extent of the simulation domain is wrong.");
    document["simulationDomain"]["z"].read<T, 2>(z);
    plbIOError(util::lessEqual(z[1], z[0]), "The z-extent of the simulation domain is wrong.");

    fullDomain = Cuboid<T>(x[0], x[1], y[0], y[1], z[0], z[1]);

    document["gridDensityFunctionFile"].read(gridDensityFunctionFile);

    xPeriodic = false;
    try {
        document["xPeriodic"].read(xPeriodic);
    } catch (PlbIOException &exception) {
    }

    yPeriodic = false;
    try {
        document["yPeriodic"].read(yPeriodic);
    } catch (PlbIOException &exception) {
    }

    zPeriodic = false;
    try {
        document["zPeriodic"].read(zPeriodic);
    } catch (PlbIOException &exception) {
    }

    gridDensityScaleFactor = (T)1;
    try {
        document["gridDensityScaleFactor"].read(gridDensityScaleFactor);
        plbIOError(
            util::lessEqual(gridDensityScaleFactor, (T)0),
            "The grid density scale factor must be greater than 0.");
    } catch (PlbIOException &exception) {
    }

    document["minLevel"].read(minLeafLevel);
    plbIOError(minLeafLevel < 0, "The minimum grid level must be non-negative.");
    document["maxLevel"].read(maxLeafLevel);
    plbIOError(
        maxLeafLevel < minLeafLevel,
        "The maximum grid level must be greater than or equal to the minimum grid level.");

    useSamples = false;
    try {
        document["useSamples"].read(useSamples);
    } catch (PlbIOException &exception) {
    }

    numSamples = -1;
    if (useSamples) {
        try {
            document["numSamples"].read(numSamples);
        } catch (PlbIOException &exception) {
        }
    }

    maxIter = 100;
    try {
        document["maxIter"].read(maxIter);
        plbIOError(maxIter < 0, "The maximum number of iterations must be non-negative.");
    } catch (PlbIOException &exception) {
    }

    removeBlocks = true;
    try {
        document["removeBlocks"].read(removeBlocks);
    } catch (PlbIOException &exception) {
    }

    document["n"].read(n);
    plbIOError(
        !(n % 2 == 0 && n >= 6),
        "The block size (n) must be an even positive integer greater than 5.");
    document["numProcesses"].read(numProcesses);
    plbIOError(numProcesses <= 0, "The number of processes must be greater than 0.");

    fineToCoarse = true;
    try {
        document["fineToCoarse"].read(fineToCoarse);
    } catch (PlbIOException &exception) {
    }

    numLevelsToGroupBlocks = 0;
    try {
        document["numLevelsToGroupBlocks"].read(numLevelsToGroupBlocks);
    } catch (PlbIOException &exception) {
    }
    plbIOError(
        numLevelsToGroupBlocks < 0, "The number of levels to group blocks must be non-negative.");

    // We have taken the decision that the overlaps are grouped by default
    // in a level one less that the grouping level of the leaf blocks.
    numLevelsToGroupOverlaps = (numLevelsToGroupBlocks > 0 ? numLevelsToGroupBlocks - 1 : 0);
    try {
        document["numLevelsToGroupOverlaps"].read(numLevelsToGroupOverlaps);
    } catch (PlbIOException &exception) {
    }
    plbIOError(
        numLevelsToGroupOverlaps < 0,
        "The number of levels to group overlaps must be non-negative.");

    strongGrouping = false;
    try {
        document["strongGrouping"].read(strongGrouping);
    } catch (PlbIOException &exception) {
    }

    // The strong grouping (single overlap processId determination) is only
    // possible when the overlaps are grouped at a level greater than or equal to
    // numLevelsToGroupBlocks-1. If they are grouped at a level less than
    // numLevelsToGroupBlocks-1 then we are able to achieve strong coupling
    // with its current meaning.
    if (numLevelsToGroupOverlaps < numLevelsToGroupBlocks - 1) {
        strongGrouping = false;
    }

    outDir = "./";
    try {
        document["output"]["directory"].read(outDir);
        if (outDir[outDir.size() - 1] != '/') {
            outDir += '/';
        }
        abortIfCannotCreateFileInDir(outDir, "plb-checkfile.txt");
    } catch (PlbIOException &exception) {
    }

    verbose = true;
    try {
        document["output"]["verbose"].read(verbose);
    } catch (PlbIOException &exception) {
    }

    stlOutput = true;
    try {
        document["output"]["stlOutput"].read(stlOutput);
    } catch (PlbIOException &exception) {
    }

    stlBaseName = "octree";
    if (stlOutput) {
        try {
            document["output"]["stlBaseName"].read(stlBaseName);
        } catch (PlbIOException &exception) {
        }
    }
}

template <typename T>
OctreeGridGenerator<T>::OctreeGridGenerator(
    Cuboid<T> const &fullDomain_, std::string gridDensityFunctionFile_, int minLeafLevel_,
    int maxLeafLevel_, plint n_, plint numProcesses_, bool xPeriodic_, bool yPeriodic_,
    bool zPeriodic_, T gridDensityScaleFactor_, bool useSamples_, plint numSamples_, plint maxIter_,
    bool removeBlocks_, bool fineToCoarse_, int numLevelsToGroupBlocks_,
    int numLevelsToGroupOverlaps_, bool strongGrouping_, std::string outDir_, bool verbose_,
    bool stlOutput_, std::string stlBaseName_) :
    fullDomain(fullDomain_),
    gridDensityFunctionFile(gridDensityFunctionFile_),
    xPeriodic(xPeriodic_),
    yPeriodic(yPeriodic_),
    zPeriodic(zPeriodic_),
    gridDensityScaleFactor(gridDensityScaleFactor_),
    minLeafLevel(minLeafLevel_),
    maxLeafLevel(maxLeafLevel_),
    useSamples(useSamples_),
    numSamples(numSamples_),
    maxIter(maxIter_),
    removeBlocks(removeBlocks_),
    n(n_),
    numProcesses(numProcesses_),
    fineToCoarse(fineToCoarse_),
    numLevelsToGroupBlocks(numLevelsToGroupBlocks_),
    numLevelsToGroupOverlaps(numLevelsToGroupOverlaps_),
    strongGrouping(strongGrouping_),
    outDir(outDir_),
    verbose(verbose_),
    stlOutput(stlOutput_),
    stlBaseName(stlBaseName_),
    outS(&pcout),
    errS(&pcerr)
{
    plbIOError(
        util::lessEqual(fullDomain.x1(), fullDomain.x0()),
        "The x-extent of the simulation domain is wrong.");
    plbIOError(
        util::lessEqual(fullDomain.y1(), fullDomain.y0()),
        "The y-extent of the simulation domain is wrong.");
    plbIOError(
        util::lessEqual(fullDomain.z1(), fullDomain.z0()),
        "The z-extent of the simulation domain is wrong.");
    plbIOError(
        util::lessEqual(gridDensityScaleFactor, (T)0),
        "The grid density scale factor must be greater than 0.");
    plbIOError(minLeafLevel < 0, "The minimum grid level must be non-negative.");
    plbIOError(
        maxLeafLevel < minLeafLevel,
        "The maximum grid level must be greater than or equal to the minimum grid level.");
    plbIOError(maxIter < 0, "The maximum number of iterations must be non-negative.");
    plbIOError(
        !(n % 2 == 0 && n >= 6),
        "The block size (n) must be an even positive integer greater than 5.");
    plbIOError(numProcesses <= 0, "The number of processes must be greater than 0.");
    plbIOError(
        numLevelsToGroupBlocks < 0, "The number of levels to group blocks must be non-negative.");

    if (numLevelsToGroupOverlaps < 0) {
        // We have taken the decision that the overlaps are grouped by default
        // in a level one less that the grouping level of the leaf blocks.
        numLevelsToGroupOverlaps = (numLevelsToGroupBlocks > 0 ? numLevelsToGroupBlocks - 1 : 0);
    }

    // The strong grouping (single overlap processId determination) is only
    // possible when the overlaps are grouped at a level greater than or equal to
    // numLevelsToGroupBlocks-1. If they are grouped at a level less than
    // numLevelsToGroupBlocks-1 then we are able to achieve strong coupling
    // with its current meaning.
    if (numLevelsToGroupOverlaps < numLevelsToGroupBlocks - 1) {
        strongGrouping = false;
    }
}

template <typename T>
OctreeGridStructure OctreeGridGenerator<T>::generateOctreeGridStructure()
{
    using namespace octreeDataUtil;

    if (verbose) {
        (*outS) << "=============================" << std::endl;
        (*outS) << "Palabos Octree Grid Generator" << std::endl;
        (*outS) << "=============================" << std::endl;
        (*outS) << std::endl;
    }

    // Octree initialization.

    if (verbose) {
        (*outS) << "Creating the octree data structure... " << std::flush;
    }
    OctreeData<T> *rootOctreeData = new OctreeData<T>(fullDomain);
    OctreeNode<OctreeData<T> > *root = new OctreeNode<OctreeData<T> >(0, rootOctreeData);
    if (verbose) {
        (*outS) << "done." << std::endl;
        (*outS) << std::endl;
    }

    // Octree refinement.

    if (verbose) {
        (*outS) << "Refining the octree." << std::endl;
    }
    bool boundFromBelow = true;  // In the beginning we act like grid density is positive, and all
                                 // blocks are allocated.
    GridDensityFunction<T> *boundedGridDensityFunction =
        new GridDensityFunction<T>(gridDensityFunctionFile, boundFromBelow);
    bool gridDensityHasNegativeValues = boundedGridDensityFunction->hasNegativeValues();
    RefineOctree<T> refineOctree(
        *boundedGridDensityFunction, gridDensityScaleFactor, minLeafLevel, maxLeafLevel, useSamples,
        numSamples);
    int iter = 0;
    do {
        if (verbose) {
            (*outS) << "    Octree refinement iteration: " << iter << std::endl;
        }
        refineOctree.reset();
        processOctreePostOrder(root, refineOctree);
        iter++;
    } while (refineOctree.refinedOctree() && iter < maxIter);

    if (refineOctree.refinedOctree()) {
        (*errS) << std::endl;
        (*errS) << "Octree refinement failed after: " << maxIter << " iterations." << std::endl;
        exit(1);
    }
    if (verbose) {
        (*outS) << std::endl;
    }

    delete boundedGridDensityFunction;
    boundedGridDensityFunction = 0;

    CountOctreeLeafNodes<T> countOctreeLeafNodes;
    if (verbose) {
        processOctreePostOrder(root, countOctreeLeafNodes);
        (*outS) << "Number of leaf nodes before calibration: "
                << countOctreeLeafNodes.getNumLeafNodes() << std::endl;
        (*outS) << std::endl;
    }

    GatherOctreeLeafCuboids<T> gatherOctreeLeafCuboids;
    TriangleSet<T> *triangleSet = 0;
    if (stlOutput) {
        if (verbose) {
            (*outS) << "Saving the octree before calibration... " << std::flush;
        }
        processOctreePostOrder(root, gatherOctreeLeafCuboids);
        TriangleSet<T> *triangleSet = new TriangleSet<T>(DBL);
        for (plint iCuboid = 0; iCuboid < (plint)gatherOctreeLeafCuboids.getCuboids().size();
             iCuboid++) {
            Array<plint, 3> nSegments(1, 1, 1);
            TriangleSet<T> cuboidSet =
                constructCuboid<T>(gatherOctreeLeafCuboids.getCuboids()[iCuboid], nSegments);
            triangleSet->append(cuboidSet);
        }
        triangleSet->writeBinarySTL(outDir + stlBaseName + "_before_calibration.stl");
        delete triangleSet;
        triangleSet = 0;
        if (verbose) {
            (*outS) << "done." << std::endl;
            (*outS) << std::endl;
        }
    }

    // Octree calibration (two adjacent tree levels must not differ by more than 1).

    if (verbose) {
        (*outS) << "Calibrating the octree." << std::endl;
    }
    CalibrateOctree<T> calibrateOctree;
    OctreePeriodicExtension<OctreeData<T> > *octreePeriodicExtension =
        new OctreePeriodicExtension<OctreeData<T> >(root, xPeriodic, yPeriodic, zPeriodic);
    root = octreePeriodicExtension->get();
    iter = 0;
    do {
        if (verbose) {
            (*outS) << "    Octree calibration iteration: " << iter << std::endl;
        }
        calibrateOctree.reset();
        processOctreePostOrder(root, calibrateOctree);
        iter++;
    } while (calibrateOctree.calibratedOctree() && iter < maxIter);
    root = octreePeriodicExtension->release();
    delete octreePeriodicExtension;
    octreePeriodicExtension = 0;

    if (calibrateOctree.calibratedOctree()) {
        (*errS) << std::endl;
        (*errS) << "Octree calibration failed after: " << maxIter << " iterations." << std::endl;
        exit(1);
    }
    if (verbose) {
        (*outS) << std::endl;
    }

    if (verbose) {
        countOctreeLeafNodes.reset();
        processOctreePostOrder(root, countOctreeLeafNodes);
        (*outS) << "Number of leaf nodes after calibration: "
                << countOctreeLeafNodes.getNumLeafNodes() << std::endl;
        (*outS) << std::endl;
    }

    int newMinLeafLevel = std::numeric_limits<int>::max();
    int newMaxLeafLevel = -1;
    getAllocatedMinMaxOctreeLeafNodeLevels(root, newMinLeafLevel, newMaxLeafLevel);
    if (verbose) {
        (*outS) << "Minimum leaf level: " << newMinLeafLevel << std::endl;
        (*outS) << "Maximum leaf level: " << newMaxLeafLevel << std::endl;
        (*outS) << std::endl;
    }

    minLeafLevel = newMinLeafLevel;
    maxLeafLevel = newMaxLeafLevel;

    if (stlOutput) {
        if (verbose) {
            (*outS) << "Saving the octree after calibration... " << std::flush;
        }
        triangleSet = new TriangleSet<T>(DBL);
        gatherOctreeLeafCuboids.reset();
        processOctreePostOrder(root, gatherOctreeLeafCuboids);
        for (plint iCuboid = 0; iCuboid < (plint)gatherOctreeLeafCuboids.getCuboids().size();
             iCuboid++) {
            Array<plint, 3> nSegments(1, 1, 1);
            TriangleSet<T> cuboidSet =
                constructCuboid<T>(gatherOctreeLeafCuboids.getCuboids()[iCuboid], nSegments);
            triangleSet->append(cuboidSet);
        }
        triangleSet->writeBinarySTL(outDir + stlBaseName + "_after_calibration.stl");
        delete triangleSet;
        triangleSet = 0;
        if (verbose) {
            (*outS) << "done." << std::endl;
            (*outS) << std::endl;
        }
    }

    // Octree node removal based on negative grid density values.

    if (gridDensityHasNegativeValues && removeBlocks) {
        if (verbose) {
            (*outS) << "Removing octree nodes based on negative grid density values." << std::endl;
        }
        bool boundFromBelow = false;  // Now we want to retain negative values.
        GridDensityFunction<T> *gridDensityFunction =
            new GridDensityFunction<T>(gridDensityFunctionFile, boundFromBelow);
        DeallocateOctreeLeafNodes<T> deallocateOctreeLeafNodes(
            *gridDensityFunction, maxLeafLevel, useSamples, numSamples);
        processOctreePostOrder(root, deallocateOctreeLeafNodes);
        if (verbose) {
            (*outS) << std::endl;
        }

        delete gridDensityFunction;
        gridDensityFunction = 0;

        if (verbose) {
            countOctreeLeafNodes.reset();
            processOctreePostOrder(root, countOctreeLeafNodes);
            (*outS) << "Number of leaf nodes after node removal: "
                    << countOctreeLeafNodes.getNumLeafNodes() << std::endl;
            (*outS) << std::endl;
        }

        int newMinLeafLevel = std::numeric_limits<int>::max();
        int newMaxLeafLevel = -1;
        getAllocatedMinMaxOctreeLeafNodeLevels(root, newMinLeafLevel, newMaxLeafLevel);
        if (verbose) {
            (*outS) << "Minimum leaf level: " << newMinLeafLevel << std::endl;
            (*outS) << "Maximum leaf level: " << newMaxLeafLevel << std::endl;
            (*outS) << std::endl;
        }

        minLeafLevel = newMinLeafLevel;
        maxLeafLevel = newMaxLeafLevel;

        if (stlOutput) {
            if (verbose) {
                (*outS) << "Saving the octree after node removal... " << std::flush;
            }
            triangleSet = new TriangleSet<T>(DBL);
            gatherOctreeLeafCuboids.reset();
            processOctreePostOrder(root, gatherOctreeLeafCuboids);
            for (plint iCuboid = 0; iCuboid < (plint)gatherOctreeLeafCuboids.getCuboids().size();
                 iCuboid++) {
                Array<plint, 3> nSegments(1, 1, 1);
                TriangleSet<T> cuboidSet =
                    constructCuboid<T>(gatherOctreeLeafCuboids.getCuboids()[iCuboid], nSegments);
                triangleSet->append(cuboidSet);
            }
            triangleSet->writeBinarySTL(outDir + stlBaseName + "_after_node_removal.stl");
            delete triangleSet;
            triangleSet = 0;
            if (verbose) {
                (*outS) << "done." << std::endl;
                (*outS) << std::endl;
            }
        }
    }

    // From physical to lattice units.

    if (verbose) {
        (*outS) << "Passing to lattice units... " << std::flush;
    }

    T Lx = fullDomain.x1() - fullDomain.x0();
    T Ly = fullDomain.y1() - fullDomain.y0();
    T Lz = fullDomain.z1() - fullDomain.z0();

    T Lmin = std::min(std::min(Lx, Ly), Lz);

    T lx = Lx / (T)util::intTwoToThePower(maxLeafLevel);
    T ly = Ly / (T)util::intTwoToThePower(maxLeafLevel);
    T lz = Lz / (T)util::intTwoToThePower(maxLeafLevel);

    T lmin = Lmin / (T)util::intTwoToThePower(maxLeafLevel);

    dxFinest = lmin / (T)n;
    T eps = (T)0.01 * dxFinest;

    T floatNx = (T)n * lx / lmin;
    plint intNx = util::roundToInt(floatNx);
    T diffNx = std::fabs(floatNx - (T)intNx);
    bool problemWithX = !util::isZero(diffNx, eps) || intNx % 2 != 0;

    T floatNy = (T)n * ly / lmin;
    plint intNy = util::roundToInt(floatNy);
    T diffNy = std::fabs(floatNy - (T)intNy);
    bool problemWithY = !util::isZero(diffNy, eps) || intNy % 2 != 0;

    T floatNz = (T)n * lz / lmin;
    plint intNz = util::roundToInt(floatNz);
    T diffNz = std::fabs(floatNz - (T)intNz);
    bool problemWithZ = !util::isZero(diffNz, eps) || intNz % 2 != 0;

    if (problemWithX || problemWithY || problemWithZ) {
        if (verbose) {
            (*errS) << std::endl;
            (*errS) << std::endl;
            (*errS) << "The Octree grid cannot be constructed: Problem with the simulation domain "
                       "and the block size."
                    << std::endl;
            (*errS) << "    Current n (block size): " << n << std::endl;
            if (problemWithX) {
                plint newIntNx = intNx;
                if (newIntNx % 2 != 0) {
                    newIntNx++;
                }
                T newLx = (T)newIntNx * Lmin / (T)n;
                (*errS)
                    << "    Current Lx (length of the simulation domain in the x-axis direction): "
                    << Lx << std::endl;
                (*errS) << "    Use, for example: Lx = " << newLx << std::endl;
            }
            if (problemWithY) {
                plint newIntNy = intNy;
                if (newIntNy % 2 != 0) {
                    newIntNy++;
                }
                T newLy = (T)newIntNy * Lmin / (T)n;
                (*errS)
                    << "    Current Ly (length of the simulation domain in the y-axis direction): "
                    << Ly << std::endl;
                (*errS) << "    Use, for example: Ly = " << newLy << std::endl;
            }
            if (problemWithZ) {
                plint newIntNz = intNz;
                if (newIntNz % 2 != 0) {
                    newIntNz++;
                }
                T newLz = (T)newIntNz * Lmin / (T)n;
                (*errS)
                    << "    Current Lz (length of the simulation domain in the z-axis direction): "
                    << Lz << std::endl;
                (*errS) << "    Use, for example: Lz = " << newLz << std::endl;
            }
        } else {
            (*errS) << std::endl;
            (*errS) << std::endl;
            (*errS) << "The Octree grid cannot be constructed: Problem with the simulation domain "
                       "and the block size."
                    << std::endl;
        }

        exit(1);
    }

    nx = intNx;
    ny = intNy;
    nz = intNz;

    if (verbose) {
        (*outS) << "done." << std::endl;
    }

    // Find and tag the overlaps (first we touch all nodes of the tree to optimize).

    if (verbose) {
        (*outS) << "Defining overlaps... " << std::flush;
    }
    TouchAllOctreeNodes<T> touchAllOctreeNodes;
    processOctreePostOrder(root, touchAllOctreeNodes);

    AssignOctreeOverlapNodes<T> assignOctreeOverlapNodes;
    octreePeriodicExtension =
        new OctreePeriodicExtension<OctreeData<T> >(root, xPeriodic, yPeriodic, zPeriodic);
    root = octreePeriodicExtension->get();
    processOctreePostOrder(root, assignOctreeOverlapNodes);
    root = octreePeriodicExtension->release();
    delete octreePeriodicExtension;
    octreePeriodicExtension = 0;

    CheckOctreeOverlapNodes<T> checkOctreeOverlapNodes;
    octreePeriodicExtension =
        new OctreePeriodicExtension<OctreeData<T> >(root, xPeriodic, yPeriodic, zPeriodic);
    root = octreePeriodicExtension->get();
    processOctreePreOrder(root, checkOctreeOverlapNodes);
    root = octreePeriodicExtension->release();
    delete octreePeriodicExtension;
    octreePeriodicExtension = 0;
    if (verbose) {
        (*outS) << "done." << std::endl;
    }

    // Compute bulks for nodes to be allocated.

    if (verbose) {
        (*outS) << "Computing bulks... " << std::flush;
    }
    ComputeOctreeBulks<T> computeOctreeBulks(dxFinest, maxLeafLevel, fullDomain, nx, ny, nz);
    processOctreePostOrder(root, computeOctreeBulks);
    if (verbose) {
        (*outS) << "done." << std::endl;
    }

    // Assign leaf and overlap node ids.

    if (verbose) {
        (*outS) << "Assigning block ids... " << std::flush;
    }
    AssignOctreeLeafAndOverlapIds<T> assignOctreeLeafAndOverlapIds;
    processOctreePostOrder(root, assignOctreeLeafAndOverlapIds);
    if (verbose) {
        (*outS) << "done." << std::endl;
    }

    // Assign process ids.

    plint numCellsPerBlock = nx * ny * nz;

    OctreeProcessLoads octreeProcessLoads(numProcesses, numCellsPerBlock, minLeafLevel);

    std::string directionPhrase =
        (fineToCoarse ? "from fine to coarse levels" : "from coarse to fine levels");

    int numBlocksPerGroup = std::pow(8, numLevelsToGroupOverlaps);
    std::string blockWord = (numBlocksPerGroup == 1 ? "block" : "blocks");

    if (verbose) {
        (*outS) << "Assigning overlap process ids " << directionPhrase << ", in groups of "
                << numBlocksPerGroup << " " << blockWord << "... " << std::flush;
    }
    if (fineToCoarse) {
        for (int level = maxLeafLevel - 1; level >= minLeafLevel; level--) {
            AssignOctreeOverlapProcessIds<T> assignOctreeOverlapProcessIds(
                level, numLevelsToGroupOverlaps, octreeProcessLoads);
            processOctreePreOrder(root, assignOctreeOverlapProcessIds);
        }
    } else {
        for (int level = minLeafLevel; level < maxLeafLevel; level++) {
            AssignOctreeOverlapProcessIds<T> assignOctreeOverlapProcessIds(
                level, numLevelsToGroupOverlaps, octreeProcessLoads);
            processOctreePreOrder(root, assignOctreeOverlapProcessIds);
        }
    }
    if (verbose) {
        (*outS) << "done." << std::endl;
    }

    if (verbose) {
        (*outS) << "Enforcing process assignment constraints... " << std::flush;
    }
    EnforceOctreeProcessAssignmentConstraints<T> enforceOctreeProcessAssignmentConstraints(
        octreeProcessLoads);
    octreePeriodicExtension =
        new OctreePeriodicExtension<OctreeData<T> >(root, xPeriodic, yPeriodic, zPeriodic);
    root = octreePeriodicExtension->get();
    processOctreePreOrder(root, enforceOctreeProcessAssignmentConstraints);
    root = octreePeriodicExtension->release();
    delete octreePeriodicExtension;
    octreePeriodicExtension = 0;
    if (verbose) {
        (*outS) << "done." << std::endl;
    }

    numBlocksPerGroup = std::pow(8, numLevelsToGroupBlocks);
    blockWord = (numBlocksPerGroup == 1 ? "block" : "blocks");

    if (verbose) {
        (*outS) << "Assigning normal block process ids " << directionPhrase << ", in groups of "
                << numBlocksPerGroup << " " << blockWord << "... " << std::flush;
    }
    if (fineToCoarse) {
        for (int level = maxLeafLevel; level >= minLeafLevel; level--) {
            AssignOctreeLeafProcessIds<T> assignLeafLeafProcessIds(
                level, numLevelsToGroupBlocks, strongGrouping, octreeProcessLoads);
            processOctreePreOrder(root, assignLeafLeafProcessIds);
        }
    } else {
        for (int level = minLeafLevel; level <= maxLeafLevel; level++) {
            AssignOctreeLeafProcessIds<T> assignLeafLeafProcessIds(
                level, numLevelsToGroupBlocks, strongGrouping, octreeProcessLoads);
            processOctreePreOrder(root, assignLeafLeafProcessIds);
        }
    }
    if (verbose) {
        (*outS) << "done." << std::endl;
    }

    // Visualize the leaf and overlap bulks in the finest grid level lattice units.

    if (stlOutput) {
        if (verbose) {
            (*outS) << "Saving the grid in lattice units... " << std::flush;
        }
        bool gatherOverlaps = false;  // We are interested in leaf nodes at the moment.
        int whichLevel = -1;          // All levels.
        plint whichProcess = -1;      // All processes.
        GatherOctreeCuboidsInLatticeUnits<T> gatherOctreeCuboidsInLatticeUnits(
            gatherOverlaps, maxLeafLevel, whichLevel, whichProcess);
        processOctreePostOrder(root, gatherOctreeCuboidsInLatticeUnits);
        bool highResolution = false;
        triangleSet = new TriangleSet<T>(DBL);
        for (plint iCuboid = 0;
             iCuboid < (plint)gatherOctreeCuboidsInLatticeUnits.getCuboids().size(); iCuboid++)
        {
            Array<plint, 3> nSegments =
                (highResolution ? Array<plint, 3>(nx - 1, ny - 1, nz - 1)
                                : Array<plint, 3>(1, 1, 1));
            TriangleSet<T> cuboidSet = constructCuboid<T>(
                gatherOctreeCuboidsInLatticeUnits.getCuboids()[iCuboid], nSegments);
            triangleSet->append(cuboidSet);
        }
        triangleSet->writeBinarySTL(outDir + stlBaseName + "_in_finest_lattice_units.stl");
        delete triangleSet;
        triangleSet = 0;

        highResolution = false;
        for (plint iProcess = 0; iProcess < numProcesses; iProcess++) {
            bool gatherOverlaps = false;  // We are interested in leaf nodes at the moment.
            int whichLevel = -1;          // All levels.
            plint whichProcess = iProcess;
            GatherOctreeCuboidsInLatticeUnits<T> gatherOctreeCuboidsInLatticeUnits(
                gatherOverlaps, maxLeafLevel, whichLevel, whichProcess);
            processOctreePostOrder(root, gatherOctreeCuboidsInLatticeUnits);

            triangleSet = new TriangleSet<T>(DBL);
            for (plint iCuboid = 0;
                 iCuboid < (plint)gatherOctreeCuboidsInLatticeUnits.getCuboids().size(); iCuboid++)
            {
                Array<plint, 3> nSegments =
                    (highResolution ? Array<plint, 3>(nx - 1, ny - 1, nz - 1)
                                    : Array<plint, 3>(1, 1, 1));
                TriangleSet<T> cuboidSet = constructCuboid<T>(
                    gatherOctreeCuboidsInLatticeUnits.getCuboids()[iCuboid], nSegments);
                triangleSet->append(cuboidSet);
            }
            std::string fname =
                createFileName(
                    outDir + stlBaseName + "_in_finest_lattice_units_process_", iProcess, 8)
                + ".stl";
            triangleSet->writeBinarySTL(fname);
            delete triangleSet;
            triangleSet = 0;
        }

        highResolution = false;
        for (int level = minLeafLevel; level <= maxLeafLevel; level++) {
            bool gatherOverlaps = false;  // We are interested in leaf nodes at the moment.
            int whichLevel = level;
            plint whichProcess = -1;  // All processes.
            GatherOctreeCuboidsInLatticeUnits<T> gatherOctreeCuboidsInLatticeUnits(
                gatherOverlaps, maxLeafLevel, whichLevel, whichProcess);
            processOctreePostOrder(root, gatherOctreeCuboidsInLatticeUnits);

            triangleSet = new TriangleSet<T>(DBL);
            for (plint iCuboid = 0;
                 iCuboid < (plint)gatherOctreeCuboidsInLatticeUnits.getCuboids().size(); iCuboid++)
            {
                Array<plint, 3> nSegments =
                    (highResolution ? Array<plint, 3>(nx - 1, ny - 1, nz - 1)
                                    : Array<plint, 3>(1, 1, 1));
                TriangleSet<T> cuboidSet = constructCuboid<T>(
                    gatherOctreeCuboidsInLatticeUnits.getCuboids()[iCuboid], nSegments);
                triangleSet->append(cuboidSet);
            }
            std::string fname =
                createFileName(outDir + stlBaseName + "_in_finest_lattice_units_level_", level, 2)
                + ".stl";
            triangleSet->writeBinarySTL(fname);
            delete triangleSet;
            triangleSet = 0;
        }

        highResolution = false;
        for (int level = minLeafLevel; level < maxLeafLevel; level++) {
            bool gatherOverlaps = true;  // We are interested in overlaps at the moment.
            int whichLevel = level;
            plint whichProcess = -1;  // All processes.
            GatherOctreeCuboidsInLatticeUnits<T> gatherOctreeCuboidsInLatticeUnits(
                gatherOverlaps, maxLeafLevel, whichLevel, whichProcess);
            processOctreePostOrder(root, gatherOctreeCuboidsInLatticeUnits);

            triangleSet = new TriangleSet<T>(DBL);
            for (plint iCuboid = 0;
                 iCuboid < (plint)gatherOctreeCuboidsInLatticeUnits.getCuboids().size(); iCuboid++)
            {
                Array<plint, 3> nSegments =
                    (highResolution ? Array<plint, 3>(nx - 1, ny - 1, nz - 1)
                                    : Array<plint, 3>(1, 1, 1));
                TriangleSet<T> cuboidSet = constructCuboid<T>(
                    gatherOctreeCuboidsInLatticeUnits.getCuboids()[iCuboid], nSegments);
                triangleSet->append(cuboidSet);
            }
            std::string fname =
                createFileName(
                    outDir + stlBaseName + "_overlaps_in_finest_lattice_units_level_", level, 2)
                + ".stl";
            triangleSet->writeBinarySTL(fname);
            delete triangleSet;
            triangleSet = 0;
        }

        if (verbose) {
            (*outS) << "done." << std::endl;
            (*outS) << std::endl;
        }
    }

    // Write grid statistics.

    ComputeOctreeGridStatistics<T> computeOctreeGridStatistics(
        minLeafLevel, maxLeafLevel, numProcesses, numCellsPerBlock);
    processOctreePostOrder(root, computeOctreeGridStatistics);

    if (verbose) {
        (*outS) << "----------------------" << std::endl;
        (*outS) << "Octree Grid Statistics" << std::endl;
        (*outS) << "----------------------" << std::endl;
        (*outS) << std::endl;

        plint numLeafBlocks = 0;
        plint numOverlapBlocks = 0;
        plint numBlocks = 0;
        for (int level = minLeafLevel; level <= maxLeafLevel; level++) {
            plint numLeafBlocksOnLevel = computeOctreeGridStatistics.getNumLeafBlocksOnLevel(level);
            plint numOverlapBlocksOnLevel =
                computeOctreeGridStatistics.getNumOverlapBlocksOnLevel(level);
            plint numBlocksOnLevel = numLeafBlocksOnLevel + numOverlapBlocksOnLevel;

            numLeafBlocks += numLeafBlocksOnLevel;
            numOverlapBlocks += numOverlapBlocksOnLevel;
            numBlocks += numBlocksOnLevel;

            (*outS) << "On level: " << level << std::endl;
            (*outS) << "    Number of blocks (not overlaps)          : " << numLeafBlocksOnLevel
                    << std::endl;
            (*outS) << "    Number of overlap blocks                 : " << numOverlapBlocksOnLevel
                    << std::endl;
            (*outS) << "    Total number of blocks                   : " << numBlocksOnLevel
                    << std::endl;
            (*outS) << "    Number of lattice nodes (not in overlaps): "
                    << numLeafBlocksOnLevel * numCellsPerBlock << std::endl;
            (*outS) << "    Number of overlap lattice nodes          : "
                    << numOverlapBlocksOnLevel * numCellsPerBlock << std::endl;
            (*outS) << "    Total number of lattice nodes            : "
                    << numBlocksOnLevel * numCellsPerBlock << std::endl;
            (*outS) << std::endl;
            for (plint iProcess = 0; iProcess < numProcesses; iProcess++) {
                plint numBlocksOnLevelOnProcess =
                    computeOctreeGridStatistics.getNumBlocksOnLevelOnProcess(level, iProcess);
                (*outS) << "    On process: " << std::right << std::setw(3) << iProcess
                        << "; Total number of blocks: " << std::setw(3) << numBlocksOnLevelOnProcess
                        << "; Relative load: "
                        << computeOctreeGridStatistics.getRelativeLoadOnLevelOnProcess(
                               level, iProcess)
                        << std::endl;
            }
            (*outS) << std::endl;
            (*outS) << "    Mean of relative load per process on this level              : " << 1
                    << std::endl;
            (*outS) << "    Standard deviation of relative load per process on this level: "
                    << computeOctreeGridStatistics.getStandardDeviationOfRelativeLoadOnLevel(level)
                    << std::endl;
            (*outS) << "    Range of relative load per process (max - min) on this level : "
                    << computeOctreeGridStatistics.getRangeOfRelativeLoadOnLevel(level)
                    << std::endl;
            (*outS) << std::endl;
        }

        for (plint iProcess = 0; iProcess < numProcesses; iProcess++) {
            plint numBlocksOnProcess = computeOctreeGridStatistics.getNumBlocksOnProcess(iProcess);
            (*outS) << "On process: " << iProcess << std::endl;
            (*outS) << "    Total number of blocks               : " << numBlocksOnProcess
                    << std::endl;
            (*outS) << "    Total number of cells                : "
                    << numBlocksOnProcess * numCellsPerBlock << std::endl;
            (*outS) << "    Relative load (all levels considered): "
                    << computeOctreeGridStatistics.getRelativeLoadOnProcess(iProcess) << std::endl;
        }
        (*outS) << std::endl;

        (*outS) << "All levels considered: " << std::endl;
        (*outS) << "    Number of blocks (not overlaps)          : " << numLeafBlocks << std::endl;
        (*outS) << "    Number of overlap blocks                 : " << numOverlapBlocks
                << std::endl;
        (*outS) << "    Total number of blocks                   : " << numBlocks << std::endl;
        (*outS) << "    Number of lattice nodes (not in overlaps): "
                << numLeafBlocks * numCellsPerBlock << std::endl;
        (*outS) << "    Number of overlap lattice nodes          : "
                << numOverlapBlocks * numCellsPerBlock << std::endl;
        (*outS) << "    Total number of lattice nodes            : " << numBlocks * numCellsPerBlock
                << std::endl;
        (*outS) << std::endl;

        (*outS) << "All processes considered: " << std::endl;
        (*outS) << "    Mean of relative load per process                                      : "
                << 1 << std::endl;
        (*outS) << "    Standard deviation of relative load per process (all levels considered): "
                << computeOctreeGridStatistics.getStandardDeviationOfRelativeLoad() << std::endl;
        (*outS) << "    Range of relative load per process (max - min)  (all levels considered): "
                << computeOctreeGridStatistics.getRangeOfRelativeLoad() << std::endl;

        (*outS) << std::endl;

        (*outS) << "At finest level  : dx = " << dxFinest << std::endl;
        (*outS) << "At coarsest level: dx = "
                << dxFinest * util::intTwoToThePower(maxLeafLevel - minLeafLevel) << std::endl;
        (*outS) << std::endl;

        (*outS) << "Uniform block size: " << nx << " x " << ny << " x " << nz << std::endl;
        (*outS) << std::endl;
    }

    bool allProcessesHaveBlocks = true;
    for (plint iProcess = 0; iProcess < numProcesses; iProcess++) {
        plint numBlocksOnProcess = computeOctreeGridStatistics.getNumBlocksOnProcess(iProcess);
        if (numBlocksOnProcess == 0) {
            allProcessesHaveBlocks = false;
            break;
        }
    }
    if (!allProcessesHaveBlocks) {
        (*errS) << std::endl;
        (*errS)
            << "The Octree grid cannot be used by Palabos: there are processes that have no blocks."
            << std::endl;
        exit(1);
    }

    // Construct the OctreeGridStructure.

    if (verbose) {
        (*outS) << "Constructing the octree grid structure... " << std::flush;
    }
    OctreeGridStructure octreeGridStructure;
    ConstructOctreeGridStructure<T> constructOctreeGridStructure(octreeGridStructure, minLeafLevel);
    octreePeriodicExtension =
        new OctreePeriodicExtension<OctreeData<T> >(root, xPeriodic, yPeriodic, zPeriodic);
    root = octreePeriodicExtension->get();
    processOctreePostOrder(root, constructOctreeGridStructure);
    root = octreePeriodicExtension->release();
    delete octreePeriodicExtension;
    octreePeriodicExtension = 0;
    if (verbose) {
        (*outS) << "done." << std::endl;
        (*outS) << std::endl;

        (*outS) << "Level bounding boxes (in lattice units of the respective level):" << std::endl;
        for (int level = minLeafLevel; level <= maxLeafLevel; level++) {
            (*outS) << "    Octree level " << level << ": ";
            Box3D bb = octreeGridStructure.getBoundingBox(level - minLeafLevel);
            (*outS) << " bounding box = [" << bb.x0 << ", " << bb.x1 << "; " << bb.y0 << ", "
                    << bb.y1 << "; " << bb.z0 << ", " << bb.z1 << "] " << std::endl;
        }
        (*outS) << std::endl;
    }

    freeOctree(root);

    return (octreeGridStructure);
}

template <typename T>
Cuboid<T> const &OctreeGridGenerator<T>::getFullDomain() const
{
    return (fullDomain);
}

template <typename T>
int OctreeGridGenerator<T>::getMinLevel() const
{
    return (minLeafLevel);
}

template <typename T>
int OctreeGridGenerator<T>::getMaxLevel() const
{
    return (maxLeafLevel);
}

template <typename T>
plint OctreeGridGenerator<T>::getNumProcesses() const
{
    return (numProcesses);
}

template <typename T>
plint OctreeGridGenerator<T>::getBlockNx() const
{
    return (nx);
}

template <typename T>
plint OctreeGridGenerator<T>::getBlockNy() const
{
    return (ny);
}

template <typename T>
plint OctreeGridGenerator<T>::getBlockNz() const
{
    return (nz);
}

template <typename T>
T OctreeGridGenerator<T>::getDxFinestLevel() const
{
    return (dxFinest);
}

template <typename T>
std::string OctreeGridGenerator<T>::getOutDir() const
{
    return (outDir);
}

template <typename T>
void OctreeGridGenerator<T>::setOutputAndErrorStreams(
    Parallel_ostream &outStream, Parallel_ostream &errStream)
{
    outS = &outStream;
    errS = &errStream;
}

}  // namespace plb

#endif  // OCTREE_GRID_GENERATOR_HH
