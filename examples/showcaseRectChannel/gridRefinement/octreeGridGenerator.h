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

#ifndef OCTREE_GRID_GENERATOR_H
#define OCTREE_GRID_GENERATOR_H

#include <vector>

#include "core/array.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "gridRefinement/octree.h"
#include "gridRefinement/octreeGridStructure.h"
#include "io/parallelIO.h"
#include "libraryInterfaces/TINYXML_xmlIO.h"
#include "multiBlock/multiDataField3D.h"
#include "multiBlock/multiDataField3D.hh"

namespace plb {

template <typename T>
struct OctreeData {
    OctreeData(Cuboid<T> const &cuboid_);

    Cuboid<T> cuboid;
    bool touch;  // This is a variable to be used for optimizing
                 // the certain processes (like refinement and
                 // finding the overlap blocks). If its value is
                 // true, this means that a processing function
                 // may "touch" the node, and do something with
                 // its data.

    Box3D bulk;       // This is the bulk of the Palabos blocks. The bulk is
                      // not needed during the process of the tree creation
                      // (this is why it does not appear in the constructor).
    plint blockId;    // After the full tree is constructed, we need to
                      // have a unique id for each leaf and overlap node,
                      // that will become blocks for Palabos.
    bool isOverlap;   // A boolean to determine if a non-leaf cell is an
                      // overlap needed for the coupling with the next
                      // grid level.
    bool allocate;    // A boolean for each cell to determine if it will
                      // be allocated and included in the grid structure.
                      // This variable has a meaning only for the leaf
                      // nodes and is set to true for all non-leaves.
                      // The user will define regions where nodes will not
                      // be allocated, by providing negative grid density
                      // values. If a block is "fully contained" in a
                      // region of space with negative grid density,
                      // then it will not be allocated. In the majority
                      // of the following algorithms, a leaf node is
                      // considered only if it is allocated (for example
                      // only the allocated leaves are refined, counted,
                      // etc).
    plint processId;  // MPI process id.
};

namespace octreeDataUtil {

template <typename T>
inline bool isAllocatableLeaf(OctreeNode<OctreeData<T> > *node);

template <typename T>
inline bool isOverlap(OctreeNode<OctreeData<T> > *node);

template <typename T>
inline bool isAllocatable(OctreeNode<OctreeData<T> > *node);

template <typename T>
Cuboid<T> getOctant(Cuboid<T> const &cuboid, int i);

template <typename T>
void getAllocatedMinMaxOctreeLeafNodeLevels(
    OctreeNode<OctreeData<T> > *root, int &minLevel, int &maxLevel);

}  // namespace octreeDataUtil

template <typename T>
struct AssignOctreeData {
    void operator()(OctreeNode<OctreeData<T> > *node);
};

template <typename T>
class GridDensityFunction {
public:
    GridDensityFunction(std::string fileName, bool boundFromBelow);
    GridDensityFunction(
        Cuboid<T> const &cuboid_, T dx_, ScalarField3D<T> const &gridDensity_, bool boundFromBelow);
    GridDensityFunction(GridDensityFunction<T> const &rhs);
    void swap(GridDensityFunction<T> &rhs);
    GridDensityFunction<T> &operator=(GridDensityFunction<T> const &rhs);
    GridDensityFunction<T> *clone() const;
    ~GridDensityFunction();
    bool hasNegativeValues() const;
    T operator()(Array<T, 3> const &pos) const;
    T operator()(Cuboid<T> const &c) const;

private:
    Cuboid<T> cuboid;
    T dx;
    plint nx, ny, nz;
    ScalarField3D<T> *gridDensity;
    bool existNegativeValues;
};

template <typename T>
class RefineOctree {
public:
    RefineOctree(
        GridDensityFunction<T> const &gridDensityFunction_, T gridDensityScaleFactor_,
        int minLeafLevel_, int maxLeafLevel_, bool useSamples_, plint numSamples_);
    void reset();
    bool refinedOctree() const;
    void operator()(OctreeNode<OctreeData<T> > *node);

private:
    T getMaxScaledGridDensityFromPoints(Cuboid<T> const &cuboid, plint n) const;
    T getMaxScaledGridDensityFromVolumes(Cuboid<T> const &cuboid) const;

private:
    AssignOctreeData<T> assignOctreeData;
    GridDensityFunction<T> const &gridDensityFunction;
    T gridDensityScaleFactor;
    int minLeafLevel;
    int maxLeafLevel;
    bool useSamples;
    plint numSamples;
    bool refined;
};

template <typename T>
class CalibrateOctree {
public:
    CalibrateOctree();
    void reset();
    bool calibratedOctree() const;
    void operator()(OctreeNode<OctreeData<T> > *node);

private:
    AssignOctreeData<T> assignOctreeData;
    bool calibrated;
};

template <typename T>
class DeallocateOctreeLeafNodes {
public:
    DeallocateOctreeLeafNodes(
        GridDensityFunction<T> const &gridDensityFunction_, int maxLeafLevel_, bool useSamples_,
        plint numSamples_);
    void operator()(OctreeNode<OctreeData<T> > *node);

private:
    bool isMaxGridDensityNegativeOnPoints(Cuboid<T> const &cuboid, plint n) const;
    bool isMaxGridDensityNegativeInVolumes(Cuboid<T> const &cuboid) const;

private:
    GridDensityFunction<T> const &gridDensityFunction;
    int maxLeafLevel;
    bool useSamples;
    plint numSamples;
};

template <typename T>
class CountOctreeLeafNodes {
public:
    CountOctreeLeafNodes();
    void reset();
    plint getNumLeafNodes() const;
    void operator()(OctreeNode<OctreeData<T> > *node);

private:
    plint numLeaves;
};

template <typename T>
class GatherOctreeLeafCuboids {
public:
    void reset();
    std::vector<Cuboid<T> > const &getCuboids() const;
    void operator()(OctreeNode<OctreeData<T> > *node);

private:
    std::vector<Cuboid<T> > cuboids;
};

template <typename T>
struct TouchAllOctreeNodes {
    void operator()(OctreeNode<OctreeData<T> > *node);
};

template <typename T>
struct AssignOctreeOverlapNodes {
    void operator()(OctreeNode<OctreeData<T> > *node);
};

template <typename T>
struct CheckOctreeOverlapNodes {
    void operator()(OctreeNode<OctreeData<T> > *node);
};

template <typename T>
class ComputeOctreeBulks {
public:
    ComputeOctreeBulks(
        T dxFinest_, int maxLeafLevel_, Cuboid<T> const &fullDomain_, plint nx_, plint ny_,
        plint nz_);
    void operator()(OctreeNode<OctreeData<T> > *node);

private:
    T dxFinest;
    int maxLeafLevel;
    Cuboid<T> fullDomain;
    plint nx, ny, nz;
};

template <typename T>
class AssignOctreeLeafAndOverlapIds {
public:
    AssignOctreeLeafAndOverlapIds();
    void operator()(OctreeNode<OctreeData<T> > *node);

private:
    plint blockId;
};

// We consider overlap blocks to have the same "weight" as normal blocks, although there
// exist optimizations to make the overlap blocks smaller. Due to this fact there is an
// error in the estimation of the process loads.
class OctreeProcessLoads {
public:
    OctreeProcessLoads(plint numProcesses_, plint numCellsPerBlock_, int minLeafLevel_);
    void clearLoads();
    void clearLoad(plint processId);
    plint getNumProcesses() const;
    plint getNumCellsPerBlock() const;
    int getMinLeafLevel() const;
    plint getLoad(plint processId) const;
    void addLoad(plint processId, int level);
    plint getProcessWithMinLoad() const;

private:
    plint numProcesses;
    plint numCellsPerBlock;
    int minLeafLevel;
    std::vector<plint> loads;
};

template <typename T>
class AssignOctreeOverlapProcessIds {
public:
    AssignOctreeOverlapProcessIds(
        int level_, int numLevelsToGroupOverlaps_, OctreeProcessLoads &octreeProcessLoads_);
    void operator()(OctreeNode<OctreeData<T> > *node);

private:
    bool needsProcessId(OctreeNode<OctreeData<T> > *node);
    void assignProcessIdsPreOrder(OctreeNode<OctreeData<T> > *root);

private:
    int level;
    int numLevelsToGroupOverlaps;
    OctreeProcessLoads &octreeProcessLoads;
    int groupingLevel;
    plint singleProcessId;
};

template <typename T>
class EnforceOctreeProcessAssignmentConstraints {
public:
    EnforceOctreeProcessAssignmentConstraints(OctreeProcessLoads &octreeProcessLoads_);
    void operator()(OctreeNode<OctreeData<T> > *node);

private:
    OctreeProcessLoads &octreeProcessLoads;
};

template <typename T>
class AssignOctreeLeafProcessIds {
public:
    AssignOctreeLeafProcessIds(
        int level_, int numLevelsToGroupBlocks_, bool strongGrouping_,
        OctreeProcessLoads &octreeProcessLoads_);
    void operator()(OctreeNode<OctreeData<T> > *node);

private:
    bool needsProcessId(OctreeNode<OctreeData<T> > *node);
    void determineSingleOverlapProcessIdPreOrder(OctreeNode<OctreeData<T> > *root);
    void assignProcessIdsPreOrder(OctreeNode<OctreeData<T> > *root);

private:
    int level;
    int numLevelsToGroupBlocks;
    bool strongGrouping;
    OctreeProcessLoads &octreeProcessLoads;
    int overlapLevel;
    int groupingLevel;
    plint singleProcessId;
};

template <typename T>
class GatherOctreeCuboidsInLatticeUnits {
public:
    GatherOctreeCuboidsInLatticeUnits(
        bool gatherOverlaps_, int maxLeafLevel_, int whichLevel_, plint whichProcess_);
    std::vector<Cuboid<T> > const &getCuboids() const;
    void operator()(OctreeNode<OctreeData<T> > *node);

private:
    bool gatherOverlaps;
    int maxLeafLevel, whichLevel;
    plint whichProcess;
    std::vector<Cuboid<T> > cuboids;
};

// Load Balancing.
//
// We define the load on a process to be the number of cells (of each level) that belong to this
// process, multiplied by the number of times "collideAndStream" is called for the each specific
// level inside one lattice Boltzmann cycle. Then, the relative load on a process is the actual load
// over the average load, so perfect load balancing corresponds to a relative load of 1.
//
// In the octree grid generation, the number of cells per block is uniform over all blocks and
// levels, so it factors out when the relative load is computed. We keep it in the calculations of
// the load, because one day we might be interested in the load itself (not its relative value).
//
// Different levels might have different MultiBlocks with different envelope widths, and also
// different data processors associated. For instance, when the immersed boundary method is used,
// the finest level has rhoBar and j with envelope width of 4, plus a particle field with an
// envelope width of 2. In such a case, the finest level has also more data processors integrated,
// to perform the immersed boundary iterations. For this reason, load balancing is required "per
// level", as well as taking all levels under consideration (with their relative "temporal" weight).
//
// In the grid statistics, we do not consider the envelope widths, because they are not part of the
// grid generation and process attribution proceedures. As noted above, there might be many
// different MultiBlocks with many different envelopes, but the grid generator has no control over
// this. Because in the octree domain decomposition, all blocks have the same number of cells, it
// would be good that the user does not choose small bulk sizes when large envelope widths are
// expected, since this will have a penalty on the global memory footprint during execution.
//
// In the grid statistics, we do not consider the "communication cost", meaning the cost of the
// duplication of the envelopes. If the bulk sizes are small then the communication cost will be
// high relative to the "computation cost". Optimizations to "merge adjacent blocks" take place in
// the manipulation of the octree grid structure. In the grid statistics, we do not consider the
// optimization to make the overlap blocks smaller either.
template <typename T>
class ComputeOctreeGridStatistics {
public:
    ComputeOctreeGridStatistics(
        int minLeafLevel_, int maxLeafLevel_, plint numProcesses_, plint numCellsPerBlock_);
    plint getNumLeafBlocksOnLevel(int level) const;
    plint getNumOverlapBlocksOnLevel(int level) const;
    plint getNumBlocksOnLevelOnProcess(int level, plint processId) const;
    plint getNumBlocksOnProcess(plint processId) const;
    T getRelativeLoadOnLevelOnProcess(int level, plint processId) const;
    T getRelativeLoadOnProcess(plint processId) const;
    T getStandardDeviationOfRelativeLoadOnLevel(int level) const;
    T getStandardDeviationOfRelativeLoad() const;
    T getRangeOfRelativeLoadOnLevel(int level) const;
    T getRangeOfRelativeLoad() const;
    void operator()(OctreeNode<OctreeData<T> > *node);

private:
    int minLeafLevel, maxLeafLevel;
    plint numProcesses, numCellsPerBlock;
    std::vector<plint> numLeafBlocksOnLevel;
    std::vector<plint> numOverlapBlocksOnLevel;
    std::vector<std::vector<plint> > numBlocksOnLevelOnProcess;
    std::vector<std::vector<plint> > loadOnLevelOnProcess;
    std::vector<plint> totalLoadOnLevel;
    plint totalLoad;
};

template <typename T>
class ConstructOctreeGridStructure {
public:
    ConstructOctreeGridStructure(OctreeGridStructure &octreeGridStructure_, int minLeafLevel_);
    void operator()(OctreeNode<OctreeData<T> > *node);

private:
    OctreeGridStructure &octreeGridStructure;
    int minLeafLevel;
};

template <typename T>
class OctreeGridGenerator {
public:
    OctreeGridGenerator(std::string xmlFileName);
    OctreeGridGenerator(
        Cuboid<T> const &fullDomain_, std::string gridDensityFunctionFile_, int minLeafLevel_,
        int maxLeafLevel_, plint n_, plint numProcesses_, bool xPeriodic_ = false,
        bool yPeriodic_ = false, bool zPeriodic_ = false, T gridDensityScaleFactor_ = (T)1,
        bool useSamples_ = false, plint numSamples_ = -1, plint maxIter_ = 100,
        bool removeBlocks_ = true, bool fineToCoarse_ = true, int numLevelsToGroupBlocks_ = 0,
        int numLevelsToGroupOverlaps_ = -1, bool strongGrouping_ = false,
        std::string outDir_ = "./", bool verbose_ = true, bool stlOutput_ = true,
        std::string stlBaseName_ = "octree");

    OctreeGridStructure generateOctreeGridStructure();

    Cuboid<T> const &getFullDomain() const;
    int getMinLevel() const;
    int getMaxLevel() const;
    plint getNumProcesses() const;
    plint getBlockNx() const;
    plint getBlockNy() const;
    plint getBlockNz() const;
    T getDxFinestLevel() const;
    std::string getOutDir() const;

    void setOutputAndErrorStreams(Parallel_ostream &outStream, Parallel_ostream &errStream);

private:
    // Geometry and grid density function.

    Cuboid<T> fullDomain;
    std::string gridDensityFunctionFile;
    bool xPeriodic, yPeriodic, zPeriodic;
    T gridDensityScaleFactor;

    // Octree construction parameters.

    int minLeafLevel;
    int maxLeafLevel;
    bool useSamples;
    plint numSamples;
    plint maxIter;
    bool removeBlocks;

    // Block size and number of processes.

    plint n;
    plint numProcesses;

    // Process assignment parameters.

    // Assign processIds from the finest to the coarsest level or the other way around.
    // We choose that both overlaps and allocatable leaves are grouped in the same
    // direction (fine -> coarse, or coarse -> fine).
    bool fineToCoarse;

    // The octree nodes can be grouped or not. By grouping nodes we mean
    // that we assign adjacent nodes to the same process (they have
    // the same processId).
    // If numLevelsToGroup* = 0, this means that no grouping takes
    // place.
    // If numLevelsToGroup* = 1, this means that the adjacent
    // nodes are grouped "per level", so in groups of 8 nodes at
    // most.
    // If numLevelsToGroup* = 2, this means that the groups have
    // 64 nodes at most.
    // In general, if numLevelsToGroup* = l, then the respective
    // groups have at most 8^l nodes.
    int numLevelsToGroupBlocks;    // For allocatable leaves.
    int numLevelsToGroupOverlaps;  // For overlaps.

    // If strongGrouping is true, then the overlap blocks are grouped with
    // all their leaf children.
    // If strongGrouping is false, then the overlap blocks are grouped
    // only with their leaf children to which they are coupled.
    bool strongGrouping;

    // Output parameters.

    std::string outDir;
    bool verbose;
    bool stlOutput;
    std::string stlBaseName;

    // Derived parameters (not user-provided).

    plint nx, ny, nz;
    T dxFinest;

    // Output and error parallel streams.

    Parallel_ostream *outS;
    Parallel_ostream *errS;
};

}  // namespace plb

#endif  // OCTREE_GRID_GENERATOR_H
