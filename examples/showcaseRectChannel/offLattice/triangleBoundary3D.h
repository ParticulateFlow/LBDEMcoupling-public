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

#ifndef TRIANGLE_BOUNDARY_3D_H
#define TRIANGLE_BOUNDARY_3D_H

#include <stack>

#include "atomicBlock/atomicContainerBlock3D.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "multiBlock/multiContainerBlock3D.h"
#include "multiBlock/multiDataField3D.h"
#include "multiBlock/redistribution3D.h"
#include "offLattice/boundaryShapes3D.h"
#include "offLattice/offLatticeBoundaryProfiles3D.h"
#include "offLattice/triangleSet.h"
#include "offLattice/triangularSurfaceMesh.h"
#include "particles/multiParticleField3D.h"

namespace plb {

template <typename T>
class DEFscaledMesh {
public:
    DEFscaledMesh(TriangleSet<T> const &triangleSet_);
    DEFscaledMesh(
        TriangleSet<T> const &triangleSet_, plint resolution_, plint referenceDirection_,
        plint margin_, plint extraLayer);
    DEFscaledMesh(
        TriangleSet<T> const &triangleSet_, plint resolution_, plint referenceDirection_,
        plint margin_, Dot3D location);
    DEFscaledMesh(DEFscaledMesh<T> const &rhs);
    ~DEFscaledMesh();
    DEFscaledMesh<T> &operator=(DEFscaledMesh<T> const &rhs);
    void swap(DEFscaledMesh<T> &rhs);

public:  // Mesh usage interface.
    /// Get a reference to the currently active mesh.
    TriangularSurfaceMesh<T> &getMesh();
    /// Get a const reference to the currently active mesh.
    TriangularSurfaceMesh<T> const &getMesh() const;
    plint getMargin() const;
    void setMargin(plint margin_)
    {
        margin = margin_;
    }
    std::vector<Array<T, 3> > const &getVertexList() const
    {
        return vertexList;
    }
    std::vector<plint> const &getEmanatingEdgeList() const
    {
        return emanatingEdgeList;
    }
    std::vector<Edge> const &getEdgeList() const
    {
        return edgeList;
    }
    Array<T, 3> getPhysicalLocation() const
    {
        return physicalLocation;
    }
    T getDx() const
    {
        return dx;
    }
    void setPhysicalLocation(Array<T, 3> physicalLocation_)
    {
        physicalLocation = physicalLocation_;
    }
    void setDx(T dx_)
    {
        dx = dx_;
    }

private:
    void initialize(
        TriangleSet<T> const &triangleSet_, plint resolution_, plint referenceDirection_,
        Dot3D location);

private:
    std::vector<Array<T, 3> > vertexList;
    /// Each vertex has exactly one emanating edge. This is a structural
    ///   information.
    std::vector<plint> emanatingEdgeList;
    /// Edges are a structural information.
    std::vector<Edge> edgeList;
    TriangularSurfaceMesh<T> *mesh;
    plint margin;
    Array<T, 3> physicalLocation;
    T dx;
};

template <typename T>
struct VertexProperty3D {
    virtual ~VertexProperty3D() { }
    // Wall portions which should not follow the laws of elasticity are rigid.
    virtual bool isRigid() const = 0;
    // Inlet/outlet nodes are for example not part of the wall.
    virtual bool isWall() const = 0;
    virtual VertexProperty3D<T> *clone() const = 0;
};

template <typename T>
struct RigidWallProperty3D : public VertexProperty3D<T> {
    virtual bool isRigid() const
    {
        return true;
    }
    virtual bool isWall() const
    {
        return true;
    }
    virtual RigidWallProperty3D<T> *clone() const
    {
        return new RigidWallProperty3D<T>(*this);
    }
};

template <typename T>
struct InletOutletProperty3D : public VertexProperty3D<T> {
    virtual bool isRigid() const
    {
        return true;
    }
    virtual bool isWall() const
    {
        return false;
    }
    virtual InletOutletProperty3D<T> *clone() const
    {
        return new InletOutletProperty3D<T>(*this);
    }
};

template <typename T>
bool isRigid(VertexProperty3D<T> const *property)
{
    if (property) {
        return property->isRigid();
    } else {
        return false;
    }
}

template <typename T>
bool isWall(VertexProperty3D<T> const *property)
{
    if (property) {
        return property->isWall();
    } else {
        return true;
    }
}

template <typename T>
class TriangleBoundary3D;

template <typename T, class SurfaceData>
class BoundaryProfiles3D {
public:
    BoundaryProfiles3D();
    ~BoundaryProfiles3D();
    BoundaryProfiles3D(BoundaryProfiles3D<T, SurfaceData> const &rhs);
    BoundaryProfiles3D<T, SurfaceData> &operator=(BoundaryProfiles3D<T, SurfaceData> const &rhs);
    void swap(BoundaryProfiles3D<T, SurfaceData> &rhs);
    void defineProfile(plint tag, BoundaryProfile3D<T, SurfaceData> *profile);
    void resetProfiles(std::map<plint, BoundaryProfile3D<T, SurfaceData> *> profiles_);
    void defineInletOutletTags(TriangleBoundary3D<T> const &boundary, plint sortDirection);
    void setWallProfile(BoundaryProfile3D<T, SurfaceData> *wallProfile_);
    void setInletOutlet(std::vector<BoundaryProfile3D<T, SurfaceData> *> inletOutlets);
    void setInletOutlet(
        BoundaryProfile3D<T, SurfaceData> *profile1, BoundaryProfile3D<T, SurfaceData> *profile2);
    void setInletOutlet(
        BoundaryProfile3D<T, SurfaceData> *profile1, BoundaryProfile3D<T, SurfaceData> *profile2,
        BoundaryProfile3D<T, SurfaceData> *profile3);
    void setInletOutlet(
        BoundaryProfile3D<T, SurfaceData> *profile1, BoundaryProfile3D<T, SurfaceData> *profile2,
        BoundaryProfile3D<T, SurfaceData> *profile3, BoundaryProfile3D<T, SurfaceData> *profile4);
    void adjustInletOutlet(TriangleBoundary3D<T> const &boundary, plint sortDirection);
    BoundaryProfile3D<T, SurfaceData> const &getProfile(
        TriangleBoundary3D<T> const &boundary, plint iTriangle) const;

private:
    void replaceProfile(plint id, BoundaryProfile3D<T, SurfaceData> *newProfile);
    void clearProfiles();

private:
    BoundaryProfile3D<T, SurfaceData> *wallProfile;
    std::map<plint, BoundaryProfile3D<T, SurfaceData> *> profiles;
    std::vector<plint> inletOutletIds;
    std::vector<Array<T, 3> > lidNormal;
    std::vector<Array<T, 3> > lidCenter;
    std::vector<T> lidRadius;
};

template <typename T>
class TriangleBoundary3D {
public:
    template <typename TMesh>
    TriangleBoundary3D(DEFscaledMesh<TMesh> const &defMesh, bool automaticCloseHoles = true);
    TriangleBoundary3D(DEFscaledMesh<T> const &defMesh, bool automaticCloseHoles = true);
    ~TriangleBoundary3D();
    TriangleBoundary3D(TriangleBoundary3D<T> const &rhs);
    TriangleBoundary3D<T> &operator=(TriangleBoundary3D<T> const &rhs);
    void swap(TriangleBoundary3D<T> &rhs);

public:  // Mesh usage interface.
    /// Select the mesh which is subsequently being referred to in calls to
    ///   the methods of this class.
    TriangleBoundary3D<T> const &select(plint whichTopology, plint whichVertices) const;
    /// Select the mesh which is subsequently being referred to in calls to
    ///   the methods of this class. Save the previous selection, which
    ///   can be recovered through a call to popSelect.
    TriangleBoundary3D<T> const &pushSelect(plint whichTopology, plint whichVertices) const;
    /// Recover the previous selection of the mesh, stored through a
    ///   call to pushSelect.
    TriangleBoundary3D<T> const &popSelect() const;
    void getSelection(plint &whichTopology, plint &whichVertices) const;
    /// Get a reference to the currently active mesh.
    TriangularSurfaceMesh<T> &getMesh();
    /// Get a const reference to the currently active mesh.
    TriangularSurfaceMesh<T> const &getMesh() const;
    /// Get the material property (for example the elasticity constants)
    ///   implemented on a given vertex (the answer is independent of the
    ///   currently active mesh).
    VertexProperty3D<T> const *getVertexProperty(plint iVertex) const;
    /// Get intersection between a line segment (fromPoint,fromPoint+direction)
    ///   and a given triangle in the currently active mesh; return true incase
    ///   of success.
    bool intersectSegment(
        plint iTriangle, AtomicBlock3D *boundaryArg, Array<T, 3> const &fromPoint,
        Array<T, 3> const &direction, Array<T, 3> &locatedPoint, T &distance,
        Array<T, 3> &wallNormal) const;
    /// Given a point p on the surface of the shape, determine its "continuous normal".
    ///   If the shape is for example piecewise linear, the normal is adjusted to vary
    ///   continuously over the surface.
    Array<T, 3> computeContinuousNormal(
        Array<T, 3> const &p, plint iTriangle, bool isAreaWeighted = false) const;
    /// Create a new set of vertices, and an associated open and closed mesh.
    void cloneVertexSet(plint whichVertexSet);
    /// Coordinates of the lower-left corner in physical units.
    Array<T, 3> getPhysicalLocation() const
    {
        return physicalLocation;
    }
    /// Size of a grid spacing.
    T getDx() const
    {
        return dx;
    }

public:  // Mesh preparation interface.
    /// Get a list of all inlets and outlets (no specific sorting order).
    std::vector<Lid> const &getInletOutlet() const;
    /// Get a list of all inlets and outlets, sorted along a given direction.
    /** This information can be used as a hint to select the boundary condition
     *    associated to each inlet/outlet. Access the variable lid.baryCenter
     *    to know the location of a given inlet/outlet.
     **/
    std::vector<Lid> getInletOutlet(plint sortDirection) const;
    template <typename DomainFunctional>
    plint tagDomain(DomainFunctional functional);
    template <typename DomainFunctional>
    plint tagDomain(
        DomainFunctional functional, Array<T, 3> normal, T angleTolerance, plint previousTag = -1);
    /// Tag all lids whose barycenters are inside the given cuboid and return the integer tag.
    plint tagLids(Cuboid<T> const &c);
    std::vector<plint> getInletOutletIds(plint sortDirection) const;
    void getLidProperties(
        plint sortDirection, std::vector<Array<T, 3> > &normal, std::vector<Array<T, 3> > &center,
        std::vector<T> &radius) const;
    /// Define the material property to be used on all vertices contained in
    ///   the domain specified by the domain functional.
    /** Returns the tag which was assigned to the corresponding vertices. **/
    template <typename DomainFunctional>
    plint setVertexProperty(VertexProperty3D<T> const &property, DomainFunctional functional);
    plint getMargin() const;
    /// Return the tag (id of boundary-portion with specific boundary condition)
    ///   of a triangle.
    plint getTag(plint iTriangle) const;
    std::vector<plint> const &getTriangleTags() const
    {
        return triangleTagList;
    }
    std::vector<plint> const &getVertexTags() const
    {
        return vertexTagList;
    }

private:
    plint currentMesh() const;
    void defineMeshes();
    /// Detect holes, close them and register them as potential inlets/outlets.
    ///   They are default initialized to no-slip. Use setInletOutlet to assign
    ///   different conditions.
    void closeHoles();
    void assignLidVertexProperty();

private:
    /// Assign a new tag to all triangles corresponding to one of the provided inlets/outlets.
    ///   The tag is taken in increasing integer value according to a sorting
    ///   or the inlets/outlets along the given space direction.
    void tagInletOutlet(std::vector<Lid> &newLids);
    /// There may exist more than one set of vertices, for example in
    ///   case of a moving wall which has current vertex positions and
    ///   equilibrium vertex positions.
    std::vector<std::vector<Array<T, 3> > > vertexLists;
    /// Each vertex has exactly one emanating edge. This is a structural
    ///   information which is identical for all sets of vertices.
    std::vector<plint> emanatingEdgeLists[2];
    /// Edges are a structural information which is shared by all sets
    ///   of vertices.
    std::vector<Edge> edgeLists[2];
    /// For each set of vertices there is a mesh, with a reference to the
    ///   given vertices (individual to each mesh) and to the
    ///   emanatingEdgeList and edgeList (same for all meshes).
    std::vector<TriangularSurfaceMesh<T> > meshes;
    /// The triangle type is an indirect index which links to the boundary
    ///   condition implemented by each triangle and defined in boundaryProfiles.
    std::vector<plint> triangleTagList;
    plint currentTagNum;
    /// The vertex type is an indirect index which links to generic material
    ///   properties implemented at that vertex and defined in vertexProperties.
    std::vector<plint> vertexTagList;
    /// Vertex properties, indexed by the vertex type in vertexTagList.
    std::vector<VertexProperty3D<T> *> vertexProperties;
    /// Inlets and outlets, saved as a collection of triangles.
    std::vector<Lid> lids;
    plint margin;
    Array<T, 3> physicalLocation;
    T dx;
    mutable std::stack<plint> topology;
    mutable std::stack<plint> vertexSet;
};

template <typename T, class SurfaceData>
class TriangleFlowShape3D : public BoundaryShape3D<T, SurfaceData> {
public:
    TriangleFlowShape3D(
        TriangleBoundary3D<T> const &boundary_,
        BoundaryProfiles3D<T, SurfaceData> const &profiles_);
    virtual bool isInside(Dot3D const &location) const;
    virtual bool isOutside(Dot3D const &location) const;
    virtual bool pointOnSurface(
        Array<T, 3> const &fromPoint, Array<T, 3> const &direction, Array<T, 3> &locatedPoint,
        T &distance, Array<T, 3> &wallNormal, SurfaceData &surfaceData, OffBoundary::Type &bdType,
        plint &id) const;
    virtual Array<T, 3> computeContinuousNormal(
        Array<T, 3> const &p, plint id, bool isAreaWeighted = false) const;
    virtual bool intersectsSurface(Array<T, 3> const &p1, Array<T, 3> const &p2, plint &id) const;
    virtual plint getTag(plint id) const;
    virtual bool distanceToSurface(Array<T, 3> const &point, T &distance, bool &isBehind) const;
    virtual TriangleFlowShape3D<T, SurfaceData> *clone() const;
    /// Use this clone function to provide the meshed data to this object.
    /** The arguments are:
     *  0: The voxel flags (ScalarField3D<int>),
     *  1: The hash container (AtomicContainerBlock3D),
     *  2: The boundary argument: an additional argument needed by BoundaryProfiles
     *     in order to compute the boundary condition. In dynamic walls this is for
     *     example often a MultiParticleField, used to determined the wall velocity.
     **/
    virtual TriangleFlowShape3D<T, SurfaceData> *clone(std::vector<AtomicBlock3D *> args) const;

private:
    TriangleBoundary3D<T> const &boundary;
    BoundaryProfiles3D<T, SurfaceData> const &profiles;
    /// Data from previous voxelization.
    ScalarField3D<int> *voxelFlags;
    /// Needed for fast access to the mesh.
    AtomicContainerBlock3D *hashContainer;
    AtomicBlock3D *boundaryArg;
};

template <typename T>
class VoxelizedDomain3D {
public:
    VoxelizedDomain3D(
        TriangleBoundary3D<T> const &boundary_, int flowType_, plint extraLayer_,
        plint borderWidth_, plint envelopeWidth_, plint blockSize_, plint gridLevel_ = 0,
        bool dynamicMesh_ = false);
    VoxelizedDomain3D(
        TriangleBoundary3D<T> const &boundary_, int flowType_, Box3D const &boundingBox,
        plint borderWidth_, plint envelopeWidth_, plint blockSize_, plint gridLevel_ = 0,
        bool dynamicMesh_ = false);
    // For faster results, the "seed" should contain at least one of the domain corners.
    VoxelizedDomain3D(
        TriangleBoundary3D<T> const &boundary_, int flowType_, Box3D const &boundingBox,
        plint borderWidth_, plint envelopeWidth_, plint blockSize_, Box3D const &seed,
        plint gridLevel_ = 0, bool dynamicMesh_ = false);
    VoxelizedDomain3D(VoxelizedDomain3D<T> const &rhs);
    ~VoxelizedDomain3D();
    MultiScalarField3D<int> &getVoxelMatrix();
    MultiScalarField3D<int> const &getVoxelMatrix() const;
    MultiContainerBlock3D &getTriangleHash();
    MultiBlockManagement3D const &getMultiBlockManagement() const;
    template <class ParticleFieldT>
    void adjustVoxelization(MultiParticleField3D<ParticleFieldT> &particles, bool dynamicMesh);
    void reparallelize(MultiBlockRedistribute3D const &redistribute);
    void reparallelize(MultiBlockManagement3D const &newManagement);
    TriangleBoundary3D<T> const &getBoundary() const
    {
        return boundary;
    }
    int getFlowType() const
    {
        return flowType;
    }

private:
    VoxelizedDomain3D<T> &operator=(VoxelizedDomain3D<T> const &rhs) { }
    void createSparseVoxelMatrix(
        MultiScalarField3D<int> &fullVoxelMatrix, plint blockSize_, plint envelopeWidth_);
    void computeSparseVoxelMatrix(
        MultiScalarField3D<int> &fullVoxelMatrix, plint blockSize, plint envelopeWidth);
    void extendEnvelopeWidth(MultiScalarField3D<int> &fullVoxelMatrix, plint envelopeWidth);
    void createTriangleHash();
    template <class ParticleFieldT>
    void reCreateTriangleHash(MultiParticleField3D<ParticleFieldT> &particles);
    void computeOuterMask();

private:
    int flowType;
    plint borderWidth;
    TriangleBoundary3D<T> const &boundary;
    MultiScalarField3D<int> *voxelMatrix;
    MultiContainerBlock3D *triangleHash;
};

template <typename T>
void addLayer(MultiScalarField3D<T> &matrix, Box3D const &domain, T previousLayer);

template <typename T>
class AddLayerFunctional3D : public BoxProcessingFunctional3D_S<T> {
public:
    AddLayerFunctional3D(T previousLayer_);
    virtual void process(Box3D domain, ScalarField3D<T> &voxels);
    virtual AddLayerFunctional3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    T previousLayer;
};

}  // namespace plb

#endif  // TRIANGLE_BOUNDARY_3D_H
