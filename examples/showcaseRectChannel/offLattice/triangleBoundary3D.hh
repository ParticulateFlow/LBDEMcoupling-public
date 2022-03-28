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

#ifndef TRIANGLE_BOUNDARY_3D_HH
#define TRIANGLE_BOUNDARY_3D_HH

#include <cmath>
#include <limits>

#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "dataProcessors/dataAnalysisWrapper3D.h"
#include "dataProcessors/dataInitializerWrapper3D.h"
#include "multiBlock/defaultMultiBlockPolicy3D.h"
#include "multiBlock/multiBlockGenerator3D.h"
#include "multiBlock/nonLocalTransfer3D.h"
#include "offLattice/makeSparse3D.h"
#include "offLattice/offLatticeBoundaryProfiles3D.h"
#include "offLattice/triangleBoundary3D.h"
#include "offLattice/triangleToDef.h"
#include "offLattice/triangularSurfaceMesh.h"
#include "offLattice/voxelizer.h"

namespace plb {

/******** class BoundaryProfiles3D *****************************************/

template <typename T, class SurfaceData>
BoundaryProfiles3D<T, SurfaceData>::BoundaryProfiles3D()
{
    wallProfile = generateDefaultWallProfile3D<T, SurfaceData>();
    replaceProfile(0, wallProfile->clone());
}

template <typename T, class SurfaceData>
BoundaryProfiles3D<T, SurfaceData>::~BoundaryProfiles3D()
{
    clearProfiles();
    delete wallProfile;
}

template <typename T, class SurfaceData>
BoundaryProfiles3D<T, SurfaceData>::BoundaryProfiles3D(
    BoundaryProfiles3D<T, SurfaceData> const &rhs)
{
    wallProfile = rhs.wallProfile->clone();
    typename std::map<plint, BoundaryProfile3D<T, SurfaceData> *>::const_iterator it =
        rhs.profiles.begin();
    for (; it != rhs.profiles.end(); ++it) {
        profiles.insert(
            std::pair<plint, BoundaryProfile3D<T, SurfaceData> *>(it->first, it->second->clone()));
    }
    inletOutletIds = rhs.inletOutletIds;
    lidNormal = rhs.lidNormal;
    lidCenter = rhs.lidCenter;
    lidRadius = rhs.lidRadius;
}

template <typename T, class SurfaceData>
BoundaryProfiles3D<T, SurfaceData> &BoundaryProfiles3D<T, SurfaceData>::operator=(
    BoundaryProfiles3D<T, SurfaceData> const &rhs)
{
    BoundaryProfiles3D<T, SurfaceData>(rhs).swap(*this);
    return *this;
}

template <typename T, class SurfaceData>
void BoundaryProfiles3D<T, SurfaceData>::swap(BoundaryProfiles3D<T, SurfaceData> &rhs)
{
    std::swap(wallProfile, rhs.wallProfile);
    profiles.swap(rhs.profiles);
    inletOutletIds.swap(rhs.inletOutletIds);
    lidNormal.swap(rhs.lidNormal);
    lidCenter.swap(rhs.lidCenter);
    lidRadius.swap(rhs.lidRadius);
}

template <typename T, class SurfaceData>
void BoundaryProfiles3D<T, SurfaceData>::setWallProfile(
    BoundaryProfile3D<T, SurfaceData> *wallProfile_)
{
    PLB_PRECONDITION(wallProfile_);
    delete wallProfile;
    wallProfile = wallProfile_;
    replaceProfile(0, wallProfile->clone());
}

template <typename T, class SurfaceData>
void BoundaryProfiles3D<T, SurfaceData>::resetProfiles(
    std::map<plint, BoundaryProfile3D<T, SurfaceData> *> profiles_)
{
    PLB_ASSERT(!profiles_.empty());
    clearProfiles();
    profiles = profiles_;
}

template <typename T, class SurfaceData>
void BoundaryProfiles3D<T, SurfaceData>::defineInletOutletTags(
    TriangleBoundary3D<T> const &boundary, plint sortDirection)
{
    inletOutletIds = boundary.getInletOutletIds(sortDirection);
    boundary.getLidProperties(sortDirection, lidNormal, lidCenter, lidRadius);
}

template <typename T, class SurfaceData>
void BoundaryProfiles3D<T, SurfaceData>::adjustInletOutlet(
    TriangleBoundary3D<T> const &boundary, plint sortDirection)
{
    boundary.getLidProperties(sortDirection, lidNormal, lidCenter, lidRadius);
    for (pluint iProfile = 0; iProfile < inletOutletIds.size(); ++iProfile) {
        plint id = inletOutletIds[iProfile];
        typename std::map<plint, BoundaryProfile3D<T, SurfaceData> *>::iterator it =
            profiles.find(id);
        PLB_ASSERT(it != profiles.end());
        it->second->setNormal(lidNormal[iProfile]);
        it->second->defineCircularShape(lidCenter[iProfile], lidRadius[iProfile]);
    }
}

template <typename T, class SurfaceData>
void BoundaryProfiles3D<T, SurfaceData>::setInletOutlet(
    std::vector<BoundaryProfile3D<T, SurfaceData> *> inletOutlets)
{
    PLB_ASSERT(inletOutletIds.size() == inletOutlets.size());
    std::map<plint, BoundaryProfile3D<T, SurfaceData> *> newProfiles;
    newProfiles[0] = wallProfile->clone();
    for (pluint iProfile = 0; iProfile < inletOutlets.size(); ++iProfile) {
        plint id = inletOutletIds[iProfile];
#ifdef PLB_DEBUG
        typename std::map<plint, BoundaryProfile3D<T, SurfaceData> *>::const_iterator it =
            newProfiles.find(id);
        PLB_ASSERT(it == newProfiles.end());
#endif
        newProfiles[id] = inletOutlets[iProfile];
        newProfiles[id]->setNormal(lidNormal[iProfile]);
        newProfiles[id]->defineCircularShape(lidCenter[iProfile], lidRadius[iProfile]);
    }
    resetProfiles(newProfiles);
}

template <typename T, class SurfaceData>
void BoundaryProfiles3D<T, SurfaceData>::setInletOutlet(
    BoundaryProfile3D<T, SurfaceData> *profile1, BoundaryProfile3D<T, SurfaceData> *profile2)
{
    std::vector<BoundaryProfile3D<T, SurfaceData> *> inletOutlets(2);
    inletOutlets[0] = profile1;
    inletOutlets[1] = profile2;
    setInletOutlet(inletOutlets);
}

template <typename T, class SurfaceData>
void BoundaryProfiles3D<T, SurfaceData>::setInletOutlet(
    BoundaryProfile3D<T, SurfaceData> *profile1, BoundaryProfile3D<T, SurfaceData> *profile2,
    BoundaryProfile3D<T, SurfaceData> *profile3)
{
    std::vector<BoundaryProfile3D<T, SurfaceData> *> inletOutlets(3);
    inletOutlets[0] = profile1;
    inletOutlets[1] = profile2;
    inletOutlets[2] = profile3;
    setInletOutlet(inletOutlets);
}

template <typename T, class SurfaceData>
void BoundaryProfiles3D<T, SurfaceData>::setInletOutlet(
    BoundaryProfile3D<T, SurfaceData> *profile1, BoundaryProfile3D<T, SurfaceData> *profile2,
    BoundaryProfile3D<T, SurfaceData> *profile3, BoundaryProfile3D<T, SurfaceData> *profile4)
{
    std::vector<BoundaryProfile3D<T, SurfaceData> *> inletOutlets(4);
    inletOutlets[0] = profile1;
    inletOutlets[1] = profile2;
    inletOutlets[2] = profile3;
    inletOutlets[3] = profile4;
    setInletOutlet(inletOutlets);
}

template <typename T, class SurfaceData>
void BoundaryProfiles3D<T, SurfaceData>::defineProfile(
    plint tag, BoundaryProfile3D<T, SurfaceData> *profile)
{
    typename std::map<plint, BoundaryProfile3D<T, SurfaceData> *>::iterator it = profiles.find(tag);
    if (it != profiles.end()) {
        delete it->second;
    }
    profiles[tag] = profile;
}

template <typename T, class SurfaceData>
BoundaryProfile3D<T, SurfaceData> const &BoundaryProfiles3D<T, SurfaceData>::getProfile(
    TriangleBoundary3D<T> const &boundary, plint iTriangle) const
{
    PLB_ASSERT(iTriangle < (plint)boundary.getTriangleTags().size());
    plint triangleTag = boundary.getTriangleTags()[iTriangle];
    typename std::map<plint, BoundaryProfile3D<T, SurfaceData> *>::const_iterator it =
        profiles.find(triangleTag);
    if (it == profiles.end()) {
        PLB_ASSERT(wallProfile);
        return *wallProfile;
    } else {
        PLB_ASSERT(it->second);
        return *(it->second);
    }
}

template <typename T, class SurfaceData>
void BoundaryProfiles3D<T, SurfaceData>::replaceProfile(
    plint id, BoundaryProfile3D<T, SurfaceData> *newProfile)
{
    typename std::map<plint, BoundaryProfile3D<T, SurfaceData> *>::iterator it = profiles.find(id);
    if (it != profiles.end()) {
        delete it->second;
    }
    profiles[id] = newProfile;
}

template <typename T, class SurfaceData>
void BoundaryProfiles3D<T, SurfaceData>::clearProfiles()
{
    typename std::map<plint, BoundaryProfile3D<T, SurfaceData> *>::iterator it = profiles.begin();
    for (; it != profiles.end(); ++it) {
        delete it->second;
    }
    profiles.clear();
}

/******** class DEFscaledMesh *****************************************/

template <typename T>
DEFscaledMesh<T>::DEFscaledMesh(TriangleSet<T> const &triangleSet_) : margin(0)
{
    initialize(triangleSet_, 0, 0, Dot3D());
}

template <typename T>
DEFscaledMesh<T>::DEFscaledMesh(
    TriangleSet<T> const &triangleSet_, plint resolution_, plint referenceDirection_, plint margin_,
    Dot3D location) :
    margin(margin_)
{
    initialize(triangleSet_, resolution_, referenceDirection_, location);
}

template <typename T>
DEFscaledMesh<T>::DEFscaledMesh(
    TriangleSet<T> const &triangleSet_, plint resolution_, plint referenceDirection_, plint margin_,
    plint extraLayer) :
    margin(margin_)
{
    plint layer = margin + extraLayer;
    initialize(triangleSet_, resolution_, referenceDirection_, Dot3D(layer, layer, layer));
}

template <typename T>
DEFscaledMesh<T>::~DEFscaledMesh()
{
    delete mesh;
}

template <typename T>
void DEFscaledMesh<T>::initialize(
    TriangleSet<T> const &triangleSet_, plint resolution_, plint referenceDirection_,
    Dot3D location)
{
    T eps = triangleSet_.getEpsilon();

    constructSurfaceMesh<T>(
        triangleSet_.getTriangles(), vertexList, emanatingEdgeList, edgeList, eps);
    mesh = new TriangularSurfaceMesh<T>(vertexList, emanatingEdgeList, edgeList);

    if (resolution_ != 0) {
        // Convert the mesh to lattice units.
        toLatticeUnits<T>(getMesh(), resolution_, referenceDirection_, physicalLocation, dx);

        Array<T, 3> luOffset(location.x, location.y, location.z);
        getMesh().translate(luOffset);
        physicalLocation -= luOffset * dx;
    }
}

template <typename T>
DEFscaledMesh<T>::DEFscaledMesh(DEFscaledMesh<T> const &rhs) :
    vertexList(rhs.vertexList),
    emanatingEdgeList(rhs.emanatingEdgeList),
    edgeList(rhs.edgeList),
    margin(rhs.margin),
    physicalLocation(rhs.physicalLocation),
    dx(rhs.dx)
{
    mesh = new TriangularSurfaceMesh<T>(vertexList, emanatingEdgeList, edgeList);
}

template <typename T>
DEFscaledMesh<T> &DEFscaledMesh<T>::operator=(DEFscaledMesh<T> const &rhs)
{
    DEFscaledMesh<T>(rhs).swap(*this);
    return *this;
}

template <typename T>
void DEFscaledMesh<T>::swap(DEFscaledMesh<T> &rhs)
{
    vertexList.swap(rhs.vertexList);
    emanatingEdgeList.swap(rhs.emanatingEdgeList);
    edgeList.swap(rhs.edgeList);
    std::swap(margin, rhs.margin);
    std::swap(physicalLocation, rhs.physicalLocation);
    std::swap(dx, rhs.dx);
    delete mesh;
    mesh = new TriangularSurfaceMesh<T>(vertexList, emanatingEdgeList, edgeList);
    delete rhs.mesh;
    rhs.mesh = new TriangularSurfaceMesh<T>(rhs.vertexList, rhs.emanatingEdgeList, rhs.edgeList);
}

template <typename T>
TriangularSurfaceMesh<T> const &DEFscaledMesh<T>::getMesh() const
{
    return *mesh;
}

template <typename T>
TriangularSurfaceMesh<T> &DEFscaledMesh<T>::getMesh()
{
    return *mesh;
}

template <typename T>
plint DEFscaledMesh<T>::getMargin() const
{
    return margin;
}

/******** class TriangleBoundary3D *****************************************/

template <typename T>
TriangleBoundary3D<T>::TriangleBoundary3D(
    DEFscaledMesh<T> const &defMesh, bool automaticCloseHoles) :
    currentTagNum(0),
    margin(defMesh.getMargin()),
    physicalLocation(defMesh.getPhysicalLocation()),
    dx(defMesh.getDx())
{
    topology.push(1);   // By default, closed mesh.
    vertexSet.push(0);  // By default, static mesh.

    vertexLists.reserve(3);  // Vertex lists are expensive to copy. Better
                             //   pre-allocate a slot for three of them.
    vertexLists.push_back(defMesh.getVertexList());
    emanatingEdgeLists[0] = defMesh.getEmanatingEdgeList();
    edgeLists[0] = defMesh.getEdgeList();

    emanatingEdgeLists[1] = emanatingEdgeLists[0];
    edgeLists[1] = edgeLists[0];

    meshes.push_back(TriangularSurfaceMesh<T>(vertexLists[0], emanatingEdgeLists[0], edgeLists[0]));
    meshes.push_back(TriangularSurfaceMesh<T>(vertexLists[0], emanatingEdgeLists[1], edgeLists[1]));

    // Prepare the vector "triangle type", which later on will inform on
    //   the type of boundary condition implemented by a given triangle.
    //   The default, 0, stands for no-slip.
    triangleTagList.resize(meshes[1].getNumTriangles());
    std::fill(triangleTagList.begin(), triangleTagList.end(), 0);

    if (automaticCloseHoles) {
        closeHoles();
    }
}

template <typename T>
template <typename TMesh>
TriangleBoundary3D<T>::TriangleBoundary3D(
    DEFscaledMesh<TMesh> const &defMesh, bool automaticCloseHoles) :
    currentTagNum(0),
    margin(defMesh.getMargin()),
    physicalLocation(defMesh.getPhysicalLocation()),
    dx(defMesh.getDx())
{
    topology.push(1);   // By default, closed mesh.
    vertexSet.push(0);  // By default, static mesh.

    vertexLists.reserve(3);  // Vertex lists are expensive to copy. Better
                             //   pre-allocate a slot for three of them.
    std::vector<Array<TMesh, 3> > const &vertexList = defMesh.getVertexList();
    std::vector<Array<T, 3> > newVertexList(vertexList.size());
    for (pluint i = 0; i < vertexList.size(); ++i) {
        newVertexList[i] = Array<T, 3>(vertexList[i]);
    }
    vertexLists.push_back(newVertexList);
    emanatingEdgeLists[0] = defMesh.getEmanatingEdgeList();
    edgeLists[0] = defMesh.getEdgeList();

    emanatingEdgeLists[1] = emanatingEdgeLists[0];
    edgeLists[1] = edgeLists[0];

    meshes.push_back(TriangularSurfaceMesh<T>(vertexLists[0], emanatingEdgeLists[0], edgeLists[0]));
    meshes.push_back(TriangularSurfaceMesh<T>(vertexLists[0], emanatingEdgeLists[1], edgeLists[1]));

    // Prepare the vector "triangle type", which later on will inform on
    //   the type of boundary condition implemented by a given triangle.
    //   The default, 0, stands for no-slip.
    triangleTagList.resize(meshes[1].getNumTriangles());
    std::fill(triangleTagList.begin(), triangleTagList.end(), 0);

    if (automaticCloseHoles) {
        closeHoles();
    }
}

template <typename T>
void TriangleBoundary3D<T>::closeHoles()
{
    // Close the holes and assign inlets/outlets.
    std::vector<Lid> newlids = meshes[1].closeHoles();

    // Prepare the vector "triangle type", which later on will inform on
    //   the type of boundary condition implemented by a given triangle.
    //   The default, 0, stands for no-slip.
    pluint oldNumTriangles = triangleTagList.size();
    triangleTagList.resize(meshes[1].getNumTriangles());
    std::fill(triangleTagList.begin() + oldNumTriangles, triangleTagList.end(), 0);

    // Assign default functions to inlet/outlets to avoid undefined state.
    tagInletOutlet(newlids);
    lids.insert(lids.end(), newlids.begin(), newlids.end());
    assignLidVertexProperty();
}

template <typename T>
TriangleBoundary3D<T>::~TriangleBoundary3D()
{
    for (pluint iProp = 0; iProp < vertexProperties.size(); ++iProp) {
        delete vertexProperties[iProp];
    }
}

template <typename T>
TriangleBoundary3D<T>::TriangleBoundary3D(TriangleBoundary3D<T> const &rhs) :
    vertexLists(rhs.vertexLists),
    triangleTagList(rhs.triangleTagList),
    currentTagNum(rhs.currentTagNum),
    vertexTagList(rhs.vertexTagList),
    vertexProperties(rhs.vertexProperties.size()),
    lids(rhs.lids),
    margin(rhs.margin),
    physicalLocation(rhs.physicalLocation),
    dx(rhs.dx),
    topology(rhs.topology),
    vertexSet(rhs.vertexSet)
{
    emanatingEdgeLists[0] = rhs.emanatingEdgeLists[0];
    emanatingEdgeLists[1] = rhs.emanatingEdgeLists[1];
    edgeLists[0] = rhs.edgeLists[0];
    edgeLists[1] = rhs.edgeLists[1];

    defineMeshes();
    for (pluint iProp = 0; iProp < vertexProperties.size(); ++iProp) {
        vertexProperties[iProp] = rhs.vertexProperties[iProp]->clone();
    }
}

template <typename T>
void TriangleBoundary3D<T>::defineMeshes()
{
    meshes.clear();
    for (pluint iVertices = 0; iVertices < vertexLists.size(); ++iVertices) {
        meshes.push_back(
            TriangularSurfaceMesh<T>(vertexLists[iVertices], emanatingEdgeLists[0], edgeLists[0]));
        meshes.push_back(
            TriangularSurfaceMesh<T>(vertexLists[iVertices], emanatingEdgeLists[1], edgeLists[1]));
    }
}

template <typename T>
TriangleBoundary3D<T> &TriangleBoundary3D<T>::operator=(TriangleBoundary3D<T> const &rhs)
{
    TriangleBoundary3D<T>(rhs).swap(*this);
    return *this;
}

template <typename T>
void TriangleBoundary3D<T>::swap(TriangleBoundary3D<T> &rhs)
{
    vertexLists.swap(rhs.vertexLists);
    emanatingEdgeLists[0].swap(rhs.emanatingEdgeLists[0]);
    emanatingEdgeLists[1].swap(rhs.emanatingEdgeLists[1]);
    edgeLists[0].swap(rhs.edgeLists[0]);
    edgeLists[1].swap(rhs.edgeLists[1]);
    triangleTagList.swap(rhs.triangleTagList);
    std::swap(currentTagNum, rhs.currentTagNum);
    vertexTagList.swap(rhs.vertexTagList);
    vertexProperties.swap(rhs.vertexProperties);
    std::swap(lids, rhs.lids);
    std::swap(margin, rhs.margin);
    std::swap(physicalLocation, rhs.physicalLocation);
    std::swap(dx, rhs.dx);
    std::swap(topology, rhs.topology);
    std::swap(vertexSet, rhs.vertexSet);
    defineMeshes();
    rhs.defineMeshes();
}

template <typename T>
TriangleBoundary3D<T> const &TriangleBoundary3D<T>::select(
    plint whichTopology, plint whichVertices) const
{
    PLB_PRECONDITION(whichTopology == 0 || whichTopology == 1);
    PLB_PRECONDITION(whichVertices >= 0 && whichVertices < (plint)vertexLists.size());
    topology.top() = whichTopology;
    vertexSet.top() = whichVertices;
    return *this;
}

template <typename T>
TriangleBoundary3D<T> const &TriangleBoundary3D<T>::pushSelect(
    plint whichTopology, plint whichVertices) const
{
    PLB_PRECONDITION(whichTopology == 0 || whichTopology == 1);
    PLB_PRECONDITION(whichVertices >= 0 && whichVertices < (plint)vertexLists.size());
    topology.push(whichTopology);
    vertexSet.push(whichVertices);
    return *this;
}

template <typename T>
TriangleBoundary3D<T> const &TriangleBoundary3D<T>::popSelect() const
{
    PLB_PRECONDITION(topology.size() >= 2);
    PLB_PRECONDITION(vertexSet.size() >= 2);
    topology.pop();
    vertexSet.pop();
    return *this;
}

template <typename T>
void TriangleBoundary3D<T>::getSelection(plint &whichTopology, plint &whichVertices) const
{
    whichTopology = topology.top();
    whichVertices = vertexSet.top();
}

template <typename T>
plint TriangleBoundary3D<T>::currentMesh() const
{
    return 2 * vertexSet.top() + topology.top();
}

template <typename T>
TriangularSurfaceMesh<T> const &TriangleBoundary3D<T>::getMesh() const
{
    return meshes[currentMesh()];
}

template <typename T>
TriangularSurfaceMesh<T> &TriangleBoundary3D<T>::getMesh()
{
    return meshes[currentMesh()];
}

template <typename T>
plint TriangleBoundary3D<T>::getMargin() const
{
    return margin;
}

template <typename T>
plint TriangleBoundary3D<T>::getTag(plint iTriangle) const
{
    PLB_ASSERT(iTriangle < (plint)triangleTagList.size());
    return triangleTagList[iTriangle];
}

template <typename T>
VertexProperty3D<T> const *TriangleBoundary3D<T>::getVertexProperty(plint iVertex) const
{
    if (vertexTagList.empty()) {
        return 0;
    }
    PLB_ASSERT(iVertex < (plint)vertexTagList.size());
    return vertexProperties[vertexTagList[iVertex]];
}

template <typename T>
bool TriangleBoundary3D<T>::intersectSegment(
    plint iTriangle, AtomicBlock3D *boundaryArg, Array<T, 3> const &fromPoint,
    Array<T, 3> const &direction, Array<T, 3> &locatedPoint, T &distance,
    Array<T, 3> &wallNormal) const
{
    int flag = 0;  // Intersection with line segment.
    Array<T, 3> point2(fromPoint + direction);
    bool doesIntersect = getMesh().pointOnTriangle(
                             fromPoint, point2, flag, iTriangle, locatedPoint, wallNormal, distance)
                         == 1;
    return doesIntersect;
}

template <typename T>
Array<T, 3> TriangleBoundary3D<T>::computeContinuousNormal(
    Array<T, 3> const &p, plint iTriangle, bool isAreaWeighted) const
{
    return getMesh().computeContinuousNormal(p, iTriangle, isAreaWeighted);
}

template <typename T>
void TriangleBoundary3D<T>::cloneVertexSet(plint whichVertexSet)
{
    PLB_PRECONDITION(whichVertexSet < (plint)vertexLists.size() && whichVertexSet >= 0);
    vertexLists.push_back(vertexLists[whichVertexSet]);
    plint newVertexSet = (plint)vertexLists.size() - 1;
    plint numVertexOpen = vertexLists[newVertexSet].size();
    for (plint iLid = 0; iLid < (plint)lids.size(); ++iLid) {
        numVertexOpen -= lids[iLid].numAddedVertices;
    }
    meshes.push_back(TriangularSurfaceMesh<T>(
        vertexLists[newVertexSet], emanatingEdgeLists[0], edgeLists[0], numVertexOpen));
    meshes.push_back(
        TriangularSurfaceMesh<T>(vertexLists[newVertexSet], emanatingEdgeLists[1], edgeLists[1]));
}

template <typename T>
std::vector<Lid> const &TriangleBoundary3D<T>::getInletOutlet() const
{
    PLB_PRECONDITION(topology.top() == 1);
    return lids;
}

template <typename T>
std::vector<Lid> TriangleBoundary3D<T>::getInletOutlet(plint sortDirection) const
{
    PLB_PRECONDITION(topology.top() == 1);
    std::vector<Lid> lids_copy(lids);
    std::sort(lids_copy.begin(), lids_copy.end(), LidLessThan<T>(sortDirection, getMesh()));
    return lids_copy;
}

template <typename T>
std::vector<plint> TriangleBoundary3D<T>::getInletOutletIds(plint sortDirection) const
{
    std::map<plint, plint> triangleToOriginalLid;
    for (pluint iLid = 0; iLid < lids.size(); ++iLid) {
        triangleToOriginalLid[lids[iLid].firstTriangle] = iLid;
    }
    std::vector<Lid> tmpLids(lids);
    std::sort(tmpLids.begin(), tmpLids.end(), LidLessThan<T>(sortDirection, getMesh()));
    std::vector<plint> ids(tmpLids.size());
    for (pluint iLid = 0; iLid < tmpLids.size(); ++iLid) {
        plint originalId = triangleToOriginalLid[tmpLids[iLid].firstTriangle];
        // It is assumed that the ID of the original wall is 0,
        // while the lids have continuous IDs starting at 1.
        ids[iLid] = originalId + 1;
    }
    return ids;
}

template <typename T>
void TriangleBoundary3D<T>::getLidProperties(
    plint sortDirection, std::vector<Array<T, 3> > &normal, std::vector<Array<T, 3> > &center,
    std::vector<T> &radius) const
{
    // Lid properties can only be computed in a closed mesh, by definition.
    PLB_PRECONDITION(topology.top() == 1);
    std::vector<Lid> tmpLids(lids);
    std::sort(tmpLids.begin(), tmpLids.end(), LidLessThan<T>(sortDirection, getMesh()));
    normal.resize(tmpLids.size());
    center.resize(tmpLids.size());
    radius.resize(tmpLids.size());
    for (pluint iLid = 0; iLid < tmpLids.size(); ++iLid) {
        normal[iLid] = computeNormal(getMesh(), tmpLids[iLid]);
        center[iLid] = computeBaryCenter(getMesh(), tmpLids[iLid]);
        radius[iLid] = (T)0.5
                       * (computeInnerRadius(getMesh(), tmpLids[iLid])
                          + computeOuterRadius(getMesh(), tmpLids[iLid]));
    }
}

template <typename T>
void TriangleBoundary3D<T>::tagInletOutlet(std::vector<Lid> &newLids)
{
    // Inlet/Outlet can only be set for a closed mesh, by definition.
    PLB_PRECONDITION(topology.top() == 1);

    for (pluint iLid = 0; iLid < newLids.size(); ++iLid) {
        ++currentTagNum;  // Tag 0 is for default wall portions.
        newLids[iLid].tag = currentTagNum;
        plint firstTriangle = newLids[iLid].firstTriangle;
        plint numTriangles = newLids[iLid].numTriangles;
        for (plint iTriangle = firstTriangle; iTriangle < firstTriangle + numTriangles; ++iTriangle)
        {
            triangleTagList[iTriangle] = currentTagNum;
        }
    }
}

template <typename T>
template <typename DomainFunctional>
plint TriangleBoundary3D<T>::tagDomain(DomainFunctional functional)
{
    // Make sure we're working with the closed mesh.
    PLB_PRECONDITION(topology.top() == 1);
    ++currentTagNum;
    plint newTag = currentTagNum;
    for (plint iTriangle = 0; iTriangle < getMesh().getNumTriangles(); ++iTriangle) {
        bool isInside = true;
        for (plint iVertex = 0; iVertex < 3; ++iVertex) {
            Array<T, 3> vertex = getMesh().getVertex(iTriangle, iVertex);
            if (!functional(vertex)) {
                isInside = false;
                break;
            }
        }
        if (isInside) {
            triangleTagList[iTriangle] = newTag;
        }
    }
    return newTag;
}

template <typename T>
template <typename DomainFunctional>
plint TriangleBoundary3D<T>::tagDomain(
    DomainFunctional functional, Array<T, 3> normal, T angleTolerance, plint previousTag)
{
    // Make sure we're working with the closed mesh.
    PLB_PRECONDITION(topology.top() == 1);
    ++currentTagNum;
    plint newTag = currentTagNum;
    for (plint iTriangle = 0; iTriangle < getMesh().getNumTriangles(); ++iTriangle) {
        if (previousTag < 0 || triangleTagList[iTriangle] == previousTag) {
            bool isInside = true;
            for (plint iVertex = 0; iVertex < 3; ++iVertex) {
                Array<T, 3> vertex = getMesh().getVertex(iTriangle, iVertex);
                if (!functional(vertex)) {
                    isInside = false;
                    break;
                }
            }
            if (isInside) {
                Array<T, 3> triangleNormal = getMesh().computeTriangleNormal(iTriangle);
                T normN = norm(triangleNormal);  // We decide not to tag zero-area triangles.
                if (!util::isZero(normN)
                    && std::fabs(angleBetweenVectors(normal, triangleNormal) < angleTolerance)) {
                    triangleTagList[iTriangle] = newTag;
                }
            }
        }
    }
    return newTag;
}

template <typename T>
plint TriangleBoundary3D<T>::tagLids(Cuboid<T> const &c)
{
    // Make sure we're working with the closed mesh.
    PLB_PRECONDITION(topology.top() == 1);
    ++currentTagNum;
    plint newTag = currentTagNum;

    for (pluint iLid = 0; iLid < lids.size(); iLid++) {
        Array<T, 3> const &bc = getMesh().getVertex(lids[iLid].centerVertex);
        if (contained(bc, c)) {
            plint firstTriangle = lids[iLid].firstTriangle;
            plint numTriangles = lids[iLid].numTriangles;
            for (plint iTriangle = firstTriangle; iTriangle < firstTriangle + numTriangles;
                 ++iTriangle) {
                triangleTagList[iTriangle] = newTag;
            }
        }
    }

    return newTag;
}

template <typename T>
template <typename DomainFunctional>
plint TriangleBoundary3D<T>::setVertexProperty(
    VertexProperty3D<T> const &property, DomainFunctional functional)
{
    // Make sure we're working with the closed mesh.
    PLB_PRECONDITION(topology.top() == 1);
    if (vertexTagList.empty()) {
        PLB_ASSERT(vertexProperties.empty());
        vertexTagList.resize(getMesh().getNumVertices());
        std::fill(vertexTagList.begin(), vertexTagList.end(), 0);
        vertexProperties.push_back(0);
    }
    vertexProperties.push_back(property.clone());
    plint nextTag = (plint)vertexProperties.size() - 1;
    for (plint iVertex = 0; iVertex < getMesh().getNumVertices(); ++iVertex) {
        Array<T, 3> const &vertex = getMesh().getVertex(iVertex);
        bool isInside = functional(vertex);
        if (isInside) {
            vertexTagList[iVertex] = nextTag;
        }
    }
    return nextTag;
}

template <typename T>
void TriangleBoundary3D<T>::assignLidVertexProperty()
{
    // Make sure we're working with the closed mesh.
    PLB_PRECONDITION(topology.top() == 1);
    std::vector<Lid> const &lids = getInletOutlet();
    if (lids.empty())
        return;

    if (vertexTagList.empty()) {
        PLB_ASSERT(vertexProperties.empty());
        vertexTagList.resize(getMesh().getNumVertices());
        std::fill(vertexTagList.begin(), vertexTagList.end(), 0);
        vertexProperties.push_back(0);
    }

    vertexProperties.push_back(new InletOutletProperty3D<T>);
    plint lidVertexTag = (plint)vertexProperties.size() - 1;
    for (pluint iLid = 0; iLid < lids.size(); ++iLid) {
        plint centerVertex = lids[iLid].centerVertex;
        vertexTagList[centerVertex] = lidVertexTag;
    }
}

/******** class TriangleFlowShape3D ****************************************/

template <typename T, class SurfaceData>
TriangleFlowShape3D<T, SurfaceData>::TriangleFlowShape3D(
    TriangleBoundary3D<T> const &boundary_, BoundaryProfiles3D<T, SurfaceData> const &profiles_) :
    boundary(boundary_), profiles(profiles_), voxelFlags(0), hashContainer(0), boundaryArg(0)
{ }

template <typename T, class SurfaceData>
bool TriangleFlowShape3D<T, SurfaceData>::isInside(Dot3D const &location) const
{
    PLB_PRECONDITION(voxelFlags);
    Dot3D localPos = location - voxelFlags->getLocation();
    return voxelFlag::insideFlag(voxelFlags->get(localPos.x, localPos.y, localPos.z));
}

template <typename T, class SurfaceData>
bool TriangleFlowShape3D<T, SurfaceData>::isOutside(Dot3D const &location) const
{
    PLB_PRECONDITION(voxelFlags);
    Dot3D localPos = location - voxelFlags->getLocation();
    return voxelFlag::outsideFlag(voxelFlags->get(localPos.x, localPos.y, localPos.z));
}

template <typename T, class SurfaceData>
bool TriangleFlowShape3D<T, SurfaceData>::pointOnSurface(
    Array<T, 3> const &fromPoint, Array<T, 3> const &direction, Array<T, 3> &locatedPoint,
    T &distance, Array<T, 3> &wallNormal, SurfaceData &surfaceData, OffBoundary::Type &bdType,
    plint &id) const
{
    PLB_PRECONDITION(hashContainer);  // Make sure these arguments have
    PLB_PRECONDITION(boundaryArg);    //   been provided by the user through
                                      //   the clone function.
    static const T maxDistance = std::sqrt((T)3);
    Array<T, 2> xRange(fromPoint[0] - maxDistance, fromPoint[0] + maxDistance);
    Array<T, 2> yRange(fromPoint[1] - maxDistance, fromPoint[1] + maxDistance);
    Array<T, 2> zRange(fromPoint[2] - maxDistance, fromPoint[2] + maxDistance);
    TriangleHash<T> triangleHash(*hashContainer);
    std::vector<plint> possibleTriangles;
    if (id >= 0 && id < boundary.getMesh().getNumTriangles()) {
        possibleTriangles.push_back(id);
    } else {
        triangleHash.getTriangles(xRange, yRange, zRange, possibleTriangles);
    }

    Array<T, 3> tmpLocatedPoint;
    T tmpDistance;
    Array<T, 3> tmpNormal;
    T shortestDistance = T();
    plint locatedTriangle = -1;

    for (pluint iPossible = 0; iPossible < possibleTriangles.size(); ++iPossible) {
        plint iTriangle = possibleTriangles[iPossible];
        if (boundary.intersectSegment(
                iTriangle, boundaryArg, fromPoint, direction, tmpLocatedPoint, tmpDistance,
                tmpNormal))
        {
            if (locatedTriangle == -1 || tmpDistance < shortestDistance) {
                shortestDistance = tmpDistance;
                locatedTriangle = iTriangle;
                locatedPoint = tmpLocatedPoint;
                distance = tmpDistance;
                wallNormal = tmpNormal;
                profiles.getProfile(boundary, iTriangle)
                    .getData(locatedPoint, iTriangle, boundaryArg, surfaceData, bdType);
            }
        }
    }
    if (locatedTriangle != -1) {
        id = locatedTriangle;
        return true;
    } else {
        return false;
    }
}

template <typename T, class SurfaceData>
Array<T, 3> TriangleFlowShape3D<T, SurfaceData>::computeContinuousNormal(
    Array<T, 3> const &p, plint id, bool isAreaWeighted) const
{
    return boundary.computeContinuousNormal(p, id, isAreaWeighted);
}

template <typename T, class SurfaceData>
bool TriangleFlowShape3D<T, SurfaceData>::intersectsSurface(
    Array<T, 3> const &p1, Array<T, 3> const &p2, plint &id) const
{
    PLB_PRECONDITION(hashContainer);  // Make sure these arguments have
    PLB_PRECONDITION(boundaryArg);    //   been provided by the user through
                                      //   the clone function.
    static const T maxDistance = std::sqrt((T)3);
    Array<T, 2> xRange(p1[0] - maxDistance, p1[0] + maxDistance);
    Array<T, 2> yRange(p1[1] - maxDistance, p1[1] + maxDistance);
    Array<T, 2> zRange(p1[2] - maxDistance, p1[2] + maxDistance);
    TriangleHash<T> triangleHash(*hashContainer);
    std::vector<plint> possibleTriangles;
    if (id >= 0 && id < boundary.getMesh().getNumTriangles()) {
        possibleTriangles.push_back(id);
    } else {
        triangleHash.getTriangles(xRange, yRange, zRange, possibleTriangles);
    }

    std::vector<plint> selection;
    for (pluint iPossible = 0; iPossible < possibleTriangles.size(); ++iPossible) {
        plint iTriangle = possibleTriangles[iPossible];
        int flag = 0;
        Array<T, 3> intersection, normal;
        T distance;
        if (boundary.getMesh().pointOnTriangle(
                p1, p2, flag, iTriangle, intersection, normal, distance)) {
            selection.push_back(iTriangle);
        }
    }
    if (selection.empty()) {
        return false;
    } else if (selection.size() == 1) {
        id = selection[0];
        return true;
    } else {
        Array<T, 3> locatedPoint;
        T distance;
        Array<T, 3> wallNormal;
        SurfaceData surfaceData;
        OffBoundary::Type bdType;
        return pointOnSurface(
            p1, p2 - p1, locatedPoint, distance, wallNormal, surfaceData, bdType, id);
    }
}

template <typename T, class SurfaceData>
plint TriangleFlowShape3D<T, SurfaceData>::getTag(plint id) const
{
    return boundary.getTag(id);
}

template <typename T, class SurfaceData>
bool TriangleFlowShape3D<T, SurfaceData>::distanceToSurface(
    Array<T, 3> const &point, T &distance, bool &isBehind) const
{
    PLB_PRECONDITION(hashContainer);  // Make sure these arguments have
    PLB_PRECONDITION(boundaryArg);    //   been provided by the user through
                                      //   the clone function.
    T maxDistance = std::sqrt((T)3);
    Array<T, 2> xRange(point[0] - maxDistance, point[0] + maxDistance);
    Array<T, 2> yRange(point[1] - maxDistance, point[1] + maxDistance);
    Array<T, 2> zRange(point[2] - maxDistance, point[2] + maxDistance);
    TriangleHash<T> triangleHash(*hashContainer);
    std::vector<plint> possibleTriangles;
    triangleHash.getTriangles(xRange, yRange, zRange, possibleTriangles);

    T tmpDistance;
    bool tmpIsBehind;
    bool triangleFound = false;

    for (pluint iPossible = 0; iPossible < possibleTriangles.size(); ++iPossible) {
        plint iTriangle = possibleTriangles[iPossible];
        boundary.getMesh().distanceToTriangle(point, iTriangle, tmpDistance, tmpIsBehind);
        if (!triangleFound || tmpDistance < distance) {
            distance = tmpDistance;
            isBehind = tmpIsBehind;
            triangleFound = true;
        }
    }
    return triangleFound;
}

template <typename T, class SurfaceData>
TriangleFlowShape3D<T, SurfaceData> *TriangleFlowShape3D<T, SurfaceData>::clone() const
{
    return new TriangleFlowShape3D<T, SurfaceData>(*this);
}

template <typename T, class SurfaceData>
TriangleFlowShape3D<T, SurfaceData> *TriangleFlowShape3D<T, SurfaceData>::clone(
    std::vector<AtomicBlock3D *> args) const
{
    PLB_PRECONDITION(args.size() == 3);
    TriangleFlowShape3D<T, SurfaceData> *newShape = new TriangleFlowShape3D<T, SurfaceData>(*this);
    newShape->voxelFlags = dynamic_cast<ScalarField3D<int> *>(args[0]);
    newShape->hashContainer = dynamic_cast<AtomicContainerBlock3D *>(args[1]);
    newShape->boundaryArg = args[2];
    PLB_ASSERT(newShape->voxelFlags);
    PLB_ASSERT(newShape->hashContainer);
    PLB_ASSERT(newShape->boundaryArg);
    return newShape;
}

/******** class VoxelizedDomain3D *****************************************/

template <typename T>
VoxelizedDomain3D<T>::VoxelizedDomain3D(
    TriangleBoundary3D<T> const &boundary_, int flowType_, plint extraLayer_, plint borderWidth_,
    plint envelopeWidth_, plint blockSize_, plint gridLevel_, bool dynamicMesh_) :
    flowType(flowType_), borderWidth(borderWidth_), boundary(boundary_)
{
    PLB_ASSERT(flowType == voxelFlag::inside || flowType == voxelFlag::outside);
    PLB_ASSERT(boundary.getMargin() >= borderWidth);
    if (dynamicMesh_) {
        boundary.pushSelect(1, 1);  // Closed, Dynamic.
    } else {
        boundary.pushSelect(1, 0);  // Closed, Static.
    }
    std::unique_ptr<MultiScalarField3D<int> > fullVoxelMatrix(
        voxelize(boundary.getMesh(), boundary.getMargin() + extraLayer_, borderWidth));
    fullVoxelMatrix->setRefinementLevel(gridLevel_);
    createSparseVoxelMatrix(*fullVoxelMatrix, blockSize_, envelopeWidth_);
    createTriangleHash();
    boundary.popSelect();
}

template <typename T>
VoxelizedDomain3D<T>::VoxelizedDomain3D(
    TriangleBoundary3D<T> const &boundary_, int flowType_, Box3D const &boundingBox,
    plint borderWidth_, plint envelopeWidth_, plint blockSize_, plint gridLevel_,
    bool dynamicMesh_) :
    flowType(flowType_), borderWidth(borderWidth_), boundary(boundary_)
{
    PLB_ASSERT(flowType == voxelFlag::inside || flowType == voxelFlag::outside);
    PLB_ASSERT(boundary.getMargin() >= borderWidth);
    if (dynamicMesh_) {
        boundary.pushSelect(1, 1);  // Closed, Dynamic.
    } else {
        boundary.pushSelect(1, 0);  // Closed, Static.
    }
    std::unique_ptr<MultiScalarField3D<int> > fullVoxelMatrix(
        voxelize(boundary.getMesh(), boundingBox, borderWidth));
    fullVoxelMatrix->setRefinementLevel(gridLevel_);
    createSparseVoxelMatrix(*fullVoxelMatrix, blockSize_, envelopeWidth_);
    createTriangleHash();
    boundary.popSelect();
}

template <typename T>
VoxelizedDomain3D<T>::VoxelizedDomain3D(
    TriangleBoundary3D<T> const &boundary_, int flowType_, Box3D const &boundingBox,
    plint borderWidth_, plint envelopeWidth_, plint blockSize_, Box3D const &seed, plint gridLevel_,
    bool dynamicMesh_) :
    flowType(flowType_), borderWidth(borderWidth_), boundary(boundary_)
{
    PLB_ASSERT(flowType == voxelFlag::inside || flowType == voxelFlag::outside);
    PLB_ASSERT(boundary.getMargin() >= borderWidth);
    if (dynamicMesh_) {
        boundary.pushSelect(1, 1);  // Closed, Dynamic.
    } else {
        boundary.pushSelect(1, 0);  // Closed, Static.
    }
    std::unique_ptr<MultiScalarField3D<int> > fullVoxelMatrix(
        voxelize(boundary.getMesh(), boundingBox, borderWidth, seed));
    fullVoxelMatrix->setRefinementLevel(gridLevel_);
    createSparseVoxelMatrix(*fullVoxelMatrix, blockSize_, envelopeWidth_);
    createTriangleHash();
    boundary.popSelect();
}

template <typename T>
void VoxelizedDomain3D<T>::createSparseVoxelMatrix(
    MultiScalarField3D<int> &fullVoxelMatrix, plint blockSize_, plint envelopeWidth_)
{
    if (blockSize_ > 0) {
        computeSparseVoxelMatrix(fullVoxelMatrix, blockSize_, envelopeWidth_);
    } else {
        // blockSize=0 means: don't create a sparse representation of the
        //   multi-block. The envelope of the voxel-matrix needs to be extended
        //   for future usage, though.
        extendEnvelopeWidth(fullVoxelMatrix, envelopeWidth_);
    }
}

template <typename T>
VoxelizedDomain3D<T>::VoxelizedDomain3D(VoxelizedDomain3D<T> const &rhs) :
    boundary(rhs.boundary),
    voxelMatrix(new MultiScalarField3D<int>(*rhs.voxelMatrix)),
    triangleHash(new MultiContainerBlock3D(*rhs.triangleHash))
{ }

template <typename T>
VoxelizedDomain3D<T>::~VoxelizedDomain3D()
{
    delete voxelMatrix;
    delete triangleHash;
}

template <typename T>
MultiScalarField3D<int> &VoxelizedDomain3D<T>::getVoxelMatrix()
{
    return *voxelMatrix;
}

template <typename T>
MultiScalarField3D<int> const &VoxelizedDomain3D<T>::getVoxelMatrix() const
{
    return *voxelMatrix;
}

template <typename T>
MultiContainerBlock3D &VoxelizedDomain3D<T>::getTriangleHash()
{
    return *triangleHash;
}

template <typename T>
template <class ParticleFieldT>
void VoxelizedDomain3D<T>::adjustVoxelization(
    MultiParticleField3D<ParticleFieldT> &particles, bool dynamicMesh)
{
    if (dynamicMesh) {
        boundary.pushSelect(1, 1);  // Closed, Dynamic.
    } else {
        boundary.pushSelect(1, 0);  // Closed, Static.
    }
    reCreateTriangleHash(particles);
    MultiScalarField3D<int> *newVoxelMatrix =
        revoxelize(boundary.getMesh(), *voxelMatrix, *triangleHash, borderWidth).release();
    std::swap(voxelMatrix, newVoxelMatrix);
    delete newVoxelMatrix;
    boundary.popSelect();
}

template <typename T>
void VoxelizedDomain3D<T>::reparallelize(MultiBlockRedistribute3D const &redistribute)
{
    MultiBlockManagement3D newManagement =
        redistribute.redistribute(voxelMatrix->getMultiBlockManagement());
    MultiScalarField3D<int> *newVoxelMatrix = new MultiScalarField3D<int>(
        newManagement, voxelMatrix->getBlockCommunicator().clone(),
        voxelMatrix->getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiScalarAccess<int>(), 0);
    copyNonLocal(*voxelMatrix, *newVoxelMatrix, voxelMatrix->getBoundingBox());
    std::swap(voxelMatrix, newVoxelMatrix);
    delete newVoxelMatrix;
    delete triangleHash;
    createTriangleHash();
}

template <typename T>
void VoxelizedDomain3D<T>::reparallelize(MultiBlockManagement3D const &newManagement)
{
    MultiScalarField3D<int> *newVoxelMatrix = new MultiScalarField3D<int>(
        newManagement, voxelMatrix->getBlockCommunicator().clone(),
        voxelMatrix->getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiScalarAccess<int>(), 0);
    copyNonLocal(*voxelMatrix, *newVoxelMatrix, voxelMatrix->getBoundingBox());
    std::swap(voxelMatrix, newVoxelMatrix);
    delete newVoxelMatrix;
    delete triangleHash;
    createTriangleHash();
}

template <typename T>
MultiBlockManagement3D const &VoxelizedDomain3D<T>::getMultiBlockManagement() const
{
    return voxelMatrix->getMultiBlockManagement();
}

template <typename T>
void VoxelizedDomain3D<T>::computeSparseVoxelMatrix(
    MultiScalarField3D<int> &fullVoxelMatrix, plint blockSize, plint envelopeWidth)
{
    // Initialized to zero.
    MultiScalarField3D<int> domainMatrix((MultiBlock3D const &)fullVoxelMatrix);
    setToConstant(domainMatrix, fullVoxelMatrix, flowType, domainMatrix.getBoundingBox(), 1);
    setToConstant(
        domainMatrix, fullVoxelMatrix, voxelFlag::borderFlag(flowType),
        domainMatrix.getBoundingBox(), 1);
    for (int iLayer = 1; iLayer <= boundary.getMargin(); ++iLayer) {
        addLayer(domainMatrix, domainMatrix.getBoundingBox(), iLayer);
    }
    MultiBlockManagement3D sparseBlockManagement = computeSparseManagement(
        *plb::reparallelize(domainMatrix, blockSize, blockSize, blockSize), envelopeWidth);

    voxelMatrix = new MultiScalarField3D<int>(
        sparseBlockManagement, fullVoxelMatrix.getBlockCommunicator().clone(),
        fullVoxelMatrix.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiScalarAccess<int>(), voxelFlag::undetermined);
    copyNonLocal(fullVoxelMatrix, *voxelMatrix, voxelMatrix->getBoundingBox());
}

template <typename T>
void VoxelizedDomain3D<T>::extendEnvelopeWidth(
    MultiScalarField3D<int> &fullVoxelMatrix, plint envelopeWidth)
{
    MultiBlockManagement3D const &oldManagement = fullVoxelMatrix.getMultiBlockManagement();
    MultiBlockManagement3D newManagement(
        oldManagement.getSparseBlockStructure(), oldManagement.getThreadAttribution().clone(),
        envelopeWidth, oldManagement.getRefinementLevel());
    voxelMatrix = new MultiScalarField3D<int>(
        newManagement, fullVoxelMatrix.getBlockCommunicator().clone(),
        fullVoxelMatrix.getCombinedStatistics().clone(),
        defaultMultiBlockPolicy3D().getMultiScalarAccess<int>(), voxelFlag::undetermined);
    plb::copy(fullVoxelMatrix, *voxelMatrix, fullVoxelMatrix.getBoundingBox());
}

template <typename T>
void VoxelizedDomain3D<T>::createTriangleHash()
{
    triangleHash = new MultiContainerBlock3D(*voxelMatrix);
    std::vector<MultiBlock3D *> hashArg;
    hashArg.push_back(triangleHash);
    applyProcessingFunctional(
        new CreateTriangleHash<T>(boundary.getMesh()), triangleHash->getBoundingBox(), hashArg);
}

template <typename T>
template <class ParticleFieldT>
void VoxelizedDomain3D<T>::reCreateTriangleHash(MultiParticleField3D<ParticleFieldT> &particles)
{
    // The lids are non-parallel, an info which must be provided
    //   to the hash algorithm by means of a list of barycenters.
    //   This is a necessary and sufficient information because
    //   all triangles connected to the barycenters are non-parallel,
    //   and all non-parallel triangles are connected to a barycenter.
    std::vector<Lid> const &lids = boundary.getInletOutlet();
    plint numLids = (plint)lids.size();
    std::vector<plint> lidBaryCenters(numLids);
    for (plint iLid = 0; iLid < numLids; ++iLid) {
        lidBaryCenters[iLid] = lids[iLid].centerVertex;
    }

    std::vector<MultiBlock3D *> hashParticleArg;
    hashParticleArg.push_back(triangleHash);
    hashParticleArg.push_back(&particles);
    applyProcessingFunctional(
        new ReAssignTriangleHash<T, ParticleFieldT>(boundary.getMesh(), lidBaryCenters),
        triangleHash->getBoundingBox(), hashParticleArg);
}

/* ******** DetectBorderLineFunctional3D ************************************* */

template <typename T>
void addLayer(MultiScalarField3D<T> &matrix, Box3D const &domain, T previousLayer)
{
    applyProcessingFunctional(new AddLayerFunctional3D<T>(previousLayer), domain, matrix);
}

template <typename T>
AddLayerFunctional3D<T>::AddLayerFunctional3D(T previousLayer_) : previousLayer(previousLayer_)
{ }

template <typename T>
void AddLayerFunctional3D<T>::process(Box3D domain, ScalarField3D<T> &voxels)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                for (plint dx = -1; dx <= 1; ++dx)
                    for (plint dy = -1; dy <= 1; ++dy)
                        for (plint dz = -1; dz <= 1; ++dz)
                            if (!(dx == 0 && dy == 0 && dz == 0)) {
                                plint nextX = iX + dx;
                                plint nextY = iY + dy;
                                plint nextZ = iZ + dz;
                                if (voxels.get(iX, iY, iZ) == 0
                                    && voxels.get(nextX, nextY, nextZ) == previousLayer) {
                                    voxels.get(iX, iY, iZ) = previousLayer + 1;
                                }
                            }
            }
        }
    }
}

template <typename T>
AddLayerFunctional3D<T> *AddLayerFunctional3D<T>::clone() const
{
    return new AddLayerFunctional3D<T>(*this);
}

template <typename T>
void AddLayerFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT AddLayerFunctional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

}  // namespace plb

#endif  // TRIANGLE_BOUNDARY_3D_HH
