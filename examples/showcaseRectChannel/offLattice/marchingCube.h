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

#ifndef MARCHING_CUBE_H
#define MARCHING_CUBE_H

#include <vector>

#include "core/globalDefs.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "offLattice/boundaryShapes3D.h"
#include "offLattice/offLatticeBoundaryProfiles3D.h"
#include "offLattice/triangleBoundary3D.h"
#include "offLattice/triangleSet.h"
#include "offLattice/triangularSurfaceMesh.h"

namespace plb {

template <typename T>
class IsoSurfaceDefinition3D;

/// Get an iso-surface by means of the marching cube algorithms.
/** The iso-surface is defined in very generic terms by the isoSurfaceDefinition,
 * and the surfDefinitionArgs are whatever arguments the isoSurfaceDefinition
 * needs. The isoSurfaceDefinition can compute a finite amount of iso-surfaces,
 * the IDs of which are provided by the last argument. If the last argument is
 * omitted, all available iso-surfaces are computed.
 * The iso-surface is returned as a set of triangles, in the first argument.
 **/
template <typename T>
void isoSurfaceMarchingCube(
    std::vector<typename TriangleSet<T>::Triangle> &triangles,
    std::vector<MultiBlock3D *> surfDefinitionArgs, IsoSurfaceDefinition3D<T> *isoSurfaceDefinition,
    Box3D const &domain, std::vector<plint> surfaceIds = std::vector<plint>());

/// This wrapper call to the marching-cube algorithm remeshes the surface of a voxelized domain.
template <typename T>
void isoSurfaceMarchingCube(
    std::vector<typename TriangleSet<T>::Triangle> &triangles,
    VoxelizedDomain3D<T> &voxelizedDomain, Box3D const &domain);

/// This wrapper call to the marching-cube algorithm computes iso-surfaces from a scalar-field.
template <typename T>
void isoSurfaceMarchingCube(
    std::vector<typename TriangleSet<T>::Triangle> &triangles, MultiScalarField3D<T> &scalarField,
    std::vector<T> const &isoLevels, Box3D const &domain);

/// This wrapper call to the marching-cube algorithm computes iso-surfaces from an analytical
/// description.
template <typename T, class Function>
void isoSurfaceMarchingCube(
    std::vector<typename TriangleSet<T>::Triangle> &triangles, MultiBlock3D &block,
    Function const &function, Box3D const &domain);

template <typename T, template <typename U> class Descriptor>
TriangleSet<T> vofToTriangles(MultiScalarField3D<T> &scalarField, T threshold, Box3D domain);

template <typename T, template <typename U> class Descriptor>
TriangleSet<T> vofToTriangles(MultiScalarField3D<T> &scalarField, T threshold);

template <typename T>
class IsoSurfaceDefinition3D {
public:
    virtual ~IsoSurfaceDefinition3D() { }
    virtual bool isInside(plint surfaceId, Array<plint, 3> const &position) const = 0;
    virtual bool isValid(Array<plint, 3> const &position) const
    {
        return true;
    }
    virtual Array<T, 3> getSurfacePosition(
        plint surfaceId, Array<plint, 3> const &p1, Array<plint, 3> const &p2) const = 0;
    virtual void setArguments(std::vector<AtomicBlock3D *> const &arguments) = 0;
    virtual IsoSurfaceDefinition3D<T> *clone() const = 0;
    virtual plint getNumArgs() const = 0;
    virtual std::vector<plint> getSurfaceIds() const = 0;

public:
    bool edgeIsValid(plint iX, plint iY, plint iZ, int edge) const;
};

template <typename T>
class ScalarFieldIsoSurface3D : public IsoSurfaceDefinition3D<T> {
public:
    ScalarFieldIsoSurface3D(std::vector<T> const &isoValues_);
    virtual bool isInside(plint surfaceId, Array<plint, 3> const &position) const;
    virtual Array<T, 3> getSurfacePosition(
        plint surfaceId, Array<plint, 3> const &p1, Array<plint, 3> const &p2) const;
    virtual void setArguments(std::vector<AtomicBlock3D *> const &arguments);
    virtual ScalarFieldIsoSurface3D<T> *clone() const;
    virtual plint getNumArgs() const
    {
        return 1;
    }
    virtual std::vector<plint> getSurfaceIds() const;

private:
    std::vector<T> isoValues;
    ScalarField3D<T> *scalar;
    Dot3D location;
};

template <typename T, class Function>
class AnalyticalIsoSurface3D : public IsoSurfaceDefinition3D<T> {
public:
    AnalyticalIsoSurface3D(Function const &function_) : function(function_) { }
    virtual bool isInside(plint surfaceId, Array<plint, 3> const &position) const;
    virtual Array<T, 3> getSurfacePosition(
        plint surfaceId, Array<plint, 3> const &p1, Array<plint, 3> const &p2) const;
    virtual void setArguments(std::vector<AtomicBlock3D *> const &arguments) { }
    virtual AnalyticalIsoSurface3D<T, Function> *clone() const;
    virtual plint getNumArgs() const
    {
        return 0;
    }
    virtual std::vector<plint> getSurfaceIds() const;

private:
    class WrappedIsInside {
    public:
        WrappedIsInside(Array<T, 3> const &p1_, Array<T, 3> const &p2_, Function const &function_) :
            p1(p1_), p2(p2_), function(function_)
        { }
        T operator()(T position) const
        {
            if (function.floatIsInside(p1 + position * (p2 - p1))) {
                return (T)1;
            } else {
                return (T)-1;
            }
        }

    private:
        Array<T, 3> p1, p2;
        Function function;
    };

private:
    Function function;
};

template <typename T, class SurfaceData>
class BoundaryShapeIsoSurface3D : public IsoSurfaceDefinition3D<T> {
public:
    BoundaryShapeIsoSurface3D(BoundaryShape3D<T, SurfaceData> *shape_);
    virtual ~BoundaryShapeIsoSurface3D();
    BoundaryShapeIsoSurface3D(BoundaryShapeIsoSurface3D<T, SurfaceData> const &rhs);
    BoundaryShapeIsoSurface3D<T, SurfaceData> &operator=(
        BoundaryShapeIsoSurface3D<T, SurfaceData> const &rhs);
    void swap(BoundaryShapeIsoSurface3D<T, SurfaceData> &rhs);
    virtual bool isInside(plint surfaceId, Array<plint, 3> const &position) const;
    virtual Array<T, 3> getSurfacePosition(
        plint surfaceId, Array<plint, 3> const &p1, Array<plint, 3> const &p2) const;
    /// Arguments are:
    /// 1. voxelizedDomain.getVoxelMatrix()
    /// 2. voxelizedDomain.getTriangleHash()
    /// 3. Argument needed by the boundary profiles.
    virtual void setArguments(std::vector<AtomicBlock3D *> const &arguments);
    virtual BoundaryShapeIsoSurface3D<T, SurfaceData> *clone() const;
    virtual plint getNumArgs() const
    {
        return 3;
    }
    virtual std::vector<plint> getSurfaceIds() const;

private:
    BoundaryShape3D<T, SurfaceData> *shape;
};

template <typename T>
class MarchingCubeSurfaces3D : public BoxProcessingFunctional3D {
public:
    typedef typename TriangleSet<T>::Triangle Triangle;

public:
    MarchingCubeSurfaces3D(
        std::vector<plint> surfaceIds_, IsoSurfaceDefinition3D<T> *isoSurface_,
        bool edgeOrientedData_ = false);
    ~MarchingCubeSurfaces3D();
    MarchingCubeSurfaces3D(MarchingCubeSurfaces3D<T> const &rhs);
    MarchingCubeSurfaces3D<T> &operator=(MarchingCubeSurfaces3D<T> const &rhs);
    void swap(MarchingCubeSurfaces3D<T> &rhs);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual void defaultImplementation(Box3D domain, AtomicContainerBlock3D *triangleContainer);
    virtual void edgeOriented(Box3D domain, AtomicContainerBlock3D *triangleContainer);
    virtual MarchingCubeSurfaces3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    void setEdgeOrientedEnvelope(plint edgeOrientedEnvelope_)
    {
        edgeOrientedEnvelope = edgeOrientedEnvelope_;
    }

public:
    class TriangleSetData : public ContainerBlockData {
    public:
        std::vector<Triangle> triangles;
        virtual TriangleSetData *clone() const
        {
            return new TriangleSetData(*this);
        }
    };
    class EdgeOrientedTriangleSetData : public ContainerBlockData {
    public:
        /// Holds the topology of a triangle, by stating on which edges its
        /// vertices are situated.
        struct OnEdgeTriangle {
            /** The edge-attribution is stored for two vertices only, because
             * the first vertex is, by definition, stored on the current edge.
             * An edge is defined by four coordinates: 3 coordinates for the cell,
             * and an identifier (0,1, or 2) for the edge.
             */
            Array<plint, 4> vertex1, vertex2;
        };

    public:
        EdgeOrientedTriangleSetData(plint nx, plint ny, plint nz) : data(nx, ny, nz) { }
        virtual EdgeOrientedTriangleSetData *clone() const
        {
            return new EdgeOrientedTriangleSetData(*this);
        }
        void addTriangle(plint iX, plint iY, plint iZ, int iEdge, OnEdgeTriangle const &triangle)
        {
            switch (iEdge) {
            case 0:
                data.get(iX, iY, iZ).edge1triangles.push_back(triangle);
                break;
            case 1:
                data.get(iX, iY, iZ).edge2triangles.push_back(triangle);
                break;
            case 2:
                data.get(iX, iY, iZ).edge3triangles.push_back(triangle);
                break;
            default:
                PLB_ASSERT(false);
            }
        }
        std::vector<OnEdgeTriangle> const &getTriangles(
            plint iX, plint iY, plint iZ, int iEdge) const
        {
            switch (iEdge) {
            case 0:
                return data.get(iX, iY, iZ).edge1triangles;
            case 1:
                return data.get(iX, iY, iZ).edge2triangles;
            case 2:
                return data.get(iX, iY, iZ).edge3triangles;
            default:
                PLB_ASSERT(false);
            }
        }
        void setVertex(plint iX, plint iY, plint iZ, int iEdge, Array<T, 3> const &vertex)
        {
            switch (iEdge) {
            case 0:
                data.get(iX, iY, iZ).edge1Vertex = vertex;
                data.get(iX, iY, iZ).edge1VertexDefined = true;
                break;
            case 1:
                data.get(iX, iY, iZ).edge2Vertex = vertex;
                data.get(iX, iY, iZ).edge2VertexDefined = true;
                break;
            case 2:
                data.get(iX, iY, iZ).edge3Vertex = vertex;
                data.get(iX, iY, iZ).edge3VertexDefined = true;
                break;
            default:
                PLB_ASSERT(false);
            }
        }
        bool getVertex(plint iX, plint iY, plint iZ, int iEdge, Array<T, 3> &vertex) const
        {
            switch (iEdge) {
            case 0:
                vertex = data.get(iX, iY, iZ).edge1Vertex;
                return data.get(iX, iY, iZ).edge1VertexDefined;
            case 1:
                vertex = data.get(iX, iY, iZ).edge2Vertex;
                return data.get(iX, iY, iZ).edge2VertexDefined;
            case 2:
                vertex = data.get(iX, iY, iZ).edge3Vertex;
                return data.get(iX, iY, iZ).edge3VertexDefined;
            default:
                PLB_ASSERT(false);
            }
        }
        std::vector<Array<Array<T, 3>, 3> > reconstructTriangles(
            plint iX, plint iY, plint iZ, plint iEdge) const
        {
            std::vector<Array<Array<T, 3>, 3> > triangles;
            EdgeData &localData = data.get(iX, iY, iZ);
            switch (iEdge) {
            case 0:
                if (!localData.edge1VertexDefined)
                    break;
                for (pluint i = 0; i < localData.edge1triangles; ++i) {
                    Array<T, 3> vertex1 = localData.edge1Vertex;
                    Array<plint, 4> v2info = localData.edge1triangles[i].vertex1;
                    Array<T, 3> vertex2 =
                        data.get(v2info[0], v2info[1], v2info[2]).getVertex(v2info[3]);
                    Array<plint, 4> v3info = localData.edge1triangles[i].vertex2;
                    Array<T, 3> vertex3 =
                        data.get(v3info[0], v3info[1], v3info[2]).getVertex(v3info[3]);
                    triangles.push_back(Array<Array<T, 3>, 3>(vertex1, vertex2, vertex3));
                }
                break;
            case 1:
                if (!localData.edge2VertexDefined)
                    break;
                for (pluint i = 0; i < localData.edge2triangles; ++i) {
                    Array<T, 3> vertex1 = localData.edge1Vertex;
                    Array<plint, 4> v2info = localData.edge2triangles[i].vertex1;
                    Array<T, 3> vertex2 =
                        data.get(v2info[0], v2info[1], v2info[2]).getVertex(v2info[3]);
                    Array<plint, 4> v3info = localData.edge2triangles[i].vertex2;
                    Array<T, 3> vertex3 =
                        data.get(v3info[0], v3info[1], v3info[2]).getVertex(v3info[3]);
                    triangles.push_back(Array<Array<T, 3>, 3>(vertex1, vertex2, vertex3));
                }
                break;
            case 2:
                if (!localData.edge3VertexDefined)
                    break;
                for (pluint i = 0; i < localData.edge3triangles; ++i) {
                    Array<T, 3> vertex1 = localData.edge1Vertex;
                    Array<plint, 4> v2info = localData.edge3triangles[i].vertex1;
                    Array<T, 3> vertex2 =
                        data.get(v2info[0], v2info[1], v2info[2]).getVertex(v2info[3]);
                    Array<plint, 4> v3info = localData.edge3triangles[i].vertex2;
                    Array<T, 3> vertex3 =
                        data.get(v3info[0], v3info[1], v3info[2]).getVertex(v3info[3]);
                    triangles.push_back(Array<Array<T, 3>, 3>(vertex1, vertex2, vertex3));
                }
                break;
            }
            return triangles;
        }
        T getVertexArea(plint iX, plint iY, plint iZ, int iEdge) const
        {
            EdgeData &localData = data.get(iX, iY, iZ);
            switch (iEdge) {
            case 0:
                if (!localData.edge1VertexDefined) {
                    return -1.0;
                } else {
                    std::vector<Array<Array<T, 3>, 3> > triangles =
                        reconstructTriangles(iX, iY, iZ, iEdge);
                    T area = T();
                    for (pluint i = 0; i < triangles.size(); ++i) {
                        Array<Array<T, 3>, 3> const &triangle = triangles[i];
                        T nextArea = computeTriangleArea(triangle[0], triangle[1], triangle[2]);
                        area += nextArea;
                    }
                    return area / (T)3.0;
                }
            case 1:
                if (!localData.edge2VertexDefined) {
                    return -1.0;
                } else {
                    std::vector<Array<Array<T, 3>, 3> > triangles =
                        reconstructTriangles(iX, iY, iZ, iEdge);
                    T area = T();
                    for (pluint i = 0; i < triangles.size(); ++i) {
                        Array<Array<T, 3>, 3> const &triangle = triangles[i];
                        T nextArea = computeTriangleArea(triangle[0], triangle[1], triangle[2]);
                        area += nextArea;
                    }
                    return area / (T)3.0;
                }
            case 2:
                if (!localData.edge3VertexDefined) {
                    return -1.0;
                } else {
                    std::vector<Array<Array<T, 3>, 3> > triangles =
                        reconstructTriangles(iX, iY, iZ, iEdge);
                    T area = T();
                    for (pluint i = 0; i < triangles.size(); ++i) {
                        Array<Array<T, 3>, 3> const &triangle = triangles[i];
                        T nextArea = computeTriangleArea(triangle[0], triangle[1], triangle[2]);
                        area += nextArea;
                    }
                    return area / (T)3.0;
                }
            default:
                PLB_ASSERT(false);
            }
        }
        std::vector<T> const &getScalars(plint iX, plint iY, plint iZ, int iEdge) const
        {
            switch (iEdge) {
            case 0:
                return data.get(iX, iY, iZ).scalars1;
            case 1:
                return data.get(iX, iY, iZ).scalars2;
            case 2:
                return data.get(iX, iY, iZ).scalars3;
            default:
                PLB_ASSERT(false);
            }
        }
        std::vector<T> &getScalars(plint iX, plint iY, plint iZ, int iEdge)
        {
            switch (iEdge) {
            case 0:
                return data.get(iX, iY, iZ).scalars1;
            case 1:
                return data.get(iX, iY, iZ).scalars2;
            case 2:
                return data.get(iX, iY, iZ).scalars3;
            default:
                PLB_ASSERT(false);
            }
        }
        bool isEdgeVertexDefined(plint iX, plint iY, plint iZ, int iEdge)
        {
            switch (iEdge) {
            case 0:
                return data.get(iX, iY, iZ).edge1VertexDefined;
            case 1:
                return data.get(iX, iY, iZ).edge2VertexDefined;
            case 2:
                return data.get(iX, iY, iZ).edge3VertexDefined;
            default:
                PLB_ASSERT(false);
            }
        }
        plint getNx() const
        {
            return data.getNx();
        }
        plint getNy() const
        {
            return data.getNy();
        }
        plint getNz() const
        {
            return data.getNz();
        }
        Box3D getBoundingBox() const
        {
            return data.getBoundingBox();
        }

    private:
        struct EdgeData {
            EdgeData() :
                edge1Vertex(Array<T, 3>::zero()),
                edge2Vertex(Array<T, 3>::zero()),
                edge3Vertex(Array<T, 3>::zero()),
                edge1VertexDefined(false),
                edge2VertexDefined(false),
                edge3VertexDefined(false)
            { }
            Array<T, 3> const &getVertex(plint iEdge) const
            {
                switch (iEdge) {
                case 0:
                    return edge1Vertex;
                case 1:
                    return edge2Vertex;
                case 2:
                    return edge3Vertex;
                }
            }
            // Each element in the three following vectors defines the topology of
            // a triangle.
            std::vector<OnEdgeTriangle> edge1triangles;
            std::vector<OnEdgeTriangle> edge2triangles;
            std::vector<OnEdgeTriangle> edge3triangles;
            // Every edge has at most one vertex which is shared by all
            // triangles that have a vertex on this edge.
            Array<T, 3> edge1Vertex, edge2Vertex, edge3Vertex;
            std::vector<T> scalars1, scalars2, scalars3;
            bool edge1VertexDefined, edge2VertexDefined, edge3VertexDefined;
        };

    private:
        ScalarField3D<EdgeData> data;
    };

private:
    void marchingCubeImpl(
        plint iX, plint iY, plint iZ, plint surfaceId, std::vector<Triangle> &triangles,
        int &cubeIndex, std::vector<Array<T, 3> > &vertlist);
    void polygonize(
        plint iX, plint iY, plint iZ, plint surfaceId, std::vector<Triangle> &triangles);
    /// Edge attribution contains three integers to label the cell ID,
    /// and one integer to label one of the three edges assigned to this cell.
    void polygonize(
        plint iX, plint iY, plint iZ, plint surfaceId, std::vector<Triangle> &triangles,
        std::vector<Array<plint, 4> > &edgeAttributions);
    static void removeFromVertex(
        Array<T, 3> const &p0, Array<T, 3> const &p1, Array<T, 3> &intersection);

private:
    std::vector<plint> surfaceIds;
    IsoSurfaceDefinition3D<T> *isoSurface;
    bool edgeOrientedData;
    plint edgeOrientedEnvelope;
};

struct MarchingCubeConstants {
    static const int edgeTable[256];
    static const int triTable[256][16];
    static const int edgeNeighb[12][3];
    static const int edgeOrient[12];
};

}  // namespace plb

#endif  // MARCHING_CUBE_H
